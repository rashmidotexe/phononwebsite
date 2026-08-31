#!/usr/bin/env python3
"""
generate_zumba_jsons.py
------------------------
Build the 27 Porto-notation zumbadb/ZMO_pol_<pol>_ei_<ei>_es_<es>.json files
for the Zn2Mo3O8 (ZMO) Raman visualizer from the underlying QE dynamical
matrix (zmo.dyn, carrying the DFPT Raman tensor) and the four LO-TO-split
eigenvector files (zmo_0.eig, zmo_x.eig, zmo_y.eig, zmo_z.eig).

Two things this fixes relative to the first version of this script:

1. Symmetry labels (mode_symmetry[]) no longer mislabel Raman-silent
   ("dark") modes. See classify_modes()/mode_symmetry.py for why the old
   fixed absolute tolerance let numerical noise through, and why a
   doubly-degenerate check on E_1/E_2 is needed on top of that.
2. The "vectors" field (used for the 3D animation) is now taken from the
   SAME per-polarization eigenvector file used to compute that geometry's
   frequencies/Raman activities/symmetries, instead of being left as a
   fixed baseline that doesn't match most of the 27 geometries.

Usage:
    python3 generate_zumba_jsons.py --base zmo.json --dyn zmo.dyn \
        --modes-0 zmo_0.eig --modes-x zmo_x.eig --modes-y zmo_y.eig --modes-z zmo_z.eig \
        --out-dir zumbadb
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import re
import sys

import numpy as np
from scipy.constants import c, h, k

from phononweb.scripts.zumba_symmetry import classify_modes

AXES = ["x", "y", "z"]

AXIS_VEC = {
    "x": np.array([1.0, 0.0, 0.0]),
    "y": np.array([0.0, 1.0, 0.0]),
    "z": np.array([0.0, 0.0, 1.0]),
}

_SPECIES_LINE_RE = re.compile(r"^\s*(\d+)\s+'[^']*'\s+([\d.eEdD+-]+)")
_ATOM_POL_RE = re.compile(r"atom\s*#\s*(\d+)\s+pol\.\s+(\d+)")
_FREQ_LINE_RE = re.compile(r"\s*freq\s*\(\s*\d+\s*\)\s*=.*\[cm-1\]")
_FREQ_VALUE_RE = re.compile(r"=\s*([\-\d.]+)\s*\[cm-1\]")


def parse_dyn_for_raman(fname):
    """Parse a QE .dyn file's per-type masses and DFPT Raman tensor block."""
    with open(fname) as fh:
        lines = fh.readlines()

    header = lines[2].split()
    ntyp, nat = int(header[0]), int(header[1])

    type_mass = {}
    for line in lines[3:3 + ntyp]:
        m = _SPECIES_LINE_RE.match(line)
        if m:
            type_mass[int(m.group(1))] = float(m.group(2))

    masses = np.zeros(nat)
    for line in lines[3 + ntyp: 3 + ntyp + nat]:
        parts = line.split()
        iatom, itype = int(parts[0]) - 1, int(parts[1])
        masses[iatom] = type_mass[itype]

    raman_tens = np.zeros((nat, 3, 3, 3))
    raman_start = None
    for i, line in enumerate(lines):
        if "Raman tensor" in line:
            raman_start = i + 1
            break

    if raman_start is not None:
        i = raman_start
        while i < len(lines):
            m = _ATOM_POL_RE.match(lines[i].strip())
            if m:
                iatom, ipol = int(m.group(1)) - 1, int(m.group(2)) - 1
                mat = []
                for _ in range(3):
                    i += 1
                    mat.append([float(v) for v in lines[i].split()])
                raman_tens[iatom, ipol] = np.array(mat)
            i += 1

    return masses, raman_tens


def parse_modes_for_eigvecs(fname, nat):
    """Parse a matdyn/dynmat-style .eig file into (freqs, eigenvectors)."""
    with open(fname) as fh:
        lines = fh.readlines()

    freqs, eigvecs = [], []
    i = 0
    while i < len(lines):
        if _FREQ_LINE_RE.match(lines[i]):
            freqs.append(float(_FREQ_VALUE_RE.search(lines[i]).group(1)))
            evec = np.zeros((nat, 3), dtype=complex)
            for s in range(nat):
                i += 1
                nums = re.findall(r"[\-\d.Ee+\-]+", lines[i])
                evec[s, 0] = complex(float(nums[0]), float(nums[1]))
                evec[s, 1] = complex(float(nums[2]), float(nums[3]))
                evec[s, 2] = complex(float(nums[4]), float(nums[5]))
            eigvecs.append(evec)
        i += 1

    return np.array(freqs), np.array(eigvecs)


def stokes_factor(freqs_cm1, temperature_k):
    f = np.asarray(freqs_cm1, dtype=float) * 1e2 * c
    factor = np.zeros_like(f)
    mask = freqs_cm1 > 1.0
    exponent = h * f[mask] / (k * temperature_k)
    n_v = 1.0 / np.expm1(exponent)
    factor[mask] = (n_v + 1.0) / f[mask]
    return factor


def compute_intensities(ei, es, R_modes, freqs, temperature_k):
    e_i, e_s = AXIS_VEC[ei], AXIS_VEC[es]
    intensities = np.array([abs(e_s @ R_modes[n] @ e_i) ** 2 for n in range(len(freqs))])
    intensities *= stokes_factor(freqs, temperature_k)
    return intensities


def eigvecs_to_json(eigvecs):
    """(nmodes, natoms, 3) complex -> nested [mode][atom][axis][re, im] lists."""
    return [
        [[[float(c.real), float(c.imag)] for c in atom] for atom in mode]
        for mode in eigvecs
    ]


def generate(base_json, dyn_file, modes_files, out_dir, temperature_k=300.0):
    with open(base_json) as f:
        base_data = json.load(f)

    gamma_index = None
    for i, qpt in enumerate(base_data.get("qpoints", [])):
        if all(abs(x) < 1e-5 for x in qpt):
            gamma_index = i
            break
    if gamma_index is None:
        raise ValueError("Gamma point not found in base JSON qpoints")

    masses, raman_tens = parse_dyn_for_raman(dyn_file)
    nat = len(masses)
    sqrt_m = np.sqrt(masses)

    os.makedirs(out_dir, exist_ok=True)
    written = []

    for pol in AXES:
        eig_key = pol if pol in modes_files else "0"
        freqs, eigvecs = parse_modes_for_eigvecs(modes_files[eig_key], nat)

        n_base = len(base_data["eigenvalues"][gamma_index])
        if len(freqs) != n_base:
            print(
                f"  WARNING: {modes_files[eig_key]} has {len(freqs)} modes, "
                f"base JSON has {n_base}",
                file=sys.stderr,
            )

        mw_evec = eigvecs / sqrt_m[np.newaxis, :, np.newaxis]
        R_modes = np.einsum("sabg,nsa->nbg", raman_tens, mw_evec)
        sym_labels = classify_modes(R_modes, freqs)
        vectors_json = eigvecs_to_json(eigvecs)

        for ei in AXES:
            for es in AXES:
                label = f"{pol}({ei},{es}){pol}"
                out_name = f"ZMO_pol_{pol}_ei_{ei}_es_{es}.json"
                out_path = os.path.join(out_dir, out_name)

                intensities = compute_intensities(ei, es, R_modes, freqs, temperature_k)

                data = copy.deepcopy(base_data)
                data["eigenvalues"][gamma_index] = freqs.tolist()
                data["vectors"][gamma_index] = vectors_json
                data["raman_intensities"] = intensities.tolist()
                data["gamma_index"] = gamma_index
                data["mode_symmetry"] = sym_labels
                data["porto_label"] = label
                data["porto_pol"] = pol
                data["porto_ei"] = ei
                data["porto_es"] = es

                with open(out_path, "w") as f:
                    json.dump(data, f, indent=1)

                written.append(out_name)

    return written


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, help="base site json, e.g. zmo.json")
    parser.add_argument("--dyn", required=True, help="QE .dyn file with the Raman tensor")
    parser.add_argument("--modes-0", required=True)
    parser.add_argument("--modes-x", required=True)
    parser.add_argument("--modes-y", required=True)
    parser.add_argument("--modes-z", required=True)
    parser.add_argument("--out-dir", default="zumbadb")
    parser.add_argument("--temperature-k", type=float, default=300.0)
    args = parser.parse_args()

    modes_files = {"0": args.modes_0, "x": args.modes_x, "y": args.modes_y, "z": args.modes_z}
    written = generate(args.base, args.dyn, modes_files, args.out_dir, args.temperature_k)
    print(f"wrote {len(written)} files to {args.out_dir}/")


if __name__ == "__main__":
    main()

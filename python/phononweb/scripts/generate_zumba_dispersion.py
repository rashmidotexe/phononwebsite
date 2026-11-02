#!/usr/bin/env python3
"""
generate_zumba_dispersion.py
----------------------------
Build zumbadb/ZMO_dispersion.json, the full phonon dispersion of Zn2Mo3O8
(ZMO) shown in the Dispersion tab, from QE matdyn.x output along a
high-symmetry path (matdyn.modes + the matdyn.in that defined the path).

The file is shared by all 27 scattering-geometry JSONs: the page prepends
each geometry's own Gamma point (with its LO-TO frequencies, Raman
intensities and symmetry labels) to this path when it loads a geometry, so
the 27 small Gamma-only files stay as they are.

matdyn.modes holds normalised displacements (e / sqrt(m), renormalised), the
same convention as the Gamma-point vectors of the geometry JSONs, so they are
used as they are. Every --stride-th q-point is kept to keep the file small.

The force constants carry Z* and eps_inf (ph.x computes them at Gamma in a
ldisp run), so matdyn adds the LO-TO term, which depends on the direction q
approaches Gamma. At a Gamma inside the path (K-Gamma-A) matdyn gives one set
of frequencies, the in-plane limit, so the line is broken just after it rather
than drawn up to the out-of-plane LO modes of the Gamma-A segment.

Usage:
    python3 generate_zumba_dispersion.py \
        --modes matdyn.modes --matdyn-in matdyn.in \
        --base zumbadb/ZMO_pol_x_ei_x_es_x.json --out zumbadb/ZMO_dispersion.json
"""

from __future__ import annotations

import argparse
import base64
import json
import re

import numpy as np

_Q_RE = re.compile(r"^\s*q\s*=\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)")
_FREQ_RE = re.compile(r"freq\s*\(\s*\d+\s*\)\s*=.*=\s*([-\d.Ee+]+)\s*\[cm-1\]")
_NUM_RE = re.compile(r"[-+]?\d+\.\d*(?:[Ee][-+]?\d+)?")


def read_matdyn_modes(fname, nat):
    """matdyn.modes -> (qpoints_cart [nq,3] in 2pi/alat, freqs [nq,nmodes], disp [nq,nmodes,nat,3])."""
    with open(fname) as fh:
        lines = fh.readlines()

    qpts, freqs, disps = [], [], []
    i = 0
    while i < len(lines):
        m = _Q_RE.match(lines[i])
        if m:
            qpts.append([float(x) for x in m.groups()])
            freqs.append([])
            disps.append([])
        else:
            m = _FREQ_RE.search(lines[i])
            if m:
                freqs[-1].append(float(m.group(1)))
                vec = np.zeros((nat, 3), dtype=complex)
                for s in range(nat):
                    i += 1
                    v = [float(x) for x in _NUM_RE.findall(lines[i])]
                    vec[s] = [complex(v[0], v[1]), complex(v[2], v[3]), complex(v[4], v[5])]
                disps[-1].append(vec)
        i += 1

    return np.array(qpts), np.array(freqs), np.array(disps)


def read_matdyn_path(fname):
    """matdyn.in -> [(label, npoints), ...] of the band path."""
    with open(fname) as fh:
        text = fh.read()

    tail = text.split("/", 1)[1].split("\n")
    tail = [l for l in tail if l.strip()]
    npath = int(tail[0].split()[0])
    path = []
    for line in tail[1:1 + npath]:
        body, _, label = line.partition("!")
        path.append((label.strip(), int(body.split()[3])))
    return path


def encode_q11_int16(vectors):
    """[nq,nmodes,nat,3] complex -> the page's "vectors_compressed" q11-int16-base64 payload."""
    arr = np.stack([vectors.real, vectors.imag], axis=-1).astype(np.float32)
    quantized = np.clip(np.rint(arr * 2047.0), -2047, 2047).astype("<i2")
    return {
        "format": "q11-int16-base64",
        "scale": 2047,
        "shape": list(quantized.shape),
        "data": base64.b64encode(quantized.tobytes(order="C")).decode("ascii"),
    }


def build(modes_file, matdyn_in, base_json, stride=2):
    with open(base_json) as fh:
        base = json.load(fh)
    lattice = np.array(base["lattice"])
    nat = base["natoms"]

    path = read_matdyn_path(matdyn_in)

    qcart, freqs, disp = read_matdyn_modes(modes_file, nat)
    nq = len(qcart)

    # high-symmetry indices along the full path (matdyn band form: n points per segment)
    hs = [0]
    for _, n in path[:-1]:
        hs.append(hs[-1] + n)
    if hs[-1] != nq - 1:
        raise ValueError(f"{matdyn_in} path has {hs[-1] + 1} points, {modes_file} has {nq}")
    if any(n % stride for _, n in path[:-1]):
        raise ValueError(f"stride {stride} does not divide every path segment")

    keep = np.arange(0, nq, stride)
    labels = [label for label, _ in path]

    # cartesian 2pi/alat -> reduced: q_red_i = q_cart . a_i / alat
    alat = np.linalg.norm(lattice[0])
    qred = qcart @ lattice.T / alat

    # path length in 1/Angstrom (2pi/alat units scaled by 2pi/alat)
    step = np.linalg.norm(np.diff(qcart, axis=0), axis=1) * 2 * np.pi / alat
    dist = np.concatenate([[0.0], np.cumsum(step)])


    gamma_inside = [i // stride for i, l in zip(hs[1:-1], labels[1:-1]) if l.upper() == "GAMMA"]
    starts = [0] + [g + 1 for g in gamma_inside]
    line_breaks = [[a, b] for a, b in zip(starts, starts[1:] + [len(keep)])]

    return {
        "name": base["name"],
        "source": "QE matdyn.x, LDA, 3x3x2 q-grid (05_phonon_raman/02_ph/dispersion), LO-TO from Z* and eps_inf",
        "qpoints": (np.round(qred[keep], 8) + 0.0).tolist(),  # + 0.0: no -0.0
        "distances": np.round(dist[keep], 6).tolist(),
        "highsym_qpts": [[int(i // stride), "GAMMA" if l.upper() == "GAMMA" else l] for i, l in zip(hs, labels)],
        "eigenvalues": np.round(freqs[keep], 4).tolist(),
        "line_breaks": line_breaks,
        "vectors_compressed": encode_q11_int16(disp[keep]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--modes", required=True, help="matdyn.modes along the band path")
    parser.add_argument("--matdyn-in", required=True, help="matdyn.in that defined the path")
    parser.add_argument("--base", required=True, help="a ZMO geometry JSON (lattice, atom order)")
    parser.add_argument("--out", default="zumbadb/ZMO_dispersion.json")
    parser.add_argument("--stride", type=int, default=2, help="keep every n-th q-point")
    args = parser.parse_args()

    data = build(args.modes, args.matdyn_in, args.base, args.stride)
    with open(args.out, "w") as fh:
        json.dump(data, fh, separators=(",", ":"))
    print(f"wrote {args.out}: {len(data['qpoints'])} q-points, {len(data['eigenvalues'][0])} modes")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
generate_zumba_resonance.py
---------------------------
Build zumbadb/resonance/ZMO_resonance.json for the "Resonance" view of the
Zn2Mo3O8 (ZMO) Raman visualizer from an exciton-phonon resonant Raman run
(LetzElPhC + Yambo GW-BSE + PhdScripts raman_driver).

Inputs:
    --tensor     Raman_tensor_e0_iR.npy  [N_wL, n_modes, 3, 3] complex, raw
                 PhdScripts units, computed on np.linspace(*--el-grid)
    --freqs      ph_freq_e0.npy          [n_modes] phonon energies (Hartree)
    --excitons   BS_eigs.npy             [n_exc] exciton energies (Hartree)
    --dipoles    ex_dip_iR.npy           [3, n_exc] exciton dipoles (same gauge)
    --experiment exp446_perp.dat         optional measured resonance profile (E_eV, I)

Intensity convention (same as ref/analysis/plot_excph_raman.py):
    alpha = R * sqrt(2 w_nu) * sqrt(w_L / |w_L - w_nu|)
        (undoes the 1/sqrt(2 w_nu) and sqrt(|w_L - w_nu| / w_L) that
        raman.compute_Raman_oneph_exc folds into R, i.e. raman_qe_units'
        length-gauge conversion up to a constant)
    I = |e_i . alpha . e_s|^2 * (n(w_nu) + 1) / w_nu        (Stokes, 300 K;
        alpha[..., incoming, outgoing], see geometry_intensity)
summed over the partners of degenerate E modes. Overall constants cancel
because the page only shows normalised intensities and ratios.

Usage:
    python3 generate_zumba_resonance.py \
        --tensor Raman_tensor_e0_iR.npy --freqs ph_freq_e0.npy \
        --excitons BS_eigs.npy --dipoles ex_dip_iR.npy \
        --experiment exp446_perp.dat --out zumbadb/resonance/ZMO_resonance.json

The enhancement factor reference ("offres") is the same calculation at
1.0 eV, the off-resonant energy ref/08_excph_raman/05_raman_gw_selectivity
uses; each geometry's "ref_peak" (the internal standard for selectivity
ratios) is the 442.3 cm-1 A_1 line where it is allowed, as in that folder's
"rel. to 442" tables, and otherwise the strongest peak at 1.0 eV.
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from phononweb.scripts.generate_zumba_jsons import AXIS_VEC, stokes_factor

HA2EV = 27.211386245988
HA2CM = 219474.6313705
T_KELVIN = 300.0
DEGENERACY_CM = 0.5
SILENT_FRACTION = 1e-4

# Backscattering geometries offered on the page, as (propagation, e_i, e_s),
# with the C6v species each one allows (same rules as PORTO_RULES in
# src/zumbamain.js).
GEOMETRIES = [
    ("x", "y", "y", ["A_1", "E_2"]),
    ("x", "y", "z", ["E_1"]),
    ("x", "z", "z", ["A_1"]),
    ("y", "x", "x", ["A_1", "E_2"]),
    ("y", "x", "z", ["E_1"]),
    ("y", "z", "z", ["A_1"]),
    ("z", "x", "y", ["E_2"]),
]


def to_polarizability(raman_tensor, energies_ev, freqs_ha):
    """Raw PhdScripts tensor -> polarizability derivative (arbitrary units)."""
    w_l = np.asarray(energies_ev, dtype=float)[:, None] / HA2EV
    w = np.asarray(freqs_ha, dtype=float)[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.sqrt(2.0 * w) * np.sqrt(w_l / np.abs(w_l - w))
    return raman_tensor * scale[:, :, None, None]


def group_peaks(alpha, freqs_cm, stokes):
    """
    Group modes into Raman peaks (degenerate partners together) and label
    them, the same way as classify() in ref/analysis/plot_excph_raman.py:
    singlets are A_1, pairs are E_1 when the xz/yz components carry most of
    the weight and E_2 otherwise. Acoustic modes and peaks weaker than
    SILENT_FRACTION of the strongest one (B_1/B_2, silent) are dropped.
    """
    intensity = np.abs(alpha) ** 2 * stokes[None, :, None, None]
    peak = intensity.sum(axis=(2, 3)).max(axis=0)
    optical = range(3, len(freqs_cm))
    top = max(peak[m] for m in optical)

    used, peaks = set(), []
    for m in optical:
        if m in used or peak[m] < SILENT_FRACTION * top:
            continue
        group = [n for n in optical
                 if n not in used and abs(freqs_cm[n] - freqs_cm[m]) < DEGENERACY_CM]
        used.update(group)
        if len(group) > 2:
            # a near-degenerate singlet swept into an E pair is only harmless
            # if it is negligible next to the pair it is summed with
            strength = np.sort(peak[group])[::-1]
            if strength[2:].max() >= 1e-3 * strength[0]:
                raise ValueError(f"modes {group} at {freqs_cm[m]:.1f} cm-1 are not one degenerate peak")
        P = intensity[:, group].max(axis=0).sum(axis=0)
        e1_weight = P[0, 2] + P[2, 0] + P[1, 2] + P[2, 1]
        if len(group) == 1:
            sym = "A_1"
        else:
            sym = "E_1" if e1_weight > 0.5 * P.sum() else "E_2"
        peaks.append({"freq": float(np.mean(freqs_cm[group])), "sym": sym, "modes": group})
    return peaks


def geometry_intensity(alpha, stokes, modes, ei, es):
    """
    I(w_L) of one peak in one geometry, summed over degenerate partners.
    The resonant tensor is not symmetric: in raman.compute_Raman_oneph_exc
    the first light index carries the incoming 1/(w_L - E_S) denominator and
    the second the outgoing 1/(w_L - E_S - w_nu) one, so the amplitude is
    e_i . alpha . e_s.
    """
    e_i, e_s = AXIS_VEC[ei], AXIS_VEC[es]
    return sum(np.abs(np.einsum("a,wab,b->w", e_i, alpha[:, n], e_s)) ** 2 * stokes[n] for n in modes)


def absorption(exc_ev, osc, energies_ev, fwhm):
    """Lorentzian-broadened exciton absorption along x, y, z (arbitrary units)."""
    hw = fwhm / 2.0
    lor = hw / np.pi / ((energies_ev[:, None] - exc_ev[None, :]) ** 2 + hw ** 2)
    return lor @ osc.T                                   # [n_E, 3]


def merge_excitons(exc_ev, osc, tol=1e-3):
    """
    Merge degenerate excitons (all within tol eV of the group's first
    member) into one stick: oscillator strengths add, the energy is their
    strength-weighted mean. Measuring from the first member keeps the
    sub-meV quasi-continuum above ~2.8 eV from chaining into wide blocks.
    """
    order = np.argsort(exc_ev)
    groups, current = [], [order[0]]
    for b in order[1:]:
        if exc_ev[b] - exc_ev[current[0]] > tol:
            groups.append(current)
            current = []
        current.append(b)
    groups.append(current)

    energies, strengths = [], []
    for g in groups:
        s = osc[:, g]
        w = s.sum(axis=0)
        energies.append(np.average(exc_ev[g], weights=w) if w.sum() > 0 else exc_ev[g].mean())
        strengths.append(s.sum(axis=1))
    return np.array(energies), np.array(strengths).T     # [n], [3, n]


def read_experiment(fname):
    """Two-column (E_eV, I) text file with '#' comments."""
    return np.loadtxt(fname, comments="#", ndmin=2)[:, :2]


def sig(values, digits=4):
    return [float(f"{v:.{digits}g}") for v in np.asarray(values, dtype=float)]


def porto_label(prop, ei, es):
    """Backscattering Porto notation as the page shows it, e.g. -X(YY)X."""
    return f"-{prop.upper()}({ei.upper()}{es.upper()}){prop.upper()}"


def build_from_arrays(raw, el_grid, freqs_ha, exc_ha, dipoles, experiment=None, emin=2.3, emax=3.0,
                      offres_energy=1.0, ref_freq=442.3, broadening=0.2, bright=1e-3, description="",
                      exciton_count=None):
    el_all = np.linspace(el_grid[0], el_grid[1], raw.shape[0])
    finite = np.isfinite(raw).all(axis=(1, 2, 3))       # the driver returns inf at w_L = 0
    raw, el_all = raw[finite], el_all[finite]
    freqs_ha = np.asarray(freqs_ha, dtype=float)
    freqs_cm = freqs_ha * HA2CM
    stokes = stokes_factor(freqs_cm, T_KELVIN)
    alpha_all = to_polarizability(raw.astype(np.complex128), el_all, freqs_ha)

    peaks = group_peaks(alpha_all, freqs_cm, stokes)

    window = (el_all >= emin - 1e-9) & (el_all <= emax + 1e-9)
    energies = el_all[window]
    alpha = alpha_all[window]
    # Off-resonant reference for the enhancement factor: a laser energy well
    # below the first exciton but far above every phonon energy (the w_L -> 0
    # end of the grid only holds the BSE-window part of the tensor)
    i_off = int(np.argmin(np.abs(el_all - offres_energy)))
    if abs(el_all[i_off] - offres_energy) > 1e-6:
        raise ValueError(f"no grid point at {offres_energy} eV (nearest {el_all[i_off]:.6f})")
    alpha_off = alpha_all[i_off:i_off + 1]

    geometries = []
    for prop, ei, es, allowed in GEOMETRIES:
        curves = [geometry_intensity(alpha, stokes, p["modes"], ei, es) for p in peaks]
        top = max(cv.max() for cv in curves)
        kept, rows, offres = [], [], []
        for i, (p, cv) in enumerate(zip(peaks, curves)):
            # symmetry-forbidden peaks only carry numerical leakage (<= 1 %)
            if p["sym"] not in allowed or cv.max() < SILENT_FRACTION * top:
                continue
            kept.append(i)
            rows.append(sig(cv))
            offres.append(float(geometry_intensity(alpha_off, stokes, p["modes"], ei, es)[0]))
        # internal standard for selectivity ratios: the user's reference line
        # where its species is allowed, otherwise the strongest off-resonant peak
        ref_candidates = [i for i in kept if abs(peaks[i]["freq"] - ref_freq) < 1.0]
        ref_peak = ref_candidates[0] if ref_candidates else kept[int(np.argmax(offres))]
        geometries.append({
            "id": f"{prop}({ei}{es}){prop}",
            "label": porto_label(prop, ei, es),
            "prop": prop, "ei": ei, "es": es,
            "allowed": allowed,
            "peaks": kept,
            "intensity": rows,
            "offres": sig(offres),
            "ref_peak": ref_peak,
        })

    exc_ev = np.asarray(exc_ha).real.astype(float) * HA2EV
    osc = np.abs(np.asarray(dipoles)[:3]) ** 2
    abs_e = np.round(np.arange(emin - 0.1, emax + 0.1 + 1e-9, 0.002), 4)
    eps = absorption(exc_ev, osc, abs_e, broadening)
    m_ev, m_osc = merge_excitons(exc_ev, osc)
    keep = (m_ev >= abs_e[0]) & (m_ev <= abs_e[-1]) & (m_osc.max(axis=0) >= bright * m_osc.max())

    data = {
        "material": "Zn2Mo3O8",
        "description": description,
        "temperature": T_KELVIN,
        "broadening": broadening,
        "exciton_count": exciton_count,
        "window": [emin, emax],
        "grid_step": float(el_all[1] - el_all[0]),
        "offres_energy": offres_energy,
        "energies": sig(energies, 6),
        "peaks": [{"freq": round(p["freq"], 2), "sym": p["sym"], "modes": p["modes"]} for p in peaks],
        "geometries": geometries,
        "absorption": {
            "energies": sig(abs_e, 6),
            "x": sig(eps[:, 0]), "y": sig(eps[:, 1]), "z": sig(eps[:, 2]),
        },
        "excitons": {
            "energies": sig(m_ev[keep], 6),
            "x": sig(m_osc[0, keep]), "y": sig(m_osc[1, keep]), "z": sig(m_osc[2, keep]),
        },
    }
    if experiment is not None:
        a1 = [i for i, p in enumerate(peaks) if p["sym"] == "A_1"]
        data["experiment"] = {
            "label": "Gudima et al., Fig. 4b — 446 cm⁻¹, ⊥c (-X(YY)X)",
            "geometries": ["x(yy)x", "y(xx)y"],
            # the measured 446 cm-1 line is the calculated 442.3 A_1 peak
            # (same mapping as ref/08_excph_raman/05_raman_gw_selectivity/vs_exp.py)
            "calc_peak": min(a1, key=lambda i: abs(peaks[i]["freq"] - 442.3)),
            "err_E": 0.03, "err_I": 0.15,
            "note": "hand-digitised, CaF₂-normalised",
            "points": [[float(e), float(v)] for e, v in experiment],
        }
    return data


def build(args):
    return build_from_arrays(
        np.load(args.tensor, allow_pickle=False),
        args.el_grid,
        np.load(args.freqs, allow_pickle=False),
        np.load(args.excitons, allow_pickle=False),
        np.load(args.dipoles, allow_pickle=False),
        experiment=read_experiment(args.experiment) if args.experiment else None,
        emin=args.emin, emax=args.emax, offres_energy=args.offres, ref_freq=args.ref_peak,
        broadening=args.broadening, bright=args.bright,
        description=args.description, exciton_count=args.exciton_count,
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tensor", required=True)
    ap.add_argument("--freqs", required=True)
    ap.add_argument("--excitons", required=True)
    ap.add_argument("--dipoles", required=True)
    ap.add_argument("--el-grid", type=float, nargs=2, default=[0.0, 4.0],
                    help="laser-energy range (eV) the tensor was computed on")
    ap.add_argument("--emin", type=float, default=2.3)
    ap.add_argument("--emax", type=float, default=3.0)
    ap.add_argument("--broadening", type=float, default=0.2,
                    help="exciton FWHM (eV) of the Raman run, reused for the absorption")
    ap.add_argument("--bright", type=float, default=1e-3,
                    help="keep excitons with oscillator strength >= this fraction of the brightest")
    ap.add_argument("--offres", type=float, default=1.0,
                    help="off-resonant laser energy (eV) the enhancement factor is taken against")
    ap.add_argument("--ref-peak", type=float, default=442.3,
                    help="peak (cm-1) selectivity ratios are taken against wherever it is allowed")
    ap.add_argument("--experiment", default=None,
                    help="optional measured resonance profile (E_eV, I), e.g. exp446_perp.dat")
    ap.add_argument("--exciton-count", type=int, default=17280,
                    help="number of excitons in the Raman run (metadata only)")
    ap.add_argument("--description",
                    default="LDA, GW-BSE, DIP_iR (length gauge), 17280 excitons, exciton FWHM 0.2 eV, 300 K Stokes")
    ap.add_argument("--out", default="zumbadb/resonance/ZMO_resonance.json")
    args = ap.parse_args(argv)

    data = build(args)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(data, fh, separators=(",", ":"))
        fh.write("\n")
    print(f"wrote {args.out}: {len(data['energies'])} laser energies, {len(data['peaks'])} peaks, "
          f"{len(data['geometries'])} geometries")


if __name__ == "__main__":
    main()

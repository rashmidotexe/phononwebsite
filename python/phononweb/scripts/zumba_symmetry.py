#!/usr/bin/env python3
"""
Assign a C6v Raman symmetry label (A_1, E_1, E_2, or "dark"/silent) to each
phonon mode from its Raman tensor, projected onto the mode's eigendisplacement.

Two corrections against numerical noise on top of the tensor-shape heuristic
that decides A_1 vs E_1 vs E_2 -- see classify_modes() docstring for why
both are needed; a fixed absolute "is this dark" tolerance let modes with
real tensor magnitude ~1e-7 (7-8 orders of magnitude below genuinely active
modes, ~1e-2 to 1) get treated as significant, because the noise-floor scale
is calculation-specific, not a universal constant.
"""

from __future__ import annotations

import numpy as np


def classify_tensor_shape(tensor_real, relative_off_diagonal_floor=0.15):
    """
    Classify a single mode's (already magnitude-gated) Raman tensor by its
    shape, per the standard C6v Raman tensor forms:
      A_1(z):  diag(a, a, b)
      E_1(x,y): dominant xz/zx or yz/zy off-diagonals
      E_2:      dominant xy/yx off-diagonal, or (xx, yy) with opposite sign
    """
    max_val = np.max(np.abs(tensor_real))
    normalized = tensor_real / max_val
    magnitude = np.abs(normalized)

    e1_score = magnitude[0, 2] + magnitude[2, 0] + magnitude[1, 2] + magnitude[2, 1]
    if e1_score > relative_off_diagonal_floor:
        return "E_1"

    e2_off_score = magnitude[0, 1] + magnitude[1, 0]
    if e2_off_score > relative_off_diagonal_floor:
        return "E_2"

    if normalized[0, 0] * normalized[1, 1] < -relative_off_diagonal_floor:
        return "E_2"
    return "A_1"


def classify_modes(
    R_modes,
    freqs,
    acoustic_freq_tol_cm1=3.0,
    magnitude_rel_floor=0.01,
):
    """
    Assign a symmetry label to every mode, or "dark" if it isn't Raman
    active. R_modes[n] is the mode's Raman tensor (3x3, possibly complex --
    only the real part is physical for a q=0 mode without an imaginary
    dielectric response). freqs[n] is that mode's frequency in cm^-1.

    Two corrections against numerical noise, applied in order:

    1. Acoustic modes (|freq| <= acoustic_freq_tol_cm1) are always dark.
       Pure translations cannot change the polarizability tensor; any
       nonzero projection here is an ASR/numerical residual, not signal,
       regardless of how large it looks in absolute terms.
    2. For optical modes, "dark" is decided by comparing each mode's tensor
       magnitude (max |entry|) to the LARGEST magnitude among optical modes
       in *this* calculation (a relative threshold), not a fixed absolute
       tolerance -- the real noise floor depends on how tightly this
       specific DFPT run converged. This is what catches modes with real
       tensor magnitude ~1e-7 to 1e-8 (7-8 orders of magnitude below
       genuinely active modes) that a fixed tolerance like 1e-8 lets
       through, since it can't tell "small but real" from "pure noise" in
       absolute terms.

    An earlier version of this function also demoted E_1/E_2 modes lacking
    an exactly frequency-degenerate partner to "dark", on the assumption
    that E-symmetry modes must be doubly degenerate. That assumption does
    NOT hold here: validated against independently-computed reference
    spectra (Lorentzian-broadened at T=300K, FWHM~2cm-1, for x(z,z)x,
    x(y,z)x, z(x,x)z, z(x,y)z), several genuinely active, non-degenerate
    modes were being incorrectly killed by that check -- e.g. modes at
    262.9, 366.5, 452.6, 484.2, 530.6, 710.7 cm-1 all show real, correctly
    predicted peaks in the reference x(y,z)x spectrum despite having no
    close-frequency partner. Degeneracy is not enforced.

    Returns a list of labels, length len(freqs), values in
    {"A_1", "E_1", "E_2", "dark"}.
    """
    freqs = np.asarray(freqs, dtype=float)
    n = len(freqs)
    is_acoustic = np.abs(freqs) <= acoustic_freq_tol_cm1
    magnitudes = np.array([np.max(np.abs(np.asarray(R_modes[i]).real)) for i in range(n)])

    optical_magnitudes = magnitudes[~is_acoustic]
    max_optical_magnitude = optical_magnitudes.max() if len(optical_magnitudes) else 0.0
    dark_threshold = max_optical_magnitude * magnitude_rel_floor

    labels = []
    for i in range(n):
        if is_acoustic[i] or magnitudes[i] < dark_threshold:
            labels.append("dark")
        else:
            labels.append(classify_tensor_shape(np.asarray(R_modes[i]).real))

    return labels

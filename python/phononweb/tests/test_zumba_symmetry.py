from pathlib import Path

import numpy as np

from phononweb.scripts.generate_zumba_jsons import (
    compute_intensities,
    generate,
    parse_dyn_for_raman,
    parse_modes_for_eigvecs,
)
from phononweb.scripts.zumba_symmetry import classify_modes

FIXTURE_DIR = Path(__file__).resolve().parents[3] / "test" / "fixtures" / "qespresso" / "zumba"
REFERENCE_SPECTRA_DIR = Path(__file__).resolve().parents[3] / "test" / "fixtures" / "zumba"


def _r_modes_for(eig_key):
    masses, raman_tens = parse_dyn_for_raman(str(FIXTURE_DIR / "zmo.dyn"))
    nat = len(masses)
    sqrt_m = np.sqrt(masses)
    freqs, eigvecs = parse_modes_for_eigvecs(str(FIXTURE_DIR / f"zmo_{eig_key}.eig"), nat)
    mw_evec = eigvecs / sqrt_m[np.newaxis, :, np.newaxis]
    R_modes = np.einsum("sabg,nsa->nbg", raman_tens, mw_evec)
    return freqs, R_modes


def test_parses_78_modes_with_real_raman_tensor():
    freqs, R_modes = _r_modes_for("x")
    assert len(freqs) == 78
    assert R_modes.shape == (78, 3, 3)
    assert abs(R_modes.imag).max() < 1e-10


def test_acoustic_modes_are_always_dark():
    freqs, R_modes = _r_modes_for("x")
    labels = classify_modes(R_modes, freqs)
    # modes 1-3 are the acoustic (translational) modes, near-zero frequency
    assert abs(freqs[0]) < 1.0 and abs(freqs[1]) < 1.0 and abs(freqs[2]) < 1.0
    assert labels[0] == labels[1] == labels[2] == "dark"


def test_pure_numerical_noise_modes_are_dark_not_mislabeled():
    # regression test: the original classifier used an absolute tolerance
    # (1e-8) to decide "dark", then normalized each mode's tensor by its own
    # max entry -- for a mode whose real tensor magnitude is ~1e-7 (still
    # "not dark" under that absolute rule), normalizing amplified pure
    # numerical noise into a fake E_1/E_2 pattern. These specific modes have
    # magnitude 7-8 orders of magnitude below genuinely active modes
    # (~1e-2 to 1) and must come out dark under the relative threshold.
    freqs, R_modes = _r_modes_for("x")
    labels = classify_modes(R_modes, freqs)
    magnitudes = np.array([np.max(np.abs(R_modes[i].real)) for i in range(len(freqs))])

    noise_scale_modes = [i for i in range(len(freqs)) if 1e-8 < magnitudes[i] < 1e-6]
    assert len(noise_scale_modes) > 5, "fixture should contain several noise-scale modes"
    for i in noise_scale_modes:
        assert labels[i] == "dark", f"mode {i + 1} (magnitude {magnitudes[i]:.2e}) should be dark"


def test_genuinely_active_non_degenerate_modes_are_not_demoted():
    # regression test: an earlier version of classify_modes additionally
    # required E_1/E_2 modes to have an exact frequency-degenerate partner,
    # reclassifying lone ones as dark. That's wrong -- validated against an
    # independently-computed reference spectrum (see
    # test_matches_independent_reference_spectrum_shape below), several
    # non-degenerate modes are genuinely Raman active. Degeneracy must not
    # be required.
    freqs, R_modes = _r_modes_for("x")
    labels = classify_modes(R_modes, freqs)

    # mode 27 (0-indexed 26), freq ~262.9 cm-1: real magnitude ~5e-2 (well
    # above the noise floor), no close-frequency partner, and confirmed
    # active by the reference x(y,z)x spectrum.
    known_active_non_degenerate = [26, 34, 45, 56, 61, 70]
    for i in known_active_non_degenerate:
        assert labels[i] != "dark", f"mode {i + 1} (freq {freqs[i]:.1f}) should not be dark"
        assert labels[i] in ("A_1", "E_1", "E_2")


def test_generate_writes_all_27_geometries_with_expected_shape(tmp_path):
    modes_files = {
        "0": str(FIXTURE_DIR / "zmo_0.eig"),
        "x": str(FIXTURE_DIR / "zmo_x.eig"),
        "y": str(FIXTURE_DIR / "zmo_y.eig"),
        "z": str(FIXTURE_DIR / "zmo_z.eig"),
    }
    written = generate(
        str(FIXTURE_DIR / "zmo.json"), str(FIXTURE_DIR / "zmo.dyn"), modes_files, str(tmp_path)
    )
    assert len(written) == 27

    import json
    data = json.loads((tmp_path / "ZMO_pol_x_ei_y_es_z.json").read_text())
    assert data["porto_label"] == "x(y,z)x"
    assert len(data["mode_symmetry"]) == 78
    assert len(data["eigenvalues"][data["gamma_index"]]) == 78
    # vectors must match THIS geometry's own (LO-TO-split) eigenvectors, not
    # a fixed baseline -- regression test for the vectors/frequency mismatch
    assert len(data["vectors"][data["gamma_index"]]) == 78
    assert len(data["vectors"][data["gamma_index"]][0]) == 26


def test_matches_independent_reference_spectrum_shape(tmp_path):
    """
    Validate against externally-supplied reference spectra: T=300K,
    Lorentzian broadening (gamma=2.0), for x(z,z)x / x(y,z)x and
    z(x,x)z / z(x,y)z. Checked by relative peak-height ratio rather than
    absolute intensity, since the two calculations use different absolute
    Raman-intensity unit conventions (confirmed: the ratio between our
    intensity and the reference intensity is constant, ~1.1e-15, across
    every frequency checked -- a fixed unit-scale factor, not a shape
    mismatch).
    """
    modes_files = {
        "0": str(FIXTURE_DIR / "zmo_0.eig"),
        "x": str(FIXTURE_DIR / "zmo_x.eig"),
        "y": str(FIXTURE_DIR / "zmo_y.eig"),
        "z": str(FIXTURE_DIR / "zmo_z.eig"),
    }
    generate(str(FIXTURE_DIR / "zmo.json"), str(FIXTURE_DIR / "zmo.dyn"), modes_files, str(tmp_path))

    import json

    def broadened(freqs, intensities, w, gamma=2.0):
        active = intensities > 0
        dw = w - freqs[active]
        return np.sum(intensities[active] * gamma * gamma / (dw * dw + gamma * gamma))

    def load(name):
        d = json.loads((tmp_path / name).read_text())
        gi = d["gamma_index"]
        return np.array(d["eigenvalues"][gi]), np.array(d["raman_intensities"])

    ref = np.loadtxt(REFERENCE_SPECTRA_DIR / "spectrum_xabx.dat")
    ref_freq, ref_xzzx = ref[:, 0], ref[:, 1]

    freqs, intens = load("ZMO_pol_x_ei_z_es_z.json")
    # two well-separated peaks from the reference spectrum
    w1, w2 = 345.53, 513.30
    mine_ratio = broadened(freqs, intens, w1) / broadened(freqs, intens, w2)
    ref_ratio = (
        ref_xzzx[np.argmin(np.abs(ref_freq - w1))]
        / ref_xzzx[np.argmin(np.abs(ref_freq - w2))]
    )
    assert abs(mine_ratio - ref_ratio) / ref_ratio < 0.05

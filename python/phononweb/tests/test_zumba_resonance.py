import numpy as np
import pytest

from phononweb.scripts.generate_zumba_resonance import (
    HA2CM,
    HA2EV,
    build_from_arrays,
    merge_excitons,
    to_polarizability,
)

EL_GRID = (0.0, 4.0)
N_EL = 3001


def _a1(a, b, leak=0.0):
    return np.array([[a, 0, leak], [0, a, 0], [leak, 0, b]], dtype=complex)


def _e2_pair(c):
    return [np.array([[c, 0, 0], [0, -c, 0], [0, 0, 0]], dtype=complex),
            np.array([[0, c, 0], [c, 0, 0], [0, 0, 0]], dtype=complex)]


def _e1_pair(d):
    return [np.array([[0, 0, d], [0, 0, 0], [d, 0, 0]], dtype=complex),
            np.array([[0, 0, 0], [0, 0, d], [0, d, 0]], dtype=complex)]


def _e1_pair_incoming_inplane(d):
    """A resonant E_1 pair that only scatters in-plane -> z light: R[in, out]."""
    return [np.array([[0, 0, 0], [0, 0, d], [0, 0, 0]], dtype=complex),
            np.array([[0, 0, d], [0, 0, 0], [0, 0, 0]], dtype=complex)]


def _synthetic(extra_mode=None, e1=None):
    """
    Raw PhdScripts-style tensor for 3 acoustic modes, an A_1 at 442.3 (with
    a symmetry-forbidden xz leak), an E_1 pair at 515.5, an E_2 pair at
    534.6 and an A_1 at 782.5 cm-1, all resonant at 2.4 eV.
    """
    freqs_cm = [5.0, 5.0, 5.0, 442.3, 515.5, 515.5, 534.6, 534.6, 782.5]
    e1 = e1 if e1 is not None else _e1_pair(0.2)
    cores = [np.zeros((3, 3), complex)] * 3 + [_a1(1.0, 0.3, leak=0.1)] + e1 + _e2_pair(0.5) + [_a1(0.4, 0.8)]
    if extra_mode is not None:
        freqs_cm.append(extra_mode[0])
        cores.append(extra_mode[1])
    freqs_ha = np.array(freqs_cm) / HA2CM
    el = np.linspace(*EL_GRID, N_EL)
    resonance = 1.0 / (el - 2.4 - 0.1j)
    alpha = resonance[:, None, None, None] * np.array(cores)[None]
    w_l = el[:, None] / HA2EV
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.sqrt(2.0 * freqs_ha[None, :]) * np.sqrt(w_l / np.abs(w_l - freqs_ha[None, :]))
        raw = alpha / scale[:, :, None, None]
    raw[0] = np.inf                                       # the driver's w_L = 0 row
    return raw, freqs_ha, alpha


def _excitons():
    exc_ev = np.array([2.4040, 2.4045, 2.6000, 3.5000])
    dip = np.array([[2.0, 1.0, 0.5, 9.0], [0.0, 1.0, 0.5, 9.0], [0.0, 0.0, 0.1, 9.0]], dtype=complex)
    return exc_ev / HA2EV, dip


def _build(**kw):
    raw, freqs_ha, _ = _synthetic(kw.pop("extra_mode", None), kw.pop("e1", None))
    exc_ha, dip = _excitons()
    return build_from_arrays(raw, EL_GRID, freqs_ha, exc_ha, dip, **kw)


def _geometry(data, gid):
    return next(g for g in data["geometries"] if g["id"] == gid)


def test_polarizability_undoes_the_drivers_scale():
    raw, freqs_ha, alpha = _synthetic()
    el = np.linspace(*EL_GRID, N_EL)
    recovered = to_polarizability(raw[1:], el[1:], freqs_ha)
    assert np.allclose(recovered[:, 3:], alpha[1:, 3:])


def test_peaks_are_grouped_and_labelled():
    data = _build()
    labels = [(round(p["freq"], 1), p["sym"], len(p["modes"])) for p in data["peaks"]]
    assert labels == [(442.3, "A_1", 1), (515.5, "E_1", 2), (534.6, "E_2", 2), (782.5, "A_1", 1)]


def test_forbidden_leakage_is_dropped_per_geometry():
    data = _build()
    syms = {g["id"]: {data["peaks"][i]["sym"] for i in g["peaks"]} for g in data["geometries"]}
    # the A_1 xz leak is visible in x(yz)x / y(xz)y but A_1 is not allowed there
    assert syms["y(xz)y"] == {"E_1"}
    assert syms["x(yz)x"] == {"E_1"}
    assert syms["x(yy)x"] == {"A_1", "E_2"}
    assert syms["x(zz)x"] == {"A_1"}
    assert syms["z(xy)z"] == {"E_2"}


def test_crossed_geometries_put_e_i_on_the_incoming_index():
    # alpha[y, z] = in y, out z: seen in -X(YZ)X (e_i = y, e_s = z) and
    # -Y(XZ)Y, which would both be dark with the light indices swapped
    data = _build(e1=_e1_pair_incoming_inplane(0.2))
    for gid in ("x(yz)x", "y(xz)y"):
        g = _geometry(data, gid)
        assert g["peaks"] == [1]
        assert max(g["intensity"][0]) > 0


def test_intensity_is_ei_alpha_es_summed_over_partners():
    data = _build()
    g = _geometry(data, "x(yy)x")
    energies = np.array(data["energies"])
    k = g["peaks"].index(2)                                # the E_2 pair
    # |e_y . alpha . e_y|^2 = c^2 |f|^2 from partner 1 only, times the Stokes factor
    from phononweb.scripts.generate_zumba_jsons import stokes_factor
    expected = 0.25 * np.abs(1.0 / (energies - 2.4 - 0.1j)) ** 2 * stokes_factor(np.array([534.6]), 300.0)[0]
    assert np.allclose(g["intensity"][k], expected, rtol=1e-3)


def test_window_and_offres_reference():
    data = _build()
    energies = np.array(data["energies"])
    assert energies[0] == pytest.approx(2.3) and energies[-1] == pytest.approx(3.0)
    assert len(energies) == 526
    g = _geometry(data, "x(zz)x")
    k = g["peaks"].index(3)                                # 782.5 A_1, zz = 0.8
    from phononweb.scripts.generate_zumba_jsons import stokes_factor
    expected = 0.64 * abs(1.0 / (1.0 - 2.4 - 0.1j)) ** 2 * stokes_factor(np.array([782.5]), 300.0)[0]
    assert g["offres"][k] == pytest.approx(expected, rel=1e-3)
    assert data["offres_energy"] == 1.0
    assert "static" not in g


def test_offres_needs_a_grid_point():
    with pytest.raises(ValueError, match="no grid point"):
        _build(offres_energy=1.0005)


def test_ref_peak_is_442_where_allowed_else_strongest_off_resonance():
    data = _build()
    assert _geometry(data, "x(yy)x")["ref_peak"] == 0
    # zz: 782.5 is stronger at 1 eV (0.8 vs 0.3) but 442.3 is still the reference
    assert _geometry(data, "x(zz)x")["ref_peak"] == 0
    assert _geometry(data, "x(yz)x")["ref_peak"] == 1
    assert _geometry(data, "z(xy)z")["ref_peak"] == 2


def test_labels_are_backscattering_porto_strings():
    data = _build()
    assert [g["label"] for g in data["geometries"]] == [
        "-X(YY)X", "-X(YZ)X", "-X(ZZ)X", "-Y(XX)Y", "-Y(XZ)Y", "-Y(ZZ)Y", "-Z(XY)Z"]


def test_excitons_within_1_meV_merge():
    exc_ha, dip = _excitons()
    energies, osc = merge_excitons(exc_ha * HA2EV, np.abs(dip) ** 2)
    assert len(energies) == 3
    # weights |d|^2 = 4 and 2 -> (4 * 2.4040 + 2 * 2.4045) / 6
    assert energies[0] == pytest.approx((4 * 2.4040 + 2 * 2.4045) / 6)
    assert osc[:, 0] == pytest.approx([5.0, 1.0, 0.0])


def test_a_quasi_continuum_does_not_chain_into_one_exciton():
    ladder = 2.9 + 0.0006 * np.arange(20)                  # 0.6 meV apart, 11.4 meV wide
    energies, _ = merge_excitons(ladder, np.ones((3, 20)))
    assert len(energies) > 1
    # every stick only gathers states within 1 meV of each other
    assert np.all(np.diff(energies) > 0.0005)


def test_excitons_outside_the_plot_range_are_dropped():
    data = _build()
    assert data["excitons"]["energies"] == pytest.approx([2.40417, 2.6], abs=1e-5)


def test_degenerate_group_with_a_strong_third_mode_is_rejected():
    with pytest.raises(ValueError, match="not one degenerate peak"):
        _build(extra_mode=(534.7, _a1(0.5, 0.5)))


def test_experiment_block_maps_to_the_442_peak():
    points = np.array([[2.30, 1.0], [2.43, 4.4]])
    data = _build(experiment=points)
    exp = data["experiment"]
    assert exp["calc_peak"] == 0
    assert exp["points"] == [[2.30, 1.0], [2.43, 4.4]]
    assert set(exp["geometries"]) == {"x(yy)x", "y(xx)y"}

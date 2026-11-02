import base64
import json

import numpy as np

from phononweb.scripts.generate_zumba_dispersion import build

NAT = 1
NMODES = 3
LATTICE = [[2.0, 0.0, 0.0], [-1.0, 3 ** 0.5, 0.0], [0.0, 0.0, 4.0]]

# Gamma -M -K -Gamma -A, 2 points per segment: 9 q-points (cartesian, 2pi/alat)
CRYST = [(0, 0, 0), (0.25, 0, 0), (0.5, 0, 0), (5 / 12, 1 / 6, 0), (1 / 3, 1 / 3, 0),
         (1 / 6, 1 / 6, 0), (0, 0, 0), (0, 0, 0.25), (0, 0, 0.5)]

MATDYN_IN = """&input
  flfrc='x.fc', q_in_band_form=.true., q_in_cryst_coord=.true.,
/
5
0.0 0.0 0.0  2   ! Gamma
0.5 0.0 0.0  2   ! M
0.333333 0.333333 0.0  2   ! K
0.0 0.0 0.0  2   ! Gamma
0.0 0.0 0.5  1   ! A
"""


def _write_inputs(tmp_path):
    lat = np.array(LATTICE)
    alat = np.linalg.norm(lat[0])
    rec = np.linalg.inv(lat).T * alat  # reduced -> cartesian in 2pi/alat
    lines = []
    for iq, c in enumerate(CRYST):
        q = np.array(c) @ rec
        lines.append(f"     q = {q[0]:12.8f} {q[1]:12.8f} {q[2]:12.8f}")
        for n in range(NMODES):
            lines.append(f"     freq ({n + 1:5d}) = 0.0 [THz] = {10.0 * iq + n:12.6f} [cm-1]")
            v = [0.0] * 6
            v[2 * n] = 1.0
            lines.append(" ( " + "  ".join(f"{x:9.6f}" for x in v) + " )")
    (tmp_path / "matdyn.modes").write_text("\n".join(lines) + "\n")
    (tmp_path / "matdyn.in").write_text(MATDYN_IN)
    base = {"name": "t", "natoms": NAT, "lattice": LATTICE, "atom_types": ["Zn"]}
    (tmp_path / "base.json").write_text(json.dumps(base))


def _decode(payload):
    raw = np.frombuffer(base64.b64decode(payload["data"]), "<i2")
    return raw.reshape(payload["shape"]) / payload["scale"]


def test_path_labels_qpoints_and_stride(tmp_path):
    _write_inputs(tmp_path)
    d = build(tmp_path / "matdyn.modes", tmp_path / "matdyn.in", tmp_path / "base.json", stride=2)

    assert len(d["qpoints"]) == 5
    assert d["highsym_qpts"] == [[0, "GAMMA"], [1, "M"], [2, "K"], [3, "GAMMA"], [4, "A"]]
    np.testing.assert_allclose(d["qpoints"], [CRYST[i] for i in range(0, 9, 2)], atol=1e-6)
    assert d["eigenvalues"][1] == [20.0, 21.0, 22.0]
    assert all(np.diff(d["distances"]) > 0)


def test_line_break_after_inner_gamma(tmp_path):
    # matdyn's inner Gamma carries the in-plane LO-TO limit: the Gamma-A line
    # must not be drawn from it
    _write_inputs(tmp_path)
    d = build(tmp_path / "matdyn.modes", tmp_path / "matdyn.in", tmp_path / "base.json", stride=2)
    assert d["line_breaks"] == [[0, 4], [4, 5]]


def test_vectors_round_trip(tmp_path):
    _write_inputs(tmp_path)
    d = build(tmp_path / "matdyn.modes", tmp_path / "matdyn.in", tmp_path / "base.json", stride=1)
    vec = _decode(d["vectors_compressed"])
    assert vec.shape == (9, NMODES, NAT, 3, 2)
    for n in range(NMODES):
        np.testing.assert_allclose(vec[0, n, 0, :, 0], np.eye(3)[n], atol=1e-3)
    assert np.all(vec[..., 1] == 0)

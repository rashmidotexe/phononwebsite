import json
import zipfile
from collections import Counter

import numpy as np
import pytest

from phononweb.scripts.prepare_rruff import (
    build_payload,
    cells_match,
    conventional_cell,
    formula_from_name,
    load_materials,
    match_materials,
    parse_cell,
    parse_chemistry,
    parse_spectrum,
    read_samples,
    reduced_formula,
    resample,
    select_spectra,
)


def rruff_file(name, rruff_id, chemistry, cell, points=((100, 0), (144, 10), (200, 1), (640, 5), (700, 0))):
    header = [
        f"##NAMES={name}",
        f"##RRUFFID={rruff_id}",
        f"##IDEAL CHEMISTRY={chemistry}",
        f"##URL=https://www.rruff.net/odr/rruff_sample/{rruff_id}",
        f"##CELL PARAMETERS={cell}",
        "##RAMAN WAVELENGTH=532",
    ]
    return "\n".join(header) + "\n\n" + "\n".join(f"{x}, {y}" for x, y in points) + "\n"


ANATASE_CELL = "a: 3.785 b: 3.785 c: 9.514 alpha: 90 beta: 90 gamma: 90 volume: 136.3 crystal system: tetragonal"
RUTILE_CELL = "a: 4.594 b: 4.594 c: 2.959 alpha: 90 beta: 90 gamma: 90 volume: 62.4 crystal system: tetragonal"


@pytest.mark.parametrize("text, expected", [
    ("Al_2_O_3_", {"Al": 2, "O": 3}),
    ("Ca(CO_3_)", {"Ca": 1, "C": 1, "O": 3}),
    ("KFe^2+^_3_(AlSi_3_O_10_)(OH)_2_", {"K": 1, "Fe": 3, "Al": 1, "Si": 3, "O": 12, "H": 2}),
    ("CaSO_4_·2H_2_O", {"Ca": 1, "S": 1, "O": 6, "H": 4}),
    ("Ti^4+^O_2_", {"Ti": 1, "O": 2}),
    ("Na_4_[Fe^2+^_10_Ti_2_]O_4_[Si_12_O_36_]", {"Na": 4, "Fe": 10, "Ti": 2, "O": 40, "Si": 12}),
])
def test_parse_chemistry(text, expected):
    assert parse_chemistry(text) == Counter(expected)


@pytest.mark.parametrize("text", [
    "(Mg,Fe)_2_SiO_4_",
    "Na_0.5_Cl",
    "(Mg_1-x_Fe_x_)O",
    "K[box]Fe_3_",
    "REE(PO_4_)",
    "Xx_2_O",
    "",
    None,
])
def test_parse_chemistry_rejects_non_stoichiometric(text):
    assert parse_chemistry(text) is None


def test_reduced_formula_and_website_names():
    assert reduced_formula(Counter({"Be": 2, "O": 2})) == (("Be", 1), ("O", 1))
    assert reduced_formula(formula_from_name("Ba1Te1")) == reduced_formula(parse_chemistry("BaTe"))
    assert formula_from_name("Qq2") is None


def test_parse_cell():
    cell = parse_cell(ANATASE_CELL)
    assert cell["family"] == "tetragonal"
    assert cell["volume"] == 136.3
    assert parse_cell("a: 1 b: 1 c: 1 alpha: 90 beta: 90 gamma: 120 volume: 1 crystal system: trigonal")["family"] == "hexagonal"
    assert parse_cell("a: 1 b: 1") is None


def test_parse_spectrum_sorts_points():
    header, x, y = parse_spectrum("##NAMES=Test\n\n300, 2\n100, 1\n200,3\n")
    assert header["NAMES"] == "Test"
    assert list(x) == [100, 200, 300]
    assert list(y) == [1, 3, 2]


def test_resample_on_grid_and_normalizes():
    start, values = resample(np.array([101.0, 141.0, 300.0]), np.array([0.0, 10.0, 0.0]), step=2.0, xmax=200.0)
    assert start == 102.0
    assert len(values) == 50
    assert max(values) == 1000
    assert values[(140 - 102) // 2] < values[(142 - 102) // 2]
    assert resample(np.array([1.0, 2.0]), np.array([1.0, 1.0])) is None


def tetragonal_structure(a, c, positions, numbers):
    return [[a, 0, 0], [0, a, 0], [0, 0, c]], positions, numbers


RUTILE = tetragonal_structure(4.65, 2.97,
                              [[0, 0, 0], [0.5, 0.5, 0.5], [0.305, 0.305, 0], [0.695, 0.695, 0],
                               [0.805, 0.195, 0.5], [0.195, 0.805, 0.5]],
                              [22, 22, 8, 8, 8, 8])


def test_conventional_cell_and_polymorph_matching():
    ours = conventional_cell(*RUTILE)
    assert ours["family"] == "tetragonal"
    assert ours["space_group"] == 136
    assert cells_match(parse_cell(RUTILE_CELL), ours)
    assert not cells_match(parse_cell(ANATASE_CELL), ours)


def test_read_select_and_match(tmp_path):
    archive = tmp_path / "excellent_unoriented.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("Rutile__R040049-1__Raman__532__0__unoriented__Raman_Data_Processed__1.txt",
                    rruff_file("Rutile", "R040049", "TiO_2_", RUTILE_CELL))
        zf.writestr("Rutile__R040049-1__Raman__532__0__unoriented__Raman_Data_RAW__1.txt",
                    rruff_file("Rutile", "R040049", "TiO_2_", RUTILE_CELL))
        zf.writestr("Anatase__R060277-1__Raman__532__0__unoriented__Raman_Data_Processed__2.txt",
                    rruff_file("Anatase", "R060277", "TiO_2_", ANATASE_CELL))
        zf.writestr("Olivine__R1-1__Raman__532__0__unoriented__Raman_Data_Processed__3.txt",
                    rruff_file("Olivine", "R1", "(Mg,Fe)_2_SiO_4_", RUTILE_CELL))
    oriented = tmp_path / "fair_oriented.zip"
    with zipfile.ZipFile(oriented, "w") as zf:
        zf.writestr("Rutile__R050031-1__Raman__780__0__depolarized__Raman_Data_Processed__4.txt",
                    rruff_file("Rutile", "R050031", "TiO_2_", RUTILE_CELL))

    samples = read_samples([archive, oriented], log=lambda *_: None)
    assert sorted(s["rruff_id"] for s in samples) == ["R040049", "R050031", "R060277"]

    manifest = tmp_path / "models.json"
    manifest.write_text(json.dumps([{"id": "2657", "name": "Ti2O4"}, {"id": "149", "name": "Si2"}]))
    materials = load_materials([(manifest, "phonondb")])

    cells = {"mp-2657": conventional_cell(*RUTILE)}
    matches = match_materials(samples, materials, cells.get, log=lambda *_: None)
    assert sorted(matches) == ["mp-2657"]
    assert sorted(s["rruff_id"] for s in matches["mp-2657"]) == ["R040049", "R050031"]

    chosen = select_spectra(matches["mp-2657"], max_spectra=3)
    assert [s["rruff_id"] for s in chosen] == ["R040049", "R050031"]

    payload = build_payload("mp-2657", chosen)
    assert payload["material"] == "mp-2657"
    assert payload["reference"]["doi"] == "10.1515/9783110417104-003"
    first = payload["spectra"][0]
    assert first["laser_nm"] == 532 and first["quality"] == "excellent" and not first["oriented"]
    assert first["start"] == 100.0 and first["step"] == 2.0
    assert max(first["intensity"]) == 1000

    excluded = match_materials(samples, materials, cells.get, keep_material=lambda m: False,
                               log=lambda *_: None)
    assert excluded == {}

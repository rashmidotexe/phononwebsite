#!/usr/bin/env python
#
# This file is part of the phononwebsite project
#
"""
Prepare experimental Raman spectra of the RRUFF project for the phonon website.

Reads the RRUFF Raman zip archives (https://rruff.info/zipped_data_files/raman/),
keeps the samples that match a material of the website, and writes one small
json per material to data/rruff/<mp-id>.json plus data/rruff/index.json.

A sample matches a material when
  1. the ideal chemistry of the sample has the same reduced formula,
  2. both are in the same crystal family, and
  3. the conventional cell agrees: sorted lattice lengths within
     --length-tolerance and volume within --volume-tolerance.
The space group is not stored in the RRUFF Raman files, but the cell is, and it
separates polymorphs such as anatase/rutile or calcite/aragonite.

Each spectrum is resampled on a regular grid and stored as integers 0-1000.

RRUFF: B. Lafuente, R. T. Downs, H. Yang and N. Stone, "The power of databases:
the RRUFF project", in Highlights in Mineralogical Crystallography,
T. Armbruster and R. M. Danisi eds., De Gruyter (2016), pp. 1-30.
"""

import argparse
import concurrent.futures
import gzip
import json
import math
import re
import sys
import time
import urllib.error
import urllib.request
import zipfile
from collections import Counter
from functools import reduce
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from phononweb.units import chemical_symbols

PHONONDB_URL = "https://henriquemiranda.github.io/phononwebsite/data/phonondb2017/"
MATERIALS_PROJECT_URL = "https://materialsproject-parsed.s3.amazonaws.com/ph-bandstructures/dfpt/"
RAMANDB_URL = "https://ramandb.oulu.fi/figs/tensors_yaml/"

REFERENCE = {
    "text": "RRUFF project (B. Lafuente, R. T. Downs, H. Yang and N. Stone, "
            "Highlights in Mineralogical Crystallography, De Gruyter (2016) 1-30)",
    "url": "https://rruff.info/",
    "doi": "10.1515/9783110417104-003",
}

QUALITY_RANK = {"excellent": 0, "fair": 1, "unrated": 2, "poor": 3}
ELEMENTS = set(chemical_symbols[1:])

#crystal families, trigonal and hexagonal are not told apart by the cell
FAMILIES = {
    "triclinic": "triclinic", "anorthic": "triclinic",
    "monoclinic": "monoclinic",
    "orthorhombic": "orthorhombic",
    "tetragonal": "tetragonal",
    "trigonal": "hexagonal", "rhombohedral": "hexagonal", "hexagonal": "hexagonal",
    "cubic": "cubic", "isometric": "cubic",
}


def space_group_family(number):
    if number <= 2: return "triclinic"
    if number <= 15: return "monoclinic"
    if number <= 74: return "orthorhombic"
    if number <= 142: return "tetragonal"
    if number <= 194: return "hexagonal"
    return "cubic"


def parse_chemistry(text):
    """
    RRUFF ideal chemistry, e.g. 'KFe^2+^_3_(AlSi_3_O_10_)(OH)_2_' -> Counter
    None for solid solutions, vacancies, variable or non-integer occupations
    """
    if not text:
        return None
    text = re.sub(r"\^[^^]*\^", "", text.strip())
    if re.search(r"[,?]|box|REE|\bx\b", text):
        return None
    parts = re.split(r"\s*[·•]\s*|&#183;", text)
    total = Counter()
    for part in parts:
        match = re.match(r"^(\d+)(.*)$", part)
        factor, part = (int(match.group(1)), match.group(2)) if match else (1, part)
        counts = _parse_group(part)
        if counts is None:
            return None
        for element, n in counts.items():
            total[element] += n * factor
    if not total or any(abs(n - round(n)) > 1e-6 for n in total.values()):
        return None
    return Counter({element: int(round(n)) for element, n in total.items()})


def _parse_group(text):
    #RRUFF uses (), [] and {} alike for grouping
    tokens = re.findall(r"[A-Z][a-z]?|_[^_]*_|\S", text)
    stack = [Counter()]
    last = None
    for token in tokens:
        if token in "([{":
            stack.append(Counter())
            last = None
        elif token in ")]}":
            if len(stack) == 1:
                return None
            last = stack.pop()
            for element, n in last.items():
                stack[-1][element] += n
        elif token.startswith("_"):
            try:
                factor = float(token.strip("_"))
            except ValueError:
                return None
            if last is None:
                return None
            for element, n in last.items():
                stack[-1][element] += n * (factor - 1)
            last = None
        elif token in ELEMENTS:
            last = Counter({token: 1})
            stack[-1][token] += 1
        else:
            return None
    if len(stack) != 1:
        return None
    return stack[0]


def reduced_formula(counts):
    """Counter -> hashable formula with the counts divided by their gcd"""
    divisor = reduce(math.gcd, counts.values())
    return tuple(sorted((element, n // divisor) for element, n in counts.items()))


def formula_from_name(name):
    """website names such as 'Be2O2' or 'Ba1Te1' -> Counter"""
    counts = Counter()
    for element, n in re.findall(r"([A-Z][a-z]?)(\d*)", name):
        if element not in ELEMENTS:
            return None
        counts[element] += int(n) if n else 1
    return counts or None


def parse_cell(text):
    """'a: 5.465 b: 9.281 ... volume: 512.9 crystal system: monoclinic' -> dict"""
    if not text:
        return None
    values = {}
    for key in ("a", "b", "c", "alpha", "beta", "gamma", "volume"):
        match = re.search(r"\b" + key + r":\s*([-+0-9.eE]+)", text)
        if not match:
            return None
        values[key] = float(match.group(1))
    match = re.search(r"crystal system:\s*([A-Za-z]+)", text)
    family = FAMILIES.get(match.group(1).lower()) if match else None
    if family is None or values["volume"] <= 0:
        return None
    values["family"] = family
    return values


def parse_spectrum(text):
    """RRUFF txt file -> (header dict, x, y)"""
    header = {}
    xs, ys = [], []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("##"):
            key, _, value = line[2:].partition("=")
            header[key.strip().upper()] = value.strip()
        elif line:
            fields = re.split(r"[,\s]+", line)
            if len(fields) >= 2:
                try:
                    xs.append(float(fields[0]))
                    ys.append(float(fields[1]))
                except ValueError:
                    pass
    order = np.argsort(xs)
    return header, np.asarray(xs)[order], np.asarray(ys)[order]


def resample(x, y, step=2.0, xmax=1400.0):
    """
    spectrum on a regular grid from the first point to xmax, as integers 0-1000
    returns (start, values) or None when nothing is left
    """
    if len(x) < 2:
        return None
    start = math.ceil(max(x[0], 0.0) / step) * step
    stop = min(x[-1], xmax)
    if stop - start < 10 * step:
        return None
    grid = np.arange(start, stop + 1e-9, step)
    values = np.interp(grid, x, y)
    values = values - min(values.min(), 0.0)
    peak = values.max()
    if not peak > 0:
        return None
    return float(start), [int(round(v)) for v in values / peak * 1000]


def conventional_cell(lattice, positions, numbers, symprec=1e-2):
    """our structure -> dict with the crystal family, sorted lengths and volume"""
    import spglib
    dataset = spglib.get_symmetry_dataset((lattice, positions, numbers), symprec=symprec)
    if dataset is None:
        return None
    std = np.asarray(dataset.std_lattice)
    return {
        "family": space_group_family(dataset.number),
        "space_group": int(dataset.number),
        "lengths": sorted(float(v) for v in np.linalg.norm(std, axis=1)),
        "volume": float(abs(np.linalg.det(std))),
    }


def cells_match(rruff_cell, our_cell, length_tolerance=0.06, volume_tolerance=0.15):
    if rruff_cell["family"] != our_cell["family"]:
        return False
    if abs(our_cell["volume"] / rruff_cell["volume"] - 1) > volume_tolerance:
        return False
    rruff_lengths = sorted((rruff_cell["a"], rruff_cell["b"], rruff_cell["c"]))
    return all(abs(ours / theirs - 1) <= length_tolerance
               for ours, theirs in zip(our_cell["lengths"], rruff_lengths))


def describe_file(zip_name, file_name):
    """quality and orientation from the archive name, laser from the file name"""
    stem = Path(zip_name).stem.lower()
    quality = next((q for q in QUALITY_RANK if stem.startswith(q)), "unrated")
    fields = Path(file_name).name.split("__")
    laser = fields[3] if len(fields) > 3 and fields[3].isdigit() else None
    return {"quality": quality, "oriented": not stem.endswith("unoriented"), "laser": laser}


def read_samples(zip_paths, log=print):
    """processed spectra of the RRUFF archives with a usable chemistry and cell"""
    samples = []
    for zip_path in zip_paths:
        kept = 0
        with zipfile.ZipFile(zip_path) as archive:
            for info in archive.infolist():
                name = info.filename
                if info.is_dir() or "Processed" not in name or not name.endswith(".txt"):
                    continue
                text = archive.read(info).decode("utf-8", errors="replace")
                header, x, y = parse_spectrum(text)
                counts = parse_chemistry(header.get("IDEAL CHEMISTRY"))
                cell = parse_cell(header.get("CELL PARAMETERS"))
                if counts is None or cell is None or len(x) < 2:
                    continue
                sample = describe_file(zip_path, name)
                sample.update({
                    "name": header.get("NAMES", "").split(",")[0].strip(),
                    "rruff_id": header.get("RRUFFID", ""),
                    "url": header.get("URL", ""),
                    "laser": header.get("RAMAN WAVELENGTH") or sample["laser"],
                    "formula": reduced_formula(counts),
                    "has_hydrogen": "H" in counts,
                    "cell": cell,
                    "file": name,
                    "x": x,
                    "y": y,
                })
                samples.append(sample)
                kept += 1
        log(f"{Path(zip_path).name}: {kept} processed spectra with a usable chemistry and cell")
    return samples


def load_materials(manifests):
    """[(path, source)] of models.json -> {mp-id: {'formula', 'source'}}, first source wins"""
    materials = {}
    for path, source in manifests:
        for entry in json.loads(Path(path).read_text()):
            material_id = "mp-" + str(entry["id"]).replace("mp-", "")
            counts = formula_from_name(entry.get("name", ""))
            if counts and material_id not in materials:
                materials[material_id] = {"formula": reduced_formula(counts), "source": source}
    return materials


def _fetch(url, retries=3):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return None
            if attempt == retries - 1:
                raise
        except urllib.error.URLError:
            if attempt == retries - 1:
                raise
        time.sleep(2 ** attempt)


def fetch_structure(material_id, source, cache_dir):
    """(lattice, reduced positions, atomic numbers) of a website material, cached"""
    cache = Path(cache_dir) / (material_id + ".json")
    if cache.exists():
        data = json.loads(cache.read_text())
    else:
        if source == "phonondb":
            raw = _fetch(PHONONDB_URL + material_id + ".json.gz")
            if raw is None:
                return None
            payload = json.loads(gzip.decompress(raw))
            data = {
                "lattice": payload["lattice"],
                "positions": payload["atom_pos_red"],
                "numbers": payload["atom_numbers"],
            }
        else:
            raw = _fetch(MATERIALS_PROJECT_URL + material_id + ".json.gz")
            if raw is None:
                return None
            structure = json.loads(gzip.decompress(raw))["structure"]
            data = {
                "lattice": structure["lattice"]["matrix"],
                "positions": [site["abc"] for site in structure["sites"]],
                "numbers": [chemical_symbols.index(site["species"][0]["element"])
                            for site in structure["sites"]],
            }
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(data))
    return data["lattice"], data["positions"], data["numbers"]


def has_ramandb(material_id, cache_dir):
    """whether the Computational Raman Database has tensors for the material, cached"""
    cache = Path(cache_dir) / "ramandb.json"
    known = json.loads(cache.read_text()) if cache.exists() else {}
    if material_id not in known:
        known[material_id] = _fetch(RAMANDB_URL + material_id + ".yaml") is not None
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(known, indent=0, sort_keys=True))
        time.sleep(0.2)
    return known[material_id]


def select_spectra(samples, max_spectra=3):
    """best spectra first: quality, unoriented, one per RRUFF sample"""
    ranked = sorted(samples, key=lambda s: (QUALITY_RANK[s["quality"]], s["oriented"],
                                            s["rruff_id"], s["file"]))
    chosen, seen = [], set()
    for sample in ranked:
        if sample["rruff_id"] in seen:
            continue
        seen.add(sample["rruff_id"])
        chosen.append(sample)
        if len(chosen) == max_spectra:
            break
    return chosen


def build_payload(material_id, samples, step=2.0):
    spectra = []
    for sample in samples:
        xmax = 4000.0 if sample["has_hydrogen"] else 1400.0
        grid = resample(sample["x"], sample["y"], step=step, xmax=xmax)
        if grid is None:
            continue
        start, values = grid
        spectra.append({
            "name": sample["name"],
            "rruff_id": sample["rruff_id"],
            "url": sample["url"],
            "laser_nm": int(sample["laser"]) if sample["laser"] and str(sample["laser"]).isdigit() else None,
            "quality": sample["quality"],
            "oriented": sample["oriented"],
            "start": start,
            "step": step,
            "intensity": values,
        })
    if not spectra:
        return None
    return {"material": material_id, "reference": REFERENCE, "spectra": spectra}


def match_materials(samples, materials, get_cell, keep_material=lambda material_id: True,
                    length_tolerance=0.06, volume_tolerance=0.15, log=print):
    """{mp-id: [samples]} for the website materials matching RRUFF samples"""
    by_formula = {}
    for sample in samples:
        by_formula.setdefault(sample["formula"], []).append(sample)

    candidates = sorted(m for m, info in materials.items() if info["formula"] in by_formula)
    log(f"{len(candidates)} website materials share a reduced formula with a RRUFF sample")

    matches = {}
    for material_id in candidates:
        our_cell = get_cell(material_id)
        if our_cell is None:
            continue
        found = [s for s in by_formula[materials[material_id]["formula"]]
                 if cells_match(s["cell"], our_cell, length_tolerance, volume_tolerance)]
        if found and keep_material(material_id):
            matches[material_id] = found
    log(f"{len(matches)} materials match a RRUFF sample by formula and cell")
    return matches


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("zips", nargs="+", help="RRUFF Raman zip archives")
    parser.add_argument("--phonondb-manifest", default="data/phonondb2017/models.json")
    parser.add_argument("--mp-manifest", default="data/mpdb/models.json")
    parser.add_argument("-o", "--output-dir", default="data/rruff")
    parser.add_argument("--cache-dir", default=".rruff-cache",
                        help="where downloaded structures and lookups are kept between runs")
    parser.add_argument("--max-spectra", type=int, default=3)
    parser.add_argument("--step", type=float, default=2.0, help="grid spacing in cm-1")
    parser.add_argument("--length-tolerance", type=float, default=0.06)
    parser.add_argument("--volume-tolerance", type=float, default=0.15)
    parser.add_argument("--all-materials", action="store_true",
                        help="keep materials without Raman tensors in the Computational Raman Database")
    parser.add_argument("--workers", type=int, default=8)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    samples = read_samples(args.zips)
    materials = load_materials([(args.phonondb_manifest, "phonondb"), (args.mp_manifest, "mp")])

    formulas = {s["formula"] for s in samples}
    candidates = [m for m, info in materials.items() if info["formula"] in formulas]
    print(f"downloading {len(candidates)} structures (cached in {args.cache_dir})")
    cells = {}

    def load_cell(material_id):
        try:
            structure = fetch_structure(material_id, materials[material_id]["source"],
                                        Path(args.cache_dir) / "structures")
        except Exception as error:
            print(f"  {material_id}: {error}")
            return material_id, None
        return material_id, conventional_cell(*structure) if structure else None

    with concurrent.futures.ThreadPoolExecutor(args.workers) as pool:
        for material_id, cell in pool.map(load_cell, candidates):
            cells[material_id] = cell

    keep = (lambda m: True) if args.all_materials else \
        (lambda m: has_ramandb(m, args.cache_dir))
    matches = match_materials(samples, materials, cells.get, keep,
                              args.length_tolerance, args.volume_tolerance)

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for old in output.glob("mp-*.json"):
        old.unlink()
    index = []
    total = 0
    for material_id in sorted(matches):
        payload = build_payload(material_id, select_spectra(matches[material_id], args.max_spectra),
                                step=args.step)
        if payload is None:
            continue
        text = json.dumps(payload, separators=(",", ":"))
        (output / (material_id + ".json")).write_text(text)
        index.append(material_id)
        total += len(text)
    (output / "index.json").write_text(json.dumps(
        {"reference": REFERENCE, "materials": index}, indent=1))
    print(f"wrote {len(index)} materials ({total / 1e6:.2f} MB) to {output}")


if __name__ == "__main__":
    main()

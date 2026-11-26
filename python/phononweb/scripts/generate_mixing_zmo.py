#!/usr/bin/env python3
"""
Zn2Mo3O8 (ZMO) for the mixing page: the q2r force constants of the Zumba
calculation (zmo.fc, 3x3x2 grid, with Born charges and epsilon_infinity) written
in the dynamical-matrix payload of the PhononDB files (phonopy-dynamical-matrix-v1),
so that the page evaluates D at any q, as for a PhononDB material.

    D(q) = SR(q) + dd(q)

SR are the short-range force constants of zmo.fc (q2r has taken the dipole-dipole
part out), with matdyn's asr='simple', on the Wigner-Seitz images of the 3x3x2
supercell. dd is the Gonze reciprocal sum with QE's Ewald parameter (rgd_blk:
the alpha written in zmo.fc, in (2 pi / alat)^2, i.e. lambda = sqrt(alpha) / alat),
so that SR + dd is matdyn's D (the antisymmetric part q2r leaves in the short-range
force constants of the Mo atoms then cancels).
The band path is the one of the Zumba dispersion (zumbadb/ZMO_dispersion.json).

usage: generate_mixing_zmo.py zmo.fc zumbadb/ZMO_dispersion.json out.json.gz [--sign +1|-1]
"""
import argparse
import gzip
import json
import math

import numpy as np

BOHR = 0.529177210903          # Angstrom
RY = 13.605693122994           # eV
AMU_RY = 911.44424310865645    # QE mass unit (Ry atomic units) per amu
PHONOPY_THZ = 15.633302300230191
NAC_UNIT = 14.39965172592227   # e^2 / (4 pi eps0) in eV Angstrom


def parse_fc(path):
    with open(path) as handle:
        tokens = handle.read().split('\n')
    lines = iter(tokens)
    head = next(lines).split()
    ntyp, nat, ibrav = int(head[0]), int(head[1]), int(head[2])
    celldm = [float(v) for v in head[3:9]]
    if ibrav != 4:
        raise ValueError('expected a hexagonal cell (ibrav = 4)')
    species = {}
    for _ in range(ntyp):
        line = next(lines)
        parts = line.split("'")
        species[int(parts[0])] = (parts[1].strip(), float(parts[2]) / AMU_RY)
    types, tau = [], []
    for _ in range(nat):
        parts = next(lines).split()
        types.append(species[int(parts[1])][0])
        tau.append([float(v) for v in parts[2:5]])
    masses = []
    for name in types:
        masses.append(next(mass for key, (label, mass) in species.items() if label == name))
    flag_line = next(lines).split()
    if flag_line[0] != 'T':
        raise ValueError('zmo.fc has no Born charges / dielectric tensor')
    # QE writes its Ewald parameter (in (2 pi / alat)^2) after the flag; older versions
    # used alpha = 1
    alpha = float(flag_line[1]) if len(flag_line) > 1 else 1.0
    eps = [[float(v) for v in next(lines).split()] for _ in range(3)]
    born = []
    for _ in range(nat):
        next(lines)
        born.append([[float(v) for v in next(lines).split()] for _ in range(3)])
    nr = [int(v) for v in next(lines).split()]
    frc = np.zeros((nr[0], nr[1], nr[2], 3, 3, nat, nat))
    for _ in range(9 * nat * nat):
        i, j, na, nb = [int(v) for v in next(lines).split()]
        for _ in range(nr[0] * nr[1] * nr[2]):
            m1, m2, m3, value = next(lines).split()
            frc[int(m1) - 1, int(m2) - 1, int(m3) - 1, i - 1, j - 1, na - 1, nb - 1] = float(value)
    alat = celldm[0] * BOHR
    lattice = alat * np.array([[1, 0, 0], [-0.5, math.sqrt(3) / 2, 0], [0, 0, celldm[2]]])
    tau = alat * np.array(tau)
    return dict(types=types, masses=masses, lattice=lattice, positions_car=tau, eps=np.array(eps),
                born=np.array(born), nr=nr, frc=frc, alat=alat, alpha=alpha)


def impose_simple_asr(frc):
    # matdyn asr='simple': the self term R = 0 takes minus the sum of all the others
    nat = frc.shape[-1]
    for na in range(nat):
        total = frc[:, :, :, :, :, na, :].sum(axis=(0, 1, 2, 5))
        frc[0, 0, 0, :, :, na, na] -= total
    return frc


def payload_from_fc(data, sign):
    lattice = data['lattice']
    nr = data['nr']
    nat = len(data['types'])
    inverse = np.linalg.inv(lattice)
    red = data['positions_car'] @ inverse
    super_lattice = np.diag(nr) @ lattice
    images = [np.array([a, b, c]) for a in range(-2, 3) for b in range(-2, 3) for c in range(-2, 3)]
    fc_to_ev = RY / BOHR ** 2
    cells = [(m1, m2, m3) for m3 in range(nr[2]) for m2 in range(nr[1]) for m1 in range(nr[0])]
    s2pp, fc_compact, svecs, multiplicity = [], [[None] * (len(cells) * nat) for _ in range(nat)], [], []
    super_index = 0
    for cell in cells:
        for nb in range(nat):
            s2pp.append(nb)
            row = []
            for na in range(nat):
                # vectors from atom na (cell 0) to atom nb in cell R (and its supercell images)
                base = np.array(cell, dtype=float) + red[nb] - red[na]
                candidates = [base + image * np.array(nr) for image in images]
                lengths = [np.linalg.norm(v @ lattice) for v in candidates]
                shortest = min(lengths)
                chosen = [v for v, length in zip(candidates, lengths) if length < shortest + 1e-4]
                row.append([len(chosen), len(svecs)])
                svecs.extend(v.tolist() for v in chosen)
                # the force constant of the pair (na at 0, nb at R): QE stores it under +R or
                # under -R depending on its phase convention (checked against matdyn)
                index = tuple((sign * c) % n for c, n in zip(cell, nr))
                fc_compact[na][super_index] = (data['frc'][index][:, :, na, nb] * fc_to_ev).tolist()
            multiplicity.append(row)
            super_index += 1
    volume = abs(np.linalg.det(lattice))
    lam = math.sqrt(data['alpha']) / data['alat']
    trace = np.trace(data['eps'])
    g_cutoff = lam * math.sqrt(-12 * math.log(1e-10) / trace)
    reciprocal = np.linalg.inv(lattice).T
    radius = math.ceil(g_cutoff * max(np.linalg.norm(lattice, axis=1)))
    g_list = []
    for a in range(-radius, radius + 1):
        for b in range(-radius, radius + 1):
            for c in range(-radius, radius + 1):
                g = a * reciprocal[0] + b * reciprocal[1] + c * reciprocal[2]
                if g @ g < g_cutoff ** 2:
                    g_list.append(g.tolist())
    born = data['born']
    pos = data['positions_car']
    dd_real = np.zeros((nat, 3, 3))
    dd_imag = np.zeros((nat, 3, 3))
    for g in np.array(g_list):
        if np.linalg.norm(g) < 1e-5:
            continue
        geg = g @ data['eps'] @ g
        pref = math.exp(-geg / (4 * lam * lam)) / geg
        zg = np.einsum('k,ika->ia', g, born)  # Z_i^T G
        for i in range(nat):
            phase = 2 * math.pi * (pos[i] - pos) @ g
            outer = np.einsum('a,jb->jab', zg[i], zg) * pref
            dd_real[i] += np.einsum('jab,j->ab', outer, np.cos(phase))
            dd_imag[i] += np.einsum('jab,j->ab', outer, np.sin(phase))
    # QE subtracts this q = 0 sum as it is (rgd_blk), not its hermitian part as phonopy
    # does: with Mo's non-symmetric Born tensors it has an antisymmetric part, which
    # cancels the one q2r left in the short-range force constants
    return {
        'format': 'phonopy-dynamical-matrix-v1',
        'acoustic_sum_rule': 'off',
        'primitive_natoms': nat,
        'supercell_natoms': len(s2pp),
        'primitive_lattice': lattice.tolist(),
        'frequency_conversion_factor': PHONOPY_THZ,
        'masses': data['masses'],
        'force_constants_compact': fc_compact,
        'shortest_vectors': svecs,
        'multiplicity': multiplicity,
        's2pp_map': s2pp,
        'nac': {
            'method': 'gonze',
            'born': born.tolist(),
            'dielectric': data['eps'].tolist(),
            'nac_factor': NAC_UNIT * 4 * math.pi / volume,
            'positions_car': pos.tolist(),
            'g_list': g_list,
            'dd_q0': {'real': dd_real.tolist(), 'imag': dd_imag.tolist()},
            'lambda': lam,
            # QE's Ewald parameter, not phonopy's rule: the short-range force constants
            # were split with it
            'lambda_fixed': True,
            'q_direction_tolerance': 1e-5,
        },
    }


def straight_segments(path):
    # matdyn's band path is one polyline (Gamma-M-K-Gamma-A-L-H-A); the mixing page
    # samples straight segments, one per pair of labels, each with its two ends
    marks = path['highsym_qpts']
    qpoints, distances, line_breaks, labels = [], [], [], []
    for (start, label_a), (end, label_b) in zip(marks[:-1], marks[1:]):
        first = len(qpoints)
        qpoints.extend(path['qpoints'][start:end + 1])
        distances.extend(path['distances'][start:end + 1])
        line_breaks.append([first, len(qpoints)])
        if not labels:
            labels.append([first, label_a])
        labels.append([len(qpoints) - 1, label_b])
    return dict(qpoints=qpoints, distances=distances, line_breaks=line_breaks, highsym_qpts=labels)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fc')
    parser.add_argument('dispersion')
    parser.add_argument('output')
    # -1: QE stores the pair (na at R, nb at 0) under R (matdyn phase exp(-i q.R)); checked
    # against the matdyn dispersion to 0.1 cm-1
    parser.add_argument('--sign', type=int, default=-1)
    args = parser.parse_args()
    data = parse_fc(args.fc)
    data['frc'] = impose_simple_asr(data['frc'])
    with open(args.dispersion) as handle:
        path = straight_segments(json.load(handle))
    lattice = data['lattice']
    red = data['positions_car'] @ np.linalg.inv(lattice)
    out = {
        'name': 'Zn2Mo3O8',
        'natoms': len(data['types']),
        'lattice': lattice.tolist(),
        'atom_types': data['types'],
        'atom_numbers': [{'Zn': 30, 'Mo': 42, 'O': 8}[t] for t in data['types']],
        'formula': 'Zn4Mo6O16',
        'qpoints': path['qpoints'],
        'repetitions': [2, 2, 1],
        'atom_pos_red': red.tolist(),
        'atom_pos_car': data['positions_car'].tolist(),
        'masses': data['masses'],
        'distances': path['distances'],
        'highsym_qpts': path['highsym_qpts'],
        'line_breaks': path['line_breaks'],
        'average_mass': float(np.mean(data['masses'])),
        'mode_amplitude_convention': 'avg-mass-normalized',
        'dynamical_matrix': payload_from_fc(data, args.sign),
    }
    with gzip.open(args.output, 'wt') as handle:
        json.dump(out, handle)


if __name__ == '__main__':
    main()

/*
 * mixingphonons.js — virtual-crystal mixing of two materials of the same
 * structure family on the mixing page: 2-atom zinc blende / diamond, 2-atom
 * rocksalt (Fm-3m), 3-atom fluorite / antifluorite (Fm-3m), 4-atom wurtzite /
 * lonsdaleite (P6_3mc), 8-atom chalcopyrite (I-42d), a ternary wurtzite
 * (beta-NaFeO2 type Pna2_1, enargite type Pmn2_1, Cu2SiS3 type Cmc2_1;
 * WURTZITE_SUPERSTRUCTURES), a cation-ordered ternary rocksalt (alpha-NaFeO2
 * type R-3m, gamma-LiFeO2 type I4_1/amd, Li2ZrO3 type C2/c, ...;
 * ROCKSALT_SUPERSTRUCTURES) or fluorite / antifluorite (scheelite, CaUO4 type,
 * MgAgAs type LiMgN, ...; FLUORITE_SUPERSTRUCTURES), or an A B X3 perovskite /
 * antiperovskite: the 5-atom cell (cubic or polar) or a supercell of the
 * octahedral tilts (Pnma, R-3c, ...; PEROVSKITE_SUPERSTRUCTURES). A binary (a
 * 5-atom perovskite) with force constants can also be folded into the cell of
 * an ordered ternary (a tilted perovskite) on its lattice. Transition-metal
 * dichalcogenides mix within their polytype (2H, 1T, 3R; TMD_SUPERSTRUCTURES).
 * Rutiles (A X2, with the CaCl2-type distortion) mix with each other and fold into
 * the cation-ordered rutiles: trirutile A B2 X6, A B X4, A B3 X8
 * (RUTILE_SUPERSTRUCTURES). Corundums (A2 X3) mix with each other and fold into
 * ilmenite (A B X3) and A3 B X6 on the same cell (CORUNDUM_SUPERSTRUCTURES). CsCl-type
 * A X mix with each other; the zinc-blende-derived (stuffed tetrahedral) cells with
 * all four fcc sites filled, Li3Bi-type A3 X, fold into the full and inverse Heuslers
 * A2 B X (BCC_SUPERSTRUCTURES; the half-Heuslers, three sites, are the MgAgAs type
 * of FLUORITE_SUPERSTRUCTURES). Materials of different structures on the same
 * lattice mix too, site by site on a common cell (LATTICE_WORLDS): zinc blende,
 * chalcopyrite, rocksalt, fluorite, antifluorite, the Heuslers, CsCl and their
 * orderings on the fcc sites, perovskite, antiperovskite and CsCl on the simple
 * cubic ones; a site only one of them fills is a vacancy in the other, occupied in
 * proportion to the mixture (computeUnionMixedPhonon). Force constants
 * and site masses are mixed separately:
 *
 *     Phi_k(q)   = M_k^1/2 D_k(q) M_k^1/2
 *     M_x        = (1 - x) M_1 + x M_2
 *     D_alloy(q) = M_x^-1/2 [(1 - x) Phi_1(q) + x Phi_2(q)] M_x^-1/2
 *
 * D is the mass-weighted dynamical matrix in THz^2. Mixing D directly would
 * break the acoustic sum rule whenever the two materials have different mass
 * ratios (the rigid translation (sqrt(m_1), sqrt(m_2)) is a zero mode of D_1 or
 * D_2 but not of their average), lifting the acoustic branches off zero at
 * Gamma. The sum rule is linear in Phi, so mixing Phi keeps it for every x.
 *
 * Gauge: what is mixed are the force constants Phi_ij(R), i.e. D in the cell
 * gauge exp(2 pi i q.R) with the same lattice vector R for corresponding atoms.
 * D is kept in the atom gauge exp(2 pi i q.(R + tau_j - tau_i)), which depends
 * on the internal positions; when the two end-members' positions differ (u,
 * ordered anions, octahedral tilts) both are brought to the gauge of the mixed
 * positions before mixing (changeGauge). Mixing the atom-gauge matrices directly
 * would add matrices in two different gauges: exact at Gamma and at the ends,
 * but not periodic in q (omega(q + G) != omega(q)).
 *
 * Long-range (LO-TO) part: the full Ewald dipole-dipole sum of each end-member
 * (reciprocal part as in phonopy, plus the real-space part phonopy leaves in its
 * short-range force constants; gonze.js getGonzeDipoleDipole) is subtracted, the
 * remainders are mixed, and the full sum of the virtual crystal is rebuilt from
 * the linearly mixed Born charges and dielectric tensor on the mixed lattice and
 * positions. The full sum does not depend on the Ewald parameter, so the split
 * is the same in any cell (a binary folded into a supercell has the folded
 * short-range part of its own primitive cell). phonopy's reciprocal part alone
 * depends on its Ewald parameter, which is set per cell from the volume and
 * dielectric tensor: subtracting only that part made the mixture depend on the
 * cell and on the two materials' parameters. For Materials Project the Born
 * charges and dielectric tensor come from the MP DFPT collection
 * (mpdielectric.js). If they are not available the long-range term is mixed
 * along with the rest of Phi.
 *
 * Before mixing, both end-members are brought to a common representation:
 *
 *   gauge  : "atom" Bloch phase exp(2 pi i q.(R + tau_j - tau_i)), which is what
 *            phonopy force constants give (PhononDB). Materials Project
 *            eigendisplacements are in the "cell" gauge exp(2 pi i q.R) and are
 *            converted with the atomic positions of the file.
 *   sites  : zinc blende: canonical order [cation, anion] with the anion at
 *            +(1/4,1/4,1/4) from the cation; rocksalt the same with the anion
 *            at (1/2,1/2,1/2). Fluorite: [fcc site, tetrahedral site at
 *            +(1/4,1/4,1/4), the one at -(1/4,1/4,1/4)], antifluorite (anion on
 *            the fcc site: Li2O, Mg2Si) as a separate family. Chalcopyrite: the 8 sites of the
 *            I-42d template (lower-valence cation on 4a, the other on 4b, anions
 *            on 8d), matched with any origin. Wurtzite: [cation, cation, anion,
 *            anion], cations at (1/3,2/3,0) and (2/3,1/3,1/2) with each anion
 *            u c above its cation, matched with any origin after rotating the
 *            cell onto a1 || x. A structure stored as the inversion image of
 *            the template (the opposite polarity in wurtzite) is handled in the
 *            atom gauge by D(q) -> D(q)^*; every other setting of the same
 *            lattice differs from these by a symmetry of the crystal.
 *            Ternary wurtzite / rocksalt / (anti)fluorite: the cell is read as a supercell of the
 *            parent and its sites matched, role by role, onto a template on that
 *            supercell (any basis, origin and orientation; an improper rotation
 *            is again handled as the inversion image). The cation roles follow
 *            the number of atoms, then the formal valence that balances the
 *            anion charge (Tl+ in DyTlS2, Tl3+ in NaTlO2).
 *   q-path : the high-symmetry segments of the end-member with the fewer
 *            segments, kept where the other end-member has a segment on the
 *            same line (in either direction) and cut to the stretch both cover.
 *            For zinc blende and wurtzite the segments coincide (the hexagonal
 *            special points do not move with c/a); for chalcopyrite the zone-
 *            boundary points Sigma, Y, ... move with c/a, so segments ending
 *            there are cut at the nearer of the two and segments whose line
 *            itself moves are dropped. Ternary wurtzites with force constants
 *            are evaluated directly on the standard path of the canonical cell,
 *            since the files run their paths along different but equivalent
 *            axes. Both end-members are sampled at the same
 *            reduced q. All files of a family share the same primitive vectors,
 *            so reduced coordinates compare directly.
 *
 * Swapping the end-members and x -> 1-x gives the same D_alloy, since the path
 * choice, the sampling and the lattice interpolation are all symmetric.
 */

import { buildDynamicalMatrixBlocks } from './dynamicalmatrix.js';
import { getGonzeReciprocalCorrection, getGonzeRealSpaceRemainder, getGonzeDipoleDipole } from './gonze.js';
import { solveComplexHermitianWithEigenWasm } from './eigenwasm.js';
import * as atomic_data from './atomic_data.js';
import * as utils from './utils.js';
import * as mat from './mat.js';

const QTOL = 1e-5;
const THZ_TO_CM1 = 33.35641;
// phonopy units: sqrt(eV/Angstrom^2/amu) -> THz, and Hartree*Bohr in eV*Angstrom
// (the NAC unit conversion stored in every PhononDB file)
const PHONOPY_THZ = 15.633302300230191;
const PHONOPY_NAC_UNIT = 14.39965172592227;

function sameQ(a, b) {
    return Math.abs(a[0] - b[0]) < QTOL &&
           Math.abs(a[1] - b[1]) < QTOL &&
           Math.abs(a[2] - b[2]) < QTOL;
}

function isGamma(q) {
    return Math.abs(q[0]) < QTOL && Math.abs(q[1]) < QTOL && Math.abs(q[2]) < QTOL;
}

function lerpQ(qA, qB, t) {
    return [
        qA[0] + t * (qB[0] - qA[0]),
        qA[1] + t * (qB[1] - qA[1]),
        qA[2] + t * (qB[2] - qA[2]),
    ];
}

function zeroMatrix(size) {
    let real = [];
    let imag = [];
    for (let i = 0; i < size; i++) {
        real.push(new Array(size).fill(0));
        imag.push(new Array(size).fill(0));
    }
    return { real: real, imag: imag };
}

function scaleMatrix(matrix, factor) {
    let size = matrix.real.length;
    let out = zeroMatrix(size);
    for (let i = 0; i < size; i++) {
        for (let j = 0; j < size; j++) {
            out.real[i][j] = matrix.real[i][j] * factor;
            out.imag[i][j] = matrix.imag[i][j] * factor;
        }
    }
    return out;
}

function mixForceConstantMatrices(d1, masses1, d2, masses2, massesMixed, x) {
    /*
    D_x = M_x^-1/2 [(1-x) M_1^1/2 D_1 M_1^1/2 + x M_2^1/2 D_2 M_2^1/2] M_x^-1/2
    with the masses given per site, each covering 3 rows/columns
    */
    let size = d1.real.length;
    let out = zeroMatrix(size);
    for (let i = 0; i < size; i++) {
        let si = Math.floor(i / 3);
        for (let j = 0; j < size; j++) {
            let sj = Math.floor(j / 3);
            let w1 = (1 - x) * Math.sqrt(masses1[si] * masses1[sj]);
            let w2 = x * Math.sqrt(masses2[si] * masses2[sj]);
            let norm = 1 / Math.sqrt(massesMixed[si] * massesMixed[sj]);
            out.real[i][j] = (w1 * d1.real[i][j] + w2 * d2.real[i][j]) * norm;
            out.imag[i][j] = (w1 * d1.imag[i][j] + w2 * d2.imag[i][j]) * norm;
        }
    }
    return out;
}

function lerpMatrix(a, b, x) {
    let size = a.real.length;
    let out = zeroMatrix(size);
    for (let i = 0; i < size; i++) {
        for (let j = 0; j < size; j++) {
            out.real[i][j] = (1 - x) * a.real[i][j] + x * b.real[i][j];
            out.imag[i][j] = (1 - x) * a.imag[i][j] + x * b.imag[i][j];
        }
    }
    return out;
}

function normalizeLabel(label) {
    let text = String(label).replace(/\$/g, '').trim();
    if (text === '\\Gamma' || text.toUpperCase() === 'GAMMA' || text === 'Γ' || text === 'G') {
        return 'GAMMA';
    }
    return text;
}

function parseComplex(value) {
    if (typeof value === 'number') {
        return [value, 0];
    }
    if (Array.isArray(value)) {
        return [Number(value[0]), Number(value[1])];
    }
    let text = String(value).trim().replace(/^\(/, '').replace(/\)$/, '').replace(/\s+/g, '');
    if (!text.endsWith('j')) {
        return [Number(text), 0];
    }
    let splitIndex = -1;
    for (let i = 1; i < text.length; i++) {
        if ((text[i] === '+' || text[i] === '-') && text[i - 1] !== 'e' && text[i - 1] !== 'E') {
            splitIndex = i;
        }
    }
    if (splitIndex === -1) {
        return [0, Number(text.slice(0, -1))];
    }
    return [Number(text.slice(0, splitIndex)), Number(text.slice(splitIndex, -1))];
}

// formal charge of the usual anions (the magnitude)
const ANION_CHARGE = {
    O: 2, S: 2, Se: 2, Te: 2,
    N: 3, P: 3, As: 3,
    F: 1, Cl: 1, Br: 1, I: 1, H: 1,
};

function getCationSymbol(name, atomTypes) {
    // in a binary the species that is not a usual anion (PbTe is listed as TePb)
    let species = atomTypes.filter((type, i) => atomTypes.indexOf(type) === i);
    if (species.length === 2) {
        let anions = species.filter((type) => ANION_CHARGE[type]);
        if (anions.length === 1) {
            return species.find((type) => type !== anions[0]);
        }
    }
    let match = /^([A-Z][a-z]?)/.exec(name || '');
    if (match && atomTypes.indexOf(match[1]) !== -1) {
        return match[1];
    }
    return atomTypes[0];
}

function isFccLattice(lattice) {
    // three primitive vectors of equal length at 60 degrees
    let length = lattice.map((row) => Math.sqrt(mat.vec_dot(row, row)));
    return [[0, 1], [0, 2], [1, 2]].every(([i, j]) =>
        Math.abs(length[i] - length[j]) < 1e-3 * length[0] &&
        Math.abs(mat.vec_dot(lattice[i], lattice[j]) / (length[i] * length[j]) - 0.5) < 1e-3
    );
}

function getCanonicalSites(atomTypes, positionsRed, cationSymbol, lattice = null) {
    /*
    returns the original atom indices in [cation, anion] order, whether the
    structure has to be inverted to put the anion at +(1/4,1/4,1/4) and the
    family: zinc blende, or rocksalt with the anion at (1/2,1/2,1/2)
    */
    let isElemental = atomTypes[0] === atomTypes[1];
    let cation = isElemental ? 0 : atomTypes.indexOf(cationSymbol);
    let anion = 1 - cation;

    let offset = [0, 1, 2].map((i) => {
        let d = positionsRed[anion][i] - positionsRed[cation][i];
        return d - Math.floor(d + 1e-6);
    });
    if (!isElemental && offset.every((d) => Math.abs(d - 0.5) < 1e-3) && (!lattice || isFccLattice(lattice))) {
        // rocksalt is centrosymmetric: no inversion needed
        return { order: [cation, anion], invert: false, family: 'rocksalt' };
    }
    let isPlus = offset.every((d) => Math.abs(d - 0.25) < 1e-3);
    let isMinus = offset.every((d) => Math.abs(d - 0.75) < 1e-3);
    if (!isPlus && !isMinus) {
        throw new Error('Not a 2-atom zinc blende / diamond or rocksalt structure: anion offset ' + JSON.stringify(offset));
    }
    if (isMinus && isElemental) {
        // both sites are equivalent, relabel instead of inverting
        return { order: [anion, cation], invert: false, family: 'zincblende' };
    }
    return { order: [cation, anion], invert: isMinus, family: 'zincblende' };
}

// species on the tetrahedral sites of a fluorite (CaF2, CoSi2); on the fcc site
// they make an antifluorite (Li2O, Mg2Si, Be2C)
const ANION_LIKE = new Set(Object.keys(ANION_CHARGE).concat(['C', 'Si', 'Ge', 'Sn', 'Pb', 'Sb', 'Bi']));

function isFluoriteCell(endpoint) {
    try {
        getFluoriteSites(endpoint.atomTypes, endpoint.positionsRed, endpoint.lattice);
        return true;
    } catch (error) {
        return false;
    }
}

function getFluoriteSites(atomTypes, positionsRed, lattice) {
    /*
    a 3-atom fluorite or antifluorite cell: [fcc site, tetrahedral site at
    +(1/4,1/4,1/4), the one at -(1/4,1/4,1/4)]. The structure is
    centrosymmetric, labelling the two tetrahedral sites by their offset needs no
    inversion. Fluorite if the tetrahedral species is an anion, else antifluorite.
    */
    if (atomTypes.length !== 3 || !isFccLattice(lattice)) {
        throw new Error('Not a 3-atom fluorite cell');
    }
    for (let center = 0; center < 3; center++) {
        let others = [0, 1, 2].filter((i) => i !== center);
        if (atomTypes[others[0]] !== atomTypes[others[1]] || atomTypes[center] === atomTypes[others[0]]) {
            continue;
        }
        let offset = (i) => positionsRed[i].map((v, k) => {
            let d = v - positionsRed[center][k];
            return d - Math.floor(d + 1e-6);
        });
        let isAt = (d, value) => d.every((v) => Math.abs(v - value) < 1e-3);
        let plus = others.find((i) => isAt(offset(i), 0.25));
        let minus = others.find((i) => isAt(offset(i), 0.75));
        if (plus === undefined || minus === undefined) {
            continue;
        }
        return {
            order: [center, plus, minus],
            invert: false,
            family: ANION_LIKE.has(atomTypes[plus]) ? 'fluorite' : 'antifluorite',
        };
    }
    throw new Error('Not a 3-atom fluorite / antifluorite structure');
}

// usual cation valence, used to put the lower-valence cation of a chalcopyrite
// (Cu in CuGaS2, Zn in ZnGeP2) on the 4a site and the other one on 4b
const CATION_VALENCE = {
    Li: 1, Na: 1, K: 1, Cu: 1, Ag: 1, Au: 1,
    Be: 2, Mg: 2, Ca: 2, Zn: 2, Cd: 2, Hg: 2,
    B: 3, Al: 3, Ga: 3, In: 3, Tl: 3,
    C: 4, Si: 4, Ge: 4, Sn: 4, Pb: 4,
    Rb: 1, Cs: 1, Sr: 2, Ba: 2, Sc: 3, Y: 3, La: 3, Pr: 3, Nd: 3, Sm: 3, Tb: 3, Dy: 3,
    Ho: 3, Er: 3, Tm: 3, Lu: 3, Cr: 3, Fe: 3, Rh: 3, Bi: 3, Ti: 4, Zr: 4, Hf: 4, Pd: 4,
    Nb: 5, Ta: 5, P: 5, Sb: 5, U: 6, Re: 7,
};

// chalcopyrite I-42d (No. 122), conventional fractional coordinates: A on 4a,
// B on 4b, anions on 8d (x, 1/4, 1/8). Each entry lists the two fixed
// coordinates and the axis/sign of the free parameter u; with the I centring
// these are the 8 atoms of the primitive cell, in canonical order.
const CHALCOPYRITE_SITES = [
    { role: 'A', fixed: [0, 0, 0] },
    { role: 'A', fixed: [0, 0.5, 0.25] },
    { role: 'B', fixed: [0, 0, 0.5] },
    { role: 'B', fixed: [0, 0.5, 0.75] },
    { role: 'X', fixed: [null, 0.25, 0.125], uAxis: 0, uSign: 1 },
    { role: 'X', fixed: [null, 0.75, 0.125], uAxis: 0, uSign: -1 },
    { role: 'X', fixed: [0.75, null, 0.875], uAxis: 1, uSign: 1 },
    { role: 'X', fixed: [0.25, null, 0.875], uAxis: 1, uSign: -1 },
];

function wrapFraction(value) {
    return value - Math.round(value);
}

function getBodyCenteredTetragonalFrame(lattice) {
    /*
    conventional axes of a body-centred tetragonal primitive cell stored as
    (-a/2, a/2, c/2), (a/2, -a/2, c/2), (a/2, a/2, -c/2) in some cartesian frame:
    a2 + a3 = a x, a1 + a3 = a y, a1 + a2 = c z
    */
    let sum = (u, v) => [u[0] + v[0], u[1] + v[1], u[2] + v[2]];
    let ex = sum(lattice[1], lattice[2]);
    let ey = sum(lattice[0], lattice[2]);
    let ez = sum(lattice[0], lattice[1]);
    let a = Math.sqrt(mat.vec_dot(ex, ex));
    let b = Math.sqrt(mat.vec_dot(ey, ey));
    let c = Math.sqrt(mat.vec_dot(ez, ez));
    let orthogonal = Math.abs(mat.vec_dot(ex, ey)) + Math.abs(mat.vec_dot(ex, ez)) + Math.abs(mat.vec_dot(ey, ez));
    if (Math.abs(a - b) > 1e-3 * a || orthogonal > 1e-3 * a * c) {
        throw new Error('Unsupported setting: not the standard body-centred tetragonal primitive cell');
    }
    return { a: a, c: c, axes: [mat.vec_scale(ex, 1 / a), mat.vec_scale(ey, 1 / a), mat.vec_scale(ez, 1 / c)] };
}

export function chalcopyriteLattice(a, c) {
    return [[-a / 2, a / 2, c / 2], [a / 2, -a / 2, c / 2], [a / 2, a / 2, -c / 2]];
}

export function chalcopyritePositions(a, c, u) {
    /*
    reduced coordinates of the canonical 8 sites in chalcopyriteLattice(a, c)
    */
    let inverse = mat.matrix_inverse(chalcopyriteLattice(a, c));
    return CHALCOPYRITE_SITES.map((site) => {
        let f = site.fixed.slice();
        if (site.role === 'X') {
            f[site.uAxis] = site.uSign * u;
        }
        let cartesian = [f[0] * a, f[1] * a, f[2] * c];
        return [0, 1, 2].map((j) => cartesian[0] * inverse[0][j] + cartesian[1] * inverse[1][j] + cartesian[2] * inverse[2][j]);
    });
}

function getChalcopyriteRoles(atomTypes) {
    /*
    species -> role: the species with 4 atoms is the anion, the cation with the
    lower valence (then lower atomic number) sits on 4a
    */
    let counts = {};
    atomTypes.forEach((type) => { counts[type] = (counts[type] || 0) + 1; });
    let species = Object.keys(counts);
    let anions = species.filter((type) => counts[type] === 4);
    let cations = species.filter((type) => counts[type] === 2);
    if (atomTypes.length !== 8 || anions.length !== 1 || cations.length !== 2) {
        throw new Error('Not an A B X2 chalcopyrite cell: ' + JSON.stringify(counts));
    }
    let key = (type) => [CATION_VALENCE[type] || 9, atomic_data.atomic_number[type]];
    cations.sort((p, q) => key(p)[0] - key(q)[0] || key(p)[1] - key(q)[1]);
    let roles = {};
    roles[cations[0]] = 'A';
    roles[cations[1]] = 'B';
    roles[anions[0]] = 'X';
    return roles;
}

function matchChalcopyriteSites(atomTypes, positionsRed, lattice) {
    /*
    map the atoms of a chalcopyrite cell onto CHALCOPYRITE_SITES. The cell may be
    stored with any origin, or as the inversion image of the template (every
    other setting of the same lattice differs from these by a symmetry of the
    crystal); the inversion image is handled like in zinc blende, D -> D^*.
    Returns the file atom index of every canonical site, the inversion flag, the
    rotation to the conventional axes and the structural parameters a, c, u.
    */
    let roles = getChalcopyriteRoles(atomTypes);
    let frame = getBodyCenteredTetragonalFrame(lattice);
    let cartesian = utils.red_car_list(positionsRed, lattice);
    let toConventional = (d) => [
        mat.vec_dot(d, frame.axes[0]) / frame.a,
        mat.vec_dot(d, frame.axes[1]) / frame.a,
        mat.vec_dot(d, frame.axes[2]) / frame.c,
    ];
    const tol = 2e-2;
    let fixedMatch = (f, fixed) => fixed.every((v, i) => v === null || Math.abs(wrapFraction(f[i] - v)) < tol);

    for (let sign of [1, -1]) {
        for (let origin = 0; origin < atomTypes.length; origin++) {
            if (roles[atomTypes[origin]] !== 'A') {
                continue;
            }
            let order = new Array(CHALCOPYRITE_SITES.length).fill(-1);
            let uValues = [];
            let ok = true;
            for (let i = 0; i < atomTypes.length && ok; i++) {
                let d = [0, 1, 2].map((k) => sign * (cartesian[i][k] - cartesian[origin][k]));
                let f = toConventional(d);
                let found = -1;
                for (let s = 0; s < CHALCOPYRITE_SITES.length && found < 0; s++) {
                    let site = CHALCOPYRITE_SITES[s];
                    if (site.role !== roles[atomTypes[i]] || order[s] >= 0) {
                        continue;
                    }
                    // the I centring (1/2, 1/2, 1/2) maps a site onto its partner
                    for (let shift of [0, 0.5]) {
                        let g = f.map((v) => v - shift);
                        if (!fixedMatch(g, site.fixed)) {
                            continue;
                        }
                        if (site.role === 'X') {
                            let u = wrapFraction(site.uSign * g[site.uAxis]);
                            if (u < 0.15 || u > 0.35) {
                                continue;
                            }
                            uValues.push(u);
                        }
                        found = s;
                        break;
                    }
                }
                if (found < 0) {
                    ok = false;
                } else {
                    order[found] = i;
                }
            }
            let uMean = uValues.reduce((p, v) => p + v, 0) / (uValues.length || 1);
            if (ok && uValues.every((u) => Math.abs(u - uMean) < 1e-2)) {
                let identity = frame.axes.every((axis, k) => Math.abs(axis[k] - 1) < 1e-6);
                return {
                    order: order,
                    invert: sign < 0,
                    rotation: identity ? null : frame.axes,
                    structure: { a: frame.a, c: frame.c, u: uMean },
                };
            }
        }
    }
    throw new Error('Could not match the cell onto the chalcopyrite (I-42d) template');
}

function getChalcopyriteSpecialPoints(a, c) {
    /*
    body-centred tetragonal (c > a) special points of Setyawan & Curtarolo, in
    reduced coordinates of chalcopyriteLattice(a, c). Plain unicode names: the
    chart only converts the first TeX name of a merged label such as "Σ1|Σ"
    */
    let eta = (1 + (a * a) / (c * c)) / 4;
    let zeta = (a * a) / (2 * c * c);
    return [
        { label: 'GAMMA', q: [0, 0, 0] },
        { label: 'X', q: [0, 0, 0.5] },
        { label: 'N', q: [0, 0.5, 0] },
        { label: 'P', q: [0.25, 0.25, 0.25] },
        { label: 'Z', q: [0.5, 0.5, -0.5] },
        { label: 'Σ', q: [-eta, eta, eta] },
        { label: 'Σ1', q: [eta, 1 - eta, -eta] },
        { label: 'Y', q: [-zeta, zeta, 0.5] },
        { label: 'Y1', q: [0.5, 0.5, -zeta] },
    ];
}

// wurtzite P6_3mc (No. 186) in the hexagonal cell wurtziteLattice(a, c): cations
// at (1/3, 2/3, 0) and (2/3, 1/3, 1/2), each anion u c above its cation. Listed
// relative to the first cation, in canonical order.
const WURTZITE_SITES = [
    { role: 'C', red: [0, 0, 0] },
    { role: 'C', red: [1 / 3, -1 / 3, 0.5] },
    { role: 'X', red: [0, 0, 0] },
    { role: 'X', red: [1 / 3, -1 / 3, 0.5] },
];

export function wurtziteLattice(a, c) {
    return [[a, 0, 0], [-a / 2, (a * Math.sqrt(3)) / 2, 0], [0, 0, c]];
}

export function wurtzitePositions(u) {
    /*
    reduced coordinates of the canonical 4 sites in wurtziteLattice(a, c)
    */
    return WURTZITE_SITES.map((site) => {
        let f = [site.red[0] + 1 / 3, site.red[1] + 2 / 3, site.red[2] + (site.role === 'X' ? u : 0)];
        return f.map((v) => v - Math.floor(v + 1e-9));
    });
}

function getHexagonalFrame(lattice) {
    /*
    rotation onto the axes of wurtziteLattice: a1 along x, a3 along z and a2 at
    120 degrees from a1 in a right-handed frame (the setting of PhononDB and of
    the Materials Project, up to a rotation about z)
    */
    let a = Math.sqrt(mat.vec_dot(lattice[0], lattice[0]));
    let c = Math.sqrt(mat.vec_dot(lattice[2], lattice[2]));
    let ex = mat.vec_scale(lattice[0], 1 / a);
    let ez = mat.vec_scale(lattice[2], 1 / c);
    let ey = [
        ez[1] * ex[2] - ez[2] * ex[1],
        ez[2] * ex[0] - ez[0] * ex[2],
        ez[0] * ex[1] - ez[1] * ex[0],
    ];
    let expected = [0, 1, 2].map((k) => a * (-0.5 * ex[k] + (Math.sqrt(3) / 2) * ey[k]));
    let mismatch = Math.max.apply(null, [0, 1, 2].map((k) => Math.abs(lattice[1][k] - expected[k])));
    if (mismatch > 1e-3 * a || Math.abs(mat.vec_dot(ex, ez)) > 1e-4) {
        throw new Error('Unsupported setting: not the standard hexagonal cell with a1 and a2 at 120 degrees');
    }
    return { a: a, c: c, axes: [ex, ey, ez] };
}

function matchWurtziteSites(atomTypes, positionsRed, lattice, cationSymbol) {
    /*
    map the atoms of a 4-atom wurtzite cell onto WURTZITE_SITES, with any origin
    and either polarity; a cell of the opposite polarity is the inversion image of
    the template, handled like in zinc blende by D -> D^*. The other cation as
    origin covers the other stacking in the plane. Returns the file atom index of
    every canonical site, the inversion flag, the rotation to the axes of
    wurtziteLattice and the structural parameters a, c, u.
    */
    let elemental = atomTypes.every((type) => type === atomTypes[0]);
    let role = (type) => (elemental ? null : type === cationSymbol ? 'C' : 'X');
    if (atomTypes.length !== 4 || (!elemental && atomTypes.filter((type) => type === cationSymbol).length !== 2)) {
        throw new Error('Not a 4-atom wurtzite cell: ' + atomTypes.join(' '));
    }
    let frame = getHexagonalFrame(lattice);
    const tol = 2e-2;

    for (let sign of [1, -1]) {
        for (let origin = 0; origin < 4; origin++) {
            if (role(atomTypes[origin]) === 'X') {
                continue;
            }
            // the rotation keeps the lattice vectors, so reduced coordinates carry over
            let order = new Array(4).fill(-1);
            let uValues = [];
            for (let i = 0; i < 4; i++) {
                let f = [0, 1, 2].map((k) => sign * (positionsRed[i][k] - positionsRed[origin][k]));
                for (let s = 0; s < 4; s++) {
                    let site = WURTZITE_SITES[s];
                    let siteRole = role(atomTypes[i]);
                    if (order[s] >= 0 || (siteRole && siteRole !== site.role)) {
                        continue;
                    }
                    let d = f.map((v, k) => wrapFraction(v - site.red[k]));
                    if (Math.abs(d[0]) > tol || Math.abs(d[1]) > tol) {
                        continue;
                    }
                    if (site.role === 'C' ? Math.abs(d[2]) > tol : d[2] < 0.3 || d[2] > 0.45) {
                        continue;
                    }
                    if (site.role === 'X') {
                        uValues.push(d[2]);
                    }
                    order[s] = i;
                    break;
                }
            }
            if (order.every((i) => i >= 0) && Math.abs(uValues[0] - uValues[1]) < 1e-2) {
                let identity = frame.axes.every((axis, k) => Math.abs(axis[k] - 1) < 1e-6);
                return {
                    order: order,
                    invert: sign < 0,
                    rotation: identity ? null : frame.axes,
                    structure: { a: frame.a, c: frame.c, u: (uValues[0] + uValues[1]) / 2 },
                };
            }
        }
    }
    throw new Error('Could not match the cell onto the wurtzite (P6_3mc) template');
}

// hexagonal special points, in reduced coordinates of wurtziteLattice(a, c);
// unlike the body-centred tetragonal ones they do not depend on c/a
const HEXAGONAL_SPECIAL_POINTS = [
    { label: 'GAMMA', q: [0, 0, 0] },
    { label: 'M', q: [0.5, 0, 0] },
    { label: 'K', q: [1 / 3, 1 / 3, 0] },
    { label: 'A', q: [0, 0, 0.5] },
    { label: 'L', q: [0.5, 0, 0.5] },
    { label: 'H', q: [1 / 3, 1 / 3, 0.5] },
];

// ---------------------------------------------------------------------------
// ternary wurtzite derivatives: cation orderings on a supercell of wurtzite
// ---------------------------------------------------------------------------

// simple orthorhombic special points, reduced coordinates of an orthogonal cell
const ORTHORHOMBIC_SPECIAL_POINTS = [
    { label: 'GAMMA', q: [0, 0, 0] },
    { label: 'X', q: [0.5, 0, 0] },
    { label: 'Y', q: [0, 0.5, 0] },
    { label: 'Z', q: [0, 0, 0.5] },
    { label: 'S', q: [0.5, 0.5, 0] },
    { label: 'U', q: [0.5, 0, 0.5] },
    { label: 'T', q: [0, 0.5, 0.5] },
    { label: 'R', q: [0.5, 0.5, 0.5] },
];

// seekpath path of a simple orthorhombic cell
const ORTHORHOMBIC_PATH = [
    ['GAMMA', 'X'], ['X', 'S'], ['S', 'Y'], ['Y', 'GAMMA'], ['GAMMA', 'Z'], ['Z', 'U'],
    ['U', 'R'], ['R', 'T'], ['T', 'Z'], ['Y', 'T'], ['U', 'X'], ['S', 'R'],
];

// Each structure is a supercell of the wurtzite cell wurtziteLattice(a, c):
// `supercell` holds its lattice vectors in the wurtzite basis, and every site sits
// on wurtzite sublattice `sub` (WURTZITE_SITES) translated by the wurtzite lattice
// vector `n`. Roles: A = the cation listed first (more atoms, then lower valence),
// B = the other cation, X = anion. Derived from PhononDB mp-2979 (ZnGeN2),
// mp-3934 (Cu3PS4) and mp-9248 (Cu2SiS3) with analyzeWurtziteDerivative.
const superstructureSites = (list) => list.map((site) => ({ role: site[0], sub: site[1], n: site.slice(2) }));
export const WURTZITE_SUPERSTRUCTURES = {
    'wurtzite-ABX2': {
        parent: 'wurtzite',
        title: 'beta-NaFeO2 type A B X2 (Pna2_1)',
        supercell: [[2, 0, 0], [1, 2, 0], [0, 0, 1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 0, 2, 1, 0], ['A', 1, 0, 0, 0], ['A', 1, 1, 1, 0],
            ['B', 0, 1, 0, 0], ['B', 0, 1, 1, 0], ['B', 1, 1, 0, 0], ['B', 1, 2, 1, 0],
            ['X', 2, 0, 0, 0], ['X', 2, 1, 0, 0], ['X', 2, 1, 1, 0], ['X', 2, 2, 1, 0],
            ['X', 3, 0, 0, 0], ['X', 3, 1, 0, 0], ['X', 3, 1, 1, 0], ['X', 3, 2, 1, 0],
        ]),
        specialPoints: ORTHORHOMBIC_SPECIAL_POINTS,
        path: ORTHORHOMBIC_PATH,
    },
    'wurtzite-A3BX4': {
        parent: 'wurtzite',
        title: 'enargite type A3 B X4 (Pmn2_1)',
        supercell: [[2, 0, 0], [1, 2, 0], [0, 0, 1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 0, 1, 0, 0], ['A', 0, 2, 1, 0], ['A', 1, 0, 0, 0], ['A', 1, 1, 1, 0], ['A', 1, 2, 1, 0],
            ['B', 0, 1, 1, 0], ['B', 1, 1, 0, 0],
            ['X', 2, 0, 0, 0], ['X', 2, 1, 0, 0], ['X', 2, 1, 1, 0], ['X', 2, 2, 1, 0],
            ['X', 3, 0, 0, 0], ['X', 3, 1, 0, 0], ['X', 3, 1, 1, 0], ['X', 3, 2, 1, 0],
        ]),
        specialPoints: ORTHORHOMBIC_SPECIAL_POINTS,
        path: ORTHORHOMBIC_PATH,
    },
    'wurtzite-A2BX3': {
        parent: 'wurtzite',
        title: 'A2 B X3 (Cmc2_1)',
        // C-centred: the primitive cell is not orthogonal and its special points
        // depend on the cell shape, the labels of the file are kept
        supercell: [[1, -1, 0], [2, 1, 0], [0, 0, 1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 0, 1, 0, 0], ['A', 1, 1, 0, 0], ['A', 1, 2, 0, 0],
            ['B', 0, 2, 0, 0], ['B', 1, 0, 0, 0],
            ['X', 2, 0, 0, 0], ['X', 2, 1, 0, 0], ['X', 2, 2, 0, 0],
            ['X', 3, 0, 0, 0], ['X', 3, 1, 0, 0], ['X', 3, 2, 0, 0],
        ]),
        specialPoints: null,
        path: null,
    },
};

// cation orderings on a supercell of rocksalt (rocksaltLattice(a), sub 0 =
// cation, 1 = anion), same encoding as WURTZITE_SUPERSTRUCTURES. Derived with
// analyzeSuperstructure from PhononDB mp-15786 (KHoS2), mp-676360 (NaNdS2),
// mp-11151 (NaPrTe2), mp-978844 (SrCa3O4), mp-755287
// (Ca3CdO4), mp-4156 (Li2ZrO3), mp-7219 (Na2ZrSe3) and mp-3540 (Li2SnO3).
// NaBiS2 (Fdd2) has the cation order of the gamma-LiFeO2 type, distorted. The
// special points of these cells depend on their shape: the path of the file is
// kept. referencePoints: where the Materials Project files of the family put the
// special points fixed by symmetry (from mp-29678, mp-33526, mp-755287, mp-4156
// and mp-7219), see referencePointScore.
export const ROCKSALT_SUPERSTRUCTURES = {
    'rocksalt-ABX2-R3m': {
        parent: 'rocksalt',
        title: 'alpha-NaFeO2 type A B X2 (R-3m)',
        label: 'R-3m',
        supercell: [[0, -1, 2], [-1, 0, 2], [-1, -1, 2]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['B', 0, -1, -1, 3], ['X', 1, 0, 0, 0], ['X', 1, -1, -1, 3],
        ]),
        referencePoints: { L: [0.5, 0, 0], Z: [0.5, 0.5, 0.5], F: [0.5, 0.5, 0] },
        specialPoints: null,
        path: null,
    },
    'rocksalt-ABX2-I41amd': {
        parent: 'rocksalt',
        title: 'gamma-LiFeO2 type A B X2 (I4_1/amd)',
        label: 'I4_1/amd',
        supercell: [[0, -2, 1], [-2, 0, 1], [1, 1, -2]],
        sites: superstructureSites([
            ['A', 0, 0, -1, 0], ['A', 0, 0, 0, 0], ['B', 0, -1, 0, 0], ['B', 0, -1, -1, 1], ['X', 1, 0, -1, 0], ['X', 1, -1, 0, 0],
            ['X', 1, 0, 0, 0], ['X', 1, -1, -1, 1],
        ]),
        referencePoints: { X: [0.5, 0.5, 0], Z: [0.5, 0.5, 0.5], N: [0, 0.5, 0], P: [0.25, 0.25, 0.25] },
        specialPoints: null,
        path: null,
    },
    'rocksalt-ABX2-Fd3m': {
        parent: 'rocksalt',
        title: 'A B X2 (Fd-3m)',
        label: 'Fd-3m',
        supercell: [[0, 0, 2], [-2, 0, 2], [0, -2, 2]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 0, 0, 0, 1], ['A', 0, 0, -1, 2], ['A', 0, -1, 0, 2], ['B', 0, 0, -1, 1], ['B', 0, -1, 0, 1],
            ['B', 0, -1, -1, 2], ['B', 0, -1, -1, 3], ['X', 1, 0, 0, 0], ['X', 1, 0, -1, 1], ['X', 1, -1, 0, 1], ['X', 1, 0, 0, 1],
            ['X', 1, -1, -1, 2], ['X', 1, 0, -1, 2], ['X', 1, -1, 0, 2], ['X', 1, -1, -1, 3],
        ]),
        specialPoints: null,
        path: null,
    },
    'rocksalt-A3BX4-Pm3m': {
        parent: 'rocksalt',
        title: 'A3 B X4 (Pm-3m)',
        label: 'Pm-3m',
        supercell: [[1, -1, -1], [-1, 1, -1], [1, 1, -1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, -1], ['A', 0, 1, 0, -1], ['A', 0, 0, 0, 0], ['B', 0, 0, 1, -1], ['X', 1, 0, 0, -1], ['X', 1, 1, 0, -1],
            ['X', 1, 0, 1, -1], ['X', 1, 0, 0, 0],
        ]),
        specialPoints: null,
        path: null,
    },
    'rocksalt-A3BX4-Cmmm': {
        parent: 'rocksalt',
        title: 'A3 B X4 (Cmmm)',
        label: 'Cmmm',
        supercell: [[0, -1, -1], [-2, 1, 1], [0, 1, -1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, -1], ['A', 0, 0, 0, 0], ['A', 0, -1, 1, 0], ['B', 0, -1, 0, 0], ['X', 1, 0, 0, -1], ['X', 1, -1, 0, 0],
            ['X', 1, 0, 0, 0], ['X', 1, -1, 1, 0],
        ]),
        referencePoints: { S: [0, 0.5, 0], R: [0, 0.5, 0.5], Z: [0, 0, 0.5], Y: [0.5, 0.5, 0], T: [0.5, 0.5, 0.5] },
        specialPoints: null,
        path: null,
    },
    'rocksalt-A2BX3-C2c': {
        parent: 'rocksalt',
        title: 'Li2ZrO3 type A2 B X3 (C2/c)',
        label: 'C2/c',
        supercell: [[1, 0, -2], [1, 0, 1], [0, -2, 1]],
        sites: superstructureSites([
            ['A', 0, 1, -1, -1], ['A', 0, 1, 0, -1], ['A', 0, 1, -1, 0], ['A', 0, 0, 0, 0], ['B', 0, 1, 0, 0], ['B', 0, 1, -1, 1],
            ['X', 1, 1, -1, -1], ['X', 1, 1, 0, -1], ['X', 1, 1, -1, 0], ['X', 1, 0, 0, 0], ['X', 1, 1, 0, 0], ['X', 1, 1, -1, 1],
        ]),
        referencePoints: { Z: [0, 0, 0.5], I: [0.5, 0.5, 0.5], X: [0.5, 0.5, 0], N: [0, 0.5, 0], M: [0, 0.5, 0.5] },
        specialPoints: null,
        path: null,
    },
    'rocksalt-A2BX3-C2m': {
        parent: 'rocksalt',
        title: 'Li2MnO3 type A2 B X3 (C2/m)',
        label: 'C2/m',
        supercell: [[-1, -1, 0], [2, -1, 0], [-1, 0, 2]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 0, -1, -1, 1], ['A', 0, 0, -1, 1], ['A', 0, 1, -1, 1], ['B', 0, 0, -1, 0], ['B', 0, 1, -1, 0],
            ['X', 1, 0, -1, 0], ['X', 1, 1, -1, 0], ['X', 1, 0, 0, 0], ['X', 1, -1, -1, 1], ['X', 1, 0, -1, 1], ['X', 1, 1, -1, 1],
        ]),
        referencePoints: { Z: [0, 0, 0.5], I: [0.5, 0.5, 0.5], X: [0.5, 0.5, 0], N: [0, 0.5, 0], M: [0, 0.5, 0.5] },
        specialPoints: null,
        path: null,
    },
    'rocksalt-A2BX3-C2c24': {
        parent: 'rocksalt',
        title: 'Li2SnO3 type A2 B X3 (C2/c)',
        label: 'C2/c',
        supercell: [[0, -1, 2], [0, 2, -1], [-4, 1, 1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 0, 0, 1, 0], ['A', 0, -2, 2, 0], ['A', 0, -1, 0, 1], ['A', 0, 0, 0, 1], ['A', 0, -3, 1, 1],
            ['A', 0, -2, 1, 1], ['A', 0, -2, 0, 2], ['B', 0, -1, 1, 0], ['B', 0, -1, 1, 1], ['B', 0, -3, 2, 1], ['B', 0, -3, 1, 2],
            ['X', 1, 0, 0, 0], ['X', 1, -1, 1, 0], ['X', 1, 0, 1, 0], ['X', 1, -2, 2, 0], ['X', 1, -1, 0, 1], ['X', 1, 0, 0, 1],
            ['X', 1, -3, 1, 1], ['X', 1, -2, 1, 1], ['X', 1, -1, 1, 1], ['X', 1, -3, 2, 1], ['X', 1, -2, 0, 2], ['X', 1, -3, 1, 2],
        ]),
        specialPoints: null,
        path: null,
    },
};

// cation orderings on a supercell of fluorite (sub 0: fcc site, 1 and 2: the
// tetrahedral sites at +-(1/4,1/4,1/4)): on the fcc site in fluorite (anion X on
// the tetrahedral sites), on the tetrahedral sites in antifluorite (X on the fcc
// site). Derived with analyzeSuperstructure from PhononDB mp-37906 (LiMgN),
// mp-685991 (Li5SiP3), mp-5104 (CaZnF4), mp-756235 (Ba3SrI8), mp-3960 (CaUO4),
// mp-752798 (CeTh2O6), mp-676074 (Na2ThF6) and mp-581857 (KBiF4). Distorted
// variants share the cation order: fergusonite (C2/c) and NaYF4 (I4_1/amd) that
// of scheelite, Zr3GeO8 (I-42m) that of Ba3SrI8, Li5TiP3 (C2/c) that of Li5SiP3.
export const FLUORITE_SUPERSTRUCTURES = {
    'antifluorite-ABX-F43m': {
        parent: 'antifluorite',
        title: 'MgAgAs type A B X (F-43m)',
        label: 'F-43m',
        supercell: [[0, 0, 1], [-1, 0, 1], [0, -1, 1]],
        sites: superstructureSites([
            ['A', 1, 0, 0, 0], ['B', 2, 0, 0, 0], ['X', 0, 0, 0, 0],
        ]),
        referencePoints: { X: [0.5, 0, 0.5], W: [0.5, 0.25, -0.25], L: [0.5, 0.5, 0.5] },
        specialPoints: null,
        path: null,
    },
    'antifluorite-A5BX3-Cc': {
        parent: 'antifluorite',
        title: 'A5 B X3 (Cc)',
        label: 'Cc',
        supercell: [[2, -1, -1], [1, 0, -2], [0, 2, 0]],
        sites: superstructureSites([
            ['A', 1, 2, 0, -2], ['A', 1, 2, 1, -2], ['A', 1, 1, 0, -1], ['A', 1, 1, 1, -1], ['A', 1, 0, 1, 0], ['A', 2, 2, 0, -2],
            ['A', 2, 2, 1, -2], ['A', 2, 1, 0, -1], ['A', 2, 1, 1, -1], ['A', 2, 0, 1, 0], ['B', 1, 0, 0, 0], ['B', 2, 0, 0, 0],
            ['X', 0, 2, 0, -2], ['X', 0, 2, 1, -2], ['X', 0, 1, 0, -1], ['X', 0, 1, 1, -1], ['X', 0, 0, 0, 0], ['X', 0, 0, 1, 0],
        ]),
        specialPoints: null,
        path: null,
    },
    'fluorite-ABX4-I41a': {
        parent: 'fluorite',
        title: 'scheelite type A B X4 (I4_1/a; fergusonite C2/c, I4_1/amd)',
        label: 'I4_1/a',
        supercell: [[1, 1, -2], [1, 1, 0], [0, -2, 1]],
        sites: superstructureSites([
            ['A', 0, 1, 0, -1], ['A', 0, 0, 0, 0], ['B', 0, 1, 1, -1], ['B', 0, 1, 0, 0], ['X', 1, 1, 0, -1], ['X', 1, 1, 1, -1],
            ['X', 1, 0, 0, 0], ['X', 1, 1, 0, 0], ['X', 2, 1, 0, -1], ['X', 2, 1, 1, -1], ['X', 2, 0, 0, 0], ['X', 2, 1, 0, 0],
        ]),
        specialPoints: null,
        path: null,
    },
    'fluorite-A3BX8-I4mmm': {
        parent: 'fluorite',
        title: 'A3 B X8 (I4/mmm)',
        label: 'I4/mmm',
        supercell: [[0, -2, 1], [-2, 0, 1], [1, 1, -2]],
        sites: superstructureSites([
            ['A', 0, -1, 0, 0], ['A', 0, 0, 0, 0], ['A', 0, -1, -1, 1], ['B', 0, 0, -1, 0], ['X', 1, 0, -1, 0], ['X', 1, -1, 0, 0],
            ['X', 1, 0, 0, 0], ['X', 1, -1, -1, 1], ['X', 2, 0, -1, 0], ['X', 2, -1, 0, 0], ['X', 2, 0, 0, 0], ['X', 2, -1, -1, 1],
        ]),
        specialPoints: null,
        path: null,
    },
    'fluorite-ABX4-R3m': {
        parent: 'fluorite',
        title: 'CaUO4 type A B X4 (R-3m)',
        label: 'R-3m',
        supercell: [[-1, -1, 0], [-1, 0, -1], [0, -1, -1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['B', 0, -1, -1, -1], ['X', 1, -1, -1, -1], ['X', 1, 0, 0, 0], ['X', 2, -1, -1, -1], ['X', 2, 0, 0, 0],
        ]),
        specialPoints: null,
        path: null,
    },
    'fluorite-A2BX6-Immm': {
        parent: 'fluorite',
        title: 'A2 B X6 (Immm)',
        label: 'Immm',
        supercell: [[1, 0, 1], [2, 0, -1], [-2, 1, 0]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 0, 1, 0, 0], ['B', 0, 2, 0, 0], ['X', 1, 0, 0, 0], ['X', 1, 1, 0, 0], ['X', 1, 2, 0, 0],
            ['X', 2, 0, 0, 0], ['X', 2, 1, 0, 0], ['X', 2, 2, 0, 0],
        ]),
        specialPoints: null,
        path: null,
    },
    'fluorite-A2BX6-C2c': {
        parent: 'fluorite',
        title: 'A2 B X6 (C2/c)',
        label: 'C2/c',
        supercell: [[-1, 0, 2], [2, 0, -1], [-1, 2, -1]],
        sites: superstructureSites([
            ['A', 0, 1, 1, -1], ['A', 0, 0, 0, 0], ['A', 0, 1, 0, 0], ['A', 0, -1, 1, 1], ['B', 0, 0, 1, 0], ['B', 0, 0, 0, 1],
            ['X', 1, 1, 1, -1], ['X', 1, 0, 0, 0], ['X', 1, 1, 0, 0], ['X', 1, 0, 1, 0], ['X', 1, 0, 0, 1], ['X', 1, -1, 1, 1],
            ['X', 2, 1, 1, -1], ['X', 2, 0, 0, 0], ['X', 2, 1, 0, 0], ['X', 2, 0, 1, 0], ['X', 2, 0, 0, 1], ['X', 2, -1, 1, 1],
        ]),
        specialPoints: null,
        path: null,
    },
    'fluorite-ABX4-Fd3m': {
        parent: 'fluorite',
        title: 'A B X4 (Fd-3m)',
        label: 'Fd-3m',
        supercell: [[-2, 0, 0], [0, 0, -2], [0, -2, 0]],
        sites: superstructureSites([
            ['A', 0, 0, 0, -1], ['A', 0, 0, -1, 0], ['A', 0, -1, 0, 0], ['A', 0, 0, 0, 0], ['B', 0, -1, -1, -1], ['B', 0, 0, -1, -1],
            ['B', 0, -1, 0, -1], ['B', 0, -1, -1, 0], ['X', 1, -1, -1, -1], ['X', 1, 0, -1, -1], ['X', 1, -1, 0, -1], ['X', 1, 0, 0, -1],
            ['X', 1, -1, -1, 0], ['X', 1, 0, -1, 0], ['X', 1, -1, 0, 0], ['X', 1, 0, 0, 0], ['X', 2, -1, -1, -1], ['X', 2, 0, -1, -1],
            ['X', 2, -1, 0, -1], ['X', 2, 0, 0, -1], ['X', 2, -1, -1, 0], ['X', 2, 0, -1, 0], ['X', 2, -1, 0, 0], ['X', 2, 0, 0, 0],
        ]),
        specialPoints: null,
        path: null,
    },
};

// perovskites: the 5-atom cell itself ('perovskite', 'antiperovskite': sub 0 =
// A, 1 = B, 2-4 = X) and the supercells of the octahedral tilts, one per
// supercell lattice (the sites carry no ordering, so every tilt system on the
// same supercell shares the template). Derived with analyzeSuperstructure from
// PhononDB mp-5777 (BaTiO3), mp-3466 (PrAlO3), mp-4438 (CaSnO3), mp-616378
// (CsSnI3), mp-669414 (PbHfO3), mp-3671 (NaNbO3), mp-985585 (Li3OCl) and mp-4192
// (Ca3AsN); referencePoints from Materials Project mp-30056 (CsCaBr3) and
// mp-2920 (LaAlO3). With no ordering on the sites every orientation of a tilted
// cell matches, but they put the tilts (or the polarization) in different senses:
// referenceStructure, the canonical cell (at unit volume) and positions of the
// material it was derived from, picks the match with the same distortion, so two
// materials are mixed with their tilts aligned (see matchScore).
export const PEROVSKITE_SUPERSTRUCTURES = {
    'perovskite': {
        parent: 'perovskite',
        title: 'A B X3 perovskite (Pm-3m, polar P4mm, R3m, Amm2)',
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['B', 1, 0, 0, 0], ['X', 2, 0, 0, 0], ['X', 3, 0, 0, 0], ['X', 4, 0, 0, 0],
        ]),
        referencePoints: { X: [0.5, 0, 0], M: [0.5, 0.5, 0], R: [0.5, 0.5, 0.5] },
        referenceStructure: {
            lattice: [[1.0055, 0.002, 0], [0.002, 1.0055, 0], [0, 0, 0.9892]],
            positions: [
                [0, 0, 0], [0.5138, 0.5138, 0.5], [0.4845, 0.4845, 0], [0.4872, -0.0229, 0.5],
                [-0.0229, 0.4872, 0.5],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'perovskite-R': {
        parent: 'perovskite',
        title: 'anti-phase tilts a-a-a-, a0a0c-, a-a-c0 (R-3c, R3c, I4/mcm, Imma)',
        supercell: [[1, 1, 0], [1, -1, 0], [-1, 0, -1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 0, 1, 0, 0], ['B', 1, 0, 0, 0], ['B', 1, 1, 0, 0], ['X', 2, 0, 0, 0], ['X', 2, 1, 0, 0],
            ['X', 3, 0, 0, 0], ['X', 3, 1, 0, 0], ['X', 4, 0, 0, 0], ['X', 4, 1, 0, 0],
        ]),
        referencePoints: { L: [0.5, 0, 0], Z: [0.5, 0, 0.5], F: [0.5, 0.5, 0.5] },
        referenceStructure: {
            lattice: [[0.7964, 0.7973, 0.0007], [0.7911, -0.7902, 0.0025], [-0.7984, 0.0018, -0.7953]],
            positions: [
                [0, 0, 0], [0.5, 0.5, 0], [0.25, -0.25, -0.5], [0.75, 0.25, -0.5],
                [0.4538, 0, 0], [1.0462, 0.5, 0], [0, 0, -0.5462], [0.5, 0.5, -0.4538],
                [-0.0462, -0.5, -0.5462], [0.5462, 0, -0.4538],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'perovskite-Pnma': {
        parent: 'perovskite',
        title: 'GdFeO3 type tilts a-b+a- (Pnma)',
        supercell: [[1, 1, 0], [0, 0, -2], [-1, 1, 0]],
        sites: superstructureSites([
            ['A', 0, 0, 0, -1], ['A', 0, 0, 1, -1], ['A', 0, 0, 0, 0], ['A', 0, 0, 1, 0], ['B', 1, 0, 0, -1], ['B', 1, 0, 1, -1],
            ['B', 1, 0, 0, 0], ['B', 1, 0, 1, 0], ['X', 2, 0, 0, -1], ['X', 2, 0, 1, -1], ['X', 2, 0, 0, 0], ['X', 2, 0, 1, 0],
            ['X', 3, 0, 0, -1], ['X', 3, 0, 1, -1], ['X', 3, 0, 0, 0], ['X', 3, 0, 1, 0], ['X', 4, 0, 0, -1], ['X', 4, 0, 1, -1],
            ['X', 4, 0, 0, 0], ['X', 4, 0, 1, 0],
        ]),
        referenceStructure: {
            lattice: [[0.6411, 0.6411, 0], [0, 0, -1.2574], [-0.6203, 0.6203, 0]],
            positions: [
                [-0.1107, 0.5, 0.0303], [0.3893, 0.5, 0.5], [0, 0, 0], [0.5, 0, 0.5303],
                [0.4447, 0.25, 0.0152], [0.9447, 0.25, 0.5152], [0.4447, -0.25, 0.0152], [0.9447, -0.25, 0.5152],
                [0.4872, 0.5, -0.0935], [0.9872, 0.5, 0.6239], [0.4021, 0, 0.1239], [0.9021, 0, 0.4065],
                [0.2449, 0.1931, -0.2889], [0.6444, 0.3069, 0.3192], [0.2449, -0.1931, -0.2889], [0.6444, -0.3069, 0.3192],
                [0.1444, 0.3069, 0.2111], [0.7449, 0.1931, 0.8192], [0.1444, -0.3069, 0.2111], [0.7449, -0.1931, 0.8192],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'perovskite-P4mbm': {
        parent: 'perovskite',
        title: 'in-phase tilt a0a0c+ (P4/mbm)',
        supercell: [[0, -1, -1], [0, -1, 1], [-1, 0, 0]],
        sites: superstructureSites([
            ['A', 0, 0, -1, 0], ['A', 0, 0, 0, 0], ['B', 1, 0, -1, 0], ['B', 1, 0, 0, 0], ['X', 2, 0, -1, 0], ['X', 2, 0, 0, 0],
            ['X', 3, 0, -1, 0], ['X', 3, 0, 0, 0], ['X', 4, 0, -1, 0], ['X', 4, 0, 0, 0],
        ]),
        referenceStructure: {
            lattice: [[0, -0.7899, -0.7899], [0, -0.7899, 0.7899], [-0.8013, 0, 0]],
            positions: [
                [0.5, 0.5, 0], [0, 0, 0], [0, 0.5, -0.5], [-0.5, 0, -0.5],
                [0.3079, 0.3079, -0.5], [-0.3079, -0.3079, -0.5], [0.1921, 0.8079, -0.5], [-0.1921, 0.1921, -0.5],
                [0, 0.5, 0], [-0.5, 0, 0],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'perovskite-Pbam': {
        parent: 'perovskite',
        title: 'antiferroelectric PbZrO3 type (Pbam)',
        supercell: [[-1, -1, 0], [2, -2, 0], [0, 0, 2]],
        sites: superstructureSites([
            ['A', 0, 1, -2, 0], ['A', 0, 0, -1, 0], ['A', 0, 1, -1, 0], ['A', 0, 0, 0, 0], ['A', 0, 1, -2, 1], ['A', 0, 0, -1, 1],
            ['A', 0, 1, -1, 1], ['A', 0, 0, 0, 1], ['B', 1, 1, -2, 0], ['B', 1, 0, -1, 0], ['B', 1, 1, -1, 0], ['B', 1, 0, 0, 0],
            ['B', 1, 1, -2, 1], ['B', 1, 0, -1, 1], ['B', 1, 1, -1, 1], ['B', 1, 0, 0, 1], ['X', 2, 1, -2, 0], ['X', 2, 0, -1, 0],
            ['X', 2, 1, -1, 0], ['X', 2, 0, 0, 0], ['X', 2, 1, -2, 1], ['X', 2, 0, -1, 1], ['X', 2, 1, -1, 1], ['X', 2, 0, 0, 1],
            ['X', 3, 1, -2, 0], ['X', 3, 0, -1, 0], ['X', 3, 1, -1, 0], ['X', 3, 0, 0, 0], ['X', 3, 1, -2, 1], ['X', 3, 0, -1, 1],
            ['X', 3, 1, -1, 1], ['X', 3, 0, 0, 1], ['X', 4, 1, -2, 0], ['X', 4, 0, -1, 0], ['X', 4, 1, -1, 0], ['X', 4, 0, 0, 0],
            ['X', 4, 1, -2, 1], ['X', 4, 0, -1, 1], ['X', 4, 1, -1, 1], ['X', 4, 0, 0, 1],
        ]),
        referenceStructure: {
            lattice: [[-0.5014, -0.5014, 0], [1.0054, -1.0054, 0], [0, 0, 0.9918]],
            positions: [
                [0.4168, 0.7527, 0], [0.5, 0.2527, 0], [-0.0832, 0.5, 0], [0, 0, 0],
                [0.4224, 0.7477, 0.5], [0.4944, 0.2477, 0.5], [-0.0776, 0.505, 0.5], [-0.0056, 0.005, 0.5],
                [-0.0478, 0.752, 0.2497], [-0.0353, 0.252, 0.2497], [-0.5478, 0.5007, 0.2497], [-0.5353, 0.0007, 0.2497],
                [-0.0478, 0.752, 0.7503], [-0.0353, 0.252, 0.7503], [-0.5478, 0.5007, 0.7503], [-0.5353, 0.0007, 0.7503],
                [-0.0209, 0.7211, 0], [-0.0623, 0.2211, 0], [-0.5209, 0.5316, 0], [-0.5623, 0.0316, 0],
                [0.0018, 0.7804, 0.5], [-0.0849, 0.2804, 0.5], [-0.4982, 0.4724, 0.5], [-0.5849, -0.0276, 0.5],
                [0.2084, 0.8764, 0.2067], [0.2084, 0.3764, 0.2291], [-0.265, 0.6163, 0.2797], [-0.3181, 0.1364, 0.2797],
                [0.2084, 0.8764, 0.7932], [0.2084, 0.3764, 0.7709], [-0.2651, 0.6163, 0.7203], [-0.3181, 0.1364, 0.7203],
                [0.235, 0.6364, 0.2797], [0.1819, 0.1163, 0.2797], [-0.2916, 0.3764, 0.2067], [-0.2916, -0.1236, 0.2291],
                [0.2349, 0.6364, 0.7203], [0.1819, 0.1163, 0.7203], [-0.2916, 0.3764, 0.7932], [-0.2916, -0.1236, 0.7709],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'perovskite-Pbcm': {
        parent: 'perovskite',
        title: 'NaNbO3 P phase type (Pbcm, Pca2_1)',
        supercell: [[-1, 1, 0], [1, 1, 0], [0, 0, -4]],
        sites: superstructureSites([
            ['A', 0, 0, 0, -3], ['A', 0, 0, 1, -3], ['A', 0, 0, 0, -2], ['A', 0, 0, 1, -2], ['A', 0, 0, 0, -1], ['A', 0, 0, 1, -1],
            ['A', 0, 0, 0, 0], ['A', 0, 0, 1, 0], ['B', 1, 0, 0, -3], ['B', 1, 0, 1, -3], ['B', 1, 0, 0, -2], ['B', 1, 0, 1, -2],
            ['B', 1, 0, 0, -1], ['B', 1, 0, 1, -1], ['B', 1, 0, 0, 0], ['B', 1, 0, 1, 0], ['X', 2, 0, 0, -3], ['X', 2, 0, 1, -3],
            ['X', 2, 0, 0, -2], ['X', 2, 0, 1, -2], ['X', 2, 0, 0, -1], ['X', 2, 0, 1, -1], ['X', 2, 0, 0, 0], ['X', 2, 0, 1, 0],
            ['X', 3, 0, 0, -3], ['X', 3, 0, 1, -3], ['X', 3, 0, 0, -2], ['X', 3, 0, 1, -2], ['X', 3, 0, 0, -1], ['X', 3, 0, 1, -1],
            ['X', 3, 0, 0, 0], ['X', 3, 0, 1, 0], ['X', 4, 0, 0, -3], ['X', 4, 0, 1, -3], ['X', 4, 0, 0, -2], ['X', 4, 0, 1, -2],
            ['X', 4, 0, 0, -1], ['X', 4, 0, 1, -1], ['X', 4, 0, 0, 0], ['X', 4, 0, 1, 0],
        ]),
        referenceStructure: {
            lattice: [[-0.4985, 0.4985, 0], [0.5053, 0.5053, 0], [0, 0, -1.9849]],
            positions: [
                [-0.002, 0.045, 0.75], [0.5162, 0.545, 0.75], [0, 0, 0.5], [0.5142, 0.5, 0.5],
                [-0.002, -0.045, 0.25], [0.5162, 0.455, 0.25], [0, 0, 0], [0.5142, 0.5, 0],
                [0.0127, 0.5242, 0.6252], [0.5015, 1.0242, 0.6252], [0.0127, 0.4758, 0.3748], [0.5015, 0.9758, 0.3748],
                [0.0127, 0.4758, 0.1252], [0.5015, 0.9758, 0.1252], [0.0127, 0.5242, -0.1252], [0.5015, 1.0242, -0.1252],
                [0.0766, 0.4747, 0.75], [0.4376, 0.9747, 0.75], [-0.0614, 0.5, 0.5], [0.5755, 1, 0.5],
                [0.0766, 0.5253, 0.25], [0.4376, 1.0253, 0.25], [-0.0614, 0.5, 0], [0.5755, 1, 0],
                [-0.2839, 0.2854, 0.643], [0.2992, 0.7066, 0.6077], [-0.2019, 0.2146, 0.357], [0.215, 0.7934, 0.3923],
                [-0.2019, 0.2146, 0.143], [0.215, 0.7934, 0.1077], [-0.2839, 0.2854, -0.143], [0.2992, 0.7066, -0.1077],
                [0.215, 0.2066, 0.6077], [0.7981, 0.7854, 0.643], [0.2992, 0.2934, 0.3923], [0.7161, 0.7146, 0.357],
                [0.2992, 0.2934, 0.1077], [0.7161, 0.7146, 0.143], [0.215, 0.2066, -0.1077], [0.7981, 0.7854, -0.143],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'antiperovskite': {
        parent: 'antiperovskite',
        title: 'X3 B A antiperovskite (Pm-3m)',
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['B', 1, 0, 0, 0], ['X', 2, 0, 0, 0], ['X', 3, 0, 0, 0], ['X', 4, 0, 0, 0],
        ]),
        referencePoints: { X: [0.5, 0, 0], M: [0.5, 0.5, 0], R: [0.5, 0.5, 0.5] },
        referenceStructure: {
            lattice: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            positions: [
                [0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5],
                [0, 0.5, 0.5],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'antiperovskite-Pnma': {
        parent: 'antiperovskite',
        title: 'GdFeO3 type tilts (Pnma)',
        supercell: [[0, -1, -1], [-2, 0, 0], [0, 1, -1]],
        sites: superstructureSites([
            ['A', 0, -1, 0, -1], ['A', 0, 0, 0, -1], ['A', 0, -1, 0, 0], ['A', 0, 0, 0, 0], ['B', 1, -1, 0, -1], ['B', 1, 0, 0, -1],
            ['B', 1, -1, 0, 0], ['B', 1, 0, 0, 0], ['X', 2, -1, 0, -1], ['X', 2, 0, 0, -1], ['X', 2, -1, 0, 0], ['X', 2, 0, 0, 0],
            ['X', 3, -1, 0, -1], ['X', 3, 0, 0, -1], ['X', 3, -1, 0, 0], ['X', 3, 0, 0, 0], ['X', 4, -1, 0, -1], ['X', 4, 0, 0, -1],
            ['X', 4, -1, 0, 0], ['X', 4, 0, 0, 0],
        ]),
        referenceStructure: {
            lattice: [[0, -0.629, -0.629], [-1.2639, 0, 0], [0, 0.6289, -0.6289]],
            positions: [
                [0.4668, 0.5, 0.5], [0.5, 0, 0.5043], [-0.0332, 0.5, 0.0043], [0, 0, 0],
                [-0.0166, 0.25, 0.5022], [-0.0166, -0.25, 0.5022], [-0.5166, 0.25, 0.0022], [-0.5166, -0.25, 0.0022],
                [0.266, 0.231, 0.7196], [0.266, -0.231, 0.7196], [-0.2992, 0.269, 0.2847], [-0.2992, -0.269, 0.2847],
                [0.2008, 0.269, 0.2196], [0.2008, -0.269, 0.2196], [-0.234, 0.231, -0.2153], [-0.234, -0.231, -0.2153],
                [-0.0112, 0.5, 0.539], [-0.022, 0, 0.4654], [-0.5112, 0.5, -0.0346], [-0.522, 0, 0.039],
            ],
        },
        specialPoints: null,
        path: null,
    },
};

// layered MX2 (transition-metal dichalcogenides): every polytype is its own
// family, matched as a single cell. Sites relative to the first metal, in the
// hexagonal cell wurtziteLattice(a, c) (1T: CdI2 type P-3m1; 2H: MoS2 type
// P6_3/mmc, two layers) or the rhombohedral primitive cell rhombohedralLattice(a,
// c) (3R: R3m, three layers, every atom on the 3-fold axis); the heights of the
// chalcogens are the usual ones, a cell is matched to the nearest site.
const LAYERED_POLYTYPES = {
    'tmd-1T': {
        cell: 'hexagonal',
        sites: [
            { role: 'A', red: [0, 0, 0] },
            { role: 'X', red: [1 / 3, 2 / 3, 0.25] },
            { role: 'X', red: [2 / 3, 1 / 3, -0.25] },
        ],
    },
    'tmd-2H': {
        cell: 'hexagonal',
        sites: [
            { role: 'A', red: [0, 0, 0] },
            { role: 'A', red: [1 / 3, -1 / 3, 0.5] },
            { role: 'X', red: [1 / 3, -1 / 3, -0.125] },
            { role: 'X', red: [1 / 3, -1 / 3, 0.125] },
            { role: 'X', red: [0, 0, 0.375] },
            { role: 'X', red: [0, 0, 0.625] },
        ],
    },
    // corundum (Al2O3, R-3c): the rhombohedral primitive cell, cations at (0,0,+-z),
    // (0,0,1/2+-z) and anions at (x,0,1/4) and equivalents in hexagonal axes, z = 0.352,
    // x = 0.306, here relative to the first cation in the rhombohedral basis
    corundum: {
        cell: 'rhombohedral',
        sites: [
            { role: 'A', red: [0, 0, 0] },
            { role: 'A', red: [0.296, 0.296, 0.296] },
            { role: 'A', red: [0.5, 0.5, 0.5] },
            { role: 'A', red: [-0.204, -0.204, -0.204] },
            { role: 'X', red: [0.204, -0.408, -0.102] },
            { role: 'X', red: [-0.102, 0.204, -0.408] },
            { role: 'X', red: [-0.408, -0.102, 0.204] },
            { role: 'X', red: [0.092, -0.296, 0.398] },
            { role: 'X', red: [0.398, 0.092, -0.296] },
            { role: 'X', red: [-0.296, 0.398, 0.092] },
        ],
    },
    'tmd-3R': {
        cell: 'rhombohedral',
        sites: [
            { role: 'A', red: [0, 0, 0] },
            { role: 'X', red: [0.25, 0.25, 0.25] },
            { role: 'X', red: [0.417, 0.417, 0.417] },
        ],
    },
};

// transition-metal dichalcogenides: the cell of each polytype (LAYERED_POLYTYPES)
// is its own template; no ordered ternaries are in the databases
const layeredTemplate = (name, title, extra = {}) => Object.assign({
    parent: name,
    title: title,
    supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    sites: superstructureSites(LAYERED_POLYTYPES[name].sites.map((site, sub) => [site.role, sub, 0, 0, 0])),
    specialPoints: null,
    path: null,
}, extra);
export const TMD_SUPERSTRUCTURES = {
    'tmd-1T': layeredTemplate('tmd-1T', '1T (CdI2 type, P-3m1)', {
        specialPoints: HEXAGONAL_SPECIAL_POINTS,
        referencePoints: { M: [0.5, 0, 0], K: [1 / 3, 1 / 3, 0], L: [0.5, 0, 0.5], H: [1 / 3, 1 / 3, 0.5] },
    }),
    'tmd-2H': layeredTemplate('tmd-2H', '2H (MoS2 type, P6_3/mmc)', {
        specialPoints: HEXAGONAL_SPECIAL_POINTS,
        referencePoints: { M: [0.5, 0, 0], K: [1 / 3, 1 / 3, 0], L: [0.5, 0, 0.5], H: [1 / 3, 1 / 3, 0.5] },
    }),
    'tmd-3R': layeredTemplate('tmd-3R', '3R (R3m)'),
};

// rutile: the 6-atom cell itself ('rutile', sub 0-1 cations, 2-5 anions) and cation
// orderings on supercells of it: trirutile (c tripled, A B B along the chains), and
// A B X4 / A B3 X8 orderings. Derived with analyzeSuperstructure from PhononDB
// mp-1873 (ZnF2), mp-3653 (MgSb2O6), mp-754478 (HfZrO4), mp-756638 (NbRhO4) and
// mp-29892 (LiSb3O8); referencePoints from the tetragonal Materials Project files
// (mp-571166 CaBr2, mp-776532 HfO2). As for the tilted perovskites, an orthorhombic
// (CaCl2 type) distortion fits the tetragonal template in two senses: each template
// keeps the distortion of a reference (for 'rutile' CaCl2 mp-23214, Pnnm) and the
// closest match is used (matchScore).
export const RUTILE_SUPERSTRUCTURES = {
    'rutile': {
        parent: 'rutile',
        title: 'A X2 rutile (P4_2/mnm, CaCl2 type Pnnm)',
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 1, 0, 0, 0], ['X', 2, 0, 0, 0], ['X', 3, 0, 0, 0], ['X', 4, 0, 0, 0], ['X', 5, 0, 0, 0],
        ]),
        referencePoints: { X: [0, 0.5, 0], M: [0.5, 0.5, 0], Z: [0, 0, 0.5], R: [0, 0.5, 0.5], A: [0.5, 0.5, 0.5] },
        referenceStructure: {
            lattice: [[1.173, 0, 0], [0, 1.1316, 0], [0, 0, 0.7534]],
            positions: [
                [0, 0, 0], [0.5, 0.5, 0.5], [0.3341, 0.2715, 0], [-0.3341, -0.2715, 0],
                [0.8341, 0.2285, 0.5], [0.1659, 0.7715, 0.5],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'rutile-AB2X6': {
        parent: 'rutile',
        title: 'trirutile A B2 X6 (P4_2/mnm)',
        supercell: [[-1, 0, 0], [0, -1, 0], [0, 0, 3]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 0, 0, 0, 1], ['A', 1, 0, 0, 1], ['A', 1, 0, 0, 2], ['B', 0, 0, 0, 2], ['B', 1, 0, 0, 0],
            ['X', 2, 0, 0, 0], ['X', 2, 0, 0, 1], ['X', 2, 0, 0, 2], ['X', 3, 0, 0, 0], ['X', 3, 0, 0, 1], ['X', 3, 0, 0, 2],
            ['X', 4, 0, 0, 0], ['X', 4, 0, 0, 1], ['X', 4, 0, 0, 2], ['X', 5, 0, 0, 0], ['X', 5, 0, 0, 1], ['X', 5, 0, 0, 2],
        ]),
        referenceStructure: {
            lattice: [[-0.7948, 0, 0], [0, -0.7948, 0], [0, 0, 1.5828]],
            positions: [
                [0, 0, 0], [0, 0, 0.3349], [-0.5, -0.5, 0.5], [-0.5, -0.5, 0.8349],
                [0, 0, 0.6675], [-0.5, -0.5, 0.1675], [-0.3024, -0.3024, -0.0038], [-0.3024, -0.3024, 0.3387],
                [-0.3094, -0.3094, 0.6675], [0.3024, 0.3024, -0.0038], [0.3024, 0.3024, 0.3387], [0.3094, 0.3094, 0.6675],
                [-0.8094, -0.1906, 0.1675], [-0.8024, -0.1976, 0.4962], [-0.8024, -0.1976, 0.8387], [-0.1906, -0.8094, 0.1675],
                [-0.1976, -0.8024, 0.4962], [-0.1976, -0.8024, 0.8387],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'rutile-ABX4-Cmmm': {
        parent: 'rutile',
        title: 'A B X4 (Cmmm)',
        supercell: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['B', 1, 0, 0, 0], ['X', 2, 0, 0, 0], ['X', 3, 0, 0, 0], ['X', 4, 0, 0, 0], ['X', 5, 0, 0, 0],
        ]),
        referenceStructure: {
            lattice: [[-1.1449, -0.0012, 0], [0.0012, 1.1449, 0], [0, 0, -0.7629]],
            positions: [
                [0, 0, 0], [-0.5, 0.5, -0.5], [-0.307, 0.307, 0], [0.307, -0.307, 0],
                [-0.8038, 0.1962, -0.5], [-0.1962, 0.8038, -0.5],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'rutile-ABX4-I41md': {
        parent: 'rutile',
        title: 'A B X4 (I4_1md)',
        supercell: [[-1, 0, 1], [1, 0, 1], [0, 1, -1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 1, 0, 0, 0], ['B', 0, 0, 0, 1], ['B', 1, 0, 0, 1], ['X', 2, 0, 0, 0], ['X', 2, 0, 0, 1],
            ['X', 3, 0, 0, 0], ['X', 3, 0, 0, 1], ['X', 4, 0, 0, 0], ['X', 4, 0, 0, 1], ['X', 5, 0, 0, 0], ['X', 5, 0, 0, 1],
        ]),
        referenceStructure: {
            lattice: [[-0.9291, 0, 0.5793], [0.9291, 0, 0.5793], [0, 0.9291, -0.5793]],
            positions: [
                [0, 0, 0], [0.25, 0.75, 0.5], [0.5216, 0.5216, 0], [0.7716, 1.2716, 0.5],
                [-0.0064, 0.2844, 0.2907], [0.5051, 0.8105, 0.3054], [-0.0064, -0.2971, -0.2907], [0.5051, 0.1997, -0.3054],
                [-0.0471, 0.7436, 0.2093], [0.4497, 1.2551, 0.1946], [0.5344, 0.7436, 0.7907], [1.0605, 1.2551, 0.8054],
            ],
        },
        specialPoints: null,
        path: null,
    },
    'rutile-AB3X8': {
        parent: 'rutile',
        title: 'A B3 X8 (P2_1/c)',
        supercell: [[-1, 0, -1], [0, -1, 0], [-1, 0, 3]],
        sites: superstructureSites([
            ['A', 0, -1, 0, 0], ['A', 0, 0, 0, 0], ['A', 0, -1, 0, 2], ['A', 1, -1, 0, 0], ['A', 1, -1, 0, 1], ['A', 1, -1, 0, 2],
            ['B', 0, -1, 0, 1], ['B', 1, 0, 0, 0], ['X', 2, -1, 0, 0], ['X', 2, 0, 0, 0], ['X', 2, -1, 0, 1], ['X', 2, -1, 0, 2],
            ['X', 3, -1, 0, 0], ['X', 3, 0, 0, 0], ['X', 3, -1, 0, 1], ['X', 3, -1, 0, 2], ['X', 4, -1, 0, 0], ['X', 4, 0, 0, 0],
            ['X', 4, -1, 0, 1], ['X', 4, -1, 0, 2], ['X', 5, -1, 0, 0], ['X', 5, 0, 0, 0], ['X', 5, -1, 0, 1], ['X', 5, -1, 0, 2],
        ]),
        referenceStructure: {
            lattice: [[-0.7243, 0, -0.4789], [0, -0.7277, 0], [-0.7211, 0, 1.4204]],
            positions: [
                [0.7367, -0.0215, 0.251], [0, 0, 0], [0.2633, 0.0215, 0.749], [0.2633, -0.5215, 0.249],
                [0, -0.5, 0.5], [-0.2633, -0.4785, 0.751], [0.5, 0, 0.5], [-0.5, -0.5, 0],
                [0.5126, -0.3094, 0.1743], [-0.2211, -0.2938, -0.0787], [0.2649, -0.3276, 0.4185], [0.0207, -0.2821, 0.6739],
                [0.9793, 0.2821, 0.3261], [0.2211, 0.2938, 0.0787], [0.7351, 0.3276, 0.5815], [0.4874, 0.3094, 0.8257],
                [0.0207, -0.2179, 0.1739], [-0.7351, -0.1724, -0.0815], [-0.2211, -0.2062, 0.4213], [-0.4874, -0.1906, 0.6743],
                [0.4874, -0.8094, 0.3257], [-0.2649, -0.8276, 0.0815], [0.2211, -0.7938, 0.5787], [-0.0207, -0.7821, 0.8261],
            ],
        },
        specialPoints: null,
        path: null,
    },
};

// corundum: the 10-atom rhombohedral cell itself ('corundum', sub 0-3 cations, 4-9
// anions, LAYERED_POLYTYPES.corundum) and the cation orderings on it: ilmenite
// (A and B in alternate layers along c) and A3 B X6. The LiNbO3 type (R3c) is a third
// ordering of the same cell; it is listed with the tilted perovskites. Derived with
// analyzeSuperstructure from PhononDB mp-1143 (Al2O3), mp-3771 (MgTiO3) and
// mp-675104 (Na3YCl6); referencePoints from Materials Project mp-974381 (Rh2S3).
export const CORUNDUM_SUPERSTRUCTURES = {
    'corundum': {
        parent: 'corundum',
        title: 'A2 X3 corundum (R-3c)',
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 1, 0, 0, 0], ['A', 2, 0, 0, 0], ['A', 3, 0, 0, 0], ['X', 4, 0, 0, 0], ['X', 5, 0, 0, 0], ['X', 6, 0, 0, 0], ['X', 7, 0, 0, 0], ['X', 8, 0, 0, 0], ['X', 9, 0, 0, 0],
        ]),
        referencePoints: { L: [0, 0, 0.5], Z: [0.5, 0.5, 0.5], F: [0.5, 0, 0.5] },
        specialPoints: null,
        path: null,
    },
    'corundum-ABX3': {
        parent: 'corundum',
        title: 'ilmenite A B X3 (R-3)',
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 1, 0, 0, 0], ['B', 2, 0, 0, 0], ['B', 3, 0, 0, 0], ['X', 4, 0, 0, 0], ['X', 5, 0, 0, 0], ['X', 6, 0, 0, 0], ['X', 7, 0, 0, 0], ['X', 8, 0, 0, 0], ['X', 9, 0, 0, 0],
        ]),
        referencePoints: { L: [0, 0, 0.5], Z: [0.5, 0.5, 0.5], F: [0.5, 0, 0.5] },
        specialPoints: null,
        path: null,
    },
    'corundum-A3BX6': {
        parent: 'corundum',
        title: 'A3 B X6 (R-3)',
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        sites: superstructureSites([
            ['A', 0, 0, 0, 0], ['A', 1, 0, 0, 0], ['A', 3, 0, 0, 0], ['B', 2, 0, 0, 0], ['X', 4, 0, 0, 0], ['X', 5, 0, 0, 0], ['X', 6, 0, 0, 0], ['X', 7, 0, 0, 0], ['X', 8, 0, 0, 0], ['X', 9, 0, 0, 0],
        ]),
        referencePoints: { L: [0, 0, 0.5], Z: [0.5, 0.5, 0.5], F: [0.5, 0, 0.5] },
        specialPoints: null,
        path: null,
    },
};

// the CsCl type A X on its 2-atom simple cubic cell, and the zinc-blende-derived
// cells with all four sites of the fcc cell filled (BCC_LAYOUTS): the Li3Bi (BiF3,
// D0_3) type A3 X and the Heusler orderings A2 B X on it, B on the octahedral site
// (full Heusler, Cu2MnAl type: K2CsSb, Li2NaSb, Ba2AuSb) or on a tetrahedral one
// (inverse Heusler, Hg2CuTi type: Li2CuSb). The half-Heuslers fill three of the four
// sites; they are ordered fluorites ('antifluorite-ABX-F43m'). Derived with analyzeSuperstructure from PhononDB mp-22865
// (CsCl), mp-2074 (Li3Sb), mp-581024 (K2CsSb) and Materials Project mp-15988
// (Li2CuSb). referencePoints: the standard simple cubic and fcc points.
const SIMPLE_CUBIC_POINTS = { X: [0, 0.5, 0], M: [0.5, 0.5, 0], R: [0.5, 0.5, 0.5] };
const FCC_POINTS = { X: [0.5, 0, 0.5], L: [0.5, 0.5, 0.5], W: [0.5, 0.25, 0.75], K: [0.375, 0.375, 0.75], U: [0.625, 0.25, 0.625] };
export const BCC_SUPERSTRUCTURES = {
    'cscl': {
        parent: 'cscl',
        title: 'CsCl type A X (Pm-3m)',
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        sites: superstructureSites([['A', 0, 0, 0, 0], ['X', 1, 0, 0, 0]]),
        referencePoints: SIMPLE_CUBIC_POINTS,
        specialPoints: null,
        path: null,
    },
    'heusler': {
        parent: 'heusler',
        title: 'Li3Bi type A3 X (Fm-3m)',
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        sites: superstructureSites([['A', 1, 0, 0, 0], ['A', 2, 0, 0, 0], ['A', 3, 0, 0, 0], ['X', 0, 0, 0, 0]]),
        referencePoints: FCC_POINTS,
        specialPoints: null,
        path: null,
    },
    'heusler-A2BX': {
        parent: 'heusler',
        title: 'full Heusler A2 B X (Cu2MnAl type, Fm-3m)',
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        sites: superstructureSites([['A', 2, 0, 0, 0], ['A', 3, 0, 0, 0], ['B', 1, 0, 0, 0], ['X', 0, 0, 0, 0]]),
        referencePoints: FCC_POINTS,
        specialPoints: null,
        path: null,
    },
    'heusler-A2BX-inverse': {
        parent: 'heusler',
        title: 'inverse Heusler A2 B X (Hg2CuTi type, F-43m)',
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        sites: superstructureSites([['A', 1, 0, 0, 0], ['A', 3, 0, 0, 0], ['B', 2, 0, 0, 0], ['X', 0, 0, 0, 0]]),
        referencePoints: FCC_POINTS,
        specialPoints: null,
        path: null,
    },
};

export const SUPERSTRUCTURES = Object.assign({}, WURTZITE_SUPERSTRUCTURES, ROCKSALT_SUPERSTRUCTURES,
    FLUORITE_SUPERSTRUCTURES, PEROVSKITE_SUPERSTRUCTURES, TMD_SUPERSTRUCTURES, RUTILE_SUPERSTRUCTURES,
    CORUNDUM_SUPERSTRUCTURES, BCC_SUPERSTRUCTURES);

function rowTimesMatrix(v, m) {
    return [0, 1, 2].map((j) => v[0] * m[0][j] + v[1] * m[1][j] + v[2] * m[2][j]);
}

function matrixTimesVector(m, v) {
    return m.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
}

// cations whose usual valence depends on the compound (Tl+ in DyTlS2, Tl3+ in NaTlO2)
const VARIABLE_VALENCE = { Tl: [1, 3], Pb: [2, 4], Sn: [2, 4], Cu: [1, 2], Ce: [3, 4] };

function getCationValences(counts, cations, anion) {
    /*
    formal valence of the two cations of a ternary: the choice that balances the
    anion charge, else the usual one
    */
    let options = cations.map((type) => VARIABLE_VALENCE[type] || [CATION_VALENCE[type] || 9]);
    let charge = ANION_CHARGE[anion] ? counts[anion] * ANION_CHARGE[anion] : null;
    let chosen = options.map((list) => list[0]);
    if (charge !== null) {
        search:
        for (let v0 of options[0]) {
            for (let v1 of options[1]) {
                if (counts[cations[0]] * v0 + counts[cations[1]] * v1 === charge) {
                    chosen = [v0, v1];
                    break search;
                }
            }
        }
    }
    let valence = {};
    cations.forEach((type, i) => { valence[type] = chosen[i]; });
    return valence;
}

function getSuperstructureRoles(atomTypes, xFraction = 0.5) {
    /*
    species -> role for a ternary cell: X is the anion, alone on its sublattice
    (half of the atoms in wurtzite and rocksalt, 2/3 in fluorite, 1/3 in
    antifluorite; a usual anion if several species qualify), the two cations on
    the ordered sublattice are A and B, ordered by number of atoms, then valence,
    then atomic number
    */
    let counts = {};
    atomTypes.forEach((type) => { counts[type] = (counts[type] || 0) + 1; });
    let species = Object.keys(counts);
    let isX = (type) => Math.abs(counts[type] - xFraction * atomTypes.length) < 1e-9;
    let anions = species.filter(isX);
    if (anions.length > 1) {
        anions = anions.filter((type) => ANION_CHARGE[type]);
    }
    let cations = species.filter((type) => type !== anions[0]);
    if (anions.length !== 1 || cations.length !== 2) {
        throw new Error('Not a ternary A B X cell: ' + JSON.stringify(counts));
    }
    let valence = getCationValences(counts, cations, anions[0]);
    let key = (type) => [-counts[type], valence[type], atomic_data.atomic_number[type]];
    cations.sort((p, q) => {
        let kp = key(p);
        let kq = key(q);
        return kp[0] - kq[0] || kp[1] - kq[1] || kp[2] - kq[2];
    });
    let roles = {};
    roles[cations[0]] = 'A';
    roles[cations[1]] = 'B';
    roles[anions[0]] = 'X';
    return roles;
}

function getNearestBonds(origin, targets, positionsCar, lattice, count) {
    // the `count` shortest vectors from atom `origin` to atoms in `targets` (any image)
    let bonds = [];
    for (let j of targets) {
        for (let n0 = -1; n0 <= 1; n0++) {
            for (let n1 = -1; n1 <= 1; n1++) {
                for (let n2 = -1; n2 <= 1; n2++) {
                    let v = [0, 1, 2].map((d) => positionsCar[j][d] + n0 * lattice[0][d] + n1 * lattice[1][d] + n2 * lattice[2][d] - positionsCar[origin][d]);
                    if (mat.vec_dot(v, v) > 1e-8) {
                        bonds.push(v);
                    }
                }
            }
        }
    }
    bonds.sort((p, q) => mat.vec_dot(p, p) - mat.vec_dot(q, q));
    return bonds.slice(0, count);
}

function* wurtziteDescriptions(atomTypes, positionsRed, lattice, isAnion) {
    /*
    every way of reading the cell as a supercell of an ideal wurtzite: an
    orthogonal frame T (rows: canonical x, y, z in the cartesian frame of the file,
    z along the cation -> anion bond parallel to c, a basal bond at -30 degrees
    from x as in wurtziteLattice) and a cation as origin. Yields T, the origin, the
    wurtzite lattice H, the supercell matrix M (file lattice = M H in the frame T)
    and, for every atom, its wurtzite sublattice and lattice vector.
    */
    let natoms = atomTypes.length;
    let positionsCar = utils.red_car_list(positionsRed, lattice);
    let cations = [];
    let anions = [];
    for (let i = 0; i < natoms; i++) {
        (isAnion(atomTypes[i]) ? anions : cations).push(i);
    }
    if (2 * cations.length !== natoms) {
        return;
    }
    let unit = (v) => mat.vec_scale(v, 1 / Math.sqrt(mat.vec_dot(v, v)));
    let bonds0 = getNearestBonds(cations[0], anions, positionsCar, lattice, 4).map(unit);
    let cosTol = Math.cos(20 * Math.PI / 180);
    let polar = bonds0.findIndex((v) => cations.every((i) =>
        getNearestBonds(i, anions, positionsCar, lattice, 4).some((b) => mat.vec_dot(unit(b), v) > cosTol)
    ));
    if (polar < 0) {
        return;
    }
    let ez = bonds0[polar];

    // c of the parent: the shortest lattice vector along z; a from the volume
    let c = Infinity;
    for (let n0 = -2; n0 <= 2; n0++) {
        for (let n1 = -2; n1 <= 2; n1++) {
            for (let n2 = -2; n2 <= 2; n2++) {
                let v = [0, 1, 2].map((d) => n0 * lattice[0][d] + n1 * lattice[1][d] + n2 * lattice[2][d]);
                let length = Math.sqrt(mat.vec_dot(v, v));
                if (length > 1e-6 && mat.vec_dot(v, ez) > length * Math.cos(5 * Math.PI / 180)) {
                    c = Math.min(c, length);
                }
            }
        }
    }
    if (!isFinite(c)) {
        return;
    }
    let volume = Math.abs(mat.matrix_determinant(lattice));
    let a = Math.sqrt((2 * volume) / (cations.length * (Math.sqrt(3) / 2) * c));
    let H = wurtziteLattice(a, c);
    let Hinverse = mat.matrix_inverse(H);
    let wurtziteOffsets = WURTZITE_SITES.map((site) => [site.red[0], site.red[1], site.red[2] + (site.role === 'X' ? 0.375 : 0)]);

    for (let b = 0; b < 4; b++) {
        if (b === polar) {
            continue;
        }
        let basal = [0, 1, 2].map((d) => bonds0[b][d] - mat.vec_dot(bonds0[b], ez) * ez[d]);
        basal = unit(basal);
        let cross = [ez[1] * basal[2] - ez[2] * basal[1], ez[2] * basal[0] - ez[0] * basal[2], ez[0] * basal[1] - ez[1] * basal[0]];
        for (let handedness of [1, -1]) {
            // the basal bond lies at -30 degrees from x: x = cos30 basal + sin30 (z x basal)
            let ex = [0, 1, 2].map((d) => (Math.sqrt(3) / 2) * basal[d] + 0.5 * handedness * cross[d]);
            let ey = [ez[1] * ex[2] - ez[2] * ex[1], ez[2] * ex[0] - ez[0] * ex[2], ez[0] * ex[1] - ez[1] * ex[0]];
            ey = mat.vec_scale(ey, handedness);
            let T = [ex, ey, ez];
            let rotated = lattice.map((row) => matrixTimesVector(T, row));
            let Mraw = rotated.map((row) => rowTimesMatrix(row, Hinverse));
            let M = Mraw.map((row) => row.map(Math.round));
            if (Mraw.some((row, i) => row.some((v, j) => Math.abs(v - M[i][j]) > 0.2)) ||
                Math.round(Math.abs(mat.matrix_determinant(M))) * 2 !== cations.length) {
                continue;
            }
            for (let origin of cations) {
                let sites = [];
                for (let i = 0; i < natoms; i++) {
                    let d = [0, 1, 2].map((k) => positionsCar[i][k] - positionsCar[origin][k]);
                    let f = rowTimesMatrix(matrixTimesVector(T, d), Hinverse);
                    let found = null;
                    for (let sub = 0; sub < 4 && !found; sub++) {
                        if ((WURTZITE_SITES[sub].role === 'X') !== isAnion(atomTypes[i])) {
                            continue;
                        }
                        let r = f.map((v, k) => v - wurtziteOffsets[sub][k]);
                        let n = r.map(Math.round);
                        if (Math.abs(r[0] - n[0]) < 0.15 && Math.abs(r[1] - n[1]) < 0.15 && Math.abs(r[2] - n[2]) < 0.12) {
                            found = { sub: sub, n: n };
                        }
                    }
                    if (!found) {
                        break;
                    }
                    sites.push(found);
                }
                if (sites.length === natoms) {
                    yield { T: T, origin: origin, H: H, M: M, sites: sites, a: a, c: c };
                }
            }
        }
    }
}

// parents on the fcc lattice rocksaltLattice(a): sublattice 0 is the fcc site at
// the origin, the others are offset by `red`. Rocksalt (Fm-3m): cation on the fcc
// site, anion at (1/2,1/2,1/2). Fluorite (CaF2, Fm-3m): cation on the fcc site,
// anions on both tetrahedral sites (1/4,1/4,1/4) and (3/4,3/4,3/4); antifluorite
// (Li2O, Mg2Si): the same with the anion on the fcc site. `centerIsX`: the anion
// (role X) sits on the fcc site; `axesFromOuter`: the cubic axes come from the
// nearest-neighbour vectors among the outer atoms (the simple cubic lattice of
// tetrahedral sites, spacing a/2) rather than the fcc site -> outer bonds;
// `xFraction`: the share of the atoms on the anion sublattice.
const ROCKSALT_SITES = [
    { role: 'C', red: [0, 0, 0] },
    { role: 'X', red: [0.5, 0.5, 0.5] },
];
const FLUORITE_SITES = [
    { red: [0, 0, 0] },
    { red: [0.25, 0.25, 0.25] },
    { red: [0.75, 0.75, 0.75] },
];
const FCC_LAYOUTS = {
    rocksalt: { sites: ROCKSALT_SITES, centerIsX: false, axesFromOuter: false, xFraction: 1 / 2 },
    fluorite: { sites: FLUORITE_SITES, centerIsX: false, axesFromOuter: true, xFraction: 2 / 3 },
    antifluorite: { sites: FLUORITE_SITES, centerIsX: true, axesFromOuter: true, xFraction: 1 / 3 },
};

export function rocksaltLattice(a) {
    return [[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]];
}

function* fccDescriptions(layout, atomTypes, positionsRed, lattice, isAnion) {
    /*
    every way of reading the cell as a supercell of a parent on the fcc lattice
    (FCC_LAYOUTS), as in wurtziteDescriptions. The cubic axes are three bond
    directions of length a/2 (rocksalt: cation -> anion of the octahedra;
    fluorite: between neighbouring tetrahedral sites), averaged over the cell;
    the cell is read in the fcc basis they span (e_y + e_z, e_x + e_z, e_x + e_y),
    so a homogeneous strain of the parent (the trigonal stretch of the layered
    alpha-NaFeO2 type reaches c/a = 5.9 instead of 4.9) drops out and only the
    local differences between the bonds remain. Every assignment of the three
    axes to x, y, z with either sign is tried; T is the orthogonal frame closest
    to it. The origin is an atom on the fcc site.
    */
    let natoms = atomTypes.length;
    let positionsCar = utils.red_car_list(positionsRed, lattice);
    let centers = [];
    let outer = [];
    for (let i = 0; i < natoms; i++) {
        (isAnion(atomTypes[i]) === layout.centerIsX ? centers : outer).push(i);
    }
    if (!centers.length || layout.sites.length * centers.length !== natoms) {
        return;
    }
    let sources = layout.axesFromOuter ? outer : centers;
    let cosine = (p, q) => mat.vec_dot(p, q) / Math.sqrt(mat.vec_dot(p, p) * mat.vec_dot(q, q));
    let axes = [];
    for (let bond of getNearestBonds(sources[0], outer, positionsCar, lattice, 6)) {
        if (!axes.some((axis) => Math.abs(cosine(axis, bond)) > 0.7)) {
            axes.push(bond);
        }
    }
    if (axes.length !== 3) {
        return;
    }
    let sums = axes.map(() => [0, 0, 0]);
    for (let source of sources) {
        for (let bond of getNearestBonds(source, outer, positionsCar, lattice, 6)) {
            let cosines = axes.map((axis) => cosine(axis, bond));
            let k = [0, 1, 2].reduce((best, j) => (Math.abs(cosines[j]) > Math.abs(cosines[best]) ? j : best), 0);
            let sign = cosines[k] > 0 ? 1 : -1;
            sums[k] = sums[k].map((v, d) => v + sign * bond[d]);
        }
    }
    let average = sums.map((v) => mat.vec_scale(v, 1 / (2 * sources.length)));
    let volume = Math.abs(mat.matrix_determinant(lattice));
    let a = Math.cbrt((4 * volume) / centers.length);
    let H = rocksaltLattice(a);
    let permutations = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];

    for (let permutation of permutations) {
        for (let signs = 0; signs < 8; signs++) {
            let e = permutation.map((k, i) => mat.vec_scale(average[k], (signs >> i) & 1 ? -1 : 1));
            let add = (p, q) => p.map((v, d) => v + q[d]);
            let basis = [add(e[1], e[2]), add(e[0], e[2]), add(e[0], e[1])];
            let inverse = mat.matrix_inverse(basis);
            let Mraw = lattice.map((row) => rowTimesMatrix(row, inverse));
            let M = Mraw.map((row) => row.map(Math.round));
            if (Mraw.some((row, r) => row.some((v, k) => Math.abs(v - M[r][k]) > 0.1)) ||
                Math.round(Math.abs(mat.matrix_determinant(M))) !== centers.length) {
                continue;
            }
            let T = orthogonalPolarFactor(e.map((v) => mat.vec_scale(v, 1 / Math.sqrt(mat.vec_dot(v, v)))));
            for (let origin of centers) {
                let sites = [];
                for (let n = 0; n < natoms; n++) {
                    let d = [0, 1, 2].map((k) => positionsCar[n][k] - positionsCar[origin][k]);
                    let f = rowTimesMatrix(d, inverse);
                    let candidates = isAnion(atomTypes[n]) === layout.centerIsX ? [0] : layout.sites.map((_, i) => i).slice(1);
                    let found = null;
                    for (let sub of candidates) {
                        let r = f.map((v, k) => v - layout.sites[sub].red[k]);
                        let rounded = r.map(Math.round);
                        if (r.every((v, k) => Math.abs(v - rounded[k]) < 0.25)) {
                            found = { sub: sub, n: rounded };
                        }
                    }
                    if (!found) {
                        break;
                    }
                    sites.push(found);
                }
                if (sites.length === natoms) {
                    yield { T: T, origin: origin, H: H, M: M, sites: sites, a: a };
                }
            }
        }
    }
}

// CsCl (B2, Pm-3m): the simple cubic cell b I, cation at the origin, anion at
// (1/2,1/2,1/2). Li3Bi / BiF3 (D0_3, Fm-3m) and the Heuslers ordered on it: the fcc
// cell rocksaltLattice(a) with the anion on the fcc site, then the octahedral site
// (1/2,1/2,1/2) and the two tetrahedral sites (1/4,1/4,1/4), (3/4,3/4,3/4). Zinc
// blende fills the fcc site and one tetrahedral site, a half-Heusler adds the
// octahedral one (an ordered fluorite, FLUORITE_SUPERSTRUCTURES), and these cells
// fill all four; a full Heusler A2 B X (Cu2MnAl type, K2CsSb) has B on the octahedral
// site, an inverse Heusler (Hg2CuTi type, Li2CuSb) on a tetrahedral one. A zinc blende
// itself cannot be folded into them: the two empty sites would be vacancies, which
// the virtual crystal (masses on every site) does not describe. In both cells every
// site of a body-centred cubic lattice of edge b (= a / 2 for the fcc cell) is
// occupied, which is how the cell is read (bccDescriptions).
const CSCL_SITES = [
    { red: [0, 0, 0] },
    { red: [0.5, 0.5, 0.5] },
];
const HEUSLER_SITES = [
    { red: [0, 0, 0] },
    { red: [0.5, 0.5, 0.5] },
    { red: [0.25, 0.25, 0.25] },
    { red: [0.75, 0.75, 0.75] },
];
const BCC_LAYOUTS = {
    cscl: { sites: CSCL_SITES, fcc: false },
    heusler: { sites: HEUSLER_SITES, fcc: true },
};

function getBccCubicAxes(positionsCar, lattice) {
    /*
    the cube edges (length b) of a cell whose atoms fill a body-centred cubic
    lattice: the eight nearest neighbours of an atom lie along the body diagonals,
    +-d_1 ... +-d_4 with the d_k towards the corners of a tetrahedron, and
    d_1 + d_k (k = 2, 3, 4) are the edges. The diagonals are averaged over the
    neighbours of every atom, so a homogeneous strain drops out as in
    fccDescriptions. Null when the neighbours are not along the diagonals.
    */
    let all = positionsCar.map((_, i) => i);
    let cosine = (p, q) => mat.vec_dot(p, q) / Math.sqrt(mat.vec_dot(p, p) * mat.vec_dot(q, q));
    let diagonals = [];
    for (let bond of getNearestBonds(0, all, positionsCar, lattice, 8)) {
        // the corners of a tetrahedron: cos = -1/3 between any two
        if (!diagonals.length || (diagonals.length < 4 && diagonals.every((d) => {
            let c = cosine(d, bond);
            return c < -0.1 && c > -0.6;
        }))) {
            diagonals.push(bond);
        }
    }
    if (diagonals.length !== 4) {
        return null;
    }
    let sums = diagonals.map(() => [0, 0, 0]);
    for (let i of all) {
        for (let bond of getNearestBonds(i, all, positionsCar, lattice, 8)) {
            let cosines = diagonals.map((d) => cosine(d, bond));
            let k = [0, 1, 2, 3].reduce((best, j) => (Math.abs(cosines[j]) > Math.abs(cosines[best]) ? j : best), 0);
            if (Math.abs(cosines[k]) < 0.85) {
                return null;
            }
            let sign = cosines[k] > 0 ? 1 : -1;
            sums[k] = sums[k].map((v, d) => v + sign * bond[d]);
        }
    }
    let average = sums.map((v) => mat.vec_scale(v, 1 / (2 * all.length)));
    return [1, 2, 3].map((k) => average[0].map((v, d) => v + average[k][d]));
}

function* bccDescriptions(layout, atomTypes, positionsRed, lattice) {
    /*
    every way of reading the cell as a supercell of a parent on the
    body-centred cubic lattice (BCC_LAYOUTS), as in fccDescriptions: the cube
    edges from getBccCubicAxes, assigned to x, y, z in every order and with either
    sign. Atoms are placed on the sites by position alone (the species decide the
    roles); any atom can be the origin.
    */
    let natoms = atomTypes.length;
    let positionsCar = utils.red_car_list(positionsRed, lattice);
    let axes = getBccCubicAxes(positionsCar, lattice);
    if (!axes) {
        return;
    }
    let volume = Math.abs(mat.matrix_determinant(lattice));
    let edge = Math.cbrt((2 * volume) / natoms);
    let H = layout.fcc ? rocksaltLattice(2 * edge) : [[edge, 0, 0], [0, edge, 0], [0, 0, edge]];
    let permutations = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let add = (p, q) => p.map((v, d) => v + q[d]);
    for (let permutation of permutations) {
        for (let signs = 0; signs < 8; signs++) {
            let e = permutation.map((k, i) => mat.vec_scale(axes[k], (signs >> i) & 1 ? -1 : 1));
            let basis = layout.fcc ? [add(e[1], e[2]), add(e[0], e[2]), add(e[0], e[1])] : e;
            let inverse = mat.matrix_inverse(basis);
            let Mraw = lattice.map((row) => rowTimesMatrix(row, inverse));
            let M = Mraw.map((row) => row.map(Math.round));
            if (Mraw.some((row, r) => row.some((v, k) => Math.abs(v - M[r][k]) > 0.1)) ||
                Math.round(Math.abs(mat.matrix_determinant(M))) * layout.sites.length !== natoms) {
                continue;
            }
            let T = orthogonalPolarFactor(e.map((v) => mat.vec_scale(v, 1 / Math.sqrt(mat.vec_dot(v, v)))));
            for (let origin = 0; origin < natoms; origin++) {
                let sites = [];
                for (let n = 0; n < natoms; n++) {
                    let d = [0, 1, 2].map((k) => positionsCar[n][k] - positionsCar[origin][k]);
                    let f = rowTimesMatrix(d, inverse);
                    // the sites are 1/4 apart (Heusler) or 1/2 (CsCl)
                    let found = null;
                    let best = 0.1;
                    layout.sites.forEach((site, sub) => {
                        let r = f.map((v, k) => v - site.red[k]);
                        let rounded = r.map(Math.round);
                        let deviation = Math.max.apply(null, r.map((v, k) => Math.abs(v - rounded[k])));
                        if (deviation < best) {
                            best = deviation;
                            found = { sub: sub, n: rounded };
                        }
                    });
                    if (!found) {
                        break;
                    }
                    sites.push(found);
                }
                if (sites.length === natoms) {
                    yield { T: T, origin: origin, H: H, M: M, sites: sites, a: layout.fcc ? 2 * edge : edge };
                }
            }
        }
    }
}

// how strongly a species takes the anion site of a CsCl or Heusler cell: a usual
// anion, then a p-block element (Sb in K3Sb, Sn in Sr2HgSn), then gold (the auride
// in CsAu, Ba2AuSb)
const bccAnionRank = (type) => (ANION_CHARGE[type] ? 3 : ANION_LIKE.has(type) ? 2 : type === 'Au' ? 1 : 0);

function getCsClRoles(atomTypes) {
    /*
    species -> role of a CsCl-type A X cell: X the anion (bccAnionRank), A the
    cation; in an intermetallic the heavier species takes X
    */
    let species = atomTypes.filter((type, i) => atomTypes.indexOf(type) === i);
    let count = (type) => atomTypes.filter((t) => t === type).length;
    if (species.length !== 2 || count(species[0]) !== count(species[1])) {
        throw new Error('Not an A X cell');
    }
    // an intermetallic (no anion): the heavier species on the X site
    let rank = (type) => bccAnionRank(type) * 1000 + atomic_data.atomic_number[type];
    let x = rank(species[0]) > rank(species[1]) ? species[0] : species[1];
    let roles = {};
    species.forEach((type) => { roles[type] = type === x ? 'X' : 'A'; });
    return roles;
}

function getHeuslerRoles(atomTypes) {
    /*
    species -> role of an A3 X (Li3Bi type) or A2 B X (Heusler) cell: X the
    anion with a quarter of the atoms (bccAnionRank; Sb rather than Au in Ba2AuSb,
    Cs in K2CsSb), A the cation with three quarters or half of them, B the other
    */
    let counts = {};
    atomTypes.forEach((type) => { counts[type] = (counts[type] || 0) + 1; });
    let species = Object.keys(counts);
    let quarter = species.filter((type) => 4 * counts[type] === atomTypes.length);
    // an intermetallic (no anion, CoOsAl2): the heavier species on the X site
    let rank = (type) => bccAnionRank(type) * 1000 + atomic_data.atomic_number[type];
    let x = quarter.slice().sort((p, q) => rank(q) - rank(p));
    if (!x.length) {
        throw new Error('No species with a quarter of the atoms: ' + JSON.stringify(counts));
    }
    let roles = {};
    roles[x[0]] = 'X';
    let cations = species.filter((type) => type !== x[0]);
    if (cations.length === 1 && 4 * counts[cations[0]] === 3 * atomTypes.length) {
        roles[cations[0]] = 'A';
    } else if (cations.length === 2 && cations.some((type) => 2 * counts[type] === atomTypes.length)) {
        cations.forEach((type) => { roles[type] = 2 * counts[type] === atomTypes.length ? 'A' : 'B'; });
    } else {
        throw new Error('Not an A3 X or A2 B X cell: ' + JSON.stringify(counts));
    }
    return roles;
}

// the parent lattices of the ordered ternaries: how to read a cell as a
// supercell of the parent, where a template site sits in the parent basis, and
// how far the frame from the bonds may be from the refined one
// perovskite ABX3 (Pm-3m): the simple cubic cell a I, A at the origin, B at
// (1/2,1/2,1/2), X at the face centres; antiperovskite (Li3OCl, Ca3PN) has the
// same sites with the cation on X. A and B are told apart by geometry: B is the
// centre of the X octahedra.
const PEROVSKITE_SITES = [
    { red: [0, 0, 0] },
    { red: [0.5, 0.5, 0.5] },
    { red: [0.5, 0.5, 0] },
    { red: [0.5, 0, 0.5] },
    { red: [0, 0.5, 0.5] },
];

function getPerovskiteRoles(atomTypes, positionsRed, lattice) {
    /*
    species -> role of an A B X3 cell: X has 3/5 of the atoms; B is the cation
    with six X neighbours forming an octahedron (the six within a factor 1.45 of
    each other, as in tetragonal PbTiO3, the seventh farther), A the other one;
    chlorates and iodates (three short bonds) are not perovskites
    */
    let counts = {};
    atomTypes.forEach((type) => { counts[type] = (counts[type] || 0) + 1; });
    let species = Object.keys(counts);
    let x = species.filter((type) => 5 * counts[type] === 3 * atomTypes.length);
    let cations = species.filter((type) => 5 * counts[type] === atomTypes.length);
    if (species.length !== 3 || x.length !== 1 || cations.length !== 2) {
        throw new Error('Not an A B X3 cell: ' + JSON.stringify(counts));
    }
    let positionsCar = utils.red_car_list(positionsRed, lattice);
    let xAtoms = atomTypes.map((type, i) => i).filter((i) => atomTypes[i] === x[0]);
    let octahedral = (type) => atomTypes.every((t, i) => {
        if (t !== type) {
            return true;
        }
        let d = getNearestBonds(i, xAtoms, positionsCar, lattice, 7).map((v) => Math.sqrt(mat.vec_dot(v, v)));
        return d[5] / d[0] < 1.45 && d[6] / d[5] > 1.1;
    });
    let mean = (type) => {
        let total = 0;
        atomTypes.forEach((t, i) => {
            if (t === type) {
                getNearestBonds(i, xAtoms, positionsCar, lattice, 6).forEach((v) => { total += Math.sqrt(mat.vec_dot(v, v)); });
            }
        });
        return total;
    };
    let b = cations.filter(octahedral);
    if (b.length === 2) {
        b = [mean(cations[0]) <= mean(cations[1]) ? cations[0] : cations[1]];
    }
    if (b.length !== 1) {
        throw new Error('No B X6 octahedra: not a perovskite');
    }
    let roles = {};
    roles[cations.find((type) => type !== b[0])] = 'A';
    roles[b[0]] = 'B';
    roles[x[0]] = 'X';
    return roles;
}

function* perovskiteDescriptions(atomTypes, positionsRed, lattice, isAnion, roles) {
    /*
    every way of reading an A B X3 cell as a supercell of the cubic perovskite, as
    in fccDescriptions: the cubic axes are the B -> X bonds of the octahedra
    (length a/2) averaged over all B, so tilts of the octahedra, which alternate
    between neighbours, average out; the cell is read in the basis of twice these
    vectors. The origin is an A atom.
    */
    let natoms = atomTypes.length;
    let positionsCar = utils.red_car_list(positionsRed, lattice);
    let role = (i) => roles[atomTypes[i]];
    let aAtoms = [];
    let bAtoms = [];
    let xAtoms = [];
    for (let i = 0; i < natoms; i++) {
        (role(i) === 'A' ? aAtoms : role(i) === 'B' ? bAtoms : xAtoms).push(i);
    }
    if (!bAtoms.length || aAtoms.length !== bAtoms.length || xAtoms.length !== 3 * bAtoms.length) {
        return;
    }
    let cosine = (p, q) => mat.vec_dot(p, q) / Math.sqrt(mat.vec_dot(p, p) * mat.vec_dot(q, q));
    let axes = [];
    for (let bond of getNearestBonds(bAtoms[0], xAtoms, positionsCar, lattice, 6)) {
        if (!axes.some((axis) => Math.abs(cosine(axis, bond)) > 0.7)) {
            axes.push(bond);
        }
    }
    if (axes.length !== 3) {
        return;
    }
    let sums = axes.map(() => [0, 0, 0]);
    for (let source of bAtoms) {
        for (let bond of getNearestBonds(source, xAtoms, positionsCar, lattice, 6)) {
            let cosines = axes.map((axis) => cosine(axis, bond));
            let k = [0, 1, 2].reduce((best, j) => (Math.abs(cosines[j]) > Math.abs(cosines[best]) ? j : best), 0);
            let sign = cosines[k] > 0 ? 1 : -1;
            sums[k] = sums[k].map((v, d) => v + sign * bond[d]);
        }
    }
    // twice the mean bond: the cubic lattice vectors. Tilts shorten the mean by
    // cos(tilt) along every axis (6 % at 20 degrees): scale to the volume per B
    let volume = Math.abs(mat.matrix_determinant(lattice));
    let average = sums.map((v) => mat.vec_scale(v, 1 / bAtoms.length));
    let scale = Math.cbrt(volume / (bAtoms.length * Math.abs(mat.matrix_determinant(average))));
    average = average.map((v) => mat.vec_scale(v, scale));
    let a = Math.cbrt(volume / bAtoms.length);
    let H = [[a, 0, 0], [0, a, 0], [0, 0, a]];
    let permutations = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];

    for (let permutation of permutations) {
        for (let signs = 0; signs < 8; signs++) {
            let basis = permutation.map((k, i) => mat.vec_scale(average[k], (signs >> i) & 1 ? -1 : 1));
            let inverse = mat.matrix_inverse(basis);
            let Mraw = lattice.map((row) => rowTimesMatrix(row, inverse));
            let M = Mraw.map((row) => row.map(Math.round));
            if (Mraw.some((row, r) => row.some((v, k) => Math.abs(v - M[r][k]) > 0.1)) ||
                Math.round(Math.abs(mat.matrix_determinant(M))) !== bAtoms.length) {
                continue;
            }
            let T = orthogonalPolarFactor(basis.map((v) => mat.vec_scale(v, 1 / Math.sqrt(mat.vec_dot(v, v)))));
            for (let origin of aAtoms) {
                let sites = [];
                for (let n = 0; n < natoms; n++) {
                    let d = [0, 1, 2].map((k) => positionsCar[n][k] - positionsCar[origin][k]);
                    let f = rowTimesMatrix(d, inverse);
                    let candidates = role(n) === 'A' ? [0] : role(n) === 'B' ? [1] : [2, 3, 4];
                    // the nearest of the candidate sites; strong tilts (NdRhO3, YScO3)
                    // move atoms up to ~0.27 a from the ideal site
                    let found = null;
                    let best = 0.35;
                    for (let sub of candidates) {
                        let r = f.map((v, k) => v - PEROVSKITE_SITES[sub].red[k]);
                        let rounded = r.map(Math.round);
                        let deviation = Math.max.apply(null, r.map((v, k) => Math.abs(v - rounded[k])));
                        if (deviation < best) {
                            best = deviation;
                            found = { sub: sub, n: rounded };
                        }
                    }
                    if (!found) {
                        break;
                    }
                    sites.push(found);
                }
                if (sites.length === natoms) {
                    yield { T: T, origin: origin, H: H, M: M, sites: sites, a: a };
                }
            }
        }
    }
}

export function rhombohedralLattice(a, c) {
    // primitive vectors of the rhombohedral lattice with hexagonal axes a, c
    let r = a / Math.sqrt(3);
    return [[a / 2, r / 2, c / 3], [-a / 2, r / 2, c / 3], [0, -r, c / 3]];
}

function getLayeredRoles(atomTypes) {
    // species -> role of an M X2 cell: X the chalcogen (2/3 of the atoms), A the metal
    let counts = {};
    atomTypes.forEach((type) => { counts[type] = (counts[type] || 0) + 1; });
    let species = Object.keys(counts);
    let x = species.filter((type) => ['S', 'Se', 'Te'].indexOf(type) !== -1 && 3 * counts[type] === 2 * atomTypes.length);
    if (species.length !== 2 || x.length !== 1) {
        throw new Error('Not an M X2 dichalcogenide cell: ' + JSON.stringify(counts));
    }
    let roles = {};
    roles[species.find((type) => type !== x[0])] = 'A';
    roles[x[0]] = 'X';
    return roles;
}

function* layeredDescriptions(polytype, atomTypes, positionsRed, lattice, isAnion) {
    /*
    every way of reading a layered M X2 cell as the polytype's cell: the in-plane
    axis x along one of the six shortest lattice vectors (in the layers), z along
    either normal, either handedness; a and c from the cell (c from the volume, so
    the interlayer distance is the file's own). The origin is a metal atom.
    */
    let natoms = atomTypes.length;
    if (natoms !== polytype.sites.length) {
        return;
    }
    let positionsCar = utils.red_car_list(positionsRed, lattice);
    let cross = (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];
    let unit = (v) => mat.vec_scale(v, 1 / Math.sqrt(mat.vec_dot(v, v)));
    let vectors = [];
    for (let n0 = -2; n0 <= 2; n0++) {
        for (let n1 = -2; n1 <= 2; n1++) {
            for (let n2 = -2; n2 <= 2; n2++) {
                let v = [0, 1, 2].map((d) => n0 * lattice[0][d] + n1 * lattice[1][d] + n2 * lattice[2][d]);
                if (mat.vec_dot(v, v) > 1e-6) {
                    vectors.push(v);
                }
            }
        }
    }
    vectors.sort((p, q) => mat.vec_dot(p, p) - mat.vec_dot(q, q));
    let a = Math.sqrt(mat.vec_dot(vectors[0], vectors[0]));
    let inPlane = vectors.filter((v) => Math.sqrt(mat.vec_dot(v, v)) < 1.02 * a);
    if (inPlane.length !== 6) {
        return;
    }
    let other = inPlane.find((v) => Math.abs(mat.vec_dot(unit(v), unit(inPlane[0]))) < 0.9);
    let normal = unit(cross(inPlane[0], other));
    let volume = Math.abs(mat.matrix_determinant(lattice));
    let hexagonalVolume = (a * a * Math.sqrt(3)) / 2;
    let c = polytype.cell === 'hexagonal' ? volume / hexagonalVolume : (3 * volume) / hexagonalVolume;
    let H = polytype.cell === 'hexagonal' ? wurtziteLattice(a, c) : rhombohedralLattice(a, c);
    let Hinverse = mat.matrix_inverse(H);
    // every full reading, offered closest first: a site pattern that is nearly its own
    // mirror image (the corundum cation column: -0.204, 0, 0.296, 0.5) also fits
    // mirrored, with other roles, within the tolerance (ilmenite: 0.09 against 0.01)
    let readings = [];

    for (let axis of inPlane) {
        for (let sign of [1, -1]) {
            for (let handedness of [1, -1]) {
                let ex = unit(axis);
                let ez = mat.vec_scale(normal, sign);
                let ey = mat.vec_scale(cross(ez, ex), handedness);
                let T = [ex, ey, ez];
                let rotated = lattice.map((row) => matrixTimesVector(T, row));
                let Mraw = rotated.map((row) => rowTimesMatrix(row, Hinverse));
                let M = Mraw.map((row) => row.map(Math.round));
                if (Mraw.some((row, r) => row.some((v, k) => Math.abs(v - M[r][k]) > 0.1)) ||
                    Math.round(Math.abs(mat.matrix_determinant(M))) !== 1) {
                    continue;
                }
                for (let origin = 0; origin < natoms; origin++) {
                    if (isAnion(atomTypes[origin])) {
                        continue;
                    }
                    let sites = [];
                    for (let n = 0; n < natoms; n++) {
                        let d = [0, 1, 2].map((k) => positionsCar[n][k] - positionsCar[origin][k]);
                        let f = rowTimesMatrix(matrixTimesVector(T, d), Hinverse);
                        let found = null;
                        let best = 0.12;
                        polytype.sites.forEach((site, sub) => {
                            if ((site.role === 'X') !== isAnion(atomTypes[n])) {
                                return;
                            }
                            let r = f.map((v, k) => v - site.red[k]);
                            let rounded = r.map(Math.round);
                            let deviation = Math.max.apply(null, r.map((v, k) => Math.abs(v - rounded[k])));
                            if (deviation < best) {
                                best = deviation;
                                found = { sub: sub, n: rounded };
                            }
                        });
                        if (!found) {
                            break;
                        }
                        found.deviation = best;
                        sites.push(found);
                    }
                    if (sites.length === natoms && new Set(sites.map((site) => site.sub)).size === natoms) {
                        readings.push({ deviation: Math.max.apply(null, sites.map((site) => site.deviation)),
                            reading: { T: T, origin: origin, H: H, M: M, sites: sites, a: a, c: c } });
                    }
                }
            }
        }
    }
    readings.sort((p, q) => p.deviation - q.deviation);
    for (let entry of readings) {
        yield entry.reading;
    }
}

function getCorundumRoles(atomTypes, positionsRed, lattice) {
    /*
    species -> role of a corundum-derived cell: X the anion (3/5 of the atoms, a
    usual anion: anti-corundum Ca3N2 is left out), A the cation, or A and B for an
    ordered ternary (ilmenite). Every cation must sit in an octahedron of six
    anions (within a factor 1.35, the seventh farther): the calcite-type carbonates,
    borates and nitrates (planar CO3 / BO3 / NO3) and the P2S6 / Si2Te6 layers have
    the same site pattern but not this chemistry
    */
    let counts = {};
    atomTypes.forEach((type) => { counts[type] = (counts[type] || 0) + 1; });
    let species = Object.keys(counts);
    let x = species.filter((type) => 5 * counts[type] === 3 * atomTypes.length && ANION_LIKE.has(type));
    if (x.length !== 1 || species.length < 2 || species.length > 3) {
        throw new Error('Not an A2 X3 corundum cell: ' + JSON.stringify(counts));
    }
    if (positionsRed && lattice) {
        let positionsCar = utils.red_car_list(positionsRed, lattice);
        let anions = atomTypes.map((type, i) => i).filter((i) => atomTypes[i] === x[0]);
        atomTypes.forEach((type, i) => {
            if (type === x[0]) {
                return;
            }
            let d = getNearestBonds(i, anions, positionsCar, lattice, 7).map((v) => Math.sqrt(mat.vec_dot(v, v)));
            if (!(d[5] / d[0] < 1.35 && d[6] / d[5] > 1.1)) {
                throw new Error('No ' + type + ' ' + x[0] + '6 octahedra: not a corundum-type cell');
            }
        });
    }
    if (species.length === 3) {
        return getSuperstructureRoles(atomTypes, 3 / 5);
    }
    let roles = {};
    roles[species.find((type) => type !== x[0])] = 'A';
    roles[x[0]] = 'X';
    return roles;
}

const layeredParent = (name) => ({
    describe: (...args) => layeredDescriptions(LAYERED_POLYTYPES[name], ...args),
    sitePosition: (site) => site.n.map((v, k) => v + LAYERED_POLYTYPES[name].sites[site.sub].red[k]),
    frameTolerance: 15,
    xFraction: 2 / 3,
    originRole: 'A',
    getRoles: getLayeredRoles,
});

// rutile (TiO2, P4_2/mnm): the tetragonal cell diag(a, a, c), cations at the
// origin and at the body centre, anions at +-(u, u, 0) and (1/2 +- u, 1/2 -+ u,
// 1/2) with u ~ 0.305 (sites are matched to the nearest, u varies a little)
const RUTILE_U = 0.305;
const RUTILE_SITES = [
    { role: 'C', red: [0, 0, 0] },
    { role: 'C', red: [0.5, 0.5, 0.5] },
    { role: 'X', red: [RUTILE_U, RUTILE_U, 0] },
    { role: 'X', red: [-RUTILE_U, -RUTILE_U, 0] },
    { role: 'X', red: [0.5 + RUTILE_U, 0.5 - RUTILE_U, 0.5] },
    { role: 'X', red: [0.5 - RUTILE_U, 0.5 + RUTILE_U, 0.5] },
];

function getRutileRoles(atomTypes) {
    /*
    species -> role of a rutile-derived A X2 / ternary cell: X the anion (2/3 of the
    atoms, a usual anion), the cations A (and B for an ordered ternary, by number
    of atoms, then valence)
    */
    let counts = {};
    atomTypes.forEach((type) => { counts[type] = (counts[type] || 0) + 1; });
    let species = Object.keys(counts);
    let x = species.filter((type) => 3 * counts[type] === 2 * atomTypes.length && ANION_LIKE.has(type));
    if (x.length !== 1 || species.length < 2 || species.length > 3) {
        throw new Error('Not an A X2 rutile cell: ' + JSON.stringify(counts));
    }
    if (species.length === 3) {
        return getSuperstructureRoles(atomTypes, 2 / 3);
    }
    let roles = {};
    roles[species.find((type) => type !== x[0])] = 'A';
    roles[x[0]] = 'X';
    return roles;
}

function* rutileDescriptions(atomTypes, positionsRed, lattice, isAnion) {
    /*
    every way of reading the cell as a supercell of rutile: c along the shortest
    cation -> cation vectors (the chains of edge-sharing octahedra, +-c), averaged
    over all cations; the a axes from the next cation neighbours (+-a/2, +-a/2,
    +-c/2) projected on the plane, as the sums of two at right angles; a from the
    volume. Every in-plane axis, either direction of c and either handedness is
    tried. The origin is a cation.
    */
    let natoms = atomTypes.length;
    let positionsCar = utils.red_car_list(positionsRed, lattice);
    let cations = [];
    let anions = [];
    for (let i = 0; i < natoms; i++) {
        (isAnion(atomTypes[i]) ? anions : cations).push(i);
    }
    if (cations.length < 2 || cations.length % 2 || 2 * cations.length !== anions.length) {
        return;
    }
    let unit = (v) => mat.vec_scale(v, 1 / Math.sqrt(mat.vec_dot(v, v)));
    let cross = (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];
    let near = getNearestBonds(cations[0], cations, positionsCar, lattice, 10);
    let axis = unit(near[0]);
    // c: the chain vector, averaged over all cations with a common sign
    let sum = [0, 0, 0];
    for (let cation of cations) {
        for (let bond of getNearestBonds(cation, cations, positionsCar, lattice, 2)) {
            let sign = mat.vec_dot(bond, axis) > 0 ? 1 : -1;
            sum = sum.map((v, d) => v + sign * bond[d]);
        }
    }
    let cVector = mat.vec_scale(sum, 1 / (2 * cations.length));
    let c = Math.sqrt(mat.vec_dot(cVector, cVector));
    let ezBase = unit(cVector);
    if (near.slice(0, 2).some((bond) => Math.abs(mat.vec_dot(unit(bond), ezBase)) < 0.95)) {
        return;
    }
    // in-plane directions of the next shell, and the a axes between them
    let planar = near.slice(2).map((bond) => unit(bond.map((v, d) => v - mat.vec_dot(bond, ezBase) * ezBase[d])));
    let directions = [];
    for (let p of planar) {
        for (let q of planar) {
            if (Math.abs(mat.vec_dot(p, q)) < 0.25) {
                let candidate = unit(p.map((v, d) => v + q[d]));
                if (!directions.some((other) => mat.vec_dot(other, candidate) > 0.95)) {
                    directions.push(candidate);
                }
            }
        }
    }
    if (!directions.length) {
        return;
    }
    let volume = Math.abs(mat.matrix_determinant(lattice));
    let a = Math.sqrt((2 * volume) / (cations.length * c));
    let H = [[a, 0, 0], [0, a, 0], [0, 0, c]];
    let Hinverse = mat.matrix_inverse(H);

    for (let ex of directions) {
        for (let sign of [1, -1]) {
            for (let handedness of [1, -1]) {
                let ez = mat.vec_scale(ezBase, sign);
                let ey = mat.vec_scale(cross(ez, ex), handedness);
                let T = [ex, ey, ez];
                let rotated = lattice.map((row) => matrixTimesVector(T, row));
                let Mraw = rotated.map((row) => rowTimesMatrix(row, Hinverse));
                let M = Mraw.map((row) => row.map(Math.round));
                if (Mraw.some((row, r) => row.some((v, k) => Math.abs(v - M[r][k]) > 0.1)) ||
                    Math.round(Math.abs(mat.matrix_determinant(M))) * 2 !== cations.length) {
                    continue;
                }
                for (let origin of cations) {
                    let sites = [];
                    for (let n = 0; n < natoms; n++) {
                        let d = [0, 1, 2].map((k) => positionsCar[n][k] - positionsCar[origin][k]);
                        let f = rowTimesMatrix(matrixTimesVector(T, d), Hinverse);
                        // rutile and its CaCl2-type distortions stay within 0.06 of the
                        // ideal sites; marcasite (FeS2, OsAs2: anion pairs) is at 0.1 and more
                        let found = null;
                        let best = 0.09;
                        RUTILE_SITES.forEach((site, sub) => {
                            if ((site.role === 'X') !== isAnion(atomTypes[n])) {
                                return;
                            }
                            let r = f.map((v, k) => v - site.red[k]);
                            let rounded = r.map(Math.round);
                            let deviation = Math.max.apply(null, r.map((v, k) => Math.abs(v - rounded[k])));
                            if (deviation < best) {
                                best = deviation;
                                found = { sub: sub, n: rounded };
                            }
                        });
                        if (!found) {
                            break;
                        }
                        sites.push(found);
                    }
                    if (sites.length === natoms) {
                        yield { T: T, origin: origin, H: H, M: M, sites: sites, a: a, c: c };
                    }
                }
            }
        }
    }
}

const perovskiteParent = (anti) => ({
    describe: perovskiteDescriptions,
    sitePosition: (site) => site.n.map((v, k) => v + PEROVSKITE_SITES[site.sub].red[k]),
    frameTolerance: 25,
    xFraction: 3 / 5,
    originRole: 'A',
    getRoles: getPerovskiteRoles,
    // the X species is an anion in perovskite, a cation in antiperovskite
    acceptsX: (type) => ANION_LIKE.has(type) !== anti,
});

const fccParent = (name) => ({
    describe: (...args) => fccDescriptions(FCC_LAYOUTS[name], ...args),
    sitePosition: (site) => site.n.map((v, k) => v + FCC_LAYOUTS[name].sites[site.sub].red[k]),
    frameTolerance: 25,
    // fraction of the atoms on the anion sublattice, and the role on the fcc site
    xFraction: FCC_LAYOUTS[name].xFraction,
    originRole: FCC_LAYOUTS[name].centerIsX ? 'X' : 'A',
    lattice: true,
});
const bccParent = (name, xFraction, originRole, getRoles) => ({
    describe: (atomTypes, positionsRed, lattice) => bccDescriptions(BCC_LAYOUTS[name], atomTypes, positionsRed, lattice),
    sitePosition: (site) => site.n.map((v, k) => v + BCC_LAYOUTS[name].sites[site.sub].red[k]),
    frameTolerance: 25,
    xFraction: xFraction,
    originRole: originRole,
    getRoles: getRoles,
});
const PARENT_STRUCTURES = {
    wurtzite: {
        describe: wurtziteDescriptions,
        sitePosition: (site) => templateParentPosition(site, 0.375),
        frameTolerance: 15,
        xFraction: 0.5,
        originRole: 'A',
    },
    rocksalt: fccParent('rocksalt'),
    fluorite: fccParent('fluorite'),
    antifluorite: fccParent('antifluorite'),
    perovskite: perovskiteParent(false),
    antiperovskite: perovskiteParent(true),
    rutile: {
        describe: rutileDescriptions,
        sitePosition: (site) => site.n.map((v, k) => v + RUTILE_SITES[site.sub].red[k]),
        frameTolerance: 15,
        xFraction: 2 / 3,
        originRole: 'A',
        getRoles: getRutileRoles,
    },
    'tmd-1T': layeredParent('tmd-1T'),
    'tmd-2H': layeredParent('tmd-2H'),
    'tmd-3R': layeredParent('tmd-3R'),
    corundum: Object.assign(layeredParent('corundum'), { xFraction: 3 / 5, getRoles: getCorundumRoles }),
    cscl: bccParent('cscl', 1 / 2, 'A', getCsClRoles),
    heusler: bccParent('heusler', 1 / 4, 'X', getHeuslerRoles),
};

function getParentRoles(parent, atomTypes, positionsRed, lattice) {
    // species -> role (A, B, X) of a ternary cell read on the given parent
    let roles = parent.getRoles
        ? parent.getRoles(atomTypes, positionsRed, lattice)
        : getSuperstructureRoles(atomTypes, parent.xFraction);
    if (parent.acceptsX) {
        let x = Object.keys(roles).find((type) => roles[type] === 'X');
        if (!parent.acceptsX(x)) {
            throw new Error('The X species ' + x + ' does not fit this parent');
        }
    }
    return roles;
}

export function analyzeWurtziteDerivative(atomTypes, positionsRed, lattice, supercell = null) {
    return analyzeSuperstructure('wurtzite', atomTypes, positionsRed, lattice, supercell);
}

export function analyzeSuperstructure(parent, atomTypes, positionsRed, lattice, supercell = null) {
    /*
    template of a ternary cell, as in SUPERSTRUCTURES: a reading as a supercell
    of the parent (wurtzite, rocksalt, fluorite or antifluorite) with an atom
    of the origin role (A, X in antifluorite) as origin (and, if
    given, the supercell in that basis), the sites with their roles, lattice
    vectors reduced into the supercell, ordered by role, sublattice and position
    */
    let parentStructure = PARENT_STRUCTURES[parent];
    let roles = getParentRoles(parentStructure, atomTypes, positionsRed, lattice);
    for (let d of parentStructure.describe(atomTypes, positionsRed, lattice, (type) => roles[type] === 'X', roles)) {
        if (roles[atomTypes[d.origin]] !== parentStructure.originRole) {
            continue;
        }
        let M = d.M;
        if (supercell) {
            let U = mat.matrix_multiply(M, mat.matrix_inverse(supercell));
            if (!isIntegerMatrix(U) || Math.abs(Math.abs(mat.matrix_determinant(U)) - 1) > 1e-6) {
                continue;
            }
            M = supercell;
        }
        let inverse = mat.matrix_inverse(M);
        let sites = d.sites.map((site, i) => {
            // reduce n into the supercell: n - floor(n M^-1) M
            let f = rowTimesMatrix(site.n, inverse).map((v) => Math.floor(v + 1e-6));
            let back = rowTimesMatrix(f, M);
            return { role: roles[atomTypes[i]], sub: site.sub, n: site.n.map((v, k) => Math.round(v - back[k])) };
        });
        let rank = { A: 0, B: 1, X: 2 };
        sites.sort((p, q) => rank[p.role] - rank[q.role] || p.sub - q.sub ||
            p.n[2] - q.n[2] || p.n[1] - q.n[1] || p.n[0] - q.n[0]);
        return { supercell: M, sites: sites, a: d.a, c: d.c };
    }
    return null;
}

function orthogonalPolarFactor(K) {
    /*
    the orthogonal factor Q of K = Q P (P symmetric positive), by Newton's
    iteration X <- (X + X^-T) / 2; it keeps the sign of det K
    */
    let X = K.map((row) => row.slice());
    for (let iteration = 0; iteration < 50; iteration++) {
        let inverseT = mat.matrix_transpose(mat.matrix_inverse(X));
        let next = X.map((row, i) => row.map((v, j) => (v + inverseT[i][j]) / 2));
        let change = Math.max.apply(null, next.map((row, i) => Math.max.apply(null, row.map((v, j) => Math.abs(v - X[i][j])))));
        X = next;
        if (change < 1e-14) {
            break;
        }
    }
    return X;
}

function isIntegerMatrix(m, tol = 1e-6) {
    return m.every((row) => row.every((v) => Math.abs(v - Math.round(v)) < tol));
}

function matchSuperstructure(atomTypes, positionsRed, lattice, template, score = null) {
    /*
    map a ternary cell onto the sites of a SUPERSTRUCTURES template: the same
    supercell of the parent (any basis) with the same roles on the same
    sublattice sites modulo the supercell. The frame may need an improper
    rotation; then the canonical cell is the inversion image, D -> D^*.
    Returns the file atom of every site, the inversion flag, the proper rotation
    to the canonical axes, the canonical lattice and reduced positions, and the
    map of reduced q from the file basis to the canonical one.
    A cell matches in several equivalent ways (symmetries of the template);
    with `score` the best scoring match is returned, else the first.
    */
    let best = null;
    let bestScore = -Infinity;
    let parent = PARENT_STRUCTURES[template.parent];
    let roles = getParentRoles(parent, atomTypes, positionsRed, lattice);
    let templateInverse = mat.matrix_inverse(template.supercell);
    let reciprocal = utils.rec_lat(lattice);
    let positionsCar = utils.red_car_list(positionsRed, lattice);
    let originSite = template.sites.find((site) => site.sub === 0 && site.n.every((v) => v === 0));
    for (let d of parent.describe(atomTypes, positionsRed, lattice, (type) => roles[type] === 'X', roles)) {
        if (roles[atomTypes[d.origin]] !== originSite.role) {
            continue;
        }
        let U = mat.matrix_multiply(d.M, templateInverse);
        if (!isIntegerMatrix(U) || Math.abs(Math.abs(mat.matrix_determinant(U)) - 1) > 1e-6) {
            continue;
        }
        let order = template.sites.map((site) => {
            for (let i = 0; i < atomTypes.length; i++) {
                if (d.sites[i].sub !== site.sub || roles[atomTypes[i]] !== site.role) {
                    continue;
                }
                let shift = [0, 1, 2].map((k) => d.sites[i].n[k] - site.n[k]);
                if (rowTimesMatrix(shift, templateInverse).every((v) => Math.abs(v - Math.round(v)) < 1e-6)) {
                    return i;
                }
            }
            return -1;
        });
        if (order.some((i) => i < 0) || new Set(order).size !== order.length) {
            continue;
        }

        // refine the frame: the orthogonal T (same handedness as the one from the
        // bonds) that best maps the cell, in the template basis, onto the ideal
        // supercell of the parent (orthogonal Procrustes), so that e.g. an
        // orthorhombic cell lies exactly along the canonical axes
        let inverseU = mat.matrix_inverse(U).map((row) => row.map(Math.round));
        let basis = mat.matrix_multiply(inverseU, lattice);
        let ideal = mat.matrix_multiply(template.supercell, d.H);
        let T = orthogonalPolarFactor(mat.matrix_multiply(mat.matrix_transpose(ideal), basis));
        if (T.some((row, i) => mat.vec_dot(row, d.T[i]) < Math.cos(parent.frameTolerance * Math.PI / 180))) {
            continue;
        }
        let invert = mat.matrix_determinant(T) < 0;
        let R = invert ? T.map((row) => row.map((v) => -v)) : T;
        let rotated = lattice.map((row) => matrixTimesVector(T, row));
        let canonicalLattice = mat.matrix_multiply(inverseU, rotated);
        let canonicalInverse = mat.matrix_inverse(canonicalLattice);
        let positions = template.sites.map((site, s) => {
            let i = order[s];
            let cartesian = matrixTimesVector(T, [0, 1, 2].map((k) => positionsCar[i][k] - positionsCar[d.origin][k]));
            let reduced = rowTimesMatrix(cartesian, canonicalInverse);
            let ideal = rowTimesMatrix(parent.sitePosition(site), templateInverse);
            return reduced.map((v, k) => v - Math.round(v - ideal[k]));
        });
        let identity = R.every((row, k) => Math.abs(row[k] - 1) < 1e-9);
        // reduced q: file basis -> cartesian -> rotated by R -> canonical basis
        let qMap = (q) => {
            let qCart = rowTimesMatrix(q, reciprocal);
            let rotatedQ = matrixTimesVector(R, qCart);
            return canonicalLattice.map((row) => mat.vec_dot(rotatedQ, row));
        };
        // the map is linear, q_canonical = q_file Q: rows of Q are the images of the unit vectors
        let qMapInverse = mat.matrix_inverse([[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(qMap));
        let result = {
            order: order,
            invert: invert,
            rotation: identity ? null : R,
            structure: { lattice: canonicalLattice, positionsRed: positions },
            qMap: qMap,
            qMapInverse: qMapInverse,
        };
        if (!score) {
            return result;
        }
        let value = score(result);
        if (value > bestScore) {
            best = result;
            bestScore = value;
        }
    }
    if (best) {
        return best;
    }
    throw new Error('Could not match the cell onto the ' + template.title + ' template');
}

function referencePointScore(segments, template) {
    /*
    files without force constants can only be sampled on their own path, and two
    files may run the same path in orientations related by a symmetry of the
    crystal (L at (0,-1/2,0) in one, (0,0,1/2) in the other). Of the equivalent
    matches, prefer the one that puts the labelled points of the file where the
    template has them (modulo the reciprocal lattice), so the paths line up. On a
    tie, prefer the points exactly there (signs included): two files whose points
    agree only modulo the reciprocal lattice can still run their segments on lines
    related by q -> -q (Z-R from (0,0,1/2) in one, from (0,0,-1/2) in the other).
    */
    if (!template.referencePoints || !segments.length) {
        return null;
    }
    return (match) => {
        let count = 0;
        let exact = 0;
        for (let segment of segments) {
            for (let [label, q] of [[segment.labelA, segment.qA], [segment.labelB, segment.qB]]) {
                let reference = template.referencePoints[normalizeLabel(label)];
                if (!reference) {
                    continue;
                }
                let mapped = match.qMap(q);
                if (mapped.every((v, k) => Math.abs(v - reference[k] - Math.round(v - reference[k])) < 1e-3)) {
                    count += 1;
                    if (mapped.every((v, k) => Math.abs(v - reference[k]) < 1e-3)) {
                        exact += 1;
                    }
                }
            }
        }
        // exact hits only break ties: fewer than 1000 labelled points per path
        return count + exact / 1000;
    };
}

function matchScore(segments, template) {
    /*
    how to choose among the equivalent matches of a cell (matchSuperstructure):
    first the distortion closest to the template's referenceStructure (positions
    in parent units and the cell shape at unit volume), then, among matches equal
    up to a symmetry of the crystal, the labelled points (referencePointScore)
    */
    let pointScore = referencePointScore(segments, template);
    let reference = template.referenceStructure;
    if (!reference) {
        return pointScore;
    }
    return (match) => {
        let lattice = match.structure.lattice;
        let scale = Math.cbrt(Math.abs(mat.matrix_determinant(lattice)));
        let deviation = 0;
        lattice.forEach((row, i) => row.forEach((v, j) => {
            deviation += (v / scale - reference.lattice[i][j]) ** 2;
        }));
        match.structure.positionsRed.forEach((position, s) => {
            let d = rowTimesMatrix(position.map((v, k) => v - reference.positions[s][k]), template.supercell);
            deviation += mat.vec_dot(d, d);
        });
        return -Math.round(deviation * 1e4) * 1e6 + (pointScore ? pointScore(match) : 0);
    };
}

function getSuperstructureFamilies(atomTypes, positionsRed, lattice) {
    // the SUPERSTRUCTURES entries with the same number of atoms per role
    let count = (list) => ['A', 'B', 'X'].map((role) => list.filter((r) => r === role).length).join(',');
    return Object.keys(SUPERSTRUCTURES).filter((family) => {
        let template = SUPERSTRUCTURES[family];
        let roles;
        try {
            roles = getParentRoles(PARENT_STRUCTURES[template.parent], atomTypes, positionsRed, lattice);
        } catch (error) {
            return false;
        }
        return count(atomTypes.map((type) => roles[type])) === count(template.sites.map((site) => site.role));
    });
}

function templateParentPosition(site, u) {
    // wurtzite-basis coordinates of a template site for the anion parameter u
    let base = WURTZITE_SITES[site.sub].red;
    return [0, 1, 2].map((k) => site.n[k] + base[k] + (k === 2 && site.role === 'X' ? u : 0));
}

function relabelSegments(segments, specialPoints) {
    /*
    name segment ends by their coordinates; the PhononDB files carry seekpath
    names and occasionally a misplaced label
    */
    let name = (q, fallback) => {
        for (let point of specialPoints) {
            if (point.q.every((v, i) => Math.abs(v - q[i]) < 2e-3)) {
                return point.label;
            }
        }
        // an equivalent point: the same modulo the reciprocal lattice, or its inverse
        for (let point of specialPoints) {
            for (let sign of [1, -1]) {
                if (point.q.every((v, i) => Math.abs(wrapFraction(v - sign * q[i])) < 2e-3)) {
                    return point.label;
                }
            }
        }
        return fallback;
    };
    segments.forEach((segment) => {
        segment.labelA = name(segment.qA, segment.labelA);
        segment.labelB = name(segment.qB, segment.labelB);
    });
}

function getSegmentsFromLabels(qpoints, labelPoints) {
    /*
    split a path into high-symmetry segments. A segment runs between two
    consecutive q-points lying on labelled points that are not adjacent in the
    list (adjacent labelled points are either a duplicated corner or a jump).
    */
    let labelled = [];
    for (let i = 0; i < qpoints.length; i++) {
        for (let j = 0; j < labelPoints.length; j++) {
            if (sameQ(qpoints[i], labelPoints[j].q)) {
                labelled.push({ index: i, label: labelPoints[j].label });
                break;
            }
        }
    }

    let segments = [];
    for (let k = 0; k + 1 < labelled.length; k++) {
        let start = labelled[k];
        let end = labelled[k + 1];
        if (end.index - start.index < 2) {
            continue;
        }
        segments.push({
            start: start.index,
            end: end.index + 1,
            labelA: start.label,
            labelB: end.label,
        });
    }
    return segments;
}

function makeSegment(qpoints, segment, sampleAt) {
    let qA = qpoints[segment.start];
    let qB = qpoints[segment.end - 1];
    return {
        qA: qA,
        qB: qB,
        labelA: segment.labelA,
        labelB: segment.labelB,
        npoints: segment.end - segment.start,
        sample: sampleAt,
    };
}

function buildPhononDBEndpoint(raw, entry) {
    let payload = raw.dynamical_matrix;
    if (!payload.primitive_lattice) {
        payload.primitive_lattice = raw.lattice;
    }
    let factor = Number(payload.frequency_conversion_factor) || 1.0;
    let factorSq = factor * factor;
    let qpoints = raw.qpoints;

    let nac = null;
    let shortRangePayload = null;
    if (payload.nac && payload.nac.method === 'gonze') {
        let volume = Math.abs(mat.matrix_determinant(payload.primitive_lattice));
        nac = {
            born: payload.nac.born,
            dielectric: payload.nac.dielectric,
            // phonopy: nac_factor = unit_conversion * 4 pi / volume
            unitConversion: payload.nac.nac_factor * volume / (4 * Math.PI),
            tolerance: payload.nac.q_direction_tolerance || 1e-5,
        };
        shortRangePayload = Object.assign({}, payload, { nac: null });
    }

    // each line break is one straight segment; the labels are only stored at the
    // segment ends, the start of a segment carries the label of the previous end
    let labelAtIndex = {};
    for (let i = 0; i < raw.highsym_qpts.length; i++) {
        labelAtIndex[raw.highsym_qpts[i][0]] = normalizeLabel(raw.highsym_qpts[i][1]);
    }

    let segments = [];
    let lineBreaks = raw.line_breaks || [[0, qpoints.length]];
    for (let i = 0; i < lineBreaks.length; i++) {
        let segment = {
            start: lineBreaks[i][0],
            end: lineBreaks[i][1],
        };
        segment.labelA = labelAtIndex[segment.start] !== undefined
            ? labelAtIndex[segment.start]
            : labelAtIndex[segment.start - 1];
        segment.labelB = labelAtIndex[segment.end - 1];
        let qA = qpoints[segment.start];
        let qB = qpoints[segment.end - 1];
        let madeSegment = makeSegment(qpoints, segment, function(t) {
            let q = lerpQ(qA, qB, t);
            // non-analytic term at Gamma: approach along the segment
            let qDirection = isGamma(q) ? [qB[0] - qA[0], qB[1] - qA[1], qB[2] - qA[2]] : null;
            return scaleMatrix(buildDynamicalMatrixBlocks(payload, q, qDirection), factorSq);
        });
        if (shortRangePayload) {
            // the stored force constants are the Gonze short-range part: without
            // the dipole-dipole term they give the analytic remainder of D(q)
            madeSegment.sampleShortRange = function(t) {
                return scaleMatrix(buildDynamicalMatrixBlocks(shortRangePayload, lerpQ(qA, qB, t)), factorSq);
            };
        }
        segments.push(madeSegment);
    }

    return {
        atomTypes: raw.atom_types.slice(),
        positionsRed: raw.atom_pos_red,
        lattice: raw.lattice,
        masses: payload.masses.slice(),
        segments: segments,
        nac: nac,
        frequencyFactor: factor,
        // force constants give D at any q (reduced, file frame), not only on the path;
        // at q equivalent to Gamma the non-analytic term needs the approach direction
        evaluate: function(q, direction) {
            let atGamma = q.every((v) => Math.abs(v - Math.round(v)) < QTOL);
            return scaleMatrix(buildDynamicalMatrixBlocks(payload, q, atGamma ? direction : null), factorSq);
        },
        evaluateShortRange: shortRangePayload ? function(q) {
            return scaleMatrix(buildDynamicalMatrixBlocks(shortRangePayload, q), factorSq);
        } : null,
    };
}

function buildMaterialsProjectEndpoint(raw, dielectric) {
    let structure = raw.structure;
    let sites = structure.sites;
    let atomTypes = sites.map((site) => site.label);
    let positionsRed = sites.map((site) => site.abc);
    let masses = atomTypes.map((type) => atomic_data.atomic_mass[atomic_data.atomic_number[type]]);
    let natoms = sites.length;
    let size = natoms * 3;
    let qpoints = raw.qpoints;
    let frequencies = raw.frequencies;
    let eigendisplacements = raw.eigendisplacements;
    let nbands = frequencies.length;

    // D(q) at the q-points of the file, in the atom gauge
    let matrices = new Array(qpoints.length).fill(null);
    let getMatrix = function(qIndex) {
        if (matrices[qIndex]) {
            return matrices[qIndex];
        }
        let matrix = zeroMatrix(size);
        let q = qpoints[qIndex];
        for (let n = 0; n < nbands; n++) {
            // eigendisplacements are u = e / sqrt(m) up to a constant: recover e
            let vector = [];
            let normSq = 0;
            for (let a = 0; a < natoms; a++) {
                let sqrtMass = Math.sqrt(masses[a]);
                for (let d = 0; d < 3; d++) {
                    let c = parseComplex(eigendisplacements[n][qIndex][a][d]);
                    let re = c[0] * sqrtMass;
                    let im = c[1] * sqrtMass;
                    vector.push([re, im]);
                    normSq += re * re + im * im;
                }
            }
            if (!(normSq > 0)) {
                continue;
            }
            let nu = frequencies[n][qIndex];
            let weight = (nu < 0 ? -nu * nu : nu * nu) / normSq;
            for (let i = 0; i < size; i++) {
                for (let j = 0; j < size; j++) {
                    // e_i conj(e_j)
                    let re = vector[i][0] * vector[j][0] + vector[i][1] * vector[j][1];
                    let im = vector[i][1] * vector[j][0] - vector[i][0] * vector[j][1];
                    matrix.real[i][j] += weight * re;
                    matrix.imag[i][j] += weight * im;
                }
            }
        }

        // cell gauge -> atom gauge: D_ij *= exp(2 pi i q.(tau_j - tau_i))
        for (let a = 0; a < natoms; a++) {
            for (let b = 0; b < natoms; b++) {
                let phase = 2 * Math.PI * (
                    q[0] * (positionsRed[b][0] - positionsRed[a][0]) +
                    q[1] * (positionsRed[b][1] - positionsRed[a][1]) +
                    q[2] * (positionsRed[b][2] - positionsRed[a][2])
                );
                let cos = Math.cos(phase);
                let sin = Math.sin(phase);
                for (let da = 0; da < 3; da++) {
                    for (let db = 0; db < 3; db++) {
                        let i = a * 3 + da;
                        let j = b * 3 + db;
                        let re = matrix.real[i][j];
                        let im = matrix.imag[i][j];
                        matrix.real[i][j] = re * cos - im * sin;
                        matrix.imag[i][j] = re * sin + im * cos;
                    }
                }
            }
        }

        matrices[qIndex] = matrix;
        return matrix;
    };

    // with the Born charges and dielectric tensor of the same calculation, split
    // D into the dipole-dipole term and the analytic (short-range) remainder
    let lattice = structure.lattice.matrix;
    let nac = null;
    let getShortRangeMatrix = null;
    // the collection lists the sites in the order of the band-structure file; atoms
    // of the same element can carry different (rotated) Born tensors, so match by
    // index, and by element only when every element occurs once
    let sameSites = dielectric && dielectric.labels.length === atomTypes.length &&
        dielectric.labels.every((label, i) => label === atomTypes[i]);
    let uniqueSpecies = new Set(atomTypes).size === atomTypes.length;
    if (dielectric && !sameSites && !uniqueSpecies) {
        console.log("Materials Project Born charges listed in another site order, interpolating the LO-TO term");
        dielectric = null;
    }
    if (dielectric) {
        let born = atomTypes.map((type, i) => {
            let index = sameSites ? i : dielectric.labels.indexOf(type);
            if (index < 0) {
                throw new Error('No Born charge for ' + type);
            }
            return dielectric.born[index];
        });
        nac = {
            born: born,
            dielectric: dielectric.dielectric,
            unitConversion: PHONOPY_NAC_UNIT,
            tolerance: 1e-5,
        };
        let ownPayload = buildNacPayload(
            born, dielectric.dielectric, PHONOPY_NAC_UNIT,
            lattice, masses, utils.red_car_list(positionsRed, lattice), nac.tolerance
        );
        let shortRange = new Array(qpoints.length).fill(null);
        getShortRangeMatrix = function(qIndex, direction) {
            if (!shortRange[qIndex]) {
                let q = qpoints[qIndex];
                let matrix = scaleMatrix(getMatrix(qIndex), 1);
                // at Gamma the file holds the limit taken along this segment
                let dd = getGonzeReciprocalCorrection(ownPayload, q, isGamma(q) ? direction : null);
                addBlockTensor(matrix, dd, -PHONOPY_THZ * PHONOPY_THZ);
                shortRange[qIndex] = matrix;
            }
            return shortRange[qIndex];
        };
    }

    let labelPoints = Object.keys(raw.labels_dict).map((label) => ({
        q: raw.labels_dict[label],
        label: normalizeLabel(label),
    }));

    let segments = getSegmentsFromLabels(qpoints, labelPoints).map(function(segment) {
        let qA = qpoints[segment.start];
        let qB = qpoints[segment.end - 1];
        let length = mat.distance(qA, qB);
        let fractions = [];
        for (let k = segment.start; k < segment.end; k++) {
            fractions.push(mat.distance(qpoints[k], qA) / length);
        }
        let direction = [qB[0] - qA[0], qB[1] - qA[1], qB[2] - qA[2]];
        // linear interpolation of D between the sampled q-points of the file
        let interpolate = function(getter, t) {
            let k = 0;
            while (k < fractions.length - 2 && fractions[k + 1] < t) {
                k++;
            }
            let t0 = fractions[k];
            let t1 = fractions[k + 1];
            let w = t1 > t0 ? Math.min(1, Math.max(0, (t - t0) / (t1 - t0))) : 0;
            return lerpMatrix(getter(segment.start + k, direction), getter(segment.start + k + 1, direction), w);
        };
        let madeSegment = makeSegment(qpoints, segment, (t) => interpolate(getMatrix, t));
        if (getShortRangeMatrix) {
            madeSegment.sampleShortRange = (t) => interpolate(getShortRangeMatrix, t);
        }
        return madeSegment;
    });

    return {
        atomTypes: atomTypes,
        positionsRed: positionsRed,
        lattice: lattice,
        masses: masses,
        segments: segments,
        nac: nac,
        frequencyFactor: PHONOPY_THZ,
    };
}

export function buildMixingEndpoint(raw, entry, dielectric = null) {
    /*
    raw is the parsed json of a PhononDB (internal format with a dynamical_matrix)
    or a Materials Project OpenData phonon file; entry is its menu entry.
    dielectric ({ born, dielectric, labels }, see mpdielectric.js) supplies the
    Born charges and dielectric tensor of a Materials Project material.
    */
    let endpoint;
    if (raw.dynamical_matrix) {
        endpoint = buildPhononDBEndpoint(raw, entry);
    } else if (raw.eigendisplacements && raw.structure) {
        endpoint = buildMaterialsProjectEndpoint(raw, dielectric);
    } else {
        throw new Error('No dynamical matrix available for ' + (entry && entry.name));
    }
    let canonical;
    let natoms = endpoint.atomTypes.length;
    let nspecies = new Set(endpoint.atomTypes).size;
    // the menu entry names the family of an ordered ternary; several share the
    // composition and cell size (e.g. 8-atom chalcopyrite and gamma-LiFeO2 type)
    let hint = entry && SUPERSTRUCTURES[entry.family] ? entry.family : null;
    let chalcopyriteError = null;
    if (natoms === 8 && !hint) {
        try {
            canonical = matchChalcopyriteSites(endpoint.atomTypes, endpoint.positionsRed, endpoint.lattice);
        } catch (error) {
            chalcopyriteError = error;
        }
    }
    // a 2-atom cell that is neither zinc blende nor rocksalt (CsCl), a 4-atom one
    // that is not wurtzite (Li3Sb) is read on the superstructure templates
    let twoAtom = null;
    let fourAtom = null;
    if (natoms === 2 && !hint) {
        try {
            twoAtom = getCanonicalSites(endpoint.atomTypes, endpoint.positionsRed,
                getCationSymbol(entry && entry.name, endpoint.atomTypes), endpoint.lattice);
        } catch (error) {
            chalcopyriteError = error;
        }
    }
    if (natoms === 4 && nspecies <= 2 && !hint) {
        try {
            fourAtom = matchWurtziteSites(endpoint.atomTypes, endpoint.positionsRed, endpoint.lattice,
                getCationSymbol(entry && entry.name, endpoint.atomTypes));
        } catch (error) {
            chalcopyriteError = error;
        }
    }
    if (twoAtom) {
        canonical = twoAtom;
        endpoint.family = canonical.family;
        endpoint.rotation = null;
        if (canonical.family === 'rocksalt') {
            endpoint.structure = { a: Math.SQRT2 * Math.sqrt(mat.vec_dot(endpoint.lattice[0], endpoint.lattice[0])) };
        }
    } else if (natoms === 3 && nspecies === 2 && !hint && isFluoriteCell(endpoint)) {
        canonical = getFluoriteSites(endpoint.atomTypes, endpoint.positionsRed, endpoint.lattice);
        endpoint.family = canonical.family;
        endpoint.rotation = null;
        endpoint.structure = { a: Math.SQRT2 * Math.sqrt(mat.vec_dot(endpoint.lattice[0], endpoint.lattice[0])) };
    } else if (fourAtom) {
        canonical = fourAtom;
        endpoint.family = 'wurtzite';
        endpoint.rotation = canonical.rotation;
        endpoint.structure = canonical.structure;
        relabelSegments(endpoint.segments, HEXAGONAL_SPECIAL_POINTS);
    } else if (canonical) {
        endpoint.family = 'chalcopyrite';
        endpoint.rotation = canonical.rotation;
        endpoint.structure = canonical.structure;
        relabelSegments(endpoint.segments, getChalcopyriteSpecialPoints(canonical.structure.a, canonical.structure.c));
    } else {
        let families = hint ? [hint] : getSuperstructureFamilies(endpoint.atomTypes, endpoint.positionsRed, endpoint.lattice);
        let family = null;
        let lastError = chalcopyriteError;
        for (let candidate of families) {
            try {
                canonical = matchSuperstructure(endpoint.atomTypes, endpoint.positionsRed, endpoint.lattice,
                    SUPERSTRUCTURES[candidate], matchScore(endpoint.segments, SUPERSTRUCTURES[candidate]));
                family = candidate;
                break;
            } catch (error) {
                lastError = error;
            }
        }
        if (!family) {
            throw lastError || new Error('Mixing supports 2-atom zinc blende, rocksalt and CsCl, 3-atom fluorite, ' +
                '4-atom wurtzite and Li3Bi type, 5-atom perovskite, 8-atom chalcopyrite and ternary wurtzite, rocksalt, ' +
                'fluorite, perovskite and Heusler cells');
        }
        let template = SUPERSTRUCTURES[family];
        endpoint.family = family;
        endpoint.rotation = canonical.rotation;
        endpoint.structure = canonical.structure;
        // the canonical cell may use another basis than the file: path in canonical reduced q
        endpoint.segments.forEach((segment) => {
            segment.qA = canonical.qMap(segment.qA);
            segment.qB = canonical.qMap(segment.qB);
        });
        if (template.specialPoints) {
            relabelSegments(endpoint.segments, template.specialPoints);
        }
        if (endpoint.evaluate) {
            // force constants: D (file frame, see toCanonical) at any canonical reduced q,
            // so that both materials are sampled on the same canonical path; the files
            // run their paths along different but equivalent directions
            let toFile = (q) => rowTimesMatrix(q, canonical.qMapInverse);
            let evaluate = endpoint.evaluate;
            let evaluateShortRange = endpoint.evaluateShortRange;
            endpoint.sampleAtQ = (q, direction) => evaluate(toFile(q), direction ? toFile(direction) : null);
            if (evaluateShortRange) {
                endpoint.sampleShortRangeAtQ = (q) => evaluateShortRange(toFile(q));
            }
        }
    }
    // some files put a label on the wrong point (PhononDB mp-5840 has "X" at q = 0
    // and "G" off Gamma): a reciprocal lattice vector is Gamma, and Gamma is
    // nowhere else
    let isLatticeVector = (q) => q.every((v) => Math.abs(v - Math.round(v)) < QTOL);
    let fixLabel = (label, q) => (isLatticeVector(q) ? 'GAMMA' : label === 'GAMMA' ? '' : label);
    endpoint.segments.forEach((segment) => {
        segment.labelA = fixLabel(segment.labelA, segment.qA);
        segment.labelB = fixLabel(segment.labelB, segment.qB);
    });
    endpoint.name = entry && entry.name ? entry.name : utils.get_formula(endpoint.atomTypes);
    // how the material is drawn on its own (the mixture at x = 0 is drawn the same way)
    endpoint.repetitions = raw.repetitions || [3, 3, 3];
    endpoint.order = canonical.order;
    endpoint.invert = canonical.invert;
    endpoint.siteTypes = canonical.order.map((i) => endpoint.atomTypes[i]);
    endpoint.siteMasses = canonical.order.map((i) => endpoint.masses[i]);
    if (endpoint.nac) {
        // Born charges and the dielectric tensor are rank-2 tensors: unchanged by the
        // inversion, rotated to the conventional axes like the blocks of D
        endpoint.nac.siteBorn = canonical.order.map((i) => rotateTensor(endpoint.rotation, endpoint.nac.born[i]));
        endpoint.nac.siteDielectric = rotateTensor(endpoint.rotation, endpoint.nac.dielectric);
    }
    return endpoint;
}

function rotateTensor(rotation, tensor) {
    /*
    R T R^T for a 3x3 tensor; rotation null means identity
    */
    if (!rotation) {
        return tensor;
    }
    return [0, 1, 2].map((a) => [0, 1, 2].map((b) => {
        let value = 0;
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                value += rotation[a][i] * tensor[i][j] * rotation[b][j];
            }
        }
        return value;
    }));
}

function toCanonical(endpoint, matrix) {
    /*
    D of the file -> canonical site order, conventional axes, and the complex
    conjugate for a cell stored as the inversion image of the template
    */
    let size = matrix.real.length;
    let index = [];
    for (let s = 0; s < endpoint.order.length; s++) {
        for (let d = 0; d < 3; d++) {
            index.push(endpoint.order[s] * 3 + d);
        }
    }
    let sign = endpoint.invert ? -1 : 1;
    let out = zeroMatrix(size);
    for (let i = 0; i < size; i++) {
        for (let j = 0; j < size; j++) {
            out.real[i][j] = matrix.real[index[i]][index[j]];
            out.imag[i][j] = sign * matrix.imag[index[i]][index[j]];
        }
    }
    let rotation = endpoint.rotation;
    if (rotation) {
        let natoms = size / 3;
        for (let p = 0; p < natoms; p++) {
            for (let r = 0; r < natoms; r++) {
                for (let part of ['real', 'imag']) {
                    let block = [0, 1, 2].map((a) => [0, 1, 2].map((b) => out[part][p * 3 + a][r * 3 + b]));
                    let rotated = rotateTensor(rotation, block);
                    for (let a = 0; a < 3; a++) {
                        for (let b = 0; b < 3; b++) {
                            out[part][p * 3 + a][r * 3 + b] = rotated[a][b];
                        }
                    }
                }
            }
        }
    }
    return out;
}

function foldIntoSupercell(primitive, family, superLattice, superPositionsRed, sublattice, structure) {
    /*
    D of a supercell of `primitive` by band folding, so that a binary can be
    mixed with an ordered ternary on the same lattice:

        D^S_IJ(Q) = 1/N sum_g D_ij(Q + g) exp(-2 pi i g.(r_J - r_I))

    over the N supercell reciprocal vectors g that differ modulo the primitive
    reciprocal lattice (i, j: sublattices of the supercell sites I, J; atom gauge,
    D_ij(k) = sum_R Phi_ij(R) exp(2 pi i k.(R + tau_j - tau_i))). With this sign
    the sum does not depend on which representative g of each class is used.
    The primitive cell is taken in its canonical frame (no rotation): the file
    lattice is the canonical one. Q + g leaves the band path, so this needs
    force constants (PhononDB).
    */
    let natoms = sublattice.length;
    let size = natoms * 3;
    let count = Math.round(natoms / primitive.siteTypes.length);
    let primitiveLattice = primitive.lattice;
    let positionsCar = utils.red_car_list(superPositionsRed, superLattice);
    let superReciprocal = utils.rec_lat(superLattice);

    // folding vectors: supercell reciprocal vectors, distinct modulo the primitive
    // reciprocal lattice; the shortest of each class, so that Q = 0 maps onto k = 0 exactly
    let candidates = [];
    for (let n0 = -3; n0 <= 3; n0++) {
        for (let n1 = -3; n1 <= 3; n1++) {
            for (let n2 = -3; n2 <= 3; n2++) {
                candidates.push([0, 1, 2].map((d) => n0 * superReciprocal[0][d] + n1 * superReciprocal[1][d] + n2 * superReciprocal[2][d]));
            }
        }
    }
    candidates.sort((g, h) => mat.vec_dot(g, g) - mat.vec_dot(h, h));
    let folding = [];
    let seen = new Set();
    for (let g of candidates) {
        let key = primitiveLattice.map((p) => {
            let f = mat.vec_dot(g, p);
            return Math.round((f - Math.floor(f + 1e-9)) * 1e6) % 1000000;
        }).join(',');
        if (!seen.has(key)) {
            seen.add(key);
            folding.push(g);
        }
    }
    if (folding.length !== count) {
        throw new Error('Expected ' + count + ' folding vectors, found ' + folding.length);
    }

    let toCartesian = (q) => [0, 1, 2].map((d) => q[0] * superReciprocal[0][d] + q[1] * superReciprocal[1][d] + q[2] * superReciprocal[2][d]);
    let toPrimitive = (k) => primitiveLattice.map((p) => mat.vec_dot(k, p));
    let sampleAtQ = function(q, direction) {
        let qCart = toCartesian(q);
        let directionPrimitive = direction ? toPrimitive(toCartesian(direction)) : null;
        let out = zeroMatrix(size);
        for (let g of folding) {
            let k = [0, 1, 2].map((d) => qCart[d] + g[d]);
            let block = toCanonical(primitive, primitive.evaluate(toPrimitive(k), directionPrimitive));
            for (let I = 0; I < natoms; I++) {
                for (let J = 0; J < natoms; J++) {
                    let phase = -2 * Math.PI * (
                        g[0] * (positionsCar[J][0] - positionsCar[I][0]) +
                        g[1] * (positionsCar[J][1] - positionsCar[I][1]) +
                        g[2] * (positionsCar[J][2] - positionsCar[I][2])
                    );
                    let cos = Math.cos(phase) / count;
                    let sin = Math.sin(phase) / count;
                    for (let p = 0; p < 3; p++) {
                        for (let r = 0; r < 3; r++) {
                            let re = block.real[sublattice[I] * 3 + p][sublattice[J] * 3 + r];
                            let im = block.imag[sublattice[I] * 3 + p][sublattice[J] * 3 + r];
                            out.real[I * 3 + p][J * 3 + r] += re * cos - im * sin;
                            out.imag[I * 3 + p][J * 3 + r] += re * sin + im * cos;
                        }
                    }
                }
            }
        }
        return out;
    };

    let siteTypes = sublattice.map((s) => primitive.siteTypes[s]);
    let siteMasses = sublattice.map((s) => primitive.siteMasses[s]);
    let nac = null;
    let sampleShortRangeAtQ = null;
    if (primitive.nac) {
        nac = {
            siteBorn: sublattice.map((s) => primitive.nac.siteBorn[s]),
            siteDielectric: primitive.nac.siteDielectric,
            unitConversion: primitive.nac.unitConversion,
            tolerance: primitive.nac.tolerance,
        };
        // short-range part: folded D minus the dipole-dipole term of the supercell built
        // with the same rules as the mixture's, so the end-member is reproduced exactly
        let ownPayload = buildNacPayload(nac.siteBorn, nac.siteDielectric, nac.unitConversion,
            superLattice, siteMasses, positionsCar, nac.tolerance);
        let factorSq = primitive.frequencyFactor * primitive.frequencyFactor;
        sampleShortRangeAtQ = function(q, direction) {
            let matrix = sampleAtQ(q, direction);
            let dd = getGonzeReciprocalCorrection(ownPayload, q, isGamma(q) ? direction : null);
            addBlockTensor(matrix, dd, -factorSq);
            return matrix;
        };
    }

    return {
        family: family,
        anyQ: true,
        // the binary it was folded from and the sublattice of every site, to draw the
        // supercell in the binary's own setting (presentInFileFrame)
        foldedFrom: primitive,
        sublattice: sublattice,
        name: primitive.name,
        lattice: superLattice,
        structure: structure,
        siteTypes: siteTypes,
        siteMasses: siteMasses,
        order: siteTypes.map((_, i) => i),
        invert: false,
        rotation: null,
        nac: nac,
        frequencyFactor: primitive.frequencyFactor,
        segments: [],
        sampleAtQ: sampleAtQ,
        sampleShortRangeAtQ: sampleShortRangeAtQ,
    };
}

export function zincblendeAsChalcopyrite(zb) {
    /*
    A zinc blende expressed in the 8-atom chalcopyrite cell, so it can be mixed
    with a chalcopyrite. With c = 2a and u = 1/4 the chalcopyrite template is
    zinc blende: the A and B sites lie on the cation sublattice, the anion sites
    on the anion sublattice, and the body-centred tetragonal primitive cell holds
    4 zinc blende primitive cells (see foldIntoSupercell).
    */
    if (!zb.evaluate) {
        throw new Error('Mixing a zinc blende with a chalcopyrite needs its force constants (PhononDB); ' +
            'the Materials Project files only cover the band path');
    }
    let primitive = zb.lattice;
    let a = Math.SQRT2 * Math.sqrt(mat.vec_dot(primitive[0], primitive[0]));
    let fcc = [[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]];
    if (!primitive.every((row, i) => row.every((v, j) => Math.abs(v - fcc[i][j]) < 1e-6 * a))) {
        throw new Error('Unsupported zinc blende setting: not the standard fcc primitive cell');
    }
    let c = 2 * a;
    let u = 0.25;
    // template sites 0-3 (A, B) are cations, 4-7 anions
    let sublattice = CHALCOPYRITE_SITES.map((site) => (site.role === 'X' ? 1 : 0));
    return foldIntoSupercell(zb, 'chalcopyrite', chalcopyriteLattice(a, c), chalcopyritePositions(a, c, u),
        sublattice, { a: a, c: c, u: u });
}

export function wurtziteAsSuperstructure(wz, family) {
    return binaryAsSuperstructure(wz, family);
}

export function binaryAsSuperstructure(binary, family) {
    /*
    A binary wurtzite, rocksalt, fluorite or antifluorite expressed in the cell of an ordered ternary on
    its lattice (SUPERSTRUCTURES), every site on its own sublattice of the
    parent, so it can be mixed with that ternary (see foldIntoSupercell)
    */
    let template = SUPERSTRUCTURES[family];
    if (!binary.evaluate) {
        throw new Error('Mixing a ' + template.parent + ' with a ternary ' + template.parent +
            ' needs its force constants (PhononDB); the Materials Project files only cover the band path');
    }
    if (SUPERSTRUCTURES[binary.family] && binary.structure && binary.structure.lattice) {
        // a 5-atom perovskite read by matchSuperstructure: fold its canonical cell (the
        // file may be rotated or in another setting), evaluated at canonical q; its
        // sites are the parent sublattices in order
        let cell = binary.structure.lattice;
        let primitive = Object.assign({}, binary, { lattice: cell, evaluate: binary.sampleAtQ });
        let lattice = mat.matrix_multiply(template.supercell, cell);
        let inverse = mat.matrix_inverse(template.supercell);
        let positionsRed = template.sites.map((site) => rowTimesMatrix(
            binary.structure.positionsRed[site.sub].map((v, k) => v + site.n[k]), inverse));
        let folded = foldIntoSupercell(primitive, family, lattice, positionsRed, template.sites.map((site) => site.sub),
            { lattice: lattice, positionsRed: positionsRed });
        // the proxy holds the canonical cell; the file setting is the binary's
        folded.foldedFrom = binary;
        return folded;
    }
    if (FCC_LAYOUTS[template.parent]) {
        // rocksalt, fluorite, antifluorite: the binary sites are the sublattices in order
        let a = binary.structure.a;
        let fcc = rocksaltLattice(a);
        if (binary.rotation || !binary.lattice.every((row, i) => row.every((v, j) => Math.abs(v - fcc[i][j]) < 1e-6 * a))) {
            throw new Error('Unsupported ' + template.parent + ' setting: not the standard fcc primitive cell');
        }
        let lattice = mat.matrix_multiply(template.supercell, fcc);
        let inverse = mat.matrix_inverse(template.supercell);
        let parent = PARENT_STRUCTURES[template.parent];
        let positionsRed = template.sites.map((site) => rowTimesMatrix(parent.sitePosition(site), inverse));
        return foldIntoSupercell(binary, family, lattice, positionsRed, template.sites.map((site) => site.sub),
            { lattice: lattice, positionsRed: positionsRed });
    }
    let wz = binary;
    let { a, c, u } = wz.structure;
    let hexagonal = wurtziteLattice(a, c);
    if (wz.rotation || !wz.lattice.every((row, i) => row.every((v, j) => Math.abs(v - hexagonal[i][j]) < 1e-6 * a))) {
        throw new Error('Unsupported wurtzite setting: not the standard hexagonal cell');
    }
    let lattice = mat.matrix_multiply(template.supercell, hexagonal);
    let inverse = mat.matrix_inverse(template.supercell);
    let positionsRed = template.sites.map((site) => rowTimesMatrix(templateParentPosition(site, u), inverse));
    return foldIntoSupercell(wz, family, lattice, positionsRed, template.sites.map((site) => site.sub),
        { lattice: lattice, positionsRed: positionsRed });
}

// Lattice worlds: one parent cell whose sites contain those of several families, so
// that two materials on the same lattice that fill different sites, or the same
// sites with other roles, can still be mixed site by site. A site filled in only one
// end-member is a vacancy in the other (computeUnionMixedPhonon).
// fcc: rocksaltLattice(a) with the fcc site, the octahedral site and the two
//      tetrahedral sites (HEUSLER_SITES 0, 1/2, 1/4, 3/4): zinc blende and
//      chalcopyrite fill 0 and 1/4, rocksalt 0 and 1/2, fluorite and antifluorite 0,
//      1/4 and 3/4, the Li3Bi type and the Heuslers all four, CsCl (folded into the
//      doubled cell) all four.
// cubic: the simple cubic cell a I with the corner, the centre and the three face
//      centres (PEROVSKITE_SITES): perovskite and antiperovskite fill all five,
//      CsCl the corner and the centre.
// For each member: the world sublattice of each of its own sublattices (binaries:
// its sites in canonical order; templates: the parent's sublattices), or how its
// cell is placed in the world ('cscl-fold', 'chalcopyrite').
const LATTICE_WORLDS = {
    fcc: {
        sites: HEUSLER_SITES,
        members: {
            zincblende: [0, 2], rocksalt: [0, 1], fluorite: [0, 2, 3], antifluorite: [0, 2, 3],
            heusler: [0, 1, 2, 3], cscl: 'cscl-fold', chalcopyrite: 'chalcopyrite',
        },
        // the sites in the order of a shift by (1/4,1/4,1/4) in the fcc basis: shifting
        // material 2 by a multiple of it maps the site set onto itself
        cycle: [0, 2, 1, 3],
    },
    cubic: {
        sites: PEROVSKITE_SITES,
        members: { perovskite: [0, 1, 2, 3, 4], antiperovskite: [0, 1, 2, 3, 4], cscl: [0, 1] },
        cycle: null,
    },
};
// the doubled CsCl cell of the fcc world: rows of rocksaltLattice(2 b) in the basis b I
const CSCL_TO_FCC = [[0, 1, 1], [1, 0, 1], [1, 1, 0]];

function worldMemberKey(family) {
    let template = SUPERSTRUCTURES[family];
    return template ? template.parent : (family || 'zincblende');
}

export function getLatticeWorlds(family) {
    // the lattice worlds a family belongs to (LATTICE_WORLDS)
    let key = worldMemberKey(family);
    return Object.keys(LATTICE_WORLDS).filter((world) => LATTICE_WORLDS[world].members[key]);
}

function chalcopyriteWorldSupercell() {
    // the body-centred tetragonal chalcopyrite cell (c = 2 a) in the fcc basis
    return mat.matrix_multiply(chalcopyriteLattice(1, 2), mat.matrix_inverse(rocksaltLattice(1))).map((row) => row.map(Math.round));
}

const IDENTITY3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

function cubicOperations() {
    // the 48 orthogonal maps of the cube onto itself (signed permutations), identity first
    let permutations = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let out = [];
    for (let permutation of permutations) {
        for (let signs = 0; signs < 8; signs++) {
            out.push([0, 1, 2].map((i) => [0, 1, 2].map((j) => (j === permutation[i] ? ((signs >> i) & 1 ? -1 : 1) : 0))));
        }
    }
    return out;
}

let worldOperationCache = {};
function worldOperations(world) {
    /*
    the point operations of the world's parent: O (cartesian) and P, the same map
    in the parent basis (H O^T = P H, rows of H the parent lattice vectors)
    */
    if (!worldOperationCache[world]) {
        let H = world === 'fcc' ? rocksaltLattice(1) : IDENTITY3;
        let Hinverse = mat.matrix_inverse(H);
        worldOperationCache[world] = cubicOperations().map((O) => ({
            O: O,
            P: mat.matrix_multiply(mat.matrix_multiply(H, mat.matrix_transpose(O)), Hinverse).map((row) => row.map(Math.round)),
        }));
    }
    return worldOperationCache[world];
}

function transformWorldSite(site, P, world) {
    // a world site (sublattice, n) under the operation P of the parent basis
    let sites = LATTICE_WORLDS[world].sites;
    let v = rowTimesMatrix(sites[site.sub].red.map((x, k) => x + site.n[k]), P);
    let sub = sites.findIndex((candidate) => v.every((x, k) => Math.abs(x - candidate.red[k] - Math.round(x - candidate.red[k])) < 1e-6));
    return { sub: sub, n: v.map((x, k) => Math.round(x - sites[sub].red[k])) };
}

function csclFccSites() {
    // the world sites of CsCl folded into the doubled cell (csclInFccWorld), in its order
    let inverse = mat.matrix_inverse(CSCL_TO_FCC);
    let sites = [];
    for (let m of cellTranslations(CSCL_TO_FCC)) {
        CSCL_SITES.forEach((site) => {
            let p = rowTimesMatrix(site.red.map((v, k) => v + m[k]), inverse);
            let sub = HEUSLER_SITES.findIndex((h) => p.every((v, k) => Math.abs(v - h.red[k] - Math.round(v - h.red[k])) < 1e-6));
            sites.push({ sub: sub, n: p.map((v, k) => Math.round(v - HEUSLER_SITES[sub].red[k])) });
        });
    }
    return sites;
}

function familyWorldCell(family, world) {
    /*
    the cell of a family in the world's parent basis and the world sites of its
    canonical sites, as worldCell gives them for an end-member; CsCl in the fcc world
    is already folded into the doubled cell (prefolded: needs force constants)
    */
    let member = LATTICE_WORLDS[world].members[worldMemberKey(family)];
    if (member === 'cscl-fold') {
        return { supercell: IDENTITY3, sites: csclFccSites(), prefolded: true };
    }
    if (member === 'chalcopyrite') {
        let S = chalcopyriteWorldSupercell();
        return {
            supercell: S,
            sites: chalcopyritePositions(1, 2, 0.25).map((p, s) => {
                let sub = CHALCOPYRITE_SITES[s].role === 'X' ? 2 : 0;
                let parent = rowTimesMatrix(p, S);
                return { sub: sub, n: parent.map((v, k) => Math.round(v - HEUSLER_SITES[sub].red[k])) };
            }),
            prefolded: false,
        };
    }
    let template = SUPERSTRUCTURES[family];
    if (template) {
        return { supercell: template.supercell, sites: template.sites.map((site) => ({ sub: member[site.sub], n: site.n.slice() })), prefolded: false };
    }
    return { supercell: IDENTITY3, sites: member.map((sub) => ({ sub: sub, n: [0, 0, 0] })), prefolded: false };
}

function sameIntegerMatrix(p, q) {
    return p.every((row, i) => row.every((v, j) => Math.abs(v - q[i][j]) < 1e-9));
}

function supercellRatio(big, small) {
    // R with big = R small, if integer (big is a supercell of small), else null
    let R = mat.matrix_multiply(big, mat.matrix_inverse(small));
    return isIntegerMatrix(R) ? R.map((row) => row.map(Math.round)) : null;
}

function determinant(S) {
    return Math.abs(Math.round(mat.matrix_determinant(S)));
}

function commonSupercell(S1, S2) {
    /*
    the smallest cell that is a supercell of both (the intersection of the two
    lattices): one of the two if it contains the other (S1 when they are the same
    lattice), else the Hermite normal form of the smallest such lattice
    */
    if (supercellRatio(S2, S1)) {
        return supercellRatio(S1, S2) ? S1 : S2;
    }
    if (supercellRatio(S1, S2)) {
        return S1;
    }
    let d1 = determinant(S1);
    let d2 = determinant(S2);
    let gcd = (a, b) => (b ? gcd(b, a % b) : a);
    let lcm = d1 / gcd(d1, d2) * d2;
    for (let N = lcm; N <= d1 * d2; N += lcm) {
        for (let h11 = 1; h11 <= N; h11++) {
            if (N % h11) {
                continue;
            }
            for (let h22 = 1; h22 <= N / h11; h22++) {
                if ((N / h11) % h22) {
                    continue;
                }
                let h33 = N / h11 / h22;
                for (let h12 = 0; h12 < h22; h12++) {
                    for (let h13 = 0; h13 < h33; h13++) {
                        for (let h23 = 0; h23 < h33; h23++) {
                            let S = [[h11, h12, h13], [0, h22, h23], [0, 0, h33]];
                            if (supercellRatio(S, S1) && supercellRatio(S, S2)) {
                                return S;
                            }
                        }
                    }
                }
            }
        }
    }
    return null;
}

function foldWorldSites(sites, S, target) {
    // the world sites of a cell S repeated over the supercell target (as foldIntoWorld)
    let R = supercellRatio(target, S);
    let out = [];
    for (let m of cellTranslations(R)) {
        let shift = rowTimesMatrix(m, S);
        sites.forEach((site) => out.push({ sub: site.sub, n: site.n.map((v, k) => v + shift[k]) }));
    }
    return out;
}

function worldSiteKey(site, S) {
    return site.sub + ':' + reduceIntoCell(site.n, S).join(',');
}

function shiftWorldSite(site, j, world) {
    // a site moved by j quarters of the body diagonal (the fcc world's cycle)
    let cycle = LATTICE_WORLDS[world].cycle;
    if (!j) {
        return site;
    }
    let p = cycle.indexOf(site.sub) + j;
    let carry = Math.floor(p / cycle.length);
    return { sub: cycle[p % cycle.length], n: site.n.map((v) => v + carry) };
}

let worldPlanCache = new Map();
export function getWorldMixing(family1, family2) {
    /*
    how two families of different structures mix in a common lattice world, or null.
    Material 2 may be turned by any of the 48 operations of the cube (the same tilted
    cell along other axes) and both are placed in the smallest common supercell S;
    a material whose cell is smaller than S is folded into it, which evaluates D off
    the band path and needs its force constants (a cell of the same size in another
    basis is only relabelled). Returns the world and the options, best first: the
    operation, S, whether each material is folded, and how many sites coincide
    (after the best shift). The options run from the smallest common cell and the
    fewest folds to the most shared sites.
    */
    let cacheKey = family1 + '|' + family2;
    if (worldPlanCache.has(cacheKey)) {
        return worldPlanCache.get(cacheKey);
    }
    let worlds = getLatticeWorlds(family1).filter((world) => getLatticeWorlds(family2).indexOf(world) !== -1);
    let plan = null;
    if (worlds.length) {
        // CsCl with a perovskite: the cubic world, no fold
        let world = worlds.indexOf('cubic') !== -1 ? 'cubic' : worlds[0];
        let cell1 = familyWorldCell(family1, world);
        let cell2 = familyWorldCell(family2, world);
        let shifts = LATTICE_WORLDS[world].cycle ? [0, 1, 2, 3] : [0];
        let options = [];
        if (!(cell1.prefolded && cell2.prefolded)) {
            worldOperations(world).forEach((operation, index) => {
                let S2 = mat.matrix_multiply(cell2.supercell, operation.P).map((row) => row.map(Math.round));
                let S = commonSupercell(cell1.supercell, S2);
                if (!S) {
                    return;
                }
                let det = determinant(S);
                let fold1 = cell1.prefolded || det > determinant(cell1.supercell);
                let fold2 = cell2.prefolded || det > determinant(S2);
                let keys1 = new Set(foldWorldSites(cell1.sites, cell1.supercell, S).map((site) => worldSiteKey(site, S)));
                let sites2 = foldWorldSites(cell2.sites.map((site) => transformWorldSite(site, operation.P, world)), S2, S);
                // the most sites in common over the translations prepareWorldPair tries
                // (quarters of the body diagonal and lattice vectors modulo S), so that
                // the options are the same whichever material comes first
                let shared = 0;
                for (let j of shifts) {
                    let moved = sites2.map((site) => shiftWorldSite(site, j, world));
                    for (let n of cellTranslations(S)) {
                        let count = moved.filter((site) => keys1.has(worldSiteKey({ sub: site.sub, n: site.n.map((v, k) => v + n[k]) }, S))).length;
                        shared = Math.max(shared, count);
                    }
                }
                options.push({ index: index, O: operation.O, P: operation.P, supercell: S, det: det, fold1: fold1, fold2: fold2, shared: shared });
            });
        }
        options.sort((p, q) => p.det - q.det || (p.fold1 + p.fold2) - (q.fold1 + q.fold2) ||
            q.shared - p.shared || p.index - q.index);
        if (options.length) {
            plan = { world: world, options: options };
        }
    }
    worldPlanCache.set(cacheKey, plan);
    return plan;
}

function worldCell(endpoint, world) {
    /*
    an end-member in a lattice world: its supercell S of the world parent (canonical
    lattice = S H), canonical lattice and reduced positions, and the world site
    (sublattice, lattice vector n in the parent basis) of each of its sites
    */
    let member = LATTICE_WORLDS[world].members[worldMemberKey(endpoint.family)];
    let template = SUPERSTRUCTURES[endpoint.family];
    let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    if (member === 'chalcopyrite') {
        // the sites from the ideal cell (c = 2a, u = 1/4): cations on the fcc site,
        // anions on the tetrahedral site 1/4
        let S = chalcopyriteWorldSupercell();
        let ideal = chalcopyritePositions(1, 2, 0.25);
        return {
            supercell: S,
            cell: endpointLattice(endpoint),
            positions: endpointPositions(endpoint),
            sites: ideal.map((p, s) => {
                let sub = CHALCOPYRITE_SITES[s].role === 'X' ? 2 : 0;
                let parent = rowTimesMatrix(p, S);
                return { sub: sub, n: parent.map((v, k) => Math.round(v - HEUSLER_SITES[sub].red[k])) };
            }),
        };
    }
    if (template) {
        return {
            supercell: template.supercell,
            cell: endpoint.structure.lattice,
            positions: endpoint.structure.positionsRed,
            sites: template.sites.map((site) => ({ sub: member[site.sub], n: site.n.slice() })),
        };
    }
    // a binary in the standard fcc primitive cell (the files all are)
    let a = Math.SQRT2 * Math.sqrt(mat.vec_dot(endpoint.lattice[0], endpoint.lattice[0]));
    let fcc = rocksaltLattice(a);
    if (endpoint.rotation || !endpoint.lattice.every((row, i) => row.every((v, j) => Math.abs(v - fcc[i][j]) < 1e-6 * a))) {
        throw new Error('Unsupported ' + endpoint.family + ' setting: not the standard fcc primitive cell');
    }
    return {
        supercell: identity,
        cell: endpoint.lattice,
        positions: endpointPositions(endpoint),
        sites: member.map((sub) => ({ sub: sub, n: [0, 0, 0] })),
    };
}

function cellTranslations(R) {
    // the lattice vectors m (integer, basis of the small cell) distinct modulo the supercell R
    let count = Math.abs(Math.round(mat.matrix_determinant(R)));
    let inverse = mat.matrix_inverse(R);
    let reach = Math.max.apply(null, R.map((row) => row.reduce((s, v) => s + Math.abs(v), 0)));
    let found = new Map();
    for (let m0 = -reach; m0 <= reach; m0++) {
        for (let m1 = -reach; m1 <= reach; m1++) {
            for (let m2 = -reach; m2 <= reach; m2++) {
                let m = [m0, m1, m2];
                let f = rowTimesMatrix(m, inverse).map((v) => Math.floor(v + 1e-9));
                let back = rowTimesMatrix(f, R);
                let reduced = m.map((v, k) => Math.round(v - back[k]));
                found.set(reduced.join(','), reduced);
            }
        }
    }
    let list = Array.from(found.values());
    if (list.length !== count) {
        throw new Error('Expected ' + count + ' translations of the supercell, found ' + list.length);
    }
    list.sort((p, q) => Math.abs(p[0]) + Math.abs(p[1]) + Math.abs(p[2]) - Math.abs(q[0]) - Math.abs(q[1]) - Math.abs(q[2]));
    return list;
}

function foldIntoWorld(endpoint, info, R, worldSitesOf = null) {
    /*
    an end-member folded into the supercell R of its own cell (canonical lattice
    R C), with the world sites of the folded cell: each of its sites once per lattice
    vector m of its cell inside the supercell (n + m S in the parent basis, or
    worldSitesOf for CsCl, whose cell is not one of the world's)
    */
    let evaluate = endpoint.sampleAtQ || endpoint.evaluate;
    if (endpoint.family === 'chalcopyrite' && !chalcopyriteInCanonicalBasis(endpoint)) {
        throw new Error('Unsupported chalcopyrite setting: not the standard body-centred tetragonal cell');
    }
    if (!evaluate) {
        throw new Error('Mixing ' + endpoint.name + ' with a material of another structure needs its force constants ' +
            '(PhononDB); the Materials Project files only cover the band path');
    }
    let primitive = endpoint.sampleAtQ
        ? Object.assign({}, endpoint, { lattice: info.cell, evaluate: endpoint.sampleAtQ })
        : endpoint;
    let superLattice = mat.matrix_multiply(R, info.cell);
    let inverse = mat.matrix_inverse(R);
    let positions = [];
    let sublattice = [];
    let sites = [];
    for (let m of cellTranslations(R)) {
        info.positions.forEach((p, s) => {
            positions.push(rowTimesMatrix(p.map((v, k) => v + m[k]), inverse));
            sublattice.push(s);
            if (!worldSitesOf) {
                let shift = rowTimesMatrix(m, info.supercell);
                sites.push({ sub: info.sites[s].sub, n: info.sites[s].n.map((v, k) => v + shift[k]) });
            }
        });
    }
    let folded = foldIntoSupercell(primitive, endpoint.family, superLattice, positions, sublattice,
        { lattice: superLattice, positionsRed: positions });
    folded.foldedFrom = endpoint;
    // its band path in the supercell's reduced q (q R^T), in case the other material
    // is folded too and has no path of its own
    let toSuper = (q) => rowTimesMatrix(q, mat.matrix_transpose(R));
    folded.pathSegments = (endpoint.pathSegments || endpoint.segments).map((segment) => ({
        qA: toSuper(segment.qA),
        qB: toSuper(segment.qB),
        labelA: segment.labelA,
        labelB: segment.labelB,
        npoints: segment.npoints,
    }));
    return {
        endpoint: folded,
        info: {
            supercell: mat.matrix_multiply(R, info.supercell).map((row) => row.map(Math.round)),
            cell: superLattice,
            positions: positions,
            sites: worldSitesOf ? worldSitesOf(positions) : sites,
        },
    };
}

function csclInFccWorld(endpoint) {
    // CsCl folded into the doubled cell, rocksaltLattice(2 b): its sites from the
    // ideal positions (corner -> fcc and octahedral sites, centre -> tetrahedral ones)
    let info = {
        supercell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        cell: endpoint.structure.lattice,
        positions: endpoint.structure.positionsRed,
        sites: null,
    };
    let inverse = mat.matrix_inverse(CSCL_TO_FCC);
    let ideal = [];
    for (let m of cellTranslations(CSCL_TO_FCC)) {
        CSCL_SITES.forEach((site) => ideal.push(rowTimesMatrix(site.red.map((v, k) => v + m[k]), inverse)));
    }
    let placed = foldIntoWorld(endpoint, info, CSCL_TO_FCC, () => ideal.map((p) => {
        let sub = HEUSLER_SITES.findIndex((site) => p.every((v, k) => Math.abs(v - site.red[k] - Math.round(v - site.red[k])) < 1e-6));
        return { sub: sub, n: p.map((v, k) => Math.round(v - HEUSLER_SITES[sub].red[k])) };
    }));
    // the doubled cell is the fcc world's primitive cell
    placed.info.supercell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    return placed;
}

function reduceIntoCell(n, S) {
    let f = rowTimesMatrix(n, mat.matrix_inverse(S)).map((v) => Math.floor(v + 1e-6));
    let back = rowTimesMatrix(f, S);
    return n.map((v, k) => Math.round(v - back[k]));
}

function hermiteRows(S) {
    // the Hermite normal form of an integer lattice basis (rows): upper triangular,
    // positive diagonal, entries above the diagonal reduced; the same for any basis
    let A = S.map((row) => row.map(Math.round));
    let r = 0;
    for (let c = 0; c < 3 && r < 3; c++) {
        for (;;) {
            let pivot = -1;
            for (let i = r; i < 3; i++) {
                if (A[i][c] !== 0 && (pivot < 0 || Math.abs(A[i][c]) < Math.abs(A[pivot][c]))) {
                    pivot = i;
                }
            }
            if (pivot < 0) {
                break;
            }
            [A[r], A[pivot]] = [A[pivot], A[r]];
            let done = true;
            for (let i = r + 1; i < 3; i++) {
                let q = Math.floor(A[i][c] / A[r][c]);
                A[i] = A[i].map((v, k) => v - q * A[r][k]);
                if (A[i][c] !== 0) {
                    done = false;
                }
            }
            if (done) {
                break;
            }
        }
        if (A[r][c] < 0) {
            A[r] = A[r].map((v) => -v);
        }
        for (let i = 0; i < r; i++) {
            let q = Math.floor(A[i][c] / A[r][c]);
            A[i] = A[i].map((v, k) => v - q * A[r][k]);
        }
        r += 1;
    }
    return A;
}

function isAnionSite(type) {
    return bccAnionRank(type) > 0;
}

function reframeEndpoint(endpoint, info, O, M) {
    /*
    an end-member turned by the cube operation O (cartesian) and expressed in the
    basis M C O^T of its cell (M integer, det +-1): blocks of D and the Born and
    dielectric tensors rotated, O D O^T; reduced positions p M^-1 and reduced q
    q M^T. The reduced coordinates do not change under the rotation itself.
    */
    if (sameIntegerMatrix(O, IDENTITY3) && sameIntegerMatrix(M, IDENTITY3)) {
        return { endpoint: endpoint, info: info };
    }
    let rotateBlocks = (matrix) => {
        let n = matrix.real.length / 3;
        let out = zeroMatrix(3 * n);
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                for (let part of ['real', 'imag']) {
                    let block = [0, 1, 2].map((a) => [0, 1, 2].map((b) => matrix[part][3 * i + a][3 * j + b]));
                    let rotated = rotateTensor(O, block);
                    for (let a = 0; a < 3; a++) {
                        for (let b = 0; b < 3; b++) {
                            out[part][3 * i + a][3 * j + b] = rotated[a][b];
                        }
                    }
                }
            }
        }
        return out;
    };
    let canonical = (matrix) => rotateBlocks(toCanonical(endpoint, matrix));
    let back = mat.matrix_transpose(mat.matrix_inverse(M));
    let toOld = (q) => q && rowTimesMatrix(q, back);
    let toNew = (q) => rowTimesMatrix(q, mat.matrix_transpose(M));
    let evaluate = endpoint.sampleAtQ || endpoint.evaluate;
    let evaluateShortRange = endpoint.sampleShortRangeAtQ || endpoint.evaluateShortRange;
    let sampleAtQ = evaluate ? (q, direction) => canonical(evaluate(toOld(q), toOld(direction))) : null;
    let sampleShortRangeAtQ = evaluateShortRange ? (q, direction) => canonical(evaluateShortRange(toOld(q), toOld(direction))) : null;
    let lattice = mat.matrix_multiply(M, info.cell.map((row) => matrixTimesVector(O, row)));
    let Minverse = mat.matrix_inverse(M);
    let positions = info.positions.map((p) => rowTimesMatrix(p, Minverse));
    let nac = endpoint.nac ? Object.assign({}, endpoint.nac, {
        siteBorn: endpoint.nac.siteBorn.map((born) => rotateTensor(O, born)),
        siteDielectric: rotateTensor(O, endpoint.nac.siteDielectric),
    }) : null;
    let reframed = {
        family: endpoint.family,
        name: endpoint.name,
        anyQ: endpoint.anyQ,
        siteTypes: endpoint.siteTypes,
        siteMasses: endpoint.siteMasses,
        order: endpoint.siteTypes.map((_, i) => i),
        invert: false,
        rotation: null,
        lattice: lattice,
        positionsRed: positions,
        structure: { lattice: lattice, positionsRed: positions },
        nac: nac,
        frequencyFactor: endpoint.frequencyFactor,
        segments: endpoint.segments.map((segment) => Object.assign({}, segment, {
            qA: toNew(segment.qA),
            qB: toNew(segment.qB),
            sample: (t) => canonical(segment.sample(t)),
            sampleShortRange: segment.sampleShortRange ? (t) => canonical(segment.sampleShortRange(t)) : null,
        })),
        pathSegments: endpoint.pathSegments ? endpoint.pathSegments.map((segment) => Object.assign({}, segment, {
            qA: toNew(segment.qA),
            qB: toNew(segment.qB),
        })) : null,
        sampleAtQ: sampleAtQ,
        sampleShortRangeAtQ: sampleShortRangeAtQ,
        evaluate: sampleAtQ,
        evaluateShortRange: sampleShortRangeAtQ,
    };
    return {
        endpoint: reframed,
        info: {
            supercell: mat.matrix_multiply(M, info.supercell).map((row) => row.map(Math.round)),
            cell: lattice,
            positions: positions,
            sites: info.sites,
        },
    };
}

function chalcopyriteInCanonicalBasis(endpoint) {
    /*
    a chalcopyrite file evaluated at canonical reduced q: its cell is the canonical
    body-centred tetragonal one, turned (as zincblendeAsChalcopyrite assumes)
    */
    let T = endpoint.rotation || IDENTITY3;
    let rotated = endpoint.lattice.map((row) => matrixTimesVector(T, row));
    let U = mat.matrix_multiply(rotated, mat.matrix_inverse(endpointLattice(endpoint)));
    return U.every((row, i) => row.every((v, j) => Math.abs(v - (i === j ? 1 : 0)) < 1e-3));
}

function canFold(endpoint) {
    return !!(endpoint.sampleAtQ || endpoint.evaluate) &&
        (endpoint.family !== 'chalcopyrite' || chalcopyriteInCanonicalBasis(endpoint));
}

function prepareWorldPair(endpoint1, endpoint2) {
    /*
    two end-members of different structures in their common lattice world
    (getWorldMixing): material 2 turned by the chosen cube operation, both folded
    into the common supercell where their cells are smaller, and the union of their
    sites. Material 2 may also be shifted by a translation that maps the world's
    site set onto itself (a quarter of the body diagonal in the fcc world). Of the
    options that need only the force constants available and share the most sites,
    the one that puts anions on anions; on a tie the same relative arrangement
    whichever material comes first. Returns the end-members, and for every union site
    its index in each (-1: vacancy).
    */
    let plan = getWorldMixing(endpoint1.family, endpoint2.family);
    if (!plan) {
        return null;
    }
    let world = plan.world;
    let feasible = plan.options.filter((option) => (!option.fold1 || canFold(endpoint1)) && (!option.fold2 || canFold(endpoint2)));
    if (!feasible.length) {
        let missing = !canFold(endpoint1) ? endpoint1 : endpoint2;
        throw new Error('Mixing ' + missing.name + ' with ' + (missing === endpoint1 ? endpoint2 : endpoint1).name +
            ' folds it into a larger cell, which needs its force constants (PhononDB); the Materials Project files only cover the band path');
    }
    let top = feasible.filter((option) => option.det === feasible[0].det &&
        option.fold1 + option.fold2 === feasible[0].fold1 + feasible[0].fold2 && option.shared === feasible[0].shared);
    let place = (endpoint) => (LATTICE_WORLDS[world].members[worldMemberKey(endpoint.family)] === 'cscl-fold'
        ? csclInFccWorld(endpoint)
        : { endpoint: endpoint, info: worldCell(endpoint, world) });
    let into = (placed, S) => {
        // folded into S, or (the same lattice in another basis) relabelled
        let R = supercellRatio(S, placed.info.supercell);
        if (sameIntegerMatrix(R, IDENTITY3)) {
            return placed;
        }
        if (Math.abs(determinant(R)) === 1) {
            return reframeEndpoint(placed.endpoint, placed.info, IDENTITY3, R);
        }
        return foldIntoWorld(placed.endpoint, placed.info, R);
    };
    let base1 = place(endpoint1);
    let base2 = place(endpoint2);
    let shifts = LATTICE_WORLDS[world].cycle ? [0, 1, 2, 3] : [0];
    let order = (endpoint1.name || '') <= (endpoint2.name || '');
    let best = null;
    // ties: the same arrangement whichever material comes first. Material 2 turned by
    // O and translated by t is, seen from material 2, material 1 turned by O^-1 and
    // translated by -O^-1 t; both are described from the material first by name, with
    // t split into quarters of the body diagonal and a lattice vector of the parent,
    // reduced modulo the common cell (in its Hermite normal form, so that the same
    // coset is described the same way in any basis)
    let operations = worldOperations(world);
    let quarter = [0.25, 0.25, 0.25];
    let isInteger = (v) => v.every((x) => Math.abs(x - Math.round(x)) < 1e-6);
    let tieKey = (option, t, S) => {
        let index = option.index;
        let cell = S;
        if (!order) {
            index = operations.findIndex((operation) => sameIntegerMatrix(operation.O, mat.matrix_transpose(option.O)));
            let P = operations[index].P;
            t = rowTimesMatrix(t, P).map((v) => -v);
            cell = mat.matrix_multiply(S, P).map((row) => row.map(Math.round));
        }
        let j = shifts.find((jj) => isInteger(t.map((v, k) => v - jj * quarter[k])));
        let n = reduceIntoCell(t.map((v, k) => Math.round(v - j * quarter[k])), hermiteRows(cell));
        return [index, j].concat(n);
    };
    let before = (p, q) => {
        for (let i = 0; i < p.length; i++) {
            if (p[i] !== q[i]) {
                return p[i] < q[i];
            }
        }
        return false;
    };
    for (let option of top) {
        let S = option.supercell;
        let p1 = into(base1, S);
        let turned = reframeEndpoint(base2.endpoint, Object.assign({}, base2.info, {
            sites: base2.info.sites.map((site) => transformWorldSite(site, option.P, world)),
            supercell: mat.matrix_multiply(base2.info.supercell, option.P).map((row) => row.map(Math.round)),
        }), option.O, IDENTITY3);
        // reframeEndpoint keeps supercell M S; here the cell is S2 P after the rotation
        turned.info.supercell = mat.matrix_multiply(base2.info.supercell, option.P).map((row) => row.map(Math.round));
        let p2 = into(turned, S);
        let keys1 = p1.info.sites.map((site) => worldSiteKey(site, S));
        let index1 = new Map(keys1.map((k, i) => [k, i]));
        // translations of material 2: quarters of the body diagonal (the fcc world's
        // site cycle) and lattice vectors of the parent modulo the common cell (an
        // ordering is not invariant under them)
        let lattice = cellTranslations(S);
        for (let j of shifts) {
            for (let n of lattice) {
                let keys2 = p2.info.sites.map((site) => {
                    let moved = shiftWorldSite(site, j, world);
                    return worldSiteKey({ sub: moved.sub, n: moved.n.map((v, k) => v + n[k]) }, S);
                });
                let shared = 0;
                let agree = 0;
                keys2.forEach((k, i) => {
                    if (index1.has(k)) {
                        shared += 1;
                        if (isAnionSite(p1.endpoint.siteTypes[index1.get(k)]) === isAnionSite(p2.endpoint.siteTypes[i])) {
                            agree += 1;
                        }
                    }
                });
                let score = shared * 1000 + agree;
                let t = quarter.map((v, k) => v * j + n[k]);
                if (best && score < best.score) {
                    continue;
                }
                let tie = tieKey(option, t, S);
                if (!best || score > best.score || before(tie, best.tie)) {
                    best = { option: option, score: score, tie: tie, shift: j, translation: n, t: t,
                        keys1: keys1, keys2: keys2, index1: index1, p1: p1, p2: p2, S: S };
                }
            }
        }
    }
    let { p1, p2, S, keys1, index1 } = best;
    // the translation in reduced coordinates of the common cell
    let shiftRed = rowTimesMatrix(best.t, mat.matrix_inverse(S));
    let positions2 = p2.info.positions.map((p) => p.map((v, k) => v + shiftRed[k]));
    let map1 = keys1.map((_, i) => i);
    let map2 = keys1.map(() => -1);
    best.keys2.forEach((k, i) => {
        if (index1.has(k)) {
            let j = index1.get(k);
            map2[j] = i;
            // the same site may lie a lattice vector away in the other cell: put it on
            // material 1's copy, so that the mixture does not move it across the cell
            let target = p1.info.positions[j];
            positions2[i] = positions2[i].map((v, d) => v - Math.round(v - target[d]));
        } else {
            map1.push(-1);
            map2.push(i);
        }
    });
    return {
        endpoint1: p1.endpoint,
        endpoint2: p2.endpoint,
        union: {
            world: world,
            map1: map1,
            map2: map2,
            cell1: p1.info.cell,
            cell2: p2.info.cell,
            positions1: p1.info.positions,
            positions2: positions2,
            // the positions the sampled D of each end-member refers to (its gauge)
            ownPositions1: p1.info.positions,
            ownPositions2: p2.info.positions,
            shift2: best.shift,
            // the chosen arrangement (operation index, shift) and its tie key, for checks
            arrangement: { operation: best.option.index, shift: best.shift, translation: best.translation, tie: best.tie, score: best.score },
        },
    };
}

function lineParameter(segment, q) {
    /*
    t with q = qA + t (qB - qA), or null when q is not on the line of the segment
    */
    let d = [0, 1, 2].map((i) => segment.qB[i] - segment.qA[i]);
    let t = [0, 1, 2].reduce((sum, i) => sum + (q[i] - segment.qA[i]) * d[i], 0) / mat.vec_dot(d, d);
    let residual = Math.max.apply(null, [0, 1, 2].map((i) => Math.abs(segment.qA[i] + t * d[i] - q[i])));
    return residual < 1e-4 ? t : null;
}

function findOverlap(segment, others) {
    /*
    the other segment lying on the same line with the longest common stretch.
    lo/hi delimit that stretch on `segment`, mapOther turns a parameter of
    `segment` into one of the other segment
    */
    let best = null;
    // time reversal: D(-q) = D(q)^* (real force constants, either gauge), so the
    // image -q of another segment, sampled with the complex conjugate, covers the
    // same physics; two files may run the same line in opposite halves (Gamma -> Z
    // towards (1/2,1/2,1/2) in one, (-1/2,-1/2,-1/2) in the other, e.g. ilmenites,
    // whose R-3 symmetry has inversion but not the two-fold axes of corundum)
    let conjugate = (sampler) => sampler && ((t) => {
        let matrix = sampler(t);
        return { real: matrix.real, imag: matrix.imag.map((row) => row.map((v) => -v)) };
    });
    let reversed = others.map((other) => Object.assign({}, other, {
        qA: other.qA.map((v) => -v),
        qB: other.qB.map((v) => -v),
        sample: conjugate(other.sample),
        sampleShortRange: conjugate(other.sampleShortRange),
    }));
    // the segments themselves first: on a tie the direct one is kept
    for (let other of others.concat(reversed)) {
        let tA = lineParameter(segment, other.qA);
        let tB = lineParameter(segment, other.qB);
        if (tA === null || tB === null) {
            continue;
        }
        let lo = Math.max(0, Math.min(tA, tB));
        let hi = Math.min(1, Math.max(tA, tB));
        if (hi - lo < 0.05 || (best && hi - lo <= best.hi - best.lo)) {
            continue;
        }
        let near = (a, b) => Math.abs(a - b) < 1e-6;
        best = {
            other: other,
            lo: near(lo, 0) ? 0 : lo,
            hi: near(hi, 1) ? 1 : hi,
            mapOther: (t) => (near(tA, 0) && near(tB, 1) ? t : near(tA, 1) && near(tB, 0) ? 1 - t : (t - tA) / (tB - tA)),
            // label of the other segment at a parameter of `segment` where it ends
            otherLabel: (t) => (near(t, tA) ? other.labelA : near(t, tB) ? other.labelB : null),
        };
    }
    return best;
}

function buildTemplatePath(endpoint1, endpoint2, template) {
    /*
    the standard path of a ternary wurtzite cell, both end-members evaluated at
    the same canonical q; points per segment from its length in the average of
    the two cells, so the path is the same whichever material comes first
    */
    let lattice = endpoint1.structure.lattice.map((row, i) => row.map((v, j) => (v + endpoint2.structure.lattice[i][j]) / 2));
    let reciprocal = utils.rec_lat(lattice);
    let point = (label) => template.specialPoints.find((p) => p.label === label).q;
    return template.path.map(([labelA, labelB]) => {
        let qA = point(labelA);
        let qB = point(labelB);
        let direction = [0, 1, 2].map((d) => qB[d] - qA[d]);
        let length = Math.sqrt(mat.vec_dot(rowTimesMatrix(direction, reciprocal), rowTimesMatrix(direction, reciprocal)));
        let atQ = (sampler) => sampler && ((t) => sampler(lerpQ(qA, qB, t), direction));
        return {
            qA: qA,
            qB: qB,
            labelA: labelA,
            labelB: labelB,
            npoints: Math.max(5, Math.round(length * 250) + 1),
            sample1: atQ(endpoint1.sampleAtQ),
            sample2: atQ(endpoint2.sampleAtQ),
            sampleShortRange1: atQ(endpoint1.sampleShortRangeAtQ),
            sampleShortRange2: atQ(endpoint2.sampleShortRangeAtQ),
        };
    });
}

export function buildCommonPath(endpoint1, endpoint2) {
    /*
    the segments of the end-member with the fewer segments (on a tie the PhononDB
    one, then endpoint1),
    each kept where the other end-member has a segment on the same line, in
    either direction, and cut to the stretch both cover. For zinc blende the
    segments coincide; for chalcopyrite the zone-boundary points Sigma, Y, ...
    move with c/a, so a segment ending there is cut at the nearer of the two.
    */
    let template = SUPERSTRUCTURES[endpoint1.family];
    if (endpoint1.sampleAtQ && endpoint2.sampleAtQ && template && template.path) {
        // both known at any q: the standard path of the canonical cell
        return buildTemplatePath(endpoint1, endpoint2, template);
    }
    if (endpoint1.sampleAtQ || endpoint2.sampleAtQ) {
        // one end-member can be evaluated at any q (a binary folded into the cell of
        // the ternary, or force constants): the path is the whole path of the other
        let real;
        if (!endpoint1.sampleAtQ) {
            real = endpoint1;
        } else if (!endpoint2.sampleAtQ) {
            real = endpoint2;
        } else {
            // both known at any q, no standard path: the path of the file of one of
            // them, the same one whichever material comes first
            let hasPath = [endpoint1, endpoint2].filter((endpoint) => endpoint.segments.length);
            if (!hasPath.length) {
                // both folded into a common cell (two orderings): the path of one of
                // them, the same one whichever material comes first, evaluated for both
                let source = [endpoint1, endpoint2].filter((endpoint) => endpoint.pathSegments && endpoint.pathSegments.length)
                    .sort((p, q) => (p.name < q.name ? -1 : p.name > q.name ? 1 : 0))[0];
                if (source) {
                    return source.pathSegments.map((segment) => {
                        let direction = [0, 1, 2].map((d) => segment.qB[d] - segment.qA[d]);
                        let atQ = (sampler) => sampler && ((t) => sampler(lerpQ(segment.qA, segment.qB, t), direction));
                        return {
                            qA: segment.qA,
                            qB: segment.qB,
                            labelA: segment.labelA,
                            labelB: segment.labelB,
                            npoints: segment.npoints,
                            sample1: atQ(endpoint1.sampleAtQ),
                            sample2: atQ(endpoint2.sampleAtQ),
                            sampleShortRange1: atQ(endpoint1.sampleShortRangeAtQ),
                            sampleShortRange2: atQ(endpoint2.sampleShortRangeAtQ),
                        };
                    });
                }
            }
            real = hasPath.length === 2
                ? (endpoint2.name < endpoint1.name ? endpoint2 : endpoint1)
                : hasPath[0] || endpoint1;
        }
        let folded = real === endpoint1 ? endpoint2 : endpoint1;
        return real.segments.map((segment) => {
            let direction = [0, 1, 2].map((d) => segment.qB[d] - segment.qA[d]);
            let atQ = (sampler) => sampler && ((t) => sampler(lerpQ(segment.qA, segment.qB, t), direction));
            let sampleFolded = atQ(folded.sampleAtQ);
            let shortRangeFolded = atQ(folded.sampleShortRangeAtQ);
            return {
                qA: segment.qA,
                qB: segment.qB,
                labelA: segment.labelA,
                labelB: segment.labelB,
                npoints: segment.npoints,
                sample1: endpoint1 === real ? segment.sample : sampleFolded,
                sample2: endpoint2 === real ? segment.sample : sampleFolded,
                sampleShortRange1: endpoint1 === real ? segment.sampleShortRange : shortRangeFolded,
                sampleShortRange2: endpoint2 === real ? segment.sampleShortRange : shortRangeFolded,
            };
        });
    }

    // on a tie the path must not depend on which material is material 1 (segments
    // can run in opposite directions, e.g. H-K and K-H): prefer the PhononDB one
    let lesser;
    if (endpoint1.segments.length !== endpoint2.segments.length) {
        lesser = endpoint2.segments.length < endpoint1.segments.length ? endpoint2 : endpoint1;
    } else {
        lesser = !endpoint1.evaluate && endpoint2.evaluate ? endpoint2 : endpoint1;
    }
    let other = lesser === endpoint1 ? endpoint2 : endpoint1;

    let path = [];
    for (let i = 0; i < lesser.segments.length; i++) {
        let segment = lesser.segments[i];
        let match = findOverlap(segment, other.segments);
        if (!match) {
            continue;
        }
        let lo = match.lo;
        let hi = match.hi;
        let oLo = match.mapOther(lo);
        let oHi = match.mapOther(hi);
        let onLesser = (sampler) => sampler && ((t) => sampler(lo + t * (hi - lo)));
        let onOther = (sampler) => sampler && ((t) => sampler(oLo + t * (oHi - oLo)));
        let sampleLesser = onLesser(segment.sample);
        let sampleOther = onOther(match.other.sample);
        let shortRangeLesser = onLesser(segment.sampleShortRange);
        let shortRangeOther = onOther(match.other.sampleShortRange);
        path.push({
            qA: lerpQ(segment.qA, segment.qB, lo),
            qB: lerpQ(segment.qA, segment.qB, hi),
            labelA: lo === 0 ? segment.labelA : (match.otherLabel(lo) || ''),
            labelB: hi === 1 ? segment.labelB : (match.otherLabel(hi) || ''),
            npoints: Math.max(
                Math.round((hi - lo) * (segment.npoints - 1)) + 1,
                Math.round(Math.abs(oHi - oLo) * (match.other.npoints - 1)) + 1,
                2
            ),
            sample1: lesser === endpoint1 ? sampleLesser : sampleOther,
            sample2: lesser === endpoint1 ? sampleOther : sampleLesser,
            sampleShortRange1: lesser === endpoint1 ? shortRangeLesser : shortRangeOther,
            sampleShortRange2: lesser === endpoint1 ? shortRangeOther : shortRangeLesser,
        });
    }
    return path;
}

export function getParentFamily(family) {
    // the binary family an ordered ternary family is a supercell of
    if (family === 'chalcopyrite') {
        return 'zincblende';
    }
    // (the 5-atom perovskite is a template of its own parent, not a superstructure)
    let template = SUPERSTRUCTURES[family];
    return template && template.parent !== family ? template.parent : null;
}

export function sampleMixingPair(endpoint1, endpoint2) {
    /*
    canonical D_1 and D_2 at every q-point of the common path; independent of x.
    When both end-members carry Born charges and dielectric tensors (PhononDB),
    their short-range parts are kept too, so the dipole-dipole term can be
    recomputed for the mixed crystal instead of being interpolated.
    */
    let union = null;
    if (endpoint1.family !== endpoint2.family) {
        // a binary goes into the cell of an ordered ternary on the same lattice
        let fold = (binary, family) => (family === 'chalcopyrite'
            ? zincblendeAsChalcopyrite(binary)
            : binaryAsSuperstructure(binary, family));
        if (getParentFamily(endpoint2.family) === endpoint1.family) {
            endpoint1 = fold(endpoint1, endpoint2.family);
        } else if (getParentFamily(endpoint1.family) === endpoint2.family) {
            endpoint2 = fold(endpoint2, endpoint1.family);
        } else {
            // other structures on a common lattice: site by site, with vacancies
            let world = prepareWorldPair(endpoint1, endpoint2);
            if (!world) {
                throw new Error('Cannot mix a ' + endpoint1.family + ' cell with a ' + endpoint2.family + ' cell');
            }
            endpoint1 = world.endpoint1;
            endpoint2 = world.endpoint2;
            union = world.union;
        }
    }
    let recomputeNac = !!(endpoint1.nac && endpoint2.nac);
    let shortRange1 = [];
    let shortRange2 = [];
    let qDirections = [];
    let path = buildCommonPath(endpoint1, endpoint2);
    if (!path.length) {
        // two files known only on their band paths, of different types (a cubic and a
        // polar perovskite): nothing to sample both on
        throw new Error('The band paths of ' + endpoint1.name + ' and ' + endpoint2.name +
            ' have no segment in common; a PhononDB file for one of them can be evaluated on the other\'s path');
    }
    let qpoints = [];
    let matrices1 = [];
    let matrices2 = [];
    let lineBreaks = [];
    let labels = [];

    for (let s = 0; s < path.length; s++) {
        let segment = path[s];
        let start = qpoints.length;
        for (let k = 0; k < segment.npoints; k++) {
            let t = k / (segment.npoints - 1);
            let q = lerpQ(segment.qA, segment.qB, t);
            qpoints.push(q);
            matrices1.push(toCanonical(endpoint1, segment.sample1(t)));
            matrices2.push(toCanonical(endpoint2, segment.sample2(t)));
            if (recomputeNac) {
                shortRange1.push(toCanonical(endpoint1, segment.sampleShortRange1(t)));
                shortRange2.push(toCanonical(endpoint2, segment.sampleShortRange2(t)));
                qDirections.push(isGamma(q) ? [
                    segment.qB[0] - segment.qA[0],
                    segment.qB[1] - segment.qA[1],
                    segment.qB[2] - segment.qA[2],
                ] : null);
            }
        }
        lineBreaks.push([start, qpoints.length]);
        labels.push([segment.labelA, segment.labelB]);
    }

    if (recomputeNac) {
        shortRange1 = subtractRealSpaceDipole(endpoint1, shortRange1, qpoints);
        shortRange2 = subtractRealSpaceDipole(endpoint2, shortRange2, qpoints);
    }

    return {
        endpoint1: endpoint1,
        endpoint2: endpoint2,
        qpoints: qpoints,
        matrices1: matrices1,
        matrices2: matrices2,
        recomputeNac: recomputeNac,
        shortRange1: shortRange1,
        shortRange2: shortRange2,
        qDirections: qDirections,
        lineBreaks: lineBreaks,
        labels: labels,
        // sites of two different structures (prepareWorldPair), null within a family
        union: union,
    };
}

function eigenvalueToCm1(value) {
    let magnitude = Math.sqrt(Math.abs(value)) * THZ_TO_CM1;
    return value < 0 ? -magnitude : magnitude;
}

async function diagonalize(matrix) {
    let solution = await solveComplexHermitianWithEigenWasm(matrix.real, matrix.imag);
    return {
        frequencies: solution.values.map(eigenvalueToCm1),
        vectors: solution.vectors,
    };
}

export async function computeReferenceEigenvalues(sampled) {
    let eigenvalues = [];
    for (let k = 0; k < sampled.qpoints.length; k++) {
        eigenvalues.push((await diagonalize(sampled.matrices1[k])).frequencies);
    }
    return eigenvalues;
}

function getPathDistances(qpoints, lineBreaks, lattice) {
    let rec = utils.rec_lat(lattice);
    let distances = new Array(qpoints.length).fill(0);
    let dist = 0;
    for (let s = 0; s < lineBreaks.length; s++) {
        let start = lineBreaks[s][0];
        let end = lineBreaks[s][1];
        distances[start] = dist;
        for (let k = start + 1; k < end; k++) {
            dist += mat.distance(utils.red_car(qpoints[k - 1], rec), utils.red_car(qpoints[k], rec));
            distances[k] = dist;
        }
    }
    return distances;
}

function getHighSymmetryPoints(lineBreaks, labels) {
    /*
    [index, label] list; a segment end and the next segment start share the same
    distance, so different labels there are merged as "B|A"
    */
    let points = [];
    for (let s = 0; s < lineBreaks.length; s++) {
        let start = lineBreaks[s][0];
        let end = lineBreaks[s][1];
        let labelA = labels[s][0];
        if (points.length && points[points.length - 1][0] === start - 1) {
            let previous = points.pop();
            labelA = previous[1] === labelA ? labelA : previous[1] + '|' + labelA;
        }
        points.push([start, labelA]);
        points.push([end - 1, labels[s][1]]);
    }
    return points;
}

function addBlockTensor(matrix, blocks, factor) {
    /*
    add a [i][j][alpha][beta][re,im] block tensor, scaled by factor, to a matrix
    */
    for (let i = 0; i < blocks.length; i++) {
        for (let j = 0; j < blocks[i].length; j++) {
            for (let alpha = 0; alpha < 3; alpha++) {
                for (let beta = 0; beta < 3; beta++) {
                    matrix.real[i * 3 + alpha][j * 3 + beta] += factor * blocks[i][j][alpha][beta][0];
                    matrix.imag[i * 3 + alpha][j * 3 + beta] += factor * blocks[i][j][alpha][beta][1];
                }
            }
        }
    }
}

function getGonzeGList(lattice, reciprocal, gCutoff) {
    /*
    reciprocal lattice points (cartesian, no 2 pi) inside the G cutoff, as in
    phonopy DynamicalMatrixGL._get_G_list
    */
    // G = sum_i n_i b_i has n_i = G.a_i, so |n_i| <= gCutoff |a_i|
    let gRadius = Math.ceil(gCutoff * Math.max.apply(null, lattice.map((a) => Math.sqrt(mat.vec_dot(a, a)))));
    let gList = [];
    for (let a = -gRadius; a <= gRadius; a++) {
        for (let b = -gRadius; b <= gRadius; b++) {
            for (let c = -gRadius; c <= gRadius; c++) {
                let g = [0, 1, 2].map((d) => a * reciprocal[0][d] + b * reciprocal[1][d] + c * reciprocal[2][d]);
                if (mat.vec_dot(g, g) < gCutoff * gCutoff) {
                    gList.push(g);
                }
            }
        }
    }
    return gList;
}

function getGonzeDdQ0(gList, born, dielectric, positionsCar, lambda) {
    /*
    sum_j sum_{G != 0} Z_i^T K(G) Z_j exp(2 pi i G.(r_i - r_j)): the q = 0 term that
    phonopy subtracts from the diagonal blocks to keep the acoustic sum rule.
    phonopy stores the hermitian part of each 3x3 block; it only differs from the
    plain sum when the Born tensors are not symmetric (e.g. chalcopyrite)
    */
    let natoms = born.length;
    let ddQ0 = { real: [], imag: [] };
    for (let i = 0; i < natoms; i++) {
        ddQ0.real.push([[0, 0, 0], [0, 0, 0], [0, 0, 0]]);
        ddQ0.imag.push([[0, 0, 0], [0, 0, 0], [0, 0, 0]]);
    }
    let l2 = 4 * lambda * lambda;
    for (let g = 0; g < gList.length; g++) {
        let G = gList[g];
        if (Math.sqrt(mat.vec_dot(G, G)) < 1e-5) {
            continue;
        }
        let gEpsG = 0;
        for (let a = 0; a < 3; a++) {
            for (let b = 0; b < 3; b++) {
                gEpsG += G[a] * dielectric[a][b] * G[b];
            }
        }
        let prefactor = Math.exp(-gEpsG / l2) / gEpsG;
        for (let i = 0; i < natoms; i++) {
            for (let j = 0; j < natoms; j++) {
                let phase = 2 * Math.PI * (
                    G[0] * (positionsCar[i][0] - positionsCar[j][0]) +
                    G[1] * (positionsCar[i][1] - positionsCar[j][1]) +
                    G[2] * (positionsCar[i][2] - positionsCar[j][2])
                );
                let cos = Math.cos(phase);
                let sin = Math.sin(phase);
                // (Z_i^T G)(G^T Z_j) = outer(Z_i^T G, Z_j^T G)
                let zi = [0, 1, 2].map((alpha) => born[i][0][alpha] * G[0] + born[i][1][alpha] * G[1] + born[i][2][alpha] * G[2]);
                let zj = [0, 1, 2].map((beta) => born[j][0][beta] * G[0] + born[j][1][beta] * G[1] + born[j][2][beta] * G[2]);
                for (let alpha = 0; alpha < 3; alpha++) {
                    for (let beta = 0; beta < 3; beta++) {
                        let value = zi[alpha] * zj[beta] * prefactor;
                        ddQ0.real[i][alpha][beta] += value * cos;
                        ddQ0.imag[i][alpha][beta] += value * sin;
                    }
                }
            }
        }
    }
    for (let i = 0; i < natoms; i++) {
        let re = ddQ0.real[i];
        let im = ddQ0.imag[i];
        ddQ0.real[i] = [0, 1, 2].map((a) => [0, 1, 2].map((b) => (re[a][b] + re[b][a]) / 2));
        ddQ0.imag[i] = [0, 1, 2].map((a) => [0, 1, 2].map((b) => (im[a][b] - im[b][a]) / 2));
    }
    return ddQ0;
}

function buildMixedNacPayload(endpoint1, endpoint2, x, lattice, masses, positionsRed) {
    /*
    Gonze dipole-dipole parameters of the virtual crystal: Born charges and
    dielectric tensor mixed linearly on the mixed lattice
    */
    let nac1 = endpoint1.nac;
    let nac2 = endpoint2.nac;
    let mix = (a, b) => (1 - x) * a + x * b;
    let born = nac1.siteBorn.map((_, s) => [0, 1, 2].map((a) => [0, 1, 2].map((b) =>
        mix(nac1.siteBorn[s][a][b], nac2.siteBorn[s][a][b])
    )));
    let dielectric = [0, 1, 2].map((a) => [0, 1, 2].map((b) => mix(nac1.siteDielectric[a][b], nac2.siteDielectric[a][b])));
    return buildNacPayload(
        born, dielectric,
        mix(nac1.unitConversion, nac2.unitConversion),
        lattice, masses, utils.red_car_list(positionsRed, lattice),
        nac1.tolerance
    );
}

export function buildNacPayload(born, dielectric, unitConversion, lattice, masses, positionsCar, tolerance, lambdaOverride = null) {
    /*
    payload for getGonzeReciprocalCorrection: the G list, Ewald parameter, dd_q0
    and unit factor built for the given lattice with phonopy's own rules (or, with
    lambdaOverride, for that Ewald parameter, the G list cut where the Gaussian
    reaches 1e-10 as phonopy does)
    */
    let volume = Math.abs(mat.matrix_determinant(lattice));
    let reciprocal = utils.rec_lat(lattice);
    let trace = dielectric[0][0] + dielectric[1][1] + dielectric[2][2];
    // phonopy defaults: 300 G points, exp(-G eps G / 4 lambda^2) = 1e-10 at the cutoff
    let gCutoff = Math.cbrt(3 * 300 / (4 * Math.PI) / volume);
    let lambda = Math.sqrt(-(gCutoff * gCutoff) * trace / 3 / 4 / Math.log(1e-10));
    if (lambdaOverride) {
        lambda = lambdaOverride;
        gCutoff = lambda * Math.sqrt(-12 * Math.log(1e-10) / trace);
    }
    let gList = getGonzeGList(lattice, reciprocal, gCutoff);

    return {
        masses: masses,
        primitive_lattice: lattice,
        nac: {
            method: 'gonze',
            born: born,
            dielectric: dielectric,
            positions_car: positionsCar,
            g_list: gList,
            lambda: lambda,
            dd_q0: getGonzeDdQ0(gList, born, dielectric, positionsCar, lambda),
            nac_factor: unitConversion * 4 * Math.PI / volume,
            unit_conversion: unitConversion,
            q_direction_tolerance: tolerance,
        },
    };
}

function imposeAcousticSumRule(matrices, qpoints, masses, positionsRed, lattice = null) {
    /*
    Remove the residual acoustic-sum-rule violation of the mixed crystal.

    phonopy keeps the hermitian part of its dipole-dipole q = 0 correction, so for
    non-symmetric Born tensors that term carries a small sum-rule violation which
    the material's own short-range force constants cancel exactly. Short-range
    parts mix linearly but the dipole term is quadratic in Z, so in a mixture the
    two no longer cancel (zero at x = 0 and 1, and for symmetric Born tensors).

    With r_i = sum_j Phi_ij(Gamma) the row sums of the mixed force constants, the
    smallest hermitian correction with the same row sums is
        Delta_ij = (r_i + r_j^T) / N - (sum_k r_k) / N^2
    It is subtracted as intra-cell (R = 0) force constants, i.e. with the phase
    exp(2 pi i q.(tau_j - tau_i)) in the atom gauge, at every q. With `lattice`, on
    the shortest vector from i to j instead (averaged over equally short images), so
    that it does not depend on which lattice copy of each atom the positions name
    (a mixture of two structures, computeUnionMixedPhonon).
    */
    let gamma = qpoints.findIndex((q) => isGamma(q));
    if (gamma < 0) {
        return;
    }
    let n = masses.length;
    let rowSums = [];
    for (let i = 0; i < n; i++) {
        rowSums.push([0, 1, 2].map((a) => [0, 1, 2].map((b) => {
            let sum = 0;
            for (let j = 0; j < n; j++) {
                sum += Math.sqrt(masses[i] * masses[j]) * matrices[gamma].real[i * 3 + a][j * 3 + b];
            }
            return sum;
        })));
    }
    let total = [0, 1, 2].map((a) => [0, 1, 2].map((b) => rowSums.reduce((s, r) => s + r[a][b], 0)));
    let delta = (i, j, a, b) => (rowSums[i][a][b] + rowSums[j][b][a]) / n - total[a][b] / (n * n);

    // the lattice vectors to add to tau_j - tau_i: none, or the shortest images
    let images = [];
    for (let i = 0; i < n; i++) {
        images.push([]);
        for (let j = 0; j < n; j++) {
            if (!lattice) {
                images[i].push([[0, 0, 0]]);
                continue;
            }
            let d = positionsRed[j].map((v, k) => v - positionsRed[i][k]);
            let base = d.map((v) => -Math.round(v));
            // a skewed cell (the Hermite form of a common supercell) can put the
            // shortest image a few cells away from the rounded difference
            let candidates = [];
            for (let n0 = -3; n0 <= 3; n0++) {
                for (let n1 = -3; n1 <= 3; n1++) {
                    for (let n2 = -3; n2 <= 3; n2++) {
                        let shift = [base[0] + n0, base[1] + n1, base[2] + n2];
                        let cartesian = rowTimesMatrix(d.map((v, k) => v + shift[k]), lattice);
                        candidates.push({ shift: shift, length: mat.vec_dot(cartesian, cartesian) });
                    }
                }
            }
            let shortest = Math.min.apply(null, candidates.map((c) => c.length));
            images[i].push(candidates.filter((c) => c.length < shortest * (1 + 1e-6) + 1e-12).map((c) => c.shift));
        }
    }
    for (let k = 0; k < qpoints.length; k++) {
        let q = qpoints[k];
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                let cos = 0;
                let sin = 0;
                for (let shift of images[i][j]) {
                    let phase = 2 * Math.PI * (
                        q[0] * (positionsRed[j][0] - positionsRed[i][0] + shift[0]) +
                        q[1] * (positionsRed[j][1] - positionsRed[i][1] + shift[1]) +
                        q[2] * (positionsRed[j][2] - positionsRed[i][2] + shift[2])
                    );
                    cos += Math.cos(phase) / images[i][j].length;
                    sin += Math.sin(phase) / images[i][j].length;
                }
                let scale = 1 / Math.sqrt(masses[i] * masses[j]);
                for (let a = 0; a < 3; a++) {
                    for (let b = 0; b < 3; b++) {
                        let value = delta(i, j, a, b) * scale;
                        matrices[k].real[i * 3 + a][j * 3 + b] -= value * cos;
                        matrices[k].imag[i * 3 + a][j * 3 + b] -= value * sin;
                    }
                }
            }
        }
    }
}

function endpointLattice(endpoint) {
    // canonical lattice of an end-member, as computeMixedPhonon mixes it
    if (endpoint.structure && endpoint.structure.lattice) {
        return endpoint.structure.lattice;
    }
    if (endpoint.family === 'wurtzite') {
        return wurtziteLattice(endpoint.structure.a, endpoint.structure.c);
    }
    if (endpoint.family === 'chalcopyrite') {
        return chalcopyriteLattice(endpoint.structure.a, endpoint.structure.c);
    }
    return endpoint.lattice;
}

function subtractRealSpaceDipole(endpoint, matrices, qpoints) {
    /*
    The short-range parts so far are D minus phonopy's reciprocal-space dipole sum,
    whose Ewald parameter depends on the cell: the rest of the dipole interaction
    (the real-space part of the Ewald sum) is still in them. Take it out too, in the
    end-member's canonical cell, so that what is mixed linearly is D minus the full,
    lambda-independent dipole-dipole term (getGonzeDipoleDipole), which the mixture
    rebuilds from the mixed Born charges and dielectric tensor.
    */
    let lattice = endpointLattice(endpoint);
    let positions = endpointPositions(endpoint);
    let payload = buildNacPayload(endpoint.nac.siteBorn, endpoint.nac.siteDielectric, endpoint.nac.unitConversion,
        lattice, endpoint.siteMasses, utils.red_car_list(positions, lattice), endpoint.nac.tolerance);
    let factorSq = endpoint.frequencyFactor * endpoint.frequencyFactor;
    return matrices.map((matrix, k) => {
        let out = { real: matrix.real.map((row) => row.slice()), imag: matrix.imag.map((row) => row.slice()) };
        addBlockTensor(out, getGonzeRealSpaceRemainder(payload, qpoints[k]), -factorSq);
        return out;
    });
}

function endpointPositions(endpoint) {
    // canonical reduced positions of an end-member, as computeMixedPhonon mixes them
    if (endpoint.structure && endpoint.structure.lattice) {
        return endpoint.structure.positionsRed;
    }
    if (endpoint.family === 'wurtzite') {
        return wurtzitePositions(endpoint.structure.u);
    }
    if (endpoint.family === 'chalcopyrite') {
        return chalcopyritePositions(endpoint.structure.a, endpoint.structure.c, endpoint.structure.u);
    }
    return FCC_LAYOUTS[endpoint.family]
        ? FCC_LAYOUTS[endpoint.family].sites.map((site) => site.red.slice())
        : [[0, 0, 0], [0.25, 0.25, 0.25]];
}

function changeGauge(matrix, q, from, to) {
    /*
    D in the atom gauge of the positions `from` -> atom gauge of `to` (reduced, same
    cell and site order): D_ij * exp(2 pi i q.((to_j - to_i) - (from_j - from_i))).
    Through the cell gauge exp(2 pi i q.R) this keeps the force constants Phi_ij(R)
    and only moves the atoms.
    */
    let n = from.length;
    let shift = from.map((p, i) => 2 * Math.PI * (q[0] * (to[i][0] - p[0]) + q[1] * (to[i][1] - p[1]) + q[2] * (to[i][2] - p[2])));
    if (shift.every((v) => Math.abs(v) < 1e-12)) {
        return matrix;
    }
    let out = { real: matrix.real.map((row) => row.slice()), imag: matrix.imag.map((row) => row.slice()) };
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            let phase = shift[j] - shift[i];
            let cos = Math.cos(phase);
            let sin = Math.sin(phase);
            for (let a = 0; a < 3; a++) {
                for (let b = 0; b < 3; b++) {
                    let re = matrix.real[3 * i + a][3 * j + b];
                    let im = matrix.imag[3 * i + a][3 * j + b];
                    out.real[3 * i + a][3 * j + b] = re * cos - im * sin;
                    out.imag[3 * i + a][3 * j + b] = re * sin + im * cos;
                }
            }
        }
    }
    return out;
}

function presentInFileFrame(data, endpoint) {
    /*
    The mixture is computed in the canonical cell of material 1 (site order of the
    template, origin on its first site, rotated onto the template axes, and the
    inversion image when the file is stored in the opposite setting). Show it in
    the setting of material 1's own file instead: at x = 0 the page then draws
    exactly the material on its own, and the mixture changes continuously from it.

    With T the orthogonal map file -> canonical (the rotation, times -1 for the
    inversion image), the canonical lattice rows are C0 = U^-1 (F T^T) for the
    file rows F and an integer U. The mixed cell C_x, reduced positions p and q
    points go back as F_x = U C_x T, p U^-1 + (origin of the file), q in the
    reciprocal basis of F_x, displacements T^T e, atoms in the order of the file.
    The eigenvectors are in the atom gauge, so wrapping a position back into the
    cell needs no phase change. A binary folded into a supercell (material 1 a
    binary, material 2 an ordered ternary) is drawn in the supercell, but in the
    binary's own setting: its axes, its inversion setting, and every site on the
    position of an atom of the binary's file.
    */
    if (endpoint.anyQ) {
        let binary = endpoint.foldedFrom;
        if (!binary || binary.anyQ || !binary.positionsRed || !binary.order || !endpoint.sublattice) {
            return data;
        }
        let R = binary.rotation || [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let sign = binary.invert ? -1 : 1;
        let T = R.map((row) => row.map((v) => sign * v));
        let cell = endpointLattice(endpoint);
        let sites = endpointPositions(endpoint);
        // the supercell in the binary's cartesian frame, and the shift that puts a site
        // of sublattice 0 on the atom of the binary's file that is its site 0
        let fileCell = cell.map((row) => rowTimesMatrix(row, T));
        let s0 = endpoint.sublattice.indexOf(0);
        let fileAtom = rowTimesMatrix(binary.positionsRed[binary.order[0]], binary.lattice);
        let siteAtom = rowTimesMatrix(rowTimesMatrix(sites[s0], cell), T);
        let shift = rowTimesMatrix([0, 1, 2].map((k) => fileAtom[k] - siteAtom[k]), mat.matrix_inverse(fileCell));
        // a fold into a cell of the same size (the MgAgAs type antifluorite) is the
        // binary's own cell: drawn exactly as the file, its basis and atom order
        let sameSize = sites.length === binary.order.length;
        let order = sameSize ? endpoint.sublattice.map((sub) => binary.order[sub]) : sites.map((_, i) => i);
        let positionsRed = new Array(sites.length);
        sites.forEach((p, i) => { positionsRed[order[i]] = p.map((v, k) => v + shift[k]); });
        return presentInFileFrame(data, {
            lattice: sameSize ? binary.lattice : fileCell,
            positionsRed: sameSize ? binary.positionsRed : positionsRed,
            order: order,
            rotation: binary.rotation,
            invert: binary.invert,
            structure: { lattice: cell, positionsRed: sites },
            repetitions: data.repetitions,
        });
    }
    if (!endpoint.positionsRed || !endpoint.order) {
        return data;
    }
    let R = endpoint.rotation || [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let sign = endpoint.invert ? -1 : 1;
    let T = R.map((row) => row.map((v) => sign * v));
    let file = endpoint.lattice;
    let canonical0 = endpointLattice(endpoint);
    let positions0 = endpointPositions(endpoint);
    let rotatedFile = file.map((row) => matrixTimesVector(T, row));
    let Uraw = mat.matrix_multiply(rotatedFile, mat.matrix_inverse(canonical0));
    let U = Uraw.map((row) => row.map(Math.round));
    if (!isIntegerMatrix(Uraw, 1e-4) || Math.abs(Math.abs(mat.matrix_determinant(U)) - 1) > 1e-6) {
        return data;
    }
    let Uinverse = mat.matrix_inverse(U);
    let n = data.natoms;
    let order = endpoint.order;
    // the place of every atom in the file order; in a mixture of two structures
    // (site_index: its site in the union of both) the sites of material 1 come in the
    // order of its file, those only material 2 fills after them
    let slots = [];
    for (let s = 0; s < n; s++) {
        let site = data.site_index ? data.site_index[s] : s;
        slots.push(site < order.length ? order[site] : site);
    }
    let sortedSlots = slots.slice().sort((p, q) => p - q);
    let place = slots.map((slot) => sortedSlots.indexOf(slot));

    // rows of F_x: canonical rows in the file basis, then back to file cartesian axes
    let lattice = mat.matrix_multiply(U, data.lattice).map((row) => rowTimesMatrix(row, T));
    let toFileReduced = (p) => rowTimesMatrix(p, Uinverse);
    let origin = [0, 1, 2].map((k) => endpoint.positionsRed[order[0]][k] - toFileReduced(positions0[0])[k]);
    let positionsRed = new Array(n);
    let atomTypes = new Array(n);
    for (let s = 0; s < n; s++) {
        let p = toFileReduced(data.atom_pos_red[s]).map((v, k) => v + origin[k]);
        positionsRed[place[s]] = p.map((v) => v - Math.floor(v + 1e-9));
        atomTypes[place[s]] = data.atom_types[s];
    }

    let canonicalReciprocal = utils.rec_lat(data.lattice);
    let qpoints = data.qpoints.map((q) => {
        let cartesian = rowTimesMatrix(q, canonicalReciprocal);
        let fileCartesian = rowTimesMatrix(cartesian, T);
        return lattice.map((row) => mat.vec_dot(row, fileCartesian));
    });
    // T^T e for every atom: as a row vector, e T
    let back = (component) => [0, 1, 2].map((j) => component[0] * T[0][j] + component[1] * T[1][j] + component[2] * T[2][j]);
    let vectors = data.vectors.map((modes) => modes.map((atoms) => {
        let out = new Array(n);
        for (let s = 0; s < n; s++) {
            let re = back(atoms[s].map((c) => c[0]));
            let im = back(atoms[s].map((c) => c[1]));
            out[place[s]] = [0, 1, 2].map((k) => [re[k], im[k]]);
        }
        return out;
    }));

    return Object.assign({}, data, {
        atom_types: atomTypes,
        atom_numbers: atomTypes.map((type) => atomic_data.atomic_number[type]),
        atom_pos_red: positionsRed,
        atom_pos_car: utils.red_car_list(positionsRed, lattice),
        lattice: lattice,
        qpoints: qpoints,
        vectors: vectors,
        repetitions: endpoint.repetitions || data.repetitions,
    });
}

export async function computeMixedPhonon(sampled, x, massFraction = x) {
    /*
    diagonalize the virtual-crystal D_alloy (mixed force constants and masses)
    on the common path and return it in the internal json format of PhononJson.

    x mixes the bonding: force constants, Born charges, dielectric tensor and
    lattice. massFraction mixes the site masses, M = (1-y) M_1 + y M_2, and
    defaults to x (plain virtual crystal). The masses only divide the mixed
    force constants, so any choice keeps the acoustic sum rule.
    */
    if (sampled.union) {
        return computeUnionMixedPhonon(sampled, x, massFraction);
    }
    let endpoint1 = sampled.endpoint1;
    let endpoint2 = sampled.endpoint2;
    let natoms = endpoint1.siteTypes.length;

    let lattice;
    let positionsRed;
    if (endpoint1.structure && endpoint1.structure.lattice) {
        // ternary wurtzite or rocksalt: canonical lattice and reduced positions mixed linearly
        let s1 = endpoint1.structure;
        let s2 = endpoint2.structure;
        lattice = s1.lattice.map((row, i) => row.map((v, j) => (1 - x) * v + x * s2.lattice[i][j]));
        positionsRed = s1.positionsRed.map((row, i) => row.map((v, j) => (1 - x) * v + x * s2.positionsRed[i][j]));
    } else if (endpoint1.family === 'chalcopyrite' || endpoint1.family === 'wurtzite') {
        // virtual crystal: a, c and the anion parameter u mixed linearly
        let s1 = endpoint1.structure;
        let s2 = endpoint2.structure;
        let a = (1 - x) * s1.a + x * s2.a;
        let c = (1 - x) * s1.c + x * s2.c;
        let u = (1 - x) * s1.u + x * s2.u;
        if (endpoint1.family === 'wurtzite') {
            lattice = wurtziteLattice(a, c);
            positionsRed = wurtzitePositions(u);
        } else {
            lattice = chalcopyriteLattice(a, c);
            positionsRed = chalcopyritePositions(a, c, u);
        }
    } else {
        lattice = [0, 1, 2].map((i) => [0, 1, 2].map((j) =>
            (1 - x) * endpoint1.lattice[i][j] + x * endpoint2.lattice[i][j]
        ));
        positionsRed = FCC_LAYOUTS[endpoint1.family]
            ? FCC_LAYOUTS[endpoint1.family].sites.map((site) => site.red.slice())
            : [[0, 0, 0], [0.25, 0.25, 0.25]];
    }
    let masses = endpoint1.siteMasses.map((_, s) => (1 - massFraction) * endpoint1.siteMasses[s] + massFraction * endpoint2.siteMasses[s]);
    let atomTypes = x < 0.5 ? endpoint1.siteTypes.slice() : endpoint2.siteTypes.slice();

    let nacPayload = sampled.recomputeNac
        ? buildMixedNacPayload(endpoint1, endpoint2, x, lattice, masses, positionsRed)
        : null;
    let frequencyFactorSq = sampled.recomputeNac ? endpoint1.frequencyFactor * endpoint1.frequencyFactor : 1;

    // the virtual crystal mixes the force constants Phi(R), not the atom-gauge D(q):
    // each end-member is brought to the gauge of the mixed positions first (the two
    // differ by phases 2 pi q.(tau_1 - tau_2) when the internal positions differ:
    // u, tilts, ordered anions). Unchanged at Gamma.
    let positions1 = endpointPositions(endpoint1);
    let positions2 = endpointPositions(endpoint2);
    let gauge1 = (matrix, q) => changeGauge(matrix, q, positions1, positionsRed);
    let gauge2 = (matrix, q) => changeGauge(matrix, q, positions2, positionsRed);

    let matrices = sampled.qpoints.map((q, k) => {
        if (nacPayload) {
            // mixed short-range force constants + dipole-dipole term of the mixed crystal
            let matrix = mixForceConstantMatrices(
                gauge1(sampled.shortRange1[k], q), endpoint1.siteMasses,
                gauge2(sampled.shortRange2[k], q), endpoint2.siteMasses,
                masses, x
            );
            // the full Ewald sum, matching what was subtracted from the end-members
            let dd = getGonzeDipoleDipole(nacPayload, q, sampled.qDirections[k]);
            addBlockTensor(matrix, dd, frequencyFactorSq);
            return matrix;
        }
        return mixForceConstantMatrices(
            gauge1(sampled.matrices1[k], q), endpoint1.siteMasses,
            gauge2(sampled.matrices2[k], q), endpoint2.siteMasses,
            masses, x
        );
    });
    imposeAcousticSumRule(matrices, sampled.qpoints, masses, positionsRed);

    let eigenvalues = [];
    let vectors = [];
    for (let k = 0; k < sampled.qpoints.length; k++) {
        let solution = await diagonalize(matrices[k]);
        eigenvalues.push(solution.frequencies);

        // mass-weighted eigenvectors -> displacements with the mixed site masses
        let modes = [];
        for (let n = 0; n < solution.vectors.length; n++) {
            let vector = solution.vectors[n];
            let atoms = [];
            for (let a = 0; a < natoms; a++) {
                let invSqrtMass = 1 / Math.sqrt(masses[a]);
                atoms.push([0, 1, 2].map((d) => [
                    vector[a * 3 + d][0] * invSqrtMass,
                    vector[a * 3 + d][1] * invSqrtMass,
                ]));
            }
            modes.push(atoms);
        }
        vectors.push(modes);
    }

    // displacements in the convention of the files (PhononDB 'avg-mass-normalized':
    // e_i sqrt(<m> / m_i)), so the mixture is drawn like material 1 on its own
    let averageMass = masses.reduce((sum, mass) => sum + mass, 0) / natoms;
    let amplitude = Math.sqrt(averageMass);
    vectors.forEach((modes) => modes.forEach((atoms) => atoms.forEach((axes) => axes.forEach((component) => {
        component[0] *= amplitude;
        component[1] *= amplitude;
    }))));

    let data = {
        name: utils.get_formula(atomTypes),
        natoms: natoms,
        atom_types: atomTypes,
        atom_numbers: atomTypes.map((type) => atomic_data.atomic_number[type]),
        atom_pos_red: positionsRed,
        atom_pos_car: utils.red_car_list(positionsRed, lattice),
        lattice: lattice,
        formula: utils.get_formula(atomTypes),
        qpoints: sampled.qpoints,
        distances: getPathDistances(sampled.qpoints, sampled.lineBreaks, lattice),
        line_breaks: sampled.lineBreaks,
        highsym_qpts: getHighSymmetryPoints(sampled.lineBreaks, sampled.labels),
        eigenvalues: eigenvalues,
        vectors: vectors,
        repetitions: natoms > 4 ? [2, 2, 2] : [3, 3, 3],
        average_mass: averageMass,
        mode_amplitude_convention: 'avg-mass-normalized',
        family: endpoint1.family,
        // 'recomputed': dipole-dipole term from mixed Born charges and dielectric tensor
        // 'interpolated': the long-range term is mixed along with the force constants
        nac_mode: sampled.recomputeNac ? 'recomputed' : 'interpolated',
    };
    return presentInFileFrame(data, endpoint1);
}

async function computeUnionMixedPhonon(sampled, x, massFraction) {
    /*
    the virtual crystal of two structures on a common lattice (prepareWorldPair):
    every site either end-member fills, with a vacancy where only the other one
    does. A vacancy has no mass, force constants or charge, so a site filled only by
    material 2 carries, with its occupation x, the mass x M_2, force constants x Phi_2
    and Born charge x Z_2 (material 1: 1 - x). Its own block of D stays that of
    material 2 (a local mode of the atom in the lattice), its coupling to the other
    sites goes as sqrt(x), and it drops out at x = 0; the acoustic sum rule and charge
    neutrality hold for every x. Sites both fill mix as in computeMixedPhonon (masses
    with massFraction). At x = 0 and 1 the end-members are reproduced exactly.
    */
    let endpoint1 = sampled.endpoint1;
    let endpoint2 = sampled.endpoint2;
    let union = sampled.union;
    let map1 = union.map1;
    let map2 = union.map2;
    // positions of every union site in each end-member: a vacancy at the other's site
    let p1 = map1.map((i, s) => (i >= 0 ? union.positions1[i] : union.positions2[map2[s]]));
    let p2 = map2.map((i, s) => (i >= 0 ? union.positions2[i] : union.positions1[map1[s]]));
    let lattice = union.cell1.map((row, i) => row.map((v, j) => (1 - x) * v + x * union.cell2[i][j]));
    let mixedPositions = p1.map((p, s) => p.map((v, k) => (1 - x) * v + x * p2[s][k]));
    let masses = map1.map((i1, s) => {
        let i2 = map2[s];
        if (i1 >= 0 && i2 >= 0) {
            return (1 - massFraction) * endpoint1.siteMasses[i1] + massFraction * endpoint2.siteMasses[i2];
        }
        return i1 >= 0 ? (1 - x) * endpoint1.siteMasses[i1] : x * endpoint2.siteMasses[i2];
    });
    let active = [];
    masses.forEach((mass, s) => {
        if (mass > 0) {
            active.push(s);
        }
    });
    let natoms = active.length;
    let activeMasses = active.map((s) => masses[s]);
    let activePositions = active.map((s) => mixedPositions[s]);
    // the end-members' own sites at the mixed positions, for the gauge
    let ownTargets = (map, count) => {
        let out = new Array(count);
        map.forEach((i, s) => {
            if (i >= 0) {
                out[i] = mixedPositions[s];
            }
        });
        return out;
    };
    let to1 = ownTargets(map1, endpoint1.siteMasses.length);
    let to2 = ownTargets(map2, endpoint2.siteMasses.length);
    // D in the atom gauge does not change when an atom is relabelled by a lattice
    // vector (the dipole-dipole term neither): the gauge change moves each atom only
    // by its displacement less whole lattice vectors, taken relative to the common
    // shift of material 2 (prepareWorldPair), so that all atoms agree
    let gaugeFrom = (own, targets) => {
        let common = targets[0].map((v, k) => v - own[0][k]);
        return own.map((p, i) => p.map((v, k) => {
            let d = targets[i][k] - v - common[k];
            return targets[i][k] - common[k] - (d - Math.round(d));
        }));
    };
    let from1 = gaugeFrom(union.ownPositions1, to1);
    let from2 = gaugeFrom(union.ownPositions2, to2);

    let recompute = sampled.recomputeNac;
    let nacPayload = null;
    if (recompute) {
        let nac1 = endpoint1.nac;
        let nac2 = endpoint2.nac;
        let born = active.map((s) => [0, 1, 2].map((a) => [0, 1, 2].map((b) =>
            (map1[s] >= 0 ? (1 - x) * nac1.siteBorn[map1[s]][a][b] : 0) +
            (map2[s] >= 0 ? x * nac2.siteBorn[map2[s]][a][b] : 0)
        )));
        let dielectric = [0, 1, 2].map((a) => [0, 1, 2].map((b) => (1 - x) * nac1.siteDielectric[a][b] + x * nac2.siteDielectric[a][b]));
        nacPayload = buildNacPayload(born, dielectric, (1 - x) * nac1.unitConversion + x * nac2.unitConversion,
            lattice, activeMasses, utils.red_car_list(activePositions, lattice), nac1.tolerance);
    }
    let factorSq = recompute ? endpoint1.frequencyFactor * endpoint1.frequencyFactor : 1;

    let matrices = sampled.qpoints.map((q, k) => {
        let d1 = changeGauge(recompute ? sampled.shortRange1[k] : sampled.matrices1[k], q, from1, to1);
        let d2 = changeGauge(recompute ? sampled.shortRange2[k] : sampled.matrices2[k], q, from2, to2);
        let out = zeroMatrix(3 * natoms);
        let parts = [[map1, endpoint1, d1, 1 - x], [map2, endpoint2, d2, x]];
        for (let a = 0; a < natoms; a++) {
            for (let b = 0; b < natoms; b++) {
                let norm = 1 / Math.sqrt(masses[active[a]] * masses[active[b]]);
                for (let [map, endpoint, d, weight] of parts) {
                    let i = map[active[a]];
                    let j = map[active[b]];
                    if (i < 0 || j < 0 || weight === 0) {
                        continue;
                    }
                    let w = weight * Math.sqrt(endpoint.siteMasses[i] * endpoint.siteMasses[j]) * norm;
                    for (let alpha = 0; alpha < 3; alpha++) {
                        for (let beta = 0; beta < 3; beta++) {
                            out.real[3 * a + alpha][3 * b + beta] += w * d.real[3 * i + alpha][3 * j + beta];
                            out.imag[3 * a + alpha][3 * b + beta] += w * d.imag[3 * i + alpha][3 * j + beta];
                        }
                    }
                }
            }
        }
        if (nacPayload) {
            addBlockTensor(out, getGonzeDipoleDipole(nacPayload, q, sampled.qDirections[k]), factorSq);
        }
        return out;
    });
    imposeAcousticSumRule(matrices, sampled.qpoints, activeMasses, activePositions, lattice);

    let averageMass = activeMasses.reduce((sum, mass) => sum + mass, 0) / natoms;
    let amplitude = Math.sqrt(averageMass);
    let eigenvalues = [];
    let vectors = [];
    for (let k = 0; k < sampled.qpoints.length; k++) {
        let solution = await diagonalize(matrices[k]);
        eigenvalues.push(solution.frequencies);
        // displacements, in the convention of the files (see computeMixedPhonon)
        vectors.push(solution.vectors.map((vector) => activeMasses.map((mass, a) => [0, 1, 2].map((d) => [
            vector[a * 3 + d][0] * amplitude / Math.sqrt(mass),
            vector[a * 3 + d][1] * amplitude / Math.sqrt(mass),
        ]))));
    }
    // the species drawn on a site: the one there, on a site both fill the majority one
    let atomTypes = active.map((s) => {
        let i1 = map1[s];
        let i2 = map2[s];
        if (i1 >= 0 && (i2 < 0 || x < 0.5)) {
            return endpoint1.siteTypes[i1];
        }
        return endpoint2.siteTypes[i2];
    });
    let data = {
        name: utils.get_formula(atomTypes),
        natoms: natoms,
        atom_types: atomTypes,
        atom_numbers: atomTypes.map((type) => atomic_data.atomic_number[type]),
        atom_pos_red: activePositions,
        atom_pos_car: utils.red_car_list(activePositions, lattice),
        lattice: lattice,
        formula: utils.get_formula(atomTypes),
        qpoints: sampled.qpoints,
        distances: getPathDistances(sampled.qpoints, sampled.lineBreaks, lattice),
        line_breaks: sampled.lineBreaks,
        highsym_qpts: getHighSymmetryPoints(sampled.lineBreaks, sampled.labels),
        eigenvalues: eigenvalues,
        vectors: vectors,
        repetitions: natoms > 4 ? [2, 2, 2] : [3, 3, 3],
        average_mass: averageMass,
        mode_amplitude_convention: 'avg-mass-normalized',
        family: endpoint1.family,
        nac_mode: recompute ? 'recomputed' : 'interpolated',
        // the union site of every atom: material 1's come first, in its order
        site_index: active,
        // sites filled in only one of the two materials at this x
        vacancy_sites: active.filter((s) => map1[s] < 0 || map2[s] < 0).length,
    };
    return presentInFileFrame(data, endpoint1);
}

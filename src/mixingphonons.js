/*
 * mixingphonons.js — virtual-crystal mixing of two materials of the same
 * structure family on the mixing page: 2-atom zinc blende / diamond, or 8-atom
 * chalcopyrite (I-42d). Force constants and site masses are mixed separately:
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
 * Long-range (LO-TO) part: only the short-range force constants are mixed, and
 * the dipole-dipole term is recomputed for the virtual crystal from the
 * linearly mixed Born charges and dielectric tensor on the mixed lattice. For
 * PhononDB the stored (Gonze) force constants already are the short-range part.
 * For Materials Project the Born charges and dielectric tensor come from the MP
 * DFPT collection (mpdielectric.js), and the material's own dipole-dipole term
 * is subtracted from D rebuilt at the q-points of the file. If they are not
 * available the long-range term is mixed along with the rest of Phi.
 *
 * Before mixing, both end-members are brought to a common representation:
 *
 *   gauge  : "atom" Bloch phase exp(2 pi i q.(R + tau_j - tau_i)), which is what
 *            phonopy force constants give (PhononDB). Materials Project
 *            eigendisplacements are in the "cell" gauge exp(2 pi i q.R) and are
 *            converted with the atomic positions of the file.
 *   sites  : zinc blende: canonical order [cation, anion] with the anion at
 *            +(1/4,1/4,1/4) from the cation. Chalcopyrite: the 8 sites of the
 *            I-42d template (lower-valence cation on 4a, the other on 4b, anions
 *            on 8d), matched with any origin. A structure stored as the
 *            inversion image of the template is handled in the atom gauge by
 *            D(q) -> D(q)^*; every other setting of the same lattice differs
 *            from these by a symmetry of the crystal.
 *   q-path : the high-symmetry segments of the end-member with the fewer
 *            segments, kept where the other end-member has a segment on the
 *            same line (in either direction) and cut to the stretch both cover.
 *            For zinc blende the segments coincide; for chalcopyrite the zone-
 *            boundary points Sigma, Y, ... move with c/a, so segments ending
 *            there are cut at the nearer of the two and segments whose line
 *            itself moves are dropped. Both end-members are sampled at the same
 *            reduced q. All files of a family share the same primitive vectors,
 *            so reduced coordinates compare directly.
 *
 * Swapping the end-members and x -> 1-x gives the same D_alloy, since the path
 * choice, the sampling and the lattice interpolation are all symmetric.
 */

import { buildDynamicalMatrixBlocks } from './dynamicalmatrix.js';
import { getGonzeReciprocalCorrection } from './gonze.js';
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

function getCationSymbol(name, atomTypes) {
    let match = /^([A-Z][a-z]?)/.exec(name || '');
    if (match && atomTypes.indexOf(match[1]) !== -1) {
        return match[1];
    }
    return atomTypes[0];
}

function getCanonicalSites(atomTypes, positionsRed, cationSymbol) {
    /*
    returns the original atom indices in [cation, anion] order and whether the
    structure has to be inverted to put the anion at +(1/4,1/4,1/4)
    */
    let isElemental = atomTypes[0] === atomTypes[1];
    let cation = isElemental ? 0 : atomTypes.indexOf(cationSymbol);
    let anion = 1 - cation;

    let offset = [0, 1, 2].map((i) => {
        let d = positionsRed[anion][i] - positionsRed[cation][i];
        return d - Math.floor(d + 1e-6);
    });
    let isPlus = offset.every((d) => Math.abs(d - 0.25) < 1e-3);
    let isMinus = offset.every((d) => Math.abs(d - 0.75) < 1e-3);
    if (!isPlus && !isMinus) {
        throw new Error('Not a 2-atom zinc blende / diamond structure: anion offset ' + JSON.stringify(offset));
    }
    if (isMinus && isElemental) {
        // both sites are equivalent, relabel instead of inverting
        return { order: [anion, cation], invert: false };
    }
    return { order: [cation, anion], invert: isMinus };
}

// usual cation valence, used to put the lower-valence cation of a chalcopyrite
// (Cu in CuGaS2, Zn in ZnGeP2) on the 4a site and the other one on 4b
const CATION_VALENCE = {
    Li: 1, Na: 1, K: 1, Cu: 1, Ag: 1, Au: 1,
    Be: 2, Mg: 2, Ca: 2, Zn: 2, Cd: 2, Hg: 2,
    B: 3, Al: 3, Ga: 3, In: 3, Tl: 3,
    C: 4, Si: 4, Ge: 4, Sn: 4, Pb: 4,
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
    if (dielectric && !sameSites && !uniqueSpecies && atomTypes[0] !== atomTypes[1]) {
        console.log("Materials Project Born charges listed in another site order, interpolating the LO-TO term");
        dielectric = null;
    }
    if (dielectric) {
        let born = atomTypes.map((type, i) => {
            let index = sameSites || atomTypes[0] === atomTypes[1] ? i : dielectric.labels.indexOf(type);
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
    if (endpoint.atomTypes.length === 2) {
        let cationSymbol = getCationSymbol(entry && entry.name, endpoint.atomTypes);
        canonical = getCanonicalSites(endpoint.atomTypes, endpoint.positionsRed, cationSymbol);
        endpoint.family = 'zincblende';
        endpoint.rotation = null;
    } else if (endpoint.atomTypes.length === 8) {
        canonical = matchChalcopyriteSites(endpoint.atomTypes, endpoint.positionsRed, endpoint.lattice);
        endpoint.family = 'chalcopyrite';
        endpoint.rotation = canonical.rotation;
        endpoint.structure = canonical.structure;
        relabelSegments(endpoint.segments, getChalcopyriteSpecialPoints(canonical.structure.a, canonical.structure.c));
    } else {
        throw new Error('Mixing supports 2-atom zinc blende and 8-atom chalcopyrite cells');
    }
    endpoint.name = entry && entry.name ? entry.name : utils.get_formula(endpoint.atomTypes);
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

export function zincblendeAsChalcopyrite(zb) {
    /*
    A zinc blende expressed in the 8-atom chalcopyrite cell, so it can be mixed
    with a chalcopyrite. With c = 2a and u = 1/4 the chalcopyrite template is
    zinc blende: the A and B sites lie on the cation sublattice, the anion sites
    on the anion sublattice, and the body-centred tetragonal primitive cell holds
    4 zinc blende primitive cells. D of this supercell follows from band folding,

        D^S_IJ(Q) = 1/4 sum_g D_ij(Q + g) exp(-2 pi i g.(r_J - r_I))

    over the 4 supercell reciprocal vectors g that differ modulo the zinc blende
    reciprocal lattice (i, j: sublattices of the supercell sites I, J; atom gauge,
    D_ij(k) = sum_R Phi_ij(R) exp(2 pi i k.(R + tau_j - tau_i))). With this sign
    the sum does not depend on which representative g of each class is used.
    Q + g leaves the band path, so this needs force constants (PhononDB).
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
    let lattice = chalcopyriteLattice(a, c);
    let positionsCar = utils.red_car_list(chalcopyritePositions(a, c, u), lattice);
    // template sites 0-3 (A, B) are cations, 4-7 anions
    let sublattice = CHALCOPYRITE_SITES.map((site) => (site.role === 'X' ? 1 : 0));
    let superReciprocal = utils.rec_lat(lattice);

    // folding vectors: supercell reciprocal vectors, distinct modulo the fcc reciprocal
    // lattice; the shortest of each class, so that Q = 0 maps onto k = 0 exactly
    let candidates = [];
    for (let n0 = -2; n0 <= 2; n0++) {
        for (let n1 = -2; n1 <= 2; n1++) {
            for (let n2 = -2; n2 <= 2; n2++) {
                candidates.push([0, 1, 2].map((d) => n0 * superReciprocal[0][d] + n1 * superReciprocal[1][d] + n2 * superReciprocal[2][d]));
            }
        }
    }
    candidates.sort((g, h) => mat.vec_dot(g, g) - mat.vec_dot(h, h));
    let folding = [];
    let seen = new Set();
    for (let g of candidates) {
        let key = primitive.map((p) => {
            let f = mat.vec_dot(g, p);
            return Math.round((f - Math.floor(f + 1e-9)) * 1e6) % 1000000;
        }).join(',');
        if (!seen.has(key)) {
            seen.add(key);
            folding.push(g);
        }
    }
    if (folding.length !== 4) {
        throw new Error('Expected 4 folding vectors, found ' + folding.length);
    }

    let toCartesian = (q) => [0, 1, 2].map((d) => q[0] * superReciprocal[0][d] + q[1] * superReciprocal[1][d] + q[2] * superReciprocal[2][d]);
    let toPrimitive = (k) => primitive.map((p) => mat.vec_dot(k, p));

    let sampleAtQ = function(q, direction) {
        let qCart = toCartesian(q);
        let directionPrimitive = direction ? toPrimitive(toCartesian(direction)) : null;
        let out = zeroMatrix(24);
        for (let g of folding) {
            let k = [0, 1, 2].map((d) => qCart[d] + g[d]);
            let block = toCanonical(zb, zb.evaluate(toPrimitive(k), directionPrimitive));
            for (let I = 0; I < 8; I++) {
                for (let J = 0; J < 8; J++) {
                    let phase = -2 * Math.PI * (
                        g[0] * (positionsCar[J][0] - positionsCar[I][0]) +
                        g[1] * (positionsCar[J][1] - positionsCar[I][1]) +
                        g[2] * (positionsCar[J][2] - positionsCar[I][2])
                    );
                    let cos = Math.cos(phase) / 4;
                    let sin = Math.sin(phase) / 4;
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

    let siteTypes = sublattice.map((s) => zb.siteTypes[s]);
    let siteMasses = sublattice.map((s) => zb.siteMasses[s]);
    let nac = null;
    let sampleShortRangeAtQ = null;
    if (zb.nac) {
        nac = {
            siteBorn: sublattice.map((s) => zb.nac.siteBorn[s]),
            siteDielectric: zb.nac.siteDielectric,
            unitConversion: zb.nac.unitConversion,
            tolerance: zb.nac.tolerance,
        };
        // short-range part: folded D minus the dipole-dipole term of the supercell built
        // with the same rules as the mixture's, so the end-member is reproduced exactly
        let ownPayload = buildNacPayload(nac.siteBorn, nac.siteDielectric, nac.unitConversion,
            lattice, siteMasses, positionsCar, nac.tolerance);
        let factorSq = zb.frequencyFactor * zb.frequencyFactor;
        sampleShortRangeAtQ = function(q, direction) {
            let matrix = sampleAtQ(q, direction);
            let dd = getGonzeReciprocalCorrection(ownPayload, q, isGamma(q) ? direction : null);
            addBlockTensor(matrix, dd, -factorSq);
            return matrix;
        };
    }

    return {
        family: 'chalcopyrite',
        anyQ: true,
        name: zb.name,
        lattice: lattice,
        structure: { a: a, c: c, u: u },
        siteTypes: siteTypes,
        siteMasses: siteMasses,
        order: siteTypes.map((_, i) => i),
        invert: false,
        rotation: null,
        nac: nac,
        frequencyFactor: zb.frequencyFactor,
        segments: [],
        sampleAtQ: sampleAtQ,
        sampleShortRangeAtQ: sampleShortRangeAtQ,
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
    for (let other of others) {
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

export function buildCommonPath(endpoint1, endpoint2) {
    /*
    the segments of the end-member with the fewer segments (endpoint1 on a tie),
    each kept where the other end-member has a segment on the same line, in
    either direction, and cut to the stretch both cover. For zinc blende the
    segments coincide; for chalcopyrite the zone-boundary points Sigma, Y, ...
    move with c/a, so a segment ending there is cut at the nearer of the two.
    */
    if (endpoint1.anyQ || endpoint2.anyQ) {
        // a zinc blende folded into the chalcopyrite cell can be evaluated at any q:
        // the path is the whole path of the chalcopyrite
        let real = endpoint1.anyQ ? endpoint2 : endpoint1;
        let folded = endpoint1.anyQ ? endpoint1 : endpoint2;
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

    let lesser = endpoint2.segments.length < endpoint1.segments.length ? endpoint2 : endpoint1;
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

export function sampleMixingPair(endpoint1, endpoint2) {
    /*
    canonical D_1 and D_2 at every q-point of the common path; independent of x.
    When both end-members carry Born charges and dielectric tensors (PhononDB),
    their short-range parts are kept too, so the dipole-dipole term can be
    recomputed for the mixed crystal instead of being interpolated.
    */
    if (endpoint1.family !== endpoint2.family) {
        let families = [endpoint1.family, endpoint2.family].sort().join('+');
        if (families !== 'chalcopyrite+zincblende') {
            throw new Error('Cannot mix a ' + endpoint1.family + ' cell with a ' + endpoint2.family + ' cell');
        }
        // the zinc blende goes into the chalcopyrite cell
        if (endpoint1.family === 'zincblende') {
            endpoint1 = zincblendeAsChalcopyrite(endpoint1);
        } else {
            endpoint2 = zincblendeAsChalcopyrite(endpoint2);
        }
    }
    let recomputeNac = !!(endpoint1.nac && endpoint2.nac);
    let shortRange1 = [];
    let shortRange2 = [];
    let qDirections = [];
    let path = buildCommonPath(endpoint1, endpoint2);
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

export function buildNacPayload(born, dielectric, unitConversion, lattice, masses, positionsCar, tolerance) {
    /*
    payload for getGonzeReciprocalCorrection: the G list, Ewald parameter, dd_q0
    and unit factor built for the given lattice with phonopy's own rules
    */
    let volume = Math.abs(mat.matrix_determinant(lattice));
    let reciprocal = utils.rec_lat(lattice);
    // phonopy defaults: 300 G points, exp(-G eps G / 4 lambda^2) = 1e-10 at the cutoff
    let gCutoff = Math.cbrt(3 * 300 / (4 * Math.PI) / volume);
    let gList = getGonzeGList(lattice, reciprocal, gCutoff);
    let trace = dielectric[0][0] + dielectric[1][1] + dielectric[2][2];
    let lambda = Math.sqrt(-(gCutoff * gCutoff) * trace / 3 / 4 / Math.log(1e-10));

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
            q_direction_tolerance: tolerance,
        },
    };
}

function imposeAcousticSumRule(matrices, qpoints, masses, positionsRed) {
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
    exp(2 pi i q.(tau_j - tau_i)) in the atom gauge, at every q.
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

    for (let k = 0; k < qpoints.length; k++) {
        let q = qpoints[k];
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                let phase = 2 * Math.PI * (
                    q[0] * (positionsRed[j][0] - positionsRed[i][0]) +
                    q[1] * (positionsRed[j][1] - positionsRed[i][1]) +
                    q[2] * (positionsRed[j][2] - positionsRed[i][2])
                );
                let cos = Math.cos(phase);
                let sin = Math.sin(phase);
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

export async function computeMixedPhonon(sampled, x, massFraction = x) {
    /*
    diagonalize the virtual-crystal D_alloy (mixed force constants and masses)
    on the common path and return it in the internal json format of PhononJson.

    x mixes the bonding: force constants, Born charges, dielectric tensor and
    lattice. massFraction mixes the site masses, M = (1-y) M_1 + y M_2, and
    defaults to x (plain virtual crystal). The masses only divide the mixed
    force constants, so any choice keeps the acoustic sum rule.
    */
    let endpoint1 = sampled.endpoint1;
    let endpoint2 = sampled.endpoint2;
    let natoms = endpoint1.siteTypes.length;

    let lattice;
    let positionsRed;
    if (endpoint1.family === 'chalcopyrite') {
        // virtual crystal: a, c and the anion parameter u mixed linearly
        let s1 = endpoint1.structure;
        let s2 = endpoint2.structure;
        let a = (1 - x) * s1.a + x * s2.a;
        let c = (1 - x) * s1.c + x * s2.c;
        lattice = chalcopyriteLattice(a, c);
        positionsRed = chalcopyritePositions(a, c, (1 - x) * s1.u + x * s2.u);
    } else {
        lattice = [0, 1, 2].map((i) => [0, 1, 2].map((j) =>
            (1 - x) * endpoint1.lattice[i][j] + x * endpoint2.lattice[i][j]
        ));
        positionsRed = [[0, 0, 0], [0.25, 0.25, 0.25]];
    }
    let masses = endpoint1.siteMasses.map((_, s) => (1 - massFraction) * endpoint1.siteMasses[s] + massFraction * endpoint2.siteMasses[s]);
    let atomTypes = x < 0.5 ? endpoint1.siteTypes.slice() : endpoint2.siteTypes.slice();

    let nacPayload = sampled.recomputeNac
        ? buildMixedNacPayload(endpoint1, endpoint2, x, lattice, masses, positionsRed)
        : null;
    let frequencyFactorSq = sampled.recomputeNac ? endpoint1.frequencyFactor * endpoint1.frequencyFactor : 1;

    let matrices = sampled.qpoints.map((q, k) => {
        if (nacPayload) {
            // mixed short-range force constants + dipole-dipole term of the mixed crystal
            let matrix = mixForceConstantMatrices(
                sampled.shortRange1[k], endpoint1.siteMasses,
                sampled.shortRange2[k], endpoint2.siteMasses,
                masses, x
            );
            let dd = getGonzeReciprocalCorrection(nacPayload, q, sampled.qDirections[k]);
            addBlockTensor(matrix, dd, frequencyFactorSq);
            return matrix;
        }
        return mixForceConstantMatrices(
            sampled.matrices1[k], endpoint1.siteMasses,
            sampled.matrices2[k], endpoint2.siteMasses,
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

    return {
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
        repetitions: natoms > 2 ? [2, 2, 2] : [3, 3, 3],
        family: endpoint1.family,
        // 'recomputed': dipole-dipole term from mixed Born charges and dielectric tensor
        // 'interpolated': the long-range term is mixed along with the force constants
        nac_mode: sampled.recomputeNac ? 'recomputed' : 'interpolated',
    };
}

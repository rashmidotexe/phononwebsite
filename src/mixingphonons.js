/*
 * mixingphonons.js — virtual-crystal mixing of two 2-atom zinc blende / diamond
 * materials on the mixing page. Force constants and site masses are mixed
 * separately:
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
 * Long-range (LO-TO) part: for two PhononDB materials only the short-range
 * (Gonze) force constants are mixed, and the dipole-dipole term is recomputed
 * for the virtual crystal from the linearly mixed Born charges and dielectric
 * tensor on the mixed lattice. Materials Project files have no Born charges or
 * dielectric tensor, so with one of them the long-range term is mixed along
 * with the rest of Phi.
 *
 * Before mixing, both end-members are brought to a common representation:
 *
 *   gauge  : "atom" Bloch phase exp(2 pi i q.(R + tau_j - tau_i)), which is what
 *            phonopy force constants give (PhononDB). Materials Project
 *            eigendisplacements are in the "cell" gauge exp(2 pi i q.R) and are
 *            converted with the atomic positions of the file.
 *   sites  : canonical order [cation, anion] with the anion at +(1/4,1/4,1/4)
 *            from the cation. A structure stored with the anion at -(1/4,1/4,1/4)
 *            is the inversion image; in the atom gauge inversion maps
 *            D(q) -> D(q)^*.
 *   q-path : only the high-symmetry segments of the end-member with the fewer
 *            segments that are also present (in either direction) in the other
 *            end-member. Every segment is resampled at the same fractions t in
 *            [0,1] for both end-members, so D_1 and D_2 are always taken at the
 *            same reduced q. All 2-atom fcc files share the same primitive
 *            vectors, so reduced coordinates compare directly.
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
    };
}

function buildMaterialsProjectEndpoint(raw) {
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
        return makeSegment(qpoints, segment, function(t) {
            // linear interpolation of D between the sampled q-points of the file
            let k = 0;
            while (k < fractions.length - 2 && fractions[k + 1] < t) {
                k++;
            }
            let t0 = fractions[k];
            let t1 = fractions[k + 1];
            let w = t1 > t0 ? Math.min(1, Math.max(0, (t - t0) / (t1 - t0))) : 0;
            return lerpMatrix(getMatrix(segment.start + k), getMatrix(segment.start + k + 1), w);
        });
    });

    return {
        atomTypes: atomTypes,
        positionsRed: positionsRed,
        lattice: structure.lattice.matrix,
        masses: masses,
        segments: segments,
    };
}

export function buildMixingEndpoint(raw, entry) {
    /*
    raw is the parsed json of a PhononDB (internal format with a dynamical_matrix)
    or a Materials Project OpenData phonon file; entry is its menu entry
    */
    let endpoint;
    if (raw.dynamical_matrix) {
        endpoint = buildPhononDBEndpoint(raw, entry);
    } else if (raw.eigendisplacements && raw.structure) {
        endpoint = buildMaterialsProjectEndpoint(raw);
    } else {
        throw new Error('No dynamical matrix available for ' + (entry && entry.name));
    }
    if (endpoint.atomTypes.length !== 2) {
        throw new Error('Mixing needs a 2-atom primitive cell');
    }

    let cationSymbol = getCationSymbol(entry && entry.name, endpoint.atomTypes);
    let canonical = getCanonicalSites(endpoint.atomTypes, endpoint.positionsRed, cationSymbol);
    endpoint.name = entry && entry.name ? entry.name : utils.get_formula(endpoint.atomTypes);
    endpoint.order = canonical.order;
    endpoint.invert = canonical.invert;
    endpoint.siteTypes = canonical.order.map((i) => endpoint.atomTypes[i]);
    endpoint.siteMasses = canonical.order.map((i) => endpoint.masses[i]);
    if (endpoint.nac) {
        // Born charges are rank-2 tensors, unchanged by the inversion
        endpoint.nac.siteBorn = canonical.order.map((i) => endpoint.nac.born[i]);
    }
    return endpoint;
}

function toCanonical(endpoint, matrix) {
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
    return out;
}

function findSegment(segments, qA, qB) {
    for (let i = 0; i < segments.length; i++) {
        if (sameQ(segments[i].qA, qA) && sameQ(segments[i].qB, qB)) {
            return { segment: segments[i], reversed: false };
        }
    }
    for (let i = 0; i < segments.length; i++) {
        if (sameQ(segments[i].qA, qB) && sameQ(segments[i].qB, qA)) {
            return { segment: segments[i], reversed: true };
        }
    }
    return null;
}

export function buildCommonPath(endpoint1, endpoint2) {
    /*
    segments of the end-member with the fewer segments (endpoint1 on a tie)
    that the other end-member also has, in either direction
    */
    let lesser = endpoint2.segments.length < endpoint1.segments.length ? endpoint2 : endpoint1;
    let other = lesser === endpoint1 ? endpoint2 : endpoint1;

    let path = [];
    for (let i = 0; i < lesser.segments.length; i++) {
        let segment = lesser.segments[i];
        let match = findSegment(other.segments, segment.qA, segment.qB);
        if (!match) {
            continue;
        }
        let orient = (sampler) => (sampler && match.reversed ? (t) => sampler(1 - t) : sampler);
        let sampleLesser = segment.sample;
        let sampleOther = orient(match.segment.sample);
        let shortRangeLesser = segment.sampleShortRange;
        let shortRangeOther = orient(match.segment.sampleShortRange);
        path.push({
            qA: segment.qA,
            qB: segment.qB,
            labelA: segment.labelA,
            labelB: segment.labelB,
            npoints: Math.max(segment.npoints, match.segment.npoints),
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
    phonopy subtracts from the diagonal blocks to keep the acoustic sum rule
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
    return ddQ0;
}

function buildMixedNacPayload(endpoint1, endpoint2, x, lattice, masses, positionsRed) {
    /*
    Gonze dipole-dipole parameters of the virtual crystal: Born charges and
    dielectric tensor mixed linearly, the G list, Ewald parameter, dd_q0 and
    unit factor rebuilt for the mixed lattice with phonopy's own rules
    */
    let nac1 = endpoint1.nac;
    let nac2 = endpoint2.nac;
    let mix = (a, b) => (1 - x) * a + x * b;
    let born = [0, 1].map((s) => [0, 1, 2].map((a) => [0, 1, 2].map((b) =>
        mix(nac1.siteBorn[s][a][b], nac2.siteBorn[s][a][b])
    )));
    let dielectric = [0, 1, 2].map((a) => [0, 1, 2].map((b) => mix(nac1.dielectric[a][b], nac2.dielectric[a][b])));
    let positionsCar = utils.red_car_list(positionsRed, lattice);

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
            nac_factor: mix(nac1.unitConversion, nac2.unitConversion) * 4 * Math.PI / volume,
            q_direction_tolerance: nac1.tolerance,
        },
    };
}

export async function computeMixedPhonon(sampled, x) {
    /*
    diagonalize the virtual-crystal D_alloy (mixed force constants and masses)
    on the common path and return it in the internal json format of PhononJson
    */
    let endpoint1 = sampled.endpoint1;
    let endpoint2 = sampled.endpoint2;
    let natoms = 2;

    let lattice = [0, 1, 2].map((i) => [0, 1, 2].map((j) =>
        (1 - x) * endpoint1.lattice[i][j] + x * endpoint2.lattice[i][j]
    ));
    let masses = [0, 1].map((s) => (1 - x) * endpoint1.siteMasses[s] + x * endpoint2.siteMasses[s]);
    let atomTypes = x < 0.5 ? endpoint1.siteTypes.slice() : endpoint2.siteTypes.slice();
    let positionsRed = [[0, 0, 0], [0.25, 0.25, 0.25]];

    let nacPayload = sampled.recomputeNac
        ? buildMixedNacPayload(endpoint1, endpoint2, x, lattice, masses, positionsRed)
        : null;
    let frequencyFactorSq = sampled.recomputeNac ? endpoint1.frequencyFactor * endpoint1.frequencyFactor : 1;

    let eigenvalues = [];
    let vectors = [];
    for (let k = 0; k < sampled.qpoints.length; k++) {
        let matrix;
        if (nacPayload) {
            // mixed short-range force constants + dipole-dipole term of the mixed crystal
            matrix = mixForceConstantMatrices(
                sampled.shortRange1[k], endpoint1.siteMasses,
                sampled.shortRange2[k], endpoint2.siteMasses,
                masses, x
            );
            let dd = getGonzeReciprocalCorrection(nacPayload, sampled.qpoints[k], sampled.qDirections[k]);
            addBlockTensor(matrix, dd, frequencyFactorSq);
        } else {
            matrix = mixForceConstantMatrices(
                sampled.matrices1[k], endpoint1.siteMasses,
                sampled.matrices2[k], endpoint2.siteMasses,
                masses, x
            );
        }
        let solution = await diagonalize(matrix);
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
        repetitions: [3, 3, 3],
        // 'recomputed': dipole-dipole term from mixed Born charges and dielectric tensor
        // 'interpolated': the long-range term is mixed along with the force constants
        nac_mode: sampled.recomputeNac ? 'recomputed' : 'interpolated',
    };
}

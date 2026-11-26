import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, it } from 'node:test';

import {
    analyzeSuperstructure,
    buildMixingEndpoint,
    sampleMixingPair,
    computeMixedPhonon,
    computeReferenceEigenvalues,
    binaryAsSuperstructure,
    buildNacPayload,
} from '../src/mixingphonons.js';
import { getGonzeDipoleDipole, getGonzeReciprocalCorrection } from '../src/gonze.js';
import { solveComplexHermitianWithEigenWasm } from '../src/eigenwasm.js';
import * as utils from '../src/utils.js';

// perovskite, dichalcogenide and Materials Project path tests, and the physics checks
// of the mixing (dipole-dipole split, gauge); split from mixing-phonons.test.mjs so
// that each file stays well within the test runner's time per file

const __dirname = dirname(fileURLToPath(import.meta.url));

function loadFixture(name) {
    return JSON.parse(gunzipSync(readFileSync(join(__dirname, 'fixtures', 'mixing', name))));
}

function endpoint(file, name) {
    return buildMixingEndpoint(loadFixture(file), { name: name });
}

function maxAbsDiff(a, b) {
    let worst = 0;
    for (let k = 0; k < a.length; k++) {
        for (let n = 0; n < a[k].length; n++) {
            worst = Math.max(worst, Math.abs(a[k][n] - b[k][n]));
        }
    }
    return worst;
}

describe('perovskite dynamical-matrix mixing', () => {
    let gammaIndices = (sampled) => sampled.qpoints
        .map((q, k) => (q.every((v) => Math.abs(v - Math.round(v)) < 1e-9) ? k : -1))
        .filter((k) => k >= 0);
    let opticalDiff = (a, b) => {
        let worst = 0;
        a.forEach((row, k) => row.forEach((v, n) => {
            if (Math.abs(v) > 1) {
                worst = Math.max(worst, Math.abs(v - b[k][n]));
            }
        }));
        return worst;
    };
    // the three modes closest to zero at Gamma: cubic perovskites have soft
    // (imaginary) modes there besides the acoustic ones
    let acoustic = (phonon, k) => phonon.eigenvalues[k].map(Math.abs).sort((p, q) => p - q)[2];
    let eigen = async (m) => (await solveComplexHermitianWithEigenWasm(m.real, m.imag)).values.sort((p, q) => p - q);

    it('tells A and B apart by the octahedra, in cubic and polar cells', () => {
        let srtio3 = endpoint('phonondb-SrTiO3-mp-5229.json.gz', 'SrTiO3');
        assert.equal(srtio3.family, 'perovskite');
        assert.deepEqual(srtio3.siteTypes, ['Sr', 'Ti', 'O', 'O', 'O']);
        // tetragonal PbTiO3: Ti off-centre in its octahedron, still B
        let pbtio3 = endpoint('phonondb-PbTiO3-mp-20459.json.gz', 'PbTiO3');
        assert.equal(pbtio3.family, 'perovskite');
        assert.deepEqual(pbtio3.siteTypes, ['Pb', 'Ti', 'O', 'O', 'O']);
        // antiperovskite: O in the Li6 octahedra, Cl on the A site
        let li3clo = endpoint('phonondb-Li3ClO-mp-985585.json.gz', 'Li3ClO');
        assert.equal(li3clo.family, 'antiperovskite');
        assert.deepEqual(li3clo.siteTypes, ['Cl', 'O', 'Li', 'Li', 'Li']);
        let catio3 = endpoint('phonondb-CaTiO3-mp-4019.json.gz', 'CaTiO3');
        assert.equal(catio3.family, 'perovskite-Pnma');
        assert.equal(catio3.siteTypes.length, 20);
    });

    it('aligns the octahedral tilts whatever the setting of the file', () => {
        // the same CaTiO3 with the origin on another Ca and the atoms reordered: every
        // orientation of the tilted cell matches the template, the tilts must come out
        // the same
        let raw = loadFixture('phonondb-CaTiO3-mp-4019.json.gz');
        let reference = buildMixingEndpoint(raw, { name: 'CaTiO3' }).structure;
        let shifted = loadFixture('phonondb-CaTiO3-mp-4019.json.gz');
        let ca = shifted.atom_types.map((type, i) => i).filter((i) => shifted.atom_types[i] === 'Ca');
        let shift = shifted.atom_pos_red[ca[1]].map((v, k) => v - shifted.atom_pos_red[ca[0]][k]);
        let order = shifted.atom_types.map((_, i) => i).reverse();
        shifted.atom_types = order.map((i) => raw.atom_types[i]);
        shifted.atom_pos_red = order.map((i) => raw.atom_pos_red[i].map((v, k) => v - shift[k]));
        let structure = buildMixingEndpoint(shifted, { name: 'CaTiO3' }).structure;
        let worst = 0;
        structure.positionsRed.forEach((p, s) => p.forEach((v, k) => {
            worst = Math.max(worst, Math.abs(v - reference.positionsRed[s][k]));
        }));
        assert.ok(worst < 1e-6, 'positions differ by ' + worst);
    });

    it('folds a 5-atom perovskite into the Pnma cell', async () => {
        // D of the 20-atom supercell at Q has the eigenvalues of SrTiO3 at the 4 points Q + g
        let srtio3 = endpoint('phonondb-SrTiO3-mp-5229.json.gz', 'SrTiO3');
        let folded = binaryAsSuperstructure(srtio3, 'perovskite-Pnma');
        assert.equal(folded.siteTypes.filter((type) => type === 'O').length, 12);
        let cell = srtio3.structure.lattice;
        let superReciprocal = utils.rec_lat(folded.lattice);
        let toCart = (q) => [0, 1, 2].map((d) => q[0] * superReciprocal[0][d] + q[1] * superReciprocal[1][d] + q[2] * superReciprocal[2][d]);
        let toPrimitive = (k) => cell.map((p) => p[0] * k[0] + p[1] * k[1] + p[2] * k[2]);
        for (let Q of [[0.13, -0.07, 0.21], [0.5, 0, 0]]) {
            let qCart = toCart(Q);
            let candidates = [];
            for (let a = -2; a <= 2; a++) {
                for (let b = -2; b <= 2; b++) {
                    for (let c = -2; c <= 2; c++) {
                        candidates.push(toCart([a, b, c]));
                    }
                }
            }
            candidates.sort((g, h) => g.reduce((t, v) => t + v * v, 0) - h.reduce((t, v) => t + v * v, 0));
            let seen = new Set();
            let expected = [];
            for (let g of candidates) {
                let key = toPrimitive(g).map((v) => Math.round((v - Math.floor(v + 1e-9)) * 1e6) % 1000000).join(',');
                if (seen.has(key)) {
                    continue;
                }
                seen.add(key);
                expected.push(...await eigen(srtio3.sampleAtQ(toPrimitive([0, 1, 2].map((d) => qCart[d] + g[d])), null)));
            }
            assert.equal(seen.size, 4);
            expected.sort((p, q) => p - q);
            let actual = await eigen(folded.sampleAtQ(Q, null));
            actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-8 * (1 + Math.abs(expected[i]))));
        }
    });

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        let srtio3 = endpoint('phonondb-SrTiO3-mp-5229.json.gz', 'SrTiO3');
        let kmgf3 = endpoint('phonondb-KMgF3-mp-3448.json.gz', 'KMgF3');
        let catio3 = endpoint('phonondb-CaTiO3-mp-4019.json.gz', 'CaTiO3');
        let namgh3 = endpoint('phonondb-NaMgH3-mp-23730.json.gz', 'NaMgH3');
        for (let [a, b] of [[srtio3, kmgf3], [srtio3, catio3], [catio3, namgh3]]) {
            let sampled = sampleMixingPair(a, b);
            let reversed = sampleMixingPair(b, a);
            assert.equal(sampled.recomputeNac, true);
            assert.deepEqual(sampled.labels, reversed.labels);
            let end0 = await computeMixedPhonon(sampled, 0);
            assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-4);
            let end1 = await computeMixedPhonon(sampled, 1);
            assert.ok(opticalDiff(end1.eigenvalues, await computeReferenceEigenvalues(reversed)) < 1e-4);
            let mixed = await computeMixedPhonon(sampled, 0.5, 0.3);
            for (let k of gammaIndices(sampled)) {
                assert.ok(acoustic(mixed, k) < 0.05, `${a.name} + ${b.name}: ${acoustic(mixed, k)}`);
            }
            let mirrored = await computeMixedPhonon(reversed, 0.5, 0.7);
            assert.ok(opticalDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-6);
        }
    });

    it('mixes a perovskite with an antiperovskite site by site', () => {
        // the same five sites with cations and anions swapped (see mixing-phonons-3)
        let srtio3 = endpoint('phonondb-SrTiO3-mp-5229.json.gz', 'SrTiO3');
        let li3clo = endpoint('phonondb-Li3ClO-mp-985585.json.gz', 'Li3ClO');
        let catio3 = endpoint('phonondb-CaTiO3-mp-4019.json.gz', 'CaTiO3');
        assert.equal(sampleMixingPair(srtio3, li3clo).union.map1.length, 5);
        // the 5-atom antiperovskite folded into the 20-atom tilted cell
        assert.equal(sampleMixingPair(li3clo, catio3).union.map1.length, 20);
    });
});

describe('Materials Project path alignment', () => {
    it('lines up two MgAgAs type antifluorite files', () => {
        // the two files put L at different, equivalent points; the reference points of
        // the template choose matches in which they coincide
        let limgn = buildMixingEndpoint(loadFixture('mpdb-LiMgN-mp-37906-trimmed.json.gz'), { name: 'LiMgN' });
        let liznas = buildMixingEndpoint(loadFixture('mpdb-LiZnAs-mp-9124-trimmed.json.gz'), { name: 'LiZnAs' });
        let points = (e) => e.segments.reduce((total, segment) => total + segment.npoints, 0);
        let sampled = sampleMixingPair(limgn, liznas);
        assert.ok(sampled.qpoints.length >= 0.9 * Math.min(points(limgn), points(liznas)));
    });
});

describe('transition-metal dichalcogenide mixing', () => {
    let gammaIndices = (sampled) => sampled.qpoints
        .map((q, k) => (q.every((v) => Math.abs(v - Math.round(v)) < 1e-9) ? k : -1))
        .filter((k) => k >= 0);
    let opticalDiff = (a, b) => {
        let worst = 0;
        a.forEach((row, k) => row.forEach((v, n) => {
            if (Math.abs(v) > 1) {
                worst = Math.max(worst, Math.abs(v - b[k][n]));
            }
        }));
        return worst;
    };

    it('matches the 2H, 1T and 3R polytypes', () => {
        let ws2 = endpoint('phonondb-WS2-2H-mp-224.json.gz', 'WS2');
        assert.equal(ws2.family, 'tmd-2H');
        assert.deepEqual(ws2.siteTypes, ['W', 'W', 'S', 'S', 'S', 'S']);
        // the hexagonal cell along the canonical axes, the interlayer distance kept
        let lattice = ws2.structure.lattice;
        assert.ok(Math.abs(lattice[1][0] / lattice[0][0] + 0.5) < 1e-6 && Math.abs(lattice[2][0]) < 1e-6);
        assert.ok(lattice[2][2] > 3 * lattice[0][0]);
        // 1T: the chalcogens above and below the metal in either setting
        let hfs2 = endpoint('phonondb-HfS2-1T-mp-985829.json.gz', 'HfS2');
        let zrs2 = endpoint('phonondb-ZrS2-1T-mp-1186.json.gz', 'ZrS2');
        assert.equal(hfs2.family, 'tmd-1T');
        assert.equal(zrs2.family, 'tmd-1T');
        assert.deepEqual(hfs2.siteTypes, ['Hf', 'S', 'S']);
        // 3R: 3 atoms like 1T, but rhombohedral
        let mos2 = endpoint('phonondb-MoS2-3R-mp-1434.json.gz', 'MoS2');
        assert.equal(mos2.family, 'tmd-3R');
        assert.deepEqual(mos2.siteTypes, ['Mo', 'S', 'S']);
    });

    it('maps a Materials Project 3R cell (other axes, polar) onto the PhononDB one', () => {
        let pdb = endpoint('phonondb-MoSe2-3R-mp-7581.json.gz', 'MoSe2');
        let mp = buildMixingEndpoint(loadFixture('mpdb-MoSe2-3R-mp-7581-trimmed.json.gz'), { name: 'MoSe2' });
        assert.equal(mp.family, 'tmd-3R');
        let relative = (a, b) => {
            let n = 0;
            let d = 0;
            a.real.forEach((row, i) => row.forEach((v, j) => {
                n += (v - b.real[i][j]) ** 2 + (a.imag[i][j] - b.imag[i][j]) ** 2;
                d += v * v + a.imag[i][j] ** 2;
            }));
            return Math.sqrt(n / d);
        };
        let median = (sampled) => {
            let r = sampled.matrices1.map((m, k) => relative(m, sampled.matrices2[k])).sort((p, q) => p - q);
            return r[r.length >> 1];
        };
        assert.ok(median(sampleMixingPair(pdb, mp)) < 0.02);
        mp.invert = !mp.invert;
        assert.ok(median(sampleMixingPair(pdb, mp)) > 0.1);
    });

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        let ws2 = endpoint('phonondb-WS2-2H-mp-224.json.gz', 'WS2');
        let wse2 = endpoint('phonondb-WSe2-2H-mp-1821.json.gz', 'WSe2');
        let hfs2 = endpoint('phonondb-HfS2-1T-mp-985829.json.gz', 'HfS2');
        let zrs2 = endpoint('phonondb-ZrS2-1T-mp-1186.json.gz', 'ZrS2');
        let mos2 = endpoint('phonondb-MoS2-3R-mp-1434.json.gz', 'MoS2');
        let mose2 = endpoint('phonondb-MoSe2-3R-mp-7581.json.gz', 'MoSe2');
        for (let [a, b] of [[ws2, wse2], [hfs2, zrs2], [mos2, mose2]]) {
            let sampled = sampleMixingPair(a, b);
            let reversed = sampleMixingPair(b, a);
            assert.equal(sampled.recomputeNac, true);
            assert.deepEqual(sampled.labels, reversed.labels);
            let end0 = await computeMixedPhonon(sampled, 0);
            assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-3);
            let end1 = await computeMixedPhonon(sampled, 1);
            assert.ok(opticalDiff(end1.eigenvalues, await computeReferenceEigenvalues(reversed)) < 1e-3);
            let mixed = await computeMixedPhonon(sampled, 0.5, 0.3);
            for (let k of gammaIndices(sampled)) {
                for (let n = 0; n < 3; n++) {
                    assert.ok(Math.abs(mixed.eigenvalues[k][n]) < 0.05, `${a.name} + ${b.name}: ${mixed.eigenvalues[k][n]}`);
                }
            }
            let mirrored = await computeMixedPhonon(reversed, 0.5, 0.7);
            assert.ok(opticalDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-6);
        }
    });

    it('mixes different polytypes on the close-packed grid', () => {
        // the layers of 2H, 1T and 3R on one grid: a common cell with vacancies
        let ws2 = endpoint('phonondb-WS2-2H-mp-224.json.gz', 'WS2');
        let hfs2 = endpoint('phonondb-HfS2-1T-mp-985829.json.gz', 'HfS2');
        let mos2 = endpoint('phonondb-MoS2-3R-mp-1434.json.gz', 'MoS2');
        assert.equal(sampleMixingPair(ws2, mos2).union.world, 'hex');
        assert.equal(sampleMixingPair(hfs2, mos2).union.world, 'hex');
        assert.equal(sampleMixingPair(ws2, hfs2).union.map1.length, 10);
    });
});

describe('mixing physics: dipole-dipole split and gauge', () => {
    // D of an end-member in its canonical frame (no rotation: site order and inversion only)
    let canonicalD = (e, q) => {
        let file = e.evaluate(q, null);
        let n = e.order.length;
        let out = { real: [], imag: [] };
        for (let s = 0; s < n; s++) {
            for (let a = 0; a < 3; a++) {
                let re = [];
                let im = [];
                for (let t = 0; t < n; t++) {
                    for (let b = 0; b < 3; b++) {
                        re.push(file.real[e.order[s] * 3 + a][e.order[t] * 3 + b]);
                        im.push((e.invert ? -1 : 1) * file.imag[e.order[s] * 3 + a][e.order[t] * 3 + b]);
                    }
                }
                out.real.push(re);
                out.imag.push(im);
            }
        }
        return out;
    };
    let addDipole = (matrix, payload, q, factor) => {
        let out = { real: matrix.real.map((r) => r.slice()), imag: matrix.imag.map((r) => r.slice()) };
        let dd = getGonzeDipoleDipole(payload, q, null);
        dd.forEach((row, i) => row.forEach((block, j) => block.forEach((line, a) => line.forEach((v, b) => {
            out.real[3 * i + a][3 * j + b] += factor * v[0];
            out.imag[3 * i + a][3 * j + b] += factor * v[1];
        }))));
        return out;
    };
    let ownPayload = (e, lattice, positionsRed) => buildNacPayload(e.nac.siteBorn, e.nac.siteDielectric, e.nac.unitConversion,
        lattice, e.siteMasses, utils.red_car_list(positionsRed, lattice), 1e-5);
    // a "sampled" pair at chosen q-points (no Gamma), for computeMixedPhonon
    let customPair = (e1, e2, qpoints, sample1, sample2, payload1, payload2) => {
        let f = e1.frequencyFactor * e1.frequencyFactor;
        let m1 = qpoints.map(sample1);
        let m2 = qpoints.map(sample2);
        return {
            endpoint1: e1,
            endpoint2: e2,
            qpoints: qpoints,
            matrices1: m1,
            matrices2: m2,
            recomputeNac: true,
            shortRange1: m1.map((m, k) => addDipole(m, payload1, qpoints[k], -f)),
            shortRange2: m2.map((m, k) => addDipole(m, payload2, qpoints[k], -f)),
            qDirections: qpoints.map(() => null),
            lineBreaks: [[0, qpoints.length]],
            labels: [['', '']],
        };
    };
    let foldingVectors = (superLattice, primitiveLattice) => {
        let reciprocal = utils.rec_lat(superLattice);
        let candidates = [];
        for (let a = -2; a <= 2; a++) {
            for (let b = -2; b <= 2; b++) {
                for (let c = -2; c <= 2; c++) {
                    candidates.push([a, b, c]);
                }
            }
        }
        // the shortest representative of each class, as foldIntoSupercell takes (the
        // truncated reciprocal sum makes D(k) periodic only to ~1e-5)
        let vectors = candidates.map((n) => [0, 1, 2].map((d) => n[0] * reciprocal[0][d] + n[1] * reciprocal[1][d] + n[2] * reciprocal[2][d]));
        vectors.sort((p, q) => (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]));
        let seen = new Set();
        let out = [];
        for (let g of vectors) {
            let k = primitiveLattice.map((p) => p[0] * g[0] + p[1] * g[1] + p[2] * g[2]);
            let key = k.map((v) => Math.round((v - Math.floor(v + 1e-9)) * 1e6) % 1000000).join(',');
            if (!seen.has(key)) {
                seen.add(key);
                out.push(g);
            }
        }
        return out;
    };

    it('uses a full dipole-dipole term that does not depend on the Ewald parameter', () => {
        for (let [file, name] of [['phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN'], ['phonondb-CuGaS2-mp-5238.json.gz', 'CuGaS2']]) {
            let e = endpoint(file, name);
            let positions = utils.red_car_list(e.positionsRed, e.lattice);
            let payload = (lambda) => buildNacPayload(e.nac.born, e.nac.dielectric, e.nac.unitConversion, e.lattice, e.masses, positions, 1e-5, lambda);
            let lambda = payload(null).nac.lambda;
            for (let q of [[0.11, -0.07, 0.23], [0.5, 0, 0]]) {
                let a = getGonzeDipoleDipole(payload(lambda), q).flat(3);
                let b = getGonzeDipoleDipole(payload(1.6 * lambda), q).flat(3);
                let reciprocalOnly = getGonzeReciprocalCorrection(payload(1.6 * lambda), q).flat(3);
                let scale = Math.max(...a.map((v) => Math.hypot(v[0], v[1])));
                let diff = Math.max(...a.map((v, k) => Math.hypot(v[0] - b[k][0], v[1] - b[k][1])));
                let diffReciprocal = Math.max(...a.map((v, k) => Math.hypot(v[0] - reciprocalOnly[k][0], v[1] - reciprocalOnly[k][1])));
                assert.ok(diff < 1e-6 * scale, name + ': ' + diff / scale);
                // phonopy's reciprocal part alone does depend on it
                assert.ok(diffReciprocal > 1e-2 * scale);
            }
        }
    });

    it('splits off the same short-range part in the primitive cell and in a supercell', () => {
        for (let [file, name, family] of [['phonondb-NaCl-mp-22862.json.gz', 'NaCl', 'rocksalt-ABX2-R3m'], ['phonondb-CaF2-mp-2741.json.gz', 'CaF2', 'fluorite-ABX4-I41a']]) {
            let e = endpoint(file, name);
            let f = e.frequencyFactor * e.frequencyFactor;
            let positions = e.family === 'rocksalt' ? [[0, 0, 0], [0.5, 0.5, 0.5]] : [[0, 0, 0], [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]];
            let primitivePayload = ownPayload(e, e.lattice, positions);
            let folded = binaryAsSuperstructure(e, family);
            let superPayload = ownPayload(folded, folded.lattice, folded.structure.positionsRed);
            // short range of the primitive cell, folded with the same machinery
            let primitiveShortRange = Object.assign({}, e, {
                order: e.order.map((_, s) => s),
                invert: false,
                evaluate: (q) => addDipole(canonicalD(e, q), primitivePayload, q, -f),
            });
            let foldedShortRange = binaryAsSuperstructure(primitiveShortRange, family);
            for (let Q of [[0.1, 0.05, -0.07], [0.5, 0, 0], [0, 0.5, 0.5]]) {
                let direct = addDipole(folded.sampleAtQ(Q, null), superPayload, Q, -f);
                let viaPrimitive = foldedShortRange.sampleAtQ(Q, null);
                let scale = 0;
                let diff = 0;
                direct.real.forEach((row, i) => row.forEach((v, j) => {
                    scale = Math.max(scale, Math.hypot(v, direct.imag[i][j]));
                    diff = Math.max(diff, Math.hypot(v - viaPrimitive.real[i][j], direct.imag[i][j] - viaPrimitive.imag[i][j]));
                }));
                assert.ok(diff < 1e-5 * scale, name + ': ' + diff / scale);
            }
        }
    });

    it('commutes mixing with folding', async () => {
        // mixing two binaries and folding the mixture = folding both and mixing in the
        // supercell; with GaN / AlN the internal parameter u differs (gauge)
        let x = 0.4;
        for (let [f1, n1, f2, n2, family] of [
            ['phonondb-NaCl-mp-22862.json.gz', 'NaCl', 'phonondb-KCl-mp-23193.json.gz', 'KCl', 'rocksalt-ABX2-R3m'],
            ['phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN', 'phonondb-AlN-wurtzite-mp-661.json.gz', 'AlN', 'wurtzite-ABX2'],
        ]) {
            let e1 = endpoint(f1, n1);
            let e2 = endpoint(f2, n2);
            let s1 = binaryAsSuperstructure(e1, family);
            let s2 = binaryAsSuperstructure(e2, family);
            let Qs = [[0.13, -0.07, 0.21], [0.31, 0.12, -0.18]];
            let superPair = customPair(s1, s2, Qs, (q) => s1.sampleAtQ(q, null), (q) => s2.sampleAtQ(q, null),
                ownPayload(s1, s1.structure.lattice, s1.structure.positionsRed),
                ownPayload(s2, s2.structure.lattice, s2.structure.positionsRed));
            let superMixed = await computeMixedPhonon(superPair, x);

            let positions = (e) => e.family === 'rocksalt' ? [[0, 0, 0], [0.5, 0.5, 0.5]] : wurtzitePositionsOf(e.structure.u);
            let latticeOf = (e) => e.family === 'rocksalt' ? e.lattice : e.lattice;
            let gs = foldingVectors(s1.structure.lattice, latticeOf(e1));
            let superReciprocal = utils.rec_lat(s1.structure.lattice);
            for (let k = 0; k < Qs.length; k++) {
                let Q = Qs[k];
                let qCart = [0, 1, 2].map((d) => Q[0] * superReciprocal[0][d] + Q[1] * superReciprocal[1][d] + Q[2] * superReciprocal[2][d]);
                let ks = gs.map((g) => latticeOf(e1).map((p) => p[0] * (qCart[0] + g[0]) + p[1] * (qCart[1] + g[1]) + p[2] * (qCart[2] + g[2])));
                let primitivePair = customPair(e1, e2, ks, (q) => canonicalD(e1, q), (q) => canonicalD(e2, q),
                    ownPayload(e1, latticeOf(e1), positions(e1)), ownPayload(e2, latticeOf(e2), positions(e2)));
                let primitiveMixed = await computeMixedPhonon(primitivePair, x);
                let expected = primitiveMixed.eigenvalues.flat().sort((p, q) => p - q);
                let actual = superMixed.eigenvalues[k].slice().sort((p, q) => p - q);
                assert.equal(actual.length, expected.length);
                actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-3, `${n1} + ${n2}: ${v} vs ${expected[i]}`));
            }
        }
    });
});

describe('mixing physics: periodicity', () => {
    it('gives frequencies periodic in the reciprocal lattice when the internal positions differ', async () => {
        // force constants Phi(R) mixed (not atom-gauge D's of different positions): the
        // mixed crystal is a crystal, omega(q + G) = omega(q). GaN and AlN differ in u;
        // mixing the atom-gauge D's directly breaks this by ~0.1 cm-1 here.
        let canonical = (e, q) => {
            let file = e.evaluate(q, null);
            let n = e.order.length;
            let out = { real: [], imag: [] };
            for (let s = 0; s < n; s++) {
                for (let a = 0; a < 3; a++) {
                    let re = [];
                    let im = [];
                    for (let t = 0; t < n; t++) {
                        for (let b = 0; b < 3; b++) {
                            re.push(file.real[e.order[s] * 3 + a][e.order[t] * 3 + b]);
                            im.push((e.invert ? -1 : 1) * file.imag[e.order[s] * 3 + a][e.order[t] * 3 + b]);
                        }
                    }
                    out.real.push(re);
                    out.imag.push(im);
                }
            }
            return out;
        };
        let gan = endpoint('phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN');
        let aln = endpoint('phonondb-AlN-wurtzite-mp-661.json.gz', 'AlN');
        // shifts along c, where the positions differ (u); the files themselves are
        // periodic there to ~1e-5 cm-1 (farther out, phonopy's truncated G sum is not)
        let q = [0.17, 0.08, 0.21];
        let qpoints = [q, [q[0], q[1], q[2] + 1], [q[0], q[1], q[2] - 1]];
        let pair = {
            endpoint1: gan,
            endpoint2: aln,
            qpoints: qpoints,
            matrices1: qpoints.map((k) => canonical(gan, k)),
            matrices2: qpoints.map((k) => canonical(aln, k)),
            recomputeNac: false,
            lineBreaks: [[0, qpoints.length]],
            labels: [['', '']],
        };
        let mixed = await computeMixedPhonon(pair, 0.5);
        for (let k = 1; k < qpoints.length; k++) {
            mixed.eigenvalues[k].forEach((v, n) => assert.ok(Math.abs(v - mixed.eigenvalues[0][n]) < 1e-3,
                `q + G ${k}: ${v} vs ${mixed.eigenvalues[0][n]}`));
        }
    });
});

function wurtzitePositionsOf(u) {
    // canonical wurtzite positions, as in mixingphonons.wurtzitePositions
    return [[1 / 3, 2 / 3, 0], [2 / 3, 1 / 3, 0.5], [1 / 3, 2 / 3, u], [2 / 3, 1 / 3, 0.5 + u]]
        .map((p) => p.map((v) => v - Math.floor(v + 1e-9)));
}

describe('mixing: presentation in the setting of material 1', () => {
    let fileStructure = (raw) => raw.structure
        ? { lattice: raw.structure.lattice.matrix, positions: raw.structure.sites.map((s) => s.abc), types: raw.structure.sites.map((s) => s.label) }
        : { lattice: raw.lattice, positions: raw.atom_pos_red, types: raw.atom_types };

    it('draws the mixture at x = 0 exactly as material 1 on its own', async () => {
        // inverted wurtzite (GaN), tilted perovskite (CaTiO3), and a Materials Project 3R
        // cell stored with c along x: the mixture is computed in a canonical cell (other
        // origin, atom order, axes, or the inversion image) but drawn in the file's setting
        let cases = [
            ['phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN', 'phonondb-AlN-wurtzite-mp-661.json.gz', 'AlN'],
            ['phonondb-CaTiO3-mp-4019.json.gz', 'CaTiO3', 'phonondb-NaMgH3-mp-23730.json.gz', 'NaMgH3'],
            ['mpdb-MoSe2-3R-mp-7581-trimmed.json.gz', 'MoSe2', 'phonondb-MoS2-3R-mp-1434.json.gz', 'MoS2'],
        ];
        for (let [f1, n1, f2, n2] of cases) {
            let raw = loadFixture(f1);
            let e1 = buildMixingEndpoint(raw, { name: n1 });
            let sampled = sampleMixingPair(e1, endpoint(f2, n2));
            let mixed = await computeMixedPhonon(sampled, 0);
            let file = fileStructure(raw);
            assert.deepEqual(mixed.atom_types, file.types, n1);
            mixed.lattice.forEach((row, i) => row.forEach((v, j) => assert.ok(Math.abs(v - file.lattice[i][j]) < 1e-6, n1)));
            mixed.atom_pos_red.forEach((p, i) => p.forEach((v, k) => {
                let d = v - file.positions[i][k];
                assert.ok(Math.abs(d - Math.round(d)) < 1e-6, n1 + ' atom ' + i);
            }));
            assert.deepEqual(mixed.repetitions, raw.repetitions || [3, 3, 3]);
            // the modes are the file's own: same eigenvectors at the same (file) q
            if (e1.evaluate) {
                let k = Math.floor(mixed.qpoints.length * 0.37);
                let D = e1.evaluate(mixed.qpoints[k], null);
                let solution = await solveComplexHermitianWithEigenWasm(D.real, D.imag);
                let order = solution.values.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]);
                let n = mixed.natoms;
                for (let mode = 0; mode < 3 * n; mode++) {
                    let [value, index] = order[mode];
                    let isolated = [mode - 1, mode + 1].every((other) => other < 0 || other >= 3 * n ||
                        Math.abs(order[other][0] - value) > 1e-3 * Math.abs(value));
                    if (!isolated) {
                        continue;
                    }
                    let re = 0;
                    let im = 0;
                    let normFile = 0;
                    let normMixed = 0;
                    for (let a = 0; a < n; a++) {
                        for (let c = 0; c < 3; c++) {
                            let f = solution.vectors[index][3 * a + c];
                            let g = mixed.vectors[k][mode][a][c].map((v) => v * Math.sqrt(e1.masses[a]));
                            re += f[0] * g[0] + f[1] * g[1];
                            im += f[0] * g[1] - f[1] * g[0];
                            normFile += f[0] * f[0] + f[1] * f[1];
                            normMixed += g[0] * g[0] + g[1] * g[1];
                        }
                    }
                    assert.ok(Math.abs(Math.hypot(re, im) / Math.sqrt(normFile * normMixed) - 1) < 1e-8, n1 + ' mode ' + mode);
                }
            }
        }
    });

    it('draws a binary folded into the ternary cell in the binary\'s own setting', async () => {
        // material 1 a binary, material 2 an ordered ternary: at x = 0 the supercell holds
        // the binary, every atom on an atom of its file, the cell vectors lattice vectors
        // of the file (inverted GaN, rotated SrTiO3), and a fold into a cell of the same
        // size (Li2O in the MgAgAs type cell) is the file itself
        let inverse = (M) => {
            let d = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
            return [
                [(M[1][1] * M[2][2] - M[1][2] * M[2][1]) / d, (M[0][2] * M[2][1] - M[0][1] * M[2][2]) / d, (M[0][1] * M[1][2] - M[0][2] * M[1][1]) / d],
                [(M[1][2] * M[2][0] - M[1][0] * M[2][2]) / d, (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / d, (M[0][2] * M[1][0] - M[0][0] * M[1][2]) / d],
                [(M[1][0] * M[2][1] - M[1][1] * M[2][0]) / d, (M[0][1] * M[2][0] - M[0][0] * M[2][1]) / d, (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / d],
            ];
        };
        let row = (v, M) => [0, 1, 2].map((j) => v[0] * M[0][j] + v[1] * M[1][j] + v[2] * M[2][j]);
        let frac = (v) => Math.abs(v - Math.round(v));
        let cases = [
            ['phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN', 'phonondb-ZnGeN2-mp-2979.json.gz', 'ZnGeN2'],
            ['phonondb-NaCl-mp-22862.json.gz', 'NaCl', 'phonondb-NaInS2-mp-20289.json.gz', 'NaInS2'],
            ['phonondb-SrTiO3-mp-5229.json.gz', 'SrTiO3', 'phonondb-CaTiO3-mp-4019.json.gz', 'CaTiO3'],
            ['phonondb-Li2O-mp-1960.json.gz', 'Li2O', 'phonondb-LiMgN-mp-37906.json.gz', 'LiMgN'],
            ['phonondb-TiO2-mp-2657.json.gz', 'TiO2', 'phonondb-ZnSb2O6-mp-3188.json.gz', 'ZnSb2O6'],
        ];
        for (let [f1, n1, f2, n2] of cases) {
            let raw = loadFixture(f1);
            let mixed = await computeMixedPhonon(sampleMixingPair(endpoint(f1, n1), endpoint(f2, n2)), 0);
            let fileInverse = inverse(raw.lattice);
            utils.red_car_list(mixed.atom_pos_red, mixed.lattice).forEach((cartesian, i) => {
                let reduced = row(cartesian, fileInverse);
                let onSite = raw.atom_pos_red.some((p, j) => raw.atom_types[j] === mixed.atom_types[i] &&
                    reduced.every((v, k) => frac(v - p[k]) < 1e-6));
                assert.ok(onSite, `${n1}: atom ${i} (${mixed.atom_types[i]}) is not on a site of the file`);
            });
            mixed.lattice.forEach((vector) => row(vector, fileInverse).forEach((v) => assert.ok(frac(v) < 1e-6, n1)));
            if (mixed.natoms === raw.natoms) {
                assert.deepEqual(mixed.atom_types, raw.atom_types);
                mixed.lattice.forEach((r, i) => r.forEach((v, j) => assert.ok(Math.abs(v - raw.lattice[i][j]) < 1e-6)));
            }
        }
    });
});

describe('rutile mixing', () => {
    let gammaIndices = (sampled) => sampled.qpoints
        .map((q, k) => (q.every((v) => Math.abs(v - Math.round(v)) < 1e-9) ? k : -1))
        .filter((k) => k >= 0);
    let opticalDiff = (a, b) => {
        let worst = 0;
        a.forEach((row, k) => row.forEach((v, n) => {
            if (Math.abs(v) > 1) {
                worst = Math.max(worst, Math.abs(v - b[k][n]));
            }
        }));
        return worst;
    };
    let eigen = async (m) => (await solveComplexHermitianWithEigenWasm(m.real, m.imag)).values.sort((p, q) => p - q);

    it('matches rutile and its CaCl2-type distortion, not marcasite', () => {
        let tio2 = endpoint('phonondb-TiO2-mp-2657.json.gz', 'TiO2');
        assert.equal(tio2.family, 'rutile');
        assert.deepEqual(tio2.siteTypes, ['Ti', 'Ti', 'O', 'O', 'O', 'O']);
        let cacl2 = endpoint('phonondb-CaCl2-Pnnm-mp-23214.json.gz', 'CaCl2');
        assert.equal(cacl2.family, 'rutile');
        // marcasite OsAs2 (As-As pairs): the anions are 0.15 of the cell off the rutile sites
        assert.throws(() => endpoint('phonondb-OsAs2-mp-2455.json.gz', 'OsAs2'));
        let znsb2o6 = endpoint('phonondb-ZnSb2O6-mp-3188.json.gz', 'ZnSb2O6');
        assert.equal(znsb2o6.family, 'rutile-AB2X6');
        assert.deepEqual(znsb2o6.siteTypes.slice(0, 6), ['Sb', 'Sb', 'Sb', 'Sb', 'Zn', 'Zn']);
    });

    it('aligns the CaCl2-type distortion whatever the setting of the file', () => {
        // the orthorhombic distortion fits the tetragonal template in two senses, and
        // which one a file lands in depends on its cell basis (the PhononDB and Materials
        // Project files of CaCl2 differed); the reference distortion of the template picks
        // the same one. Here the same crystal with the basis relabelled a, b, c -> b, c, a
        // and the atoms in reverse order (without the reference: 0.06 apart)
        let raw = loadFixture('phonondb-CaCl2-Pnnm-mp-23214.json.gz');
        let reference = buildMixingEndpoint(raw, { name: 'CaCl2' }).structure;
        let relabelled = loadFixture('phonondb-CaCl2-Pnnm-mp-23214.json.gz');
        relabelled.lattice = [raw.lattice[1], raw.lattice[2], raw.lattice[0]];
        relabelled.atom_pos_red = raw.atom_pos_red.map((p) => [p[1], p[2], p[0]]).reverse();
        relabelled.atom_types = raw.atom_types.slice().reverse();
        let structure = buildMixingEndpoint(relabelled, { name: 'CaCl2' }).structure;
        structure.positionsRed.forEach((p, s) => p.forEach((v, k) => assert.ok(Math.abs(v - reference.positionsRed[s][k]) < 1e-6)));
    });

    it('folds a rutile into the trirutile cell', async () => {
        // D of the 18-atom trirutile cell at Q has the eigenvalues of TiO2 at the 3 points Q + g
        let tio2 = endpoint('phonondb-TiO2-mp-2657.json.gz', 'TiO2');
        let folded = binaryAsSuperstructure(tio2, 'rutile-AB2X6');
        let cell = tio2.structure.lattice;
        let reciprocal = utils.rec_lat(folded.structure.lattice);
        let Q = [0.13, -0.07, 0.21];
        let qCart = [0, 1, 2].map((d) => Q[0] * reciprocal[0][d] + Q[1] * reciprocal[1][d] + Q[2] * reciprocal[2][d]);
        let candidates = [];
        for (let a = -2; a <= 2; a++) {
            for (let b = -2; b <= 2; b++) {
                for (let c = -2; c <= 2; c++) {
                    candidates.push([0, 1, 2].map((d) => a * reciprocal[0][d] + b * reciprocal[1][d] + c * reciprocal[2][d]));
                }
            }
        }
        candidates.sort((g, h) => g.reduce((t, v) => t + v * v, 0) - h.reduce((t, v) => t + v * v, 0));
        let seen = new Set();
        let expected = [];
        for (let g of candidates) {
            let key = cell.map((p) => {
                let t = p[0] * g[0] + p[1] * g[1] + p[2] * g[2];
                return Math.round((t - Math.floor(t + 1e-9)) * 1e6) % 1000000;
            }).join(',');
            if (seen.has(key)) {
                continue;
            }
            seen.add(key);
            let k = cell.map((p) => p[0] * (qCart[0] + g[0]) + p[1] * (qCart[1] + g[1]) + p[2] * (qCart[2] + g[2]));
            expected.push(...await eigen(tio2.sampleAtQ(k, null)));
        }
        assert.equal(seen.size, 3);
        expected.sort((p, q) => p - q);
        let actual = await eigen(folded.sampleAtQ(Q, null));
        actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-8 * (1 + Math.abs(expected[i]))));
    });

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        let tio2 = endpoint('phonondb-TiO2-mp-2657.json.gz', 'TiO2');
        let sno2 = endpoint('phonondb-SnO2-mp-856.json.gz', 'SnO2');
        let znsb2o6 = endpoint('phonondb-ZnSb2O6-mp-3188.json.gz', 'ZnSb2O6');
        for (let [a, b] of [[sno2, tio2], [tio2, znsb2o6]]) {
            let sampled = sampleMixingPair(a, b);
            let reversed = sampleMixingPair(b, a);
            assert.equal(sampled.recomputeNac, true);
            let end0 = await computeMixedPhonon(sampled, 0);
            assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-3);
            let end1 = await computeMixedPhonon(sampled, 1);
            assert.ok(opticalDiff(end1.eigenvalues, await computeReferenceEigenvalues(reversed)) < 1e-3);
            let mixed = await computeMixedPhonon(sampled, 0.5, 0.3);
            for (let k of gammaIndices(sampled)) {
                let near = mixed.eigenvalues[k].map(Math.abs).sort((p, q) => p - q);
                assert.ok(near[2] < 0.05, `${a.name} + ${b.name}: ${near[2]}`);
            }
            let mirrored = await computeMixedPhonon(reversed, 0.5, 0.7);
            assert.ok(opticalDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-6);
        }
    });
});

describe('corundum mixing', () => {
    let gammaIndices = (sampled) => sampled.qpoints
        .map((q, k) => (q.every((v) => Math.abs(v - Math.round(v)) < 1e-9) ? k : -1))
        .filter((k) => k >= 0);
    let opticalDiff = (a, b) => {
        let worst = 0;
        a.forEach((row, k) => row.forEach((v, n) => {
            if (Math.abs(v) > 1) {
                worst = Math.max(worst, Math.abs(v - b[k][n]));
            }
        }));
        return worst;
    };
    let eigen = async (m) => (await solveComplexHermitianWithEigenWasm(m.real, m.imag)).values.sort((p, q) => p - q);

    it('matches corundum and ilmenite, not calcite', () => {
        let al2o3 = endpoint('phonondb-Al2O3-mp-1143.json.gz', 'Al2O3');
        assert.equal(al2o3.family, 'corundum');
        assert.deepEqual(al2o3.siteTypes, ['Al', 'Al', 'Al', 'Al', 'O', 'O', 'O', 'O', 'O', 'O']);
        // MgTiO3 and CdTiO3 (the face-sharing pair longer than the in-plane distance):
        // ilmenite, the A and B layers alternating
        for (let [file, name, a, b] of [['phonondb-MgTiO3-mp-3771.json.gz', 'MgTiO3', 'Mg', 'Ti'], ['phonondb-CdTiO3-mp-13641.json.gz', 'CdTiO3', 'Cd', 'Ti']]) {
            let e = endpoint(file, name);
            assert.equal(e.family, 'corundum-ABX3', name);
            assert.deepEqual(e.siteTypes.slice(0, 4), [a, a, b, b], name);
            // ilmenite: across every shared face an A and a B (the nearest cation of site 0 is B)
            let cartesian = utils.red_car_list(e.structure.positionsRed, e.structure.lattice);
            let lattice = e.structure.lattice;
            let nearest = null;
            for (let j = 1; j < 4; j++) {
                for (let n0 = -1; n0 <= 1; n0++) {
                    for (let n1 = -1; n1 <= 1; n1++) {
                        for (let n2 = -1; n2 <= 1; n2++) {
                            let d = [0, 1, 2].map((k) => cartesian[j][k] + n0 * lattice[0][k] + n1 * lattice[1][k] + n2 * lattice[2][k] - cartesian[0][k]);
                            let r = Math.hypot(...d);
                            // along the 3-fold axis: the face-sharing partner
                            let along = Math.abs(d[0]) < 1e-3 && Math.abs(d[1]) < 1e-3;
                            if (along && (!nearest || r < nearest.r)) {
                                nearest = { r: r, type: e.siteTypes[j] };
                            }
                        }
                    }
                }
            }
            assert.equal(nearest.type, b, name);
        }
        // calcite has the corundum site pattern but planar CO3 groups, not octahedra
        let calcite = JSON.parse(readFileSync(join(__dirname, 'fixtures', 'mixing', 'structure-CaCO3-calcite-mp-3953.json')));
        assert.throws(() => analyzeSuperstructure('corundum', calcite.atom_types, calcite.atom_pos_red, calcite.lattice), /octahedra/);
    });

    it('folds a corundum into the ilmenite cell of the same size', async () => {
        let al2o3 = endpoint('phonondb-Al2O3-mp-1143.json.gz', 'Al2O3');
        let folded = binaryAsSuperstructure(al2o3, 'corundum-ABX3');
        for (let q of [[0.13, -0.07, 0.21], [0.5, 0, 0]]) {
            let expected = await eigen(al2o3.sampleAtQ(q, null));
            let actual = await eigen(folded.sampleAtQ(q, null));
            actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-8 * (1 + Math.abs(v))));
        }
    });

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        let al2o3 = endpoint('phonondb-Al2O3-mp-1143.json.gz', 'Al2O3');
        let ga2o3 = endpoint('phonondb-Ga2O3-mp-1243.json.gz', 'Ga2O3');
        let mgtio3 = endpoint('phonondb-MgTiO3-mp-3771.json.gz', 'MgTiO3');
        let cdtio3 = endpoint('phonondb-CdTiO3-mp-13641.json.gz', 'CdTiO3');
        for (let [a, b] of [[al2o3, ga2o3], [al2o3, mgtio3], [mgtio3, cdtio3]]) {
            let sampled = sampleMixingPair(a, b);
            let reversed = sampleMixingPair(b, a);
            assert.equal(sampled.recomputeNac, true);
            let end0 = await computeMixedPhonon(sampled, 0);
            assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-3);
            let end1 = await computeMixedPhonon(sampled, 1);
            assert.ok(opticalDiff(end1.eigenvalues, await computeReferenceEigenvalues(reversed)) < 1e-3);
            let mixed = await computeMixedPhonon(sampled, 0.5, 0.3);
            for (let k of gammaIndices(sampled)) {
                let near = mixed.eigenvalues[k].map(Math.abs).sort((p, q) => p - q);
                assert.ok(near[2] < 0.05, `${a.name} + ${b.name}: ${near[2]}`);
            }
            let mirrored = await computeMixedPhonon(reversed, 0.5, 0.7);
            assert.ok(opticalDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-6);
        }
    });
});

describe('Materials Project paths related by time reversal', () => {
    it('mixes two files whose paths run in opposite halves (q -> -q)', async () => {
        // a Materials Project file and its own image under q -> -q (D(-q) = D(q)^*): before,
        // they shared no segment; now they mix, and give the material's own frequencies
        let raw = loadFixture('mpdb-NaInS2-mp-20289-trimmed.json.gz');
        let e1 = buildMixingEndpoint(raw, { name: 'NaInS2' });
        let e2 = buildMixingEndpoint(loadFixture('mpdb-NaInS2-mp-20289-trimmed.json.gz'), { name: 'NaInS2' });
        e2.segments = e2.segments.map((segment) => Object.assign({}, segment, {
            qA: segment.qA.map((v) => -v),
            qB: segment.qB.map((v) => -v),
            sample: (t) => {
                let matrix = segment.sample(t);
                return { real: matrix.real, imag: matrix.imag.map((row) => row.map((v) => -v)) };
            },
        }));
        let sampled = sampleMixingPair(e1, e2);
        let points = e1.segments.reduce((total, segment) => total + segment.npoints, 0);
        assert.ok(sampled.qpoints.length >= 0.9 * points);
        let mixed = await computeMixedPhonon(sampled, 0.5);
        let own = await computeReferenceEigenvalues(sampled);
        // up to the acoustic-sum-rule correction the mixture applies (the file's own
        // small violation at Gamma, removed at every q): ~1e-5 relative
        mixed.eigenvalues.forEach((row, k) => row.forEach((v, n) => {
            if (Math.abs(v) > 1) {
                assert.ok(Math.abs(v - own[k][n]) < 1e-2);
            }
        }));
    });
});

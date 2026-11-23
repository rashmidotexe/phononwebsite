import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, it } from 'node:test';

import {
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

    it('refuses a perovskite with an antiperovskite', () => {
        let srtio3 = endpoint('phonondb-SrTiO3-mp-5229.json.gz', 'SrTiO3');
        let li3clo = endpoint('phonondb-Li3ClO-mp-985585.json.gz', 'Li3ClO');
        let catio3 = endpoint('phonondb-CaTiO3-mp-4019.json.gz', 'CaTiO3');
        assert.throws(() => sampleMixingPair(srtio3, li3clo), /Cannot mix/);
        assert.throws(() => sampleMixingPair(li3clo, catio3), /Cannot mix/);
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

    it('does not mix different polytypes', () => {
        let ws2 = endpoint('phonondb-WS2-2H-mp-224.json.gz', 'WS2');
        let hfs2 = endpoint('phonondb-HfS2-1T-mp-985829.json.gz', 'HfS2');
        let mos2 = endpoint('phonondb-MoS2-3R-mp-1434.json.gz', 'MoS2');
        assert.throws(() => sampleMixingPair(ws2, mos2), /Cannot mix/);
        assert.throws(() => sampleMixingPair(hfs2, mos2), /Cannot mix/);
        assert.throws(() => sampleMixingPair(ws2, hfs2), /Cannot mix/);
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

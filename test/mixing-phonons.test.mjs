import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, it } from 'node:test';

import {
    buildMixingEndpoint,
    buildCommonPath,
    sampleMixingPair,
    computeMixedPhonon,
    computeReferenceEigenvalues,
    zincblendeAsChalcopyrite,
} from '../src/mixingphonons.js';
import { solveComplexHermitianWithEigenWasm } from '../src/eigenwasm.js';
import * as utils from '../src/utils.js';
import { alphaIdToNumber, indexDielectricRows } from '../src/mpdielectric.js';

const __dirname = dirname(fileURLToPath(import.meta.url));

function loadFixture(name) {
    return JSON.parse(gunzipSync(readFileSync(join(__dirname, 'fixtures', 'mixing', name))));
}

function mpDielectric() {
    // GaAs and AlN rows of the MP DFPT phonon collection, as read from its parquet file
    let rows = JSON.parse(readFileSync(join(__dirname, 'fixtures', 'mixing', 'mp-dfpt-dielectric-GaAs-AlN.json')));
    return indexDielectricRows(rows);
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

describe('zinc blende dynamical-matrix mixing', () => {
    it('restricts a PhononDB / Materials Project pair to the shorter common path', () => {
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        let gaas = endpoint('mpdb-GaAs-mp-2534.json.gz', 'GaAs');
        assert.equal(gap.segments.length, 6);
        assert.equal(gaas.segments.length, 10);

        let labels = (path) => path.map((s) => s.labelA + '-' + s.labelB).join(' ');
        let expected = 'GAMMA-X X-K K-GAMMA GAMMA-L L-W W-X';
        assert.equal(labels(buildCommonPath(gap, gaas)), expected);
        assert.equal(labels(buildCommonPath(gaas, gap)), expected);
    });

    it('is symmetric under swapping the materials and x -> 1-x', async () => {
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        let gaas = endpoint('mpdb-GaAs-mp-2534.json.gz', 'GaAs');

        let forward = await computeMixedPhonon(sampleMixingPair(gap, gaas), 0.3);
        let backward = await computeMixedPhonon(sampleMixingPair(gaas, gap), 0.7);

        assert.deepEqual(forward.highsym_qpts, backward.highsym_qpts);
        // acoustic modes at Gamma are sqrt(~1e-16 roundoff) ~ 1e-5 cm-1, the rest agree to ~1e-12
        assert.ok(maxAbsDiff(forward.eigenvalues, backward.eigenvalues) < 1e-4);
        assert.ok(maxAbsDiff([forward.distances], [backward.distances]) < 1e-12);
    });

    it('keeps the acoustic modes at zero at Gamma for every x', async () => {
        // GaP and GaAs have very different anion/cation mass ratios, which is
        // exactly the case where averaging the mass-weighted D breaks the sum rule
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        let gaas = endpoint('mpdb-GaAs-mp-2534.json.gz', 'GaAs');
        let sampled = sampleMixingPair(gap, gaas);

        for (let x of [0, 0.25, 0.5, 0.75, 1]) {
            let mixed = await computeMixedPhonon(sampled, x);
            let gamma = mixed.highsym_qpts.filter((point) => point[1] === 'GAMMA').map((point) => point[0]);
            assert.ok(gamma.length > 0);
            for (let k of gamma) {
                for (let n = 0; n < 3; n++) {
                    assert.ok(Math.abs(mixed.eigenvalues[k][n]) < 0.5, `x=${x} q=${k} mode ${n}: ${mixed.eigenvalues[k][n]}`);
                }
            }
        }
    });

    it('reproduces the Materials Project frequencies at x = 0', async () => {
        let raw = loadFixture('mpdb-GaAs-mp-2534.json.gz');
        let gaas = buildMixingEndpoint(raw, { name: 'GaAs' });
        let aln = endpoint('mpdb-AlN-mp-1700.json.gz', 'AlN');
        let sampled = sampleMixingPair(gaas, aln);
        let mixed = await computeMixedPhonon(sampled, 0);
        let reference = await computeReferenceEigenvalues(sampled);

        let matched = 0;
        for (let k = 0; k < mixed.qpoints.length; k++) {
            let i = raw.qpoints.findIndex((q) => q.every((v, d) => Math.abs(v - mixed.qpoints[k][d]) < 1e-6));
            if (i < 0) {
                continue;
            }
            matched += 1;
            let expected = raw.frequencies.map((band) => band[i] * 33.35641).sort((a, b) => a - b);
            for (let n = 0; n < expected.length; n++) {
                assert.ok(Math.abs(expected[n] - mixed.eigenvalues[k][n]) < 1e-3);
            }
        }
        // AlN is sampled more densely, so only the segment ends are file q-points of GaAs
        assert.ok(matched >= 2 * sampled.lineBreaks.length);
        assert.ok(maxAbsDiff(mixed.eigenvalues, reference) < 1e-4);
    });

    it('recomputes the dipole-dipole term for PhononDB pairs and keeps the end-members exact', async () => {
        // GaP (P at -1/4, inverted) and AlN (N at +1/4) both carry Born charges
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        let aln = endpoint('phonondb-AlN-mp-1700.json.gz', 'AlN');
        let sampled = sampleMixingPair(gap, aln);
        assert.equal(sampled.recomputeNac, true);
        let interpolated = Object.assign({}, sampled, { recomputeNac: false });

        let opticalDiff = (a, b) => {
            let worst = 0;
            a.forEach((row, k) => row.forEach((v, n) => {
                if (Math.abs(v) > 1) {
                    worst = Math.max(worst, Math.abs(v - b[k][n]));
                }
            }));
            return worst;
        };

        for (let x of [0, 1]) {
            let recomputed = await computeMixedPhonon(sampled, x);
            assert.equal(recomputed.nac_mode, 'recomputed');
            // short-range part + rebuilt dipole term == the original dynamical matrix
            assert.ok(opticalDiff(recomputed.eigenvalues, (await computeMixedPhonon(interpolated, x)).eigenvalues) < 1e-6);
        }

        let half = await computeMixedPhonon(sampled, 0.5);
        let gamma = half.highsym_qpts.filter((point) => point[1] === 'GAMMA').map((point) => point[0]);
        for (let k of gamma) {
            for (let n = 0; n < 3; n++) {
                assert.ok(Math.abs(half.eigenvalues[k][n]) < 0.5);
            }
        }
        // AlN and GaP have very different dielectric constants: the LO mode at Gamma
        // must differ from simply interpolating the long-range term
        let halfInterpolated = await computeMixedPhonon(interpolated, 0.5);
        assert.ok(Math.abs(half.eigenvalues[gamma[0]][5] - halfInterpolated.eigenvalues[gamma[0]][5]) > 5);

        let mirrored = await computeMixedPhonon(sampleMixingPair(aln, gap), 0.5);
        assert.ok(opticalDiff(half.eigenvalues, mirrored.eigenvalues) < 1e-6);
    });

    it('interpolates the long-range term when Materials Project Born charges are missing', async () => {
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        let gaas = endpoint('mpdb-GaAs-mp-2534.json.gz', 'GaAs');
        let mixed = await computeMixedPhonon(sampleMixingPair(gap, gaas), 0.5);
        assert.equal(mixed.nac_mode, 'interpolated');
    });

    it('decodes Materials Project AlphaIDs', () => {
        assert.equal(alphaIdToNumber('aaaaaaaa'), 0);
        assert.equal(alphaIdToNumber('aaaaadtm'), 2534);
        assert.equal(alphaIdToNumber('aaaaacnk'), 1700);
    });

    it('splits Materials Project D into a direction-independent short-range part', () => {
        let table = mpDielectric();
        let raw = loadFixture('mpdb-GaAs-mp-2534.json.gz');
        let gaas = buildMixingEndpoint(raw, { name: 'GaAs' }, table.get(2534));
        assert.ok(gaas.nac);

        // the full D at Gamma depends on the approach direction (LO-TO), the
        // short-range remainder must not
        let atGamma = [];
        for (let segment of gaas.segments) {
            if (segment.qA.every((v) => Math.abs(v) < 1e-8)) atGamma.push([segment, 0]);
            if (segment.qB.every((v) => Math.abs(v) < 1e-8)) atGamma.push([segment, 1]);
        }
        assert.ok(atGamma.length >= 3);
        let spread = (matrices) => {
            let worst = 0;
            let scale = 0;
            for (let m of matrices) {
                m.real.flat().forEach((v, i) => {
                    scale = Math.max(scale, Math.abs(v));
                    worst = Math.max(worst, Math.abs(v - matrices[0].real.flat()[i]));
                });
            }
            return worst / scale;
        };
        assert.ok(spread(atGamma.map(([s, t]) => s.sample(t))) > 1e-2);
        assert.ok(spread(atGamma.map(([s, t]) => s.sampleShortRange(t))) < 1e-5);
    });

    it('recomputes the dipole-dipole term for Materials Project pairs', async () => {
        let table = mpDielectric();
        let gaasRaw = loadFixture('mpdb-GaAs-mp-2534.json.gz');
        let gaas = buildMixingEndpoint(gaasRaw, { name: 'GaAs' }, table.get(2534));
        let aln = buildMixingEndpoint(loadFixture('mpdb-AlN-mp-1700.json.gz'), { name: 'AlN' }, table.get(1700));
        let sampled = sampleMixingPair(gaas, aln);
        assert.equal(sampled.recomputeNac, true);

        // at x = 0, on the q-points of the GaAs file, the file frequencies come back
        let mixed = await computeMixedPhonon(sampled, 0);
        assert.equal(mixed.nac_mode, 'recomputed');
        let matched = 0;
        for (let k = 0; k < mixed.qpoints.length; k++) {
            let i = gaasRaw.qpoints.findIndex((q) => q.every((v, d) => Math.abs(v - mixed.qpoints[k][d]) < 1e-6));
            if (i < 0) {
                continue;
            }
            matched += 1;
            let expected = gaasRaw.frequencies.map((band) => band[i] * 33.35641).sort((a, b) => a - b);
            for (let n = 3; n < expected.length; n++) {
                assert.ok(Math.abs(expected[n] - mixed.eigenvalues[k][n]) < 1e-3);
            }
        }
        assert.ok(matched >= 2 * sampled.lineBreaks.length);

        let half = await computeMixedPhonon(sampled, 0.5);
        let gamma = half.highsym_qpts.filter((point) => point[1] === 'GAMMA').map((point) => point[0]);
        for (let k of gamma) {
            for (let n = 0; n < 3; n++) {
                assert.ok(Math.abs(half.eigenvalues[k][n]) < 0.05);
            }
        }
        let mirrored = await computeMixedPhonon(sampleMixingPair(aln, gaas), 0.5);
        assert.ok(maxAbsDiff(half.eigenvalues, mirrored.eigenvalues) < 1e-4);
    });

    it('tunes the site masses independently of the character', async () => {
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        let aln = endpoint('phonondb-AlN-mp-1700.json.gz', 'AlN');
        let sampled = sampleMixingPair(gap, aln);

        // by default the masses follow x
        let linked = await computeMixedPhonon(sampled, 0.5);
        assert.ok(maxAbsDiff(linked.eigenvalues, (await computeMixedPhonon(sampled, 0.5, 0.5)).eigenvalues) < 1e-12);

        let heavy = await computeMixedPhonon(sampled, 0.5, 0);
        let light = await computeMixedPhonon(sampled, 0.5, 1);
        let gamma = heavy.highsym_qpts.filter((point) => point[1] === 'GAMMA').map((point) => point[0]);
        for (let mixed of [heavy, light]) {
            // same force constants, only the masses differ: the sum rule still holds
            for (let k of gamma) {
                for (let n = 0; n < 3; n++) {
                    assert.ok(Math.abs(mixed.eigenvalues[k][n]) < 0.5);
                }
            }
            // the bonding (lattice) still follows x
            assert.deepEqual(mixed.lattice, linked.lattice);
        }
        // GaP masses (Ga, P) are heavier than AlN masses (Al, N): lower optical modes
        assert.ok(heavy.eigenvalues[gamma[0]][5] < linked.eigenvalues[gamma[0]][5]);
        assert.ok(linked.eigenvalues[gamma[0]][5] < light.eigenvalues[gamma[0]][5]);

        // a -> b at (x, y) is b -> a at (1 - x, 1 - y)
        let mirrored = await computeMixedPhonon(sampleMixingPair(aln, gap), 0.7, 0.8);
        let forward = await computeMixedPhonon(sampled, 0.3, 0.2);
        let optical = (a, b) => {
            let worst = 0;
            a.forEach((row, k) => row.forEach((v, n) => {
                if (Math.abs(v) > 1) {
                    worst = Math.max(worst, Math.abs(v - b[k][n]));
                }
            }));
            return worst;
        };
        assert.ok(optical(forward.eigenvalues, mirrored.eigenvalues) < 1e-6);
    });

    it('maps an inverted setting of the same compound onto the other one', async () => {
        // AlN is stored with N at +1/4 in PhononDB and at -1/4 in Materials Project
        let pdb = endpoint('phonondb-AlN-mp-1700.json.gz', 'AlN');
        let mp = endpoint('mpdb-AlN-mp-1700.json.gz', 'AlN');
        assert.equal(pdb.invert, false);
        assert.equal(mp.invert, true);
        assert.deepEqual(pdb.siteTypes, ['Al', 'N']);
        assert.deepEqual(mp.siteTypes, ['Al', 'N']);

        let sampled = sampleMixingPair(pdb, mp);
        let end0 = await computeMixedPhonon(sampled, 0);
        let end1 = await computeMixedPhonon(sampled, 1);
        let half = await computeMixedPhonon(sampled, 0.5);
        let mean = end0.eigenvalues.map((row, k) => row.map((v, n) => (v + end1.eigenvalues[k][n]) / 2));

        // the two DFT data sets differ by ~13 cm-1, the mixture must sit in between
        assert.ok(maxAbsDiff(half.eigenvalues, mean) < 2);
    });
});

describe('chalcopyrite dynamical-matrix mixing', () => {
    let gammaIndices = (phonon) => phonon.highsym_qpts.filter((point) => point[1] === 'GAMMA').map((point) => point[0]);
    let opticalDiff = (a, b) => {
        let worst = 0;
        a.forEach((row, k) => row.forEach((v, n) => {
            if (Math.abs(v) > 1) {
                worst = Math.max(worst, Math.abs(v - b[k][n]));
            }
        }));
        return worst;
    };

    function zngep2Dielectric() {
        let rows = JSON.parse(readFileSync(join(__dirname, 'fixtures', 'mixing', 'mp-dfpt-dielectric-ZnGeP2.json')));
        return indexDielectricRows(rows).get(4524);
    }

    it('matches 8-atom cells onto the I-42d template, in either orientation', () => {
        let cugas2 = endpoint('phonondb-CuGaS2-mp-5238.json.gz', 'CuGaS2');
        let cualse2 = endpoint('phonondb-CuAlSe2-mp-8016.json.gz', 'CuAlSe2');
        assert.equal(cugas2.family, 'chalcopyrite');
        assert.deepEqual(cugas2.siteTypes, ['Cu', 'Cu', 'Ga', 'Ga', 'S', 'S', 'S', 'S']);
        assert.deepEqual(cualse2.siteTypes, ['Cu', 'Cu', 'Al', 'Al', 'Se', 'Se', 'Se', 'Se']);
        // CuAlSe2 is stored as the inversion image of the template
        assert.equal(cugas2.invert, false);
        assert.equal(cualse2.invert, true);
        assert.ok(Math.abs(cugas2.structure.u - 0.2459) < 1e-3);
        // seekpath names of the file are replaced by coordinate-based ones
        let labels = cugas2.segments.map((s) => s.labelA + '-' + s.labelB);
        assert.ok(labels.includes('Σ-GAMMA'));
        assert.ok(!labels.some((label) => label.includes('S_0')));
    });

    it('folds a zinc blende into the chalcopyrite cell', async () => {
        // D of the 8-atom supercell at Q has the eigenvalues of the zinc blende at the
        // 4 points Q + g, also at Q = 0 for a polar crystal
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        let folded = zincblendeAsChalcopyrite(gap);
        assert.equal(folded.family, 'chalcopyrite');
        assert.deepEqual(folded.siteTypes, ['Ga', 'Ga', 'Ga', 'Ga', 'P', 'P', 'P', 'P']);
        assert.ok(Math.abs(folded.structure.c - 2 * folded.structure.a) < 1e-12);

        let eigen = async (m) => (await solveComplexHermitianWithEigenWasm(m.real, m.imag)).values.sort((p, q) => p - q);
        let superReciprocal = utils.rec_lat(folded.lattice);
        let candidates = [];
        for (let a = -2; a <= 2; a++) {
            for (let b = -2; b <= 2; b++) {
                for (let c = -2; c <= 2; c++) {
                    candidates.push([0, 1, 2].map((d) => a * superReciprocal[0][d] + b * superReciprocal[1][d] + c * superReciprocal[2][d]));
                }
            }
        }
        candidates.sort((g, h) => g.reduce((s, v) => s + v * v, 0) - h.reduce((s, v) => s + v * v, 0));
        let direction = [0, 0, 1];
        for (let Q of [[0, 0, 0], [0.13, -0.07, 0.21], [0.5, 0.5, -0.5]]) {
            let qCart = [0, 1, 2].map((d) => Q[0] * superReciprocal[0][d] + Q[1] * superReciprocal[1][d] + Q[2] * superReciprocal[2][d]);
            let seen = new Set();
            let expected = [];
            for (let g of candidates) {
                let key = gap.lattice.map((p) => {
                    let f = g[0] * p[0] + g[1] * p[1] + g[2] * p[2];
                    return Math.round((f - Math.floor(f + 1e-9)) * 1e6) % 1000000;
                }).join(',');
                if (seen.has(key)) {
                    continue;
                }
                seen.add(key);
                let k = gap.lattice.map((p) => p[0] * (qCart[0] + g[0]) + p[1] * (qCart[1] + g[1]) + p[2] * (qCart[2] + g[2]));
                let directionCart = [0, 1, 2].map((d) => direction[0] * superReciprocal[0][d] + direction[1] * superReciprocal[1][d] + direction[2] * superReciprocal[2][d]);
                let directionPrimitive = gap.lattice.map((p) => p[0] * directionCart[0] + p[1] * directionCart[1] + p[2] * directionCart[2]);
                expected.push(...await eigen(gap.evaluate(k, directionPrimitive)));
            }
            expected.sort((p, q) => p - q);
            let actual = await eigen(folded.sampleAtQ(Q, direction));
            actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-9 * (1 + Math.abs(expected[i]))));
        }
    });

    it('mixes a PhononDB zinc blende with a chalcopyrite', async () => {
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        let cugas2 = endpoint('phonondb-CuGaS2-mp-5238.json.gz', 'CuGaS2');
        let sampled = sampleMixingPair(gap, cugas2);
        let reverse = sampleMixingPair(cugas2, gap);
        assert.equal(sampled.recomputeNac, true);
        // the folded zinc blende can be evaluated anywhere: the whole chalcopyrite path
        assert.deepEqual(sampled.labels, cugas2.segments.map((s) => [s.labelA, s.labelB]));

        let zincblendeEnd = await computeMixedPhonon(sampled, 0);
        assert.ok(opticalDiff(zincblendeEnd.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-6);
        let chalcopyriteEnd = await computeMixedPhonon(reverse, 0);
        assert.ok(opticalDiff(chalcopyriteEnd.eigenvalues, await computeReferenceEigenvalues(reverse)) < 1e-6);

        let half = await computeMixedPhonon(sampled, 0.5, 0.3);
        assert.equal(half.natoms, 8);
        for (let k of gammaIndices(half)) {
            for (let n = 0; n < 3; n++) {
                assert.ok(Math.abs(half.eigenvalues[k][n]) < 0.05);
            }
        }
        let mirrored = await computeMixedPhonon(reverse, 0.5, 0.7);
        assert.ok(opticalDiff(half.eigenvalues, mirrored.eigenvalues) < 1e-6);
    });

    it('refuses a Materials Project zinc blende with a chalcopyrite', () => {
        // its D is only known on the band path, folding needs it at Q + g
        let gaas = endpoint('mpdb-GaAs-mp-2534.json.gz', 'GaAs');
        let cugas2 = endpoint('phonondb-CuGaS2-mp-5238.json.gz', 'CuGaS2');
        assert.throws(() => sampleMixingPair(gaas, cugas2), /force constants/);
    });

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        let cugas2 = endpoint('phonondb-CuGaS2-mp-5238.json.gz', 'CuGaS2');
        let cualse2 = endpoint('phonondb-CuAlSe2-mp-8016.json.gz', 'CuAlSe2');
        let sampled = sampleMixingPair(cugas2, cualse2);
        assert.equal(sampled.recomputeNac, true);

        // short-range part + rebuilt dipole term (with phonopy's hermitian dd_q0)
        // gives back each material's own dynamical matrix
        let end0 = await computeMixedPhonon(sampled, 0);
        assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-6);
        let end1 = await computeMixedPhonon(sampled, 1);
        let mirroredEnd = await computeMixedPhonon(sampleMixingPair(cualse2, cugas2), 0);
        assert.ok(opticalDiff(end1.eigenvalues, mirroredEnd.eigenvalues) < 1e-6);

        for (let x of [0.25, 0.5, 0.75]) {
            let mixed = await computeMixedPhonon(sampled, x, 1 - x);
            assert.equal(mixed.natoms, 8);
            for (let k of gammaIndices(mixed)) {
                for (let n = 0; n < 3; n++) {
                    assert.ok(Math.abs(mixed.eigenvalues[k][n]) < 0.05, `x=${x} mode ${n}: ${mixed.eigenvalues[k][n]}`);
                }
            }
        }

        let forward = await computeMixedPhonon(sampled, 0.3, 0.6);
        let mirrored = await computeMixedPhonon(sampleMixingPair(cualse2, cugas2), 0.7, 0.4);
        assert.ok(opticalDiff(forward.eigenvalues, mirrored.eigenvalues) < 1e-6);
    });

    it('splits Materials Project chalcopyrite D with per-site Born charges', () => {
        let zngep2 = buildMixingEndpoint(loadFixture('mpdb-ZnGeP2-mp-4524-trimmed.json.gz'), { name: 'ZnGeP2' }, zngep2Dielectric());
        assert.ok(zngep2.nac);
        // the two Zn (and the four P) carry differently oriented Born tensors
        assert.notDeepEqual(zngep2.nac.born[0], zngep2.nac.born[1]);

        // the short-range remainder is the same at Gamma from every direction,
        // while the full D carries the direction-dependent LO-TO splitting
        let atGamma = [];
        let fullAtGamma = [];
        for (let segment of zngep2.segments) {
            for (let t of [0, 1]) {
                let q = t ? segment.qB : segment.qA;
                if (q.every((v) => Math.abs(v) < 1e-9)) {
                    atGamma.push(segment.sampleShortRange(t));
                    fullAtGamma.push(segment.sample(t));
                }
            }
        }
        assert.equal(atGamma.length, 3);
        let relative = (a, b) => {
            let num = 0;
            let den = 0;
            for (let i = 0; i < 24; i++) {
                for (let j = 0; j < 24; j++) {
                    num += (a.real[i][j] - b.real[i][j]) ** 2 + (a.imag[i][j] - b.imag[i][j]) ** 2;
                    den += a.real[i][j] ** 2;
                }
            }
            return Math.sqrt(num / den);
        };
        for (let i = 1; i < atGamma.length; i++) {
            assert.ok(relative(atGamma[0], atGamma[i]) < 1e-4);
            assert.ok(relative(fullAtGamma[0], fullAtGamma[i]) > 1e-2);
        }
    });

    it('cuts segments ending at c/a-dependent points at the nearer one', () => {
        let cugas2 = endpoint('phonondb-CuGaS2-mp-5238.json.gz', 'CuGaS2');
        let zngep2 = buildMixingEndpoint(loadFixture('mpdb-ZnGeP2-mp-4524-trimmed.json.gz'), { name: 'ZnGeP2' }, zngep2Dielectric());
        let path = buildCommonPath(cugas2, zngep2);
        assert.deepEqual(path.map((s) => s.labelA + '-' + s.labelB), ['GAMMA-X', 'Σ-GAMMA', 'GAMMA-Z']);

        // Sigma = (-eta, eta, eta) with eta = (1 + a^2/c^2)/4: the larger c/a gives
        // the point nearer to Gamma, where the common stretch ends
        let eta = (s) => (1 + (s.a * s.a) / (s.c * s.c)) / 4;
        let nearer = Math.min(eta(cugas2.structure), eta(zngep2.structure));
        let sigmaGamma = path[1];
        assert.ok(Math.abs(sigmaGamma.qA[1] - nearer) < 1e-6);

        // a -> b and b -> a use the same q-points
        let reverse = buildCommonPath(zngep2, cugas2);
        assert.deepEqual(reverse.map((s) => s.labelA + '-' + s.labelB), path.map((s) => s.labelA + '-' + s.labelB));
    });
});

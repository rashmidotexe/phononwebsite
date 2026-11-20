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
    wurtziteAsSuperstructure,
    binaryAsSuperstructure,
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

describe('wurtzite dynamical-matrix mixing', () => {
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
    let relative = (a, b) => {
        let num = 0;
        let den = 0;
        for (let i = 0; i < a.real.length; i++) {
            for (let j = 0; j < a.real.length; j++) {
                num += (a.real[i][j] - b.real[i][j]) ** 2 + (a.imag[i][j] - b.imag[i][j]) ** 2;
                den += a.real[i][j] ** 2 + a.imag[i][j] ** 2;
            }
        }
        return Math.sqrt(num / den);
    };
    let median = (values) => values.slice().sort((p, q) => p - q)[values.length >> 1];

    function mpAlN() {
        let rows = JSON.parse(readFileSync(join(__dirname, 'fixtures', 'mixing', 'mp-dfpt-dielectric-AlN-wurtzite.json')));
        let dielectric = indexDielectricRows(rows).get(661);
        return buildMixingEndpoint(loadFixture('mpdb-AlN-wurtzite-mp-661-trimmed.json.gz'), { name: 'AlN' }, dielectric);
    }

    it('matches 4-atom cells onto the P6_3mc template, with either polarity', () => {
        let gan = endpoint('phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN');
        let aln = endpoint('phonondb-AlN-wurtzite-mp-661.json.gz', 'AlN');
        let alnMP = mpAlN();
        for (let e of [gan, aln, alnMP]) {
            assert.equal(e.family, 'wurtzite');
            assert.equal(e.siteTypes[0], e.siteTypes[1]);
            assert.equal(e.siteTypes[2], e.siteTypes[3]);
        }
        // the GaN file lists N first and has the opposite polarity
        assert.deepEqual(gan.siteTypes, ['Ga', 'Ga', 'N', 'N']);
        assert.equal(gan.invert, true);
        assert.equal(aln.invert, false);
        assert.ok(Math.abs(aln.structure.u - 0.382) < 1e-3);
        assert.ok(Math.abs(aln.structure.c / aln.structure.a - 1.600) < 1e-3);
        // Materials Project stores a1 at -60 degrees from x: rotated onto the PhononDB axes
        assert.equal(aln.rotation, null);
        assert.ok(alnMP.rotation);
        assert.ok(alnMP.nac);
        let labels = aln.segments.map((s) => s.labelA + '-' + s.labelB);
        assert.deepEqual(labels, ['GAMMA-M', 'M-K', 'K-GAMMA', 'GAMMA-A', 'A-L', 'L-H', 'H-A', 'L-M', 'H-K']);
    });

    it('maps the Materials Project setting of a compound onto the PhononDB one', () => {
        let pdb = endpoint('phonondb-AlN-wurtzite-mp-661.json.gz', 'AlN');
        let mp = mpAlN();
        let differences = (sampled) => sampled.matrices1.map((m, k) => relative(m, sampled.matrices2[k]));

        // two DFT data sets of the same AlN agree to a few per cent
        let good = differences(sampleMixingPair(pdb, mp));
        assert.ok(median(good) < 0.03, 'median ' + median(good));
        // the opposite polarity is a different crystal orientation and does not
        mp.invert = !mp.invert;
        assert.ok(median(differences(sampleMixingPair(pdb, mp))) > 0.3);
    });

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        let gan = endpoint('phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN');
        let aln = endpoint('phonondb-AlN-wurtzite-mp-661.json.gz', 'AlN');
        let sampled = sampleMixingPair(gan, aln);
        assert.equal(sampled.recomputeNac, true);

        let end0 = await computeMixedPhonon(sampled, 0);
        assert.equal(end0.natoms, 4);
        assert.equal(end0.family, 'wurtzite');
        assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-4);
        let end1 = await computeMixedPhonon(sampled, 1);
        let mirroredEnd = await computeMixedPhonon(sampleMixingPair(aln, gan), 0);
        assert.ok(opticalDiff(end1.eigenvalues, mirroredEnd.eigenvalues) < 1e-4);

        for (let x of [0.25, 0.5, 0.75]) {
            let mixed = await computeMixedPhonon(sampled, x, 1 - x);
            for (let k of gammaIndices(mixed)) {
                for (let n = 0; n < 3; n++) {
                    assert.ok(Math.abs(mixed.eigenvalues[k][n]) < 0.05, `x=${x} mode ${n}: ${mixed.eigenvalues[k][n]}`);
                }
            }
        }

        let forward = await computeMixedPhonon(sampled, 0.3, 0.6);
        let mirrored = await computeMixedPhonon(sampleMixingPair(aln, gan), 0.7, 0.4);
        assert.ok(opticalDiff(forward.eigenvalues, mirrored.eigenvalues) < 1e-6);
    });

    it('takes the same path whichever material is material 1', async () => {
        // equal segment counts: the path must not follow material 1, the files run
        // their last segment in opposite directions
        let gan = endpoint('phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN');
        let alnMP = mpAlN();
        let labels = (path) => path.map((s) => s.labelA + '-' + s.labelB);
        assert.deepEqual(labels(buildCommonPath(gan, alnMP)), labels(buildCommonPath(alnMP, gan)));

        let forward = await computeMixedPhonon(sampleMixingPair(gan, alnMP), 0.4);
        let mirrored = await computeMixedPhonon(sampleMixingPair(alnMP, gan), 0.6);
        assert.ok(opticalDiff(forward.eigenvalues, mirrored.eigenvalues) < 1e-6);
    });

    it('does not mix a wurtzite with a zinc blende', () => {
        let gan = endpoint('phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN');
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        assert.throws(() => sampleMixingPair(gan, gap), /Cannot mix/);
        assert.throws(() => sampleMixingPair(gap, gan), /Cannot mix/);
    });
});

describe('ternary wurtzite dynamical-matrix mixing', () => {
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
    let eigen = async (m) => (await solveComplexHermitianWithEigenWasm(m.real, m.imag)).values.sort((p, q) => p - q);
    const ORTHORHOMBIC_PATH = ['GAMMA-X', 'X-S', 'S-Y', 'Y-GAMMA', 'GAMMA-Z', 'Z-U', 'U-R', 'R-T', 'T-Z', 'Y-T', 'U-X', 'S-R'];

    it('matches beta-NaFeO2 type cells onto the Pna2_1 template', () => {
        let zngen2 = endpoint('phonondb-ZnGeN2-mp-2979.json.gz', 'ZnGeN2');
        let ligas2 = endpoint('phonondb-LiGaS2-mp-3647.json.gz', 'LiGaS2');
        for (let e of [zngen2, ligas2]) {
            assert.equal(e.family, 'wurtzite-ABX2');
            assert.ok(e.sampleAtQ);
            // the canonical cell lies along the axes: 2a x sqrt(3)a x c of the wurtzite
            let lattice = e.structure.lattice;
            lattice.forEach((row, i) => row.forEach((v, j) => {
                if (i !== j) {
                    assert.ok(Math.abs(v) < 1e-9);
                }
            }));
            assert.ok(Math.abs(lattice[0][0] / lattice[1][1] - 2 / Math.sqrt(3)) < 0.08);
        }
        assert.deepEqual(zngen2.siteTypes, ['Zn', 'Zn', 'Zn', 'Zn', 'Ge', 'Ge', 'Ge', 'Ge', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N']);
        assert.deepEqual(ligas2.siteTypes.slice(0, 8), ['Li', 'Li', 'Li', 'Li', 'Ga', 'Ga', 'Ga', 'Ga']);
        // one of the two needs an improper rotation onto the template: its image under inversion
        assert.notEqual(zngen2.invert, ligas2.invert);
    });

    it('folds a binary wurtzite into the ternary cell', async () => {
        // D of the 16-atom supercell at Q has the eigenvalues of GaN at the 4 points Q + g
        let gan = endpoint('phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN');
        let folded = wurtziteAsSuperstructure(gan, 'wurtzite-ABX2');
        assert.equal(folded.siteTypes.length, 16);
        assert.equal(folded.siteTypes.filter((type) => type === 'Ga').length, 8);

        let superReciprocal = utils.rec_lat(folded.lattice);
        let toCart = (q) => [0, 1, 2].map((d) => q[0] * superReciprocal[0][d] + q[1] * superReciprocal[1][d] + q[2] * superReciprocal[2][d]);
        let toPrimitive = (k) => gan.lattice.map((p) => p[0] * k[0] + p[1] * k[1] + p[2] * k[2]);
        for (let Q of [[0.13, -0.07, 0.21], [0.5, 0.5, 0]]) {
            let qCart = toCart(Q);
            let expected = [];
            // the supercell is 4 wurtzite cells: one g per class modulo the wurtzite
            // reciprocal lattice, the shortest as in the folding (the truncated
            // dipole-dipole sum makes D(k) periodic only to ~1e-5)
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
            for (let g of candidates) {
                let key = toPrimitive(g).map((v) => Math.round((v - Math.floor(v + 1e-9)) * 1e6) % 1000000).join(',');
                if (seen.has(key)) {
                    continue;
                }
                seen.add(key);
                let k = toPrimitive([0, 1, 2].map((d) => qCart[d] + g[d]));
                expected.push(...await eigen(gan.evaluate(k, null)));
            }
            assert.equal(seen.size, 4);
            expected.sort((p, q) => p - q);
            let actual = await eigen(folded.sampleAtQ(Q, null));
            actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-9 * (1 + Math.abs(expected[i]))));
        }
    });

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        // the two files run their paths along different (equivalent) axes: both are
        // evaluated on the standard path of the canonical cell
        let zngen2 = endpoint('phonondb-ZnGeN2-mp-2979.json.gz', 'ZnGeN2');
        let ligas2 = endpoint('phonondb-LiGaS2-mp-3647.json.gz', 'LiGaS2');
        let sampled = sampleMixingPair(zngen2, ligas2);
        assert.equal(sampled.recomputeNac, true);
        assert.deepEqual(sampled.labels.map((l) => l.join('-')), ORTHORHOMBIC_PATH);

        let end0 = await computeMixedPhonon(sampled, 0);
        assert.equal(end0.natoms, 16);
        assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-4);
        let end1 = await computeMixedPhonon(sampled, 1);
        let mirroredEnd = await computeMixedPhonon(sampleMixingPair(ligas2, zngen2), 0);
        assert.ok(opticalDiff(end1.eigenvalues, mirroredEnd.eigenvalues) < 1e-4);

        let mixed = await computeMixedPhonon(sampled, 0.5, 0.3);
        for (let k of gammaIndices(mixed)) {
            for (let n = 0; n < 3; n++) {
                assert.ok(Math.abs(mixed.eigenvalues[k][n]) < 0.05, `mode ${n}: ${mixed.eigenvalues[k][n]}`);
            }
        }
        let mirrored = await computeMixedPhonon(sampleMixingPair(ligas2, zngen2), 0.5, 0.7);
        assert.ok(opticalDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-6);
    });

    it('mixes a PhononDB wurtzite with a ternary wurtzite', async () => {
        let gan = endpoint('phonondb-GaN-wurtzite-mp-804.json.gz', 'GaN');
        let zngen2 = endpoint('phonondb-ZnGeN2-mp-2979.json.gz', 'ZnGeN2');
        let sampled = sampleMixingPair(gan, zngen2);
        assert.equal(sampled.endpoint1.family, 'wurtzite-ABX2');
        assert.deepEqual(sampled.labels.map((l) => l.join('-')), ORTHORHOMBIC_PATH);

        // GaN end: the folded GaN; ZnGeN2 end: ZnGeN2 itself
        let end0 = await computeMixedPhonon(sampled, 0);
        assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-4);
        let end1 = await computeMixedPhonon(sampled, 1);
        let own = await computeReferenceEigenvalues(sampleMixingPair(zngen2, gan));
        assert.ok(opticalDiff(end1.eigenvalues, own) < 1e-4);

        let mixed = await computeMixedPhonon(sampled, 0.5);
        for (let k of gammaIndices(mixed)) {
            for (let n = 0; n < 3; n++) {
                assert.ok(Math.abs(mixed.eigenvalues[k][n]) < 0.05);
            }
        }
    });

    it('refuses pairs on different lattices', () => {
        let zngen2 = endpoint('phonondb-ZnGeN2-mp-2979.json.gz', 'ZnGeN2');
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        let cugas2 = endpoint('phonondb-CuGaS2-mp-5238.json.gz', 'CuGaS2');
        assert.throws(() => sampleMixingPair(zngen2, gap), /Cannot mix/);
        assert.throws(() => sampleMixingPair(cugas2, zngen2), /Cannot mix/);
        // a Materials Project wurtzite is only known on its band path
        let rows = JSON.parse(readFileSync(join(__dirname, 'fixtures', 'mixing', 'mp-dfpt-dielectric-AlN-wurtzite.json')));
        let alnMP = buildMixingEndpoint(loadFixture('mpdb-AlN-wurtzite-mp-661-trimmed.json.gz'), { name: 'AlN' }, indexDielectricRows(rows).get(661));
        assert.throws(() => sampleMixingPair(alnMP, zngen2), /force constants/);
    });
});

describe('rocksalt dynamical-matrix mixing', () => {
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
    let caoMP = () => {
        let rows = JSON.parse(readFileSync(join(__dirname, 'fixtures', 'mixing', 'mp-dfpt-dielectric-CaO.json')));
        return buildMixingEndpoint(loadFixture('mpdb-CaO-mp-2605-trimmed.json.gz'), { name: 'CaO' }, indexDielectricRows(rows).get(2605));
    };

    it('matches 2-atom cells with the anion at (1/2,1/2,1/2) as rocksalt', () => {
        let nacl = endpoint('phonondb-NaCl-mp-22862.json.gz', 'NaCl');
        assert.equal(nacl.family, 'rocksalt');
        assert.deepEqual(nacl.siteTypes, ['Na', 'Cl']);
        assert.equal(nacl.invert, false);
        // the cation is the species that is not a usual anion, whatever the name order (TePb)
        let named = buildMixingEndpoint(loadFixture('phonondb-NaCl-mp-22862.json.gz'), { name: 'ClNa' });
        assert.deepEqual(named.siteTypes, ['Na', 'Cl']);
        let cao = caoMP();
        assert.equal(cao.family, 'rocksalt');
        assert.deepEqual(cao.siteTypes, ['Ca', 'O']);
        assert.ok(cao.nac);
        // zinc blende stays zinc blende
        assert.equal(endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP').family, 'zincblende');
    });

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        let nacl = endpoint('phonondb-NaCl-mp-22862.json.gz', 'NaCl');
        let kcl = endpoint('phonondb-KCl-mp-23193.json.gz', 'KCl');
        let sampled = sampleMixingPair(nacl, kcl);
        assert.equal(sampled.recomputeNac, true);
        let end0 = await computeMixedPhonon(sampled, 0);
        assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-6);
        let end1 = await computeMixedPhonon(sampled, 1);
        assert.ok(opticalDiff(end1.eigenvalues, await computeReferenceEigenvalues(sampleMixingPair(kcl, nacl))) < 1e-6);

        let mixed = await computeMixedPhonon(sampled, 0.4, 0.7);
        let gamma = gammaIndices(sampled);
        assert.ok(gamma.length >= 2);
        for (let k of gamma) {
            for (let n = 0; n < 3; n++) {
                assert.ok(Math.abs(mixed.eigenvalues[k][n]) < 0.05);
            }
        }
        let mirrored = await computeMixedPhonon(sampleMixingPair(kcl, nacl), 0.6, 0.3);
        assert.ok(maxAbsDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-4);
    });

    it('recomputes the dipole-dipole term with a Materials Project rocksalt', async () => {
        let nacl = endpoint('phonondb-NaCl-mp-22862.json.gz', 'NaCl');
        let caoRaw = loadFixture('mpdb-CaO-mp-2605-trimmed.json.gz');
        let sampled = sampleMixingPair(nacl, caoMP());
        assert.equal(sampled.recomputeNac, true);

        // at x = 1, on the q-points of the CaO file away from Gamma, the file frequencies come back
        let mixed = await computeMixedPhonon(sampled, 1);
        let matched = 0;
        for (let k = 0; k < mixed.qpoints.length; k++) {
            if (gammaIndices(sampled).includes(k)) {
                continue;
            }
            let i = caoRaw.qpoints.findIndex((q) => q.every((v, d) => Math.abs(v - mixed.qpoints[k][d]) < 1e-6));
            if (i < 0) {
                continue;
            }
            matched += 1;
            let expected = caoRaw.frequencies.map((band) => band[i] * 33.35641).sort((a, b) => a - b);
            for (let n = 0; n < expected.length; n++) {
                assert.ok(Math.abs(expected[n] - mixed.eigenvalues[k][n]) < 1e-2, `${expected[n]} vs ${mixed.eigenvalues[k][n]}`);
            }
        }
        assert.ok(matched >= sampled.lineBreaks.length);

        let half = await computeMixedPhonon(sampled, 0.5);
        for (let k of gammaIndices(sampled)) {
            for (let n = 0; n < 3; n++) {
                assert.ok(Math.abs(half.eigenvalues[k][n]) < 0.05);
            }
        }
    });

    it('does not mix a rocksalt with a zinc blende', () => {
        let nacl = endpoint('phonondb-NaCl-mp-22862.json.gz', 'NaCl');
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        assert.throws(() => sampleMixingPair(nacl, gap), /Cannot mix/);
    });
});

describe('ternary rocksalt dynamical-matrix mixing', () => {
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

    it('matches cation-ordered cells onto their template', () => {
        let nains2 = endpoint('phonondb-NaInS2-mp-20289.json.gz', 'NaInS2');
        assert.equal(nains2.family, 'rocksalt-ABX2-R3m');
        assert.deepEqual(nains2.siteTypes, ['Na', 'In', 'S', 'S']);
        assert.ok(nains2.sampleAtQ);
        // charge balance: Tl+ next to Dy3+ takes the A site of Na+
        let dytls2 = endpoint('phonondb-DyTlS2-mp-31166.json.gz', 'DyTlS2');
        assert.equal(dytls2.family, 'rocksalt-ABX2-R3m');
        assert.deepEqual(dytls2.siteTypes, ['Tl', 'Dy', 'S', 'S']);
        // same composition and cell size as a chalcopyrite: the gamma-LiFeO2 ordering
        let lisco2 = endpoint('phonondb-LiScO2-mp-5840.json.gz', 'LiScO2');
        assert.equal(lisco2.family, 'rocksalt-ABX2-I41amd');
        assert.deepEqual(lisco2.siteTypes, ['Li', 'Li', 'Sc', 'Sc', 'O', 'O', 'O', 'O']);
        // the menu entry may name the family; the Materials Project file of the same compound
        let nains2MP = buildMixingEndpoint(loadFixture('mpdb-NaInS2-mp-20289-trimmed.json.gz'),
            { name: 'NaInS2', family: 'rocksalt-ABX2-R3m' });
        assert.equal(nains2MP.family, 'rocksalt-ABX2-R3m');
        assert.ok(!nains2MP.sampleAtQ);
        assert.throws(() => buildMixingEndpoint(loadFixture('phonondb-NaInS2-mp-20289.json.gz'),
            { name: 'NaInS2', family: 'rocksalt-ABX2-I41amd' }), /Could not match/);
    });

    it('labels Gamma only at reciprocal lattice vectors', () => {
        // the LiScO2 file has "X" at q = 0 and "G" off Gamma
        let lisco2 = endpoint('phonondb-LiScO2-mp-5840.json.gz', 'LiScO2');
        let isLatticeVector = (q) => q.every((v) => Math.abs(v - Math.round(v)) < 1e-5);
        for (let segment of lisco2.segments) {
            for (let [label, q] of [[segment.labelA, segment.qA], [segment.labelB, segment.qB]]) {
                assert.equal(label === 'GAMMA', isLatticeVector(q), label + ' ' + q);
            }
        }
    });

    it('folds a binary rocksalt into the ternary cell', async () => {
        // D of the 4-atom supercell at Q has the eigenvalues of NaCl at the 2 points Q + g
        let nacl = endpoint('phonondb-NaCl-mp-22862.json.gz', 'NaCl');
        let folded = binaryAsSuperstructure(nacl, 'rocksalt-ABX2-R3m');
        assert.deepEqual(folded.siteTypes, ['Na', 'Na', 'Cl', 'Cl']);

        let superReciprocal = utils.rec_lat(folded.lattice);
        let toCart = (q) => [0, 1, 2].map((d) => q[0] * superReciprocal[0][d] + q[1] * superReciprocal[1][d] + q[2] * superReciprocal[2][d]);
        let toPrimitive = (k) => nacl.lattice.map((p) => p[0] * k[0] + p[1] * k[1] + p[2] * k[2]);
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
                expected.push(...await eigen(nacl.evaluate(toPrimitive([0, 1, 2].map((d) => qCart[d] + g[d])), null)));
            }
            assert.equal(seen.size, 2);
            expected.sort((p, q) => p - q);
            let actual = await eigen(folded.sampleAtQ(Q, null));
            actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-9 * (1 + Math.abs(expected[i]))));
        }
    });

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        let nacl = endpoint('phonondb-NaCl-mp-22862.json.gz', 'NaCl');
        let nains2 = endpoint('phonondb-NaInS2-mp-20289.json.gz', 'NaInS2');
        let dytls2 = endpoint('phonondb-DyTlS2-mp-31166.json.gz', 'DyTlS2');
        for (let [a, b] of [[nacl, nains2], [nains2, dytls2]]) {
            let sampled = sampleMixingPair(a, b);
            let reversed = sampleMixingPair(b, a);
            assert.equal(sampled.recomputeNac, true);
            // the path does not depend on which material comes first
            assert.deepEqual(sampled.labels, reversed.labels);
            let end0 = await computeMixedPhonon(sampled, 0);
            assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-4);
            let end1 = await computeMixedPhonon(sampled, 1);
            assert.ok(opticalDiff(end1.eigenvalues, await computeReferenceEigenvalues(reversed)) < 1e-4);

            let mixed = await computeMixedPhonon(sampled, 0.5, 0.3);
            for (let k of gammaIndices(sampled)) {
                for (let n = 0; n < 3; n++) {
                    assert.ok(Math.abs(mixed.eigenvalues[k][n]) < 0.05, `mode ${n}: ${mixed.eigenvalues[k][n]}`);
                }
            }
            let mirrored = await computeMixedPhonon(reversed, 0.5, 0.7);
            assert.ok(opticalDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-6);
        }
    });

    it('lines up the paths of two Materials Project files', () => {
        // the two files run the rhombohedral path in orientations related by a
        // symmetry of the crystal; the canonical cells are chosen so they coincide
        let nains2 = buildMixingEndpoint(loadFixture('mpdb-NaInS2-mp-20289-trimmed.json.gz'), { name: 'NaInS2' });
        let klas2 = buildMixingEndpoint(loadFixture('mpdb-KLaS2-mp-15781-trimmed.json.gz'), { name: 'KLaS2' });
        let points = (e) => e.segments.reduce((total, segment) => total + segment.npoints, 0);
        let sampled = sampleMixingPair(nains2, klas2);
        assert.ok(sampled.qpoints.length >= 0.9 * Math.min(points(nains2), points(klas2)));
        for (let label of ['L', 'Z', 'F']) {
            assert.ok(sampled.labels.some((pair) => pair.includes(label)), label);
        }
    });

    it('refuses pairs on different lattices or without force constants', () => {
        let nains2 = endpoint('phonondb-NaInS2-mp-20289.json.gz', 'NaInS2');
        let lisco2 = endpoint('phonondb-LiScO2-mp-5840.json.gz', 'LiScO2');
        let zngen2 = endpoint('phonondb-ZnGeN2-mp-2979.json.gz', 'ZnGeN2');
        let gap = endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP');
        assert.throws(() => sampleMixingPair(nains2, lisco2), /Cannot mix/);
        assert.throws(() => sampleMixingPair(zngen2, nains2), /Cannot mix/);
        assert.throws(() => sampleMixingPair(gap, nains2), /Cannot mix/);
        let rows = JSON.parse(readFileSync(join(__dirname, 'fixtures', 'mixing', 'mp-dfpt-dielectric-CaO.json')));
        let cao = buildMixingEndpoint(loadFixture('mpdb-CaO-mp-2605-trimmed.json.gz'), { name: 'CaO' }, indexDielectricRows(rows).get(2605));
        assert.throws(() => sampleMixingPair(cao, nains2), /force constants/);
    });
});

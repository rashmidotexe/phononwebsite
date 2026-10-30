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
} from '../src/mixingphonons.js';
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

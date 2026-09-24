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
        assert.ok(maxAbsDiff(forward.eigenvalues, backward.eigenvalues) < 1e-8);
        assert.ok(maxAbsDiff([forward.distances], [backward.distances]) < 1e-12);
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
        assert.ok(maxAbsDiff(mixed.eigenvalues, reference) < 1e-8);
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

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
    getWorldMixing,
} from '../src/mixingphonons.js';
import { solveComplexHermitianWithEigenWasm } from '../src/eigenwasm.js';

// the body-centred cubic families: CsCl-type binaries, and Li3Bi-type binaries
// with the full and inverse Heuslers on their cell

const __dirname = dirname(fileURLToPath(import.meta.url));

function loadFixture(name) {
    return JSON.parse(gunzipSync(readFileSync(join(__dirname, 'fixtures', 'mixing', name))));
}

function endpoint(file, name) {
    return buildMixingEndpoint(loadFixture(file), { name: name });
}

const wrap = (v) => v - Math.floor(v + 1e-6);

describe('CsCl and Heusler mixing', () => {
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

    it('reads CsCl, Li3Bi-type and full Heusler cells', () => {
        let cscl = endpoint('phonondb-CsCl-mp-22865.json.gz', 'CsCl');
        assert.equal(cscl.family, 'cscl');
        assert.deepEqual(cscl.siteTypes, ['Cs', 'Cl']);
        let k3sb = endpoint('phonondb-K3Sb-mp-10159.json.gz', 'K3Sb');
        assert.equal(k3sb.family, 'heusler');
        assert.deepEqual(k3sb.siteTypes, ['K', 'K', 'K', 'Sb']);
        let k2cssb = endpoint('phonondb-K2CsSb-mp-581024.json.gz', 'K2CsSb');
        assert.equal(k2cssb.family, 'heusler-A2BX');
        assert.deepEqual(k2cssb.siteTypes, ['K', 'K', 'Cs', 'Sb']);
        // Cs on the octahedral site: (1/2,1/2,1/2) from Sb in the fcc cell
        let positions = k2cssb.structure.positionsRed;
        positions[2].forEach((v, k) => assert.ok(Math.abs(wrap(v - positions[3][k]) - 0.5) < 1e-6));
    });

    it('reads an inverse Heusler and intermetallics', () => {
        // K2CsSb with Cs and a K swapped: Cs on a tetrahedral site (Hg2CuTi type)
        let raw = loadFixture('phonondb-K2CsSb-mp-581024.json.gz');
        let types = raw.atom_types.slice();
        let cs = types.indexOf('Cs');
        let sb = types.indexOf('Sb');
        let tetrahedral = types.findIndex((type, i) => type === 'K' && raw.atom_pos_red[i].every((v, k) =>
            Math.abs(wrap(v - raw.atom_pos_red[sb][k]) - 0.25) < 1e-3 || Math.abs(wrap(v - raw.atom_pos_red[sb][k]) - 0.75) < 1e-3));
        [types[cs], types[tetrahedral]] = [types[tetrahedral], types[cs]];
        let inverse = analyzeSuperstructure('heusler', types, raw.atom_pos_red, raw.lattice);
        let b = inverse.sites.find((site) => site.role === 'B');
        assert.ok(b.sub === 2 || b.sub === 3);
        assert.equal(inverse.sites.find((site) => site.role === 'X').sub, 0);
        // NiAl (B2): no anion, the heavier Ni takes the X site
        let cubic = [[2.88, 0, 0], [0, 2.88, 0], [0, 0, 2.88]];
        let nial = analyzeSuperstructure('cscl', ['Ni', 'Al'], [[0, 0, 0], [0.5, 0.5, 0.5]], cubic);
        assert.deepEqual(nial.sites.map((site) => site.role), ['A', 'X']);
        // a rocksalt cell is not read as CsCl
        let nacl = analyzeSuperstructure('cscl', ['Na', 'Cl'], [[0, 0, 0], [0.5, 0.5, 0.5]],
            [[0, 2.8, 2.8], [2.8, 0, 2.8], [2.8, 2.8, 0]]);
        assert.equal(nacl, null);
    });

    it('folds a Li3Bi-type binary into the full Heusler cell', async () => {
        let k3sb = endpoint('phonondb-K3Sb-mp-10159.json.gz', 'K3Sb');
        let folded = binaryAsSuperstructure(k3sb, 'heusler-A2BX');
        for (let q of [[0.13, -0.07, 0.21], [0.5, 0, 0.5]]) {
            let expected = await eigen(k3sb.sampleAtQ(q, null));
            let actual = await eigen(folded.sampleAtQ(q, null));
            actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-8 * (1 + Math.abs(v))));
        }
    });

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        let cscl = endpoint('phonondb-CsCl-mp-22865.json.gz', 'CsCl');
        let csbr = endpoint('phonondb-CsBr-mp-22906.json.gz', 'CsBr');
        let k3sb = endpoint('phonondb-K3Sb-mp-10159.json.gz', 'K3Sb');
        let cs3sb = endpoint('phonondb-Cs3Sb-mp-10378.json.gz', 'Cs3Sb');
        let k2cssb = endpoint('phonondb-K2CsSb-mp-581024.json.gz', 'K2CsSb');
        for (let [a, b] of [[cscl, csbr], [k3sb, cs3sb], [k3sb, k2cssb], [k2cssb, cs3sb]]) {
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
        // a CsCl-type cell folds into the Li3Bi-type cell: the same four sites
        let union = sampleMixingPair(cscl, k3sb).union;
        assert.equal(union.map1.length, 4);
        assert.ok(union.map1.every((i) => i >= 0) && union.map2.every((i) => i >= 0));
    });
});

describe('mixing two structures on a common lattice (vacancies)', () => {
    let opticalDiff = (a, b) => {
        let worst = 0;
        a.forEach((row, k) => row.forEach((v, n) => {
            if (Math.abs(v) > 1) {
                worst = Math.max(worst, Math.abs(v - b[k][n]));
            }
        }));
        return worst;
    };
    let pairs = [
        ['phonondb-NaCl-mp-22862.json.gz', 'NaCl', 'phonondb-GaP-mp-2490.json.gz', 'GaP', 3],
        ['phonondb-NaCl-mp-22862.json.gz', 'NaCl', 'phonondb-CaF2-mp-2741.json.gz', 'CaF2', 3],
        ['phonondb-CaF2-mp-2741.json.gz', 'CaF2', 'phonondb-Li2O-mp-1960.json.gz', 'Li2O', 3],
        ['phonondb-GaP-mp-2490.json.gz', 'GaP', 'phonondb-K3Sb-mp-10159.json.gz', 'K3Sb', 4],
        ['phonondb-CsCl-mp-22865.json.gz', 'CsCl', 'phonondb-NaCl-mp-22862.json.gz', 'NaCl', 4],
        ['phonondb-SrTiO3-mp-5229.json.gz', 'SrTiO3', 'phonondb-Li3ClO-mp-985585.json.gz', 'Li3ClO', 5],
        ['phonondb-CsCl-mp-22865.json.gz', 'CsCl', 'phonondb-SrTiO3-mp-5229.json.gz', 'SrTiO3', 5],
    ];

    it('keeps end-members exact, acoustic modes at zero and mirror symmetry', async () => {
        for (let [file1, name1, file2, name2, sites] of pairs) {
            let label = name1 + ' + ' + name2;
            let sampled = sampleMixingPair(endpoint(file1, name1), endpoint(file2, name2));
            let reversed = sampleMixingPair(endpoint(file2, name2), endpoint(file1, name1));
            assert.equal(sampled.union.map1.length, sites, label);
            let end0 = await computeMixedPhonon(sampled, 0);
            // material 1 alone (a CsCl folded into the fcc cell shows that cell)
            assert.equal(end0.natoms, sampled.union.map1.filter((i) => i >= 0).length, label);
            assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-4, label);
            let end1 = await computeMixedPhonon(sampled, 1);
            assert.ok(opticalDiff(end1.eigenvalues, await computeReferenceEigenvalues(reversed)) < 1e-4, label);
            let mixed = await computeMixedPhonon(sampled, 0.4, 0.6);
            assert.equal(mixed.natoms, sites, label);
            sampled.qpoints.forEach((q, k) => {
                if (q.every((v) => Math.abs(v - Math.round(v)) < 1e-9)) {
                    let near = mixed.eigenvalues[k].map(Math.abs).sort((p, r) => p - r);
                    assert.ok(near[2] < 0.05, label + ': ' + near[2]);
                }
            });
            let mirrored = await computeMixedPhonon(reversed, 0.6, 0.4);
            assert.ok(opticalDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-6, label);
        }
    });

    it('adds the vacancy sites continuously: material 1 unchanged at small x', async () => {
        // a vacancy site couples to the rest as sqrt(x): a local mode near one of
        // NaCl's bands (the LO mode at 235 cm-1) splits it by about 1 cm-1 at x = 1e-4
        let sampled = sampleMixingPair(endpoint('phonondb-NaCl-mp-22862.json.gz', 'NaCl'),
            endpoint('phonondb-CaF2-mp-2741.json.gz', 'CaF2'));
        let reference = await computeReferenceEigenvalues(sampled);
        let small = await computeMixedPhonon(sampled, 1e-8);
        assert.equal(small.natoms, 3);
        reference.forEach((row, k) => row.forEach((v) => {
            assert.ok(small.eigenvalues[k].some((w) => Math.abs(w - v) < 0.1), `q ${k}: ${v}`);
        }));
    });

    it('tells which pairs mix and which material is folded', () => {
        let best = (f1, f2) => {
            let plan = getWorldMixing(f1, f2);
            let option = plan.options[0];
            return [plan.world, option.det, option.fold1, option.fold2];
        };
        assert.deepEqual(best('rocksalt', 'zincblende'), ['fcc', 1, false, false]);
        assert.deepEqual(best('perovskite', 'antiperovskite'), ['cubic', 1, false, false]);
        assert.deepEqual(best('cscl', 'perovskite'), ['cubic', 1, false, false]);
        assert.deepEqual(best('cscl', 'heusler'), ['fcc', 1, true, false]);
        assert.deepEqual(best('zincblende', 'rocksalt-ABX2-R3m'), ['fcc', 2, true, false]);
        // the same tilted cell along other axes: turned, relabelled, no fold
        assert.deepEqual(best('perovskite-Pnma', 'antiperovskite-Pnma'), ['cubic', 4, false, false]);
        // two orderings on different cells: the smallest common supercell, both folded
        let rocksalts = getWorldMixing('rocksalt-ABX2-R3m', 'rocksalt-ABX2-I41amd').options[0];
        assert.ok(rocksalts.fold1 && rocksalts.det % 2 === 0);
        // rutile (a buckled hcp anion array) and rocksalt on the close-packed grid
        assert.equal(getWorldMixing('rutile', 'rocksalt').world, 'hex');
        assert.equal(getWorldMixing('wurtzite', 'rocksalt').world, 'hex');
    });
});

describe('mixing orderings on different cells or along other axes', () => {
    let opticalDiff = (a, b) => {
        let worst = 0;
        a.forEach((row, k) => row.forEach((v, n) => {
            if (Math.abs(v) > 1) {
                worst = Math.max(worst, Math.abs(v - b[k][n]));
            }
        }));
        return worst;
    };
    let atZero = (sampled, mixed) => {
        let worst = 0;
        sampled.qpoints.forEach((q, k) => {
            if (q.every((v) => Math.abs(v) < 1e-9)) {
                worst = Math.max(worst, mixed.eigenvalues[k].map(Math.abs).sort((p, r) => p - r)[2]);
            }
        });
        return worst;
    };

    it('mixes two rocksalt orderings in their common supercell', async () => {
        // R-3m (2 cations) and I4_1/amd (4): both folded into the 16-site common cell,
        // the path of one of them carried into it
        let nains2 = () => endpoint('phonondb-NaInS2-mp-20289.json.gz', 'NaInS2');
        let lisco2 = () => endpoint('phonondb-LiScO2-mp-5840.json.gz', 'LiScO2');
        let sampled = sampleMixingPair(nains2(), lisco2());
        let reversed = sampleMixingPair(lisco2(), nains2());
        assert.equal(sampled.union.map1.length, 16);
        assert.ok(sampled.qpoints.length > 50);
        let end0 = await computeMixedPhonon(sampled, 0);
        assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-3);
        let mixed = await computeMixedPhonon(sampled, 0.4, 0.6);
        assert.ok(atZero(sampled, mixed) < 0.05);
        let mirrored = await computeMixedPhonon(reversed, 0.6, 0.4);
        assert.ok(opticalDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-6);
    });

    it('turns a tilted cell onto the axes of the other', () => {
        // perovskite-Pnma and antiperovskite-Pnma: the same cell along other axes
        let plan = getWorldMixing('perovskite-Pnma', 'antiperovskite-Pnma');
        let option = plan.options[0];
        assert.equal(option.det, 4);
        assert.ok(!option.fold1 && !option.fold2);
        assert.ok(option.index > 0);
    });
});

describe('close-packed structures on one grid (wurtzite, zinc blende, rocksalt, ...)', () => {
    let opticalDiff = (a, b) => {
        let worst = 0;
        a.forEach((row, k) => row.forEach((v, n) => {
            if (Math.abs(v) > 1) {
                worst = Math.max(worst, Math.abs(v - b[k][n]));
            }
        }));
        return worst;
    };
    let zmo = () => buildMixingEndpoint(JSON.parse(gunzipSync(readFileSync(join(__dirname, '..', 'data', 'mixing', 'zumba', 'ZMO.json.gz')))),
        { name: 'Zn2Mo3O8', family: 'zmo' });

    it('reads the Zumba Zn2Mo3O8 force constants as matdyn does', async () => {
        // the q2r force constants with QE's Ewald parameter: the matdyn dispersion of
        // the Zumba page, away from Gamma (where the LO direction is a choice)
        let sampled = sampleMixingPair(zmo(), zmo());
        let reference = await computeReferenceEigenvalues(sampled);
        let matdyn = JSON.parse(readFileSync(join(__dirname, '..', 'zumbadb', 'ZMO_dispersion.json')));
        let worst = 0;
        sampled.qpoints.forEach((q, k) => {
            let j = matdyn.qpoints.findIndex((p) => p.every((v, i) => Math.abs(v - q[i]) < 1e-6));
            if (j < 0 || q.every((v) => Math.abs(v) < 1e-9)) {
                return;
            }
            let ours = reference[k].slice().sort((a, b) => a - b);
            matdyn.eigenvalues[j].slice().sort((a, b) => a - b).forEach((v, i) => { worst = Math.max(worst, Math.abs(v - ours[i])); });
        });
        assert.ok(worst < 0.01, String(worst));
        let half = await computeMixedPhonon(sampled, 0);
        assert.ok(opticalDiff(half.eigenvalues, reference) < 1e-6);
    });

    it('mixes wurtzite with zinc blende, rocksalt, corundum, 1T and perovskite', async () => {
        let aln = () => endpoint('phonondb-AlN-wurtzite-mp-661.json.gz', 'AlN');
        let pairs = [
            [aln, () => endpoint('phonondb-GaP-mp-2490.json.gz', 'GaP'), 20],
            [aln, () => endpoint('phonondb-NaCl-mp-22862.json.gz', 'NaCl'), 22],
            [aln, () => endpoint('phonondb-Al2O3-mp-1143.json.gz', 'Al2O3'), 16],
            [() => endpoint('phonondb-ZrS2-1T-mp-1186.json.gz', 'ZrS2'), aln, 5],
            [() => endpoint('phonondb-SrTiO3-mp-5229.json.gz', 'SrTiO3'), () => endpoint('phonondb-NaCl-mp-22862.json.gz', 'NaCl'), 8],
        ];
        for (let [first, second, sites] of pairs) {
            let sampled = sampleMixingPair(first(), second());
            let reversed = sampleMixingPair(second(), first());
            let label = sampled.endpoint1.name + ' + ' + sampled.endpoint2.name;
            assert.equal(sampled.union.world, 'hex', label);
            assert.equal(sampled.union.map1.length, sites, label);
            let end0 = await computeMixedPhonon(sampled, 0);
            assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-4, label);
            let end1 = await computeMixedPhonon(sampled, 1);
            assert.ok(opticalDiff(end1.eigenvalues, await computeReferenceEigenvalues(reversed)) < 1e-4, label);
            let mixed = await computeMixedPhonon(sampled, 0.4, 0.6);
            let mirrored = await computeMixedPhonon(reversed, 0.6, 0.4);
            assert.ok(opticalDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-4, label);
            sampled.qpoints.forEach((q, k) => {
                if (q.every((v) => Math.abs(v) < 1e-9)) {
                    assert.ok(mixed.eigenvalues[k].map(Math.abs).sort((p, r) => p - r)[2] < 0.05, label);
                }
            });
        }
    });

    it('mixes rutile and trirutile with 1T, wurtzite and corundum', async () => {
        // rutile: the hcp anions of 1T / wurtzite, the cations in half the octahedral
        // holes; one cation of the 2-cell 1T moves (CdI2 fills alternate layers fully,
        // rutile half of every layer)
        let tio2 = () => endpoint('phonondb-TiO2-mp-2657.json.gz', 'TiO2');
        let pairs = [
            [tio2, () => endpoint('phonondb-ZrS2-1T-mp-1186.json.gz', 'ZrS2'), 7],
            [tio2, () => endpoint('phonondb-AlN-wurtzite-mp-661.json.gz', 'AlN'), 10],
            [() => endpoint('phonondb-ZnSb2O6-mp-3188.json.gz', 'ZnSb2O6'), () => endpoint('phonondb-ZrS2-1T-mp-1186.json.gz', 'ZrS2'), null],
        ];
        for (let [first, second, sites] of pairs) {
            let sampled = sampleMixingPair(first(), second());
            let reversed = sampleMixingPair(second(), first());
            let label = sampled.endpoint1.name + ' + ' + sampled.endpoint2.name;
            assert.equal(sampled.union.world, 'hex', label);
            if (sites) {
                assert.equal(sampled.union.map1.length, sites, label);
            }
            let end0 = await computeMixedPhonon(sampled, 0);
            assert.ok(opticalDiff(end0.eigenvalues, await computeReferenceEigenvalues(sampled)) < 1e-4, label);
            let end1 = await computeMixedPhonon(sampled, 1);
            assert.ok(opticalDiff(end1.eigenvalues, await computeReferenceEigenvalues(reversed)) < 1e-4, label);
            let mixed = await computeMixedPhonon(sampled, 0.4, 0.6);
            let mirrored = await computeMixedPhonon(reversed, 0.6, 0.4);
            assert.ok(opticalDiff(mixed.eigenvalues, mirrored.eigenvalues) < 1e-4, label);
        }
    });

    it('folds a wurtzite into the Zn2Mo3O8 cell', async () => {
        let sampled = sampleMixingPair(endpoint('phonondb-AlN-wurtzite-mp-661.json.gz', 'AlN'), zmo());
        assert.equal(sampled.union.world, 'hex');
        // AlN folded 8 times (2 x 2 in the layer, 4 layers): 32 sites, 13 on Zn2Mo3O8 sites
        assert.equal(sampled.union.map1.length, 45);
        let end1 = await computeMixedPhonon(sampled, 1);
        let alone = await computeReferenceEigenvalues(sampleMixingPair(zmo(), zmo()));
        assert.ok(opticalDiff(end1.eigenvalues, alone) < 1e-4);
    });

    it('caps the size of the common cell', () => {
        // zinc blende (3 layers) with Zn2Mo3O8 (4 layers, 2 x 2): 12 layers, too many sites
        assert.equal(getWorldMixing('zincblende', 'zmo'), null);
        assert.equal(getWorldMixing('wurtzite', 'zmo').world, 'hex');
    });
});

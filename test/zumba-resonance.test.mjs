import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

import {
  buildShapes,
  computeNormalisers,
  dominantSegments,
  efCell,
  evToNm,
  formatEF,
  nmToEV,
  normalise,
  parseEnergyInput,
  snapIndex,
  spectrumAt,
  wavelengthColor,
  xGrid,
} from '../src/zumbaresonance.js';

const energies = Array.from({ length: 526 }, (_, i) => 2.3 + i * 0.7 / 525);

describe('resonance energy helpers', () => {
  it('snaps to the grid and clamps to the window', () => {
    assert.equal(snapIndex(energies, 2.3), 0);
    assert.equal(snapIndex(energies, 3.0), 525);
    assert.equal(snapIndex(energies, 2.0), 0);
    assert.equal(snapIndex(energies, 3.5), 525);
    assert.equal(snapIndex(energies, 2.42), 90);
    for (const i of [0, 1, 90, 300, 524, 525]) assert.equal(snapIndex(energies, energies[i]), i);
  });

  it('converts between eV and nm', () => {
    assert.ok(Math.abs(evToNm(2.331) - 531.9) < 0.1);
    assert.ok(Math.abs(nmToEV(488) - 2.5407) < 1e-3);
    assert.ok(Math.abs(nmToEV(evToNm(2.7)) - 2.7) < 1e-12);
  });

  it('reads energies typed in eV or nm', () => {
    assert.equal(parseEnergyInput('2.54'), 2.54);
    assert.equal(parseEnergyInput(' 2.54 eV '), 2.54);
    assert.ok(Math.abs(parseEnergyInput('488nm') - nmToEV(488)) < 1e-12);
    assert.ok(Math.abs(parseEnergyInput('488 nm') - nmToEV(488)) < 1e-12);
    assert.ok(Math.abs(parseEnergyInput('488') - nmToEV(488)) < 1e-12);
    assert.equal(parseEnergyInput('abc'), null);
    assert.equal(parseEnergyInput(''), null);
    assert.equal(parseEnergyInput('0'), null);
    // out-of-range values are returned as typed; the view clamps them
    assert.equal(snapIndex(energies, parseEnergyInput('5')), 525);
  });

  it('colours energies like the laser wavelength', () => {
    for (let nm = 413; nm <= 540; nm += 1) assert.match(wavelengthColor(nm), /^rgb\(\d{1,3},\d{1,3},\d{1,3}\)$/);
    // same values as wavelength_rgb() of ref/analysis/plot_excph_raman.py, x 0.75
    assert.equal(wavelengthColor(532), 'rgb(60,191,0)');
    assert.equal(wavelengthColor(488), 'rgb(0,184,191)');
    assert.equal(wavelengthColor(1239.84 / 3.0), 'rgb(75,0,169)');
  });
});

describe('resonance normalisation', () => {
  // three peaks on four laser energies
  const intensity = [
    [1, 4, 2, 0.5],
    [2, 1, 1, 0.25],
    [0.5, 0.5, 8, 0.1],
  ];
  const norm = computeNormalisers(intensity);

  it('finds the strongest resonance, the strongest peak per energy and each peak maximum', () => {
    assert.equal(norm.M, 8);
    assert.equal(norm.mPeak, 2);
    assert.equal(norm.mIndex, 2);
    assert.deepEqual(Array.from(norm.s), [2, 4, 8, 0.5]);
    assert.deepEqual(Array.from(norm.P), [4, 2, 8]);
    assert.deepEqual(Array.from(norm.Pidx), [1, 0, 2]);
  });

  it("'strongest' puts the whole window on one scale", () => {
    const all = intensity.flatMap((row, k) => row.map((_, i) => normalise(intensity, norm, k, i, 'strongest')));
    assert.equal(Math.max(...all), 1);
    assert.equal(normalise(intensity, norm, 0, 3, 'strongest'), 0.5 / 8);
  });

  it("'each' scales every spectrum to its own strongest peak", () => {
    for (let i = 0; i < 4; i++) {
      const col = intensity.map((_, k) => normalise(intensity, norm, k, i, 'each'));
      assert.equal(Math.max(...col), 1);
    }
  });

  it('merges runs of the dominant peak', () => {
    assert.deepEqual(dominantSegments(intensity), [
      { k: 1, i0: 0, i1: 0 },
      { k: 0, i0: 1, i1: 1 },
      { k: 2, i0: 2, i1: 2 },
      { k: 0, i0: 3, i1: 3 },
    ]);
    assert.deepEqual(dominantSegments([[1, 1, 1], [0, 0, 2]]), [{ k: 0, i0: 0, i1: 1 }, { k: 1, i0: 2, i1: 2 }]);
  });

  it('builds height-1 Lorentzians whose sum is the broadened spectrum', () => {
    const xs = xGrid([100, 104.3], 0.5, 200);
    assert.ok(xs.includes(104.3));
    assert.ok(xs.every((x, j) => j === 0 || x > xs[j - 1]));
    const shapes = buildShapes(xs, [100, 104.3], 2);
    const at = (x) => xs.indexOf(x);
    assert.equal(shapes[0][at(100)], 1);
    assert.equal(shapes[0][at(102)], 0.5);           // HWHM = 2
    const y = spectrumAt(shapes, [0.5, 0.25]);
    assert.ok(Math.abs(y[at(100)] - (0.5 + 0.25 * 4 / (4.3 ** 2 + 4))) < 1e-12);
  });
});

describe('resonance enhancement factor', () => {
  it('is the ratio to the reference energy', () => {
    assert.deepEqual(efCell(930, 1, 1, true), { value: 930, state: 'ok' });
  });

  it("flags peaks too weak at the reference as 'res.-only'", () => {
    assert.equal(efCell(5, 0.0099, 1, true).state, 'resonly');
    assert.equal(efCell(5, 0.01, 1, true).state, 'ok');
    assert.equal(efCell(5, 0, 1, true).state, 'resonly');
  });

  it("greys out peaks that are barely visible at E_L", () => {
    assert.deepEqual(efCell(2, 1, 1, false), { value: 2, state: 'weak' });
  });

  it('formats with two significant figures', () => {
    assert.equal(formatEF(2.4), '×2.4');
    assert.equal(formatEF(0.414), '×0.41');
    assert.equal(formatEF(930.4), '×930');
    assert.equal(formatEF(3843), '×3.8·10³');
    assert.equal(formatEF(9990), '×1.0·10⁴');
    assert.equal(formatEF(1.2e5), '×1.2·10⁵');
    assert.equal(formatEF(NaN), '—');
  });
});

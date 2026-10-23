import assert from 'node:assert/strict';
import { afterEach, beforeEach, describe, it } from 'node:test';

import {
  setupLegacyTestEnv,
  teardownLegacyTestEnv,
} from './helpers/legacy-test-env.mjs';

import { fetchRruffSpectra, getRruffSeries } from '../src/rruff.js';

const PAYLOAD = {
  material: 'mp-2657',
  reference: { text: 'RRUFF project', url: 'https://rruff.info/', doi: '10.1515/9783110417104-003' },
  spectra: [
    { name: 'Rutile', rruff_id: 'R040049', url: 'https://rruff.example/R040049', laser_nm: 532,
      start: 100, step: 2, intensity: [0, 500, 1000, 250] },
    { name: 'Rutile', rruff_id: 'R050031', url: 'https://rruff.example/R050031', laser_nm: null,
      start: 100, step: 2, intensity: [1000, 0] },
  ],
};

function mockFetch(files) {
  const calls = [];
  globalThis.fetch = async (url) => {
    calls.push(url);
    if (!(url in files)) { return new Response('', { status: 404 }); }
    return new Response(JSON.stringify(files[url]), { status: 200 });
  };
  return calls;
}

describe('rruff', () => {
  let originalFetch;

  beforeEach(() => {
    originalFetch = globalThis.fetch;
  });

  afterEach(() => {
    globalThis.fetch = originalFetch;
  });

  it('only requests materials listed in the index, once each', async () => {
    const calls = mockFetch({
      'base1/index.json': { materials: ['mp-2657'] },
      'base1/mp-2657.json': PAYLOAD,
    });

    assert.deepEqual(await fetchRruffSpectra('mp-2657', 'base1/'), PAYLOAD);
    assert.deepEqual(await fetchRruffSpectra('mp-2657', 'base1/'), PAYLOAD);
    assert.equal(await fetchRruffSpectra('mp-149', 'base1/'), null);
    assert.deepEqual(calls, ['base1/index.json', 'base1/mp-2657.json']);
  });

  it('returns nothing without an index or an id', async () => {
    const calls = mockFetch({});
    assert.equal(await fetchRruffSpectra('mp-2657', 'base2/'), null);
    assert.equal(await fetchRruffSpectra(null, 'base2/'), null);
    assert.deepEqual(calls, ['base2/index.json']);
  });

  it('scales the spectra to the computed peak and shows only the first', () => {
    const series = getRruffSeries(PAYLOAD, 2.0, 104);

    assert.equal(series.length, 2);
    assert.deepEqual(series[0].data, [[100, 0], [102, 1], [104, 2]]);
    assert.equal(series[0].name, 'Rutile R040049 · 532 nm');
    assert.equal(series[1].name, 'Rutile R050031');
    assert.equal(series[0].visible, true);
    assert.equal(series[1].visible, false);
    assert.deepEqual(getRruffSeries(null, 1, 100), []);
  });
});

describe('PhononWebpage.loadExperimentalRaman', () => {
  let dom;
  let PhononWebpage;
  let originalFetch;

  beforeEach(async () => {
    ({ dom } = setupLegacyTestEnv());
    ({ PhononWebpage } = await import('../src/phononwebpage.js'));
    originalFetch = globalThis.fetch;
    mockFetch({
      'page/index.json': { materials: ['mp-2657'] },
      'page/mp-2657.json': PAYLOAD,
    });
  });

  afterEach(() => {
    globalThis.fetch = originalFetch;
    teardownLegacyTestEnv(dom);
  });

  function makePage(phonon) {
    return {
      phonon: phonon,
      ramanDbId: 'mp-2657',
      rruffDataUrl: 'page/',
      plotted: 0,
      plotRaman() { this.plotted++; },
    };
  }

  it('redraws the Raman spectrum with the measured data', async () => {
    const page = makePage({ raman_intensities: [0, 1] });
    await PhononWebpage.prototype.loadExperimentalRaman.call(page);

    assert.deepEqual(page.phonon.experimental_raman, PAYLOAD);
    assert.equal(page.plotted, 1);
  });

  it('keeps the data for later when there is no computed spectrum yet', async () => {
    const page = makePage({ raman_intensities: null });
    await PhononWebpage.prototype.loadExperimentalRaman.call(page);

    assert.deepEqual(page.phonon.experimental_raman, PAYLOAD);
    assert.equal(page.plotted, 0);
  });

  it('ignores spectra of a material that is no longer shown', async () => {
    const first = { raman_intensities: [0, 1] };
    const page = makePage(first);
    const pending = PhononWebpage.prototype.loadExperimentalRaman.call(page);
    page.phonon = { raman_intensities: [0, 1] };
    await pending;

    assert.equal(first.experimental_raman, undefined);
    assert.equal(page.phonon.experimental_raman, undefined);
    assert.equal(page.plotted, 0);
  });
});

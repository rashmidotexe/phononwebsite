import assert from 'node:assert/strict';
import { afterEach, beforeEach, describe, it } from 'node:test';

import {
  setupLegacyTestEnv,
  teardownLegacyTestEnv,
} from './helpers/legacy-test-env.mjs';

describe('PhononWebpage plot view', () => {
  let dom;
  let PhononWebpage;

  beforeEach(async () => {
    ({ dom } = setupLegacyTestEnv());
    document.body.insertAdjacentHTML('beforeend', `
      <div class="flex-highcharts">
        <div id="plot-view-controls" style="display:none">
          <button data-plot-view="dispersion"></button>
          <button data-plot-view="raman"></button>
        </div>
        <div id="raman-spectrum"><div id="raman-spectrum-chart"></div><div id="raman-table-container"></div></div>
      </div>`);
    ({ PhononWebpage } = await import('../src/phononwebpage.js'));
  });

  afterEach(() => {
    teardownLegacyTestEnv(dom);
  });

  function makePage() {
    const page = Object.create(PhononWebpage.prototype);
    page.activePlotView = null;
    page.defaultPlotView = 'dispersion';
    page.phonon = {
      eigenvalues: [[0, 0, 0, 320, 690, 700]],
      raman_intensities: [0, 0, 0, 0.1, 1, 0.5],
      gamma_index: 0,
    };
    return page;
  }

  function selectedButton() {
    return document.querySelector('[data-plot-view][aria-selected="true"]').dataset.plotView;
  }

  it('opens materials with Raman data on the dispersion', () => {
    const page = makePage();
    page.plotRaman();

    assert.equal(page.activePlotView, 'dispersion');
    assert.equal(selectedButton(), 'dispersion');
    assert.equal(document.getElementById('plot-view-controls').style.display, 'block');
    assert.equal(document.getElementById('highcharts').style.display, '');
    assert.equal(document.getElementById('raman-spectrum').style.display, 'none');
  });

  it('keeps the view the user picked while it stays available', () => {
    const page = makePage();
    page.plotRaman();
    page.setPlotView('raman');
    page.plotRaman();

    assert.equal(page.activePlotView, 'raman');
    assert.equal(selectedButton(), 'raman');
  });

  it('offers the Raman view for materials with only measured spectra', () => {
    const page = makePage();
    page.phonon.raman_intensities = null;
    page.phonon.kpoints = [[0.5, 0, 0], [0, 0, 0]];
    page.phonon.eigenvalues = [[1, 2, 3, 4, 5, 6], [0, 0, 0, 320, 690, 700]];
    page.phonon.experimental_raman = {
      reference: { doi: '10.1515/9783110417104-003' },
      spectra: [{ name: 'Rutile', rruff_id: 'R040049', start: 100, step: 2, intensity: [0, 1000, 0] }],
    };
    page.plotRaman();

    assert.equal(document.getElementById('plot-view-controls').style.display, 'block');
    assert.ok(document.getElementById('raman-measured-only'));
    //the Gamma modes are still listed so they can be animated
    assert.match(document.getElementById('raman-table').textContent, /690\.00/);
  });

  it('hides the Raman view when there is no Raman data at all', () => {
    const page = makePage();
    page.phonon.raman_intensities = null;
    page.plotRaman();

    assert.equal(document.getElementById('plot-view-controls').style.display, 'none');
  });

  it('lets a page open on the Raman spectrum instead', () => {
    const page = makePage();
    page.defaultPlotView = 'raman';
    page.plotRaman();

    assert.equal(selectedButton(), 'raman');
  });
});

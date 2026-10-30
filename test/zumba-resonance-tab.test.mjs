import assert from 'node:assert/strict';
import { afterEach, beforeEach, describe, it } from 'node:test';

import {
  setupLegacyTestEnv,
  teardownLegacyTestEnv,
} from './helpers/legacy-test-env.mjs';

describe('Zumba resonance tab', () => {
  let dom;
  let PhononWebpage;
  let installResonanceTab;

  beforeEach(async () => {
    ({ dom } = setupLegacyTestEnv());
    window.requestAnimationFrame = (callback) => callback();
    document.getElementById('highcharts').remove();
    document.body.insertAdjacentHTML('beforeend', `
      <div class="flex-container">
        <div class="flex-menu">
          <div id="resonance-menu" style="display:none"></div>
          <div id="panel1"></div>
          <div id="panel2"></div>
        </div>
        <div class="flex-highcharts">
          <div id="plot-view-controls" style="display:none">
            <button data-plot-view="dispersion"></button>
            <button data-plot-view="raman"></button>
            <button data-plot-view="resonance"></button>
          </div>
          <div id="plot-chart-area">
            <div id="highcharts"></div>
            <div id="raman-spectrum"><div id="raman-spectrum-chart"></div><div id="raman-table-container"></div></div>
            <div id="resonance-view"></div>
          </div>
        </div>
      </div>`);
    ({ PhononWebpage } = await import('../src/phononwebpage.js'));
    ({ installResonanceTab } = await import('../src/zumbaresonancetab.js'));
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
    const calls = [];
    const view = {
      activate: () => calls.push('activate'),
      deactivate: () => calls.push('deactivate'),
      reflow: () => calls.push('reflow'),
    };
    installResonanceTab(page, view);
    return { page, calls };
  }

  const display = (id) => document.getElementById(id).style.display;
  const selected = () => document.querySelector('[data-plot-view][aria-selected="true"]').dataset.plotView;

  it('shows the resonance view and its menu in place of the geometry panels', () => {
    const { page, calls } = makePage();
    page.plotRaman();
    page.setPlotView('resonance');

    assert.equal(page.activePlotView, 'resonance');
    assert.equal(selected(), 'resonance');
    assert.equal(display('resonance-view'), 'block');
    assert.equal(display('highcharts'), 'none');
    assert.equal(display('raman-spectrum'), 'none');
    assert.equal(display('resonance-menu'), '');
    assert.equal(display('panel1'), 'none');
    assert.equal(display('panel2'), 'none');
    assert.ok(document.querySelector('.flex-container').classList.contains('resonance-active'));
    assert.deepEqual(calls, ['deactivate', 'activate', 'reflow']);
  });

  it('restores the Raman and dispersion views and the geometry panels', () => {
    const { page, calls } = makePage();
    page.plotRaman();
    page.setPlotView('resonance');
    page.setPlotView('raman');

    assert.equal(page.activePlotView, 'raman');
    assert.equal(selected(), 'raman');
    assert.equal(display('resonance-view'), 'none');
    assert.equal(display('raman-spectrum'), 'flex');
    assert.equal(display('resonance-menu'), 'none');
    assert.equal(display('panel1'), '');
    assert.equal(display('panel2'), '');
    assert.ok(!document.querySelector('.flex-container').classList.contains('resonance-active'));

    page.setPlotView('dispersion');
    assert.equal(selected(), 'dispersion');
    assert.equal(display('highcharts'), '');
    assert.equal(display('raman-spectrum'), 'none');
    assert.deepEqual(calls, ['deactivate', 'activate', 'reflow', 'deactivate', 'deactivate']);
  });

  it('keeps the resonance view when a geometry reload replots the Raman data', () => {
    const { page } = makePage();
    page.plotRaman();
    page.setPlotView('resonance');
    page.plotRaman();

    assert.equal(page.activePlotView, 'resonance');
    assert.equal(selected(), 'resonance');
    assert.equal(display('resonance-view'), 'block');
    assert.equal(display('raman-spectrum'), 'none');
    assert.equal(display('highcharts'), 'none');
  });

  it('opens on the dispersion as before', () => {
    const { page, calls } = makePage();
    page.plotRaman();

    assert.equal(page.activePlotView, 'dispersion');
    assert.equal(display('resonance-view'), 'none');
    assert.equal(display('panel1'), '');
    assert.deepEqual(calls, ['deactivate']);
  });
});

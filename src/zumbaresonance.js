/*
 * zumbaresonance.js  —  "Resonance" view of zumba.html
 *
 * Exciton-phonon resonant Raman of Zn2Mo3O8 against the excitation energy
 * E_L (zumbadb/resonance/ZMO_resonance.json, written by
 * python/phononweb/scripts/generate_zumba_resonance.py):
 *   - a spectrum at the slider energy, with pinned overlays, or a waterfall
 *   - resonance profiles of every peak over the GW-BSE absorption
 *   - a peak table with enhancement factors
 *
 * Every plotted value derives from the stick intensities I_m(E) of the
 * selected geometry; see normalise() for the two normalisations.
 */

import Highcharts from 'highcharts';

export const HC_EV_NM = 1239.84;
export const LASER_LINES = [532, 514.5, 488, 473, 457.9, 441.6];   // nm
const PALETTE = ['#D55E00', '#E69F00', '#CC79A7', '#8C510A', '#B2182B', '#F781BF'];
const LIVE_COLOR = '#2c3e50';
const FAINT_COLOR = '#bdc3c7';
const GAMMA = 2.0;            // Lorentzian HWHM (cm-1), same as the Raman tab
const X_MAX = 850;
const MAX_PINS = 8;
const MAX_SELECTED = 6;
const RESONLY = 0.01;         // EF denominator floor (fraction of the strongest peak at E_ref)
const WEAK = 0.03;            // visibility cut at E_L (the cut of 05_raman_gw_selectivity)
const LOG_FLOOR = 1e-4;
const PLAY_SECONDS = 9;
const SYMBOLS = { A_1: 'circle', E_1: 'triangle', E_2: 'diamond' };
const SUB = { A_1: 'A₁', E_1: 'E₁', E_2: 'E₂' };
const EQUIVALENT = { 'x(yy)x': 'y(xx)y', 'y(xx)y': 'x(yy)x', 'x(yz)x': 'y(xz)y',
                     'y(xz)y': 'x(yz)x', 'x(zz)x': 'y(zz)y', 'y(zz)y': 'x(zz)x' };

/* ---------- pure helpers (no DOM, unit tested) ---------- */

export function evToNm(ev) { return HC_EV_NM / ev; }
export function nmToEV(nm) { return HC_EV_NM / nm; }

export function snapIndex(energies, ev) {
    if (!(ev > energies[0])) return 0;
    const n = energies.length;
    if (ev >= energies[n - 1]) return n - 1;
    const step = (energies[n - 1] - energies[0]) / (n - 1);
    return Math.max(0, Math.min(n - 1, Math.round((ev - energies[0]) / step)));
}

export function indexToEV(energies, i) { return energies[i]; }

// '2.54', '2.54 eV', '488', '488nm' -> energy in eV (unclamped), or null
export function parseEnergyInput(text) {
    const m = String(text).trim().toLowerCase().match(/^([0-9]*\.?[0-9]+)\s*(ev|nm)?$/);
    if (!m) return null;
    const v = parseFloat(m[1]);
    if (!(v > 0)) return null;
    const nm = m[2] === 'nm' || (m[2] !== 'ev' && v > 100);
    return nm ? nmToEV(v) : v;
}

// wavelength_rgb() of ref/analysis/plot_excph_raman.py, darkened for a white background
export function wavelengthColor(nm, darken = 0.75) {
    let r, g, b;
    if (nm < 380) { r = 0.35; g = 0; b = 0.45; }
    else if (nm > 780) { r = 0.45; g = 0; b = 0; }
    else {
        if (nm < 440)      { r = -(nm - 440) / 60; g = 0; b = 1; }
        else if (nm < 490) { r = 0; g = (nm - 440) / 50; b = 1; }
        else if (nm < 510) { r = 0; g = 1; b = -(nm - 510) / 20; }
        else if (nm < 580) { r = (nm - 510) / 70; g = 1; b = 0; }
        else if (nm < 645) { r = 1; g = -(nm - 645) / 65; b = 0; }
        else               { r = 1; g = 0; b = 0; }
        const f = nm < 420 ? 0.3 + 0.7 * (nm - 380) / 40 : (nm > 700 ? 0.3 + 0.7 * (780 - nm) / 80 : 1);
        r *= f; g *= f; b *= f;
    }
    const c = (v) => Math.round(255 * Math.max(0, Math.min(1, v * darken)));
    return `rgb(${c(r)},${c(g)},${c(b)})`;
}

const SUPERSCRIPT = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };

export function formatEF(v) {
    if (!Number.isFinite(v)) return '—';
    if (v >= 1000) {
        const e = Math.floor(Math.log10(v));
        let mant = (v / 10 ** e).toFixed(1);
        let exp = e;
        if (mant === '10.0') { mant = '1.0'; exp += 1; }
        return `×${mant}·10${String(exp).split('').map((ch) => SUPERSCRIPT[ch]).join('')}`;
    }
    if (v >= 10) return `×${Math.round(v)}`;
    return `×${Number(v.toPrecision(2))}`;
}

// EF of one peak; 'resonly' when its reference intensity is too small for
// the ratio to mean anything, 'weak' when it is barely visible at E_L
export function efCell(iNow, iRef, maxAtRef, visibleNow) {
    if (!(iRef >= RESONLY * maxAtRef) || !(iRef > 0)) return { value: NaN, state: 'resonly' };
    return { value: iNow / iRef, state: visibleNow ? 'ok' : 'weak' };
}

// per-geometry constants of the normalisations
export function computeNormalisers(intensity) {
    const nPeaks = intensity.length;
    const nE = intensity[0].length;
    const s = new Float64Array(nE);
    const P = new Float64Array(nPeaks);
    const Pidx = new Int32Array(nPeaks);
    let M = 0, mPeak = 0, mIndex = 0;
    for (let k = 0; k < nPeaks; k++) {
        const row = intensity[k];
        for (let i = 0; i < nE; i++) {
            const v = row[i];
            if (v > s[i]) s[i] = v;
            if (v > P[k]) { P[k] = v; Pidx[k] = i; }
            if (v > M) { M = v; mPeak = k; mIndex = i; }
        }
    }
    return { M, mPeak, mIndex, s, P, Pidx };
}

// n_m(E): 'strongest' -> I/M (one scale for the window), 'each' -> I/s(E)
export function normalise(intensity, norm, k, i, mode) {
    const v = intensity[k][i];
    return mode === 'each' ? (norm.s[i] > 0 ? v / norm.s[i] : 0) : v / norm.M;
}

// height-1 Lorentzians of every peak on the x grid
export function buildShapes(xs, freqs, gamma = GAMMA) {
    const g2 = gamma * gamma;
    return freqs.map((f) => {
        const out = new Float64Array(xs.length);
        for (let j = 0; j < xs.length; j++) {
            const d = xs[j] - f;
            out[j] = g2 / (d * d + g2);
        }
        return out;
    });
}

export function spectrumAt(shapes, weights) {
    const out = new Float64Array(shapes.length ? shapes[0].length : 0);
    for (let k = 0; k < shapes.length; k++) {
        const w = weights[k];
        if (!w) continue;
        const sh = shapes[k];
        for (let j = 0; j < out.length; j++) out[j] += w * sh[j];
    }
    return out;
}

// runs of the strongest peak along E: [{k, i0, i1}]
export function dominantSegments(intensity) {
    const nE = intensity[0].length;
    const segs = [];
    for (let i = 0; i < nE; i++) {
        let best = 0;
        for (let k = 1; k < intensity.length; k++) if (intensity[k][i] > intensity[best][i]) best = k;
        const last = segs[segs.length - 1];
        if (last && last.k === best) last.i1 = i;
        else segs.push({ k: best, i0: i, i1: i });
    }
    return segs;
}

export function xGrid(freqs, step = 0.5, xmax = X_MAX) {
    const xs = [];
    for (let x = 0; x <= xmax + 1e-9; x += step) xs.push(Math.round(x * 1000) / 1000);
    for (const f of freqs) if (f >= 0 && f <= xmax) xs.push(f);
    return Float64Array.from(new Set(xs)).sort();
}

export function symLabel(sym) { return SUB[sym] || sym; }

function fmtI(v) {
    if (!(v > 0)) return '0';
    return v >= 0.1 ? v.toFixed(2) : Number(v.toPrecision(2)).toString();
}

function fmtEV(ev) { return `${ev.toFixed(3)} eV (${evToNm(ev).toFixed(1)} nm)`; }

function escapeHTML(text) {
    return String(text).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
}

/* ---------- the view ---------- */

export class ResonanceView {

    constructor({ view, menu, dataUrl, ramanFileFor }) {
        this.view = view;
        this.menu = menu;
        this.dataUrl = dataUrl;
        this.ramanFileFor = ramanFileFor;

        this.loading = null;
        this.built = false;
        this.charts = {};

        this.geom = 0;
        this.ie = 0;
        this.pins = [];
        this.normMode = 'strongest';
        this.chartMode = 'spectrum';
        this.profScale = 'log';
        this.selected = [];                 // global peak ids, oldest first
        this.colorOf = new Map();
        this.efRef = 'offres';
        this.showEF = true;
        this.showDFPT = false;
        this.showExp = false;
        this.wfTraces = 29;
        this.wfHeight = 3;
        this.sort = { key: 'I', dir: -1 };
        this.playing = false;
        this.pending = null;
        this.dfptCache = new Map();
    }

    /* ----- lifecycle ----- */

    activate() {
        if (this.built) { this.reflow(); return; }
        if (!this.loading) {
            this.view.innerHTML = '<p class="res-status">Loading resonance data…</p>';
            this.loading = fetch(this.dataUrl)
                .then((res) => {
                    if (!res.ok) throw new Error(`HTTP ${res.status}`);
                    return res.json();
                })
                .then((data) => { this.init(data); })
                .catch((err) => {
                    this.loading = null;
                    this.view.innerHTML = `<p class="res-status">Could not load resonance data (${escapeHTML(err.message)}).
                        <button type="button" class="res-retry">Retry</button></p>`;
                    this.view.querySelector('.res-retry').addEventListener('click', () => this.activate());
                });
        }
    }

    deactivate() { this.stopPlay(); }

    reflow() {
        Object.values(this.charts).forEach((chart) => {
            if (chart && chart.renderTo && chart.renderTo.offsetParent !== null) chart.reflow();
        });
    }

    init(data) {
        this.data = data;
        this.energies = data.energies;
        this.nE = data.energies.length;
        this.laserIdx = LASER_LINES.map((nm) => snapIndex(this.energies, nmToEV(nm)));
        const yy = data.geometries[0];
        const yyNorm = computeNormalisers(yy.intensity);
        this.yyMax = yyNorm.M;
        this.ie = yyNorm.mIndex;
        this.buildMenu();
        this.buildView();
        this.built = true;
        this.setGeometry(0);
    }

    /* ----- DOM ----- */

    buildMenu() {
        const options = this.data.geometries.map((g, i) => {
            const title = `light along ${g.prop}, e_i ∥ ${g.ei}, e_s ∥ ${g.es}`;
            return `<option value="${i}" title="${title}">${escapeHTML(g.label)}</option>`;
        }).join('');
        this.menu.innerHTML = `
            <h3>Resonance geometry</h3>
            <select id="res-geom">${options}</select>
            <div class="res-menu-line" id="res-species"></div>
            <div class="res-menu-line res-menu-note" id="res-equiv"></div>
            <div class="res-menu-line" id="res-strongest"></div>
            <div class="res-menu-line res-menu-small">${escapeHTML(this.provenance())}</div>
            <div class="res-menu-line res-menu-small">3D view shows the Raman-tab calculation; it is not linked to these peaks.</div>`;
        this.menu.querySelector('#res-geom').addEventListener('change', (e) => {
            this.stopPlay();
            this.setGeometry(parseInt(e.target.value, 10));
        });
    }

    provenance() {
        const d = this.data;
        const n = d.exciton_count ? `${d.exciton_count.toLocaleString('en').replace(',', ' ')} excitons` : '';
        return ['LDA · GW-BSE · DIP_iR (length gauge)', n, `exciton FWHM ${d.broadening} eV`,
                `Stokes ${d.temperature} K`].filter(Boolean).join(' · ');
    }

    buildView() {
        const lasers = LASER_LINES.map((nm, j) =>
            `<button type="button" class="res-chip res-laser" data-laser="${j}" title="${fmtEV(this.energies[this.laserIdx[j]])}">
                <span class="res-swatch" style="background:${wavelengthColor(nm)}"></span>${nm}</button>`).join('');
        const ticks = this.laserIdx.map((i) => `<option value="${i}"></option>`).join('');
        const efOptions = [`<option value="offres">${this.data.offres_energy.toFixed(1)} eV (off-resonant, exc-ph)</option>`]
            .concat(LASER_LINES.map((nm, j) => `<option value="${j}">${nm} nm (${this.energies[this.laserIdx[j]].toFixed(3)} eV)</option>`))
            .join('');

        this.view.innerHTML = `
            <div class="res-bar">
                <div class="res-row">
                    <button type="button" class="res-play" title="Sweep E_L across the window">▶</button>
                    <input type="range" class="res-slider" min="0" max="${this.nE - 1}" step="1" list="res-ticks" aria-label="Excitation energy">
                    <datalist id="res-ticks">${ticks}</datalist>
                    <input type="text" class="res-energy" size="7" title="eV, or a wavelength such as 488nm">
                    <span class="res-readout"></span>
                    <span class="res-badge"></span>
                </div>
                <div class="res-row">
                    <span class="res-label">Laser</span>${lasers}
                    <button type="button" class="res-chip res-pin" title="Overlay the spectrum at this energy">Pin</button>
                    <select class="res-overlay" aria-label="Overlay presets">
                        <option value="">Overlay ▾</option>
                        <option value="lasers">Pin all laser lines</option>
                        <option value="even">Pin N evenly spaced</option>
                        <option value="clear">Clear pins</option>
                    </select>
                    <input type="number" class="res-even" min="2" max="${MAX_PINS}" value="5" title="N for 'evenly spaced'">
                    <span class="res-pins"></span>
                </div>
                <div class="res-row">
                    <span class="res-label">Normalise</span>
                    <span class="res-seg" data-control="norm">
                        <button type="button" data-value="strongest" title="One scale for the whole window: 1 = the strongest peak at any E_L">Strongest resonance</button>
                        <button type="button" data-value="each" title="Each spectrum (and each profile) scaled to its own maximum">Each curve</button>
                    </span>
                    <span class="res-norm-info"></span>
                </div>
            </div>

            <div class="res-section">
                <div class="res-head">
                    <span class="res-seg" data-control="chartMode">
                        <button type="button" data-value="spectrum">Spectrum</button>
                        <button type="button" data-value="waterfall">Waterfall</button>
                    </span>
                    <span class="res-opts-spectrum">
                        <label><input type="checkbox" class="res-show-ef" checked> EF bars</label>
                        <label title="The DFPT spectrum of the Raman tab, own maximum = 1"><input type="checkbox" class="res-show-dfpt"> Non-resonant DFPT</label>
                    </span>
                    <span class="res-opts-waterfall">
                        <label>Traces <select class="res-wf-traces">
                            <option>8</option><option>15</option><option selected>29</option><option>48</option></select></label>
                        <label>Ridge height <input type="range" class="res-wf-height" min="1" max="10" step="1" value="${this.wfHeight}"></label>
                    </span>
                    <span class="res-hint">drag to zoom · shift-drag to pan</span>
                </div>
                <div class="res-chart res-chart-a" id="res-chart-spectrum"></div>
                <div class="res-chart res-chart-a" id="res-chart-waterfall"></div>
            </div>

            <div class="res-section">
                <div class="res-head">
                    <strong>Resonance profiles</strong>
                    <span class="res-seg" data-control="profScale">
                        <button type="button" data-value="log">Log</button>
                        <button type="button" data-value="lin">Lin</button>
                    </span>
                    <label class="res-exp-toggle" title="Gudima et al. Fig. 4b, hand-digitised; shape only">
                        <input type="checkbox" class="res-show-exp"> Experiment (446 cm⁻¹ ⊥c)</label>
                    <span class="res-hint">click a line to select it · click the plot to set E_L</span>
                </div>
                <div class="res-chart" id="res-chart-profiles"></div>
                <p class="res-foot res-exp-caption"></p>
                <p class="res-foot">Top strip: strongest peak at each E_L. Bottom: GW-BSE absorption for the light polarisations of
                    this geometry and the bright excitons (one common scale). The ~2.7–2.8 eV ⊥c feature is real (Biswas 2017).</p>
            </div>

            <div class="res-section">
                <div class="res-head">
                    <strong class="res-table-title"></strong>
                    <label>EF reference <select class="res-efref">${efOptions}</select></label>
                </div>
                <div class="res-table-wrap"><table class="res-table">
                    <thead><tr>
                        <th data-sort="sel" title="Selected peaks are highlighted in every chart"></th>
                        <th data-sort="freq">ω (cm⁻¹)</th>
                        <th data-sort="sym">Sym</th>
                        <th data-sort="I">I (norm.)</th>
                        <th data-sort="ef" class="res-th-ef"></th>
                        <th data-sort="rel" class="res-th-rel"></th>
                        <th data-sort="peak">Peak at</th>
                    </tr></thead>
                    <tbody></tbody>
                </table></div>
                <p class="res-foot res-caveat"></p>
            </div>

            <details class="res-section res-about">
                <summary>About this calculation — note: the selective enhancement seen in experiment is not reproduced</summary>
                <ul>
                    <li>${escapeHTML(this.data.description)}; phonons at q = 0.</li>
                    <li>I = |ê<sub>s</sub>·α·ê<sub>i</sub>|²(n+1)/ω, summed over degenerate partners, with no ω⁴ prefactor and no
                        absorption, penetration-depth or reflectivity correction: compare with CaF₂-normalised data, not raw counts.</li>
                    <li>A constant exciton linewidth of ${this.data.broadening} eV (FWHM) is used, so features narrower than about 0.05 eV are not resolved.</li>
                    <li>The exciton–phonon tensor only holds the BSE band window, 10–40 % of the full DFPT tensor off resonance:
                        off-resonant ratios and absolute enhancement factors are unreliable, and antiresonance dips may be too deep.</li>
                    <li>TO phonons only (no Fröhlich/LO term). The GW gap (2.79 eV) is below the published ~3.15 eV, so exciton energies
                        are uncertain by about 0.1–0.3 eV.</li>
                    <li>These phonon frequencies come from the exciton–phonon DFPT run and differ by up to ~10 cm⁻¹ from the Raman tab.</li>
                    <li>Notations related by the hexagonal symmetry give the same result up to numerical noise.
                        -X(ZZ)X and -X(YZ)X are still rising at 3.0 eV: their resonance lies above the window.</li>
                </ul>
            </details>`;

        const q = (sel) => this.view.querySelector(sel);
        this.el = {
            play: q('.res-play'), slider: q('.res-slider'), energy: q('.res-energy'), readout: q('.res-readout'),
            badge: q('.res-badge'), pins: q('.res-pins'), pin: q('.res-pin'), overlay: q('.res-overlay'),
            even: q('.res-even'), normInfo: q('.res-norm-info'), tbody: q('.res-table tbody'),
            tableTitle: q('.res-table-title'), caveat: q('.res-caveat'), efref: q('.res-efref'),
            thEF: q('.res-th-ef'), thRel: q('.res-th-rel'), expToggle: q('.res-exp-toggle'),
            expCaption: q('.res-exp-caption'), optsSpectrum: q('.res-opts-spectrum'),
            optsWaterfall: q('.res-opts-waterfall'),
        };
        this.wireControls();
        this.createCharts();
    }

    wireControls() {
        const el = this.el;
        el.slider.addEventListener('input', () => {
            this.stopPlay();
            this.requestIndex(parseInt(el.slider.value, 10));
        });
        el.slider.addEventListener('keydown', (e) => {
            const jump = { PageUp: 8, PageDown: -8 }[e.key];
            if (jump) { e.preventDefault(); this.stopPlay(); this.setIndex(this.ie + jump); }
        });
        const commitEnergy = () => {
            const ev = parseEnergyInput(el.energy.value);
            if (ev === null) { this.flash(el.energy); this.updateReadouts(); return; }
            const lo = this.energies[0], hi = this.energies[this.nE - 1];
            if (ev < lo - 1e-9 || ev > hi + 1e-9) this.flash(el.energy);
            this.stopPlay();
            this.setIndex(snapIndex(this.energies, ev));
        };
        el.energy.addEventListener('keydown', (e) => { if (e.key === 'Enter') commitEnergy(); });
        el.energy.addEventListener('change', commitEnergy);
        el.play.addEventListener('click', () => (this.playing ? this.stopPlay() : this.startPlay()));

        this.view.querySelectorAll('.res-laser').forEach((chip) => {
            chip.addEventListener('click', () => {
                this.stopPlay();
                this.setIndex(this.laserIdx[parseInt(chip.dataset.laser, 10)]);
            });
        });
        el.pin.addEventListener('click', () => this.addPins([this.ie]));
        el.overlay.addEventListener('change', () => {
            const action = el.overlay.value;
            el.overlay.value = '';
            if (action === 'lasers') this.addPins(this.laserIdx);
            if (action === 'even') {
                const n = Math.max(2, Math.min(MAX_PINS, parseInt(el.even.value, 10) || 5));
                el.even.value = n;
                this.pins = [];
                this.addPins(Array.from({ length: n }, (_, j) => Math.round(j * (this.nE - 1) / (n - 1))));
            }
            if (action === 'clear') { this.pins = []; this.renderPins(); }
        });
        el.pins.addEventListener('click', (e) => {
            const chip = e.target.closest('[data-pin]');
            if (!chip) return;
            const i = parseInt(chip.dataset.pin, 10);
            if (e.target.closest('.res-pin-x')) {
                this.pins = this.pins.filter((p) => p !== i);
                this.renderPins();
            } else {
                this.stopPlay();
                this.setIndex(i);
            }
        });

        this.view.querySelectorAll('.res-seg').forEach((seg) => {
            seg.addEventListener('click', (e) => {
                const button = e.target.closest('button[data-value]');
                if (!button) return;
                const control = seg.dataset.control;
                if (this[control === 'norm' ? 'normMode' : control] === button.dataset.value) return;
                if (control === 'norm') { this.normMode = button.dataset.value; this.refreshAll(); }
                if (control === 'chartMode') this.setChartMode(button.dataset.value);
                if (control === 'profScale') { this.profScale = button.dataset.value; this.updateProfiles(); }
                this.syncSegments();
            });
        });
        this.view.querySelector('.res-show-ef').addEventListener('change', (e) => {
            this.showEF = e.target.checked;
            this.layoutSpectrumAxes();
            this.updateSpectrum();
        });
        this.view.querySelector('.res-show-dfpt').addEventListener('change', (e) => {
            this.showDFPT = e.target.checked;
            this.updateDFPT();
        });
        this.view.querySelector('.res-show-exp').addEventListener('change', (e) => {
            this.showExp = e.target.checked;
            const expPeak = this.data.experiment && this.data.experiment.calc_peak;
            if (this.showExp && expPeak !== undefined && !this.selected.includes(expPeak)) this.toggleSelect(expPeak);
            else this.updateProfiles();
        });
        this.view.querySelector('.res-wf-traces').addEventListener('change', (e) => {
            this.wfTraces = parseInt(e.target.value, 10);
            this.updateWaterfall(true);
        });
        this.view.querySelector('.res-wf-height').addEventListener('input', (e) => {
            this.wfHeight = parseInt(e.target.value, 10);
            this.updateWaterfall(true);
        });
        el.efref.addEventListener('change', () => {
            this.efRef = el.efref.value;
            this.computeRefs();
            this.layoutSpectrumAxes();
            this.updateSpectrum();
            this.buildTable();
        });

        const thead = this.view.querySelector('.res-table thead');
        thead.addEventListener('click', (e) => {
            const th = e.target.closest('th[data-sort]');
            if (!th) return;
            const key = th.dataset.sort;
            this.sort = { key, dir: this.sort.key === key ? -this.sort.dir : (key === 'freq' || key === 'sym' ? 1 : -1) };
            this.buildTable();
        });
        el.tbody.addEventListener('click', (e) => {
            const row = e.target.closest('tr[data-peak]');
            if (row) this.toggleSelect(parseInt(row.dataset.peak, 10));
        });
        el.tbody.addEventListener('mouseover', (e) => {
            const row = e.target.closest('tr[data-peak]');
            this.hover(row ? parseInt(row.dataset.peak, 10) : null);
        });
        el.tbody.addEventListener('mouseleave', () => this.hover(null));
        this.view.addEventListener('keydown', (e) => { if (e.key === 'Escape') this.hover(null); });
    }

    flash(input) {
        input.classList.add('res-invalid');
        setTimeout(() => input.classList.remove('res-invalid'), 900);
    }

    syncSegments() {
        const value = { norm: this.normMode, chartMode: this.chartMode, profScale: this.profScale };
        this.view.querySelectorAll('.res-seg').forEach((seg) => {
            seg.querySelectorAll('button').forEach((b) => {
                b.setAttribute('aria-pressed', String(b.dataset.value === value[seg.dataset.control]));
            });
        });
        this.el.optsSpectrum.style.display = this.chartMode === 'spectrum' ? '' : 'none';
        this.el.optsWaterfall.style.display = this.chartMode === 'waterfall' ? '' : 'none';
    }

    /* ----- charts ----- */

    baseOptions(zoomType) {
        return {
            chart: {
                animation: false,
                zooming: {
                    type: zoomType,
                    resetButton: { position: { align: 'right', verticalAlign: 'top', x: -10, y: 30 } },
                },
                panning: { enabled: true, type: zoomType },
                panKey: 'shift',
                pinchType: zoomType,
                spacingTop: 8,
                style: { fontFamily: 'inherit' },
            },
            credits: { enabled: false },
            accessibility: { enabled: false },
            plotOptions: {
                series: {
                    animation: false,
                    states: { inactive: { enabled: false }, hover: { lineWidthPlus: 0 } },
                    marker: { enabled: false },
                },
            },
        };
    }

    createCharts() {
        const self = this;
        const spectrumEl = this.view.querySelector('#res-chart-spectrum');
        const waterfallEl = this.view.querySelector('#res-chart-waterfall');
        const profilesEl = this.view.querySelector('#res-chart-profiles');

        const base = this.baseOptions('xy');
        this.charts.spectrum = Highcharts.chart(spectrumEl, Highcharts.merge(base, {
            // the legend grows down the right side as pins are added
            chart: { zooming: { resetButton: { position: { align: 'left', verticalAlign: 'top', x: 10, y: 10 } } } },
            title: { text: '', align: 'left', style: { fontSize: '14px' } },
            subtitle: { text: '', align: 'left' },
            xAxis: { title: { text: 'Raman shift (cm⁻¹)' }, min: 0, max: X_MAX },
            yAxis: [
                { id: 'y-spec', title: { text: '' }, min: 0, startOnTick: false, endOnTick: false, top: '0%', height: '72%', offset: 0 },
                { id: 'y-ef', title: { text: '' }, type: 'logarithmic', top: '78%', height: '22%', offset: 0,
                  zoomEnabled: false, startOnTick: false, endOnTick: false, gridLineDashStyle: 'Dot' },
            ],
            legend: { align: 'right', verticalAlign: 'top', layout: 'vertical', floating: true, y: 24,
                      itemStyle: { fontSize: '11px', fontWeight: 'normal' }, maxHeight: 140,
                      backgroundColor: 'rgba(255,255,255,0.85)',
                      // visibility is driven by the pins and checkboxes, not the legend
                      events: { itemClick: () => false } },
            tooltip: { useHTML: true, formatter() { return self.pointTooltip(this); } },
            series: [
                { id: 'dfpt', type: 'line', name: 'Non-resonant DFPT (Raman tab), own max = 1; positions differ ≤ 10 cm⁻¹',
                  data: [], color: '#95a5a6', dashStyle: 'Dash', lineWidth: 1.2, zIndex: 1, enableMouseTracking: false,
                  visible: false, showInLegend: false, yAxis: 0 },
                { id: 'live', type: 'line', name: 'current E_L', data: [], color: LIVE_COLOR, lineWidth: 2, zIndex: 3,
                  enableMouseTracking: false, yAxis: 0 },
                { id: 'markers', type: 'scatter', name: 'peaks', data: [], yAxis: 0, zIndex: 4, showInLegend: false,
                  cursor: 'pointer', stickyTracking: false, marker: { enabled: true },
                  dataLabels: { enabled: true, allowOverlap: false, style: { fontSize: '10px', fontWeight: 'normal', textOutline: '2px white' },
                                formatter() { return this.point.options.label || null; } },
                  point: { events: { click() { self.toggleSelect(this.options.peak); } } } },
                { id: 'ef', type: 'column', name: 'EF', data: [], yAxis: 1, zIndex: 2, showInLegend: false,
                  pointWidth: 4, grouping: false, threshold: null, borderWidth: 0, cursor: 'pointer',
                  point: { events: { click() { self.toggleSelect(this.options.peak); } } } },
            ],
        }));

        this.charts.waterfall = Highcharts.chart(waterfallEl, Highcharts.merge(base, {
            chart: { events: { click(e) { self.onWaterfallClick(e); } } },
            title: { text: '', align: 'left', style: { fontSize: '14px' } },
            subtitle: { text: '', align: 'left' },
            legend: { enabled: false },
            xAxis: { title: { text: 'Raman shift (cm⁻¹)' }, min: 0, max: X_MAX },
            yAxis: [
                { title: { text: 'E_L (eV)' }, startOnTick: false, endOnTick: false, gridLineWidth: 0 },
                { title: { text: 'λ (nm)' }, linkedTo: 0, opposite: true, gridLineWidth: 0,
                  startOnTick: false, endOnTick: false,
                  labels: { formatter() { return evToNm(this.value).toFixed(0); } } },
            ],
            tooltip: { useHTML: true, formatter() {
                return `<b>${this.x.toFixed(1)} cm⁻¹</b><br>${fmtEV(self.energies[self.ie])}<br>I = ${fmtI(this.point.options.norm)}`;
            } },
            series: [],
        }));
        waterfallEl.style.display = 'none';

        this.charts.profiles = Highcharts.chart(profilesEl, Highcharts.merge(this.baseOptions('x'), {
            chart: { events: { click(e) { self.setIndexFromClick(e.xAxis[0].value); } } },
            title: { text: '' },
            legend: { enabled: false },
            xAxis: [
                { title: { text: 'Excitation energy E_L (eV)' }, min: this.energies[0], max: this.energies[this.nE - 1],
                  startOnTick: false, endOnTick: false },
                { linkedTo: 0, opposite: true, title: { text: null },
                  labels: { formatter() { return `${evToNm(this.value).toFixed(0)} nm`; }, style: { fontSize: '10px' } } },
            ],
            yAxis: [
                { id: 'y-strip', top: '0%', height: '5%', min: 0, max: 1, title: { text: null }, labels: { enabled: false },
                  gridLineWidth: 0, zoomEnabled: false, offset: 0 },
                { id: 'y-prof', top: '8%', height: '60%', offset: 0, title: { text: '' }, zoomEnabled: false,
                  startOnTick: false, endOnTick: false },
                { id: 'y-abs', top: '72%', height: '28%', offset: 0, min: 0, max: 1.05, title: { text: 'absorption' },
                  zoomEnabled: false, startOnTick: false, endOnTick: false, gridLineWidth: 0 },
                { id: 'y-exc', top: '72%', height: '28%', offset: 0, min: 0, max: 1.1, opposite: true, title: { text: null },
                  labels: { enabled: false }, gridLineWidth: 0, zoomEnabled: false },
            ],
            // narrow line trackers, so clicks between the faint profiles reach
            // the plot (set E_L) instead of selecting the nearest peak
            tooltip: { useHTML: true, snap: 1, formatter() { return self.profileTooltip(this); } },
            plotOptions: {
                series: {
                    events: { click(e) {
                        if (this.options.peak !== undefined) { self.toggleSelect(this.options.peak); return false; }
                        return undefined;
                    } },
                },
            },
            series: [],
        }));
    }

    layoutSpectrumAxes() {
        const chart = this.charts.spectrum;
        const show = this.showEF && this.chartMode === 'spectrum';
        chart.yAxis[0].update({ height: show ? '72%' : '100%' }, false);
        const [lo, hi] = this.efRange;
        const label = this.efRef === 'offres' ? `${this.data.offres_energy.toFixed(1)} eV` : `${LASER_LINES[this.efRef]} nm`;
        chart.yAxis[1].update({
            visible: show, min: lo, max: hi,
            title: { text: `×EF vs ${label}`, style: { fontSize: '10px' } },
            plotLines: this.efRef === 'offres' ? [] : [{ value: 1, color: '#7f8c8d', width: 1, dashStyle: 'Dash' }],
        }, false);
        chart.get('ef').setVisible(show, false);
    }

    /* ----- state changes ----- */

    setGeometry(index) {
        this.geom = index;
        const g = this.data.geometries[index];
        this.g = g;
        this.peakIds = g.peaks;
        this.rowOf = new Map(g.peaks.map((id, k) => [id, k]));
        this.freqs = g.peaks.map((id) => this.data.peaks[id].freq);
        this.norm = computeNormalisers(g.intensity);
        this.segments = dominantSegments(g.intensity);

        this.xs = xGrid(this.freqs);
        this.shapes = buildShapes(this.xs, this.freqs);
        this.centre = this.freqs.map((f) => this.xs.indexOf(f));
        this.wfXs = xGrid(this.freqs, 1.0);
        this.wfShapes = buildShapes(this.wfXs, this.freqs);
        this.spectrumTop = { strongest: this.curveMax('strongest'), each: this.curveMax('each') };

        const kept = this.selected.filter((id) => this.rowOf.has(id));
        if (!kept.length) {
            this.selected = [];
            this.colorOf.clear();
            Array.from(this.norm.P.keys()).sort((a, b) => this.norm.P[b] - this.norm.P[a]).slice(0, 4)
                .forEach((k) => this.addSelected(g.peaks[k]));
        } else {
            this.selected = kept;
            for (const id of Array.from(this.colorOf.keys())) if (!kept.includes(id)) this.colorOf.delete(id);
        }

        const sel = this.menu.querySelector('#res-geom');
        if (sel) sel.value = String(index);
        const hasExp = !!(this.data.experiment && this.data.experiment.geometries.includes(g.id));
        this.el.expToggle.style.display = hasExp ? '' : 'none';
        if (!hasExp) {
            this.showExp = false;
            this.el.expToggle.querySelector('input').checked = false;
        }

        this.computeRefs();
        this.updateMenu();
        this.syncSegments();
        this.layoutSpectrumAxes();
        this.refreshAll();
        if (this.showDFPT) this.updateDFPT();
    }

    // max of the broadened curve over the window, so the y axis never rescales while sliding
    curveMax(mode) {
        let top = 0;
        for (let i = 0; i < this.nE; i++) {
            for (let a = 0; a < this.freqs.length; a++) {
                let v = 0;
                const j = this.centre[a];
                for (let k = 0; k < this.freqs.length; k++) v += normalise(this.g.intensity, this.norm, k, i, mode) * this.shapes[k][j];
                if (v > top) top = v;
            }
        }
        return Math.max(top, 1);
    }

    computeRefs() {
        const g = this.g;
        const ref = g.peaks.map((_, k) => (this.efRef === 'offres' ? g.offres[k] : g.intensity[k][this.laserIdx[this.efRef]]));
        this.ref = ref;
        this.refMax = Math.max(...ref);
        this.refSum = ref.reduce((a, b) => a + b, 0);
        this.refRow = this.rowOf.get(g.ref_peak);
        // fixed log range of the EF bars over the whole window
        let lo = Infinity, hi = 0;
        g.intensity.forEach((row, k) => {
            if (!(ref[k] >= RESONLY * this.refMax)) return;
            for (const v of row) {
                const ef = v / ref[k];
                if (ef > 0) { lo = Math.min(lo, ef); hi = Math.max(hi, ef); }
            }
        });
        if (!(hi > 0)) { lo = 0.1; hi = 10; }
        this.efRange = [Math.max(1e-2, 10 ** Math.floor(Math.log10(lo))), 10 ** Math.ceil(Math.log10(hi))];
    }

    updateMenu() {
        const g = this.g, P = this.data.peaks;
        const counts = g.allowed.map((s) => symLabel(s)).join(' + ');
        this.menu.querySelector('#res-species').textContent = `${counts} · ${g.peaks.length} peaks`;

        const partnerId = EQUIVALENT[g.id];
        const partner = this.data.geometries.find((x) => x.id === partnerId);
        let equiv = '';
        if (partner) {
            let diff = 0;
            g.peaks.forEach((id, k) => {
                const kk = partner.peaks.indexOf(id);
                const a = g.intensity[k], b = kk >= 0 ? partner.intensity[kk] : null;
                for (let i = 0; i < this.nE; i++) diff = Math.max(diff, Math.abs(a[i] - (b ? b[i] : 0)));
            });
            const rel = diff / this.norm.M;
            const how = rel < 1e-6 ? 'identical' : `computed curves differ < ${Math.max(0.01, Math.ceil(rel * 1000) / 10)} %`;
            equiv = `≡ ${partner.label} by hexagonal symmetry (${how})`;
        }
        this.menu.querySelector('#res-equiv').textContent = equiv;

        const m = this.norm, top = P[g.peaks[m.mPeak]];
        const pct = 100 * m.M / this.yyMax;
        const edge = m.mIndex === this.nE - 1 ? ' — at the window edge; the resonance lies above 3.0 eV' : '';
        this.menu.querySelector('#res-strongest').textContent =
            `Strongest: ${top.freq.toFixed(1)} ${symLabel(top.sym)} at ${fmtEV(this.energies[m.mIndex])} · ` +
            `${pct >= 10 ? pct.toFixed(0) : pct.toFixed(1)} % of ${this.data.geometries[0].label}${edge}`;
    }

    requestIndex(i) {
        this.pending = i;
        if (this.frame) return;
        this.frame = requestAnimationFrame(() => {
            this.frame = null;
            if (this.pending !== null) { const p = this.pending; this.pending = null; this.setIndex(p); }
        });
    }

    setIndex(i) {
        this.ie = Math.max(0, Math.min(this.nE - 1, i));
        this.updateReadouts();
        this.updateSpectrum();
        this.updateWaterfallCursor();
        this.moveCursor();
        this.updateTableValues();
    }

    setIndexFromClick(ev) {
        if (!Number.isFinite(ev)) return;
        this.stopPlay();
        this.setIndex(snapIndex(this.energies, ev));
    }

    onWaterfallClick(e) {
        const ev = e.yAxis[0].value;
        const x = e.xAxis[0].value;
        let near = -1, best = 6;
        this.freqs.forEach((f, k) => { if (Math.abs(f - x) < best) { best = Math.abs(f - x); near = k; } });
        if (near >= 0) this.toggleSelect(this.peakIds[near]);
        this.setIndexFromClick(ev);
    }

    refreshAll() {
        this.updateReadouts();
        this.renderPins(false);
        this.updateSpectrum(true);
        this.updateWaterfall(true);
        this.updateProfiles();
        this.buildTable();
    }

    setChartMode(mode) {
        const spectrum = this.charts.spectrum, waterfall = this.charts.waterfall;
        const from = mode === 'waterfall' ? spectrum : waterfall;
        const to = mode === 'waterfall' ? waterfall : spectrum;
        const ext = from.xAxis[0].getExtremes();
        this.chartMode = mode;
        spectrum.renderTo.style.display = mode === 'spectrum' ? '' : 'none';
        waterfall.renderTo.style.display = mode === 'waterfall' ? '' : 'none';
        this.syncSegments();
        to.reflow();
        if (mode === 'waterfall') this.updateWaterfall(true);
        else { this.layoutSpectrumAxes(); this.updateSpectrum(); }
        if (ext.userMin !== undefined || ext.userMax !== undefined) {
            to.xAxis[0].setExtremes(ext.min, ext.max);
            if (!to.resetZoomButton) to.showResetZoom();
        } else if (to.resetZoomButton || to.xAxis[0].getExtremes().userMin !== undefined) {
            to.zoomOut();
        }
    }

    /* ----- selection ----- */

    addSelected(id) {
        if (this.selected.length >= MAX_SELECTED) {
            const old = this.selected.shift();
            this.colorOf.delete(old);
        }
        const used = new Set(this.colorOf.values());
        this.colorOf.set(id, PALETTE.find((c) => !used.has(c)) || PALETTE[0]);
        this.selected.push(id);
    }

    toggleSelect(id) {
        if (this.selected.includes(id)) {
            this.selected = this.selected.filter((x) => x !== id);
            this.colorOf.delete(id);
        } else {
            this.addSelected(id);
        }
        this.updateSpectrum();
        this.updateWaterfallGuides();
        this.updateProfiles();
        this.buildTable();
        const prof = this.view.querySelector('#res-chart-profiles');
        const box = prof.getBoundingClientRect(), vbox = this.view.getBoundingClientRect();
        if (box.bottom < vbox.top + 60 || box.top > vbox.bottom - 60) prof.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
    }

    // emphasise one peak's profile while its table row is hovered
    hover(id) {
        if (this.hovered === id) return;
        const prev = this.hovered;
        this.hovered = id;
        this.el.tbody.querySelectorAll('tr').forEach((tr) => {
            tr.classList.toggle('res-hover', id !== null && parseInt(tr.dataset.peak, 10) === id);
        });
        const chart = this.charts.profiles;
        chart.series.forEach((s) => {
            const o = s.options;
            if (o.peak === undefined) return;
            if (o.peak === id) s.update({ lineWidth: 3.5, zIndex: 10 }, false);
            else if (o.peak === prev) s.update({ lineWidth: o.baseWidth, zIndex: o.baseZ }, false);
        });
        chart.redraw(false);
    }

    /* ----- readouts and pins ----- */

    updateReadouts() {
        const ev = this.energies[this.ie];
        const nm = evToNm(ev);
        this.el.slider.value = String(this.ie);
        this.el.energy.value = ev.toFixed(3);
        this.el.readout.innerHTML = `= ${nm.toFixed(1)} nm <span class="res-swatch" style="background:${wavelengthColor(nm)}"></span>`;
        const near = LASER_LINES.findIndex((l) => Math.abs(nmToEV(l) - ev) <= 0.004);
        this.el.badge.textContent = near >= 0 ? `≈ ${LASER_LINES[near]} nm line` : '';
        this.view.querySelectorAll('.res-laser').forEach((chip) => {
            chip.classList.toggle('res-active', this.laserIdx[parseInt(chip.dataset.laser, 10)] === this.ie);
        });
        this.el.pin.disabled = this.pins.length >= MAX_PINS || this.pins.includes(this.ie);
        this.el.pin.title = this.pins.length >= MAX_PINS ? `max ${MAX_PINS} — remove one` : 'Overlay the spectrum at this energy';
        const m = this.norm, top = this.data.peaks[this.g.peaks[m.mPeak]];
        this.el.normInfo.textContent = this.normMode === 'strongest'
            ? `1 = ${top.freq.toFixed(1)} ${symLabel(top.sym)} at ${this.energies[m.mIndex].toFixed(3)} eV`
            : '1 = strongest peak of each spectrum';
    }

    addPins(indices) {
        for (const i of indices) {
            if (this.pins.length >= MAX_PINS) break;
            if (!this.pins.includes(i)) this.pins.push(i);
        }
        this.pins.sort((a, b) => a - b);
        this.renderPins();
    }

    // f(E) = s(E)/M: what an 'each curve' spectrum was scaled by
    hiddenFactor(i) { return this.norm.s[i] / this.norm.M; }

    renderPins(redraw = true) {
        this.el.pins.innerHTML = this.pins.map((i) => {
            const ev = this.energies[i], nm = evToNm(ev);
            const f = this.normMode === 'each' ? ` <span class="res-factor">×${Number(this.hiddenFactor(i).toPrecision(2))}</span>` : '';
            return `<span class="res-chip res-pinchip" data-pin="${i}" title="Move E_L here">
                <span class="res-swatch" style="background:${wavelengthColor(nm)}"></span>${ev.toFixed(3)} eV · ${nm.toFixed(1)} nm${f}
                <button type="button" class="res-pin-x" title="Remove">×</button></span>`;
        }).join('');

        const chart = this.charts.spectrum;
        chart.series.filter((s) => s.options.isPin).forEach((s) => s.remove(false));
        this.pins.forEach((i, n) => {
            const nm = evToNm(this.energies[i]);
            const f = this.normMode === 'each' ? `, ×${Number(this.hiddenFactor(i).toPrecision(2))}` : '';
            chart.addSeries({
                isPin: true, pinIndex: i, type: 'line', yAxis: 0, zIndex: 2,
                name: `${this.energies[i].toFixed(3)} eV · ${nm.toFixed(1)} nm${f}`,
                color: wavelengthColor(nm), lineWidth: 1.2, dashStyle: n % 2 ? 'ShortDash' : 'Solid',
                enableMouseTracking: false, data: this.curveData(i),
            }, false);
        });
        this.updateReadouts();
        this.drawEnergyLines();
        if (redraw) chart.redraw(false);
    }

    /* ----- spectrum ----- */

    weights(i) {
        return this.g.peaks.map((_, k) => normalise(this.g.intensity, this.norm, k, i, this.normMode));
    }

    curveData(i) {
        const y = spectrumAt(this.shapes, this.weights(i));
        const out = new Array(this.xs.length);
        for (let j = 0; j < this.xs.length; j++) out[j] = [this.xs[j], y[j]];
        return out;
    }

    updateSpectrum(axes = false) {
        const chart = this.charts.spectrum;
        if (!chart) return;
        const i = this.ie, g = this.g, P = this.data.peaks;
        const w = this.weights(i);
        const y = spectrumAt(this.shapes, w);
        const live = new Array(this.xs.length);
        for (let j = 0; j < this.xs.length; j++) live[j] = [this.xs[j], y[j]];
        const ev = this.energies[i];
        chart.get('live').setData(live, false, false, false);

        const order = Array.from(w.keys()).sort((a, b) => y[this.centre[b]] - y[this.centre[a]]);
        const labelled = new Set(order.slice(0, 5));
        const strongest = this.norm.s[i];
        const markers = g.peaks.map((id, k) => {
            const selected = this.colorOf.get(id);
            const visible = g.intensity[k][i] >= WEAK * strongest;
            const fill = selected || '#95a5a6';
            return {
                x: this.freqs[k], y: y[this.centre[k]], peak: id,
                label: labelled.has(k) ? `${Math.round(this.freqs[k])} ${symLabel(P[id].sym)}` : null,
                marker: { symbol: SYMBOLS[P[id].sym] || 'circle', radius: selected ? 5.5 : 4,
                          fillColor: visible || selected ? fill : Highcharts.color(fill).setOpacity(0.5).get(),
                          lineColor: selected ? '#333' : '#7f8c8d', lineWidth: 1 },
            };
        });
        chart.get('markers').setData(markers, false, false, false);

        const efPoints = [];
        if (this.showEF) {
            g.peaks.forEach((id, k) => {
                const cell = efCell(g.intensity[k][i], this.ref[k], this.refMax, g.intensity[k][i] >= WEAK * strongest);
                if (cell.state === 'resonly') return;
                const color = this.colorOf.get(id) || '#95a5a6';
                efPoints.push({ x: this.freqs[k], y: cell.value, peak: id,
                                color: cell.state === 'weak' ? Highcharts.color(color).setOpacity(0.3).get() : color });
            });
        }
        chart.get('ef').setData(efPoints, false, false, false);

        const dominant = P[g.peaks[order[0]]];
        const sMax = strongest / this.norm.M;
        const tail = this.normMode === 'each'
            ? `hidden factor ×${Number(this.hiddenFactor(i).toPrecision(2))} (spectrum max = ${(100 * sMax).toFixed(sMax < 0.1 ? 1 : 0)} % of the strongest resonance)`
            : `spectrum max = ${(100 * sMax).toFixed(sMax < 0.1 ? 1 : 0)} % of the strongest resonance`;
        chart.setTitle(
            { text: `Raman spectrum — ${escapeHTML(g.label)}` },
            { text: `E_L = ${fmtEV(ev)} · dominant ${dominant.freq.toFixed(1)} ${symLabel(dominant.sym)} · ${tail}` }, false);

        if (axes) {
            // fixed for the whole window, so resonance shows up as height while sliding
            const top = this.spectrumTop[this.normMode];
            chart.yAxis[0].update({
                max: 1.05 * top,
                title: { text: this.normMode === 'each' ? 'I / I_max (this spectrum)' : 'I / I_max (window)' },
            }, false);
            chart.yAxis[0].setExtremes(0, 1.05 * top, false);
        }
        chart.redraw(false);
    }

    updateDFPT() {
        const chart = this.charts.spectrum;
        const series = chart.get('dfpt');
        const box = this.view.querySelector('.res-show-dfpt');
        if (!this.showDFPT) { series.update({ visible: false, showInLegend: false }); return; }
        const g = this.g;
        const url = this.ramanFileFor(g.prop, g.ei, g.es);
        let entry = this.dfptCache.get(url);
        if (!entry) {
            entry = fetch(url).then((res) => {
                if (!res.ok) throw new Error(`HTTP ${res.status}`);
                return res.json();
            }).then((ph) => {
                const gi = ph.gamma_index || 0;
                const f = ph.eigenvalues[gi], inten = ph.raman_intensities || [];
                const xs = [];
                for (let x = 0; x <= X_MAX; x += 0.5) xs.push(x);
                const y = xs.map((x) => {
                    let v = 0;
                    for (let n = 0; n < f.length; n++) {
                        if (inten[n] > 0) { const d = x - f[n]; v += inten[n] * GAMMA * GAMMA / (d * d + GAMMA * GAMMA); }
                    }
                    return v;
                });
                const top = Math.max(...y) || 1;
                return xs.map((x, j) => [x, y[j] / top]);
            });
            this.dfptCache.set(url, entry);
        }
        entry.then((data) => {
            if (this.g !== g || !this.showDFPT) return;
            series.setData(data, false);
            series.update({ visible: true, showInLegend: true });
        }).catch(() => {
            // forget the failure so ticking the box again retries
            this.dfptCache.delete(url);
            if (this.g !== g) return;
            box.checked = false;
            box.parentElement.title = 'Could not load the Raman-tab spectrum for this geometry; tick again to retry';
            this.showDFPT = false;
            series.update({ visible: false, showInLegend: false });
        });
    }

    pointTooltip(ctx) {
        const id = ctx.point.options.peak;
        if (id === undefined) return false;
        return this.peakSummary(id, this.ie);
    }

    // profile = true describes the value plotted on a resonance profile, which
    // in 'each' mode is scaled to the peak's own window maximum
    peakSummary(id, i, profile = false) {
        const g = this.g, k = this.rowOf.get(id), p = this.data.peaks[id];
        const I = g.intensity[k][i];
        const s = this.norm.s[i];
        const eachProfile = profile && this.normMode === 'each';
        const n = eachProfile ? I / this.norm.P[k] : normalise(g.intensity, this.norm, k, i, this.normMode);
        const how = eachProfile ? '÷ this peak\'s window maximum'
            : (this.normMode === 'each' ? '÷ strongest of this spectrum' : '÷ strongest resonance');
        const maxLine = eachProfile
            ? `<br>its maximum = ${this.strengthPct(k)} % of the strongest resonance, at ${fmtEV(this.energies[this.norm.Pidx[k]])}` : '';
        let total = 0;
        for (const row of g.intensity) total += row[i];
        const ef = efCell(I, this.ref[k], this.refMax, I >= WEAK * s);
        const rel = this.relCell(k, i);
        const refLabel = this.efRef === 'offres' ? `${this.data.offres_energy.toFixed(1)} eV` : `${LASER_LINES[this.efRef]} nm`;
        const efText = ef.state === 'resonly' ? 'res.-only' : formatEF(ef.value);
        return `<b>${p.freq.toFixed(1)} cm⁻¹ ${symLabel(p.sym)}</b> · ${fmtEV(this.energies[i])}<br>
            I = ${fmtI(n)} (${how})${maxLine}<br>
            ${(100 * I / s).toFixed(0)} % of the strongest peak here · share ${(100 * I / total).toFixed(0)} %<br>
            ×EF ${efText} vs ${refLabel} · ×rel ${rel.state === 'resonly' ? 'res.-only' : formatEF(rel.value)}`;
    }

    // a peak's window maximum as a percentage of the strongest resonance
    strengthPct(k) {
        const pct = 100 * this.norm.P[k] / this.norm.M;
        return pct >= 10 ? pct.toFixed(0) : pct.toPrecision(2);
    }

    // selectivity: EF of this peak over EF of the reference peak
    relCell(k, i) {
        const g = this.g, r = this.refRow;
        const mine = efCell(g.intensity[k][i], this.ref[k], this.refMax, true);
        const theirs = efCell(g.intensity[r][i], this.ref[r], this.refMax, true);
        if (mine.state === 'resonly' || theirs.state === 'resonly') return { value: NaN, state: 'resonly' };
        return { value: mine.value / theirs.value, state: 'ok' };
    }

    /* ----- waterfall ----- */

    waterfallIndices() {
        const n = this.wfTraces;
        return Array.from({ length: n }, (_, j) => Math.round(j * (this.nE - 1) / (n - 1)));
    }

    updateWaterfall(rebuild = false) {
        const chart = this.charts.waterfall;
        if (!chart || this.chartMode !== 'waterfall') return;
        if (rebuild) {
            while (chart.series.length) chart.series[0].remove(false);
            const idx = this.waterfallIndices();
            const dE = (this.energies[this.nE - 1] - this.energies[0]) / (idx.length - 1);
            this.wfAmp = this.wfHeight * dE;
            idx.forEach((i, n) => {
                const E = this.energies[i];
                const y = spectrumAt(this.wfShapes, this.weights(i));
                const data = new Array(this.wfXs.length);
                for (let j = 0; j < data.length; j++) data[j] = [this.wfXs[j], E + this.wfAmp * y[j]];
                const label = this.normMode === 'each' && n % Math.ceil(idx.length / 10) === 0
                    ? `×${Number(this.hiddenFactor(i).toPrecision(2))}` : null;
                chart.addSeries({
                    type: 'area', data, threshold: E, color: wavelengthColor(evToNm(E)), fillColor: '#ffffff',
                    lineWidth: 1.2, zIndex: idx.length - n, enableMouseTracking: false,
                    dataLabels: label ? { enabled: true, formatter() { return this.point.index === this.series.data.length - 1 ? label : null; },
                                          align: 'left', x: 2, style: { fontSize: '9px', fontWeight: 'normal' } } : { enabled: false },
                }, false);
            });
            chart.addSeries({ id: 'wf-current', type: 'line', color: '#000', lineWidth: 2, zIndex: idx.length + 1, data: [] }, false);
            const top = this.spectrumTop[this.normMode];
            chart.yAxis[0].update({ min: this.energies[0] - 0.01, max: this.energies[this.nE - 1] + this.wfAmp * top * 1.05 }, false);
            chart.setTitle({ text: `Waterfall — ${escapeHTML(this.g.label)}` },
                { text: `${idx.length} spectra, colour = laser wavelength · ${this.normMode === 'each'
                    ? 'each spectrum scaled to its own maximum (×factor at the right)' : 'one scale for the whole window'}` }, false);
            this.updateWaterfallGuides(false);
        }
        this.updateWaterfallCursor(false);
        chart.redraw(false);
    }

    updateWaterfallCursor(redraw = true) {
        const chart = this.charts.waterfall;
        if (!chart || this.chartMode !== 'waterfall') return;
        const current = chart.get('wf-current');
        if (!current) return;
        const E = this.energies[this.ie];
        const w = this.weights(this.ie);
        const y = spectrumAt(this.wfShapes, w);
        const data = new Array(this.wfXs.length);
        for (let j = 0; j < data.length; j++) data[j] = { x: this.wfXs[j], y: E + this.wfAmp * y[j], norm: y[j] };
        current.setData(data, false, false, false);
        if (redraw) chart.redraw(false);
    }

    updateWaterfallGuides(redraw = true) {
        const chart = this.charts.waterfall;
        if (!chart) return;
        const axis = chart.xAxis[0];
        (this.wfGuides || []).forEach((id) => axis.removePlotLine(id));
        this.wfGuides = this.selected.map((id) => {
            const pid = `wf-guide-${id}`;
            axis.addPlotLine({ id: pid, value: this.data.peaks[id].freq, color: this.colorOf.get(id), width: 1, dashStyle: 'Dash', zIndex: 50 });
            return pid;
        });
        if (redraw && this.chartMode === 'waterfall') chart.redraw(false);
    }

    /* ----- profiles ----- */

    updateProfiles() {
        const chart = this.charts.profiles;
        if (!chart) return;
        const g = this.g, P = this.data.peaks, E = this.energies;
        const log = this.profScale === 'log';
        this.hovered = null;
        while (chart.series.length) chart.series[0].remove(false);

        // dominant-peak strip
        const shade = ['#d5d8dc', '#aeb6bf'];
        const strip = [];
        const stripLabels = [];
        this.segments.forEach((seg, n) => {
            const id = g.peaks[seg.k];
            const color = this.colorOf.get(id) || shade[n % 2];
            for (let i = seg.i0; i <= seg.i1; i++) strip.push({ x: E[i], y: 1, color, seg: n });
            if (E[seg.i1] - E[seg.i0] >= 0.015) {
                stripLabels.push({ x: (E[seg.i0] + E[seg.i1]) / 2, y: 0.5, seg: n, label: `${Math.round(P[id].freq)}` });
            }
        });
        chart.addSeries({ type: 'column', yAxis: 0, data: strip, pointPadding: 0, groupPadding: 0, borderWidth: 0,
                          pointRange: E[1] - E[0], turboThreshold: 0, isStrip: true }, false);
        chart.addSeries({ type: 'scatter', yAxis: 0, data: stripLabels, marker: { enabled: false }, isStrip: true,
                          dataLabels: { enabled: true, formatter() { return this.point.options.label; }, y: 8,
                                        style: { fontSize: '9px', fontWeight: 'normal', textOutline: 'none' } } }, false);

        // profiles: everything faint, the selection on top in its colours
        const floor = (v) => (log ? Math.max(v, LOG_FLOOR) : v);
        const value = (k, i) => {
            const v = g.intensity[k][i];
            return this.normMode === 'each' ? v / this.norm.P[k] : v / this.norm.M;
        };
        g.peaks.forEach((id, k) => {
            const color = this.colorOf.get(id);
            const data = new Array(this.nE);
            for (let i = 0; i < this.nE; i++) data[i] = { x: E[i], y: floor(value(k, i)), raw: value(k, i) };
            const pk = this.norm.Pidx[k];
            data[pk] = { ...data[pk], marker: color ? { enabled: true, symbol: 'circle', radius: 4, fillColor: '#fff', lineColor: color, lineWidth: 1.5 } : undefined };
            if (color) {
                data[this.nE - 1] = { ...data[this.nE - 1],
                    dataLabels: { enabled: true, format: `${Math.round(P[id].freq)} ${symLabel(P[id].sym)}` +
                                      `${this.normMode === 'each' ? ` · max ${this.strengthPct(k)} %` : ''}${pk === this.nE - 1 ? ' (edge)' : ''}`,
                                  align: 'left', x: 4, crop: false, overflow: 'allow', style: { color, fontSize: '10px', textOutline: '2px white' } } };
            }
            chart.addSeries({
                type: 'line', yAxis: 1, data, peak: id, turboThreshold: 0, cursor: 'pointer',
                color: color || FAINT_COLOR, lineWidth: color ? 2.5 : 1, opacity: color ? 1 : 0.7,
                zIndex: color ? 5 : 1, baseWidth: color ? 2.5 : 1, baseZ: color ? 5 : 1,
            }, false);
        });
        if (this.normMode === 'strongest') {
            const data = E.map((e, i) => {
                let t = 0;
                for (const row of g.intensity) t += row[i];
                return { x: e, y: floor(t / this.norm.M), raw: t / this.norm.M };
            });
            chart.addSeries({ type: 'line', yAxis: 1, data, name: 'Σ all peaks', color: LIVE_COLOR, dashStyle: 'Dash',
                              lineWidth: 1.5, zIndex: 4, isTotal: true, turboThreshold: 0 }, false);
        }
        this.addExperiment(value, floor);
        this.addAbsorption();

        // room for the Σ curve, which can exceed the strongest single peak
        const top = this.normMode === 'strongest' ? Math.max(1, this.totalMax()) : 1;
        const title = { text: this.normMode === 'each' ? 'I / I_max (this peak)' : 'I / I_max (window)' };
        chart.yAxis[1].update(log
            ? { type: 'logarithmic', min: LOG_FLOOR, max: 1.5 * top, title }
            : { type: 'linear', min: 0, max: 1.05 * top, title }, false);
        chart.redraw(false);
    }

    totalMax() {
        let top = 0;
        for (let i = 0; i < this.nE; i++) {
            let t = 0;
            for (const row of this.g.intensity) t += row[i];
            top = Math.max(top, t / this.norm.M);
        }
        return top;
    }

    addExperiment(value, floor) {
        const chart = this.charts.profiles;
        const exp = this.data.experiment;
        const axis = chart.xAxis[0];
        axis.removePlotBand('res-exp-band');
        this.el.expCaption.textContent = '';
        if (!this.showExp || !exp || !exp.geometries.includes(this.g.id)) return;
        const k = this.rowOf.get(exp.calc_peak);
        const E = this.energies;
        const lo = E[0], hi = 2.68;
        const pts = exp.points.filter(([e]) => e >= lo - 1e-9 && e <= E[this.nE - 1] + 1e-9);
        if (k === undefined || !pts.length) return;
        let calcMax = 0;
        for (let i = 0; i < this.nE; i++) if (E[i] <= hi) calcMax = Math.max(calcMax, value(k, i));
        const expMax = Math.max(...pts.filter(([e]) => e <= hi).map(([, v]) => v));
        const scale = calcMax / expMax;
        chart.addSeries({
            type: 'scatter', yAxis: 1, name: exp.label, isExp: true, zIndex: 8,
            marker: { enabled: true, symbol: 'diamond', radius: 4, fillColor: '#fff', lineColor: '#000', lineWidth: 1.2 },
            data: pts.map(([e, v]) => ({ x: e, y: floor(v * scale), raw: v * scale, measured: v })),
        }, false);
        axis.addPlotBand({ id: 'res-exp-band', from: lo, to: 2.755, color: '#fdf2e9', zIndex: 0,
                           label: { text: 'measured range (450–650 nm)', style: { fontSize: '9px', color: '#b9770e' }, y: 30 } });
        this.el.expCaption.textContent = `Diamonds: ${exp.label}, ${exp.note}, scaled once to the calculated ` +
            `${this.data.peaks[exp.calc_peak].freq.toFixed(1)} A₁ profile (shape only). The peak position agrees; the calculated ` +
            'profile falls faster above 2.6 eV, and the selective enhancement seen in experiment is not reproduced.';
    }

    addAbsorption() {
        const chart = this.charts.profiles;
        const g = this.g, abs = this.data.absorption, exc = this.data.excitons, E = this.energies;
        const lo = E[0], hi = E[this.nE - 1];
        const perp = abs.x.map((v, j) => (v + abs.y[j]) / 2);
        let scale = 0;
        abs.energies.forEach((e, j) => { if (e >= lo && e <= hi) scale = Math.max(scale, perp[j]); });
        const wantPerp = g.ei !== 'z' || g.es !== 'z';
        const wantPar = g.ei === 'z' || g.es === 'z';
        if (wantPerp) {
            chart.addSeries({ type: 'area', yAxis: 2, name: 'GW-BSE absorption ⊥c', isAbs: true, color: '#7f8c8d',
                              fillColor: 'rgba(189,195,199,0.5)', lineWidth: 1, enableMouseTracking: false,
                              data: abs.energies.map((e, j) => [e, perp[j] / scale]) }, false);
        }
        if (wantPar) {
            const both = wantPerp;
            chart.addSeries({ type: both ? 'line' : 'area', yAxis: 2, name: 'GW-BSE absorption ∥c', isAbs: true,
                              color: '#34495e', dashStyle: both ? 'Dash' : 'Solid', fillColor: 'rgba(52,73,94,0.25)',
                              lineWidth: 1, enableMouseTracking: false,
                              data: abs.energies.map((e, j) => [e, abs.z[j] / scale]) }, false);
        }
        // bright-exciton sticks on the same common scale
        const sticks = [];
        let top = 0;
        exc.energies.forEach((e, j) => { if (e >= lo && e <= hi) top = Math.max(top, (exc.x[j] + exc.y[j]) / 2); });
        const add = (proj, color, pol) => {
            let polTop = 0;
            exc.energies.forEach((e, j) => { if (e >= lo && e <= hi) polTop = Math.max(polTop, proj(j)); });
            exc.energies.forEach((e, j) => {
                if (e < lo || e > hi || proj(j) < 0.01 * polTop) return;
                sticks.push({ x: e, y: proj(j) / top, color, pol, j });
            });
        };
        if (wantPerp) add((j) => (exc.x[j] + exc.y[j]) / 2, '#7f8c8d', '⊥c');
        if (wantPar) add((j) => exc.z[j], '#34495e', '∥c');
        sticks.sort((a, b) => b.y - a.y);
        const shown = sticks.slice(0, 60);
        shown.slice(0, 5).forEach((s) => { s.label = s.x.toFixed(3); });
        chart.addSeries({ type: 'column', yAxis: 3, name: 'bright excitons', isExc: true, pointWidth: 1.5, minPointLength: 2,
                          borderWidth: 0, grouping: false, data: shown.sort((a, b) => a.x - b.x), turboThreshold: 0,
                          dataLabels: { enabled: true, formatter() { return this.point.options.label || null; }, crop: false,
                                        style: { fontSize: '9px', fontWeight: 'normal', textOutline: '2px white' } } }, false);
    }

    // laser lines and pins on the profile chart (they span all three panes)
    drawEnergyLines() {
        const chart = this.charts.profiles;
        if (!chart) return;
        const axis = chart.xAxis[0];
        (this.profLines || []).forEach((id) => axis.removePlotLine(id));
        const lines = [];
        LASER_LINES.forEach((nm, j) => {
            lines.push({ id: `res-laser-${j}`, value: this.energies[this.laserIdx[j]], color: '#95a5a6', width: 1, dashStyle: 'Dot', zIndex: 3,
                         label: { text: `${nm}`, rotation: 0, y: 10, x: 2, style: { fontSize: '9px', color: '#7f8c8d' } } });
        });
        this.pins.forEach((i) => {
            lines.push({ id: `res-pin-${i}`, value: this.energies[i], color: wavelengthColor(evToNm(this.energies[i])), width: 1, zIndex: 4 });
        });
        lines.forEach((l) => axis.addPlotLine(l));
        this.profLines = lines.map((l) => l.id);
        this.moveCursor();
    }

    moveCursor() {
        const chart = this.charts.profiles;
        if (!chart) return;
        const axis = chart.xAxis[0];
        axis.removePlotLine('res-cursor');
        axis.addPlotLine({ id: 'res-cursor', value: this.energies[this.ie], color: LIVE_COLOR, width: 1.5, zIndex: 6,
                           label: { text: 'E_L', rotation: 0, verticalAlign: 'bottom', y: -4, x: 3, style: { fontSize: '10px', color: LIVE_COLOR } } });
    }

    profileTooltip(ctx) {
        const o = ctx.point.options, s = ctx.series.options;
        if (s.isStrip) {
            const seg = this.segments[o.seg];
            if (!seg) return false;
            const p = this.data.peaks[this.g.peaks[seg.k]];
            return `${this.energies[seg.i0].toFixed(3)}–${this.energies[seg.i1].toFixed(3)} eV: <b>${p.freq.toFixed(1)} ${symLabel(p.sym)}</b> strongest`;
        }
        if (s.isExc) {
            const exc = this.data.excitons, j = o.j;
            return `exciton at <b>${exc.energies[j].toFixed(3)} eV</b> (${o.pol})<br>|d<sub>x</sub>|² ${fmtI(exc.x[j])} · |d<sub>y</sub>|² ${fmtI(exc.y[j])} · |d<sub>z</sub>|² ${fmtI(exc.z[j])}`;
        }
        if (s.isExp) {
            const exp = this.data.experiment;
            return `<b>Experiment</b> ${o.x.toFixed(3)} ± ${exp.err_E} eV · I ${o.measured} ± ${exp.err_I} (CaF₂-norm.)<br>` +
                `${escapeHTML(exp.note)} · exp. 446 ↔ calc. ${this.data.peaks[exp.calc_peak].freq.toFixed(1)} A₁ ` +
                '(446.3 E₂ lies 4 cm⁻¹ away)<br>points below 2.3 eV are not shown (≤ 2.1 eV features are defect-related)';
        }
        if (s.isTotal) return `<b>Σ all peaks</b> · ${fmtEV(o.x)}<br>${fmtI(o.raw)} of the strongest single resonance`;
        if (s.peak === undefined) return false;
        const i = snapIndex(this.energies, o.x);
        return this.peakSummary(s.peak, i, true);
    }

    /* ----- table ----- */

    buildTable() {
        const g = this.g, P = this.data.peaks;
        const refLabel = this.efRef === 'offres' ? `${this.data.offres_energy.toFixed(1)} eV` : `${LASER_LINES[this.efRef]} nm`;
        const refPeak = P[g.ref_peak];
        this.el.thEF.textContent = `×EF vs ${refLabel}`;
        this.el.thRel.textContent = `×rel ÷ ${refPeak.freq.toFixed(1)} ${symLabel(refPeak.sym)}`;
        this.el.thRel.title = 'Enhancement relative to the reference peak: 1 = enhanced alike, > 1 = gains on it';
        this.view.querySelectorAll('.res-table th[data-sort]').forEach((th) => {
            th.classList.toggle('res-sorted', th.dataset.sort === this.sort.key);
            th.dataset.dir = this.sort.dir > 0 ? '▲' : '▼';
        });

        this.el.tbody.innerHTML = g.peaks.map((id) => {
            const k = this.rowOf.get(id), p = P[id];
            const color = this.colorOf.get(id);
            const pk = this.norm.Pidx[k];
            const edge = pk === 0 || pk === this.nE - 1 ? ' (edge)' : '';
            return `<tr data-peak="${id}" class="${color ? 'res-selected' : ''}">
                <td><span class="res-sel" style="${color ? `background:${color};border-color:${color}` : ''}"></span></td>
                <td>${p.freq.toFixed(1)}</td>
                <td>${symLabel(p.sym)}</td>
                <td class="res-cell-I"><span class="res-bar-fill"></span><span class="res-bar-text"></span></td>
                <td class="res-cell-ef"></td>
                <td class="res-cell-rel"></td>
                <td>${this.energies[pk].toFixed(3)} eV / ${evToNm(this.energies[pk]).toFixed(0)} nm${edge}</td>
            </tr>`;
        }).join('');
        this.rows = new Map();
        this.el.tbody.querySelectorAll('tr').forEach((tr) => {
            this.rows.set(parseInt(tr.dataset.peak, 10), {
                tr, bar: tr.querySelector('.res-bar-fill'), text: tr.querySelector('.res-bar-text'),
                ef: tr.querySelector('.res-cell-ef'), rel: tr.querySelector('.res-cell-rel'),
            });
        });
        this.el.caveat.textContent = `EF vs ${refLabel} comes from the same exciton–phonon calculation, whose tensor holds only ` +
            'the BSE band window (10–40 % of the full DFPT tensor off resonance): absolute EFs are overestimated, unevenly between ' +
            `modes, so compare modes with ×rel. "res.-only" = below ${RESONLY * 100} % of the strongest peak at the reference energy; ` +
            `grey = below ${WEAK * 100} % of the strongest peak at E_L. Laser-line references are more robust, but 532 nm is already pre-resonant.`;
        this.updateTableValues();
    }

    updateTableValues() {
        if (!this.rows) return;
        const g = this.g, i = this.ie, s = this.norm.s[i];
        let total = 0;
        for (const row of g.intensity) total += row[i];
        const label = g.label;
        const gtot = total / this.refSum;
        this.el.tableTitle.textContent = `Peaks at E_L = ${fmtEV(this.energies[i])} · ${label} · whole spectrum ${formatEF(gtot)}`;

        const values = new Map();
        g.peaks.forEach((id, k) => {
            const I = g.intensity[k][i];
            const n = this.normMode === 'each' ? I / s : I / this.norm.M;
            const visible = I >= WEAK * s;
            const ef = efCell(I, this.ref[k], this.refMax, visible);
            const rel = this.relCell(k, i);
            values.set(id, { n, ef, rel, visible, k });
            const row = this.rows.get(id);
            row.bar.style.width = `${Math.min(100, 100 * n)}%`;
            row.text.textContent = fmtI(n);
            const efHTML = ef.state === 'resonly'
                ? '<span class="res-pill" title="Below 1 % of the strongest peak at the reference energy in this calculation; the ratio is not meaningful">res.-only</span>'
                : formatEF(ef.value);
            row.ef.innerHTML = efHTML;
            row.ef.classList.toggle('res-weak', ef.state === 'weak');
            row.rel.innerHTML = rel.state === 'resonly' ? '<span class="res-pill">res.-only</span>' : formatEF(rel.value);
            row.rel.classList.toggle('res-weak', !visible);
        });

        // keep the rows in the chosen order
        const key = this.sort.key, dir = this.sort.dir, P = this.data.peaks;
        const val = (id) => {
            const v = values.get(id);
            switch (key) {
                case 'freq': return P[id].freq;
                case 'sym': return P[id].sym + String(P[id].freq).padStart(8, '0');
                case 'sel': return this.selected.includes(id) ? 1 : 0;
                case 'ef': return v.ef.state === 'resonly' ? -Infinity : v.ef.value;
                case 'rel': return v.rel.state === 'resonly' ? -Infinity : v.rel.value;
                case 'peak': return this.norm.Pidx[v.k];
                default: return v.n;
            }
        };
        const order = g.peaks.slice().sort((a, b) => {
            const va = val(a), vb = val(b);
            return va < vb ? -dir : va > vb ? dir : 0;
        });
        const current = Array.from(this.el.tbody.children).map((tr) => parseInt(tr.dataset.peak, 10));
        if (order.some((id, n) => id !== current[n])) order.forEach((id) => this.el.tbody.appendChild(this.rows.get(id).tr));
    }

    /* ----- play ----- */

    // sweeps the window in PLAY_SECONDS whatever the frame rate
    startPlay() {
        if (this.ie >= this.nE - 1) this.setIndex(0);
        this.playing = true;
        this.el.play.textContent = '❚❚';
        this.el.play.title = 'Pause';
        const msPerStep = 1000 * PLAY_SECONDS / (this.nE - 1);
        let last = null, carry = 0;
        const step = (t) => {
            if (!this.playing) return;
            if (last !== null) carry += t - last;
            last = t;
            const n = Math.floor(carry / msPerStep);
            if (n > 0) {
                carry -= n * msPerStep;
                this.setIndex(this.ie + n);
            }
            if (this.ie >= this.nE - 1) { this.stopPlay(); return; }
            this.playFrame = requestAnimationFrame(step);
        };
        this.playFrame = requestAnimationFrame(step);
    }

    stopPlay() {
        if (!this.playing) return;
        this.playing = false;
        cancelAnimationFrame(this.playFrame);
        if (this.el) { this.el.play.textContent = '▶'; this.el.play.title = 'Sweep E_L across the window'; }
    }
}

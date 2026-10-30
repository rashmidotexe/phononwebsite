/*
 * zumbaresonancetab.js  —  adds the "Resonance" plot view to a PhononWebpage
 *
 * The shared setPlotView() only knows 'dispersion' and 'raman'. Both the tab
 * buttons and plotRaman() call it through the instance, so overriding it on
 * the instance also keeps the Resonance view up when Geometry 1 reloads
 * behind it, without touching phononwebpage.js.
 */

export function installResonanceTab(page, resonanceView, doc = document) {
    const originalSetPlotView = page.setPlotView.bind(page);
    const byId = (id) => doc.getElementById(id);
    const show = (id, visible, display = '') => {
        const el = byId(id);
        if (el) el.style.display = visible ? display : 'none';
    };

    page.setPlotView = function (view) {
        const resonance = view === 'resonance';
        const container = doc.querySelector('.flex-container');
        const widthChanged = !!container && container.classList.contains('resonance-active') !== resonance;

        if (container) container.classList.toggle('resonance-active', resonance);
        show('resonance-view', resonance, 'block');
        show('resonance-menu', resonance);
        show('panel1', !resonance);
        show('panel2', !resonance);

        if (resonance) {
            page.activePlotView = 'resonance';
            show('highcharts', false);
            show('raman-spectrum', false);
            doc.querySelectorAll('#plot-view-controls [data-plot-view]').forEach((button) => {
                button.setAttribute('aria-selected', String(button.dataset.plotView === 'resonance'));
            });
            resonanceView.activate();
        } else {
            resonanceView.deactivate();
            originalSetPlotView(view);
        }

        if (widthChanged) {
            // the 3D renderer and every chart size themselves on window resize
            const win = doc.defaultView;
            win.requestAnimationFrame(() => {
                win.dispatchEvent(new win.Event('resize'));
                if (resonance) resonanceView.reflow();
            });
        }
    };
}

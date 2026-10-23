/*
 * rruff.js — experimental Raman spectra of the RRUFF project
 * (B. Lafuente, R. T. Downs, H. Yang and N. Stone, Highlights in Mineralogical
 * Crystallography, De Gruyter (2016) 1-30) shown over the computed spectrum.
 *
 * python/phononweb/scripts/prepare_rruff.py writes one small json per matching
 * material to data/rruff/<mp-id>.json and lists them in data/rruff/index.json,
 * so only materials that have spectra are requested.
 */

export const RRUFF_DATA_URL = 'data/rruff/';

export const RRUFF_REFERENCE = {
    text: 'RRUFF (Lafuente et al. 2016)',
    url: 'https://doi.org/10.1515/9783110417104-003',
};

const indexCache = new Map();
const spectraCache = new Map();

function loadIndex(baseUrl) {
    if (!indexCache.has(baseUrl)) {
        let request = fetch(baseUrl + 'index.json')
            .then((response) => response.ok ? response.json() : null)
            .then((index) => index ? new Set(index.materials) : null)
            .catch(() => null);
        indexCache.set(baseUrl, request);
    }
    return indexCache.get(baseUrl);
}

export function fetchRruffSpectra(materialId, baseUrl = RRUFF_DATA_URL) {
    /*
    experimental spectra of a material, null if there are none
    */
    if (!materialId) {
        return Promise.resolve(null);
    }
    let key = baseUrl + materialId;
    if (!spectraCache.has(key)) {
        let request = loadIndex(baseUrl).then((materials) => {
            if (!materials || !materials.has(materialId)) { return null; }
            return fetch(baseUrl + materialId + '.json')
                .then((response) => response.ok ? response.json() : null);
        }).catch((error) => {
            spectraCache.delete(key);
            console.warn('Unable to load RRUFF spectra for ' + materialId + ':', error);
            return null;
        });
        spectraCache.set(key, request);
    }
    return spectraCache.get(key);
}

export function getRruffSeries(payload, height, maxFrequency) {
    /*
    Highcharts line series of the spectra scaled to the given height,
    only the first one visible
    */
    if (!payload || !Array.isArray(payload.spectra)) {
        return [];
    }
    return payload.spectra.map((spectrum, index) => {
        let data = [];
        for (let i = 0; i < spectrum.intensity.length; i++) {
            let x = spectrum.start + i * spectrum.step;
            if (maxFrequency != null && x > maxFrequency) { break; }
            data.push([x, spectrum.intensity[i] / 1000 * height]);
        }
        let laser = spectrum.laser_nm ? ' · ' + spectrum.laser_nm + ' nm' : '';
        return {
            name: spectrum.name + ' ' + spectrum.rruff_id + laser,
            type: 'line',
            data: data,
            color: ['#27ae60', '#8e44ad', '#d35400'][index % 3],
            lineWidth: 1,
            marker: { enabled: false },
            enableMouseTracking: false,
            visible: index === 0,
            rruffUrl: spectrum.url,
        };
    });
}

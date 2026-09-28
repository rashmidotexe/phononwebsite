/*
 * mpdielectric.js — Born effective charges and electronic dielectric tensors of
 * the Materials Project DFPT phonon database (G. Petretto et al.).
 *
 * The phonon band structures the site loads from the OpenData bucket only carry
 * frequencies and eigendisplacements. The Born charges and epsilon_infinity of
 * the same calculations are in the static "phonon" collection of the
 * materialsproject-build bucket, stored as a Delta Lake table of parquet files
 * partitioned by phonon_method. They are fetched at runtime (CORS is open on
 * the bucket) and indexed by the numeric part of the mp-id.
 *
 * The collection identifies materials with Materials Project AlphaIDs, the
 * base-26 lower-case spelling of the old integer id ("aaaaaaaa" = 0).
 */

import { parquetReadObjects } from 'hyparquet';
import { compressors } from 'hyparquet-compressors';

const BUCKET = 'https://materialsproject-build.s3.amazonaws.com';
const DFPT_PREFIX = 'static-collections/phonon/phonon_method=dfpt/';

let dielectricPromise = null;

export function alphaIdToNumber(alphaId) {
    let number = 0;
    for (let i = 0; i < alphaId.length; i++) {
        number = number * 26 + (alphaId.charCodeAt(i) - 97);
    }
    return number;
}

async function findDfptParquetUrl() {
    /*
    the parquet file name carries a random id that changes when the collection
    is rewritten, so list the partition and take the newest file
    */
    let response = await fetch(BUCKET + '/?list-type=2&prefix=' + encodeURIComponent(DFPT_PREFIX));
    if (!response.ok) {
        throw new Error('Unable to list the Materials Project phonon collection: HTTP ' + response.status);
    }
    let xml = new DOMParser().parseFromString(await response.text(), 'application/xml');
    let newest = null;
    for (let node of Array.from(xml.getElementsByTagName('Contents'))) {
        let key = node.getElementsByTagName('Key')[0].textContent;
        let modified = node.getElementsByTagName('LastModified')[0].textContent;
        if (key.endsWith('.parquet') && (!newest || modified > newest.modified)) {
            newest = { key: key, modified: modified };
        }
    }
    if (!newest) {
        throw new Error('No DFPT parquet file found in the Materials Project phonon collection');
    }
    return BUCKET + '/' + newest.key;
}

export function indexDielectricRows(rows) {
    /*
    rows with identifier, born, epsilon_electronic and structure -> Map keyed by
    the numeric mp-id with the Born charges in the order of the site labels
    */
    let byId = new Map();
    for (let i = 0; i < rows.length; i++) {
        let row = rows[i];
        if (!row.identifier || !row.born || !row.epsilon_electronic || !row.structure) {
            continue;
        }
        byId.set(alphaIdToNumber(row.identifier), {
            born: row.born,
            dielectric: row.epsilon_electronic,
            labels: row.structure.sites.map((site) => site.label),
        });
    }
    return byId;
}

async function fetchDielectricTable() {
    let url = await findDfptParquetUrl();
    let response = await fetch(url);
    if (!response.ok) {
        throw new Error('Unable to download the Materials Project DFPT dielectric data: HTTP ' + response.status);
    }
    let rows = await parquetReadObjects({
        file: await response.arrayBuffer(),
        compressors: compressors,
        columns: ['identifier', 'born', 'epsilon_electronic', 'structure'],
    });
    return indexDielectricRows(rows);
}

export function loadMaterialsProjectDielectric() {
    /*
    Map of mp-id number -> { born, dielectric, labels }, fetched once per page
    */
    if (!dielectricPromise) {
        dielectricPromise = fetchDielectricTable();
        dielectricPromise.catch(() => { dielectricPromise = null; });
    }
    return dielectricPromise;
}

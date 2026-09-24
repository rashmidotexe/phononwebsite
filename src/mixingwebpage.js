import { MaterialsProjectDB } from './mpdb.js';
import { LocalPhononDB } from './localphonondb.js';
import { PhononJson } from './phononjson.js';
import { PhononHighcharts } from './phononhighcharts.js';
import { PhononWebpage } from './phononwebpage.js';
import {
    buildMixingEndpoint,
    sampleMixingPair,
    computeMixedPhonon,
    computeReferenceEigenvalues,
} from './mixingphonons.js';
import * as utils from './utils.js';

class RawPhononJson extends PhononJson {
    /*
    reuse the PhononJson download/decompression and keep the parsed json as is
    */
    getFromJson(json, callback) {
        this.raw = json;
        callback();
    }
}

function loadRawJson(url) {
    return new Promise(function(resolve, reject) {
        let loader = new RawPhononJson();
        loader.getFromURL(url, function() { resolve(loader.raw); }, {
            onError: function(error) {
                reject(new Error(error && error.message ? error.message : 'Unable to load ' + url));
            }
        });
    });
}

export class MixingHighcharts extends PhononHighcharts {
    /*
    dispersion plot with the bands of material 1 drawn as a dotted reference
    underneath the mixed bands
    */

    getGraph(phonon) {
        super.getGraph(phonon);

        let reference = phonon.reference_eigenvalues;
        if (!reference) {
            return;
        }

        let dists = phonon.distances;
        let nbands = reference[0].length;
        for (let n = 0; n < nbands; n++) {
            for (let i = 0; i < phonon.line_breaks.length; i++) {
                let data = [];
                for (let k = phonon.line_breaks[i][0]; k < phonon.line_breaks[i][1]; k++) {
                    data.push([dists[k], reference[k][n]]);
                }
                this.highcharts.push({
                    name: 'reference',
                    isReferenceSeries: true,
                    color: '#555555',
                    dashStyle: 'Dot',
                    lineWidth: 1.5,
                    zIndex: 4,
                    enableMouseTracking: false,
                    showInLegend: false,
                    states: { inactive: { opacity: 1 } },
                    marker: { enabled: false },
                    data: data
                });
            }
        }
    }
}

export class MixingWebpage extends PhononWebpage {
    /*
    Phonon page restricted to the curated list of 2-atom zinc blende / diamond
    materials in data/mixing/models.json, drawn from the PhononDB (A. Togo)
    and Materials Project (G. Petretto et al.) databases.

    Material 1 is picked in the regular list, material 2 in the second list and
    the character slider sets x in D_alloy = (1-x) D_1 + x D_2.
    */

    constructor(visualizer, dispersion) {
        super(visualizer, dispersion);
        this.mixingList = 'data/mixing/models.json';
        // copies of the PhononDB files outside the LFS-tracked data/phonondb2017
        this.phonondbRoot = 'data/mixing/phonondb';
        this.material1 = null;
        this.material2 = null;
        this.character = 0;
        this.material2FilterQuery = '';
        this.rawCache = new Map();
        this.pairCache = null;
        this.mixToken = 0;
    }

    setMaterial2List(dom_mat2) { this.dom_mat2 = dom_mat2; }

    setMaterial2FilterInput(dom_input) {
        dom_input.on('input', () => {
            this.material2FilterQuery = dom_input.val() || '';
            this.renderMaterial2Menu();
        });
    }

    setCharacterInput(dom_range, dom_value) {
        this.dom_character_value = dom_value;
        dom_range.on('input', () => {
            this.character = Number(dom_range.val()) / 100;
            this.updateCharacterLabel();
            if (this.material2) {
                this.refreshMixing(false);
            }
        });
        this.character = Number(dom_range.val()) / 100;
        this.updateCharacterLabel();
    }

    updateCharacterLabel() {
        if (!this.dom_character_value) {
            return;
        }
        let percent = Math.round(this.character * 100);
        let name = this.material2 ? utils.format_formula_html(this.material2.name) : 'material 2';
        this.dom_character_value.html(percent + " % " + name);
    }

    updateMenu() {
        let self = this;

        this.materialsIndex = [];
        if (this.dom_mat) { this.dom_mat.empty(); }
        if (this.dom_ref) { this.dom_ref.empty(); }

        let phonondb = new LocalPhononDB();
        let mpdb = new MaterialsProjectDB();

        function getReference(db) {
            return db.author+", "+"<a href="+db.url+">"+db.name+"</a> ("+db.year+")";
        }

        function makeMaterial(entry) {
            let m = {
                id: entry.id,
                name: entry.name,
                type: "json",
                link: "https://materialsproject.org/materials/mp-"+entry.id,
            };
            if (entry.source === phonondb.name) {
                m.source = phonondb.name;
                m.reference = getReference(phonondb);
                m.url = self.phonondbRoot+"/mp-"+entry.id+".json.gz";
            } else {
                m.source = mpdb.name;
                m.reference = getReference(mpdb);
                m.url = "https://materialsproject-parsed.s3.amazonaws.com/ph-bandstructures/dfpt/mp-"+entry.id+".json.gz";
            }
            return m;
        }

        function addMaterials(entries) {
            for (let i=0; i<entries.length; i++) {
                self.materialsIndex.push(makeMaterial(entries[i]));
            }
            self.renderMaterialsMenu();
        }

        $.get(this.mixingList, function(entries) {
            addMaterials(entries.filter((entry) => entry.source === phonondb.name));

            //materials project entries are only listed if the OpenData bucket is reachable
            let mpEntries = entries.filter((entry) => entry.source === mpdb.name);
            mpdb.checkAvailability(function(isAvailable) {
                if (isAvailable) {
                    addMaterials(mpEntries);
                } else {
                    console.log("Skipping Materials Project phonons because the OpenData bucket is unreachable from this browser.");
                }
            });
        });
    }

    renderMaterialsMenu() {
        super.renderMaterialsMenu();

        let materialsHeading = document.querySelector("#material-list h3");
        if (materialsHeading && this.dom_mat) {
            materialsHeading.textContent = "Choose material 1 (" + this.dom_mat.children().length + "):";
        }
        this.markSelectedMaterials();
        this.renderMaterial2Menu();
    }

    getReferenceIndices() {
        let indices = new Map();
        let sorted = this.materialsIndex.slice().sort(this.compareMaterialsForMenu.bind(this));
        for (let i = 0; i < sorted.length; i++) {
            let key = this.getMaterialReferenceKey(sorted[i]);
            if (!indices.has(key)) {
                indices.set(key, indices.size + 1);
            }
        }
        return indices;
    }

    renderMaterial2Menu() {
        let dom_mat2 = this.dom_mat2;
        if (!dom_mat2) {
            return;
        }
        dom_mat2.empty();

        let tokens = this.material2FilterQuery
            .toLowerCase()
            .split(/[\s,]+/)
            .filter((token) => token.length > 0);
        let referenceIndices = this.getReferenceIndices();
        let materials = this.materialsIndex
            .filter((material) => this.materialMatchesFilter(material, tokens))
            .filter((material) => this.isReferenceEnabled(this.getMaterialReferenceKey(material)))
            .sort(this.compareMaterialsForMenu.bind(this));

        let addItem = (label, material) => {
            let li = document.createElement("LI");
            let a = document.createElement("A");
            a.innerHTML = label;
            a.onclick = () => { this.selectMaterial2(material); };
            if ((this.material2 ? this.material2.url : null) === (material ? material.url : null)) {
                li.className = "selected";
            }
            li.appendChild(a);
            dom_mat2.append(li);
        };

        addItem("None", null);
        for (let i = 0; i < materials.length; i++) {
            let m = materials[i];
            let index = referenceIndices.get(this.getMaterialReferenceKey(m));
            addItem(utils.format_formula_html(m.name) + " [" + index + "]", m);
        }
    }

    markSelectedMaterials() {
        if (!this.dom_mat || !this.material1) {
            return;
        }
        let url = this.material1.url;
        let sorted = this.materialsIndex
            .filter((material) => this.materialMatchesFilter(material, this.getMaterialFilterTokens()))
            .filter((material) => this.isReferenceEnabled(this.getMaterialReferenceKey(material)))
            .sort(this.compareMaterialsForMenu.bind(this));
        let items = this.dom_mat.children();
        for (let i = 0; i < items.length && i < sorted.length; i++) {
            items[i].className = sorted[i].url === url ? "selected" : "";
        }
    }

    findMaterial(url) {
        for (let i = 0; i < this.materialsIndex.length; i++) {
            if (this.materialsIndex[i].url === url) {
                return this.materialsIndex[i];
            }
        }
        return null;
    }

    loadURL(url_vars, callback) {
        /*
        clicking material 1 (or the default material) selects material 1
        */
        let url = url_vars.json;
        let entry = this.findMaterial(url);
        let name = entry ? entry.name : String(url_vars.name || '').replace(/\s*\[\d+\]\s*$/, '').replace(/<[^>]+>/g, '');
        this.material1 = {
            url: url,
            name: name,
            link: url_vars.link,
            url_vars: url_vars,
        };
        this.markSelectedMaterials();
        this.refreshMixing(true, callback);
    }

    selectMaterial2(material) {
        this.material2 = material ? { url: material.url, name: material.name } : null;
        this.renderMaterial2Menu();
        this.updateCharacterLabel();
        this.refreshMixing(true);
    }

    getRawJson(url) {
        if (!this.rawCache.has(url)) {
            let promise = loadRawJson(url);
            promise.catch(() => { this.rawCache.delete(url); });
            this.rawCache.set(url, promise);
        }
        return this.rawCache.get(url);
    }

    async getSampledPair() {
        let key = this.material1.url + "::" + this.material2.url;
        if (this.pairCache && this.pairCache.key === key) {
            return this.pairCache;
        }
        let raws = await Promise.all([this.getRawJson(this.material1.url), this.getRawJson(this.material2.url)]);
        let endpoint1 = buildMixingEndpoint(raws[0], this.material1);
        let endpoint2 = buildMixingEndpoint(raws[1], this.material2);
        let sampled = sampleMixingPair(endpoint1, endpoint2);
        if (!sampled.qpoints.length) {
            throw new Error("The two materials have no high-symmetry segment in common.");
        }
        let reference = await computeReferenceEigenvalues(sampled);
        this.pairCache = { key: key, sampled: sampled, reference: reference };
        return this.pairCache;
    }

    getMixingTitle() {
        let percent = Math.round(this.character * 100);
        return utils.format_formula_html(this.material1.name) +
            " + " + utils.format_formula_html(this.material2.name) +
            " (" + percent + " % " + utils.format_formula_html(this.material2.name) + ")";
    }

    async refreshMixing(pairChanged, callback) {
        if (!this.material1) {
            return;
        }

        if (!this.material2) {
            // plain material 1, as on the phonon page
            this.mixToken += 1;
            super.loadURL(this.material1.url_vars, callback);
            return;
        }

        let token = ++this.mixToken;
        let x = this.character;
        if (pairChanged) {
            this.startLoadingFeedback(this.material1.name + " + " + this.material2.name);
        }

        try {
            let pair = await this.getSampledPair();
            let data = await computeMixedPhonon(pair.sampled, x);
            if (token !== this.mixToken) {
                return;
            }

            let phonon = new PhononJson();
            phonon.getFromInternalJson(data, () => {});
            // alloy eigenvectors are in the atom gauge: animate with the atomic phase
            phonon.addatomphase = true;
            phonon.reference_eigenvalues = pair.reference;

            this.phonon = phonon;
            this.k = Math.min(this.k, phonon.kpoints.length - 1);
            delete this.link;
            if (pairChanged) {
                this.finishLoadingFeedback();
                this.loadCallback();
            } else {
                this.update();
            }
            this.name = this.getMixingTitle();
            this.updatePage();
            if (callback) {
                callback();
            }
        } catch (error) {
            console.error(error);
            if (token === this.mixToken) {
                this.failLoadingFeedback({ message: error.message });
            }
        }
    }
}

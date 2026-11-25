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
    getParentFamily,
    getWorldMixing,
} from './mixingphonons.js';
import { loadMaterialsProjectDielectric } from './mpdielectric.js';
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
    Phonon page restricted to the curated list of 2-atom zinc blende / diamond,
    2-atom rocksalt, 3-atom fluorite and antifluorite, 4-atom wurtzite /
    lonsdaleite, 8-atom chalcopyrite, ternary wurtzite (beta-NaFeO2, enargite and
    Cu2SiS3 type) and cation-ordered ternary rocksalt (alpha-NaFeO2,
    gamma-LiFeO2, Li2ZrO3 type, ...), fluorite (scheelite, CaUO4 type, ...) and
    antifluorite (MgAgAs type half-Heuslers LiMgN, LiZnAs; Li5SiP3) materials, and ABX3 perovskites and
    antiperovskites (the 5-atom cell, cubic or polar, and the supercells of the
    octahedral tilts: Pnma, R-3c, ...), the 2H, 1T and 3R polytypes of the
    transition-metal dichalcogenides, and rutiles (with the CaCl2-type
    distortion) with the trirutiles and other cation orderings on them,
    corundums with the ilmenites and A3BX6 ordering, CsCl-type binaries, and the
    zinc-blende-derived Li3Bi-type (D0_3) binaries with the full and inverse
    Heuslers on them (all four fcc sites filled; the half-Heuslers fill three and
    are listed with the antifluorites) in
    data/mixing/models.json, drawn from the PhononDB (A. Togo) and Materials
    Project (G. Petretto et al.) databases. Materials mix within their family; a
    binary with force constants (PhononDB) can also be mixed with an ordered
    ternary on its lattice (zinc blende with chalcopyrite, the others with their
    ternaries), folded into the cell of the ternary. Materials of different
    structures on the same lattice (the fcc-site families, or perovskite,
    antiperovskite and CsCl) mix site by site on a common cell, with a vacancy
    on every site only one of them fills (getWorldMixing).

    Material 1 is picked in the regular list, material 2 in the second list and
    the character slider sets the fraction x of material 2 in the virtual
    crystal (see mixingphonons.js).
    */

    constructor(visualizer, dispersion) {
        super(visualizer, dispersion);
        this.mixingList = 'data/mixing/models.json';
        // copies of the PhononDB files outside the LFS-tracked data/phonondb2017
        this.phonondbRoot = 'data/mixing/phonondb';
        this.material1 = null;
        this.material2 = null;
        this.character = 0;
        // site masses follow the character unless tuned separately
        this.massLinked = true;
        this.massCharacter = 0;
        this.material2FilterQuery = '';
        this.rawCache = new Map();
        this.pairCache = null;
        this.mixToken = 0;
        // a slider-triggered mixture is waiting to run
        this.mixScheduled = false;
    }

    setMaterial2List(dom_mat2, dom_note) { this.dom_mat2 = dom_mat2; this.dom_mat2_note = dom_note; }

    setMaterial2FilterInput(dom_input) {
        dom_input.on('input', () => {
            this.material2FilterQuery = dom_input.val() || '';
            this.renderMaterial2Menu();
        });
    }

    setMixingNote(dom_note) { this.dom_mixing_note = dom_note; }

    updateMixingNote(nacMode, sampled = null) {
        if (!this.dom_mixing_note) {
            return;
        }
        // the dipole-dipole term is normally recomputed from the mixed Born charges and
        // dielectric tensor; only the fallback, where it is interpolated, is flagged,
        // with its cause: a PhononDB file without Born charges, or Materials Project
        // charges that could not be loaded
        let note = "";
        let union = sampled && sampled.union;
        if (union && this.material2) {
            // two structures on a common lattice: say which sites only one of them fills
            let only1 = union.map2.filter((i) => i < 0).length;
            let only2 = union.map1.filter((i) => i < 0).length;
            let parts = [];
            if (only1) {
                parts.push(only1 + " filled only in " + sampled.endpoint1.name);
            }
            if (only2) {
                parts.push(only2 + " only in " + sampled.endpoint2.name);
            }
            note = "Mixed site by site on a common cell of " + union.map1.length + " sites" +
                (parts.length ? " (" + parts.join(", ") + "; vacancies in the other, occupied in proportion to the mixture)" : "") + ". ";
        }
        if (nacMode === 'interpolated' && this.material2) {
            let missing = sampled ? [sampled.endpoint1, sampled.endpoint2].filter((endpoint) => !endpoint.nac) : [];
            let phonondb = missing.filter((endpoint) => endpoint.evaluate || endpoint.anyQ);
            note += (phonondb.length && phonondb.length === missing.length
                ? "LO–TO splitting interpolated: no Born charges in the PhononDB file of " +
                    phonondb.map((endpoint) => endpoint.name).join(" and ") + "."
                : "LO–TO splitting interpolated: the Materials Project Born charges could not be loaded.");
        }
        this.dom_mixing_note.text(note);
    }

    setCharacterInput(dom_range, dom_value) {
        this.dom_character_value = dom_value;
        dom_range.on('input', () => {
            this.character = Number(dom_range.val()) / 100;
            this.updateCharacterLabel();
            if (this.massLinked) {
                this.setMassCharacter(this.character);
            }
            if (this.material2) {
                this.refreshMixing(false);
            }
        });
        this.character = Number(dom_range.val()) / 100;
        this.updateCharacterLabel();
    }

    setMassInput(dom_toggle, dom_range, dom_value) {
        /*
        mass slider: disabled and following the character slider by default,
        independent once "tune mass separately" is checked
        */
        this.dom_mass_range = dom_range;
        this.dom_mass_value = dom_value;
        dom_toggle.on('change', () => {
            this.massLinked = !dom_toggle.prop('checked');
            dom_range.prop('disabled', this.massLinked);
            if (this.massLinked && this.massCharacter !== this.character) {
                this.setMassCharacter(this.character);
                if (this.material2) {
                    this.refreshMixing(false);
                }
            }
            this.updateMassLabel();
        });
        dom_range.on('input', () => {
            if (this.massLinked) {
                return;
            }
            this.setMassCharacter(Number(dom_range.val()) / 100);
            if (this.material2) {
                this.refreshMixing(false);
            }
        });
        this.massLinked = !dom_toggle.prop('checked');
        dom_range.prop('disabled', this.massLinked);
        this.setMassCharacter(this.massLinked ? this.character : Number(dom_range.val()) / 100);
    }

    setMassCharacter(fraction) {
        this.massCharacter = fraction;
        if (this.dom_mass_range) {
            this.dom_mass_range.val(Math.round(fraction * 100));
        }
        this.updateMassLabel();
    }

    getMassFraction() {
        return this.massLinked ? this.character : this.massCharacter;
    }

    updateMassLabel() {
        if (!this.dom_mass_value) {
            return;
        }
        let percent = Math.round(this.getMassFraction() * 100);
        let pair = this.pairCache && this.material2 && this.pairCache.key === this.material1.url + "::" + this.material2.url
            ? this.pairCache.sampled
            : null;
        if (!pair) {
            this.dom_mass_value.html(this.massLinked ? "follows character %" : percent + " %");
            return;
        }
        // site masses of the virtual crystal, e.g. "Ga 69.72 · P/As 51.35 amu"
        let y = this.getMassFraction();
        let e1 = pair.endpoint1;
        let e2 = pair.endpoint2;
        // one entry per distinct site (2 in zinc blende, A / B / anion in chalcopyrite)
        let sites = [];
        let add = (entry) => {
            if (sites.indexOf(entry) < 0) {
                sites.push(entry);
            }
        };
        if (pair.union) {
            // two structures on a common cell: the sites of the union, a vacancy ("–")
            // where only one of them has an atom, occupied with the character x
            let x = this.character;
            pair.union.map1.forEach((i1, s) => {
                let i2 = pair.union.map2[s];
                let mass = i1 >= 0 && i2 >= 0 ? (1 - y) * e1.siteMasses[i1] + y * e2.siteMasses[i2]
                    : i1 >= 0 ? (1 - x) * e1.siteMasses[i1] : x * e2.siteMasses[i2];
                let t1 = i1 >= 0 ? e1.siteTypes[i1] : "–";
                let t2 = i2 >= 0 ? e2.siteTypes[i2] : "–";
                add((t1 === t2 ? t1 : t1 + "/" + t2) + " " + mass.toFixed(2));
            });
        } else {
            e1.siteTypes.forEach((type, s) => {
                let mass = (1 - y) * e1.siteMasses[s] + y * e2.siteMasses[s];
                let types = type === e2.siteTypes[s] ? type : type + "/" + e2.siteTypes[s];
                add(types + " " + mass.toFixed(2));
            });
        }
        let prefix = this.massLinked ? "" : percent + " % " + utils.format_formula_html(this.material2.name) + ": ";
        this.dom_mass_value.html(prefix + sites.join(" · ") + " amu");
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
                // only materials of the same structure family can be mixed
                family: entry.family || "zincblende",
                spacegroup: entry.spacegroup,
            };
            if (entry.source === phonondb.name) {
                m.source = phonondb.name;
                m.reference = getReference(phonondb);
                // zinc blende files ship with the page, the others come from the upstream site
                m.url = (m.family === "zincblende" ? self.phonondbRoot : phonondb.root)+"/mp-"+entry.id+".json.gz";
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
        // the base list shows the formula only: tag every family but zinc blende,
        // the same compound is often listed in several
        let materials = this.getMenuMaterials();
        let items = this.dom_mat ? this.dom_mat.children() : [];
        for (let i = 0; i < items.length && i < materials.length; i++) {
            let tag = this.getFamilyTag(materials[i]);
            if (tag) {
                items[i].querySelector("a").insertAdjacentHTML("beforeend", tag);
            }
        }
        this.markSelectedMaterials();
        this.renderMaterial2Menu();
    }

    getMenuMaterials() {
        // materials of the material 1 list, in the order they are shown
        return this.materialsIndex
            .filter((material) => this.materialMatchesFilter(material, this.getMaterialFilterTokens()))
            .filter((material) => this.isReferenceEnabled(this.getMaterialReferenceKey(material)))
            .sort(this.compareMaterialsForMenu.bind(this));
    }

    getFamilyTag(material) {
        /*
        the structure family; ternary rocksalts, (anti)fluorites and tilted or polar
        perovskites also name their space group,
        e.g. LiDyS2 is listed in the alpha-NaFeO2 (R-3m) and gamma-LiFeO2
        (I4_1/amd) orderings
        */
        let family = material && material.family ? material.family : "";
        let tag = "";
        if (family.indexOf("wurtzite") === 0) {
            tag = "wurtzite";
        } else if (family.indexOf("tmd-") === 0) {
            // the polytype of a transition-metal dichalcogenide: 2H, 1T, 3R
            tag = family.slice(4);
        } else if (family === "rutile-AB2X6") {
            tag = "trirutile";
        } else if (family === "corundum-ABX3") {
            tag = "ilmenite";
        } else if (family === "cscl") {
            tag = "CsCl type";
        } else if (family === "heusler") {
            tag = "Li<sub>3</sub>Bi type";
        } else if (family === "heusler-A2BX") {
            tag = "Heusler";
        } else if (family === "heusler-A2BX-inverse") {
            tag = "inverse Heusler";
        } else if (family === "antifluorite-ABX-F43m") {
            // MgAgAs type: zinc blende with the octahedral site filled
            tag = "half-Heusler";
        } else {
            let spacegroup = (material.spacegroup || "").replace(/_(\d)/g, "<sub>$1</sub>");
            for (let parent of ["rocksalt", "fluorite", "antifluorite", "perovskite", "antiperovskite", "rutile", "corundum"]) {
                if (family === parent) {
                    // the cubic cell is the default; a polar 5-atom perovskite names its group
                    let cubic = !spacegroup || spacegroup === "Fm-3m" || spacegroup === "Pm-3m" ||
                        spacegroup === "P4<sub>2</sub>/mnm" || spacegroup === "R-3c";
                    tag = parent + (cubic ? "" : " " + spacegroup);
                } else if (family.indexOf(parent + "-") === 0) {
                    tag = parent + (spacegroup ? " " + spacegroup : "");
                }
            }
        }
        return tag ? ' <span class="family-tag">' + tag + '</span>' : "";
    }

    materialMatchesFilter(material, tokens) {
        // a filter token may also name the structure family, e.g. "wurtzite", or a
        // word of its tag ("heusler" finds the full, inverse and half-Heuslers)
        let tag = this.getFamilyTag(material).replace(/<[^>]*>/g, "").toLowerCase();
        return tokens.every((token) =>
            super.materialMatchesFilter(material, [token]) ||
            (token.length > 2 && (material.family || "").indexOf(token) === 0) ||
            (token.length > 2 && tag.indexOf(token.toLowerCase()) !== -1)
        );
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
        let material1 = this.material1 ? this.findMaterial(this.material1.url) : null;
        let materials = this.materialsIndex
            .filter((material) => this.canMix(material1, material))
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
            addItem(utils.format_formula_html(m.name) + this.getFamilyTag(m) + " [" + index + "]", m);
        }
        if (this.dom_mat2_note) {
            this.dom_mat2_note.html(this.getMaterial2Note(material1, referenceIndices));
        }
    }

    getMaterial2Note(material1, referenceIndices) {
        /*
        explain the pairs canMix leaves out: folding a binary into a ternary cell
        needs its force constants, which only the PhononDB entries carry
        */
        let FAMILY_NAMES = {
            zincblende: ["zinc blende", "Chalcopyrites"],
            wurtzite: ["wurtzite", "Ternary wurtzites"],
            rocksalt: ["rocksalt", "Ternary rocksalts"],
            fluorite: ["fluorite", "Ternary fluorites"],
            antifluorite: ["antifluorite", "Ternary antifluorites"],
            perovskite: ["5-atom perovskite", "Tilted perovskites"],
            antiperovskite: ["5-atom antiperovskite", "Tilted antiperovskites"],
            rutile: ["rutile", "Ordered rutiles"],
            corundum: ["corundum", "Ilmenites and ordered corundums"],
            heusler: ["Li<sub>3</sub>Bi-type compound", "Full and inverse Heuslers"],
        };
        let family1 = material1 ? material1.family : "zincblende";
        let parent = FAMILY_NAMES[family1] ? family1 : getParentFamily(family1);
        if (!FAMILY_NAMES[parent]) {
            return "";
        }
        let reference = this.materialsIndex.find((m) => m.family === parent && m.source === "phonondb");
        let index = reference ? referenceIndices.get(this.getMaterialReferenceKey(reference)) : null;
        let phonondb = "PhononDB" + (index ? " [" + index + "]" : "");
        let [binaryName, ternaryName] = FAMILY_NAMES[parent];
        if (family1 === parent) {
            let hasTernaries = this.materialsIndex.some((m) => getParentFamily(m.family) === parent);
            if (!material1 || material1.source === "phonondb" || !hasTernaries) {
                return "";
            }
            return ternaryName + " are listed only when material 1 is a " + phonondb + " "
                + binaryName + ": mixing them requires its force constants.";
        }
        let hasMpBinaries = this.materialsIndex.some((m) => m.family === parent && m.source !== "phonondb");
        if (!hasMpBinaries) {
            return "";
        }
        return "Only " + phonondb + " " + binaryName + "s are listed: mixing with a "
            + (ternaryName.charAt(0).toLowerCase() + ternaryName.slice(1)).replace(/s$/, "") + " requires their force constants.";
    }

    markSelectedMaterials() {
        if (!this.dom_mat || !this.material1) {
            return;
        }
        let url = this.material1.url;
        let sorted = this.getMenuMaterials();
        let items = this.dom_mat.children();
        for (let i = 0; i < items.length && i < sorted.length; i++) {
            items[i].className = sorted[i].url === url ? "selected" : "";
        }
    }

    loadRamanDatabase() {
        // the mixing page shows no Raman spectra (and has no chart to draw them in)
        return Promise.resolve();
    }

    loadExperimentalRaman() {
        return Promise.resolve();
    }

    canMix(material1, material2) {
        /*
        materials of the same structure family mix directly; a binary can be
        folded into the cell of an ordered ternary on its lattice (zinc blende ->
        chalcopyrite, wurtzite -> ternary wurtzites) when its force constants are
        available (PhononDB), the Materials Project files only cover the band path
        */
        let family1 = material1 ? material1.family : "zincblende";
        if (family1 === material2.family) {
            return true;
        }
        let binary = getParentFamily(material2.family) === family1 ? material1
            : getParentFamily(family1) === material2.family ? material2
            : null;
        if (binary) {
            return binary.source === "phonondb";
        }
        // other structures on a common lattice, with vacancies on the sites only one
        // of them fills; a material folded into the other's cell needs force constants
        let world = getWorldMixing(family1 || "zincblende", material2.family || "zincblende");
        if (!world || !material1) {
            return false;
        }
        return !world.folds.length ||
            world.folds.some((which) => (which === 1 ? material1 : material2).source === "phonondb");
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
        // drop a material 2 that cannot be mixed with the new material 1
        let material2 = this.material2 ? this.findMaterial(this.material2.url) : null;
        if (material2 && !this.canMix(this.findMaterial(url), material2)) {
            this.material2 = null;
            this.updateCharacterLabel();
            this.updateMassLabel();
        }
        this.markSelectedMaterials();
        this.renderMaterial2Menu();
        this.refreshMixing(true, callback);
    }

    loadCallback() {
        super.loadCallback();
        if (!this.material2 && this.material1 && this.material1.name) {
            // the menu name (cation first: NaCl, PbTe) rather than the formula stored
            // in the file (ClNa, TePb)
            this.name = utils.format_formula_html(this.material1.name);
            this.updatePage();
        }
    }

    selectMaterial2(material) {
        this.material2 = material ? { url: material.url, name: material.name } : null;
        this.renderMaterial2Menu();
        this.updateCharacterLabel();
        this.updateMassLabel();
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

    async getMaterialsProjectDielectric(url, raw) {
        /*
        Born charges and dielectric tensor of a Materials Project material, from
        the MP DFPT phonon collection; null (long-range term interpolated) if the
        collection cannot be loaded
        */
        let match = /mp-(\d+)\.json\.gz$/.exec(url);
        if (raw.dynamical_matrix || !match) {
            return null;
        }
        try {
            let table = await loadMaterialsProjectDielectric();
            return table.get(Number(match[1])) || null;
        } catch (error) {
            console.log("Materials Project Born charges unavailable, interpolating the LO-TO term:", error);
            return null;
        }
    }

    async getSampledPair() {
        let key = this.material1.url + "::" + this.material2.url;
        if (this.pairCache && this.pairCache.key === key) {
            return this.pairCache;
        }
        let raws = await Promise.all([this.getRawJson(this.material1.url), this.getRawJson(this.material2.url)]);
        let dielectric1 = await this.getMaterialsProjectDielectric(this.material1.url, raws[0]);
        let dielectric2 = await this.getMaterialsProjectDielectric(this.material2.url, raws[1]);
        let endpoint1 = buildMixingEndpoint(raws[0], this.material1, dielectric1);
        let endpoint2 = buildMixingEndpoint(raws[1], this.material2, dielectric2);
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
        let mass = this.massLinked ? "" : ", mass " + Math.round(this.massCharacter * 100) + " %";
        let familyName = this.pairCache ? this.pairCache.sampled.endpoint1.family : "";
        let union = this.pairCache ? this.pairCache.sampled.union : null;
        let family = ["wurtzite", "rocksalt", "fluorite", "antifluorite", "perovskite", "antiperovskite", "rutile", "corundum"]
            .find((name) => familyName.indexOf(name) === 0);
        if (familyName.indexOf("tmd-") === 0) {
            family = familyName.slice(4);
        } else if (familyName === "cscl") {
            family = "CsCl type";
        } else if (familyName === "heusler") {
            family = "Li<sub>3</sub>Bi type";
        } else if (familyName.indexOf("heusler-") === 0) {
            family = "Heusler";
        } else if (familyName === "antifluorite-ABX-F43m") {
            family = "half-Heusler";
        }
        if (union) {
            // two structures on a common cell
            family = union.world === "fcc" ? "common fcc cell" : "common cubic cell";
        }
        family = family ? family + ", " : "";
        return utils.format_formula_html(this.material1.name) +
            " + " + utils.format_formula_html(this.material2.name) +
            " (" + family + percent + " % " + utils.format_formula_html(this.material2.name) + mass + ")";
    }

    async refreshMixing(pairChanged, callback, scheduled = false) {
        if (!this.material1) {
            return;
        }

        if (!this.material2) {
            // plain material 1, as on the phonon page
            this.mixToken += 1;
            this.updateMixingNote(null);
            super.loadURL(this.material1.url_vars, callback);
            return;
        }

        if (!pairChanged && !scheduled) {
            // slider moves: computing a mixture blocks the page, so the moves made
            // meanwhile queue up; run once, from a timer, with the latest values
            // instead of once per move
            if (!this.mixScheduled) {
                this.mixScheduled = true;
                setTimeout(() => {
                    this.mixScheduled = false;
                    if (this.material2) {
                        this.refreshMixing(false, null, true);
                    }
                }, 0);
            }
            return;
        }

        let token = ++this.mixToken;
        // a plain material 1 load still in flight must not replace the mixture
        this.beginLoad();
        let x = this.character;
        let massFraction = this.getMassFraction();
        if (pairChanged) {
            this.startLoadingFeedback(this.material1.name + " + " + this.material2.name);
        }

        try {
            let pair = await this.getSampledPair();
            let data = await computeMixedPhonon(pair.sampled, x, massFraction);
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
            // a mixture of two structures has fewer bands at 0 and 100 % (no vacancy sites)
            this.n = Math.min(this.n, data.eigenvalues[0].length - 1);
            delete this.link;
            if (pairChanged) {
                this.finishLoadingFeedback();
                this.loadCallback();
            } else {
                this.update();
            }
            this.name = this.getMixingTitle();
            this.updateMixingNote(data.nac_mode, pair.sampled);
            this.updateMassLabel();
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

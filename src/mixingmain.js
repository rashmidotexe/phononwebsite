import $ from 'jquery';
import { initPhononPage } from './phononpage.js';
import { MixingWebpage, MixingHighcharts } from './mixingwebpage.js';

initPhononPage(
    {
        json: "data/mixing/phonondb/mp-149.json.gz",
        name: "Si2 [1]",
        link: "https://materialsproject.org/materials/mp-149",
    },
    {
        webpageClass: MixingWebpage,
        highchartsClass: MixingHighcharts,
        setup: function(p) {
            p.setMaterial2List($('#mat2'));
            p.setMaterial2FilterInput($('#materials_filter2'));
            p.setCharacterInput($('#character_range'), $('#character_value'));
            p.setMixingNote($('#mixing_note'));
        },
    }
);

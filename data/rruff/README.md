# Experimental Raman spectra (RRUFF)

Measured Raman spectra from the [RRUFF project](https://rruff.info/), shown
over the computed spectrum on the phonon page. One file per material
(`mp-<id>.json`) plus `index.json`, which lists the materials and is read first
so only materials with spectra are requested.

Please cite: B. Lafuente, R. T. Downs, H. Yang and N. Stone, "The power of
databases: the RRUFF project", in *Highlights in Mineralogical
Crystallography*, T. Armbruster and R. M. Danisi eds., De Gruyter (2016),
pp. 1–30, [doi:10.1515/9783110417104-003](https://doi.org/10.1515/9783110417104-003).

## Regenerate

Download the Raman archives from https://rruff.info/zipped_data_files/raman/
(keep them out of the repository), then from the repository root:

    python python/phononweb/scripts/prepare_rruff.py /path/to/rruff/*.zip

A sample is kept when its ideal chemistry has the same reduced formula as a
website material, both are in the same crystal family, and the conventional
cells agree (lattice lengths within 6 %, volume within 15 %). Materials without
a computed Raman spectrum are kept too, and the page then shows the measured
spectrum on its own; `--only-ramandb` keeps only materials with Raman tensors
in the Computational Raman Database. Each spectrum is the processed RRUFF spectrum
resampled every 2 cm⁻¹ and scaled to integers 0–1000; up to three spectra per
material, best quality and unoriented first.

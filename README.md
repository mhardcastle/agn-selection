# agn-selection

Python code to select AGN from the LoTSS DR2 catalogue as described by
Hardcastle et al 2025, and to make plots for that paper.

Assumes STILTS is available on the command line, otherwise
dependencies are astropy and matplotlib.

Code is released under a CC-BY 4.0 licence, which means derivative works must credit the author. Citing the published paper constitutes sufficient credit.

## selection

Starting from the v1.1 DR2 catalogue, the Drake et al 2024 emission-line catalogue, and the SDSS DR16 quasar catalogue:

* `python selection.py`
* `python class_cm.py`
* `python absmag.py fcoz_class_DR16Q.fits`
* `python apply_cuts.py`
* `python cc_class.py` (for colour-colour classifications)

## plots

* absmag_elc.py
* absmag_lum_w2_density.py
* absmag_lum_w2_spc_density.py
* absmag_lum_w3_density.py
* absmag_lum_w3_spc_density.py
* cc_class.py
* core_alpha_cdist.py
* core_prom_lum_with3c.py
* core_vi_cdist.py
* cores_cdist.py
* fit_lf.py
* make_agn_sf_hists.py
* make_ccc_plots2.py
* make_pdd.py
* numbercounts.py
* plot_agn_mass_lum.py
* plot_agn_mass_size.py
* plot_lf.py
* plot_lf_fits.py
* plot_lf_ms.py
* plots.py
* typecount.py
* willott.py
* wisecc_class_density.py
* wisecc_class_density_elc.py


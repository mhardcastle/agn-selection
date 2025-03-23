# agn-selection

Python code to select AGN from the LoTSS DR2 catalogue as described by
Hardcastle et al 2025, and to make plots for that paper.

Assumes STILTS is available on the command line, otherwise
dependencies are astropy and matplotlib.

## selection

Starting from the v1.1 DR2 catalogue, the Drake et al 2024 emission-line catalogue, and the SDSS DR16 quasar catalogue:

* `python selection.py`
* `python class_cm.py`
* `python absmag.py fcoz_class_DR16Q.fits`
* `python apply_cuts.py`



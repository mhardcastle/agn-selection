from astropy.table import Table
import numpy as np

t=Table.read('classified.fits')
print('    Original sample&%i&\\' % len(t))
print('    Star formation exclusion&%i&%i\\' % (np.sum(t['SF_EXCLUDE']),np.sum(t['SF_EXCLUDE_BROAD'])))
print('    RQQ exclusion&%i&%i\\' % (np.sum(t['RQQ_EXCLUDE']),np.sum(t['RQQ_EXCLUDE_BROAD'])))
print('    Remaining AGN&%i&%i\\' % (np.sum(t['AGN_BROAD']),np.sum(t['AGN_NARROW'])))

import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
from plots import *

from cores_cdist import cdf

t=Table.read('large_cores_with_vi_ccc.fits')
#t=t[t['Core_3000']>0.005]
#t=t[t['Core_3000']<0.020]
setupfig()
colors={'DR16Q':ccol(3),'N':ccol(5),'RI':ccol(12),'Q':ccol(7),'RE':ccol(8)}

t=t[~np.isnan(t['Variability_index'])]
print(len(t),'with VI')

cdf(t['Variability_index'],label='All sources (%i)' % len(t),bootstrap=100,color='black',lw=2)
for label in ['RE','Q','N','RI','DR16Q']:
    st=t[t[label]]
    cdf(st['Variability_index'],label=label+' (%i)' % len(st),bootstrap=100,color=colors[label])

plt.xlim(0.5,10)
plt.ylim(0.5,1)
plt.grid()
plt.legend(loc=0)
plt.xlabel('Core variability index $\chi^2_{r,VI}$')
savefig('core_vi_cdist.pdf')
plt.show()


import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
from plots import *
import sys

from cores_cdist import cdf

try:
    first=sys.argv[1]=='FIRST'
except:
    first=False
        
t=Table.read('large_cores_ccc.fits')
#t=t[t['Core_3000']>0.005]
#t=t[t['Core_3000']<0.020]
setupfig()
colors={'DR16Q':ccol(3),'N':ccol(5),'RI':ccol(12),'Q':ccol(7),'RE':ccol(8)}

if first:
    k='Core_alpha_1400_3000'
else:
    k='Core_alpha_144_3000'

t=t[~np.isnan(t[k])]
print(len(t),'with',k)

cdf(-t[k],label='All sources (%i)' % len(t),bootstrap=100,color='black',lw=2)
for label in ['RE','Q','N','RI','DR16Q']:
    st=t[t[label]]
    cdf(-st[k],label=label+' (%i)' % len(st),bootstrap=100,color=colors[label])

#plt.xlim(-3.5,0)
plt.xlabel(k)
plt.legend(loc=0)
plt.grid()
if first:
    plt.xlabel('Core spectral index $\\alpha_{1400}^{3000}$')
    plt.xlim(-2,2)
    savefig('core_alpha_cdist_1400_3000.pdf')
else:
    plt.xlabel('Core spectral index $\\alpha_{144}^{3000}$')
    savefig('core_alpha_cdist_144_3000.pdf')
plt.show()


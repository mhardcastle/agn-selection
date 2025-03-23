from astropy.table import Table
from astropy.io import fits

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import astropy.units as u
import numpy as np
import sys

# compute the power-law spectral index based on https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html

t=Table.read(sys.argv[1])

z=t['z_best']
ld=cosmo.luminosity_distance(z).value*1e6

zero=[2.699,3.339,5.174]
fc=[309.540,171.787,31.674]

for i in range(3):
    t['flux_w%i' % (i+1)]=fc[i]*10**-((t['mag_w%i' % (i+1)]-zero[i])/2.5)

t['alpha_w1_w2']=np.log10(t['flux_w1']/t['flux_w2'])/np.log10(3.4/4.6) # usual convention, steep spectrum is more positive -- because 3.4 and 4.6 are wavelengths.
mu=2.5*np.log10(ld**2.0*(1+t['z_best'])**(t['alpha_w1_w2']-1))-5
t['mu']=np.array(mu)*u.mag
t['abs_w2']=t['mag_w2']-t['mu']

t['alpha_w2_w3']=np.log10(t['flux_w2']/t['flux_w3'])/np.log10(4.6/12.1) # usual convention, steep spectrum is more positive -- because 3.4 and 4.6 are wavelengths.
mu=2.5*np.log10(ld**2.0*(1+z)**(t['alpha_w2_w3']-1))-5
t['mu']=np.array(mu)*u.mag
t['abs_w3']=t['mag_w3']-t['mu']

t.write(sys.argv[1].replace('.fits','_absmag.fits'),overwrite=True)

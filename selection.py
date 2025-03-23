from astropy.table import Table
import numpy as np

t=Table.read('combined-release-v1.1-LM_opt_mass.fits')

filt=t['Total_flux']>1.1
tfc=np.sum(filt)
print('After flux cut:',tfc)
filt&=t['z_best']>0.01
print('After z_best filt',np.sum(filt))

print('Now there are',np.sum(t[filt]['WISE_Src']=='unwise2019'),'unwise sources')

for band in range(1,4):
    filt&=~t['mag_w'+str(band)].mask
    filt&=~np.isnan(t['mag_w'+str(band)])
    print('After band %i filt: %i' % (band,np.sum(filt)))

t=t[filt]
print(len(t))
print(100.0*len(t)/tfc)

t.write('fcoz.fits',overwrite=True)

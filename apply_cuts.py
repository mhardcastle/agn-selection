from astropy.table import Table
import numpy as np
import sys

try:
    infile=sys.argv[1]
    outroot=sys.argv[1].replace('.fits','_')
except:
    infile='fcoz_class_DR16Q_absmag.fits'
    outroot=''

t=Table.read(infile)

def divide(m):
    return -m/2.5+14

def qline(m):
    return 25.3+(-m-27)*2.0/7

logl=np.log10(t['L_144'])

sf_exclude=logl<divide(t['abs_w3'])
sf_exclude&=logl<24.8

t['SF_EXCLUDE_BROAD']=sf_exclude
t['SF_EXCLUDE_BROAD'].mask=0
sf_exclude&=~t['magerr_w3'].mask

t['SF_EXCLUDE']=sf_exclude
t['SF_EXCLUDE'].mask=0

print('SF cut excludes',np.sum(sf_exclude))
t[sf_exclude].write(outroot+'sf.fits',overwrite=True)

logl=np.log10(t['L_144'])
      
rqq=t['abs_w3']<-27
rqq&=logl>=24.8
rqq&=logl<qline(t['abs_w3'])
t['RQQ_EXCLUDE_BROAD']=rqq
rqq&=~t['magerr_w3'].mask
t['RQQ_EXCLUDE']=rqq

print('RQQ exclude',np.sum(rqq))

t[rqq].write(outroot+'rqq.fits',overwrite=True)
#t[~rqq].write('norqq.fits',overwrite=True)

t['AGN_NARROW']=~(t['SF_EXCLUDE_BROAD'] | t['RQQ_EXCLUDE_BROAD'])
t['AGN_BROAD']=~(t['SF_EXCLUDE'] | t['RQQ_EXCLUDE'])
t.write(outroot+'classified.fits',overwrite=True)
t=t[t['AGN_BROAD']]
t.write(outroot+'agn.fits',overwrite=True)

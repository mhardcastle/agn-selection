from astropy.table import Table
import os
import numpy as np

tf=Table.read('fcoz.fits')

if not os.path.isfile('classifications_stripped.fits'):
    t=Table.read('Classifications_ABD_v1.2.fits')

    filt=t['CLASS_z_WARNING']==0 # else can't believe lines
    t=t[filt]

    for k in tf.colnames:
        if k in t.colnames and k!='Source_Name':
            del(t[k])

    t.write('classifications_stripped.fits',overwrite=True)
    print(len(t))

os.system('stilts tmatch2 in1=fcoz.fits in2=classifications_stripped.fits out=fcoz_classifications.fits matcher=exact values1=Source_Name values2=Source_Name join=all1')

t=Table.read('fcoz_classifications.fits')
t['ELC']=~t['Source_Name_2'].mask
print(np.sum(t['ELC']))
del(t['Source_Name_2'])
del(t['RA_1'])
del(t['DEC_1'])
t['Source_Name_1'].name='Source_Name'
t['RA_1b'].name='RA'
t['DEC_1b'].name='DEC'
t['ra_1c'].name='ra'
t['dec_1c'].name='dec'
t['ra_1a'].name='SDSS_ELC_RA'
t['dec_1a'].name='SDSS_ELC_DEC'

t.write('fcoz_classifications.fits',overwrite=True)

os.system('stilts tmatch2 in1=fcoz_classifications.fits in2=DR16Q_v4.fits matcher=sky params=1 values1="ID_RA ID_DEC" values2="RA DEC" join=all1 out=fcoz_class_DR16Q.fits')

t=Table.read('fcoz_class_DR16Q.fits')
t['DR16Q']=~t['Separation'].mask
t['RA_1'].name='RA'
t['DEC_1'].name='DEC'
t['RA_2'].name='SDSS_RA'
t['DEC_2'].name='SDSS_DEC'
t['ra_1a'].name='ra'
t['dec_1a'].name='dec'
t['objid_1'].name='objid'
t['OBJID_2'].name='SDSS_OBJID'

t.write('fcoz_class_DR16Q.fits',overwrite=True)

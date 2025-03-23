from astropy.table import Table
import sys
import numpy as np

t=Table.read(sys.argv[1])

corrs=[2.699,3.339,5.174,6.620]
for i in range(4):
    band='mag_w%i' % (i+1)
    vband=band+'v'
    t[vband]=t[band]-corrs[i]

t['w2w3']=t['mag_w2v']-t['mag_w3v']
t['w1w2']=t['mag_w1v']-t['mag_w2v']

print('Total objects',len(t))
#t['RE']=(t['abs_w3']<-24.5) & (t['magerr_w3']>0)
#t['RE'].mask=np.zeros_like(t['RE'],dtype=bool)
#t['RI']=(t['abs_w3']>-24.5)
t['RE']=((t['w2w3']>2) & (t['magerr_w3']>0)) | (t['w1w2']>0.4)
t['RE'].mask=0
t['RI']=(t['w2w3']<2) & (t['w1w2']<0.4)
#t['RI'].mask=0
print('RE is',np.sum(t['RE']))
print('RI is',np.sum(t['RI']))
neither=(~t['RE']) & (~t['RI'])
t['U']=neither
t['U'].mask=0
print('Unclassified is',np.sum(t['U']))
t['Q']=t['RE'] & (t['w1w2']>0.75)
t['N']=t['RE'] & (t['w1w2']<=0.75)
#print(np.sum(t['RE']),np.sum(t['RI']),np.sum(t['U']),np.sum(t['Q']),np.sum(t['N']))
t.write(sys.argv[1].replace('.fits','_ccc.fits'),overwrite=True)

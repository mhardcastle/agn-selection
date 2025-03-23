import glob
import matplotlib.pyplot as plt
import numpy as np
from plots import *
from fit_lf import model
from willott import *

setupfig()

ms_agn=np.loadtxt('/home/mjh/lofar/catalogues/ms-agn.txt')
ms_sf=np.loadtxt('/home/mjh/lofar/catalogues/ms-sf.txt')
cols=[ccol(11),ccol(7)]
vals=[]
g=sorted(glob.glob('LF*0.01-0.3-W1.csv'))
for i,f in enumerate(g):
    label=f.replace('LF_','').split('-')[0]
    if label=='Total': continue
    d=np.loadtxt(f,delimiter=',')
    d[:,1]*=5700/5200
    d[:,2]*=5700/5200
    #d[:,1]=np.where(d[:,1]>0.001,np.nan,d[:,1])
    #d[:,1]=np.where(d[:,1]<1.5*d[:,2],np.nan,d[:,1])
    d=d[~np.isnan(d[:,1])]
    plt.errorbar(d[:,0],np.log10(d[:,1]),yerr=np.log10(d[:,1]+d[:,2])-np.log10(d[:,1]),label=label,color=cols[i],marker='o')

scale=0.7*np.log10(1400.0/150)

plt.errorbar(ms_agn[:,0]+scale,ms_agn[:,1]+np.log10(2.5),yerr=(ms_agn[:,2],ms_agn[:,3]),color=cols[0],ls='--',label='MS AGN')
plt.errorbar(ms_sf[:,0]+scale,ms_sf[:,1]+np.log10(2.5),yerr=(ms_sf[:,2],ms_sf[:,3]),color=cols[1],ls='--',label='MS SF')
plt.legend(loc=0,ncol=2)
ax=plt.gca()
set_logx(ax)
set_logy(ax)
plt.xlim(21,29)
plt.ylim(-9.3,-2)
plt.xlabel('$L_{144}$ (W Hz$^{-1}$)')
plt.ylabel('$\\rho$ (Mpc$^{-3}$ dex$^{-1}$)')
savefig('lf-ms.pdf')
#plt.show()

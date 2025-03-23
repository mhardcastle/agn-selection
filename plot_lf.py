import glob
import matplotlib.pyplot as plt
import numpy as np
from plots import *
from fit_lf import model
from willott import *

setupfig()
overplot_willott=False

colors=[ccol(2),ccol(3),ccol(11),ccol(4),ccol(5),ccol(6),ccol(8),ccol(7),ccol(9),ccol(10)]

ms_agn=np.loadtxt('/home/mjh/lofar/catalogues/ms-agn.txt')

vals=[]
g=sorted(glob.glob('LF*W1-rebin-fr_corr-samples.npy'))
g=[f.replace('-samples.npy','.csv') for f in g]
for i,f in enumerate(g):
    zrange=f.replace('.csv','').replace('LF_AGN-','').split('-')
    if zrange[0]=='0.0': zrange[0]='0.01'
    label='$%s < z < %s$' % (zrange[0],zrange[1])  

    d=np.loadtxt(f,delimiter=',')
    d[:,1]=np.where(d[:,1]>0.001,np.nan,d[:,1])
    d[:,1]=np.where(d[:,1]<1.5*d[:,2],np.nan,d[:,1])
    d=d[~np.isnan(d[:,1])]
    plt.errorbar(d[:,0],np.log10(d[:,1]),yerr=np.log10(d[:,1]+d[:,2])-np.log10(d[:,1]),label=label,color=colors[i],marker='o',ls='none')
    samplefile=f.replace('.csv','-samples.npy')
    samples=np.load(samplefile)
    ss=np.random.choice(len(samples),100)
    xv=np.linspace(np.min(d[:,0]),np.max(d[:,0]),100)
    if not overplot_willott:
        for j in range(len(ss)):
            if samples[ss[j],2]>0.7: continue
            plt.plot(xv,np.log10(model(xv,samples[ss[j],0],samples[ss[j],1],samples[ss[j],2],samples[ss[j],3])),color=colors[i],alpha=0.07)
    if overplot_willott:
        wlf=rho(10**xv, (float(zrange[1])+float(zrange[0]))/2, **params_model_C_Omega0)
        plt.plot(xv,np.log10(wlf),color=colors[i])

scale=0.7*np.log10(1400.0/150)

#plt.errorbar(ms_agn[:,0]+scale,ms_agn[:,1]+np.log10(2.5),yerr=(ms_agn[:,2],ms_agn[:,3]),color='black',ls='--',label='Mauch \& Sadler')
plt.legend(loc=0,ncol=2)
ax=plt.gca()
set_logx(ax)
set_logy(ax)
plt.xlim(21,29)
plt.ylim(-9.0,-3.4)
plt.xlabel('$L_{144}$ (W Hz$^{-1}$)')
plt.ylabel('$\\rho$ (Mpc$^{-3}$ dex$^{-1}$)')
if overplot_willott:
    savefig('lfwillott.pdf')
else:
    savefig('lfdata.pdf')
#plt.show()

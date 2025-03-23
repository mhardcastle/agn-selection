import glob
import matplotlib.pyplot as plt
import numpy as np
from plots import *
from astropy.cosmology import FlatLambdaCDM,z_at_value
import astropy.units as u

colors=[ccol(2),ccol(3),ccol(11),ccol(4),ccol(5),ccol(6),ccol(8),ccol(7),ccol(9),ccol(10)]

cosmo=FlatLambdaCDM(H0=70, Om0=0.3)
setupfig(size=(18,12))
g=sorted(glob.glob('LF*W1-rebin-fr_corr.csv'))
labels=['logC','$L_*$ (W Hz$^{-1}$)','$\\alpha$','$\\beta$','$\\rho_{25}$ (Mpc$^{-3}$ dex$^{-1}$)']
outd=[[],[],[],[],[]]
for f in g:
    zrange=f.replace('.csv','').replace('LF_AGN-','').split('-')
    zmid=(float(zrange[0])+float(zrange[1]))/2
    samples=np.load(f.replace('.csv','-samples.npy'))
    meanparms=np.median(samples,axis=0)
    for i in range(5):
        errors=np.percentile(samples[:,i],(16,84))
        merr=((meanparms[i]-errors[0])+(errors[1]-meanparms[i]))/2
        outd[i].append([zmid,meanparms[i],merr])
        #print(labels[i],meanparms[i],np.percentile(samples[:,i],(16,84)))

d=np.array(outd)
ages=np.array([1,2,3,4,5,6,7,8])*u.Gyr
zvals=z_at_value(cosmo.lookback_time,ages)

for n,i in enumerate([4,1,2,3]):
    plt.subplot(2,2,n+1)
    #plt.subplots_adjust(vspace=0.01)
    ax=plt.gca()
    sd=d[i]
    for j in range(len(sd)):
        plt.errorbar(sd[j,0],sd[j,1],yerr=sd[j,2],marker='o',ls='none',c=colors[j])
    if n==0 or n==1:
        ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
        ax.xaxis.set_ticks(zvals)
        ax.xaxis.set_ticklabels(ages.value)
        plt.xlabel('$t_\mathrm{lookback}$ (Gyr)',labelpad=10)
        ax.xaxis.set_label_position('top')
    else:
        plt.xlabel('$z$')
    plt.ylabel(labels[i])
    if i==2:
        plt.plot([0,1.2],[0.49,0.49],color='black',ls=':')
        plt.xlim([0,1.2])
    if i==3:
        plt.plot([0,1.2],[1.27,1.27],color='black',ls=':')
        plt.xlim([0,1.2])
    if i==1 or i==4:
        set_logy(plt.gca())
    if i==1:
        z=np.linspace(0,1.2,100)
        exponent=8
        plt.plot(z,25.5+exponent*np.log10(1+z),color=ccol(12),label=f'$(1+z)^{exponent}$')
        plt.legend(loc=0)
    if i==4:
        z=np.linspace(0,1.2,100)
        plt.plot(z,-5.5+3.48*np.log10(1+z),color=ccol(13),label='Willott+ 01')
        # Pure density evolution from Wang+ 24
        A_D=-0.77
        B_D=2.68
        alpha_D=A_D*z+B_D
        plt.plot(z,-5.4+alpha_D*np.log10(1+z),color=ccol(14),label='Wang+ 24')
        plt.legend(loc=0)

savefig('parameters.pdf')
plt.show()

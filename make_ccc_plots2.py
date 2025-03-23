from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from plots import *
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import astropy.units as u
from astropy.stats import binom_conf_interval

def errorbar(mp,k,n,**kwargs):
    mpp=mp[n>4]
    npp=n[n>4]
    kp=k[n>4]
    kwargs['markersize']=10
    ci=binom_conf_interval(kp,npp,interval='jeffreys')
    yerr=np.abs(ci-kp/npp)
    plt.errorbar(mpp,kp/npp,yerr=yerr,**kwargs)
    plt.fill_between(mpp,ci[0],ci[1],alpha=0.3,color=kwargs['color'])

setupfig(size=(22,10))
cols={'RE':ccol(5),'RI':ccol(12),'U':'gray','Q':ccol(3)}
ccols={'RE':ccol(8),'RI':ccol(13),'U':ccol(2),'Q':ccol(7)}

t=Table.read('agn_ccc.fits')
t=t[t['z_best']<1.2]
zv=np.log10(1+t['z_best'])
lv=np.log10(t['L_144'])
t['zv']=zv
t['lv']=lv

t3c=Table.read('3crr_2.txt',format='ascii')
t3c=t3c[t3c['z']<1.2]
t3c['L144']=t3c['L']*4*np.pi*(178/144)**0.7
t3c=t3c[t3c['L144']>1e24]
t3c['RI']=t3c['Type']==0
t3c['RE']=t3c['Type']>0
t3c['Q']=t3c['Type']>1
t3c['zv']=np.log10(1+t3c['z'])
t3c['lv']=np.log10(t3c['L144'])

td24=Table.read('Classifications_ABD_v1.2.fits')
td24=td24[td24['RADIO_EXCESS']>0.95]
td24['RI']=td24['CLASS_LINELERG']>0.8
td24['RE']=td24['CLASS_HERG']>0.8
td24['U']=~(td24['RI'] | td24['RE'])
td24['Q']=False
td24['zv']=np.log10(1+td24['z_best'])
td24['lv']=np.log10(td24['L_144'])

tb23=Table.read('en1_classifications_dr1.fits')
tb23=tb23[tb23['RadioAGN_final']==1]
tb23=tb23[tb23['z_best']<1.2]
tb23['RI']=tb23['Overall_class']=='LERG'
tb23['RE']=tb23['Overall_class']=='HERG'
tb23['U']=tb23['Overall_class']=='Unc'
tb23['zv']=np.log10(1+tb23['z_best'])
ld=cosmo.luminosity_distance(tb23['z_best'])
lr=4*1e-26*np.pi*tb23['S_150MHz']*ld.to(u.m)**2.0*(1+tb23['z_best'])**(-0.3)
tb23['L_144']=lr
tb23['lv']=np.log10(lr.value)

v='lv'
bins=np.linspace(20.25,29.75,20)
n,bins=np.histogram(t[v],bins=bins)
mpf=0.5*(bins[0:-1]+bins[1:])
#plt.plot(mp,n,label='All',color='black')
d={}
d['All']=n
for k in ['RE','RI','U','Q']:
    d[k],_=np.histogram(t[v][t[k]],bins=bins)

for subplot in range(1,4):

    plt.subplot(1,3,subplot)
    #    plt.plot(mp,d[k],label=k,color=cols[k])
    #plt.yscale('log')
    #plt.legend(loc=0)
    #plt.subplot(2,2,2+i*2)
    for k in ['RE','RI','U','Q'] if subplot==1 else ['RE','RI','U']:
        errorbar(mpf,d[k],n,label=k,color=cols[k],marker='o')

    # switch to coarse bins
    bins=np.linspace(20.5,29.5,10)
    mp=0.5*(bins[0:-1]+bins[1:])

    if subplot==1:
        plt.title('(a) Comparison with 3CRR')
        n3c,_=np.histogram(t3c[v],bins=bins)
        for k in ['RE','RI','Q']:
            n3cs,_=np.histogram(t3c[v][t3c[k]],bins=bins)
            errorbar(mp,n3cs,n3c,label='3C '+k,ls='--',color=ccols[k],marker='*')
    elif subplot==2:
        plt.title('(b) Comparison with D24')
        nd24,_=np.histogram(td24[v],bins=bins)
        for k in ['RE','RI','U']:
            nd24s,_=np.histogram(td24[v][td24[k]],bins=bins)
            errorbar(mp,nd24s,nd24,label='D24 '+k,ls=':',color=ccols[k],marker='v')
    elif subplot==3:
        plt.title('(c) Comparison with Best et al (2023)')
        nb23,_=np.histogram(tb23[v],bins=bins)
        for k in ['RE','RI']:
            nb23s,_=np.histogram(tb23[v][tb23[k]],bins=bins)
            errorbar(mp,nb23s,nb23,label='B23 '+k,ls='-.',color=ccols[k],marker='s')

    plt.legend(loc=0,ncols=3) # was in loop
    if subplot==1: plt.ylabel('Fraction')
    plt.xlabel('$L_{144}$ (W Hz$^{-1}$)')
    set_logx(plt.gca())
    plt.xlim(21,29)
    plt.ylim(0,1.25)
    
savefig('ccc_comparison.pdf')

    

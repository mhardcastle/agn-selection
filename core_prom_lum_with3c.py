import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
from plots import *
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
colors={'DR16Q':ccol(3),'N':ccol(5),'RI':ccol(12),'Q':ccol(3),'RE':ccol(8)}

rmsf=open('rms.txt').readlines()
rms={}
for l in rmsf:
    bits=l.rstrip().split(',')
    rms[bits[0]]=float(bits[1])

names=[]
prominence=[]
z=[]
re=[]
q=[]
n=[]
ri=[]
lum=[]

lines=open('mullin_cores.txt').readlines()
for l in lines:
    bits=l.rstrip().split()
    if bits[5].startswith('<'): continue
    names.append(bits[0])
    z.append(float(bits[1]))
    flux=float(bits[3])
    alpha=float(bits[4])
    ld=cosmo.luminosity_distance(z[-1])
    lr=4*1e-26*np.pi*flux*ld.to(u.m)**2.0*(1+z[-1])**(alpha-1)
    lr*=(178.0/144.0)**alpha
    lum.append(lr)
    ty=bits[2]
    re.append(ty in ['N','B','Q'])
    ri.append(ty=='E')
    n.append(ty=='N')
    q.append(ty in ['B','Q'])
    prominence.append(np.log10(float(bits[5])))
    

tm=Table([names,z,lum,prominence,re,ri,n,q],names=['Source_Name','z','L_144','prominence','RE','RI','N','Q'])


t=Table.read('large_cores_ccc.fits')
setupfig()
#t=t[t['Core_3000']>0]
#print(len(t))
t['prominence']=3+np.log10(np.array(t['Core_3000'])/t['Total_flux'])
plimits=[]
for r in t:
    if np.isnan(r['prominence']):
        plimits.append(3+np.log10(5*rms[r['Source_Name']]/r['Total_flux']))
    else:
        plimits.append(np.nan)
t['plimit']=plimits
t.write('temp.fits',overwrite=True)

for name,table,marker,alpha,size in [('RLAGN',t,'o',0.3,10),('3CRR',tm,'*',1.0,100)]:

    for label in ['RE','RI','Q']:
        st=table[table[label]]
        print(len(st))
        plt.scatter(np.log10(st['L_144']),st['prominence'],marker=marker,alpha=alpha,label=name+' '+label,color=colors[label],s=size)
        if 'plimit' in table.colnames:
            plt.scatter(np.log10(st['L_144']),st['plimit'],marker=r'$\downarrow$',alpha=alpha,color=colors[label],s=size*4)

set_logx(plt.gca())
plt.xlim(26,29)
plt.xlabel('Luminosity (W Hz$^{-1}$)')
plt.ylabel('log$_{10}$(Core prominence at 3 or 8 GHz)')
plt.legend(loc=0)
savefig('core_prom_lum_with3c.pdf')
plt.show()


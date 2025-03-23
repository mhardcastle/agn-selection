from plots import *
from matplotlib.patches import Ellipse
from density_contours2 import DC

t=Table.read('fcoz_class_DR16Q_absmag.fits')
filt=~np.isnan(t['abs_w3'])
filt&=~np.isnan(t['L_144'])
#filt&=~np.isnan(t['magerr_w3'])
t=t[filt]

print(len(t))
print(np.sum(t['magerr_w3'].mask))
    
setupfig()
extent=(-18,-34,20.5,29)
plt.hexbin(t['abs_w3'],np.log10(t['L_144']),mincnt=5,gridsize=100,alpha=1.0,bins='log',cmap='Greys',extent=extent,label=None)


sf=t[t['CLASS_SFG']>0.95]
filt=t['RADIO_EXCESS']>0.95
filt&=~(t['r_50']>8) # remove large
re=t[filt]
rq=t[t['CLASS_RQAGN']>0.95]
hl=t[t['L_144']>1e25]

print('RE',len(re),'SF',len(sf),'HL',len(hl))

#plt.scatter(rq['w3w3'],rq['w1w3'],color='yellow',s=2,alpha=0.2)
d=DC()
d.density_contours(sf['abs_w3'],np.log10(sf['L_144']),*extent,ccol(7),0.8,12,label='D24 SFG',legend=True,scatter=200)
d.density_contours(re['abs_w3'],np.log10(re['L_144']),*extent,ccol(11),0.8,12,label='D24 RXG',legend=True,scatter=200)
d.do_scatter()

#plt.xlim((-0.5,7))
#plt.ylim((-1,2.5))

plt.xlabel('Absolute $W3$ magnitude')
ax=plt.gca()
set_logy(ax)
plt.ylabel('144-MHz radio luminosity (W Hz$^{-1}$)')

def divide(m):
    return -m/2.5+14

plt.plot([-18,-34],divide(np.array([-18,-34])),color='black',label='SF exclusion')
ax.add_patch(Ellipse((-27,24.85),2.7,0.7,0,label="`Blob'",ls='--',edgecolor=ccol(14),facecolor='none'))

plt.ylim(extent[2:4])
plt.xlim(extent[0:2])
plt.legend(loc=0,frameon=False)

savefig('absmag_lum_w3_density.pdf')


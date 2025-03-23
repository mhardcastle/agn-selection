from plots import *
from matplotlib.patches import Ellipse
from density_contours2 import DC

t=Table.read('fcoz_class_DR16Q_absmag.fits')
filt=~np.isnan(t['abs_w2'])
filt&=~np.isnan(t['L_144'])
t=t[filt]
print(len(t))
    
setupfig()
extent=(-16,-32,20.5,29)
plt.hexbin(t['abs_w2'],np.log10(t['L_144']),mincnt=5,gridsize=100,alpha=1.0,bins='log',cmap='Greys',extent=extent,label=None)

sf=t[t['CLASS_SFG']>0.95]
filt=t['RADIO_EXCESS']>0.95
filt&=~(t['r_50']>8) # remove large
re=t[filt]
#rq=t[t['CLASS_RQAGN']>0.95]
#hl=t[t['L_144']>1e25]

print('RE',len(re),'SF',len(sf))

#plt.scatter(rq['w2w3'],rq['w1w2'],color='yellow',s=2,alpha=0.2)
d=DC()
d.density_contours(sf['abs_w2'],np.log10(sf['L_144']),*extent,ccol(7),0.8,12,label='D24 SFG',legend=True,scatter=200)
d.density_contours(re['abs_w2'],np.log10(re['L_144']),*extent,ccol(11),0.8,12,label='D24 RXG',legend=True,scatter=200)
d.do_scatter()

#plt.xlim((-0.5,7))
#plt.ylim((-1,2.5))

plt.xlabel('Absolute $W2$ magnitude')
set_logy(plt.gca())
plt.ylabel('144-MHz radio luminosity (W Hz$^{-1}$)')
plt.ylim(extent[2:4])
plt.xlim(extent[0:2])

ax=plt.gca()
ax.add_patch(Ellipse((-24.5,24.85),2.4,0.7,0,label="`Blob'",ls='--',edgecolor=ccol(14),facecolor='none'))
handles, labels = ax.get_legend_handles_labels()
newArtist1 = plt.Line2D((0,1),(0,0), color='gray', marker='h', markersize=8, alpha=0.2, linestyle='')
newlabel=['All sources']
leg=plt.legend([newArtist1]+handles,newlabel+labels,loc=0,numpoints=1,prop={'size':20}, frameon=False,markerscale=1.5)
for i,lh in enumerate(leg.legendHandles): 
    if i>0:
        lh.set_alpha(1)


savefig('absmag_lum_w2_density.pdf')


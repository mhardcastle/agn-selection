from plots import *
from matplotlib.patches import Ellipse
from density_contours2 import DC

t=Table.read('fcoz_class_DR16Q_absmag.fits')
#t=t[~np.isnan(t['magerr_w3'])]
filt=~np.isnan(t['abs_w2'])
filt&=~np.isnan(t['L_144'])
#filt&=~np.isnan(t['magerr_w3'])
t=t[filt]

q=t[t['DR16Q']]
print(len(t))
    
setupfig()
extent=(-17,-32,20.5,29)
plt.hexbin(t['abs_w2'],np.log10(t['L_144']),mincnt=5,gridsize=100,alpha=1.0,bins='log',cmap='Greys',extent=extent,label=None)


herg=t[t['CLASS_HERG']>0.80]
lerg=t[(t['CLASS_LINELERG']>0.80) & (t['r_50']<8)]

print('HERG',len(herg),'LERG',len(lerg),'Q',len(q))
d=DC()
d.density_contours(q['abs_w2'],np.log10(q['L_144']),*extent,ccol(3),0.8,12,label='DR16Q',legend=True,scatter=200)
d.density_contours(lerg['abs_w2'],np.log10(lerg['L_144']),*extent,ccol(12),0.8,12,label='D24 LERG',legend=True,scatter=200)
d.density_contours(herg['abs_w2'],np.log10(herg['L_144']),*extent,ccol(5),0.8,12,label='D24 HERG',legend=True,scatter=200)
d.do_scatter()

#nd=t[np.isnan(t['magerr_w3'])]
#plt.scatter(nd['abs_w2'],np.log10(nd['L_144']),color='yellow',s=2,alpha=0.3,label='ND')


'''
Limit code
z=np.arange(0,4,0.1)
w3s=[]
l144s=[]
for i in range(len(z)-1):
    zmin=z[i]
    zmax=z[i+1]
    filt=(nd['z_best']>zmin) & (nd['z_best']<zmax)
    c=nd[filt]
    if len(c)==0: continue
    w3=np.nanmedian(c['abs_w2'])
    l144=np.min(c['L_144'])
    print(zmin,zmax,len(c),w3,l144)
    w3s.append(w3)
    l144s.append(l144)
    
plt.scatter(w3s,np.log10(np.array(l144s)),color='orange',marker='+',label='Limit line')
'''

#plt.xlim((-0.5,7))
#plt.ylim((-1,2.5))

plt.xlabel('Absolute $W2$ magnitude')
ax=plt.gca()
set_logy(ax)
plt.ylabel('144-MHz radio luminosity (W Hz$^{-1}$)')

def divide(m):
    return -m/2.5+14


ax.add_patch(Ellipse((-24.5,24.85),2.4,0.7,0,label="`Blob'",ls='--',edgecolor=ccol(14),facecolor='none'))

plt.ylim(20.5,29)
plt.xlim(-17,-32)
ax=plt.gca()
handles, labels = ax.get_legend_handles_labels()
newArtist1 = plt.Line2D((0,1),(0,0), color='gray', marker='h', markersize=8, alpha=0.2, linestyle='')
newlabel=['All sources']
leg=plt.legend([newArtist1]+handles,newlabel+labels,loc=0,numpoints=1,prop={'size':20}, frameon=False,markerscale=1.5)
for i,lh in enumerate(leg.legendHandles): 
    if i>0:
        lh.set_alpha(1)

#plt.show()

savefig('absmag_lum_w2_spc_density.pdf')

from plots import *
from scipy import stats
from density_contours2 import DC

t=Table.read('fcoz_class_DR16Q_absmag.fits')
#t=Table.read('nosf.fits')
for i in range(3):
    band='mag_w%i' % (i+1)
    err=band.replace('mag','magerr')
    filt=~np.isnan(t[err])
    #    filt&=t[err]<0.3
    print('band',i+1,'reducing to',np.sum(filt))
    t=t[filt]

print(len(t))
    
setupfig()
corrs=[2.699,3.339,5.174,6.620]
table=t
for i in range(4):
    band='mag_w%i' % (i+1)
    vband=band+'v'
    table[vband]=table[band]-corrs[i]

t['w2w3']=t['mag_w2v']-t['mag_w3v']
t['w1w2']=t['mag_w1v']-t['mag_w2v']
extent=(0,5.5,-0.4,2.4)
    
plt.hexbin(t['w2w3'],t['w1w2'],mincnt=5,gridsize=100,alpha=1.0,bins='log',cmap='Greys',extent=extent,label=None,vmax=3000)

sf=t[t['CLASS_SFG']>0.95]
filt=t['RADIO_EXCESS']>0.95
filt&=~(t['r_50']>8) # remove large
re=t[filt]
rq=t[t['CLASS_RQAGN']>0.95]
hl=t[t['L_144']>1e26]
q=t[t['DR16Q']]

print('RE',len(re),'SF',len(sf),'HL',len(hl))

ax=plt.gca()
#plt.scatter(rq['w2w3'],rq['w1w2'],color='yellow',s=2,alpha=0.2)
d=DC()
d.density_contours(sf['w2w3'],sf['w1w2'],*extent,ccol(7),0.8,12,label='D24 SFG',legend=True,scatter=200)
d.density_contours(re['w2w3'],re['w1w2'],*extent,ccol(11),0.8,12,label='D24 RXG',legend=True,scatter=200)
d.density_contours(hl['w2w3'],hl['w1w2'],*extent,ccol(6),0.8,12,label='$L_{150}>10^{26}$ W Hz$^{-1}$',legend=True,scatter=200)
d.density_contours(q['w2w3'],q['w1w2'],*extent,ccol(3),0.8,12,label='DR16Q',legend=True,scatter=200)
d.do_scatter()

plt.xlabel('$W2-W3$ (Vega)')
plt.ylabel('$W1-W2$ (Vega)')
plt.xlim(extent[:2])
plt.ylim(extent[2:])

# shape is defined in AB mags
shape=np.array([(10,-0.8),(1.20,-0.8),(1.10,-0.55),(1.10,-0.25),(1.6,0),(2.2,0.3),(10,0.3)])
shapev=np.copy(shape)
shapev[:,0]+=corrs[2]-corrs[1]
shapev[:,1]+=corrs[1]-corrs[0]
plt.plot(shapev[:,0],shapev[:,1],color='black',label='SF region boundary')

ax=plt.gca()
handles, labels = ax.get_legend_handles_labels()
newArtist1 = plt.Line2D((0,1),(0,0), color='gray', marker='h', markersize=8, alpha=0.2, linestyle='')
newlabel=['All sources (%i)' % len(t)]
leg=plt.legend([newArtist1]+handles,newlabel+labels,loc=0,numpoints=1,prop={'size':20}, frameon=False,markerscale=1.5)
for i,lh in enumerate(leg.legendHandles): 
    if i>0:
        lh.set_alpha(1)

savefig('wisecc_class_density.pdf')

plt.show()

from plots import *
from scipy import stats
from density_contours2 import DC

t=Table.read('agn.fits')

for i in range(3):
    band='mag_w%i' % (i+1)
    err=band.replace('mag','magerr')
    filt=~np.isnan(t[err])
    print('band',i+1,'reducing to',np.sum(filt))
    t=t[filt]

print(len(t))
#t=t[t['abs_w3']<-24.5]

setupfig(size=(20,10))
corrs=[2.699,3.339,5.174,6.620]
table=t
for i in range(4):
    band='mag_w%i' % (i+1)
    vband=band+'v'
    table[vband]=table[band]-corrs[i]

t['w2w3']=t['mag_w2v']-t['mag_w3v']
t['w1w2']=t['mag_w1v']-t['mag_w2v']
extent=(0,5.0,-0.4,1.6)
plt.subplot(121)    
plt.hexbin(t['w2w3'],t['w1w2'],mincnt=5,gridsize=100,alpha=1.0,bins='log',cmap='Greys',extent=extent,label=None,vmax=1000)

herg=t[t['CLASS_HERG']>0.80]
lerg=t[(t['CLASS_LINELERG']>0.80)]
q=t[t['DR16Q']]

ax=plt.gca()
#plt.scatter(rq['w2w3'],rq['w1w2'],color='yellow',s=2,alpha=0.2)
d=DC()
d.density_contours(lerg['w2w3'],lerg['w1w2'],*extent,ccol(12),0.8,12,label='D24 LERG',legend=True,scatter=200)
d.density_contours(herg['w2w3'],herg['w1w2'],*extent,ccol(5),0.8,12,label='D24 HERG',legend=True,scatter=200)
d.density_contours(q['w2w3'],q['w1w2'],*extent,ccol(3),0.8,12,label='DR16Q',legend=True,scatter=200)
d.do_scatter()

plt.xlabel('$W2-W3$ (Vega)')
plt.ylabel('$W1-W2$ (Vega)')
plt.xlim(extent[:2])
plt.ylim(extent[2:])

ax=plt.gca()
handles, labels = ax.get_legend_handles_labels()
newArtist1 = plt.Line2D((0,1),(0,0), color='gray', marker='h', markersize=8, alpha=0.2, linestyle='')
newlabel=['All sources (%i)' % len(t)]
leg=plt.legend([newArtist1]+handles,newlabel+labels,loc=0,numpoints=1,prop={'size':20}, frameon=False,markerscale=1.5)
for i,lh in enumerate(leg.legendHandles): 
    if i>0:
        lh.set_alpha(1)

plt.subplot(122)
hq=t[np.logical_or(t['DR16Q'],t['CLASS_HERG']>0.80)]
#extent=(0,5)
extent=(-0.4,1.6)
k='w1w2'
_,bins,_=plt.hist(t[k],bins=60,range=extent,alpha=0.2,color='gray',label='All sources')
plt.hist(lerg[k],bins=bins,label='D24 LERG',alpha=0.5,color=ccol(12))
plt.hist(herg[k],bins=bins,label='D24 HERG',alpha=0.5,color=ccol(5))
plt.hist(q[k],bins=bins,label='DR16Q RL',alpha=0.5,color=ccol(3))
plt.hist(hq[k],bins=bins,label='RE',alpha=0.5,color='black',histtype='step')
plt.yscale('log')
plt.xlim(extent[0],extent[1])
plt.legend(loc=0)
plt.xlabel('$W2-W3$')
plt.ylabel('Number of sources')
        
savefig('wisecc_class_density_elc.pdf')

plt.show()

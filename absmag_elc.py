from plots import *

t=Table.read('agn.fits')
t=t[t['magerr_w3']>0]
herg=t[t['CLASS_HERG']>0.80]
lerg=t[(t['CLASS_LINELERG']>0.80) & (t['r_50']<8)]
q=t[t['DR16Q']]
hq=t[np.logical_or(t['DR16Q'],t['CLASS_HERG']>0.80)]
setupfig(size=(20,10))
extent=(-33,-19)
plt.subplot(121)
for band in [2,3]:
    k='abs_w'+str(band)
    _,bins,_=plt.hist(t[k],bins=60,range=extent,alpha=0.2,color='gray',label='All sources')
    plt.hist(lerg[k],bins=bins,label='D24 LERG',alpha=0.5,color=ccol(12))
    plt.hist(herg[k],bins=bins,label='D24 HERG',alpha=0.5,color=ccol(5))
    plt.hist(q[k],bins=bins,label='DR16Q RL',alpha=0.5,color=ccol(3))
    if k=='abs_w3':
        plt.hist(hq[k],bins=bins,label='RE',alpha=0.5,color='black',histtype='step')
    plt.yscale('log')
    plt.xlim(extent[1],extent[0])
    plt.legend(loc=0)
    plt.xlabel('Absolute $W%i$ magnitude' % band)
    plt.ylabel('Number of sources')
    plt.subplot(122)

savefig('absmag_elc_hist.pdf')
    
plt.show()



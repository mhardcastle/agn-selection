from plots import *

def bootstrap(data,function,iters=1000):
    if len(data)==0:
        return 0
    result=np.zeros(iters)
    for i in range(iters):
        x=np.random.choice(data,size=len(data),replace=True)
        result[i]=function(x)
    return result

setupfig()
nbins=30

t=Table.read('agn.fits')
filt=t['flag_mass']
filt&=t['z_best']<1.2
filt&=t['L_144']>1e25
filt&=t['Resolved']
t=t[filt]

logs=np.array(np.log10(t['Size']))
print(len(t))
lower=1.8
higher=3.8
bins=np.linspace(lower,higher,nbins+1)
indices=np.digitize(logs,bins=bins,right=True)
median=[]
yerr=[]
perc=[]
bc=(bins[:-1]+bins[1:])*0.5
for i in range(nbins):
    sample=np.array(t[indices==i]['Mass_median'])
    median.append(np.median(sample))
    yerr.append(np.percentile(bootstrap(sample,np.median)-median[-1],(16,84)))
    perc.append(np.percentile(sample,(5,95)))

yerr=np.abs(np.array(yerr).T)
perc=np.array(perc)
#plt.plot(bc,median)
plt.errorbar(bc,median,yerr=yerr)
print(median)
plt.fill_between(bc,perc[:,0],perc[:,1],alpha=0.2,label='All sources')

zbins=np.linspace(0,1.2,7)
for j in range(len(zbins)-1):
    filt=t['z_best']>zbins[j]
    filt&=t['z_best']<zbins[j+1]
    print(j,zbins[j],zbins[j+1],np.sum(filt))
    newls=logs[filt]
    newt=t[filt]
    median=[]
    indices=np.digitize(newls,bins=bins,right=True)
    for i in range(nbins):
        sample=np.array(newt[indices==i]['Mass_median'])
        median.append(np.median(sample))
    print(median)
    plt.plot(bc,median,label='$%.1f \le z < %.1f$ (%i)' % (zbins[j],zbins[j+1],len(newt)),ls=':')

ax=plt.gca()
set_logx(ax)
set_logy(ax)
plt.xlabel('Size (kpc)')
plt.ylabel('$M_*/M_\odot$')
plt.xlim(lower,higher)
plt.ylim(10,12)
savefig('agn_mass_size.pdf')
plt.show()


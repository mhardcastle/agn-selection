from plots import *
from astropy.stats import binom_conf_interval

def pltsc(numbers,area=5200,**kwargs):
    # numbers is a np.hist
    counts=numbers/area # to per sq deg
    counts/=width
    counts/=3.0462e-4 # to per sr
    #counts*=midpoint**2.5
    counts*=midpoint**2.5
    #plt.plot(centres,np.log10(counts),**kwargs)
    yerr=counts/np.sqrt(numbers)
    yerr=np.where(numbers<2,0,yerr)
    lower=np.log10(counts)-np.log10(counts-yerr)
    upper=np.log10(counts+yerr)-np.log10(counts)
    plt.errorbar(centres,np.log10(counts),yerr=[lower,upper],**kwargs)
    
    return counts


#plt.plot(centres,np.log10(hist),label='all')

setupfig(size=(20,10),ticks=False)
plt.subplot(121)
set_ticks(plt.gca())

t=Table.read('en1_classifications_dr1.fits')
#t=t[t['S_150MHz']>1.1e-3]
print(len(t))
logf=np.log10(t['S_150MHz'])

hist,bins=np.histogram(logf,bins=30)
limits=10**bins
centres=0.5*(bins[:-1]+bins[1:])
midpoint=0.5*(limits[:-1]+limits[1:])
width=limits[1:]-limits[:-1]

#plt.subplot(121)
pltsc(hist,label='EN1 all sources',color='black',ls=':',area=7.15)

cols={'RLAGN':ccol(11),'SFG':ccol(7),'RQAGN':ccol(3)}

for cut in ['RLAGN','SFG','RQAGN']:
    label=cut
    if cut!='RLAGN':
        ss=t[t['Overall_class']==cut]
    else:
        ss=t[(t['Overall_class']=='LERG') | (t['Overall_class']=='HERG')]
    print(label,len(ss))
    logf=np.log10(ss['S_150MHz'])
    nhist,_=np.histogram(logf,bins=bins)
    pltsc(nhist,area=7.15,label='EN1 '+label,ls=':',color=cols[label])
#plt.errorbar(np.nan,np.nan,yerr=np.nan,color='white',label=' ')
pltsc(nhist/10,area=7.15,label=' ',color='white')

    
#plt.gca().set_prop_cycle(None)
t=Table.read('fc.fits')
print(len(t))
logf=np.log10(t['Total_flux'])-3

hist,bins=np.histogram(logf,bins=30)
limits=10**bins
centres=0.5*(bins[:-1]+bins[1:])
midpoint=0.5*(limits[:-1]+limits[1:])
width=limits[1:]-limits[:-1]

pltsc(hist,label='DR2 all sources',color='black',area=5700)

t=Table.read('classified.fits')
logf=np.log10(t['Total_flux'])-3
hist,_=np.histogram(logf,bins=bins)
pltsc(hist,label='DR2 with IDs',color='grey')

cols={'RLAGN':ccol(11),'SF':ccol(7),'RQQ':ccol(3)}

for cut in ['AGN_BROAD','SF_EXCLUDE','RQQ_EXCLUDE']:
    label=cut.split('_')[0]
    if label=='AGN': label='RLAGN'
    ss=t[t[cut]]
    logf=np.log10(ss['Total_flux'])-3
    nhist,_=np.histogram(logf,bins=bins)
    pltsc(nhist,label='DR2 '+label,color=cols[label])

set_logx(plt.gca())
set_logy(plt.gca())
plt.xlim(-4,2)
plt.ylim(-3.2,4)
plt.xlabel('Flux density (Jy)')
plt.ylabel('$S^{5/2}n$ (sr$^{-1}$ Jy$^{1.5}$)')
plt.legend(loc=0,ncol=2)

plt.subplot(122)
set_ticks(plt.gca())
for cut in ['AGN_BROAD','SF_EXCLUDE','RQQ_EXCLUDE']:
    label=cut.split('_')[0]
    if label=='AGN': label='RLAGN'
    ss=t[t[cut]]
    logf=np.log10(ss['Total_flux'])-3
    nhist,_=np.histogram(logf,bins=bins)
    plt.plot(centres,np.log10(nhist/hist),label=label,color=cols[label])
    ci=binom_conf_interval(nhist,hist,interval='jeffreys')
    plt.fill_between(centres,np.log10(ci[0]),np.log10(ci[1]),alpha=0.3,color=cols[label])
    
    
set_logx(plt.gca())
set_logy(plt.gca())
plt.xlim(-3,2)
plt.ylim(-4,0.1)

plt.xlabel('Flux density (Jy)')
plt.ylabel('Fraction of classified sources')
plt.legend(loc=0)

savefig('numbercounts.pdf')
plt.show()

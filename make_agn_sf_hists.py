from plots import *

ta=Table.read('agn.fits')
ts=Table.read('sf.fits')
tq=Table.read('rqq.fits')

setupfig(size=(20,20))
plt.subplot(221)
_,bins,_=plt.hist(np.log10(ta['L_144']),bins=50,alpha=0.5,color=ccol(11),label='AGN')
plt.hist(np.log10(ts['L_144']),bins=bins,alpha=0.5,color=ccol(7),label='SF')
plt.hist(np.log10(tq['L_144']),bins=bins,alpha=0.5,color=ccol(3),label='RQQ')
ax=plt.gca()
xticks=ax.get_xticks()
ax.set_xticklabels(powerticks2(xticks))
plt.xlabel('$L_{144}$ (W Hz$^{-1}$)')
plt.legend(loc=0)

plt.subplot(222)
tar=ta[ta['Resolved']]
tsr=ts[ts['Resolved']]
tqr=tq[tq['Resolved']]
tar=tar[tar['Size']>0]
tsr=tsr[tsr['Size']>0]
tqr=tqr[tqr['Size']>0]

_,bins,_=plt.hist(np.log10(tar['Size']),bins=50,alpha=0.5,color=ccol(11))
plt.hist(np.log10(tsr['Size']),bins=bins,alpha=0.5,color=ccol(7))
plt.hist(np.log10(tqr['Size']),bins=bins,alpha=0.5,color=ccol(3))
ax=plt.gca()
xticks=ax.get_xticks()
ax.set_xticklabels(dologticks(xticks))
plt.xlabel('Physical size (kpc)')

plt.subplot(223)

_,bins,_=plt.hist(np.log10(1+ta['z_best']),range=(0,np.log10(4)),bins=50,alpha=0.5,color=ccol(11))
plt.hist(np.log10(1+ts['z_best']),bins=bins,alpha=0.5,color=ccol(7))
plt.hist(np.log10(1+tq['z_best']),bins=bins,alpha=0.5,color=ccol(3))
plt.xlim(0,np.log10(4))
ax=plt.gca()
xticks=ax.get_xticks()
ax.set_xticklabels(dologticks(xticks,subtract=1))
plt.xlabel('$z$')

plt.subplot(224)
_,bins,_=plt.hist(ta[ta['flag_mass']]['Mass_median'],range=(9,12),bins=50,alpha=0.5,color=ccol(11))
plt.hist(ts[ts['flag_mass']]['Mass_median'],bins=bins,alpha=0.5,color=ccol(7))
ax=plt.gca()
set_logx(ax)
plt.xlim(9,12)
plt.xlabel('$M_*$ (solar masses)')

savefig('agn_sf_hists.pdf')


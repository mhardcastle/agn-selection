from plots import *
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from solver import Evolve_RG
import glob
from synch_constants import *
from scipy import optimize
from scipy.interpolate import CubicSpline
import sys

try:
    png=sys.argv[1]=='png'
except:
    png=False

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
t3c=Table.read('3crr.txt',format='ascii')

t=Table.read('agn.fits')
#t=t[t['z_best']<1]
tr=t[t['Resolved']]

tu=t[~t['Resolved']]
print(len(tr),len(tu))
z=tu['z_best']
setupfig(fontsize=20 if png else 12)
#plt.scatter(tu['Length'],tu['L150'],alpha=0.1,label='Limits')
#plt.scatter(tr['Length'],tr['L150'],alpha=0.1,label='Resolved')
x_range=[0,4]
y_range=[21,30]
plt.hexbin(np.log10(tu['Size']),np.log10(tu['L_144']),bins='log',cmap=plt.get_cmap('Blues'),extent=(x_range[0],x_range[1],y_range[0],y_range[1]),gridsize=60,label=None)
plt.hexbin(np.log10(tr['Size']),np.log10(tr['L_144']),bins='log',cmap=plt.get_cmap('Greens'),extent=(x_range[0],x_range[1],y_range[0],y_range[1]),alpha=0.5,gridsize=60,label=None)
plt.scatter(np.log10(t3c['Size']),np.log10(t3c['L_178']*(178.0/150.0)**t3c['alpha']*np.pi*4),alpha=0.5,color='magenta',label='3CRR',zorder=3)

plt.xlabel('Total projected size (kpc)')
plt.ylabel('$L_{144}$ (W Hz$^{-1}$)')
ax = plt.gca()
set_logx(ax)
set_logy(ax)

handles, labels = ax.get_legend_handles_labels()
newArtist1 = plt.Line2D((0,1),(0,0), color='green', marker='h', markersize=8, alpha=0.5, linestyle='')
newArtist2 = plt.Line2D((0,1),(0,0), color='blue', marker='h', markersize=8, alpha=0.5, linestyle='')
newlabel=['Resolved','Size limits']
plt.legend([newArtist1,newArtist2]+handles,newlabel+labels,loc=4,numpoints=1,prop={'size':12}, frameon=False)

g=glob.glob('/home/mjh/git/analytic/dr2_agn/save*.0_13.4.pickle')
step=50*Myr
theta=np.arctan(5.7/2.5)
M=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
for f in g:
    env=Evolve_RG.load(f)
    print(env.Q)
    #plt.plot(np.log10(2*env.R/kpc),np.log10(env.synch*env.corrs[:,0]),zorder=2,ls=':',color='black')
    cs_r=CubicSpline(np.log10(env.tv),np.log10(2*env.R/kpc))
    cs_s=CubicSpline(np.log10(env.tv),np.log10(env.synch*env.corrs[:,0]))

    # Now we need to find the time value that causes the curve to cross the line

    txv=np.log10(2*env.R/kpc)-1.5
    tyv=np.log10(env.synch*env.corrs[:,0])-21
    tc=M.dot([txv,tyv])
    maxv=np.where(tc[1]<0)[0][0]-1
    print(maxv,maxv+1)
    print(env.tv[maxv],env.tv[maxv+1])

    def zero_find_function(logt):
        txv=cs_r(logt)-1.5
        tyv=cs_s(logt)-21
        tc=M.dot([txv,tyv])
        return tc[1]

    tzero=optimize.root(zero_find_function, [np.log10(env.tv[maxv]),np.log10(env.tv[maxv+1])]).x[0]
    print('cutoff is',10**tzero/Myr/1000,'Gyr')

    tv=np.linspace(np.log10(env.tv[1]),tzero,100)
    
    plt.plot(cs_r(tv),cs_s(tv),ls='--' if 'z1' not in f else ':',color='black',zorder=2)
    if 'z1' not in f:
        tv=np.arange(0,env.tv[-1]+step,step)
        tv=np.log10(tv[1:])
        tv=tv[tv<tzero]
        plt.scatter(cs_r(tv),cs_s(tv),zorder=2,marker='+',color='black')

# Line is 1.5 to 4, 21 to 26.7
    
plt.plot((1.5,x_range[1]),(y_range[0],26.7),ls='-')
    
plt.xlim(x_range)
plt.ylim(y_range)
    
#for lh in leg.legendHandles: 
#    lh.set_alpha(1)
#plt.show()
if png:
    savefig('pdd.png')
else:
    savefig('pdd.pdf')


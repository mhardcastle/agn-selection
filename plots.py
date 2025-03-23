from __future__ import division
from astropy.table import Table
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings('ignore')

# colour-blindness colours
def ccol(i):
    colours=[[0,0,0],[0,73,73],[0,146,146],[255,109,182],[255,182,219],[73,0,146],[0,109,219],[182,109,255],[109,182,255],[182,219,255],[146,0,0],[146,73,0],[219,209,0],[36,255,36],[255,255,109]]
    i-=1
    return [v/255.0 for v in colours[i]]

import os

def set_ticks(axis):
    axis.tick_params(axis='both',top=True,bottom=True,left=True,right=True,direction='in',which='both',pad=7)

def setupfig(size=(10,10),fontsize=20,ticks=True):
    rc('font',**{'family':'serif','serif':['Times'],'size':fontsize})
    rc('text', usetex=True)
    fig=plt.figure(figsize=size)
    plt.rcParams.update({'font.size':fontsize})
    if ticks:
        set_ticks(plt.gca()) # in current matplotlib does not work with subplot
        #formerly plt.gca().tick_params(axis='both',top=True,bottom=True,left=True,right=True,direction='in',which='both',pad=7)
    return fig

def savefig(filename):
    plt.tight_layout()
    plt.savefig(filename)
    if '.pdf' in filename:
        os.system('pdfcrop '+filename)

def powerticks(ticks):
    rticks=[]
    for t in ticks:
        i=int(t)
        if float(i)!=t:
            rticks.append('')
        elif i==0:
            rticks.append('1')
        elif i==-1:
            rticks.append('0.1')
        elif i==1:
            rticks.append('10')
        else:
            rticks.append('$10^{%i}$' % i)
    return rticks

def powerticks2(ticks):
    return ['$10^{%i}$' % t for t in ticks]

def dologticks(ticks,subtract=0):
    labels=[]
    for tk in ticks:
        dp=tk
        if dp>=1:
            dp=0
        else:
            dp=1+int(0.5-dp)
        labels.append(('%.'+str(dp)+'f') % ((10**tk)-subtract))
    return labels

def set_logx(ax):
    limits=ax.get_xlim()
    minpower=np.floor(limits[0])
    maxpower=np.ceil(limits[1])
    ax.set_xticks(list(range(int(minpower),int(maxpower+1))))
    labels=[]
    for v in range(int(minpower),int(maxpower+1)):
        if v>-3 and v<4:
            labels.append(str(10**v))
        else:
            labels.append("$10^{%i}$" % v)
    ax.set_xticklabels(labels)
    minorloc=[]
    for v in range(int(minpower),int(maxpower)):
        minorloc+=list(np.arange(10.0**v,10.0**(v+1),10.0**v)[1:])
    minor=ticker.FixedLocator(np.log10(minorloc))
    ax.xaxis.set_minor_locator(minor)

def set_logy(ax):
    limits=ax.get_ylim()
    minpower=np.floor(limits[0])
    maxpower=np.ceil(limits[1])
    ax.set_yticks(list(range(int(minpower),int(maxpower+1))))
    labels=[]
    for v in range(int(minpower),int(maxpower+1)):
        if v>-3 and v<4:
            labels.append(str(10**v))
        else:
            labels.append("$10^{%i}$" % v)
    ax.set_yticklabels(labels)
    minorloc=[]
    for v in range(int(minpower),int(maxpower)):
        minorloc+=list(np.arange(10.0**v,10.0**(v+1),10.0**v)[1:])
    minor=ticker.FixedLocator(np.log10(minorloc))
    ax.yaxis.set_minor_locator(minor)

if __name__=='__main__':
    setupfig()
    for i in range(1,16):
        plt.scatter(i,16,color=ccol(i),label=str(i))
    plt.legend(loc=0)
    plt.show()
    

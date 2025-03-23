import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
from plots import *

def cdf(data,bootstrap=0,alphabase=2,**kwargs):
    sdata=sorted(data)
    p=plt.step(sdata, np.linspace(0,1,len(sdata)),**kwargs)
    color=p[0].get_color()
    if 'label' in kwargs: del kwargs['label']
    if 'color' in kwargs: del kwargs['color']
    bsbins=np.linspace(min(sdata),max(sdata),80)
    bsr=[]
    for i in range(bootstrap):
        bsv=[]
        newdata=sorted(np.random.choice(data,len(data),replace=True))
        for j in bsbins:
            bsv.append(np.sum(newdata<j))
        bsr.append(bsv)
    bsr=np.array(bsr)/len(data)
    lower,upper=np.percentile(bsr,(16,84),axis=0)
    plt.fill_between(bsbins,lower,upper,alpha=alphabase/10,color=color)
        #plt.step(newdata, np.linspace(0,1,len(sdata)),color=color,alpha=alphabase/bootstrap,**kwargs)

if __name__=='__main__':
        
    t=Table.read('large_cores_ccc.fits')
    setupfig()
    colors={'DR16Q':ccol(3),'N':ccol(5),'RI':ccol(12),'Q':ccol(7),'RE':ccol(8)}

    t=t[t['Core_3000']>0]
    print(len(t),'with cores')
    t['prominence']=3+np.log10(t['Core_3000']/t['Total_flux'])

    cdf(t['prominence'],label='All sources (%i)' % len(t),bootstrap=100,color='black',lw=2)
    for label in ['RE','Q','N','RI','DR16Q']:
        st=t[t[label]]
        cdf(st['prominence'],label=label+' (%i)' % len(st),bootstrap=100,color=colors[label])

    plt.xlim(-3.5,0)
    plt.ylim(0,1)
    plt.xlabel('log$_{10}$(Core prominence at 3 GHz)')
    plt.legend(loc=0)
    plt.grid()
    savefig('cores_cdist.pdf')
    #plt.show()


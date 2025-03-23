import glob
import matplotlib.pyplot as plt
import numpy as np
import emcee
import corner

def model(x,logC,logPstar,alpha,beta):
# power-law model from Mauch and Sadler
    PstaroverP=10**(logPstar-x) # work in linear space
    return 10**logC/(PstaroverP**-alpha+PstaroverP**-beta)

def lnlike(X,x,y,yerr):
    logC,logPstar,alpha,beta=X
    return -0.5*np.sum((y - model(x,logC,logPstar,alpha,beta))**2 / yerr**2)

def lnprior(X):
    logC,logPstar,alpha,beta=X
    if -9 < logC < 0 and 20 < logPstar < 29 and 0<alpha<1 and 1<beta<2.5:
        return 0.0
    return -np.inf

def lnpost(X,x,y,yerr):
    return lnprior(X)+lnlike(X,x,y,yerr)

if __name__=='__main__':
    g=sorted(glob.glob('LF*W1-rebin-fr_corr.csv'))
    labels=['logC','logPstar','alpha','beta','phi25']
    plots=False
    outd=[[],[],[],[],[]]
    for f in g:
        print(f)
        zrange=f.replace('.csv','').replace('LF_AGN-','').split('-')
        zmid=(float(zrange[0])+float(zrange[1]))/2
        d=np.loadtxt(f,delimiter=',')
        d[:,1]=np.where(d[:,1]>0.001,np.nan,d[:,1])
        d[:,1]=np.where(d[:,1]<1.5*d[:,2],np.nan,d[:,1])
        d=d[~np.isnan(d[:,1])]
        d=d[d[:,0]>22.5]

        nwalkers=100
        ndim=4
        pos=[[-5,24,0.5,1.5]+0.05*np.random.normal(size=ndim)
             for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost,
                                    args=(d[:,0],d[:,1],d[:,2]))

        sampler.run_mcmc(pos, 2000)
        if plots:
            for i in range(ndim):
                plt.subplot(ndim,1,i+1)
                plt.plot(sampler.chain[:,:,i].transpose())
                plt.ylabel(labels[i])
            plt.xlabel('Samples')
            plt.show()
        samples=sampler.chain[:, 500:, :].reshape((-1, ndim))
        phi25=np.log10(model(25,samples[:,0],samples[:,1],samples[:,2],samples[:,3]))
        samples=(np.vstack([samples.T,phi25])).T
        with open(f.replace('.csv','-samples.npy'),'wb') as outfile:
            np.save(outfile,samples)
        meanparms=np.median(samples,axis=0)
        for i in range(5):
            errors=np.percentile(samples[:,i],(16,84))
            merr=((meanparms[i]-errors[0])+(errors[1]-meanparms[i]))/2
            outd[i].append([zmid,meanparms[i],merr])
            print(labels[i],meanparms[i],np.percentile(samples[:,i],(16,84)))

        # make triangle plot
        fig = corner.corner(samples,labels=labels)
        plt.savefig(f.replace('.csv','-corner.pdf'))
        plt.close()

        #plot data and model
        if plots:
            plt.errorbar(d[:,0],np.log10(d[:,1]),yerr=np.log10(d[:,1]+d[:,2])-np.log10(d[:,1]),label=f)
            plt.plot(d[:,0],np.log10(model(d[:,0],meanparms[0],meanparms[1],meanparms[2],meanparms[3])))
            #plt.plot(d[:,0],np.log10(model(d[:,0],-5.5,24.59,-0.49,-1.27)))
            plt.show()

    d=np.array(outd)
    for i in range(5):
        plt.subplot(2,3,i+1)
        sd=d[i]
        plt.errorbar(sd[:,0],sd[:,1],yerr=sd[:,2],marker='o',ls='none')
        plt.xlabel('z')
        plt.ylabel(labels[i])

    plt.show()

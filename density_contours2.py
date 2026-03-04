from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
from matplotlib import rc
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

# class version where the dots are stored and plotted at the end

class DC(object):
    def __init__(self):
        self.scatter=[]

    def density_contours(self,x,y,xmin,xmax,ymin,ymax,color,alpha,n_bin=10,label=None,scatter=None,sc_size=8,legend=False,level='sqrt'):
        '''
        make 2D density plots from data points
        originally by Miranda Jarvis with some modifications by Martin

        parameters
        -------------------
        x : (N,) array_like, the x axis values of the data points you want a density plot from 
        y : (N,) array_like, the y axis values of the data points you want a density plot from
        xmin : number, the minimum x axis value to calculate the density plot on 
        xmax : number, the maximum x axis value to calculate the density plot on 
        ymin : number, the minimum y axis value to calculate the density plot on 
        ymax : number, the maximum y axis value to calculate the density plot on 
        color : the colour to use for the density plot contours (fades from white to color)
        alpha : the alpha value (opaqueness 1 = opaque 0 = transparent) to use for the density plot
        n_bin : how many steps to use in creating white to color colour scale (impacts the final color and speed of colour chance of density contours)
        label : If not None, overplot this text at the peak of the density contours
        scatter : If not None or 0, plot this many example points from the distribution over the contours.
        sc_size : The size of points to plot
        legend : If true populate the legend rather than using the label
        level : a string describing the level creation algorithm ('linear', 'sqrt' and 'logarithmic' are supported) 

        returns: 
        None ; adds density contours to axis (ax) 
        '''
        #create grid to calculate density contours on 
        X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([X.ravel(), Y.ravel()])

        #reformat values 
        values = np.vstack([x,y])

        #calculate density 
        kernel = stats.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)

        #create colour scale 

        #cm = LinearSegmentedColormap.from_list('test', colors, N=n_bin)
        alpha=np.linspace(0,alpha,n_bin)
        clist=[]
        for a in alpha:
            print(a)
            clist.append(list(color[0:3])+[a,])
        print(clist)
        cm=ListedColormap(clist)

        #define what levels want to put contours at (these were arbitrarily chosen to look good, very subjective!!!
        if level=='linear':
            levels=np.linspace(np.amax(Z)/n_bin,np.amax(Z),n_bin)
        elif level=='sqrt':
            levels=np.linspace(np.sqrt(np.max(Z))/n_bin,np.sqrt(np.max(Z)),n_bin)**2
        elif level=='log':
            print(np.log10(np.max(Z)/n_bin),np.log10(np.max(Z)),n_bin)
            levels=np.logspace(np.log10(np.max(Z)/n_bin),np.log10(np.max(Z)),n_bin)
        else:
            raise NotImplementedError('Not implemented: '+level)
        print(levels)
        #do contours
        xp=np.rot90(X)
        yp=np.rot90(Y)
        zp=np.rot90(Z)
        plt.contourf(xp,yp,zp,levels,cmap=cm, origin='upper',label=(label+' (%i)' % len(x)) if legend else None)
        if label is not None and not legend:
            print(zp.shape)
            pmax=np.unravel_index(zp.argmax(), zp.shape)
            print(pmax)
            xmax=xp[pmax]
            ymax=yp[pmax]
            plt.text(xmax,ymax,label,horizontalalignment='center',verticalalignment='center')
        if scatter:
            if len(x)<scatter:
                pvx=x
                pvy=y
            else:
                pv=np.random.choice(len(x),scatter,replace=False)
                pvx=x[pv]
                pvy=y[pv]
        self.scatter.append((pvx,pvy,color,sc_size,(label+' (%i)' % len(x))))
        #plt.scatter(pvx,pvy,color=color,s=sc_size,label=(label+' (%i)' % len(x)))

    def do_scatter(self):
        for s in self.scatter:
            plt.scatter(s[0],s[1],color=s[2],s=s[3],label=s[4])
            

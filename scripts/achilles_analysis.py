import numpy as np
import h5py
import os.path
from math import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import ticker
au=1.496e13
parsec=3.086e+18
lsun=3.99e33
msun=1.99e33
rsun=6.96e10
rjupiter=7.140e9
yr=3.154e7
amu=1.660539040e-24
arad=7.5646e-15
kb=1.380658e-16
c_light=2.9979e10
def plot2d(shape,data,filenumber,x_range,y_range,zmin,zmax,time,fname='file',title='color',scale='linear',figuresize=[10,8]):
    fs=20   #fontsize
    nx=shape[1]
    ny=shape[0]
    x=np.linspace(x_range[0],x_range[1],nx+1)
    y=np.linspace(y_range[0],y_range[1],ny+1)
    Z=np.asarray(data)
    X,Y=np.meshgrid(x,y)
    plt.figure(figsize=(figuresize[0],figuresize[1]))
    if (scale=='log'):
        p=plt.pcolor(X,Y,Z,norm=colors.LogNorm(vmin=zmin,vmax=zmax),cmap='viridis')
        cb=plt.colorbar(p)
        cb.ax.tick_params(labelsize=fs)
        cb.set_ticks(np.logspace(log10(zmin),log10(zmax),num=5))
    else:
        p=plt.pcolor(X,Y,Z,vmin=zmin,vmax=zmax,cmap='viridis',)
        cb=plt.colorbar(p)
        cb.ax.tick_params(labelsize=fs)
        cb.set_ticks(np.linspace(zmin,zmax,num=5))
    plt.yticks([])
    plt.xticks([])
    plt.axes().set_aspect(ny/nx)
    plt.title(title+' t='+str(time)+' s',fontsize=fs)
    plt.tight_layout()
    #plt.show()
    plt.savefig(fname+filenumber+'.png',dpi=200,bbox_inches='tight')
    plt.close()

def plot1d(shape,data,filenumber,x_range,ymin,ymax,fname='file',scale='linear',xlabel='',ylabel='',title='',figuresize=[8,8],
        figmarker='',xscale=1,yscale=1,data2=None):
    fs=22   #fontsize
    nx=shape[0]
    x=np.linspace(x_range[0]/xscale,x_range[1]/xscale,nx)
    y=np.asarray(data)/yscale
    plt.figure(figsize=(figuresize[0],figuresize[1]))
    plt.plot(x,y,color='k',marker=figmarker)
    if data2 is None:
        data2=None
    else:
        y2=np.asarray(data2)/yscale
        plt.plot(x,y2,color='r',marker=figmarker)
    plt.ylim((ymin/yscale,ymax/yscale))
    plt.xlim((x_range[0]/xscale,x_range[1]/xscale))
    plt.xlabel(xlabel,fontsize=fs)
    plt.ylabel(ylabel,fontsize=fs)
    plt.title(title,fontsize=fs)
    if (scale=='log'):
        plt.yscale('log')
    if (fname=='tau'):
        plt.axhline(y=0.6667,color='r',linestyle='-')
    plt.tight_layout()
    plt.savefig(fname+filenumber+'.png',dpi=200)
    #plt.show()
    plt.close()

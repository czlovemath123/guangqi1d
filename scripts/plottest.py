import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
from numpy import *
import math as m
import os.path
matplotlib.rc('lines',linewidth=2)
matplotlib.rcParams.update({'font.size': 24, 'font.family': 'serif'})
for i in range(1,2):
    data=[]
    filename1='extrm'+str(i)+'.dat'
    filename2='extrm'+str(i)+'.png'
    print filename1,filename2
    outdir='/Users/zchen/bitbucket/rmhd/modules/problem/out'
#   need the absolute directory
    f=open(os.path.join(outdir,filename1),'r')
    for j in range(100):
        line=f.readline()
        a=[float(s) for s in line.split()]
        data.append(a)
    data=transpose(asarray(data))
    rhomin=min(data[0])*0.9
    rhomax=max(data[0])*1.1
    umin=min(data[1])-0.1*abs(min(data[1]))
    umax=max(data[1])+0.1*abs(max(data[1]))
    pmin=min(data[4])*0.9
    pmax=max(data[4])*1.1
    einternal=zeros(100)
    for i in range(100):
        einternal[i]=data[4,i]/data[0,i]/(data[5,i]-1)
    einternal=asarray(einternal)
    emin=min(einternal)*0.9
    emax=max(einternal)*1.1
    x=0.01*arange(100)
    plt.figure(figsize=(15,8))
    ax=plt.subplot(2,2,1)
    plt.plot(x,data[0])
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\rho$')
    plt.xlim(0,0.99)
    plt.ylim(rhomin,rhomax)
    ax=plt.subplot(2,2,2)
    plt.plot(x,data[1])
    plt.xlabel(r'$x$')
    plt.ylabel(r'$velocity$')
    plt.xlim(0,0.99)
    plt.ylim(umin,umax)
    plt.subplot(2,2,3)
    ax=plt.plot(x,data[4])
    plt.xlabel(r'$x$')
    plt.ylabel(r'$pressure$')
    plt.xlim(0,0.99)
    plt.ylim(pmin,pmax)
    plt.subplot(2,2,4)
    ax=plt.plot(x,einternal)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$e_{internal}$')
    plt.xlim(0,0.99)
    plt.ylim(emin,emax)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,filename2))
    plt.close()
    f.close()

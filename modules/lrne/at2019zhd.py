import matplotlib.pyplot as plt
import sys
import matplotlib
import numpy as np
from math import *
matplotlib.rcParams.update({'font.size': 20, 'font.family': 'serif'})

def read_file(filename,n):
    f=open(filename,'r')
    line=f.readline()
    data=[]
    i=0
    while line:
        d=[float(s) for s in line.split()]
        data.append(d)
        line=f.readline()
        i=i+1
        if (i>n):
            data.pop(0)
    f.close()
    return data

def results():
    fig,axes=plt.subplots(figsize=(8,6))
    axes.plot(times[0],luminosities[0],color='r',linestyle='-',markersize=1,label='shock')
    axes.plot(times[1],luminosities[1],color='b',linestyle='-',markersize=1,label='shock-free')
    axes.plot(t0,l0,'k.-',label='AT2019zhd')
    #axes.plot(t1,l1,'r.-',markersize=1,label='no_arad')
    #axes.plot(t2,l2,'b.-',markersize=1,label='arad')
    axes.set_yscale('log')
    axes.set_xlim(t[0],40)
    axes.set_ylim(1e38,4e39)
    axes.legend()
    axes.set_xlabel('t [day]')
    axes.set_ylabel('L [erg/s]')
    plt.tight_layout()
    plt.savefig('at2019zhd.png',bbox_inches='tight',dpi=300)
    plt.show()


narg=len(sys.argv)
simfiles=[]
if (narg>=2):
    for i in range(1,narg):
        simfile='model'+str(sys.argv[i]).zfill(2)+'/history.data'
        simfiles.append(simfile)
else:
    simfile='history.data'
    simfiles.append(simfile)
filename='at2019zhd.dat'
d=np.transpose(np.asarray(read_file(filename,100)))
t0=d[0]
l0=d[3]*1e37
times=[]
luminosities=[]
for k in range(narg-1):
    d=np.transpose(np.asarray(read_file(simfiles[k],100000)))
    t=(d[0]-4.176e14/3e10)/86400
    l=d[1]
    for i in range(len(t)):
        if (t[i]>0):
            n=i+10
            break
    for i in range(n):
        t=np.delete(t,0)
        l=np.delete(l,0)
    tmax=np.where(l==np.amax(l))
    t=t-t[tmax]
    times.append(t)
    luminosities.append(l)
#simfile='best/history.data'
#d1=np.transpose(np.asarray(read_file(simfile,100000)))
#simfile='best_arad/history.data'
#d2=np.transpose(np.asarray(read_file(simfile,100000)))
#t1=(d1[0]-4.176e14/3e10)/86400
#t2=(d2[0]-4.176e14/3e10)/86400
#for i in range(len(t1)):
#    if (t1[i]>0):
#        n=i+10
#        break
#l1=d1[1]
#for i in range(n):
#    t1=np.delete(t1,0)
#    l1=np.delete(l1,0)
#tmax=np.where(l1==np.amax(l1))
#t1=t1-t1[tmax]
#for i in range(len(t2)):
#    if (t2[i]>0):
#        n=i+10
#        break
#l2=d2[1]
#for i in range(n):
#    t2=np.delete(t2,0)
#    l2=np.delete(l2,0)
#tmax=np.where(l2==np.amax(l2))
#t2=t2-t2[tmax]

colors=['r','b','g','c']
labels=['fit',r'$\gamma=1.1$',r'$\gamma=1.4$',r'$\gamma=1.667$']
labels=['fit',r'$a_{\rm{rad}}=0$']
labels=[r'$E_{r}/e_{g}=21.5$',r'$E_{r}/e_{g}=2.15$',r'$E_{r}/e_{g}=0.215$']

fig,axes=plt.subplots(figsize=(8,6))
axes.plot(t0,l0,'k.-',label='AT2019zhd')
for i in range(len(times)):
    axes.plot(times[i],luminosities[i],color=colors[i],linestyle='-.',markersize=1,label=sys.argv[i+1])
    #axes.plot(times[i],luminosities[i],color=colors[i],linestyle='-.',markersize=1,label=labels[i])
#axes.plot(t1,l1,'r.-',markersize=1,label='no_arad')
#axes.plot(t2,l2,'b.-',markersize=1,label='arad')
axes.set_yscale('log')
axes.set_xlim(t[0],40)
axes.set_ylim(1e38,4e39)
axes.legend()
axes.set_xlabel('t [day]')
axes.set_ylabel('L [erg/s]')
plt.tight_layout()
plt.savefig('at2019zhd.png',dpi=200)
plt.show()

#results()

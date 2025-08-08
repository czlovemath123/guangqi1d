import sys
sys.path.insert(1,'../../scripts')
from assemble_2d_data import *
import numpy as np
from math import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 20, 'font.family': 'serif'})


def read_file(filename):
    f=open(filename,'r')
    line=f.readline()
    data=[]
    i=0
    while line:
        d=[float(s) for s in line.split()]
        data.append(d)
        line=f.readline()
        i=i+1
    f.close()
    return data

filename='history.data'
d=read_file(filename)
d=np.transpose(np.asarray(d))
t=d[0]
t=(t-t[0])/day
l=d[1]*2

fig,axes=plt.subplots(1,1,figsize=(8,6),sharex=True,squeeze=True)
axes.plot(t,l,'k-')
axes.set_xlabel('days')
#axes.set_xlim(t[0],t[-1])
axes.set_xlim(t[0],4e6/86400)
axes.set_ylabel(r'$L_{\odot}$')
axes.set_yscale('log')
axes.set_ylim(min(l)/2,max(l)*2)
plt.tight_layout()
plt.savefig(filename+'.png')
plt.show()

from phy_const import *
import sys
import numpy as np
import h5py
import f90nml
from math import *
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('agg')
#au=1.496e13
#parsec=3.086e+18
#G=6.67259e-8
#lsun=3.99e33
#msun=1.99e33
#rsun=6.96e10
#rjupiter=7.140e9
#mjupiter=1.899e30
#yr=3.154e7
#amu=1.660539040e-24
#arad=7.5646e-15
#kb=1.380658e-16
#h_planck=6.6260755e-27
#c_light=2.9979e10
#mearth=5.976e27
#day=86400
#sigma_sb=5.67051e-5
#NA=6.0221367e23
def read_attr(filename):
    f=h5py.File(filename,'r')
    time=f.attrs['Time'][0]
    nblocks=f.attrs['nblocks'][0]
    blk_size_nx=f.attrs['blk_size_nx'][0]
    return f,time,nblocks,blk_size_nx

def glb_level(filename):
    f,time,nblock,nx=read_attr(filename)
    level=[]
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        level_blk=np.ones(nx)*f[blkname].attrs['level'][0]
        if (i==0):
            level.append(level_blk)
            level=np.asarray(level)
        else:
            level=np.vstack((level,level_blk))
    level=np.ravel(level)
    return level

def glb_x_center(filename):
    f,time,nblock,nx=read_attr(filename)
    var=[]
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        var_blk=f[blkname]['x_center']
        if (i==0):
            var.append(var_blk)
            var=np.asarray(var)
        else:
            var=np.vstack((var,var_blk))
    var=np.ravel(var)
    return var

def glb_x_inter(filename):
    f,time,nblock,nx=read_attr(filename)
    var=[]
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        var_blk=f[blkname]['mesh_x']
        if (i==0):
            var.append(var_blk)
            var=np.asarray(var)
        else:
            var=np.vstack((var,var_blk))
    v=var[0,0]
    var=np.delete(var,0,1)
    var=np.ravel(var)
    var=np.insert(var,0,v)
    return var

def glb_cell_var(filename,vname):
    f,time,nblock,nx=read_attr(filename)
    var=[]
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        var_blk=f[blkname][vname][0]
        if (i==0):
            var.append(var_blk)
            var=np.asarray(var)
        else:
            var=np.vstack((var,var_blk))
    var=np.ravel(var)
    return var

def glb_inter_var(filename,vname):
    f,time,nblock,nx=read_attr(filename)
    var=[]
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        var_blk=f[blkname][vname][0]
        #print(var_blk)
        if (i==0):
            var.append(var_blk)
            var=np.asarray(var)
        else:
            var=np.vstack((var,var_blk))
    v=var[0,0]
    var=np.delete(var,0,1)
    var=np.ravel(var)
    var=np.insert(var,0,v)
    #exit()
    return var

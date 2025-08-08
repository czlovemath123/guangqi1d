import sys
import numpy as np
import h5py
import f90nml
from math import *
from phy_const import *
from scipy import interpolate
from scipy.interpolate import griddata
from scipy.interpolate import interpn
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
#mh=1.67e-24

cbname_dict={'rho':r'$\rho\ [\rm{g}\cdot\rm{cm}^{-3}]$','vx':r'$v_{x}\ [cm/s]$','vy':r'$v_{y}\ [cm/s]$',        \
    'viscous_heat':r'$\Lambda\ [\rm{erg}\cdot\rm{cm}^{-3}]$','temp':r'$T_{\rm{g}}$'+' [K]','pres':'p [dyn]',    \
    'Rosseland':r'$\kappa_{\rm{R}}\ [\rm{cm}^2\cdot\rm{g}^{-1}]$', 'Erad':r'$\rm{erg}\cdot\rm{cm}^3$',          \
    'entropy':'kb per particle','irrad':r'$erg/s/cm^{2}$'}

def read_attr(filename):
    f=h5py.File('out/'+filename+'.h5','r')
    time=f.attrs['Time'][0]
    nblocks=f.attrs['nblocks'][0]
    blk_size_nx=f.attrs['blk_size_nx'][0]
    blk_size_ny=f.attrs['blk_size_ny'][0]
    return f,time,nblocks,blk_size_nx,blk_size_ny

def read_levels(filename):
    f=h5py.File('out/'+filename+'.h5','r')
    levels=f['levels']
    return np.asarray(levels)

def blk_mesh(f,blkname):
    mesh_x=np.asarray(f[blkname]['mesh_x'])
    mesh_y=np.asarray(f[blkname]['mesh_y'])
    return mesh_x,mesh_y

def blk_center_polar(f,blkname):
    x=np.asarray(f[blkname]['x_center'])
    y=np.asarray(f[blkname]['y_center'])
    return x,y

def polar_center(f,blkname):
    blk_size_nx=f.attrs['blk_size_nx'][0]
    blk_size_ny=f.attrs['blk_size_ny'][0]
    r=np.asarray(f[blkname]['x_center'])
    theta=np.asarray(f[blkname]['y_center'])
    x=np.ndarray(shape=(blk_size_ny,blk_size_nx))
    y=np.ndarray(shape=(blk_size_ny,blk_size_nx))
    for i in range(blk_size_ny):
        for j in range(blk_size_nx):
            x[i][j]=r[j]*cos(theta[i])
            y[i][j]=r[j]*sin(theta[i])
    return x,y

def cartesian_mesh(f,blkname):
    blk_size_nx=f.attrs['blk_size_nx'][0]
    blk_size_ny=f.attrs['blk_size_ny'][0]
    x=f[blkname]['mesh_x']
    y=f[blkname]['mesh_y']
    mx,my=np.meshgrid(x,y,indexing='ij')

def polar_mesh(f,blkname):
    blk_size_nx=f.attrs['blk_size_nx'][0]
    blk_size_ny=f.attrs['blk_size_ny'][0]
    r=f[blkname]['mesh_x']
    theta=f[blkname]['mesh_y']
    x=np.ndarray(shape=(blk_size_ny+1,blk_size_nx+1))
    y=np.ndarray(shape=(blk_size_ny+1,blk_size_nx+1))
    for i in range(blk_size_ny+1):
        for j in range(blk_size_nx+1):
            x[i][j]=r[j]*cos(theta[i])
            y[i][j]=r[j]*sin(theta[i])
    return x,y

def cartesian_v(f,blkname):
    blk_size_nx=f.attrs['blk_size_nx'][0]
    blk_size_ny=f.attrs['blk_size_ny'][0]
    vr=f[blkname]['vx']
    vtheta=f[blkname]['vy']
    vx=np.ndarray(shape=(blk_size_ny,blk_size_nx))
    vy=np.ndarray(shape=(blk_size_ny,blk_size_nx))
    #r=np.asarray(f[blkname]['x_center'])
    theta=np.asarray(f[blkname]['y_center'])
    for i in range(blk_size_ny):
        for j in range(blk_size_nx):
            vx[i][j]=sin(theta[i])*vr[i,j]+cos(theta[i])*vtheta[i,j]
            vy[i][j]=cos(theta[i])*vr[i,j]-sin(theta[i])*vtheta[i,j]
    return vx,vy

def blk_var(f,blkname,vname):
    var=f[blkname][vname]
    return var

def cartesian_to_polar(x,y):
    r=sqrt(x**2+y**2)
    if y==0:
        if (x>0):
            theta=pi/2
        else:
            theta=-pi/2
    else:
        theta=atan(x/y)
    return r,theta

def polar_to_cartesian(r,theta):
    x=np.ndarray(shape=(len(theta),len(r)))
    y=np.ndarray(shape=(len(theta),len(r)))
    for i in range(len(theta)):
        for j in range(len(r)):
            x[i][j]=r[j]*cos(theta[i])
            y[i][j]=r[j]*sin(theta[i])
    x=x.reshape(len(r)*len(theta))
    y=y.reshape(len(r)*len(theta))
    return x,y

def glb_level_max(filename):
    f,time,nblock,nx,ny=read_attr(filename)
    level_max=0
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        level_max=max(level_max,f[blkname].attrs['level'][0])
    return level_max

def glb_cell_fix_r(filename,vname,r):
    f,time,nblock,nx,ny=read_attr(filename)
    theta=[]
    profile=[]
    cell=[]
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        x_center,y_center=blk_center_polar(f,blkname)
        x_mesh,y_mesh=blk_mesh(f,blkname)
        if (r>=x_mesh[0] and r<x_mesh[-1]):
            var=blk_var(f,blkname,vname)
            var_interpolator=interpolate.interp2d(x_center,y_center,var,kind='linear')
            var_interp=np.zeros(ny)
            for j in range(ny):
                var_interp[j]=var_interpolator(r,y_center[j])
            theta.append(y_center)
            profile.append(var_interp)
            cell.append(y_mesh)
    theta=np.asarray(theta)
    profile=np.asarray(profile)
    cell=np.asarray(cell)
    theta,profile=array_sort(theta,profile)
    cell=np.unique(cell)
    theta=np.ravel(theta)
    profile=np.ravel(profile)
    cell=np.ravel(cell)
    return theta,profile,cell

def glb_cell_fix_R(filename,vname,R):
    f,time,nblock,nx,ny=read_attr(filename)
    z=[]
    profile=[]
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        x_center,y_center=blk_center_polar(f,blkname)
        x_mesh,y_mesh=blk_mesh(f,blkname)
        r_min=x_mesh[0]
        r_max=x_mesh[nx]
        theta_min=y_mesh[0]
        theta_max=y_mesh[ny]
        R1=r_min*sin(theta_min)
        R2=r_min*sin(theta_max)
        R3=r_max*sin(theta_min)
        R4=r_max*sin(theta_max)
        if (R>=R1 and R<R4):    #has intersection
            var=blk_var(f,blkname,vname)
            var_interpolator=interpolate.interp2d(x_center,y_center,var,kind='linear')
            if (R<R2):
                z_min=sqrt(r_min**2-R**2)
            else:
                if (abs(theta_max-pi/2)<1e-4):
                    z_min=0.0
                else:
                    z_min=R/tan(theta_max)
            if (R<R3):
                z_max=R/tan(theta_min)
            else:
                z_max=sqrt(r_max**2-R**2)
            n=max(nx,ny)+1
            z_blk=np.linspace(z_min,z_max,n)
            var_interp=np.zeros(n)
            line_polar=np.zeros((2,n))
            for j in range(n):
                line_polar[0][j],line_polar[1][j]=cartesian_to_polar(R,z_blk[j])
            for j in range(n):
                var_interp[j]=var_interpolator(line_polar[0][j],line_polar[1][j])
            z.append(z_blk)
            profile.append(var_interp)
    z=np.asarray(z)
    profile=np.asarray(profile)
    z,profile=array_sort(z,profile)
    z=np.ravel(z)
    profile=np.ravel(profile)
    return z,profile

def glb_sum(filename,vname):
    f,time,nblock,nx,ny=read_attr(filename)
    sum=0
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        x_mesh,y_mesh=blk_mesh(f,blkname)
        var=f[blkname][vname]
        var=np.transpose(var)
        for j in range(nx):
            for k in range(ny):
                r1=x_mesh[j]
                r2=x_mesh[j+1]
                theta1=y_mesh[k]
                theta2=y_mesh[k+1]
                vol=1.0/3.0*(r2**3-r1**3)*(cos(theta1)-cos(theta2))*2*pi
                sum=sum+var[j,k]*vol
    return sum

def array_sort(a,b):
    n=len(a)
    for j in range(n-1,0,-1):
        for i in range(j):
            if (a[i][0]>a[i+1][0]):
                a[[i,i+1],:]=a[[i+1,i],:]
                b[[i,i+1],:]=b[[i+1,i],:]
    return a,b

def glb_cell_fix_theta(filename,vname,theta):
    f,time,nblock,nx,ny=read_attr(filename)
    r=[]
    profile=[]
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        x_center,y_center=blk_center_polar(f,blkname)
        x_mesh,y_mesh=blk_mesh(f,blkname)
        if (theta>=y_mesh[0] and theta<=y_mesh[-1]):
            var=blk_var(f,blkname,vname)
            var_interpolator=interpolate.interp2d(x_center,y_center,var,kind='linear')
            var_interp=np.zeros(nx)
            for j in range(nx):
                var_interp[j]=var_interpolator(x_center[j],theta)
            r.append(x_center)
            profile.append(var_interp)
    r=np.asarray(np.ravel(r))
    profile=np.asarray(np.ravel(profile))
    return r,profile

def glb_inter_fix_theta(filename,vname,theta):
    f,time,nblock,nx,ny=read_attr(filename)
    r=[]
    profile=[]
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        x_center,y_center=blk_center_polar(f,blkname)
        x_mesh,y_mesh=blk_mesh(f,blkname)
        if (theta>=y_mesh[0] and theta<=y_mesh[-1]):
            var=blk_var(f,blkname,vname)
            var_interpolator=interpolate.interp2d(x_mesh,y_center,var,kind='linear')
            var_interp=np.zeros(nx+1)
            for j in range(nx+1):
                var_interp[j]=var_interpolator(x_mesh[j],theta)
            if (len(profile)==0):
                r.append(x_mesh)
                profile.append(var_interp)
                r=np.asarray(r)
                profile=np.asarray(profile)
            else:
                r=np.vstack((r,x_mesh))
                profile=np.vstack((profile,var_interp))
    r0=r[0,0]
    p0=profile[0,0]
    r=np.delete(r,0,1)
    profile=np.delete(profile,0,1)
    r,profile=array_sort(r,profile)
    r=np.ravel(r)
    profile=np.ravel(profile)
    r=np.insert(r,0,r0)
    profile=np.insert(profile,0,p0)
    return r,profile

def glb_surf_fix_r(filename,vname,r):
    f,time,nblock,nx,ny=read_attr(filename)
    theta=[]
    profile=[]
    cell=[]
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        mesh_x,mesh_y=blk_mesh(f,blkname)
        y_center=f[blkname]['y_center']
        if (r>=mesh_x[0] and r<mesh_x[-1]):
            var=blk_var(f,blkname,vname)
            var_interpolator=interpolate.interp2d(mesh_x,y_center,var,kind='linear')
            var_interp=np.zeros(ny)
            for j in range(ny):
                var_interp[j]=var_interpolator(r,y_center[j])
            theta.append(y_center)
            profile.append(var_interp)
            cell.append(mesh_y)
    theta=np.asarray(theta)
    profile=np.asarray(profile)
    cell=np.asarray(cell)
    theta,profile=array_sort(theta,profile)
    cell=np.unique(cell)
    theta=np.ravel(theta)
    cell=np.ravel(cell)
    profile=np.ravel(profile)
    return theta,profile,cell

#def glb_surf_fix_theta(filename,vname,theta):
#    f,time,nblock,nx,ny=read_attr(filename)
#    r=[]
#    profile=[]
#    for i in range(nblock):
#        blkname='blk'+str(i+1).zfill(5)
#        mesh_y=f[blkname]['mesh_y']
#        x_center=f[blkname]['x_center']
#        if (theta>=mesh_y[0] and theta<mesh_y[-1]):
#            var=blk_var(f,blkname,vname)
#            var_interpolator=interpolate.interp2d(x_center,mesh_y,var,kind='linear')
#            var_interp=np.zeros(ny)
#            for j in range(ny):
#                var_interp[j]=var_interpolator(r,y_center[j])
#            theta.append(y_center)
#            profile.append(var_interp)
#    theta=np.asarray(np.ravel(theta))
#    profile=np.asarray(np.ravel(profile))
#    return theta,profile

def glb_area_fix_r(filename,r):
    f,time,nblock,nx,ny=read_attr(filename)
    theta=[]
    surf=[]
    ds=np.zeros(ny)
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        y_center=f[blkname]['y_center']
        mesh_x=f[blkname]['mesh_x']
        if (r>=mesh_x[0] and r<mesh_x[-1]):
            mesh_y=np.asarray(f[blkname]['mesh_y'])
            for j in range(ny):
                ds[j]=r*r*2*pi*(cos(mesh_y[j])-cos(mesh_y[j+1]))
            if (len(surf)==0):
                theta.append(y_center)
                theta=np.asarray(theta)
                surf.append(ds)
                surf=np.asarray(surf)
            else:
                theta=np.vstack((theta,y_center))
                surf=np.vstack((surf,ds))
    theta,surf=array_sort(theta,surf)
    surf=np.ravel(surf)
    return surf

def glb_area_fix_theta(filename,theta):
    f,time,nblock,nx,ny=read_attr(filename)
    surf=[]
    ds=np.zeros(nx)
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        x_center,y_center=blk_center_polar(f,blkname)
        if (theta>y_center[0] and theta<y_center[-1]):
            mesh_x=np.asarray(f[blkname]['mesh_x'])
            for j in range(ny):
                ds[j]=0.5*(mesh_x[j+1]**2-mesh_x[j])*2*pi*sin(theta)
            surf.append(ds)
    surf=np.asarray(np.ravel(surf))
    return surf

def xlabel_sys(lscale):
    if (lscale=='norm'):
        xlabel=r'$L/L_{x}$'
    else:
        if ((abs(lscale-rsun)/rsun)<0.01):
            xlabel=r'$R_{\odot}$'
        elif (abs(lscale-rjupiter)/rjupiter<0.01):
            xlabel=r'$R_{j}$'
        elif (abs(lscale-au)/au<0.01):
            xlabel='AU'
        else:
            xlabel='cm'
    return xlabel

def ylabel_sys(lscale):
    if (lscale=='norm'):
        ylabel=r'$L/L_{y}$'
    else:
        if ((abs(lscale-rsun)/rsun)<0.01):
            ylabel=r'$R_{\odot}$'
        elif (abs(lscale-rjupiter)/rjupiter<0.01):
            ylabel=r'$R_{j}$'
        elif (abs(lscale-au)/au<0.01):
            ylabel='AU'
        else:
            ylabel='cm'
    return ylabel

def timesys(timescale):
    if (timescale==86400):
        value=' day'
    elif (abs((timescale-31536000)/31536000)<0.01):
        value=' year'
    else:
        value=''
    return value

def meshsys(x_mesh,y_mesh,domain,lscale):
    if (lscale=='norm'):
        x=x_mesh/(domain[1]-domain[0])
        y=y_mesh/(domain[3]-domain[2])
    else:
        x=x_mesh/lscale
        y=y_mesh/lscale
    return x,y

def scalesys(ax,lscale,domain):
    if (lscale=='norm'):
        ax.set_xlim(domain[0]/(domain[1]-domain[0]),domain[1]/(domain[1]-domain[0]))
        ax.set_ylim(domain[2]/(domain[3]-domain[2]),domain[3]/(domain[3]-domain[2]))
    else:
        ax.set_xlim(domain[0],domain[1])
        ax.set_ylim(domain[2],domain[3])

def velocitysys(vx,vy):
    mag=np.sqrt(vx**2+vy**2)
    minv=np.min(mag)
    length=np.log(mag/minv)
    vx=length*vx/mag
    vy=length*vy/mag
    return vx,vy

def velo_norm(vx,vy):
    mag=np.sqrt(vx**2+vy**2)
    vx=vx/mag
    vy=vy/mag
    return vx,vy,mag

def getcbname(cbname,var,vscale=''):
    if (cbname==''):
        if (vscale=='log'):
            value=r'$\log_{10}$'+cbname_dict.get(var)
        else:
            value=cbname_dict.get(var)
    else:
        value=cbname
    return value

def getfigname(figname,var):
    if (figname==''):
        value=var
    else:
        value=figname
    return value

def subplot_vars(filename,layout,plotsize,var_list,domain_list,lscale_list,vrange_list,vscale_list='',cbname_list='',figname='combine',timescale=1,coord='cartesian',   \
    func_list='',format='.png',suptitle=False,extra_title='',wdspace=0.14,htspace=0.06,timeunit='s',figdir='',sharey=False,extra_plot=False,dpi=100):
    fig,axes=plt.subplots(nrows=layout[0],ncols=layout[1],figsize=(plotsize[0],plotsize[1]),sharey=sharey)
    ifunc=0
    if (layout[0]==1 and layout[1]==1):
        if (var_list[0]=='block'):
            subplot_block(fig,axes,filename,domain_list[0],lscale_list[0],coord=coord)
        elif (var_list[0]=='func'):
            subplot_func(fig,axes,filename,func_list[0],domain_list[0],lscale_list[0],vrange_list[0],vscale_list[0],coord=coord)
            ifunc=ifunc+1
        else:
            subplot_var(fig,axes,filename,var_list[0],domain_list[0],lscale_list[0],vrange_list[0],vscale_list[0],coord=coord)
    elif (layout[0]==1):
        for j in range(layout[1]):
            ilist=j
            if (var_list[ilist]=='block'):
                subplot_block(fig,axes[j],filename,domain_list[ilist],lscale_list[ilist],coord=coord)
            elif (var_list[ilist]=='func'):
                subplot_func(fig,axes[j],filename,func_list[ifunc],domain_list[ilist],lscale_list[ilist],vrange_list[ilist],vscale_list[ilist],coord=coord)
                ifunc=ifunc+1
            else:
                subplot_var(fig,axes[j],filename,var_list[ilist],domain_list[ilist],lscale_list[ilist],vrange_list[ilist],vscale_list[ilist],coord=coord)
    elif (layout[1]==1):
        for i in range(layout[0]):
            ilist=i
            if (var_list[ilist]=='block'):
                subplot_block(fig,axes[i],filename,domain_list[ilist],lscale_list[ilist],coord=coord)
            elif (var_list[ilist]=='func'):
                subplot_func(fig,axes[i],filename,func_list[ifunc],domain_list[ilist],lscale_list[ilist],vrange_list[ilist],vscale_list[ilist],coord=coord)
                ifunc=ifunc+1
            else:
                subplot_var(fig,axes[i],filename,var_list[ilist],domain_list[ilist],lscale_list[ilist],vrange_list[ilist],vscale_list[ilist],coord=coord)
    else:
        for i in range(layout[0]):
            for j in range(layout[1]):
                ilist=i*(layout[1])+j
                if (var_list[ilist]=='block'):
                    subplot_block(fig,axes[i,j],filename,domain_list[ilist],lscale_list[ilist],coord=coord)
                elif (var_list[ilist]=='func'):
                    subplot_func(fig,axes[i,j],filename,func_list[ifunc],domain_list[ilist],lscale_list[ilist],vrange_list[ilist],vscale_list[ilist],coord=coord)
                    ifunc=ifunc+1
                elif (var_list[ilist]==''):
                    pass
                else:
                    subplot_var(fig,axes[i,j],filename,var_list[ilist],domain_list[ilist],lscale_list[ilist],vrange_list[ilist],vscale_list[ilist],coord=coord)
    f,time,nblock,nx,ny=read_attr(filename)
    f.close()
    if (suptitle==True):
        simtime=' Time='+"{:.3e}".format(time/timescale)+timeunit
        fig.suptitle(extra_title+simtime,fontsize =22,y=0.96)
    fig.tight_layout()
    fig.subplots_adjust(wspace=wdspace,hspace=htspace)
    if (extra_plot==True):
        return fig,axes
    else:
        fig.savefig(figdir+figname+'_'+filename+format,dpi=dpi)
        plt.close(fig)
        plt.close()

def subplot_var(fig,ax,filename,var,domain,lscale,vrange,vscale='log',cbname='',figname='',coord='cartesian'):
    f,time,nblock,nx,ny=read_attr(filename)
    maxv=-1e50
    minv=1e50
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        if (coord=='spherical'):
            x_mesh,y_mesh=polar_mesh(f,blkname)
        else:
            x_mesh,y_mesh=blk_mesh(f,blkname)
        x_mesh,y_mesh=meshsys(x_mesh,y_mesh,domain,lscale)
        vblk=f[blkname][var]
        if (vscale=='log'):
            vblk=np.log10(np.maximum(np.asarray(vblk),10**(vrange[0]-2)))
        maxv=max(maxv,np.max(vblk))
        minv=min(minv,np.min(vblk))
        if (coord=='spherical'):
            pcm=ax.pcolormesh(y_mesh,x_mesh,vblk,cmap='jet',vmin=vrange[0],vmax=vrange[1],shading='auto')
        else:
            pcm=ax.pcolormesh(x_mesh,y_mesh,vblk,cmap='jet',vmin=vrange[0],vmax=vrange[1],shading='auto')
    ax.set_box_aspect((domain[3]-domain[2])/(domain[1]-domain[0]))
    ax.set_xlabel(xlabel_sys(lscale))
    if (lscale!='norm'):
        if (coord=='cartesian'):
            ax.set_ylabel(ylabel_sys(lscale))
    scalesys(ax,lscale,domain)
    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size='2%', pad=0.4, pack_start = True)
    fig.add_axes(cax)
    fig.colorbar(pcm,cax=cax,label=getcbname(cbname,var,vscale),orientation="horizontal",shrink=0.8)
    f.close()
    #plt.close(fig)

def subplot_func(fig,ax,filename,func,domain,lscale,vrange,vscale='log',coord='cartesian'):
    f,time,nblock,nx,ny=read_attr(filename)
    maxv=-1e50
    minv=1e50
    for i in range(nblock):
        blkname='blk'+str(i+1).zfill(5)
        if (coord=='spherical'):
            x_mesh,y_mesh=polar_mesh(f,blkname)
        else:
            x_mesh,y_mesh=blk_mesh(f,blkname)
        x_mesh,y_mesh=meshsys(x_mesh,y_mesh,domain,lscale)
        cbname,vblk=func(f,blkname)
        if (vscale=='log'):
            vblk=np.log10(np.maximum(np.asarray(vblk),10**(vrange[0]-2)))
        maxv=max(maxv,np.max(vblk))
        minv=min(minv,np.min(vblk))
        if (coord=='spherical'):
            pcm=ax.pcolormesh(y_mesh,x_mesh,vblk,cmap='jet',vmin=vrange[0],vmax=vrange[1],shading='auto')
        else:
            pcm=ax.pcolormesh(x_mesh,y_mesh,vblk,cmap='jet',vmin=vrange[0],vmax=vrange[1],shading='auto')
    ax.set_box_aspect((domain[3]-domain[2])/(domain[1]-domain[0]))
    ax.set_xlabel(xlabel_sys(lscale))
    if (lscale!='norm'):
        if (coord=='cartesian'):
            ax.set_ylabel(ylabel_sys(lscale))
    scalesys(ax,lscale,domain)
    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size='2%', pad=0.4, pack_start = True)
    fig.add_axes(cax)
    fig.colorbar(pcm,cax=cax,label=cbname,orientation="horizontal",shrink=0.8)
    f.close()
    #plt.close(fig)

def reshuffle_cpuid(cpu_id):
    n=cpu_id/5
    cpu_id=cpu_id*4+n
    return cpu_id

def subplot_block(fig,ax,filename,domain,lscale,coord='cartesian'):
    f,time,nblock,nx,ny=read_attr(filename)
    cmap=matplotlib.colormaps['tab20c']
    if (coord=='cartesian'):
        for i in range(nblock):
            blkname='blk'+str(i+1).zfill(5)
            cpu_id=f[blkname]['cpu_id'][0]
            x_mesh,y_mesh=blk_mesh(f,blkname)
            if (lscale=='norm'):
                x_mesh=x_mesh/(domain[1]-domain[0])
                y_mesh=y_mesh/(domain[3]-domain[2])
            xlyl=[x_mesh[0],y_mesh[0]]
            xlen=x_mesh[nx]-x_mesh[0]
            ylen=y_mesh[ny]-y_mesh[0]
            y1=[y_mesh[0],y_mesh[-1]]
            x1=[x_mesh[0],x_mesh[-1]]
            cpu_id=reshuffle_cpuid(cpu_id)
            ax.add_patch(matplotlib.patches.Rectangle((xlyl[0],xlyl[1]),xlen,ylen,fill=True,color=cmap(0.5/20.0+(cpu_id%20)/20.0)))
            ax.add_patch(matplotlib.patches.Rectangle((xlyl[0],xlyl[1]),xlen,ylen,edgecolor='black',linewidth=1,fill=False))
            #ax.set_ylabel(ylabel_sys(lscale))
            ax.set_box_aspect((domain[3]-domain[2])/(domain[1]-domain[0]))
    elif (coord=='spherical'):
        for i in range(nblock):
            blkname='blk'+str(i+1).zfill(5)
            cpu_id=f[blkname]['cpu_id'][0]
            x_mesh,y_mesh=blk_mesh(f,blkname)
            x_mesh=x_mesh/lscale
            r_min=x_mesh[0]
            r_max=x_mesh[-1]
            theta_min=pi/2-y_mesh[-1]
            theta_max=pi/2-y_mesh[0]
            cpu_id=reshuffle_cpuid(cpu_id)
            ax.add_patch(matplotlib.patches.Wedge((0,0),r_max,degrees(theta_min),degrees(theta_max),width=r_max-r_min,fill=True,color=cmap(0.5/20.0+(cpu_id%20)/20.0)))
            ax.add_patch(matplotlib.patches.Wedge((0,0),r_max,degrees(theta_min),degrees(theta_max),width=r_max-r_min,edgecolor='black',linewidth=1,fill=False))
            ax.set_xlim(domain[0],domain[1])
            ax.set_ylim(domain[2],domain[3])
            ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(xlabel_sys(lscale))
    scalesys(ax,lscale,domain)
    f.close()
    #plt.close(fig)

def get_velocity(start,end,inner,domain,lscale,stream=True):
    u=[]
    v=[]
    for i in range(start,end):
        filename='cpd'+str(i).zfill(5)
        print(filename)
        f,time,nblock,nx,ny=read_attr(filename)
        for k in range(nblock):
            blkname='blk'+str(k+1).zfill(5)
            if k==0:
                x,y=polar_center(f,blkname)
                vx,vy=cartesian_v(f,blkname)
            else:
                xx,yy=polar_center(f,blkname)
                vxvx,vyvy=cartesian_v(f,blkname)
                x,y,vx,vy=np.vstack((x,xx)),np.vstack((y,yy)),np.vstack((vx,vxvx)),np.vstack((vy,vyvy))
        x,y,vx,vy=np.ravel(x),np.ravel(y),np.ravel(vx),np.ravel(vy)
        u.append(vx)
        v.append(vy)
        f.close()
    u,v=np.mean(u,axis=0),np.mean(v,axis=0)
    if stream:
        xx=np.linspace(0,domain[1]*rjupiter,100)
        yy=np.linspace(0,domain[1]*rjupiter,100)
        xx,yy=np.meshgrid(xx,yy)
        xx,yy=np.ravel(xx),np.ravel(yy)
    else:
        rr=np.linspace(inner*rjupiter,domain[1]*1.2*rjupiter,60)
        thetatheta=np.linspace(0,pi/2+pi/10,40)
        xx,yy=polar_to_cartesian(rr,thetatheta)
    points=np.vstack((y,x)).T
    xi=np.vstack((xx,yy)).T
    uu=griddata(points,u,xi,method='linear')
    vv=griddata(points,v,xi,method='linear')
    xx,yy=meshsys(xx,yy,domain,lscale)
    return xx,yy,uu,vv

def subplot_velocity(start,end,inner,domain,lscale,coord='cartesian',figname='velocity',format='.png',timescale=1,timeunit='s',extra_title=''):
    fig,ax=plt.subplots(figsize=(10,8))
    filename='cpd'+str(end-1).zfill(5)
    f,time,nblock,nx,ny=read_attr(filename)
    x,y,u,v=get_velocity(start,end,inner,domain,lscale,stream=False)
    u,v,mag=velo_norm(u,v)
    subplot_block(fig,ax,filename,domain,lscale,coord='cartesian')
    quiver=ax.quiver(x,y,u,v,mag,scale=30,scale_units='width',pivot='tip',cmap='viridis')
    cbar=fig.colorbar(quiver,ax=ax,orientation='horizontal',shrink=0.8)
    cbar.set_label('2D velocity '+r'[$cm\cdot s^{-1}$]')
    ax.set_xlabel(xlabel_sys(lscale))
    ax.set_ylabel(ylabel_sys(lscale))
    simtime=' Time='+"{:.3e}".format(time/timescale)+timeunit
    fig.suptitle(extra_title+simtime,fontsize =22,y=0.96)
    fig.tight_layout()
    fig.savefig(figname+'_'+str(domain[1])+'_'+filename+format,dpi=300)
    plt.close(fig)
    f.close()


def subplot_streamline(start,end,inner,domain,lscale,coord='cartesian',figname='streamline',format='.png',timescale=1,timeunit='s',extra_title=''):
    fig,ax=plt.subplots(figsize=(10,8))
    filename='cpd'+str(end-1).zfill(5)
    f,time,nblock,nx,ny=read_attr(filename)
    x,y,u,v=get_velocity(start,end,inner,domain,lscale,stream=True)
    x,y,u,v=x.reshape(100,100),y.reshape(100,100),u.reshape(100,100),v.reshape(100,100)
    u,v,mag=velo_norm(u,v)
    subplot_block(fig,ax,filename,domain,lscale,coord='cartesian')
    stream=ax.streamplot(x,y,u,v,color=mag,density=5,cmap='viridis',arrowsize=1)
    cbar=fig.colorbar(stream.lines,ax=ax,orientation='horizontal',shrink=0.8)
    cbar.set_label('2D velocity '+r'[$cm\cdot s^{-1}$]')
    ax.set_xlabel(xlabel_sys(lscale))
    ax.set_ylabel(ylabel_sys(lscale))
    simtime=' Time='+"{:.3e}".format(time/timescale)+timeunit
    fig.suptitle(extra_title+simtime,fontsize =22,y=0.96)
    fig.tight_layout()
    fig.savefig(figname+'_'+str(domain[1])+'_'+filename+format,dpi=300)
    plt.close(fig)
    f.close()

def extra_subplot_streamline(fig,axes,iframe,inner,domain,lscale,coord='cartesian',figdir='',figname='streamline',format='.png',timescale=1,timeunit='s',extra_title='',dpi=100):
    filename='cpd'+str(iframe).zfill(5)
    f,time,nblock,nx,ny=read_attr(filename)
    x,y,u,v=get_velocity(iframe,iframe+1,inner,domain,lscale,stream=True)
    x,y,u,v=x.reshape(100,100),y.reshape(100,100),u.reshape(100,100),v.reshape(100,100)
    u,v,mag=velo_norm(u,v)
    stream=axes.streamplot(x,y,u,v,color='grey',density=2,arrowsize=1)
    fig.savefig(figdir+figname+format,dpi=dpi)
    plt.close(fig)
    f.close()

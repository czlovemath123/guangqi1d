import sys
sys.path.insert(1,'../../scripts')
from assemble_1d_data import *
from math import *
import matplotlib
import matplotlib.pyplot as plt
import os.path
from scipy.interpolate import interp1d
plt.rcParams.update({'font.size': 18})
np.set_printoptions(threshold=sys.maxsize)

def read_history(filename):
    nl=len(open(filename).readlines(  ))
    f=open(filename,'r')
    line=f.readline()
    a=[float(s) for s in line.split()]
    da=len(a)
    history=np.zeros((nl,da))
    history[0]=a
    for i in range(nl-1):
        line=f.readline()
        a=[float(s) for s in line.split()]
        history[i+1]=a
    return history,nl


def single_model_evolution(workdir,figname='evolution.png'):
    print('single model kippenhahn:'+workdir)
    os.chdir(workdir)
    f_var=open('output_var_info.dat','r')
    n_var=int(f_var.readline())
    var_list=[]
    for i in range(n_var):
        line=f_var.readline()
        var_list.append(line.rstrip())
    line=f_var.readline()
    f_var.close()
    filehead=line.rstrip()
    nml=f90nml.read('global.data')
    xmin=nml['meshinfo']['n_domain'][0]
    xmax=nml['meshinfo']['n_domain'][1]
    xmin_rsun=xmin/rsun
    xmax_rsun=xmax/rsun
    nframe=nml['meshinfo']['nframe']
    tfinal=nml['meshinfo']['tfinal']
    X=nml['phyinfo']['h_ratio']
    dt=tfinal/nframe
    nml2=f90nml.read('problem.data')
    ms=nml2['parameters_1d']['m_star']
    cwd=os.getcwd()
    data_dir='/out'
    fig_dir=''
    os.chdir(cwd+data_dir)
    nblk=[]
    k=1
    filenumber=str(k).zfill(5)
    filename=filehead+filenumber+'.h5'
    r=glb_x_center(filename)
    ri=glb_x_inter(filename)
    t,rho,temp,arad,kr,mu,mupres,frad,xi,presgrad,vr=([] for _ in range(11))
    while os.path.isfile(filename):
        f,time,nblocks,nx=read_attr(filename)
        t.append(time)
        rho.append(glb_cell_var(filename,'rho'))
        temp.append(glb_cell_var(filename,'temp'))
        Erad=glb_cell_var(filename,'Erad')
        Frad=glb_inter_var(filename,'Fradx')
        aradx=glb_inter_var(filename,'aradx')
        pres=glb_cell_var(filename,'pres')
        mass=glb_cell_var(filename,'rho')
        nh2=mass*glb_cell_var(filename,'H2')/mh2
        nhi=mass*glb_cell_var(filename,'HI')/mh
        nhii=mass*2*glb_cell_var(filename,'HII')/mh
        nhei=mass*glb_cell_var(filename,'HeI')/mhe
        nheii=mass*2*glb_cell_var(filename,'HeII')/mhe
        nheiii=mass*3*glb_cell_var(filename,'HeIII')/mhe
        nparticle=X*(nh2+nhi+nhii)+(1-X)*(nhei+nheii+nheii+nheiii)
        mu_profile=mass/nparticle/mh
        pgrad=-np.gradient(pres,r)
        g=mass*ms*msun*G/r**2
        presgrad.append(pgrad/g)
        mu.append(mass/nparticle/mh)
        phi=ms*msun*G/r
        vx=glb_cell_var(filename,'vx')
        vr.append(vx)
        xi.append(vx**2/2/phi)
        eg=glb_cell_var(filename,'egv')
        Frad=Frad*4*pi*ri**2
        frad.append(Frad/1e37)
        kr.append(glb_cell_var(filename,'Rosseland'))
        erad_temp=np.power(Erad/a_rad,0.25)
        g=ms*msun*G/ri**2
        aradg=aradx/g
        arad.append(aradg)
        k=k+1
        filenumber=str(k).zfill(5)
        filename=filehead+filenumber+'.h5'
    t=np.asarray(t)
    t=t/day
    r=r/rsun/1000
    ri=ri/rsun/1000
    rmin=0.01
    rmax=2
    tmin=0
    tmax=50
    rho,temp,kr,arad,mu,presgrad,frad,xi,vr=[np.asarray(z) for z in [rho,temp,kr,arad,mu,presgrad,frad,xi,vr]]
    arrays=[rho,temp,kr]
    arrays=[np.log10(arr) for arr in arrays]
    rho,temp,kr=arrays

    X, Y = np.meshgrid(t, r)
    Xi, Yi= np.meshgrid(t, ri)
    fig, axes = plt.subplots(
        nrows=2,
        ncols=4,
        figsize=(18, 16),  # Unchanged
        sharey=True,
        gridspec_kw={'top': 0.95, 'wspace': 0.02, 'hspace': 0.02}
    )
    im = axes[0, 0].pcolormesh(X, Y, rho.T, shading='auto', cmap='jet')
    cbar = fig.colorbar(im, ax=axes[0,0], orientation='horizontal',shrink=0.9, pad=0.12, aspect=40)
    cbar.set_label(r'$\log_{10}\rho\ (\rm{g}\cdot\rm{cm}^{-3})$', labelpad=5)
    axes[0,0].set_ylim(rmin,rmax)
    axes[0,0].set_xlim(tmin,tmax)
    axes[0,0].set_ylabel(r'$r\ (\times1000R_{\odot})$')
    axes[0,0].set_xlabel(r'$t\ (day)$')

    contour_plot=axes[0,0].contour(
    X, Y, vr.T,  # vr should have the same shape as X and Y
    levels=[0],   # Only plot vr=0
    colors='black',
    linewidths=2,  # Adjust thickness
    linestyles='-',  # Solid line
)
    label_positions = [
    (t[80], r[1])  # Example: Use actual (t, r) values where the contour is visible
]
    axes[0,0].clabel(
    contour_plot, 
    inline=True,          # Place label inline with the contour
    fmt=r'$v=0$',            # Text to display
    fontsize=18,          # Adjust font size
    colors='black',       # Label color
    manual=label_positions          # Let matplotlib choose label positions
)

    im = axes[0,1].pcolormesh(Xi, Yi, frad.T, shading='auto', cmap='jet')
    cbar = fig.colorbar(im, ax=axes[0,1], orientation='horizontal',shrink=0.9, pad=0.12, aspect=40)
    cbar.set_label(r'$L_{37}\ (\times10^{37}\rm{erg}\cdot\rm{s}^{-1})$', labelpad=5)
    axes[0,1].set_ylim(rmin,rmax)
    axes[0,1].set_xlim(tmin,tmax)
    axes[0,1].set_xlabel(r'$t\ (day)$')

    im = axes[0,2].pcolormesh(X, Y, kr.T, shading='auto', cmap='jet',vmin=-2,vmax=1)
    cbar = fig.colorbar(im, ax=axes[0,2], orientation='horizontal',shrink=0.9, pad=0.12, aspect=40)
    cbar.set_label(r'$log_{10}\kappa_{\rm{R}}\ (\rm{cm}^{2}\cdot\rm{g}^{-1})$', labelpad=5)
    axes[0,2].set_ylim(rmin,rmax)
    axes[0,2].set_xlim(tmin,tmax)
    axes[0,2].set_xlabel(r'$t\ (day)$')

    im = axes[0,3].pcolormesh(Xi, Yi, arad.T, shading='auto', cmap='nipy_spectral',vmin=0,vmax=2)
    cbar = fig.colorbar(im, ax=axes[0,3], orientation='horizontal',shrink=0.9, pad=0.12, aspect=40)
    cbar.set_label(r'$a_{\rm{rad}}/g$', labelpad=5)
    axes[0,3].set_ylim(rmin,rmax)
    axes[0,3].set_xlim(tmin,tmax)
    axes[0,3].set_xlabel(r'$t\ (day)$')

    im = axes[1, 0].pcolormesh(X, Y, temp.T, shading='auto', cmap='jet')
    cbar = fig.colorbar(im, ax=axes[1,0], orientation='horizontal',shrink=0.9, pad=0.12, aspect=40)
    cbar.set_label(r'$\log_{10}T_{g}\ (\rm{K})$', labelpad=5)
    axes[1,0].set_ylim(rmin,rmax)
    axes[1,0].set_xlim(tmin,tmax)
    axes[1,0].set_ylabel(r'$r\ (\times1000R_{\odot})$')
    axes[1,0].set_xlabel(r'$t\ (day)$')

    im = axes[1,1].pcolormesh(X, Y, mu.T, shading='auto', cmap='jet')
    cbar = fig.colorbar(im, ax=axes[1,1], orientation='horizontal',shrink=0.9, pad=0.12, aspect=40)
    cbar.set_label(r'$\mu$', labelpad=5)
    axes[1,1].set_ylim(rmin,rmax)
    axes[1,1].set_xlim(tmin,tmax)
    axes[1,1].set_xlabel(r'$t\ (day)$')

    im = axes[1,2].pcolormesh(X, Y, xi.T, shading='auto', cmap='seismic',vmin=0,vmax=2)
    cbar = fig.colorbar(im, ax=axes[1,2], orientation='horizontal',shrink=0.9, pad=0.12, aspect=40)
    cbar.set_label(r'$\xi$', labelpad=5)
    axes[1,2].set_ylim(rmin,rmax)
    axes[1,2].set_xlim(tmin,tmax)
    axes[1,2].set_xlabel(r'$t\ (day)$')

    im = axes[1,3].pcolormesh(X, Y, presgrad.T, shading='auto', cmap='nipy_spectral',vmin=0,vmax=2)
    cbar = fig.colorbar(im, ax=axes[1,3], orientation='horizontal',shrink=0.9, pad=0.12, aspect=40)
    cbar.set_label(r'$-\nabla p/g$', labelpad=5)
    axes[1,3].set_ylim(rmin,rmax)
    axes[1,3].set_xlim(tmin,tmax)
    axes[1,3].set_xlabel(r'$t\ (day)$')
    for i in range(2):
        for j in range(4):
            axes[i,j].text(
            0.02,  # x position in axes coordinates (0=left, 1=right)
            0.98,  # y position in axes coordinates (0=bottom, 1=top)
            'p'+str(j+i*4),  # the text
            transform=axes[i,j].transAxes,  # use axes coordinates
            color='white',  # text color
            fontsize=18,  # adjust font size as needed
            verticalalignment='top',  # align to top
            bbox=dict(boxstyle='round', facecolor='black', alpha=1)  # optional: adds a semi-transparent black background
        )

    fig.subplots_adjust(bottom=0, top=0.98, left=0.05, right=0.95,wspace=0.25, hspace=0.2)
    plt.tight_layout()
    plt.savefig('../'+figname,bbox_inches='tight', pad_inches=0,dpi=300)
    plt.show()

def plot_evolution():
    imodel=99
    workdir='model'+str(imodel).zfill(2)
    figname='single_evolution.png'
    single_model_evolution(workdir,figname)

plot_evolution()

import os
import sys
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
from assemble_1d_data import *
from math import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
from matplotlib.lines import Line2D
import os.path
from scipy.interpolate import interp1d
plt.rcParams.update({
    'font.size': 18,
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman'],
})
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
    fig_dir='/pictures/'
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
        mu.append(mass/nparticle/mh)
        pgrad=-np.gradient(pres,r)
        g=mass*ms*msun*G/r**2
        presgrad.append(pgrad/g)
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
    rmin=0.45
    rmax=3
    tmin=0
    tmax=t[-1]
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
    axes[0,0].set_xlabel(r'$Time\ [day]$')

    contour_plot=axes[0,0].contour(
    X, Y, vr.T,  # vr should have the same shape as X and Y
    levels=[0],   # Only plot vr=0
    colors='black',
    linewidths=2,  # Adjust thickness
    linestyles='-',  # Solid line
)
    label_positions = [
    (t[len(t)//2], r[1])  # place label mid-run, robust to the number of frames
]
    axes[0,0].clabel(
    contour_plot, 
    inline=True,          # Place label inline with the contour
    fmt=r'$v_{r}=0$',            # Text to display
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
    cbar.set_label(r'$-(\partial p/\partial r)/g$', labelpad=5)
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
    plt.savefig('..'+fig_dir+figname,bbox_inches='tight', pad_inches=0,dpi=300)
    plt.show()
    os.chdir('../..')

def single_model_evo2(workdir, figname='rhoL.png'):
    print('single model rho & L: ' + workdir)
    os.chdir(workdir)
    f_var = open('output_var_info.dat', 'r')
    n_var = int(f_var.readline())
    var_list = []
    for i in range(n_var):
        line = f_var.readline()
        var_list.append(line.rstrip())
    line = f_var.readline()
    f_var.close()
    filehead = line.rstrip()
    nml = f90nml.read('global.data')
    xmin = nml['meshinfo']['n_domain'][0]
    xmax = nml['meshinfo']['n_domain'][1]
    xmin_rsun = xmin / rsun
    xmax_rsun = xmax / rsun
    nframe = nml['meshinfo']['nframe']
    tfinal = nml['meshinfo']['tfinal']
    X = nml['phyinfo']['h_ratio']
    dt = tfinal / nframe
    nml2 = f90nml.read('problem.data')
    ms = nml2['parameters_1d']['m_star']
    cwd = os.getcwd()
    data_dir = '/out'
    fig_dir = '/pictures/'
    os.chdir(cwd + data_dir)
    
    k = 1
    filenumber = str(k).zfill(5)
    filename = filehead + filenumber + '.h5'
    r = glb_x_center(filename)
    ri = glb_x_inter(filename)
    r_cm = r.copy()
    t, rho, temp, frad, kr = ([] for _ in range(5))

    while os.path.isfile(filename):
        f, time, nblocks, nx = read_attr(filename)
        t.append(time)
        rho.append(glb_cell_var(filename, 'rho'))
        temp.append(glb_cell_var(filename, 'temp'))
        kr.append(glb_cell_var(filename, 'Rosseland'))
        Frad = glb_inter_var(filename, 'Fradx')
        Frad = Frad * 4 * pi * ri**2
        frad.append(Frad / 1e37)
        k = k + 1
        filenumber = str(k).zfill(5)
        filename = filehead + filenumber + '.h5'

    # Shift the simulation to the observation's phase frame (t=0 = peak):
    # subtract the best-fit phase shift. The observation data is NOT shifted.
    RTL_SHIFT = 23.7  # fallback if optimize_best.json is absent
    try:
        import json
        _bj = os.path.join(os.path.dirname(cwd), 'optimize_best.json')
        _s = json.load(open(_bj)).get('shift')
        if _s is not None:
            RTL_SHIFT = float(_s)
    except (OSError, ValueError):
        pass
    t = np.asarray(t) / day - RTL_SHIFT
    r = r / rsun / 1000
    ri = ri / rsun / 1000
    rmin = 0.5
    rmax = 5
    tmin = t[0]
    tmax = t[-1]

    rho, temp, frad, kr = [np.asarray(z) for z in [rho, temp, frad, kr]]
    log10_rho = np.log10(rho)
    log10_temp = np.log10(temp)

    # Rosseland optical depth tau, integrated inward from the outer
    # boundary (tau=1e-5) using tau' = -rho*kappa_R along r_cm (cm).
    tau_outer = 1e-5
    integrand = rho * kr
    dr_cm = np.diff(r_cm)
    seg = 0.5 * (integrand[:, :-1] + integrand[:, 1:]) * dr_cm[np.newaxis, :]
    cumseg = np.cumsum(seg[:, ::-1], axis=1)[:, ::-1]
    tau = np.empty_like(integrand)
    tau[:, :-1] = tau_outer + cumseg
    tau[:, -1] = tau_outer
    log10_tau = np.log10(tau)

    # Local blackbody luminosity L_bb = 4*pi*r^2*sigma_sb*T^4 (r in cm),
    # compared against the actual luminosity at the outer boundary L_outer(t)
    # (frad evaluated at the last, outermost interface).
    L_outer = frad[:, -1]
    L_bb = 4 * pi * r_cm[np.newaxis, :]**2 * sigma_sb * temp**4 / 1e37
    log10_Lbb_ratio = np.log10(L_bb) - np.log10(L_outer)[:, np.newaxis]
    log10_Lbb_ratio[:, r > 2] = np.nan

    X, Y = np.meshgrid(t, r)
    Xi, Yi = np.meshgrid(t, ri)
    
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(16.5, 7.2),
        sharey=True,
        gridspec_kw={'top': 0.95, 'wspace': 0.05}
    )

    # Plot log10(rho)
    im0 = axes[0].pcolormesh(X, Y, log10_rho.T, shading='auto', cmap='inferno', vmin=-16, vmax=-8)
    cbar0 = fig.colorbar(im0, ax=axes[0], orientation='horizontal', shrink=0.9, pad=0.14, aspect=30)
    cbar0.set_label(r'$\log_{10}\rho\ (\rm{g}\cdot\rm{cm}^{-3})$', labelpad=5)
    axes[0].set_ylim(rmin, rmax)
    axes[0].set_xlim(tmin, tmax)
    axes[0].set_ylabel(r'$r\ [\times1000R_{\odot}]$')
    axes[0].set_xlabel('Time [day]')
    axes[0].set_xticks([0, 50, 100])

    # White annotations with a thin black halo read well over 'inferno'.
    HALO = [patheffects.withStroke(linewidth=3, foreground='black')]

    # Temperature contours (6000 K solid, 3000 K dashed, 1500 K dotted),
    # labelled by a single legend in the black upper-center region.
    for _level, _ls in [(6000, '-'), (3000, '--'), (1500, ':')]:
        contour_plot = axes[0].contour(
            X, Y, log10_temp.T,
            levels=[np.log10(_level)],
            colors='white',
            linewidths=2,
            linestyles=_ls
        )
        contour_plot.set_path_effects(HALO)
    temp_handles = [
        Line2D([], [], color='white', linewidth=2, linestyle=':',
               label=r'$T=1500\,\rm{K}$'),
        Line2D([], [], color='white', linewidth=2, linestyle='--',
               label=r'$T=3000\,\rm{K}$'),
        Line2D([], [], color='white', linewidth=2, linestyle='-',
               label=r'$T=6000\,\rm{K}$'),
    ]
    _legtemp = axes[0].legend(handles=temp_handles, loc='upper center',
                              frameon=True, framealpha=0.6,
                              facecolor='black', edgecolor='none',
                              fontsize=21.6)
    for _txt in _legtemp.get_texts():
        _txt.set_color('white')
        _txt.set_path_effects(HALO)
    for _h in _legtemp.legend_handles:
        _h.set_path_effects(HALO)

    axes[0].text(50, 1.8, r'$\rm{H}^{-}$', color='white', fontsize=23.76,
                 ha='center', va='top', path_effects=HALO)
    axes[0].text(100, 3, 'molecule', color='white', fontsize=23.76,
                 ha='center', va='top', path_effects=HALO)

    # Plot log10(T)
    im1 = axes[1].pcolormesh(X, Y, log10_temp.T, shading='auto', cmap='inferno', vmin=2.9, vmax=4.8)
    cbar1 = fig.colorbar(im1, ax=axes[1], orientation='horizontal', shrink=0.9, pad=0.14, aspect=30)
    cbar1.set_label(r'$\log_{10}T_{g}\ (\rm{K})$', labelpad=5)
    axes[1].set_ylim(rmin, rmax)
    axes[1].set_xlim(tmin, tmax)
    axes[1].set_xlabel('Time [day]')

    # Contour where local blackbody luminosity equals L at the outer boundary
    contour_plot_Lbb = axes[1].contour(
        X, Y, log10_Lbb_ratio.T,
        levels=[0],
        colors='white',
        linewidths=2,
        linestyles='-'
    )
    contour_plot_Lbb.set_path_effects(HALO)

    # Contour of tau_Rosseland = 2/3 (photosphere)
    contour_plot_tau = axes[1].contour(
        X, Y, log10_tau.T,
        levels=[np.log10(2/3)],
        colors='white',
        linewidths=2,
        linestyles='--'
    )
    contour_plot_tau.set_path_effects(HALO)

    # Observed blackbody radius Rbb (post_UVOIR_RTL_fin.txt); plotted on the
    # same phase axis as the (already-shifted) simulation. Radius in 1000 R_sun.
    _obs_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'post_UVOIR_RTL_fin.txt')
    if os.path.isfile(_obs_path):
        _obs = np.loadtxt(_obs_path)
        _ot = _obs[:, 0]
        _orb = _obs[:, 1] / rsun / 1000
        _orbe = _obs[:, 2] / rsun / 1000
        eb = axes[1].errorbar(_ot, _orb, yerr=_orbe, fmt='o', color='white',
                              markerfacecolor='white', markeredgecolor='black',
                              markeredgewidth=1.0, markersize=5,
                              ecolor='white', elinewidth=1.3, capsize=2.5,
                              linestyle='none', zorder=6, label=r'$R_{\rm bb}$')
        for _art in [eb[0]] + list(eb[1]) + list(eb[2]):
            _art.set_path_effects(HALO)

    # Legend for the r_em and tau=2/3 contours plus the observed Rbb points
    # (with their error-bar symbol), in the upper right corner.
    b_handles = [
        Line2D([], [], color='white', linewidth=2, linestyle='-',
               label=r'$r_{\rm em}$'),
        Line2D([], [], color='white', linewidth=2, linestyle='--',
               label=r'$\tau_{\rm R}=2/3$'),
    ]
    if os.path.isfile(_obs_path):
        b_handles.append(eb)
    _legrbb = axes[1].legend(handles=b_handles, loc='upper right',
                             frameon=True, framealpha=1.0,
                             facecolor='black', edgecolor='none',
                             fontsize=21.6)
    for _txt in _legrbb.get_texts():
        _txt.set_color('white')
        _txt.set_path_effects(HALO)
    for _h in _legrbb.legend_handles:
        _h.set_path_effects(HALO)

    # Plot L
    im2 = axes[2].pcolormesh(Xi, Yi, frad.T, shading='auto', cmap='inferno', vmin=0, vmax=500)
    cbar2 = fig.colorbar(im2, ax=axes[2], orientation='horizontal', shrink=0.9, pad=0.14, aspect=30)
    cbar2.set_label(r'$L_{37}\ [\times10^{37}\rm{erg}\cdot\rm{s}^{-1}]$', labelpad=5)
    axes[2].set_ylim(rmin, rmax)
    axes[2].set_xlim(tmin, tmax)
    axes[2].set_xlabel('Time [day]')

    # 0.1 * Eddington luminosity contour. L_Edd = 4*pi*G*M*c/kappa_es with
    # kappa_es = 0.34 cm^2/g (solar-composition electron scattering) and
    # M = m_star. frad is luminosity in units of 1e37 erg/s.
    _ledd = 4 * pi * G * (ms * msun) * c_light / 0.34
    _level_01ledd = _ledd / 1e37
    contour_plot_ledd = axes[2].contour(
        Xi, Yi, frad.T,
        levels=[_level_01ledd],
        colors='white',
        linewidths=2,
        linestyles='-'
    )
    contour_plot_ledd.set_path_effects(HALO)

    # Panel text labels
    for j in range(3):
        axes[j].text(
            0.04, 0.97,
            r'\textbf{' + chr(ord('a') + j) + '}',
            transform=axes[j].transAxes,
            color='white',
            fontsize=26,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='black', alpha=1)
        )

    # Common time ticks (phase, days) on all three panels.
    for _ax in axes:
        _ax.set_xticks([0, 50, 100])

    plt.tight_layout()
    png_path = '..' + fig_dir + figname
    plt.savefig(png_path, bbox_inches='tight', pad_inches=0, dpi=300)
    plt.show()
    os.chdir('../..')

def best_model_name():
    """Return the best-fit model folder name from optimize_best.json."""
    import json
    here=os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here,'optimize_best.json')) as f:
        return json.load(f)['name']

def plot_evolution(imodel=None):
    if imodel is None:
        workdir=best_model_name()
    else:
        workdir='model'+str(int(imodel)).zfill(5)
    figname = 'single_evolution.png'
    figname2 = 'rhoL.png'
    #single_model_evolution(workdir, figname)
    single_model_evo2(workdir, figname2)

def multiple_model_evolution(workdirs,figname='arad.png',subtitles=0,uylims=[0.3,10],figsize=(18,8),tmax=500,choice='arad'):
    print('multiple:', workdirs)

    # Create subplots (keep original figsize)
    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(workdirs),
        figsize=figsize,  # Unchanged
        sharey=True,
        gridspec_kw={'top': 0.78, 'wspace': 0.02}
    )

    filehead = '1dlrne'
    all_z = []

    for i in range(len(workdirs)):
        os.chdir(workdirs[i])
        nml2 = f90nml.read('problem.data')
        ms = nml2['parameters_1d']['m_star']
        os.chdir('out')

        # Read data files (unchanged)
        k = 1
        filenumber = str(k).zfill(5)
        filename = filehead + filenumber + '.h5'
        r = glb_x_center(filename)
        ri = glb_x_inter(filename)
        t, z = [], []

        while os.path.isfile(filename):
            f, time, nblocks, nx = read_attr(filename)
            t.append(time)
            if (choice=='arad'):
                aradx = glb_inter_var(filename, 'aradx')
                g = ms * msun * G / ri**2
                aradg = aradx / g
                z.append(aradg)
            elif (choice=='pres'):
                rho=glb_cell_var(filename,'rho')
                pres=glb_cell_var(filename,'pres')
                pgrad=-np.gradient(pres,r)
                g=rho*ms*msun*G/r**2
                z.append(pgrad/g)
            elif (choice=='all'):
                aradx = glb_inter_var(filename, 'aradx')
                g = ms * msun * G / ri**2
                aradg = aradx / g
                rho=glb_cell_var(filename,'rho')
                pres=glb_cell_var(filename,'pres')
                pgrad=-np.gradient(pres,r)
                g=rho*ms*msun*G/r**2
                pgradg=pgrad/g
                interp_func = interp1d(
                    ri,           # Original grid (interface)
                    aradg,        # Data to interpolate
                    axis=0,       # Interpolate along radial dimension
                    bounds_error=False,
                    fill_value="extrapolate"
                )
                aradg_interp = interp_func(r)  # Interpolated to cell centers
                atotal=aradg_interp+pgradg
                z.append(atotal)
            elif (choice=='prad'):
                rho=glb_cell_var(filename,'rho')
                temp=glb_cell_var(filename,'temp')
                prad=a_rad*temp**4/3
                pgrad=-np.gradient(prad,r)
                g=rho*ms*msun*G/r**2
                z.append(pgrad/g)
            k += 1
            filenumber = str(k).zfill(5)
            filename = filehead + filenumber + '.h5'

        # Process data (unchanged)
        t = np.asarray(t) / day
        r = r / rsun / 1000
        ri = ri / rsun / 1000
        z = np.asarray(z)
        all_z.append(z)

        # Plot (unchanged)
        if (choice=='arad'):
            Xi, Yi = np.meshgrid(t, ri)
        elif (choice=='pres' or choice=='all' or choice=='prad'):
            Xi, Yi = np.meshgrid(t, r)
        im = axes[i].pcolormesh(
            Xi, Yi, z.T,
            shading='auto',
            cmap='nipy_spectral',
            vmin=0,
            vmax=2
        )

        # Set axis limits and labels (unchanged)
        axes[i].set_ylim(uylims[0],uylims[1])
        axes[i].set_xlim(0.1, tmax)
        axes[i].set_xlabel(r'$t$'+' (day)')

        # Add subtitle (NEW)
        if (subtitles!=0):
            axes[i].set_title(subtitles[i],pad=10)  # y=1.02 shifts title up

        if i == 0:
            axes[i].set_ylabel(r'$r\ [\times1000R_{\odot}]$')

        os.chdir('../..')

    # Colorbar (unchanged)
    cbar_ax = fig.add_axes([0.15, 0.9, 0.7, 0.02])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.ax.xaxis.set_label_position('top')
    if (choice=='arad'):
        cbar.set_label(r'$a_{\rm{rad}}/g$', labelpad=10)
    elif (choice=='pres'):
        cbar.set_label(r'$-(\partial p/\partial r)/g$', labelpad=10)
    elif (choice=='all'):
        cbar.set_label(r'$[a_{\rm{rad}}-(\partial p/\partial r)]/g$', labelpad=10)
    elif (choice=='prad'):
        cbar.set_label(r'$p_{\rm{rad}}/p$', labelpad=10)

    # Adjust layout to avoid overlap (NEW: reduced top margin)
    plt.tight_layout(rect=[0, 0, 1, 0.9])  # Original: [0, 0, 1, 1]
    plt.savefig('pictures/'+figname, bbox_inches='tight', pad_inches=0, dpi=300)
    plt.show()


# Usage:
#   python kippenhahn2.py        -> plot the best-fit model (from optimize_best.json)
#   python kippenhahn2.py NNN    -> plot modelNNN (NNN zero-padded to 5 digits)
if len(sys.argv) > 1:
    plot_evolution(sys.argv[1])
else:
    plot_evolution()


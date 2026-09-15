from math import *
import sys
sys.path.insert(1,'../../scripts')
from phy_const import *
from eos import *
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import f90nml
import copy
import shutil
plt.rcParams.update({'font.size': 20})

class Config:

    x = 0.74
    rho_floor = 1e-17
    
    # Simulation parameters
    n = 10000
    ms = 8                    # Msun
    rin = 500                   # Inner boundary (Rsun) 
    rout = 10000
    tfinal = 160
    nframe = 100
    nx = 512*2
    xratio = 80
    level = 1
    ejt = 32                   # Total ejection time (days)
    gamma_eos = 1.4
    iradiation = 4
    
    # Model parameters
    # (model folder is auto-numbered, see next_model_directory())
    bound = 'reflective'
    model_eos = 'real'


cfg = Config()
amu = 1.66053886e-24
mh2 = 3.3466E-24
mh = 1.6733e-24
me = 9.1093897e-28
mhion = mh-me
mhe = 6.646481526e-24
ionh = 2.18e-11
dish = 7.17e-12
ionhe1 = 3.94e-11
ionhe2 = 8.72e-11
kb = 1.380658E-16
h = 6.6260755E-27
a_rad = 7.56e-15
day = 86400
yr = 3.154e7
G = 6.67259e-8
msun = 1.99e33
rsun = 6.96e10
year = 31536000
x = 0.74
n = 10000




# ---------------------------------------------------------------------------
# Equation of state
# ---------------------------------------------------------------------------
# All EOS functions (zh*, e_internal, p_internal, solve_species_n,
# *_species_state*, adiabatic_cs, ...) come from eos.py via `from eos import *`
# at the top of this file -- that is the authoritative implementation and must
# not be duplicated here.
#
# The adiabatic-index table (gamma.txt, ~18 MB) is ONLY needed to compute the
# gas sound speed for the Mach-number plot, so it is loaded lazily on first use
# and cached. Model generation (the hot path in the optimizer) never touches
# it, so no per-run cost is incurred.
GAMMA_TABLE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'gamma.txt')
_GAMMA_TABLE = None   # cached (rho_axis, t_axis, table) after first load


def _get_gamma_table():
    global _GAMMA_TABLE
    if _GAMMA_TABLE is None:
        _GAMMA_TABLE = load_gamma_table(GAMMA_TABLE_FILE)
    return _GAMMA_TABLE


def sound_speed(rho, t, x=x):
    """Adiabatic gas sound speed [cm/s] from eos.adiabatic_cs + the gamma table.

    Loads gamma.txt on first call only (used by the Mach-number plot).
    """
    rho_axis, t_axis, table = _get_gamma_table()
    return adiabatic_cs(rho, t, rho_axis, t_axis, table)



 
def parametricTemp(t, Tpeak=0.75e5 - 0.015e5):

    t1 = 4.1
    t2 = 25  #25

    T0     = 0.15e5
    Tmid   = 0.86e5 - 0.05e5 + 0.04e5
    Tfinal = 0.30e5  
    tau2   = 3 
    delta1 = 1.5
    delta2 = 1.5

    m1 = (Tpeak - T0) / t1
    m2 = (Tmid - Tpeak) / (t2 - t1)
    m3 = (Tfinal - Tmid) / tau2   

    if t < t1 - delta1:
        return T0 + m1 * t

    elif t < t1 + delta1:
        tL = t1 - delta1
        tR = t1 + delta1
        
        yL = T0 + m1 * tL
        yR = Tpeak + m2 * (tR - t1)

        s = (t - tL) / (tR - tL)
        h00 = 2*s**3 - 3*s**2 + 1
        h10 = s**3 - 2*s**2 + s
        h01 = -2*s**3 + 3*s**2
        h11 = s**3 - s**2

        return (h00*yL +
                h10*(tR - tL)*m1 +
                h01*yR +
                h11*(tR - tL)*m2)

    elif t < t2 - delta2:
        return Tpeak + m2 * (t - t1)

    elif t < t2 + delta2:
        tL = t2 - delta2
        tR = t2 + delta2
        
        yL = Tpeak + m2 * (tL - t1)
        yR = Tmid + (Tfinal - Tmid) * (1 - np.exp(-(tR - t2)/tau2))
        s = (t - tL) / (tR - tL)
        h00 = 2*s**3 - 3*s**2 + 1
        h10 = s**3 - 2*s**2 + s
        h01 = -2*s**3 + 3*s**2
        h11 = s**3 - s**2

        return (h00*yL +
                h10*(tR - tL)*m2 +
                h01*yR +
                h11*(tR - tL)*m3)

    else:
        return Tmid + (Tfinal - Tmid) * (1 - np.exp(-(t - t2) / tau2))
    



def vej(t, mod):
   
    if t < mod.t1:
        if t < mod.x2 * mod.t1:
            value = mod.v1_floor1 + (mod.v1 - mod.v1_floor1) * (t / (mod.x2 * mod.t1))
        else:
            value = mod.v1_floor2 + (mod.v1 - mod.v1_floor2) * exp(-(mod.x1 * (t - mod.x2 * mod.t1) / mod.t1)**2)
    else:
        vinitial = mod.v1_floor2 + (mod.v1 - mod.v1_floor2) * exp(-(mod.x1 * (1.0 - mod.x2))**2)
        value = mod.v2 + (vinitial - mod.v2) * exp(-(t - mod.t1) / mod.dt2) - mod.v2decline * (t - mod.t1) / mod.t2
    
    value = value * sqrt(2 * mod.ms * G / mod.r)
    return value * mod.vej_factor  # overall velocity scaling (was hardcoded 1.45)



def massloss(t, mod):
    dt = mod.tm1 / 30
    dm = mod.mdot2 - mod.mdot1 - mod.mdot1rise

    if t < mod.tm1:
        xVal = (t - mod.tm1) / dt
        value = mod.mdot1 + dm * exp(xVal) / (exp(xVal) + 1) + mod.mdot1rise * (t / mod.tm1)**2
        return value
    
    elif t < mod.tm2:
        xVal = (mod.tm1 - t) / dt
        value = mod.mdot2 - dm * exp(xVal) / (exp(xVal) + 1) - mod.mdot2decline * (t - mod.tm1) / mod.tm2
        return value
    
    else:
        t_mid = mod.tm2
        xVal_mid = (mod.tm1 - t_mid) / dt
        value_mid = mod.mdot2 - dm * exp(xVal_mid) / (exp(xVal_mid) + 1) - mod.mdot2decline * (t_mid - mod.tm1) / mod.tm2
        decay_timescale = mod.tm3 
        mdot_final = mod.mdot3      
        value = mdot_final + (value_mid - mdot_final) * exp(-(t - mod.tm2) / decay_timescale)
        return value


# Model Class
class model_class:
    def __init__(self, ms, rin, time, n, bound, model_eos):
        self.ms = ms
        self.rin = rin
        self.time = time
        self.t = np.linspace(0, time, n)
        self.dt = time/n*day
        self.bound = bound
        self.diffused_rho = cfg.rho_floor
        self.eos = model_eos
        self.header = ['ms', 'rin', 'time', 'v', 'mdot', 'eratio', 'bound', 'diffused_rho', 'model_eos']
        self.value = [self.ms, self.rin, self.time, 0.8, 0.5, 0.8, self.bound, self.diffused_rho, self.eos]

def init_model_csv(header):
    file_path = 'models.csv'
    df = pd.DataFrame(columns=header)
    df.to_csv(file_path, index=False)
    print(f"Initialized new {file_path} with header: {header}")


class fittingModel:
    # Class-level overrides hook: if set (dict), every fittingModel() built
    # anywhere in this module picks them up. The optimizer sets this so that
    # save_model()'s internally-constructed fittingModel() is also tuned.
    PARAM_OVERRIDES = None

    def __init__(self, overrides=None):
        if overrides is None:
            overrides = fittingModel.PARAM_OVERRIDES
        # Use parameters from configuration
        self.ms = cfg.ms * msun  # Convert to grams
        self.r = cfg.rin * rsun  # Use unified rin, convert to cm
        
        '''
        self.vstart = 1.21
        self.k = 0.007
        '''
        
        self.t1 = 10     
        self.t2 = 20.0
       

        self.dt2 = 2.1  
        self.v1 = 1.0   

        vv = 0.0  
        self.v1_floor1 = 0.97 - 0.15  
        self.v1_floor2 = 0.90
        self.v2 = 0.92    
        self.x1 = 6  
        self.x2 = 0.55 + 0.44 
        self.v2decline =  0.32  
        

        self.tm1 = 8.0 - 1
        self.tm2 = 0.8 * cfg.ejt   # plateau->decay turnover, tied to injection window
        self.tm3 = 1.0

        self.mdot1 = 3.0
        self.mdot2 = 8
        self.mdot1rise = 0.0001
        self.mdot2decline = -10
        self.mdot3 = 2.5

        # Overall scaling factors (promoted from hardcoded literals so they
        # can be tuned). vej_factor multiplies the ejecta velocity; temp_factor
        # multiplies parametricTemp; Tpeak is the temperature peak in K.
        self.vej_factor = 1.45
        self.temp_factor = 0.85
        self.Tpeak = 0.75e5 - 0.015e5

        # Optional parameter overrides (used by the optimizer). Any attribute
        # set above can be overridden by passing PARAM_OVERRIDES, an env-driven
        # dict, or by mutating the instance before generation.
        for key, val in (overrides or {}).items():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise AttributeError(f"fittingModel has no parameter '{key}'")

        self.v_esc = sqrt(2 * G * self.ms / self.r)
        print(f"Escape velocity: {self.v_esc:.2e} cm/s")
        print(f"Stellar radius: {self.r/rsun:.1f} Rsun")
    
    def vej(self, t):
        
        return vej(t, self)
    
    def massloss(self, t):
        
        return massloss(t, self) * msun / year  # Convert to g/s
    
    def density(self, t):
        
        try:
            Mdot = self.massloss(t)
            v_ej = self.vej(t)
            
            if v_ej <= 0:
                v_ej = 1e5
                
            rho = Mdot / (4 * pi * self.r**2 * v_ej)
            rho = max(1e-20, min(rho, 1e-5))
            return rho
        except:
            return 1e-17
    
    def temperature(self, t):

        try:
            return parametricTemp(t, self.Tpeak) * self.temp_factor
        except:
            return 1e4
    
    def eradEgasRatio(self, t):
        
        try:
            rho = self.density(t)
            T = self.temperature(t)
            e_gas = e_internal(rho, T, x)
            if e_gas <= 0:
                return 0.0
                
            e_rad = a_rad * T**4
            ratio = e_rad / e_gas
            return max(0.0, min(ratio, 100.0))
        except:
            return 0.1


def save_model(model, modeldir):
    file_path = 'models.csv'
    header = model.header
    new_model = model.value

    if not os.path.isfile(file_path):
        df = pd.DataFrame(columns=header)
        df.to_csv(file_path, index=False)

    df = pd.DataFrame([new_model], columns=header)
    df.to_csv(file_path, mode='a', header=False, index=False)
    print(f"model appended: {new_model}")

    os.makedirs(modeldir, exist_ok=True)
    os.makedirs(modeldir+'/out', exist_ok=True)
    os.makedirs(modeldir+'/pictures', exist_ok=True)
    source_file = 'output_var_info.dat'
    destination = modeldir+'/output_var_info.dat'
    shutil.copy(source_file, destination)
    source_file = 'guangqi'
    destination = modeldir+'/'+modeldir
    shutil.copy(source_file, destination)
    print(f"Directory '{modeldir}' created successfully!")
    
    formatted_string = df.to_string(col_space=len(header), index=False)
    with open(modeldir+'/formatted_data.txt', 'w') as file:
        file.write(formatted_string)
    
    n = len(model.t)
    rho = np.zeros(n)
    v_arr = np.zeros(n)
    temp = np.zeros(n)
    mdot_arr = np.zeros(n)
    eratio_arr = np.zeros(n)
    asym = np.zeros(n)
    
    fitting_model = fittingModel()
    
    print("Using fitting model to calculate physical quantities...")
    
    for i in range(n):
        t_val = model.t[i]
        # Use fitting model for calculations
        v_arr[i] = fitting_model.vej(t_val)
        mdot_arr[i] = fitting_model.massloss(t_val)
        rho[i] = fitting_model.density(t_val)
        temp[i] = fitting_model.temperature(t_val)
        eratio_arr[i] = fitting_model.eradEgasRatio(t_val)
        asym[i] = 0.0
        
    # Calculate cumulative mass
    culmdot = np.zeros(n)
    for i in range(n):
        if i == 0:
            culmdot[i] = mdot_arr[i] * model.dt / msun  # Convert to solar mass
        else:
            culmdot[i] = culmdot[i-1] + mdot_arr[i] * model.dt / msun
    
    print(f"Cumulative mass range: {np.min(culmdot):.6f} - {np.max(culmdot):.6f} M☉")
    
    # Check results
    print(f"Final statistics:")
    print(f"Velocity range: {np.min(v_arr):.2e} - {np.max(v_arr):.2e} cm/s")
    print(f"Mass loss rate range: {np.min(mdot_arr):.2e} - {np.max(mdot_arr):.2e} g/s")
    print(f"Density range: {np.min(rho):.2e} - {np.max(rho):.2e} g/cm³")
    print(f"Temperature range: {np.min(temp):.2e} - {np.max(temp):.2e} K")
    print(f"eratio range: {np.min(eratio_arr):.2e} - {np.max(eratio_arr):.2e}")
    
    # Generate bcinput
    data = [model.t * day, v_arr, rho, temp, eratio_arr, asym]
    data = np.transpose(data)
    mm = np.shape(data)[0]
    nn = np.shape(data)[1]
    
    headers = ['time', 'v', 'rho', 'temp', 'eratio', 'asym']
    
    print(f"Generated {nn} columns of data: {headers}")
    
    with open(modeldir+'/bcinput.dat', 'w') as file:
        for header in headers:
            file.write(f"{header:>16}")
        file.write("\n")
        for i in range(mm):
            for j in range(nn):
                file.write(f"{data[i, j]:16.8E}")
            file.write("\n")
    
    print(f"Successfully generated bcinput.dat file with {mm} rows")
    return model.t, v_arr, rho, temp, eratio_arr, mdot_arr, culmdot

def compute_bc_derived(t, v, rho, temp, mdot):
    """Compute the derived boundary-condition quantities plotted by
    plot_bcinput (Mach number, energy ratios, energy outflow, ...).

    Returns a dict with keys: mach, erad, eg, prad_pgas, cumE, vesc.
    Factored out of plot_bcinput so several models can be overlaid on one
    figure without duplicating the analytic EOS below.
    """
    # ---------------------------------------------------------------------------
    # Vectorized analytic EOS for the ejecta temperature range (2–10 × 10^4 K).
    # This avoids the scalar Python loops and sympy polynomial solvers in eos.py.
    #
    # H:  T > ~500 K so H2 is dissociated; use HI + HII (Saha quadratic).
    # He: T in range so use HeI + HeII (Saha quadratic).
    # Both branches are analytic — no sympy, no bisection.
    # ---------------------------------------------------------------------------
    xH  = x          # hydrogen mass fraction  (0.74)
    xHe = 1.0 - x    # helium  mass fraction   (0.26)

    rhoH  = xH  * rho   # [g/cm^3]
    rhoHe = xHe * rho

    _ztr_v = lambda T, m: (2.0 * np.pi * m * kb * T)**1.5 / h**3

    # K_H = ztr(mhion) * exp(-(dish+2*ionh)/2kT) * ztr(me) / (ztr(mh)*exp(-dish/2kT))
    #      = (mhion/mh)^1.5 * ztr(me) * exp(-ionh/kT)
    qion_H  = (mhion / mh)**1.5 * _ztr_v(temp, me) * np.exp(-ionh  / (kb * temp))

    # K_He = ztr(mhe)*exp(-ionhe1/kT) * ztr(me) / ztr(mhe)  =  ztr(me)*exp(-ionhe1/kT)
    qion_He = _ztr_v(temp, me) * np.exp(-ionhe1 / (kb * temp))

    # H: quadratic Saha  nhII^2 + qion*nhII - qion*nhtot = 0
    nhtot   = rhoH  / mh
    nhII    = (-qion_H + np.sqrt(qion_H**2 + 4.0 * qion_H * nhtot)) / 2.0
    nhII    = np.maximum(nhII, 0.0)
    nhI     = np.maximum(nhtot - nhII, 0.0)
    nhelec  = nhII.copy()

    # He: quadratic Saha  nheII^2 + qion*nheII - qion*nhetot = 0
    nhetot  = rhoHe / mhe
    nheII   = (-qion_He + np.sqrt(qion_He**2 + 4.0 * qion_He * nhetot)) / 2.0
    nheII   = np.maximum(nheII, 0.0)
    nheI    = np.maximum(nhetot - nheII, 0.0)
    nheelec = nheII.copy()

    # Gas pressure: p = n_total * kB * T
    ntot = nhI + nhII + nhelec + nheI + nheII + nheelec
    pgas = ntot * kb * temp

    # Internal energy density: e = sum(n_i * epsilon_i)
    eg = (nhI  * (1.5 * kb * temp + dish / 2.0)
        + nhII * (1.5 * kb * temp + dish / 2.0 + ionh)
        + nhelec * 1.5 * kb * temp
        + nheI  * 1.5 * kb * temp
        + nheII * (1.5 * kb * temp + ionhe1)
        + nheelec * 1.5 * kb * temp)

    # Radiation energy density and pressure
    erad = a_rad * temp**4
    prad = erad / 3.0

    # Sound speed: use gamma table (vectorized interpolation)
    rho_ax, t_ax, gtable = _get_gamma_table()
    from scipy.interpolate import RegularGridInterpolator as _RGI
    _ginterp = _RGI((t_ax, rho_ax), gtable, bounds_error=False, fill_value=None)
    logrho = np.log10(np.maximum(rho,  1e-40))
    logt   = np.log10(np.maximum(temp, 1.0))
    pts = np.column_stack([logt, logrho])
    gamma_arr = _ginterp(pts)
    cs = np.sqrt(gamma_arr * pgas / rho)

    mach      = v / np.maximum(cs, 1.0)
    prad_pgas = prad / np.maximum(pgas, 1e-30)

    fitting_model = fittingModel()
    vesc = sqrt(2 * fitting_model.ms * G / fitting_model.r)

    ekin_dens = 0.5 * rho * v**2
    etot_dens = ekin_dens + eg + erad          # rad + gas internal + kinetic (no gravity)

    # Total energy outflow rate [erg/s]: mass-loss rate x specific total energy.
    # etot_dens/rho is the specific total energy [erg/g]; mdot is dM/dt [g/s].
    edot_tot = mdot * (etot_dens / rho)

    # Cumulative energy outflow [erg]: integrate the rate over time.
    # t is in days -> convert to seconds for the integral.
    t_sec = t * 86400.0
    cumE = np.concatenate([[0.0], np.cumsum(0.5 * (edot_tot[1:] + edot_tot[:-1])
                                            * np.diff(t_sec))])

    return dict(mach=mach, erad=erad, eg=eg, prad_pgas=prad_pgas,
                cumE=cumE, vesc=vesc)


def draw_bcinput_axes(axes, t, v, rho, temp, mdot, culmdot, derived,
                      color='k', label=None, guides=True):
    """Draw one model's boundary-condition curves onto an existing 3x3 axes grid.

    Layout:
      row 0:  rho          |  v_ej            |  T
      row 1:  Mach (v/c_s) |  erad/eg         |  prad/pgas
      row 2:  Edot_tot     |  mass-loss dM/dt |  cumulative dM

    `derived` is the dict returned by compute_bc_derived. `guides` draws the
    reference lines (v_esc, Mach=1, prad=pgas); set False to avoid redrawing
    them when overlaying additional models.
    """
    ls = '-'
    axes[0,0].plot(t, rho, color=color, linestyle=ls, label=label)
    axes[0,0].set_yscale('log')
    axes[0,0].set_ylabel(r'$\rho\ [\rm{g\,cm^{-3}}]$')

    axes[0,1].plot(t, v/1e5, color=color, linestyle=ls)
    if guides:
        axes[0,1].hlines(derived['vesc']/1e5, xmin=0, xmax=140, color='black',
                         linestyle='solid', linewidth=2, label=r'$v_{\rm{esc}}$')
        axes[0,1].legend(fontsize=18)
    axes[0,1].set_ylim(30, 100)
    axes[0,1].set_ylabel(r'$v_{\rm{ej}}$' + ' [km/s]')

    axes[0,2].plot(t, temp/1e4, color=color, linestyle=ls)
    axes[0,2].set_ylim(2, 10)
    axes[0,2].set_ylabel(r'$T\ [\times10^{4}K]$')

    axes[1,0].plot(t, derived['mach'], color=color, linestyle=ls)
    if guides:
        axes[1,0].axhline(1.0, color='black', linestyle='--', linewidth=1.5,
                          label=r'$\mathcal{M}=1$')
        axes[1,0].legend(fontsize=18)
    axes[1,0].set_ylim(0, max(derived['mach'].max() * 1.2, 2.0))
    axes[1,0].set_ylabel(r'$\mathcal{M} = v_{\rm{ej}}/c_s$')

    axes[1,1].plot(t, derived['erad'] / np.maximum(derived['eg'], 1e-30),
                   color=color, linestyle=ls)
    axes[1,1].set_ylim(0.01, 5.1)
    axes[1,1].set_ylabel(r'$E_{\rm{r}}/e_{\rm{g}}$')

    axes[1,2].plot(t, derived['prad_pgas'], color=color, linestyle=ls)
    if guides:
        axes[1,2].axhline(1.0, color='black', linestyle='--', linewidth=1.5,
                          label=r'$p_{\rm{r}}=p_{\rm{gas}}$')
        axes[1,2].legend(fontsize=18)
    axes[1,2].set_yscale('log')
    axes[1,2].set_ylabel(r'$p_{\rm{r}}/p_{\rm{g}}$')

    axes[2,0].plot(t, derived['cumE'], color=color, linestyle=ls)
    axes[2,0].set_ylabel(r'$E_{\rm{tot}}\ [\rm{erg}]$')
    axes[2,0].set_xlabel('t [day]')

    axes[2,1].plot(t, mdot/(msun/year), color=color, linestyle=ls)
    axes[2,1].set_ylabel(r'$\dot{M}\ [M_{\odot}\cdot\rm{yr}^{-1}]$')
    axes[2,1].set_xlabel('t [day]')

    axes[2,2].plot(t, culmdot, color=color, linestyle=ls)
    axes[2,2].set_ylabel(r'$\Delta M\ [M_{\odot}]$')
    axes[2,2].set_xlabel('t [day]')


def plot_bcinput(t, v, rho, temp, mdot, culmdot, modeldir):
    """Plot bcinput.dat diagnostics in a 3x3 grid.

    Layout:
      row 0:  rho          |  v_ej          |  T
      row 1:  Mach (v/c_s) |  erad/eg       |  prad/pgas
      row 2:  Edot_tot     |  mass-loss dM/dt |  cumulative dM
    """
    derived = compute_bc_derived(t, v, rho, temp, mdot)

    fig, axes = plt.subplots(3, 3, figsize=(20, 11), sharex=True, squeeze=True)
    draw_bcinput_axes(axes, t, v, rho, temp, mdot, culmdot, derived, color='k')
    axes[2,0].set_xlim(0, 39)

    plt.tight_layout()
    plt.subplots_adjust(hspace=.0)
    plt.savefig(f'{modeldir}/boundcond.png', bbox_inches='tight', dpi=300)
    #plt.show()

    print(f"Generated boundary condition plot: {modeldir}/boundcond.png")

def save_problem(model, modeldir):
    probnml = f90nml.read('problem.data')
    probnml['parameters_1d']['m_star'] = cfg.ms  # Use cfg.ms
    probnml['parameters_1d']['record_length'] = 11
    probnml['parameters_1d']['petsc_iter'] = 6
    probnml['parameters_1d']['petsc_qratio'] = 1.4
    probnml['parameters_1d']['floor_tauR'] = 1
    probnml['parameters_1d']['rho_floor'] = cfg.rho_floor
    probnml['parameters_1d']['post_ej_bound'] = model.bound
    
    glbnml = f90nml.read('global.data')
    glbnml['meshinfo']['n_domain'] = [cfg.rin, cfg.rout, 0, 0]  # Use cfg.rin and cfg.rout
    glbnml['meshinfo']['lengthscale'] = 6.96e10
    glbnml['meshinfo']['tfinal'] = cfg.tfinal
    glbnml['meshinfo']['timescale'] = 86400
    glbnml['meshinfo']['CFL'] = 0.5 # 0.7
    glbnml['meshinfo']['nframe'] = cfg.nframe
    glbnml['meshinfo']['nrefine_region'] = 1
    glbnml['phyinfo']['llnx'] = True
    glbnml['phyinfo']['xgeo_h'] = cfg.xratio
    glbnml['phyinfo']['llny'] = False
    
    if cfg.iradiation == 0:
        probnml['parameters_1d']['lfld_mom'] = False
        probnml['parameters_1d']['larad'] = False
        probnml['parameters_1d']['lpradgradv'] = False
        glbnml['phyinfo']['lrad_adv'] = False
    elif cfg.iradiation == 4:
        probnml['parameters_1d']['lfld_mom'] = True
        probnml['parameters_1d']['larad'] = True
        probnml['parameters_1d']['lpradgradv'] = True
        glbnml['phyinfo']['lrad_adv'] = True
    
    glbnml['phyinfo']['lam_con'] = False
    glbnml['global_parameters']['nd'] = 1
    glbnml['global_parameters']['nx'] = cfg.nx
    glbnml['global_parameters']['ny'] = 1
    glbnml['global_parameters']['blk_size_nx'] = 64
    glbnml['global_parameters']['blk_size_ny'] = 1
    glbnml['global_parameters']['gamma_gas'] = cfg.gamma_eos
    glbnml['global_parameters']['iradiation'] = cfg.iradiation
    
    glbnml['refinement'][0]['refine_xmin'] = 0
    glbnml['refinement'][0]['refine_xmax'] = cfg.rin * 1.1  
    glbnml['refinement'][0]['refine_ymin'] = 0
    glbnml['refinement'][0]['refine_ymax'] = 0
    glbnml['refinement'][0]['level'] = cfg.level
    
    probnml.write(modeldir+'/problem.data', force=True)
    glbnml.write(modeldir+'/global.data', force=True)
    print(f"Successfully generated problem.data and global.data files")


def next_model_directory(width=5):
    """Return the next free model folder name, e.g. 'model00010'.

    Scans the current directory for existing model<number> folders and picks
    max+1, zero-padded to `width` digits (00000-99999). The copied executable
    is given the same name (see save_model), so folder and executable always
    match.
    """
    max_num = 0
    for name in os.listdir('.'):
        if os.path.isdir(name) and name.startswith('model'):
            digits = name[len('model'):]
            if digits.isdigit():
                max_num = max(max_num, int(digits))
    return 'model' + str(max_num + 1).zfill(width)


# Main Program
def main():
    model_template = model_class(ms=cfg.ms, rin=cfg.rin, time=cfg.ejt, n=cfg.n,
                               bound=cfg.bound, model_eos=cfg.model_eos)
    init_model_csv(header=model_template.header)

    model_directory = next_model_directory()

    model = model_class(ms=cfg.ms, rin=cfg.rin, time=cfg.ejt, n=cfg.n,
                       bound=cfg.bound, model_eos=cfg.model_eos)

    t, v, rho, temp, eratio, mdot, culmdot = save_model(model, modeldir=model_directory)
    save_problem(model, modeldir=model_directory)

    plot_bcinput(t, v, rho, temp, mdot, culmdot, model_directory)

    print(f"Successfully created model: {model_directory}")

if __name__ == "__main__":
    main()

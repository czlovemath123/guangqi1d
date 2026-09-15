"""abao.py -- low/high-speed lrne example driver for guangqi1d.

Lives in modules/lrne and manages two complete run directories,
lowspeed/ and highspeed/ (created next to this script), from setup to
post-processing:

    init     create the two folders with the low/high-speed bcinput and the
             model03000 run configuration
    run      build guangqi for the lrne module (if needed) and execute it in
             both folders
    bc       overlay the two boundary conditions on the 3x3 diagnostic grid
             (best_boundcond.png, like guangqi2d's plot_best.py)
    lc       light curves of both cases vs the AT2025abao observations
             (best_lc_comparison.png, like guangqi2d's plot_best_lc.py)
    rhol     per-folder rhoL.png (density / temperature / luminosity maps,
             like guangqi2d's kippenhahn2.single_model_evo2) into pictures/
    vhist    mass-weighted velocity histograms at t = 0/50/100 d
             (best_vhist.png, like guangqi2d's plot_best_vhist.py)
    escmass  escaping-mass evolution for both cases (escape_mass.png, like
             guangqi2d's escape_mass.py)
    all      init -> run -> bc -> lc -> rhol -> vhist -> escmass

Typical session:
    <venv>/bin/python abao.py init
    <venv>/bin/python abao.py run
    <venv>/bin/python abao.py bc
    <venv>/bin/python abao.py lc
    <venv>/bin/python abao.py rhol
    <venv>/bin/python abao.py vhist
    <venv>/bin/python abao.py escmass
or simply:
    <venv>/bin/python abao.py all
(<venv> is any Python 3 environment with numpy and matplotlib.)
"""
import argparse
import importlib.util
import types
import os
import shutil
import string
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, '..', '..'))       # this repository
SCRIPTS = os.path.join(REPO, 'scripts')
MODULE_LRNE = os.path.join(REPO, 'modules', 'lrne')
DEFAULT_BINARY = os.path.join(REPO, 'guangqi')

# phy_const, assemble_1d_data and eos are imported by the helper modules
# (model_generater.py, kippenhahn2.py); make them importable regardless of
# how/where this script is invoked.
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

# (folder, bcinput in modules/lrne, reference model dir under --src, label, color,
#  hardcoded light-curve phase shift [days])
CASES = [
    dict(folder='lowspeed',  bcinput='bcinput_lowspeed.dat',  ref='m8sub',
         label='low-speed',  color='#003f5c', shift=23.0),
    dict(folder='highspeed', bcinput='bcinput_highspeed.dat', ref='m8sup',
         label='high-speed', color='#d62728', shift=23.6),
]

# Run configuration, mirroring the archived m8sub/m8sup model03000 campaign
# (identical for both cases; only bcinput.dat differs). path_tables points at
# this repository's own opacity tables.
DEFAULT_TABLES = os.path.join(REPO, 'tables', 'opacity')
GLOBAL_DATA = {
    'meshinfo': {
        'n_domain': [500, 10000, 0, 0],
        'lengthscale': 6.96e10,
        'tfinal': 160,
        'timescale': 86400,
        'cfl': 0.5,
        'v_boundary': [0, 0, 0, 0],
        'hydro_bound_type': [9, 4, 2, 2],
        'rad_bound_type': [8, 1, 2, 2],
        'passive_bound_type': [9, 4, 2, 2],
        'nframe': 100,
        'refine_type': 'static',
        'nrefine_region': 1,
        'nderefine_region': 0,
        'max_refine_level': 10,
        'restart': False,
        'restart_iframe': 60,
    },
    'phyinfo': {
        'h_ratio': 0.74,
        'he_ratio': 0.26,
        'petsc_rtol': 1e-10,
        'llnx': True,
        'xgeo_h': 80,
        'llny': False,
        'ygeo_h': 0.5,
        'lrad_adv': True,
        'lam_con': False,
    },
    'global_parameters': {
        'igravity': 1,
        'igeometry': 2,
        'icooling': 0,
        'iradiation': 4,
        'nd': 1,
        'nx': 512,
        'ny': 1,
        'blk_size_nx': 64,
        'blk_size_ny': 1,
        'maw': 1.0,
        'gamma_gas': 1.4,
        'path_tables': DEFAULT_TABLES,
    },
    'refinement': {
        'refine_xmin': 0,
        'refine_xmax': 550.0,
        'refine_ymin': 0,
        'refine_ymax': 0,
        'level': 2,
    },
}
PROBLEM_DATA = {
    'parameters_1d': {
        'm_star': 8,
        'lfld_mom': True,
        'larad': True,
        'dt_record': 1000.0,
        'record_length': 11,
        'petsc_iter': 6,
        'petsc_qratio': 1.4,
        'floor_taur': 1,
        'inner_floor_temp': 300.0,
        'rho_floor': 1e-17,
        'post_ej_bound': 'reflective',
        't_ram': 0.5,
        'lpradgradv': True,
        'rho0': 1e-17,
        'temp0': 5000.0,
        'ejection_rad_bc_type': 8,
    },
}
OUTPUT_VAR_INFO = '5\nrho\nvx\npres\ntemp\negv\nlrne\n'


# --------------------------------------------------------------------- init
def _nml_value(v):
    if isinstance(v, bool):
        return '.true.' if v else '.false.'
    if isinstance(v, str):
        return f"'{v}'"
    if isinstance(v, (list, tuple)):
        return ', '.join(_nml_value(x) for x in v)
    return repr(v)


def write_namelist(path, groups):
    """Write groups in the given order -- the engine reads the namelists
    sequentially (meshinfo -> phyinfo -> global_parameters), so ordering
    matters. `groups` is a list of (group name, {key: value})."""
    lines = []
    for name, kv in groups:
        lines.append(f'&{name}')
        for k, v in kv.items():
            lines.append(f'    {k} = {_nml_value(v)}')
        lines.append('/')
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


def cmd_init(args):
    src = args.src
    for case in CASES:
        d = os.path.join(HERE, case['folder'])
        os.makedirs(os.path.join(d, 'out'), exist_ok=True)
        os.makedirs(os.path.join(d, 'pictures'), exist_ok=True)

        shutil.copy(os.path.join(MODULE_LRNE, case['bcinput']),
                    os.path.join(d, 'bcinput.dat'))

        # Configs: verbatim copy from the reference model dir if available,
        # otherwise generate from the baked model03000 values.
        ref_cfg = os.path.join(src, case['ref'], 'model03000')
        if os.path.isfile(os.path.join(ref_cfg, 'global.data')):
            for f in ('global.data', 'problem.data', 'output_var_info.dat'):
                shutil.copy(os.path.join(ref_cfg, f), os.path.join(d, f))
            print(f"- {case['folder']}: configs copied from {ref_cfg}")
        else:
            gd = dict(GLOBAL_DATA)
            gd['global_parameters']['path_tables'] = args.tables
            write_namelist(os.path.join(d, 'global.data'),
                           [('meshinfo', gd['meshinfo']),
                            ('phyinfo', gd['phyinfo']),
                            ('global_parameters', gd['global_parameters']),
                            ('refinement', gd['refinement'])])
            write_namelist(os.path.join(d, 'problem.data'),
                           [('parameters_1d', PROBLEM_DATA['parameters_1d'])])
            with open(os.path.join(d, 'output_var_info.dat'), 'w') as fh:
                fh.write(OUTPUT_VAR_INFO)
            print(f"- {case['folder']}: configs generated (model03000 values)")

        print(f"- {case['folder']}: bcinput <- {case['bcinput']}")
    print("init done (use 'abao.py run' to execute both cases)")


# ---------------------------------------------------------------------- run
def ensure_lrne_build(force=False, binary=DEFAULT_BINARY):
    """Point modules/problem at lrne and (re)build the repo binary if needed."""
    link = os.path.join(REPO, 'modules', 'problem')
    target = os.path.realpath(link) if os.path.islink(link) else None
    need = force or target != os.path.join(REPO, 'modules', 'lrne') \
        or not os.path.isfile(binary)
    if not need:
        return
    if os.path.islink(link) or os.path.exists(link):
        os.remove(link)
    os.symlink('lrne', link)
    print('building guangqi for the lrne module ...')
    subprocess.run(['make', 'clean'], cwd=REPO, check=True,
                   stdout=subprocess.DEVNULL)
    subprocess.run(['make'], cwd=REPO, check=True, stdout=subprocess.DEVNULL)
    print('build ok')


def cmd_run(args):
    ensure_lrne_build(force=args.rebuild, binary=args.binary)
    env = dict(os.environ)
    # Shared-library roots; override via the environment just like the Makefile
    # (PETSC_DIR, HDF5_DIR, OPENMPI_DIR). Non-existent roots are skipped.
    lib_roots = (os.environ.get('HDF5_DIR', '/usr/local/hdf5'),
                 os.environ.get('PETSC_DIR', '/usr/local/petsc'),
                 os.environ.get('OPENMPI_DIR', '/usr/local/openmpi'))
    extra_lib = [os.path.join(p, 'lib') for p in lib_roots
                 if os.path.isdir(os.path.join(p, 'lib'))]
    env['LD_LIBRARY_PATH'] = os.pathsep.join(
        extra_lib + [env.get('LD_LIBRARY_PATH', '')]).rstrip(os.pathsep)
    for case in CASES:
        d = os.path.join(HERE, case['folder'])
        if not os.path.isfile(os.path.join(d, 'bcinput.dat')):
            raise SystemExit(f"{case['folder']}: not initialized -- run 'abao.py init' first")
        print(f"=== running {case['folder']} ({case['label']}) ===")
        with open(os.path.join(d, 'run.log'), 'w') as log:
            rc = subprocess.run(['mpiexec', '-np', str(args.np), args.binary],
                                cwd=d, env=env, stdin=subprocess.DEVNULL,
                                stdout=log, stderr=subprocess.STDOUT).returncode
        last = os.path.join(d, 'last.dat')
        if rc != 0 or not os.path.isfile(last):
            raise SystemExit(f"{case['folder']}: run failed (rc={rc}), see run.log")
        print(f"- {case['folder']}: finished, last frame = "
              f"{open(last).read().strip()} (run.log in the folder)")
    print('run done (use abao.py bc / rhol / escmass for the figures)')


# ----------------------------------------------------------------------- bc
def load_generator():
    """Load model_generater.py (for compute_bc_derived / draw_bcinput_axes)."""
    spec = importlib.util.spec_from_file_location(
        'model_generater', os.path.join(HERE, 'model_generater.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_bcinput(path):
    """Read a bcinput table -> (t [day], v, rho, temp [K]); time is seconds."""
    raw = np.loadtxt(path, skiprows=1)
    return raw[:, 0] / 86400.0, raw[:, 1], raw[:, 2], raw[:, 3]


def cmd_bc(args):
    """Overlay the two boundary conditions on the 3x3 diagnostic grid
    (same figure as guangqi2d's plot_best.py)."""
    gen = load_generator()
    fm = gen.fittingModel()

    FSCALE = 1.3
    plt.rcParams.update({
        k: (v * FSCALE if isinstance(v, (int, float)) else v)
        for k, v in plt.rcParams.items()
        if k in ('font.size', 'axes.labelsize', 'xtick.labelsize',
                 'ytick.labelsize', 'legend.fontsize')
    })

    fig, axes = plt.subplots(3, 3, figsize=(20, 11), sharex=True, squeeze=True)

    print("=== Overlaying low/high-speed bcinput boundary conditions ===")
    drawn = 0
    for case in CASES:
        path = os.path.join(MODULE_LRNE, case['bcinput'])
        if not os.path.isfile(path):
            print(f"  ✗ {path} not found, skipping")
            continue
        print(f"- {case['bcinput']}: {case['label']}")
        t, v, rho, temp = read_bcinput(path)
        mdot = 4.0 * np.pi * fm.r**2 * rho * v
        t_sec = t * 86400.0
        culmdot = np.concatenate(
            [[0.0], np.cumsum(0.5 * (mdot[1:] + mdot[:-1]) * np.diff(t_sec))]) / gen.msun
        derived = gen.compute_bc_derived(t, v, rho, temp, mdot)
        # Draw reference guide lines only once (with the first case).
        gen.draw_bcinput_axes(axes, t, v, rho, temp, mdot, culmdot, derived,
                              color=case['color'], label=case['label'],
                              guides=(drawn == 0))
        drawn += 1

    if drawn == 0:
        raise SystemExit("No bcinput files drawn.")

    # Scale legends drawn by draw_bcinput_axes (fontsize=18) by 1.3x.
    for ax in axes.flat:
        lg = ax.get_legend()
        if lg is not None:
            for txt in lg.get_texts():
                txt.set_fontsize(18 * FSCALE)

    axes[0, 0].legend(fontsize=16 * FSCALE, loc='lower center')
    axes[2, 0].set_xlim(0, args.xmax)

    # Make the escape-velocity guide a dashed line (same width as Mach=1).
    for coll in axes[0, 1].collections:
        if coll.get_label() == r'$v_{\rm{esc}}$':
            coll.set_linestyles('--')
            coll.set_linewidth([1.5])

    # Label panels a..i (row 1 = a,b,c; row 2 = d,e,f; row 3 = g,h,i).
    for ax, tag in zip(axes.flat, string.ascii_lowercase[:9]):
        ax.text(0.015, 0.88, tag, transform=ax.transAxes,
                fontsize=30, fontweight='bold', va='top', ha='left')

    # Shared x-axis label.
    for ax in axes[2]:
        ax.set_xlabel('Time [day]')

    # Thicken every plotted line by 2x.
    for ax in axes.flat:
        for line in ax.get_lines():
            line.set_linewidth(line.get_linewidth() * 2)
        for coll in ax.collections:
            if coll.get_label() == r'$v_{\rm{esc}}$':
                coll.set_linewidth(np.asarray(coll.get_linewidth()) * 2)

    # E_tot panel: fold the 1e47 scale into the axis, ticks in 1e47 units.
    ETOT_SCALE = 1e47
    axes[2, 0].yaxis.set_major_formatter(
        FuncFormatter(lambda val, _: f'{val/ETOT_SCALE:g}'))
    axes[2, 0].set_ylabel(r'$E_{\rm{tot}}\ [10^{47}\,\rm{erg}]$')
    axes[2, 0].axhline(1.25 * ETOT_SCALE, color='black', linestyle='--',
                       linewidth=3.0, label=r'$1.25\times10^{47}\,\rm{erg}$')
    axes[2, 0].legend(loc='lower right', fontsize=16 * FSCALE)

    # Delta M panel: 0.92 Msun reference.
    axes[2, 2].axhline(0.92, color='black', linestyle='--',
                       linewidth=3.0, label=r'$0.92\,M_{\odot}$')
    axes[2, 2].legend(loc='lower right', fontsize=16 * FSCALE)

    # Density panel: linear axis in units of 1e-9 (yticks 2, 4, 10).
    RHO_SCALE = 1e-9
    for line in axes[0, 0].get_lines():
        if line.get_label() in [c['label'] for c in CASES]:
            line.set_ydata(line.get_ydata() / RHO_SCALE)
    axes[0, 0].set_yscale('linear')
    axes[0, 0].yaxis.set_major_locator(FixedLocator([2, 4, 10]))
    axes[0, 0].yaxis.set_major_formatter(FuncFormatter(lambda val, _: f'{val:g}'))
    axes[0, 0].set_ylabel(r'$\rho\ [10^{-9}\,\rm{g\,cm^{-3}}]$')
    axes[0, 0].relim()
    axes[0, 0].autoscale_view()

    plt.tight_layout()
    plt.subplots_adjust(hspace=.0, wspace=.2)
    out = os.path.join(HERE, args.out)
    plt.savefig(f'{out}.png', bbox_inches='tight', dpi=300)
    print(f"\n✓ saved '{out}.png'")


# --------------------------------------------------------------------- rhol
def cmd_rhol(args):
    # kippenhahn2.py ends with a bare "if len(sys.argv) > 1: plot_evolution(...)"
    # block (no __main__ guard), so strip it before executing the module.
    src_txt = open(os.path.join(HERE, 'kippenhahn2.py')).read()
    src_txt = src_txt[:src_txt.index('if len(sys.argv) > 1:')]
    k2 = types.ModuleType('kippenhahn2')
    k2.__file__ = os.path.join(HERE, 'kippenhahn2.py')
    exec(compile(src_txt, k2.__file__, 'exec'), k2.__dict__)
    cwd = os.getcwd()
    try:
        for case in CASES:
            d = os.path.join(HERE, case['folder'])
            if not os.path.isdir(os.path.join(d, 'out')):
                print(f"  ✗ {case['folder']}/out not found, skipping (run first)")
                continue
            os.makedirs(os.path.join(d, 'pictures'), exist_ok=True)
            print(f"- {case['folder']}: rhoL.png")
            k2.single_model_evo2(d, figname='rhoL.png')
    finally:
        os.chdir(cwd)
    print('rhol done (pictures/rhoL.png in each folder)')


# ------------------------------------------------------------------ escmass
MSUN = 1.99e33
G = 6.67259e-8
M_CENTRAL = 8.5 * MSUN


def frame_list(out_dir, filehead='lrne'):
    idx, k, misses = [], 0, 0
    while misses < 2:
        if os.path.isfile(os.path.join(out_dir, f'{filehead}{k:05d}.h5')):
            idx.append(k)
            misses = 0
        else:
            misses += 1
        k += 1
    return idx


def escaping_mass(ad, fn):
    xi = ad.glb_x_inter(fn)
    rho = ad.glb_cell_var(fn, 'rho')
    vx = ad.glb_cell_var(fn, 'vx')
    rc = 0.5 * (xi[:-1] + xi[1:])
    dr = np.diff(xi)
    dM = rho * 4.0 * np.pi * rc**2 * dr / MSUN
    esc = (vx > 0.0) & (0.5 * vx**2 > G * M_CENTRAL / rc)
    return float(np.sum(dM[esc])), float(np.sum(dM))


def cmd_escmass(args):
    import assemble_1d_data as ad
    plt.rcParams.update({'font.size': 20, 'font.family': 'serif'})

    fig, ax = plt.subplots(1, 1, figsize=(5.5, 4.5), squeeze=True)
    print("=== Escaping mass (0.5 v^2 > G*8.5Msun/r, vx>0) ===")
    for case in CASES:
        out_dir = os.path.join(HERE, case['folder'], 'out')
        if not os.path.isdir(out_dir):
            print(f"  ✗ {case['folder']}/out not found, skipping (run first)")
            continue
        shift = case['shift']
        days, mesc, mtot = [], [], []
        for k in frame_list(out_dir):
            fn = os.path.join(out_dir, f'lrne{k:05d}.h5')
            _, t, _, _ = ad.read_attr(fn)
            me, mt = escaping_mass(ad, fn)
            days.append(t / 86400.0 - shift)
            mesc.append(me)
            mtot.append(mt)
        days, mesc, mtot = map(np.asarray, (days, mesc, mtot))
        if len(days) == 0:
            print(f"  ✗ no frames for {case['folder']}")
            continue
        y = mesc / np.maximum(mtot, 1e-30) if args.frac else mesc
        ax.plot(days, y, color=case['color'], linewidth=3.0, marker='o',
                markersize=4, alpha=0.95, label=case['label'])
        print(f"- {case['folder']}: final escaping mass = {mesc[-1]:.4f} Msun "
              f"({100*mesc[-1]/max(mtot[-1],1e-30):.1f}% of {mtot[-1]:.3f} Msun "
              f"at {days[-1]:.1f} d)")

    ax.set_xlabel('Time [day]', fontsize=18)
    ax.set_ylabel('Escaping mass fraction' if args.frac
                  else r'Escaping mass [$M_\odot$]', fontsize=18)
    ax.text(0.05, 0.92, 'd', fontweight='bold', transform=ax.transAxes,
            ha='center', color='black')
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.grid(True, alpha=0.3)
    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.legend(fontsize=18, loc='best')
    plt.tight_layout()
    out = os.path.join(HERE, args.out)
    plt.savefig(f'{out}.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ saved '{out}.png'")


# ----------------------------------------------------------------------- lc
# history.data columns (0-based), same convention as plot_best_lc.py/fit_score.
COL_TIME = 0
COL_LUM = 3
OBS_FILE = os.path.join(HERE, 'post_UVOIR_RTL_fin.txt')


def load_history(history_path):
    """Return (t_days, L) sorted by time from a history.data file
    (column 0 = time [s], column 3 = luminosity [erg/s])."""
    d = np.loadtxt(history_path, max_rows=200000)
    t = d[:, COL_TIME] / 86400.0
    L = d[:, COL_LUM]
    order = np.argsort(t)
    return t[order], L[order]


def read_obs(path=OBS_FILE):
    """Return (phase, L, sigL) -- columns 0, 5, 6 of the observation file
    (same convention as fit_score.read_obs)."""
    d = np.loadtxt(path, max_rows=1000)
    return d[:, 0], d[:, 5], d[:, 6]


def cmd_lc(args):
    plt.rcParams.update({'font.size': 20, 'font.family': 'serif'})
    fig, ax = plt.subplots(1, 1, figsize=(6.25, 6), squeeze=True)

    obs = None
    if os.path.isfile(OBS_FILE):
        obs = read_obs()
    elif not args.no_obs:
        raise SystemExit(f"Observation file '{OBS_FILE}' not found; "
                         f"use --no-obs to plot the models only.")

    print("=== Plotting low/high-speed light curves ===")
    plotted = 0
    for case in CASES:
        hp = os.path.join(HERE, case['folder'], 'history.data')
        if not os.path.isfile(hp):
            print(f"  ✗ {case['folder']}/history.data not found, skipping (run first)")
            continue
        shift = args.shift if args.shift is not None else case['shift']
        print(f"- {case['folder']} ({case['label']}, shift {shift:.1f} d)")
        t, L = load_history(hp)
        t = t - shift
        mask = (t >= args.xmin) & (t <= args.xmax)
        if not np.any(mask):
            print(f"  ✗ no data in range for {case['folder']}")
            continue
        ax.plot(t[mask], L[mask] / 1e39, color=case['color'], linewidth=3.2,
                alpha=0.95, zorder=10)
        print(f"  ✓ plotted {len(t[mask])} points")
        plotted += 1

    if plotted == 0:
        raise SystemExit("No model curves were drawn.")

    if obs is not None and not args.no_obs:
        phase, L_obs, sigL = obs
        ax.errorbar(phase, L_obs / 1e39, yerr=sigL / 1e39, fmt='o',
                    color='darkorange', ecolor='black', elinewidth=1.5,
                    capsize=4, capthick=1.5, markersize=8,
                    markerfacecolor='lightyellow', markeredgecolor='black',
                    markeredgewidth=1.5, linewidth=1.5, label='AT2025abao',
                    alpha=0.8, zorder=6)

    ax.set_xlim(args.xmin, args.xmax)
    ax.set_ylim(args.ymin, args.ymax)
    ax.legend(fontsize=18, loc='upper right')
    ax.set_xlabel('Time [day]', fontsize=18)
    ax.set_ylabel('Luminosity [$10^{39}$ erg/s]', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.grid(True, alpha=0.3, zorder=1)
    ax.yaxis.get_offset_text().set_fontsize(18)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, zorder=1)
    ax.text(0.05, 0.92, 'a', fontweight='bold', transform=ax.transAxes,
            ha='center', color='black')

    plt.tight_layout()
    out = os.path.join(HERE, args.out)
    plt.savefig(f'{out}.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ saved '{out}.png'")


# -------------------------------------------------------------------- vhist
V_MIN, V_MAX, V_STEP = -30.0, 190.0, 10.0
MASS_UNIT = 0.1   # histogram mass weights are in units of 0.1 Msun
TARGET_DAYS = [
    (0,   '#003f5c'),
    (50,  '#ffa600'),
    (100, '#d62728'),
]


def frame_times(ad, out_dir, filehead='lrne'):
    idx, times = [], []
    k, misses = 0, 0
    while misses < 2:
        fn = os.path.join(out_dir, f'{filehead}{k:05d}.h5')
        if os.path.isfile(fn):
            _, t, _, _ = ad.read_attr(fn)
            idx.append(k)
            times.append(t)
            misses = 0
        else:
            misses += 1
        k += 1
    return np.asarray(idx), np.asarray(times)


def shell_mass_and_v(ad, fn):
    xi = ad.glb_x_inter(fn)
    rho = ad.glb_cell_var(fn, 'rho')
    vx = ad.glb_cell_var(fn, 'vx')
    rc = 0.5 * (xi[:-1] + xi[1:])
    dr = np.diff(xi)
    dM = rho * 4.0 * np.pi * rc**2 * dr / MSUN
    return vx / 1e5, dM


def cmd_vhist(args):
    import assemble_1d_data as ad
    plt.rcParams.update({'font.size': 18, 'font.family': 'serif'})
    bins = np.arange(V_MIN, V_MAX + 0.5 * V_STEP, V_STEP)

    fig, axes = plt.subplots(1, len(CASES), figsize=(11, 4.5), sharey=True,
                             squeeze=False)
    axes = axes[0]

    for col, case in enumerate(CASES):
        ax = axes[col]
        out_dir = os.path.join(HERE, case['folder'], 'out')
        if not os.path.isdir(out_dir):
            print(f"  ✗ {case['folder']}/out not found, skipping (run first)")
            continue

        shift = case['shift']
        idx, times = frame_times(ad, out_dir)
        if len(idx) == 0:
            print(f"  ✗ no frames in {out_dir}, skipping")
            continue
        days = times / 86400.0 - shift

        print(f"- {case['folder']} (shift {shift:.1f} d, {len(idx)} frames, "
              f"phase {days.min():.1f}..{days.max():.1f} d)")

        for target, color in TARGET_DAYS:
            j = int(np.argmin(np.abs(days - target)))
            fn = os.path.join(out_dir, f'lrne{idx[j]:05d}.h5')
            v, dM = shell_mass_and_v(ad, fn)
            w = dM / MASS_UNIT if args.weight == 'mass' else None
            ax.hist(v, bins=bins, weights=w, histtype='step', linewidth=2.4,
                    color=color, alpha=0.9, label=fr'$t={target}$ d')
            print(f"    t={target:>3} d -> frame {idx[j]:>3} (actual {days[j]:.1f} d)")

        ax.text(0.3, 0.98, case['label'], transform=ax.transAxes,
                ha='center', va='top', fontsize=20, zorder=5)
        ax.set_xlabel('Velocity [km/s]', fontsize=20)
        ax.set_xlim(V_MIN, V_MAX)
        if args.logy:
            ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=14, loc='upper right')

    for ax, tag in zip(axes, ['e']):
        ax.text(0.02, 0.95, tag, transform=ax.transAxes,
                ha='left', va='top', fontsize=18, fontweight='bold', zorder=5)

    ylabel = (r'Mass [$0.1\,M_\odot$]' if args.weight == 'mass' else 'Cell count')
    axes[0].set_ylabel(ylabel, fontsize=20)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.0)
    out = os.path.join(HERE, args.out)
    plt.savefig(f'{out}.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ saved '{out}.png'")


# ---------------------------------------------------------------------- all
def cmd_all(args):
    """init -> run -> bc -> lc -> rhol -> vhist -> escmass."""
    cmd_init(argparse.Namespace(src=args.src, tables=args.tables))
    cmd_run(argparse.Namespace(binary=args.binary, np=args.np,
                               rebuild=args.rebuild))
    cmd_bc(argparse.Namespace(out='best_boundcond', xmax=39))
    cmd_lc(argparse.Namespace(shift=None, no_obs=False,
                              xmin=-18, xmax=120, ymin=0.5, ymax=5.6,
                              out='best_lc_comparison'))
    cmd_rhol(argparse.Namespace())
    cmd_vhist(argparse.Namespace(out='best_vhist', weight='mass', logy=False))
    cmd_escmass(argparse.Namespace(out='escape_mass', frac=False))
    print("\nall done.")


# --------------------------------------------------------------------- main
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    sp = sub.add_parser('init', help='create the lowspeed/ and highspeed/ folders')
    sp.add_argument('--src', default=HERE,
                    help="folder holding <m8sub,m8sup>/model03000 reference "
                         "configs (default: this folder; falls back to the "
                         "baked model03000 values)")
    sp.add_argument('--tables', default=DEFAULT_TABLES,
                    help="path_tables written into global.data when configs "
                         "are generated (default: the merged repo's tables)")
    sp.set_defaults(func=cmd_init)

    sp = sub.add_parser('run', help='execute guangqi in both folders')
    sp.add_argument('--binary', default=DEFAULT_BINARY)
    sp.add_argument('--np', type=int, default=1)
    sp.add_argument('--rebuild', action='store_true',
                    help='force make clean && make for the lrne module')
    sp.set_defaults(func=cmd_run)

    sp = sub.add_parser('bc', help='best_boundcond.png: both boundary conditions')
    sp.add_argument('--out', default='best_boundcond')
    sp.add_argument('--xmax', type=float, default=39)
    sp.set_defaults(func=cmd_bc)

    sp = sub.add_parser('rhol', help='rhoL.png per folder into pictures/')
    sp.set_defaults(func=cmd_rhol)

    sp = sub.add_parser('escmass', help='escape_mass.png for both cases')
    sp.add_argument('--out', default='escape_mass')
    sp.add_argument('--frac', action='store_true',
                    help='plot escaping-mass fraction instead of mass')
    sp.set_defaults(func=cmd_escmass)

    sp = sub.add_parser('lc', help='best_lc_comparison.png: light curves vs AT2025abao')
    sp.add_argument('--shift', type=float, default=None,
                    help='force this shift (days) for both models '
                         '(default: hardcoded 23.0 low-speed, 23.6 high-speed)')
    sp.add_argument('--no-obs', action='store_true', dest='no_obs',
                    help='do not draw the AT2025abao observations')
    sp.add_argument('--xmin', type=float, default=-18)
    sp.add_argument('--xmax', type=float, default=120)
    sp.add_argument('--ymin', type=float, default=0.5)
    sp.add_argument('--ymax', type=float, default=5.6)
    sp.add_argument('--out', default='best_lc_comparison')
    sp.set_defaults(func=cmd_lc)

    sp = sub.add_parser('vhist', help='best_vhist.png: velocity histograms')
    sp.add_argument('--out', default='best_vhist')
    sp.add_argument('--weight', choices=['mass', 'count'], default='mass',
                    help='histogram weight: shell mass (default) or cell count')
    sp.add_argument('--logy', action='store_true', help='use a log y-axis')
    sp.set_defaults(func=cmd_vhist)

    sp = sub.add_parser('all', help='init -> run -> bc -> lc -> rhol -> vhist -> escmass')
    sp.add_argument('--src', default=HERE,
                    help="folder holding <m8sub,m8sup>/model03000 reference "
                         "configs (default: this folder; falls back to the "
                         "baked model03000 values)")
    sp.add_argument('--tables', default=DEFAULT_TABLES,
                    help='path_tables written into generated global.data')
    sp.add_argument('--binary', default=DEFAULT_BINARY)
    sp.add_argument('--np', type=int, default=1)
    sp.add_argument('--rebuild', action='store_true',
                    help='force make clean && make for the lrne module')
    sp.set_defaults(func=cmd_all)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == '__main__':
    main()

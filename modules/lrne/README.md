# modules/lrne — luminous red nova (LRNe) ejecta module

This module simulates the ejecta of a common-envelope event, as observed in
luminous red novae (LRNe), with the 1D (spherically symmetric) version of
Guangqi: gas is injected through the inner boundary according to a
time-dependent ejecta table, expands through the envelope of the giant, and
the resulting bolometric light curve is measured at the outer boundary.

The module was used for the LRNe/CEE studies of the code and, in its
`abao.py` campaign, for the AT 2025abao analysis of Chen et al. 2026,
submitted to ApJL.

## The model (`problem.f90`)

**Geometry and grid.** Purely radial 1D polar geometry (`igeometry=2`,
`ny=1`): the domain runs from an inner radius `n_domain(1)` to an outer
radius `n_domain(2)` in units of `lengthscale` (Rsun). The module runs on a
**single MPI rank only** (`mpirun -np 1`); static or SMR grids only
(adaptive refinement is rejected).

**Gravity.** Point-mass potential of a central star of mass `m_star`
(`problem.data`).

**Initial condition.** Static (`v=0`) envelope with

- `rho(r)   = rho0  * (rin/r)^1.5`
- `T(r)     = temp0 * (r/rin)^-1`
- radiation energy in equilibrium with the gas, plus density floors
  (`rho_floor`, an `1/r`-like floor) to keep the outer cells well-behaved.

**Inner boundary — the ejecta table.** The inner ghost cells are set by
linear interpolation in time of `bcinput.dat` (**required**; the run aborts
without it). The table has six columns:

| col | quantity |
|-----|----------|
| 1   | time [s] |
| 2   | ejecta velocity `vej(t)` |
| 3   | ejecta density `rhoej(t)` |
| 4   | ejecta temperature `tempej(t)` (mid-plane value) |
| 5   | fixed `Erad/Egas` ratio of the injected gas |
| 6   | asymmetry opening angle in units of pi |

The temperature of the injected gas is reconstructed from the
`Erad/Egas` ratio through the equation of state (bisection between the
molecular and fully ionized branches, 500 K – 10^6 K, for the realistic
H/He EOS `ieos=2`).

**Injection phases** (`time_dependent_bound_type`):

1. `t < t_ram` (days): injection is held off (`hydro_bound_type(1)=9`,
   radiation inner boundary type 8) so that overly hot ejecta does not
   enter the domain immediately.
2. `t_ram < t < tmax`, with `tmax` the last time in `bcinput.dat`: ejecta
   injection; the inner radiation boundary is `ejection_rad_bc_type`
   (see `problem.f90`).
3. `t > tmax`: ejecta phase finished; the inner hydrodynamic boundary
   becomes `post_ej_bound` (`reflective`, `transmissive`, or
   `constant_pres`) and the inner radiation boundary a zero-gradient-like
   specified condition.

The outer radiation boundary falls off as `Erad * (r1/r)^2`.

**Radiation transport.** Flux-limited diffusion (`iradiation=4`) solved
implicitly with PETSc (`petsc_qratio`, `petsc_iter`, `petsc_rtol`).
Opacities come from the tables in `path_tables` (Rosseland and Planck),
clamped by `floor_tauR` (Rosseland optical-depth floor) and Planck-opacity
ceilings, and are valid for `200 K <= T <= 3e6 K`,
`1e-18 <= rho <= 1e-4 g/cm^3` (`opacity_gas_*` in `problem.f90`).

At start-up the code prints the ejecta Mach number and the escape velocity
of the central star.

## Input decks

- `global.data` — global namelists (`meshinfo`, `phyinfo`,
  `global_parameters`, optional `refinement`); remember to point
  `path_tables` at this repository's `tables/opacity`.
- `problem.data` — module namelist `&parameters_1d`: `m_star`, `rho0`,
  `temp0`, `t_ram`, `post_ej_bound`, `ejection_rad_bc_type`, floors,
  PETSc controls, `dt_record` (history cadence), and related switches.
- `bcinput.dat` — the ejecta table above. Templates for the two AT 2025abao
  realizations are provided: `bcinput_lowspeed.dat` and
  `bcinput_highspeed.dat` (copy one to `bcinput.dat`).
- `output_var_info.dat` — names of the variables stored in the HDF5 frames;
  prepare from `output_var_info.dat.empty`.

## Outputs

- `history.data` — one row per `dt_record`: time [s] (column 0), inner and
  outer bolometric luminosities `4*pi*r^2*F_rad` (columns 1 and 3, erg/s),
  inner and outer mass fluxes (columns 2 and 5, Msun/yr) and radiation /
  total energy fluxes (see `assemble_record_array` in `problem.f90`).
- `out/lrneNNNNN.h5` (+ `.xdmf`) — snapshot frames for ParaView.
- `last.dat` — latest frame index (restart support via `restart=.true.` and
  `restart_iframe`).

## Running

```bash
ln -sfn lrne modules/problem      # from the repository root
make clean && make                # single-rank build
cd modules/lrne
cp bcinput_highspeed.dat bcinput.dat      # or the lowspeed table
mpiexec -np 1 ../../guangqi       # single rank is mandatory in this module
```

## `abao.py` — the AT 2025abao campaign

`abao.py` is a self-contained driver that reproduces the AT 2025abao
light-curve fits of Chen et al. 2026, submitted to ApJL. It manages two
complete run directories next to itself,

- `lowspeed/`  from `bcinput_lowspeed.dat` (label "low-speed"),
- `highspeed/` from `bcinput_highspeed.dat` (label "high-speed"),

both using an `8 Msun` central star and the `model03000` configuration
(domain 500–10000 Rsun, `tfinal=160` d, 512 cells, one level-2 refinement
region below 550 Rsun). `init` copies the reference configurations from
`<m8sub|m8sup>/model03000` if available and otherwise writes the baked-in
values; only `bcinput.dat` differs between the two cases.

Subcommands (in pipeline order; `all` runs everything):

| command | what it produces |
|---------|------------------|
| `init`  | create `lowspeed/` and `highspeed/` with `bcinput.dat`, `global.data`, `problem.data`, `output_var_info.dat`, `out/` |
| `run`   | point `modules/problem` at `lrne`, rebuild `guangqi` if needed, `mpiexec -np 1` in each folder (`run.log` kept per folder) |
| `bc`    | `best_boundcond.png` — both ejecta tables on a 3x3 diagnostic grid (v, rho, T, Mdot, cumulative ejected mass, E_tot with `1.25e47 erg` and `0.92 Msun` guides) |
| `lc`    | `best_lc_comparison.png` — model light curves (outer-boundary luminosity from `history.data`) phase-shifted by 23.0 d (low-speed) / 23.6 d (high-speed), overplotted on the AT 2025abao UVOIR observations in `post_UVOIR_RTL_fin.txt` |
| `rhol`  | `pictures/rhoL.png` per folder — density/temperature/luminosity Kippenhahn-style maps (via `kippenhahn2.py`) |
| `vhist` | `best_vhist.png` — mass-weighted velocity histograms at t = 0, 50, 100 d |
| `escmass` | `escape_mass.png` — escaping mass (`0.5 v^2 > G*M/r`, `M = 8.5 Msun`) and its evolution |

Typical session (any Python 3 environment with numpy and matplotlib):

```bash
python abao.py init
python abao.py run
python abao.py all        # or the individual steps: bc, lc, rhol, vhist, escmass
```

Useful options: `run --rebuild` (force `make clean && make`), `lc --shift`
(override the phase alignment), `lc --no-obs` (models only), `escmass
--frac` (escaping-mass fraction).

## Other scripts in this folder

- `model_generater.py` — ejecta-table design and derived-quantity helpers
  (used by `abao.py bc` for the diagnostic panels).
- `kippenhahn2.py` — per-model evolution figures (`single_model_evo2`).
- `gamma.txt`, `wfst.dat`, `post_UVOIR_RTL_fin.txt` — module data and the
  AT 2025abao observed UVOIR light curve used by `abao.py lc`.

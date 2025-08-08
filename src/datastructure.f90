module datastructure
#include <petsc/finclude/petscksp.h>
use mpi
use petscksp
implicit none

real(8), parameter, public ::    &
pi=3.141592653589793238462643383279502884197169399375105d0,  &
e_exp=2.718281828459d0

real(8), parameter, public ::    &
gr=6.67259d-8,              &   !gravitational constant
NA=6.02214076d-23,          &   !Avogadro Constant
kb=1.380658d-16,            &   !boltzmann constant
rsun=6.96d10,               &   !sun radius
rjupiter=7.14d9,            &   !jupiter radius
au=1.496d13,                &
km=1d5,                     &
year=3.15569d7,             &
day=86400d0,                &
t_hubble=4.55d17,           &   !hubble time
me=9.1093897d-28,           &   !electron mass
mp=1.6726231d-24,           &   !proton mass
c_light=2.9979d10,          &   !light speed
msun=1.99d33,               &   !sun mass
lsun=3.9d33,                &   !sun luminosity
mearth=5.976d27,            &   !earth mass
mjupiter=1.899d30,          &   !jupiter mass
ljupiter=3.463798d45,       &   !estimated jupiter angular momentum
h_planck=6.6260755d-27,     &   !planck constant
amu=1.660539040d-24,        &   !atomic mass unit
mh=1.6733d-24,              &   !hydrogen atom mass
mh2=3.3466d-24,             &   !hydrogen molecule mass
mhe=6.646481526d-24,        &   !helium mass
muh2=2.0158816d0,           &   !atomic weight of hydrogen molecule
muh=1.0076849d0,            &   !atomic weight of hydrogen
muhe=4.0026d0,              &   !atomic weight of helium
sigma_sb=5.67051d-5,        &   !Stefan-Boltzmann constant
ionh=2.18d-11,              &   !atomic hydrogen binding energy
dish=7.17d-12,              &   !molecular hydrogen dissociation energy
ionhe1=3.94d-11,            &   !helium first ionization
ionhe2=8.72d-11,            &   !helium second ionization
a_rad=7.5646d-15            !radiation energy density constant

real(8), parameter :: zero=0d0, fourth=0.25d0, half=0.5d0
real(8), parameter :: one=1d0, two=2d0, three=3d0, four=4d0
complex(8), parameter :: czero=(0.d0,0.d0)


interface initialize_global_parameters
    module procedure initialize_global_parameters_int,initialize_global_parameters_real,  &
    initialize_global_parameters_logical
end interface initialize_global_parameters

interface allocate_cell_data_block
    module procedure allocate_cell_data_block_1d,allocate_cell_data_block_nd,allocate_cell_data_block_1d_int,allocate_cell_data_block_nd_int
end interface allocate_cell_data_block

interface allocate_xsurface_data_block
    module procedure allocate_xsurface_data_block_1d,allocate_xsurface_data_block_nd
end interface allocate_xsurface_data_block

interface allocate_ysurface_data_block
    module procedure allocate_ysurface_data_block_1d,allocate_ysurface_data_block_nd
end interface allocate_ysurface_data_block

interface sync_bcast
    module procedure sync_bcast_logical1,sync_bcast_logical2,sync_bcast_int1,sync_bcast_int2
end interface sync_bcast

interface update_protected
    module procedure update_protected_logical,update_protected_real,update_protected_int
end interface update_protected

interface blk_traversal
    module procedure blk_traversal1,blk_traversal2
end interface blk_traversal

type axbsystem
    Vec x
    Vec b
    Mat A
    KSP ksp
end type axbsystem

type passive_obj
    integer :: passive_type     !1-dependent,2-independent,3-implicit
    character(len=16) :: passivename
    type(passive_obj), pointer :: next=>null()
    procedure(blk_operator), pointer, nopass :: sub=>null()
end type passive_obj

type processor_variables
    logical :: processor_logical
    logical :: processor_nan
    logical :: petsc_decrease_rtol
    logical, dimension(:), allocatable :: mpi_logical,global_logical
    integer :: amr_ngrow,smr_ngrow,mr_ngrow
    integer :: amr_nderefine,smr_nderefine,mr_nderefine,mr_ntrim
    integer :: n_fulledge,n_halfedge,n_corner
    integer :: hydro_edge_corner(5),fld_edge_corner(5)      !1-xfull edge, 2-yfull edge, 3-xhalf edge, 4-yhalf edge, 5-corner
    integer :: hydro_edge_corner_non_prolongation(5),hydro_edge_corner_prolongation(3)
    integer :: hydro_flux_edge(2)                           !only need when from high level to low level
    integer, dimension(:), allocatable :: mpi_processor_oper
    integer, dimension(:), allocatable :: processor_operation_sum
    integer, dimension(:,:), allocatable :: grow_keys
    !each processor stores the sequence of the keys to a this list
    integer, dimension(:,:), allocatable :: derefine_keys
    integer, dimension(:,:), allocatable :: mpi_world_keys
    !rank0 processor stack all the keys to this list, then bcast it to the other processors
    real(8) :: collective_processor_result
    real(8) :: collective_result
    real(8), dimension(:), allocatable :: collective_array
end type processor_variables

type timedef
    real(8) :: t                            !current time
    real(8) :: t_output_initial             !the time of frame 1
    real(8) :: t_final                      !end program when t=tfinal
    real(8) :: dt_hydro                     !hydrodynamical timescale
    real(8) :: dt_processor                 !timestep of this processor
    real(8) :: dt_source                    !minimum of dt_thermal,dt_gravity,dt_radiation, and dt_hydro
    real(8) :: dt_thermal                   !thermal timescale
    real(8) :: dt_gravity                   !gravitational timescale
    real(8) :: dt_radiation                 !radiative timescale
    real(8) :: dt_frame                     !time interval between each frames
    real(8) :: t_record
    real(8) :: dt_record                    !time interval to write a record
    real(8) :: dt_radhydro
    real(8) :: t_next                       !next time of the output frame
    integer :: ntimestep
end type timedef

type scaledef
    !the default scales are cgs unit
    real(8) :: lscale                       !length scale
    real(8) :: timescale                    !timescale
    real(8) :: rscale                       !density scale
    real(8) :: tempscale                    !temperature scale
    real(8) :: vscale                       !velocity scale
    real(8) :: mscale                       !mass scale
end type scaledef

integer, parameter :: irho=1,ivx=2,ivy=3,ivz=4,ipres=5,itemp=6,iegv=7
integer, parameter :: imass=1,imomx=2,imomy=3,imomz=4,iE=5
integer, parameter :: pyro80=1,pyro70=2,pyro60=3,pyro50=4
integer, parameter :: ku_hlle=3,kl_hlle=3,ktotal=10
integer, parameter :: ku_lapack=2,kl_lapack=2,ktotal_lapack=7
integer, parameter :: transmissive=1,reflective=2,extrapolate=3,specified=9
integer, protected :: ioformat=2
integer, protected :: igravity=0    !self gravity=1, uniform gravity=2
integer, protected :: igeometry=0   !cartesian=0, cylindrical=1, spherical=2
integer, protected :: icooling=0
integer, protected :: iradiation=0
integer, protected :: nd
integer, protected :: nx
integer, protected :: ny
integer, protected :: nz
integer, protected :: llevel_max            !maxiumum logical level
integer, protected :: nx2
integer, protected :: n_hydro_guard
integer, protected :: xlb                   !x coordinate, lower boundary
integer, protected :: xub                   !x coordinate, upper boundary
integer, protected :: ylb
integer, protected :: yub
integer, protected :: nrefine_region
integer, protected :: nderefine_region
integer, protected :: max_refine_level
integer, protected :: np                    !number of processors
integer, protected :: rank
integer, protected :: restart_iframe
integer, protected :: blk_hydro_size
integer, protected :: blk_cell_size
integer, protected :: blksize
integer, protected :: blk_interface_size
integer, protected :: blk_corner_size_hydro
integer, public :: cell_var_length
integer, public :: counter
integer, protected :: nray_static=32
logical, protected :: lrad_adv=.false.
logical, protected :: lpassive=.false.
logical :: lam_con=.false.
logical, protected :: lstiff_source
logical, protected :: restart
logical, public :: truncate_dt_hydro=.false.
logical, public :: lviscous=.false.
logical, public :: lirradiation=.false.
logical, public :: lfld_heating=.false.
logical, public :: llnx=.false.
logical, public :: llny=.false.
logical, public :: llntheta=.false.
logical, public :: wedge_polar=.false.
logical, public :: lfloor_temp=.false.
logical, public :: lceiling_vx=.false.
logical, public :: lceiling_vy=.false.
logical, public :: no_implicit=.true.
logical, public :: temp_smr=.false.
logical, public :: manual_smr=.false.
logical, public :: refine_zeldovich=.false.
logical, public :: refine_density=.false.
logical, public :: refine_erad=.false.
logical, public :: same_abs_ems=.true.
logical, public :: lfld_mom=.false.
logical, public :: larad=.false.
logical, public :: lsource=.false.
logical, public :: linternal_bound=.false.
real(8), protected :: n_domain(4)
real(8), protected :: dxyz(3)
real(8), protected :: theta_ratio0,dtheta_base0
real(8), protected :: maw
real(8), protected :: gamma_gas
real(8), protected :: CFL
real(8), protected :: temp_Erad_boundary=100d0
real(8), protected :: fld_diffusion_boundary_value=0d0      
real(8), protected :: pfloor
real(8), protected :: g_uniform
real(8), protected :: base_blk_size(2)
real(8), protected :: floor_temp=1d1
real(8), public :: ceiling_vx=1d6
real(8), public :: ceiling_vy=1d6
real(8), public :: const_opacity=1d0
real(8), public :: const_Erad=1d0
real(8), public :: weight_max
character(len=32), public :: refine_type
character(len=16), protected :: left,right,east,south,west,north,se,sw,nw,ne
character(len=16), protected :: hydro_field,fld_field,passive_field,viscous_field,erad_field
character(len=16), protected :: predictor,corrector,prolong,restrict,non
character(len=16), protected :: pn_vphi,pn_erad,pn_omega !independent passive scalars
character(len=16), protected :: pn_am,pn_ekphi   !dependent passive scalars
character(len=16), protected :: all_direction(8),corner_direction(4),edge_direction(4)
character(len=16), protected :: sdensity,stgas,strad,xcoord,ycoord,next,previous,sshear
character(len=255), protected :: path_root,path_tables
integer, protected :: blk_size_nx,blk_size_ny,blk_size_nz
integer, protected :: nx_blks,ny_blks,blk_xlb,blk_xub,blk_ylb,blk_yub
type(passive_obj), pointer :: passive_list=>null()
procedure(condition_hydro), pointer :: bound_hydro=>null(),init_hydro=>null(),bound_scalar=>null(),init_scalar=>null(), &
    bound_rad=>null(),init_rad=>null(),user_source=>null()

type blockpointer
    integer :: cpu_rank
    !this blk pointer determines the nature of the blockpointer
    type(blockdef), pointer :: blk
end type blockpointer

type blockdef
    integer :: blk_id,level,static_level
    character(len=1) :: curve
    integer :: key(3)                       !(i,j,logical_level)
    integer :: nb_level(8),nb_level_1d(2)
    integer :: refine_countdown,derefine_countdown
    integer :: idx,idy
    type(blockdef), pointer :: pre,next,pa,child_head,child_tail  !blk_xl,blk_xu,blk_xlyl,blk_xuyl,blk_xlyu,blk_xuyu
    type(blockpointer) :: children(2,2),nb(12)
    real(8) :: dxyz(3)                      !coordinate difference of a cell
    !the lower limits of all three boundaries and the numerical size in x and y coordinates
    real(8) :: pos(3)
    real(8) :: erad_vertex(4)
    real(8), allocatable :: erad_xedge(:,:),erad_yedge(:,:)
    real(8) :: Erad_pre,Erad_next,Erad_pre_int,Erad_next_int,sigma_rosseland_pre,sigma_rosseland_next,rho_pre,rho_next
    real(8) :: dx_xl,dx_xu,dx_yl,dx_yu,xl_guard,xu_guard
    real(8), dimension(:), allocatable :: dy,dr,dtheta,blk_var
    real(8), dimension(:,:), allocatable :: dy2d
    real(8) :: dt_hydro_blk,dtheta0,theta_ratio
    !Mesh points, for XDMF output
    real(8), dimension(:), allocatable :: mesh_x,mesh_y
    !Implicit boundary
    real(8), dimension(:,:), allocatable :: eradxm,eradxp,eradym,eradyp
    real(8), dimension(:,:), allocatable :: eradintxm,eradintxp,eradintym,eradintyp
    !Surface quantities, values are on the interfaces
    real(8), dimension(:,:,:,:), allocatable :: xflux,yflux
    real(8), dimension(:,:,:), allocatable :: hllc_vx,hllc_vy
    real(8), dimension(:,:,:,:), allocatable :: xpflux,ypflux
    real(8), dimension(:,:,:,:), allocatable :: gx,gy,alphax,alphay
    real(8), dimension(:,:), allocatable :: gx1d,alphax1d
    real(8), dimension(:,:,:), allocatable :: erad_xflux,erad_yflux
    real(8), dimension(:,:,:), allocatable :: kx,ky,tau,krx_inter,kry_inter
    real(8), dimension(:,:,:), allocatable :: torque_xz,torque_yz,viscous_heat_xz,viscous_heat_yz
    real(8), dimension(:,:,:), allocatable :: viscous_heat,heat1,heat2,torque_work,viscous_heat_dt
    real(8), dimension(:,:,:), allocatable :: radhydro_boost_viscous_heat
    real(8), dimension(:,:,:), allocatable :: omega_xz,omega_yz
    real(8), dimension(:,:,:), allocatable :: Fradx,Frady
    real(8), dimension(:,:,:), allocatable :: surf1,surf2,lever_r,lever_theta
    real(8), dimension(:,:,:), allocatable :: gpotential_surf
    real(8), dimension(:,:,:,:), allocatable :: fld_surfx,fld_surfy
    !Volume quantities, values are at the center of the cells
    real(8), dimension(:,:,:,:), allocatable :: w,u,source
    real(8), dimension(:,:,:,:), allocatable :: u0,u1,u_muscl
    real(8), dimension(:,:,:,:), allocatable :: w_xl,w_xr,w_yl,w_yr
    real(8), dimension(:,:,:), allocatable :: vol,cs,amz,omega0,omega1,omega_muscl,omega,fr_omega,ftheta_omega,l_omega
    real(8), dimension(:,:,:), allocatable :: gpotential,erad_adv
    real(8), dimension(:,:,:), allocatable :: vphi_2dx,vphi_2dy,ek_vphi_nomix,ek_vphi_mix,ekphi_mix_heat
    !W primitive, U conservative, S source form. x,y,z (or i,j,k) and physical quantities
    !W physical quantities: 1=rho,2=vx,3=vy,4=vz,5=p
    !U physical quantities: 1=rho,2=rho vx,3=rho vy,4=rho vz,5=E
    real(8), dimension(:,:,:), allocatable :: egv,temp,iso_temp,cv,entropy
    real(8), dimension(:,:,:,:), allocatable :: xslp,yslp
    real(8), dimension(:,:,:), allocatable :: temp_xl,temp_xr,temp_yl,temp_yr,egv_xl,egv_xr,egv_yl,egv_yr
    real(8), dimension(:,:,:), allocatable :: omega_xl,omega_xr,omega_yl,omega_yr
    real(8), dimension(:,:,:), allocatable :: H2,HI,HII,HeI,HeII,HeIII,electron
    real(8), dimension(:,:,:,:), allocatable :: s_geo1,s_geo2,s_grav
    real(8), dimension(:,:,:), allocatable :: Erad,Erad_int,Erad_source,irrad,fld_heating,fld_heating_cell,fld_heating_viscous,heat_acculm
    real(8), dimension(:,:,:), allocatable :: kappa_planck,kappa_abs,kappa_rosseland                  !unit cm^2/g
    real(8), dimension(:,:,:), allocatable :: sigma_planck,sigma_abs,sigma_rosseland,aradx,arady
    real(8), dimension(:,:,:), allocatable :: stheta,ctheta
    real(8), dimension(:,:,:), allocatable :: rcylin
    real(8), dimension(:,:,:,:), allocatable :: passive_scalar,passive_scalar0
    real(8), dimension(:), allocatable :: x_center,y_center,x_interface,y_interface,irrad0
    real(8), dimension(:), allocatable :: e_guard,s_guard,w_guard,n_guard
    real(8), dimension(2) :: ne_guard,se_guard,sw_guard,nw_guard
    real(8) :: e_guardx,s_guardy,w_guardx,n_guardy
    logical :: derefine
    logical :: amr_intrinsic
    logical :: on_processor
    logical :: fld_refine(4)
    logical :: zel_spike
    logical :: refine_direction(8)
end type blockdef

type environmentdef
    !light data, parameters about the simulation
    real(8) :: h_ratio
    real(8) :: he_ratio
    integer :: nrho
    integer :: ntemp
    real(8) :: rhomin
    real(8) :: rhomax
    real(8) :: tmin
    real(8) :: tmax
    real(8), dimension(:), allocatable :: rho_eos,t_eos
    integer :: table_dim(2)
    integer :: meta_table_dim(3)
    logical :: lenforce_field_quantities
    logical :: firststep
end type environmentdef

type table1d
    real(8), dimension(:), allocatable :: t1d,xlist
    real(8) :: dx
    integer :: nd1
end type table1d

type table2d
    real(8), dimension(:,:), allocatable :: t2d          !values of the table
    real(8), dimension(:), allocatable :: xlist,ylist    !lists of sampled points
    real(8) :: dx,dy                                     !space of sampled points
    integer :: nd1,nd2                                  !sampled points
end type table2d

type uniform_table1d
    real(8) :: xmin,xmax,dx
end type uniform_table1d

type table_group
    type(table2d), dimension(:), allocatable :: table_eos
end type table_group

type particle
    real(8) :: mass
    real(8), dimension(3) :: xyz
end type particle

type stellar
    type(particle) :: core
    real(8) :: teff
    real(8) :: radius
    real(8) :: period
    real(8) :: spin
    real(8) :: luminosity               !blackbody part
    real(8) :: lum_atm                  !atmosphere luminosity
end type stellar

type(table_group), public :: tables
type(particle) :: single_star
type(stellar) :: central_star
logical, public :: EOS_test=.false.

real(8), public :: opacity_gas_rho_min=1d-18
real(8), public :: opacity_gas_rho_max=1d-5
real(8), public :: opacity_gas_pres_min
real(8), public :: opacity_gas_pres_max
real(8), public :: opacity_gas_t_min=50
real(8), public :: opacity_gas_t_max=40000
real(8), public :: viscosity=0d0
real(8), public :: alpha_vis=1d-2
real(8), public :: ralpha=1.5d13
real(8), public :: kr_ratio=1d0
real(8), public :: kp_ratio=1d0
real(8), public :: rho_slope_refine_thresh,rho_slope_derefine_thresh,erad_slope_refine_thresh,erad_slope_derefine_thresh,gradp_refine_thresh,gradp_derefine_thresh
real(8), public :: zeldovich_refine_thresh,zeldovich_derefine_thresh
real(8), public :: erad_bound(4),frad_bound(4)
real(8), public :: inner_floor_temp=3d2
real(8), dimension(:,:), allocatable, public :: temp_division,temp_h,temp_he
real(8), dimension(:,:), allocatable, public :: record_history
real(8), dimension(:,:), allocatable, public :: refine_region
real(8), dimension(:,:), allocatable, public :: derefine_region
real(8), dimension(:,:), allocatable, public :: rad_mg
!real(8), pointer :: field_ptr(:,:)
integer, public :: rad_bound_type(4),hydro_bound_type(4),passive_bound_type(4)
integer :: record_i=0,record_chunk_i=0,record_array_size=0,record_chunk_size=100
integer, public :: metallicity=2                !1~-0.3 (metal poor), 2~0.0 (solar), 3~0.3 (metal rich)
integer, public :: nframe                       !number of frames to be generated
integer, public :: iframe=0                     !current frame number
integer, public :: record_length,j_record
integer, public :: petsc_iter=1
integer, public :: vis_iter=3
integer, public :: radhydro_boost=1
integer, dimension(:), allocatable, public :: np_nblk
type(uniform_table1d) :: tab_rho_t
logical :: division_known,division_hydrogen,division_helium
logical, public :: energy_conserve_formalism=.false.
logical, public :: rad_adv=.false.
logical, public :: change_mesh=.false.
logical, public :: llalpha=.false.
logical, public :: lchange_dust_opacity=.false.
logical, public :: lupdate_smr=.false.
logical, public :: lremove_shock=.false.
logical, public :: positivity_diagnosis=.false.
logical, public :: lrosseland_pres=.false.
logical, public :: lhydro=.true.
real(8) :: timescale=1d0
real(8) :: lengthscale=1d0
real(8), public :: tempfloorscale
logical, public :: lhis_overwrite=.true.
real(8), dimension(:), allocatable, public :: dt_processors
real(8), dimension(:), allocatable, public :: t_target_output
real(8), dimension(:), allocatable, public :: theta_center,theta_interface
real(8), public :: rho_thresh_petsc1=1d-16          !rho lower bound of solving rad-hydro thermal coupling
real(8), public :: rho_thresh_petsc2=1d-16
real(8), public :: rho_thresh_petsc3=1d-9          !rho upper bound of solving rad-hydro thermal coupling
real(8), public :: sig_rosseland_floor=1d-16
real(8), public :: kappa_planck_ceil1=5d0
real(8), public :: kappa_planck_ceil2=1d0
real(8), public :: p_rosseland_atm=1d5
real(8), public :: cv_thresh
real(8), public :: petsc_rtol
real(8), public :: petsc_qratio=2d0
real(8), public :: vis_qratio=3d0
real(8), public :: radhydro_tspan=0d0
real(8), public :: kr_floor=1d-5
real(8), public :: xgeo_a,xgeo_b,xgeo_pn,xgeo_dp,xgeo_h,ygeo_a,ygeo_b,ygeo_pn,ygeo_dp,ygeo_h
real(8), public :: time_boost_start,time_boost_end
character(len=16), protected :: op_min,op_max,op_sum
character(len=32), public :: out_dir='/out'
character(len=32), protected :: xsurf_variables(8),ysurf_variables(6),cell_variables(32)
character(len=128), public :: path_out
character(len=1), dimension(2,2), protected :: A_type,B_type,C_type,D_type
integer, dimension(4,2), protected :: A_seq,B_seq,C_seq,D_seq,Z_seq
type(timedef), public :: time_sys
type(scaledef), public :: scale_sys
type(axbsystem), public :: axb_fld,axb_viscous
!blocks related
integer, public :: nblk_total,nparocessor
integer, public :: nblocks_tree
integer, public :: directional_buffer
integer, public :: npassive=0,nindependentpassive=0,nimplicitpassive=0
type(blockdef), public, pointer :: blk_head,blk_tail,blk_root,llist_head,llist_tail,blk11
type(blockdef), public, pointer :: blk_head_xl_bound,blk_head_xu_bound,blk_head_yl_bound,blk_head_yu_bound
type(environmentdef), public :: environment
type(processor_variables), public :: processor


contains

function field_calculator(blk,i,j)
    !the prototype of field variable calculation
    type(blockdef), pointer :: blk
    real(8) :: field_calculator
    integer :: i,j
end function field_calculator

!prototype subroutines

function collective_oper_blk(blk)
    type(blockdef), pointer :: blk
    real(8) :: collective_oper_blk
end function collective_oper_blk

subroutine condition_hydro(blk,i,j)
    !this prototype is used both to specify hydro conditions inside the computational domain
    !and at the boundaries, mark is used to distinguish different boundary locations or
    !some special cases
    type(blockdef), pointer :: blk
    integer :: i,j
end subroutine condition_hydro

subroutine extractarray(blk,ap)
    !the prototype of subroutine that extract an array from a block
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: ap
end subroutine extractarray

subroutine blk_operator(blk)
    !the prototype of subroutine of blk operator, works for blk_traversal subroutine
    type(blockdef), pointer :: blk
end subroutine blk_operator

subroutine blk_operator_real(blk,x)
    type(blockdef), pointer :: blk
    real(8) :: x
end subroutine blk_operator_real

subroutine blk_operator_int(blk,i)
    type(blockdef), pointer :: blk
    integer :: i
end subroutine blk_operator_int

subroutine blk_operator_logical(blk,l)
    type(blockdef), pointer :: blk
    logical :: l
end subroutine blk_operator_logical

integer function blk_get_attr_int(blk)
    type(blockdef), pointer :: blk
end function blk_get_attr_int

real function blk_get_attr_real(blk)
    type(blockdef), pointer :: blk
end function blk_get_attr_real

logical function blk_get_attr_logical(blk)
    type(blockdef), pointer :: blk
end function blk_get_attr_logical

subroutine initialize_simulation_parameters()
    real(8) :: tfinal,v_boundary(4),h_ratio,he_ratio,rhomin,rhomax,tmin,tmax,refine_xmin,refine_xmax, &
        refine_ymin,refine_ymax,derefine_xmin,derefine_xmax,derefine_ymin,derefine_ymax,r,r1,r2
    integer :: i,j,nrho,ntemp,level,ierr,nc,errcode
    character(len=16) :: temp_name
    logical :: ex
    namelist /meshinfo/ n_domain,lengthscale,tfinal,timescale,CFL,v_boundary,hydro_bound_type,rad_bound_type,    &
        passive_bound_type,nframe,refine_type,nrefine_region,nderefine_region,max_refine_level,restart,restart_iframe
    namelist /phyinfo/ h_ratio,he_ratio,petsc_rtol,llnx,xgeo_h,llny,ygeo_h,lrad_adv,lam_con
    namelist /global_parameters/ igravity,igeometry,icooling,iradiation,nd,nx,ny,    &
        blk_size_nx,blk_size_ny,maw,gamma_gas,path_tables
    namelist /refinement/ refine_xmin,refine_xmax,refine_ymin,refine_ymax,level
    namelist /derefinement/ derefine_xmin,derefine_xmax,derefine_ymin,derefine_ymax
    call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
    call mpi_comm_size(MPI_COMM_WORLD,np,ierr)
    call getcwd(path_root)
    path_out=trim(path_root)//'/out'
    left='left';right='right';east='east';south='south';west='west';north='north';se='se';sw='sw';nw='nw';ne='ne'
    pn_vphi='vphi';pn_am='am';pn_ekphi='ekphi';pn_omega='omega';pn_erad='erad';sdensity='density';strad='trad';stgas='tgas'
    xcoord='xcoord';ycoord='ycoord';next='next';sshear='shear';prolong='prolong';restrict='restrict';non='none'
    previous='previous'
    op_min='min';op_max='max';op_sum='sum'
    all_direction=(/east,south,west,north,se,sw,nw,ne/)
    corner_direction=(/se,sw,nw,ne/)
    edge_direction=(/east,south,west,north/)
    hydro_field='hydro'
    fld_field='fld'
    passive_field='passive'
    viscous_field='viscous'
    erad_field='rad'
    predictor='predictor'
    corrector='corrector'
    A_type=reshape(['D','B','A','A'],[2,2])
    B_type=reshape(['B','A','B','C'],[2,2])
    C_type=reshape(['C','C','D','B'],[2,2])
    D_type=reshape(['A','D','C','D'],[2,2])
    A_seq=transpose(reshape([1,1,1,2,2,2,2,1],[2,4]))
    B_seq=transpose(reshape([2,2,1,2,1,1,2,1],[2,4]))
    C_seq=transpose(reshape([2,2,2,1,1,1,1,2],[2,4]))
    D_seq=transpose(reshape([1,1,2,1,2,2,1,2],[2,4]))
    Z_seq=transpose(reshape([1,1,2,1,1,2,2,2],[2,4]))
    allocate(np_nblk(np))
    if (rank==0) then
        open(unit=11,file=trim(path_root)//'/global.data',status='old',action='read')
        read(unit=11,nml=meshinfo)
        read(unit=11,nml=phyinfo)
        read(unit=11,nml=global_parameters)
        if (igeometry==1.or.igeometry==2) then
            n_domain(3:4)=n_domain(3:4)*pi
            n_domain(1:2)=n_domain(1:2)*lengthscale
            !print *,n_domain
        end if
        if (igeometry==0) then
            n_domain=n_domain*lengthscale
        end if
    end if
    nz=1;blk_size_nz=1
    !the order should be consistent with the order of global.data
    call mpi_bcast(n_domain,4,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lengthscale,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(tfinal,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(timescale,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(CFL,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(v_boundary,4,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hydro_bound_type,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rad_bound_type,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(passive_bound_type,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nframe,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(refine_type,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nrefine_region,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nderefine_region,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(max_refine_level,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(restart_iframe,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(h_ratio,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(he_ratio,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(petsc_rtol,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(llnx,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(xgeo_h,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(llny,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ygeo_h,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lrad_adv,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lam_con,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(igravity,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(igeometry,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(icooling,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(iradiation,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(blk_size_nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(blk_size_ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(maw,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(gamma_gas,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(path_tables,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    erad_bound=0d0
    frad_bound=0d0
    if (restart.and.restart_iframe==-1) then
        if (rank==0) then
            inquire(file="last.dat",exist=ex)
            if (ex.eqv..false.) then
                print *,'there is no recorded last frame number'
                stop
            else
                open(unit=89,file="last.dat",status="old",action="read")
                read(89,*) iframe
                close(89)
            end if
        end if
        call mpi_bcast(iframe,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    else if (restart) then
        iframe=restart_iframe
    else
        iframe=0
    end if
    rho_slope_refine_thresh=0.2
    rho_slope_derefine_thresh=0.1
    erad_slope_refine_thresh=0.03
    erad_slope_derefine_thresh=0.02
    gradp_refine_thresh=0.4
    gradp_derefine_thresh=0.3
    zeldovich_refine_thresh=1.0001d0
    zeldovich_derefine_thresh=1d0
    directional_buffer=1
    tempfloorscale=1
    time_sys%ntimestep=1
    time_sys%t_final=tfinal*timescale
    time_sys%dt_record=0d0
    time_sys%dt_hydro=0d0
    cell_var_length=5
    if (refine_type=='none') then                                       !no mesh refinement
        max_refine_level=1
    else if (refine_type=='static'.or.refine_type=='mixed') then        !static mesh refinement
        allocate(refine_region(nrefine_region,5))
        do i=1,nrefine_region
            if (rank==0) then
                read(11,nml=refinement)
                refine_region(i,1)=refine_xmin
                refine_region(i,2)=refine_xmax
                refine_region(i,3)=refine_ymin
                refine_region(i,4)=refine_ymax
                if (igeometry==1.or.igeometry==2) then
                    refine_region(i,3:4)=refine_region(i,3:4)*pi
                    refine_region(i,1:2)=refine_region(i,1:2)*lengthscale
                end if
                if (igeometry==0) then
                    refine_region(i,1:4)=refine_region(i,1:4)*lengthscale
                end if
                refine_region(i,5)=level
            end if
        end do
        if (np>1) then
            nc=nrefine_region*5
            call mpi_bcast(refine_region,nc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        end if
        allocate(derefine_region(nderefine_region,4))
        do i=1,nderefine_region
            if (rank==0) then
                read(11,nml=derefinement)
                derefine_region(i,1)=derefine_xmin
                derefine_region(i,2)=derefine_xmax
                derefine_region(i,3)=derefine_ymin
                derefine_region(i,4)=derefine_ymax
                if (igeometry==1.or.igeometry==2) then
                    derefine_region(i,3:4)=derefine_region(i,3:4)*pi
                    derefine_region(i,1:2)=derefine_region(i,1:2)*lengthscale
                end if
                if (igeometry==0) then
                    derefine_region(i,1:4)=derefine_region(i,1:4)*lengthscale
                end if
            end if
        end do
        if (np>1) then
            nc=nderefine_region*4
            call mpi_bcast(derefine_region,nc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        end if
    end if
    close(11)
    dxyz=(/(n_domain(2)-n_domain(1))/nx,(n_domain(4)-n_domain(3))/ny,0d0/)
    nrho=1201;ntemp=1001;rhomin=1d-20;rhomax=1d0;tmin=5d2;tmax=2d5
    environment%h_ratio=h_ratio
    environment%he_ratio=he_ratio
    environment%nrho=nrho
    environment%ntemp=ntemp
    environment%rhomin=rhomin
    environment%rhomax=rhomax
    environment%tmin=tmin
    environment%tmax=tmax
    environment%firststep=.true.
    xsurf_variables=[character(len=32) :: 'Fradx','aradx','viscous_heat_xz','torque_xz','mass_xflux','am_xflux','torque_work',   &
    'hydro_energy_xflux']
    ysurf_variables=[character(len=32) :: 'Frady','arady','viscous_heat_yz','torque_yz','mass_yflux','am_yflux']
    cell_variables=[character(len=32) :: 'rho','vx','vy','vz','pres','cs','temp','egv','entropy','vphi','omega','irrad',    &
    'H2','HI','HII','HeI','HeII','HeIII','Rosseland','sigma_Rosseland','Planck','sigma_Planck','viscous_heat','radhydro_boost_viscous_heat',      &
    'Erad','Erad_int','Erad_source','isotemp','fld_heating','fld_heating_viscous','heat2','ekphi']
    if (igeometry==2) then
        allocate(theta_center(-1:ny/2),theta_interface(-2:ny/2))
        do i=-2,ny/2
            theta_interface(i)=n_domain(3)+i*dxyz(2)
        end do
        do i=-1,ny/2
            theta_center(i)=n_domain(3)+(i-half)*dxyz(2)
        end do
        do i=-1,nx+2
            r1=n_domain(1)+(i-1)*dxyz(1)
            r2=r1+dxyz(1)
        end do
    end if
    allocate(dt_processors(np),processor%processor_operation_sum(np))
    allocate(processor%collective_array(np),processor%global_logical(np))
    call check_global_parameter_consistency()
    call initialize_geometry()
    call initialize_mesh()
end subroutine initialize_simulation_parameters

!subroutine update_time_sys(timescale,dt_record)
!    real(8) :: timescale,dt_record
!    time_sys%ntimestep=0
!    time_sys%t_final=time_sys%t_final*timescale
!    time_sys%dt_record=0d0
!    time_sys%dt_hydro=0d0
!end subroutine update_time_sys

subroutine create_table_1d(m,table)
    integer :: m
    type(table1d) :: table
    allocate(table%t1d(m),table%xlist(m))
    table%nd1=m
end subroutine create_table_1d

subroutine create_table_2d(m,n,table)
    integer :: m,n
    type(table2d) :: table
    allocate(table%t2d(m,n),table%xlist(m),table%ylist(n))
    table%nd1=m
    table%nd2=n
end subroutine create_table_2d

subroutine calculate_llevel_max()
    integer :: nxblock,nyblock
    character(len=128) :: alert
    if (mod(nx,blk_size_nx)/=0) then
        alert='block size is not a divider of domain size'
        call abort_guangqi(alert)
    end if
    nxblock=nx/blk_size_nx
    llevel_max=0
    do while (nxblock/=1)
        nxblock=ceiling(real(nxblock)/2)
        llevel_max=llevel_max+1
    end do
end subroutine calculate_llevel_max

subroutine block_key_to_pointer(key,blk)
    !given a block key, return the pointer to the key (blk)
    !if such a block does not exist, return the *closest* parent key
    type(blockdef), pointer :: blk
    integer :: key(3),i,llevel,idx,idy
    integer, allocatable :: seq(:),seqx(:),seqy(:)
    logical :: nonexist
    allocate(seq(key(3)));seq=-1
    i=1;llevel=key(3);idx=key(1)
    do while (llevel>0)
        seq(i)=2-mod(idx,2)
        i=i+1
        llevel=llevel-1
        idx=ceiling(real(idx)/2)
    end do
    blk=>blk_root
    do i=key(3),1,-1
        idx=seq(i);idy=1
        if (associated(blk%children(idx,idy)%blk)) then
            blk=>blk%children(idx,idy)%blk
        else
            exit
        end if
    end do
    deallocate(seq)
end subroutine block_key_to_pointer

function key_parent(key)
    integer :: key(3),key_parent(3)
    key_parent(1)=ceiling(real(key(1))/2)
    key_parent(2)=key(2)
    key_parent(3)=key(3)-1
end function key_parent

recursive subroutine release_all_not_on_processor(blk)
    type(blockdef), pointer :: blk
    integer :: idx,idy
    if (blk%on_processor.eqv..false.) call release_block(blk)
    if (associated(blk%children(1,1)%blk)) then
        do idy=1,2
            do idx=1,2
                if (associated(blk%children(idx,idy)%blk)) then
                    call release_all_not_on_processor(blk%children(idx,idy)%blk)
                end if
            end do
        end do
    end if
end subroutine release_all_not_on_processor

subroutine identify_on_processor()
    call negate_on_processor(blk_root)
    call blk_traversal(on_processor)
end subroutine identify_on_processor

recursive subroutine negate_on_processor(blk)
    type(blockdef), pointer :: blk
    integer :: idx,idy
    blk%on_processor=.false.
    idx=curve_seq_idx(blk%curve,1)
    idy=curve_seq_idy(blk%curve,1)
    if (associated(blk%children(idx,idy)%blk)) then
        do idy=1,2
            do idx=1,2
                if (associated(blk%children(idx,idy)%blk)) then
                    call negate_on_processor(blk%children(idx,idy)%blk)
                end if
            end do
        end do
    end if
end subroutine negate_on_processor

subroutine on_processor(blk)
    type(blockdef), pointer :: blk
    integer :: id_head,id_tail,i
    if (rank>0) then
        id_head=sum(np_nblk(1:rank))+1
        id_tail=sum(np_nblk(1:rank+1))
    else
        id_head=1;id_tail=np_nblk(1)
    end if
    if (blk%blk_id>=id_head.and.blk%blk_id<=id_tail) then
        blk%on_processor=.true.
    else
        blk%on_processor=.false.
    end if
end subroutine on_processor

subroutine check_flying_blk_data()
    type(blockdef), pointer :: blk
    integer :: i,ncount(np),ierr,ntotal
    i=0;call verify_on_processor(blk_root,i)
    call mpi_gather(i,1,MPI_INTEGER,ncount,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ncount,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    ntotal=sum(abs(ncount-np_nblk))
    if (ntotal/=0) then
        print *,'there is allocated but not used blk'
        print *,np_nblk
        print *,ncount
        stop
    end if
end subroutine check_flying_blk_data

recursive subroutine verify_on_processor(blk,n)
    type(blockdef), pointer :: blk
    integer :: i,n,idx,idy
    if (allocated(blk%w)) n=n+1
    do i=1,2**nd
        idx=curve_seq_idx(blk%curve,i)
        idy=curve_seq_idy(blk%curve,i)
        if (associated(blk%children(idx,idy)%blk)) then
            call verify_on_processor(blk%children(idx,idy)%blk,n)
        end if
    end do
    !end if
end subroutine verify_on_processor

subroutine display_llist()
    type(blockdef), pointer :: blk
    blk=>llist_head
    do while (.not.associated(blk,llist_tail))
        print *,rank,blk%key
        blk=>blk%next
    end do
    print *,rank,blk%key
end subroutine display_llist

subroutine domain_block_parameters()
    integer :: i
    real(8) :: dx,dy
    if (llnx) then
    else
        dx=n_domain(2)-n_domain(1)
    end if
    if (llntheta) then
    else
        dy=n_domain(4)-n_domain(3)
    end if
    base_blk_size=(/dx/nx_blks,dy/ny_blks/)
end subroutine domain_block_parameters

subroutine new_blank_block(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    allocate(blk)
    blk%blk_id=-1
    blk%level=-1
    blk%static_level=0
    nullify(blk%pre,blk%next,blk%pa,blk%child_head,blk%child_tail)
    do i=1,2
        do j=1,2
            nullify(blk%children(i,j)%blk)
        end do
    end do
    do i=1,12
        blk%nb(i)%cpu_rank=-1
        nullify(blk%nb(i)%blk)
    end do
    blk%derefine=.false.
    blk%amr_intrinsic=.false.
    blk%refine_countdown=0
    blk%refine_direction=.false.
    blk%fld_refine=.false.
    blk%zel_spike=.false.
    blk%derefine_countdown=0
    blk%on_processor=.false.
    blk%curve='Z'
end subroutine new_blank_block

subroutine block_mesh_data(blk)
    !from blk%key, calculate blk%pos and level
    type(blockdef), pointer :: blk
    integer :: i
    real(8) :: q,dtheta_level_base
    blk%level=blk%key(3)-llevel_max
    blk%dxyz(1)=dxyz(1)/2d0**blk%level
    if (llntheta) then
    else
        blk%dxyz(2)=dxyz(2)/2d0**blk%level
    end if
    blk%dxyz(3)=0d0
    allocate(blk%mesh_x(blk_size_nx+1))
    if (llnx) then
        do i=1,blk_size_nx+1
            blk%mesh_x(i)=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level,blk%key(1),i-1)
        end do
    else
        blk%pos(1)=n_domain(1)+(blk%key(1)-1)*blk%dxyz(1)*blk_size_nx
        blk%mesh_x=blk%pos(1)+(/(((i-1)*blk%dxyz(1)),i=1,blk_size_nx+1)/)
    end if
    if (llntheta) then
    else
        blk%pos(2)=n_domain(3)+(blk%key(2)-1)*blk%dxyz(2)*blk_size_ny
    end if
    blk%pos(3)=0
end subroutine block_mesh_data

subroutine display_block_data(blk)
    !including the buffer zones
    type(blockdef), pointer :: blk
    character(len=32) :: s1,fmt1
    integer :: nbx,nby,i
end subroutine display_block_data

subroutine check_leaves()
    type(blockdef), pointer :: blk
    real(8) :: r
    nullify(blk)
    blk=>blk_head
    do while (associated(blk))
        if (dxyz(1)/blk%dxyz(1)/=2**blk%level) then
            print *,'level wrong',blk%key,blk%level,dxyz(1)/blk%dxyz(1)
            stop
        end if
        if (associated(blk%pre)) then
            if (abs(blk%pre%level-blk%level)>1) then
                print *,'1,check linked blocks, pre',blk%key,blk%pre%key
                stop
            end if
            r=blk%pre%dxyz(1)/blk%dxyz(1)
            if (r>2.or.r<0.5) then
                print *,'2,check linked blocks, pre',blk%key,blk%pre%key,r
                stop
            end if
        end if
        if (associated(blk%next)) then
            if (abs(blk%next%level-blk%level)>1) then
                print *,'1,check linked blocks, next',blk%key,blk%next%key
                stop
            end if
            r=blk%next%dxyz(1)/blk%dxyz(1)
            if (r>2.or.r<0.5) then
                print *,'2,check linked blocks, next',blk%key,blk%next%key,r
                stop
            end if
        end if
        blk=>blk%next
    end do
end subroutine check_leaves

subroutine print_block_quantity(s)
    type(blockdef), pointer :: blk
    character(len=*) :: s
    blk=>blk_head
    print *,trim(s),time_sys%t,time_sys%ntimestep
    do while (associated(blk))
        select case (trim(s))
        case ('rho')
            write(*,'(2I8,8ES15.4E2)') blk%blk_id,blk%level,blk%w(1,1:blk_size_nx,1,1)
        case ('vx')
            write(*,'(2I8,8ES15.4E2)') blk%blk_id,blk%level,blk%w(2,1:blk_size_nx,1,1)
        case ('p')
            write(*,'(2I8,8ES15.4E2)') blk%blk_id,blk%level,blk%w(5,1:blk_size_nx,1,1)
        case ('temp')
            write(*,'(2I8,8ES15.4E2)') blk%blk_id,blk%level,blk%temp(1:blk_size_nx,1,1)
        end select
        blk=>blk%next
    end do
end subroutine print_block_quantity

subroutine decrease_static_level(blk)
    type(blockdef), pointer :: blk
    if (time_sys%ntimestep==20.and.refine_type=='mixed'.and.blk%static_level>0) then
        blk%static_level=0
    end if
end subroutine decrease_static_level

logical function block_zone_overlap(blk,zone)
    !logical :: block_zone_overlap
    type(blockdef), pointer :: blk
    real(8) :: zone(4),xmin,xmax,ymin,ymax
    xmin=blk%mesh_x(1)
    xmax=blk%mesh_x(blk_size_nx+1)
    if (xmin<=zone(2).and.xmax>=zone(1)) then
        block_zone_overlap=.true.
    else
        block_zone_overlap=.false.
    end if
end function block_zone_overlap

logical function block_zone_inside(blk,zone)
    type(blockdef), pointer :: blk
    real(8) :: zone(4),xmin,xmax,ymin,ymax
end function block_zone_inside

subroutine calculate_load(nblk)
    integer, allocatable :: nblk(:)
    integer :: i,blks_per_core,remainder
    blks_per_core=nblk_total/np
    nblk=blks_per_core
    remainder=mod(nblk_total,np)
    if (remainder>0) then
        do i=np,np-remainder+1,-1
            nblk(i)=nblk(i)+1
        end do
    end if
end subroutine calculate_load

subroutine reset_mr_processor_variables()
    processor%mr_ngrow=0
    processor%mr_nderefine=0
    if (.not.allocated(processor%mpi_logical)) allocate(processor%mpi_logical(np))
    processor%mpi_logical=.false.
    if (.not.allocated(processor%mpi_processor_oper)) allocate(processor%mpi_processor_oper(np))
    processor%mpi_processor_oper=0
    if (allocated(processor%grow_keys)) deallocate(processor%grow_keys)
    if (allocated(processor%derefine_keys)) deallocate(processor%derefine_keys)
    if (allocated(processor%mpi_world_keys)) deallocate(processor%mpi_world_keys)
end subroutine reset_mr_processor_variables

subroutine reset_smr_processor_variables()
    !upon restart, the smr structure may change
    processor%smr_ngrow=0
    processor%smr_nderefine=0
    if (.not.allocated(processor%mpi_logical)) allocate(processor%mpi_logical(np))
    processor%mpi_logical=.false.
    if (.not.allocated(processor%mpi_processor_oper)) allocate(processor%mpi_processor_oper(np))
    processor%mpi_processor_oper=0
    if (allocated(processor%grow_keys)) deallocate(processor%grow_keys)
    if (allocated(processor%derefine_keys)) deallocate(processor%derefine_keys)
    if (allocated(processor%mpi_world_keys)) deallocate(processor%mpi_world_keys)
end subroutine reset_smr_processor_variables

subroutine reset_amr_processor_variables()
    processor%amr_ngrow=0
    processor%amr_nderefine=0
    if (.not.allocated(processor%mpi_logical)) allocate(processor%mpi_logical(np))
    processor%mpi_logical=.false.
    if (.not.allocated(processor%mpi_processor_oper)) allocate(processor%mpi_processor_oper(np))
    processor%mpi_processor_oper=0
    if (allocated(processor%grow_keys)) deallocate(processor%grow_keys)
    if (allocated(processor%derefine_keys)) deallocate(processor%derefine_keys)
    if (allocated(processor%mpi_world_keys)) deallocate(processor%mpi_world_keys)
end subroutine reset_amr_processor_variables

function opposite_direction(direction)
    !serve as a dictionary of opposite directions
    character(len=16) :: direction,opposite_direction
    if (direction==east) opposite_direction=west
    if (direction==south) opposite_direction=north
    if (direction==west) opposite_direction=east
    if (direction==north) opposite_direction=south
    if (direction==se) opposite_direction=nw
    if (direction==sw) opposite_direction=ne
    if (direction==nw) opposite_direction=se
    if (direction==ne) opposite_direction=sw
    if (direction==next) opposite_direction=previous
    if (direction==previous) opposite_direction=next
end function opposite_direction

integer function direc_index(direction)
    !a dictionary from direction to index
    character(len=16) :: direction
    if (direction==east) direc_index=1
    if (direction==south) direc_index=2
    if (direction==west) direc_index=3
    if (direction==north) direc_index=4
    if (direction==se) direc_index=5
    if (direction==sw) direc_index=6
    if (direction==nw) direc_index=7
    if (direction==ne) direc_index=8
end function direc_index

integer function direc_nb_index(direction)
    !serve as a dictionary from direction to neighbor index
    !there could be 2 neighbors in the same direction
    character(len=16) :: direction
    !integer :: direc_nb_index
    if (direction==east) direc_nb_index=1
    if (direction==south) direc_nb_index=3
    if (direction==west) direc_nb_index=5
    if (direction==north) direc_nb_index=7
    if (direction==se) direc_nb_index=9
    if (direction==sw) direc_nb_index=10
    if (direction==nw) direc_nb_index=11
    if (direction==ne) direc_nb_index=12
end function direc_nb_index

function index_direc(id)
    !a dictionary of neighbors from index to direction
    character(len=16) :: index_direc
    integer :: id
    if (id==1) index_direc=east
    if (id==2) index_direc=south
    if (id==3) index_direc=west
    if (id==4) index_direc=north
    if (id==5) index_direc=se
    if (id==6) index_direc=sw
    if (id==7) index_direc=nw
    if (id==8) index_direc=ne
end function index_direc

function nb_index_direc(id)
    !serve as a dictionary of neighbors from neighbor index to direction
    !there could be 2 neighbors in the same direction
    character(len=16) :: nb_index_direc
    integer :: id
    if (id==1.or.id==2) nb_index_direc=east
    if (id==3.or.id==4) nb_index_direc=south
    if (id==5.or.id==6) nb_index_direc=west
    if (id==7.or.id==8) nb_index_direc=north
    if (id==9) nb_index_direc=se
    if (id==10) nb_index_direc=sw
    if (id==11) nb_index_direc=nw
    if (id==12) nb_index_direc=ne
end function nb_index_direc

subroutine direc_to_nb_key(blk,direction,key_nb)
    !given the blk and a direction, return the key of the neighbor
    type(blockdef), pointer :: blk
    integer :: key(3),key_nb(3),level
    character(len=16) :: direction
    key=blk%key
    level=blk%level
    key_nb=0
    if (direction==east) then
        if (key(1)==nx_blks*2**blk%level.and.hydro_bound_type(2)==3) then
            key_nb=(/1,key(2),key(3)/)
        else if (key(1)/=nx_blks*2**blk%level) then
            key_nb=(/key(1)+1,key(2),key(3)/)
        end if
    end if
    if (direction==south) then
        if (key(2)==1.and.hydro_bound_type(3)==3) then
            key_nb=(/key(1),ny_blks*2**(key(3)-llevel_max),key(3)/)
        else if (key(2)/=1) then
            key_nb=(/key(1),key(2)-1,key(3)/)
        end if
    end if
    if (direction==west) then
        if (key(1)==1.and.hydro_bound_type(1)==3) then
            key_nb=(/nx_blks*2**(key(3)-llevel_max),key(2),key(3)/)
        else if (key(1)/=1) then
            key_nb=(/key(1)-1,key(2),key(3)/)
        end if
    end if
    if (direction==north) then
        if (key(2)==ny_blks*2**blk%level.and.hydro_bound_type(4)==3) then
            key_nb=(/key(1),1,key(3)/)
        else if (key(2)/=ny_blks*2**blk%level) then
            key_nb=(/key(1),key(2)+1,key(3)/)
        end if
    end if
    if (direction==se) then
        if (key(1)==nx_blks*2**level.and.key(2)/=1) then
            if (hydro_bound_type(2)==3) then
                key_nb=(/1,key(2)-1,key(3)/)
            else if (hydro_bound_type(2)==1.or.hydro_bound_type(2)==2) then
                key_nb=(/key(1),key(2)-1,key(3)/)
            end if
        else if (key(1)/=nx_blks*2**level.and.key(2)==1) then
            if (hydro_bound_type(3)==3) then
                key_nb=(/key(1)+1,ny_blks*2**(key(3)-llevel_max),key(3)/)
            else if (hydro_bound_type(3)==1.or.hydro_bound_type(3)==2) then
                key_nb=(/key(1)+1,key(2),key(3)/)
            end if
        else if (key(1)==nx_blks*2**level.and.key(2)==1) then
            if (hydro_bound_type(2)==3.and.hydro_bound_type(3)==3) then
                key_nb=(/1,ny_blks*2**(key(3)-llevel_max),key(3)/)
            else if (hydro_bound_type(2)==3.and.hydro_bound_type(3)/=3) then
                key_nb=(/1,key(2),key(3)/)
            else if (hydro_bound_type(2)/=3.and.hydro_bound_type(3)==3) then
                key_nb=(/key(1),ny_blks*2**(key(3)-llevel_max),key(3)/)
            end if
        else
            key_nb=(/key(1)+1,key(2)-1,key(3)/)
        end if
    end if
    if (direction==sw) then
        if (key(1)==1.and.key(2)/=1) then
            if (hydro_bound_type(1)==3) then
                key_nb=(/nx_blks*2**(key(3)-llevel_max),key(2)-1,key(3)/)
            else if (hydro_bound_type(1)==1.or.hydro_bound_type(1)==2) then
                key_nb=(/key(1),key(2)-1,key(3)/)
            end if
        else if (key(1)/=1.and.key(2)==1) then
            if (hydro_bound_type(3)==3) then
                key_nb=(/key(1)-1,ny_blks*2**(key(3)-llevel_max),key(3)/)
            else if (hydro_bound_type(3)==1.or.hydro_bound_type(3)==2) then
                key_nb=(/key(1)-1,key(2),key(3)/)
            end if
        else if (key(1)==1.and.key(2)==1) then
            if (hydro_bound_type(1)==3.and.hydro_bound_type(3)==3) then
                key_nb=(/nx_blks*2**(key(3)-llevel_max),ny_blks*2**(key(3)-llevel_max),key(3)/)
            else if (hydro_bound_type(1)==3.and.hydro_bound_type(3)/=3) then
                key_nb=(/nx_blks*2**(key(3)-llevel_max),key(2),key(3)/)
            else if (hydro_bound_type(1)/=3.and.hydro_bound_type(3)==3) then
                key_nb=(/key(1),ny_blks*2**(key(3)-llevel_max),key(3)/)
            end if
        else
            key_nb=(/key(1)-1,key(2)-1,key(3)/)
        end if
    end if
    if (direction==nw) then
        if (key(1)==1.and.key(2)/=ny_blks*2**level) then
            if (hydro_bound_type(1)==3) then
                key_nb=(/nx_blks*2**(key(3)-llevel_max),key(2)+1,key(3)/)
            else if (hydro_bound_type(1)==1.or.hydro_bound_type(1)==2) then
                key_nb=(/key(1),key(2)+1,key(3)/)
            end if
        else if (key(1)/=1.and.key(2)==ny_blks*2**level) then
            if (hydro_bound_type(4)==3) then
                key_nb=(/key(1)-1,1,key(3)/)
            else if (hydro_bound_type(4)==1.or.hydro_bound_type(4)==2) then
                key_nb=(/key(1)-1,key(2),key(3)/)
            end if
        else if (key(1)==1.and.key(2)==ny_blks*2**level) then
            if (hydro_bound_type(1)==3.and.hydro_bound_type(4)==3) then
                key_nb=(/nx_blks*2**(key(3)-llevel_max),1,key(3)/)
            else if (hydro_bound_type(1)==3.and.hydro_bound_type(4)/=3) then
                key_nb=(/nx_blks*2**(key(3)-llevel_max),key(2),key(3)/)
            else if (hydro_bound_type(1)/=3.and.hydro_bound_type(4)==3) then
                key_nb=(/key(1),1,key(3)/)
            end if
        else
            key_nb=(/key(1)-1,key(2)+1,key(3)/)
        end if
    end if
    if (direction==ne) then
        if (key(1)==nx_blks*2**level.and.key(2)/=ny_blks*2**level) then
            if (hydro_bound_type(2)==3) then
                key_nb=(/1,key(2)+1,key(3)/)
            else if (hydro_bound_type(2)==1.or.hydro_bound_type(2)==2) then
                key_nb=(/key(1),key(2)+1,key(3)/)
            end if
        else if (key(1)/=nx_blks*2**level.and.key(2)==ny_blks*2**level) then
            if (hydro_bound_type(4)==3) then
                key_nb=(/key(1)+1,1,key(3)/)
            else if (hydro_bound_type(4)==1.or.hydro_bound_type(4)==2) then
                key_nb=(/key(1)+1,key(2),key(3)/)
            end if
        else if (key(1)==nx_blks*2**level.and.key(2)==ny_blks*2**level) then
            if (hydro_bound_type(2)==3.and.hydro_bound_type(4)==3) then
                key_nb=(/1,1,key(3)/)
            else if (hydro_bound_type(2)==3.and.hydro_bound_type(4)/=3) then
                key_nb=(/1,key(2),key(3)/)
            else if (hydro_bound_type(2)/=3.and.hydro_bound_type(4)==3) then
                key_nb=(/key(1),1,key(3)/)
            end if
        else
            key_nb=(/key(1)+1,key(2)+1,key(3)/)
        end if
    end if
end subroutine direc_to_nb_key

function corner_name(i)
    integer :: i
    character(len=16) :: corner_name
    if (i==1) then
        corner_name=se
    else if (i==2) then
        corner_name=sw
    else if (i==3) then
        corner_name=nw
    else if (i==4) then
        corner_name=ne
    end if
end function corner_name

function corner_type(blk,direction)
    !the corner type is with respect to the computational domain
    !1 2 3
    !4 5 6
    !7 8 9
    type(blockdef), pointer :: blk
    integer :: i,j,ix,iy,corner_type
    character(len=16) :: direction
    if (direction==se) then
        i=blk_size_nx;j=1
    else if (direction==sw) then
        i=1;j=1
    else if (direction==nw) then
        i=1;j=blk_size_ny
    else if (direction==ne) then
        i=blk_size_nx;j=blk_size_ny
    end if
    ix=(blk%key(1)-1)*blk_size_nx+i
    iy=(blk%key(2)-1)*blk_size_ny+j
    if (ix==1) then
        if (iy==1) then
            corner_type=7
        else if (iy==ny*2**blk%level) then
            corner_type=1
        else
            corner_type=4
        end if
    else if (ix==nx*2**blk%level) then
        if (iy==1) then
            corner_type=9
        else if (iy==ny*2**blk%level) then
            corner_type=3
        else
            corner_type=6
        end if
    else
        if (iy==1) then
            corner_type=8
        else if (iy==ny*2**blk%level) then
            corner_type=2
        else
            corner_type=5
        end if
    end if
end function corner_type

function corner_neighbor_direction(blk,direction,bound_type)
    !given blk and direction, return its corresponding neighbour's corner direction
    !1-se, 2-sw, 3-nw, 4-ne, 0-the neighbor is itself, or with a given boundary condition
    type(blockdef), pointer :: blk
    character(len=16) :: direction
    integer :: corner_neighbor_direction,bound_type(4)
    corner_neighbor_direction=0
    select case (corner_type(blk,direction))
    case (1)
        if (bound_type(1)==3.and.bound_type(4)==3) then
            corner_neighbor_direction=1
        else if (bound_type(1)==3) then
            corner_neighbor_direction=4
        else if (bound_type(4)==3) then
            corner_neighbor_direction=2
        end if
    case (2)
        if (direction==ne) then
            if (bound_type(4)==3) then
                corner_neighbor_direction=2
            else if (bound_type(4)==1.or.bound_type(4)==2) then
                corner_neighbor_direction=3
            end if
        else if (direction==nw) then
            if (bound_type(4)==3) then
                corner_neighbor_direction=1
            else if (bound_type(4)==1.or.bound_type(4)==2) then
                corner_neighbor_direction=4
            end if
        end if
    case (3)
        if (bound_type(2)==3.and.bound_type(4)==3) then
            corner_neighbor_direction=2
        else if (bound_type(2)==3) then
            corner_neighbor_direction=3
        else if (bound_type(4)==3) then
            corner_neighbor_direction=1
        end if
    case (4)
        if (direction==nw) then
            if (bound_type(1)==3) then
                corner_neighbor_direction=1
            else if (bound_type(1)==1.or.bound_type(1)==2) then
                corner_neighbor_direction=2
            end if
        else if (direction==sw) then
            if (bound_type(1)==3) then
                corner_neighbor_direction=4
            else if (bound_type(1)==1.or.bound_type(1)==2) then
                corner_neighbor_direction=3
            end if
        end if
    case (5)
        if (direction==se) then
            corner_neighbor_direction=3
        else if (direction==sw) then
            corner_neighbor_direction=4
        else if (direction==nw) then
            corner_neighbor_direction=1
        else if (direction==ne) then
            corner_neighbor_direction=2
        end if
    case (6)
        if (direction==se) then
            if (bound_type(2)==3) then
                corner_neighbor_direction=3
            else if (bound_type(2)==1.or.bound_type(2)==2) then
                corner_neighbor_direction=4
            end if
        else if (direction==ne) then
            if (bound_type(2)==3) then
                corner_neighbor_direction=2
            else if (bound_type(2)==1.or.bound_type(2)==2) then
                corner_neighbor_direction=1
            end if
        end if
    case (7)
        if (bound_type(1)==3.and.bound_type(3)==3) then
            corner_neighbor_direction=4
        else if (bound_type(1)==3) then
            corner_neighbor_direction=1
        else if (bound_type(3)==3) then
            corner_neighbor_direction=3
        end if
    case (8)
        if (direction==se) then
            if (bound_type(3)==3) then
                corner_neighbor_direction=3
            else if (bound_type(3)==1.or.bound_type(3)==2) then
                corner_neighbor_direction=2
            end if
        else if (direction==sw) then
            if (bound_type(3)==3) then
                corner_neighbor_direction=4
            else if (bound_type(3)==1.or.bound_type(3)==2) then
                corner_neighbor_direction=1
            end if
        end if
    case (9)
        if (bound_type(2)==3.and.bound_type(3)==3) then
            corner_neighbor_direction=3
        else if (bound_type(2)==3) then
            corner_neighbor_direction=2
        else if (bound_type(3)==3) then
            corner_neighbor_direction=4
        end if
    end select
end function corner_neighbor_direction

function curve_type(curve,idx,idy)
    character(len=1) :: curve,curve_type
    integer :: idx,idy
    if (curve=='A') then
        curve_type=A_type(idx,idy)
    else if (curve=='B') then
        curve_type=B_type(idx,idy)
    else if (curve=='C') then
        curve_type=C_type(idx,idy)
    else if (curve=='D') then
        curve_type=D_type(idx,idy)
    else if (curve=='Z') then
        curve_type='Z'
    end if
end function curve_type

function curve_seq_idx(curve,i)
    integer :: i,curve_seq_idx
    character(len=1) :: curve
    if (curve=='A') then
        curve_seq_idx=A_seq(i,1)
    else if (curve=='B') then
        curve_seq_idx=B_seq(i,1)
    else if (curve=='C') then
        curve_seq_idx=C_seq(i,1)
    else if (curve=='D') then
        curve_seq_idx=D_seq(i,1)
    else if (curve=='Z') then
        curve_seq_idx=Z_seq(i,1)
    end if
end function curve_seq_idx

function curve_seq_idy(curve,i)
    integer :: i,curve_seq_idy
    character(len=1) :: curve
    if (curve=='A') then
        curve_seq_idy=A_seq(i,2)
    else if (curve=='B') then
        curve_seq_idy=B_seq(i,2)
    else if (curve=='C') then
        curve_seq_idy=C_seq(i,2)
    else if (curve=='D') then
        curve_seq_idy=D_seq(i,2)
    else if (curve=='Z') then
        curve_seq_idy=Z_seq(i,2)
    end if
end function curve_seq_idy

subroutine allocate_block_heavy_data(blk)
    type(blockdef), pointer :: blk
    call allocate_geometry_block(blk)
    !call allocate_implicit_mesh(blk)
    call allocate_hydro_block(blk)
    call allocate_eos_block(blk)
    call allocate_radiation_block(blk)
    call allocate_source_block(blk)
    call allocate_viscosity_block(blk)
    call allocate_passive_block(blk)
end subroutine allocate_block_heavy_data

subroutine release_block(blk)
    !deallocate all heavy data
    type(blockdef), pointer :: blk
    call deallocate_geometry_block(blk)
    call deallocate_hydro_block(blk)
    !if (allocated(blk%w)) deallocate(blk%w)
    !if (allocated(blk%u)) deallocate(blk%u)
    !if (allocated(blk%temp)) deallocate(blk%temp)
    !if (allocated(blk%egv)) deallocate(blk%egv)
    !if (allocated(blk%w_xl)) deallocate(blk%w_xl)
    !if (allocated(blk%w_xr)) deallocate(blk%w_xr)
    !if (allocated(blk%w_yl)) deallocate(blk%w_yl)
    !if (allocated(blk%w_yr)) deallocate(blk%w_yr)
    !if (allocated(blk%xslp)) deallocate(blk%xslp)
    !if (allocated(blk%yslp)) deallocate(blk%yslp)
    !if (allocated(blk%temp_xl)) deallocate(blk%temp_xl)
    !if (allocated(blk%temp_xr)) deallocate(blk%temp_xr)
    !if (allocated(blk%temp_yl)) deallocate(blk%temp_yl)
    !if (allocated(blk%temp_yr)) deallocate(blk%temp_yr)
    !if (allocated(blk%egv_xl)) deallocate(blk%egv_xl)
    !if (allocated(blk%egv_xr)) deallocate(blk%egv_xr)
    !if (allocated(blk%egv_yl)) deallocate(blk%egv_yl)
    !if (allocated(blk%egv_yr)) deallocate(blk%egv_yr)
    !if (allocated(blk%xflux)) deallocate(blk%xflux)
    !if (allocated(blk%yflux)) deallocate(blk%yflux)
    !if (allocated(blk%ekphi_mix_heat)) deallocate(blk%ekphi_mix_heat)
    !if (allocated(blk%hllc_vx)) deallocate(blk%hllc_vx)
    !if (allocated(blk%hllc_vy)) deallocate(blk%hllc_vy)
    !if (allocated(blk%u0)) deallocate(blk%u0)
    !if (allocated(blk%u1)) deallocate(blk%u1)
    !if (allocated(blk%cs)) deallocate(blk%cs)
    call deallocate_eos_block(blk)
    !if (allocated(blk%H2)) deallocate(blk%H2)
    !if (allocated(blk%HI)) deallocate(blk%HI)
    !if (allocated(blk%HII)) deallocate(blk%HII)
    !if (allocated(blk%HeI)) deallocate(blk%HeI)
    !if (allocated(blk%HeII)) deallocate(blk%HeII)
    !if (allocated(blk%HeIII)) deallocate(blk%HeIII)
    !if (allocated(blk%electron)) deallocate(blk%electron)
    call deallocate_radiation_block(blk)
    call deallocate_source_block(blk)
    !if (allocated(blk%erad_xedge)) deallocate(blk%erad_xedge)
    !if (allocated(blk%erad_yedge)) deallocate(blk%erad_yedge)
    !if (allocated(blk%tau)) deallocate(blk%tau)
    !if (allocated(blk%Fradx)) deallocate(blk%Fradx)
    !if (allocated(blk%Frady)) deallocate(blk%Frady)
    !if (allocated(blk%erad_xflux)) deallocate(blk%erad_xflux)
    !if (allocated(blk%erad_yflux)) deallocate(blk%erad_yflux)
    !if (allocated(blk%erad_adv)) deallocate(blk%erad_adv)
    !if (allocated(blk%fld_surfx)) deallocate(blk%fld_surfx)
    !if (allocated(blk%fld_surfy)) deallocate(blk%fld_surfy)
    !if (allocated(blk%entropy)) deallocate(blk%entropy)
    !if (allocated(blk%Erad_int)) deallocate(blk%Erad_int)
    !if (allocated(blk%Erad_source)) deallocate(blk%Erad_source)
    !if (allocated(blk%irrad)) deallocate(blk%irrad)
    !if (allocated(blk%kappa_planck)) deallocate(blk%kappa_planck)
    !if (allocated(blk%kappa_abs)) deallocate(blk%kappa_abs)
    !if (allocated(blk%kappa_rosseland)) deallocate(blk%kappa_rosseland)
    !if (allocated(blk%sigma_planck)) deallocate(blk%sigma_planck)
    !if (allocated(blk%sigma_abs)) deallocate(blk%sigma_abs)
    !if (allocated(blk%sigma_rosseland)) deallocate(blk%sigma_rosseland)
    if (allocated(blk%krx_inter)) deallocate(blk%krx_inter)
    if (allocated(blk%kry_inter)) deallocate(blk%kry_inter)
    if (allocated(blk%aradx)) deallocate(blk%aradx)
    if (allocated(blk%arady)) deallocate(blk%arady)
    !if (allocated(blk%irrad0)) deallocate(blk%irrad0)
    if (allocated(blk%dy)) deallocate(blk%dy)
    if (allocated(blk%dy2d)) deallocate(blk%dy2d)
    if (allocated(blk%dr)) deallocate(blk%dr)
    if (allocated(blk%xpflux)) deallocate(blk%xpflux)
    if (allocated(blk%ypflux)) deallocate(blk%ypflux)
    if (allocated(blk%kx)) deallocate(blk%kx)
    if (allocated(blk%ky)) deallocate(blk%ky)
    if (allocated(blk%torque_work)) deallocate(blk%torque_work)
    if (allocated(blk%viscous_heat_xz)) deallocate(blk%viscous_heat_xz)
    if (allocated(blk%omega_xz)) deallocate(blk%omega_xz)
    if (allocated(blk%torque_xz)) deallocate(blk%torque_xz)
    if (allocated(blk%viscous_heat_yz)) deallocate(blk%viscous_heat_yz)
    if (allocated(blk%omega_yz)) deallocate(blk%omega_yz)
    if (allocated(blk%torque_yz)) deallocate(blk%torque_yz)
    if (allocated(blk%surf1)) deallocate(blk%surf1)
    if (allocated(blk%surf2)) deallocate(blk%surf2)
    if (allocated(blk%lever_r)) deallocate(blk%lever_r)
    if (allocated(blk%lever_theta)) deallocate(blk%lever_theta)
    if (allocated(blk%source)) deallocate(blk%source)
    if (allocated(blk%vol)) deallocate(blk%vol)
    if (allocated(blk%rcylin)) deallocate(blk%rcylin)
    if (allocated(blk%gx1d)) deallocate(blk%gx1d)
    if (allocated(blk%alphax1d)) deallocate(blk%alphax1d)
    if (allocated(blk%gx)) deallocate(blk%gx)
    if (allocated(blk%gy)) deallocate(blk%gy)
    if (allocated(blk%alphax)) deallocate(blk%alphax)
    if (allocated(blk%alphay)) deallocate(blk%alphay)
    if (allocated(blk%e_guard)) deallocate(blk%e_guard)
    if (allocated(blk%s_guard)) deallocate(blk%s_guard)
    if (allocated(blk%w_guard)) deallocate(blk%w_guard)
    if (allocated(blk%n_guard)) deallocate(blk%n_guard)
    if (allocated(blk%omega)) deallocate(blk%omega)
    if (allocated(blk%omega0)) deallocate(blk%omega0)
    if (allocated(blk%omega1)) deallocate(blk%omega1)
    if (allocated(blk%omega_muscl)) deallocate(blk%omega_muscl)
    if (allocated(blk%amz)) deallocate(blk%amz)
    if (allocated(blk%fr_omega)) deallocate(blk%fr_omega)
    if (allocated(blk%ftheta_omega)) deallocate(blk%ftheta_omega)
    if (allocated(blk%l_omega)) deallocate(blk%l_omega)
    if (allocated(blk%gpotential)) deallocate(blk%gpotential)
    if (allocated(blk%viscous_heat)) deallocate(blk%viscous_heat)
    if (allocated(blk%viscous_heat_dt)) deallocate(blk%viscous_heat_dt)
    if (allocated(blk%radhydro_boost_viscous_heat)) deallocate(blk%radhydro_boost_viscous_heat)
    if (allocated(blk%heat1)) deallocate(blk%heat1)
    if (allocated(blk%heat2)) deallocate(blk%heat2)
    if (allocated(blk%cv)) deallocate(blk%cv)
    if (allocated(blk%omega_xl)) deallocate(blk%omega_xl)
    if (allocated(blk%omega_xr)) deallocate(blk%omega_xr)
    if (allocated(blk%omega_yl)) deallocate(blk%omega_yl)
    if (allocated(blk%omega_yr)) deallocate(blk%omega_yr)
    if (allocated(blk%s_geo1)) deallocate(blk%s_geo1)
    if (allocated(blk%s_geo2)) deallocate(blk%s_geo2)
    if (allocated(blk%s_grav)) deallocate(blk%s_grav)
    if (allocated(blk%Erad)) deallocate(blk%Erad)
    if (allocated(blk%stheta)) deallocate(blk%stheta)
    if (allocated(blk%ctheta)) deallocate(blk%ctheta)
    if (allocated(blk%passive_scalar)) deallocate(blk%passive_scalar)
    if (allocated(blk%passive_scalar0)) deallocate(blk%passive_scalar0)
    if (allocated(blk%x_center)) deallocate(blk%x_center)
    if (allocated(blk%y_center)) deallocate(blk%y_center)
    if (allocated(blk%x_interface)) deallocate(blk%x_interface)
    if (allocated(blk%y_interface)) deallocate(blk%y_interface)
    blk%on_processor=.false.
end subroutine release_block

subroutine relocate_llist_head_tail(load)
    !given a load, relocate the local list head and tail
    type(blockdef), pointer :: blk
    integer, dimension(:), allocatable :: load
    integer :: i,j
    blk=>blk_head
    do i=1,rank+1
        llist_head=>blk
        do j=1,load(i)-1
            blk=>blk%next
        end do
        llist_tail=>blk
        blk=>blk%next
    end do
end subroutine relocate_llist_head_tail

subroutine renumber_domain_blocks()
    type(blockdef), pointer :: blk
    integer :: id
    id=1
    blk=>blk_head
    do while (associated(blk))
        blk_tail=>blk
        blk%blk_id=id
        id=id+1
        blk=>blk%next
    end do
    nblk_total=blk_tail%blk_id
    nullify(blk)
end subroutine renumber_domain_blocks

subroutine destroy_block(blk)
    type(blockdef), pointer :: blk
    call release_block(blk)
    deallocate(blk)
end subroutine destroy_block

subroutine allocate_passive_block(blk)
    type(blockdef), pointer :: blk
    call allocate_cell_data_block(blk%passive_scalar,npassive)
    call allocate_cell_data_block(blk%passive_scalar0,npassive)
    call allocate_xsurface_data_block(blk%xpflux,npassive)
end subroutine allocate_passive_block

function spherical_r_center(r1,r2)
    real(8) :: r1,r2,spherical_r_center,rm,dr
    rm=(r1+r2)/2d0
    dr=r2-r1
    spherical_r_center=rm+(2d0*rm*dr**2d0)/(12d0*rm**2d0+dr**2d0)
end function spherical_r_center

function spherical_theta_center(theta1,theta2)
    real(8) :: theta1,theta2,spherical_theta_center
    spherical_theta_center=(theta1*cos(theta1)-sin(theta1)-theta2*cos(theta2)+sin(theta2))/(cos(theta1)-cos(theta2))
end function spherical_theta_center

function spherical_fr_omega(r1,r2,theta1,theta2)
    real(8) :: spherical_fr_omega,r1,r2,theta1,theta2
    spherical_fr_omega=(r2**4d0-r1**4d0)/4d0*(cos(theta2)**3d0/3d0-cos(theta1)**3d0/3d0-cos(theta2)+cos(theta1))
end function spherical_fr_omega

function spherical_ftheta_omega(r1,r2,theta1,theta2)
    real(8) :: spherical_ftheta_omega,r1,r2,theta1,theta2
    spherical_ftheta_omega=(r2**4d0-r1**4d0)/4d0*(sin(theta2)**3d0/3d0-sin(theta1)**3d0/3d0)
end function spherical_ftheta_omega

function spherical_l_omega(r1,r2,theta1,theta2)
    real(8) :: spherical_l_omega,r1,r2,theta1,theta2
    spherical_l_omega=(r2**5d0-r1**5d0)/5d0*(cos(theta2)**3d0/3d0-cos(theta1)**3d0/3d0-cos(theta2)+cos(theta1))
end function spherical_l_omega

function geometric_a(s0,sn,h)
    real(8) :: geometric_a,s0,sn,h
    geometric_a=(sn-s0)/(h-1d0)
end function geometric_a

function geometric_b(s0,sn,h)
    real(8) :: geometric_b,s0,sn,h
    geometric_b=(h*s0-sn)/(h-1d0)
end function geometric_b

function geometric_pn(s0,sn,b)
    real(8) :: geometric_pn,s0,sn,b
    geometric_pn=log((sn-b)/(s0-b))
end function geometric_pn

function geometric_si(a,b,dp,level,iblkx,i)
    real(8) :: geometric_si,a,b,dp
    integer :: iblkx,i,level
    geometric_si=a*exp(((iblkx-1)*blk_size_nx+i)*dp/(2d0**level))+b
end function geometric_si

subroutine initialize_geometry()
    if (llnx) then
        if (xgeo_h==1) then
            if (rank==0) print *,'xgeo_h error'
        end if
        xgeo_a=geometric_a(n_domain(1),n_domain(2),xgeo_h)
        xgeo_b=geometric_b(n_domain(1),n_domain(2),xgeo_h)
        xgeo_pn=geometric_pn(n_domain(1),n_domain(2),xgeo_b)
        xgeo_dp=xgeo_pn/nx
    else
    end if
    if (llny) then
        if (ygeo_h==1) then
            if (rank==0) print *,'ygeo_h error'
        end if
        ygeo_a=geometric_a(n_domain(3),n_domain(4),ygeo_h)
        ygeo_b=geometric_b(n_domain(3),n_domain(4),ygeo_h)
        ygeo_pn=geometric_pn(n_domain(3),n_domain(4),ygeo_b)
        ygeo_dp=ygeo_pn/ny
    end if
end subroutine initialize_geometry

subroutine allocate_geometry_block(blk)
    !surf1 is the surface area (line length in 2d) that is perpendicular to the 1st coordinate, same for surf2
    type(blockdef), pointer :: blk
    real(8) :: coords(3),dx,dxyz(3),r,theta,theta1,theta2,r1,r2,v1,v2,v,pos(3),r_mean,rad,vol,dtheta0,q
    integer :: i,j,ijk(3),k
    call allocate_cell_data_block(blk%vol)
    allocate(blk%x_center(blk_xlb:blk_xub),blk%x_interface(blk_xlb-1:blk_xub),blk%surf1(blk_xlb-1:blk_xub,1,1))
    if (llnx) allocate(blk%dr(blk_xlb:blk_xub))
    coords=blk%pos
    dx=blk%dxyz(1)
    if (igeometry==0) then
        do i=blk_xlb,blk_xub
            blk%x_center(i)=coords(1)+(i-half)*dx
            blk%vol(i,1,1)=dx
        end do
        do i=blk_xlb-1,blk_xub
            blk%surf1=1d0
            blk%x_interface(i)=coords(1)+i*dx
        end do
    else if (igeometry==2) then
        if (llnx) then
            do i=blk_xlb-1,blk_xub
                blk%x_interface(i)=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level,blk%key(1),i)
            end do
            !print *,rank,blk%x_interface/n_domain(1)
        else
            blk%x_interface=coords(1)+(/((i*dx),i=blk_xlb-1,blk_xub)/)
        end if
        do i=blk_xlb-1,blk_xub
            r=blk%x_interface(i)
            blk%surf1(i,1,1)=r*r
        end do
        do i=blk_xlb,blk_xub
            r1=blk%x_interface(i-1)
            r2=blk%x_interface(i)
            blk%x_center(i)=spherical_r_center(r1,r2)
            blk%vol(i,1,1)=(r2**3d0-r1**3d0)/3d0
        end do
        if (llnx) blk%dr=blk%x_interface(blk_xlb:blk_xub)-blk%x_interface(blk_xlb-1:blk_xub-1)
    end if
end subroutine allocate_geometry_block

subroutine deallocate_geometry_block(blk)
    type(blockdef), pointer :: blk
end subroutine deallocate_geometry_block

subroutine allocate_implicit_mesh(blk)
    type(blockdef), pointer :: blk
    real(8), dimension(2) :: gxl,gxu,gyl,gyu,alpha_xl,alpha_xu,alpha_yl,alpha_yu
    real(8) :: gxl1d,gxu1d,alpha_xl1d,alpha_xu1d
    integer :: i,j
    if (iradiation==4) then
        if (allocated(blk%gx1d)) deallocate(blk%gx1d)
        if (allocated(blk%alphax1d)) deallocate(blk%alphax1d)
        allocate(blk%gx1d(2,blk_size_nx),blk%alphax1d(2,blk_size_nx))
        call calculate_guard_coord(blk)
        do i=1,blk_size_nx
            call calculate_g_alpha_1d(blk,1d0,i,gxl1d,gxu1d,alpha_xl1d,alpha_xu1d)
            blk%gx1d(1,i)=gxl1d
            blk%gx1d(2,i)=gxu1d
            blk%alphax1d(1,i)=alpha_xl1d
            blk%alphax1d(2,i)=alpha_xu1d
        end do
    end if
end subroutine allocate_implicit_mesh

subroutine calculate_guard_coord(blk)
    type(blockdef), pointer :: blk
    call calculate_xl_guard_coord(blk)
    call calculate_xu_guard_coord(blk)
end subroutine calculate_guard_coord

subroutine calculate_xl_guard_coord(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,x,dx
    integer :: nb_level(2),iblkx
    nb_level=blk%nb_level_1d
    dx=blk%dxyz(1)
    if (llnx) then
        iblkx=int((blk%key(1)-1)*2d0**nb_level(1)+1)
        x1=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(1),iblkx,-1)
        x2=blk%x_interface(0)
    else
        x1=blk%x_interface(0)-dx*2d0**(-nb_level(1))
        x2=blk%x_interface(0)
    end if
    x=spherical_r_center(x1,x2)
    blk%xl_guard=x
    !print *,'xl',x,blk%x_center(0),x1,blk%x_interface(-1),blk%key(1),iblkx,nb_level(1)
end subroutine calculate_xl_guard_coord

subroutine calculate_xu_guard_coord(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,x,dx
    integer :: nb_level(2),iblkx
    nb_level=blk%nb_level_1d
    dx=blk%dxyz(1)
    if (llnx) then
        iblkx=int(blk%key(1)*2d0**nb_level(2))
        x1=blk%x_interface(blk_size_nx)
        x2=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(2),iblkx,blk_size_nx+1)
    else
        x1=blk%x_interface(blk_size_nx)
        x2=blk%x_interface(blk_size_nx)+dx*2d0**(-nb_level(2))
    end if
    x=spherical_r_center(x1,x2)
    blk%xu_guard=x
    !print *,'xu',x,blk%x_center(blk_size_nx+1),x2,blk%x_interface(blk_size_nx+1),blk%key(1),iblkx,nb_level(2)
end subroutine calculate_xu_guard_coord

subroutine calculate_e_guard_coord(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,x,dx,dx_nb,y1,y2,y,dy,dy_nb,dtheta0,dtheta
    integer :: i,j,nb_level(8),iblkx,iblky
    nb_level=blk%nb_level
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==2) then
        if (llnx) then
            iblkx=int(blk%key(1)*2d0**nb_level(1))
            x1=blk%x_interface(blk_size_nx)
            x2=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(1),iblkx,blk_size_nx+1)
        else
            x1=blk%x_interface(blk_size_nx)
            x2=blk%x_interface(blk_size_nx)+dx*2d0**(-nb_level(1))
        end if
        x=spherical_r_center(x1,x2)
        blk%e_guardx=x
        if (nb_level(1)==0) then
            blk%e_guard(1:blk_size_ny)=blk%y_center(1:blk_size_ny)
        else if (nb_level(1)==1) then
            if (llny) then
                iblky=int((blk%key(2)-1)*2d0**nb_level(1)+1)
                do j=1,int(blk_size_ny*2d0**nb_level(1))
                    y1=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+nb_level(1),iblky,j-1)
                    y2=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+nb_level(1),iblky,j)
                    blk%e_guard(j)=spherical_theta_center(y1,y2)
                end do
            else
                dy_nb=dy*2d0**(-nb_level(1))
                do j=1,blk_size_ny*2
                    y1=blk%y_interface(0)+dy_nb*(j-1)
                    y2=blk%y_interface(0)+dy_nb*j
                    blk%e_guard(j)=spherical_theta_center(y1,y2)
                end do
            end if
        else if (nb_level(1)==-1) then
            do j=1,blk_size_ny/2
                y1=blk%y_interface(2*j-2)
                y2=blk%y_interface(2*j)
                blk%e_guard(j)=spherical_theta_center(y1,y2)
            end do
        end if
    else if (igeometry==0) then
        dx_nb=dx*2d0**(-nb_level(1))
        blk%e_guardx=blk%x_interface(blk_size_nx)+dx_nb/2d0
        if (nb_level(1)==0) then
            blk%e_guard(1:blk_size_ny)=blk%y_center(1:blk_size_ny)
        else if (nb_level(1)==1) then
            dy_nb=dy/2d0
            do j=1,blk_size_ny*2
                blk%e_guard(j)=blk%y_interface(0)+dy_nb*(j-0.5d0)
            end do
        else if (nb_level(1)==-1) then
            do j=1,blk_size_ny/2
                blk%e_guard(j)=blk%y_interface(2*j-1)
            end do
        end if
    end if
end subroutine calculate_e_guard_coord

subroutine calculate_s_guard_coord(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,x,dx,dx_nb,y1,y2,y,dy,dy_nb,dtheta0,dtheta
    integer :: i,j,nb_level(8),iblkx,iblky
    nb_level=blk%nb_level
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==2) then
        if (llny) then
            iblky=int((blk%key(2)-1)*2d0**nb_level(2)+1)
            y1=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+nb_level(2),iblky,-1)
            y2=blk%y_interface(0)
        else
            dtheta=dy*2d0**(-nb_level(2))
            y1=blk%y_interface(0)-dtheta
            y2=blk%y_interface(0)
        end if
        y=spherical_theta_center(y1,y2)
        blk%s_guardy=y
        if (llnx) then
            iblkx=int((blk%key(1)-1)*2d0**nb_level(2)+1)
            do i=1,int(blk_size_nx*2d0**nb_level(2))
                x1=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(2),iblkx,i-1)
                x2=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(2),iblkx,i)
                blk%s_guard(i)=spherical_r_center(x1,x2)
            end do
        else
            dx_nb=dx*2d0**(-nb_level(2))
            do i=1,int(blk_size_nx*2d0**nb_level(2))
                x1=blk%x_interface(0)+dx_nb*(i-1)
                x2=blk%x_interface(0)+dx_nb*i
                blk%s_guard(i)=spherical_r_center(x1,x2)
            end do
        end if
    else if (igeometry==0) then
        dy_nb=dy*2d0**(-nb_level(2))
        blk%s_guardy=blk%y_interface(0)-dy_nb/2d0
        if (nb_level(2)==0) then
            blk%s_guard(1:blk_size_nx)=blk%x_center(1:blk_size_nx)
        else if (nb_level(2)==1) then
            dx_nb=dx/2d0
            do i=1,blk_size_nx*2
                blk%s_guard(i)=blk%x_interface(0)+dx_nb*(i-0.5d0)
            end do
        else if (nb_level(2)==-1) then
            do i=1,blk_size_nx/2
                blk%s_guard(i)=blk%x_interface(2*i-1)
            end do
        end if
    end if
end subroutine calculate_s_guard_coord

subroutine calculate_w_guard_coord(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,x,dx,dx_nb,y1,y2,y,dy,dy_nb,dtheta0,dtheta
    integer :: i,j,nb_level(8),iblkx,iblky
    nb_level=blk%nb_level
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==2) then
        if (llnx) then
            iblkx=int((blk%key(1)-1)*2d0**nb_level(3)+1)
            x1=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(3),iblkx,-1)
            x2=blk%x_interface(0)
        else
            x1=blk%x_interface(0)-dx*2d0**(-nb_level(3))
            x2=blk%x_interface(0)
        end if
        x=spherical_r_center(x1,x2)
        blk%w_guardx=x
        if (nb_level(3)==0) then
            blk%w_guard(1:blk_size_ny)=blk%y_center(1:blk_size_ny)
        else if (nb_level(3)==1) then
            if (llny) then
                iblky=int((blk%key(2)-1)*2d0**nb_level(3)+1)
                do j=1,int(blk_size_ny*2d0**nb_level(3))
                    y1=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+nb_level(3),iblky,j-1)
                    y2=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+nb_level(3),iblky,j)
                    blk%w_guard(j)=spherical_theta_center(y1,y2)
                end do
            else
                dy_nb=dy*2d0**(-nb_level(3))
                do j=1,blk_size_ny*2
                    y1=blk%y_interface(0)+dy_nb*(j-1)
                    y2=blk%y_interface(0)+dy_nb*j
                    blk%w_guard(j)=spherical_theta_center(y1,y2)
                end do
            end if
        else if (nb_level(3)==-1) then
            do j=1,blk_size_ny/2
                y1=blk%y_interface(2*j-2)
                y2=blk%y_interface(2*j)
                blk%w_guard(j)=spherical_theta_center(y1,y2)
            end do
        end if
    else if (igeometry==0) then
        dx_nb=dx*2d0**(-nb_level(3))
        blk%w_guardx=blk%x_interface(0)-dx_nb/2d0
        if (nb_level(3)==0) then
            blk%w_guard(1:blk_size_ny)=blk%y_center(1:blk_size_ny)
        else if (nb_level(3)==1) then
            dy_nb=dy/2d0
            do j=1,blk_size_ny*2
                blk%w_guard(j)=blk%y_interface(0)+dy_nb*(j-0.5d0)
            end do
        else if (nb_level(3)==-1) then
            do j=1,blk_size_ny/2
                blk%w_guard(j)=blk%y_interface(2*j-1)
            end do
        end if
    end if
end subroutine calculate_w_guard_coord

subroutine calculate_n_guard_coord(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,x,dx,dx_nb,y1,y2,y,dy,dy_nb,dtheta0,dtheta
    integer :: i,j,nb_level(8),iblkx,iblky
    nb_level=blk%nb_level
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==2) then
        if (llny) then
            iblky=int(blk%key(2)*2d0**nb_level(4))
            y1=blk%y_interface(blk_size_ny)
            y2=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+nb_level(4),iblky,blk_size_ny+1)
        else
            dtheta=dy*2d0**(-nb_level(4))
            y1=blk%y_interface(blk_size_ny)
            y2=blk%y_interface(blk_size_ny)+dtheta
        end if
        y=spherical_theta_center(y1,y2)
        blk%n_guardy=y
        if (llnx) then
            iblkx=int((blk%key(1)-1)*2d0**nb_level(4)+1)
            do i=1,int(blk_size_nx*2d0**nb_level(4))
                x1=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(4),iblkx,i-1)
                x2=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(4),iblkx,i)
                blk%n_guard(i)=spherical_r_center(x1,x2)
            end do
        else
            dx_nb=dx*2d0**(-nb_level(4))
            do i=1,int(blk_size_nx*2d0**nb_level(4))
                x1=blk%x_interface(0)+dx_nb*(i-1)
                x2=blk%x_interface(0)+dx_nb*i
                blk%n_guard(i)=spherical_r_center(x1,x2)
            end do
        end if
    else if (igeometry==0) then
        dy_nb=dy*2d0**(-nb_level(4))
        blk%n_guardy=blk%y_interface(blk_size_ny)+dy_nb/2d0
        if (nb_level(4)==0) then
            blk%n_guard(1:blk_size_nx)=blk%x_center(1:blk_size_nx)
        else if (nb_level(4)==1) then
            dx_nb=dx/2d0
            do i=1,blk_size_nx*2
                blk%n_guard(i)=blk%x_interface(0)+dx_nb*(i-0.5d0)
            end do
        else if (nb_level(4)==-1) then
            do i=1,blk_size_nx/2
                blk%n_guard(i)=blk%x_interface(2*i-1)
            end do
        end if
    end if
end subroutine calculate_n_guard_coord

subroutine calculate_se_guard_coord(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,x,dx,dx_nb,y1,y2,y,dy,dy_nb,dtheta0,dtheta
    integer :: i,j,nb_level(8),iblkx,iblky
    nb_level=blk%nb_level
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==2) then
        if (llnx) then
            iblkx=int(blk%key(1)*2d0**nb_level(5))
            x1=blk%x_interface(blk_size_nx)
            x2=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(5),iblkx,blk_size_nx+1)
        else
            x1=blk%x_interface(blk_size_nx)
            x2=blk%x_interface(blk_size_nx)+dx*2d0**(-nb_level(5))
        end if
        blk%se_guard(1)=spherical_r_center(x1,x2)
        if (llny) then
            iblky=int((blk%key(2)-1)*2d0**nb_level(5)+1)
            y1=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+nb_level(5),iblky,-1)
            y2=blk%y_interface(0)
        else
            dtheta=dy*2d0**(-nb_level(5))
            y1=blk%y_interface(0)-dtheta
            y2=blk%y_interface(0)
        end if
        blk%se_guard(2)=spherical_theta_center(y1,y2)
    else if (igeometry==0) then
        blk%se_guard(1)=blk%x_interface(blk_size_nx)+dx*2d0**(-nb_level(5))/2d0
        blk%se_guard(2)=blk%y_interface(0)-dy*2d0**(-nb_level(5))/2d0
    end if
end subroutine calculate_se_guard_coord

subroutine calculate_sw_guard_coord(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,x,dx,dx_nb,y1,y2,y,dy,dy_nb,dtheta0,dtheta
    integer :: i,j,nb_level(8),iblkx,iblky
    nb_level=blk%nb_level
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==2) then
        if (llnx) then
            iblkx=int((blk%key(1)-1)*2d0**nb_level(6)+1)
            x1=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(6),iblkx,-1)
            x2=blk%x_interface(0)
        else
            x1=blk%x_interface(0)-dx*2d0**(-nb_level(6))
            x2=blk%x_interface(0)
        end if
        blk%sw_guard(1)=spherical_r_center(x1,x2)
        if (llny) then
            iblky=int((blk%key(2)-1)*2d0**nb_level(6)+1)
            y1=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+nb_level(6),iblky,-1)
            y2=blk%y_interface(0)
        else
            dtheta=dy*2d0**(-nb_level(6))
            y1=blk%y_interface(0)-dtheta
            y2=blk%y_interface(0)
        end if
        blk%sw_guard(2)=spherical_theta_center(y1,y2)
    else if (igeometry==0) then
        blk%sw_guard(1)=blk%x_interface(0)-dx*2d0**(-nb_level(6))/2d0
        blk%sw_guard(2)=blk%y_interface(0)-dy*2d0**(-nb_level(6))/2d0
    end if
end subroutine calculate_sw_guard_coord

subroutine calculate_nw_guard_coord(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,x,dx,dx_nb,y1,y2,y,dy,dy_nb,dtheta0,dtheta
    integer :: i,j,nb_level(8),iblkx,iblky
    nb_level=blk%nb_level
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==2) then
        if (llnx) then
            iblkx=int((blk%key(1)-1)*2d0**nb_level(7)+1)
            x1=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(7),iblkx,-1)
            x2=blk%x_interface(0)
        else
            x1=blk%x_interface(0)-dx*2d0**(-nb_level(7))
            x2=blk%x_interface(0)
        end if
        blk%sw_guard(1)=spherical_r_center(x1,x2)
        if (llny) then
            iblky=int(blk%key(2)*2d0**nb_level(7))
            y1=blk%y_interface(blk_size_ny)
            y2=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+nb_level(7),iblky,blk_size_ny+1)
        else
            dtheta=dy*2d0**(-nb_level(7))
            y1=blk%y_interface(blk_size_ny)
            y2=blk%y_interface(blk_size_ny)+dtheta
        end if
        blk%nw_guard(2)=spherical_theta_center(y1,y2)
    else if (igeometry==0) then
        blk%nw_guard(1)=blk%x_interface(0)-dx*2d0**(-nb_level(7))/2d0
        blk%nw_guard(2)=blk%y_interface(blk_size_ny)+dy*2d0**(-nb_level(7))/2d0
    end if
end subroutine calculate_nw_guard_coord

subroutine calculate_ne_guard_coord(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,x,dx,dx_nb,y1,y2,y,dy,dy_nb,dtheta0,dtheta
    integer :: i,j,nb_level(8),iblkx,iblky
    nb_level=blk%nb_level
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==2) then
        if (llnx) then
            iblkx=int(blk%key(1)*2d0**nb_level(8))
            x1=blk%x_interface(blk_size_nx)
            x2=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+nb_level(8),iblkx,blk_size_nx+1)
        else
            x1=blk%x_interface(blk_size_nx)
            x2=blk%x_interface(blk_size_nx)+dx*2d0**(-nb_level(8))
        end if
        blk%ne_guard(1)=spherical_r_center(x1,x2)
        if (llny) then
            iblky=int(blk%key(2)*2d0**nb_level(8))
            y1=blk%y_interface(blk_size_ny)
            y2=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+nb_level(8),iblky,blk_size_ny+1)
        else
            dtheta=dy*2d0**(-nb_level(8))
            y1=blk%y_interface(blk_size_ny)
            y2=blk%y_interface(blk_size_ny)+dtheta
        end if
        blk%ne_guard(2)=spherical_theta_center(y1,y2)
    else if (igeometry==0) then
        blk%ne_guard(1)=blk%x_interface(blk_size_nx)+dx*2d0**(-nb_level(8))/2d0
        blk%ne_guard(2)=blk%y_interface(blk_size_ny)+dy*2d0**(-nb_level(8))/2d0
    end if
end subroutine calculate_ne_guard_coord

subroutine recalculate_g_alpha(blk)
    type(blockdef), pointer :: blk
    real(8), dimension(2) :: sxl,sxu,syl,syu,gxl,gxu,gyl,gyu,alpha_xl,alpha_xu,alpha_yl,alpha_yu
    real(8) :: gxl1d,gxu1d,alphaxl1d,alphaxu1d
    integer :: i,j,key(3)
    call xl_g_alpha(blk,1d0,1,gxl1d,alphaxl1d)
    call xu_g_alpha(blk,1d0,blk_size_nx,gxu1d,alphaxu1d)
    blk%gx1d(1,1)=gxl1d
    blk%gx1d(2,blk_size_nx)=gxu1d
    blk%alphax1d(1,1)=alphaxl1d
    blk%alphax1d(2,blk_size_nx)=alphaxu1d
end subroutine recalculate_g_alpha

subroutine calculate_g_alpha_1d(blk,dt,i,gxl,gxu,alpha_xl,alpha_xu)
    type(blockdef), pointer :: blk
    real(8) :: dt,gxl,gxu,alpha_xl,alpha_xu
    integer :: i
    call xl_g_alpha(blk,dt,i,gxl,alpha_xl)
    call xu_g_alpha(blk,dt,i,gxu,alpha_xu)
end subroutine calculate_g_alpha_1d

subroutine xl_g_alpha(blk,dt,i,g,alpha)
    type(blockdef), pointer :: blk
    real(8) :: dt,g,alpha,dx
    integer :: i,nb_level(2)
    nb_level=blk%nb_level_1d
    dx=blk%dxyz(1)
    if (igeometry==0) then
        g=1d0/dx
        if (i==1) then
            alpha=dt/(dx*(0.5d0+0.5d0*2d0**(-nb_level(1))))
        else
            alpha=dt/dx
        end if
    else if (igeometry==2) then
        g=blk%surf1(i-1,1,1)/blk%vol(i,1,1)
        if (i==1) then
            alpha=dt/(blk%x_center(i)-blk%xl_guard)
        else
            alpha=dt/(blk%x_center(i)-blk%x_center(i-1))
        end if
    end if
end subroutine xl_g_alpha

subroutine xu_g_alpha(blk,dt,i,g,alpha)
    type(blockdef), pointer :: blk
    real(8) :: dt,g,alpha,dx
    integer :: i,nb_level(2)
    nb_level=blk%nb_level_1d
    dx=blk%dxyz(1)
    if (igeometry==0) then
        g=1d0/dx
        if (i==blk_size_nx) then
            alpha=dt/(dx*(0.5d0+0.5d0*2d0**(-nb_level(2))))
        else
            alpha=dt/dx
        end if
    else if (igeometry==2) then
        g=blk%surf1(i,1,1)/blk%vol(i,1,1)
        if (i==blk_size_nx) then
            alpha=dt/(blk%xu_guard-blk%x_center(i))
        else
            alpha=dt/(blk%x_center(i+1)-blk%x_center(i))
        end if
    end if
end subroutine xu_g_alpha

subroutine calculate_g_alpha(blk,dt,i,j,gxl,gxu,gyl,gyu,alpha_xl,alpha_xu,alpha_yl,alpha_yu)
    type(blockdef), pointer :: blk
    real(8) :: dt
    real(8), dimension(2) :: gxl,gxu,gyl,gyu,alpha_xl,alpha_xu,alpha_yl,alpha_yu
    integer :: i,j
    call east_g_alpha(blk,dt,i,j,gxu,alpha_xu)
    call south_g_alpha(blk,dt,i,j,gyl,alpha_yl)
    call west_g_alpha(blk,dt,i,j,gxl,alpha_xl)
    call north_g_alpha(blk,dt,i,j,gyu,alpha_yu)
end subroutine calculate_g_alpha

subroutine east_g_alpha(blk,dt,i,j,gxu,alpha_xu)
    type(blockdef), pointer :: blk
    integer :: i,j,iblky
    real(8) :: dt,sxu(2),gxu(2),alpha_xu(2),rr(2),r1,r2,theta1,theta2,thetam,dtheta,dx,dy,vol
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==0) then
        if (i==blk_size_nx.and.blk%nb_level(1)/=0) then
            if (blk%nb_level(1)==1) then
                gxu=1d0/dx/2d0;alpha_xu=dt/(dx*0.75d0)
            else if (blk%nb_level(1)==-1) then
                gxu=1d0/dx;alpha_xu=dt/(dx*1.5d0)
            end if
        else
            gxu=1d0/dx;alpha_xu=dt/dx
        end if
    else if (igeometry==2) then
        if (i==blk_size_nx) then
            if (blk%nb_level(1)==1) then
                r2=blk%x_interface(i)
                theta1=blk%y_interface(j-1)
                theta2=blk%y_interface(j)
                if (llny) then
                    !dtheta=blk%dtheta0*dy**(j-1)
                    !thetam=theta1+dtheta/(1+sqrt(dy))
                    iblky=(blk%key(2)-1)*2+1
                    thetam=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+1,iblky,2*j-1)
                else
                    thetam=(theta1+theta2)/2
                end if
                sxu(1)=r2**2d0*(cos(theta1)-cos(thetam))
                sxu(2)=r2**2d0*(cos(thetam)-cos(theta2))
                vol=blk%vol(i,j,1)
                gxu=sxu/vol
            else
                gxu=blk%surf1(i,j,1)/blk%vol(i,j,1)
            end if
            alpha_xu=dt/(blk%e_guardx-blk%x_center(i))
        else
            gxu=blk%surf1(i,j,1)/blk%vol(i,j,1)
            alpha_xu=dt/(blk%x_center(i+1)-blk%x_center(i))
        end if
    end if
end subroutine east_g_alpha

subroutine south_g_alpha(blk,dt,i,j,gyl,alpha_yl)
    type(blockdef), pointer :: blk
    integer :: i,j,iblkx
    real(8) :: dt,syl(2),gyl(2),alpha_yl(2),r,r1,r2,rm,theta1,theta2,theta,dtheta0,dx,dy,dtheta,vol
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==0) then
        if (j==1.and.blk%nb_level(2)/=0) then
            if (blk%nb_level(2)==1) then
                gyl=1d0/dy/2d0;alpha_yl=dt/(dy*0.75d0)
            else if (blk%nb_level(2)==-1) then
                gyl=1d0/dy;alpha_yl=dt/(dy*1.5d0)
            end if
        else
            gyl=1d0/dy;alpha_yl=dt/dy
        end if
    else if (igeometry==2) then
        if (j==1) then
            if (blk%nb_level(2)==1) then
                r1=blk%x_interface(i-1)
                r2=blk%x_interface(i)
                if (llnx) then
                    iblkx=(blk%key(1)-1)*2+1
                    rm=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+1,iblkx,2*i-1)
                else
                    rm=(r1+r2)/2
                end if
                theta2=blk%y_interface(j-1)
                syl(1)=0.5d0*(rm**2d0-r1**2d0)*sin(theta2)
                syl(2)=0.5d0*(r2**2d0-rm**2d0)*sin(theta2)
                vol=blk%vol(i,j,1)
                gyl=syl/vol
            else
                gyl=blk%surf2(i,j-1,1)/blk%vol(i,j,1)
            end if
            r=blk%x_center(i)
            alpha_yl=dt/r/(blk%y_center(j)-blk%s_guardy)
        else
            gyl=blk%surf2(i,j-1,1)/blk%vol(i,j,1)
            r=blk%x_center(i)
            alpha_yl=dt/r/(blk%y_center(j)-blk%y_center(j-1))
        end if
    end if
end subroutine south_g_alpha

subroutine west_g_alpha(blk,dt,i,j,gxl,alpha_xl)
    type(blockdef), pointer :: blk
    integer :: i,j,iblky,key(3)
    real(8) :: dt,sxl(2),gxl(2),alpha_xl(2),rl(2),r1,r2,theta1,theta2,thetam,dtheta,dx,dy,vol
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==0) then
        if (i==1.and.blk%nb_level(3)/=0) then
            if (blk%nb_level(3)==1) then
                gxl=1d0/dx/2d0;alpha_xl=dt/(dx*0.75d0)
            else if (blk%nb_level(3)==-1) then
                gxl=1d0/dx;alpha_xl=dt/(dx*1.5d0)
            end if
        else
            gxl=1d0/dx;alpha_xl=dt/dx
        end if
    else if (igeometry==2) then
        if (i==1) then
            if (blk%nb_level(3)==1) then
                r1=blk%x_interface(i-1)
                theta1=blk%y_interface(j-1)
                theta2=blk%y_interface(j)
                if (llny) then
                    iblky=(blk%key(2)-1)*2+1
                    thetam=geometric_si(ygeo_a,ygeo_b,ygeo_dp,blk%level+1,iblky,2*j-1)
                    !dtheta=blk%dtheta0*dy**(j-1)
                    !thetam=theta1+dtheta/(1+sqrt(dy))
                else
                    thetam=(theta1+theta2)/2
                end if
                sxl(1)=r1**2d0*(cos(theta1)-cos(thetam))
                sxl(2)=r1**2d0*(cos(thetam)-cos(theta2))
                vol=blk%vol(i,j,1)
                gxl=sxl/vol
            else
                gxl=blk%surf1(i-1,j,1)/blk%vol(i,j,1)
            end if
            alpha_xl=dt/(blk%x_center(i)-blk%w_guardx)
            key=blk%key
        else
            gxl=blk%surf1(i-1,j,1)/blk%vol(i,j,1)
            alpha_xl=dt/(blk%x_center(i)-blk%x_center(i-1))
        end if
    end if
end subroutine west_g_alpha

subroutine north_g_alpha(blk,dt,i,j,gyu,alpha_yu)
    type(blockdef), pointer :: blk
    integer :: i,j,iblkx
    real(8) :: dt,syu(2),gyu(2),alpha_yu(2),r,r1,r2,rm,theta1,theta2,theta,dtheta0,dx,dy,dtheta,vol
    dx=blk%dxyz(1);dy=blk%dxyz(2)
    if (igeometry==0) then
        if (j==blk_size_ny.and.blk%nb_level(4)/=0) then
            if (blk%nb_level(4)==1) then
                gyu=1d0/dy/2d0;alpha_yu=dt/(dy*0.75d0)
            else if (blk%nb_level(4)==-1) then
                gyu=1d0/dy;alpha_yu=dt/(dy*1.5d0)
            end if
        else
            gyu=1d0/dy;alpha_yu=dt/dy
        end if
    else if (igeometry==2) then
        if (j==blk_size_ny) then
            if (blk%nb_level(4)==1) then
                r1=blk%x_interface(i-1)
                r2=blk%x_interface(i)
                if (llnx) then
                    iblkx=(blk%key(1)-1)*2+1
                    rm=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+1,iblkx,2*i-1)
                else
                    rm=(r1+r2)/2
                end if
                theta1=blk%y_interface(j)
                syu(1)=0.5d0*(rm**2d0-r1**2d0)*sin(theta1)
                syu(2)=0.5d0*(r2**2d0-rm**2d0)*sin(theta1)
                vol=blk%vol(i,j,1)
                gyu=syu/vol
            else
                gyu=blk%surf2(i,j,1)/blk%vol(i,j,1)
            end if
            r=blk%x_center(i)
            alpha_yu=dt/r/(blk%n_guardy-blk%y_center(j))
        else
            gyu=blk%surf2(i,j,1)/blk%vol(i,j,1)
            r=blk%x_center(i)
            alpha_yu=dt/r/(blk%y_center(j+1)-blk%y_center(j))
        end if
    end if
end subroutine north_g_alpha

subroutine calculate_r_center(r,dxyz,r_val)
    !r is the non-weighted center, r_val is the weighted center
    real(8) :: r,dxyz(3),r_val
    r_val=r+(2d0*r*dxyz(1)**2)/(12d0*r**2+dxyz(1)**2)
end subroutine calculate_r_center

subroutine allocate_hydro_block(blk)
    type(blockdef), pointer :: blk
    call allocate_cell_data_block(blk%w,5)
    call allocate_cell_data_block(blk%u,5)
    call allocate_cell_data_block(blk%temp)
    call allocate_cell_data_block(blk%egv)
    call allocate_cell_data_block(blk%w_xl,5)
    call allocate_cell_data_block(blk%w_xr,5)
    call allocate_cell_data_block(blk%xslp,5)
    call allocate_xsurface_data_block(blk%xflux,5)
    call allocate_cell_data_block(blk%u_muscl,5)
    call allocate_cell_data_block(blk%u0,5)
    call allocate_cell_data_block(blk%u1,5)
    call allocate_xsurface_data_block(blk%hllc_vx)
#if     ieos==2
    call allocate_cell_data_block(blk%cs)
    call allocate_cell_data_block(blk%temp_xl)
    call allocate_cell_data_block(blk%temp_xr)
    call allocate_cell_data_block(blk%egv_xl)
    call allocate_cell_data_block(blk%egv_xr)
#endif
end subroutine allocate_hydro_block

subroutine deallocate_hydro_block(blk)
    type(blockdef), pointer :: blk
    deallocate(blk%w,blk%u,blk%u0,blk%u1,blk%u_muscl,blk%temp,blk%egv,blk%w_xl,blk%w_xr,blk%xslp,blk%xflux,blk%hllc_vx)
#if     ieos==2
    deallocate(blk%cs,blk%temp_xl,blk%temp_xr,blk%egv_xl,blk%egv_xr)
#endif
end subroutine deallocate_hydro_block

subroutine allocate_eos_block(blk)
    type(blockdef), pointer :: blk
#if     ieos==2
#if     ieosmodule==1
    call allocate_cell_data_block(blk%H2)
    call allocate_cell_data_block(blk%HI)
    call allocate_cell_data_block(blk%HII)
    call allocate_cell_data_block(blk%electron)
#elif   ieosmodule==2
    call allocate_cell_data_block(blk%HI)
    call allocate_cell_data_block(blk%HII)
#elif   ieosmodule==3
    call allocate_cell_data_block(blk%H2)
    call allocate_cell_data_block(blk%HI)
    call allocate_cell_data_block(blk%HII)
    call allocate_cell_data_block(blk%HeI)
    call allocate_cell_data_block(blk%HeII)
    call allocate_cell_data_block(blk%HeIII)
    call allocate_cell_data_block(blk%electron)
#elif   ieosmodule==4
    call allocate_cell_data_block(blk%H2)
    call allocate_cell_data_block(blk%HI)
#endif
#endif
end subroutine allocate_eos_block

subroutine deallocate_eos_block(blk)
    type(blockdef), pointer :: blk
#if     ieos==2
#if     ieosmodule==1
    deallocate(blk%H2,blk%HI,blk%HII,blk%electron)
#elif   ieosmodule==2
    deallocate(blk%HI,blk%HII)
#elif   ieosmodule==3
    deallocate(blk%H2,blk%HI,blk%HII,blk%HeI,blk%HeII,blk%HeIII,blk%electron)
#elif   ieosmodule==4
    deallocate(blk%H2,blk%HI)
#endif
#endif
end subroutine deallocate_eos_block

subroutine allocate_radiation_block(blk)
    !fully implicit fld radiation transfer
    type(blockdef), pointer :: blk
    if (iradiation/=0) then
        if (lirradiation) call allocate_cell_data_block(blk%irrad)
        if (lirradiation) allocate(blk%irrad0(blk_size_ny))
        if (lfld_heating) call allocate_cell_data_block(blk%fld_heating)
        if (lfld_heating) call allocate_cell_data_block(blk%fld_heating_cell)
        if (lfld_heating.and.lsource) call allocate_cell_data_block(blk%heat_acculm)
        if (lfld_heating.and.lsource.and.lviscous) call allocate_cell_data_block(blk%fld_heating_viscous)
        call allocate_xsurface_data_block(blk%kx)
        call allocate_xsurface_data_block(blk%Fradx)
        if (lrad_adv) call allocate_xsurface_data_block(blk%erad_xflux)
        if (lrad_adv) call allocate_cell_data_block(blk%erad_adv)
        call allocate_xsurface_data_block(blk%tau)
        call allocate_cell_data_block(blk%cv)
        call allocate_cell_data_block(blk%Erad)
        call allocate_cell_data_block(blk%Erad_int)
        call allocate_cell_data_block(blk%Erad_source)
        call allocate_cell_data_block(blk%entropy)
        call allocate_cell_data_block(blk%kappa_planck)
        call allocate_cell_data_block(blk%kappa_abs)
        call allocate_cell_data_block(blk%kappa_rosseland)
        call allocate_cell_data_block(blk%sigma_planck)
        call allocate_cell_data_block(blk%sigma_abs)
        call allocate_cell_data_block(blk%sigma_rosseland)
        call allocate_xsurface_data_block(blk%krx_inter)
        call allocate_xsurface_data_block(blk%aradx)
    end if
end subroutine allocate_radiation_block

subroutine deallocate_radiation_block(blk)
    type(blockdef), pointer :: blk
    if (iradiation/=0) then
        if (lirradiation) deallocate(blk%irrad,blk%irrad0)
        if (lfld_heating) deallocate(blk%fld_heating)
        if (lfld_heating) deallocate(blk%fld_heating_cell)
        if (lfld_heating.and.lsource) deallocate(blk%heat_acculm)
        if (lfld_heating.and.lsource.and.lviscous) deallocate(blk%fld_heating_viscous)
        if (lrad_adv) deallocate(blk%erad_xflux,blk%erad_adv)
        deallocate(blk%kx,blk%Fradx,blk%tau,blk%cv,blk%Erad,blk%Erad_int,blk%Erad_source,blk%entropy,blk%kappa_planck,  &
            blk%kappa_abs,blk%kappa_rosseland,blk%sigma_planck,blk%sigma_abs,blk%sigma_rosseland,blk%krx_inter,blk%aradx)
    end if
end subroutine deallocate_radiation_block

subroutine allocate_viscosity_block(blk)
    type(blockdef), pointer :: blk
    if (lviscous) then
    end if
end subroutine allocate_viscosity_block

subroutine allocate_source_block(blk)
    type(blockdef), pointer :: blk
    real(8) :: r,pos(3),dxyz(3),m
    integer :: i
#if     usersource
    call allocate_cell_data_block(blk%iso_temp)
#endif
end subroutine allocate_source_block

subroutine deallocate_source_block(blk)
    type(blockdef), pointer :: blk
#if     usersource
    deallocate(blk%iso_temp)
#endif
end subroutine deallocate_source_block

subroutine allocate_cell_data_block_1d(dataarray)
    real(8), allocatable :: dataarray(:,:,:)
    allocate(dataarray(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    dataarray=0d0
end subroutine allocate_cell_data_block_1d

subroutine allocate_cell_data_block_1d_int(dataarray)
    integer, allocatable :: dataarray(:,:,:)
    allocate(dataarray(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    dataarray=0d0
end subroutine allocate_cell_data_block_1d_int

subroutine allocate_cell_data_block_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    allocate(dataarray(arraydim,blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    dataarray=0d0
end subroutine allocate_cell_data_block_nd

subroutine allocate_cell_data_block_nd_int(dataarray,arraydim)
    integer, allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    allocate(dataarray(arraydim,blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    dataarray=0d0
end subroutine allocate_cell_data_block_nd_int

subroutine allocate_xsurface_data_block_1d(dataarray)
    real(8), allocatable :: dataarray(:,:,:)
    allocate(dataarray(blk_xlb:blk_xub-1,blk_ylb:blk_yub,1))
    dataarray=0d0
end subroutine allocate_xsurface_data_block_1d

subroutine allocate_xsurface_data_block_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    allocate(dataarray(arraydim,blk_xlb:blk_xub-1,blk_ylb:blk_yub,1))
    dataarray=0d0
end subroutine allocate_xsurface_data_block_nd

subroutine allocate_ysurface_data_block_1d(dataarray)
    real(8), allocatable :: dataarray(:,:,:)
    allocate(dataarray(blk_xlb:blk_xub,blk_ylb:blk_yub-1,1))
    dataarray=0d0
end subroutine allocate_ysurface_data_block_1d

subroutine allocate_ysurface_data_block_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    allocate(dataarray(arraydim,blk_xlb:blk_xub,blk_ylb:blk_yub-1,1))
    dataarray=0d0
end subroutine allocate_ysurface_data_block_nd

function ipassive(passivename)
    character(len=16) :: passivename
    type(passive_obj), pointer :: po
    integer :: i,ipassive
    po=>passive_list
    do i=1,npassive
        if (po%passivename==passivename) then
            ipassive=i
            exit
        end if
        if (i/=npassive) po=>po%next
    end do
    if (po%passivename/=passivename) then
        print *,'no such passivename'
        stop
    end if
end function ipassive

function find_the_slope(diff_left,diff_right,diff_central)
    real(8) :: diff_left,diff_right,diff_central,find_the_slope,r
#if     isolver==1
    if (diff_central/=zero) then
        if (diff_left*diff_right<=0d0) then
            find_the_slope=0d0
        else
            find_the_slope=abs(diff_central)/diff_central*min(abs(diff_left),abs(diff_right))
        end if
    else
        find_the_slope=0d0
    end if
#elif   isolver==2
    if (diff_central/=zero) then
        if (diff_left*diff_right<=0) then
            find_the_slope=0d0
        else
            find_the_slope=abs(diff_central)/diff_central*0.5d0*min(abs(diff_left),abs(diff_right))
        end if
    else
        find_the_slope=0d0
    end if
#endif
end function find_the_slope

subroutine calculate_slp_cell(blk,i,j,slp)
    !calculate 1d and 2d slope
    !hydro uses the primitive variables, passive uses the specific passive scalars
    type(blockdef), pointer :: blk
    character(len=16) :: obj
    integer :: i,j,l,k
    real(8), allocatable :: v(:,:),slp_cell(:,:),slp(:,:)
    real(8) :: xc(3),yc(3),dx,dy
    l=cell_var_length
    if (iradiation/=0) l=l+1
    allocate(v(l,3),slp_cell(l,3))
    v(1:5,1:3)=blk%w(1:5,i-1:i+1,1,1)
    if (iradiation/=0) v(l,1:3)=blk%Erad(i-1:i+1,1,1)/blk%w(1,i-1:i+1,1,1)
    if (lpassive) then
        v(6:5+npassive,1)=blk%passive_scalar(1:npassive,i-1,1,1)/blk%w(1,i-1,1,1)
        v(6:5+npassive,2)=blk%passive_scalar(1:npassive,i,1,1)/blk%w(1,i,1,1)
        v(6:5+npassive,3)=blk%passive_scalar(1:npassive,i+1,1,1)/blk%w(1,i+1,1,1)
    end if
    if (igeometry==0) then
    else if (igeometry==2) then
        xc=blk%x_center(i-1:i+1)
        slp_cell(1:l,1)=(v(1:l,2)-v(1:l,1))/(xc(2)-xc(1))
        slp_cell(1:l,2)=(v(1:l,3)-v(1:l,2))/(xc(3)-xc(2))
        slp_cell(1:l,3)=(v(1:l,3)-v(1:l,1))/(xc(3)-xc(1))
        do k=1,l
            slp(k,1)=find_the_slope(slp_cell(k,1),slp_cell(k,2),slp_cell(k,3))
        end do
    end if
    deallocate(v,slp_cell)
end subroutine calculate_slp_cell

subroutine initialize_global_parameters_int(pvalue,p,printout)
    !give pvalue to p, p is protected
    !only this subroutine can pass values to global parameters
    integer :: pvalue,p
    character(len=*) :: printout
    p=pvalue
    if (rank==0) write(*,'(A32,I18)') printout,p
end subroutine initialize_global_parameters_int

subroutine initialize_global_parameters_real(pvalue,p,printout)
    !give pvalue to p, p is protected
    !only this subroutine can pass values to global parameters
    real(8) :: pvalue,p
    character(len=*) :: printout
    p=pvalue
    if (rank==0) write(*,'(A32,ES18.8E3)') printout,p
end subroutine initialize_global_parameters_real

subroutine initialize_global_parameters_logical(pvalue,p,printout)
    !give pvalue to p, p is protected
    !only this subroutine can pass values to global parameters
    logical :: pvalue,p
    character(len=*) :: printout
    p=pvalue
    if (rank==0) write(*,'(A32,L18)') printout,p
end subroutine initialize_global_parameters_logical

subroutine initialize_mesh()
    n_hydro_guard=2
    nx_blks=nx/blk_size_nx
    ny_blks=ny/blk_size_ny
    blk_xlb=1-n_hydro_guard
    blk_xub=blk_size_nx+n_hydro_guard
    blk_ylb=1;blk_yub=1
    blk_hydro_size=5*(blk_size_nx+2*n_hydro_guard)
    blk_cell_size=blk_size_nx+2*n_hydro_guard
    blk_interface_size=blk_cell_size-1
    !blksize=2*blk_hydro_size+2*blk_cell_size
    blksize=blk_hydro_size
    if (iradiation/=0) blksize=blksize+blk_cell_size
    if (lam_con) blksize=blksize+blk_cell_size
    if (rank==0) then
        write(*,'(A20,I15)') 'n_hydro_guard=',n_hydro_guard
    end if
end subroutine initialize_mesh

subroutine blk_traversal1(sub)
    type(blockdef), pointer :: blk
    procedure(blk_operator) :: sub
    integer :: i
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        call sub(blk)
        if (.not.associated(blk,llist_tail)) then
            blk=>blk%next
        end if
    end do
end subroutine blk_traversal1

subroutine blk_traversal2(sub,x)
    type(blockdef), pointer :: blk
    procedure(blk_operator_real) :: sub
    real(8) :: x
    integer :: i
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        call sub(blk,x)
        if (.not.associated(blk,llist_tail)) then
            blk=>blk%next
        end if
    end do
end subroutine blk_traversal2

subroutine blk_traversal_ind(sub)
    type(blockdef), pointer :: blk
    procedure(blk_operator) :: sub
    blk=>llist_head
    do while (.not.associated(blk,llist_tail))
        call sub(blk)
        blk=>blk%next
    end do
    call sub(blk)
end subroutine blk_traversal_ind

subroutine blk_global_traversal(sub)
    type(blockdef), pointer :: blk
    procedure(blk_operator) :: sub
    blk=>blk_head
    do while (associated(blk))
        call sub(blk)
        blk=>blk%next
    end do
    nullify(blk)
end subroutine blk_global_traversal

subroutine field_calculation(blk,sub,ap)
    type(blockdef), pointer :: blk
    procedure(field_calculator), pointer :: sub
    real(8), dimension(:,:,:), allocatable :: ap
    integer :: i,j
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            ap(i,j,1)=sub(blk,i,j)
        end do
    end do
end subroutine field_calculation

subroutine sync_bcast_logical1(lvalue,gvalue)
    !each processor has local value and shall be synced by a global value
    logical, allocatable :: lvalue(:),gvalue(:)
    integer :: i,j,reqs,ierr,stat(MPI_STATUS_SIZE)
    call mpi_isend(lvalue,np_nblk(rank+1),MPI_LOGICAL,0,1,MPI_COMM_WORLD,reqs,ierr)
    if (rank==0) then
        j=1
        do i=0,np-1
            if (i>0) j=j+np_nblk(i)
            call mpi_recv(gvalue(j),np_nblk(i+1),MPI_LOGICAL,i,1,MPI_COMM_WORLD,stat,ierr)
        end do
    end if
    call mpi_bcast(gvalue,nblk_total,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_request_free(reqs,ierr)
end subroutine sync_bcast_logical1

subroutine sync_bcast_logical2(lvalue,gvalue)
    logical, allocatable :: lvalue(:,:),gvalue(:,:)
    integer :: i,j,d(2),reqs,ierr,stat(MPI_STATUS_SIZE),send_count
    d=shape(gvalue)
    send_count=np_nblk(rank+1)*d(1)
    call mpi_isend(lvalue,send_count,MPI_LOGICAL,0,1,MPI_COMM_WORLD,reqs,ierr)
    if (rank==0) then
        j=1
        do i=0,np-1
            if (i>0) j=j+np_nblk(i)
            call mpi_recv(gvalue(1,j),np_nblk(i+1)*d(1),MPI_LOGICAL,i,1,MPI_COMM_WORLD,stat,ierr)
        end do
    end if
    call mpi_bcast(gvalue,nblk_total*d(1),MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_request_free(reqs,ierr)
end subroutine sync_bcast_logical2

subroutine sync_bcast_int1(lvalue,gvalue)
    integer, allocatable :: lvalue(:),gvalue(:)
    integer :: i,j,reqs,ierr,stat(MPI_STATUS_SIZE)
    call mpi_isend(lvalue,np_nblk(rank+1),MPI_INTEGER,0,1,MPI_COMM_WORLD,reqs,ierr)
    if (rank==0) then
        j=1
        do i=0,np-1
            if (i>0) j=j+np_nblk(i)
            call mpi_recv(gvalue(j),np_nblk(i+1),MPI_INTEGER,i,1,MPI_COMM_WORLD,stat,ierr)
        end do
    end if
    call mpi_bcast(gvalue,nblk_total,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_request_free(reqs,ierr)
end subroutine sync_bcast_int1

subroutine sync_bcast_int2(lvalue,gvalue,n_oper)
    integer, allocatable :: lvalue(:,:),gvalue(:,:),n_oper(:)
    integer :: i,j,d(2),reqs,ierr,stat(MPI_STATUS_SIZE),send_count,recv_count
    d=shape(gvalue)
    send_count=n_oper(rank+1)*d(1)
    if (send_count>0) call mpi_isend(lvalue,send_count,MPI_INTEGER,0,1,MPI_COMM_WORLD,reqs,ierr)
    if (rank==0) then
        do i=0,np-1
            if (i>0) then
                j=sum(n_oper(1:i))+1
            else
                j=1
            end if
            recv_count=n_oper(i+1)*d(1)
            if (recv_count>0) call mpi_recv(gvalue(1,j),recv_count,MPI_INTEGER,i,1,MPI_COMM_WORLD,stat,ierr)
        end do
    end if
    send_count=sum(n_oper)*d(1)
    call mpi_bcast(gvalue,send_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_request_free(reqs,ierr)
end subroutine sync_bcast_int2

subroutine collective_sub(op_blk,operation1,operation2,outcome)
    !the subroutine that carries out the collective operation
    !op_blk operates on each block of every processor
    !op_proc operates on each processor, save the result to processor%collective_processor_result
    type(blockdef), pointer :: blk
    procedure(collective_oper_blk) :: op_blk
    character(len=16) :: operation1,operation2
    real(8) :: blk_result,outcome
    real(8), dimension(:), allocatable :: blk_result_array
    integer :: i,ierr
    processor%collective_processor_result=0d0
    processor%collective_array=0d0
    processor%collective_result=0d0
    allocate(blk_result_array(np_nblk(rank+1)))
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        blk_result_array(i)=op_blk(blk)
        blk=>blk%next
    end do
    if (operation1==op_min) then
        processor%collective_processor_result=minval(blk_result_array)
    else if (operation1==op_max) then
        processor%collective_processor_result=maxval(blk_result_array)
    else if (operation1==op_sum) then
        processor%collective_processor_result=sum(blk_result_array)
    end if
    call mpi_gather(processor%collective_processor_result,1,MPI_REAL8,processor%collective_array,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    if (operation2==op_min) then
        processor%collective_result=minval(processor%collective_array)
    else if (operation2==op_max) then
        processor%collective_result=maxval(processor%collective_array)
    else if (operation2==op_sum) then
        processor%collective_result=sum(processor%collective_array)
    end if
    call mpi_bcast(processor%collective_result,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    outcome=processor%collective_result
    deallocate(blk_result_array)
end subroutine collective_sub

subroutine check_global_parameter_consistency()
    integer :: max_static_level
    character(len=256) :: alert
#if     ieos==2
    if (opacity_gas_rho_min<=environment%rhomin) then
        alert='ieos=2 and opacity_gas_rho_min<=environment%rhomin'
        call abort_guangqi(alert)
    end if
    if (opacity_gas_rho_max>=environment%rhomax) then
        alert='ieos=2 and opacity_gas_rho_max>=environment%rhomax'
        call abort_guangqi(alert)
    end if
#endif
    if (refine_type=='static'.or.refine_type=='mixed') then
        max_static_level=maxval(refine_region(:,5))
        if (max_refine_level<max_static_level) then
            alert='max_refine_level<max_static_level'
            call abort_guangqi(alert)
        end if
    end if
    if (nd==1) then
        if (n_domain(4)/=n_domain(3)) then
            alert='in 1d, y direction should be 0'
            call abort_guangqi(alert)
        end if
    end if
    if (nd/=1) then
        if (abs(sum(hydro_bound_type-passive_bound_type))/=0) then
            print *,'passive_bound needs to follow hydro_bound'
            call abort_guangqi(alert)
        end if
        if (hydro_bound_type(1)==3.or.hydro_bound_type(2)==3) then
            if (hydro_bound_type(1)/=3.or.hydro_bound_type(2)/=3) then
                alert='xl and xu must be periodic at the same time'
                call abort_guangqi(alert)
            end if
        end if
        if (hydro_bound_type(3)==3.or.hydro_bound_type(4)==3) then
            if (hydro_bound_type(3)/=3.or.hydro_bound_type(4)/=3) then
                alert='yl and yu must be periodic at the same time'
                call abort_guangqi(alert)
            end if
        end if
        if (igeometry==1) then
            alert='using cylindrical coordinate'
            call print_guangqi(alert)
            if (n_domain(1)<=0d0) then
                alert='r inner domain <= 0'
                call abort_guangqi(alert)
            end if
            if (n_domain(3)<0d0) then
                alert='phi lower domain < 0.'
                call abort_guangqi(alert)
            end if
            if (n_domain(4)>2d0*pi) then
                alert='phi upper domain > 2pi.'
                call abort_guangqi(alert)
            end if
        else if (igeometry==2) then
            alert='using polar coordinate'
            call print_guangqi(alert)
            if (n_domain(1)<=0d0) then
                alert='r inner domain <= 0'
                call abort_guangqi(alert)
            end if
            if (n_domain(4)>=pi) then
                alert='theta upper domain >= pi.'
                call abort_guangqi(alert)
            end if
        end if
        if (iradiation==0) then
            if (lrad_adv) then
                alert='no radiation, thus no radiation advection'
                call abort_guangqi(alert)
            end if
        end if
#if     ieosmodule==3
        if (environment%h_ratio>=1.or.environment%h_ratio<=0) then
            alert='H and He must both exist to use h and he mixture eos'
            call abort_guangqi(alert)
        end if
#endif
    end if
end subroutine check_global_parameter_consistency

subroutine update_protected_logical(a,b)
    logical :: a,b
    a=b
end subroutine update_protected_logical

subroutine update_protected_real(a,b)
    real(8) :: a,b
    a=b
end subroutine update_protected_real

subroutine update_protected_int(a,b)
    integer :: a,b
    a=b
end subroutine update_protected_int

!subroutine intp_interface_fix_x(blk,r,sub,var,belong)
!    type(blockdef), pointer :: blk
!    procedure(blk_profiling), pointer :: sub
!    logical :: belong
!    real(8) :: r
!    real(8), dimension(:), allocatable :: var
!    logical :: l1,l2,l3
!    belong=.false.
!    l1=(r==n_domain(1).and.blk%key(1)==1)
!    l2=(r==n_domain(2).and.blk%key(1)==nx_blks*2**blk%level)
!    l3=(r<blk%x_interface(blk_size_nx).and.r>=blk%x_interface(0))
!    if (l1.or.l2.or.l3) then
!        allocate(blk%blk_var(blk_size_ny))
!        call sub(blk,r)
!        var=blk%blk_var
!        deallocate(blk%blk_var)
!        belong=.true.
!    end if
!end subroutine intp_interface_fix_x

!subroutine intp_cell_fix_r(blk,r,sub,var,belong)
!    !2d polar coordinate
!    type(blockdef), pointer :: blk
!    procedure(blk_profiling) :: sub
!    logical :: belong
!    real(8) :: r
!    real(8), dimension(:), allocatable :: var
!    belong=.false.
!    if (r<blk%x_interface(blk_size_nx).and.r>=blk%x_interface(0)) then
!        allocate(blk%blk_var(blk_size_ny))
!        !call sub(blk,r)
!        var=blk%blk_var
!        deallocate(blk%blk_var)
!        belong=.true.
!    end if
!end subroutine intp_cell_fix_r
!
!subroutine intp_cell_fix_theta(blk,theta,sub,var,belong)
!    !2d polar coordinate
!    type(blockdef), pointer :: blk
!    procedure(blk_operator) :: sub
!    logical :: belong
!    real(8) :: theta
!    real(8), dimension(:), allocatable :: var
!    belong=.false.
!    if (theta<blk%y_interface(blk_size_ny).and.theta>=blk%y_interface(0)) then
!        allocate(blk%blk_var(blk_size_nx))
!        call sub(blk)
!        var=blk%blk_var
!        deallocate(blk%blk_var)
!        belong=.true.
!    end if
!end subroutine intp_cell_fix_theta

subroutine diagnostics()
end subroutine diagnostics

function blk_select(key1,key2)
    integer :: key1(3),key2(3)
    logical :: blk_select
    blk_select=.false.
    if (sum(abs(key1-key2))==0) blk_select=.true.
end function blk_select

subroutine check_nan()
    type(blockdef), pointer :: blk
    integer :: i,j
    character(len=128) :: str
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        do j=blk_xlb,blk_xub
            if (isnan(blk%w(1,j,1,1)).or.blk%w(1,j,1,1)<=0) then
                print *,'density',time_sys%ntimestep,rank,i,j,blk%key,blk%level
                str='nan or non-positive'
                call abort_guangqi(str)
            end if
            if (isnan(blk%temp(j,1,1))) then
                print *,'temp',time_sys%ntimestep,rank,i,j,blk%key,blk%level
                str='nan'
                call abort_guangqi(str)
            end if
            if (isnan(blk%egv(j,1,1))) then
                print *,'egv',time_sys%ntimestep,rank,i,j,blk%key,blk%level
                str='nan'
                call abort_guangqi(str)
            end if
            if (blk%w(5,j,1,1)<0) then
                print *,'pres',time_sys%ntimestep,rank,i,j,blk%blk_id,blk%next%blk_id
                print *,blk%w(5,j,1,1)
                print *,blk%next%w(5,:,1,1)
                print *,blk%level,blk%next%level
                str='negative'
                call abort_guangqi(str)
            end if
        end do
        blk=>blk%next
    end do
end subroutine check_nan

subroutine block_id_to_processor_rank(load,blk_id,rank_id)
    !given a load and a blk_id of a block, calculate the rank_id of the block
    integer, dimension(:), allocatable :: load
    integer :: blk_id,rank_id,i,istart,iend
    istart=1;iend=load(1)
    do i=1,np
        if (blk_id>=istart.and.blk_id<=iend) then
            rank_id=i-1
            exit
        end if
        if (i/=np) then
            istart=istart+load(i)
            iend=iend+load(i+1)
        end if
    end do
end subroutine block_id_to_processor_rank

subroutine abort_guangqi(str)
    integer :: errorcode,ierr
    character(len=128), intent(in) :: str
    if (rank==0) print *,trim(str)
    call PetscFinalize(ierr)
    stop
    !call mpi_finalize(ierr)
    !call mpi_abort(MPI_COMM_WORLD,errorcode,ierr)
    !stop
end subroutine abort_guangqi

subroutine print_guangqi(str)
    integer :: ierr
    character(len=128), intent(in) :: str
    if (rank==0) print *,trim(str)
end subroutine print_guangqi

end module datastructure

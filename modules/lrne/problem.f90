!the 1d module is for single processor only
module problem
use phylib
use datastructure
use mathlib
use eos
implicit none
real(8), protected :: m_star,v_esc,rho_floor,ms,lct,t_ram,floor_tauR,disk_rho0,rin,rout,tmax,rho0,temp0,const_kr,const_kp
logical, protected :: ex
real(8), dimension(:), allocatable :: time,rho,v,temp,asym_outflow,eratio
real(8), dimension(:,:), allocatable :: boundcond
integer, protected :: ejection_rad_bc_type
character(len=32), protected :: post_ej_bound
integer, dimension(2) :: diminput 


contains

subroutine initialize_problem()
    !calculate w, u, temp, and egv for all location
    character(len=128) :: str
    real(8) :: p_wind,p_clump,xlen,ylen,dt_record
    real(8) :: eg,erad,cs,t,mach
    integer :: i,j,m,n,o,error,ncols
    character(len=30) :: subdir1,subdir2
    character(len=5) :: str1,str2,str3,str4
    namelist /parameters_1d/ m_star,lfld_mom,larad,dt_record,record_length,             &
        petsc_qratio,petsc_iter,floor_tauR,inner_floor_temp,rho_floor,post_ej_bound,    &
        t_ram,lpradgradv,rho0,temp0,ejection_rad_bc_type
    open(unit=11,file=trim(path_root)//'/problem.data',status='old',action='read')
    read(unit=11,nml=parameters_1d)
    close(11)
    if (restart) lhis_overwrite=.false.
    if (refine_type=='adaptive'.or.refine_type=='mixed') then
        print *,'this module only accept SMR or fixed grid'
        stop
    end if
    path_out=trim(path_root)//trim(out_dir)
    inquire(file="bcinput.dat",exist=ex)
    if (ex.eqv..false.) then
        str='need boundary input file, bcinput.dat'
        call abort_guangqi(str)
    end if
    ncols=6
    if (rank==0) then
        call read_file('bcinput.dat',boundcond,ncols)
        diminput=shape(boundcond)
    end if
    call mpi_barrier(MPI_COMM_WORLD,error)
    call mpi_bcast(diminput,2,MPI_INTEGER,0,MPI_COMM_WORLD,error)
    if (.not.allocated(boundcond)) allocate(boundcond(diminput(1),diminput(2)))
    call mpi_bcast(boundcond,size(boundcond),MPI_REAL8,0,MPI_COMM_WORLD,error)
    allocate(time(diminput(1)),rho(diminput(1)),v(diminput(1)),temp(diminput(1)),eratio(diminput(1)),asym_outflow(diminput(1)))
    !the ejecta is time and could be angle dependent, we assume fixed erad/eg ratio, temp and rho are mid-plane quantities
    time=boundcond(:,1)
    v=boundcond(:,2)
    rho=boundcond(:,3)
    temp=boundcond(:,4)
    eratio=boundcond(:,5)
    asym_outflow=boundcond(:,6)*pi
    post_ej_bound=trim(post_ej_bound)
    tmax=maxval(time)
    t_ram=t_ram*day
    wedge_polar=.true.
    lfloor_temp=.true.
    rho_thresh_petsc1=1e-16
    rho_thresh_petsc2=1e-16
    kappa_planck_ceil1=1d0
    kappa_planck_ceil2=1d0
    opacity_gas_rho_min=1d-18
    opacity_gas_rho_max=1d-4
    opacity_gas_t_min=200
    opacity_gas_t_max=3d6
    !rin=n_domain(1)*rin
    !rout=n_domain(2)*rout
    tempfloorscale=1
    ms=m_star*msun
    sig_rosseland_floor=floor_tauR/n_domain(2)
    central_star%core%mass=m_star*msun
    v_esc=sqrt(2*gr*m_star*msun/n_domain(1))
    time_sys%dt_record=dt_record
    lct=n_domain(2)/c_light         !light crossing time
    cs=adiabatic_cs(rho(1),temp(1))
    mach=v(1)/cs
    !initialize all the cells including boundary cells
    if (rank==0) then
        write(*,'(A32,ES18.6E2)') 'ejecta Mach number',mach
        write(*,'(A32,ES18.6E2)') 'escape velocity',v_esc
        write(*,'(A32,ES18.6E2)') 'dt_record',time_sys%dt_record
    end if
    if (np>1) then
        print *,'use single processor'
        stop
    end if
end subroutine initialize_problem

subroutine read_file(filename,boundcond,ncols)
    character(len=*), intent(in) :: filename            ! Input file name
    real(8), dimension(:,:), allocatable :: boundcond   ! Array to store lines
    integer :: num_rows                                 ! Number of rows in the file
    integer :: io_status,i,unit_number,ncols
    real(8), allocatable :: line(:)                       ! Temporary buffer for reading lines
    character(len=32), allocatable :: header(:)
    character(len=32) :: form1,form2
    character(len=1) strncols
    unit_number=666
    open(newunit=unit_number, file=filename, status='old', action='read')
    allocate(line(ncols),header(ncols))
    write(strncols, '(i1)' ) ncols
    form1='('//strncols//'A16)'
    form2='('//strncols//'ES16.8E2)'
    ! First pass: Count the number of rows
    read(unit_number,trim(form1)) header
    num_rows = 0
    do
      read(unit_number,trim(form2),iostat=io_status) line
      if (io_status /= 0) exit
      num_rows = num_rows + 1
    end do
    ! Allocate the array to hold the lines
    allocate(boundcond(num_rows,ncols))
    ! Rewind the file to read the lines again
    rewind(unit_number)
    read(unit_number,trim(form1)) header
    ! Second pass: Read the lines into the array
    do i = 1, num_rows
      read(unit_number,trim(form2)) boundcond(i,:)
    end do
    deallocate(line,header)
    close(unit_number)
end subroutine read_file

!*****************************************************initial conditions*************************************************

subroutine initial_hydro(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: pos(3),t,w(5),u(5),temp,egv,r0,rr,r
    r=blk%x_center(i)
    w(1)=max(rho_profile(r),initial_floor_rho(r),rho_floor)
    w(2)=vx_profile(r)
    temp=temp_profile(r)
    w(3:4)=0d0
    w(5)=prhot(w(1),temp)
    call assign_w_to_u_cell(blk,w,i,j)
end subroutine initial_hydro

subroutine initial_rad(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: temp
    if (iradiation==4) then
        temp=blk%temp(i,j,1)
        blk%Erad(i,j,1)=a_rad*temp**4
    end if
end subroutine initial_rad

subroutine initial_scalar(blk,i,j)
    type(blockdef), pointer :: blk
    real(8) :: r,rho,l_omega
    integer :: i,j,ip_am
end subroutine initial_scalar

function initial_floor_rho(r)
    real(8) :: x,r,initial_floor_rho
    x=r/n_domain(1)
    initial_floor_rho=1e-15/x
end function initial_floor_rho

function rho_profile(r)
    real(8) :: r,rho_profile,r0
    r0=n_domain(1)
    rho_profile=rho0*(r0/r)**1.5d0
end function rho_profile

function vx_profile(r)
    real(8) :: r,vx_profile
    vx_profile=0d0
end function vx_profile

function temp_profile(r)
    real(8) :: r,temp_profile,r0
    r0=n_domain(1)
    temp_profile=temp0*(r/r0)**(-1d0)
end function temp_profile

!*****************************************************boundary conditions*************************************************

subroutine boundary_hydro(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: r,t,w(5),u(5),temp,egv
    t=time_sys%t
    r=blk%x_center(i)
    if (r<n_domain(1)) then
        w(2)=vej(t)
        w(1)=rhoej(t)
        w(3:4)=0d0
        temp=tempej(t)
        w(5)=prhot(w(1),temp)
        call assign_w_to_u_cell(blk,w,i,j)
    end if
end subroutine boundary_hydro

subroutine boundary_rad(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: v,r,r1,temp,t,rho
    if (iradiation==4) then
        r=blk%x_center(i)
        t=time_sys%t
        if (r<n_domain(1)) then
            !inner boundary
            temp=tempej(t)
            blk%Erad(i,j,1)=a_rad*temp**4d0
        end if
        if (r>n_domain(2)) then
            !outer boundary
            r1=blk%x_center(blk_size_nx)
            blk%Erad(i,j,1)=blk%Erad(blk_size_nx,j,1)*(r1/r)**2
        end if
    end if
end subroutine boundary_rad

function erad_to_egv_ratio(x)
    real(8), allocatable :: x(:)
    real(8) :: rho,temp,erad,egv,eratio,erad_to_egv_ratio
    temp=x(1)
    rho=x(2)
    eratio=x(3)
    erad=a_rad*temp**4
    egv=egvrhot(rho,temp)
    erad_to_egv_ratio=eratio-erad/egv
end function erad_to_egv_ratio

function solve_for_temp(rho,eratio)
    real(8) :: rho,eratio,lowT,highT,convergence,solve_for_temp
    integer :: conv_mode
    procedure(fun), pointer :: ptr
    real(8), allocatable :: x(:),root(:)
#if     ieos==2
    conv_mode=2
    convergence=1e-5
    lowT=5e2        !pure molecular
    highT=1e6       !fully ionized
    allocate(x(3),root(3))
    x=(/1d4,rho,eratio/)
    ptr=>erad_to_egv_ratio
    call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
    solve_for_temp=root(1)
    deallocate(x,root)
#elif   ieos==1
    solve_for_temp=pow(eratio*rho/maw/mh/(gamma_gas-1)/a_rad*kb,1d0/3d0)
#endif
end function solve_for_temp

subroutine boundary_scalar(blk,i,j)
    type(blockdef), pointer :: blk
    real(8) :: r,rho,l_omega
    integer :: i,j,ip_am
end subroutine boundary_scalar

!*******************time dependent boundary conditions**************************************************************

function vej(t)
    real(8) :: vej,t
    call interpolation_linear(t,vej,time,v)
end function vej

function rhoej(t)
    real(8) :: rhoej,t
    call interpolation_linear(t,rhoej,time,rho)
end function rhoej

function tempej(t)
    real(8) :: tempej,t
    call interpolation_linear(t,tempej,time,temp)
end function tempej

!end of time dependent boundary conditions

subroutine time_dependent_bound_type()
    real(8) :: t,temp_inner,v
    integer :: i
    t=time_sys%t
    if (t>tmax) then
        !when the ejecting phase is done, the inner hydro boundary condition changes to free
        !the inner radiation boundary condition changes to zero-gradient
        if (post_ej_bound=='reflective') then
            hydro_bound_type(1)=2
        else if (post_ej_bound=='transmissive') then
            hydro_bound_type(1)=4
        else if (post_ej_bound=='constant_pres') then
            hydro_bound_type(1)=9
        else
            print *,'unacceptable boundary condition'
            stop
        end if
        rad_bound_type(1)=8
    else if (t<t_ram) then
        !cannot allow too hot ejecta come enters the computational domain immediately
        hydro_bound_type(1)=9
        rad_bound_type(1)=8
    else
        !injecting gas into the computational domain
        hydro_bound_type(1)=9
        rad_bound_type(1)=ejection_rad_bc_type
    end if
end subroutine time_dependent_bound_type

subroutine assemble_record_array(record_array)
    type(blockdef), pointer :: blk
    real(8), allocatable :: record_array(:)
    real(8) :: sample_r(2),v,arad_work,prad_gradv,sum_record(8)
    procedure(field_calculator), pointer :: g1,g2,g3,g4
    integer :: i,error,j,iloc
    record_array=0
    record_array(1)=time_sys%t
    if (rank==0) then
        blk=>llist_head
        if (iradiation==4) then
            record_array(2)=4*pi*blk%x_interface(0)**2*blk%Fradx(0,1,1)
        else
            record_array(2)=0d0
        end if
        record_array(3)=4*pi*blk%x_interface(0)**2*blk%xflux(1,0,1,1)/msun*year
        if (iradiation==4) record_array(6)=4*pi*blk%x_interface(0)**2*blk%erad_xflux(0,1,1)
        record_array(7)=4*pi*blk%x_interface(0)**2*blk%xflux(5,0,1,1)
        if (abs(record_array(3))<1d-10) then
            record_array(3)=0d0
            record_array(7)=0d0
        end if
    end if
    if (rank==np-1) then
        blk=>blk_tail
        if (iradiation==4) then
            record_array(4)=4d0*pi*blk%x_interface(blk_size_nx)**2*blk%Fradx(blk_size_nx,1,1)
        else
            record_array(4)=0d0
        end if
        record_array(5)=4*pi*blk%x_interface(blk_size_nx)**2*blk%xflux(1,blk_size_nx,1,1)/msun*year
        if (iradiation==4) record_array(8)=4*pi*blk%x_interface(blk_size_nx)**2*blk%erad_xflux(blk_size_nx,1,1)
        record_array(9)=4*pi*blk%x_interface(blk_size_nx)**2*blk%xflux(5,blk_size_nx,1,1)
    end if
end subroutine assemble_record_array

subroutine problem_oper()
end subroutine problem_oper

subroutine finalize_problem()
end subroutine finalize_problem

function user_kr(rho,tgas)
    real(8) :: user_kr,rho,tgas
    user_kr=const_kr
end function user_kr

function user_kp(rho,Erad,tgas)
    real(8) :: user_kp,rho,Erad,tgas
    user_kp=const_kp
end function user_kp

end module problem

!the time dependent ejecta is injected to the computational domain from the inner boundary
!the ejecta information, including time, rho, vx, and temp are read in from bcinput.dat
module problem
use phylib
use datastructure
use mathlib
use eos
implicit none
real(8), dimension(5), protected :: w0,u0
real(8), protected :: egv0,E0,m_star,v_esc,rho_floor,      &
    p0,ms,sd_theta,lct,t_ram,floor_tauR,rin,rout,tmax,diffused_rho0,diffused_temp0
logical, protected :: loutflow_focusing,lam,ex
real(8), dimension(:), allocatable :: time,rho,v,temp,asym_outflow
real(8), dimension(:,:), allocatable :: boundcond
integer, dimension(2) :: diminput 


contains

subroutine initialize_problem()
    !calculate w, u, temp, and egv for all location
    character(len=128) :: str
    real(8) :: p_wind,p_clump,xlen,ylen,dt_record
    real(8) :: eg,erad,theta,dtheta,t
    integer :: i,j,m,n,o,error,ncols
    character(len=30) :: subdir1,subdir2
    character(len=5) :: str1,str2,str3,str4
    namelist /parameters_1d/ m_star,lfld_mom,larad,dt_record,record_length,   &
        petsc_qratio,petsc_iter,floor_tauR,inner_floor_temp,diffused_rho0,diffused_temp0,rho_floor
    open(unit=11,file=trim(path_root)//'/problem.data',status='old',action='read')
    read(unit=11,nml=parameters_1d)
    close(11)
    if (.not.restart) lhis_overwrite=.true.
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
    ncols=4
    if (rank==0) then
        call read_file('bcinput.dat',boundcond,ncols)
        diminput=shape(boundcond)
    end if
    call mpi_barrier(MPI_COMM_WORLD,error)
    call mpi_bcast(diminput,2,MPI_INTEGER,0,MPI_COMM_WORLD,error)
    if (.not.allocated(boundcond)) allocate(boundcond(diminput(1),diminput(2)))
    call mpi_bcast(boundcond,size(boundcond),MPI_REAL8,0,MPI_COMM_WORLD,error)
    allocate(time(diminput(1)),rho(diminput(1)),v(diminput(1)),temp(diminput(1)),asym_outflow(diminput(1)))
    time=boundcond(:,1)
    v=boundcond(:,2)
    rho=boundcond(:,3)
    temp=boundcond(:,4)
    tmax=maxval(time)
    t_ram=0.01*day          !the first t_ram day will enforce gradErad=0
    lfloor_temp=.true.
    rho_thresh_petsc1=rho_floor
    rho_thresh_petsc2=rho_floor
    kappa_planck_ceil1=1d0
    kappa_planck_ceil2=1d0
    opacity_gas_rho_min=1d-18
    opacity_gas_rho_max=1d-4
    opacity_gas_t_min=200
    opacity_gas_t_max=3d6
    rin=n_domain(1)*rin
    rout=n_domain(2)*rout
    tempfloorscale=1
    ms=m_star*msun
    !diffused_focusing=diffused_focusing*pi
    sig_rosseland_floor=floor_tauR/n_domain(2)
    central_star%core%mass=m_star*msun
    v_esc=sqrt(2*gr*m_star*msun/n_domain(1))
    !sd_theta=sd_theta*pi
    !call update_time_sys(timescale,dt_record)
    time_sys%dt_record=dt_record
    lct=n_domain(2)/c_light         !light crossing time
    !initialize all the cells including boundary cells
    if (rank==0) then
        write(*,'(A32,ES18.6E2)') 'escape velocity',v_esc
        write(*,'(A32,ES18.6E2)') 'dt_record',time_sys%dt_record
    end if
end subroutine initialize_problem

subroutine read_file(filename,boundcond,ncols)
    character(len=*), intent(in) :: filename            ! input file name
    real(8), dimension(:,:), allocatable :: boundcond   ! array to store lines
    integer :: num_rows                                 ! number of rows in the file
    integer :: io_status,i,unit_number,ncols
    real(8), allocatable :: line(:)                       ! temporary buffer for reading lines
    character(len=32), allocatable :: header(:)
    character(len=32) :: form1,form2
    character(len=1) strncols
    unit_number=666
    open(newunit=unit_number, file=filename, status='old', action='read')
    allocate(line(ncols),header(ncols))
    write(strncols, '(i1)' ) ncols
    form1='('//strncols//'a16)'
    form2='('//strncols//'es16.8e2)'
    ! first pass: count the number of rows
    read(unit_number,trim(form1)) header
    num_rows = 0
    do
      read(unit_number,trim(form2),iostat=io_status) line
      if (io_status /= 0) exit
      num_rows = num_rows + 1
    end do
    ! allocate the array to hold the lines
    allocate(boundcond(num_rows,ncols))
    ! rewind the file to read the lines again
    rewind(unit_number)
    read(unit_number,trim(form1)) header
    ! second pass: read the lines into the array
    do i = 1, num_rows
      read(unit_number,trim(form2)) boundcond(i,:)
    end do
    deallocate(line,header)
    close(unit_number)
end subroutine read_file

!*****************************************************initial conditions*************************************************

subroutine initial_hydro(blk,i,j)
    !initial hydro condition
    !w(1:5) are rho, vx, vy, vz, and Etotal, respectively
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: pos(3),t,w(5),u(5),temp,egv,r0,rr,r
    r=blk%x_center(i)
    w(1)=max(diffused_rho(r),rho_floor)
    w(2)=diffused_vx(r)
    w(3:4)=0d0
    temp=diffused_temp(r)
    w(5)=prhot(w(1),temp)
    call assign_w_to_u_cell(blk,w,i,j)
end subroutine initial_hydro

subroutine initial_rad(blk,i,j)
    !initial radiation condition, it is simply Erad=aT^4
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: temp
    temp=blk%temp(i,j,1)
    blk%Erad(i,j,1)=a_rad*temp**4
end subroutine initial_rad

subroutine initial_scalar(blk,i,j)
    type(blockdef), pointer :: blk
    real(8) :: r,theta,rho,l_omega
    integer :: i,j,ip_am
end subroutine initial_scalar

function diffused_rho(r)
    !initial condition density
    real(8) :: diffused_rho,r,r0
    r0=n_domain(1)
    diffused_rho=diffused_rho0*(r/r0)**(-1.5)
end function diffused_rho

function diffused_vx(r)
    !initial condition vr
    real(8) :: diffused_vx,r
    diffused_vx=sqrt(2*ms*gr/r)
end function diffused_vx

function diffused_temp(r)
    !initial condition temperature
    real(8) :: diffused_temp,r,r0
    r0=n_domain(1)
    diffused_temp=diffused_temp0*(r/r0)**(-1d0)
end function diffused_temp

!*****************************************************boundary conditions*************************************************

subroutine boundary_hydro(blk,i,j)
    !hydro boundary condition, including inner and outer)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: r,theta,t,w(5),u(5),temp,egv
    t=time_sys%t
    r=blk%x_center(i)
    !the default is spherical
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
    !radiation boundary condition, including inner and outer
    !only effective when the corresponding rad_bound_type=9
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: v,r,r1,theta,temp,t
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
end subroutine boundary_rad

subroutine boundary_scalar(blk,i,j)
    !passive scalar, this function in this version of the code is no complete
    type(blockdef), pointer :: blk
    real(8) :: r,theta,rho,l_omega
    integer :: i,j,ip_am
end subroutine boundary_scalar

!*******************time dependent boundary conditions**************************************************************

function vej(t)
    !time dependent velocity of the ejecta
    real(8) :: vej,t
    call interpolation_linear(t,vej,time,v)
end function vej

function rhoej(t)
    !time dependent density of the ejecta
    real(8) :: rhoej,t
    call interpolation_linear(t,rhoej,time,rho)
end function rhoej

function tempej(t)
    !time dependent temperature of the ejecta
    real(8) :: tempej,t
    call interpolation_linear(t,tempej,time,temp)
end function tempej

!end of time dependent boundary conditions

subroutine time_dependent_bound_type()
    !called in src/boundary.f90, can be used to modify the boundary condition in computation
    real(8) :: t,temp_inner,v
    integer :: i
    t=time_sys%t
    if (t>tmax) then
        !when the ejecting phase is done, the inner hydro boundary condition changes to free
        !the inner radiation boundary condition changes to zero-gradient
        hydro_bound_type(1)=1
        rad_bound_type(1)=2
    else if (t<t_ram) then
        !t_ram is a buffer
        !cannot allow too hot ejecta come enters the computational domain immediately
        hydro_bound_type(1)=9
        rad_bound_type(1)=2
    else
        !injecting gas into the computational domain
        !during this ejecta launching phase, the inner boundary also has high temperature
        !radiation can diffusion from the inner boundary to the computational domain
        hydro_bound_type(1)=9
        rad_bound_type(1)=9
    end if
end subroutine time_dependent_bound_type

subroutine assemble_record_array(record_array)
    !called in rmhd.f90, it is used to record high frequency data, the time interval is dt_record
    !the length of record_array is record_length, and it must be larger than the actual data length
    type(blockdef), pointer :: blk
    real(8), allocatable :: record_array(:)
    real(8) :: sample_r(2),v
    procedure(field_calculator), pointer :: g1,g2
    integer :: i,error
    record_array=0
    record_array(1)=time_sys%t
    if (rank==np-1) then
        blk=>blk_tail
        if (iradiation==4) then
            record_array(2)=4d0*pi*blk%x_interface(blk_size_nx)**2*blk%Fradx(blk_size_nx,1,1)
        else
            record_array(2)=0d0
        end if
    end if
end subroutine assemble_record_array

subroutine problem_oper()
    !called in rmhd.f90, after each hydro and radiation transport step,
    !it can be used to modify data
end subroutine problem_oper

subroutine finalize_problem()
end subroutine finalize_problem

function user_kp(rho,Erad,tgas)
    !user defined planck mean opacity, when iopacity=3
    real(8) :: user_kp,rho,Erad,tgas
end function user_kp

function user_kr(rho,tgas)
    !user defined Rosseland mean opacity, when iopacity=3
    real(8) :: user_kr,rho,tgas
end function user_kr

end module problem

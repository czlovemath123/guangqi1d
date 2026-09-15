module problem
use phylib
use datastructure
use mathlib
use eos
implicit none
real(8), dimension(5), protected :: w_in,u_in,wl
real(8), protected :: rho,acc_rate,v(3),temp_in,p,egv_in,E0,m_planet,m_atm,t_atm,egvl,templ,p_thresh,bottom_state(4)
real(8), protected :: rho_disk,temp_disk,s_disk,tff,t_countdown,t_record,l_internal,kh_timescale,p_bottom_ratio,xtff
logical, protected :: keep_growing,thresh_reached


contains

subroutine initialize_problem()
    !calculate w, u, temp, and egv for all location
    real(8) :: vff,acc,rho_incoming,v_incoming,a(4),cs,pram,lacc
    integer :: isave,ierr
    logical :: file_exists
    namelist /rhd_quantities/ acc_rate,m_planet,kh_timescale,petsc_iter,petsc_qratio,   &
        p_bottom_ratio,xtff,lfld_heating,lchange_dust_opacity,kp_ratio,kr_ratio
    !test: rho left, u left, temp left, rho right, u right, temp right, final time
    open(unit=11,file=trim(path_root)//'/problem.data',status='old',action='read')
    read(unit=11,nml=rhd_quantities)
    close(11)
    path_out=trim(path_root)//trim(out_dir)
    lfld_mom=.true.                             !will calculate radiation flux
    temp_in=400d0                               !initial condition: ambient temperature
    kh_timescale=kh_timescale*year
    m_atm=1d-6                                  !initial condition: mass of the atmosphere in earth mass
    t_atm=2d3                                   !initial condition: atmosphere temperature
    m_planet=m_planet*mjupiter                  !plane mass in jupiter mass
    tff=pi/2*sqrt(n_domain(2)**3/gr/m_planet)   !free fall timescale at the outer boundary
    vff=sqrt(2*m_planet*gr/n_domain(2))         !free fall speed at the outer boundary
    acc_rate=acc_rate*mearth/year
    sig_rosseland_floor=1d-2/n_domain(2)        !lowerbound for the product of Rosseland mean opacity and density
    opacity_gas_rho_min=1d-13                   !min density of the tailored opacity table
    opacity_gas_rho_max=1d-4                    !max density of the tailored opacity table
    v_incoming=vff
    rho_incoming=acc_rate/4/pi/n_domain(2)**2/v_incoming
    w_in(1)=rho_incoming
    w_in(2)=-v_incoming
    w_in(3:4)=0d0
    p=prhot(rho_incoming,temp_in)               !gas pressure at the outer boundary
    w_in(5)=p
    call eos_wtou(w_in,u_in,temp_in,egv_in)
    E0=a_rad*temp_in**4                         !initial radiation energy density
    central_star%core%mass=m_planet
    m_atm=m_atm*mearth
    rho_disk=1d-11
    temp_disk=200d0
    s_disk=srhot_per_particle(rho_disk,temp_disk)
    pram=acc_rate/4d0/pi/n_domain(1)**2d0*sqrt(2d0*gr*m_planet/n_domain(1))
    p_thresh=pram*p_bottom_ratio
    if (lfld_heating) then
        l_internal=gr*m_planet**2d0/n_domain(1)/kh_timescale                  !impose a radiation flux at the inner boundary
    else
        l_internal=0d0
    end if
    lacc=acc_rate*gr*m_planet/n_domain(1)
    if (l_internal/lacc>0.3) then
        print *,'the planet internal luminosity is too large'
        stop
    end if
    !initialize all the cells including boundary cells
    if (rank==0) then
        write(*,'(A32,7ES18.6E2)') 'boundary velocity kms',v_incoming/1e5
        write(*,'(A32,7ES18.6E2)') 'incoming sound speed kms',adiabatic_cs(rho_incoming,temp_in)/1e5
        write(*,'(A32,7ES18.6E2)') 'boundary density',rho_incoming
        write(*,'(A32,7ES18.6E2)') 'primitive quantities',w_in
        write(*,'(A32,7ES18.6E2)') 'conserved quantities',u_in
        write(*,'(A32,7ES18.6E2)') 'entropy in disk',s_disk
        write(*,'(A32,7ES18.6E2)') 'estimated pram in bar',pram/1d6
        write(*,'(A32,7ES18.6E2)') 'bottom pressure in bar',p_thresh/1d6
        write(*,'(A32,7ES18.6E2)') 'internal luminosity in lsun',l_internal/lsun
        write(*,'(A32,7ES18.6E2)') 'accretion luminosity in lsun',lacc/lsun
    end if
    keep_growing=.true.
    if (restart) then
        call read_left_boundary(isave,bottom_state,file_exists)
        if (file_exists) then
            if (isave<=iframe) then
                keep_growing=.false.
                wl(1)=bottom_state(1)
                wl(2:4)=0d0
                wl(5)=bottom_state(2)
                templ=bottom_state(3)
                egvl=bottom_state(4)
                hydro_bound_type(1)=9
            else
                print *,'not supported now'
                stop
            end if
        else
            print *,'disallowed'
            stop
        end if
        if (rank==0) print *,'free fall timescale=',tff
    end if
end subroutine initialize_problem

subroutine time_dependent_bound_type()
    type(blockdef), pointer :: blk
    real(8) :: t,v
    integer :: ierr
    character(len=128) :: str
    t=time_sys%t
    if (keep_growing) then
        blk=>llist_head
        if (associated(blk,blk_head)) then
            if (blk%w(5,1,1,1)>p_thresh) then
                bottom_state(1)=blk%w(1,1,1,1)
                bottom_state(2)=blk%w(5,1,1,1)
                bottom_state(3)=blk%temp(1,1,1)
                bottom_state(4)=blk%egv(1,1,1)
                call save_left_boundary()
                wl(1)=bottom_state(1)
                wl(2:4)=0d0
                wl(5)=bottom_state(2)
                templ=bottom_state(3)
                egvl=bottom_state(4)
                keep_growing=.false.
                t_record=time_sys%t
                t_countdown=max(xtff*tff,t_record*0.5)
                print *,'t_ff=', tff
                print *,'t_record=',t_record
                print *,'t_countdown=',t_countdown
                print *,'left pressure threshold reached, start time countdown to relax the atmosphere'
                hydro_bound_type(1)=9
            end if
        end if
        nullify(blk)
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call mpi_bcast(t_countdown,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(t_record,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(keep_growing,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    else
        t_countdown=max(xtff*tff,t_record*0.5)
        if (time_sys%t-t_record>t_countdown) then
            print *,time_sys%t,t_record,t_countdown
            str='countdown finished, the postshock state should be in equilibrium now'
            call abort_guangqi(str)
        end if
    end if
    if (t>0.1*tff.and.lfld_heating) then
        blk=>llist_head
        if (blk%key(1)==1) then
            v=sum(blk%vol(1:3,1,1))*4*pi
            blk%fld_heating(1:3,1,1)=l_internal/v
        end if
    end if
end subroutine time_dependent_bound_type

subroutine initial_hydro(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: r,t,w(5),u(5),temp,egv,rho_atm,p_atm,r_in,r_out,v_atm,w_atm(5),u_atm(5),egv_atm
    !spread the initial atmosphere mass within bottom 0.1 jupiter radius
    r=blk%x_center(i)
    r_in=n_domain(1)
    r_out=r_in+1d-1*rjupiter
    v_atm=4d0/3d0*pi*(r_out**3-r_in**3)
    rho_atm=m_atm/v_atm
    if (r<r_out) then
        w_atm(1)=rho_atm
        w_atm(2:4)=0d0
        w_atm(5)=prhot(w_atm(1),t_atm)
        w=w_atm
    else
        w=w_in
    end if
    call assign_w_to_u_cell(blk,w,i,j)
end subroutine initial_hydro

subroutine initial_rad(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: temp
    if (iradiation==4) then
        temp=blk%temp(i,j,1)
        blk%Erad(i,j,1)=a_rad*temp**4d0
    end if
end subroutine initial_rad

subroutine save_left_boundary()
    if (rank==0) then
        open(unit=16,file=trim(path_root)//'/temp.data',status='replace',action='write')
        write(unit=16,fmt='(I5,4ES18.8E2)') iframe,bottom_state
        close(16)
    end if
end subroutine save_left_boundary

subroutine read_left_boundary(isave,a,file_exists)
    real(8) :: a(4)
    integer :: isave,ierr
    logical :: file_exists
    inquire(file=trim(path_root)//"/temp.data", exist=file_exists)
    if (file_exists) then
        open(unit=16,file=trim(path_root)//'/temp.data',status='old',action='read')
        read(unit=16,fmt='(I5,4ES18.8E2)') isave,a
        close(16)
    end if
end subroutine read_left_boundary

subroutine boundary_hydro(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: r,t,w(5),u(5),temp,egv
    r=blk%x_center(i)
    if (r<n_domain(1)) then
        w=wl
        temp=templ
    else if (r>n_domain(2)) then
        w(1:4)=w_in(1:4)
        temp=blk_tail%temp(blk_size_nx,1,1)
    end if
    w(5)=prhot(w(1),temp)
    call assign_w_to_u_cell(blk,w,i,j)
end subroutine boundary_hydro

subroutine boundary_rad(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
end subroutine boundary_rad

subroutine assemble_record_array(record_array)
    real(8), allocatable :: record_array(:)
end subroutine assemble_record_array

subroutine problem_oper()
end subroutine problem_oper

subroutine finalize_problem()
end subroutine finalize_problem

subroutine initial_scalar(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
end subroutine initial_scalar

subroutine boundary_scalar(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
end subroutine boundary_scalar

end module problem

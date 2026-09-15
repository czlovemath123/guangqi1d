module radiation_common_functions
use phylib
use mathlib
use communication
use datastructure
use eos
implicit none

interface fld_conductivity
    module procedure fld_conductivity_1d
end interface fld_conductivity

contains

subroutine generate_rosseland_planck_kap_table()
    call readin_gas_opacity_tables()
    call readin_dust_opacity()
end subroutine generate_rosseland_planck_kap_table

subroutine readin_gas_opacity_tables()
    real(8), dimension(:,:,:), allocatable :: kappa_rosseland,kappa_planck
    real(8), dimension(:,:), allocatable :: kappa_R_trim,kappa_P_trim
    real(8), dimension(:), allocatable :: tgas_table,pres_table
    real(8) :: metal,pres,rho,kR,kP,metal_list(3),tgas_list(94),pres_list(126)
    character(len=64) :: fmt1
    integer :: i,j,k,nt,npres,nmetal,tgas,trad,ipres_min,ipres_max,it_min,it_max,nt_table,npres_table
    nt=94
    npres=126
    nmetal=3
    fmt1='(1F5.1,1I8,2E10.4,2E13.6)'
    opacity_gas_pres_min=prhot(opacity_gas_rho_min,opacity_gas_t_min)
    opacity_gas_pres_max=prhot(opacity_gas_rho_max,opacity_gas_t_max)
    allocate(kappa_rosseland(nmetal,nt,npres),kappa_planck(nmetal,nt,npres))
    open(unit=11,file=trim(path_tables)//'/table1T.dat',status='old',action='read')
    do i=1,nmetal
        do j=1,nt
            do k=1,npres
                read(11,fmt=fmt1) metal,tgas,pres,rho,kR,kP
                kappa_rosseland(i,j,k)=kR
                kappa_planck(i,j,k)=kP
                if (i==metallicity.and.j==1) then
                    pres_list(k)=pres
                end if
            end do
            if (i==metallicity) then
                tgas_list(j)=tgas
            end if
        end do
        metal_list(i)=metal
    end do
    close(11)
    do i=1,npres-1
        if (pres_list(i+1)>opacity_gas_pres_min) then
            ipres_min=i
            exit
        end if
    end do
    do i=npres,2,-1
        if (pres_list(i-1)<opacity_gas_pres_max) then
            ipres_max=i
            exit
        end if
    end do
    npres_table=ipres_max-ipres_min+1
    do i=1,nt-1
        if (tgas_list(i+1)>opacity_gas_t_min) then
            it_min=i
            exit
        end if
    end do
    do i=nt,2,-1
        if (tgas_list(i-1)<opacity_gas_t_max) then
            it_max=i
            exit
        end if
    end do
    nt_table=it_max-it_min+1
    call create_table_2d(nt_table,npres_table,rosseland_gas_opacity_table)
    call create_table_2d(nt_table,npres_table,planck_gas_opacity_table)
    allocate(pres_table(npres_table),tgas_table(nt_table))
    allocate(kappa_R_trim(nt_table,npres_table),kappa_P_trim(nt_table,npres_table))
    do i=ipres_min,ipres_max
        pres_table(i-ipres_min+1)=pres_list(i)
    end do
    do i=it_min,it_max
        tgas_table(i-it_min+1)=tgas_list(i)
    end do
    do i=ipres_min,ipres_max
        do j=it_min,it_max
            kappa_R_trim(j-it_min+1,i-ipres_min+1)=max(kappa_rosseland(2,j,i),kr_floor)
            kappa_P_trim(j-it_min+1,i-ipres_min+1)=kappa_planck(2,j,i)
        end do
    end do
    kappa_P_trim=log10(kappa_P_trim)
    kappa_R_trim=log10(kappa_R_trim)
    pres_table=log10(pres_table)
    tgas_table=log10(tgas_table)
    rosseland_gas_opacity_table%t2d=kappa_R_trim
    rosseland_gas_opacity_table%xlist=tgas_table
    rosseland_gas_opacity_table%ylist=pres_table
    rosseland_gas_opacity_table%nd1=nt_table
    rosseland_gas_opacity_table%nd2=npres_table
    planck_gas_opacity_table%t2d=kappa_P_trim
    planck_gas_opacity_table%xlist=tgas_table
    planck_gas_opacity_table%ylist=pres_table
    planck_gas_opacity_table%nd1=nt_table
    planck_gas_opacity_table%nd1=npres_table
    deallocate(pres_table,tgas_table,kappa_R_trim,kappa_P_trim)
end subroutine readin_gas_opacity_tables

subroutine readin_dust_opacity()
    real(8), dimension(:,:), allocatable :: kappa_rosseland,kappa_planck
    real(8) :: rho,t,kp,kr,t_table(101),rho_table(61)
    integer :: nt,nrho,i,j,k
    character(len=64) :: fmt1
    nt=101
    nrho=61
    fmt1='(3ES18.5E2)'
    allocate(kappa_rosseland(nt,nrho),kappa_planck(nt,nrho))
    open(unit=11,file=trim(path_tables)//'/kP.dat',status='old',action='read')
    open(unit=12,file=trim(path_tables)//'/kR.dat',status='old',action='read')
    do i=1,nrho
        do j=1,nt
            read(11,fmt=fmt1) rho,t,kp
            read(12,fmt=fmt1) rho,t,kr
            if (i==1) then
                t_table(j)=t
            end if
            kappa_rosseland(j,i)=kr
            kappa_planck(j,i)=kp
        end do
        rho_table(i)=rho
    end do
    close(11)
    close(12)
    call create_table_2d(nt,nrho,rosseland_dust_opacity_table)
    call create_table_2d(nt,nrho,planck_dust_opacity_table)
    t_table=log10(t_table)
    rho_table=log10(rho_table)
    if (lchange_dust_opacity) kappa_rosseland=kappa_rosseland*kr_ratio
    kappa_rosseland=log10(kappa_rosseland)
    rosseland_dust_opacity_table%t2d=kappa_rosseland
    rosseland_dust_opacity_table%xlist=t_table
    rosseland_dust_opacity_table%ylist=rho_table
    rosseland_dust_opacity_table%nd1=nt
    rosseland_dust_opacity_table%nd2=nrho
    if (lchange_dust_opacity) kappa_planck=kappa_planck*kp_ratio
    kappa_planck=log10(kappa_planck)
    planck_dust_opacity_table%t2d=kappa_planck
    planck_dust_opacity_table%xlist=t_table
    planck_dust_opacity_table%ylist=rho_table
    planck_dust_opacity_table%nd1=nt
    planck_dust_opacity_table%nd1=nrho
end subroutine readin_dust_opacity

function fld_rosseland_opacity_interp(rho,tgas)
    real(8) :: rho,pres,tgas,fld_rosseland_opacity_interp,logrho,logpres,logtgas,kappa1,kappa2,kappa
#if     iopacity==1
    if (tgas<=200) then
        fld_rosseland_opacity_interp=1.87151d0
    else if (tgas>=1d6) then
        fld_rosseland_opacity_interp=1.085812E+01 
    else
        logtgas=log10(tgas)
        pres=prhot(rho,tgas)
        logpres=log10(pres)
        logrho=log10(rho)
        if (tgas<=1500d0) then
            !use dust table
            logrho=min(max(logrho,-18d0),-8d0)
            call interpolation2d_linear(logtgas,logrho,kappa1,rosseland_dust_opacity_table%xlist,rosseland_dust_opacity_table%ylist,  &
                rosseland_dust_opacity_table%t2d)
        else
            kappa1=-10d0
        end if
        if (tgas>=1100d0) then
            !use gas table
            logpres=max(logpres,-9d0)
            call interpolation2d_linear(logtgas,logpres,kappa2,rosseland_gas_opacity_table%xlist,rosseland_gas_opacity_table%ylist,  &
                rosseland_gas_opacity_table%t2d)
        else
            kappa2=-10d0
        end if
        kappa=max(kappa1,kappa2)
        fld_rosseland_opacity_interp=10d0**kappa
    end if
#elif   iopacity==2
    fld_rosseland_opacity_interp=const_opacity
#elif   iopacity==3
    fld_rosseland_opacity_interp=user_kr(rho,tgas)
#endif
end function fld_rosseland_opacity_interp

function fld_planck_opacity_interp_trad(rho,Erad,tgas)
    real(8) :: rho,pres,Erad,trad,tgas,fld_planck_opacity_interp_trad,logrho,logpres,logtrad,kappa,kappa1,kappa2
#if     iopacity==1
    trad=(Erad/a_rad)**0.25d0
    if (trad<=200) then
        fld_planck_opacity_interp_trad=2.77749d0
    else if (trad>=1d6) then
        fld_planck_opacity_interp_trad=7.693802E+00
    else
        logtrad=log10(trad)
        pres=prhot(rho,tgas)
        logpres=log10(pres)
        logrho=log10(rho)
        if (trad<1500d0) then
            !use dust table
            logrho=min(max(logrho,-18d0),-8d0)
            call interpolation2d_linear(logtrad,logrho,kappa1,planck_dust_opacity_table%xlist,planck_dust_opacity_table%ylist,  &
                planck_dust_opacity_table%t2d)
        else
            kappa1=-10d0
        end if
        if (trad>1100d0) then
            !use gas table
            logpres=max(logpres,-9d0)
            call interpolation2d_linear(logtrad,logpres,kappa2,planck_gas_opacity_table%xlist,planck_gas_opacity_table%ylist,  &
                planck_gas_opacity_table%t2d)
        else
            kappa2=-10d0
        end if
        kappa=max(kappa1,kappa2)
        kappa=min(kappa,2d0)
        fld_planck_opacity_interp_trad=min(10d0**kappa,kappa_planck_ceil1)
        if (pres>5d4) fld_planck_opacity_interp_trad=min(max(fld_planck_opacity_interp_trad,fld_rosseland_opacity_interp(rho,trad)),kappa_planck_ceil2)
    end if
#elif   iopacity==2
    fld_planck_opacity_interp_trad=const_opacity
#elif   iopacity==3
    fld_planck_opacity_interp_trad=user_kp(rho,Erad,tgas)
#endif
end function fld_planck_opacity_interp_trad

function fld_abs_opacity_interp_trad(rho,Erad,tgas)
    real(8) :: rho,pres,Erad,trad,tgas,fld_abs_opacity_interp_trad,logrho,logpres,logtrad,kappa
    trad=(Erad/a_rad)**0.25d0
    fld_abs_opacity_interp_trad=0d0
end function fld_abs_opacity_interp_trad

function fld_flux_limiter(r)
    real(8) :: fld_flux_limiter,r
    fld_flux_limiter=(2d0+r)/(6d0+3d0*r+r**2)
end function fld_flux_limiter

subroutine calculate_planck_rosseland_opacity_block(blk)
    type(blockdef), pointer :: blk
    real(8) :: Erad_cell,trad,tgas,rho,kappa_planck_trad,kappa_planck_tgas
    integer :: i,j,key(3),level
    blk%kappa_planck=0d0;blk%kappa_rosseland=0d0
    do i=blk_xlb,blk_xub
        Erad_cell=blk%Erad(i,1,1)
        tgas=blk%temp(i,1,1)
        trad=(Erad_cell/a_rad)**0.25d0
        rho=blk%w(1,i,1,1)
        kappa_planck_trad=fld_planck_opacity_interp_trad(rho,Erad_cell,tgas)
        blk%kappa_planck(i,1,1)=kappa_planck_trad
        blk%kappa_rosseland(i,1,1)=fld_rosseland_opacity_interp(rho,trad)
    end do
    if (associated(blk,blk_head)) then
        !left boundary block
        if (rad_bound_type(1)==9) then
            Erad_cell=blk%Erad(0,1,1)
            tgas=blk%temp(0,1,1)
            trad=(Erad_cell/a_rad)**0.25d0
            rho=blk%w(1,0,1,1)
            kappa_planck_trad=fld_planck_opacity_interp_trad(rho,Erad_cell,tgas)
            blk%kappa_planck(0,1,1)=kappa_planck_trad
            blk%kappa_rosseland(0,1,1)=fld_rosseland_opacity_interp(rho,trad)
        else
            blk%kappa_planck(0,1,1)=blk%kappa_planck(1,1,1)
            blk%kappa_rosseland(0,1,1)=blk%kappa_rosseland(1,1,1)
        end if
    end if
    if (associated(blk,blk_tail)) then
        !right boundary block
        if (rad_bound_type(2)==9) then
            Erad_cell=blk%Erad(blk_size_nx+1,1,1)
            tgas=blk%temp(blk_size_nx+1,1,1)
            trad=(Erad_cell/a_rad)**0.25d0
            rho=blk%w(1,blk_size_nx+1,1,1)
            kappa_planck_trad=fld_planck_opacity_interp_trad(rho,Erad_cell,tgas)
            blk%kappa_planck(blk_size_nx+1,1,1)=kappa_planck_trad
            blk%kappa_rosseland(blk_size_nx+1,1,1)=fld_rosseland_opacity_interp(rho,trad)
        else
            blk%kappa_planck(blk_size_nx+1,1,1)=blk%kappa_planck(blk_size_nx,1,1)
            blk%kappa_rosseland(blk_size_nx+1,1,1)=blk%kappa_rosseland(blk_size_nx,1,1)
        end if
    end if
#if     iopacity==2
    blk%kappa_planck=const_opacity
    blk%kappa_rosseland=const_opacity
#endif
end subroutine calculate_planck_rosseland_opacity_block

subroutine calculate_fld_mfp_sigma_block(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dr
    do j=blk_ylb,blk_yub
        do i=blk_xlb,blk_xub
            blk%sigma_planck(i,j,1)=blk%kappa_planck(i,j,1)*blk%w(1,i,j,1)
            blk%sigma_rosseland(i,j,1)=max(blk%kappa_rosseland(i,j,1)*max(blk%w(1,i,j,1),rho_thresh_petsc2),sig_rosseland_floor)
        end do
    end do
    !outer boundary for the spherical coordinate
#if     iopacity==1
    if (igeometry==2) then
        if (blk%key(1)==nx_blks*(2**blk%level)) then
            do j=blk_ylb,blk_yub
                blk%sigma_rosseland(blk_size_nx:blk_size_nx+1,j,1)=min(blk%sigma_rosseland(blk_size_nx,j,1),3d-1/n_domain(2))
            end do
        end if
    end if
#endif
    !blk%sigma_rosseland=1d7
end subroutine calculate_fld_mfp_sigma_block

subroutine applyboundconds_rad()
    call applyboundconds_rad_1d()
end subroutine applyboundconds_rad

subroutine applyboundconds_rad_1d()
    type(blockdef), pointer :: blk
    integer :: i
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        if (associated(blk,blk_head)) then
            call coor1lower_rad_boundary(blk)
        end if
        if (associated(blk,blk_tail)) then
            call coor1upper_rad_boundary(blk)
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%next
        end if
    end do
end subroutine applyboundconds_rad_1d

subroutine coor1lower_rad_boundary(blk)
    type(blockdef), pointer :: blk
    real(8) :: pos(3),t,Frad(3),feddington
    real(8), dimension(:), allocatable :: Erad
    integer :: ijk(3),i,j
    character(len=128) :: alert
    select case(rad_bound_type(1))
    case(9)     !user specified boundary condition
        call bound_rad(blk,0,1)
    case(8)     !user specified boundary condition
        call bound_rad(blk,0,1)
    case(1)     !transmissive boundary condition
        alert="1d transimissive lower radiation boundary not implemented"
        call abort_guangqi(alert)
    case(2)     !reflective boundary condition
        blk%Erad(0,1,1)=blk%Erad(1,1,1)
    case default
    end select
end subroutine coor1lower_rad_boundary

subroutine coor1upper_rad_boundary(blk)
    type(blockdef), pointer :: blk
    real(8) :: pos(3),t,Frad(3),feddington
    real(8), dimension(:), allocatable :: Erad
    integer :: ijk(3),i,j
    character(len=128) :: alert
    select case(rad_bound_type(2))
    case(9)     !user specified boundary condition
        call bound_rad(blk,blk_size_nx+1,1)
    case(8)     !user specified boundary condition
        call bound_rad(blk,blk_size_nx+1,1)
    case(1)     !transmissive boundary condition
        !defined implicitly in radiation_common_functions.f90 calculate_fld_conductivity
        blk%Erad_next=blk%Erad(blk_size_nx,1,1)
    case(2)     !reflective boundary condition
        blk%Erad(blk_size_nx+1,1,1)=blk%Erad(blk_size_nx,1,1)
    case default
    end select
end subroutine coor1upper_rad_boundary

subroutine calculate_fld_conductivity()
    type(blockdef), pointer :: blk
    integer :: i
    call communicate_fld()
    call applyboundconds_rad()
    call blk_traversal(calculate_fld_conductivity_1d)
    !call check_conductivity()
end subroutine calculate_fld_conductivity

subroutine calculate_fld_conductivity_1d(blk)
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,xi,dx,Erad_l,Erad_u,sigma_l,sigma_u
    integer :: i,j
    do i=1,np_nblk(rank+1)
        call fld_conductivity_left(blk)
        call fld_conductivity_interior(blk)
        call fld_conductivity_right(blk)
    end do
end subroutine calculate_fld_conductivity_1d

subroutine fld_conductivity_left(blk)
    !first cell of this block
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,xi,Erad1,Erad2,sg1,sg2,kr1,kr2,gradErad,sigma,Erad,dx,dx_nb
    if (associated(blk,blk_head)) then
        if (rad_bound_type(1)==2.or.rad_bound_type(1)==8) then
            blk%kx(0,1,1)=c_light/blk%sigma_rosseland(1,1,1)/3d0
            blk%krx_inter(0,1,1)=0d0
        else if (rad_bound_type(1)==9) then
            x1=blk%x_center(0)
            x2=blk%x_center(1)
            xi=blk%x_interface(0)
            Erad1=blk%Erad(0,1,1)
            Erad2=blk%Erad(1,1,1)
            Erad=Erad1+(xi-x1)/(x2-x1)*(Erad2-Erad1)
            gradErad=(Erad2-Erad1)/(x2-x1)
            sg1=blk%sigma_rosseland(0,1,1)
            sg2=blk%sigma_rosseland(1,1,1)
            sigma=2d0*sg1*sg2/(sg1+sg2)
            blk%kx(0,1,1)=fld_conductivity(gradErad,sigma,Erad)
            kr1=blk%kappa_rosseland(0,1,1)
            kr2=blk%kappa_rosseland(1,1,1)
            blk%krx_inter(0,1,1)=2d0*kr1*kr2/(kr1+kr2)
        end if
    else
        xi=blk%x_interface(0)
        x1=blk%xl_guard
        x2=blk%x_center(1)
        Erad1=blk%Erad(0,1,1)
        Erad2=blk%Erad(1,1,1)
        Erad=Erad1+(xi-x1)/(x2-x1)*(Erad2-Erad1)
        gradErad=(Erad2-Erad1)/(x2-x1)
        sg1=blk%sigma_rosseland(0,1,1)
        sg2=blk%sigma_rosseland(1,1,1)
        sigma=2d0*sg1*sg2/(sg1+sg2)
        blk%kx(0,1,1)=fld_conductivity(gradErad,sigma,Erad)
        kr1=blk%kappa_rosseland(0,1,1)
        kr2=blk%kappa_rosseland(1,1,1)
        blk%krx_inter(0,1,1)=2d0*kr1*kr2/(kr1+kr2)
    end if
end subroutine fld_conductivity_left

subroutine fld_conductivity_interior(blk)
    !interior cells
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,xi,Erad1,Erad2,sg1,sg2,kr1,kr2,gradErad,sigma,Erad
    integer :: j
    do j=1,blk_size_nx-1
        x1=blk%x_center(j)
        x2=blk%x_center(j+1)
        xi=blk%x_interface(j)
        Erad1=blk%Erad(j,1,1)
        Erad2=blk%Erad(j+1,1,1)
        Erad=Erad1+(xi-x1)/(x2-x1)*(Erad2-Erad1)
        gradErad=(Erad2-Erad1)/(x2-x1)
        sg1=blk%sigma_rosseland(j,1,1)
        sg2=blk%sigma_rosseland(j+1,1,1)
        sigma=2d0*sg1*sg2/(sg1+sg2)
        blk%kx(j,1,1)=fld_conductivity(gradErad,sigma,Erad)
        kr1=blk%kappa_rosseland(j,1,1)
        kr2=blk%kappa_rosseland(j+1,1,1)
        blk%krx_inter(j,1,1)=2d0*kr1*kr2/(kr1+kr2)
    end do
end subroutine fld_conductivity_interior

subroutine fld_conductivity_right(blk)
    !last cell of this block
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,xi,Erad1,Erad2,sg1,sg2,kr1,kr2,gradErad,sigma,Erad,dx,dx_nb
    if (associated(blk,blk_tail)) then
        if (rad_bound_type(2)==2) then
            blk%kx(blk_size_nx,1,1)=c_light/blk%sigma_rosseland(blk_size_nx,1,1)/3d0
            blk%krx_inter(blk_size_nx,1,1)=0d0
        else if (rad_bound_type(2)==9) then
        else if (rad_bound_type(2)==1) then
            x1=blk%x_center(blk_size_nx)
            x2=blk%x_center(blk_size_nx+1)
            xi=blk%x_interface(blk_size_nx)
            Erad1=blk%Erad(blk_size_nx,1,1)
            if (igeometry==0) then
                Erad2=blk%Erad(blk_size_nx,1,1)/(1d0+3d0*blk%sigma_rosseland(blk_size_nx,1,1)*dx/2d0)
                blk%Erad(blk_size_nx+1,1,1)=Erad2
            else if (igeometry==2) then
                Erad2=blk%Erad(blk_size_nx,1,1)*x1**2/x2**2
                blk%Erad(blk_size_nx+1,1,1)=Erad2
            end if
            Erad=Erad1+(xi-x1)/(x2-x1)*(Erad2-Erad1)
            gradErad=(Erad2-Erad1)/(x2-x1)
            sg1=blk%sigma_rosseland(blk_size_nx,1,1)
            sg2=blk%sigma_rosseland(blk_size_nx,1,1)
            sigma=min((sg1+sg2)/2d0,2d0*sg1*sg2/(sg1+sg2))
            blk%kx(blk_size_nx,1,1)=fld_conductivity(gradErad,sigma,Erad)
            kr1=blk%kappa_rosseland(blk_size_nx,1,1)
            kr2=blk%kappa_rosseland(blk_size_nx+1,1,1)
            blk%krx_inter(blk_size_nx,1,1)=2d0*kr1*kr2/(kr1+kr2)
        end if
    else
        x1=blk%x_center(blk_size_nx)
        x2=blk%xu_guard
        xi=blk%x_interface(blk_size_nx)
        Erad1=blk%Erad(blk_size_nx,1,1)
        Erad2=blk%Erad(blk_size_nx+1,1,1)
        Erad=Erad1+(xi-x1)/(x2-x1)*(Erad2-Erad1)
        gradErad=(Erad2-Erad1)/(x2-x1)
        sg1=blk%sigma_rosseland(blk_size_nx,1,1)
        sg2=blk%sigma_rosseland(blk_size_nx+1,1,1)
        sigma=min((sg1+sg2)/2d0,2d0*sg1*sg2/(sg1+sg2))
        blk%kx(blk_size_nx,1,1)=fld_conductivity(gradErad,sigma,Erad)
        kr1=blk%kappa_rosseland(blk_size_nx,1,1)
        kr2=blk%kappa_rosseland(blk_size_nx+1,1,1)
        blk%krx_inter(blk_size_nx,1,1)=2d0*kr1*kr2/(kr1+kr2)
    end if
end subroutine fld_conductivity_right

function fld_conductivity_1d(gradErad,sigma,Erad)
    real(8) :: fld_conductivity_1d,gradErad,sigma,Erad,r
    r=abs(gradErad)/sigma/Erad
    fld_conductivity_1d=c_light*fld_flux_limiter(r)/sigma
end function fld_conductivity_1d

subroutine fld_irradiation_1d(blk)
    !calcualte the intensity induced by the irradiation
    !used to calculate the radiation force and thermal heating
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: tau,rho,kappa,r
    do i=1,blk_size_nx
        rho=blk%w(1,i,1,1)
        kappa=blk%kappa_planck(i,1,1)
        tau=(blk%tau(i-1,1,1)+blk%tau(i,1,1))/2
        r=blk%x_center(i)
        blk%irrad(i,1,1)=kappa*rho*blk%irrad0(1)*exp(-tau)*(n_domain(1)/r)**2d0
    end do
end subroutine fld_irradiation_1d

subroutine fld_radhydro_boost_heating()
    type(blockdef), pointer :: blk
    integer :: i,j,iblk,error
    if (mod((time_sys%ntimestep+1),radhydro_boost)==0) then
        blk=>llist_head
        do iblk=1,np_nblk(rank+1)
            blk%fld_heating=0d0
            if (lradhydro_boost) blk%fld_heating=blk%fld_heating+blk%compres_heat_acculm/time_sys%dt_radhydro
            blk=>blk%next
        end do
    end if
end subroutine fld_radhydro_boost_heating

subroutine rad_source(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: s1,s2,flux1,flux2,vol,dx,dy,sx1,sx2,sy1,sy2,xflux1,xflux2,yflux1,yflux2
    blk%Erad_source=0d0
    if (lirradiation) blk%Erad_source=blk%Erad_source+blk%irrad
    if (lfld_heating) blk%Erad_source=blk%Erad_source+blk%fld_heating
    if (lother_heating) blk%Erad_source=blk%Erad_source+blk%other_heating
    if (lrad_adv) then
        blk%erad_adv=0d0
        if (igeometry==0) then
            dx=blk%dxyz(1)
            do i=1,blk_size_nx
                flux1=blk%erad_xflux(i-1,1,1)
                flux2=blk%erad_xflux(i,1,1)
                blk%erad_adv(i,1,1)=(flux1-flux2)/dx
            end do
        else if (igeometry==2) then
            do i=1,blk_size_nx
                s1=blk%surf1(i-1,1,1)
                s2=blk%surf1(i,1,1)
                flux1=blk%erad_xflux(i-1,1,1)
                flux2=blk%erad_xflux(i,1,1)
                vol=blk%vol(i,1,1)
                blk%erad_adv(i,1,1)=(s1*flux1-s2*flux2)/vol
            end do
        end if
    end if
end subroutine rad_source

subroutine assign_beta_temp_D_1d(blk,i,dt,beta,temp,D_xl,D_xu)
    type(blockdef), pointer :: blk
    integer :: i
    real(8) :: dt,temp,beta,D_xl,D_xu
    if (blk%w(1,i,1,1)<rho_thresh_petsc1) then
        beta=0d0
    else
        beta=blk%sigma_planck(i,1,1)*c_light*dt
    end if
    temp=blk%temp(i,1,1)
    D_xl=blk%kx(i-1,1,1)
    D_xu=blk%kx(i,1,1)
end subroutine assign_beta_temp_D_1d

subroutine fld_pradgradv_1d(blk)
    type(blockdef), pointer :: blk
    integer :: i
    real(8) :: dt,pradgradv,gdv(3),f,lambda,r,Erad,gradErad,sigma,x1,x2,x3,vel1,vel2,vel3,e1,e2,e3,chi
    if (igeometry==2) then
        i=1
        x1=blk%xl_guard
        x2=blk%x_center(i)
        x3=blk%x_center(i+1)
        gdv(1)=(blk%w(2,i+1,1,1)-blk%w(2,i-1,1,1))/(x3-x1)
        gdv(2)=blk%w(2,i,1,1)/x2
        gdv(3)=blk%w(2,i,1,1)/x2
        if (blk%key(1)==1) then
            if (rad_bound_type(1)==8.or.rad_bound_type(1)==2) then
                gradErad=(blk%Erad(i+1,1,1)-blk%Erad(i,1,1))/(x3-x1)
            else if (rad_bound_type(1)==9) then
                gradErad=(blk%Erad(i+1,1,1)-blk%Erad(i-1,1,1))/(x3-x1)
            end if
        else
            gradErad=(blk%Erad(i+1,1,1)-blk%Erad(i-1,1,1))/(x3-x1)
        end if
        Erad=blk%Erad(i,1,1)
        sigma=blk%sigma_rosseland(i,1,1)
        r=abs(gradErad)/sigma/Erad
        lambda=fld_flux_limiter(r)
        f=lambda+(lambda*r)**2
        pradgradv=(f*gdv(1)+(1-f)/2*gdv(2)+(1-f)/2*gdv(3))
        blk%fld_pradgradv(i,1,1)=pradgradv
        do i=2,blk_size_nx-1
            x1=blk%x_center(i-1)
            x2=blk%x_center(i)
            x3=blk%x_center(i+1)
            gdv(1)=(blk%w(2,i+1,1,1)-blk%w(2,i-1,1,1))/(x3-x1)
            gdv(2)=blk%w(2,i,1,1)/x2
            gdv(3)=blk%w(2,i,1,1)/x2
            !direct approach
            gradErad=(blk%Erad(i+1,1,1)-blk%Erad(i-1,1,1))/(x3-x1)
            Erad=blk%Erad(i,1,1)
            sigma=blk%sigma_rosseland(i,1,1)
            r=abs(gradErad)/sigma/Erad
            lambda=fld_flux_limiter(r)
            f=lambda+(lambda*r)**2
            pradgradv=(f*gdv(1)+(1-f)/2*gdv(2)+(1-f)/2*gdv(3))
            blk%fld_pradgradv(i,1,1)=pradgradv
        end do
        i=blk_size_nx
        x1=blk%x_center(i-1)
        x2=blk%x_center(i)
        x3=blk%xu_guard
        gdv(1)=(blk%w(2,i+1,1,1)-blk%w(2,i-1,1,1))/(x3-x1)
        gdv(2)=blk%w(2,i,1,1)/x2
        gdv(3)=blk%w(2,i,1,1)/x2
        gradErad=(blk%Erad(i+1,1,1)-blk%Erad(i-1,1,1))/(x3-x1)
        Erad=blk%Erad(i,1,1)
        sigma=blk%sigma_rosseland(i,1,1)
        r=abs(gradErad)/sigma/Erad
        lambda=fld_flux_limiter(r)
        f=lambda+(lambda*r)**2
        pradgradv=(f*gdv(1)+(1-f)/2*gdv(2)+(1-f)/2*gdv(3))
        blk%fld_pradgradv(i,1,1)=pradgradv
    else
    end if
end subroutine fld_pradgradv_1d

subroutine calculate_fld_pradgradv(blk)
    type(blockdef), pointer :: blk
    call fld_pradgradv_1d(blk)
end subroutine calculate_fld_pradgradv

subroutine assign_g_alpha_1d(blk,dt,i,gxl,gxu,alpha_xl,alpha_xu)
    type(blockdef), pointer :: blk
    real(8) :: dt,gxl,gxu,alpha_xl,alpha_xu
    integer :: i
    gxl=blk%gx1d(1,i)
    gxu=blk%gx1d(2,i)
    alpha_xl=dt*blk%alphax1d(1,i)
    alpha_xu=dt*blk%alphax1d(2,i)
end subroutine assign_g_alpha_1d

subroutine erad_int_on_boundaries(blk)
    !calculate erad_int on boundaries, prepare for the frad calculations
    type(blockdef), pointer :: blk,blk_temp
    real(8) :: q,x1,x2
    integer :: i,j
    if (blk%key(1)==1) then
        if (rad_bound_type(1)==2) blk%Erad_int(0,1,1)=blk%Erad_int(1,1,1)
        if (rad_bound_type(1)==8) blk%Erad_int(0,1,1)=blk%Erad_int(1,1,1)
        if (rad_bound_type(1)==9) blk%Erad_int(0,1,1)=blk%Erad_int(1,1,1)
    end if
    if (blk%key(1)==nx_blks*2**blk%level) then
        if (rad_bound_type(2)==2) blk%Erad_int(blk_size_nx+1,1,1)=blk%Erad_int(blk_size_nx,1,1)
        if (rad_bound_type(2)==1) then
            x1=blk%x_center(blk_size_nx)
            x2=blk%x_center(blk_size_nx+1)
            if (igeometry==0) then
                q=1.5d0*blk%sigma_rosseland(blk_size_nx,1,1)*(x2-x1)
                blk%Erad_int(blk_size_nx+1,1,1)=blk%Erad_int(blk_size_nx,1,1)/(1+q)
            else if (igeometry==2) then
                blk%Erad_int(blk_size_nx+1,1,1)=blk%Erad_int(blk_size_nx,1,1)*x1**2d0/x2**2d0
            end if
        end if
    end if
end subroutine erad_int_on_boundaries

subroutine frad_for_momentum(blk)
    !calculate fradx and frady on edges and inside
    type(blockdef), pointer :: blk,blk_temp
    real(8) :: x1,x2,y1,y2,gradEx
    integer :: i,j
    x1=blk%xl_guard
    x2=blk%x_center(1)
    gradEx=(blk%Erad_int(1,1,1)-blk%Erad_int(0,1,1))/(x2-x1)
    blk%Fradx(0,1,1)=-blk%kx(0,1,1)*gradEx
    do i=1,blk_size_nx-1
        x1=blk%x_center(i)
        x2=blk%x_center(i+1)
        gradEx=(blk%Erad_int(i+1,1,1)-blk%Erad_int(i,1,1))/(x2-x1)
        blk%Fradx(i,1,1)=-blk%kx(i,1,1)*gradEx
    end do
    x1=blk%x_center(blk_size_nx)
    x2=blk%xu_guard
    gradEx=(blk%Erad_int(blk_size_nx+1,1,1)-blk%Erad_int(blk_size_nx,1,1))/(x2-x1)
    blk%Fradx(blk_size_nx,1,1)=-blk%kx(blk_size_nx,1,1)*gradEx
end subroutine frad_for_momentum

subroutine fld_aradx_arady(blk)
    !radiation acceleration is optional on top of the radiation flux calculation
    type(blockdef), pointer :: blk
    real(8) :: dx,gradEx,x1,x2,y1,y2,q,g,u(5),r,egv,dt,de,aradx,arady
    integer :: i,j,ijk(3)
    real(8), allocatable :: arad(:,:)
    do i=0,blk_size_nx
        blk%aradx(i,1,1)=blk%krx_inter(i,1,1)*blk%Fradx(i,1,1)/c_light
    end do
    dt=time_sys%dt_radhydro
    do i=1,blk_size_nx
        u=blk%u(1:5,i,1,1)
        r=blk%x_center(i)
        !aradx=(blk%aradx(i-1,1,1)+blk%aradx(i,1,1))/2
        aradx=blk%sigma_rosseland(i,1,1)/blk%w(1,i,1,1)*(blk%Fradx(i-1,1,1)+blk%Fradx(i,1,1))/2/c_light
        blk%u(1:5,i,1,1)=u
        blk%u(2,i,1,1)=u(2)+u(1)*aradx*dt
        blk%u(5,i,1,1)=u(5)+u(2)*aradx*dt
        blk%w(2,i,1,1)=blk%w(2,i,1,1)+aradx*dt
        !blk%fld_aradwork(i,1,1)=u(2)*aradx*dt
    end do
    !call convert_u_to_w_block(blk)
end subroutine fld_aradx_arady

end module radiation_common_functions

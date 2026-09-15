module muscl
use gravity
use hydro
use eos
use datastructure
use communication
use recon_evolve
use boundary
use source_control
!use viscous
implicit none

real(8), protected :: t_muscl_max

contains

function estimate_block_dt_hydro(blk)
    !estimate dt_hydro on all interfaces, save the minimum value to estimate_block_dt_hydro
    type(blockdef), pointer :: blk,blk_temp
    real(8), dimension(:,:,:), allocatable :: varray,varray_x,varray_y,vmaxarray
    real(8), allocatable :: dty(:),dtx(:),dty2d(:,:)
    real(8) :: vmax,cs,rho,temp,dt(nd),vblockmax,vlocalmax,temp2(2),egv2(2),flux(5),wl(5),wr(5)
    character(len=1) :: dir
    integer :: i,j,k,ierr
    real(8) :: estimate_block_dt_hydro
    call allocate_cell_data_block(varray_x)
    do i=1,blk_size_nx
        varray_x(i,1,1)=max(blk%hllc_vx(i-1,1,1),blk%hllc_vx(i,1,1))
    end do
    if (igeometry==0) then
        dt(1)=blk%dxyz(1)*CFL/maxval(varray_x)
        estimate_block_dt_hydro=dt(1)
    else if (igeometry==1.or.igeometry==2) then
        if (llnx) then
            allocate(dtx(blk_size_nx))
            do i=1,blk_size_nx
                dtx(i)=blk%dr(i)*CFL/varray_x(i,1,1)
            end do
            dt(1)=minval(dtx)
            deallocate(dtx)
        else
            dt(1)=blk%dxyz(1)*CFL/maxval(varray_x)
        end if
        estimate_block_dt_hydro=dt(1)
    end if
    deallocate(varray_x)
end function estimate_block_dt_hydro

subroutine vanleer_hydro_unsplit()
    !based on Stone & Gardiner 2009 "A simple unsplit Godunov method for multidimensional MHD"
    type(blockdef), pointer :: blk
    real(8) :: dt_temp,w(5)
    integer :: iblk,i,j,k,key(3),ierr
    if (lradhydro_boost) call blk_traversal(initialize_radhydro_boost)
    call applyboundconds()
    call communicate_hydro()
    call blk_traversal(muscl_blk_initialize)
    call blk_traversal(reconstruct_hydro)
    call blk_traversal(hllc_muscl)
    call communicate_flux()
    if (lirradiation) call calculate_irrad_tau()
    call collective_sub(estimate_block_dt_hydro,op_min,op_min,time_sys%dt_hydro)
    call blk_traversal(vanleer_predictor)
    call applyboundconds()
    call communicate_hydro()
    if (lrad_adv) call communicate_fld()
    call blk_traversal(reconstruct_hydro)
    call blk_traversal(hllc_muscl)
    call communicate_flux()
    if (lirradiation) call calculate_irrad_tau()
    call blk_traversal(vanleer_corrector)
    if (lradhydro_boost) call blk_traversal(radhydro_boost_heating)
end subroutine vanleer_hydro_unsplit

subroutine muscl_blk_initialize(blk)
    type(blockdef), pointer :: blk
    blk%u0=blk%u
    if (lpassive) blk%passive_scalar0=blk%passive_scalar
    if (lam_con) blk%omega0=blk%omega
end subroutine muscl_blk_initialize

subroutine vanleer_predictor(blk)
    !predictor integrator
    type(blockdef), pointer :: blk
    real(8) :: dt
    integer :: i,j,k,ip_am,key(3)
    dt=time_sys%dt_hydro
    call hydro_conservation_law(blk,dt*half)
    call geometric_source(blk,dt*half,predictor)
    call gravity_source(blk,dt*half,predictor)
    call convert_u_to_w_block(blk)
    blk%u1=blk%u
end subroutine vanleer_predictor

subroutine vanleer_corrector(blk)
    !correction integrator
    type(blockdef), pointer :: blk
    real(8) :: dt
    integer :: i,j,k,ip_am
    dt=time_sys%dt_hydro
    call hydro_conservation_law(blk,dt)
    call geometric_source(blk,dt,corrector)
    call gravity_source(blk,dt,corrector)
    call convert_u_to_w_block(blk)
end subroutine vanleer_corrector

subroutine hydro_conservation_law(blk, dt)
    type(blockdef), pointer :: blk 
    real(8), intent(in) :: dt
    real(8) :: sx1, sx2, sy1, sy2, vol, inv_vol
    real(8), dimension(5) :: u, xflux1, xflux2, yflux1, yflux2
    integer :: i, j
    
    ! Precompute pointers to frequently accessed arrays
    real(8), dimension(:,:,:), pointer :: surf1, surf2, vol_arr
    real(8), dimension(:,:,:,:), pointer :: u0, u_arr, xflux, yflux
    
    surf1 => blk%surf1
    surf2 => blk%surf2
    vol_arr => blk%vol
    u0 => blk%u0
    u_arr => blk%u
    xflux => blk%xflux
    yflux => blk%yflux
    
    do i = 1, blk_size_nx
        u = u0(1:5, i, 1, 1)
        sx1 = surf1(i-1, 1, 1)
        sx2 = surf1(i, 1, 1)
        xflux1 = xflux(1:5, i-1, 1, 1)
        xflux2 = xflux(1:5, i, 1, 1)
        vol = vol_arr(i, 1, 1)
        
        ! Precompute inverse volume to avoid division in loop
        inv_vol = dt / vol
        u = u + (xflux1 * sx1 - xflux2 * sx2) * inv_vol
        
        u_arr(1:5, i, 1, 1) = u
    end do
end subroutine hydro_conservation_law

subroutine geometric_source(blk, dt, step)
    type(blockdef), pointer :: blk
    real(8), intent(in) :: dt
    real(8) :: temp, egv, p, rho, r, vr, vtheta, vol, cot, ftheta_omega, fr_omega, omega, theta, fromega, fthetaomega, &
               mplanet, fr, ftheta, egv1, pxl, pxr, pyl, pyr, ril, rir, thetail, thetair, dpx, dpy, slpx, slpy, prs, pthetas, rhox, &
               omegax, vrx, sx1, sx2, sy1, sy2, inv_vol, r_sq, sin_theta, cos_theta, omega_sq, r_omega_sq, &
               ril_sq, rir_sq, ril_cu, rir_cu  ! Added missing variable declarations
    real(8), dimension(5) :: u, u1, u2, w
    integer :: i, j, ip_omega
    character(len=16), intent(in) :: step

    ! Precompute pointers to frequently accessed arrays
    real(8), dimension(:,:,:,:), pointer :: u_arr, u_muscl_arr, u0_arr, u1_arr
    real(8), dimension(:,:,:), pointer :: omega_muscl_arr, omega0_arr, omega1_arr
    real(8), dimension(:), pointer :: x_center, x_interface, y_center, y_interface
    real(8), dimension(:,:,:), pointer :: vol_arr, fr_omega_arr, ftheta_omega_arr
    real(8), dimension(:,:,:,:), pointer :: xslp_arr, yslp_arr

    ! Assign pointers
    u_arr => blk%u
    u0_arr => blk%u0
    u1_arr => blk%u1
    u_muscl_arr => blk%u_muscl
    x_center => blk%x_center
    x_interface => blk%x_interface
    vol_arr => blk%vol
    xslp_arr => blk%xslp

    if (igeometry == 2) then
        ! Set u_muscl based on step
        if (step == predictor) then
            u_muscl_arr = u0_arr
        else if (step == corrector) then
            u_muscl_arr = u1_arr
        end if

        do i = 1, blk_size_nx
            u1 = u_arr(1:5, i, 1, 1)
            u = u_muscl_arr(1:5, i, 1, 1)
            call eos_utow(u, w, temp, egv)
            p = w(5)
            slpx = xslp_arr(5, i, 1, 1)
            r = x_center(i)
            ril = x_interface(i-1)
            rir = x_interface(i)
            vol = vol_arr(i, 1, 1)
            inv_vol = 1.0d0 / vol

            ! Precompute common terms
            ril_sq = ril * ril
            rir_sq = rir * rir
            ril_cu = ril * ril_sq
            rir_cu = rir * rir_sq
            r_sq = r * r

            fr = (p * (-ril_sq + rir_sq) + ((-2.0d0 * ril_cu + 2.0d0 * rir_cu + 3.0d0 * r * (ril_sq - rir_sq)) * slpx) / 3.0d0) * inv_vol

            u1(2) = u1(2) + fr * dt
            u_arr(1:5, i, 1, 1) = u1
        end do
    end if
end subroutine geometric_source


subroutine initialize_radhydro_boost(blk)
    !reset the iso_temp and compres_heat_acculm
    type(blockdef), pointer :: blk
    integer :: i,j
    if (mod(time_sys%ntimestep,radhydro_boost)==0) then
        blk%iso_temp=blk%temp
        blk%compres_heat_acculm=0d0
    end if
end subroutine initialize_radhydro_boost

subroutine radhydro_boost_heating(blk)
    !accumulate the compressional heating source
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: temp,rho,w(5),egv,egv_iso_temp
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            temp=blk%iso_temp(i,j,1)
            w=blk%w(1:5,i,j,1)
            rho=w(1)
            egv=blk%egv(i,j,1)
            egv_iso_temp=egvrhot(w(1),temp)
            if (egv>egv_iso_temp) then
                blk%compres_heat_acculm(i,j,1)=blk%compres_heat_acculm(i,j,1)+egv-egv_iso_temp
                w(5)=prhot(rho,temp)
                call assign_w_to_u_cell(blk,w,i,j)
            end if
        end do
    end do
end subroutine radhydro_boost_heating

end module muscl

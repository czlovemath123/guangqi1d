module muscl
use gravity
use hydro
use eos
use datastructure
use communication
use recon_evolve
use boundary
use source_control
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
    call applyboundconds()
    call communicate_hydro()
    call blk_traversal(muscl_blk_initialize)
    call blk_traversal(reconstruct_hydro)
    call blk_traversal(hllc_muscl)
    call communicate_flux()
    call collective_sub(estimate_block_dt_hydro,op_min,op_min,time_sys%dt_hydro)
    call blk_traversal(vanleer_predictor)
    call applyboundconds()
    call communicate_hydro()
    if (lrad_adv) call communicate_fld()
    call blk_traversal(reconstruct_hydro)
    call blk_traversal(hllc_muscl)
    call communicate_flux()
    call blk_traversal(vanleer_corrector)
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

subroutine hydro_conservation_law(blk,dt)
    type(blockdef), pointer :: blk
    real(8) :: temp,egv,p,rho,r,dt,sx1,sx2,sy1,sy2,vol,egv_floor
    real(8), dimension(5) :: u,xflux1,xflux2,yflux1,yflux2
    integer :: i,j
    do i=1,blk_size_nx
        u=blk%u0(1:5,i,1,1)
        sx1=blk%surf1(i-1,1,1)
        sx2=blk%surf1(i,1,1)
        xflux1=blk%xflux(1:5,i-1,1,1)
        xflux2=blk%xflux(1:5,i,1,1)
        vol=blk%vol(i,1,1)
        u=u+(xflux1*sx1-xflux2*sx2)*dt/vol
        blk%u(1:5,i,1,1)=u
    end do
end subroutine hydro_conservation_law

subroutine geometric_source(blk,dt,step)
    type(blockdef), pointer :: blk
    real(8) :: temp,egv,p,rho,r,dt,vr,vtheta,sx1,sx2,sy1,sy2,vol,cot,ftheta_omega,fr_omega,omega,theta,mplanet,fr,  &
        ftheta,egv1,pxl,pxr,ril,rir,slpx
    real(8), dimension(5) :: u,u1,u2,w
    integer :: i,j,ip_omega
    character(len=16) :: step
    if (igeometry==1) then
    else if (igeometry==2) then
        if (step==predictor) then
            blk%u_muscl=blk%u0
        else if (step==corrector) then
            blk%u_muscl=blk%u1
        end if
        do i=1,blk_size_nx
            u1=blk%u(1:5,i,1,1)
            u=blk%u_muscl(1:5,i,1,1)
            call eos_utow(u,w,temp,egv)
            p=w(5)
            slpx=blk%xslp(5,i,1,1)
            !pxl=blk%w_xl(5,i,1,1)
            !pxr=blk%w_xr(5,i,1,1)
            r=blk%x_center(i)
            ril=blk%x_interface(i-1)
            rir=blk%x_interface(i)
            !sx1=blk%surf1(i-1,1,1)
            !sx2=blk%surf1(i,1,1)
            vol=blk%vol(i,1,1)
            !fr=(p+pxr-pxl)*(sx2-sx1)/vol-2*(pxr-pxl)*r*(rir-ril)/vol
            fr=(p*(-ril**2 + rir**2) + ((-2*ril**3 + 2*rir**3 + 3*r*(ril**2 - rir**2))*slpx)/3d0)/vol
            u1(2)=u1(2)+fr*dt
            blk%u(1:5,i,1,1)=u1
        end do
    end if
end subroutine geometric_source

end module muscl

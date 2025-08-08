module recon_evolve
use mathlib
use datastructure
use eos
use phylib
use source_control

implicit none

contains

subroutine reconstruct_hydro(blk)
    type(blockdef), pointer :: blk
    !real(8), dimension(:,:,:,:), allocatable :: w,wl,wr
    !real(8), dimension(:,:,:), allocatable :: temp
    real(8), dimension(5) :: v1,v2,v3,slp_left,slp_right,slp_central,slp
    real(8) :: xl,xc,xr,xl_face,xr_face,yl,yc,yr,yl_face,yr_face,vv1,vv2,vv3,slp1,slp2,slp3,slpp
    integer :: i,j,k,kk
    do i=blk_xlb+1,blk_xub-1
        v1=blk%w(1:5,i-1,1,1)
        v2=blk%w(1:5,i,1,1)
        v3=blk%w(1:5,i+1,1,1)
        xl=blk%x_center(i-1)
        xc=blk%x_center(i)
        xr=blk%x_center(i+1)
        slp_left=(v2-v1)/(xc-xl)
        slp_right=(v3-v2)/(xr-xc)
        slp_central=(v3-v1)/(xr-xl)
        do k=1,5
            slp(k)=find_the_slope(slp_left(k),slp_right(k),slp_central(k))
        end do
        xl_face=blk%x_interface(i-1)
        xr_face=blk%x_interface(i)
        blk%xslp(1:5,i,1,1)=slp
        blk%w_xl(1:5,i,1,1)=blk%w(1:5,i,1,1)+(xl_face-xc)*slp
        blk%w_xr(1:5,i,1,1)=blk%w(1:5,i,1,1)+(xr_face-xc)*slp
#if     ieos==2
        blk%temp_xl(i,1,1)=solvetp(blk%w_xl(5,i,1,1),blk%w_xl(1,i,1,1))
        blk%temp_xr(i,1,1)=solvetp(blk%w_xr(5,i,1,1),blk%w_xr(1,i,1,1))
        blk%egv_xl(i,1,1)=egvrhot(blk%w_xl(1,i,1,1),blk%temp_xl(i,1,1))
        blk%egv_xr(i,1,1)=egvrhot(blk%w_xr(1,i,1,1),blk%temp_xr(i,1,1))
#endif
    end do
end subroutine reconstruct_hydro

end module recon_evolve

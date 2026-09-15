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
    real(8), dimension(5) :: v1, v2, v3, slp_left, slp_right, slp_central, slp
    real(8) :: xl, xc, xr, xl_face, xr_face, yl, yc, yr, yl_face, yr_face
    real(8) :: vv1, vv2, vv3, slp1, slp2, slp3, slpp,f
    integer :: i, j, k

    ! Precompute constants and arrays
    integer :: x_lb, x_ub, y_lb, y_ub
    !real(8), dimension(:), allocatable :: x_center, x_interface, y_center, y_interface
    real(8), dimension(:), pointer :: x_center, x_interface, y_center, y_interface

    x_lb = blk_xlb + 1
    x_ub = blk_xub - 1
    y_lb = blk_ylb + 1
    y_ub = blk_yub - 1

    ! Cache frequently accessed arrays
    x_center => blk%x_center
    x_interface => blk%x_interface
    y_center => blk%y_center
    y_interface => blk%y_interface
    ! Precompute differences
    do i = x_lb, x_ub
        xl = x_center(i-1)
        xc = x_center(i)
        xr = x_center(i+1)
        xl_face = x_interface(i-1)
        xr_face = x_interface(i)

        v1 = blk%w(1:5, i-1, 1, 1)
        v2 = blk%w(1:5, i, 1, 1)
        v3 = blk%w(1:5, i+1, 1, 1)

        ! Vectorized slope calculations
        slp_left = (v2 - v1) / (xc - xl)
        slp_right = (v3 - v2) / (xr - xc)
        slp_central = (v3 - v1) / (xr - xl)

        ! Vectorized slope limiting
        do k = 1, 5
            slp(k) = find_the_slope(slp_left(k), slp_right(k), slp_central(k))
        end do

        blk%xslp(1:5, i, 1, 1) = slp
        blk%w_xl(1:5, i, 1, 1) = v2 + (xl_face - xc) * slp
        blk%w_xr(1:5, i, 1, 1) = v2 + (xr_face - xc) * slp

#if ieos == 2
        ! Only compute EOS if needed
        blk%temp_xl(i, 1, 1) = solvetp(blk%w_xl(5, i, 1, 1), blk%w_xl(1, i, 1, 1))
        blk%temp_xr(i, 1, 1) = solvetp(blk%w_xr(5, i, 1, 1), blk%w_xr(1, i, 1, 1))
        blk%egv_xl(i, 1, 1) = egvrhot(blk%w_xl(1, i, 1, 1), blk%temp_xl(i, 1, 1))
        blk%egv_xr(i, 1, 1) = egvrhot(blk%w_xr(1, i, 1, 1), blk%temp_xr(i, 1, 1))
#endif
    end do
end subroutine reconstruct_hydro

function point_grav_potential(m,r)
    real(8) :: point_grav_potential,m,r
    point_grav_potential=-gr*m/r
end function point_grav_potential

end module recon_evolve

module source_control
use datastructure
use mathlib
use eos
use phylib
use problem
use gravity
use cooling
implicit none

contains

subroutine initialize_source()
#if    usersource
    user_source=>problem_source
#endif
end subroutine initialize_source

logical function source_applicability()
    source_applicability=.true.
    !if (iradiation==4.and.lsource) then
        !if (mod(time_sys%ntimestep,radhydro_boost)==0) then
        !    source_applicability=.false.
        !end if
        !if (time_sys%ntimestep==1) source_applicability=.false.
    !end if
end function source_applicability

subroutine source_apply(blk)
    !currently, only cooling, geometry and gravity are calculated by source_apply
    type(blockdef), pointer :: blk
    character(len=128) :: alert
    integer :: ijk(3),i,j,k
    real(8) :: t,dt,dedt,de,t_source,t_source_final,lum_atm,rthetaphi(3),ml_rate,rho
    real(8), dimension(5) :: u_temp,source,dsource,source_out,dsource_stiff,dsourcedt
    !need to calculate additional quantities when specific physics is use, e.g. divv and opacity
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            call user_source(blk,i,j)
        end do
    end do
end subroutine source_apply

end module source_control

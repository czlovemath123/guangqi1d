module boundary
use datastructure
use communication
use passivescalars
use problem
implicit none

contains

subroutine initialize_hydro()
    bound_hydro=>boundary_hydro
    init_hydro=>initial_hydro
end subroutine initialize_hydro

subroutine applyboundconds()
    integer :: bound1(4),bound2(4)
    bound1=hydro_bound_type
    bound2=rad_bound_type
    call time_dependent_bound_type()
    if (sum(abs(bound1-hydro_bound_type))/=0.or.sum(abs(bound2-rad_bound_type))/=0) then
        call relink_neighbors_tree()
        call assemble_communication_pattern()
    end if
    call blk_traversal(apply_explicit_bound)
end subroutine applyboundconds

subroutine apply_explicit_bound(blk)
    type(blockdef), pointer :: blk
    call applyboundconds_hydro_1d(blk)
end subroutine apply_explicit_bound

subroutine applyboundconds_hydro_1d(blk)
    type(blockdef), pointer :: blk
    integer :: i
    if (associated(blk,blk_head)) call coor1lower_hydro_internal(blk)
    if (associated(blk,blk_tail)) call coor1upper_hydro_internal(blk)
end subroutine applyboundconds_hydro_1d

!hydro internal communication step

subroutine coor1lower_hydro_internal(blk)
    !the lower boundary condition of the first coordinate
    type(blockdef), pointer :: blk,blk_temp
    integer :: i,j
    character(len=128) :: alert
    character(len=1) :: mirror
    logical :: reflect
    select case(hydro_bound_type(1))
    case(1)     !transmissive on the lower boundary of 1st coordinate
        do i=0,blk_xlb,-1
            blk%w(1:5,i,1,1)=blk%w(1:5,1-i,1,1)
            blk%u(1:5,i,1,1)=blk%u(1:5,1-i,1,1)
            blk%temp(i,1,1)=blk%temp(1-i,1,1)
            blk%egv(i,1,1)=blk%egv(1-i,1,1)
        end do
    case(2)     !reflective on the lower boundary of 1st coordinate
        do i=0,blk_xlb,-1
            blk%w(1:5,i,1,1)=blk%w(1:5,1-i,1,1)
            blk%w(2,i,1,1)=-blk%w(2,1-i,1,1)        !overwrite
            blk%u(1:5,i,1,1)=blk%u(1:5,1-i,1,1)
            blk%u(2,i,1,1)=-blk%u(2,1-i,1,1)        !overwrite
            blk%temp(i,1,1)=blk%temp(1-i,1,1)
            blk%egv(i,1,1)=blk%egv(1-i,1,1)
        end do
    case(3)     !periodic, not common
    case(9)     !user specified boundary condition
        do i=-1,0
            call bound_hydro(blk,i,1)
        end do
    end select
end subroutine coor1lower_hydro_internal

subroutine coor1upper_hydro_internal(blk)
    !the upper boundary condition of the first coordinate
    type(blockdef), pointer :: blk,blk_temp
    integer :: ijk(3),i,j
    character(len=128) :: alert
    character(len=1) :: mirror
    logical :: reflect
    select case(hydro_bound_type(2))
    case(1)     !transmissive on the upper boundary of 1st coordinate
        do i=blk_size_nx+1,blk_xub
            blk%w(1:5,i,1,1)=blk%w(1:5,2*blk_size_nx-i+1,1,1)
            blk%u(1:5,i,1,1)=blk%u(1:5,2*blk_size_nx-i+1,1,1)
            blk%temp(i,1,1)=blk%temp(2*blk_size_nx-i+1,1,1)
            blk%egv(i,1,1)=blk%egv(2*blk_size_nx-i+1,1,1)
        end do
    case(2)     !reflective on the upper boundary of 1st coordinate
        do i=blk_size_nx+1,blk_xub
            blk%w(1:5,i,1,1)=blk%w(1:5,2*blk_size_nx-i+1,1,1)
            blk%w(2,i,1,1)=-blk%w(2,2*blk_size_nx-i+1,1,1)          !overwrite
            blk%u(1:5,i,1,1)=blk%u(1:5,2*blk_size_nx-i+1,1,1)
            blk%u(2,i,1,1)=-blk%u(2,2*blk_size_nx-i+1,1,1)          !overwrite
            blk%temp(i,1,1)=blk%temp(2*blk_size_nx-i+1,1,1)
            blk%egv(i,1,1)=blk%egv(2*blk_size_nx-i+1,1,1)
        end do
    case(3)     !periodic, not common
    case(9)     !user specified boundary condition
        do i=blk_size_nx+1,blk_xub
            call bound_hydro(blk,i,1)
        end do
    end select
end subroutine coor1upper_hydro_internal

end module boundary

module passivescalars
use datastructure
use communication
use problem
implicit none

contains

subroutine initialize_passive_scalar()
    !independent first, dependent follows
    if (lam_con) then
        if (igeometry==1) then
            call add_passive_scalar(pn_am,angular_momentum)
        else if (igeometry==2.or.igeometry==3) then
            call add_passive_scalar(pn_am,angular_momentum)
        end if
    end if
    if (npassive>0) then
        bound_scalar=>boundary_scalar
        init_scalar=>initial_scalar
    end if
    cell_var_length=cell_var_length+npassive
end subroutine initialize_passive_scalar

subroutine add_passive_scalar(passivename,sub)
    type(passive_obj), pointer :: po
    character(len=16) :: passivename
    procedure(blk_operator), optional :: sub
    if (npassive==0) then
        allocate(passive_list)
        passive_list%passivename=passivename
        if (present(sub)) passive_list%sub=>sub
        passive_list%passive_type=1
        npassive=1
        call update_protected(lpassive,.true.)
    else
        po=>passive_list
        do while (associated(po%next))
            po=>po%next
        end do
        allocate(po%next)
        po=>po%next
        po%passivename=passivename
        if (present(sub)) po%sub=>sub
        po%passive_type=1
        npassive=npassive+1
    end if
end subroutine add_passive_scalar

subroutine apply_passive_init_condition(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    if (npassive>0) then
        do j=blk_ylb,blk_yub
            do i=blk_xlb,blk_xub
                call init_scalar(blk,i,j)
            end do
        end do
    end if
end subroutine apply_passive_init_condition

subroutine integrate_ipassive(passivename,blk,dt)
    type(blockdef), pointer :: blk
    character(len=16) :: passivename
    integer :: ip,i,j
    real(8) :: dt,sx1,sx2,sy1,sy2,vol,xpflux1,xpflux2,ypflux1,ypflux2,x
    real(8) :: mxflux1,mxflux2,myflux1,myflux2
    ip=ipassive(passivename)
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            sx1=blk%surf1(i-1,j,1)
            sx2=blk%surf1(i,j,1)
            sy1=blk%surf2(i,j-1,1)
            sy2=blk%surf2(i,j,1)
            vol=blk%vol(i,j,1)
            xpflux1=blk%xpflux(ip,i-1,j,1)
            xpflux2=blk%xpflux(ip,i,j,1)
            ypflux1=blk%ypflux(ip,i,j-1,1)
            ypflux2=blk%ypflux(ip,i,j,1)
            blk%passive_scalar(ip,i,j,1)=blk%passive_scalar0(ip,i,j,1)+(xpflux1*sx1-xpflux2*sx2+ypflux1*sy1-ypflux2*sy2)*dt/vol
        end do
    end do
end subroutine integrate_ipassive

subroutine integrate_all_dependent_passive(blk,dt)
    type(blockdef), pointer :: blk
    real(8) :: dt
    integer :: i,j
end subroutine integrate_all_dependent_passive

subroutine bound_edge_passive(blk,direction)
    type(blockdef), pointer :: blk
    character(len=16) :: direction
    integer :: i,j,ip_omega,ip_am
    if (direction==east) then
        if (passive_bound_type(2)==1.or.passive_bound_type(2)==2) then
            do i=blk_size_nx+1,blk_size_nx+2
                blk%passive_scalar(1:npassive,i,1:blk_size_ny,1) &
                =blk%passive_scalar(1:npassive,blk_size_nx*2+1-i,1:blk_size_ny,1)
            end do
        else if (passive_bound_type(2)==9) then
            do j=-1,blk_size_ny+2
                do i=blk_size_nx+1,blk_size_nx+2
                    call bound_scalar(blk,i,j)
                end do
            end do
        end if
    else if (direction==south) then
        if (passive_bound_type(3)==1.or.passive_bound_type(3)==2) then
            do j=-1,0
                blk%passive_scalar(1:npassive,1:blk_size_nx,j,1) &
                =blk%passive_scalar(1:npassive,1:blk_size_nx,1-j,1)
            end do
            if (wedge_polar) then
                ip_am=ipassive(pn_am)
                do j=-1,0
                    blk%passive_scalar(ip_am,1:blk_size_nx,j,1)=-blk%passive_scalar(ip_am,1:blk_size_nx,1-j,1)
                    blk%omega(1:blk_size_nx,j,1)=-blk%omega(1:blk_size_nx,1-j,1)
                end do
            end if
        else if (passive_bound_type(3)==9) then
            do j=-1,0
                do i=-1,blk_size_nx+2
                    call bound_scalar(blk,i,j)
                end do
            end do
        end if
    else if (direction==west) then
        if (passive_bound_type(1)==1.or.passive_bound_type(1)==2) then
            do i=-1,0
                blk%passive_scalar(1:npassive,i,1:blk_size_ny,1) &
                =blk%passive_scalar(1:npassive,1-i,1:blk_size_ny,1)
            end do
        else if (passive_bound_type(1)==9) then
            do j=-1,blk_size_ny+2
                do i=-1,0
                    call bound_scalar(blk,i,j)
                end do
            end do
        end if
    else if (direction==north) then
        if (passive_bound_type(4)==1.or.passive_bound_type(4)==2) then
            do j=blk_size_ny+1,blk_size_ny+2
                blk%passive_scalar(1:npassive,1:blk_size_nx,j,1) &
                =blk%passive_scalar(1:npassive,1:blk_size_nx,2*blk_size_ny+1-j,1)
            end do
        else if (passive_bound_type(4)==9) then
            do j=blk_size_ny+1,blk_size_ny+2
                do i=-1,blk_size_nx+2
                    call bound_scalar(blk,i,j)
                end do
            end do
        end if
    end if
end subroutine bound_edge_passive

subroutine update_passive_condition_all()
    !dependent passive scalar may change in other operators.
    type(passive_obj), pointer :: po
    procedure(blk_operator), pointer :: sub
    po=>passive_list
    do while (associated(po))
        if (po%passive_type==1) then
            sub=>po%sub
            call blk_traversal(sub)
        end if
        po=>po%next
    end do
end subroutine update_passive_condition_all

subroutine upwind_passive(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    integer :: ip_am
    real(8) :: r,theta,s
    do i=0,blk_size_nx
        if (blk%xflux(1,i,1,1)>0) then
            blk%xpflux(1:npassive,i,1,1)=blk%passive_scalar(1:npassive,i,1,1)/blk%w(1,i,1,1)*blk%xflux(1,i,1,1)
        else
            blk%xpflux(1:npassive,i,1,1)=blk%passive_scalar(1:npassive,i+1,1,1)/blk%w(1,i+1,1,1)*blk%xflux(1,i,1,1)
        end if
    end do
end subroutine upwind_passive

end module passivescalars

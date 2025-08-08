module communication_1d
#include <petsc/finclude/petscsys.h>
use petscsys
use datastructure
use tree_module
use phylib
use eos
implicit none

contains

subroutine communicate_flux_internal_1d(blk)
    type(blockdef), pointer :: blk
    call communicate_flux_left_internal(blk)
    call communicate_flux_right_internal(blk)
end subroutine communicate_flux_internal_1d

subroutine communicate_flux_external_1d(blk)
    type(blockdef), pointer :: blk
    call communicate_flux_left_external(blk)
    call communicate_flux_right_external(blk)
end subroutine communicate_flux_external_1d

subroutine communicate_flux_left_internal(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    integer :: i,nb_level,level,neighbour_rank
    if (associated(blk%nb(1)%blk)) then
        blk_nb=blk%nb(1)
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank==rank) then
            blk_temp=>blk_nb%blk
            if (blk%nb_level_1d(1)>0) then
                blk%xflux(1:5,0,1,1)=blk_temp%xflux(1:5,blk_size_nx,1,1)
                if (lrad_adv) blk%erad_xflux(0,1,1)=blk_temp%erad_xflux(blk_size_nx,1,1)
            end if
        end if
    end if
end subroutine communicate_flux_left_internal

subroutine communicate_flux_right_internal(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    integer :: i,nb_level,level,neighbour_rank
    if (associated(blk%nb(2)%blk)) then
        blk_nb=blk%nb(2)
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank==rank) then
            blk_temp=>blk_nb%blk
            if (blk%nb_level_1d(2)>0) then
                blk%xflux(1:5,blk_size_nx,1,1)=blk_temp%xflux(1:5,0,1,1)
                if (lrad_adv) blk%erad_xflux(blk_size_nx,1,1)=blk_temp%erad_xflux(0,1,1)
            end if
        end if
    end if
end subroutine communicate_flux_right_internal

subroutine communicate_flux_left_external(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    integer :: i,l,nb_level,level,neighbour_rank,tag,req,ierr,stat(MPI_STATUS_SIZE)
    real(8) :: flux(6)
    l=5
    if (lrad_adv) l=l+1
    if (associated(blk%nb(1)%blk)) then
        blk_nb=blk%nb(1)
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank/=rank) then
            blk_temp=>blk_nb%blk
            if (blk%nb_level_1d(1)<0) then
                flux(1:5)=blk%xflux(1:5,0,1,1)
                if (lrad_adv) flux(6)=blk%erad_xflux(0,1,1)
                tag=1;call mpi_isend(flux,l,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
            else if (blk%nb_level_1d(1)>0) then
                tag=1;call mpi_recv(flux,l,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
                blk%xflux(1:5,0,1,1)=flux(1:5)
                if (lrad_adv) blk%erad_xflux(0,1,1)=flux(6)
            end if
        end if
    end if
end subroutine communicate_flux_left_external

subroutine communicate_flux_right_external(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    integer :: i,l,nb_level,level,neighbour_rank,tag,req,ierr,stat(MPI_STATUS_SIZE)
    real(8) :: flux(6)
    l=5
    if (lrad_adv) l=l+1
    if (associated(blk%nb(2)%blk)) then
        blk_nb=blk%nb(2)
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank/=rank) then
            blk_temp=>blk_nb%blk
            if (blk%nb_level_1d(2)<0) then
                flux(1:5)=blk%xflux(1:5,blk_size_nx,1,1)
                if (lrad_adv) flux(6)=blk%erad_xflux(blk_size_nx,1,1)
                tag=1;call mpi_isend(flux,l,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
            else if (blk%nb_level_1d(2)>0) then
                tag=1;call mpi_recv(flux,l,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
                blk%xflux(1:5,blk_size_nx,1,1)=flux(1:5)
                if (lrad_adv) blk%erad_xflux(0,1,1)=flux(6)
            end if
        end if
    end if
end subroutine communicate_flux_right_external

subroutine boundary_restriction_1d(blk,ibound,w,u,egv,temp)
    !do restriction on one boundary size of blk, the result is stored in w,u,egv,temp
    type(blockdef), pointer, intent(in) :: blk
    real(8) :: w(5,2),u(5,2),temp(2),egv(2),vol1,vol2
    integer :: ibound,i,j
    u=0d0;w=0d0;temp=0d0;egv=0d0
    select case(ibound)
    case(1)         !xl boundary
        do i=1,n_hydro_guard
            if (igeometry==0) then
                u(1:5,i)=(blk%u(1:5,2*i-1,1,1)+blk%u(1:5,2*i,1,1))/2
            else if (igeometry==2) then
                u(1:5,i)=(blk%u(1:5,2*i-1,1,1)*blk%vol(2*i-1,1,1)+blk%u(1:5,2*i,1,1)*blk%vol(2*i,1,1))/(blk%vol(2*i-1,1,1)+blk%vol(2*i,1,1))
            end if
            call eos_utow(u(1:5,i),w(1:5,i),temp(i),egv(i))
        end do
    case(2)         !xu boundary
        do i=1,n_hydro_guard
            if (igeometry==0) then
                u(1:5,i)=(blk%u(1:5,blk_size_nx-2*n_hydro_guard+2*i-1,1,1)   &
                    +blk%u(1:5,blk_size_nx-2*n_hydro_guard+2*i,1,1))/2
            else if (igeometry==2) then
                vol1=blk%vol(blk_size_nx-2*n_hydro_guard+2*i-1,1,1)
                vol2=blk%vol(blk_size_nx-2*n_hydro_guard+2*i,1,1)
                u(1:5,i)=(blk%u(1:5,blk_size_nx-2*n_hydro_guard+2*i-1,1,1)*vol1   &
                    +blk%u(1:5,blk_size_nx-2*n_hydro_guard+2*i,1,1)*vol2)/(vol1+vol2)
            end if
            call eos_utow(u(1:5,i),w(1:5,i),temp(i),egv(i))
        end do
    case(3)         !yl boundary
    case(4)         !yu boundary
    end select
end subroutine boundary_restriction_1d

subroutine boundary_prolongation_1d(blk,ibound,w,u,egv,temp)
    !do prolongation on one boundary size of blk, the result is stored in w,u,egv,temp
    type(blockdef), pointer, intent(in) :: blk
    real(8) :: w(5,2),u(5,2),temp(2),egv(2)
    real(8) :: diff_l,diff_r,diff_c,diff
    integer :: ibound,i,j
    u=0d0;w=0d0;temp=0d0;egv=0d0
    select case(ibound)
    case(1)         !xl boundary
        if (n_hydro_guard/=2) then
            print *,'not implemented'
            stop
        end if
        do j=1,5
            diff_l=blk%w(j,1,1,1)-blk%w(j,0,1,1)
            diff_r=blk%w(j,2,1,1)-blk%w(j,1,1,1)
            diff_c=(blk%w(j,2,1,1)-blk%w(j,0,1,1))/2
            if (diff_c/=zero) then
                diff=abs(diff_c)/diff_c*min(abs(diff_l),abs(diff_r),abs(diff_c))
            else
                diff=0d0
            end if
            diff=0d0
            w(j,1)=blk%w(j,1,1,1)-half*diff
            w(j,2)=blk%w(j,1,1,1)+half*diff
        end do
        do i=1,2
            call eos_wtou(w(1:5,i),u(1:5,i),temp(i),egv(i))
        end do
    case(2)         !xu boundary
        if (n_hydro_guard/=2) then
            print *,'not implemented'
            stop
        end if
        do j=1,5
            diff_l=blk%w(j,blk_size_nx,1,1)-blk%w(j,blk_size_nx-1,1,1)
            diff_r=blk%w(j,blk_size_nx+1,1,1)-blk%w(j,blk_size_nx,1,1)
            diff_c=(blk%w(j,blk_size_nx+1,1,1)-blk%w(j,blk_size_nx-1,1,1))/2
            if (diff_c/=zero) then
                diff=abs(diff_c)/diff_c*min(abs(diff_l),abs(diff_r),abs(diff_c))
            else
                diff=0d0
            end if
            diff=0d0
            w(j,1)=blk%w(j,blk_size_nx,1,1)-half*diff
            w(j,2)=blk%w(j,blk_size_nx,1,1)+half*diff
        end do
        do i=1,2
            call eos_wtou(w(1:5,i),u(1:5,i),temp(i),egv(i))
        end do
    case(3)         !yl boundary
    case(4)         !yu boundary
    end select
end subroutine boundary_prolongation_1d

subroutine communicate_blocks_hydro_internal_1d(blk)
    type(blockdef), pointer :: blk
    call communicate_hydro_left_internal(blk)
    call communicate_hydro_right_internal(blk)
end subroutine communicate_blocks_hydro_internal_1d

subroutine communicate_hydro_left_internal(blk)
    !restrict or prolongate or copy the left neighbour's boundary cells and fill blk's guard cells
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    real(8) :: w(5,2),u(5,2),temp(2),egv(2)
    integer :: i,nb_level,level,neighbour_rank
    if (associated(blk%nb(1)%blk)) then
        blk_nb=blk%nb(1)
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank==rank) then
            blk_temp=>blk_nb%blk
            level=blk%level
            nb_level=blk_temp%level
            if (nb_level>level) then
                call boundary_restriction_1d(blk_temp,2,w,u,egv,temp)
            else if (nb_level<level) then
                call boundary_prolongation_1d(blk_temp,2,w,u,egv,temp)
            else
                !same level
                w=blk_temp%w(1:5,blk_size_nx-1:blk_size_nx,1,1)
                u=blk_temp%u(1:5,blk_size_nx-1:blk_size_nx,1,1)
                egv=blk_temp%egv(blk_size_nx-1:blk_size_nx,1,1)
                temp=blk_temp%temp(blk_size_nx-1:blk_size_nx,1,1)
            end if
            blk%w(1:5,-1:0,1,1)=w
            blk%u(1:5,-1:0,1,1)=u
            blk%temp(-1:0,1,1)=temp
            blk%egv(-1:0,1,1)=egv
        end if
        nullify(blk_temp)
    end if
end subroutine communicate_hydro_left_internal

subroutine communicate_hydro_right_internal(blk)
    !restrict or prolongate or copy the right neighbour's boundary cells and fill blk's guard cells
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    real(8) :: w(5,2),u(5,2),temp(2),egv(2)
    integer :: i,nb_level,level,neighbour_rank
    if (associated(blk%nb(2)%blk)) then
        blk_nb=blk%nb(2)
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank==rank) then
            blk_temp=>blk_nb%blk
            level=blk%level
            nb_level=blk_temp%level
            if (nb_level>level) then
                call boundary_restriction_1d(blk_temp,1,w,u,egv,temp)
            else if (nb_level<level) then
                call boundary_prolongation_1d(blk_temp,1,w,u,egv,temp)
            else
                !same level
                w=blk_temp%w(1:5,1:2,1,1)
                u=blk_temp%u(1:5,1:2,1,1)
                egv=blk_temp%egv(1:2,1,1)
                temp=blk_temp%temp(1:2,1,1)
            end if
            blk%w(1:5,blk_size_nx+1:blk_size_nx+2,1,1)=w
            blk%u(1:5,blk_size_nx+1:blk_size_nx+2,1,1)=u
            blk%temp(blk_size_nx+1:blk_size_nx+2,1,1)=temp
            blk%egv(blk_size_nx+1:blk_size_nx+2,1,1)=egv
        end if
        nullify(blk_temp)
    end if
end subroutine communicate_hydro_right_internal

subroutine communicate_blocks_hydro_external_1d(blk)
    type(blockdef), pointer :: blk
    call communicate_left_neighbour_external(blk)
    call communicate_right_neighbour_external(blk)
end subroutine communicate_blocks_hydro_external_1d

subroutine communicate_left_neighbour_external(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    real(8) :: w(5,2),u(5,2),temp(2),egv(2),msg_out(12,2),msg_in(12,2)
    integer :: i,ierr,nb_level,level,neighbour_rank,send_count,recv_count,tag,req,stat(MPI_STATUS_SIZE)
    if (associated(blk%nb(1)%blk)) then
        blk_nb=blk%nb(1)
        blk_temp=>blk_nb%blk
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank/=rank) then
            level=blk%level
            nb_level=blk_temp%level
            if (nb_level>level) then
                call boundary_prolongation_1d(blk,1,w,u,egv,temp)
                msg_out(1:5,1:2)=w
                msg_out(6:10,1:2)=u
                msg_out(11,1:2)=egv
                msg_out(12,1:2)=temp
            else if (nb_level<level) then
                call boundary_restriction_1d(blk,1,w,u,egv,temp)
                msg_out(1:5,1:2)=w
                msg_out(6:10,1:2)=u
                msg_out(11,1:2)=egv
                msg_out(12,1:2)=temp
            else
                !same level
                msg_out(1:5,1:2)=blk%w(1:5,1:2,1,1)
                msg_out(6:10,1:2)=blk%u(1:5,1:2,1,1)
                msg_out(11,1:2)=blk%egv(1:2,1,1)
                msg_out(12,1:2)=blk%temp(1:2,1,1)
            end if
            send_count=size(msg_out)
            tag=1;call mpi_isend(msg_out,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
            call mpi_recv(msg_in,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
            blk%w(1:5,-1:0,1,1)=msg_in(1:5,1:2)
            blk%u(1:5,-1:0,1,1)=msg_in(6:10,1:2)
            blk%egv(-1:0,1,1)=msg_in(11,1:2)
            blk%temp(-1:0,1,1)=msg_in(12,1:2)
        end if
    end if
end subroutine communicate_left_neighbour_external

subroutine communicate_right_neighbour_external(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    real(8) :: w(5,2),u(5,2),temp(2),egv(2),msg_out(12,2),msg_in(12,2)
    integer :: i,ierr,nb_level,level,neighbour_rank,send_count,recv_count,tag,req,stat(MPI_STATUS_SIZE)
    if (associated(blk%nb(2)%blk)) then
        blk_nb=blk%nb(2)
        blk_temp=>blk_nb%blk
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank/=rank) then
            level=blk%level
            nb_level=blk_temp%level
            if (nb_level>level) then
                call boundary_prolongation_1d(blk,2,w,u,egv,temp)
                msg_out(1:5,1:2)=w
                msg_out(6:10,1:2)=u
                msg_out(11,1:2)=egv
                msg_out(12,1:2)=temp
            else if (nb_level<level) then
                call boundary_restriction_1d(blk,2,w,u,egv,temp)
                msg_out(1:5,1:2)=w
                msg_out(6:10,1:2)=u
                msg_out(11,1:2)=egv
                msg_out(12,1:2)=temp
            else
                !same level
                msg_out(1:5,1:2)=blk%w(1:5,blk_size_nx-1:blk_size_nx,1,1)
                msg_out(6:10,1:2)=blk%u(1:5,blk_size_nx-1:blk_size_nx,1,1)
                msg_out(11,1:2)=blk%egv(blk_size_nx-1:blk_size_nx,1,1)
                msg_out(12,1:2)=blk%temp(blk_size_nx-1:blk_size_nx,1,1)
            end if
            send_count=size(msg_out)
            tag=1;call mpi_isend(msg_out,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
            call mpi_recv(msg_in,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
            blk%w(1:5,blk_size_nx+1:blk_size_nx+2,1,1)=msg_in(1:5,1:2)
            blk%u(1:5,blk_size_nx+1:blk_size_nx+2,1,1)=msg_in(6:10,1:2)
            blk%egv(blk_size_nx+1:blk_size_nx+2,1,1)=msg_in(11,1:2)
            blk%temp(blk_size_nx+1:blk_size_nx+2,1,1)=msg_in(12,1:2)
        end if
    end if
end subroutine communicate_right_neighbour_external

subroutine communicate_fld_internal_1d()
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    integer :: i,neighbour_rank
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        if (associated(blk%nb(1)%blk)) then
            blk_nb=blk%nb(1)
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank==rank) then
                blk_temp=>blk_nb%blk
                blk%Erad(0,1,1)=blk_temp%Erad(blk_size_nx,1,1)
                blk%Erad_int(0,1,1)=blk_temp%Erad_int(blk_size_nx,1,1)
                blk%sigma_rosseland(0,1,1)=blk_temp%sigma_rosseland(blk_size_nx,1,1)
            end if
        end if
        if (associated(blk%nb(2)%blk)) then
            blk_nb=blk%nb(2)
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank==rank) then
                blk_temp=>blk_nb%blk
                blk%Erad(blk_size_nx+1,1,1)=blk_temp%Erad(1,1,1)
                blk%Erad_int(blk_size_nx+1,1,1)=blk_temp%Erad_int(1,1,1)
                blk%sigma_rosseland(blk_size_nx+1,1,1)=blk_temp%sigma_rosseland(1,1,1)
            end if
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%next
        end if
    end do
    nullify(blk)
end subroutine communicate_fld_internal_1d

subroutine communicate_fld_external_1d()
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    real(8) :: msg_out(3),msg_in(3)
    integer :: i,neighbour_rank,req,ierr,stat(MPI_STATUS_SIZE)
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        if (associated(blk%nb(1)%blk)) then
            blk_nb=blk%nb(1)
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank/=rank) then
                msg_out(1)=blk%Erad(1,1,1)
                msg_out(2)=blk%Erad_int(1,1,1)
                msg_out(3)=blk%sigma_rosseland(1,1,1)
                call mpi_isend(msg_out,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,req,ierr)
                call mpi_recv(msg_in,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,stat,ierr)
                blk%Erad(0,1,1)=msg_in(1)
                blk%Erad_int(0,1,1)=msg_in(2)
                blk%sigma_rosseland(0,1,1)=msg_in(3)
            end if
        end if
        if (associated(blk%nb(2)%blk)) then
            blk_nb=blk%nb(2)
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank/=rank) then
                msg_out(1)=blk%Erad(blk_size_nx,1,1)
                msg_out(2)=blk%Erad_int(blk_size_nx,1,1)
                msg_out(3)=blk%sigma_rosseland(blk_size_nx,1,1)
                call mpi_isend(msg_out,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,req,ierr)
                call mpi_recv(msg_in,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,stat,ierr)
                blk%Erad(blk_size_nx+1,1,1)=msg_in(1)
                blk%Erad_int(blk_size_nx+1,1,1)=msg_in(2)
                blk%sigma_rosseland(blk_size_nx+1,1,1)=msg_in(3)
            end if
        end if
        blk=>blk%next
    end do
    nullify(blk)
end subroutine communicate_fld_external_1d

subroutine communicate_fld2_internal_1d()
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    integer :: i,neighbour_rank
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        if (associated(blk%nb(1)%blk)) then
            blk_nb=blk%nb(1)
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank==rank) then
                blk_temp=>blk_nb%blk
                blk%Erad(0,1,1)=blk_temp%Erad(blk_size_nx,1,1)
                blk%Erad_int(0,1,1)=blk_temp%Erad_int(blk_size_nx,1,1)
                blk%sigma_rosseland(0,1,1)=blk_temp%sigma_rosseland(blk_size_nx,1,1)
            end if
        end if
        if (associated(blk%nb(2)%blk)) then
            blk_nb=blk%nb(2)
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank==rank) then
                blk_temp=>blk_nb%blk
                blk%Erad(blk_size_nx+1,1,1)=blk_temp%Erad(1,1,1)
                blk%Erad_int(blk_size_nx+1,1,1)=blk_temp%Erad_int(1,1,1)
                blk%sigma_rosseland(blk_size_nx+1,1,1)=blk_temp%sigma_rosseland(1,1,1)
            end if
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%next
        end if
    end do
    nullify(blk)
end subroutine communicate_fld2_internal_1d

subroutine communicate_fld2_external_1d()
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    real(8) :: msg_out(3),msg_in(3)
    integer :: i,neighbour_rank,req,ierr,stat(MPI_STATUS_SIZE)
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        if (associated(blk%nb(1)%blk)) then
            blk_nb=blk%nb(1)
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank/=rank) then
                msg_out(1)=blk%Erad(1,1,1)
                msg_out(2)=blk%Erad_int(1,1,1)
                msg_out(3)=blk%sigma_rosseland(1,1,1)
                call mpi_isend(msg_out,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,req,ierr)
                call mpi_recv(msg_in,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,stat,ierr)
                blk%Erad(0,1,1)=msg_in(1)
                blk%Erad_int(0,1,1)=msg_in(2)
                blk%sigma_rosseland(0,1,1)=msg_in(3)
            end if
        end if
        if (associated(blk%nb(2)%blk)) then
            blk_nb=blk%nb(2)
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank/=rank) then
                msg_out(1)=blk%Erad(blk_size_nx,1,1)
                msg_out(2)=blk%Erad_int(blk_size_nx,1,1)
                msg_out(3)=blk%sigma_rosseland(blk_size_nx,1,1)
                call mpi_isend(msg_out,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,req,ierr)
                call mpi_recv(msg_in,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,stat,ierr)
                blk%Erad(blk_size_nx+1,1,1)=msg_in(1)
                blk%Erad_int(blk_size_nx+1,1,1)=msg_in(2)
                blk%sigma_rosseland(blk_size_nx+1,1,1)=msg_in(3)
            end if
        end if
        blk=>blk%next
    end do
    nullify(blk)
end subroutine communicate_fld2_external_1d

subroutine prolongation_cell_1d(blk,i,prolonged)
    type(blockdef), pointer :: blk
    integer :: i,ii,l,iblkx
    real(8), allocatable :: slp(:,:),v(:),v_temp(:),prolonged(:,:)
    real(8) :: dx(2),dx_cell,r1,r2,r
    l=cell_var_length
    if (iradiation/=0) l=l+1
    allocate(slp(l,1),v(l),v_temp(l))
    call calculate_slp_cell(blk,i,1,slp)
    slp=0d0
    if (igeometry==0) then
        dx=blk%dxyz(1)/4d0;dx(1)=-dx(2)
    else if (igeometry==2) then
        if (llnx) then
            r1=blk%x_interface(i-1)
            iblkx=(blk%key(1)-1)*2+1
            r2=geometric_si(xgeo_a,xgeo_b,xgeo_dp,blk%level+1,iblkx,2*i-1)
            r=spherical_r_center(r1,r2)
            dx(1)=-(blk%x_center(i)-r)
            r1=r2
            r2=blk%x_interface(i)
            r=spherical_r_center(r1,r2)
            dx(2)=r-blk%x_center(i)
        else
            dx_cell=blk%dxyz(1)
            r1=blk%x_interface(i-1)
            r2=r1+dx_cell/2d0
            r=spherical_r_center(r1,r2)
            dx(1)=-(blk%x_center(i)-r)
            r1=r2
            r2=blk%x_interface(i)
            r=spherical_r_center(r1,r2)
            dx(2)=r-blk%x_center(i)
        end if
    end if
    v(1:5)=blk%w(1:5,i,1,1)
    if (iradiation/=0) v(l)=blk%Erad(i,1,1)/blk%w(1,i,1,1)
    if (lpassive) v(6:5+npassive)=blk%passive_scalar(1:npassive,i,1,1)/blk%w(1,i,1,1)
    do ii=1,2
        v_temp=v+slp(1:l,1)*dx(ii)
        prolonged(1:5,ii)=v_temp(1:5)
        if (iradiation/=0) prolonged(l,ii)=v_temp(l)*v_temp(1)
        if (lpassive) prolonged(6:5+npassive,ii)=v_temp(6:5+npassive)*v_temp(1)
    end do
    deallocate(slp,v,v_temp)
end subroutine prolongation_cell_1d

end module communication_1d

module communication
use datastructure
use phylib
use eos
use communication_1d
use problem
implicit none

interface prolongation_cell
    module procedure prolongation_cell_1d
end interface prolongation_cell

contains

subroutine communicate_flux()
    call blk_traversal(communicate_flux_internal_1d)
    call blk_traversal(communicate_flux_external_1d)
end subroutine communicate_flux

subroutine communicate_hydro()
    call blk_traversal(communicate_blocks_hydro_internal_1d)
    call blk_traversal(communicate_blocks_hydro_external_1d)
end subroutine communicate_hydro

subroutine communicate_fld()
    call communicate_fld_internal_1d()
    call communicate_fld_external_1d()
end subroutine communicate_fld

subroutine communicate_fld2()
    call communicate_fld2_internal_1d()
    call communicate_fld2_external_1d()
end subroutine communicate_fld2

subroutine assemble_communication_pattern()
    type(blockdef), pointer :: blk
    integer :: i,j,k
end subroutine assemble_communication_pattern

subroutine block_restriction(blk)
    type(blockdef), pointer :: blk,blk1,blk2,blk3,blk4,blk_temp
    real(8) :: u(5)
    real(8), allocatable :: a(:)
    integer :: i,j,k,l,ilocal,jlocal,key(3),ii,jj,iloc,jloc,idx,idy
    if (allocated(blk%w)) call release_block(blk)
    call allocate_block_heavy_data(blk)
    l=cell_var_length
    !print *,'restriction',blk%key,time_sys%ntimestep
    blk%on_processor=.true.
    blk1=>blk%children(1,1)%blk
    blk2=>blk%children(2,1)%blk
    do i=0,blk_size_nx/2
        u=(blk1%u(1:5,2*i-1,1,1)+blk1%u(1:5,2*i,1,1))/2
        call assign_u_to_w_cell(blk,u,i,1)
        if (iradiation/=0) blk%Erad(i,1,1)=(blk1%Erad(2*i-1,1,1)+blk1%Erad(2*i,1,1))/2
    end do
    do i=blk_size_nx/2+1,blk_size_nx+1
        ilocal=i-blk_size_nx/2
        u=(blk2%u(1:5,2*ilocal-1,1,1)+blk2%u(1:5,2*ilocal,1,1))/2
        call assign_u_to_w_cell(blk,u,i,1)
        if (iradiation/=0) blk%Erad(i,1,1)=(blk2%Erad(2*ilocal-1,1,1)+blk2%Erad(2*ilocal,1,1))/2
    end do
    call destroy_block(blk%children(1,1)%blk)
    call destroy_block(blk%children(2,1)%blk)
    !print *,time_sys%ntimestep,'restriction',blk%key
end subroutine block_restriction

subroutine block_prolongation(blk)
    type(blockdef), pointer :: blk,blk_temp
    real(8) :: diff_l(5),diff_r(5),diff_c(5),diff(5),omega_prolonged(2,2)
    real(8), allocatable :: prolonged(:,:,:),prolonged_1d(:,:)
    integer :: i,j,k,ilocal,jlocal,key(3),ii,jj,l,idx,idy
    logical :: guard
    guard=.true.
    l=cell_var_length
    !print *,'prolongation',blk%key,time_sys%ntimestep
    blk%refine_countdown=0
    if (iradiation/=0) l=l+1
    allocate(prolonged_1d(l,2))
    do idx=1,2
        blk_temp=>blk%children(idx,1)%blk
        blk_temp%on_processor=.true.
        blk_temp%refine_countdown=10
        call allocate_block_heavy_data(blk_temp)
        ii=1-mod(blk_temp%key(1),2)
        do i=0,blk_size_nx/2+1
            ilocal=ii*blk_size_nx/2+i
            call prolongation_cell(blk,ilocal,prolonged_1d)
            blk_temp%w(1:5,2*i-1:2*i,1,1)=prolonged_1d(1:5,1:2)
            if (iradiation/=0) blk_temp%Erad(2*i-1:2*i,1,1)=prolonged_1d(l,1:2)
            if (lpassive) blk_temp%passive_scalar(1:npassive,2*i-1:2*i,1,1)=prolonged_1d(6:5+npassive,1:2)
        end do
        call convert_w_to_u_block(blk_temp,guard)
        !print *,blk%key,blk_temp%key
        !print *,blk_temp%temp
    end do
    deallocate(prolonged_1d)
    call release_block(blk)
end subroutine block_prolongation

subroutine record_processor_amr_grow_order(blk)
    !record the keys of the blocks that are refined
    type(blockdef), pointer :: blk
    integer, dimension(:,:), allocatable :: temp_keys
    processor%amr_ngrow=processor%amr_ngrow+1
    if (processor%amr_ngrow==1) then
        !the first grown key
        allocate(processor%grow_keys(3,1))
        processor%grow_keys(1:3,1)=blk%key
    else
        !add another key to the existing list
        allocate(temp_keys(3,processor%amr_ngrow))
        temp_keys(1:3,1:processor%amr_ngrow-1)=processor%grow_keys
        temp_keys(1:3,processor%amr_ngrow)=blk%key
        deallocate(processor%grow_keys)
        allocate(processor%grow_keys(3,processor%amr_ngrow))
        processor%grow_keys=temp_keys
        deallocate(temp_keys)
    end if
end subroutine record_processor_amr_grow_order

subroutine record_processor_amr_derefine_order(blk)
    !record the keys of the blocks that are derefined
    type(blockdef), pointer :: blk
    integer, dimension(:,:), allocatable :: temp_keys
    processor%amr_nderefine=processor%amr_nderefine+1
    if (processor%amr_nderefine==1) then
        allocate(processor%derefine_keys(3,1))
        processor%derefine_keys(1:3,1)=blk%key
    else
        allocate(temp_keys(3,processor%amr_nderefine))
        temp_keys(1:3,1:processor%amr_nderefine-1)=processor%derefine_keys
        temp_keys(1:3,processor%amr_nderefine)=blk%key
        deallocate(processor%derefine_keys)
        allocate(processor%derefine_keys(3,processor%amr_nderefine))
        processor%derefine_keys=temp_keys
        deallocate(temp_keys)
    end if
end subroutine record_processor_amr_derefine_order

subroutine send_block(blk,receiver_rank,istart,msg,reqs,ireq)
    !save the block data to msg and release the block to send
    type(blockdef), pointer :: blk
    integer :: receiver_rank,tag,ierr,req,istart,iend,imsg,i,j,k,ireq
    integer, allocatable :: reqs(:)
    real(8), dimension(:), allocatable :: msg
    real(8) :: w(5)
    imsg=istart
    msg(imsg:imsg+blk_hydro_size-1)=reshape(blk%w,(/blk_hydro_size/));imsg=imsg+blk_hydro_size
    if (iradiation/=0) then
        msg(imsg:imsg+blk_cell_size-1)=reshape(blk%Erad,(/blk_cell_size/))
        imsg=imsg+blk_cell_size
    end if
    if (lam_con) then
        msg(imsg:imsg+blk_cell_size-1)=reshape(blk%omega,(/blk_cell_size/))
        imsg=imsg+blk_cell_size
    end if
    iend=imsg-1   !istart+blksize-1
    tag=blk%blk_id;call mpi_isend(msg(istart:iend),blksize,MPI_REAL8,receiver_rank,tag,MPI_COMM_WORLD,reqs(ireq),ierr)
    ireq=ireq+1
    istart=iend+1
    call release_block(blk)
end subroutine send_block

subroutine receive_block(blk,former_rank)
    type(blockdef), pointer :: blk
    integer :: former_rank,tag,ierr,stat(MPI_STATUS_SIZE),imsg,ip_am,i,j,k
    real(8), dimension(:), allocatable :: msg_in
    real(8) :: omega,l_omega,rho,w(5)
    if (.not.allocated(blk%w)) call allocate_block_heavy_data(blk)
    allocate(msg_in(blksize))
    tag=blk%blk_id;call mpi_recv(msg_in,blksize,MPI_REAL8,former_rank,tag,MPI_COMM_WORLD,stat,ierr)
    imsg=1
    blk%w=reshape(msg_in(imsg:imsg+blk_hydro_size-1),shape(blk%w));imsg=imsg+blk_hydro_size
    do j=blk_ylb,blk_yub
        do i=blk_xlb,blk_xub
            w=blk%w(1:5,i,j,1)
            if (w(1)/=0) call assign_w_to_u_cell(blk,w,i,j)
        end do
    end do
    if (iradiation/=0) then
        blk%Erad=reshape(msg_in(imsg:imsg+blk_cell_size-1),shape(blk%Erad));imsg=imsg+blk_cell_size
    end if
    if (lam_con) then
        blk%omega=reshape(msg_in(imsg:imsg+blk_cell_size-1),shape(blk%omega));imsg=imsg+blk_cell_size
        call angular_momentum(blk)
    end if
    deallocate(msg_in)
end subroutine receive_block

subroutine sync_load()
    !assume each processor has correct llist_head and llist_tail position
    !sync the np_nblk across the processors
    type(blockdef), pointer :: blk
    integer :: nblk,ierr
    call renumber_domain_blocks()
    nblk=1
    blk=>llist_head
    do while (.not.associated(blk,llist_tail))
        nblk=nblk+1
        blk=>blk%next
    end do
    np_nblk(rank+1)=nblk
    call mpi_gather(np_nblk(rank+1),1,MPI_INTEGER,np_nblk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(np_nblk,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine sync_load

subroutine sync_blk_attribute_int(get_attr_int,assign_attr_int)
    type(blockdef), pointer :: blk
    procedure(blk_operator_int) :: assign_attr_int
    procedure(blk_get_attr_int) :: get_attr_int
    integer, dimension(:), allocatable :: local_int,global_int
    integer :: i
    allocate(local_int(np_nblk(rank+1)),global_int(nblk_total))
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        local_int(i)=get_attr_int(blk)
        blk=>blk%next
    end do
    call sync_bcast(local_int,global_int)
    blk=>blk_head
    do i=1,nblk_total
        call assign_attr_int(blk,global_int(i))
        blk=>blk%next
    end do
    deallocate(local_int,global_int)
end subroutine sync_blk_attribute_int

subroutine sync_blk_attribute_logical(get_attr_logical,assign_attr_logical)
    type(blockdef), pointer :: blk
    procedure(blk_operator_logical) :: assign_attr_logical
    procedure(blk_get_attr_logical) :: get_attr_logical
    logical, dimension(:), allocatable :: local_logical,global_logical
    integer :: i
    allocate(local_logical(np_nblk(rank+1)),global_logical(nblk_total))
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        local_logical(i)=get_attr_logical(blk)
        blk=>blk%next
    end do
    call sync_bcast(local_logical,global_logical)
    blk=>blk_head
    do i=1,nblk_total
        call assign_attr_logical(blk,global_logical(i))
        blk=>blk%next
    end do
    deallocate(local_logical,global_logical)
end subroutine sync_blk_attribute_logical

subroutine sync_directional_refine()
    type(blockdef), pointer :: blk
    integer :: i,reqs,ierr,stat(MPI_STATUS_SIZE)
    logical, dimension(:,:), allocatable :: grefine,lrefine
    allocate(lrefine(8,np_nblk(rank+1)),grefine(8,nblk_total))
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        lrefine(:,i)=blk%refine_direction
        blk=>blk%next
    end do
    call sync_bcast(lrefine,grefine)
    blk=>blk_head
    do i=1,nblk_total
        blk%refine_direction=grefine(:,i)
        blk=>blk%next
    end do
    deallocate(lrefine,grefine)
end subroutine sync_directional_refine

subroutine sync_the_tree_refine()
    !sync by doing the same action, not by copying the final state of the tree
    type(blockdef), pointer :: blk
    integer :: ierr,i,j,n_oper,key(3)
    integer, dimension(:,:), allocatable :: keys
    call mpi_gather(processor%amr_ngrow,1,MPI_INTEGER,processor%mpi_processor_oper,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(processor%mpi_processor_oper,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    np_nblk=np_nblk+processor%mpi_processor_oper*(2**nd-1)
    n_oper=sum(processor%mpi_processor_oper)
    if (n_oper>0) then
        change_mesh=.true.
        allocate(processor%mpi_world_keys(3,n_oper))
        call sync_bcast(processor%grow_keys,processor%mpi_world_keys,processor%mpi_processor_oper)
        do i=1,n_oper
            key=processor%mpi_world_keys(1:3,i)
            call block_key_to_pointer(key,blk)
            call grow_tree_node(blk)
        end do
        call renumber_domain_blocks()
        call reset_amr_processor_variables()
    end if
    call sync_blk_attribute_int(get_refine_countdown,assign_refine_countdown)
end subroutine sync_the_tree_refine

integer function get_refine_countdown(blk)
    type(blockdef), pointer :: blk
    get_refine_countdown=blk%refine_countdown
end function get_refine_countdown

subroutine assign_refine_countdown(blk,i)
    type(blockdef), pointer :: blk
    integer :: i
    blk%refine_countdown=max(i-1,0)
end subroutine assign_refine_countdown

subroutine sync_the_tree_derefine()
    !sync by doing the same action, not by copying the final state of the tree
    type(blockdef), pointer :: blk
    integer :: ierr,i,j,n_oper,p_start,send_count,receive_count,key(3)
    integer, dimension(:,:), allocatable :: keys
    integer, dimension(:), allocatable :: lderefine_countdown,gderefine_countdown
    call mpi_gather(processor%amr_nderefine,1,MPI_INTEGER,processor%mpi_processor_oper,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(processor%mpi_processor_oper,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    np_nblk=np_nblk-processor%mpi_processor_oper*(2**nd-1)
    n_oper=sum(processor%mpi_processor_oper)
    if (n_oper>0) then
        change_mesh=.true.
        allocate(processor%mpi_world_keys(3,n_oper))
        call sync_bcast(processor%derefine_keys,processor%mpi_world_keys,processor%mpi_processor_oper)
        do i=1,n_oper
            key=processor%mpi_world_keys(1:3,i)
            call block_key_to_pointer(key,blk)
            call trim_tree_node(blk)
        end do
        call renumber_domain_blocks()
        call reset_amr_processor_variables()
    end if
    call sync_blk_attribute_int(get_derefine_countdown,assign_derefine_countdown)
end subroutine sync_the_tree_derefine

integer function get_derefine_countdown(blk)
    type(blockdef), pointer :: blk
    get_derefine_countdown=blk%derefine_countdown
end function get_derefine_countdown

subroutine assign_derefine_countdown(blk,i)
    type(blockdef), pointer :: blk
    integer :: i
    blk%derefine_countdown=max(i-1,0)
end subroutine assign_derefine_countdown

subroutine send_recv_cross_processor_derefine()
    type(blockdef), pointer :: pa,blk,blk_temp
    integer :: i,key(3),level,levels1d(2),levels2d(4),nsend,nrecv,istart,ierr,idx,idy
    integer, allocatable :: nsend_all(:),nrecv_all(:),srank_key(:,:),rrank_key(:,:),rank_temp(:,:)
    logical :: lon_processor(2,2),lderefine(2,2)
    nsend=0;nrecv=0
    allocate(rank_temp(4,50));rank_temp=0
    if (rank/=0) then
        blk=>llist_head
        if (mod(blk%key(1),2)/=1) then
            if (blk%derefine.and.blk%pre%derefine) then
                level=blk%level;levels1d=0
                pa=>blk%pa
                call deepest_child_level(pa%children(1,1)%blk,levels1d(1))
                call deepest_child_level(pa%children(2,1)%blk,levels1d(2))
                if (sum(abs(levels1d-level))==0) then
                    nsend=nsend+1
                    rank_temp(1:4,nsend)=(/rank-1,llist_head%key/)
                    !srank_key(1:4,nsend)=(/rank-1,llist_head%key/)
                end if
            end if
        end if
    end if
    if (nsend>0) then
        allocate(srank_key(4,nsend))
        srank_key(1:4,1:nsend)=rank_temp(1:4,1:nsend)
    end if
    if (rank/=np-1) then
        blk=>llist_tail
        if (mod(blk%key(1),2)/=0) then
            if (blk%derefine.and.blk%next%derefine) then
                level=blk%level;levels1d=0
                pa=>blk%pa
                call deepest_child_level(pa%children(1,1)%blk,levels1d(1))
                call deepest_child_level(pa%children(2,1)%blk,levels1d(2))
                if (sum(abs(levels1d-level))==0) then
                    nrecv=nrecv+1
                    rank_temp(1:4,nrecv)=(/rank+1,llist_tail%next%key/)
                    !rrank_key(1:4,nrecv)=(/rank+1,llist_tail%next%key/)
                end if
            end if
        end if
    end if
    if (nrecv>0) then
        allocate(rrank_key(4,nrecv))
        rrank_key(1:4,1:nrecv)=rank_temp(1:4,1:nrecv)
    end if
    call send_recv_with_sequence(srank_key,rrank_key,nsend,nrecv)
    allocate(nsend_all(np),nrecv_all(np))
    call mpi_gather(nsend,1,MPI_INTEGER,nsend_all,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_gather(nrecv,1,MPI_INTEGER,nrecv_all,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nsend_all,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nrecv_all,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    np_nblk=np_nblk+nrecv_all-nsend_all
    call relocate_llist_head_tail(np_nblk)
    call identify_on_processor()
    deallocate(nsend_all,nrecv_all,rank_temp)
end subroutine send_recv_cross_processor_derefine

subroutine send_recv_with_sequence(srank_key,rrank_key,nsend,nrecv)
    !send/recv with a given sequence of send/recv keys
    type(blockdef), pointer :: blk
    integer :: nsend,nrecv,i,rank1,rank2,istart,key(3),ierr,ireq
    integer, allocatable :: reqs(:)
    integer, dimension(:,:), allocatable :: srank_key,rrank_key
    real(8), dimension(:), allocatable :: msg
    if (nsend>0) then
        allocate(msg(nsend*(blksize)),reqs(nsend))
        istart=1;ireq=1
        do i=1,nsend
            rank1=srank_key(1,i)
            key=srank_key(2:4,i)
            call block_key_to_pointer(key,blk)
            call send_block(blk,rank1,istart,msg,reqs,ireq)
        end do
    end if
    if (nrecv>0) then
        do i=1,nrecv
            rank2=rrank_key(1,i)
            key=rrank_key(2:4,i)
            call block_key_to_pointer(key,blk)
            call receive_block(blk,rank2)
        end do
    end if
    if (allocated(msg)) then
        call mpi_waitall(nsend,reqs,MPI_STATUSES_IGNORE,ierr)
        deallocate(msg)
    end if
end subroutine send_recv_with_sequence

subroutine start_mpi()
    integer :: ierr
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    if (ierr .ne. 0) then
        print*,'Unable to initialize PETSc'
        stop
    end if
end subroutine start_mpi

subroutine finalize_mpi()
    integer :: ierr
    call PetscFinalize(ierr)
end subroutine finalize_mpi

end module communication

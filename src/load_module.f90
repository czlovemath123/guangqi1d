module load_module
use datastructure
use tree_module
use communication
implicit none

logical, protected :: load_rebalanced=.false.

contains

subroutine load_balancer_static()
    !only used when there is no block data move
    call renumber_domain_blocks()
    call calculate_load(np_nblk)
    call relocate_llist_head_tail(np_nblk)
    !call display_llist()
    !stop
    call identify_on_processor()
end subroutine load_balancer_static

subroutine load_balancer_dynamic()
    type(blockdef), pointer :: blk
    integer :: iblk,i,j,k,blk_id,former_rank,new_rank,ierr,istart,nsend,nrecv,key(3)
    integer, dimension(:,:), allocatable :: srank_key,rrank_key
    real(8) :: w(5),temp
    call renumber_domain_blocks()
    call assemble_send_recv_sequence(srank_key,rrank_key,nsend,nrecv)
    call send_recv_with_sequence(srank_key,rrank_key,nsend,nrecv)
    call identify_on_processor()
    call reset_refine_derefine_countdown(blk_root)
    call relink_neighbors_tree()
end subroutine load_balancer_dynamic

subroutine assemble_send_recv_sequence(srank_key,rrank_key,nsend,nrecv)
    type(blockdef), pointer :: blk
    integer :: nsend,nrecv,i,new_rank,old_rank
    integer, allocatable :: pre_load(:),new_load(:),srank_key(:,:),rrank_key(:,:),rank_temp(:,:),op_load(:)
    real(8) :: maxload,minload
    nsend=0;nrecv=0
    maxload=maxval(np_nblk)
    minload=minval(np_nblk)
    if (maxload/minload>1.1.and.real(nblk_total)/real(np)>3*2.0**nd) then
    !if (.true.) then
        if (rank==0) print *,'load optimization',time_sys%ntimestep
        allocate(pre_load(np),new_load(np),op_load(np))
        allocate(rank_temp(4,1000));rank_temp=0
        pre_load=np_nblk
        call calculate_load(new_load)
        call load_optimization(new_load,op_load)
        if ((sum(np_nblk)-sum(op_load))/=0) then
            print *,'nblk not conserved after load balancing step'
            print *,sum(np_nblk)-sum(op_load),sum(np_nblk)-sum(new_load)
            stop
        end if
        if (sum(abs(np_nblk-op_load))/=0) then
            !assemble to send block info
            change_mesh=.true.
            blk=>llist_head
            do i=1,pre_load(rank+1)
                call block_id_to_processor_rank(op_load,blk%blk_id,new_rank)
                if (rank/=new_rank) then
                    nsend=nsend+1
                    rank_temp(1:4,nsend)=(/new_rank,blk%key/)
                end if
                blk=>blk%next
            end do
            if (nsend>0) then
                allocate(srank_key(4,nsend))
                srank_key(1:4,1:nsend)=rank_temp(1:4,1:nsend)
            end if
            rank_temp=0
            call relocate_llist_head_tail(op_load)
            blk=>llist_head
            do i=1,op_load(rank+1)
                call block_id_to_processor_rank(pre_load,blk%blk_id,old_rank)
                if (rank/=old_rank) then
                    nrecv=nrecv+1
                    rank_temp(1:4,nrecv)=(/old_rank,blk%key/)
                end if
                blk=>blk%next
            end do
            if (nrecv>0) then
                allocate(rrank_key(4,nrecv))
                rrank_key(1:4,1:nrecv)=rank_temp(1:4,1:nrecv)
            end if
            np_nblk=op_load
        end if
        deallocate(pre_load,new_load,rank_temp,op_load)
    end if
end subroutine assemble_send_recv_sequence

subroutine load_optimization(load,optimized_load)
    !let the children blocks on the same core.
    type(blockdef), pointer :: blk,blk_temp
    integer, allocatable :: load(:),optimized_load(:)
    integer :: i,idx,idy,nblk,ierr
    optimized_load=load
    !call relocate_llist_head_tail(load)
    !blk=>llist_head
    !blk_temp=>blk%pa
    !if (.not.associated(blk_temp%children(1,1)%blk,blk)) then
    !    llist_head=>blk_temp%children(1,1)%blk
    !end if
    !blk=>llist_tail
    !blk_temp=>blk%pa
    !if (.not.associated(blk_temp%children(2,2)%blk,blk)) then
    !    llist_tail=>blk_temp%children(1,1)%blk%pre
    !end if
    !nblk=llist_tail%blk_id-llist_head%blk_id+1
    !call mpi_gather(nblk,1,MPI_INTEGER,optimized_load,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    !call mpi_bcast(optimized_load,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    !call relocate_llist_head_tail(np_nblk)
end subroutine load_optimization

recursive subroutine reset_refine_derefine_countdown(blk)
    !reset all non leaf blocks countdown to be 0
    type(blockdef), pointer :: blk
    integer :: idx,idy
    idx=curve_seq_idx(blk%curve,1)
    idy=curve_seq_idy(blk%curve,1)
    if (blk%on_processor.eqv..false..and.associated(blk%children(idx,idy)%blk)) then
        blk%refine_countdown=0
        blk%derefine_countdown=0
        do idy=1,2
            do idx=1,2
                if (associated(blk%children(idx,idy)%blk)) call reset_refine_derefine_countdown(blk%children(idx,idy)%blk)
            end do
        end do
    end if
end subroutine reset_refine_derefine_countdown

end module load_module

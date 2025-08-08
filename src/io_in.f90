module io_in
use tree_module
use communication
use amr_module
use io_out
use mathlib
use phylib
use boundary
use eos
use problem
use HDF5
use mpi
implicit none

contains

subroutine hdf_read()
    character(len=32) :: filename_hdf,blockname,attname
    integer(HID_T) :: file_id,att_id
    integer(HSIZE_T) :: att_dim(1)
    integer :: error,ierr,lastframe
    character(len=16) :: fmt1,s1
    if (rank==0) then
        fmt1='(I5.5)'
        write(s1,fmt1) iframe
        filename_hdf=trim(filename_head)//trim(s1)//'.h5'
        call chdir(path_out)
        att_dim=1
        call h5open_f(error)
        call h5fopen_f(filename_hdf,H5F_ACC_RDWR_F,file_id,error)
        call h5aopen_f(file_id,'nblocks_tree',att_id,error)
        call h5aread_f(att_id,h5t_native_integer,nblocks_tree,att_dim,error)
        call h5aclose_f(att_id,error)
        call h5aopen_f(file_id,'Time',att_id,error)
        call h5aread_f(att_id,h5t_native_double,time_sys%t,att_dim,error)
        call h5aclose_f(att_id,error)
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(nblocks_tree,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(time_sys%t,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call read_tree(file_id)
    if (rank==0) then
        call h5fclose_f(file_id,error)
        call h5close_f(error)
        call chdir(path_root)
    end if
end subroutine hdf_read

subroutine read_tree(file_id)
    !read in the keys and rebuild the tree
    type(blockdef), pointer :: blk
    integer(HID_T) :: file_id,dset_id
    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer :: i,error,ierr,ii,jj,kk
    integer, dimension(:,:), allocatable :: keys
    allocate(dims(2),keys(nblocks_tree,6))
    if (rank==0) then
        dims=(/nblocks_tree,6/)
        call h5dopen_f(file_id,'block_tree',dset_id,error)
        call h5dread_f(dset_id,h5t_native_integer,keys,dims,error)
        call h5dclose_f(dset_id,error)
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(keys,nblocks_tree*6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call rebuild_tree_list(keys)
    if (rank==0) then
        call reallocate_blocks(blk_root)
        call read_heavy_data(file_id)
    end if
    call distribute_heavy_data()
    call restart_complete_all_data()
    call applyboundconds()
    call restart_static_level()
#if     usersource
    call initialize_iso_temp()
#endif
    !call clean_tree_heavy_data(blk_root)
    if (refine_type=='static'.and.rank==0) then
        print *,'total block number=',nblk_total,'block per core=',real(nblk_total)/real(np)
    end if
    deallocate(dims,keys)
end subroutine read_tree

subroutine initialize_iso_temp()
    type(blockdef), pointer :: blk
    integer :: iblk
    blk=>llist_head
    do iblk=1,np_nblk(rank+1)
        blk%iso_temp=blk%temp
        blk=>blk%next
    end do
end subroutine initialize_iso_temp

subroutine rebuild_tree_list(keys)
    !preorder traversal
    !logical level have level=-1, domain blocks have level>=0
    !rebuild the links within a node. No inter nodes links
    type(blockdef), pointer :: blk,blk_pre
    integer, dimension(:,:), allocatable :: keys
    integer :: i,j,key_p(3),key(3),level,it,iloc,jloc,idx,idy
    real(8) :: coords(2)
    call domain_block_parameters()
    call new_blank_block(blk_root)
    blk_root%key=keys(1,4:6)
    blk_root%curve='Z'
    blk=>blk_root
    if (nblocks_tree==1) then
        call block_mesh_data(blk)
    else
        do i=2,nblocks_tree
            key=keys(i,4:6)
            key_p=keys(i,1:3)
            if (key(3)>blk%key(3)) then
                idx=2-mod(key(1),2)
                idy=2-mod(key(2),2)
                call new_children(blk,idx,idy)
            else
                it=1+blk%key(3)-key(3)
                do j=1,it
                    blk=>blk%pa
                end do
                idx=2-mod(key(1),2)
                idy=2-mod(key(2),2)
                call new_children(blk,idx,idy)
            end if
            blk=>blk%children(idx,idy)%blk
            if (key(3)>=llevel_max) then
                call block_mesh_data(blk)
            else
                blk%level=-1
            end if
        end do
    end if
    call link_children(blk_root)
    call build_linked_list()
    call load_balancer_static()
    call blk_global_traversal(link_neighbors)
end subroutine rebuild_tree_list

subroutine restart_static_level()
    !only modify refinement region based on refinement zones
    type(blockdef), pointer :: blk
    integer :: i,static_level,ierr,j
    real(8) :: zone(4)
    call identify_on_processor()
    if (refine_type=='mixed'.or.refine_type=='static') then
        do i=1,nrefine_region
            zone=refine_region(i,1:4)
            static_level=int(refine_region(i,5))
            call restart_static_refine(blk_root,zone,static_level)
        end do
        call sync_load()
        call relink_neighbors_tree()
        call reset_amr_processor_variables()
        call add_buffer_zone()
        call relink_neighbors_tree()
        call restart_static_derefine()
        call load_balancer_dynamic()
        !call relink_neighbors_tree()
        call blk_traversal(allocate_implicit_mesh)
        call assemble_communication_pattern()
        call communicate_hydro()
        if (iradiation==4) call communicate_fld()
    else if (refine_type=='none') then
        call blk_traversal(allocate_implicit_mesh)
    end if
end subroutine restart_static_level

recursive subroutine restart_static_refine(blk,zone,static_level)
    type(blockdef), pointer :: blk
    real(8) :: zone(4)
    integer :: static_level,idx,idy
    if (associated(blk%children(1,1)%blk)) then
        blk%static_level=blk%level
        do idy=1,2**(nd-1)
            do idx=1,2
                if (associated(blk%children(idx,idy)%blk)) then
                    call restart_static_refine(blk%children(idx,idy)%blk,zone,static_level)
                end if
            end do
        end do
    else
        blk%static_level=static_level
        if (block_zone_overlap(blk,zone).and.blk%level<static_level) then
            blk%static_level=blk%level
            call grow_smr_node(blk,zone,static_level)
            call restart_block_prolongation(blk)
        end if
    end if
end subroutine restart_static_refine

recursive subroutine restart_block_prolongation(blk)
    type(blockdef), pointer :: blk,blk_temp
    integer :: idx,idy
    if (blk%on_processor) then
        if (associated(blk%children(1,1)%blk)) then
            call block_prolongation(blk)
            blk_temp=>blk%children(1,1)%blk
            do idy=1,2**(nd-1)
                do idx=1,2
                    call restart_block_prolongation(blk%children(idx,idy)%blk)
                end do
            end do
        end if
    end if
end subroutine restart_block_prolongation

subroutine restart_static_derefine()
    !use a single core to traverse all the blks and derefine 1 level where necessary
    type(blockdef), pointer :: blk
    integer :: i
    call blk_traversal(examine_static_derefine)
    call sync_blk_attribute_logical(get_derefine,assign_derefine)
    call send_recv_cross_processor_derefine()
    call blk_traversal(derefine_amr_node)
    call identify_on_processor()
    call sync_the_tree_derefine()
    call relink_neighbors_tree()
end subroutine restart_static_derefine

subroutine examine_static_derefine(blk)
    type(blockdef), pointer :: blk
    integer :: i
    real(8) :: zone(4)
    blk%derefine=.false.
    do i=1,nderefine_region
        zone=derefine_region(i,1:4)
        if (block_zone_inside(blk,zone).and.blk%level>blk%static_level) then
            if (.not.level_diff_refine(blk)) then
                blk%derefine=.true.
            end if
        end if
    end do
end subroutine examine_static_derefine

recursive subroutine reallocate_blocks(blk)
    type(blockdef), pointer :: blk
    integer :: idx,idy
    logical :: child
    child=.false.
    outer: do idy=1,2**(nd-1)
        do idx=1,2
            if (associated(blk%children(idx,idy)%blk)) then
                call reallocate_blocks(blk%children(idx,idy)%blk)
                child=.true.
            !else
            !    call allocate_block_heavy_data(blk)
            !    exit outer
            end if
        end do
    end do outer
    if (child.eqv..false.) then
        call allocate_block_heavy_data(blk)
    end if
end subroutine reallocate_blocks

recursive subroutine clean_tree_heavy_data(blk)
    type(blockdef), pointer :: blk
    integer :: idx,idy
    if (associated(blk%children(1,1)%blk)) then
        if (allocated(blk%w)) call release_block(blk)
        do idy=1,2**(nd-1)
            do idx=1,2
                if (associated(blk%children(idx,idy)%blk)) call clean_tree_heavy_data(blk%children(idx,idy)%blk)
            end do
        end do
    end if
end subroutine clean_tree_heavy_data

subroutine read_heavy_data(file_id)
    !master rank read all the data and store it to its blocks
    type(blockdef), pointer :: blk
    integer(HSIZE_T) :: dims(2)
    integer(HID_T) :: file_id,dset_id
    integer :: i,j,error,ndset,ierr
    character(len=16) :: fmt1,s1
    character(len=32) :: groupname,dsetname(100)
    blk=>blk_head
    dims=(/blk_size_nx,blk_size_ny/)
    fmt1='(I5.5)'
    ndset=4
    dsetname(1)='rho';dsetname(2)='vx';dsetname(3)='pres';dsetname(4)='temp'
    if (iradiation==4) then
        dsetname(ndset+1)='Erad'
        ndset=ndset+1
    end if
    do i=1,nblk_total
        write(s1,fmt1) blk%blk_id
        groupname='blk'//trim(s1)
        do j=1,ndset
            call read_a_data(file_id,groupname,dsetname(j),blk,dims)
        end do
#if         imodify==1
        call restart_modify_data(blk)
#endif
        blk=>blk%next
    end do
    nullify(blk)
end subroutine read_heavy_data

subroutine distribute_heavy_data()
    !rank0 distribute the blocks to other blocks according to the load balancer
    type(blockdef), pointer :: blk
    integer :: iblk,blk_id,recv_rank,i,ierr,istart,ireq,ndist
    integer, allocatable :: reqs(:)
    real(8), dimension(:), allocatable :: msg
    if (np>1) then
        if (rank==0) then
            istart=1
            allocate(msg(blksize*(nblk_total-llist_tail%blk_id)))
            blk=>llist_tail%next
            !do while (associated(blk))
            ndist=nblk_total-llist_tail%blk_id
            ireq=1;allocate(reqs(ndist))
            do iblk=1,ndist
                blk_id=blk%blk_id
                call block_id_to_processor_rank(np_nblk,blk_id,recv_rank)
                call send_block(blk,recv_rank,istart,msg,reqs,ireq)
                blk=>blk%next
            end do
        end if
        if (rank/=0) then
            blk=>llist_head
            do i=1,np_nblk(rank+1)
                call receive_block(blk,0)
                blk=>blk%next
            end do
        end if
        if (allocated(msg)) then
            call mpi_waitall(ndist,reqs,MPI_STATUSES_IGNORE,ierr)
            deallocate(msg)
        end if
    end if
end subroutine distribute_heavy_data

subroutine read_a_data(file_id,groupname,dsetname,blk,dims)
    type(blockdef), pointer :: blk
    integer(HID_T) :: file_id,group_id,dset_id
    integer(HSIZE_T) :: dims(2)
    character(len=32) :: groupname,dsetname
    integer :: error
    real(4), allocatable :: x(:,:)
    allocate(x(blk_size_nx,blk_size_ny))
    call h5gopen_f(file_id,groupname,group_id,error)
    call h5dopen_f(group_id,trim(dsetname),dset_id,error)
    call h5dread_f(dset_id,h5t_native_real,x,dims,error)
    select case (dsetname)
    case('rho')
        blk%w(irho,1:blk_size_nx,1:blk_size_ny,1)=x
    case('vx')
        blk%w(ivx,1:blk_size_nx,1:blk_size_ny,1)=x
    case('vy')
        blk%w(ivy,1:blk_size_nx,1:blk_size_ny,1)=x
    case('pres')
        blk%w(ipres,1:blk_size_nx,1:blk_size_ny,1)=x
    case('temp')
        blk%temp(1:blk_size_nx,1:blk_size_ny,1)=x
    case('Erad')
        blk%Erad(1:blk_size_nx,1:blk_size_ny,1)=x
    case('omega')
        blk%omega(1:blk_size_nx,1:blk_size_ny,1)=x
    end select
    call h5dclose_f(dset_id,error)
    call h5gclose_f(group_id,error)
    deallocate(x)
end subroutine read_a_data

subroutine restart_complete_all_data()
    !calculate the derived physical quantities
    !fill all guard cells of all blocks
    type(blockdef), pointer :: blk
    integer :: i,j,iblk,ierr,ip_am
    real(8) :: rho,temp,egv,w(5),u(5),omega
    if (lam_con) then
        ip_am=ipassive(pn_am)
    end if
    blk=>llist_head
    do iblk=1,np_nblk(rank+1)
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                rho=blk%w(1,i,j,1)
                temp=blk%temp(i,j,1)
                if (lam_con) then
                    omega=blk%omega(i,j,1)
                    blk%passive_scalar(ip_am,i,j,1)=blk%l_omega(i,j,1)*omega*rho
                end if
                w=blk%w(1:5,i,j,1)
                call eos_wtou(w,u,temp,egv)
                blk%egv(i,j,1)=egv
                blk%u(1:5,i,j,1)=u
            end do
        end do
        blk=>blk%next
    end do
    nullify(blk)
    call assemble_communication_pattern()
    call communicate_hydro()
    if (iradiation==4) call communicate_fld()
end subroutine restart_complete_all_data

end module io_in

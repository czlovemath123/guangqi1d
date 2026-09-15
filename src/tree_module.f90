module tree_module
!tree and linked list module
use datastructure

implicit none

contains

subroutine build_logical_tree()
    !decompose the computational domain into many small blocks
    call new_blank_block(blk_root)
    blk_root%key=(/1,1,0/)
    blk_root%curve='Z'
    call base_block_binary(blk_root)
    call link_children(blk_root)
end subroutine build_logical_tree

subroutine new_children(blk,idx,idy,curve)
    !given a parent blk, and the location of a child (idx,idy)
    !create the child block
    type(blockdef), pointer :: blk,blk_child
    type(blockpointer), pointer :: blkptr
    integer :: idx,idy
    character(len=1), optional :: curve
    call new_blank_block(blk%children(idx,idy)%blk)
    if (present(curve)) then
        blk%children(idx,idy)%blk%curve=curve
    else
        blk%children(idx,idy)%blk%curve=curve_type(blk%curve,idx,idy)
    end if
    blk_child=>blk%children(idx,idy)%blk
    blk_child%idx=idx;blk_child%idy=idy
    blk_child%key=(/2*blk%key(1)-2+idx,2*blk%key(2)-2+idy,blk%key(3)+1/)
    blk_child%pa=>blk
end subroutine new_children

recursive subroutine base_block_binary(blk)
    !calculate blk%key, link the local child blocks
    type(blockdef), pointer :: blk
    integer :: llevel,level_size,level_size_next,idx,idy
    llevel=blk%key(3)
    level_size=blk_size_nx*2**(llevel_max-llevel)
    idy=1
    if (llevel<llevel_max) then
        level_size_next=level_size/2
        if (nx>(blk%key(1)*level_size-level_size_next)) then
            do idx=1,2
                call new_children(blk,idx,idy)
                call base_block_binary(blk%children(idx,idy)%blk)
            end do
        else
            call new_children(blk,1,1)
            call base_block_binary(blk%children(1,1)%blk)
        end if
    end if
end subroutine base_block_binary

recursive subroutine base_block_quad(blk)
    type(blockdef), pointer :: blk
    integer :: i,llevel,level_size(2),level_size_next(2),idx,idy
    llevel=blk%key(3)
    level_size(1)=blk_size_nx*2**(llevel_max-llevel)
    level_size(2)=blk_size_ny*2**(llevel_max-llevel)
    if (llevel<llevel_max) then
        level_size_next=level_size/2
        if (nx>(blk%key(1)*level_size(1)-level_size_next(1)).and.ny>(blk%key(2)*level_size(2)-level_size_next(2))) then
            do i=1,4
                idx=curve_seq_idx(blk%curve,i)
                idy=curve_seq_idy(blk%curve,i)
                call new_children(blk,idx,idy)
                call base_block_quad(blk%children(idx,idy)%blk)
            end do
        else if (nx>(blk%key(1)*level_size(1)-level_size_next(1)).and.ny<=(blk%key(2)*level_size(2)-level_size_next(2))) then
            do i=1,4
                idx=curve_seq_idx(blk%curve,i)
                idy=curve_seq_idy(blk%curve,i)
                if (idy==1) then        !same as the parent curve type
                    call new_children(blk,idx,idy,blk%curve)
                    call base_block_quad(blk%children(idx,idy)%blk)
                end if
            end do
        else if (nx<=(blk%key(1)*level_size(1)-level_size_next(1)).and.ny>(blk%key(2)*level_size(2)-level_size_next(2))) then
            do i=1,4
                idx=curve_seq_idx(blk%curve,i)
                idy=curve_seq_idy(blk%curve,i)
                if (idx==1) then        !same as the parent curve type
                    call new_children(blk,idx,idy,blk%curve)
                    call base_block_quad(blk%children(idx,idy)%blk)
                end if
            end do
        else if (nx<=(blk%key(1)*level_size(1)-level_size_next(1)).and.ny<=(blk%key(2)*level_size(2)-level_size_next(2))) then
            call new_children(blk,1,1)
            call base_block_quad(blk%children(1,1)%blk)
        end if
    end if
end subroutine base_block_quad

recursive subroutine reassign_curve_type(blk)
    type(blockdef), pointer :: blk
    integer :: i,llevel,level_size(2),level_size_next(2),idx,idy
    llevel=blk%key(3)
    level_size(1)=blk_size_nx*2**(llevel_max-llevel)
    level_size(2)=blk_size_ny*2**(llevel_max-llevel)
    if (llevel<llevel_max) then
        level_size_next=level_size/2
        if (nx>(blk%key(1)*level_size(1)-level_size_next(1)).and.ny>(blk%key(2)*level_size(2)-level_size_next(2))) then
            do i=1,4
                idx=curve_seq_idx(blk%curve,i)
                idy=curve_seq_idy(blk%curve,i)
                blk%children(idx,idy)%blk%curve=curve_type(blk%curve,idx,idy)
                call reassign_curve_type(blk%children(idx,idy)%blk)
            end do
        else if (nx>(blk%key(1)*level_size(1)-level_size_next(1)).and.ny<=(blk%key(2)*level_size(2)-level_size_next(2))) then
            do i=1,4
                idx=curve_seq_idx(blk%curve,i)
                idy=curve_seq_idy(blk%curve,i)
                if (idy==1) then        !same as the parent curve type
                    blk%children(idx,idy)%blk%curve=blk%curve
                    call reassign_curve_type(blk%children(idx,idy)%blk)
                end if
            end do
        else if (nx<=(blk%key(1)*level_size(1)-level_size_next(1)).and.ny>(blk%key(2)*level_size(2)-level_size_next(2))) then
            do i=1,4
                idx=curve_seq_idx(blk%curve,i)
                idy=curve_seq_idy(blk%curve,i)
                if (idx==1) then        !same as the parent curve type
                    blk%children(idx,idy)%blk%curve=blk%curve
                    call reassign_curve_type(blk%children(idx,idy)%blk)
                end if
            end do
        else if (nx<=(blk%key(1)*level_size(1)-level_size_next(1)).and.ny<=(blk%key(2)*level_size(2)-level_size_next(2))) then
            blk%children(1,1)%blk%curve=blk%curve
            call reassign_curve_type(blk%children(1,1)%blk)
        end if
    else if (llevel>=llevel_max) then
        do idy=1,2**(nd-1)
            do idx=1,2
                if (associated(blk%children(idx,idy)%blk)) then
                    blk%children(idx,idy)%blk%curve=curve_type(blk%curve,idx,idy)
                    call reassign_curve_type(blk%children(idx,idy)%blk)
                end if
            end do
        end do
    end if
end subroutine reassign_curve_type

recursive subroutine unlink_children(blk)
    !nullify all the children same level links
    type(blockdef), pointer :: blk
    integer :: i,idx,idy
    do i=1,4
        idx=curve_seq_idx(blk%curve,i)
        idy=curve_seq_idy(blk%curve,i)
        if (associated(blk%children(idx,idy)%blk)) then
            nullify(blk%children(idx,idy)%blk%next,blk%children(idx,idy)%blk%pre)
            call unlink_children(blk%children(idx,idy)%blk)
        end if
    end do
end subroutine unlink_children

recursive subroutine link_children(blk)
    !link the children within the same one-level parent node
    type(blockdef), pointer :: blk,blk1,blk2
    integer :: i,idx,idy
    nullify(blk1,blk2)
    do i=1,4
        idx=curve_seq_idx(blk%curve,i)
        idy=curve_seq_idy(blk%curve,i)
        if (associated(blk%children(idx,idy)%blk)) then
            if (associated(blk1)) then
                blk2=>blk%children(idx,idy)%blk
            else
                blk1=>blk%children(idx,idy)%blk
            end if
            if (associated(blk1).and.associated(blk2)) then
                blk1%next=>blk2
                blk2%pre=>blk1
                blk1=>blk2
            end if
        end if
    end do
    !need to modify here, (1,2) and (2,1) situation same children and parent type
    do i=1,4
        idx=curve_seq_idx(blk%curve,i)
        idy=curve_seq_idy(blk%curve,i)
        if (associated(blk%children(idx,idy)%blk)) then
            call link_children(blk%children(idx,idy)%blk)
        end if
    end do
end subroutine link_children

subroutine build_linked_list()
    !build the arbitrary linked list
    type(blockdef), pointer :: blk_builder1,blk_builder2,blk
    integer :: i
    call find_the_first_node(blk_root,blk_head)
    blk_builder1=>blk_head
    call find_the_successor_node(blk_builder1,blk_builder2)
    do while (.not.associated(blk_builder2,blk_root))
        blk_builder1%next=>blk_builder2
        blk_builder2%pre=>blk_builder1
        blk_builder1=>blk_builder2
        call find_the_successor_node(blk_builder1,blk_builder2)
    end do
    call find_the_last_node()
    call renumber_domain_blocks()
    nullify(blk_builder1,blk_builder2)
end subroutine build_linked_list

recursive subroutine find_the_first_node(blk1,blk2)
    !find the first node of branch blk1, store in blk2
    !it is based on the sequenced children relation
    type(blockdef), pointer :: blk1,blk2
    integer :: i,idx,idy
    logical :: found
    found=.false.
    do i=1,4
        idx=curve_seq_idx(blk1%curve,i)
        idy=curve_seq_idy(blk1%curve,i)
        if (associated(blk1%children(idx,idy)%blk)) then
            call find_the_first_node(blk1%children(idx,idy)%blk,blk2)
            found=.true.
            exit
        end if
    end do
    if (.not.found) blk2=>blk1      !no children node
end subroutine find_the_first_node

subroutine find_the_last_node()
    !next until none
    type(blockdef), pointer :: blk
    blk=>blk_head
    do while(associated(blk%next))
        blk=>blk%next
    end do
    blk_tail=>blk
end subroutine find_the_last_node

subroutine find_the_successor_node(blk1,blk2)
    type(blockdef), pointer :: blk1,blk2,blk_temp
    integer :: key(3)
    if (associated(blk1%next)) then                     !the next logical block is in this level
        call find_the_first_node(blk1%next,blk2)
    else
        if (associated(blk1,blk_root)) then
            blk2=>blk1
        else
            blk2=>blk1%pa                               !go to the parent levels to find the next logical block
            do while (.not.associated(blk2%next))
                if (associated(blk2,blk_root)) exit
                blk2=>blk2%pa
            end do
            if (.not.associated(blk2,blk_root)) call find_the_first_node(blk2%next,blk2)
        end if
    end if
end subroutine find_the_successor_node

subroutine linked_list_mesh_data()
    !calculate blk%coords, and blk%mesh
    type(blockdef), pointer :: blk
    call domain_block_parameters()
    call blk_global_traversal(block_mesh_data)
end subroutine linked_list_mesh_data

subroutine smr_full_domain()
    !every core keeps a copy of the full tree structure
    type(blockdef), pointer :: blk
    real(8) :: zone(4)
    integer :: i,static_level
    if (refine_type=='static'.or.refine_type=='mixed') then
        do i=1,nrefine_region
            zone=refine_region(i,1:4)
            static_level=int(refine_region(i,5))
            blk=>blk_head
            do while (associated(blk))
                if (block_zone_overlap(blk,zone).and.blk%level<static_level) then
                    call grow_smr_node(blk,zone,static_level)
                end if
                blk=>blk%next
            end do
        end do
        call renumber_domain_blocks()
        do i=1,max_refine_level-1
            call blk_global_traversal(smr_buffer_blocks)
        end do
        call renumber_domain_blocks()
    end if
end subroutine smr_full_domain

recursive subroutine grow_smr_node(blk,zone,static_level)
    !recursively create deeper levels of blocks on all processors
    !does not allocate heavy data
    type(blockdef), pointer :: blk
    real(8) :: coords(2),zone(4)
    integer :: id,level,static_level,idx,idy
    call grow_tree_node(blk)
    level=blk%level+1
    do idy=1,2**(nd-1)
        do idx=1,2
            blk%children(idx,idy)%blk%static_level=level
            if (level<static_level) then
                if (block_zone_overlap(blk%children(idx,idy)%blk,zone)) then
                    call grow_smr_node(blk%children(idx,idy)%blk,zone,static_level)
                end if
            end if
        end do
    end do
end subroutine grow_smr_node

subroutine maintain_blk_llist_grow(blk)
    !after growing a tree node, maintain the order and structure of the children blocks and linked list
    type(blockdef), pointer :: blk
    integer :: idx,idy
    idx=curve_seq_idx(blk%curve,1)
    idy=curve_seq_idy(blk%curve,1)
    if (associated(blk,llist_head)) then
        llist_head=>blk%children(idx,idy)%blk
    end if
    if (associated(blk%pre)) then
        blk%pre%next=>blk%children(idx,idy)%blk
        blk%children(idx,idy)%blk%pre=>blk%pre
    else            !this is the head block
        blk_head=>blk%children(idx,idy)%blk
    end if
    idx=curve_seq_idx(blk%curve,2**nd)
    idy=curve_seq_idy(blk%curve,2**nd)
    if (associated(blk,llist_tail)) then
        llist_tail=>blk%children(idx,idy)%blk
    end if
    if (associated(blk%next)) then
        blk%next%pre=>blk%children(idx,idy)%blk
        blk%children(idx,idy)%blk%next=>blk%next
    else            !this is the tail block
        blk_tail=>blk%children(idx,idy)%blk
    end if
    call link_children(blk)
end subroutine maintain_blk_llist_grow

subroutine maintain_blk_llist_trim(blk) 
    !blk is the parent level block
    type(blockdef), pointer :: blk,blk_pre,blk_next
    integer :: idx,idy
    idx=curve_seq_idx(blk%curve,1)
    idy=curve_seq_idy(blk%curve,1)
    if (associated(blk%children(idx,idy)%blk,llist_head)) llist_head=>blk
    if (associated(blk%children(idx,idy)%blk%pre)) then
        blk_pre=>blk%children(idx,idy)%blk%pre
        blk_pre%next=>blk
        blk%pre=>blk_pre
    else            !this is the first block
        blk_head=>blk
    end if
    idx=curve_seq_idx(blk%curve,2**nd)
    idy=curve_seq_idy(blk%curve,2**nd)
    if (associated(blk%children(idx,idy)%blk,llist_tail)) llist_tail=>blk
    if (associated(blk%children(idx,idy)%blk%next)) then
        blk_next=>blk%children(idx,idy)%blk%next
        blk%next=>blk_next
        blk_next%pre=>blk
    else            !this is the last block
        blk_tail=>blk
    end if
end subroutine maintain_blk_llist_trim

recursive subroutine unlink_neighbor_blocks(blk)
    !remove the links between blocks
    type(blockdef), pointer :: blk
    integer :: i,idx,idy
    do i=1,12
        blk%nb(i)%cpu_rank=-1
        nullify(blk%nb(i)%blk)
    end do
    do idy=1,2**(nd-1)
        do idx=1,2
            if (associated(blk%children(idx,idy)%blk)) call unlink_neighbor_blocks(blk%children(idx,idy)%blk)
        end do
    end do
end subroutine unlink_neighbor_blocks

subroutine link_neighbors(blk)
    type(blockdef), pointer :: blk
    type(blockpointer) :: blk_nb
    character(len=16) :: direction
    integer :: i
    blk%nb_level_1d=0
    call find_left_neighbor(blk)
    call find_right_neighbor(blk)
end subroutine link_neighbors

subroutine relink_neighbors_tree()
    call unlink_neighbor_blocks(blk_root)
    call blk_global_traversal(link_neighbors)
end subroutine relink_neighbors_tree

subroutine grow_tree_node(blk)
    type(blockdef), pointer :: blk
    integer :: idx,idy
    idx=curve_seq_idx(blk%curve,1)
    idy=curve_seq_idy(blk%curve,1)
    if (.not.associated(blk%children(idx,idy)%blk)) then
        do idy=1,2**(nd-1)
            do idx=1,2
                call new_children(blk,idx,idy)
                blk%children(idx,idy)%blk%on_processor=blk%on_processor
                call block_mesh_data(blk%children(idx,idy)%blk)
            end do
        end do
        call maintain_blk_llist_grow(blk)
    end if
end subroutine grow_tree_node

subroutine trim_tree_node(blk)
    !trim when blk's heavy data is not on this processor
    type(blockdef), pointer :: blk
    integer :: idx,idy
    idx=curve_seq_idx(blk%curve,1)
    idy=curve_seq_idy(blk%curve,1)
    if (associated(blk%children(idx,idy)%blk)) then
        call maintain_blk_llist_trim(blk)
        do idy=1,2**(nd-1)
            do idx=1,2
                nullify(blk%children(idx,idy)%blk)
            end do
        end do
    end if
end subroutine trim_tree_node

subroutine smr_buffer_blocks(blk)
    !block meta data, smr only
    type(blockdef), pointer :: blk,blk_temp
    integer :: level,i
    level=blk%level
    if (associated(blk%pre)) then
        blk_temp=>blk%pre
        if (blk_temp%level+1<level) then
            call grow_tree_node(blk_temp)
        end if
    end if
    if (associated(blk%next)) then
        blk_temp=>blk%next
        if (blk_temp%level+1<level) then
            call grow_tree_node(blk_temp)
        end if
    end if
end subroutine smr_buffer_blocks

subroutine find_left_neighbor(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    if (blk%key(1)/=1) then
        blk_temp=>blk%pre
        blk_nb=blk%nb(1)
        call link_nb(blk%nb(1),blk_temp)
        blk%nb_level_1d(1)=blk_temp%level-blk%level
    else
        blk%nb_level_1d(1)=0
    end if
end subroutine find_left_neighbor

subroutine find_right_neighbor(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    if (blk%key(1)/=nx_blks*2**blk%level) then
        blk_temp=>blk%next
        blk_nb=blk%nb(2)
        call link_nb(blk%nb(2),blk_temp)
        blk%nb_level_1d(2)=blk_temp%level-blk%level
    else
        blk%nb_level_1d(2)=0
    end if
end subroutine find_right_neighbor

subroutine link_nb(blk_nb,blk)
    type(blockdef), pointer :: blk
    type(blockpointer) :: blk_nb
    integer :: blk_id
    if (associated(blk_nb%blk)) nullify(blk_nb%blk)
    blk_nb%blk=>blk
    blk_id=blk%blk_id
    call block_id_to_processor_rank(np_nblk,blk_id,blk_nb%cpu_rank)
end subroutine link_nb

recursive subroutine deepest_directional_child_level(blk,direction,dp_level)
    !return the deepest children block level within blk and the given direction
    type(blockdef), pointer :: blk
    integer :: dp_level
    character(len=16) :: direction
    if (blk%level>dp_level) dp_level=blk%level
    if (direction==next) then
        if (associated(blk%children(2,1)%blk)) call deepest_directional_child_level(blk%children(2,1)%blk,direction,dp_level)
    else if (direction==previous) then
        if (associated(blk%children(1,1)%blk)) call deepest_directional_child_level(blk%children(1,1)%blk,direction,dp_level)
    end if
end subroutine deepest_directional_child_level

recursive subroutine deepest_child_level(blk,dp_level)
    !return the deepest children block level within blk
    type(blockdef), pointer :: blk
    integer :: dp_level,idx,idy
    if (blk%level>dp_level) dp_level=blk%level
    do idy=1,2**(nd-1)
        do idx=1,2
            if (associated(blk%children(idx,idy)%blk)) call deepest_child_level(blk%children(idx,idy)%blk,dp_level)
        end do
    end do
end subroutine deepest_child_level

recursive subroutine display_tree(blk)
    type(blockdef), pointer :: blk
    logical :: l1,l2,l3,l4
    integer :: idx,idy
    print *,rank,blk%key
    do idy=1,2**(nd-1)
        do idx=1,2
            if (associated(blk%children(idx,idy)%blk)) call display_tree(blk%children(idx,idy)%blk)
        end do
    end do
end subroutine display_tree

end module tree_module

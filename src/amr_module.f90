module amr_module
use load_module
use communication
use boundary
implicit none

contains

subroutine grow_amr_tree()
    !all cores have a copy of the whole tree
    !we use action_countdown to avoid a refinement immediately after a derefine
    type(blockdef), pointer :: blk
    integer :: i,j,ierr
    if (refine_type=='adaptive'.or.refine_type=='mixed') then
        change_mesh=.false.
        call communicate_hydro()
        call intrinsic_amr()
        call add_buffer_zone()
        !call add_directional_refine()
        !call add_buffer_zone()
        call derefine()
        call load_balancer_dynamic()
        call check_flying_blk_data()
        !call release_all_not_on_processor(blk_root)
        if (change_mesh) then
            if (iradiation==4) call blk_traversal(allocate_implicit_mesh)
            call assemble_communication_pattern()
            call communicate_hydro()
        end if
    end if
end subroutine grow_amr_tree

!refine the intrinsic nodes

subroutine intrinsic_amr()
    call reset_amr_processor_variables()
    call blk_traversal(grow_amr_node)
    call sync_the_tree_refine()
    call sync_blk_attribute_logical(get_amr_intrinsic,assign_amr_intrinsic)
end subroutine intrinsic_amr

recursive subroutine grow_amr_node(blk)
    !recursively grow amr node, refine the intrinsic refine zone
    type(blockdef), pointer :: blk
    integer :: idx,idy
    if (blk%level<max_refine_level) then
        if (refinement_criterion(blk)) then
            call grow_tree_node(blk)
            call record_processor_amr_grow_order(blk)
            call block_prolongation(blk)
            do idy=1,2**(nd-1)
                do idx=1,2
                    call grow_amr_node(blk%children(idx,idy)%blk)
                end do
            end do
        end if
    end if
end subroutine grow_amr_node

!add buffer nodes until it converges (extrinsic refinement)

subroutine add_buffer_zone()
    integer :: i
    do i=1,max_refine_level-1
        call add_buffer_onelevel()
    end do
    call relink_neighbors_tree()
end subroutine add_buffer_zone

subroutine add_buffer_onelevel()
    call blk_traversal(add_amr_buffer_blocks)
    call sync_the_tree_refine()
end subroutine add_buffer_onelevel

subroutine add_amr_buffer_blocks(blk)
    !see if blk's neighbours has more then 2 levels higher refinement
    !if so, refine this blk
    type(blockdef), pointer :: blk
    integer :: i,level,dp_level,key(3),idx,idy
    character(len=16) :: dir(2)
    dir=(/next,previous/)
    level=blk%level
    do i=1,2
        call nb_deepest_level(blk,dir(i),dp_level)
        if (dp_level>level+1) then
            call grow_amr_tree_node_buffer(blk)
            blk=>blk%children(2,1)%blk
            exit
        end if
    end do
end subroutine add_amr_buffer_blocks

subroutine grow_amr_tree_node_buffer(blk)
    type(blockdef), pointer :: blk
    integer :: id,level
    call grow_tree_node(blk)
    call record_processor_amr_grow_order(blk)
    call block_prolongation(blk)
end subroutine grow_amr_tree_node_buffer

!refine the nodes that are due to wave pheonomenon (wave transported to the edge of the blocks)

subroutine add_directional_refine()
    call blk_traversal(label_directional_refine)
    call sync_directional_refine()
    call blk_traversal_ind(grow_directional_refine)
    call sync_the_tree_refine()
    call add_buffer_zone()
end subroutine add_directional_refine

subroutine label_directional_refine(blk)
    type(blockdef), pointer :: blk
    integer :: i,j,key(3)
    logical :: l1,l2,l3
    real(8) :: tg,Erad,trad
    blk%refine_direction=.false.
    do i=1,directional_buffer
        blk%refine_direction(1)=directional_refine(blk,i,1).or.blk%refine_direction(1)
        blk%refine_direction(2)=directional_refine(blk,blk_size_nx-i+1,1).or.blk%refine_direction(2)
    end do
end subroutine label_directional_refine

function directional_refine(blk,i,j)
    type(blockdef), pointer :: blk
    logical :: l1,l2,l3,directional_refine
    integer :: i,j
    l1=.false.;l2=.false.;l3=.false.
    l1=shock_refine(blk,xcoord,i,j,gradp_refine_thresh)
    if (refine_zeldovich.and.iradiation==4) l2=zeldovich_spike(blk,i,j)
    if (iradiation==4) l3=slope_refine(blk,strad,xcoord,i,j,erad_slope_refine_thresh)
    directional_refine=l1.or.l2.or.l3
end function directional_refine

subroutine grow_directional_refine(blk)
    type(blockdef), pointer :: blk
    type(blockpointer) :: blk_nb
    integer :: i,ct
    character(len=16) :: direction
    do i=1,12
        direction=nb_index_direc(i)
        blk_nb=blk%nb(i)
        if (associated(blk_nb%blk)) then
            if (blk_nb%blk%level>blk%level) then
                if (direction==east) then
                    if (blk_nb%blk%refine_direction(3)) then
                        call grow_amr_tree_node_buffer(blk)
                        exit
                    end if
                else if (direction==south) then
                    if (blk_nb%blk%refine_direction(4)) then
                        call grow_amr_tree_node_buffer(blk)
                        exit
                    end if
                else if (direction==west) then
                    if (blk_nb%blk%refine_direction(1)) then
                        call grow_amr_tree_node_buffer(blk)
                        exit
                    end if
                else if (direction==north) then
                    if (blk_nb%blk%refine_direction(2)) then
                        call grow_amr_tree_node_buffer(blk)
                        exit
                    end if
                else
                    ct=corner_type(blk,direction)
                    select case(ct)
                    case (1)
                        if (hydro_bound_type(1)==3.and.hydro_bound_type(4)==3) then
                            if (blk_nb%blk%refine_direction(5)) then
                                call grow_amr_tree_node_buffer(blk)
                                exit
                            end if
                        end if
                    case (2)
                        if (hydro_bound_type(4)==3) then
                            if (direction==ne) then
                                if (blk_nb%blk%refine_direction(6)) then
                                    call grow_amr_tree_node_buffer(blk)
                                    exit
                                end if
                            else if (direction==nw) then
                                if (blk_nb%blk%refine_direction(5)) then
                                    call grow_amr_tree_node_buffer(blk)
                                    exit
                                end if
                            end if
                        end if
                    case (3)
                        if (hydro_bound_type(2)==3.and.hydro_bound_type(4)==3) then
                            if (blk_nb%blk%refine_direction(6)) then
                                call grow_amr_tree_node_buffer(blk)
                                exit
                            end if
                        end if
                    case (4)
                        if (hydro_bound_type(1)==3) then
                            if (direction==nw) then
                                if (blk_nb%blk%refine_direction(5)) then
                                    call grow_amr_tree_node_buffer(blk)
                                    exit
                                end if
                            else if (direction==sw) then
                                if (blk_nb%blk%refine_direction(8)) then
                                    call grow_amr_tree_node_buffer(blk)
                                    exit
                                end if
                            end if
                        end if
                    case (5)
                        if (direction==se) then
                            if (blk_nb%blk%refine_direction(7)) then
                                call grow_amr_tree_node_buffer(blk)
                                exit
                            end if
                        else if (direction==sw) then
                            if (blk_nb%blk%refine_direction(8)) then
                                call grow_amr_tree_node_buffer(blk)
                                exit
                            end if
                        else if (direction==nw) then
                            if (blk_nb%blk%refine_direction(5)) then
                                call grow_amr_tree_node_buffer(blk)
                                exit
                            end if
                        else if (direction==ne) then
                            if (blk_nb%blk%refine_direction(6)) then
                                call grow_amr_tree_node_buffer(blk)
                                exit
                            end if
                        end if
                    case (6)
                        if (hydro_bound_type(2)==3) then
                            if (direction==ne) then
                                if (blk_nb%blk%refine_direction(6)) then
                                    call grow_amr_tree_node_buffer(blk)
                                    exit
                                end if
                            else if (direction==se) then
                                if (blk_nb%blk%refine_direction(7)) then
                                    call grow_amr_tree_node_buffer(blk)
                                    exit
                                end if
                            end if
                        end if
                    case (7)
                        if (hydro_bound_type(1)==3.and.hydro_bound_type(3)==3) then
                            if (blk_nb%blk%refine_direction(8)) then
                                call grow_amr_tree_node_buffer(blk)
                                exit
                            end if
                        end if
                    case (8)
                        if (hydro_bound_type(3)==3) then
                            if (direction==se) then
                                if (blk_nb%blk%refine_direction(7)) then
                                    call grow_amr_tree_node_buffer(blk)
                                    exit
                                end if
                            else if (direction==sw) then
                                if (blk_nb%blk%refine_direction(8)) then
                                    call grow_amr_tree_node_buffer(blk)
                                    exit
                                end if
                            end if
                        end if
                    case (9)
                        if (hydro_bound_type(2)==3.and.hydro_bound_type(3)==3) then
                            if (blk_nb%blk%refine_direction(7)) then
                                call grow_amr_tree_node_buffer(blk)
                                exit
                            end if
                        end if
                    end select
                end if
            end if
        end if
    end do
end subroutine grow_directional_refine

!derefine the nodes

subroutine derefine()
    call examine_derefine()
    call sync_blk_attribute_logical(get_derefine,assign_derefine)
    call send_recv_cross_processor_derefine()
    call blk_traversal(derefine_amr_node)
    call identify_on_processor()
    call sync_the_tree_derefine()
    call relink_neighbors_tree()
end subroutine derefine

subroutine derefine_amr_node(blk)
    !derefine if all its cousins are on processor and to be derefined
    type(blockdef), pointer :: blk,pa
    logical :: l1,l2,l3,l4,l(2,2)
    integer :: i,j,key(3),idy,idx
    l=.true.
    pa=>blk%pa
    do idy=1,2**(nd-1)
        do idx=1,2
            if (associated(pa%children(idx,idy)%blk)) then
                l(idx,idy)=(pa%children(idx,idy)%blk%on_processor.and.pa%children(idx,idy)%blk%derefine)
            end if
        end do
    end do
    if (all(l)) then                !derefine
        call maintain_blk_llist_trim(pa)
        call block_restriction(pa)
        call record_processor_amr_derefine_order(pa)
        blk=>pa
    end if
end subroutine derefine_amr_node

function refinement_criterion(blk)
    !intrinsic refinement, only depend on local data
    type(blockdef), pointer :: blk
    logical :: refinement_criterion,lshock,ldensity,ltrad,l1,l2,luser
    real(8) :: pl,pr,vxl,vxr,egv,Etotal,gradp,trad,xl,xc,xr
    integer :: i,j,k,key(3)
    lshock=.false.;ldensity=.false.;ltrad=.false.;luser=.false.
#if     user_amr
    do i=1,blk_size_nx
        luser=amr_criterion(blk,i,1)
        if (luser) exit
    end do
#else
    do i=1,blk_size_nx
        lshock=shock_refine(blk,xcoord,i,1,gradp_refine_thresh)
        ldensity=slope_refine(blk,sdensity,xcoord,i,1,rho_slope_refine_thresh)
        !if (iradiation==4) l1=zeldovich_spike(blk,i,1)
        if (iradiation==4) ltrad=slope_refine(blk,strad,xcoord,i,1,erad_slope_refine_thresh)
        if (lshock.or.ldensity.or.ltrad) exit
    end do
#endif
    if ((lshock.or.ldensity.or.ltrad.or.luser)) then
        if (blk%derefine_countdown<=0) then
            refinement_criterion=.true.
        else
            refinement_criterion=.false.
        end if
        blk%amr_intrinsic=.true.
    else
        refinement_criterion=.false.
        blk%amr_intrinsic=.false.
    end if
end function refinement_criterion

function shock_refine(blk,coord,i,j,thresh)
    type(blockdef), pointer :: blk
    character(len=16) :: coord
    integer :: i,j
    real(8) :: pl,pr,egv,Etotal,gradp,thresh
    logical :: shock_refine
    shock_refine=.false.
    if (coord==xcoord) then
        pl=blk%w(5,i-1,j,1)
        pr=blk%w(5,i+1,j,1)
    else if (coord==ycoord) then
        pl=blk%w(5,i,j-1,1)
        pr=blk%w(5,i,j+1,1)
    end if
    egv=blk%egv(i,j,1)
    Etotal=blk%u(5,i,j,1)
    gradp=abs(pr-pl)/min(pl,pr)
    if (gradp>thresh.and.egv/Etotal>0.1) shock_refine=.true.
end function shock_refine

function zeldovich_spike(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: tg,thresh,trad
    logical :: zeldovich_spike
    zeldovich_spike=.false.
    tg=blk%temp(i,j,1)
    if (tg>1.2d3) then
        trad=(blk%Erad(i,j,1)/a_rad)**0.25d0
        if (tg>=trad*1.1) zeldovich_spike=.true.
    end if
end function zeldovich_spike

function zeldovich_spike_derefine(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: tg,thresh,trad
    logical :: zeldovich_spike_derefine
    zeldovich_spike_derefine=.false.
    tg=blk%temp(i,j,1)
    if (tg>1.2d3) then
        trad=(blk%Erad(i,j,1)/a_rad)**0.25d0
        if (tg>=trad*1.02) zeldovich_spike_derefine=.true.
    end if
end function zeldovich_spike_derefine

function slope_refine(blk,field,coord,i,j,thresh)
    type(blockdef), pointer :: blk
    character(len=16) :: field,coord
    integer :: i,j
    logical :: slope_refine
    real(8) :: q(3),thresh
    if (field==sdensity) then
        if (coord==xcoord) then
            q=blk%w(1,i-1:i+1,j,1)
        else if (coord==ycoord) then
            q=blk%w(1,i,j-1:j+1,1)
        end if
    else if (field==strad) then
        if (coord==xcoord) then
            q=blk%temp(i-1:i+1,j,1)
        else if (coord==ycoord) then
            q=blk%temp(i,j-1:j+1,1)
        end if
    else if (field==stgas) then
        if (coord==xcoord) then
            q=blk%temp(i-1:i+1,j,1)
        else if (coord==ycoord) then
            q=blk%temp(i,j-1:j+1,1)
        end if
    else if (field==sshear) then
        if (coord==xcoord) then
            q=blk%w(2,i,j-1:j+1,1)*blk%dxyz(1)
        else if (coord==ycoord) then
            q=blk%w(3,i-1:i+1,j,1)*blk%dxyz(1)
        end if
    end if
    if (abs(q(3)-q(1))/2/q(2)>thresh) then
        slope_refine=.true.
    else
        slope_refine=.false.
    end if
end function slope_refine

function level_diff_refine(blk)
    !determine whether blk should be refined based on the level difference condition
    type(blockdef), pointer :: blk,blk_temp
    type(blockpointer) :: blk_nb
    integer :: i,level
    logical :: level_diff_refine
    level_diff_refine=.false.
    do i=1,12
        blk_nb=blk%nb(i)
        if (associated(blk_nb%blk)) then
            blk_temp=>blk_nb%blk
            if ((blk_temp%level-blk%level)>0) then
                level_diff_refine=.true.
                exit
            end if
        end if
    end do
end function level_diff_refine

logical function refine_neighbor(blk)
    type(blockdef), pointer :: blk,blk_temp
    character(len=16) :: direction
    type(blockpointer) :: blk_nb
    integer :: i,ct,rd
    refine_neighbor=.false.
end function refine_neighbor

function derefinement_criterion(blk)
    !do not derefine if:
    !the block contains discontinuities
    !any of the block's neighbour has 1 level higher AMR already
    type(blockdef), pointer :: blk,blk_temp
    integer :: level,i,j,key(3)
    real(8) :: trad
    logical :: discon,level_diff,directional_refine
    logical :: derefinement_criterion,ldensity,l1,l2,l3,lzeldovich
    derefinement_criterion=.true.;discon=.false.;level_diff=.false.;directional_refine=.false.
    l1=.true.;l2=.true.;l3=.true.
    if (blk%level>blk%static_level.and.blk%refine_countdown<=0) then
        level_diff=level_diff_refine(blk)
        discon=blk%amr_intrinsic
        directional_refine=.false.!refine_neighbour(blk)
        if (level_diff.or.discon.or.directional_refine) then
            derefinement_criterion=.false.
            blk%refine_countdown=10
        else
#if     user_amr
            derefinement_criterion=amr_derefine(blk)
#else
            derefinement_criterion=.true.
            if (refine_zeldovich.and.iradiation==4) l1=derefine_zeldovich(blk)
            if (refine_density) l2=derefine_density(blk)
            if (refine_erad.and.iradiation==4) l3=derefine_erad(blk)
            derefinement_criterion=(l1.and.l2.and.l3)
#endif
        end if
    else
        derefinement_criterion=.false.
    end if
end function derefinement_criterion

function derefine_zeldovich(blk)
    type(blockdef), pointer :: blk
    logical :: lzeldovich,derefine_zeldovich
    integer :: i,j
    lzeldovich=.false.
    derefine_zeldovich=.true.
    outer: do j=1,blk_size_ny
        do i=1,blk_size_nx
            lzeldovich=lzeldovich.or.zeldovich_spike_derefine(blk,i,j)
            if (lzeldovich) then
                derefine_zeldovich=.false.
                exit outer
            end if
        end do
    end do outer
end function derefine_zeldovich

function derefine_density(blk)
    type(blockdef), pointer :: blk
    logical :: l1,l2,ldensity,derefine_density
    integer :: i,j
    ldensity=.false.
    derefine_density=.true.
    outer: do j=1,blk_size_ny
        do i=1,blk_size_nx
            l1=slope_refine(blk,sdensity,xcoord,i,j,rho_slope_derefine_thresh)
            l2=slope_refine(blk,sdensity,ycoord,i,j,rho_slope_derefine_thresh)
            ldensity=ldensity.or.l1.or.l2
            if (ldensity) then
                derefine_density=.false.
                exit outer
            end if
        end do
    end do outer
end function derefine_density

function derefine_erad(blk)
    type(blockdef), pointer :: blk
    logical :: l1,l2,lerad,derefine_erad
    integer :: i,j
    lerad=.false.
    derefine_erad=.true.
    outer: do j=1,blk_size_ny
        do i=1,blk_size_nx
            l1=slope_refine(blk,strad,xcoord,i,j,erad_slope_derefine_thresh)
            lerad=lerad.or.l1
            if (lerad) then
                derefine_erad=.false.
                exit outer
            end if
        end do
    end do outer
end function derefine_erad

subroutine examine_derefine()
    !gather all the blk%derefine information on one processor
    !check which blk should derefine and bcast their key to all others
    type(blockdef), pointer :: blk
    integer :: i,j,reqs,ierr,stat(MPI_STATUS_SIZE)
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        blk%derefine=derefinement_criterion(blk)
        blk=>blk%next
    end do
end subroutine examine_derefine

logical function get_amr_intrinsic(blk)
    type(blockdef), pointer :: blk
    get_amr_intrinsic=blk%amr_intrinsic
end function get_amr_intrinsic

subroutine assign_amr_intrinsic(blk,amr_intrinsic)
    type(blockdef), pointer :: blk
    logical :: amr_intrinsic
    blk%amr_intrinsic=amr_intrinsic
end subroutine assign_amr_intrinsic

logical function get_derefine(blk)
    type(blockdef), pointer :: blk
    get_derefine=blk%derefine
end function get_derefine

subroutine assign_derefine(blk,derefine)
    type(blockdef), pointer :: blk
    logical :: derefine
    blk%derefine=derefine
end subroutine assign_derefine

end module amr_module

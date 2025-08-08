program rmhd
use hdf5
use mathlib
use datastructure
use eos
use communication
use amr_module
use phylib
use boundary
use io_in
use io_out
use gravity
use hydroscheme
use problem
use radiation
use cooling
use source_control
use passivescalars
implicit none

integer :: c1,c2,count_rate,count_max
integer :: ierr,errcode
real(8) :: dt_phy,t1,t2,dt_cpu
character(len=32) :: settings(3)
call problem_gen()
call problem_start()
call cpu_time(t1)
call system_clock(c1,count_rate,count_max)
call advan()
call system_clock(c2,count_rate,count_max)
call cpu_time(t2)
dt_phy=real(c2-c1,kind=8)/count_rate
dt_cpu=t2-t1
if (rank==0) then
    print *,'physical time is:',dt_phy
    print *,'cpu time is:',dt_cpu
end if
call finalize()

contains

subroutine problem_gen()
    integer :: i,ierr,errcode
    type(blockdef), pointer :: blk
    call start_mpi()
    call initialize_simulation_parameters()
    call initialize_eos_environment()
    call calculate_llevel_max()
    call initialize_problem()
    call initialize_hydro()
    call initialize_source()
    call initialize_passive_scalar()
    call dsetname_filename()
    call simulation_config()
    if (restart) then
        call hdf_read()
        time_sys%t_output_initial=time_sys%t
        time_sys%dt_frame=(time_sys%t_final-time_sys%t_output_initial)/(nframe-iframe)
        allocate(t_target_output(nframe-iframe))
        do i=1,nframe-iframe
            t_target_output(i)=time_sys%t_output_initial+time_sys%dt_frame*i
        end do
        time_sys%t_next=max(time_sys%t_output_initial+time_sys%dt_frame,time_sys%t_final)
    else
        call build_logical_tree()
        call build_linked_list()
        call linked_list_mesh_data()
        call smr_full_domain()
        call load_balancer_static()
        call blk_traversal(link_neighbors)
        call assemble_communication_pattern()
        time_sys%t=zero
        time_sys%t_output_initial=zero
        time_sys%dt_frame=(time_sys%t_final-time_sys%t_output_initial)/nframe
        allocate(t_target_output(nframe))
        do i=1,nframe
            t_target_output(i)=time_sys%dt_frame*i
        end do
        time_sys%t_next=max(time_sys%t_output_initial+time_sys%dt_frame,time_sys%t_final)
        call blk_traversal(allocate_block_heavy_data)
        if (iradiation==4) call blk_traversal(allocate_implicit_mesh)
    end if
    if (rank==0) print *,'total block number=',nblk_total,'block per core=',real(nblk_total)/real(np)
    call initialize_radiation()
end subroutine problem_gen

subroutine problem_start()
    if (.not.restart) then
        call blk_traversal(apply_hydro_condition)
        call blk_traversal(apply_rad_condition)
        call blk_traversal(apply_passive_init_condition)
        call output_blocks_parallel()
    end if
end subroutine problem_start

subroutine first_step()
end subroutine first_step

!evolve with time. Hydro, radiation, magnetic field, geometry, and gravity are doing their job here.
subroutine advan()
    integer :: i,j,k,ierr,iupdate
    real(8) :: t1,t2,tsample,t3,t4
    real(8), allocatable :: record_array(:)
    character(len=30) :: fmt1
    character(len=5) :: mark
    character(len=128) :: alert
    logical ::converged,file_exists
    i=0;j=0;t1=0;t2=0
    tsample=time_sys%t_output_initial
    iupdate=0
#if     irecord==1
    allocate(record_array(record_length))
    inquire(file="history.data", exist=file_exists)
    if (rank==0) then
        if (file_exists) then
            if (lhis_overwrite) then
                open(unit=88,file='history.data',status='replace',action='write')
                close(88)
            end if
        else
            open(unit=88,file='history.data',status='new',action='write')
            close(88)
        end if
    end if
    time_sys%t_record=time_sys%t
#endif
    call status_quo_columns()
    do while(time_sys%t<time_sys%t_final)
        iframe=iframe+1
        i=i+1
        time_sys%t_next=t_target_output(i)
        call status_quo(t1,t2)
        call cpu_time(t1)
        do while(time_sys%t<time_sys%t_next)
            if (temp_smr) call blk_traversal(decrease_static_level)
            call applyboundconds()
            call grow_amr_tree()
            if (lhydro) call hydro_step()
            if (lsource.and.source_applicability()) then
                call applyboundconds()
                call communicate_hydro()
                call blk_traversal(source_apply)
            end if
            call radiation_transfer()
            call problem_oper()
#if         iconverge==1
            call convergence_test(converged,tsample)
#endif
            time_sys%t=time_sys%t+time_sys%dt_hydro
#if     irecord==1
            if (time_sys%dt_record==0) then
                alert='dt_record=0'
                call abort_guangqi(alert)
            else
                if (time_sys%t-time_sys%t_record>time_sys%dt_record.and.time_sys%ntimestep>0) then
                    time_sys%t_record=time_sys%t
                    call assemble_record_array(record_array)
                    if (rank==0) then
                        open(unit=88,file='history.data',status='old',action='write',POSITION="append")
                        write(unit=88,fmt='(100ES16.8E2)') record_array
                        close(88)
                    end if
                end if
            end if
#endif
            time_sys%ntimestep=time_sys%ntimestep+1
        end do
        call output_blocks_parallel()
        call cpu_time(t2)
    end do
#if     irecord==1
    deallocate(record_array)
#endif
end subroutine advan

subroutine finalize()
    !deallocate all the arrays and do some final output
    integer :: ierr
    call finalize_problem()
    call finalize_radiation()
    call finalize_mpi()
end subroutine finalize

end program rmhd

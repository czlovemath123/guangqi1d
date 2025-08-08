module radiation
#include <petsc/finclude/petscksp.h>
use phylib
use mathlib
use datastructure
use problem
use eos
use petsc_fld
use radiation_common_functions
implicit none

contains

subroutine radiation_transfer()
    integer :: ierr
    if (iradiation==4) then
        !time_sys%tradhydro_end=time_sys%tradhydro_end+time_sys%dt_hydro
        radhydro_tspan=radhydro_tspan+time_sys%dt_hydro
        if (mod(time_sys%ntimestep,radhydro_boost)==0) then
            !petsc flux-limited diffusion
            time_sys%dt_radhydro=radhydro_tspan
            if (refine_type=='static'.or.refine_type=='none') then
                call petsc_fld_solve(axb_fld,ierr)
            else if (refine_type=='adaptive'.or.refine_type=='mixed') then
                if (change_mesh) then
                    call petsc_fld_recycle(axb_fld,ierr)
                    call petsc_fld_initialize(axb_fld,ierr)
                    call blk_traversal(calculate_guard_coord)
                    call blk_traversal(recalculate_g_alpha)
                end if
                call petsc_fld_solve(axb_fld,ierr)
            else if (refine_type=='updated') then
                call petsc_fld_recycle(axb_fld,ierr)
                call petsc_fld_initialize(axb_fld,ierr)
                call petsc_fld_solve(axb_fld,ierr)
                refine_type='static'
            end if
            radhydro_tspan=0d0
        end if
    end if
end subroutine radiation_transfer

subroutine initialize_radiation()
    integer :: ierr
    bound_rad=>boundary_rad
    init_rad=>initial_rad
    if (iradiation==1) then
    else if (iradiation==2) then
    else if (iradiation==3) then
        !call jacobi_init()
    else if (iradiation==4) then
#if     iopacity==3
#else
        call generate_rosseland_planck_kap_table()
#endif
        call petsc_fld_initialize(axb_fld,ierr)
        directional_buffer=2
    else if (iradiation==5) then
    end if
end subroutine initialize_radiation

subroutine finalize_radiation()
    integer :: ierr
    if (iradiation==1) then
    else if (iradiation==2) then
    else if (iradiation==3) then
    else if (iradiation==4) then
        call petsc_fld_finalize(axb_fld,ierr)
    else if (iradiation==5) then
    end if
end subroutine finalize_radiation

end module radiation

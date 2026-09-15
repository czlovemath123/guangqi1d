module petsc_fld
#include <petsc/finclude/petscksp.h>
use datastructure
use eos
use io_out
use communication
use phylib
use boundary
use radiation_common_functions
use petsc_fld_1d
use petscksp
implicit none

logical, protected :: petsc_not_converged

contains 

subroutine petsc_fld_initialize(userctx,ierr)
    PetscErrorCode ierr
    type(axbsystem) userctx
    Mat     A
    Vec     b,x
    KSP    ksp
    PetscInt Ntot,submat_nrow,n_local_rows
    PetscInt, allocatable :: nc(:),d_nz(:),d_nnz(:),o_nz(:),o_nnz(:)
    integer :: i
    cv_thresh=rho_thresh_petsc1*kb/(gamma_gas-1)/maw/amu
    Ntot=2*nblk_total*blk_size_nx
    n_local_rows=2*np_nblk(rank+1)*blk_size_nx
    call MatCreateAIJ(PETSC_COMM_WORLD,n_local_rows,n_local_rows,Ntot,Ntot,4,PETSC_NULL_INTEGER,1,PETSC_NULL_INTEGER,A,ierr)!;petsccall(ierr)
    call VecCreateMPI(PETSC_COMM_WORLD,n_local_rows,Ntot,b,ierr)!;petsccall(ierr)
    call VecDuplicate(b,x,ierr)!;petsccall(ierr)
    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)!;petsccall(ierr)
    call KSPSetOperators(ksp,A,A,ierr)!;petsccall(ierr)
    userctx%x = x
    userctx%b = b
    userctx%A = A
    userctx%ksp = ksp
end subroutine petsc_fld_initialize

subroutine diagonal_and_off_diagonal_number(d_nnz,o_nnz)
    type(blockdef), pointer :: blk
    PetscInt, allocatable :: d_nnz(:),o_nnz(:)
    integer :: i,j,iblk,iblk0,ii_E,ii_T
    blk=>llist_head
    iblk0=blk%blk_id
    do iblk=1,np_nblk(rank+1)
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                ii_T=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2+1
                d_nnz(ii_T)=2;o_nnz(ii_T)=0
            end do
        end do
        do j=2,blk_size_ny-1
            do i=2,blk_size_nx-1
                ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
                d_nnz(ii_E)=6;o_nnz(ii_E)=0
            end do
        end do
        j=1
        if (blk%nb_level(2)==1) then
            do i=2,blk_size_nx-1
                ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
                d_nnz(ii_E)=7;o_nnz(ii_E)=2
            end do
        else
            do i=2,blk_size_nx-1
                ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
                d_nnz(ii_E)=6;o_nnz(ii_E)=2
            end do
        end if
        i=1
        ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
        d_nnz(ii_E)=8;o_nnz(ii_E)=4
        i=blk_size_nx
        ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
        d_nnz(ii_E)=8;o_nnz(ii_E)=4
        j=blk_size_ny
        if (blk%nb_level(4)==1) then
            do i=2,blk_size_nx-1
                ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
                d_nnz(ii_E)=7;o_nnz(ii_E)=2
            end do
        else
            do i=2,blk_size_nx-1
                ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
                d_nnz(ii_E)=6;o_nnz(ii_E)=2
            end do
        end if
        i=1
        ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
        d_nnz(ii_E)=8;o_nnz(ii_E)=4
        i=blk_size_nx
        ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
        d_nnz(ii_E)=8;o_nnz(ii_E)=4
        i=1
        if (blk%nb_level(3)==1) then
            do j=2,blk_size_ny-1
                ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
                d_nnz(ii_E)=7;o_nnz(ii_E)=2
            end do
        else
            do j=2,blk_size_ny-1
                ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
                d_nnz(ii_E)=6;o_nnz(ii_E)=2
            end do
        end if
        i=blk_size_nx
        if (blk%nb_level(1)==1) then
            do j=2,blk_size_ny-1
                ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
                d_nnz(ii_E)=7;o_nnz(ii_E)=2
            end do
        else
            do j=2,blk_size_ny-1
                ii_E=((blk%blk_id-iblk0)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
                d_nnz(ii_E)=6;o_nnz(ii_E)=2
            end do
        end if
        blk=>blk%next
    end do
end subroutine diagonal_and_off_diagonal_number

subroutine petsc_fld_solve(userctx,ierr)
    type(blockdef), pointer :: blk
    PetscErrorCode ierr
    type(axbsystem) userctx
    PetscScalar, pointer :: g(:)
    PetscInt, allocatable, dimension(:) :: mat_blks
    KSP, allocatable, dimension(:) :: subksp
    real(8) :: dt,dx,v,petsc_rtol_temp,dtemp,ddt,q,ddt0,vsump,vsumn,vreb,vcheck,vsump_mpi,vsumn_mpi
    PC   pc,subpc
    KSP  ksp
    Vec  b,x
    Mat  A
    PetscViewer viewer
    integer :: Ntot,i,j,k,II,JJ,nlocal,first,nblocks,nsubcycle,blk_id,m,n,iter,iblk,error
    x    = userctx%x
    b    = userctx%b
    A    = userctx%A
    ksp  = userctx%ksp
    call communicate_hydro()
    call communicate_fld()
    call applyboundconds()
    call applyboundconds_rad()
    dt=time_sys%dt_radhydro
    call blk_traversal(get_cv_block)
    call blk_traversal(calculate_planck_rosseland_opacity_block)
    call blk_traversal(calculate_fld_mfp_sigma_block)
    call calculate_fld_conductivity()
    if (lpradgradv) call blk_traversal(calculate_fld_pradgradv)
    if (lirradiation) then
        call blk_traversal(fld_irradiation_1d)
    end if
    if (lradhydro_boost) call fld_radhydro_boost_heating()
    call blk_traversal(rad_source)
    call blk_traversal(petsc_frad_init)
    q=petsc_qratio
    if (petsc_iter>1) then
        ddt0=dt*(1d0-q)/(1d0-q**petsc_iter)
    else
        ddt0=dt
    end if
    do i=1,petsc_iter
        ddt=ddt0*q**(i-1)
#if         ieos==2
        if (i>1) call blk_traversal(get_cv_block)
#endif
        call petsc_assemble_mat_vec(ddt,A,b)
        call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)!;petsccall(ierr)
        call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)!;petsccall(ierr)
        call VecAssemblyBegin(b,ierr)!;petsccall(ierr)
        call VecAssemblyEnd(b,ierr)!;petsccall(ierr)
        call petsc_set_ksp(ksp,petsc_rtol)

        call kspsolve(ksp,b,x,ierr)!;petsccall(ierr)
        call petsc_solve_convergence(b,x,ksp)
        call petsc_accept_result(x,ddt)
    end do
    if (lfld_mom) then
        call communicate_fld2()
        call blk_traversal(petsc_fld_momentum)
    end if
    call blk_traversal(convert_rho_temp_to_u_w_block)
    if (lradhydro_boost) then
        call communicate_hydro()
        call communicate_fld()
    end if
end subroutine petsc_fld_solve

subroutine petsc_assemble_mat_vec(dt,A,b)
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: dt
    integer :: i
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        if (nd==1) call petsc_Ab_block_1d(blk,dt,A,b)
        blk=>blk%next
    end do
    nullify(blk)
end subroutine petsc_assemble_mat_vec

subroutine save_Ab(A,b)
    Vec  b
    Mat  A
    PetscViewer viewer
    integer :: ierr
    character(len=128) :: alert
    call PetscViewerASCIIOpen(MPI_COMM_WORLD,'mat.xml',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_XML,ierr)
    call PetscViewerDestroy(viewer,ierr)
    call PetscViewerASCIIOpen(MPI_COMM_WORLD,'vec.xml',viewer,ierr)
    call Vecview(b,viewer,ierr)
    call PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_XML,ierr)
    call PetscViewerDestroy(viewer,ierr)
    alert='save matrix and vector'
    call abort_guangqi(alert)
end subroutine save_Ab

subroutine petsc_solve_convergence(b,x,ksp)
    integer :: ierr
    real(8) :: petsc_rtol_temp
    KSP  ksp
    Vec  b,x
    call petsc_fld_check_positivity(x)
    petsc_rtol_temp=petsc_rtol
    do while (processor%petsc_decrease_rtol)
        if (petsc_rtol_temp<1d-12) then
            print *,'petsc does not converge',petsc_rtol_temp,time_sys%ntimestep
            print *,'timestep=',time_sys%ntimestep
            call output_blocks_parallel()
            stop
        end if
        petsc_rtol_temp=petsc_rtol_temp/10
        call petsc_set_ksp(ksp,petsc_rtol_temp)
        call KSPSolve(ksp,b,x,ierr)!;petsccall(ierr)
        call petsc_fld_check_positivity(x)
    end do
end subroutine petsc_solve_convergence

subroutine petsc_fld_check_positivity(x,A)
    type(blockdef), pointer :: blk
    Vec x
    Mat, optional :: A
    PetscScalar, pointer :: g(:)
    integer :: i,j,k,ierr,blk_id,II
    processor%processor_logical=.false.
    processor%global_logical=.false.
    processor%processor_nan=.false.
    call VecGetArrayReadF90(x,g,ierr)
    blk=>llist_head
    outer2:  do k=1,np_nblk(rank+1)
        blk_id=blk%blk_id-llist_head%blk_id+1
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                II=2*((blk_id-1)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)+1
                if (g(II)<=0.or.g(II+1)<=0) then
                    processor%processor_logical=.false.
                    if (positivity_diagnosis) then
                        print *,time_sys%ntimestep,rank,i,j
                        print *,'positivity',II,g(II),g(II+1),blk%temp(i,j,1),blk%erad(i,j,1)
                    end if
                    exit outer2
                else if (isnan(g(II)).or.g(II)>huge(g(II)).or.isnan(g(II+1)).or.g(II+1)>huge(g(II+1))) then
                    processor%processor_logical=.false.
                    processor%processor_nan=.true.
                    if (rank==0) print *,'Nan or Infinity in petsc',rank,i,j,k,g(II),g(II+1),time_sys%ntimestep
                    exit outer2
                else
                    processor%processor_logical=.true.
                end if
            end do
        end do
        blk=>blk%next
    end do outer2
    call VecRestoreArrayReadF90(x,g,ierr)
    call mpi_gather(processor%processor_logical,1,MPI_LOGICAL,processor%global_logical,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (rank==0) processor%petsc_decrease_rtol=.not.(all(processor%global_logical))
    call mpi_bcast(processor%petsc_decrease_rtol,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
end subroutine petsc_fld_check_positivity

subroutine petsc_accept_result(x,ddt)
    type(blockdef), pointer :: blk
    Vec x
    PetscScalar, pointer :: g(:)
    integer :: i,j,k,ierr,blk_id,II,key(3)
    real(8) :: maxtemp_pre,maxtemp_after,r1,r2,q,dx,dt,ddt
    dt=time_sys%dt_radhydro
    call VecGetArrayReadF90(x,g,ierr)
    blk=>llist_head
    do k=1,np_nblk(rank+1)
        blk_id=blk%blk_id-llist_head%blk_id+1
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                II=2*((blk_id-1)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)+1
                blk%Erad(i,j,1)=g(II)
                blk%Erad_int(i,j,1)=blk%Erad_int(i,j,1)+g(II)*ddt/dt
                if (blk%w(1,i,j,1)<rho_thresh_petsc1) then
                    blk%temp(i,j,1)=(g(II)/a_rad)**0.25d0
                else
                    blk%temp(i,j,1)=g(II+1)
                end if
            end do
        end do
        blk=>blk%next
    end do
    nullify(blk)
    call VecRestoreArrayReadF90(x,g,ierr)
end subroutine petsc_accept_result

subroutine petsc_set_ksp(ksp,rtol)
    KSP ksp
    KSP, allocatable, dimension(:) :: subksp
    PC pc,subpc
    real(8) :: rtol
    integer :: ierr,i,nlocal,first
    if (np==1) then
        call KSPSetType(ksp,KSPGMRES,ierr)!;petsccall(ierr)
        call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)!;petsccall(ierr)
        call KSPGetPC(ksp,pc,ierr)!;petsccall(ierr)
        call PCSetType(pc,PCLU,ierr)!;petsccall(ierr)
        !call PCFactorSetLevels(pc,10,ierr);CHKERRA(ierr)
    else
        call KSPSetType(ksp,KSPGMRES,ierr)!;petsccall(ierr)
        !call KSPSetType(ksp,KSPBCGS,ierr)!;petsccall(ierr)
        call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)!;petsccall(ierr)
        call KSPGetPC(ksp,pc,ierr)!;petsccall(ierr)
        call PCSetType(pc,PCGASM,ierr)
        call PCGASMSetType(pc,PC_ASM_RESTRICT,ierr)
        !call PCGASMSetType(pc,PC_ASM_INTERPOLATE,ierr)
        call PCGASMSetOverlap(pc,1,ierr)
        call KSPsetup(ksp,ierr)
        nlocal=1;allocate(subksp(nlocal))
        call PCGASMGetSubKSP(pc,nlocal,first,subksp,ierr)
        do i=1,nlocal
            call KSPGetPC(subksp(i),subpc,ierr)
            call PCSetType(subpc,PCLU,ierr); CHKERRA(ierr)
            !call PCSetType(subpc,PCILU,ierr); CHKERRA(ierr)
            !call PCFactorSetLevels(subpc,5,ierr);CHKERRA(ierr)
            call KSPSetType(subksp(i),KSPGMRES,ierr); CHKERRA(ierr)
            !call KSPSetType(subksp(i),KSPBCGS,ierr); CHKERRA(ierr)
            call KSPSetTolerances(subksp(i),rtol/np,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)!;petsccall(ierr)
        end do
        deallocate(subksp)
    end if
end subroutine petsc_set_ksp

subroutine petsc_frad_init(blk)
    type(blockdef), pointer :: blk
    blk%Fradx=0d0
    blk%Erad_int=0d0
    blk%aradx=0d0
    blk%fld_aradwork=0d0
    if (lpradgradv) then
        !blk%fld_pradgradv=0d0
        blk%fld_pradgradv_dt=0d0
    end if
end subroutine petsc_frad_init

subroutine petsc_fld_momentum(blk)
    !to calculate the aradx and arady
    type(blockdef), pointer :: blk
    call erad_int_on_boundaries(blk)
    call frad_for_momentum(blk)
    if (larad) call fld_aradx_arady(blk)
end subroutine petsc_fld_momentum

subroutine petsc_frad_consistency()
    type(blockdef), pointer :: blk,next
    real(8) :: frad1,frad2,kx1,kx2,Erad1,Erad2,sigma1,sigma2
    blk=>llist_head
    print *,'head',rank,blk%Fradx(0,1,1)
    blk=>llist_tail
    print *,'tail',rank,blk%Fradx(blk_size_nx,1,1)
end subroutine petsc_frad_consistency

subroutine petsc_fld_recycle(userctx,ierr)
    PetscErrorCode ierr
    type(axbsystem) userctx
    call VecDestroy(userctx%x,ierr)!;petsccall(ierr)
    call VecDestroy(userctx%b,ierr)!;petsccall(ierr)
    call MatDestroy(userctx%A,ierr)!;petsccall(ierr)
    call KSPDestroy(userctx%ksp,ierr)!;petsccall(ierr)
end subroutine petsc_fld_recycle

subroutine petsc_fld_finalize(userctx,ierr)
    PetscErrorCode ierr
    type(axbsystem) userctx
    if (refine_type=='static'.or.refine_type=='none') then
        call VecDestroy(userctx%x,ierr)!;petsccall(ierr)
        call VecDestroy(userctx%b,ierr)!;petsccall(ierr)
        call MatDestroy(userctx%A,ierr)!;petsccall(ierr)
        call KSPDestroy(userctx%ksp,ierr)!;petsccall(ierr)
    else if (refine_type=='adaptive'.or.refine_type=='mixed') then
    end if
end subroutine petsc_fld_finalize

end module petsc_fld

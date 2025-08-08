module petsc_fld_1d
#include <petsc/finclude/petscksp.h>
use datastructure
use eos
use communication
use phylib
use radiation_common_functions
use petscksp
implicit none

contains

subroutine petsc_Ab_block_1d(blk,dt,A,b)
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: dt
    integer :: i,iglb,level,itop
    level=blk%level
    itop=nx_blks*2**level*blk_size_nx
    do i=1,blk_size_nx
        iglb=(blk%key(1)-1)*blk_size_nx+i
        if (iglb==1) then
            call petsc_Ab_1d_lb(dt,A,b)
        else if (iglb==itop) then
            call petsc_Ab_1d_rb(dt,A,b)
        else
            call petsc_Ab_1d_interior(blk,i,dt,A,b)
        end if
    end do
end subroutine petsc_Ab_block_1d

subroutine petsc_Ab_1d_interior(blk,i,dt,A,b)
    !assemble the ith cell's corresponding matrix
    type(blockdef), pointer :: blk
    real(8) :: dx,dt,beta(2),alpha_xl,alpha_xu,gxl,gxu,D_xl,D_xu,cv,er,temp,v,xl,xc,xr,sxl,sxr,vol,rho
    integer :: i,iblk,ii_E,ii_T,ierr,jj
    Mat A
    Vec b
    iblk=blk%blk_id
    dx=blk%dxyz(1)
    rho=blk%w(1,i,1,1)
    call assign_beta_temp_D_1d(blk,i,dt,beta,temp,D_xl,D_xu)
    cv=blk%cv(i,1,1)
    er=blk%Erad(i,1,1)
    temp=blk%temp(i,1,1)
    call assign_g_alpha_1d(blk,dt,i,gxl,gxu,alpha_xl,alpha_xu)
    ii_E=2*((iblk-1)*blk_size_nx+i-1)
    ii_T=ii_E+1
    jj=ii_E-2
    v=-alpha_xl*gxl*D_xl*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    jj=ii_E
    v=1+(alpha_xl*gxl*D_xl+alpha_xu*gxu*D_xu+beta(2))*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    jj=ii_E+2
    v=-alpha_xu*gxu*D_xu*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    v=-4*beta(1)*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    if (rho<rho_thresh_petsc1) then
        v=cv_thresh
    else
        v=4*beta(1)*a_rad*temp**3+cv
    end if
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    v=-beta(2)
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    if (lrad_adv) then
        v=er-(3*beta(1)*a_rad*temp**4)*const_Erad+blk%erad_adv(i,1,1)*dt
    else
        v=er-(3*beta(1)*a_rad*temp**4)*const_Erad
    end if
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    if (rho<rho_thresh_petsc1) then
        v=cv_thresh*temp
    else
        v=cv*temp+3*beta(1)*a_rad*temp**4+blk%erad_source(i,1,1)*dt
    end if
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr)!;petsccall(ierr)
end subroutine petsc_Ab_1d_interior

subroutine petsc_Ab_1d_lb(dt,A,b)
    type(blockdef), pointer :: blk
    real(8) :: dt,beta(2),alpha_xl,alpha_xu,gxl,gxu,D_xl,D_xu,cv,er,temp,v,dx,xl,xc,xr,sxl,sxr,vol,rho
    integer :: i,ii_E,ii_T,jj,ierr
    Mat A
    Vec b
    blk=>blk_head
    ii_E=0
    ii_T=1
    rho=blk%w(1,1,1,1)
    i=1
    call assign_beta_temp_D_1d(blk,i,dt,beta,temp,D_xl,D_xu)
    cv=blk%cv(1,1,1)
    er=blk%Erad(1,1,1)
    temp=blk%temp(1,1,1)
    dx=blk%dxyz(1)
    call assign_g_alpha_1d(blk,dt,i,gxl,gxu,alpha_xl,alpha_xu)
    jj=ii_E
    if (rad_bound_type(1)==2) then
        v=1+(alpha_xu*gxu*D_xu+beta(2))*const_Erad
    else if (rad_bound_type(1)==9) then
        v=1+(alpha_xl*gxl*D_xl+alpha_xu*gxu*D_xu+beta(2))*const_Erad
    else if (rad_bound_type(1)==8) then
        v=1+(alpha_xu*gxu*D_xu+beta(2))*const_Erad
    end if
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    jj=ii_E+2
    v=-alpha_xu*gxu*D_xu*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    v=-4*beta(1)*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    if (rho<rho_thresh_petsc1) then
        v=cv_thresh
    else
        v=4*beta(1)*a_rad*temp**3+cv
    end if
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    v=-beta(2)
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    if (rad_bound_type(1)==2) then
        if (lrad_adv) then
            v=er-(3*beta(1)*a_rad*temp**4)*const_Erad+blk%erad_adv(1,1,1)*dt
        else
            v=er-(3*beta(1)*a_rad*temp**4)*const_Erad
        end if
    else if (rad_bound_type(1)==9) then
        if (lrad_adv) then
            !v=er-(3*beta(1)*a_rad*temp**4-alpha_xl*gxl*D_xl*erad_bound(1))*const_Erad+blk%erad_adv(1,1,1)*dt
            v=er-(3*beta(1)*a_rad*temp**4-alpha_xl*gxl*D_xl*blk%Erad(0,1,1))*const_Erad+blk%erad_adv(1,1,1)*dt
        else
            !v=er-(3*beta(1)*a_rad*temp**4-alpha_xl*gxl*D_xl*erad_bound(1))*const_Erad
            v=er-(3*beta(1)*a_rad*temp**4-alpha_xl*gxl*D_xl*blk%Erad(0,1,1))*const_Erad
        end if
    else if (rad_bound_type(1)==8) then     !fixed Frad
        v=gxl*frad_bound(1)*dt+er-(3*beta(1)*a_rad*temp**4)*const_Erad
    end if
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    if (rho<rho_thresh_petsc1) then
        v=cv_thresh*temp
    else
        v=cv*temp+3*beta(1)*a_rad*temp**4+blk%erad_source(1,1,1)*dt
    end if
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    nullify(blk)
end subroutine petsc_Ab_1d_lb

subroutine petsc_Ab_1d_rb(dt,A,b)
    type(blockdef), pointer :: blk
    real(8) :: dt,beta(2),alpha_xl,alpha_xu,gxl,gxu,D_xl,D_xu,cv,er,temp,v,dx,q,xl,xc,xr,sxl,sxr,vol,rho
    integer :: i,ii_E,ii_T,jj,ierr
    Mat A
    Vec b
    blk=>blk_tail
    ii_E=2*((nblk_total-1)*blk_size_nx+blk_size_nx-1)
    ii_T=ii_E+1
    rho=blk%w(1,blk_size_nx,1,1)
    i=blk_size_nx
    call assign_beta_temp_D_1d(blk,i,dt,beta,temp,D_xl,D_xu)
    cv=blk%cv(blk_size_nx,1,1)
    er=blk%Erad(blk_size_nx,1,1)
    temp=blk%temp(blk_size_nx,1,1)
    dx=blk%dxyz(1)
    call assign_g_alpha_1d(blk,dt,i,gxl,gxu,alpha_xl,alpha_xu)
    if (igeometry==0) then
        alpha_xl=dt/dx
        alpha_xu=dt/dx
        gxl=1d0/dx
        gxu=1d0/dx
    else if (igeometry==2) then
        xl=blk%x_center(blk_size_nx-1)
        xc=blk%x_center(blk_size_nx)
        xr=blk%x_center(blk_size_nx+1)
        alpha_xl=dt/(xc-xl)
        alpha_xu=dt/(xr-xc)
        sxl=blk%surf1(blk_size_nx-1,1,1)
        sxr=blk%surf1(blk_size_nx,1,1)
        vol=blk%vol(blk_size_nx,1,1)
        gxl=sxl/vol
        gxu=sxr/vol
    end if
    if (igeometry==0) then
        if (rad_bound_type(2)==1) then
            q=1.5d0*blk%sigma_rosseland(blk_size_nx,1,1)*dx
            v=1+(beta(2)+alpha_xl*gxl*D_xl+q*alpha_xu*gxu*D_xu/(1d0+q))*const_Erad
        else if (rad_bound_type(2)==2) then
            v=1+(beta(2)+alpha_xl*gxl*D_xl)*const_Erad
        end if
    else if (igeometry==2) then
        if (rad_bound_type(2)==1) then
            v=1+(alpha_xl*gxl*D_xl+alpha_xu*gxu*D_xu*(xr**2-xc**2)/xr**2+beta(2))
        else if (rad_bound_type(2)==2) then
            v=1+(alpha_xl*gxl*D_xl+beta(2))*const_Erad
        end if
    end if
    jj=ii_E
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    jj=ii_E-2
    v=-alpha_xl*gxl*D_xl*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    v=-4*beta(1)*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    if (rho<rho_thresh_petsc1) then
        v=cv_thresh
    else
        v=4*beta(1)*a_rad*temp**3+cv
    end if
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    v=-beta(2)
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    if (lrad_adv) then
        v=er-(3*beta(1)*a_rad*temp**4)*const_Erad+blk%erad_adv(blk_size_nx,1,1)*dt
    else
        v=er-(3*beta(1)*a_rad*temp**4)*const_Erad
    end if
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    if (rho<rho_thresh_petsc1) then
        v=cv_thresh*temp
    else
        v=cv*temp+3*beta(1)*a_rad*temp**4+blk%erad_source(blk_size_nx,1,1)*dt
    end if
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr)!;petsccall(ierr)
    nullify(blk)
end subroutine petsc_Ab_1d_rb

end module petsc_fld_1d

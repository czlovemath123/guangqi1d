module eos
use datastructure
use mathlib
use phylib
#if ieos==1
    use constant_gamma_eos
#elif ieos==2
    use eos_analytic
#endif
implicit none

contains


subroutine initialize_eos_environment()
#if     ieos==2
    !eos_hllc_analytic
    call initialize_x_ratio()
#if     ieosmodule==3
    allocate(temp_h(nrho_h_he,4),temp_he(nrho_h_he,4))
#else
    allocate(temp_division(environment%nrho,4))
#endif
    call generate_rho_t_eos()
#if     ieosmodule==3
    call generate_temperature_division_hydrogen(environment%rho_eos,environment%t_eos,temp_h)
    call generate_temperature_division_helium(environment%rho_eos,environment%t_eos,temp_he)
#else
    call generate_temperature_division(environment%rho_eos,environment%t_eos,temp_division)
#endif
    if (rank==0) print *,'using analytic eos'
#endif
end subroutine initialize_eos_environment

subroutine convert_u_to_w_block(blk,guard)
    !given u, calculate w, temp and, egv
    type(blockdef), pointer :: blk
    integer :: i,j,xliml,xlimu,yliml,ylimu
    real(8) :: u(5),w(5),temp,egv,p,rho
    character(len=128) :: alert
    logical, optional :: guard
    if (present(guard)) then
        xliml=blk_xlb;xlimu=blk_xub
        yliml=1;ylimu=1
    else
        xliml=1;xlimu=blk_size_nx
        yliml=1;ylimu=1
    end if
    do j=yliml,ylimu
        do i=xliml,xlimu
            u=blk%u(1:5,i,j,1)
            call eos_utow(u,w,temp,egv)
            blk%w(1:5,i,j,1)=w
            blk%temp(i,j,1)=temp
            blk%egv(i,j,1)=egv
        end do
    end do
end subroutine convert_u_to_w_block

subroutine convert_w_to_u_block(blk,guard)
    type(blockdef), pointer :: blk
    integer :: i,j,xliml,xlimu,yliml,ylimu
    real(8) :: u(5),w(5),temp,egv
    logical, optional :: guard
    if (present(guard)) then
        xliml=blk_xlb;xlimu=blk_xub
        yliml=1;ylimu=1
    else
        xliml=1;xlimu=blk_size_nx
        yliml=1;ylimu=1
    end if
    do j=blk_ylb,blk_yub
        do i=blk_xlb,blk_xub
            w=blk%w(1:5,i,j,1)
            call eos_wtou(w,u,temp,egv)
            blk%u(1:5,i,j,1)=u
            blk%temp(i,j,1)=temp
            blk%egv(i,j,1)=egv
        end do
    end do
end subroutine convert_w_to_u_block

subroutine calculate_entropy_block(blk)
    !used only in output step
    !erg/g
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: rho,temp,ntot,species4(4)
    real(8), allocatable :: species5(:)
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            rho=blk%w(1,i,j,1)
            temp=blk%temp(i,j,1)
            blk%entropy(i,j,1)=srhot(rho,temp)
        end do
    end do
end subroutine calculate_entropy_block

subroutine get_temp_from_egv_rho(egv,w,temp)
    real(8), allocatable, dimension(:,:,:) :: egv,temp
    real(8), allocatable, dimension(:,:,:,:) :: w
end subroutine get_temp_from_egv_rho

subroutine convert_rho_temp_to_u_w_block(blk)
    type(blockdef), pointer :: blk
    real(8), allocatable, dimension(:,:,:) :: egv_old
    real(8) :: rho,temp
    integer :: i,j
    call allocate_cell_data_block(egv_old)
    egv_old=blk%egv
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            rho=blk%w(1,i,j,1)
            temp=blk%temp(i,j,1)
            blk%egv(i,j,1)=egvrhot(rho,temp)
            blk%u(5,i,j,1)=blk%u(5,i,j,1)-egv_old(i,j,1)+blk%egv(i,j,1)
            blk%w(5,i,j,1)=prhot(rho,temp)
        end do
    end do
    deallocate(egv_old)
end subroutine convert_rho_temp_to_u_w_block

subroutine get_cv_block(blk)
    type(blockdef), pointer :: blk
    real(8) :: rho,temp
    integer :: i,j
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            rho=blk%w(1,i,j,1)
            temp=blk%temp(i,j,1)
            blk%cv(i,j,1)=eos_cv(rho,temp)
        end do
    end do
end subroutine get_cv_block

subroutine assign_u_to_w_cell(blk,u,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: w(5),u(5),temp,egv
    call eos_utow(u,w,temp,egv)
    blk%w(1:5,i,j,1)=w
    blk%u(1:5,i,j,1)=u
    blk%temp(i,j,1)=temp
    blk%egv(i,j,1)=egv
end subroutine assign_u_to_w_cell

subroutine assign_w_to_u_cell(blk,w,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: w(5),u(5),temp,egv
    call eos_wtou(w,u,temp,egv)
    blk%w(1:5,i,j,1)=w
    blk%u(1:5,i,j,1)=u
    blk%temp(i,j,1)=temp
    blk%egv(i,j,1)=egv
end subroutine assign_w_to_u_cell

subroutine examine_temperature(temp)
    real(8), dimension(:,:,:), allocatable :: temp
    integer :: i,j,k
    do i=xlb,xub
        if (temp(i,1,1)<=0d0) then
            print *,'negative temperature found at index', i
            stop
        end if
    end do
end subroutine examine_temperature

function flooregv(rho,x)
    real(8) :: flooregv,rho
    real(8), optional :: x
    if (present(x)) then
        flooregv=egvrhot(rho,inner_floor_temp/x)
    else
        flooregv=egvrhot(rho,floor_temp)
    end if
end function flooregv

function pu(u)
    real(8) :: u(5),pu,ek,egv,temp,rho
    ek=0.5d0*u(1)*((u(2)/u(1))**2+(u(3)/u(1))**2+(u(4)/u(1))**2)
#if     ieos==1
    pu=(gamma_gas-one)*(u(5)-ek)
#elif   ieos==2
    egv=u(5)-ek
    rho=u(1)
    temp=solvetegv(egv,rho)
    pu=prhot(rho,temp)
#endif
end function pu


#if     ieos==3
!subroutine cal_species_h_he(info)
!    !calculate all species (H2,HI,HII,HeI,HeII,HeIII,e-)
!    type(infodef) :: info
!    real(8) :: x,y,mu,species(7),mu1,mu2,mu3,mu4,mu5
!    integer :: i,j,k,chunk
!    x=info%h_ratio
!    y=info%he_ratio
!    mu1=1d0/(x/muh2+y/muhe)
!    mu2=1d0/(x/muh+y/muhe)
!    mu3=1d0/(2d0*x/muh+y/muhe)
!    mu4=1d0/(2d0*x/muh+2d0*y/muhe)
!    mu5=1d0/(2d0*x/muh+3d0*y/muhe)
!    do i=1,nx
!        do j=1,ny
!            do k=1,nz
!                info%mu(i,j,k)=eos_tabulated_maw_2d(info%w(i,j,k,1),info%temp(i,j,k))
!                mu=info%mu(i,j,k)
!                if (mu.le.mu1*1.0001.and.mu.gt.mu1) then
!                    call h2_hI_heI(x,y,mu1,species)
!                else if (mu.le.mu1.and.mu.gt.mu2) then
!                    call h2_hI_heI(x,y,mu,species)
!                else if (mu.le.mu2.and.mu.gt.mu3) then
!                    call hI_hII_heI_e(x,y,mu,species)
!                else if (mu.le.mu3.and.mu.gt.mu4) then
!                    call hII_heI_heII_e(x,y,mu,species)
!                else if (mu.le.mu4.and.mu.gt.mu5) then
!                    call hII_heII_heIII_e(x,y,mu,species)
!                else if (mu.le.mu5.and.mu.gt.mu5/1.0001) then
!                    call hII_heII_heIII_e(x,y,mu5,species)
!                else
!                    print *, 'mu out of bound',mu,mu1,mu2,mu3,mu4,mu5,x,y
!                    stop
!                end if
!                info%species(i,j,k,:)=species
!            end do
!        end do
!    end do
!end subroutine cal_species_h_he
#endif

end module eos

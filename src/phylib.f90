module phylib
use datastructure
use mathlib
implicit none

!**************interp method, EOS Riemann solver related tables and quantities***********
!type(table2d), public :: cs_trho,eg_trho,t_prho,p_trho,gamma_trho,eg_rhop,rho_egp,maw_trho,t_egrho
!type(table1d), public :: gamma_t,egm_t,g_t,t_g,cs_t,maw_t,g_egm,egm_g,t_egm,cs_g
!real(8), public :: rho_min_bound,rho_max_bound,p_min_bound,p_max_bound,cs_min_bound,cs_max_bound,  &
!    t_min_bound,t_max_bound,gamma_max,gamma_min,g_min_bound,g_max_bound
!**************interp method, EOS Riemann solver related tables and quantities***********

real(8) :: dust_beta,opacity_gas,opacity_dust,dens_threshold_opacity,dens_cooling_threshold,  &
    dust_formation_length_scale,I_dust_threshold,r_dust_formation,dust_radcool(4),dust_gas_ratio
integer :: idustform
public :: dust_beta,opacity_gas,opacity_dust,dens_threshold_opacity,dens_cooling_threshold,  &
    dust_formation_length_scale,I_dust_threshold,r_dust_formation,dust_radcool,idustform,dust_gas_ratio

!**************radiation transfer related tables and quantities**************************
type(table1d), public :: NK_h2o_l0
type(table2d), public :: NK_h2o_llte,NK_h2o_n0_5,NK_h2o_alpha
type(table2d), public :: kappa_lowT,kappa_highT
type(table2d), public :: rosseland_gas_opacity_table,planck_gas_opacity_table
type(table2d), public :: rosseland_dust_opacity_table,planck_dust_opacity_table
real(8), public :: critical_dtao
real(8), public :: t_min_kappa,t_max_kappa,r_min_kappa,r_max_kappa  !t=log10, r=log10(rho)-3log10(t)+18
!**************radiation transfer related tables and quantities**************************

contains

!*************boundary specification. the quantities needed for specification include****
!*********w(1:5), temp, egv, and u(1:5) at the boundary for all ieos*********************
!*********object oriented subroutines. all boundaries are objects************************
    
subroutine apply_hydro_condition(blk)
    !assign value to all cells including guard cells in all blocks
    type(blockdef), pointer :: blk
    !procedure(condition_hydro) :: sub
    real(8) :: t
    integer :: ijk(3),i,j,k,iblk
    do j=blk_ylb,blk_yub
        do i=blk_xlb,blk_xub
            call init_hydro(blk,i,j)
        end do
    end do
end subroutine apply_hydro_condition

subroutine apply_rad_condition(blk)
    type(blockdef), pointer :: blk
    integer :: ijk(3),i,j,k,iblk
    do j=blk_ylb,blk_yub
        do i=blk_xlb,blk_xub
            call init_rad(blk,i,j)
        end do
    end do
end subroutine apply_rad_condition

function blackbody_rad_power(temp)
    !Stefan-Boltzmann law, the intensity has a 1/pi relation
    real(8) :: blackbody_rad_power,temp
    blackbody_rad_power=sigma_sb*temp**4
end function blackbody_rad_power

function planck_function(temp)
    !rad intensity integrated over frequency and assume blackbody
    real(8) :: planck_function,temp
    planck_function=a_rad*c_light/4d0/pi*temp**4
end function planck_function

function planck_law_wavelength(temp,x)
    !x in cm
    real(8) :: temp,x,planck_law_wavelength
    planck_law_wavelength=2*h_planck*c_light**2/x**5/(exp(h_planck*c_light/x/kb/temp)-1)
end function planck_law_wavelength

function planck_function_peak_frequency(temp)
    real(8) :: planck_function_peak_frequency,temp
    planck_function_peak_frequency=5.879d10*temp
end function planck_function_peak_frequency

function planck_law_frequency_dlnnu(temp,nu)
    !nu in s^{-1}
    real(8) :: temp,nu,planck_law_frequency_dlnnu,s
    s=h_planck*nu/kb/temp
    if (s<1d-6) then
        planck_law_frequency_dlnnu=2d0*nu**3/c_light**2*kb*temp
    else if (s>20d0) then
        planck_law_frequency_dlnnu=2d0*h_planck*nu**4/c_light**2/exp(s)
    else
        planck_law_frequency_dlnnu=2d0*h_planck*nu**4/c_light**2/(exp(s)-1d0)
    end if
end function planck_law_frequency_dlnnu

function planck_law_dfrequency_dlnnu(temp,nu)
    real(8) :: temp,nu,planck_law_dfrequency_dlnnu,s
    s=h_planck*nu/kb/temp
    if (s<1d-6) then
        planck_law_dfrequency_dlnnu=2d0*nu**3/c_light**2*kb
    else if (s>20d0) then
        planck_law_dfrequency_dlnnu=2d0*h_planck*nu**4/c_light**2/exp(s)*s/temp
    else
        planck_law_dfrequency_dlnnu=2d0*h_planck*nu**4/c_light**2/(exp(s)-1d0)**2*exp(s)*s/temp
    end if
end function planck_law_dfrequency_dlnnu

function rad_energy_density(temp)
    !radiation energy density assume blackbody
    real(8) :: rad_energy_density,temp
    rad_energy_density=a_rad*temp**4
end function rad_energy_density

function sum_eg()
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: sum_eg
    sum_eg=0d0
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            sum_eg=sum_eg+blk%vol(j,1,1)*blk%egv(j,1,1)
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%next
        end if
    end do
end function sum_eg

function sum_ek()
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: sum_ek
    sum_ek=0d0
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            sum_ek=sum_ek+blk%vol(j,1,1)*(blk%u(5,j,1,1)-blk%egv(j,1,1))
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%next
        end if
    end do
end function sum_ek

function sum_er()
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: sum_er
    sum_er=0d0
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            sum_er=sum_er+blk%vol(j,1,1)*blk%Erad(j,1,1)
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%next
        end if
    end do
end function sum_er

function sum_gpotential()
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: sum_gpotential
    sum_gpotential=0d0
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            sum_gpotential=sum_gpotential+blk%vol(j,1,1)*blk%w(1,j,1,1)*blk%gpotential(j,1,1)
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%next
        end if
    end do
end function sum_gpotential

!collective functions

function mass_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: mass_sum
    mass_sum=sum(blk%w(1,1:blk_size_nx,1:blk_size_ny,1)*blk%vol(1:blk_size_nx,1:blk_size_ny,1))
    if (igeometry==2) mass_sum=mass_sum*2d0*pi
end function mass_sum

function angular_momentum_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j,ip_am
    real(8) :: angular_momentum_sum
    ip_am=ipassive(pn_am)
    angular_momentum_sum=sum(blk%passive_scalar(ip_am,1:blk_size_nx,1:blk_size_ny,1)*blk%vol(1:blk_size_nx,1:blk_size_ny,1))
    angular_momentum_sum=angular_momentum_sum*2d0*pi
end function angular_momentum_sum

function Egas_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: Egas_sum
    Egas_sum=sum(blk%u(5,1:blk_size_nx,1:blk_size_ny,1)*blk%vol(1:blk_size_nx,1:blk_size_ny,1))
    if (igeometry==2) Egas_sum=Egas_sum*2d0*pi
end function Egas_sum

function E2dpolar_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: E2dpolar_sum,r,theta,s,rho,omega
    E2dpolar_sum=sum(blk%u(5,1:blk_size_nx,1:blk_size_ny,1)*blk%vol(1:blk_size_nx,1:blk_size_ny,1))
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            rho=blk%w(1,i,j,1)
            omega=blk%omega(i,j,1)
            r=blk%x_center(i)
            theta=blk%x_center(j)
            s=r*sin(theta)
            E2dpolar_sum=E2dpolar_sum+0.5d0*rho*(omega*s)**2
        end do
    end do
    if (igeometry==2) E2dpolar_sum=E2dpolar_sum*2d0*pi
end function E2dpolar_sum

function Erad_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: Erad_sum
    Erad_sum=sum(blk%Erad(1:blk_size_nx,1:blk_size_ny,1)*blk%vol(1:blk_size_nx,1:blk_size_ny,1))
    if (igeometry==2) Erad_sum=Erad_sum*2d0*pi
end function Erad_sum

function heat1_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: heat1_sum
    heat1_sum=sum(blk%vol(1:blk_size_nx,1:blk_size_ny,1)*blk%heat1(1:blk_size_nx,1:blk_size_ny,1))
    if (igeometry==2) heat1_sum=heat1_sum*2d0*pi
end function heat1_sum

function grav_potential_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: grav_potential_sum,m,r
    real(8), allocatable :: phi(:,:,:)
    call allocate_cell_data_block(phi)
    m=central_star%core%mass
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            r=blk%x_center(i)
            phi(i,j,1)=-m*gr/r
        end do
    end do
    grav_potential_sum=sum(blk%vol(1:blk_size_nx,1:blk_size_ny,1)*blk%u(1,1:blk_size_nx,1:blk_size_ny,1)*phi(1:blk_size_nx,1:blk_size_ny,1))
    if (igeometry==2) grav_potential_sum=grav_potential_sum*2d0*pi
    deallocate(phi)
end function grav_potential_sum

!profiling functions

function massflux_x(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: massflux_x
    massflux_x=blk%xflux(1,i,j,1)
end function massflux_x

function fr_x(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: fr_x
    fr_x=blk%Fradx(i,j,1)
end function fr_x

function amflux_x(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j,ip
    real(8) :: amflux_x
    ip=ipassive(pn_am)
    amflux_x=blk%xpflux(ip,i,j,1)
end function amflux_x

function torque_x(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: torque_x
    torque_x=blk%torque_xz(i,j,1)
end function torque_x

function totalEgflux_x(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j,ip
    real(8) :: flux,l,r,theta,vphi,ekphi_flux,ephi_flux,totalEgflux_x,m
    ip=ipassive(pn_am)
    flux=blk%xflux(5,i,j,1)
    l=blk%xpflux(ip,i,j,1)/blk%xflux(1,i,j,1)
    r=blk%x_interface(i)
    theta=blk%y_center(j)
    vphi=l/r/sin(theta)
    ekphi_flux=0.5*blk%xflux(1,i,j,1)*vphi**2
    m=central_star%core%mass
    totalEgflux_x=flux+ekphi_flux
end function totalEgflux_x

function erad_adv_flux_x(blk,i,j)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: erad_adv_flux_x
    erad_adv_flux_x=blk%erad_xflux(i,j,1)
end function erad_adv_flux_x

subroutine irradiation_dtau(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: rho,kappa,dx
    do i=1,blk_size_nx
        rho=blk%w(1,i,1,1)
        kappa=blk%kappa_planck(i,1,1)
        dx=blk%x_interface(i)-blk%x_interface(i-1)
        blk%dtau(i,1,1)=rho*kappa*dx
    end do
end subroutine irradiation_dtau

subroutine irradiation_update_tau(blk)
    type(blockdef), pointer :: blk
    integer :: i
    do i=1,blk_size_nx
        blk%tau(i,1,1)=blk%tau(i-1,1,1)+blk%dtau(i,1,1)
    end do
end subroutine irradiation_update_tau

function zsimplegas(t,m)
    real(8) :: zsimplegas,t,m
    zsimplegas=(sqrt(8d0)*pi**1.5*(kb*m*t)**1.5)/h_planck**3
end function zsimplegas

end module phylib

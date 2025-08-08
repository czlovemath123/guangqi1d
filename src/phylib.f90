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

subroutine cal_species_h_he_indivi(x,y,mu,species)
    real(8) :: x,y,mu,mu1,mu2,mu3,mu4,mu5,species(7)
    mu1=1d0/(x/muh2+y/muhe)
    mu2=1d0/(x/muh+y/muhe)
    mu3=1d0/(2d0*x/muh+y/muhe)
    mu4=1d0/(2d0*x/muh+2d0*y/muhe)
    mu5=1d0/(2d0*x/muh+3d0*y/muhe)
    if (mu.le.mu1*1.0001.and.mu.gt.mu1) then
        call h2_hI_heI(x,y,mu1,species)
    else if (mu.le.mu1.and.mu.gt.mu2) then
        call h2_hI_heI(x,y,mu,species)
    else if (mu.le.mu2.and.mu.gt.mu3) then
        call hI_hII_heI_e(x,y,mu,species)
    else if (mu.le.mu3.and.mu.gt.mu4) then
        call hII_heI_heII_e(x,y,mu,species)
    else if (mu.le.mu4.and.mu.gt.mu5) then
        call hII_heII_heIII_e(x,y,mu,species)
    else if (mu.le.mu5.and.mu.gt.mu5/1.0001) then
        call hII_heII_heIII_e(x,y,mu5,species)
    else
        print *, 'mu out of bound',mu,mu1,mu2,mu3,mu4,mu5,x,y
        stop
    end if
end subroutine cal_species_h_he_indivi

subroutine h2_hI_heI(x,y,mu,species)
    !h2 disassociates to become hI
    real(8) :: x,y,a,mu,species(7)
    a=(1-mu*x/muh2-mu*y/muhe)/(mu*x*(1d0/muh-1d0/muh2))
    species(1)=(1d0-a)*x
    species(2)=a*x
    species(3)=0d0
    species(4)=y
    species(5)=0d0
    species(6)=0d0
    species(7)=0d0
end subroutine h2_hI_heI

subroutine hI_hII_heI_e(x,y,mu,species)
    !hI becomes hII and e-
    real(8) :: x,y,a,mu,species(7)
    a=(1-mu*x/muh-mu*y/muhe)/(mu*x*(2d0/muh-1d0/muh))
    species(1)=0d0
    species(2)=(1d0-a)*x
    species(3)=a*x
    species(4)=y
    species(5)=0d0
    species(6)=0d0
    species(7)=a*x
end subroutine hI_hII_heI_e

subroutine hII_heI_heII_e(x,y,mu,species)
    !heI becomes heII and e-
    real(8) :: x,y,a,mu,species(7)
    a=(1-2*mu*x/muh-mu*y/muhe)/(mu*y*(2d0/muhe-1d0/muhe))
    species(1)=0d0
    species(2)=0d0
    species(3)=x
    species(4)=(1d0-a)*y
    species(5)=a*y
    species(6)=0d0
    species(7)=x+a*y
end subroutine hII_heI_heII_e

subroutine hII_heII_heIII_e(x,y,mu,species)
    !heII becomes heIII and e-
    real(8) :: x,y,a,mu,species(7)
    a=(1d0-2d0*mu*x/muh-2d0*mu*y/muhe)/(mu*y*(3d0/muhe-2d0/muhe))
    species(1)=0d0
    species(2)=0d0
    species(3)=x
    species(4)=0d0
    species(5)=(1d0-a)*y
    species(6)=a*y
    species(7)=x+2d0*a*y+(1d0-a)*y
end subroutine hII_heII_heIII_e


!******************************************************************************
!*******************************end EOS part***********************************
!******************************************************************************

subroutine initialize_star(star,m,temp,r,period,spin)
    type(stellar) :: star
    real(8) :: m,loc(3),temp,r
    real(8), optional :: period,spin
    loc=0d0
    star%core%mass=m
    star%core%xyz=loc
    star%teff=temp
    star%radius=r
    star%luminosity=4*pi*r*r*sigma_sb*temp**4
    star%lum_atm=0d0
    if (present(period)) then
        star%period=period   !if the star has some periodic property
    else
        star%period=t_hubble
    end if
    if (present(spin)) then
        star%spin=spin
    else
        star%spin=0d0
    end if
end subroutine initialize_star

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

!subroutine set_floor_temp(u)
!    real(8) :: u(5),egv,egv_floor
!    egv=u(5)-half*(u(2)**2d0+u(3)**2d0+u(4)**2d0)/u(1)
!    egv_floor=egvrhot(u(1),floor_temp)
!    if (egv<egv_floor) u(5)=egv_floor+half*(u(2)**2d0+u(3)**2d0+u(4)**2d0)/u(1)
!end subroutine set_floor_temp

subroutine set_ceiling_vx(u)
    real(8) :: u(5),vx,egv
    vx=u(2)/u(1)
    if (abs(vx)>ceiling_vx) then
        egv=u(5)-half*(u(2)**2d0+u(3)**2d0+u(4)**2d0)/u(1)
        vx=vx/abs(vx)*ceiling_vx
        u(2)=u(1)*vx
        u(5)=half*(u(2)**2+u(3)**2+u(4)**2d0)/u(1)+egv
    end if
end subroutine set_ceiling_vx

subroutine set_ceiling_vy(u)
    real(8) :: u(5),vy,egv
    vy=u(2)/u(1)
    if (abs(vy)>ceiling_vy) then
        egv=u(5)-half*(u(2)**2d0+u(3)**2d0+u(4)**2d0)/u(1)
        vy=vy/abs(vy)*ceiling_vy
        u(2)=u(1)*vy
        u(5)=half*(u(2)**2+u(3)**2+u(4)**2d0)/u(1)+egv
    end if
end subroutine set_ceiling_vy

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

subroutine spherical2d_vphi(blk,vphi)
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: vphi
end subroutine spherical2d_vphi

subroutine extrapolate1d(x,v,direction)
    real(8) :: x(3),v(3)
    character(len=1) :: direction
    if (direction=='+') then
        v(3)=v(2)+(v(2)-v(1))/(x(2)-x(1))*(x(3)-x(2))
    else if (direction=='-') then
        v(1)=v(2)+(v(3)-v(2))/(x(3)-x(2))*(x(1)-x(2))
    end if
end subroutine extrapolate1d

subroutine angular_momentum(blk)
    type(blockdef), pointer :: blk
    real(8) :: r,rho,omega
    integer :: i,j,ip_am
    if (lam_con) then
        ip_am=ipassive(pn_am)
        do j=blk_ylb,blk_yub
            do i=blk_xlb,blk_xub
                rho=blk%w(1,i,j,1)
                omega=blk%omega(i,j,1)
                blk%passive_scalar(ip_am,i,j,1)=rho*blk%l_omega(i,j,1)*omega
            end do
        end do
    end if
end subroutine angular_momentum

subroutine angular_frequency(blk)
    type(blockdef), pointer :: blk
    real(8) :: r,rho,omega,l
    integer :: i,j,ip_am
    if (lam_con) then
        ip_am=ipassive(pn_am)
        blk%omega(blk_xlb:blk_xub,blk_ylb:blk_yub,1)=blk%passive_scalar(ip_am,blk_xlb:blk_xub,blk_ylb:blk_yub,1)/blk%u(1,blk_xlb:blk_xub,blk_ylb:blk_yub,1)/    &
            blk%l_omega(blk_xlb:blk_xub,blk_ylb:blk_yub,1)
        blk%omega1=blk%omega
    end if
end subroutine angular_frequency

!collective functions

function mass_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: mass_sum
    mass_sum=sum(blk%w(1,1:blk_size_nx,1:blk_size_ny,1)*blk%vol(1:blk_size_nx,1:blk_size_ny,1))
    if (igeometry==2) mass_sum=mass_sum*2d0*pi
end function mass_sum

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

function vis_heat_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: vis_heat_sum
    vis_heat_sum=sum(blk%vol(1:blk_size_nx,1:blk_size_ny,1)*blk%viscous_heat(1:blk_size_nx,1:blk_size_ny,1))
    if (igeometry==2) vis_heat_sum=vis_heat_sum*2d0*pi
end function vis_heat_sum

function vis_heat2_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: vis_heat2_sum
    vis_heat2_sum=sum(blk%vol(1:blk_size_nx,1:blk_size_ny,1)*blk%heat2(1:blk_size_nx,1:blk_size_ny,1))
    if (igeometry==2) vis_heat2_sum=vis_heat2_sum*2d0*pi
end function vis_heat2_sum

function heat1_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: heat1_sum
    heat1_sum=sum(blk%vol(1:blk_size_nx,1:blk_size_ny,1)*blk%heat1(1:blk_size_nx,1:blk_size_ny,1))
    if (igeometry==2) heat1_sum=heat1_sum*2d0*pi
end function heat1_sum

function ekphi_mix_heat_sum(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: ekphi_mix_heat_sum
    ekphi_mix_heat_sum=sum(blk%vol(1:blk_size_nx,1:blk_size_ny,1)*blk%ekphi_mix_heat(1:blk_size_nx,1:blk_size_ny,1))
    if (igeometry==2) ekphi_mix_heat_sum=ekphi_mix_heat_sum*2d0*pi
end function ekphi_mix_heat_sum

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

subroutine profile_const_x(x,profiler,v)
    type(blockdef), pointer :: blk
    procedure(field_calculator), pointer :: profiler
    real(8) :: x,v,proc_v,sum_v,x1,x2,v1,v2,surf1,surf2
    real(8), allocatable :: v_array(:)
    logical :: belong,l1,l2,l3
    integer :: i,j,k,error
    v=0d0
    proc_v=0d0
    allocate(v_array(blk_size_ny))
    blk=>llist_head
    do k=1,np_nblk(rank+1)
        belong=.false.
        l1=(x==n_domain(1).and.blk%key(1)==1)
        l2=(x==n_domain(2).and.blk%key(1)==nx_blks*2**blk%level)
        l3=(x<blk%x_interface(blk_size_nx).and.x>=blk%x_interface(0))
        if (l1.or.l2.or.l3) then
            if (l1) then
                i=1
            else if (l2) then
                i=blk_size_nx
            else if (l3) then
                do i=1,blk_size_nx
                    if (x<blk%x_interface(i).and.x>=blk%x_interface(i-1)) exit
                end do
            end if
            x1=blk%x_interface(i-1)
            x2=blk%x_interface(i)
            do j=1,blk_size_ny
                surf1=blk%surf1(i-1,j,1)
                surf2=blk%surf1(i,j,1)
                v1=profiler(blk,i-1,j)*surf1
                v2=profiler(blk,i,j)*surf2
                v_array(j)=(v1+(x-x1)/(x2-x1)*(v2-v1))
            end do
            belong=.true.
            if (belong) proc_v=proc_v+sum(v_array)
        end if
        blk=>blk%next
    end do
    if (igeometry==2) proc_v=proc_v*2d0*pi
    call mpi_reduce(proc_v,sum_v,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,error)
    call mpi_bcast(sum_v,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
    v=sum_v
    deallocate(v_array)
end subroutine profile_const_x

function zsimplegas(t,m)
    real(8) :: zsimplegas,t,m
    zsimplegas=(sqrt(8d0)*pi**1.5*(kb*m*t)**1.5)/h_planck**3
end function zsimplegas

end module phylib

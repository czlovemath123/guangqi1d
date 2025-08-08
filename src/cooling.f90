module cooling
use datastructure
use mathlib
use phylib
use eos
implicit none

real(8), dimension(:), allocatable, protected :: templist_log,l0list,ntildlist,templist
real(8), dimension(:), allocatable, protected :: templist_log_co,l0list_co,ntildlist_co,templist_co
real(8), dimension(:,:), allocatable, protected :: lltetable,n0_5table,alphatable
real(8), dimension(:,:), allocatable, protected :: lltetable_co,n0_5table_co,alphatable_co
real(8), dimension(:), allocatable, protected :: templist_log_h2,l0list_h2,lltelist_h2,  &
    n0_5list_h2,alphalist_h2

contains

!subroutine initialize_cooling(info)
!    type(infodef) :: info
!    if (icooling==1) then
!        call init_nkcool()
!        if (allocated(info%divv).eqv..false.) then
!            call allocate_cell_data(info%divv)
!        end if
!    end if
!end subroutine initialize_cooling
!
!subroutine cooling_cell(info,ijk,source,dsourcedt)
!    type(infodef) :: info
!    integer :: ijk(3)
!    real(8) :: rho,egv,h_ratio,divv,dedt,source(5),dsourcedt(5),dx,dedt_h,dedt_h2,tau
!    if (icooling==1) then
!        if (nd==1) then
!            !1d molecular and atomic cooling, no metal cooling
!            rho=source(1)
!            divv=info%divv(ijk(1),ijk(2),ijk(3))
!            egv=source(5)
!            h_ratio=info%h_ratio
!            dx=dxyz(1)
!            tau=1d0  !info%tau_fld(ijk(1),1,1)
!            if (source(5)<=zero) then
!                dedt=zero
!                print *,'negative internal energy'
!                stop
!            else
!                call nk_cooling_cell(rho,divv,egv,h_ratio,dx,tau,dedt,dedt_h,dedt_h2)
!            end if
!            dsourcedt=(/0d0,0d0,0d0,0d0,dedt/)
!        else
!            dsourcedt=0d0
!        end if
!    else
!        dsourcedt=0d0
!    end if
!end subroutine cooling_cell
!
!subroutine cooling_cell_dt(info,ijk,source,tau,dt_cooling)
!    type(infodef) :: info
!    integer :: ijk(3)
!    real(8) :: dt,t,t_stop,dt_gravity,dt_radiation,dt_cooling,source(5),tau
!    if (icooling==1) then
!        !molecular and atomic cooling, no metal cooling
!        call cooling_cell_dt_nk(info,ijk,source,tau,dt_cooling)
!    else
!        !not implemented yet
!        dt_cooling=t_hubble
!    end if
!end subroutine cooling_cell_dt
!
!!********************************begin NK cooling***********************************************
!!********************************based on Neufeld & Kaufman 1993********************************
!
!subroutine init_nkcool()
!    use eos_analytic
!    allocate(templist_log(6),l0list(6),ntildlist(10),lltetable(6,10),n0_5table(6,10),alphatable(6,10),templist(6))
!    open(unit=505,file='tables/cooling/NK_h2o.txt')
!    read(unit=505,fmt=*) templist,ntildlist,l0list,lltetable,n0_5table,alphatable
!    close(505)
!    allocate(templist_log_co(6),l0list_co(6),ntildlist_co(11),lltetable_co(6,11),  &
!        n0_5table_co(6,11),alphatable_co(6,11),templist_co(6))
!    open(unit=505,file='tables/cooling/NK_co.txt')
!    read(unit=505,fmt=*) templist_co,ntildlist_co,l0list_co,lltetable_co,n0_5table_co,alphatable_co
!    close(505)
!    allocate(templist_log_h2(17),l0list_h2(17),lltelist_h2(17),n0_5list_h2(17),alphalist_h2(17))
!    open(unit=505,file='tables/cooling/NK_h2.txt')
!    read(unit=505,fmt=*) templist_log_h2,l0list_h2,lltelist_h2,n0_5list_h2,alphalist_h2
!    close(505)
!    templist_log=log10(templist)
!    templist_log_co=log10(templist_co)
!    if (allocated(temp_division)) then
!        print *,'using analytic eos together with NK cooling'
!    else
!        call initialize_x_ratio()
!        allocate(temp_division(info%ntemp,4))
!        call generate_rho_t_eos()
!        call generate_temperature_division(environment%rho_eos,environment%t_eos,temp_division)
!        print *,'using NK cooling alone (ieos/=2)'
!        allocate(info%H2(0:nx+1,0:ny+1,1),info%HI(0:nx+1,0:ny+1,1),info%HII(0:nx+1,0:ny+1,1),info%electron(0:nx+1,0:ny+1,1))
!    end if
!end subroutine init_nkcool
!
!subroutine nk_cooling_cell(rho,divv,egv,h_ratio,dx,tau,dedt,dedt_h,dedt_h2)   
!    use eos_analytic
!    use constant_gamma_eos
!    !given w, divv, egv, and h_ratio calculate dedt
!    !templist_log,l0list,ntildlist,lltetable,n0_5table and
!    !alphatable need to be read in.
!    real(8), dimension(:), allocatable :: species
!    real(8) :: divv,dt,temp,ntild,alpha,nhion,nh,nelec,lambda,nh2,nh2o,a,rho,rhoh
!    real(8) :: temp_log,l0,llte,n0_5,coolingtime,r,kappagas,radbalance,q1,  &
!        alpha_co,nco,a_co,l0_co,llte_co,n0_5_co,alpha_h2,a_h2,l0_h2,llte_h2,n0_5_h2,  &
!        l_ion,l_rec,l_ex,l_bre,l_co,l_h2o,l_h2,egv,h_ratio,dedt,nhe,dx,dtau,kappa,r_kappa,log10_t,  &
!        l_kappa,dedt_h,dedt_h2,tau
!    !lambda the cooling rate
!    temp=calculate_temp_from_egv_rho(egv,rho)
!    rhoh=rho*h_ratio
!    dedt_h=0d0
!    dedt_h2=0d0
!    allocate(species(5))
!    call calspecies(rhoh,temp,species)
!    nh2=species(1)
!    nh=species(2)
!    nhion=species(3)
!    nelec=species(4)
!    nhe=species(5)
!    if (temp>3000d0.and.rho<2e-9) then
!        !collisional excitation cooling, recombination cooling, collisional
!        !ionization cooling, and Bremsstrahlung cooling
!        !Reference Cen, Renyue (1992) ApJS
!        l_ex=0d0;l_rec=0d0;l_ion=0d0;l_bre=0d0
!        l_ex=7.5d-19*nh*nelec/(1d0+sqrt(temp/1d5))*exp(-118348d0/temp)    !collisional excitation cooling
!        l_rec=8.7d-27*nelec*nhion*sqrt(temp)*(temp/1d3)**(-0.2)/(1+(temp/1d6)**0.7)
!        l_ion=1.27d-21*nh*nelec*sqrt(temp)/(1d0+sqrt(temp/1d5))*exp(-157809d0/temp)
!        l_bre=1.42d-27*nhion*nelec*1.5d0*sqrt(temp)
!        r_kappa=log10(rho)-3d0*log10(temp)+18d0
!        log10_t=log10(temp)
!        lambda=(l_ex+l_rec+l_ion+l_bre)   !dtau
!    else
!        lambda=0d0
!    end if
!    dedt_h=-lambda
!    if (temp>100d0.and.temp<1800d0.and.rho<2d-9) then
!        !molecular cooling. coolant include h2o, co, and h2
!        !Reference Neufeld & Kaufmann (1993) ApJ
!        nh2=max(nh2,rhoh/mh2*1d-4)      !set a floor nh2 density, otherwise may never cool
!        l_h2o=0d0;l_co=0d0;l_h2=0d0
!        nh2o=rhoh/mh2*1d-4    !2E-4*nh2
!        nco=rhoh/mh2*1d-4   !2E-4*nh2    !JH Lacy, R Knacke, TR Geballe & Tokunaga, A. T. ApJ 1994
!        !******************H2O and H2 cooling***************************
!        ntild=nh2o/max(1d-16,abs(divv)/1d5)  !divv is in cgs, ntild needs divv in km/s
!        ntild=log10(ntild)
!        temp_log=log10(temp)
!        call interpolation_linear(temp_log,l0,templist_log,l0list)
!        l0=10**(-l0)
!        call interpolation_linear(temp_log,ntild,llte,templist_log,ntildlist,lltetable)
!        llte=10**(-llte)
!        call interpolation_linear(temp_log,ntild,n0_5,templist_log,ntildlist,n0_5table)
!        n0_5=10**(n0_5)
!        call interpolation_linear(temp_log,ntild,alpha,templist_log,ntildlist,alphatable)
!        a=1d0/l0+nh2/llte+(1d0/l0)*((nh2/n0_5)**alpha)*(1d0-n0_5*l0/llte)
!        l_h2o=1d0/a
!        !*****************H2O and H2 cooling****************************
!        !*****************CO and H2 cooling******************************
!        if (temp<=2000d0) then
!            call interpolation_linear(temp_log,l0_co,templist_log_CO,l0list_CO)
!            l0_co=10**(-l0_co)
!            call interpolation_linear(temp_log,ntild,llte_co,templist_log_CO,ntildlist_CO,lltetable_CO)
!            llte_co=10**(-llte_co)
!            call interpolation_linear(temp_log,ntild,n0_5_co,templist_log_CO,ntildlist_CO,n0_5table_CO)
!            n0_5_co=10**(n0_5_co)
!            call interpolation_linear(temp_log,ntild,alpha_co,templist_log_CO,ntildlist_CO,alphatable_CO)
!            a_co=1d0/l0_co+nh2/llte_co+(1/l0_co)*((nh2/n0_5_co)**alpha_co)*(1d0-n0_5_co*l0_co/llte_co)
!            l_co=1d0/a_co
!        else
!            l_co=0d0
!        end if
!        !*****************CO and H2 cooling ends**************************
!        !*****************molecular hydrogen cooling**************
!        call interpolation_linear(temp_log,l0_h2,templist_log_H2,l0list_H2)
!        l0_h2=10**(-l0_h2)
!        call interpolation_linear(temp_log,llte_h2,templist_log_H2,lltelist_H2)
!        llte_h2=10**(-llte_h2)
!        call interpolation_linear(temp_log,n0_5_h2,templist_log_H2,n0_5list_H2)
!        n0_5_h2=10**(n0_5_h2)
!        call interpolation_linear(temp_log,alpha_h2,templist_log_H2,alphalist_H2)
!        a_h2=1d0/l0_h2+nh2/llte_h2+(1/l0_h2)*((nh2/n0_5_h2)**alpha_h2)*(1d0-n0_5_h2*l0_h2/llte_h2)
!        l_h2=1d0/a_h2
!        !****************molecular hydrogen cooling ends************
!        lambda=lambda+l_h2o*nh2*nh2o+l_co*nh2*nco+l_h2*nh2*nh2
!        dedt_h2=-(l_h2o*nh2*nh2o+l_co*nh2*nco+l_h2*nh2*nh2)
!    else
!        dedt_h2=0d0
!    end if
!    deallocate(species)
!    dedt=-lambda
!end subroutine nk_cooling_cell
!
!subroutine cooling_cell_dt_nk(info,ijk,source,tau,dt)
!    type(infodef) :: info
!    integer :: ijk(3)
!    real(8) :: dt,rho,divv,egv,h_ratio,dedt,source(5),dx,dedt_h,dedt_h2,tau
!    if (nd==1) then
!        rho=source(1)
!        divv=info%divv(ijk(1),ijk(2),ijk(3))
!        egv=source(5)
!        h_ratio=info%h_ratio
!        dx=dxyz(1)
!        call nk_cooling_cell(rho,divv,egv,h_ratio,dx,tau,dedt,dedt_h,dedt_h2)
!        if (dedt==0) then
!            dt=t_hubble
!        else
!            dt=-egv/dedt
!        end if
!    else
!        dt=t_hubble
!    end if
!end subroutine cooling_cell_dt_nk

!********************************end NK cooling*************************************************

end module cooling

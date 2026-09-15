module eos_h2_HI
use mathlib
use phylib
implicit none

real(8), protected :: x_ratio

contains

subroutine initialize_x_ratio() 
    !initialize the value of x_ratio in eos_analytic
    x_ratio=environment%h_ratio
    if (rank==0) then
        print *, 'X=',x_ratio,'   Y=',one-x_ratio
    end if
end subroutine initialize_x_ratio


subroutine gammarhot(rho,t,gamm) 
    ! Calculate the Gamma with given rho and t
    real(8) :: rho,t,gamm,p,ntot,srho,st,prho,pt,t1,t2,species(3),dlnzdlnt(3),d2lnzdlnt(3)
    real(8) :: dlnndlnrho(3),dlnndlnt(3)                                                    
    call divide_temperature_domain_known(rho,t1,t2)
    call calspecies(rho,t,species,t1,t2)
    call caldlnz(t,dlnzdlnt,d2lnzdlnt)
    call caldlnndlnrho(rho,t,species,dlnndlnrho,t1,t2)
    call caldlnndlnt(rho,t,species,dlnndlnt,t1,t2)
    call calprhopt(species,dlnndlnrho,dlnndlnt,rho,t,prho,pt)
    call calsrhost(rho,t,species,dlnndlnrho,dlnndlnt,dlnzdlnt,d2lnzdlnt,srho,st)
    ntot=sum(species)
    p=ntot*kb*t
    gamm=rho/p*(p/rho*prho-p/t*pt*srho/st)
end subroutine gammarhot

function adiabatic_cs(rho,t)
    ! calculate the adiabatic sound speed with given rho and t
    real(8) :: adiabatic_cs,rho,t,gamm,p,ntot,srho,st,prho,pt,t1,t2
    real(8), dimension(3) :: dlnndlnrho,dlnndlnt,species,dlnzdlnt,d2lnzdlnt
    call divide_temperature_domain_known(rho,t1,t2)
    call calspecies(rho,t,species,t1,t2)
    call caldlnz(t,dlnzdlnt,d2lnzdlnt)
    call caldlnndlnrho(rho,t,species,dlnndlnrho,t1,t2)
    call caldlnndlnt(rho,t,species,dlnndlnt,t1,t2)
    call calprhopt(species,dlnndlnrho,dlnndlnt,rho,t,prho,pt)
    call calsrhost(rho,t,species,dlnndlnrho,dlnndlnt,dlnzdlnt,d2lnzdlnt,srho,st)
    ntot=sum(species)
    p=ntot*kb*t
    adiabatic_cs=sqrt(p/rho*prho-p/t*pt*srho/st)
end function adiabatic_cs

function eos_cv(rho,t,temp1,temp2)
    ! d(specific energy)/dT
    real(8) :: eos_cv,rho,t,t1,t2,xi_h2,xi_h,xi_he
    real(8),dimension(3) :: xi,species,dlnndlnt,dndt
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    call calspecies(rho,t,species,t1,t2)
    call caldlnndlnt(rho,t,species,dlnndlnt,t1,t2)
    call caldndt(rho,t,species,dndt)
    xi_h2=1.5d0*kb*t
    xi_h=1.5d0*kb*t+0.5d0*dish
    xi_he=1.5d0*kb*t
    xi=(/xi_h2,xi_h,xi_he/)
    eos_cv=dot_product(dndt,xi)+sum(species)*1.5d0*kb
end function eos_cv

function solvetchi_H(chi_H,rho)
    !given chi_H and rho, find the corresponding temperature
    real(8) :: solvetchi_H,chi_H,rho,lowT,lowchi_H,highT,highchi_H,convergence,t,t1,t2,nhtot,nhetot,ntot
    real(8), dimension(:), allocatable :: x,root
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    conv_mode=2
    convergence=1e-6
    call divide_temperature_domain_known(rho,t1,t2)
    lowT=t1
    highT=t2
    lowchi_H=chi_Hrhot(rho,lowT,t1,t2)
    highchi_H=chi_Hrhot(rho,highT,t1,t2)
    if (chi_H<lowchi_H) then
        print *,'solve chi_H degenerate'
        stop
    else if (chi_H>highchi_H) then
        print *,'solve chi_H degenerate'
        stop
    else
        allocate(x(5),root(5))
        x=(/highT/2d0,rho,chi_H,t1,t2/) 
        ptr=>fun_Tchi_H
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solvetchi_H=root(1)
        deallocate(x,root)
    end if
end function solvetchi_H

function chi_Hrhot(rho,t,temp1,temp2)
    !given rho and t, calculate the ionization fraction
    real(8) :: chi_Hrhot,rho,t,t1,t2,species(3),nhtot,nhetot,ntot
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then                                         ! H2
        chi_Hrhot=0d0
    else if (t>=t2) then                                    ! H
        chi_Hrhot=1d0
    else
        call calspecies(rho,t,species,t1,t2)
        ntot=x_ratio*rho/mh
        chi_Hrhot=species(2)/ntot                           
    end if
end function chi_Hrhot

function fun_Tchi_H(x)
    real(8), dimension(:), allocatable :: x
    real(8) :: fun_Tchi_H,t,rho,chi_H,t1,t2
    t=x(1)
    rho=x(2)
    chi_H=x(3)
    t1=x(4)
    t2=x(5)
    fun_Tchi_H=chi_H-chi_Hrhot(rho,t,t1,t2)
end function fun_Tchi_H

function solvetp(p,rho)
    !given p and rho, find the corresponding temperature
    real(8) :: solvetp,p,rho,lowT,lowp,highT,highp,convergence,t,t1,t2,nhtot,nhetot,ntot
    real(8), dimension(:), allocatable :: x,root        
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    conv_mode=2
    convergence=1e-6
    call divide_temperature_domain_known(rho,t1,t2)
    lowT=t1
    highT=t2
    lowp=prhot(rho,lowT,t1,t2)
    highp=prhot(rho,highT,t1,t2)
    if (p<=lowp) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=nhtot/2+nhetot
        solvetp=p/kb/ntot
    else if (p>=highp) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=nhtot+nhetot
        solvetp=p/kb/ntot
    else
        allocate(x(5),root(5))
        x=(/highT/2d0,rho,p,t1,t2/)
        ptr=>fun_Tp
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solvetp=root(1)
        deallocate(x,root)
    end if
end function solvetp

function prhot(rho,t,temp1,temp2)
    ! calculate pressure with given rho and t
    real(8) :: prhot,rho,t,t1,t2,species(3),nhtot,nhetot,ntot
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then                                     ! H2+He
        nhtot=x_ratio*rho/mh                            
        nhetot=(1d0-x_ratio)*rho/mhe                    
        ntot=nhtot/2d0+nhetot                           ! # of particles
        prhot=ntot*kb*t
    else if (t>=t2) then                                ! H+He
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=nhtot+nhetot
        prhot=ntot*kb*t
    else                                                ! H2+H+He
        call calspecies(rho,t,species,t1,t2)
        prhot=sum(species)*kb*t
    end if
end function prhot

function fun_Tp(x)
    !x(1)=t,x(2)=rho,x(3)=p
    real(8), dimension(:), allocatable :: x
    real(8) :: fun_Tp,t,rho,p,t1,t2
    t=x(1)
    rho=x(2)
    p=x(3)
    t1=x(4)
    t2=x(5)
    fun_Tp=p-prhot(rho,t,t1,t2)
end function fun_Tp

function solvets(s,rho)
    !given specific entropy and rho, find the corresponding temperature
    real(8) :: solvets,s,rho,lowT,highT,lows,highs,convergence,t,t1,t2
    real(8) :: ah2,ah,ahe,bh2,bh,bhe,nh2tot,nhtot,nhetot
    procedure(fun), pointer :: ptr
    real(8), dimension(:), allocatable :: x,root
    integer :: conv_mode
    conv_mode=2
    call divide_temperature_domain_known(rho,t1,t2)
    lowT=t1
    highT=t2
    convergence=1d-6
    lows=srhot(rho,lowT,t1,t2)
    highs=srhot(rho,highT,t1,t2)
    if (s<=lows) then                                                                       ! Pure H2
        nh2tot=x_ratio*rho/mh2                                                              ! # of H
        nhetot=(1d0-x_ratio)*rho/mhe
        ah2=kb*nh2tot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mh2*kb/h_planck**2)-log(nh2tot))
        ahe=kb*nhetot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mhe*kb/h_planck**2)-log(nhetot))
        bh2=kb*nh2tot/rho*1.5d0
        bhe=kb*nhetot/rho*1.5d0
        solvets=exp((s-ah2-ahe)/(bh2+bhe))
    else if (s>=highs) then                                                                 ! pure H     
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ah=kb*nhtot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mh*kb/h_planck**2)-log(nhtot))
        ahe=kb*nhetot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mhe*kb/h_planck**2)-log(nhetot))
        bh=kb*nhtot/rho*1.5d0
        bhe=kb*nhetot/rho*1.5d0
        solvets=exp((s-ah-ahe)/(bh+bhe))
    else
        allocate(x(5),root(5))
        x=(/highT/2d0,rho,s,t1,t2/)
        ptr=>fun_Ts
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solvets=root(1)
        deallocate(x,root)
    end if
end function solvets

function srhot(rho,t,temp1,temp2)
    real(8) :: srhot,rho,t,t1,t2
    real(8),dimension(3) :: species,dlnzdlnt,d2lnzdlnt,mu
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    call calspecies(rho,t,species,t1,t2)
    call caldlnz(t,dlnzdlnt,d2lnzdlnt)
    call calmu(rho,t,species,mu,t1,t2)
    srhot=kb/rho*dot_product(species,1d0+dlnzdlnt-mu/kb/t)
end function srhot

function srhot_per_particle(rho,t,temp1,temp2)
    real(8) :: srhot,rho,t,t1,t2,n,srhot_per_particle
    real(8),dimension(3) :: species,dlnzdlnt,d2lnzdlnt,mu
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    call calspecies(rho,t,species,t1,t2)
    call caldlnz(t,dlnzdlnt,d2lnzdlnt)
    call calmu(rho,t,species,mu,t1,t2)
    n=sum(species)
    srhot=kb/rho*dot_product(species,1d0+dlnzdlnt-mu/kb/t)
    srhot_per_particle=srhot/kb/n*rho
end function srhot_per_particle

function fun_Ts(x)
    !x(1)=t,x(2)=rho,x(3)=s,x(4-5)=t1,t2
    real(8), dimension(:), allocatable :: x
    real(8) :: rho,t,s,fun_Ts,t1,t2,t3,t4
    t=x(1)
    rho=x(2)
    s=x(3)
    t1=x(4)
    t2=x(5)
    fun_Ts=s-srhot(rho,t,t1,t2)
end function fun_Ts

function solvetegv(egv,rho)
    !given egv and rho, find the corresponding temperature
    real(8) :: solvetegv,egv,rho,lowT,highT,lowegv,highegv,t_temp,t,t1,t2
    real(8) :: convergence,nh2tot,nhtot,nhetot,ntot
    procedure(fun), pointer :: ptr
    real(8), dimension(:), allocatable :: x,root
    integer :: conv_mode
    conv_mode=2
    convergence=1e-6
    call divide_temperature_domain_known(rho,t1,t2)
    lowT=t1
    highT=t2
    lowegv=egvrhot(rho,lowT,t1,t2)
    highegv=egvrhot(rho,highT,t1,t2)
    if (egv<=lowegv) then                           ! H2
        nh2tot=x_ratio*rho/mh2
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=nh2tot+nhetot
        solvetegv=2d0/3d0*egv/kb/ntot
    else if (egv>=highegv) then                     ! H
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=nhtot+nhetot
        solvetegv=2d0/3d0*(egv-nhtot*dish/2d0)/kb/ntot
    else
        allocate(x(5),root(5))
        x=(/highT/2d0,rho,egv,t1,t2/)
        ptr=>fun_Tegv
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solvetegv=root(1)
        deallocate(x,root)
    end if
end function solvetegv

function egvrhot(rho,t,temp1,temp2)
    real(8) :: rho,t,egvrhot,t1,t2,nh2tot,nhtot,nhetot
    real(8),dimension(3) :: energy,species
    real(8), optional :: temp1,temp2
    energy(1)=1.5d0*kb*t                                    ! H2
    energy(2)=1.5d0*kb*t+dish/2d0                           ! H
    energy(3)=1.5d0*kb*t                                    ! He
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then                                         ! H2
        nh2tot=x_ratio*rho/mh2
        nhetot=(1d0-x_ratio)*rho/mhe
        egvrhot=nh2tot*energy(1)+nhetot*energy(3)
    else if (t>=t2) then                                    ! H
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        egvrhot=nhtot*energy(2)+nhetot*energy(3)
    else
        call calspecies(rho,t,species,t1,t2)
        egvrhot=dot_product(species,energy)
    end if
end function egvrhot

function fun_Tegv(x)
    !x(1)=t,x(2)=rho,x(3)=egv,x(4-5)=t1,t2
    !calculate the root function
    real(8), dimension(:), allocatable :: x
    real(8) :: rho,t,egv,fun_Tegv,t1,t2
    t=x(1)
    rho=x(2)
    egv=x(3)
    t1=x(4)
    t2=x(5)
    fun_Tegv=egv-egvrhot(rho,t,t1,t2)
end function fun_Tegv

subroutine divide_temperature_domain_known(rho,t1,t2)   
    real(8) :: rho,t1,t2
    t1=10d0**get_t1(rho*x_ratio)                    
    t2=10d0**get_t2(rho*x_ratio)                    
end subroutine divide_temperature_domain_known

function get_t1(rho)
    real(8) :: get_t1,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=max(log10(rho),log10_rhomin)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    get_t1=temp_division(i,1)+(temp_division(i+1,1)-temp_division(i,1))*dd/dlogrho  
end function get_t1

function get_t2(rho)
    real(8) :: get_t2,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=max(log10(rho),log10_rhomin)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    get_t2=temp_division(i,2)+(temp_division(i+1,2)-temp_division(i,2))*dd/dlogrho
end function get_t2

subroutine calspecies(rho,t,species,temp1,temp2)        ! calculate nh2,nh
    real(8) :: rho,t,species(3),t1,t2,nhtot,nhetot
    real(8), optional :: temp1,temp2
    nhtot=x_ratio*rho/mh                    
    nhetot=(1d0-x_ratio)*rho/mhe                    
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then
        species(1)=nhtot/2d0                            ! H2
        species(2)=0d0                                  ! H
    else if (t>=t2) then
        species(1)=0d0
        species(2)=nhtot
    else
        species(2)=solvenhdis(rho,t)
        species(1)=species(2)**2/kdis(t)
    end if
    species(3)=nhetot                                   ! He
end subroutine calspecies

function solvenhdis(rho,t)
    !solve the number density of H at given rho and t
    real(8) :: rho,t,solvenhdis,a,nhtot
    a=kdis(t)
    nhtot=x_ratio*rho/mh
    solvenhdis=2d0*a*nhtot/(a+sqrt(a**2+8d0*a*nhtot))   
end function solvenhdis

function kdis(t)
    real(8) :: kdis,t
    kdis=zh(t)*zh(t)/zh2(t)
end function kdis

function zh(t)
    real(8) :: zh,t
    zh=zsimplegas(t,mh)/(exp(dish/(2d0*kb*t)))
end function zh

function zh2(t)
    real(8) :: zh2,t
    zh2=zsimplegas(t,mh2)
end function zh2

function zhe(t)
    real(8) :: zhe,t
     zhe=zsimplegas(t,mhe)
end function zhe

subroutine caldlnz(t,dlnzdlnt,d2lnzdlnt)
    real(8) :: t,dlnzdlnt(3),d2lnzdlnt(3)
    dlnzdlnt(1)=1.5d0
    dlnzdlnt(2)=1.5d0+dish/2d0/kb/t
    dlnzdlnt(3)=1.5d0
    d2lnzdlnt(1)=0d0
    d2lnzdlnt(2)=-dish/2d0/kb/t
    d2lnzdlnt(3)=0d0
end subroutine caldlnz

subroutine caldlnndlnrho(rho,t,species,dlnndlnrho,temp1,temp2)
    real(8) :: rho,t,dlnndlnrho(3),rhs(2),species(3),t1,t2
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then                                     ! pure H2
        dlnndlnrho(1)=x_ratio*rho/mh2/species(1)
        dlnndlnrho(2)=0d0
    else if (t>=t2) then                                ! pure H
        dlnndlnrho(1)=0d0
        dlnndlnrho(2)=x_ratio*rho/mh/species(2)
    else
        rhs(1)=x_ratio*rho/mh                           ! # of H atom
        rhs(2)=0
        dlnndlnrho(2)=rhs(1)/(species(2)+4d0*species(1))
        dlnndlnrho(1)=2d0*dlnndlnrho(2)
    end if
    if (x_ratio==1d0) then
        dlnndlnrho(3)=0d0                               ! He
    else
        dlnndlnrho(3)=(1d0-x_ratio)*rho/mhe/species(3)
    end if
end subroutine caldlnndlnrho

subroutine caldlnndlnt(rho,t,species,dlnndlnt,temp1,temp2)
    real(8) :: rho,t,dlnndlnt(3),rhs(2),species(3)
    real(8), optional :: temp1,temp2
    real(8) :: ph2,ph,t1,t2
    ph2=1.5d0                                                       ! dlnz/dlnT
    ph=1.5d0+dish/kb/t/2d0
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then                                                 ! H2                                    
        dlnndlnt=0d0
    else if (t>=t2) then                                            ! H
        dlnndlnt=0d0
    else
        rhs(1)=0d0
        rhs(2)=2d0*ph-ph2
        dlnndlnt(2)=2d0*rhs(2)*species(1)/(4d0*species(1)+species(2))
        dlnndlnt(1)=-species(2)*dlnndlnt(2)/2d0/species(1)
    end if
    dlnndlnt(3)=0d0                                                 ! He             
end subroutine caldlnndlnt

subroutine generate_temperature_division(rho,t,temp_division)
    ! to get 2 rho-t relations that divide the rho-t EOS into 3 blocks
    ! (H2), (H2,HI), (HI)
    real(8), dimension(:,:,:), allocatable :: xi,division
    real(8), dimension(:,:), allocatable :: temp_division
    real(8), dimension(:), allocatable :: rho,t
    real(8) :: log10_tmax,delta,xi_temp(2)
    integer :: i,j,k,ntemp,nrho
    nrho=size(rho)
    ntemp=size(t)
    allocate(division(ntemp,nrho,2),xi(ntemp,nrho,2))
    do i=1,nrho
        do j=1,ntemp
            call calxi_temp_division(rho(i),t(j),xi_temp)
            xi(j,i,:)=xi_temp
        end do
    end do
    log10_tmax=log10(t(ntemp))
    delta=1d-6
    division=-1d0
    where (xi(:,:,:)>delta.and.(1d0-xi(:,:,:))>delta)
        division(:,:,:)=1d0
    end where
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,2)<0.and.division(j+1,i,2)>0) then !HI emerges
                temp_division(i,1)=log10(t(j))
                exit
            else
                temp_division(i,1)=log10_tmax
            end if
        end do
    end do
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,1)>0.and.division(j+1,i,1)<0) then !H2 disappers
                temp_division(i,2)=log10(t(j))
                exit
            else
                temp_division(i,2)=log10_tmax
            end if
        end do
    end do
    division_known=.true.
    deallocate(division,xi)
end subroutine generate_temperature_division

subroutine calxi_temp_division(rho,t,xi)
    real(8) :: rho,t,xi(2),species(3)
    call calspecies_division_unk(rho,t,species)
    xi(1)=max(species(1)/x_ratio/rho*mh2,1d-16)         ! chi_h2
    xi(2)=max(species(2)/x_ratio/rho*mh,1d-16)          ! chi_H
end subroutine calxi_temp_division

subroutine calspecies_division_unk(rho,t,species)
    real(8) :: rho,t,t1,t2,species(3),nhtot,nhetot
    nhtot=x_ratio*rho/mh
    nhetot=(1d0-x_ratio)*rho/mhe
    call divide_temperature_domain(rho,t1,t2)
    if (t<=t1) then
        !fully molecular
        species(1)=nhtot/2d0
        species(2)=0d0
    else if (t>=t2) then
        !fully ionized
        species(1)=0d0
        species(2)=nhtot
    else
        !mixed
        species(2)=solvenhdis(rho,t)
        species(1)=species(2)**2/kdis(t)
    end if
    species(3)=nhetot
end subroutine calspecies_division_unk

subroutine divide_temperature_domain(rho,t1,t2)
    real(8) :: rho,t1,t2
    t1=temperature1(rho)
    t2=temperature2(rho)
end subroutine divide_temperature_domain

function temperature1(rho)  
    real(8) :: temperature1,rho
    temperature1=10**2.7
end function temperature1

function temperature2(rho)
    real(8) :: temperature2,rho,k1,k2
    k1=0.045d0
    k2=1.5d0
    if (rho>1e-10) then
        temperature2=min(10**(3.2d0+k1*10d0+k2*(log10(rho)+10d0)),1.1e6)
    else
    temperature2=10**(3.2d0+k1*(log10(rho)+20d0))
    end if
end function temperature2

subroutine calprhopt(species,dlnndlnrho,dlnndlnt,rho,t,prho,pt)
    real(8) :: species(3),rho,t,dlnndlnrho(3),dlnndlnt(3),prho,pt,ntot
    integer :: i
    ntot=sum(species)
    prho=0d0
    pt=0d0
    do i=1,3                                        
        prho=prho+species(i)/ntot*dlnndlnrho(i)
        pt=pt+species(i)/ntot*dlnndlnt(i)
    end do
    pt=pt+1                                         
end subroutine calprhopt

subroutine calsrhost(rho,t,species,dlnndlnrho,dlnndlnt,dlnzdlnt,d2lnzdlnt,srho,st)
    ! (partial(s)(rho))T, (partial(s)(T))rho
    real(8) :: rho,t,srho,st
    real(8), dimension(3) :: dsdrho,dsdt,dndrho,dndt,mu,species,dlnndlnrho,dlnndlnt,dlnzdlnt,d2lnzdlnt
    integer :: i
    call calmu(rho,t,species,mu)
    do i=1,3
        dndrho(i)=species(i)/rho*dlnndlnrho(i)
        dndt(i)=species(i)/t*dlnndlnt(i)
        dsdt(i)=kb*species(i)/rho/t*(d2lnzdlnt(i)+(1d0+dlnndlnt(i))*dlnzdlnt(i))-dndt(i)*mu(i)/rho/t
        dsdrho(i)=kb*species(i)/rho**2*((dlnndlnrho(i)-1d0)*dlnzdlnt(i)-1d0)+(species(i)/rho-dndrho(i))*mu(i)/rho/t
    end do
    srho=sum(dsdrho)
    st=sum(dsdt)
end subroutine calsrhost

subroutine calmu(rho,t,species,mu,temp1,temp2)
    real(8) :: rho,t,species(3),mu(3),t1,t2
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then                                     ! H2
        mu(1)=kb*t*log(species(1)/zh2(t))
        mu(2)=0d0
    else if (t>=t2) then                                ! H
        mu(1)=0d0
        mu(2)=kb*t*log(species(2)/zh(t))
    else
        mu(1)=kb*t*log(species(1)/zh2(t))
        mu(2)=kb*t*log(species(2)/zh(t))
    end if
    if (x_ratio==1d0) then
        mu(3)=0d0
    else
        mu(3)=kb*t*log(species(3)/zhe(t))
    end if
end subroutine calmu

subroutine caldndt(rho,t,species,dndt)
    real(8) :: rho,t,species(3),dndt(3),dlnndlnt(3)
    integer :: i
    call caldlnndlnt(rho,t,species,dlnndlnt)
    do i=1,3
        dndt(i)=species(i)/t*dlnndlnt(i)
    end do
end subroutine caldndt

subroutine eos_species(blk)
    type(blockdef), pointer :: blk
    real(8) :: rho,temp,nhtot,species(3)
    integer :: i
    do i=1,blk_size_nx
        rho=blk%w(1,i,1,1)
        temp=blk%temp(i,1,1)
        call calspecies(rho,temp,species)
        nhtot=2d0*species(1)+species(2)
        blk%H2(i,1,1)=2*species(1)/nhtot
        blk%HI(i,1,1)=species(2)/nhtot
        blk%cs(i,1,1)=adiabatic_cs(rho,temp)
    end do
end subroutine eos_species

end module eos_h2_HI

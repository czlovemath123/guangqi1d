module eos_HI_HII
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
    real(8) :: rho,t,gamm,p,ntot,srho,st,prho,pt,t1,t2,species(4),dlnzdlnt(4),d2lnzdlnt(4)
    real(8) :: dlnndlnrho(4),dlnndlnt(4)
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
    real(8) :: adiabatic_cs,rho,t,gamm,p,ntot,srho,st,prho,pt,t1,t2,species(4),dlnzdlnt(4),d2lnzdlnt(4)
    real(8) :: dlnndlnrho(4),dlnndlnt(4)
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
    real(8) :: eos_cv,rho,t,t1,t2,xi_h,xi_hion,xi_e,xi_he
    real(8) :: xi(4),species(4),dlnndlnt(4),dndt(4)
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
    xi_h=1.5d0*kb*t
    xi_hion=1.5d0*kb*t+ionh
    xi_e=1.5d0*kb*t
    xi_he=1.5d0*kb*t
    xi=(/xi_h,xi_hion,xi_e,xi_he/)
    eos_cv=dot_product(dndt,xi)+sum(species)*1.5d0*kb
end function eos_cv

function solvetchi_HII(chi_HII,rho)
    !given chi_HII and rho, find the corresponding temperature
    real(8) :: solvetchi_HII,chi_HII,rho,lowT,lowchi_HII,highT,highchi_HII,convergence,t,t1,t2,nhtot,nhetot,ntot
    real(8), dimension(:), allocatable :: x,root
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    conv_mode=2
    convergence=1e-6
    call divide_temperature_domain_known(rho,t1,t2)
    lowT=t1
    highT=t2
    lowchi_HII=chi_HIIrhot(rho,lowT,t1,t2)
    highchi_HII=chi_HIIrhot(rho,highT,t1,t2)
    if (chi_HII<lowchi_HII) then
        print *,'solve chi_HII degenerate'
        stop
    else if (chi_HII>highchi_HII) then
        print *,'solve chi_HII degenerate'
        stop
    else
        allocate(x(5),root(5))
        x=(/(lowT+highT)/2d0,rho,chi_HII,t1,t2/)
        ptr=>fun_Tchi_HII
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solvetchi_HII=root(1)
        deallocate(x,root)
    end if
end function solvetchi_HII

function chi_HIIrhot(rho,t,temp1,temp2)
    !given rho and t, calculate the ionization fraction
    real(8) :: chi_HIIrhot,rho,t,t1,t2,species(4),nhtot,nhetot,ntot
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then
        chi_HIIrhot=0d0
    else if (t>=t2) then
        chi_HIIrhot=1d0
    else
        call calspecies(rho,t,species,t1,t2)
        chi_HIIrhot=species(2)/(species(1)+species(2))
    end if
end function chi_HIIrhot

function fun_Tchi_HII(x)
    real(8), dimension(:), allocatable :: x
    real(8) :: fun_Tchi_HII,t,rho,chi_HII,t1,t2
    t=x(1)
    rho=x(2)
    chi_HII=x(3)
    t1=x(4)
    t2=x(5)
    fun_Tchi_HII=chi_HII-chi_HIIrhot(rho,t,t1,t2)
end function fun_Tchi_HII

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
        ntot=nhtot+nhetot
        solvetp=p/kb/ntot
    else if (p>=highp) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=2d0*nhtot+nhetot
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
    real(8) :: prhot,rho,t,t1,t2,species(4),nhtot,nhetot,ntot
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=nhtot+nhetot
        prhot=ntot*kb*t
    else if (t>=t2) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=2d0*nhtot+nhetot
        prhot=ntot*kb*t
    else
        call calspecies(rho,t,species,t1,t2)
        prhot=sum(species)*kb*t
    end if
end function prhot

function fun_Tp(x)
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
    real(8) :: solvets,s,rho,lowT,highT,lows,highs,convergence,t,t1,t2,ah,ae,ahe,bh,be,bhe,nhtot,nhetot,mhion
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
    if (s<=lows) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ah=kb*nhtot/rho*(2.5d0+log(four)+1.5d0*log(2d0*pi*mh*kb/h_planck**2)-log(nhtot))
        ahe=kb*nhetot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mhe*kb/h_planck**2)-log(nhetot))
        bh=kb*nhtot/rho*1.5d0
        bhe=kb*nhetot/rho*1.5d0
        solvets=exp((s-ah-ahe)/(bh+bhe))
    else if (s>=highs) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        mhion=mh-me
        ah=kb*nhtot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mhion*kb/h_planck**2)-log(nhtot))
        ae=kb*nhtot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*me*kb/h_planck**2)-log(nhtot))
        ahe=kb*nhetot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mhe*kb/h_planck**2)-log(nhetot))
        bh=kb*nhtot/rho*1.5d0
        be=kb*nhtot/rho*1.5d0
        bhe=kb*nhetot/rho*1.5d0
        solvets=exp((s-ah-ae-ahe)/(bh+be+bhe))
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
    real(8) :: srhot,rho,t,t1,t2,species(4),dlnzdlnt(4),d2lnzdlnt(4),mu(4),mu2(4)
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
    mu2=mu+kb*t*(/log(2d0),log(2d0),log(1d0),log(2d0)/)
    srhot=kb/rho*dot_product(species,1d0+dlnzdlnt-mu2/kb/t)
end function srhot

function srhot_per_particle(rho,t,temp1,temp2)
    real(8) :: srhot,rho,t,t1,t2,species(4),dlnzdlnt(4),d2lnzdlnt(4),mu(4),mu2(4),n,srhot_per_particle
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
    mu2=mu+kb*t*(/log(2d0),log(2d0),log(1d0),log(2d0)/)
    n=sum(species)
    srhot=kb/rho*dot_product(species,1d0+dlnzdlnt-mu2/kb/t)
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
    real(8) :: convergence,nhtot,nhetot,ntot
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
    if (egv<=lowegv) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=nhtot+nhetot
        solvetegv=2d0/3d0*egv/kb/ntot
    else if (egv>=highegv) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=2d0*nhtot+nhetot
        solvetegv=2d0/3d0*(egv-nhtot*ionh)/kb/ntot
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
    real(8) :: rho,t,egvrhot,t1,t2,energy(4),species(4),nhtot,nhetot
    real(8), optional :: temp1,temp2
    energy=(/1.5d0*kb*t,1.5d0*kb*t+ionh,1.5d0*kb*t,1.5d0*kb*t/)
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        egvrhot=nhtot*energy(1)+nhetot*energy(4)
    else if (t>=t2) then
        nhtot=x_ratio*rho/mh
        nhetot=(1d0-x_ratio)*rho/mhe
        egvrhot=nhtot*(energy(2)+energy(3))+nhetot*energy(4)
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

subroutine calsrhost(rho,t,species,dlnndlnrho,dlnndlnt,dlnzdlnt,d2lnzdlnt,srho,st)
    real(8) :: rho,t,dsdrho(4),dsdt(4),srho,st,species(4),mu(4),dlnndlnrho(4),dlnndlnt(4),  &
        dndrho(4),dndt(4),dlnzdlnt(4),d2lnzdlnt(4)
    integer :: i
    call calmu(rho,t,species,mu)
    do i=1,4
        dndrho(i)=species(i)/rho*dlnndlnrho(i)
        dndt(i)=species(i)/t*dlnndlnt(i)
        dsdt(i)=kb*species(i)/rho/t*(d2lnzdlnt(i)+(1d0+dlnndlnt(i))*dlnzdlnt(i))-dndt(i)*mu(i)/rho/t
        dsdrho(i)=kb*species(i)/rho**2*((dlnndlnrho(i)-1d0)*dlnzdlnt(i)-1d0)+(species(i)/rho-dndrho(i))*mu(i)/rho/t
    end do
    srho=sum(dsdrho)
    st=sum(dsdt)
end subroutine calsrhost

subroutine calprhopt(species,dlnndlnrho,dlnndlnt,rho,t,prho,pt)
    real(8) :: species(4),rho,t,dlnndlnrho(4),dlnndlnt(4),prho,pt,ntot
    integer :: i
    ntot=sum(species)
    prho=0d0
    pt=0d0
    do i=1,4
        prho=prho+species(i)/ntot*dlnndlnrho(i)
        pt=pt+species(i)/ntot*dlnndlnt(i)
    end do
    pt=pt+1
end subroutine calprhopt

subroutine caldlnz(t,dlnzdlnt,d2lnzdlnt)
    real(8) :: t,dlnzdlnt(4),d2lnzdlnt(4)
    dlnzdlnt(1)=1.5d0
    dlnzdlnt(2)=1.5d0+ionh/kb/t
    dlnzdlnt(3)=1.5d0
    dlnzdlnt(4)=1.5d0
    d2lnzdlnt(1)=0d0
    d2lnzdlnt(2)=-ionh/kb/t
    d2lnzdlnt(3)=0d0
    d2lnzdlnt(4)=0d0
end subroutine caldlnz

subroutine calmu(rho,t,species,mu,temp1,temp2)
    real(8) :: rho,t,species(4),mu(4),t1,t2,t3,t4
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then
        mu(1)=kb*t*log(species(1)/zh(t))
        mu(2)=0d0
        mu(3)=0d0
    else if (t>=t2) then
        mu(1)=0d0
        mu(2)=kb*t*log(species(2)/zhion(t))
        mu(3)=kb*t*log(species(3)/ze(t))
    else
        mu(1)=kb*t*log(species(1)/zh(t))
        mu(2)=kb*t*log(species(2)/zhion(t))
        mu(3)=kb*t*log(species(3)/ze(t))
    end if
    if (x_ratio==1d0) then
        mu(4)=0d0
    else
        mu(4)=kb*t*log(species(4)/zhe(t))
    end if
end subroutine calmu

subroutine caldlnndlnrho(rho,t,species,dlnndlnrho,temp1,temp2)
    real(8) :: rho,t,dlnndlnrho(4),rhs(2),species(4),t1,t2
    real(8), optional :: temp1,temp2
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then
        dlnndlnrho(1:3)=(/1d0,0d0,0d0/)
    else if (t>=t2) then
        dlnndlnrho(1:3)=(/0d0,1d0,1d0/)
    else
        dlnndlnrho(2)=x_ratio*rho/mh/(species(2)+2d0*species(1))
        dlnndlnrho(3)=dlnndlnrho(2)
        dlnndlnrho(1)=2d0*dlnndlnrho(2)
    end if
    if (x_ratio==1d0) then
        dlnndlnrho(4)=0d0
    else
        dlnndlnrho(4)=1d0
    end if
end subroutine caldlnndlnrho

subroutine caldndt(rho,t,species,dndt)
    real(8) :: rho,t,species(4),dndt(4),dlnndlnt(4)
    integer :: i
    call caldlnndlnt(rho,t,species,dlnndlnt)
    do i=1,4
        dndt(i)=species(i)/t*dlnndlnt(i)
    end do
end subroutine caldndt

subroutine caldlnndlnt(rho,t,species,dlnndlnt,temp1,temp2)
    real(8) :: rho,t,dlnndlnt(4),rhs(2),species(4)
    real(8), optional :: temp1,temp2
    real(8) :: phion,pe,ph,t1,t2
    ph=1.5d0
    pe=1.5d0
    phion=1.5d0+ionh/kb/t
    if (present(temp1)) then
        t1=temp1
        t2=temp2
    else
        call divide_temperature_domain_known(rho,t1,t2)
    end if
    if (t<=t1) then
        dlnndlnt=0d0
    else if (t>=t2) then
        dlnndlnt=0d0
    else
        rhs(1)=0d0
        rhs(2)=phion+pe-ph
        dlnndlnt(2)=species(1)*rhs(2)/(2d0*species(1)+species(2))
        dlnndlnt(3)=dlnndlnt(2)
        dlnndlnt(1)=-species(2)/species(1)*dlnndlnt(2)
    end if
    dlnndlnt(4)=0d0
end subroutine caldlnndlnt

subroutine generate_temperature_division(rho,t,temp_division)
    real(8), allocatable :: xi(:,:,:),division(:,:,:),temp_division(:,:),rho(:),t(:)
    real(8) :: log10_tmax,delta,xi_temp(3)
    integer :: i,j,k,ntemp,nrho
    nrho=size(rho)
    ntemp=size(t)
    allocate(division(ntemp,nrho,3),xi(ntemp,nrho,3))
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
            if (division(j,i,2)<0.and.division(j+1,i,2)>0) then !HII emerges
                temp_division(i,1)=log10(t(j))
                exit
            else
                temp_division(i,1)=log10_tmax
            end if
        end do
    end do
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,1)>0.and.division(j+1,i,1)<0) then !HI disappers
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
    real(8) :: rho,t,xi(3),species(4)
    call calspecies_division_unk(rho,t,species)
    xi(1)=species(1)/x_ratio/rho*mh
    xi(2)=species(2)/x_ratio/rho*mh
    xi(3)=xi(2)
end subroutine calxi_temp_division

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
    log10_rho=log10(rho)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    get_t1=temp_division(i,1)+(temp_division(i+1,1)-temp_division(i,1))*dd/dlogrho
end function get_t1

function get_t2(rho)
    real(8) :: get_t2,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=log10(rho)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    get_t2=temp_division(i,2)+(temp_division(i+1,2)-temp_division(i,2))*dd/dlogrho
end function get_t2

subroutine calspecies_division_unk(rho,t,species)
    real(8) :: rho,t,t1,t2,species(4),nhtot,nhetot
    nhtot=x_ratio*rho/mh
    nhetot=(1d0-x_ratio)*rho/mhe
    call divide_temperature_domain(rho,t1,t2)
    if (t<=t1) then
        !fully atomic
        species(1)=nhtot
        species(2)=0d0
        species(3)=species(2)
    else if (t>=t2) then
        !fully ionized
        species(1)=0d0
        species(2)=nhtot
        species(3)=species(2)
    else
        !mixed
        species(2)=solvenhion(rho,t)
        species(3)=species(2)
        species(1)=species(3)**2/kion(t)
    end if
    species(4)=nhetot
end subroutine calspecies_division_unk

subroutine divide_temperature_domain(rho,t1,t2)
    real(8) :: rho,t1,t2
    t1=temperature1(rho)
    t2=temperature2(rho)
end subroutine divide_temperature_domain

function temperature1(rho)
    real(8) :: temperature1,rho,k
    k=0.02d0
    temperature1=10d0**(3d0+k*(log10(rho)+20d0))
end function temperature1

function temperature2(rho)
    real(8) :: temperature2,rho,k
    k=2d0/10d0
    temperature2=min(10d0**(4d0+k*(log10(rho)+20d0)),1.1e6)
end function temperature2

subroutine calspecies(rho,t,species,temp1,temp2)
    !calculate nh,nhion and ne
    real(8) :: rho,t,species(4),t1,t2,nhtot,nhetot
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
        species(1)=nhtot
        species(2)=0d0
        species(3)=0d0
    else if (t>=t2) then
        species(1)=0d0
        species(2)=nhtot
        species(3)=species(2)
    else
        species(2)=solvenhion(rho,t)
        species(1)=species(2)**2/kion(t)
        species(3)=species(2)
    end if
    species(4)=nhetot
end subroutine calspecies

function solvenhion(rho,t)
    !solve the number density of HII at given rho and t
    real(8) :: rho,t,solvenhion,kion_temp,a,nhtot
    a=kion(t)
    nhtot=x_ratio*rho/mh
    solvenhion=2d0*a*nhtot/(a+sqrt(a**2+4d0*a*nhtot))
end function solvenhion

function kion(t)
    real(8) :: kion,t
    kion=zhion(t)*ze(t)/zh(t)
end function kion

function zh(t)
    real(8) :: zh,t,ztr,zspin,zelec
    zh=zsimplegas(t,mh)
end function zh

function ze(t)
    real(8) :: ze,t,ztr,zspin
    ze=zsimplegas(t,me)
end function ze

function zhion(t)
    real(8) :: zhion,t,ztr,zspin,zelec,mhion
    zhion=zsimplegas(t,mh)/(exp((ionh)/(kb*t)))
end function zhion

function zhe(t)
    real(8) :: zhe,t
    zhe=zsimplegas(t,mhe)
end function zhe

subroutine eos_species(blk)
    type(blockdef), pointer :: blk
    real(8) :: rho,temp,nhtot,species(4)
    integer :: i
    do i=1,blk_size_nx
        rho=blk%w(1,i,1,1)
        temp=blk%temp(i,1,1)
        call calspecies(rho,temp,species)
        nhtot=species(1)+species(2)
        blk%HI(i,1,1)=species(1)/nhtot
        blk%HII(i,1,1)=species(2)/nhtot
        blk%cs(i,1,1)=adiabatic_cs(rho,temp)
    end do
end subroutine eos_species

end module eos_HI_HII

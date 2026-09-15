module eos_h2_HI_HII
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
    !calculate the Gamma with given rho and t
    real(8) :: rho,t,gamm,p,ntot,srho,st,prho,pt,t1,t2,t3,t4
    real(8), dimension(5) :: species,dlnzdlnt,d2lnzdlnt,dlnndlnrho,dlnndlnt
    call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    call calspecies(rho,t,species,t1,t2,t3,t4)
    call caldlnz(t,dlnzdlnt,d2lnzdlnt)
    call caldlnndlnrho(rho,t,species,dlnndlnrho,t1,t2,t3,t4)
    call caldlnndlnt(rho,t,species,dlnndlnt,t1,t2,t3,t4)
    call calprhopt(species,dlnndlnrho,dlnndlnt,rho,t,prho,pt)
    call calsrhost(rho,t,species,dlnndlnrho,dlnndlnt,dlnzdlnt,d2lnzdlnt,srho,st)
    ntot=sum(species)
    p=ntot*kb*t
    gamm=rho/p*(p/rho*prho-p/t*pt*srho/st)
end subroutine gammarhot

function adiabatic_cs(rho,t)
    !calculate the adiabatic sound speed with given rho and t
    real(8) :: adiabatic_cs,rho,t,gamm,p,ntot,srho,st,prho,pt,t1,t2,t3,t4
    real(8), dimension(5) :: species,dlnzdlnt,d2lnzdlnt,dlnndlnrho,dlnndlnt
    call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    call calspecies(rho,t,species,t1,t2,t3,t4)
    call caldlnz(t,dlnzdlnt,d2lnzdlnt)
    call caldlnndlnrho(rho,t,species,dlnndlnrho,t1,t2,t3,t4)
    call caldlnndlnt(rho,t,species,dlnndlnt,t1,t2,t3,t4)
    call calprhopt(species,dlnndlnrho,dlnndlnt,rho,t,prho,pt)
    call calsrhost(rho,t,species,dlnndlnrho,dlnndlnt,dlnzdlnt,d2lnzdlnt,srho,st)
    ntot=sum(species)
    p=ntot*kb*t
    adiabatic_cs=sqrt(p/rho*prho-p/t*pt*srho/st)
end function adiabatic_cs

function eos_cv(rho,t,temp1,temp2,temp3,temp4)
    real(8) :: eos_cv,rho,t,t1,t2,t3,t4,xi_h2,xi_h,xi_hion,xi_e,xi_he,xi(5)
    real(8), optional :: temp1,temp2,temp3,temp4
    real(8), dimension(5) :: species,dlnndlnt,dndt
    if (present(temp1)) then
        t1=temp1
        t2=temp2
        t3=temp3
        t4=temp4
    else
        call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    end if
    call calspecies(rho,t,species,t1,t2,t3,t4)
    call caldlnndlnt(rho,t,species,dlnndlnt,t1,t2,t3,t4)
    call caldndt(rho,t,species,dndt)
    xi_h2=1.5d0*kb*t
    xi_h=1.5d0*kb*t+0.5d0*dish
    xi_hion=1.5d0*kb*t+0.5d0*dish+ionh
    xi_e=1.5d0*kb*t
    xi_he=1.5d0*kb*t
    xi=(/xi_h2,xi_h,xi_hion,xi_e,xi_he/)
    eos_cv=dot_product(dndt,xi)+sum(species)*1.5d0*kb
end function eos_cv

function solvetp(p,rho)
    !given p and rho, find the corresponding temperature
    real(8) :: solvetp,p,rho,lowT,lowp,highT,highp,convergence,t,nhetot,nhtot,ntot,t1,t2,t3,t4
    real(8), dimension(:), allocatable :: x,root
    procedure(fun), pointer :: ptr1,ptr2,ptr
    integer :: conv_mode
    conv_mode=2
    call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    lowT=t1
    highT=t4
    lowp=prhot(rho,lowT,t1,t2,t3,t4)                        !all molecular hydrogen
    highp=prhot(rho,highT,t1,t2,t3,t4)                      !all ionized hydrogen
    convergence=1d-5
    if (p<=lowp) then
        nhtot=x_ratio*rho/mh2   !all molecular
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=nhtot+nhetot
        solvetp=p/kb/ntot
    else if (p>=highp) then
        nhtot=x_ratio*rho/mh    !fully ionized
        nhetot=(1d0-x_ratio)*rho/mhe
        ntot=2*nhtot+nhetot
        solvetp=p/kb/ntot
    else if (p>lowp.and.p<highp) then
        allocate(x(7),root(7))
        x=(/highT/2d0,rho,p,t1,t2,t3,t4/)
        ptr=>fun_Tp
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solvetp=root(1)
        deallocate(x,root)
    end if
end function solvetp

function prhot(rho,t,temp1,temp2,temp3,temp4)
    !calculate pressure with given rho and t
    real(8) :: prhot,rho,t,ntot,t1,t2,t3,t4
    real(8), optional :: temp1,temp2,temp3,temp4
    real(8), dimension(5) :: species
    if (present(temp1)) then
        t1=temp1
        t2=temp2
        t3=temp3
        t4=temp4
    else
        call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    end if
    if (t<=t1) then
        prhot=(x_ratio*rho/mh2+(1-x_ratio)*rho/mhe)*kb*t
    elseif (t>=t4) then
        prhot=(x_ratio*rho/mh*2d0+(1-x_ratio)*rho/mhe)*kb*t
    elseif (t2<t3) then
        if (t<=t2) then
            call calspecies_H2_HI(rho,t,species)
            ntot=sum(species)
            prhot=ntot*kb*t
        elseif (t>=t3) then
            call calspecies_HI_HII(rho,t,species)
            ntot=sum(species)
            prhot=ntot*kb*t
        else
            prhot=(x_ratio*rho/mh+(1-x_ratio)*rho/mhe)*kb*t
        end if
    elseif (t2>t3) then
        if (t<=t3) then
            call calspecies_H2_HI(rho,t,species)
            ntot=sum(species)
            prhot=ntot*kb*t
        elseif (t>=t2) then
            call calspecies_HI_HII(rho,t,species)
            ntot=sum(species)
            prhot=ntot*kb*t
        else
            call calspecies_H2_HI_HII(rho,t,species)
            ntot=sum(species)
            prhot=ntot*kb*t
        end if
    elseif (t2==t3) then
        if(t<=t2) then
            call calspecies_H2_HI(rho,t,species)
            ntot=sum(species)
            prhot=ntot*kb*t
        else
            call calspecies_HI_HII(rho,t,species)
            ntot=sum(species)
            prhot=ntot*kb*t
        end if
    end if
end function prhot

function fun_Tp(x)
    !x(1)=t,x(2)=rho,x(3)=p
    real(8), dimension(:), allocatable :: x
    real(8) :: fun_Tp,t,rho,p,ntot,t1,t2,t3,t4
    t=x(1)
    rho=x(2)
    p=x(3)
    t1=x(4)
    t2=x(5)
    t3=x(6)
    t4=x(7)
    fun_Tp=p-prhot(rho,t,t1,t2,t3,t4)
end function fun_Tp

function solvets(s,rho)
    !given specific entropy and rho, find the corresponding temperature
    real(8) :: solvets,s,rho,lowT,highT,lows,highs,mhion
    real(8) :: convergence,delta,nhetot,nhtot,ntot,t1,t2,t3,t4,ah2,ah,ae,ahe,bh2,bh,be,bhe
    procedure(fun), pointer :: ptr
    real(8), dimension(:), allocatable :: x,root
    integer :: conv_mode
    conv_mode=2
    call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    lowT=t1
    highT=t4
    convergence=1e-5
    lows=srhot(rho,lowT,t1,t2,t3,t4)
    highs=srhot(rho,highT,t1,t2,t3,t4)
    if (s<=lows) then  !purely molecular
        nhtot=x_ratio*rho/mh2
        nhetot=(1-x_ratio)*rho/mhe
        ah2=kb*nhtot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mh2*kb/h_planck**2)-log(nhtot))
        ahe=kb*nhetot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mhe*kb/h_planck**2)-log(nhetot))
        bh2=kb*nhtot/rho*1.5d0
        bhe=kb*nhetot/rho*1.5d0
        solvets=exp((s-ah2-ahe)/(bh2+bhe))
    else if (s>=highs) then  !fully ionized
        nhtot=x_ratio*rho/mh
        nhetot=(1-x_ratio)*rho/mhe
        mhion=mh-me
        ah=kb*nhtot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mhion*kb/h_planck**2)-log(nhtot))
        ae=kb*nhtot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*me*kb/h_planck**2)-log(nhtot))
        ahe=kb*nhetot/rho*(2.5d0+log(two)+1.5d0*log(2d0*pi*mhe*kb/h_planck**2)-log(nhetot))
        bh=kb*nhtot/rho*1.5d0
        be=kb*nhtot/rho*1.5d0
        bhe=kb*nhetot/rho*1.5d0
        solvets=exp((s-ah-ae-ahe)/(bh+be+bhe))
    else if (s<highs.and.s>lows) then
        allocate(x(7),root(7))
        x=(/highT/2d0,rho,s,t1,t2,t3,t4/)
        ptr=>fun_Ts
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solvets=root(1)
        deallocate(x,root)
    end if
end function solvets

function srhot(rho,t,temp1,temp2,temp3,temp4)
    real(8) :: srhot,rho,t,t1,t2,t3,t4
    real(8), optional :: temp1,temp2,temp3,temp4
    real(8), dimension(5) :: species,dlnzdlnt,d2lnzdlnt,mu,si
    if (present(temp1)) then
        t1=temp1
        t2=temp2
        t3=temp3
        t4=temp4
    else
        call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    end if
    call calspecies(rho,t,species,t1,t2,t3,t4)
    call caldlnz(t,dlnzdlnt,d2lnzdlnt)
    call calmu(rho,t,species,mu,t1,t2,t3,t4)
    si=1d0+dlnzdlnt-mu/kb/t
    si(1)=si(1)-log(4d0)
    si(2)=si(2)-log(4d0)
    si(3)=si(3)-log(2d0)
    si(4)=si(4)-log(2d0)
    si(5)=si(5)-log(4d0)
    srhot=kb/rho*dot_product(species,1d0+dlnzdlnt-mu/kb/t)
end function srhot

function srhot_per_particle(rho,t,temp1,temp2,temp3,temp4)
    real(8) :: srhot,rho,t,srhot_per_particle,n,t1,t2,t3,t4
    real(8), optional :: temp1,temp2,temp3,temp4
    real(8), dimension(5) :: species,dlnzdlnt,d2lnzdlnt,mu
    if (present(temp1)) then
        t1=temp1
        t2=temp2
        t3=temp3
        t4=temp4
    else
        call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    end if
    call calspecies(rho,t,species,t1,t2,t3,t4)
    call caldlnz(t,dlnzdlnt,d2lnzdlnt)
    call calmu(rho,t,species,mu,t1,t2,t3,t4)
    srhot=kb/rho*dot_product(species,1d0+dlnzdlnt-mu/kb/t)
    n=sum(species)
    srhot_per_particle=srhot/kb/n*rho
end function srhot_per_particle

function fun_Ts(x)
    !x(1)=t,x(2)=rho,x(3)=s,x(4-7)=t1,t2,t3,t4
    real(8), dimension(:), allocatable :: x,species
    real(8) :: rho,t,s,fun_Ts,t1,t2,t3,t4
    t=x(1)
    rho=x(2)
    s=x(3)
    t1=x(4)
    t2=x(5)
    t3=x(6)
    t4=x(7)
    fun_Ts=s-srhot(rho,t,t1,t2,t3,t4)
end function fun_Ts

function solvetegv(egv,rho)
    !given egv and rho, find the corresponding temperature
    real(8) :: solvetegv,egv,rho,lowT,highT,lowegv,highegv,t_temp,t
    real(8) :: convergence,delta,nhetot,nhtot,ntot,t1,t2,t3,t4
    procedure(fun), pointer :: ptr
    real(8), dimension(:), allocatable :: x,root
    integer :: conv_mode
    conv_mode=2
    convergence=1e-5
    call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    lowT=t1
    highT=t4
    lowegv=egvrhot(rho,lowT,t1,t2,t3,t4)
    highegv=egvrhot(rho,highT,t1,t2,t3,t4)
    if (egv<=lowegv) then  !purely molecular
        nhtot=x_ratio*rho/mh2
        nhetot=(1-x_ratio)*rho/mhe
        ntot=nhtot+nhetot
        solvetegv=2d0/3d0*egv/kb/ntot
    else if (egv>=highegv) then  !fully ionized
        nhtot=x_ratio*rho/mh
        nhetot=(1-x_ratio)*rho/mhe
        ntot=2*nhtot+nhetot
        solvetegv=2d0/3d0*(egv-nhtot*ionh-0.5d0*nhtot*dish)/kb/ntot
    else if (egv<highegv.and.egv>lowegv) then
        allocate(x(7),root(7))
        x=(/highT/2d0,rho,egv,t1,t2,t3,t4/)
        ptr=>fun_Tegv
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solvetegv=root(1)
        deallocate(x,root)
    end if
end function solvetegv

function egvrhot(rho,t,temp1,temp2,temp3,temp4)
    !calculate egv with given rho and t, assume rho-t division is known
    real(8) :: egvrhot,rho,t,t1,t2,t3,t4,energy(5)
    real(8), optional :: temp1,temp2,temp3,temp4
    real(8), dimension(5) :: species
    energy(1)=1.5d0*kb*t                 !H2
    energy(2)=1.5d0*kb*t+dish/2d0        !HI
    energy(3)=1.5d0*kb*t+dish/2d0+ionh   !HII
    energy(4)=1.5d0*kb*t                 !electron
    energy(5)=1.5d0*kb*t                 !He
    if (present(temp1)) then
        t1=temp1
        t2=temp2
        t3=temp3
        t4=temp4
    else
        call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    end if
    if (t<=t1) then
        egvrhot=x_ratio*rho/mh2*energy(1)+(1-x_ratio)*rho/mhe*energy(5)
    elseif (t>=t4) then
        egvrhot=x_ratio*rho/mh*(energy(3)+energy(4))+(1-x_ratio)*rho/mhe*energy(5)
    elseif (t2<t3) then
        if (t<=t2) then
            call calspecies_H2_HI(rho,t,species)
            egvrhot=dot_product(species,energy)
        elseif (t>=t3) then
            call calspecies_HI_HII(rho,t,species)
            egvrhot=dot_product(species,energy)
        else
            !egvrhot=rho/mh*energy_h
            egvrhot=x_ratio*rho/mh*energy(2)+(1-x_ratio)*rho/mhe*energy(5)
        end if
    elseif (t2>t3) then
        if (t<=t3) then
            call calspecies_H2_HI(rho,t,species)
            egvrhot=dot_product(species,energy)
        elseif (t>=t2) then
            call calspecies_HI_HII(rho,t,species)
            egvrhot=dot_product(species,energy)
        else
            call calspecies_H2_HI_HII(rho,t,species)
            egvrhot=dot_product(species,energy)
        end if
    elseif (t2==t3) then
        if(t<=t2) then
            call calspecies_H2_HI(rho,t,species)
            egvrhot=dot_product(species,energy)
        else
            call calspecies_HI_HII(rho,t,species)
            egvrhot=dot_product(species,energy)
        end if
    end if
end function egvrhot

function fun_Tegv(x)
    !x(1)=t,x(2)=rho,x(3)=egv,x(4-7)=t1,t2,t3,t4
    !calculate the root function
    real(8), dimension(:), allocatable :: x,species
    real(8) :: rho,t,egv,fun_Tegv,t1,t2,t3,t4
    t=x(1)
    rho=x(2)
    egv=x(3)
    t1=x(4)
    t2=x(5)
    t3=x(6)
    t4=x(7)
    fun_Tegv=egv-egvrhot(rho,t,t1,t2,t3,t4)
end function fun_Tegv

function solveth(h,rho)
    !given h=p+egv and rho, find the corresponding temperature
    real(8) :: solveth,h,rho,lowT,lowh,highT,highh,convergence,t,nhetot,nhtot,ntot,t1,t2,t3,t4
    real(8), dimension(:), allocatable :: x,root
    procedure(fun), pointer :: ptr1,ptr2,ptr
    integer :: conv_mode
    conv_mode=2
    call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    lowT=t1
    highT=t4
    convergence=1e-5
    lowh=hrhot(rho,lowT,t1,t2,t3,t4)
    highh=hrhot(rho,highT,t1,t2,t3,t4)
    if (h<=lowh) then
        nhtot=x_ratio*rho/mh2   !all molecular
        nhetot=(1-x_ratio)*rho/mhe
        ntot=nhtot+nhetot
        solveth=h/kb/ntot/2.5d0
    else if (h>=highh) then
        nhtot=x_ratio*rho/mh    !fully ionized
        nhetot=(1-x_ratio)*rho/mhe
        ntot=2*nhtot+nhetot
        solveth=(h-nhtot*ionh-0.5d0*nhtot*dish)/kb/ntot/2.5d0
    else if (h>lowh.and.h<highh) then
        allocate(x(7),root(7))
        x=(/highT/2d0,rho,h,t1,t2,t3,t4/)
        ptr=>fun_Th
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solveth=root(1)
        deallocate(x,root)
    end if
end function solveth

function fun_Th(x)
    !x(1)=t,x(2)=rho,x(3)=h,x(4-7)=t1,t2,t3,t4
    !calculate the root function
    real(8), dimension(:), allocatable :: x,species
    real(8) :: rho,t,h,fun_Th,t1,t2,t3,t4
    t=x(1)
    rho=x(2)
    h=x(3)
    t1=x(4)
    t2=x(5)
    t3=x(6)
    t4=x(7)
    fun_Th=h-hrhot(rho,t,t1,t2,t3,t4)
end function fun_Th

function hrhot(rho,t,temp1,temp2,temp3,temp4)
    real(8) :: hrhot,rho,t,p,egv
    real(8), optional :: temp1,temp2,temp3,temp4
    if (present(temp1)) then
        p=prhot(rho,t,temp1,temp2,temp3,temp4)
        egv=egvrhot(rho,t,temp1,temp2,temp3,temp4)
    else
        p=prhot(rho,t)
        egv=egvrhot(rho,t)
    end if
    hrhot=p+egv
end function hrhot

subroutine calsrhost(rho,t,species,dlnndlnrho,dlnndlnt,dlnzdlnt,d2lnzdlnt,srho,st)
    real(8) :: rho,t,srho,st
    real(8), dimension(5) :: dsdrho,dsdt,dndrho,dndt,mu,species,dlnndlnrho,dlnndlnt,dlnzdlnt,d2lnzdlnt
    integer :: i
    call calmu(rho,t,species,mu)
    do i=1,5
        dndrho(i)=species(i)/rho*dlnndlnrho(i)
        dndt(i)=species(i)/t*dlnndlnt(i)
        dsdt(i)=kb*species(i)/rho/t*(d2lnzdlnt(i)+(1+dlnndlnt(i))*dlnzdlnt(i))-dndt(i)*mu(i)/rho/t
        dsdrho(i)=kb*species(i)/rho**2*((dlnndlnrho(i)-1)*dlnzdlnt(i)-1)+(species(i)/rho-dndrho(i))*mu(i)/rho/t
    end do
    srho=sum(dsdrho)
    st=sum(dsdt)
end subroutine calsrhost

subroutine calprhopt(species,dlnndlnrho,dlnndlnt,rho,t,prho,pt)
    real(8) :: rho,t,prho,pt,ntot
    real(8), dimension(5) :: species,dlnndlnrho,dlnndlnt
    integer :: i
    ntot=sum(species)
    prho=0d0
    pt=0d0
    do i=1,5
        prho=prho+species(i)/ntot*dlnndlnrho(i)
        pt=pt+species(i)/ntot*dlnndlnt(i)
    end do
    pt=pt+1  !this is actually dlnpdlnt
end subroutine calprhopt

subroutine caldlnz(t,dlnzdlnt,d2lnzdlnt)
    real(8) :: t
    real(8), dimension(5) :: dlnzdlnt,d2lnzdlnt
    dlnzdlnt(1)=1.5d0
    dlnzdlnt(2)=1.5d0+dish/2d0/kb/t
    dlnzdlnt(3)=1.5d0+(dish/2d0+ionh)/kb/t
    dlnzdlnt(4)=1.5d0
    dlnzdlnt(5)=1.5d0
    d2lnzdlnt(1)=0d0
    d2lnzdlnt(2)=-dish/2d0/kb/t
    d2lnzdlnt(3)=-(dish/2d0+ionh)/kb/t
    d2lnzdlnt(4)=0d0
    d2lnzdlnt(5)=0d0
end subroutine caldlnz

subroutine calmu(rho,t,species,mu,temp1,temp2,temp3,temp4)
    real(8) :: rho,t,t1,t2,t3,t4
    real(8), optional :: temp1,temp2,temp3,temp4
    real(8), dimension(5) :: species,mu
    if (present(temp1)) then
        t1=temp1
        t2=temp2
        t3=temp3
        t4=temp4
    else
        call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    end if
    if (t<=t1) then
        mu(1)=kb*t*log(species(1)/zh2(t))
        mu(2)=0d0
        mu(3)=0d0
        mu(4)=0d0
    elseif (t>=t4) then
        mu(1)=0d0
        mu(2)=0d0
        mu(3)=kb*t*log(species(3)/zhion(t))
        mu(4)=kb*t*log(species(4)/ze(t))
    elseif (t2<t3) then
        if (t<=t2) then
            mu(1)=kb*t*log(species(1)/zh2(t))
            mu(2)=kb*t*log(species(2)/zh(t))
            mu(3)=0d0
            mu(4)=0d0
        elseif (t>=t3) then
            mu(1)=0d0
            mu(2)=kb*t*log(species(2)/zh(t))
            mu(3)=kb*t*log(species(3)/zhion(t))
            mu(4)=kb*t*log(species(4)/ze(t))
        else
            mu(1)=0d0
            mu(2)=kb*t*log(species(2)/zh(t))
            mu(3)=0d0
            mu(4)=0d0
        end if
    elseif (t2>t3) then
        if (t<=t3) then
            mu(1)=kb*t*log(species(1)/zh2(t))
            mu(2)=kb*t*log(species(2)/zh(t))
            mu(3)=0d0
            mu(4)=0d0
        elseif (t>=t2) then
            mu(1)=0d0
            mu(2)=kb*t*log(species(2)/zh(t))
            mu(3)=kb*t*log(species(3)/zhion(t))
            mu(4)=kb*t*log(species(4)/ze(t))
        else
            mu(1)=kb*t*log(species(1)/zh2(t))
            mu(2)=kb*t*log(species(2)/zh(t))
            mu(3)=kb*t*log(species(3)/zhion(t))
            mu(4)=kb*t*log(species(4)/ze(t))
        end if
    elseif (t2==t3) then
        if(t<=t2) then
            mu(1)=kb*t*log(species(1)/zh2(t))
            mu(2)=kb*t*log(species(2)/zh(t))
            mu(3)=0d0
            mu(4)=0d0
        else
            mu(1)=0d0
            mu(2)=kb*t*log(species(2)/zh(t))
            mu(3)=kb*t*log(species(3)/zhion(t))
            mu(4)=kb*t*log(species(4)/ze(t))
        end if
    end if
    if (x_ratio==1d0) then
        mu(5)=0d0
    else
        mu(5)=kb*t*log(species(5)/zhe(t))
    end if
end subroutine calmu

subroutine caldlnndlnrho(rho,t,species,dlnndlnrho,temp1,temp2,temp3,temp4)
    real(8) :: rho,t,rhs(2),a(3,3),b(3),x(3),ntot,t1,t2,t3,t4
    real(8), optional :: temp1,temp2,temp3,temp4
    real(8), dimension(5) :: species,dlnndlnrho
    if (present(temp1)) then
        t1=temp1
        t2=temp2
        t3=temp3
        t4=temp4
    else
        call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    end if
    if (t<=t1) then
        dlnndlnrho(1)=x_ratio*rho/mh2/species(1)
        dlnndlnrho(2)=0d0
        dlnndlnrho(3)=0d0
        dlnndlnrho(4)=0d0
    elseif (t>=t4) then
        dlnndlnrho(1)=0d0
        dlnndlnrho(2)=0d0
        dlnndlnrho(3)=x_ratio*rho/mh/species(3)
        dlnndlnrho(4)=dlnndlnrho(3)
    elseif (t2<t3) then
        if (t<=t2) then
            dlnndlnrho(3)=0d0
            dlnndlnrho(4)=0d0
            rhs(1)=x_ratio*rho/mh
            rhs(2)=0d0
            dlnndlnrho(2)=rhs(1)/(4d0*species(1)+species(2))
            dlnndlnrho(1)=2d0*dlnndlnrho(2)
        elseif (t>=t3) then
            dlnndlnrho(1)=0d0
            rhs(1)=x_ratio*rho/mh
            rhs(2)=0d0
            dlnndlnrho(3)=rhs(1)/(2*species(2)+species(3))
            dlnndlnrho(2)=2*dlnndlnrho(3)
            dlnndlnrho(4)=dlnndlnrho(3)
        else
            dlnndlnrho(1)=0d0
            dlnndlnrho(2)=x_ratio*rho/mh/species(2)
            dlnndlnrho(3)=0d0
            dlnndlnrho(4)=0d0
        end if
    elseif (t2>t3) then
        if (t<=t3) then
            dlnndlnrho(3)=0d0
            dlnndlnrho(4)=0d0
            rhs(1)=x_ratio*rho/mh
            rhs(2)=0d0
            dlnndlnrho(2)=rhs(1)/(4d0*species(1)+species(2))
            dlnndlnrho(1)=2d0*dlnndlnrho(2)
        elseif (t>=t2) then
            dlnndlnrho(1)=0d0
            rhs(1)=x_ratio*rho/mh
            rhs(2)=0d0
            dlnndlnrho(3)=rhs(1)/(2*species(2)+species(3))
            dlnndlnrho(2)=2*dlnndlnrho(3)
            dlnndlnrho(4)=dlnndlnrho(3)
        else
            ntot=x_ratio*rho/mh
            dlnndlnrho(3)=ntot/(8d0*species(1)+2d0*species(2)+species(3))
            dlnndlnrho(4)=dlnndlnrho(3)
            dlnndlnrho(2)=2d0*dlnndlnrho(3)
            dlnndlnrho(1)=2d0*dlnndlnrho(2)
        end if
    elseif (t2==t3) then
        if(t<=t2) then
            dlnndlnrho(3)=0d0
            dlnndlnrho(4)=0d0
            rhs(1)=x_ratio*rho/mh
            rhs(2)=0d0
            dlnndlnrho(2)=rhs(1)/(4d0*species(1)+species(2))
            dlnndlnrho(1)=2d0*dlnndlnrho(2)
        else
            dlnndlnrho(1)=0d0
            rhs(1)=x_ratio*rho/mh
            rhs(2)=0d0
            dlnndlnrho(3)=rhs(1)/(2*species(2)+species(3))
            dlnndlnrho(2)=2*dlnndlnrho(3)
            dlnndlnrho(4)=dlnndlnrho(3)
        end if
    end if
    if (x_ratio==1d0) then
        dlnndlnrho(5)=0d0
    else
        dlnndlnrho(5)=(1d0-x_ratio)*rho/mhe/species(5)
    end if
end subroutine caldlnndlnrho

subroutine caldndt(rho,t,species,dndt)
    real(8) :: rho,t
    real(8), dimension(5) :: species,dndt,dlnndlnt
    integer :: i
    call caldlnndlnt(rho,t,species,dlnndlnt)
    do i=1,5
        dndt(i)=species(i)/t*dlnndlnt(i)
    end do
end subroutine caldndt

subroutine caldlnndlnt(rho,t,species,dlnndlnt,temp1,temp2,temp3,temp4)
    real(8) :: rho,t,a(3,3),b(3),x(3),rhs(2),ph2,ph,phion,pe,t1,t2,t3,t4
    real(8), optional :: temp1,temp2,temp3,temp4
    real(8), dimension(5) :: species,dlnndlnt
    !partition of energy
    ph2=1.5d0
    ph=1.5d0+dish/kb/t/2d0
    phion=1.5d0+(dish/2d0+ionh)/kb/t
    pe=1.5d0
    if (present(temp1)) then
        t1=temp1
        t2=temp2
        t3=temp3
        t4=temp4
    else
        call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    end if
    if (t<=t1) then
        dlnndlnt=0d0
    elseif (t>=t4) then
        dlnndlnt=0d0
    elseif (t2<t3) then
        if (t<=t2) then
            dlnndlnt(3)=0d0
            dlnndlnt(4)=0d0
            rhs(1)=0d0
            rhs(2)=2d0*ph-ph2
            dlnndlnt(2)=2d0*rhs(2)*species(1)/(4d0*species(1)+species(2))
            dlnndlnt(1)=-species(2)*dlnndlnt(2)/2d0/species(1)
        elseif (t>=t3) then
            dlnndlnt(1)=0d0
            rhs(1)=0d0
            rhs(2)=phion+pe-ph
            dlnndlnt(3)=rhs(2)*species(2)/(2d0*species(2)+species(3))
            dlnndlnt(2)=-species(3)*dlnndlnt(3)/species(2)
            dlnndlnt(4)=dlnndlnt(3)
        else
            dlnndlnt=0d0
        end if
    elseif (t2>t3) then
        if (t<=t3) then
            dlnndlnt(3)=0d0
            dlnndlnt(4)=0d0
            rhs(1)=0d0
            rhs(2)=2d0*ph-ph2
            dlnndlnt(2)=2d0*rhs(2)*species(1)/(4d0*species(1)+species(2))
            dlnndlnt(1)=-species(2)*dlnndlnt(2)/2d0/species(1)
        elseif (t>=t2) then
            dlnndlnt(1)=0d0
            rhs(1)=0d0
            rhs(2)=phion+pe-ph
            dlnndlnt(3)=rhs(2)*species(2)/(2d0*species(2)+species(3))
            dlnndlnt(2)=-species(3)*dlnndlnt(3)/species(2)
            dlnndlnt(4)=dlnndlnt(3)
        else
            b=(/0d0,2d0*ph-ph2,phion+pe-ph/)
            dlnndlnt(3)=(2d0*species(1)*b(2)+(4d0*species(1)+species(2))*b(3))/(8d0*species(1)+2d0*species(2)+species(3))
            dlnndlnt(4)=dlnndlnt(3)
            dlnndlnt(2)=2d0*dlnndlnt(3)-b(3)
            dlnndlnt(1)=2d0*dlnndlnt(2)-b(2)
        end if
    elseif (t2==t3) then
        if(t<=t2) then
            dlnndlnt(3)=0d0
            dlnndlnt(4)=0d0
            rhs(1)=0d0
            rhs(2)=2d0*ph-ph2
            dlnndlnt(2)=2d0*rhs(2)*species(1)/(4d0*species(1)+species(2))
            dlnndlnt(1)=-species(2)*dlnndlnt(2)/2d0/species(1)
        else
            dlnndlnt(1)=0d0
            rhs(1)=0d0
            rhs(2)=phion+pe-ph
            dlnndlnt(3)=rhs(2)*species(2)/(2d0*species(2)+species(3))
            dlnndlnt(2)=-species(3)*dlnndlnt(3)/species(2)
            dlnndlnt(4)=dlnndlnt(3)
        end if
    end if
    dlnndlnt(5)=0d0
end subroutine caldlnndlnt

subroutine generate_temperature_division(rho,t,temp_division)
    !to get four rho-t relations that divide the rho-t EOS into 6 blocks
    !(H2), (H2,HI), (HI), (H2,HI,HII,e),(HI,HII,e), and (HII,e)
    real(8), dimension(:,:,:), allocatable :: xi,division
    real(8), dimension(:,:), allocatable :: temp_division
    real(8), dimension(:), allocatable :: rho,t
    real(8) :: log10_tmax,delta,xi_temp(4)
    integer :: i,j,k,ntemp,nrho
    nrho=size(rho)
    ntemp=size(t)
    allocate(division(ntemp,nrho,4),xi(ntemp,nrho,4))
    do i=1,nrho
        do j=1,ntemp
            call calxi_temp_division(rho(i),t(j),xi_temp)
            xi(j,i,:)=xi_temp
        end do
    end do
    log10_tmax=log10(t(ntemp))
    delta=1e-6
    division=-1d0
    where (xi(:,:,:)>delta.and.(1-xi(:,:,:))>delta)
        division(:,:,:)=1d0
    end where
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,2)<0.and.division(j+1,i,2)>0) then  !HI emerges
                temp_division(i,1)=log10(t(j))
                exit
            else
                temp_division(i,1)=log10_tmax
            end if
        end do
    end do
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,1)>0.and.division(j+1,i,1)<0) then  !H2 disappears
                temp_division(i,2)=log10(t(j))
                exit
            else
                temp_division(i,2)=log10_tmax
            end if
        end do
    end do
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,3)<0.and.division(j+1,i,3)>0) then  !HII emerges
                temp_division(i,3)=log10(t(j))
                exit
            else
                temp_division(i,3)=log10_tmax
            end if
        end do
    end do
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,3)>0.and.division(j+1,i,3)<0) then  !HI disappears
                temp_division(i,4)=log10(t(j))
                exit
            else
                temp_division(i,4)=log10_tmax
            end if
        end do
    end do
    division_known=.true.
    deallocate(division,xi)
end subroutine generate_temperature_division

subroutine calxi_temp_division(rho,t,xi)
    real(8) :: rho,t,xi(4)
    real(8), dimension(5) :: species
    call calspecies_division_unk(rho,t,species)
    xi(1)=species(1)/x_ratio/rho*mh2
    xi(2)=species(2)/x_ratio/rho*mh
    xi(3)=species(3)/x_ratio/rho*mh
    xi(4)=xi(3)
end subroutine calxi_temp_division

subroutine divide_temperature_domain_known(rho,t1,t2,t3,t4)
    !given rho, calculate the four temperatures
    real(8) :: rho,dlogrho,t1,t2,t3,t4
    t1=10d0**get_t1(rho*x_ratio)
    t2=10d0**get_t2(rho*x_ratio)
    t3=10d0**get_t3(rho*x_ratio)
    t4=10d0**get_t4(rho*x_ratio)
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
    integer :: i,n
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=log10(rho)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    get_t2=temp_division(i,2)+(temp_division(i+1,2)-temp_division(i,2))*dd/dlogrho
end function get_t2

function get_t3(rho)
    real(8) :: get_t3,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=log10(rho)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    get_t3=temp_division(i,3)+(temp_division(i+1,3)-temp_division(i,3))*dd/dlogrho
end function get_t3

function get_t4(rho)
    real(8) :: get_t4,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=log10(rho)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    get_t4=temp_division(i,4)+(temp_division(i+1,4)-temp_division(i,4))*dd/dlogrho
end function get_t4

subroutine calspecies_division_unk(rho,t,species)
    real(8) :: rho,t,t1,t2,t3,t4
    real(8), dimension(5) :: species
    character(len=128) :: alert
    call divide_temperature_domain(rho,t1,t2,t3,t4)
    if (t<=t1) then   !purely molecular
        species(1)=x_ratio*rho/mh2
        species(2)=0d0
        species(3)=0d0
        species(4)=0d0
    else if (t>=t4) then  !fully ionized
        species(1)=0d0
        species(2)=0d0
        species(3)=x_ratio*rho/mh
        species(4)=x_ratio*rho/mh
    else if (t>t1.and.t<=t2) then  !H2,HI
        species(3)=0d0
        species(4)=0d0
        species(2)=solvenh_quadratic(rho,t)
        species(1)=species(2)**2/kdis(t)
    else if (t>t2.and.t<t3) then  !four species
        species(3)=solvenhion_quartic(rho,t)
        species(4)=species(3)
        species(2)=species(3)**2/kion(t)
        species(1)=species(2)**2/kdis(t)
    else if (t>=t3.and.t<t4) then  !HI, HII, e
        species(1)=0d0
        species(3)=solvenhion_quadratic(rho,t)
        species(2)=species(3)**2/kion(t)
        species(4)=species(3)
    else
        print *, 'calspecies wrong region2'
        stop
    end if
    species(5)=(1-x_ratio)*rho/mhe
end subroutine calspecies_division_unk

subroutine calspecies_H2_HI(rho,t,species)  !only H2 and HI present
    real(8) :: rho,t
    real(8), dimension(5) :: species
    species(3)=0d0
    species(4)=0d0
    species(2)=solvenh_quadratic(rho,t)
    species(1)=species(2)**2/kdis(t)
    species(5)=(1d0-x_ratio)*rho/mhe
end subroutine calspecies_H2_HI

subroutine calspecies_HI_HII(rho,t,species)  !only HI and HII present
    real(8) :: rho,t
    real(8), dimension(5) :: species
    species(1)=0d0
    species(3)=solvenhion_quadratic(rho,t)
    species(2)=species(3)**2/kion(t)
    species(4)=species(3)
    species(5)=(1d0-x_ratio)*rho/mhe
end subroutine calspecies_HI_HII

subroutine calspecies_H2_HI_HII(rho,t,species)   !four species
    real(8) :: rho,t
    real(8), dimension(5) :: species
    species(3)=solvenhion_quartic(rho,t)
    species(4)=species(3)
    species(2)=species(3)**2/kion(t)
    species(1)=species(2)**2/kdis(t)
    species(5)=(1d0-x_ratio)*rho/mhe
end subroutine calspecies_H2_HI_HII

subroutine divide_temperature_domain(rho,t1,t2,t3,t4)
    real(8) :: rho,t1,t2,t3,t4
    t1=temperature1(rho)
    t2=temperature2(rho)
    t3=temperature3(rho)
    t4=temperature4(rho)
end subroutine divide_temperature_domain

function temperature1(rho)
    real(8) :: rho,temperature1
    temperature1=10**2.7
end function temperature1

function temperature2(rho)
    real(8) :: rho,temperature2,k
    k=0.032d0
    if (rho>1e-10) then
        temperature2=10**(3.2d0+k*10d0)
    else
        temperature2=10**(3.2d0+k*(log10(rho)+20d0))
    end if
end function temperature2

function temperature3(rho)
    real(8) :: rho,temperature3,k1,k2
    k1=0.038d0
    k2=1.25d0
    if (rho>1d-10) then
        temperature3=min(10**(3.2d0+k1*10d0+k2*(log10(rho)+10d0)),1.1e6)
    else if (rho<=1d-10) then
        temperature3=10**(3.2d0+k1*(log10(rho)+20d0))
    end if
end function temperature3

function temperature4(rho)
    real(8) :: rho,temperature4,k
    k=2d0/10d0
    temperature4=min(10**(4d0+k*(log10(rho)+20d0)),1.1e6)
end function temperature4

subroutine calspecies(rho,t,species,temp1,temp2,temp3,temp4)
    real(8) :: rho,t,t1,t2,t3,t4
    real(8), optional :: temp1,temp2,temp3,temp4
    real(8), dimension(5) :: species
    if (present(temp1)) then
        t1=temp1
        t2=temp2
        t3=temp3
        t4=temp4
    else
        call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    end if
    if (t<=t1) then
        species=(/x_ratio*rho/mh2,0d0,0d0,0d0,(1-x_ratio)*rho/mhe/)
    else if (t>=t4) then
        species=(/0d0,0d0,x_ratio*rho/mh,x_ratio*rho/mh,(1-x_ratio)*rho/mhe/)
    else if (t2<t3) then
        if (t<=t2) then
            call calspecies_H2_HI(rho,t,species)
        elseif (t>=t3) then
            call calspecies_HI_HII(rho,t,species)
        else
            species=(/0d0,x_ratio*rho/mh,0d0,0d0,(1-x_ratio)*rho/mhe/)
        end if
    else if (t2>t3) then
        if (t<=t3) then
            call calspecies_H2_HI(rho,t,species)
        elseif (t>=t2) then
            call calspecies_HI_HII(rho,t,species)
        else
            call calspecies_H2_HI_HII(rho,t,species)
        end if
    else if (t2==t3) then
        if(t<=t2) then
            call calspecies_H2_HI(rho,t,species)
        else
            call calspecies_HI_HII(rho,t,species)
        end if
    end if
end subroutine calspecies

function solvenhion_quartic(rho,t)
    !solve the number density of HII at given rho and t, all species present
    real(8) :: rho,t,solvenhion_quartic,nhtot,a(5),b,c
    complex(8) :: z(4)
    b=kion(t)
    c=kdis(t)
    a(1)=-x_ratio*rho/mh
    a(2)=1d0
    a(3)=1d0/b
    a(4)=0d0
    a(5)=2d0/b/b/c
    call quarticroots(a,z)
    solvenhion_quartic=dble(z(1))
end function solvenhion_quartic

function solvenhion_quadratic(rho,t)
    !solve the number density of HII at given rho and t, only HI, HII, and e present
    real(8) :: rho,t,solvenhion_quadratic,nhtot,a
    a=kion(t)
    nhtot=x_ratio*rho/mh
    solvenhion_quadratic=2d0*a*nhtot/(a+sqrt(a**2+4d0*a*nhtot))
end function solvenhion_quadratic

function solvenh_quadratic(rho,t)
    !solve the number density of HII at given rho and t, only H2 and HI present
    real(8) :: rho,t,solvenh_quadratic,nhtot,a
    a=kdis(t)
    nhtot=x_ratio*rho/mh
    solvenh_quadratic=2d0*a*nhtot/(a+sqrt(a**2+8d0*a*nhtot))
end function solvenh_quadratic

function kion(t)
    real(8) :: kion,t
    kion=zhion(t)*ze(t)/zh(t)
end function kion

function kdis(t)
    real(8) :: kdis,t
    kdis=zh(t)*zh(t)/zh2(t)
end function kdis

function zh(t)
    real(8) :: zh,t,ztr,zspin,zelec
    zh=zsimplegas(t,mh)/(exp(dish/(2d0*kb*t)))
end function zh

function ze(t)
    real(8) :: ze,t
    ze=zsimplegas(t,me)
end function ze

function zhion(t)
    real(8) :: zhion,t,mhion
    zhion=zsimplegas(t,mh)/(exp((dish+2d0*ionh)/(2d0*kb*t)))
end function zhion

function zh2(t)
    real(8) :: zh2,t
    zh2=zsimplegas(t,mh2)
end function zh2

function zhe(t)
    real(8) :: zhe,t
    zhe=zsimplegas(t,mhe)
end function zhe

subroutine eos_species(blk)
    type(blockdef), pointer :: blk
    real(8) :: rho,temp,nhtot,species(5)
    integer :: i,j
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            rho=blk%w(1,i,j,1)
            temp=blk%temp(i,j,1)
            call calspecies(rho,temp,species)
            nhtot=2*species(1)+species(2)+species(3)
            blk%H2(i,j,1)=2*species(1)/nhtot
            blk%HI(i,j,1)=species(2)/nhtot
            blk%HII(i,j,1)=species(3)/nhtot
            blk%cs(i,j,1)=adiabatic_cs(rho,temp)
        end do
    end do
end subroutine eos_species

end module eos_h2_HI_HII

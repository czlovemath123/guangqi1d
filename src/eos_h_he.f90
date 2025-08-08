module eos_h_he
use mathlib
use phylib
implicit none

real(8), protected :: x_ratio,y_ratio
integer, protected :: ntemp_h_he,nrho_h_he

contains

subroutine initialize_x_ratio()
    !initialize the value of x_ratio in eos_analytic
    x_ratio=environment%h_ratio
    y_ratio=1d0-x_ratio
    nrho_h_he=environment%nrho
    ntemp_h_he=environment%ntemp
    if (rank==0) then
        print *, 'X=',x_ratio,'   Y=',one-x_ratio
    end if
end subroutine initialize_x_ratio

subroutine gammarhot(rho,t,gamm)
    !calculate the Gamma with given rho and t
    real(8) :: rho,t,gamm,maw,p,ntot,srho,st,prho,pt,rho1,rho2
    real(8), dimension(4) :: s_h,s_he,dlnzdlnt_h,dlnzdlnt_he, &
        d2lnzdlnt_h,d2lnzdlnt_he,dlnndlnrho_h,dlnndlnrho_he,  &
        dlnndlnt_h,dlnndlnt_he,t_h,t_he
    rho1=rho*x_ratio
    rho2=rho*y_ratio
    call divide_domain_hydrogen(rho1,t_h)
    call divide_domain_helium(rho2,t_he)
    call calspecies_hydrogen(rho1,t,s_h,t_h)
    call calspecies_helium(rho2,t,s_he,t_he)
    call caldlnz_hydrogen(t,dlnzdlnt_h,d2lnzdlnt_h)
    call caldlnz_helium(t,dlnzdlnt_he,d2lnzdlnt_he)
    call caldlnndlnrho_hydrogen(rho1,t,s_h,dlnndlnrho_h,t_h)
    call caldlnndlnrho_helium(rho2,t,s_he,dlnndlnrho_he,t_he)
    call caldlnndlnt_hydrogen(rho1,t,s_h,dlnndlnt_h,t_h)
    call caldlnndlnt_helium(rho2,t,s_he,dlnndlnt_he,t_he)
    call calprhopt(s_h,s_he,dlnndlnrho_h,dlnndlnrho_he,dlnndlnt_h,dlnndlnt_he,prho,pt)
    call calsrhost(rho,t,s_h,s_he,dlnndlnrho_h,dlnndlnrho_he,dlnndlnt_h,dlnndlnt_he,      &
        dlnzdlnt_h,dlnzdlnt_he,d2lnzdlnt_h,d2lnzdlnt_he,srho,st)
    ntot=sum(s_h)+sum(s_he)
    p=ntot*kb*t
    gamm=rho/p*(p/rho*prho-p/t*pt*srho/st)
end subroutine gammarhot

function adiabatic_cs(rho,t)
    !calculate the adiabatic sound speed with given rho and t
    real(8) :: adiabatic_cs,rho,t,gamm,maw,p,ntot,srho,st,prho,pt,rho1,rho2
    real(8), dimension(4) :: s_h,s_he,dlnzdlnt_h,dlnzdlnt_he, &
        d2lnzdlnt_h,d2lnzdlnt_he,dlnndlnrho_h,dlnndlnrho_he,  &
        dlnndlnt_h,dlnndlnt_he,t_h,t_he
    rho1=rho*x_ratio
    rho2=rho*y_ratio
    call divide_domain_hydrogen(rho1,t_h)
    call divide_domain_helium(rho2,t_he)
    call calspecies_hydrogen(rho1,t,s_h,t_h)
    call calspecies_helium(rho2,t,s_he,t_he)
    call caldlnz_hydrogen(t,dlnzdlnt_h,d2lnzdlnt_h)
    call caldlnz_helium(t,dlnzdlnt_he,d2lnzdlnt_he)
    call caldlnndlnrho_hydrogen(rho1,t,s_h,dlnndlnrho_h,t_h)
    call caldlnndlnrho_helium(rho2,t,s_he,dlnndlnrho_he,t_he)
    call caldlnndlnt_hydrogen(rho1,t,s_h,dlnndlnt_h,t_h)
    call caldlnndlnt_helium(rho2,t,s_he,dlnndlnt_he,t_he)
    call calprhopt(s_h,s_he,dlnndlnrho_h,dlnndlnrho_he,dlnndlnt_h,dlnndlnt_he,prho,pt)
    call calsrhost(rho,t,s_h,s_he,dlnndlnrho_h,dlnndlnrho_he,dlnndlnt_h,dlnndlnt_he,      &
        dlnzdlnt_h,dlnzdlnt_he,d2lnzdlnt_h,d2lnzdlnt_he,srho,st)
    ntot=sum(s_h)+sum(s_he)
    p=ntot*kb*t
    adiabatic_cs=sqrt(p/rho*prho-p/t*pt*srho/st)
end function adiabatic_cs

function eos_cv(rho,t)
    real(8) :: eos_cv,rho,rho1,rho2,t,xi_h2,xi_hi,xi_hii,xi_hei,xi_heii,xi_heiii,xi_e,  &
        xi(7),s(7),dndt(7)
    real(8), dimension(4) :: s_h,s_he,dndt_h,dndt_he,t_h,t_he
    rho1=rho*x_ratio
    rho2=rho*y_ratio
    call divide_domain_hydrogen(rho1,t_h)
    call divide_domain_helium(rho2,t_he)
    call calspecies_hydrogen(rho1,t,s_h,t_h)
    call calspecies_helium(rho2,t,s_he,t_he)
    xi_h2=1.5d0*kb*t
    xi_hi=1.5d0*kb*t+dish/2d0
    xi_hii=1.5d0*kb*t+dish/2d0+ionh
    xi_hei=1.5d0*kb*t
    xi_heii=1.5d0*kb*t+ionhe1
    xi_heiii=1.5d0*kb*t+ionhe2
    xi_e=1.5d0*kb*t
    xi=(/xi_h2,xi_hi,xi_hii,xi_hei,xi_heii,xi_heiii,xi_e/)
    call caldndt_hydrogen(rho1,t,s_h,dndt_h,t_h)
    call caldndt_helium(rho2,t,s_he,dndt_he,t_he)
    s(1:3)=s_h(1:3)
    s(4:6)=s_he(1:3)
    s(7)=s_h(4)+s_he(4)
    dndt(1:3)=dndt_h(1:3)
    dndt(4:6)=dndt_he(1:3)
    dndt(7)=dndt_h(4)+dndt_he(4)
    eos_cv=dot_product(dndt,xi)+sum(s)*1.5d0*kb
end function eos_cv

function solvetp(p,rho)
    !given p and rho, find the corresponding temperature
    real(8) :: solvetp,p,rho,lowT,lowp,highT,highp,convergence,t,ntot,t_h(4),t_he(4),rho1,rho2
    real(8), dimension(:), allocatable :: x,root
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    conv_mode=2
    rho1=rho*x_ratio
    rho2=rho*y_ratio
    call divide_domain_hydrogen(rho1,t_h)
    call divide_domain_helium(rho2,t_he)
    lowT=t_h(1)
    highT=t_he(4)
    lowp=prhot(rho,lowT)
    highp=prhot(rho,highT)
    convergence=1d-5
    if (p<=lowp) then
        ntot=rho1/mh2+rho2/mhe
        solvetp=p/kb/ntot
    else if (p>=highp) then
        ntot=2d0*rho1/mh+3d0*rho2/mhe
        solvetp=p/kb/ntot
    else
        allocate(x(11),root(11))
        x=(/highT/2d0,rho,p,t_h(1:4),t_he(1:4)/)
        ptr=>fun_Tp_h_he
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solvetp=root(1)
        deallocate(x,root)
    end if
end function solvetp

function prhot(rho,t,t_h,t_he)
    real(8) :: prhot,rho,rho1,rho2,t,ntot,s_h(4),s_he(4)
    real(8), optional :: t_h(4),t_he(4)
    rho1=rho*x_ratio
    rho2=rho*y_ratio
    if (present(t_h)) then
        call calspecies_hydrogen(rho1,t,s_h,t_h)
        call calspecies_helium(rho2,t,s_he,t_he)
    else
        call calspecies_hydrogen(rho1,t,s_h)
        call calspecies_helium(rho2,t,s_he)
    end if
    ntot=sum(s_h)+sum(s_he)
    prhot=ntot*kb*t
end function prhot

function prhot_helium(rho,t,t_he)
    !calculate pressure with given rho and t
    real(8) :: prhot_helium,rho,t,ntot,species(4)
    real(8), optional :: t_he(4)
    if (present(t_he)) then
        call calspecies_helium(rho,t,species,t_he)
    else
        call calspecies_helium(rho,t,species)
    end if
    ntot=sum(species)
    prhot_helium=ntot*kb*t
end function prhot_helium

function prhot_hydrogen(rho,t,t_h)
    !calculate pressure with given rho and t
    real(8) :: prhot_hydrogen,rho,t,ntot,t1,t2,t3,t4,species(4)
    real(8), optional :: t_h(4)
    if (present(t_h)) then
        call calspecies_hydrogen(rho,t,species,t_h)
    else
        call calspecies_hydrogen(rho,t,species)
    end if
    ntot=sum(species)
    prhot_hydrogen=ntot*kb*t
end function prhot_hydrogen

function fun_Tp_h_he(x)
    real(8), dimension(:), allocatable :: x
    real(8) :: fun_Tp_h_he,t,rho,p,t_h(4),t_he(4)
    t=x(1)
    rho=x(2)
    p=x(3)
    t_h=x(4:7)
    t_he=x(8:11)
    fun_Tp_h_he=p-prhot(rho,t,t_h,t_he)
end function fun_Tp_h_he

function srhot(rho,t,t_h,t_he)
    real(8) :: srhot,rho,rho1,rho2,t,ntot
    real(8), dimension(4) :: dlnzdlnt_h,dlnzdlnt_he,d2lnzdlnt_h,d2lnzdlnt_he,mu_h,mu_he,s_h,s_he,si_h,si_he
    real(8), optional :: t_h(4),t_he(4)
    rho1=rho*x_ratio
    rho2=rho*y_ratio
    if (present(t_h)) then
        call calspecies_hydrogen(rho1,t,s_h,t_h)
        call calspecies_helium(rho2,t,s_he,t_he)
    else
        call calspecies_hydrogen(rho1,t,s_h)
        call calspecies_helium(rho2,t,s_he)
    end if
    call caldlnz_hydrogen(t,dlnzdlnt_h,d2lnzdlnt_h)
    call caldlnz_helium(t,dlnzdlnt_he,d2lnzdlnt_he)
    call calmu_hydrogen(rho1,t,s_h,mu_h,t_h)
    call calmu_helium(rho2,t,s_he,mu_he,t_he)
    si_h=1d0+dlnzdlnt_h-mu_h/kb/t
    si_he=1d0+dlnzdlnt_he-mu_he/kb/t
    srhot=kb/rho*(dot_product(s_h,si_h)+dot_product(s_he,si_he))
end function srhot

function srhot_per_particle(rho,t,t_h,t_he)
    real(8) :: srhot_per_particle,srhot,rho,rho1,rho2,t,ntot
    real(8), dimension(4) :: dlnzdlnt_h,dlnzdlnt_he,d2lnzdlnt_h,d2lnzdlnt_he,mu_h,mu_he,s_h,s_he,si_h,si_he
    real(8), optional :: t_h(4),t_he(4)
    rho1=rho*x_ratio
    rho2=rho*y_ratio
    if (present(t_h)) then
        call calspecies_hydrogen(rho1,t,s_h,t_h)
        call calspecies_helium(rho2,t,s_he,t_he)
    else
        call calspecies_hydrogen(rho1,t,s_h)
        call calspecies_helium(rho2,t,s_he)
    end if
    call caldlnz_hydrogen(t,dlnzdlnt_h,d2lnzdlnt_h)
    call caldlnz_helium(t,dlnzdlnt_he,d2lnzdlnt_he)
    call calmu_hydrogen(rho1,t,s_h,mu_h,t_h)
    call calmu_helium(rho2,t,s_he,mu_he,t_he)
    si_h=1d0+dlnzdlnt_h-mu_h/kb/t
    si_he=1d0+dlnzdlnt_he-mu_he/kb/t
    srhot=kb/rho*(dot_product(s_h,si_h)+dot_product(s_he,si_he))
    srhot_per_particle=srhot/kb/(sum(s_h)+sum(s_he))*rho
end function srhot_per_particle

function solvetegv(egv,rho)
    !given egv and rho, find the corresponding temperature
    real(8) :: solvetegv,egv,rho,lowT,lowegv,highT,highegv,convergence,t,ntot,nhtot,nhetot,t_h(4),t_he(4),rho1,rho2
    real(8), dimension(:), allocatable :: x,root
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    conv_mode=2
    rho1=rho*x_ratio
    rho2=rho*y_ratio
    call divide_domain_hydrogen(rho1,t_h)
    call divide_domain_helium(rho2,t_he)
    lowT=t_h(1)
    highT=t_he(4)
    lowegv=egvrhot(rho,lowT)
    highegv=egvrhot(rho,highT)
    convergence=1d-5
    if (egv<=lowegv) then
        ntot=rho1/mh2+rho2/mhe
        solvetegv=2d0/3d0*egv/kb/ntot
    else if (egv>=highegv) then
        nhtot=rho1/mh
        nhetot=rho2/mhe
        ntot=2d0*nhtot+3d0*nhetot
        solvetegv=2d0/3d0*(egv-nhtot*(ionh+0.5d0*dish)-(ionhe1+ionhe2)*nhetot)/kb/ntot
    else
        allocate(x(11),root(11))
        x=(/highT/2d0,rho,egv,t_h(1:4),t_he(1:4)/)
        ptr=>fun_Tegv_h_he
        call bisectionroot(ptr,lowT,highT,x,convergence,root,conv_mode)
        solvetegv=root(1)
        deallocate(x,root)
    end if
end function solvetegv

function egvrhot(rho,t,t_h,t_he)
    real(8) :: egvrhot,rho,rho1,rho2,t,ntot,s_h(4),s_he(4),egvrhot1,egvrhot2
    real(8), optional :: t_h(4),t_he(4)
    rho1=rho*x_ratio
    rho2=rho*y_ratio
    if (present(t_h)) then
        egvrhot1=egvrhot_hydrogen(rho1,t,t_h)
        egvrhot2=egvrhot_helium(rho2,t,t_he)
    else
        egvrhot1=egvrhot_hydrogen(rho1,t)
        egvrhot2=egvrhot_helium(rho2,t)
    end if
    egvrhot=egvrhot1+egvrhot2
end function egvrhot

function egvrhot_helium(rho,t,t_he)
    !calculate egv with given rho and t, assume rho-t division is known
    real(8) :: egvrhot_helium,rho,t,energy(4),species(4)
    real(8), optional :: t_he(4)
    energy(1)=1.5d0*kb*t
    energy(2)=1.5d0*kb*t+ionhe1
    energy(3)=1.5d0*kb*t+ionhe1+ionhe2
    energy(4)=1.5d0*kb*t
    if (present(t_he)) then
        call calspecies_helium(rho,t,species,t_he)
    else
        call calspecies_helium(rho,t,species)
    end if
    egvrhot_helium=dot_product(species,energy)
end function egvrhot_helium

function egvrhot_hydrogen(rho,t,t_h)
    !calculate egv with given rho and t, assume rho-t division is known
    real(8) :: egvrhot_hydrogen,rho,t,t1,t2,t3,t4,energy(4),species(4)
    real(8), optional :: t_h(4)
    energy(1)=1.5d0*kb*t
    energy(2)=1.5d0*kb*t+dish/2d0
    energy(3)=1.5d0*kb*t+dish/2d0+ionh
    energy(4)=1.5d0*kb*t
    if (present(t_h)) then
        call calspecies_hydrogen(rho,t,species,t_h)
    else
        call calspecies_hydrogen(rho,t,species)
    end if
    egvrhot_hydrogen=dot_product(species,energy)
end function egvrhot_hydrogen

function fun_Tegv_h_he(x)
    real(8), dimension(:), allocatable :: x
    real(8) :: fun_Tegv_h_he,t,rho,egv,t_h(4),t_he(4)
    t=x(1)
    rho=x(2)
    egv=x(3)
    t_h=x(4:7)
    t_he=x(8:11)
    fun_Tegv_h_he=egv-egvrhot(rho,t,t_h,t_he)
end function fun_Tegv_h_he

subroutine calsrhost(rho,t,s_h,s_he,dlnndlnrho_h,dlnndlnrho_he,dlnndlnt_h,dlnndlnt_he,      &
        dlnzdlnt_h,dlnzdlnt_he,d2lnzdlnt_h,d2lnzdlnt_he,srho,st)
    real(8) :: rho,t,srho,st,rho1,rho2
    real(8), dimension(4) :: s_h,s_he,dsdrho_h,dsdrho_he,dsdt_h,dsdt_he,dndrho_h,dndrho_he, &
        dndt_h,dndt_he,mu_h,mu_he,dlnndlnrho_h,dlnndlnrho_he,dlnndlnt_h,dlnndlnt_he,        &
        dlnzdlnt_h,dlnzdlnt_he,d2lnzdlnt_h,d2lnzdlnt_he
    integer :: i
    rho1=rho*x_ratio
    rho2=rho*y_ratio
    call calmu_hydrogen(rho1,t,s_h,mu_h)
    call calmu_helium(rho2,t,s_he,mu_he)
    do i=1,4
        dndrho_h(i)=s_h(i)/rho*dlnndlnrho_h(i)
        dndrho_he(i)=s_he(i)/rho*dlnndlnrho_he(i)
        dndt_h(i)=s_h(i)/t*dlnndlnt_h(i)
        dndt_he(i)=s_he(i)/t*dlnndlnt_he(i)
        dsdt_h(i)=kb*s_h(i)/rho/t*(d2lnzdlnt_h(i)+(1+dlnndlnt_h(i))*dlnzdlnt_h(i))-dndt_h(i)*mu_h(i)/rho/t
        dsdt_he(i)=kb*s_he(i)/rho/t*(d2lnzdlnt_he(i)+(1+dlnndlnt_he(i))*dlnzdlnt_he(i))-dndt_he(i)*mu_he(i)/rho/t
        dsdrho_h(i)=kb*s_h(i)/rho**2*((dlnndlnrho_h(i)-1)*dlnzdlnt_h(i)-1)+(s_h(i)/rho-dndrho_h(i))*mu_h(i)/rho/t
        dsdrho_he(i)=kb*s_he(i)/rho**2*((dlnndlnrho_he(i)-1)*dlnzdlnt_he(i)-1)+(s_he(i)/rho-dndrho_he(i))*mu_he(i)/rho/t
    end do
    srho=sum(dsdrho_h)+sum(dsdrho_he)
    st=sum(dsdt_h)+sum(dsdt_he)
end subroutine calsrhost

subroutine calprhopt(s_h,s_he,dlnndlnrho_h,dlnndlnrho_he,dlnndlnt_h,dlnndlnt_he,prho,pt)
    real(8) :: prho,pt,ntot
    real(8), dimension(4) :: s_h,s_he,dlnndlnrho_h,dlnndlnrho_he,dlnndlnt_h,dlnndlnt_he
    integer :: i
    ntot=sum(s_h)+sum(s_he)
    prho=0d0
    pt=0d0
    do i=1,4
        prho=prho+s_h(i)/ntot*dlnndlnrho_h(i)
        prho=prho+s_he(i)/ntot*dlnndlnrho_he(i)
        pt=pt+s_h(i)/ntot*dlnndlnt_h(i)
        pt=pt+s_he(i)/ntot*dlnndlnt_he(i)
    end do
    pt=pt+1
end subroutine calprhopt

subroutine caldlnz_helium(t,dlnzdlnt,d2lnzdlnt)
    real(8) :: t,dlnzdlnt(4),d2lnzdlnt(4)
    dlnzdlnt(1)=1.5d0
    dlnzdlnt(2)=1.5d0+ionhe1/kb/t
    dlnzdlnt(3)=1.5d0+(ionhe1+ionhe2)/kb/t
    dlnzdlnt(4)=1.5d0
    d2lnzdlnt(1)=0d0
    d2lnzdlnt(2)=-ionhe1/kb/t
    d2lnzdlnt(3)=-(ionhe1+ionhe2)/kb/t
    d2lnzdlnt(4)=0d0
end subroutine caldlnz_helium

subroutine caldlnz_hydrogen(t,dlnzdlnt,d2lnzdlnt)
    real(8) :: t,dlnzdlnt(4),d2lnzdlnt(4)
    dlnzdlnt(1)=1.5d0
    dlnzdlnt(2)=1.5d0+dish/2d0/kb/t
    dlnzdlnt(3)=1.5d0+(dish/2d0+ionh)/kb/t
    dlnzdlnt(4)=1.5d0
    d2lnzdlnt(1)=0d0
    d2lnzdlnt(2)=-dish/2d0/kb/t
    d2lnzdlnt(3)=-(dish/2d0+ionh)/kb/t
    d2lnzdlnt(4)=0d0
end subroutine caldlnz_hydrogen

subroutine calmu_helium(rho,t,species,mu,temp)
    real(8) :: rho,t,t1,t2,t3,t4,species(4),mu(4),t_he(4)
    real(8), optional :: temp(4)
    if (present(temp)) then
        t_he=temp
    else
        call divide_domain_helium(rho,t_he)
    end if
    t1=t_he(1)
    t2=t_he(2)
    t3=t_he(3)
    t4=t_he(4)
    if (t<=t1) then
        mu(1)=kb*t*log(species(1)/zhe(t))
        mu(2)=0d0
        mu(3)=0d0
        mu(4)=0d0
    else if (t>=t4) then
        mu(1)=0d0
        mu(2)=0d0
        mu(3)=kb*t*log(species(3)/zheion2(t))
        mu(4)=kb*t*log(species(4)/ze(t))
    else if (t2<=t3) then
        if (t<=t2) then
            mu(1)=kb*t*log(species(1)/zhe(t))
            mu(2)=kb*t*log(species(2)/zheion1(t))
            mu(3)=0d0
            mu(4)=kb*t*log(species(4)/ze(t))
        else if (t>=t3) then
            mu(1)=0d0
            mu(2)=kb*t*log(species(2)/zheion1(t))
            mu(3)=kb*t*log(species(3)/zheion2(t))
            mu(4)=kb*t*log(species(4)/ze(t))
        else
            mu(1)=0d0
            mu(2)=kb*t*log(species(2)/zheion1(t))
            mu(3)=0d0
            mu(4)=kb*t*log(species(4)/ze(t))
        end if
    else if (t2>t3) then
        if (t<=t3) then
            mu(1)=kb*t*log(species(1)/zhe(t))
            mu(2)=kb*t*log(species(2)/zheion1(t))
            mu(3)=0d0
            mu(4)=kb*t*log(species(4)/ze(t))
        else if (t>=t2) then
            mu(1)=0d0
            mu(2)=kb*t*log(species(2)/zheion1(t))
            mu(3)=kb*t*log(species(3)/zheion2(t))
            mu(4)=kb*t*log(species(4)/ze(t))
        else
            mu(1)=kb*t*log(species(1)/zhe(t))
            mu(2)=kb*t*log(species(2)/zheion1(t))
            mu(3)=kb*t*log(species(3)/zheion2(t))
            mu(4)=kb*t*log(species(4)/ze(t))
        end if
    end if
end subroutine calmu_helium

subroutine calmu_hydrogen(rho,t,species,mu,temp)
    real(8) :: rho,t,t1,t2,t3,t4,species(4),mu(4),t_h(4)
    real(8), optional :: temp(4)
    if (present(temp)) then
        t_h=temp
    else
        call divide_domain_hydrogen(rho,t_h)
    end if
    t1=t_h(1)
    t2=t_h(2)
    t3=t_h(3)
    t4=t_h(4)
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
end subroutine calmu_hydrogen

subroutine caldndt_helium(rho,t,species,dndt,temp)
    real(8) :: rho,t,species(4),dndt(4),dlnndlnt(4),t_he(4)
    integer :: i
    real(8), optional :: temp(4)
    if (present(temp)) then
        t_he=temp
    else
        call divide_domain_helium(rho,t_he)
    end if
    call caldlnndlnt_helium(rho,t,species,dlnndlnt,t_he)
    do i=1,4
        dndt(i)=species(i)/t*dlnndlnt(i)
    end do
end subroutine caldndt_helium

subroutine caldndt_hydrogen(rho,t,species,dndt,temp)
    real(8) :: rho,t,species(4),dndt(4),dlnndlnt(4),t_h(4)
    integer :: i
    real(8), optional :: temp(4)
    if (present(temp)) then
        t_h=temp
    else
        call divide_domain_hydrogen(rho,t_h)
    end if
    call caldlnndlnt_hydrogen(rho,t,species,dlnndlnt,t_h)
    do i=1,4
        dndt(i)=species(i)/t*dlnndlnt(i)
    end do
end subroutine caldndt_hydrogen

subroutine caldlnndlnrho_helium(rho,t,species,dlnndlnrho,temp)
    real(8) :: rho,t,x,ntot,t_he(4),species(4),dlnndlnrho(4),t1,t2,t3,t4
    real(8), optional :: temp(4)
    if (present(temp)) then
        t_he=temp
    else
        call divide_domain_helium(rho,t_he)
    end if
    t1=t_he(1)
    t2=t_he(2)
    t3=t_he(3)
    t4=t_he(4)
    ntot=rho/mhe
    if (t<=t1) then
        dlnndlnrho(1)=1d0
        dlnndlnrho(2)=0d0
        dlnndlnrho(3)=0d0
        dlnndlnrho(4)=0d0
    else if (t>=t4) then
        dlnndlnrho(1)=0d0
        dlnndlnrho(2)=0d0
        dlnndlnrho(3)=1d0
        dlnndlnrho(4)=1d0
    else if (t2<=t3) then
        if (t<=t2) then
            dlnndlnrho(2)=ntot/(2d0*species(1)+species(2))
            dlnndlnrho(1)=2d0*dlnndlnrho(2)
            dlnndlnrho(4)=dlnndlnrho(2)
            dlnndlnrho(3)=0d0
        else if (t>=t3) then
            x=species(4)*species(3)+species(4)*species(2)+species(2)*species(3)
            dlnndlnrho(1)=0d0
            dlnndlnrho(2)=ntot*(species(4)+2*species(3))/x
            dlnndlnrho(3)=dlnndlnrho(2)/(species(4)+2*species(3))*2*species(3)
            dlnndlnrho(4)=dlnndlnrho(2)/(species(4)+2*species(3))*species(4)
        else
            dlnndlnrho(1)=0d0
            dlnndlnrho(2)=ntot/species(2)
            dlnndlnrho(3)=0d0
            dlnndlnrho(4)=dlnndlnrho(2)
        end if
    else if (t2>t3) then
        if (t<=t3) then
            dlnndlnrho(2)=ntot/(2d0*species(1)+species(2))
            dlnndlnrho(1)=2d0*dlnndlnrho(2)
            dlnndlnrho(4)=dlnndlnrho(2)
            dlnndlnrho(3)=0d0
        else if (t>=t2) then
            x=species(4)*species(3)+species(4)*species(2)+species(2)*species(3)
            dlnndlnrho(1)=0d0
            dlnndlnrho(2)=ntot*(species(4)+2*species(3))/x
            dlnndlnrho(3)=dlnndlnrho(2)/(species(4)+2*species(3))*2*species(3)
            dlnndlnrho(4)=dlnndlnrho(2)/(species(4)+2*species(3))*species(4)
        else
            x=species(3)*species(2)+species(4)*ntot+species(1)*(4*species(3)+species(2))
            dlnndlnrho(1)=ntot*(species(2)+4*species(3)+species(4))/x
            dlnndlnrho(2)=ntot*(species(4)+2*species(3))/x
            dlnndlnrho(3)=ntot*2*species(3)/x
            dlnndlnrho(4)=ntot*species(4)/x
        end if
    end if
end subroutine caldlnndlnrho_helium

subroutine caldlnndlnrho_hydrogen(rho,t,species,dlnndlnrho,temp)
    real(8) :: rho,t,rhs(2),ntot,t_h(4),t1,t2,t3,t4,species(4),dlnndlnrho(4)
    real(8), optional :: temp(4)
    if (present(temp)) then
        t_h=temp
    else
        call divide_domain_hydrogen(rho,t_h)
    end if
    t1=t_h(1)
    t2=t_h(2)
    t3=t_h(3)
    t4=t_h(4)
    if (t<=t1) then
        dlnndlnrho(1)=1d0
        dlnndlnrho(2)=0d0
        dlnndlnrho(3)=0d0
        dlnndlnrho(4)=0d0
    elseif (t>=t4) then
        dlnndlnrho(1)=0d0
        dlnndlnrho(2)=0d0
        dlnndlnrho(3)=1d0
        dlnndlnrho(4)=1d0
    elseif (t2<=t3) then
        if (t<=t2) then
            dlnndlnrho(3)=0d0
            dlnndlnrho(4)=0d0
            rhs(1)=rho/mh
            rhs(2)=0d0
            dlnndlnrho(2)=rhs(1)/(4d0*species(1)+species(2))
            dlnndlnrho(1)=2d0*dlnndlnrho(2)
        elseif (t>=t3) then
            dlnndlnrho(1)=0d0
            rhs(1)=rho/mh
            rhs(2)=0d0
            dlnndlnrho(3)=rhs(1)/(2*species(2)+species(3))
            dlnndlnrho(2)=2*dlnndlnrho(3)
            dlnndlnrho(4)=dlnndlnrho(3)
        else
            dlnndlnrho(1)=0d0
            dlnndlnrho(2)=1d0
            dlnndlnrho(3)=0d0
            dlnndlnrho(4)=0d0
        end if
    elseif (t2>t3) then
        if (t<=t3) then
            dlnndlnrho(3)=0d0
            dlnndlnrho(4)=0d0
            rhs(1)=rho/mh
            rhs(2)=0d0
            dlnndlnrho(2)=rhs(1)/(4d0*species(1)+species(2))
            dlnndlnrho(1)=2d0*dlnndlnrho(2)
        elseif (t>=t2) then
            dlnndlnrho(1)=0d0
            rhs(1)=rho/mh
            rhs(2)=0d0
            dlnndlnrho(3)=rhs(1)/(2*species(2)+species(3))
            dlnndlnrho(2)=2*dlnndlnrho(3)
            dlnndlnrho(4)=dlnndlnrho(3)
        else
            ntot=rho/mh
            dlnndlnrho(3)=ntot/(8d0*species(1)+2d0*species(2)+species(3))
            dlnndlnrho(4)=dlnndlnrho(3)
            dlnndlnrho(2)=2d0*dlnndlnrho(3)
            dlnndlnrho(1)=2d0*dlnndlnrho(2)
        end if
    end if
end subroutine caldlnndlnrho_hydrogen

subroutine caldlnndlnt_helium(rho,t,species,dlnndlnt,temp)
    real(8) :: rho,t,x1,x2,y,pheI,pheII,pheIII,pe,ntot,t1,t2,t3,t4,species(4),dlnndlnt(4),t_he(4)
    real(8), optional :: temp(4)
    !partition of energy
    pheI=1.5d0
    pheII=1.5d0+ionhe1/kb/t
    pheIII=1.5d0+(ionhe1+ionhe2)/kb/t
    pe=1.5d0
    if (present(temp)) then
        t_he=temp
    else
        call divide_domain_helium(rho,t_he)
    end if
    t1=t_he(1)
    t2=t_he(2)
    t3=t_he(3)
    t4=t_he(4)
    if (t<=t1) then
        dlnndlnt=0d0
    elseif (t>=t4) then
        dlnndlnt=0d0
    elseif (t2<=t3) then
        if (t<=t2) then
            x1=pheII+pe-pheI
            dlnndlnt(2)=x1*species(1)/(2d0*species(1)+species(2))
            dlnndlnt(1)=-x1*species(2)/(2d0*species(1)+species(2))
            dlnndlnt(3)=0d0
            dlnndlnt(4)=dlnndlnt(2)
        elseif (t>=t3) then
            x1=pheIII+pe-pheII
            y=species(4)*species(3)+species(4)*species(2)+species(2)*species(3)
            dlnndlnt(1)=0d0
            dlnndlnt(2)=-x1*species(4)*species(3)/y
            dlnndlnt(3)=x1*species(2)*species(4)/y
            dlnndlnt(4)=x1*species(2)*species(3)/y
        else
            dlnndlnt=0d0
        end if
    elseif (t2>t3) then
        if (t<=t3) then
            x1=pheII+pe-pheI
            dlnndlnt(2)=x1*species(1)/(2d0*species(1)+species(2))
            dlnndlnt(1)=-x1*species(2)/(2d0*species(1)+species(2))
            dlnndlnt(3)=0d0
            dlnndlnt(4)=dlnndlnt(2)
        elseif (t>=t2) then
            x1=pheIII+pe-pheII
            y=species(4)*species(3)+species(4)*species(2)+species(2)*species(3)
            dlnndlnt(1)=0d0
            dlnndlnt(2)=-x1*species(4)*species(3)/y
            dlnndlnt(3)=x1*species(2)*species(4)/y
            dlnndlnt(4)=x1*species(2)*species(3)/y
        else
            x1=pheII+pe-pheI
            x2=pheIII+pe-pheII
            ntot=rho/mhe
            y=species(3)*species(2)+species(4)*ntot+species(1)*(4*species(3)+species(2))
            dlnndlnt(1)=((-x1+x2)*species(2)*species(3)-species(4)*((x1+x2)*species(3)+x1*species(2)))/y
            dlnndlnt(2)=(2*(x1-x2)*species(1)*species(3)+species(4)*(x1*species(1)-x2*species(3)))/y
            dlnndlnt(3)=((-x1+x2)*species(1)*species(2)+species(4)*((x1+x2)*species(1)+x2*species(2)))/y
            dlnndlnt(4)=(x2*species(3)*species(2)+species(1)*(2*(x1+x2)*species(3)+x1*species(2)))/y
        end if
    end if
end subroutine caldlnndlnt_helium

subroutine caldlnndlnt_hydrogen(rho,t,species,dlnndlnt,temp)
    real(8) :: rho,t,rhs(2),b(3),ph2,ph,phion,pe,t1,t2,t3,t4,species(4),dlnndlnt(4),t_h(4)
    real(8), optional :: temp(4)
    !partition of energy
    ph2=1.5d0
    ph=1.5d0+dish/kb/t/2d0
    phion=1.5d0+(dish/2d0+ionh)/kb/t
    pe=1.5d0
    if (present(temp)) then
        t_h=temp
    else
        call divide_domain_hydrogen(rho,t_h)
    end if
    t1=t_h(1)
    t2=t_h(2)
    t3=t_h(3)
    t4=t_h(4)
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
end subroutine caldlnndlnt_hydrogen

subroutine generate_temperature_division_helium(rho,t,temp_division)
    !to get four rho-t relations that divide the rho-t EOS into 6 blocks
    !(HeI), (HeI,HeII), (HeII), (HeI,HeII,HeIII,e),(HeII,HeIII,e), and (HeIII,e)
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
            call calxi_temp_division_helium(rho(i),t(j),xi_temp) 
            xi(j,i,:)=xi_temp
        end do
    end do
    log10_tmax=log10(t(ntemp))
    delta=1e-5
    division=-1d0
    where (xi(:,:,:)>delta.and.(1-xi(:,:,:))>delta)
        division(:,:,:)=1d0
    end where
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,2)<0.and.division(j+1,i,2)>0) then  !HeII emerges
                temp_division(i,1)=log10(t(j))
                exit
            else
                temp_division(i,1)=log10_tmax
            end if
        end do
    end do
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,1)>0.and.division(j+1,i,1)<0) then  !HeI disappears
                temp_division(i,2)=log10(t(j))
                exit
            else
                temp_division(i,2)=log10_tmax
            end if
        end do
    end do
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,3)<0.and.division(j+1,i,3)>0) then  !HeIII emerges
                temp_division(i,3)=log10(t(j))
                exit
            else
                temp_division(i,3)=log10_tmax
            end if
        end do
    end do
    do i=1,nrho
        do j=1,ntemp-1
            if (division(j,i,3)>0.and.division(j+1,i,3)<0) then  !HeII disappears
                temp_division(i,4)=log10(t(j))
                exit
            else
                temp_division(i,4)=log10_tmax
            end if
        end do
    end do
    division_helium=.true.
    deallocate(division,xi)
end subroutine generate_temperature_division_helium

subroutine generate_temperature_division_hydrogen(rho,t,temp_division)
    !to get four rho-t relations that divide the rho-t EOS into 6 blocks
    !(H2), (H2,HI), (HI), (H2,HI,HII,e),(HI,HII,e), and (HII,e)
    real(8), dimension(:,:,:), allocatable :: xi,division
    real(8), dimension(:,:), allocatable :: gamm,eg,maw,p,temp_division
    real(8), dimension(:), allocatable :: rho,t
    real(8) :: log10_tmax,delta,xi_temp(4)
    integer :: i,j,k,ntemp,nrho
    nrho=size(rho)
    ntemp=size(t)
    allocate(division(ntemp,nrho,4),xi(ntemp,nrho,4))
    do i=1,nrho
        do j=1,ntemp
            call calxi_temp_division_hydrogen(rho(i),t(j),xi_temp) 
            xi(j,i,:)=xi_temp
        end do
    end do
    log10_tmax=log10(t(ntemp))
    delta=1e-5
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
    division_hydrogen=.true.
    deallocate(division,xi)
end subroutine generate_temperature_division_hydrogen

subroutine calxi_temp_division_helium(rho,t,xi)
    real(8) :: rho,t,xi(3),species(4)
    call calspecies_division_helium_unk(rho,t,species)
    xi(1)=max(species(1)/rho*mhe,1d-16)
    xi(2)=max(species(2)/rho*mhe,1d-16)
    xi(3)=max(species(3)/rho*mhe,1d-16)
end subroutine calxi_temp_division_helium

subroutine calxi_temp_division_hydrogen(rho,t,xi)
    real(8) :: rho,t,xi(3),species(4)
    call calspecies_division_hydrogen_unk(rho,t,species)
    xi(1)=max(species(1)/rho*mh2,1d-16)
    xi(2)=max(species(2)/rho*mh,1d-16)
    xi(3)=max(species(3)/rho*mh,1d-16)
end subroutine calxi_temp_division_hydrogen

subroutine divide_domain_helium(rho,t)
    real(8) :: rho,dlogrho,t(4)
    t(1)=10**get_helium_t1(rho)
    t(2)=10**get_helium_t2(rho)
    t(3)=10**get_helium_t3(rho)
    t(4)=10**get_helium_t4(rho)
end subroutine divide_domain_helium

subroutine divide_domain_hydrogen(rho,t)
    real(8) :: rho,dlogrho,t(4)
    t(1)=10**get_hydrogen_t1(rho)
    t(2)=10**get_hydrogen_t2(rho)
    t(3)=10**get_hydrogen_t3(rho)
    t(4)=10**get_hydrogen_t4(rho)
end subroutine divide_domain_hydrogen

function get_helium_t1(rho)
    real(8) :: get_helium_t1,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=max(log10(rho),log10_rhomin)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    if (i/=ntemp_h_he) then
        get_helium_t1=temp_he(i,1)+(temp_he(i+1,1)-temp_he(i,1))*dd/dlogrho
    else
        get_helium_t1=temp_he(i,1)
    end if
end function get_helium_t1

function get_helium_t2(rho)
    real(8) :: get_helium_t2,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i,n
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=max(log10(rho),log10_rhomin)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    if (i/=ntemp_h_he) then
        get_helium_t2=temp_he(i,2)+(temp_he(i+1,2)-temp_he(i,2))*dd/dlogrho
    else
        get_helium_t2=temp_he(i,2)
    end if
end function get_helium_t2

function get_helium_t3(rho)
    real(8) :: get_helium_t3,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=max(log10(rho),log10_rhomin)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    if (i/=ntemp_h_he) then
        get_helium_t3=temp_he(i,3)+(temp_he(i+1,3)-temp_he(i,3))*dd/dlogrho
    else
        get_helium_t3=temp_he(i,3)
    end if
end function get_helium_t3

function get_helium_t4(rho)
    real(8) :: get_helium_t4,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=max(log10(rho),log10_rhomin)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    if (i/=ntemp_h_he) then
        get_helium_t4=temp_he(i,4)+(temp_he(i+1,4)-temp_he(i,4))*dd/dlogrho
    else
        get_helium_t4=temp_he(i,4)
    end if
end function get_helium_t4

function get_hydrogen_t1(rho)
    real(8) :: get_hydrogen_t1,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=max(log10(rho),log10_rhomin)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    if (i/=ntemp_h_he) then
        get_hydrogen_t1=temp_h(i,1)+(temp_h(i+1,1)-temp_h(i,1))*dd/dlogrho
    else
        get_hydrogen_t1=temp_h(i,1)
    end if
end function get_hydrogen_t1

function get_hydrogen_t2(rho)
    real(8) :: get_hydrogen_t2,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i,n
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=max(log10(rho),log10_rhomin)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    if (i/=ntemp_h_he) then
        get_hydrogen_t2=temp_h(i,2)+(temp_h(i+1,2)-temp_h(i,2))*dd/dlogrho
    else
        get_hydrogen_t2=temp_h(i,2)
    end if
end function get_hydrogen_t2

function get_hydrogen_t3(rho)
    real(8) :: get_hydrogen_t3,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=max(log10(rho),log10_rhomin)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    if (i/=ntemp_h_he) then
        get_hydrogen_t3=temp_h(i,3)+(temp_h(i+1,3)-temp_h(i,3))*dd/dlogrho
    else
        get_hydrogen_t3=temp_h(i,3)
    end if
end function get_hydrogen_t3

function get_hydrogen_t4(rho)
    real(8) :: get_hydrogen_t4,rho,log10_rho,dlogrho,log10_rhomin,dd
    integer :: i
    log10_rhomin=tab_rho_t%xmin
    dlogrho=tab_rho_t%dx
    log10_rho=max(log10(rho),log10_rhomin)
    i=floor((log10_rho-log10_rhomin)/dlogrho+1)
    dd=log10_rho-log10_rhomin-dlogrho*(i-1)
    if (i/=ntemp_h_he) then
        get_hydrogen_t4=temp_h(i,4)+(temp_h(i+1,4)-temp_h(i,4))*dd/dlogrho
    else
        get_hydrogen_t4=temp_h(i,4)
    end if
end function get_hydrogen_t4

subroutine calspecies_division_helium_unk(rho,t,species)
    real(8) :: rho,t,t1,t2,t3,t4,species(4)
    call divide_domain_helium_unk(rho,t1,t2,t3,t4)
    if (t<=t1) then   !purely HeI
        species(1)=rho/mhe
        species(2)=0d0
        species(3)=0d0
        species(4)=0d0
    else if (t>=t4) then  !purely HeIII
        species(1)=0d0
        species(2)=0d0
        species(3)=rho/mhe
        species(4)=2*rho/mhe
    else if (t>t1.and.t<=t2) then  !HeI and HeII
        call calspecies_HeI_HeII(rho,t,species)
    else
        call calspecies_HeI_HeII_HeIII(rho,t,species)
    end if
end subroutine calspecies_division_helium_unk

subroutine calspecies_division_hydrogen_unk(rho,t,species)
    real(8) :: rho,t,t1,t2,t3,t4,species(4)
    call divide_domain_hydrogen_unk(rho,t1,t2,t3,t4)
    if (t<=t1) then   !purely molecular
        species(1)=rho/mh2
        species(2)=0d0
        species(3)=0d0
        species(4)=0d0
    else if (t>=t4) then  !fully ionized
        species(1)=0d0
        species(2)=0d0
        species(3)=rho/mh
        species(4)=rho/mh
    else if (t>t1.and.t<=t2) then  !H2,HI
        species(3)=0d0
        species(4)=0d0
        species(2)=solvenh_quadratic(rho,t)
        species(1)=species(2)**2/kdis(t)
    else if (t>t2.and.t<t3) then  !four species
        species(3)=solvenhion_quartic(rho,t)
        species(4)=species(3)
        species(2)=species(3)**2/khion(t)
        species(1)=species(2)**2/kdis(t)
    else if (t>=t3.and.t<t4) then  !HI, HII, e
        species(1)=0d0
        species(3)=solvenhion_quadratic(rho,t)
        species(2)=species(3)**2/khion(t)
        species(4)=species(3)
    end if
end subroutine calspecies_division_hydrogen_unk

subroutine divide_domain_helium_unk(rho,t1,t2,t3,t4)
    real(8) :: rho,t1,t2,t3,t4
    t1=t_helium1(rho)
    t2=t_helium2(rho)
    t3=t_helium3(rho)
    t4=t_helium4(rho)
end subroutine divide_domain_helium_unk

subroutine divide_domain_hydrogen_unk(rho,t1,t2,t3,t4)
    real(8) :: rho,t1,t2,t3,t4
    t1=t_hydrogen1(rho)
    t2=t_hydrogen2(rho)
    t3=t_hydrogen3(rho)
    t4=t_hydrogen4(rho)
end subroutine divide_domain_hydrogen_unk

function t_helium1(rho)
    real(8) :: rho,t_helium1
    t_helium1=10**3.5
end function t_helium1

function t_helium2(rho)
    real(8) :: rho,t_helium2,k1,k2
    k1=0.02d0
    k2=0.03d0
    if (rho<1e-15) then
        t_helium2=10**(4.05d0+k1*(log10(rho)+20d0))
    else
        t_helium2=10**(4d0+k2*(log10(rho)+20d0))
    end if
end function t_helium2

function t_helium3(rho)
    real(8) :: rho,t_helium3,k1,k2,k3
    k1=0.04d0
    k2=0.08d0
    k3=0.8
    if (rho>=1d-5) then
        t_helium3=min(10**(9d0+k3*log10(rho)),1.1e6)
    else if (rho>1d-10.and.rho<1d-5) then
        t_helium3=min(10**(4.6d0+k2*(log10(rho)+10d0)),1.1e6)
    else if (rho<=1d-10) then
        t_helium3=10**(4.2d0+k1*(log10(rho)+20d0))
    end if
end function t_helium3

function t_helium4(rho)
    real(8) :: rho,t_helium4,k
    k=2d0/10d0
    t_helium4=min(10**(4.5d0+k*(log10(rho)+20d0)),1.1e6)
end function t_helium4

function t_hydrogen1(rho)
    real(8) :: rho,t_hydrogen1
    t_hydrogen1=10**2.7
end function t_hydrogen1

function t_hydrogen2(rho)
    real(8) :: rho,t_hydrogen2,k
    k=0.032d0
    if (rho>1e-10) then
        t_hydrogen2=10**(3.2d0+k*10d0)
    else
        t_hydrogen2=10**(3.2d0+k*(log10(rho)+20d0))
    end if
end function t_hydrogen2

function t_hydrogen3(rho)
    real(8) :: rho,t_hydrogen3,k1,k2
    k1=0.038d0
    k2=1.25d0
    if (rho>1d-10) then
        t_hydrogen3=min(10**(3.2d0+k1*10d0+k2*(log10(rho)+10d0)),1.1e6)
    else if (rho<=1d-10) then
        t_hydrogen3=10**(3.2d0+k1*(log10(rho)+20d0))
    end if
end function t_hydrogen3

function t_hydrogen4(rho)
    real(8) :: rho,t_hydrogen4,k
    k=2d0/10d0
    t_hydrogen4=min(10**(4d0+k*(log10(rho)+20d0)),1.1e6)
end function t_hydrogen4

subroutine calspecies_helium(rho,t,species,temp)
    real(8) :: rho,t,t1,t2,t3,t4,species(4),t_he(4)
    real(8), optional :: temp(4)
    if (present(temp)) then
        t_he=temp
    else
        call divide_domain_helium(rho,t_he)
    end if
    t1=t_he(1)
    t2=t_he(2)
    t3=t_he(3)
    t4=t_he(4)
    if (t<=t1) then
        species=(/rho/mhe,0d0,0d0,0d0/)
    else if (t>=t4) then
        species=(/0d0,0d0,rho/mhe,2*rho/mhe/)
    else if (t2<=t3) then
        if (t<=t2) then
            call calspecies_HeI_HeII(rho,t,species)
        elseif (t>=t3) then
            call calspecies_HeII_HeIII(rho,t,species)
        else
            species=(/0d0,rho/mhe,0d0,rho/mhe/)
        end if
    else if (t2>t3) then
        if (t<=t3) then
            call calspecies_HeI_HeII(rho,t,species)
        elseif (t>=t2) then
            call calspecies_HeII_HeIII(rho,t,species)
        else
            call calspecies_HeI_HeII_HeIII(rho,t,species)
        end if
    end if
end subroutine calspecies_helium

subroutine calspecies_hydrogen(rho,t,species,temp)
    real(8) :: rho,t,t1,t2,t3,t4,species(4),t_h(4)
    real(8), optional :: temp(4)
    if (present(temp)) then
        t_h=temp
    else
        call divide_domain_hydrogen(rho,t_h)
    end if
    t1=t_h(1)
    t2=t_h(2)
    t3=t_h(3)
    t4=t_h(4)
    if (t<=t1) then
        species=(/rho/mh2,0d0,0d0,0d0/)
    elseif (t>=t4) then
        species=(/0d0,0d0,rho/mh,rho/mh/)
    elseif (t2<=t3) then
        if (t<=t2) then
            call calspecies_H2_HI(rho,t,species)
        elseif (t>=t3) then
            call calspecies_HI_HII(rho,t,species)
        else
            species=(/0d0,rho/mh,0d0,0d0/)
        end if
    elseif (t2>t3) then
        if (t<=t3) then
            call calspecies_H2_HI(rho,t,species)
        elseif (t>=t2) then
            call calspecies_HI_HII(rho,t,species)
        else
            call calspecies_H2_HI_HII(rho,t,species)
        end if
    end if
end subroutine calspecies_hydrogen

subroutine calspecies_HeI_HeII_HeIII(rho,t,species)
    real(8) :: rho,t,species(4)
    species(3)=solvenheion2_cubic(rho,t)
    species(4)=2*species(3)/(1-species(3)/kheion2(t))
    species(2)=species(3)*species(4)/kheion2(t)
    species(1)=species(2)*species(4)/kheion1(t)
end subroutine calspecies_HeI_HeII_HeIII

subroutine calspecies_HeII_HeIII(rho,t,species)
    real(8) :: rho,t,species(4)
    species(1)=0d0
    species(3)=solvenheion2_quadratic(rho,t)
    species(4)=species(3)+rho/mhe
    species(2)=species(3)*species(4)/kheion2(t)
end subroutine calspecies_HeII_HeIII

subroutine calspecies_HeI_HeII(rho,t,species)
    real(8) :: rho,t,species(4)
    species(2)=solvenheion1_quadratic(rho,t)
    species(1)=species(2)**2d0/kheion1(t)
    species(3)=0d0
    species(4)=species(2)
end subroutine calspecies_HeI_HeII

subroutine calspecies_H2_HI(rho,t,species)  !only H2 and HI present
    real(8) :: rho,t,species(4)
    species(3)=0d0
    species(4)=0d0
    species(2)=solvenh_quadratic(rho,t)
    species(1)=species(2)**2/kdis(t)
end subroutine calspecies_H2_HI

subroutine calspecies_HI_HII(rho,t,species)  !only HI and HII present
    real(8) :: rho,t,species(4)
    species(1)=0d0
    species(3)=solvenhion_quadratic(rho,t)
    species(2)=species(3)**2/khion(t)
    species(4)=species(3)
end subroutine calspecies_HI_HII

subroutine calspecies_H2_HI_HII(rho,t,species)   !four species
    real(8) :: rho,t,species(4)
    species(3)=solvenhion_quartic(rho,t)
    species(4)=species(3)
    species(2)=species(3)**2/khion(t)
    species(1)=species(2)**2/kdis(t)
end subroutine calspecies_H2_HI_HII

function solvenheion2_cubic(rho,t)
    real(8) :: rho,t,solvenheion2_cubic,nhetot,a(4),b,c
    complex(8) :: z(3)
    b=kheion1(t)
    c=kheion2(t)
    nhetot=rho/mhe
    a(1)=b*c**2*nhetot
    a(2)=-b*c*(c+2*nhetot)
    a(3)=b*nhetot
    a(4)=b-4*c
    call cubicroots(a,z)
    solvenheion2_cubic=dble(z(1))
end function solvenheion2_cubic

function solvenheion2_quadratic(rho,t)
    real(8) :: rho,t,solvenheion2_quadratic,nhetot,a
    a=kheion2(t)
    nhetot=rho/mhe
    solvenheion2_quadratic=2d0*nhetot*a/(nhetot+a+sqrt(nhetot**2+a**2+6*a*nhetot))
end function solvenheion2_quadratic

function solvenheion1_quadratic(rho,t)
    real(8) :: rho,t,solvenheion1_quadratic,nhetot,a
    a=kheion1(t)
    nhetot=rho/mhe
    solvenheion1_quadratic=2d0*a*nhetot/(a+sqrt(a**2+4d0*a*nhetot))
end function solvenheion1_quadratic

function solvenhion_quartic(rho,t)
    !solve the number density of HII at given rho and t, all species present
    real(8) :: rho,t,solvenhion_quartic,nhtot,a(5),b,c
    complex(8) :: z(4)
    b=khion(t)
    c=kdis(t)
    a(1)=-rho/mh
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
    a=khion(t)
    nhtot=rho/mh
    solvenhion_quadratic=2d0*a*nhtot/(a+sqrt(a**2+4d0*a*nhtot))
end function solvenhion_quadratic

function solvenh_quadratic(rho,t)
    !solve the number density of HII at given rho and t, only H2 and HI present
    real(8) :: rho,t,solvenh_quadratic,nhtot,a
    a=kdis(t)
    nhtot=rho/mh
    solvenh_quadratic=2d0*a*nhtot/(a+sqrt(a**2+8d0*a*nhtot))
end function solvenh_quadratic

function kheion2(t)
    real(8) :: kheion2,t
    kheion2=zheion2(t)*ze(t)/zheion1(t)
end function kheion2

function kheion1(t)
    real(8) :: kheion1,t
    kheion1=zheion1(t)*ze(t)/zhe(t)
end function kheion1

function khion(t)
    real(8) :: khion,t
    khion=zhion(t)*ze(t)/zh(t)
end function khion

function kdis(t)
    real(8) :: kdis,t
    kdis=zh(t)*zh(t)/zh2(t)
end function kdis

function ze(t)
    real(8) :: ze,t
    ze=zsimplegas(t,me)
end function ze

function zh(t)
    real(8) :: zh,t,ztr,zspin,zelec
    zh=zsimplegas(t,mh)/(exp(dish/(2d0*kb*t)))
end function zh

function zh2(t)
    real(8) :: zh2,t
    zh2=zsimplegas(t,mh2)
end function zh2

function zhion(t)
    real(8) :: zhion,t
    zhion=zsimplegas(t,mh)/(exp((dish+2d0*ionh)/(2d0*kb*t)))
end function zhion

function zhe(t)
    real(8) :: zhe,t
    zhe=zsimplegas(t,mhe)
end function zhe

function zheion1(t)
    real(8) :: zheion1,t
    zheion1=zsimplegas(t,mhe)/exp(ionhe1/kb/t)
end function zheion1

function zheion2(t)
    real(8) :: zheion2,t
    zheion2=zsimplegas(t,mhe)/exp((ionhe1+ionhe2)/kb/t)
end function zheion2

subroutine eos_species(blk)
    type(blockdef), pointer :: blk
    real(8) :: rho,temp,nhtot,nhetot,species(4)
    integer :: i,j
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            rho=blk%w(1,i,j,1)
            temp=blk%temp(i,j,1)
            call calspecies_hydrogen(rho,temp,species)
            nhtot=2*species(1)+species(2)+species(3)
            blk%H2(i,j,1)=2*species(1)/nhtot
            blk%HI(i,j,1)=species(2)/nhtot
            blk%HII(i,j,1)=species(3)/nhtot
            call calspecies_helium(rho,temp,species)
            nhetot=species(1)+species(2)+species(3)
            blk%HeI(i,j,1)=species(1)/nhetot
            blk%HeII(i,j,1)=species(2)/nhetot
            blk%HeIII(i,j,1)=species(3)/nhetot
            blk%cs(i,j,1)=adiabatic_cs(rho,temp)
        end do
    end do
end subroutine eos_species

end module eos_h_he

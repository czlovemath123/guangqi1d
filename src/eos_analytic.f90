module eos_analytic
use mathlib
use phylib
#if     ieosmodule==1
    use eos_h2_HI_HII
#elif   ieosmodule==2
    use eos_HI_HII
#elif   ieosmodule==3
    use eos_h_he
#elif   ieosmodule==4
    use eos_h2_HI
#endif

implicit none

contains

subroutine eos_utow(u,w,temp,egv)
    real(8), dimension(5) :: u,w
    real(8) :: rho,vx,vy,vz,egv,p,temp
    rho=u(1)
    vx=u(2)/rho
    vy=u(3)/rho
    vz=u(4)/rho
    egv=(u(5)-half*rho*(vx*vx+vy*vy+vz*vz))
    temp=solvetegv(egv,rho)
    p=prhot(rho,temp)
    w(1)=rho
    w(2)=vx
    w(3)=vy
    w(4)=vz
    w(5)=p
end subroutine eos_utow

subroutine eos_wtou(w,u,temp,egv)
    real(8), dimension(5) :: u,w
    real(8) :: rho,vx,vy,vz,egv,p,temp
    rho=w(1)
    vx=w(2)
    vy=w(3)
    vz=w(4)
    p=w(5)
    temp=solvetp(p,rho)
    egv=egvrhot(rho,temp)
    u(1)=rho
    u(2)=rho*vx
    u(3)=rho*vy
    u(4)=rho*vz
    u(5)=egv+half*rho*(vx*vx+vy*vy+vz*vz)
end subroutine eos_wtou

subroutine eos_hllc_analytic_wtou(w,egv,u)
    !given w and egv, calculate u
    real(8), dimension(5) :: w,u
    real(8) :: rho,vx,vy,vz,E,egv,p
    !for w, it is rho, vx, vy, vz, p, t, egv
    !egv is energy per unit volume
    rho=w(1)
    vx=w(2)
    vy=w(3)
    vz=w(4)
    p=w(5)
    E=egv+0.5*rho*(vx*vx+vy*vy+vz*vz)
    u(1)=rho
    u(2)=rho*vx
    u(3)=rho*vy
    u(4)=rho*vz
    u(5)=E
end subroutine eos_hllc_analytic_wtou

subroutine utosource(u,source)
    real(8) :: u(5),source(5),rho,vx,vy,vz
    source(1:4)=u(1:4)
    rho=u(1)
    vx=u(2)/rho
    vy=u(3)/rho
    vz=u(4)/rho
    source(5)=u(5)-half*rho*(vx**2+vy**2+vz**2)
end subroutine utosource

subroutine sourcetou(source,u)
    real(8) :: u(5),source(5),rho,vx,vy,vz
    u(1:4)=source(1:4)
    rho=source(1)
    vx=source(2)/rho
    vy=source(3)/rho
    vz=source(4)/rho
    u(5)=source(5)+half*rho*(vx**2+vy**2+vz**2)
end subroutine sourcetou

subroutine eos_hllc_analytic_wtoxflux(w,egv,xflux)
    !given w and egv, calculate xflux
    real(8), dimension(5) :: w,xflux
    real(8) :: rho,vx,vy,vz,egv,p,t
    !egv is energy per unit volume
    rho=w(1)
    vx=w(2)
    vy=w(3)
    vz=w(4)
    p=w(5)
    xflux(1)=rho*vx
    xflux(2)=rho*vx*vx+p
    xflux(3)=rho*vx*vy
    xflux(4)=rho*vx*vz
    xflux(5)=vx*(0.5*rho*(vx*vx+vy*vy+vz*vz)+egv+p)
end subroutine eos_hllc_analytic_wtoxflux

subroutine generate_rho_t_eos()
    real(8) :: rhologspan,tlogspan,dlogrho,dlogt
    integer :: i,j
    allocate(environment%rho_eos(environment%nrho),environment%t_eos(environment%ntemp))
    rhologspan=log10(environment%rhomax/environment%rhomin)
    dlogrho=rhologspan/(environment%nrho-1)
    tlogspan=log10(environment%tmax/environment%tmin)
    dlogt=tlogspan/(environment%ntemp-1)
    do i=1,environment%ntemp
        environment%t_eos(i)=environment%tmin*10**((i-1)*dlogt)
    end do
    do j=1,environment%nrho
        environment%rho_eos(j)=environment%rhomin*10**((j-1)*dlogrho)
    end do
    tab_rho_t%xmin=log10(environment%rhomin)
    tab_rho_t%xmax=log10(environment%rhomax)
    tab_rho_t%dx=dlogrho
end subroutine generate_rho_t_eos

end module eos_analytic

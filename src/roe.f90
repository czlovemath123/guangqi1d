module roe
use datastructure
use mathlib
use phylib
implicit none
contains
subroutine calxflux(wl,wr,xflux,vmax)
    real(8) :: xflux(5),xlflux(5),xrflux(5),wl(5),wr(5),ul(5),ur(5),du(5),  &
                alpha(5),lambda(5),k(5,5)
    real(8) :: rhot,ut,vt,wt,ht,at,speedt,du5bar
    real(8) :: hl,hr
    real(8) :: pstar,ustar,rhostarl,rhostarr,astarl,astarr
    real(8) :: al,ar,lambda1l,lambda1r,lambda5l,lambda5r,vmax
    integer :: i
    hl=(0.5*wl(1)*(wl(2)**2+wl(3)**2+wl(4)**2)+wl(5)/(gamma_gas-1)+wl(5))/wl(1)
    hr=(0.5*wr(1)*(wr(2)**2+wr(3)**2+wr(4)**2)+wr(5)/(gamma_gas-1)+wr(5))/wr(1)
    rhot=sqrt(wl(1)*wr(1))
    ut=(sqrt(wl(1))*wl(2)+sqrt(wr(1))*wr(2))/(sqrt(wl(1))+sqrt(wr(1)))
    vt=(sqrt(wl(1))*wl(3)+sqrt(wr(1))*wr(3))/(sqrt(wl(1))+sqrt(wr(1)))
    wt=(sqrt(wl(1))*wl(4)+sqrt(wr(1))*wr(4))/(sqrt(wl(1))+sqrt(wr(1)))
    speedt=sqrt(ut**2+vt**2+wt**2)
    ht=(sqrt(wl(1))*hl+sqrt(wr(1))*hr)/(sqrt(wl(1))+sqrt(wr(1)))
    at=sqrt((gamma_gas-1)*(ht-0.5*speedt**2))
    k(1,1)=1
    k(2,1)=ut-at
    k(3,1)=vt
    k(4,1)=wt
    k(5,1)=ht-ut*at
    k(1,2)=1
    k(2,2)=ut
    k(3,2)=vt
    k(4,2)=wt
    k(5,2)=0.5*speedt**2
    k(1,3)=0
    k(2,3)=0
    k(3,3)=1
    k(4,3)=0
    k(5,3)=vt
    k(1,4)=0
    k(2,4)=0
    k(3,4)=0
    k(4,4)=1
    k(5,4)=wt
    k(1,5)=1
    k(2,5)=ut+at
    k(3,5)=vt
    k(4,5)=wt
    k(5,5)=ht+ut*at
    lambda(1)=ut-at
    lambda(2:4)=ut
    lambda(5)=ut+at
    call wtou(wl,ul)
    call wtou(wr,ur)
    du=ur-ul
    du5bar=du(5)-(du(3)-vt*du(1))*vt-(du(4)-wt*du(1))*wt
    alpha(2)=(gamma_gas-1)/at/at*(du(1)*(ht-ut*ut)+ut*du(2)-du5bar)
    alpha(1)=0.5/at*(du(1)*(ut+at)-du(2)-at*alpha(2))
    alpha(3)=du(3)-vt*du(1)
    alpha(4)=du(4)-wt*du(1)
    alpha(5)=du(1)-(alpha(1)+alpha(2))
    !identify rarefaction wave and impose entropy fix. if no rarefaction wave, no need
    !to add entropy fix
    !adaptive non-iterative riemann solver
    call anrsx(wl,wr,pstar,ustar,rhostarl,rhostarr)
    !roe averaged riemann solver
    !rhostarl=wl(1)+alpha(1)
    !ustar=(wl(1)*wl(2)+alpha(1)*(ut-at))/(wl(1)+alpha(1))
    !pstar=(gamma_gas-1)*(ul(1)+alpha(1)*(ht-ut*at)-0.5d0*rhostarl*ustar**2)
    if (pstar<wl(5)) then
        al=sqrt(gamma_gas*wl(5)/wl(1))
        astarl=sqrt(gamma_gas*pstar/rhostarl)
        lambda1l=wl(2)-al
        lambda1r=ustar-astarl
        if (lambda1l<0 .and. lambda1r>0) then
            lambda(1)=lambda1l*(lambda1r-lambda(1))/(lambda1r-lambda1l)
        end if
    end if
    if (pstar<wr(5)) then
        ar=sqrt(gamma_gas*wr(5)/wr(1))
        astarr=sqrt(gamma_gas*pstar/rhostarr)
        lambda5l=ustar+astarr
        lambda5r=wr(2)+ar
        if (lambda5l<0 .and. lambda5r>0) then
            lambda(5)=lambda5r*(lambda(5)-lambda5l)/(lambda5r-lambda5l)
        end if
    end if
    call wtoxflux(wl,xlflux)
    call wtoxflux(wr,xrflux)
    xflux=0.5*(xlflux+xrflux)
    do i=1,5
        xflux=xflux-0.5*alpha(i)*abs(lambda(i))*k(:,i)
    end do
    vmax=max(abs(lambda(1)),abs(lambda(5)))
end subroutine calxflux

subroutine roe_averaged_rs(wl,wr,pstar)
    real(8) :: wl(5),wr(5),ul(5),ur(5),alpha(2),du(5),pstar
    real(8) :: rhot,ut,vt,wt,ht,at,speedt,du5bar
    real(8) :: hl,hr,rhostarl,ustar
    hl=(0.5*wl(1)*(wl(2)**2+wl(3)**2+wl(4)**2)+wl(5)/(gamma_gas-1)+wl(5))/wl(1)
    hr=(0.5*wr(1)*(wr(2)**2+wr(3)**2+wr(4)**2)+wr(5)/(gamma_gas-1)+wr(5))/wr(1)
    rhot=sqrt(wl(1)*wr(1))
    ut=(sqrt(wl(1))*wl(2)+sqrt(wr(1))*wr(2))/(sqrt(wl(1))+sqrt(wr(1)))
    vt=(sqrt(wl(1))*wl(3)+sqrt(wr(1))*wr(3))/(sqrt(wl(1))+sqrt(wr(1)))
    wt=(sqrt(wl(1))*wl(4)+sqrt(wr(1))*wr(4))/(sqrt(wl(1))+sqrt(wr(1)))
    speedt=sqrt(ut**2+vt**2+wt**2)
    ht=(sqrt(wl(1))*hl+sqrt(wr(1))*hr)/(sqrt(wl(1))+sqrt(wr(1)))
    at=sqrt((gamma_gas-1)*(ht-0.5*speedt**2))
    call wtou(wl,ul)
    call wtou(wr,ur)
    du=ur-ul
    alpha(2)=(gamma_gas-1)/at/at*(du(1)*(ht-ut*ut)+ut*du(2)-du5bar)
    alpha(1)=0.5/at*(du(1)*(ut+at)-du(2)-at*alpha(2))
    rhostarl=wl(1)+alpha(1)
    ustar=(wl(1)*wl(2)+alpha(1)*(ut-at))/(wl(1)+alpha(1))
    pstar=(gamma_gas-1)*(ul(5)+alpha(1)*(ht-ut*at)-0.5d0*rhostarl*ustar**2)
end subroutine roe_averaged_rs

subroutine anrsx(wl,wr,pstar,ustar,rhostarl,rhostarr)
    real(8) :: wl(5),wr(5)
    real(8) :: pstar,ustar,rhostarl,rhostarr
    real(8) :: ppvrs,pmin,pmax,qanrs,cl,cr,al,ar,z,fl,fr,gl,gr,p0
    real(8) :: aal,aar,bbl,bbr
    !primitive variable riemann solver
    qanrs=2.0
    al=sqrt(gamma_gas*wl(5)/wl(1))
    ar=sqrt(gamma_gas*wr(5)/wr(1))
    cl=wl(1)*al
    cr=wr(1)*ar
    ppvrs=1/(cl+cr)*(cr*wl(5)+cl*wr(5)+cl*cr*(wl(2)-wr(2)))
    pmin=min(wl(5),wr(5))
    pmax=max(wl(5),wr(5))
    if (pmax/pmin<qanrs .and. ppvrs<pmax .and. pmin<ppvrs) then
        pstar=ppvrs
        ustar=1/(cl+cr)*(cl*wl(2)+cr*wr(2)+(wl(5)-wr(5)))
        rhostarl=wl(1)+(pstar-wl(5))/al/al
        rhostarr=wr(1)+(pstar-wr(5))/ar/ar
    !two shock riemann solver
    else if (pmin<ppvrs) then
        p0=max(0.0,ppvrs)
        aal=1/(gamma_gas+1)/wl(1)
        aar=1/(gamma_gas+1)/wr(1)
        bbl=(gamma_gas-1)/(gamma_gas+1)*wl(5)
        bbr=(gamma_gas-1)/(gamma_gas+1)*wr(5)
        gl=sqrt(aal/(p0+bbl))
        gr=sqrt(aar/(p0+bbr))
        pstar=(gl*wl(5)+gr*wr(5)-(wr(2)-wl(2)))/(gl+gr)
        ustar=0.5*(wl(2)+wr(2))+0.5*((pstar-wr(5))*gr-(pstar-wl(5))*gl)
        rhostarl=wl(1)*(pstar/wl(5)+(gamma_gas-1)/(gamma_gas+1))/((gamma_gas-1)/(gamma_gas+1)*pstar/wl(5)+1)
        rhostarr=wr(1)*(pstar/wr(5)+(gamma_gas-1)/(gamma_gas+1))/((gamma_gas-1)/(gamma_gas+1)*pstar/wr(5)+1)
    !two rarefaction riemann solver
    else
        z=(gamma_gas-1)/gamma_gas/2
        pstar=pow((al+ar-(gamma_gas-1)/2*(wr(2)-wl(2)))/(al/pow(wl(5),z)+ar/pow(wr(5),z)),1/z)
        rhostarl=wl(1)*pow(pstar/wl(5),1/gamma_gas)
        rhostarr=wr(1)*pow(pstar/wr(5),1/gamma_gas)
        fl=2*al/(gamma_gas-1)*(pow(pstar/wl(5),z)-1)
        fr=2*ar/(gamma_gas-1)*(pow(pstar/wr(5),z)-1)
        ustar=0.5*(wl(2)+wr(2))+0.5*(fr-fl)
    end if
end subroutine anrsx

end module roe

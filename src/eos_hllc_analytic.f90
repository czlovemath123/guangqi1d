module eos_hllc_analytic
use datastructure
use phylib
use mathlib
use passivescalars
use eos_analytic
implicit none
contains

function eos_analytic_gleft(wl,p,gamm)
    real(8), dimension(5) :: wl
    real(8) :: p,al,bl,eos_analytic_gleft,gamm
    al=2/(gamm+1)/wl(1)
    bl=(gamm-1)/(gamm+1)*wl(5)
    eos_analytic_gleft=sqrt(al/(p+bl))
end function eos_analytic_gleft

function eos_analytic_gright(wr,p,gamm)
    real(8), dimension(5) :: wr
    real(8) :: p,ar,br,eos_analytic_gright,gamm
    ar=2/(gamm+1)/wr(1)
    br=(gamm-1)/(gamm+1)*wr(5)
    eos_analytic_gright=sqrt(ar/(p+br))
end function eos_analytic_gright

subroutine eos_analytic_hllc_pvrs(wl,wr,pvrs,gammal,gammar)
    real(8), dimension(5) :: wl,wr
    real(8) :: ppvrs,pvrs,rhoa,aa,al,ar,gammal,gammar
    al=sqrt(gammal*wl(5)/wl(1))
    ar=sqrt(gammar*wr(5)/wr(1))
    aa=half*(al+ar)
    rhoa=half*(WL(1)+WR(1))
    ppvrs=half*(WL(5)+WR(5))-half*(WR(2)-WL(2))*rhoa*aa
    pvrs=max(zero,ppvrs)
end subroutine eos_analytic_hllc_pvrs

subroutine eos_analytic_hllc_trrs(wl,wr,p,gammal,gammar,gamm)
    real(8), dimension(5) :: wl,wr
    real(8) :: p,z,al,ar,gammal,gammar,gamm,zeta
    al=sqrt(gammal*wl(5)/wl(1))
    ar=sqrt(gammar*wr(5)/wr(1))
    z=(gamm-1)/two/gamm
    zeta=max((al+ar-gamm*z*(wr(2)-wl(2)))/(al/pow(wl(5),z)+ar/pow(wr(5),z)),0d0)
    p=max(pow(zeta,1d0/z),pfloor)
end subroutine eos_analytic_hllc_trrs

subroutine eos_analytic_hllc_tsrs(wl,wr,p,gammal,gammar,gamm)
    real(8), dimension(5) :: wl,wr
    real(8) :: p,gl,gr,p0,pvrs,gammal,gammar,gamm,x
    call eos_analytic_hllc_pvrs(wl,wr,pvrs,gammal,gammar)
    p0=max(zero,pvrs)
    gl=eos_analytic_gleft(wl,p0,gammal)
    gr=eos_analytic_gright(wr,p0,gammar)
    p=(gl*wl(5)+gr*wr(5)-(wr(2)-wl(2)))/(gl+gr)
end subroutine eos_analytic_hllc_tsrs

!eos_analytic_hllc_pstar should do the same thing as hllc_pstar
subroutine eos_analytic_hllc_pstar(wl,wr,pstar,gammal,gammar,gamm)
    real(8), dimension(5) :: wl,wr
    real(8) :: pstar,q,qmax,pvrs,trrs,tsrs,pmax,pmin,gammal,gammar,gamm
    call eos_analytic_hllc_pvrs(wl,wr,pvrs,gammal,gammar)
    pmax=max(wl(5),wr(5))
    pmin=min(wl(5),wr(5))
    q=pmax/pmin
    qmax=2.0
    if (q<=qmax.and.pvrs<=pmax.and.pvrs>=pmin) then
        pstar=pvrs
    else
        if (pvrs<pmin) then
            call eos_analytic_hllc_trrs(wl,wr,pstar,gammal,gammar,gamm)
        else
            call eos_analytic_hllc_tsrs(wl,wr,pstar,gammal,gammar,gamm)
        end if
    end if
end subroutine eos_analytic_hllc_pstar

subroutine eos_hllc_analytic_slsrsstar(wl,wr,sl,sr,sstar,egv,gamm,gammal,gammar,pstar)
    !Hybrid estimate in "Restoration of the contact surface in the HLL-Riemann solver"
    !calculate sl, sr, and sstar
    real(8), dimension(5) :: wl,wr
    real(8) :: sl,sr,sl_temp,sr_temp,sstar,pstar,wsl,wsr,al,ar,ql,qr,gamm,gammal,gammar
    real(8) :: sstar_real(2),p_jump,temp(2),egv(2)
    call eos_analytic_hllc_pstar(wl,wr,pstar,gammal,gammar,gamm)
    !calculate sl,sr and sstar
    if (pstar>wl(5)) then
        !left shock
        al=sqrt(gammal*wl(5)/wl(1))
        ql=sqrt(1d0+(gammal+1d0)/2d0/gammal*(pstar/wl(5)-1d0))
        sl_temp=wl(2)-al*ql
    else
        !left rarefaction
        al=sqrt(gammal*wl(5)/wl(1))
        sl_temp=wl(2)-al
    end if
    if (pstar>wr(5)) then
        !right shock
        ar=sqrt(gammar*wr(5)/wr(1))
        !ar=sqrt(min(gammal,gammar)*wr(5)/wr(1))
        qr=sqrt(1d0+(gammar+1d0)/2d0/gammar*(pstar/wr(5)-1d0))
        sr_temp=wr(2)+ar*qr
    else
        !right rarefaction
        ar=sqrt(gammar*wr(5)/wr(1))
        sr_temp=wr(2)+ar
    end if
    sl=min(sl_temp,sr_temp)
    sr=max(sl_temp,sr_temp)
    sstar=(wr(5)-wl(5)+wl(1)*wl(2)*(sl-wl(2))-wr(1)*wr(2)*(sr-wr(2)))/(wl(1)*(sl-wl(2))-wr(1)*(sr-wr(2)))
end subroutine eos_hllc_analytic_slsrsstar

subroutine eos_hllc_analytic_ustar(w,egv,sstar,s,ustar)
    !given w, egv, sstar and sl or sr, calculate ustar
    real(8), dimension(5) :: w,ustar
    real(8) :: rho,t,sstar,h,egv,s,E,vx,vy,vz
    !egv is energy per unit volume (include radiation if present)
    !E is total energy
    h=w(1)*(s-w(2))/(s-sstar)
    ustar(1)=h
    ustar(2)=h*sstar
    ustar(3)=h*w(3)
    ustar(4)=h*w(4)
    rho=w(1)
    vx=w(2)
    vy=w(3)
    vz=w(4)
    E=egv+half*rho*(vx*vx+vy*vy+vz*vz)
    ustar(5)=h*(E/w(1)+(sstar-w(2))*(sstar+w(5)/w(1)/(s-w(2))))
end subroutine eos_hllc_analytic_ustar

subroutine eos_hllc_analytic_local_gamma(wl,wr,gamm,gammal,gammar,temp)
    !calculate gammas
    real(8), dimension(5) :: wl,wr
    real(8) :: gammal,gammar,p,rho,t
    real(8) :: gamm,pstar,e,rho_post,e_post,wsl,wsr,temp(2)
    integer :: i
    character(len=128) :: alert
    rho=wl(1)
    t=temp(1)
    call gammarhot(rho,t,gammal)
    rho=wr(1)
    t=temp(2)
    call gammarhot(rho,t,gammar)
    gamm=max(gammal,gammar)
end subroutine eos_hllc_analytic_local_gamma

subroutine eos_hllc_analytic_flux(wl,wr,temp,egv,flux,vlocalmax)
    !extended infomation:rho,vx,vy,vz,p,temp,egv
    !given left and right states (wl,wr), temp, and egv, calculate flux and vlocalmax
    real(8), dimension(5) :: wl,wr,flux,flux_temp,ur,ul,ustarl,ustarr
    real(8) :: sl,sr,sstar,vlocalmax,gamm,gammal,gammar,pstar,temp(2),egv(2)
    call eos_hllc_analytic_local_gamma(wl,wr,gamm,gammal,gammar,temp)
    call eos_hllc_analytic_slsrsstar(wl,wr,sl,sr,sstar,egv,gamm,gammal,gammar,pstar)
    vlocalmax=max(abs(sl),abs(sr))
    if (sl>0) then
        call eos_hllc_analytic_wtoxflux(wl,egv(1),flux)
    else if (sl<=0.and.sstar>0) then
        call eos_hllc_analytic_wtoxflux(wl,egv(1),flux_temp)
        call eos_hllc_analytic_wtou(wl,egv(1),ul)
        call eos_hllc_analytic_ustar(wl,egv(1),sstar,sl,ustarl)
        flux=flux_temp+sl*(ustarl-ul)
    else if (sstar<=0.and.sr>0) then
        call eos_hllc_analytic_wtoxflux(wr,egv(2),flux_temp)
        call eos_hllc_analytic_wtou(wr,egv(2),ur)
        call eos_hllc_analytic_ustar(wr,egv(2),sstar,sr,ustarr)
        flux=flux_temp+sr*(ustarr-ur)
    else
        call eos_hllc_analytic_wtoxflux(wr,egv(2),flux)
    end if
end subroutine eos_hllc_analytic_flux

subroutine hllc_muscl(blk)
    !use W_xl, W_xr, W_yl, W_yr to calculate fluxes
    type(blockdef), pointer :: blk
    real(8), dimension(5) :: wl,wr,flux
    real(8), dimension(:,:,:), allocatable :: vmaxarray
    real(8) :: vlocalmax,temp(2),egv(2),vblockmax(nd),dt(nd)
    integer :: i,j
    !for 1d the default is x direction
    blk%xflux=0d0
    do i=0,blk_size_nx
        wl=blk%w_xr(1:5,i,1,1)
        wr=blk%w_xl(1:5,i+1,1,1)
        temp(1)=blk%temp_xr(i,1,1)
        temp(2)=blk%temp_xl(i+1,1,1)
        egv(1)=blk%egv_xr(i,1,1)
        egv(2)=blk%egv_xl(i+1,1,1)
        call eos_hllc_analytic_flux(wl,wr,temp,egv,flux,vlocalmax)
        blk%xflux(1:5,i,1,1)=flux
        blk%hllc_vx(i,1,1)=vlocalmax
    end do
    if (lpassive) call upwind_passive(blk)
    if (lrad_adv) then
        blk%erad_xflux=0d0
        do i=0,blk_size_nx
            if (blk%xflux(1,i,1,1)>0) then
                blk%erad_xflux(i,1,1)=blk%xflux(1,i,1,1)/blk%w(1,i,1,1)*blk%erad(i,1,1)
            else
                blk%erad_xflux(i,1,1)=blk%xflux(1,i,1,1)/blk%w(1,i+1,1,1)*blk%erad(i+1,1,1)
            end if
        end do
    end if
end subroutine hllc_muscl

end module eos_hllc_analytic

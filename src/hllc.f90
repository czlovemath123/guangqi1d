module hllc
use phylib
use mathlib
use constant_gamma_eos
use passivescalars
use datastructure
implicit none
contains

!P. Batten, N. Clarke, C. Lambert, and D. M. Causon "On the choince of wavespeeds for the HLLC Riemann solver" 1997

subroutine hllc_roe_averaged(wl,wr,q_roe,c_roe)
    !use Roe average Riemann solver to find q_roe, c_roe
    real(8), dimension(5) :: wl,wr,w_roe
    real(8) :: egv(2),hl,hr,rhot,ut,vt,wt,ht,h,tt,q_roe,c_roe
    real(8) :: sqrtrhol,sqrtrhor,sqrtrho,ep
    egv=(/wl(5)/(gamma_gas-1d0),wr(5)/(gamma_gas-1d0)/)
    hl=(half*wl(1)*(wl(2)**2+wl(3)**2+wl(4)**2)+egv(1)+wl(5))/wl(1)
    hr=(half*wr(1)*(wr(2)**2+wr(3)**2+wr(4)**2)+egv(2)+wr(5))/wr(1)
    rhot=sqrt(wl(1)*wr(1))
    sqrtrhol=sqrt(wl(1))
    sqrtrhor=sqrt(wr(1))
    sqrtrho=sqrt(wl(1))+sqrt(wr(1))
    ut=(sqrtrhol*wl(2)+sqrtrhor*wr(2))/sqrtrho
    vt=(sqrtrhol*wl(3)+sqrtrhor*wr(3))/sqrtrho
    wt=(sqrtrhol*wl(4)+sqrtrhor*wr(4))/sqrtrho
    ht=(sqrtrhol*hl+sqrtrhor*hr)/sqrtrho
    ep=(ht-half*(ut**2+vt**2+wt**2))
    c_roe=sqrt((gamma_gas-1d0)*ep)
    !c_roe=sqrt((gamma_gas-1d0)*ep/gamma_gas)
    q_roe=ut
end subroutine hllc_roe_averaged

function hllc_p(wl,wr)
    real(8) :: wl(5),wr(5),hllc_p
    real(8) :: sl,sr,sm,vlocalmax,egv(2)
    real(8) :: ql,qr,q_roe,cl,cr,c_roe,dsl,dsr,dsm
    real(8) :: rhom,pm,momum,momvm,momwm,el,er,em
    call hllc_roe_averaged(wl,wr,q_roe,c_roe)
    ql=wl(2)
    qr=wr(2)
    cl=sqrt(gamma_gas*wl(5)/wl(1))
    cr=sqrt(gamma_gas*wr(5)/wr(1))
    egv(1)=wl(5)/(gamma_gas-1d0)
    egv(2)=wr(5)/(gamma_gas-1d0)
    sl=min(ql-cl,q_roe-c_roe)
    sr=max(qr+cr,q_roe+c_roe)
    dsl=sl-ql
    dsr=sr-qr
    sm=(wr(1)*qr*dsr-wl(1)*ql*dsl+wl(5)-wr(5))/(wr(1)*dsr-wl(1)*dsl)
    if (sl>zero) then
        hllc_p=wl(5)
    else if (sl<=zero.and.sm>zero) then
        dsm=sl-sm
        rhom=wl(1)*dsl/dsm
        pm=wl(5)-wl(1)*dsl*(ql-sm)
        momum=(dsl*wl(1)*wl(2)+(pm-wl(5)))/dsm
        momvm=dsl*wl(1)*wl(3)/dsm
        momwm=dsl*wl(1)*wl(4)/dsm
        el=egv(1)+half*wl(1)*(wl(2)**2+wl(3)**2+wl(4)**2)
        em=(dsl*el-ql*wl(5)+pm*sm)/dsm
        hllc_p=pm
    else if (sm<=zero.and.sr>zero) then
        dsm=sr-sm
        rhom=wr(1)*dsr/dsm
        pm=wr(5)-wr(1)*dsr*(qr-sm)
        momum=(dsr*wr(1)*wr(2)+(pm-wr(5)))/dsm
        momvm=dsr*wr(1)*wr(3)/dsm
        momwm=dsr*wr(1)*wr(4)/dsm
        er=egv(2)+half*wr(1)*(wr(2)**2+wr(3)**2+wr(4)**2)
        em=(dsr*er-qr*wr(5)+pm*sm)/dsm
        hllc_p=pm
    else
        hllc_p=wr(5)
    end if
end function hllc_p

subroutine hllc_flux(wl,wr,flux,vlocalmax)
    !given left and right states (wl,wr), calculate flux and vlocalmax
    real(8), dimension(5) :: wl,wr,flux,ul,ur,fl,fr
    real(8) :: sl,sr,sm,vlocalmax,egv(2)
    real(8) :: ql,qr,q_roe,cl,cr,c_roe,dsl,dsr,dsm,bl,br
    real(8) :: rhom,pm,momum,momvm,momwm,el,er,em,ustar(5),wstar(5),pml,pmr
    call hllc_roe_averaged(wl,wr,q_roe,c_roe)
    ql=wl(2)
    qr=wr(2)
    cl=sqrt(gamma_gas*wl(5)/wl(1))
    cr=sqrt(gamma_gas*wr(5)/wr(1))
    bl=sqrt((gamma_gas-1d0)/2d0/gamma_gas)
    br=sqrt((gamma_gas-1d0)/2d0/gamma_gas)
    egv(1)=wl(5)/(gamma_gas-1d0)
    egv(2)=wr(5)/(gamma_gas-1d0)
    sl=min(ql-bl*cl,q_roe-c_roe)
    sr=max(qr+br*cr,q_roe+c_roe)
    dsl=sl-ql
    dsr=sr-qr
    sm=(wr(1)*qr*dsr-wl(1)*ql*dsl+wl(5)-wr(5))/(wr(1)*dsr-wl(1)*dsl)
    if (sl>zero) then
        call wtoflux(wl,flux)
    else if (sl<=zero.and.sm>zero) then
        dsm=sl-sm
        rhom=wl(1)*dsl/dsm
        pm=wl(5)-wl(1)*dsl*(ql-sm)
        momum=(dsl*wl(1)*wl(2)+(pm-wl(5)))/dsm
        momvm=dsl*wl(1)*wl(3)/dsm
        momwm=dsl*wl(1)*wl(4)/dsm
        el=egv(1)+half*wl(1)*(wl(2)**2+wl(3)**2+wl(4)**2)
        em=(dsl*el-ql*wl(5)+pm*sm)/dsm
        flux(1)=rhom*sm
        flux(2)=momum*sm+pm
        flux(3)=momvm*sm
        flux(4)=momwm*sm
        flux(5)=(em+pm)*sm
    else if (sm==0d0) then
        dsm=sl-sm
        rhom=wl(1)*dsl/dsm
        pml=wl(5)-wl(1)*dsl*(ql-sm)
        dsm=sr-sm
        rhom=wr(1)*dsr/dsm
        pmr=wr(5)-wr(1)*dsr*(qr-sm)
        pm=(pml+pmr)/2d0
        flux=0d0
        flux(2)=pm
    else if (sm<zero.and.sr>=zero) then
        dsm=sr-sm
        rhom=wr(1)*dsr/dsm
        pm=wr(5)-wr(1)*dsr*(qr-sm)
        momum=(dsr*wr(1)*wr(2)+(pm-wr(5)))/dsm
        momvm=dsr*wr(1)*wr(3)/dsm
        momwm=dsr*wr(1)*wr(4)/dsm
        er=egv(2)+half*wr(1)*(wr(2)**2+wr(3)**2+wr(4)**2)
        em=(dsr*er-qr*wr(5)+pm*sm)/dsm
        flux(1)=rhom*sm
        flux(2)=momum*sm+pm
        flux(3)=momvm*sm
        flux(4)=momwm*sm
        flux(5)=(em+pm)*sm
    else
        call wtoflux(wr,flux)
    end if
    vlocalmax=max(abs(sl),abs(sr))
end subroutine hllc_flux

subroutine hllc_muscl(blk)
    !use W_xl, W_xr, W_yl, W_yr to calculate fluxes
    type(blockdef), pointer :: blk
    real(8), dimension(5) :: wl,wr,flux
    real(8), dimension(:,:,:), allocatable :: vmaxarray
    real(8) :: vlocalmax,vblockmax,dt(nd),pm
    integer :: i,j,k
    !for 1d, the default is x direction
    blk%xflux=0d0
    do i=0,blk_size_nx
        wl=blk%w_xr(1:5,i,1,1)
        wr=blk%w_xl(1:5,i+1,1,1)
        call hllc_flux(wl,wr,flux,vlocalmax)
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

end module hllc

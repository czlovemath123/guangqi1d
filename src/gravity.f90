module gravity
use datastructure
use mathlib
use phylib
use eos
implicit none

contains

subroutine gravity_source(blk,dt,step)
    !gravity source is second order, t=0 is stored in s_grav1 and t=dt is stored in s_grav2
    type(blockdef), pointer :: blk
    real(8) :: dt
    character(len=16), optional :: step
    character(len=128) :: alert
    if (igravity==1) then
        call central_gravity(blk,dt,step)
    else if (igravity==2) then
        call uniform_gravity(blk,dt)
    end if
end subroutine gravity_source

subroutine central_gravity(blk,dt,step)
    type(blockdef), pointer :: blk
    real(8) :: m,dt,u(5),u1(5),r,egv,rho,v,x,vol
    real(8) :: ril,rir,rhoxl,rhoxr,vxl,vxr,rhox,vx
    integer :: i,j
    character(len=16) :: step
    m=central_star%core%mass
    if (igeometry==1) then
    else if (igeometry==2) then
        if (step==predictor) then
            blk%u_muscl=blk%u0
        else if (step==corrector) then
            blk%u_muscl=blk%u1
        end if
        do i=1,blk_size_nx
            u1=blk%u(1:5,i,1,1)
            u=blk%u_muscl(1:5,i,1,1)
            vol=blk%vol(i,1,1)
            rho=u(1)
            v=u(2)/u(1)
            r=blk%x_center(i)
            ril=blk%x_interface(i-1)
            rir=blk%x_interface(i)
            rhoxl=blk%w_xl(1,i,1,1)
            rhoxr=blk%w_xr(1,i,1,1)
            vxl=blk%w_xl(2,i,1,1)
            vxr=blk%w_xr(2,i,1,1)
            rhox=(rhoxr-rhoxl)/(rir-ril)
            vx=(vxr-vxl)/(rir-ril)
            u1(2)=u1(2)-m*gr/vol*gravity_source_mom(rho,rhox,r,ril,rir)*dt
            u1(5)=u1(5)-m*gr/vol*gravity_source_energy(rho,v,rhox,vx,r,ril,rir)*dt
            !u1(2)=u1(2)-rho*m*gr/r**2d0*dt
            !u1(5)=u1(5)-u(2)*m*gr/r**2d0*dt
            egv=u1(5)-half*u1(2)**2d0/u1(1)
            x=r/n_domain(1)/tempfloorscale
            if (egv<flooregv(u1(1),x)) then
                egv=flooregv(u1(1),x)
            end if
            u1(5)=egv+half*u1(2)**2d0/u1(1)
            blk%u(1:5,i,1,1)=u1
        end do
    end if
end subroutine central_gravity

function gravity_source_mom(rho0,rhox,x0,x1,x2)
    real(8) :: rho0,rhox,x0,x1,x2,gravity_source_mom
    !gravity_source_mom=-(rho0*x1) + rhox*x0*x1 - (rhox*x1**2)/2. + rho0*x2 - rhox*x0*x2 + (rhox*x2**2)/2.
    gravity_source_mom=-0.5*((x1 - x2)*(2*rho0 + rhox*(-2*x0 + x1 + x2)))
end function gravity_source_mom

function gravity_source_energy(rho0,v0,rhox,vx,x0,x1,x2)
    real(8) :: rho0,v0,rhox,vx,x0,x1,x2,gravity_source_energy
    !gravity_source_energy=-((x1 - x2)*(6*rho0*v0 + 3*rhox*v0*(-2*x0 + x1 + x2) + 3*rho0*vx*(-2*x0 + x1 + x2) + 2*rhox*vx*(3*x0**2 + x1**2 + x1*x2 + x2**2 - 3*x0*(x1 + x2))))/6.0
    gravity_source_energy=-0.5*((x1 - x2)*(2*rho0*v0 + rhox*v0*(-2*x0 + x1 + x2) + rho0*vx*(-2*x0 + x1 + x2)))
end function gravity_source_energy

subroutine uniform_gravity(blk,dt)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dt,u(5),egv,r,theta
    dt=time_sys%dt_hydro
    if (igeometry==0) then
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                u=blk%u(1:5,i,j,1)
                egv=u(5)-half*(u(2)**2d0+u(3)**2d0)/u(1)
                u(3)=u(3)-u(1)*g_uniform*dt
                u(5)=egv+half*(u(2)**2d0+u(3)**2d0)/u(1)
                blk%u(1:5,i,j,1)=u
            end do
        end do
    else if (igeometry==2) then
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                u=blk%u(1:5,i,j,1)
                egv=u(5)-half*(u(2)**2d0+u(3)**2d0)/u(1)
                r=blk%x_center(i)
                theta=blk%y_center(j)
                u(2)=u(2)-u(1)*g_uniform*cos(theta)*dt
                u(3)=u(3)-u(1)*g_uniform*sin(theta)*dt
                u(5)=egv+half*(u(2)**2d0+u(3)**2d0)/u(1)
                blk%u(1:5,i,j,1)=u
            end do
        end do
    end if
end subroutine uniform_gravity

subroutine ff_timescale(m,r,dt)
    real(8) :: m,r,dt
    dt=pi/2d0*pow(r,1.5d0)/sqrt(2*gr*m)
end subroutine ff_timescale

end module gravity

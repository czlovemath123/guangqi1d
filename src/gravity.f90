module gravity
use datastructure
use mathlib
use phylib
use eos
use radiation_common_functions
implicit none

interface gravity_source_mom
    module procedure gravity_source_mom_1d
end interface gravity_source_mom

interface gravity_source_energy
    module procedure gravity_source_energy_1d
end interface gravity_source_energy

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
    real(8) :: m,dt,u(5),u1(5),r,y,egv,rho,omega,theta,x,ril,rir,yil,yir, &
        rhoxl,rhoxr,rhoyl,rhoyr,rhox,rhoy,vx,vxx,vxy,vol
    real(8), allocatable :: acc(:,:,:)
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
        !call allocate_cell_data_block(acc)
        !if (lirradiation) then
        !    call fld_irradiation_1d(blk)
        !    blk%aradx_irrad(1:blk_size_nx,1,1)=blk%irrad(1:blk_size_nx,1,1)/blk%w(1,1:blk_size_nx,1,1)/c_light
        !    acc=blk%aradx_irrad
        !end if
        do i=1,blk_size_nx
            u1=blk%u(1:5,i,1,1)
            u=blk%u_muscl(1:5,i,1,1)
            rho=u(1)
            r=blk%x_center(i)
            vol=blk%vol(i,1,1)
            ril=blk%x_interface(i-1)
            rir=blk%x_interface(i)
            rhox=blk%xslp(1,i,1,1)
            vx=u(2)/u(1)
            vxx=blk%xslp(2,i,1,1)
            u1(2)=u1(2)-m*gr*dt*gravity_source_mom(rho,rhox,r,ril,rir)/vol!+rho*acc(i,1,1)*dt
            u1(5)=u1(5)-m*gr*dt*gravity_source_energy(rho,vx,rhox,vxx,r,ril,rir)/vol!+u(2)*acc(i,1,1)*dt
            egv=u1(5)-half*u1(2)**2d0/u1(1)
            x=r/n_domain(1)/tempfloorscale
            if (egv<flooregv(u1(1),x)) then
                egv=flooregv(u1(1),x)
            end if
            u1(5)=egv+half*u1(2)**2d0/u1(1)
            blk%u(1:5,i,1,1)=u1
        end do
        !deallocate(acc)
    end if
end subroutine central_gravity

function gravity_source_mom_1d(rho0,rhox,x0,x1,x2)
    !second order source term rho=rho0+rhox(x-x0)
    real(8) :: rho0,rhox,x0,x1,x2,gravity_source_mom_1d
    gravity_source_mom_1d=-0.5*((x1 - x2)*(2*rho0 + rhox*(-2*x0 + x1 + x2)))
end function gravity_source_mom_1d

function gravity_source_energy_1d(rho0,v0,rhox,vx,x0,x1,x2)
    !second order source term rho=rho0+rhox(x-x0), v=v0+vx(x-x0)
    real(8) :: rho0,v0,rhox,vx,x0,x1,x2,gravity_source_energy_1d
    gravity_source_energy_1d=-0.5*((x1 - x2)*(2*rho0*v0 + rhox*v0*(-2*x0 + x1 + x2) + rho0*vx*(-2*x0 + x1 + x2)))
end function gravity_source_energy_1d

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

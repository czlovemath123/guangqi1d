module mathlib
use datastructure

implicit none

interface interpolation_linear
    module procedure interpolation1d_linear,interpolation2d_linear
end interface interpolation_linear

interface cspline
    module procedure cspline_single,cspline_array
end interface cspline

interface cspline_derivative1
    module procedure cspline_deri1_single,cspline_deri1_array
end interface cspline_derivative1

interface gradient
    module procedure gradient_scalar,gradient_vector
end interface gradient

type fun_array
    procedure(fun), pointer, nopass :: ptr
end type fun_array

real(8),parameter,private:: eps=epsilon(one)

character(len=*),parameter,public:: polyroots_version= "1.3 (4 jan 1999)"
integer,private:: outputcode
!    =0 degenerate equation
!    =1 one real root
!    =21 two identical real roots
!    =22 two distinct real roots
!    =23 two complex roots
!    =31 multiple real roots
!    =32 one real and two complex roots
!    =33 three distinct real roots
!    =41
!    =42 two real and two complex roots
!    =43
!    =44 four complex roots

private:: cuberoot
public:: linearroot
private:: onelargetwosmall
public:: quadraticroots
public:: cubicroots
public:: quarticroots
public:: solvepolynomial

!----------------------------------------------------------------------------


contains

function fun(x)
    real(8) :: fun
    real(8), dimension(:), allocatable :: x
end function fun

! ---------------------------------------------------------------------------
! purpose - solve for the roots of a polynomial equation with real
!   coefficients, up to quartic order. retrun a code indicating the nature
!   of the roots found.

! authors - alfred h. morris, naval surface weapons center, dahlgren,va
!           william l. davis, naval surface weapons center, dahlgren,va
!           alan miller,  csiro mathematical & information sciences
!                         clayton, victoria, australia 3169
!                         http://www.mel.dms.csiro.au/~alan
!           ralph l. carmichael, public domain aeronautical software
!                         http://www.pdas.com
!     revision history
!   date  vers person  statement of changes
!    ??    1.0 ahm&wld original coding                     
! 27feb97  1.1   am    converted to be compatible with elf90
! 12jul98  1.2   rlc   module format; numerous style changes
!  4jan99  1.3   rlc   made the tests for zero constant term exactly zero

function cuberoot(x) result(f)
! ---------------------------------------------------------------------------
! purpose - compute the cube root of a real(8) number. if the argument is
!   negative, then the cube root is also negative.

  real(8),intent(in) :: x
  real(8):: f
!----------------------------------------------------------------------------
  if (x < zero) then
    f=-exp(log(-x)/three)
  else if (x > zero) then
    f=exp(log(x)/three)
  else
    f=zero
  end if
  return
end function cuberoot   ! ---------------------------------------------------

!+
subroutine linearroot(a, z)
! ---------------------------------------------------------------------------
! purpose - computes the roots of the real polynomial
!              a(1) + a(2)*z 
!     and stores the results in z. it is assumed that a(2) is non-zero.
  real(8),intent(in),dimension(:):: a
  real(8),intent(out):: z
!----------------------------------------------------------------------------
  if (a(2)==0.0) then
    z=zero
  else
    z=-a(1)/a(2)
  end if
  return
end subroutine linearroot   ! -----------------------------------------------

!+
subroutine onelargetwosmall(a1,a2,a4,w, z)
! ---------------------------------------------------------------------------
! purpose - compute the roots of a cubic when one root, w, is known to be
!   much larger in magnitude than the other two

  real(8),intent(in):: a1,a2,a4
  real(8),intent(in):: w
  complex(8),intent(out),dimension(:):: z


  real(8),dimension(3):: aq
!----------------------------------------------------------------------------
  aq(1)=a1
  aq(2)=a2+a1/w
  aq(3)=-a4*w
  call quadraticroots(aq, z)
  z(3)=cmplx(w,zero,8)
  
  if (aimag(z(1)) == zero) return
  z(3)=z(2)
  z(2)=z(1)
  z(1)=cmplx(w,zero,8)
  return
end subroutine onelargetwosmall   ! -----------------------------------------

!+
subroutine quadraticroots(a, z)
! ---------------------------------------------------------------------------
! purpose - computes the roots of the real polynomial
!              a(1) + a(2)*z + a(3)*z**2
!     and stores the results in z.  it is assumed that a(3) is nonzero.

  real(8),intent(in),dimension(:):: a
  complex(8),intent(out),dimension(:):: z


  real(8):: d, r, w, x, y
!----------------------------------------------------------------------------
  if(a(1)==0.0) then     ! eps is a global module constant (private)
    z(1) = czero               ! one root is obviously zero
    z(2) = cmplx(-a(2)/a(3), zero,8)    ! remainder is a linear eq.
    outputcode=21   ! two identical real roots
    return
  end if

  d = a(2)*a(2) - four*a(1)*a(3)             ! the discriminant
  if (abs(d) <= two*eps*a(2)*a(2)) then
    z(1) = cmplx(-half*a(2)/a(3), zero, 8) ! discriminant is tiny
    z(2) = z(1)
    outputcode=22  ! two distinct real roots
    return
  end if

  r = sqrt(abs(d))
  if (d < zero) then
    x = -half*a(2)/a(3)        ! negative discriminant => roots are complex   
    y = abs(half*r/a(3))
    z(1) = cmplx(x, y, 8)
    z(2) = cmplx(x,-y, 8)   ! its conjugate
    outputcode=23                        !  complex roots
    return
  end if

  if (a(2) /= zero) then              ! see numerical recipes, sec. 5.5
    w = -(a(2) + sign(r,a(2)))
    z(1) = cmplx(two*a(1)/w,  zero, 8)
    z(2) = cmplx(half*w/a(3), zero, 8)
    outputcode=22           ! two real roots
    return
  end if

  x = abs(half*r/a(3))   ! a(2)=0 if you get here
  z(1) = cmplx( x, zero, 8)
  z(2) = cmplx(-x, zero, 8)
  outputcode=22
  return
end subroutine quadraticroots   ! -------------------------------------------

!+
subroutine cubicroots(a, z)
!----------------------------------------------------------------------------
! purpose - compute the roots of the real polynomial
!              a(1) + a(2)*z + a(3)*z**2 + a(4)*z**3
  real(8),intent(in),dimension(:):: a
  complex(8),intent(out),dimension(:):: z

  real(8),parameter:: rt3=1.7320508075689d0    ! (sqrt(3)
  real (8) :: aq(3), arg, c, cf, d, p, p1, q, q1
  real(8):: r, ra, rb, rq, rt
  real(8):: r1, s, sf, sq, sum, t, tol, t1, w
  real(8):: w1, w2, x, x1, x2, x3, y, y1, y2, y3

! note -   it is assumed that a(4) is non-zero. no test is made here.
!----------------------------------------------------------------------------
  if (a(1)==0.0) then
    z(1) = czero  ! one root is obviously zero
    call quadraticroots(a(2:4), z(2:3))   ! remaining 2 roots here
    return
  end if

  p = a(3)/(three*a(4))
  q = a(2)/a(4)
  r = a(1)/a(4)
  tol = four*eps

  c = zero
  t = a(2) - p*a(3)
  if (abs(t) > tol*abs(a(2))) c = t/a(4)

  t = two*p*p - q
  if (abs(t) <= tol*abs(q)) t = zero
  d = r + p*t
  if (abs(d) <= tol*abs(r)) go to 110

!           set  sq = (a(4)/s)**2 * (c**3/27 + d**2/4)

  s = max(abs(a(1)), abs(a(2)), abs(a(3)))
  p1 = a(3)/(three*s)
  q1 = a(2)/s
  r1 = a(1)/s

  t1 = q - 2.25d0*p*p
  if (abs(t1) <= tol*abs(q)) t1 = zero
  w = fourth*r1*r1
  w1 = half*p1*r1*t
  w2 = q1*q1*t1/27.0d0

  if (w1 >= zero) then
    w = w + w1
    sq = w + w2
  else if (w2 < zero) then
    sq = w + (w1 + w2)
  else
    w = w + w2
    sq = w + w1
  end if

  if (abs(sq) <= tol*w) sq = zero
  rq = abs(s/a(4))*sqrt(abs(sq))
  if (sq >= zero) go to 40

!                   all roots are real

  arg = atan2(rq, -half*d)
  cf = cos(arg/three)
  sf = sin(arg/three)
  rt = sqrt(-c/three)
  y1 = two*rt*cf
  y2 = -rt*(cf + rt3*sf)
  y3 = -(d/y1)/y2

  x1 = y1 - p
  x2 = y2 - p
  x3 = y3 - p

  if (abs(x1) > abs(x2)) call swap_single(x1,x2)
  if (abs(x2) > abs(x3)) call swap_single(x2,x3)
  if (abs(x1) > abs(x2)) call swap_single(x1,x2)

  w = x3

  if (abs(x2) < 0.1d0*abs(x3)) go to 70
  if (abs(x1) < 0.1d0*abs(x2)) x1 = - (r/x3)/x2
  z(1) = cmplx(x1, zero,8)
  z(2) = cmplx(x2, zero,8)
  z(3) = cmplx(x3, zero,8)
  return

!                  real and complex roots

40 ra =cuberoot(-half*d - sign(rq,d))
  rb = -c/(three*ra)
  t = ra + rb
  w = -p
  x = -p
  if (abs(t) <= tol*abs(ra)) go to 41
  w = t - p
  x = -half*t - p
  if (abs(x) <= tol*abs(p)) x = zero
  41 t = abs(ra - rb)
  y = half*rt3*t
  
  if (t <= tol*abs(ra)) go to 60
  if (abs(x) < abs(y)) go to 50
  s = abs(x)
  t = y/x
  go to 51
50 s = abs(y)
  t = x/y
51 if (s < 0.1d0*abs(w)) go to 70
  w1 = w/s
  sum = one + t*t
  if (w1*w1 < 0.01d0*sum) w = - ((r/sum)/s)/s
  z(1) = cmplx(w, zero,8)
  z(2) = cmplx(x, y,8)
  z(3) = cmplx(x,-y,8)
  return

!               at least two roots are equal

60 if (abs(x) < abs(w)) go to 61
  if (abs(w) < 0.1d0*abs(x)) w = - (r/x)/x
  z(1) = cmplx(w, zero,8)
  z(2) = cmplx(x, zero,8)
  z(3) = z(2)
  return
  61 if (abs(x) < 0.1d0*abs(w)) go to 70
  z(1) = cmplx(x, zero,8)
  z(2) = z(1)
  z(3) = cmplx(w, zero,8)
  return

!     here w is much larger in magnitude than the other roots.
!     as a result, the other roots may be exceedingly inaccurate
!     because of roundoff error.  to deal with this, a quadratic
!     is formed whose roots are the same as the smaller roots of
!     the cubic.  this quadratic is then solved.

!     this code was written by william l. davis (nswc).

70 aq(1) = a(1)
  aq(2) = a(2) + a(1)/w
  aq(3) = -a(4)*w
  call quadraticroots(aq, z)
  z(3) = cmplx(w, zero,8)
  
  if (aimag(z(1)) == zero) return
  z(3) = z(2)
  z(2) = z(1)
  z(1) = cmplx(w, zero,8)
  return
!-----------------------------------------------------------------------


!                   case when d = 0

110 z(1) = cmplx(-p, zero,8)
  w = sqrt(abs(c))
  if (c < zero) go to 120
  z(2) = cmplx(-p, w,8)
  z(3) = cmplx(-p,-w,8)
  return

120 if (p /= zero) go to 130
  z(2) = cmplx(w, zero,8)
  z(3) = cmplx(-w, zero,8)
  return

130 x = -(p + sign(w,p))
  z(3) = cmplx(x, zero,8)
  t = three*a(1)/(a(3)*x)
  if (abs(p) > abs(t)) go to 131
  z(2) = cmplx(t, zero,8)
  return
131 z(2) = z(1)
  z(1) = cmplx(t, zero,8)
  return
end subroutine cubicroots   ! -----------------------------------------------


!+
subroutine quarticroots(a,z)
!----------------------------------------------------------------------------
! purpose - compute the roots of the real polynomial
!               a(1) + a(2)*z + ... + a(5)*z**4

  real(8), intent(in)     :: a(5)
  complex(8), intent(out) :: z(4)

  complex(8) :: w
  real(8):: b,b2, c, d, e, h, p, q, r, t
  real(8),dimension(4):: temp
  real(8):: u, v, v1, v2, x, x1, x2, x3, y


! note - it is assumed that a(5) is non-zero. no test is made here

!----------------------------------------------------------------------------

  if (a(1)==0.0) then
    z(1) = czero    !  one root is obviously zero
    call cubicroots(a(2:), z(2:))
    return
  end if


  b = a(4)/(four*a(5))
  c = a(3)/a(5)
  d = a(2)/a(5)
  e = a(1)/a(5)
  b2 = b*b

  p = half*(c - 6.0d0*b2)
  q = d - two*b*(c - four*b2)
  r = b2*(c - three*b2) - b*d + e

! solve the resolvent cubic equation. the cubic has at least one
! nonnegative real root.  if w1, w2, w3 are the roots of the cubic
! then the roots of the originial equation are
!     z = -b + csqrt(w1) + csqrt(w2) + csqrt(w3)
! where the signs of the square roots are chosen so
! that csqrt(w1) * csqrt(w2) * csqrt(w3) = -q/8.

  temp(1) = -q*q/64.0d0
  temp(2) = 0.25d0*(p*p - r)
  temp(3) =  p
  temp(4) = one
  call cubicroots(temp,z)
  if (aimag(z(2)) /= zero) go to 60

!         the resolvent cubic has only real roots
!         reorder the roots in increasing order

  x1 = dble(z(1))
  x2 = dble(z(2))
  x3 = dble(z(3))
  if (x1 > x2) call swap_single(x1,x2)
  if (x2 > x3) call swap_single(x2,x3)
  if (x1 > x2) call swap_single(x1,x2)

  u = zero
  if (x3 > zero) u = sqrt(x3)
  if (x2 <= zero) go to 41
  if (x1 >= zero) go to 30
  if (abs(x1) > x2) go to 40
  x1 = zero

30 x1 = sqrt(x1)
  x2 = sqrt(x2)
  if (q > zero) x1 = -x1
  temp(1) = (( x1 + x2) + u) - b
  temp(2) = ((-x1 - x2) + u) - b
  temp(3) = (( x1 - x2) - u) - b
  temp(4) = ((-x1 + x2) - u) - b
  call selectsort_poly(temp)
  if (abs(temp(1)) >= 0.1d0*abs(temp(4))) go to 31
  t = temp(2)*temp(3)*temp(4)
  if (t /= zero) temp(1) = e/t
31 z(1) = cmplx(temp(1), zero,8)
  z(2) = cmplx(temp(2), zero,8)
  z(3) = cmplx(temp(3), zero,8)
  z(4) = cmplx(temp(4), zero,8)
  return

40 v1 = sqrt(abs(x1))
v2 = zero
go to 50
41 v1 = sqrt(abs(x1))
v2 = sqrt(abs(x2))
if (q < zero) u = -u

50 x = -u - b
y = v1 - v2
z(1) = cmplx(x, y,8)
z(2) = cmplx(x,-y,8)
x =  u - b
y = v1 + v2
z(3) = cmplx(x, y,8)
z(4) = cmplx(x,-y,8)
return

!                the resolvent cubic has complex roots

60 t = dble(z(1))
x = zero
if (t < zero) then
  go to 61
else if (t == zero) then
  go to 70
else
  go to 62
end if
61 h = abs(dble(z(2))) + abs(aimag(z(2)))
if (abs(t) <= h) go to 70
go to 80
62 x = sqrt(t)
if (q > zero) x = -x

70 w = sqrt(z(2))
  u = two*dble(w)
  v = two*abs(aimag(w))
  t =  x - b
  x1 = t + u
  x2 = t - u
  if (abs(x1) <= abs(x2)) go to 71
  t = x1
  x1 = x2
  x2 = t
71 u = -x - b
  h = u*u + v*v
  if (x1*x1 < 0.01d0*min(x2*x2,h)) x1 = e/(x2*h)
  z(1) = cmplx(x1, zero,8)
  z(2) = cmplx(x2, zero,8)
  z(3) = cmplx(u, v,8)
  z(4) = cmplx(u,-v,8)
  return

80 v = sqrt(abs(t))
  z(1) = cmplx(-b, v,8)
  z(2) = cmplx(-b,-v,8)
  z(3) = z(1)
  z(4) = z(2)
  return

end subroutine quarticroots

!+
subroutine selectsort_poly(a)
! ---------------------------------------------------------------------------
! purpose - reorder the elements of in increasing order.
  real(8),intent(in out),dimension(:):: a

  integer:: j
  integer,dimension(1):: k
! note - this is a n**2 method. it should only be used for small arrays. <25
!----------------------------------------------------------------------------
  do j=1,size(a)-1
    k=minloc(a(j:))
    if (j /= k(1)) call swap_single(a(k(1)),a(j))
  end do
  return
end subroutine selectsort_poly   ! -----------------------------------------------

!+
subroutine solvepolynomial(quarticcoeff, cubiccoeff, quadraticcoeff, &
  linearcoeff, constantcoeff, code, root1,root2,root3,root4)
! ---------------------------------------------------------------------------
  real(8),intent(in):: quarticcoeff
  real(8),intent(in):: cubiccoeff, quadraticcoeff
  real(8),intent(in):: linearcoeff, constantcoeff
  integer,intent(out):: code
  complex(8),intent(out):: root1,root2,root3,root4
  real(8),dimension(5):: a
  complex(8),dimension(5):: z
!----------------------------------------------------------------------------
  a(1)=constantcoeff
  a(2)=linearcoeff
  a(3)=quadraticcoeff
  a(4)=cubiccoeff
  a(5)=quarticcoeff

  if (quarticcoeff /= zero) then
    call quarticroots(a,z)  
  else if (cubiccoeff /= zero) then
    call cubicroots(a,z)
  else if (quadraticcoeff /= zero) then
    call quadraticroots(a,z)
  else if (linearcoeff /= zero) then
    z(1)=cmplx(-constantcoeff/linearcoeff, 0, 8)
    outputcode=1
  else
    outputcode=0    !  { no roots }
  end if

  code=outputcode
  if (outputcode > 0) root1=z(1)
  if (outputcode > 1) root2=z(2)
  if (outputcode > 23) root3=z(3)
  if (outputcode > 99) root4=z(4)
  return
end subroutine solvepolynomial   ! ------------------------------------------

subroutine quadratic_real(a,b,c,x,code)
    !solve ax^2+bx+c=0, a cannot be zero, discriminant cannot be too small
    !code=0 no real root, 1 two identical real roots, 2 two distinct real roots
    !x are the roots
    real(8) :: a,b,c,x(2),delta
    integer :: code
    if (a.eq.0d0) then
        print *, 'quadratic_real, a=0'
        stop
    end if
    !delta is the discriminant
    delta=b**2-4d0*a*c
    if (delta.lt.0d0) then  !complex roots
        code=0
    else if (delta.eq.0d0) then
        code=1
        x(1)=-b/2d0/a
        x(2)=-b/2d0/a
    else if (delta.gt.0d0) then
        code=2
        if (a.gt.0d0) then
            x(1)=(-b+sqrt(delta))/2d0/a
            x(2)=(-b-sqrt(delta))/2d0/a  !x(1)>x(2)
        else
            x(1)=(-b-sqrt(delta))/2d0/a
            x(2)=(-b+sqrt(delta))/2d0/a  !x(1)>x(2)
        end if
    else
        print *, 'quadratic_real wrong'
        stop
    end if
end subroutine quadratic_real

!***********************************end polynomial root***********************

!****************************vector and tensor operator***********************

subroutine divergence1d(input,div_input,dx)  !center finite difference
    !input(0:n+1), div_input(0:n+1)
    real(8), dimension(:), allocatable :: input,div_input
    integer :: i,n
    real(8) :: dx
    n=size(div_input)-2
    do i=1,n
        div_input(i)=(input(i+1)-input(i-1))/2/dx
    end do
end subroutine divergence1d

subroutine divergence1d_spherical(input,div_input,dr,r_bound)
    !div in spherical coordinate is radii dependent
    !input(0:n+1), div_input(0:n+1)
    real(8), dimension(:), allocatable :: input,div_input
    integer :: i,n
    real(8) :: r_bound,r_in,r_out,r_center,dr,v_left,v_right
    !r_bound is the inner radii of div_input
    n=size(div_input)-2
    do i=1,n
        r_in=r_bound+(i-1)*dr
        r_out=r_bound+i*dr
        r_center=r_bound+(i-half)*dr
        v_left=(input(i-1)+input(i))/2
        v_right=(input(i)+input(i+1))/2
        div_input(i)=(1d0/r_center/r_center)*(r_out*r_out*v_left-r_in*r_in*v_right)/dr
    end do
end subroutine divergence1d_spherical

!************************end vector and tensor operator***********************



!******************************************************************************************
!******************************************************************************************
!Begin general root finding section
!******************************************************************************************
!******************************************************************************************

subroutine newtonraphson(ptr1,ptr2,root,estimate,tolerance)
    !newton-raphson method to find root with an array of parameters
    !ptr1 will point to the function
    !ptr2 will point to the derivative of the function
    procedure(fun), pointer :: ptr1,ptr2
    real(8) :: root,tolerance,diff
    real(8), dimension(:), allocatable :: xold,xnew,estimate
    integer :: n,i
    n=size(estimate)
    allocate(xold(n),xnew(n))
    xold=estimate
    xnew=estimate
    xnew(1)=xold(1)-ptr1(xold)/ptr2(xold)
    diff=abs(xnew(1)-xold(1))
    i=0
    do while (diff>tolerance)
        xold(1)=xnew(1)
        xnew(1)=xold(1)-ptr1(xold)/ptr2(xold)
        diff=abs(xnew(1)-xold(1))/abs(xold(1))
        i=i+1
        if (i>100) then
            print *,'newtonraphson too many iterations'
            stop
        end if
    end do
    root=xnew(1)
    deallocate(xold,xnew)
end subroutine newtonraphson

subroutine bisectionroot(ptr,a,b,x,converge,root,conv_mode)
    real(8), dimension(:), allocatable :: x,xleft,xright,lastroot,root
    !x is the array to be fed to ptr function, x(1) will be the variable
    !root is the array of solution
    real(8) :: last,a,b,converge,temp,fa,fb,eps
    !a and b are the min and max of the possible variable, b>a absolutely
    procedure (fun), pointer :: ptr
    integer :: slope,i,n
    integer, optional :: conv_mode  !1=absolute conv  2=relative conv
    eps=2*epsilon(1d0)
    n=size(x)
    allocate(xleft(n),xright(n),lastroot(n))
    xleft=x
    xright=x
    xleft(1)=a
    xright(1)=b
    fa=ptr(xleft)
    fb=ptr(xright)
    i=0
    !check for trivial case
    if (abs(fa)<eps) then
        root=xleft
        goto 10
    end if
    if (abs(fb)<eps) then
        root=xright
        goto 10
    end if
    if (fa*fb>0) then
        write (*,*) "bisectionroot: wrong (a,b)"
        print *,a,b
        print *,x,ptr(xleft),ptr(xright)
        stop
    else if (ptr(xleft)>ptr(xright)) then
        slope=1
    else
        slope=-1
    end if
    lastroot=xleft
    last=fa
    root=x   !x(1) is the initial guess of the root
    temp=ptr(root)
    if (present(conv_mode).and.conv_mode==2) then   !relative convergence criterion
        do while (2d0*(xright(1)-xleft(1))/(xright(1)+xleft(1))>converge)
            lastroot=root
            if (temp*slope.lt.0) then
                xright=root
            else
                xleft=root
            end if
            last=temp
            root(1)=(xleft(1)+xright(1))/2d0
            temp=ptr(root)
            i=i+1
            if (i>100) then
                write(*,*) "bisection too many iterations"
                stop
            end if
            if ((xleft(1)+xright(1)).eq.0d0) then
                print *,'bisection conv_mode wrong'
                stop
            end if
        end do
    else   !absolute convergence criterion
        do while ((xright(1)-xleft(1))>converge)
            lastroot=root
            if (temp*slope.lt.0) then
                xright=root
            else
                xleft=root
            end if
            last=temp
            root(1)=(xleft(1)+xright(1))/2d0
            temp=ptr(root)
            i=i+1
            if (i>100) then
                write(*,*) "bisection too many iterations"
                stop
            end if
        end do
    end if
10  deallocate(xleft,xright,lastroot)
end subroutine bisectionroot

subroutine brent_method(ptr,a,b,x,toll,delta,root)  !need to be reprogramed
    real(8), dimension(:), allocatable :: x,xa,xb,xc,xs,root
    !x is the array to be fed to ptr function, x(1) will be the variable
    !root is the array of solution
    real(8) :: a,b,c,d,s,toll,delta,fa,fb,fc,fs
    !a and b are the min and max of the possible variable, b>a absolutely
    procedure (fun), pointer :: ptr
    integer :: slope,i,n
    logical :: flag,flag1,flag2,flag3,flag4,flag5
    d=0d0
    i=0
    n=size(x)
    allocate(xa(n),xb(n),xc(n),xs(n))
    xa=x
    xb=x
    xc=x
    xs=x
    xa(1)=a
    xb(1)=b
    fa=ptr(xa)
    fb=ptr(xb)
    if (fa.eq.0d0) then
        root=xa
        goto 20
    end if
    if (fb.eq.0d0) then
        root=xb
        goto 20
    end if
    if (fa*fb.gt.0d0) then
        print *,'brent_method: solution not bracketed'
        stop
    end if
    if (abs(fa).lt.abs(fb)) then
        call swap_array(xa,xb)
        call swap_single(fa,fb)
    end if
    xc=xa
    flag=.true.
    do while ((abs(xb(1)-xa(1)).gt.toll).and.(abs(fb).gt.delta))
        if (xc(1).ne.xa(1).and.xc(1).ne.xb(1)) then  !inverse quadratic interpolation
            !print *, 'a',fa,fb,fc,xa(1),xb(1),xc(1)
            if (fb.eq.fc) then
                root=xb
                goto 20
            else
                s=xa(1)*fb*fc/(fa-fb)/(fa-fc)+xb(1)*fa*fc/(fb-fa)/(fb-fc)+xc(1)*fa*fb/(fc-fa)/(fc-fb)
            end if
        else   !secant method
            !print *, 'b',fa,fb,fc,xa(1),xb(1),xc(1)
            s=xb(1)-fb*(xb(1)-xa(1))/(fb-fa)
        end if
        if (isnan(s)) then
            print *, '1',fa,fb,fc,xa(1),xb(1),xc(1)
            stop
        end if
        flag1=(s-(3d0*xa(1)+xb(1))/4d0)*(s-xb(1)).gt.0d0
        flag2=(flag.eqv..true.).and.(abs(s-xb(1)).ge.abs(xb(1)-xc(1))/2d0)
        flag3=(flag.eqv..false.).and.(abs(s-xb(1)).ge.abs(xc(1)-d)/2d0)
        flag4=(flag.eqv..true.).and.(abs(xb(1)-xc(1)).lt.toll)
        flag5=(flag.eqv..false.).and.(abs(xc(1)-d).lt.toll)
        !print *,flag,flag1,flag2,flag3,flag4,flag5
        if (flag1.or.flag2.or.flag3.or.flag4.or.flag5) then   !bisection method
            s=(xa(1)+xb(1))/2d0
            flag=.true.
        else
            flag=.false.
        end if
        if (isnan(s)) then
            print *, '2'
            stop
        end if
10      xs(1)=s
        fs=ptr(xs)
        fa=ptr(xa)
        d=xc(1)
        xc=xb
        if (fa*fs.le.0d0) then
            xb(1)=s
        else
            xa(1)=s
        end if
        fa=ptr(xa)
        fb=ptr(xb)
        fc=ptr(xc)
        if (abs(fa).lt.abs(fb)) then
            call swap_array(xa,xb)
            call swap_single(fa,fb)
        end if
        i=i+1
    end do
    root=xb
20  deallocate(xa,xb,xc,xs)
end subroutine brent_method

!******************************************************************************************
!******************************************************************************************
!End general root finding section
!******************************************************************************************
!******************************************************************************************

subroutine swap_single(a,b)
    !exchange the value of a and b
    real(8) :: a,b,c
    c=a
    a=b
    b=c
end subroutine swap_single

subroutine swap_array(a,b)
    !exchange array a and array b
    real(8), dimension(:), allocatable :: a,b,c
    integer :: n
    n=size(a)
    allocate(c(n))
    c=a
    a=b
    b=c
    deallocate(c)
end subroutine swap_array

function pow(a,b)
    real(8) :: a,b,pow
    pow=a**b
end function pow

!******************************************************************************************
!******************************************************************************************
!Begin interpolation section
!******************************************************************************************
!******************************************************************************************

subroutine loglogloginversefunction(x,y,z,s,t,bound)
    ! given x,y coord and s(x,y) table, calculate the inverse function x(s,y)
    ! the result will be stroed in t, as t(s,y)
    ! x,y should all be in its log form and increase monotonically
    ! s is not in log form
    ! t will be log10log10log10
    ! z will be calculated, increase monotonically
    ! bound(1) is the value assigned when the t value is below the lower bound
    ! bound(2) is the value assigned when the t value is above the higher bound
    real(8), dimension(:), allocatable :: x,y,z,z_extend
    real(8), dimension(:,:), allocatable :: s,t,s_log
    real(8) :: zmax,zmin,dlogz,zlogspan,xleft,dx
    real(8), optional :: bound(2)
    integer :: i,j,k,mark,m,n,pts
    n=size(x)
    m=size(y)
    pts=size(t,1)
    allocate(s_log(n,m),z_extend(0:pts+1))
    s_log=log10(s)
    zmax=maxval(s_log)
    zmin=minval(s_log)
    !range and domain of the first variable in x(s,y)
    zlogspan=zmax-zmin
    dlogz=zlogspan/(pts-1)
    do i=1,pts
        z(i)=zmin+(i-1)*dlogz
    end do
    do i=0,pts+1
        z_extend(i)=zmin+(i-1)*dlogz
    end do
    do j=1,m
        do i=1,pts
            do k=1,n-1
                dx=x(k+1)-x(k)
                if (z_extend(i).le.s_log(k+1,j).and.z_extend(i).ge.s_log(k,j)) then
                    xleft=x(k)
                    mark=k
                    t(i,j)=xleft+dx*(z_extend(i)-s_log(mark,j))/(s_log(mark+1,j)-s_log(mark,j))
                    exit
                else if (z_extend(i).ge.s_log(n,j).and.z_extend(i-1).le.s_log(n,j)) then
                    xleft=x(n)
                    mark=n
                    t(i,j)=xleft+dx*(z_extend(i)-s_log(mark,j))/(s_log(mark,j)-s_log(mark-1,j))
                    exit
                else if (z_extend(i).le.s_log(1,j).and.z_extend(i+1).ge.s_log(1,j)) then
                    xleft=x(1)
                    mark=1
                    t(i,j)=xleft+dx*(z_extend(i)-s_log(mark,j))/(s_log(mark+1,j)-s_log(mark,j))
                    exit
                else if (z_extend(i).lt.s_log(1,j)) then
                    ! at such y(j), this z_extend(i) is lower than the lowest possible x
                    if (present(bound)) then
                        t(i,j)=bound(1)
                    else
                        t(i,j)=x(1)-dx
                    end if
                else if (z_extend(i).gt.s_log(n,j)) then
                    ! at such y(j), this z_extend(i) is higher than the highest possible x
                    if (present(bound)) then
                        t(i,j)=bound(2)
                    else
                        t(i,j)=x(n)+dx
                    end if
                else
                end if
            end do
        end do
    end do
    deallocate(s_log,z_extend)
end subroutine loglogloginversefunction

subroutine interpolation1d_linear(x,val,xlist,table)
    !piecewise linear interpolation
    !x is the interpolated point
    integer :: i,x1,x2,xdim
    real(8) :: x,val
    real(8), allocatable :: xlist(:),table(:)
    xdim=size(xlist)
    do i=1,xdim
        if (x > xlist(i) .AND. x <= xlist(i+1)) then
            x1=i
            x2=i+1
            exit
        else if (x <= xlist(1)) then
            x1=1
            x2=2
            exit
        else if (x > xlist(xdim)) then
            x1=xdim-1
            x2=xdim
            exit
        end if
    end do
    val=(table(x1)*(xlist(x2)-x)+table(x2)*(x-xlist(x1)))/(xlist(x2)-xlist(x1))
end subroutine interpolation1d_linear

function quadratic_interpolation(x,xx,yy)
    real(8) :: quadratic_interpolation,x,xx(3),yy(3),l(3)
    l(1)=(x-xx(2))*(x-xx(3))/(xx(1)-xx(2))/(xx(1)-xx(3))
    l(2)=(x-xx(1))*(x-xx(3))/(xx(2)-xx(1))/(xx(2)-xx(3))
    l(3)=(x-xx(1))*(x-xx(2))/(xx(3)-xx(1))/(xx(3)-xx(2))
    quadratic_interpolation=dot_product(l,yy)
end function quadratic_interpolation

subroutine interpolation2d_linear(x,y,val,xlist,ylist,table)
    !bi-linear interpolation
    !x,y,(val) are the interpolated points (value)
    !xlist and ylist are from small to large
    integer :: i,x1,x2,y1,y2,xdim,ydim
    real(8) :: x,y,val
    real(8), allocatable :: table(:,:),xlist(:),ylist(:)
    xdim=size(xlist)
    ydim=size(ylist)
    do i=1,xdim
        if (x > xlist(i) .AND. x <= xlist(i+1)) then
            x1=i
            x2=i+1
            exit
        else if (x <= xlist(1)) then  !if x is smaller than the smallest
            x1=1
            x2=2
            exit
        else if (x > xlist(xdim)) then  !if x is greater than the largest
            x1=xdim-1
            x2=xdim
            exit
        end if
    end do
    do i=1,ydim
        if (y > ylist(i) .AND. y <= ylist(i+1)) then
            y1=i
            y2=i+1
            exit
        else if (y <= ylist(1)) then  !if y is smaller than the smallest
            y1=1
            y2=2
            exit
        else if (y > ylist(ydim)) then  !if y is greater than the largest
            y1=ydim-1
            y2=ydim
            exit
        end if
    end do
    val=(table(x1,y1)*(xlist(x2)-x)*(ylist(y2)-y)+&
    table(x2,y1)*(x-xlist(x1))*(ylist(y2)-y)+&
    table(x1,y2)*(xlist(x2)-x)*(y-ylist(y1))+&
    table(x2,y2)*(x-xlist(x1))*(y-ylist(y1)))/&
    ((xlist(x2)-xlist(x1))*(ylist(y2)-ylist(y1)))
end subroutine interpolation2d_linear

subroutine rect_bilinear_interp(p1,p2,p3,p4,p,v)
    real(8) :: p1(3),p2(3),p3(3),p4(3),p(2),v,x,y,x1,x2,y1,y2
    x=p(1);y=p(2);x1=p1(1);x2=p2(1);y1=p1(2);y2=p4(2)
    v=(x2-x)*(y2-y)/(x2-x1)/(y2-y1)*p1(3)+(x2-x)*(y-y1)/(x2-x1)/(y2-y1)*p2(3)   &
        +(x-x1)*(y2-y)/(x2-x1)/(y2-y1)*p3(3)+(x-x1)*(y-y1)/(x2-x1)/(y2-y1)*p4(3)
end subroutine rect_bilinear_interp

function unitsq_bilinear_interp(p,x,y)
    real(8) :: x,y,p(4),unitsq_bilinear_interp
    unitsq_bilinear_interp=(1d0-x)*(1d0-y)*p(1)+x*(1d0-y)*p(2)+(1d0-x)*y*p(3)+x*y*p(4)
end function unitsq_bilinear_interp

subroutine map_quadrilateral_to_unitsq(p1,p2,p3,p4,p,psq)
    real(8), dimension(2) :: p1,p2,p3,p4,p,psq
    real(8) :: a(4),b(4),aa,bb,cc,m,l,det,m1,m2
    a(1)=p1(1)
    a(2)=-p1(1)+p2(1)
    a(3)=-p1(1)+p3(1)
    a(4)=p1(1)-p2(1)-p3(1)+p4(1)
    b(1)=p1(2)
    b(2)=-p1(2)+p2(2)
    b(3)=-p1(2)+p3(2)
    b(4)=p1(2)-p2(2)-p3(2)+p4(2)
    aa=a(4)*b(3)-a(3)*b(4)
    bb=a(4)*b(1)-a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+p(1)*b(4)-p(2)*a(4)
    cc=a(2)*b(1)-a(1)*b(2)+p(1)*b(2)-p(2)*a(2)
    if (aa==0) then
        if (bb==0) then
            print *,'quadrilateral problem'
            stop
        else
            m=-cc/bb
        end if
    else
        det=sqrt(bb*bb-4d0*aa*cc)
        m1=(-bb+det)/2d0/aa
        m2=(-bb-det)/2d0/aa
        if (m1>0.and.m1<1) then
            m=m1
        else
            m=m2
        end if
    end if
    l=(p(1)-a(1)-a(3)*m)/(a(2)+a(4)*m)
    psq(1)=l;psq(2)=m
end subroutine map_quadrilateral_to_unitsq

subroutine interpolation2d_nonuniform(x,y,xlist,ylist,table)
    integer :: i,x1,x2,y1,y2,xdim,ydim
    real(8) :: x,y,val
    real(8), allocatable :: table(:,:),xlist(:),ylist(:)
    xdim=size(xlist)
    ydim=size(ylist)
end subroutine interpolation2d_nonuniform

subroutine cspline_integ_array(x,y,t,x_cspline,v,v0)
    !calculate v(i)=int_{x(1)}^{x(i)} ydx+v0 with cubic spline
    real(8), dimension(:), allocatable :: x,y,t,x_cspline,v,v_temp
    real(8), optional :: v0
    real(8) :: v_init
    integer :: i,j,n,m,n_total,k
    if (present(v0)) then
        v_init=v0
    else
        v_init=0d0
    end if
    n=size(x)
    n_total=size(x_cspline)
    m=(n_total-1)/(n-1)
    allocate(v_temp(n))
    v_temp=0d0
    v_temp(1)=v_init
    do i=1,n-1
        do j=1,m+1
            k=(i-1)*m+j
            call cspline_integ_single(x(i),x(i+1),y(i),y(i+1),t(i),t(i+1),  &
                x_cspline(k),v(k),v_temp(i))
        end do
        v_temp(i+1)=v(i*m+1)
    end do
    deallocate(v_temp)
end subroutine cspline_integ_array

subroutine cspline_integ_single(x1,x2,y1,y2,t1,t2,x,v,v0)
    !calculate v=int_{x(1)}^{x} ydx+v0 with cubic spline
    real(8) :: x1,x2,y1,y2,t1,t2,x,z,g,v,h00,h01,h10,h11,v_init
    real(8), optional :: v0
    if (present(v0)) then
        v_init=v0
    else
        v_init=0d0
    end if
    z=(x-x1)/(x2-x1)
    g=x-x1
    h00=z**3*g/2d0-z**2*g+g
    h10=z**3*g/4d0-2d0*z**2*g/3d0+z*g/2d0
    h01=-z**3*g/2d0+z**2*g
    h11=z**3*g/4d0-z**2*g/3d0
    print *,h00,h01,h00+h01
    v=h00*y1+h10*(x2-x1)*t1+h01*y2+h11*(x2-x1)*t2+v_init
end subroutine cspline_integ_single

subroutine cspline_deri1_array(x,y,t,x_cspline,v)
    real(8), dimension(:), allocatable :: x,y,t,x_cspline,v
    integer :: i,j,n,m,n_total,k
    n=size(x)
    n_total=size(x_cspline)
    m=(n_total-1)/(n-1)
    do i=1,n-1
        do j=1,m
            k=(i-1)*m+j
            call cspline_deri1_single(x(i),x(i+1),y(i),y(i+1),t(i),t(i+1),  &
                x_cspline(k),v(k))
        end do
    end do
    v(n_total)=t(n)
end subroutine cspline_deri1_array

subroutine cspline_deri1_single(x1,x2,y1,y2,t1,t2,x,v)
    !the first order derivative of cubic spline
    !x1<x<x2, find the interpolated value of v
    !see subroutine tangent for the meaning of x, y, and t
    real(8) :: x1,x2,y1,y2,t1,t2,x,z,v,h00,h01,h10,h11,dzdx
    z=(x-x1)/(x2-x1)
    dzdx=1d0/(x2-x1)
    h00=6d0*z**2-6d0*z
    h10=3d0*z**2-4d0*z+1d0
    h01=-6d0*z**2+6d0*z
    h11=3d0*z**2-2d0*z
    v=(h00*y1+h10*(x2-x1)*t1+h01*y2+h11*(x2-x1)*t2)*dzdx
end subroutine cspline_deri1_single

subroutine cspline_init(x,n_insert,x_spline,u,v,w,p,q)
    !uniformly insert n_insert points between every two points in x
    !the inserted points (including points in x) are stored in x_spline
    !u,v,w,p,q are the arrays that are being created at the same time
    !x does not need to be uniform
    real(8), allocatable :: x(:),x_spline(:),u(:)
    real(8), allocatable, optional :: v(:),w(:),p(:),q(:)
    real(8) :: dx
    integer :: i,j,n,n_insert,n_total
    n=size(x)
    if (allocated(x_spline)) then
        n_total=size(x_spline)
        allocate(u(n_total))
    else
        n_total=(n-1)*(n_insert+1)+1
        allocate(x_spline(n_total),u(n_total))
    end if
    if (present(v)) allocate(v(n_total))
    if (present(w)) allocate(w(n_total))
    if (present(p)) allocate(p(n_total))
    if (present(q)) allocate(q(n_total))
    do i=1,n-1
        dx=(x(i+1)-x(i))/(n_insert+1)
        do j=1,n_insert+1
            x_spline((i-1)*(n_insert+1)+j)=x(i)+(j-1)*dx
        end do
    end do
    x_spline(n_total)=x(n)
end subroutine cspline_init

subroutine cspline_array(x,y,t,x_cspline,v)
    !cubic spline for an array
    !points in x_cspline are uniformly inserted
    !x,y,t are the data points value
    !x_cspline contains the points to be interpolated
    !v is the interpolated value
    real(8), dimension(:), allocatable :: x,y,t,x_cspline,v
    integer :: i,j,n,m,n_total,k
    n=size(x)
    n_total=size(x_cspline)
    m=(n_total-1)/(n-1)
    do i=1,n-1
        do j=1,m
            k=(i-1)*m+j
            call cspline_single(x(i),x(i+1),y(i),y(i+1),t(i),t(i+1),  &
                x_cspline(k),v(k))
        end do
    end do
    v(n_total)=y(n)
end subroutine cspline_array

subroutine cspline_single(x1,x2,y1,y2,t1,t2,x,v)
    !cubic spline for a single value
    !x1<x<x2, find the interpolated value of v
    !see subroutine tangent for the meaning of x, y, and t
    real(8) :: x1,x2,y1,y2,t1,t2,x,z,v,h00,h01,h10,h11
    if (x.eq.x1) then
        v=y1
    else if (x.eq.x2) then
        v=y2
    else
        z=(x-x1)/(x2-x1)
        h00=2d0*z**3-3d0*z**2+1d0
        h10=z**3-2d0*z**2+z
        h01=-2d0*z**3+3d0*z**2
        h11=z**3-z**2
        v=h00*y1+h10*(x2-x1)*t1+h01*y2+h11*(x2-x1)*t2
    end if
end subroutine cspline_single

subroutine cspline_five_points(v,v_interp,dx)
    !uniform five points
    !**********     v1     v2   v_interp(1)   v3   v_interp(2)   v4    v5     **************
    real(8) :: v(5),v_interp(2),v_deri(3),dx,x0,x
    x0=0d0
    x=dx/2
    v_deri(1)=(v(3)-v(1))/2/dx
    v_deri(2)=(v(4)-v(2))/2/dx
    v_deri(3)=(v(5)-v(3))/2/dx
    call cspline_single(x0,dx,v(2),v(3),v_deri(1),v_deri(2),x,v_interp(1))
    call cspline_single(x0,dx,v(3),v(4),v_deri(2),v_deri(3),x,v_interp(2))
end subroutine cspline_five_points

subroutine monotone_tangent(x,y,t)
    !given points (x,y), calculate the monotone tangent of dy/dx at x
    !and store in t, delta is the slope between points
    !x does not need to be uniform
    real(8), dimension(:), allocatable :: x,y,t,delta,alpha,beta
    real(8) :: tau
    integer :: i,n
    n=size(x)
    allocate(alpha(n-1),beta(n-1))
    alpha=0d0
    beta=0d0
    call slope(x,y,delta)
    call tangent(x,y,t)
    i=1
    do i=1,n-1
        if (delta(i).eq.0d0) then
            t(i)=0d0
            t(i+1)=0d0
        else
            alpha(i)=t(i)/delta(i)
            beta(i)=t(i+1)/delta(i)
            if (i.ne.1) then
                if ((alpha(i).lt.0d0).or.(beta(i-1).lt.0d0)) then
                    t(i)=0d0
                else if ((alpha(i)**2+beta(i)**2).gt.9d0) then
                    tau=3d0/(alpha(i)**2+beta(i)**2)
                    t(i)=tau*alpha(i)*delta(i)
                    t(i+1)=tau*beta(i)*delta(i)
                end if
            end if
        end if
    end do
    deallocate(alpha,beta)
end subroutine monotone_tangent

subroutine slope(x,y,delta)
    !given points (x,y), calculate the slope between points and store in delta
    !x does not need to be uniform
    real(8), dimension(:), allocatable :: x,y,delta
    integer :: i,n
    n=size(x)
    allocate(delta(n-1))
    do i=1,n-1
        delta(i)=(y(i+1)-y(i))/(x(i+1)-x(i))
    end do
end subroutine slope

subroutine tangent(x,y,t)
    !given points (x,y), calculate the tangent of dy/dx at x and store in t
    !x does not need to be uniform
    real(8), dimension(:), allocatable :: x,y,t
    integer :: i,n
    n=size(x)
    allocate(t(n))
    t(1)=(y(2)-y(1))/(x(2)-x(1))
    t(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
    do i=2,n-1
        t(i)=(y(i+1)-y(i-1))/(x(i+1)-x(i-1))
    end do
end subroutine tangent

!******************************************************************************************
!******************************************************************************************
!End interpolation section
!******************************************************************************************
!******************************************************************************************


!******************************************************************************************
!******************************************************************************************
!Begin integration
!******************************************************************************************
!******************************************************************************************

subroutine integrate_function(ptr,x,val,v_init)
    !calculate val=int_{x(1,:)}^{x(n,:)}ptr dx+v_init with piecewise linear method
    !x(:,1) is the variable, x(:,2,3,...) are the parameters
    procedure (fun), pointer :: ptr
    real(8), dimension(:,:), allocatable :: x
    real(8), dimension(:), allocatable :: x1,x2
    real(8) :: val,v0
    real(8), optional :: v_init
    integer :: n,m,i
    n=size(x,1)
    m=size(x,2)
    allocate(x1(m),x2(m))
    if (present(v_init)) then
        v0=v_init
    else
        v0=0d0
    end if
    val=v0
    do i=2,n
        x1=x(i-1,:)
        x2=x(i,:)
        val=val+(ptr(x1)+ptr(x2))/2d0*(x2(1)-x1(1))
    end do
    deallocate(x1,x2)
end subroutine integrate_function

subroutine integrate_profile_directive(x,y,v,v_init)
    !given coordinate x and value y, x does not need to be uniform
    !assume piecewise constant y
    !calculate v(i)=v_init+int_{x(1)}^{x(i)}y dx, v is also an array
    real(8), dimension(:), allocatable :: x,y,v
    real(8) :: v0
    real(8), optional :: v_init
    integer :: n,i
    n=size(x)
    if (present(v_init)) then
        v0=v_init
    else
        v0=0d0
    end if
    v(1)=v0
    do i=2,n
        v(i)=v(i-1)+(y(i)+y(i-1))*(x(i)-x(i-1))/2d0
    end do
end subroutine integrate_profile_directive

subroutine integrate_profile_distance(x,y,v,v_init)
    !given coordinate x and value y, x does not need to be uniform
    !assume piecewise constant y
    !calculate v(i)=v_init+int_{x(1)}^{x(i)}y abs(dx), v is also an array
    real(8), dimension(:), allocatable :: x,y,v
    real(8) :: v0
    real(8), optional :: v_init
    integer :: n,i
    n=size(x)
    if (present(v_init)) then
        v0=v_init
    else
        v0=0d0
    end if
    v(1)=v0
    do i=2,n
        v(i)=v(i-1)+(y(i)+y(i-1))*abs(x(i)-x(i-1))/2d0
    end do
end subroutine integrate_profile_distance

!******************************************************************************************
!******************************************************************************************
!End integration
!******************************************************************************************
!******************************************************************************************


!******************************************************************************************
!******************************************************************************************
!Begin vector calculus
!******************************************************************************************
!******************************************************************************************

subroutine gradient_scalar(u,gradu,dx)
    !u is a scalar field, uniform grid dx
    !gradu is a vector
    !central scheme
    real(8), allocatable, intent(in) :: u(:,:,:)
    real(8), intent(in) :: dx
    real(8), allocatable, intent(inout) :: gradu(:,:,:,:)
    integer :: i,j
    if (nd==1) then
        do i=1,nx
            gradu(i,1,1,1)=(u(i+1,1,1)-u(i-1,1,1))/2/dx
            gradu(i,1,1,2)=zero
            gradu(i,1,1,3)=zero
        end do
    else if (nd==2) then
        do j=1,ny
            do i=1,nx
                gradu(i,j,1,1)=(u(i+1,j,1)-u(i-1,j,1))/2/dx
                gradu(i,j,1,2)=(u(i,j+1,1)-u(i,j-1,1))/2/dx
                gradu(i,j,1,3)=zero
            end do
        end do
    else
    end if
end subroutine gradient_scalar

subroutine gradient_vector(u,gradu,dx)
    !u is a 3d vector field (1d of which is uniform), uniform grid dx
    !gradu is rank 2 tensor field that has been flattened
    !du1/dx1    du1/dx2    du1/dx3
    !du2/dx1    du2/dx2    du2/dx3
    !du3/dx1    du3/dx2    du3/dx3
    !central scheme
    real(8), allocatable, intent(in) :: u(:,:,:,:)
    real(8), intent(in) :: dx
    real(8), allocatable, intent(inout) :: gradu(:,:,:,:)
    real(8) :: local_tensor(3,3),array(9)
    integer :: i,j
    if (nd==1) then
        do i=1,nx
            local_tensor(1,1)=(u(i+1,1,1,1)-u(i-1,1,1,1))/2/dx
            local_tensor(2,1)=zero
            local_tensor(3,1)=zero
            local_tensor(1,2)=zero
            local_tensor(2,2)=zero
            local_tensor(3,2)=zero
            local_tensor(1,3)=zero
            local_tensor(2,3)=zero
            local_tensor(3,3)=zero
            array=reshape(local_tensor,(/9/))
            gradu(i,1,1,1:9)=array
        end do
    else if (nd==2) then
        do j=1,ny
            do i=1,nx
                !fortran is column major, each row of the tensor is a gradient of the corresponding scalar
                local_tensor(1,1)=(u(i+1,j,1,1)-u(i-1,j,1,1))/2/dx
                local_tensor(2,1)=(u(i+1,j,1,2)-u(i-1,j,1,2))/2/dx
                local_tensor(3,1)=zero
                local_tensor(1,2)=(u(i,j+1,1,1)-u(i,j-1,1,1))/2/dx
                local_tensor(2,2)=(u(i,j+1,1,2)-u(i,j-1,1,2))/2/dx
                local_tensor(3,2)=zero
                local_tensor(1,3)=zero
                local_tensor(2,3)=zero
                local_tensor(3,3)=zero
                array=reshape(local_tensor,(/9/))
                gradu(i,j,1,1:9)=array
            end do
        end do
    else
    end if
end subroutine gradient_vector

function slope_limiter(r)
    real(8) :: slope_limiter,omega,r,xi_l,xi_r,rnew
    if (r<=zero) then
        slope_limiter=zero
    else
        omega=0d0
        slope_limiter=min(2/(1+r),2*r/(1+r))
    end if
end function slope_limiter

function minmod(diff_left,diff_right)
    real(8) :: diff_left,diff_right,minmod
    if (diff_left*diff_right<=0) then
        minmod=0d0
    else
        minmod=abs(diff_left)/diff_left*min(abs(diff_left),abs(diff_right))
    end if
end function minmod

subroutine divergence(u,divu,dx)
    !u is a 3d vector field (1d of which is uniform), uniform grid dx
    !central scheme
    real(8), allocatable, intent(in) :: u(:,:,:,:)
    real(8), intent(in) :: dx
    real(8), allocatable, intent(out) :: divu(:,:,:)
    integer :: i,j
    do j=1,ny
        do i=1,nx
            divu(i,j,1)=(u(i+1,j,1,1)-u(i-1,j,1,1)+u(i,j+1,1,2)-u(i,j-1,1,2))/dx/2
        end do
    end do
end subroutine divergence

!******************************************************************************************
!******************************************************************************************
!End vector calculus
!******************************************************************************************
!******************************************************************************************

!******************************************************************************************
!******************************************************************************************
!Begin ODE and PDE section
!******************************************************************************************
!******************************************************************************************

subroutine advection_1d(f,v,x,dt,output,ncores,geometry,flux)
    !solve the equation df/dt+div(fv)=0 where v is piecewise constant
    !0,1,2 geometry, uniform grid, upwind
    !f and v are dependent variables, f should include boundary
    !f of next time (dt) is stored in output
    !flux stores the flux due to advection, (erg/cm^2/s)
    real(8), dimension(:), allocatable :: f,v,output,x,deltax,f_next,flux
    real(8) :: dt,dx
    integer :: xdim,n,i,chunk,ncores,geometry
    dx=x(2)-x(1)
    xdim=size(f)
    allocate(deltax(xdim),f_next(xdim))
    f_next=f
    n=xdim-2
    chunk=ceiling(real(xdim)/ncores)
    do i=1,xdim
        deltax(i)=v(i)*dt
        if (v(i).ge.0) then
            if (i.eq.xdim) then
            !the right boundary
                f_next(i)=f_next(i)-deltax(i)/dx*f(i)
            else
                f_next(i)=f_next(i)-deltax(i)/dx*f(i)
                f_next(i+1)=f_next(i+1)+deltax(i)/dx*f(i)
            end if
        else
            if (i.eq.1) then
            !the left boundary
                f_next(i)=f_next(i)+deltax(i)/dx*f(i)
            else
                f_next(i)=f_next(i)+deltax(i)/dx*f(i)
                f_next(i-1)=f_next(i-1)-deltax(i)/dx*f(i)
            end if
        end if
    end do
    do i=1,n
        output(i)=f_next(i+1)-dt*geometry*f(i+1)*v(i)/x(i+1)
    end do
    do i=1,n+1
        flux(i)=0
        if (v(i).gt.0) then  !upwind
            flux(i)=flux(i)+(f(i)*v(i))
        end if
        if (v(i+1).lt.0) then  !upwind
            flux(i)=flux(i)+(f(i+1)*v(i+1))
        end if
    end do
    deallocate(deltax,f_next)
end subroutine advection_1d

!******************************************************************************************
!******************************************************************************************
!End ODE and PDE section
!******************************************************************************************
!******************************************************************************************

!******************************************************************************************
!******************************************************************************************
!Begin geometry section
!******************************************************************************************
!******************************************************************************************

subroutine rotate_xyz_to_yzx(a)
    !specifically for 3d euler equations, component 2-4 are xyz
    real(8) :: a(5),b(5)
    b=a
    a(1)=b(1)
    a(2)=b(3)
    a(3)=b(4)
    a(4)=b(2)
    a(5)=b(5)
end subroutine rotate_xyz_to_yzx

subroutine rotate_yzx_to_xyz(a)
    !specifically for 3d euler equations, component 2-4 are xyz
    real(8) :: a(5),b(5)
    b=a
    a(1)=b(1)
    a(2)=b(4)
    a(3)=b(2)
    a(4)=b(3)
    a(5)=b(5)
end subroutine rotate_yzx_to_xyz

subroutine ray_sphere_3d(s,e,c,r,d,code)
    !calculate the distance of intersection of a ray with a sphere
    !s=the location of the origin of the ray
    !e=the unit directional vector of the ray
    !c=the center of the sphere
    !r=radius of the sphere, d(1),d(2) are two possible intersection distance
    !code=0 no intersection, 1 one intersection, 2 two intersections
    !ray starts at the sphere does not count towards intersection
    real(8) :: s(3),e(3),c(3),r,d(2),a,b,xi
    integer :: code
    a=1d0
    b=2d0*((s(1)-c(1))*e(1)+(s(2)-c(2))*e(2)+(s(3)-c(3))*e(3))
    xi=(s(1)-c(1))**2+(s(2)-c(2))**2+(s(3)-c(3))**2-r**2
    d=0
    call quadratic_real(a,b,xi,d,code)
    if (code.eq.0) then
    else if (code.eq.1) then
        code=0   !=0 means ray instersects the sphere tangentially
    else if (code.eq.2) then      !d(1)>d(2)
        if (d(1).le.0d0) then     !=0 means ray starts at the sphere and only intersect once
            code=0
        else if (d(2).le.0d0) then   !=0 means ray starts at the sphere and intersect twice
            code=1
        else if (d(2).gt.0d0) then
            code=2
        else
            print *, 'ray_sphere_3d'
            stop
        end if
    else
        print *, 'there should be no such code, ray_sphere_3d'
        stop
    end if
end subroutine ray_sphere_3d

subroutine ray_block_single(s,e,c,r,d,blocked,dx)
    !given a ray and a sphere find if the ray is blocked by the sphere
    !if true, d return the value of distance to the sphere
    !d_temp(1)-d_temp(2) must be greater than dx or the ray is not
    !counted as blocked
    integer :: i,code
    real(8) :: s(3),e(3),c(3),r,d,d_temp(2)
    real(8), optional :: dx
    logical :: blocked
    blocked=.false.
    call ray_sphere_3d(s,e,c,r,d_temp,code)
    if (code.gt.1) then
        if (present(dx)) then
            if ((d_temp(1)-d_temp(2)).gt.dx) then
                d=d_temp(2)
                blocked=.true.
            end if
        else
            d=d_temp(2)
            blocked=.true.
        end if
    end if
end subroutine ray_block_single

subroutine lineseg_sphere_3d_group(s0,s1,c,r,d,n_seg)
    !change ray to a finite line segment compared to ray_sphere_3d_group
    !the line is from s0 to s1
    real(8) :: s0(3),s1(3),e(3),c(3),r(:),d_s0s1
    real(8), allocatable :: d(:),d_temp(:)
    integer :: n_seg,i,n,n_new,n_intersect
    d_s0s1=norm2(s1-s0)
    if (d_s0s1.eq.0d0) then
        print *,'lineseg wrong'
        stop
    end if
    e=(s1-s0)/norm2(s1-s0)
    call ray_sphere_3d_group(s0,e,c,r,d_temp,n_intersect)
    n=size(r)
    if (n_intersect/=0) then
        do i=1,n_intersect
            if (d_temp(i).gt.d_s0s1) then
                n_new=i  !line segment end here
                exit
            end if
        end do
        if (i==n_intersect+1) n_new=i  !line segment outside all the spheres
        allocate(d(n_new))
        d(1:n_new-1)=d_temp(1:n_new-1)
        d(n_new)=d_s0s1
        n_seg=n_new
    else
        n_seg=0
    end if
    if (n_seg==0) then
        print *,e, 'lineseg wrong'
        stop
    end if
    deallocate(d_temp)
end subroutine lineseg_sphere_3d_group

subroutine ray_sphere_3d_group(s,e,c,r,d,n_intersect)
    !calculate the distance of intersection of a ray with a group of concentric spheres
    !s=the location of the origin of the ray
    !e=the unit directional vector of the ray
    !c=the center of the group of spheres
    !r=radius of the group of spheres, from small to large
    !d(:) are the intersection distance from small to large, d(1).ne.0
    !n_intersect is the number of intersection
    real(8) :: s(3),e(3),c(3),r(:)
    real(8), allocatable :: d_temp(:,:),d(:)
    integer :: n_intersect,n,i,n0_intersect,n1_intersect,n2_intersect,j
    integer, allocatable :: code(:)
    n=size(r)  !how many spheres are there in the group
    n_intersect=0
    n0_intersect=0
    n1_intersect=0
    n2_intersect=0
    allocate(d_temp(2,n),code(n))
    do i=1,n
        call ray_sphere_3d(s,e,c,r(i),d_temp(:,i),code(i))
        n_intersect=n_intersect+code(i)
        if (code(i).eq.0) n0_intersect=n0_intersect+1  !# of sphere has 0 intersection
        if (code(i).eq.1) n1_intersect=n1_intersect+1  !# of sphere has 1 intersection
        if (code(i).eq.2) n2_intersect=n2_intersect+1  !# of sphere has 2 intersection
    end do
    allocate(d(n_intersect))
    j=1   !dedicated to count the index of single intersection sphere
    do i=1,n
        if (code(i).eq.1) then
            d(n2_intersect*2+j)=d_temp(1,i)
            j=j+1
        else if (code(i).eq.2) then
            d(n0_intersect+n2_intersect-i+1)=d_temp(2,i)
            d(i-n0_intersect+n2_intersect)=d_temp(1,i)
        end if
    end do
    deallocate(d_temp,code)
end subroutine ray_sphere_3d_group

subroutine ray_sphere_3d_forbid_depend(s0,e,c,r_forbid,r_depend,s1,dx)
    !forbid and dependence spheres are concentric
    !given s0,e(unit vec from s0 to s1), calculate s1
    !s1 is the farthest point on the dependence sphere of the ray
    !if the ray is blocked by the forbid sphere, s1 is the closer point
    !if the ray is neither blocked nor does intersect with the depend sphere, s1=s0
    real(8) :: s0(3),s1(3),e(3),c(3),r_depend(:),r_forbid,d_forbid
    real(8), optional :: dx   !to let the end point be farther away
    real(8), allocatable :: d(:)
    integer :: n_intersect
    logical :: blocked
    call ray_block_single(s0,e,c,r_forbid,d_forbid,blocked,dx)
    call ray_sphere_3d_group(s0,e,c,r_depend,d,n_intersect)
    if (blocked) then
        if (present(dx)) then
            s1=s0+e*(d_forbid+dx)
        else
            s1=s0+e*d_forbid
        end if
    else
        if (n_intersect==0) then
            s1=s0
        else
            if (present(dx)) then
                s1=s0+e*(d(n_intersect)+dx)
            else
                s1=s0+e*d(n_intersect)
            end if
        end if
    end if
end subroutine ray_sphere_3d_forbid_depend

subroutine polar_angle_0th_mom_integrate(I_bundle,mom0,theta)
    !integrate I over all polar angles, default is uniform dtheta
    !azimuthally symmetric, size(I_bundle) must be an even number
    real(8), dimension(:), allocatable :: I_bundle
    real(8), dimension(:), allocatable, optional :: theta
    integer :: n_polar_rays,i
    real(8) :: theta1,theta2,dtheta,mom0
    n_polar_rays=size(I_bundle)
    dtheta=2d0*pi/n_polar_rays
    mom0=0d0
    do i=0,n_polar_rays/2
        if (i.eq.0) then
            theta1=0d0
            theta2=dtheta/2d0
        else if (i.eq.n_polar_rays/2) then
            theta1=pi-dtheta/2d0
            theta2=pi
        else
            theta1=(i-0.5d0)*dtheta
            theta2=theta1+dtheta
        end if
        mom0=mom0+(cos(theta1)-cos(theta2))*2d0*pi*I_bundle(i+1)
    end do
    mom0=mom0
end subroutine polar_angle_0th_mom_integrate

subroutine solid_angle_0th_mom_integrate(I_bundle,mom0,thetadphi)
    !integrate I over all solid angles, default is uniform dtheta and dphi
    real(8), dimension(:,:), allocatable :: I_bundle
    real(8), dimension(:,:), allocatable, optional :: thetadphi
    real(8) :: theta1,theta2,dtheta,phi1,phi2,dphi,mom0
end subroutine solid_angle_0th_mom_integrate

subroutine polar_angle_1st_mom_integrate(I_bundle,mom1,theta)
    !integrate I*vec over all polar angles, default is uniform dtheta
    !azimuthally symmetric, size(I_bundle) must be an even number
    !all rays are incoming
    real(8), dimension(:), allocatable :: I_bundle
    real(8), dimension(:), allocatable, optional :: theta
    integer :: n_polar_rays,i
    real(8) :: theta1,theta2,dtheta,mom1(3),vec(3)
    n_polar_rays=size(I_bundle)
    dtheta=2d0*pi/n_polar_rays
    mom1=0d0
    do i=0,n_polar_rays/2
        if (i.eq.0) then
            theta1=0d0
            theta2=dtheta/2d0
        else if (i.eq.n_polar_rays/2) then
            theta1=pi-dtheta/2d0
            theta2=pi
        else
            theta1=(i-0.5d0)*dtheta
            theta2=theta1+dtheta
        end if
        vec=(/0d0,-cos(i*dtheta),0d0/)  !incoming rays
        mom1=mom1+(cos(theta1)-cos(theta2))*2d0*pi*I_bundle(i+1)*vec
    end do
end subroutine polar_angle_1st_mom_integrate

!******************************************************************************************
!******************************************************************************************
!End geometry section
!******************************************************************************************
!******************************************************************************************

subroutine tridiagonal_matrix(n,a,b,c,d,x)
    !see wikipedia https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    integer :: i,n
    real(8), dimension(:) :: a(n-1),b(n),c(n-1),d(n),x(n),cprime(n-1),dprime(n)
    do i=1,n-1
        if (i==1) then
            cprime(i)=c(i)/b(i)
        else
            cprime(i)=c(i)/(b(i)-a(i)*cprime(i-1))
        end if
    end do
    do i=1,n
        if (i==1) then
            dprime(i)=d(i)/b(i)
        else
            dprime(i)=(d(i)-a(i-1)*dprime(i-1))/(b(i)-a(i-1)*cprime(i-1))
        end if
    end do
    do i=n,1,-1
        if (i==n) then
            x(i)=dprime(i)
        else
            x(i)=dprime(i)-cprime(i)*x(i+1)
        end if
    end do
end subroutine tridiagonal_matrix


!******************************************************************************************
!******************************************************************************************
!Begin computer science section
!******************************************************************************************
!******************************************************************************************

subroutine sorting(a,b)
    !given an array a, sort the array from small to large and store the sorted array in b
    real(8), dimension(:) :: a,b
end subroutine sorting

function standard_deviation(a)
    real(8), dimension(:), allocatable :: a
    real(8) :: standard_deviation,s,avg
    integer :: i,n
    n=size(a)
    s=0d0
    avg=sum(a)/n
    do i=1,n
        s=s+(a(i)-avg)**2d0
    end do
    standard_deviation=sqrt(s/n)
end function standard_deviation

subroutine linear_regression(x,y,a,b)
    real(8), dimension(:), allocatable :: x,y
    real(8) :: a,b,xm,ym,s1,s2
    integer :: n,i
    n=size(x)
    xm=sum(x)/n
    ym=sum(y)/n
    s1=0d0;s2=0d0
    do i=1,n
        s1=s1+(x(i)-xm)**2d0
        s2=s2+(x(i)-xm)*y(i)
    end do
    b=s2/s1
    a=ym-b*xm
end subroutine linear_regression


!******************************************************************************************
!******************************************************************************************
!End computer science section
!******************************************************************************************
!******************************************************************************************

function logistic(L,k,x,x0)
    !logistic function has asymptotic behavior
    !L=asymptote, k=steepness, x0=midpoint
    real(8) :: logistic,L,k,x,x0
    logistic=L/(1+exp(-k*(x-x0)))
end function logistic

function normal_distribution(x,mean,sd)
    !normalized gaussian distribution
    real(8) :: normal_distribution,mean,sd,x
    normal_distribution=1/sqrt(2*pi*sd**2)*exp(-(x-mean)**2/2/sd**2)
end function normal_distribution

function gaussian(x,mu,sigma)
    !not normalized gaussian distribution
    real(8) :: gaussian,x,mu,sigma
    gaussian=exp(-(x-mu)**2/sigma**2)
end function gaussian

function sigmoid(x)
    real(8) :: sigmoid,x
    sigmoid=1d0/(1d0+exp(-x))
end function sigmoid

end module mathlib

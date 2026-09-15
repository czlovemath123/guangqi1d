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

!***********************************end polynomial root***********************

!****************************vector and tensor operator***********************

!************************end vector and tensor operator***********************



!******************************************************************************************
!******************************************************************************************
!Begin general root finding section
!******************************************************************************************
!******************************************************************************************

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

function pow(a,b)
    real(8) :: a,b,pow
    pow=a**b
end function pow

!******************************************************************************************
!******************************************************************************************
!Begin interpolation section
!******************************************************************************************
!******************************************************************************************

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

function unitsq_bilinear_interp(p,x,y)
    real(8) :: x,y,p(4),unitsq_bilinear_interp
    unitsq_bilinear_interp=(1d0-x)*(1d0-y)*p(1)+x*(1d0-y)*p(2)+(1d0-x)*y*p(3)+x*y*p(4)
end function unitsq_bilinear_interp

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

!subroutine rotate_xyz_to_yzx(a)
!    !specifically for 3d euler equations, component 2-4 are xyz
!    real(8) :: a(5),b(5)
!    b=a
!    a(1)=b(1)
!    a(2)=b(3)
!    a(3)=b(4)
!    a(4)=b(2)
!    a(5)=b(5)
!end subroutine rotate_xyz_to_yzx

!subroutine rotate_yzx_to_xyz(a)
!    !specifically for 3d euler equations, component 2-4 are xyz
!    real(8) :: a(5),b(5)
!    b=a
!    a(1)=b(1)
!    a(2)=b(4)
!    a(3)=b(2)
!    a(4)=b(3)
!    a(5)=b(5)
!end subroutine rotate_yzx_to_xyz

!******************************************************************************************
!******************************************************************************************
!End geometry section
!******************************************************************************************
!******************************************************************************************


!******************************************************************************************
!******************************************************************************************
!Begin computer science section
!******************************************************************************************
!******************************************************************************************

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

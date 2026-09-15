module hydro
use mathlib
use phylib
#if     isolver==1
    use hllc
#elif   isolver==2
    use eos_hllc_analytic
#endif
implicit none
contains

end module hydro

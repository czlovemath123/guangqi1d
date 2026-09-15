module hydroscheme
use muscl

implicit none
contains

subroutine hydro_step()
    call vanleer_hydro_unsplit()
end subroutine hydro_step

end module hydroscheme

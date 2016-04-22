
subroutine external_force(mode)
#include "control_simulation.h"
    use commons
    implicit none
    integer, intent(in) :: mode
    logical, parameter :: debug=.false.
    
    select case (mode)
    case(1) ! Elastic force  on the beads. This is an axial force
        ! taking the x axis in the position of the center of the MD box
        ! the axis goes through the points: (0,Ly/2,Lz/2) and (Lx,Ly/2,Lz/2)

        do i_part = 1, n_mon_tot
            force(2,i_part) = force(2,i_part) - k_elastic*(r0(2,i_part)-boundary(2)/2.)
            force(3,i_part) = force(3,i_part) - k_elastic*(r0(3,i_part)-boundary(3)/2.)
        end do
    end select

end subroutine external_force

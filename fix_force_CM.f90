subroutine fix_force_CM(mode)
    use commons
#include "control_simulation.h"
    implicit none
    integer, intent(in) :: mode
! This routine calculate the net force over the center of mass of the principal polymer
! chain and substract it from the position of each polymer beads.
! WARNING! First save ring force with ring_net_force(1)
    select case (mode)
	case(1)
! The force over the CM of the chain should be zero all the time
! instead, these lines can't do that because force(:,i_part) >> Fcm/N. 
! so numeric error achieve force_CM ~ 3 e-14 (!= 0)
! 
        vec_dummy(:) = 0.
        do i_part = 1 , part_init_d
            vec_dummy(:) = vec_dummy(:) + force(:,i_part)
        end do
        vec_dummy(:) = vec_dummy(:)/part_init_d ! force per bead!

        do i_part = 1 , part_init_d
            force(:,i_part) = force (:,i_part) - vec_dummy(:)
        end do
 ! debug
    !write(88,"(3g20.7)") vec_dummy(:) ! elefante force over CM chain
    fr_CM(:) = fr_CM(:) + vec_dummy(:) ! to print in the end of run
    count_CM_force = count_CM_force +1

!!!!# if RINGS == 1

    case(2)
    if( rfijo.eq.1) then
! rings are kept in X, Y & Z original position
        do i_part = part_init_d+1, n_mon_tot 
       ! y,z axes are special cases
            force(1,i_part) = 0.
            force(2,i_part) = 0.
            force(3,i_part) = 0.
        end do

     else if (rfijo .eq. 0) then
! rings are kept in Y & Z original position
        do i_part = part_init_d+1, n_mon_tot 
       ! y,z axes are special cases
            force(2,i_part) = 0.
            force(3,i_part) = 0.
        end do

!     else if (rfijo == 2) then
!        print '(/a/)', " * Rings are free. Correcting its forces too." 
!
!        vec_dummy(:) = 0.
!        do i_part = 1 , n_mon_tot 
!            vec_dummy(:) = vec_dummy(:) + force(:,i_part)
!        end do
!        vec_dummy(:) = vec_dummy(:)/n_mon_tot ! force per bead!
!
!        write(3, "(i,3f16.5)") i_time, vec_dummy(:) ! print in a file!   elefante
!!        print '(/a,3f16.5/)',"   -  Fcm (force per bead) = ", vec_dummy(:) 
!        do i_part = 1 , n_mon_tot
!              force(1,i_part) = force(1,i_part) - vec_dummy(1) 
!              force(2,i_part) = force(2,i_part) - vec_dummy(2) 
!              force(3,i_part) = force(3,i_part) - vec_dummy(3) 
!        end do
     end if
!!!!#endif
    case default
        print*, "Something is wrong with the fix_force_CM foutine. check!, itime", i_time
        stop
    end select
end subroutine fix_force_CM

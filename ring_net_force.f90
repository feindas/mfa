subroutine ring_net_force(mode)
    use commons
#include "control_simulation.h"
    implicit none
    integer, intent(in) :: mode
	logical :: file_exists

# if RINGS != 0

!print*. size(ring_force(:,dim=2)  )
! implicits variables in fortran.
!
! If the program is here, is because exist at least one ring
select case (mode)
case(-1)  ! close 
    close(26) 

case(0) ! Init file. 
!    inquire(file="ring_force.dat",exist=file_exists)
    open(unit=26, file='ring_force.dat', status='unknown')
     write(26,*) "# Fx_1 , Sx_1 , Fy_1 , Sy_1 , etc...  Fx_2 , Sx_2 , etc... "

!    if(.not.file_exists) then
!    endif
    inv_n_safe = 1./dble(n_safe)
    inv_n_obser = 1./dble(n_obser)
    
case(1) ! calculating mean forces
      ring_force(:,:) = 0.
! First Ring = 1
do i_part = 1,n_mon_d*n_chain_d
        ring_force(:,1) = ring_force(:,1) + force(:,i_part + n_mon*n_chain)
end do
     total_ring_force(:,1) = total_ring_force(:,1) + ring_force(:,1)
     total_ring_force2(:,1) = total_ring_force2(:,1) + ring_force(:,1)*ring_force(:,1)

    if(n_ring > 1) then
! Second Ring = 2
do i_part = 1,n_mon_e*n_chain_e
        ring_force(:,2) = ring_force(:,2) + force(:,i_part + n_mon*n_chain + n_mon_d*n_chain_d)
end do
     total_ring_force(:,2) = total_ring_force(:,2) + ring_force(:,2)
     total_ring_force2(:,2) = total_ring_force2(:,2) + ring_force(:,2)*ring_force(:,2)
endif

case(2) 
! saving forces with n_obser frequency ! *** WARNING *** ! be carefully between
! case(0) and case(2)

if (n_ring == 1) then
        write(26,'(6(g20.13,x))') total_ring_force(1,1)*inv_n_obser, &
sqrt( ( total_ring_force2(1,1) - total_ring_force(1,1)*total_ring_force(1,1)*inv_n_obser ) *inv_n_obser ), &
total_ring_force(2,1)*inv_n_obser, &
sqrt( ( total_ring_force2(2,1) - total_ring_force(2,1)*total_ring_force(2,1)*inv_n_obser ) *inv_n_obser ), &
total_ring_force(3,1)*inv_n_obser, &
sqrt( ( total_ring_force2(3,1) - total_ring_force(3,1)*total_ring_force(3,1)*inv_n_obser ) *inv_n_obser )
else if (n_ring == 2) then
!        write(26,'(6g16.10,6g16.10)') total_ring_force(1,1)*inv_n_obser, &
        write(26,'(12(g20.13,x))') total_ring_force(1,1)*inv_n_obser, &
sqrt( ( total_ring_force2(1,1) - total_ring_force(1,1)*total_ring_force(1,1)*inv_n_obser ) *inv_n_obser ), &
total_ring_force(2,1)*inv_n_obser, &
sqrt( ( total_ring_force2(2,1) - total_ring_force(2,1)*total_ring_force(2,1)*inv_n_obser ) *inv_n_obser ), &
total_ring_force(3,1)*inv_n_obser, &
sqrt( ( total_ring_force2(3,1) - total_ring_force(3,1)*total_ring_force(3,1)*inv_n_obser ) *inv_n_obser ), &
total_ring_force(1,2)*inv_n_obser, &
sqrt( ( total_ring_force2(1,2) - total_ring_force(1,2)*total_ring_force(1,2)*inv_n_obser ) *inv_n_obser ), &
total_ring_force(2,2)*inv_n_obser, &
sqrt( ( total_ring_force2(2,2) - total_ring_force(2,2)*total_ring_force(2,2)*inv_n_obser ) *inv_n_obser ), &
total_ring_force(3,2)*inv_n_obser, &
sqrt( ( total_ring_force2(3,2) - total_ring_force(3,2)*total_ring_force(3,2)*inv_n_obser ) *inv_n_obser )
        endif

      total_ring_force(:,:) = 0.
     total_ring_force2(:,:) = 0.

case(3)  ! SAVING FORCES (WITH n_safe FREQUENCY) ! *** WARNING *** !

if (n_ring == 1) then
        write(26,'(6(g20.13,x))') total_ring_force(1,1)*inv_n_safe, &
sqrt( ( total_ring_force2(1,1) - total_ring_force(1,1)*total_ring_force(1,1)*inv_n_safe ) *inv_n_safe), &
total_ring_force(2,1)*inv_n_safe, &
sqrt( ( total_ring_force2(2,1) - total_ring_force(2,1)*total_ring_force(2,1)*inv_n_safe ) *inv_n_safe), &
total_ring_force(3,1)*inv_n_safe, &
sqrt( ( total_ring_force2(3,1) - total_ring_force(3,1)*total_ring_force(3,1)*inv_n_safe ) *inv_n_safe)
else if (n_ring == 2) then
!        write(26,'(6g16.10,6g16.10)') total_ring_force(1,1)*inv_n_safe, &
        write(26,'(12(g20.13,x))') total_ring_force(1,1)*inv_n_safe, &
sqrt( ( total_ring_force2(1,1) - total_ring_force(1,1)*total_ring_force(1,1)*inv_n_safe ) *inv_n_safe), &
total_ring_force(2,1)*inv_n_safe, &
sqrt( ( total_ring_force2(2,1) - total_ring_force(2,1)*total_ring_force(2,1)*inv_n_safe ) *inv_n_safe), &
total_ring_force(3,1)*inv_n_safe, &
sqrt( ( total_ring_force2(3,1) - total_ring_force(3,1)*total_ring_force(3,1)*inv_n_safe ) *inv_n_safe), &
total_ring_force(1,2)*inv_n_safe, &
sqrt( ( total_ring_force2(1,2) - total_ring_force(1,2)*total_ring_force(1,2)*inv_n_safe ) *inv_n_safe), &
total_ring_force(2,2)*inv_n_safe, &
sqrt( ( total_ring_force2(2,2) - total_ring_force(2,2)*total_ring_force(2,2)*inv_n_safe ) *inv_n_safe), &
total_ring_force(3,2)*inv_n_safe, &
sqrt( ( total_ring_force2(3,2) - total_ring_force(3,2)*total_ring_force(3,2)*inv_n_safe ) *inv_n_safe)
        endif
    close(24)
   
      total_ring_force(:,:) = 0.
     total_ring_force2(:,:) = 0.

case default
     print *, "Something is wrong with the ring_net_force routine. please check!"
end select
# endif
end subroutine ring_net_force

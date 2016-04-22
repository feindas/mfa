subroutine bond_distance(mode)
#include 'control_simulation.h'
    use commons

 implicit none
 integer, intent(in) :: mode
 integer :: tot_time
 integer :: i,j
 real(kind=8)   :: tot_time_real
 include 'fftw3.f'
! include 'fftw.f'

tot_time = n_relax + n_obser
tot_time_real=tot_time

#if SYSTEM == 4

select case (mode)
case(0)
!Initializing


    allocate (bond_dist(n_mon-1), bond_dist2(n_mon-1), error_bond_dist(n_mon-1) )
    allocate (add_bond_dist(2), mean_bond_dist(2), add_bond_dist2(2), mean_bond_dist2(2))
     

! add_bond_dist(:) !cummulative value for bond distance. dim 2. 1: bonds outside the rings
															!   2: bonds between the rings
! mean_bond_dist(:) !mean value for bond distance. dim 2. 1: bonds outside the rings
														    !   2: bonds between the rings																								

! add_bond_dist2(:) !cuadratic cummulative value for bond distance. dim 2. 1: bonds outside the rings
																		!  2: bonds between the rings
! mean_bond_dist2(:) !cuadratic mean value for bond distance. dim 2. 1: bonds outside the rings
																		!  2: bonds between the rings

do j=1,n_mon-1
     bond_dist(j)=0.
     bond_dist2(j)=0.
     error_bond_dist(j)=0.
end do

do i=1,2
     add_bond_dist(i)=0.
     mean_bond_dist(i)=0.
     add_bond_dist2(i)=0.
     mean_bond_dist2(i)=0.
end do



case(2)
!Observation


do j=1,n_mon-1
     bond_dist(j)=bond_dist(j)+sqrt( (r0(1,j+1) - r0(1,j))*(r0(1,j+1) - r0(1,j)) + (r0(2,j+1) - r0(2,j))*(r0(2,j+1) - r0(2,j)) + (r0(3,j+1) - r0(3,j))*(r0(3,j+1) - r0(3,j)) )
     bond_dist2(j)=bond_dist2(j)+( (r0(1,j+1) - r0(1,j))*(r0(1,j+1) - r0(1,j)) + (r0(2,j+1) - r0(2,j))*(r0(2,j+1) - r0(2,j)) + (r0(3,j+1) - r0(3,j))*(r0(3,j+1) - r0(3,j)) ) 
end do


case(1)
!Storing

    open(unit=35, file='bond_distances.dat',status='unknown')
!    open(unit=36, file='error_bond_distances.dat',status='unknown')

    do j=1,n_mon-1
!    error_bond_dist(j)=sqrt( (bond_dist2(j)/tot_time) - (bond_dist(j)/tot_time)*(bond_dist(j)/tot_time) ) / sqrt(tot_time_real)
    write(35,*) j, bond_dist(j)/(tot_time-n_relax), bond_dist2(j)/(tot_time-n_relax)
!    write(36,*) j, error_bond_dist(j)    

    end do
end select
#endif /*SYSTEM 4 */
end subroutine bond_distance 

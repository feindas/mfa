subroutine gen_chain(chain_type)

use ziggurat !, only: rnor, uni
use commons
use util

implicit none

integer ::  n,i,j, ierror, seed=147892376
integer, intent(in) :: chain_type
!integer, parameter :: n = 60000 !number of beads
!real (kind = 8), parameter :: sigma = 1
real (kind = 8) :: r(1:3),modran,D, l=1.
real (kind = 8) :: PIii,amp,peri,peri2, amp2
!------------------------------------------
!real (kind = 8), allocatable :: M(:,:)
! In the mfa_prog M(3,N) is r0(3,N)
!------------------------------------------

!   ------   Entry parameters stage   ------   

!open(unit=11, file = 'pchain', status = 'old',iostat=ierror )
!if(ierror.eq.0) then
!        read(11,*)
!        read (11, *) n, l, D
!        print *, "# bolitas: ", n, "bolitas distance(radio en&
!plano yz): ", l,"end-to-end distance: ", D
!else 
!        print *, "Error de lectura de valores de entrada."
!end if

!close(unit=11)

	write(*,*) "Executing  gen_chain ... "


    n = n_mon*n_chain
    D = boundary(1)

    select case(chain_type)
    case(0) !chain with regular positions in x axis, random in y and z axis
! Allocations - It's not necessary now
!        allocate (M(3,n))

! first particle is clampted in x=0, y=0, z=0
! the position of each particle is stored in the M's colum
!       M(:,1)=0.
!        r0(1,1)=0.
        r0(1,1)=(boundary(1)/n)/2.0
        r0(2,1)=boundary(2)/2.0 !+ (2.0*uni()-1)*l
        r0(3,1)=boundary(3)/2.0 !+ (2.0*uni()-1)*l
! Initialized in init_param.f90
!call zigset(seed)
!
        do j=1,n_mon*n_chain
            r(2) = (2.0*uni()-1) !the first random run from -1 to 1
            r(3) = (2.0*uni()-1) !the second random run from -1 to 1
            modran = sqrt( r(2)*r(2)+r(3)*r(3) )
            r(2)=r(2)/modran
            r(3)=r(3)/modran
! equispaced x position
!        M(1,j)= (D/N)*(j-1)
!r0(1,j)= boundary(1) + (D/N)*(j-1) 
            r0(1,j)= (boundary(1)/n)*( j - 1/2.0 )
! y & z position in circle of radius l
!        M(2,j)= r(2)*l
!        M(3,j)= r(3)*l
            r0(2,j)= boundary(2)/2.0 + r(2)*l
            r0(3,j)= boundary(3)/2.0 + r(3)*l

            ! ------ the new chain is stored in file: new_chain
            !open(unit=10, file = 'new_chain', status = 'new', action = 'write', &
            !        position = 'append', iostat = ierror)
            !if (ierror.eq.0) then
            !        write(10,* ) n_mon
            !        write(10,* ) 
            !        do i_part=1,n_mon
            !!                write( 10, * ) " O ",  M(1,i), M(2,i), M(3,i)
            !                write( 10, * ) " O ",  r0(1,i_part), r0(2,i_part), r0(3,i_part)
            !        end do
            !end if
            !close(unit=10)
        end do

    case(1) !chain with harmonic positions in x axis, random in y and z axis
	
       print '(/a/)'," Generating harmonic chain "
	
! Allocations - It's not necessary now
!        allocate (M(3,n))

! first particle is clampted in x=0, y=0, z=0
! the position of each particle is stored in the M's colum
!       M(:,1)=0.
!        r0(1,1)=0.
        r0(1,1)=(boundary(1)/dble(n))/2.0
        r0(2,1)=boundary(2)/2.0 !+ (2.0*uni()-1)*l
        r0(3,1)=boundary(3)/2.0 !+ (2.0*uni()-1)*l
! Initialized in init_param.f90
!call zigset(seed)
!
	
        do j=1,n
            r(2) = (2.0*uni()-1) !the first random run from -1 to 1
            r(3) = (2.0*uni()-1) !the second random run from -1 to 1
            modran = sqrt( r(2)*r(2)+r(3)*r(3) )
! ? Romina: explain             
            r(2)=0.*r(2)/modran
            r(3)=0.*r(3)/modran
! harmonic x position
            PIii=4.D0*DATAN(1.D0)

            peri = 512.
            peri2 = 170.
            amp=1.8
            amp2=1.5
            !    print '(/a,f12.4/)', "Using  period: ",peri 

            r0(1,j)= (boundary(1)/dble(n))*( j - 1/2.0 )
            r0(2,j)= (boundary(2)/2.0) + r(2)*l
            write (204,'(3f16.8)') boundary(2), r(2)
            r0(3,j)= boundary(3)/2. + amp*sin(j*2*PIii/peri)
            r0(3,j)=r0(3,j) + amp2*sin(j*2*PIii/peri2)
            !       r0(3,j)=r0(3,j) + amp2*sin(j*2*PIii/peri2 + PIii/(peri2*6.))


            !print *,amp*sin(j*2*PIii/peri),amp2*sin(j*2*PIii/peri2)
            ! y & z position in circle of radius l
            !        M(2,j)= r(2)*l
            !        M(3,j)= r(3)*l


            ! ------ the new chain is stored in file: new_chain
            !open(unit=10, file = 'new_chain', status = 'new', action = 'write', &
            !        position = 'append', iostat = ierror)
            !if (ierror.eq.0) then
            !        write(10,* ) n_mon
            !        write(10,* ) 
            !        do i_part=1,n_mon
            !!                write( 10, * ) " O ",  M(1,i), M(2,i), M(3,i)
            !                write( 10, * ) " O ",  r0(1,i_part), r0(2,i_part), r0(3,i_part)
            !        end do
            !end if
            !close(unit=10)


            !          write (101,'(a,3f16.8)') "Cl", r0(:,j)
            !          write (201,'(3f16.8)') r0(1,j), r0(3,j)
   end do
    case(2) ! Chain with regular positions in x axis, random in y and z axis  plus  fix-ends in the faces of the MD-box  
        
! First bead:  fixed in the left wall of the box
        r0(1,1)=0.2*sigma(1,1)               ! (boundary(1)/n)/2.0
        r0(2,1)=boundary(2)/2.0 !+ (2.0*uni()-1)*l
        r0(3,1)=boundary(3)/2.0 !+ (2.0*uni()-1)*l

        do j=2,n_mon*n_chain-1
            r(2) = (2.0*uni()-1) !the first random run from -1 to 1
            r(3) = (2.0*uni()-1) !the second random run from -1 to 1
            modran = sqrt( r(2)*r(2)+r(3)*r(3) )
            r(2)=r(2)/modran
            r(3)=r(3)/modran
            r0(1,j)= (boundary(1)/n)*( j - 1/2.0 )
            r0(2,j)= boundary(2)/2.0 + r(2)*l
            r0(3,j)= boundary(3)/2.0 + r(3)*l
        end do
!   Last bead: fixed in the right wall of the box        
        r0(1,n_mon*n_chain)=boundary(1)-0.2*sigma(1,1)
        r0(2,n_mon*n_chain)=boundary(2)/2.0 
        r0(3,n_mon*n_chain)=boundary(3)/2.0 
#ifdef      RUN_2D
           r0(2,:)= boundary(2)/2.0 
#endif

end select
end subroutine gen_chain

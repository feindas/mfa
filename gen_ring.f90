subroutine gen_ring(mode)
    use commons
    use functions
!    use ziggurat
    implicit none
 
    integer, intent(in) :: mode !number of rings to be generated
    integer :: j, k
    real (kind = 8) :: radio=1.5, dist_inter_chain = 0.98  
    
    radio = ring_radio
    
	print*, "WARNING! *****************************************************"
	print*, "WARNING! * THESE RINGS HAVE radio = ", radio, " * "
	print*, "WARNING! *****************************************************"

	

select case (mode)
case(0) 
    ! no rings
    print*, " These runs don't have rings"

case(1)
!The ring will be centered in the middle of the chain
    do k= 1, n_chain_d
        do j=1, n_mon_d
        r0(1,(j+n_mon*n_chain)+(k-1)*n_mon_d) = boundary(1)/2.0 + (dist_inter_chain/2) * (2*k-n_chain_d-1)
        r0(2,(j+n_mon*n_chain)+(k-1)*n_mon_d) = boundary(2)/2.0 + radio*cos( (2*(j-1)+k-1)*pi/dble(n_mon_d) )
        r0(3,(j+n_mon*n_chain)+(k-1)*n_mon_d) = boundary(3)/2.0 + radio*sin( (2*(j-1)+k-1)*pi/dble(n_mon_d) )

!      r0(1,j) = boundary(1)/2.0 ! in the middle of the chain in X direction
!      r0(2,j) = boundary(2)/2.0 - radio*cos(6.28318*dble(j)/dble(n_mon_d)) !equispaced angle position
!      r0(3,j) = boundary(3)/2.0 - radio*sin(6.28318*dble(j)/dble(n_mon_d))
!      ! print*, 'facundo_ring_X', r0(1,j), ' facundo_ring_Y', r0(2,j), ' facundo_ring_Z', r0(3,j)
        end do
    end do
! calculate the center of mass of the ring
!do j=n_mon*n_chain+1, n_mon*n_chain + n_mon_d
!        r0(1,j) = boundary(1)/2 !r0(1,p)   !in x random position, nothig to do
!        r0(2,j) = boundary(2)/2 + radio*cos(6.28318*dble(j)/dble(n_mon_d)) !equispaced angle position
!        r0(3,j) = boundary(3)/2 + radio*sin(6.28318*dble(j)/dble(n_mon_d))
!        print*, 'facundo_ring_X', r0(1,j), ' facundo_ring_Y', r0(2,j), ' facundo_ring_Z', r0(3,j)

!        CM_ring(1,1) = CM_ring(1,1) + r0(1,j) ! x of the first ring
!        CM_ring(2,1) = CM_ring(2,1) + r0(2,j) ! y of the first ring
!        CM_ring(3,1) = CM_ring(3,1) + r0(3,j) ! z of the first ring
!enddo

!        CM_ring(:,1) = CM_ring(:,1)/n_mon_d

! ------ the new ring is stored in file: new_ring
!open(unit=10, file = 'new_ring.xyz', status = 'new', action = 'write', &
!        position = 'append')
!
!        write(10,* ) n_mon*n_chain+n_mon_d*n_chain_d
!        write(10,* ) 
!        do i_part=1,n_mon*n_chain
!                write( 10, * ) " N ",  r0(1,i_part), r0(2,i_part), r0(3,i_part)
!        end do
!        write(10,*) n_mon_d
!        write(10,*) 
!        do i_part=n_mon*n_chain+1,n_mon*n_chain + n_mon_d
!                write( 10, * ) " O ",  r0(1,i_part), r0(2,i_part), r0(3,i_part)
!        end do
!close(unit=10)
!
!stop

case(2) ! two rings
! Rings are centered in the middle of the chain minus(plus) ring_distance/2

 do k= 1, n_chain_d
    do j=n_mon*n_chain+1, n_mon*n_chain + n_mon_d
    r0(1,j+(k-1)*n_mon_d) = (boundary(1)-ring_dist + dist_inter_chain * (2*k-n_chain_d-1)  )/2.0
    r0(2,j+(k-1)*n_mon_d) = boundary(2)/2.0 + radio*cos( dble(2*(j-1) + k-1)*pi/dble(n_mon_d) )
    r0(3,j+(k-1)*n_mon_d) = boundary(3)/2.0 + radio*sin( dble(2*(j-1) + k-1)*pi/dble(n_mon_d) )
   
!    do j=n_mon*n_chain+1, n_mon*n_chain + n_mon_d
!        r0(1,j) = (boundary(1)-ring_dist)/2.0 !r0(1,p)   !in x random position, nothig to do
!        r0(2,j) = boundary(2)/2.0 - radio*cos(6.28318*dble(j)/dble(n_mon_d)) !equispaced angle position
!        r0(3,j) = boundary(3)/2.0 - radio*sin(6.28318*dble(j)/dble(n_mon_d))
    
        end do
    end do

!        do j=(n_mon*n_chain+ n_mon_d) + 1, (n_mon*n_chain + n_mon_d) + n_mon_e
!        r0(1,j) = (boundary(1)+ring_dist)/2.0 !r0(1,p)   !in x random position, nothig to do
!        r0(2,j) = boundary(2)/2.0 - radio*cos(6.28318*dble(j)/dble(n_mon_e)) !equispaced angle position
!        r0(3,j) = boundary(3)/2.0 - radio*sin(6.28318*dble(j)/dble(n_mon_e))
!        end do
    do k=1, n_chain_e
        do j=(n_mon*n_chain+ n_mon_d*n_chain_d) + 1, (n_mon*n_chain + n_mon_d*n_chain_d) + n_mon_e
    r0(1,j+(k-1)*n_mon_e) = (boundary(1)+ring_dist + dist_inter_chain * (2*k-n_chain_e-1)  )/2.0
    r0(2,j+(k-1)*n_mon_e) = boundary(2)/2.0 + radio*cos( dble(2*(j-1) + k-1)*pi/dble(n_mon_e) )
    r0(3,j+(k-1)*n_mon_e) = boundary(3)/2.0 + radio*sin( dble(2*(j-1) + k-1)*pi/dble(n_mon_e) )
!        r0(1,j+(i-1)*n_mon_e) = (boundary(1)+ring_dist + dist_inter_chain * (2*i-n_chain_e-1)  )/2.0
!        r0(2,j+(i-1)*n_mon_e) = boundary(2)/2.0 + radio*cos(  (2*j + i-1)*pi/dble(n_mon_e) )
!        r0(3,j+(i-1)*n_mon_e) = boundary(3)/2.0 - radio*sin(  (2*j + i-1)*pi/dble(n_mon_e) )
        end do
    end do
    
    case default
            print *, "  Sorry! check the variable n_ring = ", n_ring
            stop
end select

end subroutine gen_ring

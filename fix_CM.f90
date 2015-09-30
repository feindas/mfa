subroutine fix_CM()
    use commons
    implicit none
    real (kind = 8) :: x=0,y=0,z=0
! *******************************************************************************
! This routine calculate the position of center of mass of the principal polymer
! chain set it in the middle of the simulation box. 
! Then put the center of mass of the rings in the same place. <- NOT so sure
! *******************************************************************************

! Performing a rigid shift over r0_unfold. Setting the center of mass of the chain (NOT PHYSICS)
! Only for a pretty view 
    vec_dummy(:) = 0.
    r_dummy = 0.
    do i_part = 1 , n_mon
      vec_dummy(:) = vec_dummy(:) + mass(i_part)*r0_unfold(:,i_part)
      r_dummy = r_dummy + mass(i_part)
    end do
      vec_dummy(:) = vec_dummy(:)/r_dummy
!      print '(a,3f16.5)', "  S4 First - Rcm_UNFOLD (chain): ", vec_dummy(:)   

! Para mirar la pelicula es importante r0_unfold, r0 no debe cambiar (representa la dinamica)
    do i_part = 1 , n_mon_tot
      r0_unfold(:,i_part) = r0_unfold(:,i_part) - vec_dummy(:) + boundary(:)/2.
!r0 no deberia cambiar      r0(:,i_part) = r0(:,i_part) - vec_dummy(:) + boundary(:)/2. 
    end do
! Este shift rigido no deberia modificar la posicion de los anillos cuando anillos fijos. elefante

end subroutine fix_CM


subroutine dpd_forces(inv_sqrt_r_2)
    use commons
    use ziggurat, only: rnor,uni
#include "control_simulation.h"
    implicit none
    real(kind=8) , intent(in) :: inv_sqrt_r_2
    real(kind=8) :: g_rand,r_versor(3),delta_v(3),rrc
!    integer , intent(in)      :: i_type,j_type
!
!   -------- DPD forces calculation : Version inside the force caculation  
!
! WARN: check-uout if works with different masses as is. 

!    *** Gaussian random distribution. sig^2 = 1 , <x> = 0

!      g_rand(:) = (/ rnor(),rnor(),rnor() /)    

          g_rand =  rnor()

!    *** Uniform distribution, claimed to be as good as the gaussian
! (fullfills 1st and second momentum properties of the random variable)

!       g_rand = sqrt(3.)*(2.*uni() - 1.) 

!       print *,"random",g_rand(:)
!      if(debug)      print*,"random",g_rand 

! --- Weight functions 
#if DPD_WEIGHT == 0
#       ifndef DPD_CUT_OFF
             w_d = (1. - sqrt(r_2*inv_range_2(i_type,j_type)))  
#       else      
             w_d = (1. - sqrt(r_2/r_cut_dpd_2) )
#       endif      
          w_r = w_d
          w_d = w_d*w_d
#endif      

!!note:  #if DPD_WEIGHT == 1 as the functions are constants, they are defined in init_params.f90
!!               and w_d=w_r for ever in the rest of the program
!!note:        w_d = 1.
!!note:  #endif      

#if DPD_WEIGHT == 2
      w_d = sqrt( 1. - sqrt(r_2/range_2(i_type,j_type)) )  
      w_r = sqrt(w_d)
#endif      
#if DPD_WEIGHT == 3
      rrc= sqrt(r_2/range_2(i_type,j_type))
      w_d = ( 1. - rrc )**0.25  
      w_r = sqrt(w_d)
#endif      
#if DPD_WEIGHT == 4
      rrc= sqrt(r_2/range_2(i_type,j_type))
      w_d = ( ( 1. - rrc )**0.25 ) *100.* exp (-(1-rrc)**2)  
      w_r = sqrt(w_d)
#endif      

! ----  Random force computation:      

          r_versor(:) = delta_r(:)*inv_sqrt_r_2
         vec_dummy(:) = sig*w_r*g_rand*r_versor(:) 

!#       if BIN_TYPE == 1             
!             vec_dummy(:) = 0.5*vec_dummy(:)
!#       endif
#   if BIN_TYPE == 0 
      force(:,i_part) =  force(:,i_part) + vec_dummy(:)   
      force(:,j_part) =  force(:,j_part) - vec_dummy(:)   
#   elif BIN_TYPE == 1
      if (i_part>j_part) then 
      force(:,i_part) =  force(:,i_part) + vec_dummy(:)   
      force(:,j_part) =  force(:,j_part) - vec_dummy(:)   
      endif
#   endif
     
! Dissipative force computation      

      delta_v(:) = v(:,i_part) - v(:,j_part)

      r_dummy = delta_v(1)*r_versor(1) + delta_v(2)*r_versor(2) + delta_v(3)*r_versor(3)

      vec_dummy(:) = -1.*friction(1)*w_d*r_dummy*r_versor(:)
        
!#       if BIN_TYPE == 1             
!             vec_dummy(:) = 0.5*vec_dummy(:)
!#       endif


      force(:,i_part) =  force(:,i_part) + vec_dummy(:)   
#     if BIN_TYPE == 0      
      force(:,j_part) =  force(:,j_part) - vec_dummy(:)   
#     endif      



end subroutine dpd_forces

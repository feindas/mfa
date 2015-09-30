
!! lgv_forces !! 
! * Computes random and friction forces for langevin thermostat 
! * It is assumed that the routine is called inside the fluid_fluid routine

subroutine lgv_forces()
    use commons
    use ziggurat, only: rnor,uni
#include "control_simulation.h"
    implicit none
!    real(kind=8) , intent(in) :: inv_sqrt_r_2
    real(kind=8) :: g_rand,r_versor(3),delta_v(3),rrc


!    *** Uniform distribution, claimed to be as good as the gaussian
! (fullfills 1st and second momentum properties of the random variable)
! See: Paul & Duenweg 

!       g_rand = sqrt(3.)*(2.*uni() - 1.) 


! ----  Random force computation:      


      vec_dummy(1) = sig*rnor()*sqrt(mass(i_part))
      vec_dummy(2) = sig*rnor()*sqrt(mass(i_part))
      vec_dummy(3) = sig*rnor()*sqrt(mass(i_part))
        

!      force(:,i_part) =  force(:,i_part) + vec_dummy(:)*sqrt(mass(i_part))   

     
! ---- Dissipative force computation      

      vec_dummy(:) = vec_dummy(:) -friction(1)*mass(i_part)*v(:,i_part)

#   ifdef LGVX_0
            vec_dummy(1) = 0. ! No thermostat in X
#   endif        

      force(:,i_part) =  force(:,i_part) + vec_dummy(:)   


end subroutine lgv_forces 

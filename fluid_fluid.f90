      subroutine fluid_fluid()
#include "control_simulation.h"
      use commons 
!!!#ifdef DPD_EMBEDDED
!!!      use ziggurat, only: rnor,uni
!!!#endif
     implicit none
     real (kind=8) :: delta_v(3),r_versor(3),g_rand,rrc,f_ipart(3) ! Needed for DPD_EMBEDDED only
     real (kind=8) :: l_eps,r_cut2,r_61
     real (kind=8) :: f_cou_real(3)
     real(kind=8)  :: inv_r_2,inv_sqrt_r_2
! cache blocking 
     integer :: ii_part,ii_neigh
      logical, parameter :: debug=.false.
      
#       if SYSTEM == 2 || SYSTEM == 3
            v_coul = 0.
#       endif

!!#ifdef DPD_EMBEDDED
! NOTE: this is at the beggining of dpd_forces_ll when the DPD forces are not
! embededd here. This must be added if dpd_forces is never called in the program
! flow. 

! Force switch-on update  

!---  Define temperature ramp

!!!obs #       if THERMOSTAT == 0 || THERMOSTAT == 1  /* if DPD is used we do not use T ramp */
!!!obs 
!!!obs          if((r_time.lt.temp_time).and.(temp_init.ne.temp_final)) then
!!!obs           temp = temp_init + r_time*(temp_final-temp_init)/temp_time
!!!obs          else
!!!obs           temp = temp_final
!!!obs          end if
!!!obs 
!!!obs #       endif      
! 
!---  Define effective minimium distance used for force "switch on" at the beginning 
!     This is done for each time step

#ifdef FORCE_SWITCH_ON

      if(r_time.lt.r_2_min_time) then
       r_2_min = (1.-r_time/r_2_min_time)*r_2_min_init !*range_2(:,:)!ORI*r_2_min_init
      else
       r_2_min = 0.
      end if
#endif      


      v_fluid_fluid = 0.
!
! LJ fluid-fluid force and V  calculation 
!

     !ori do i_part = 1,n_mon_tot  !n_mon_tot= brushes + droplet/melt
      do i_part = 1,n_part  !n_mon_tot= brushes + droplet/melt


          i_type = a_type(i_part)
          i_dummy = ff_list(0,i_part)

        !  print *, ff_list(:,1) ; stop
          f_ipart(:) =0.

          do i_neigh = 1, i_dummy


              j_part = ff_list(i_neigh,i_part)
              j_type = a_type(j_part)

              r_cut2 = range_2(i_type,j_type)

              ! HPC

              delta_r(1) = r0(1,i_part) - r0(1,j_part)
              delta_r(2) = r0(2,i_part) - r0(2,j_part)
              delta_r(3) = r0(3,i_part) - r0(3,j_part)

              !----- PBC ----
              ! HPC

              delta_r(1) = delta_r(1) - boundary(1)*int(2.*delta_r(1)*inv_boundary(1))
              delta_r(2) = delta_r(2) - boundary(2)*int(2.*delta_r(2)*inv_boundary(2))
#       if SYMMETRY == 1
              delta_r(3) = delta_r(3) - boundary(3)*int(2.*delta_r(3)*inv_boundary(3))
#       endif

              r_2 =  delta_r(1)*delta_r(1)  + delta_r(2)*delta_r(2) +  delta_r(3)*delta_r(3)

              !-----  check whether interaction takes place

              if( r_2 .lt. r_cut2 ) then
#           ifdef FORCE_SWITCH_ON 
                  r_2 = max(r_2,r_2_min) 
#           endif                 

                  inv_r_2 = 1./r_2
                  inv_sqrt_r_2 = sqrt(1./r_2)
                  l_eps = epsil(i_type,j_type)


                  r_61= sigma_2(i_type,j_type)*inv_r_2 
                  r_6 = r_61*r_61*r_61

                  r_12 = r_6*r_6
                  pot_loc = (r_12-r_6) - e_shift(i_type,j_type)

                  v_fluid_fluid = v_fluid_fluid + l_eps*pot_loc

                  r_dummy = l_eps*(-12*r_12+6*r_6)*inv_r_2

                  force_loc(1) = r_dummy*delta_r(1)
                  force_loc(2) = r_dummy*delta_r(2)
                  force_loc(3) = r_dummy*delta_r(3)
                  ! HPC

#       if BIN_TYPE == 0 /* binning.f90 */
                  f_ipart(1) = f_ipart(1) -  force_loc(1)
                  f_ipart(2) = f_ipart(2) -  force_loc(2)
                  f_ipart(3) = f_ipart(3) -  force_loc(3)
#       endif 
                  force(1,j_part) = force(1,j_part) + force_loc(1)
                  force(2,j_part) = force(2,j_part) + force_loc(2)
                  force(3,j_part) = force(3,j_part) + force_loc(3)

                  !  --- DPD calculation if DPD has not its own cut-off  
#if THERMOSTAT == 0
                  !!#       ifdef DPD_EMBEDDED
#     ifndef DPD_CUT_OFF
#           if SYMMETRY == 0
#               if SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 3 /* if there are brushes with fix heads*/
                  !        Exclude heads from DPD calculation 
#                  ifndef FREE_HEADS
#                  endif

#           if WALL == 1 /* explicit wall*/
                  if(i_type.eq.n_type) cycle       ! Excluding thermostat for wall-fluid
                  if(j_type.eq.n_type) cycle       ! Excluding heads thermostat for wall-fluid

#           endif                    


#      endif

#                  if PINNED == 1
                  if(i_type.eq.3) cycle       !  Excluding  fixed particles
                  if(i_type.eq.4) cycle       ! 
                  if(j_type.eq.3) cycle       ! 
                  if(j_type.eq.4) cycle       ! 
#                  endif 
#               endif 

                  call dpd_forces(inv_sqrt_r_2)

#       endif /* not DPD_CUT_OFF */
#endif /*THERMOSTAT = 0 */


#       if SYSTEM == 2  || SYSTEM  == 3 /* charged systems */
                  !!! CLAUDIO NARAMBUENA

                  ! Coulomb potential and interaction in REAL space 

                  call ewald_real()

#       endif

              end if ! if the particle is inside the interaction sphere 

!
! DPD has its own cutoff radius 
!

#ifdef DPD_CUT_OFF
#       if THERMOSTAT == 0
              !!#              ifdef DPD_EMBEDDED
#                      if SYMMETRY == 0
#                              if SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 3
              !                               Exclude heads from DPD calculation 
#                                   ifndef FREE_HEADS
              if(i_type.eq.1) cycle       ! Excluding heads !!!
              if(j_type.eq.1) cycle       ! Excluding heads !!!
#                                   endif
#                                   if WALL == 1 /* explicit wall*/
              if(i_type.eq.n_type) cycle       ! Excluding thermostat for wall-fluid
              if(j_type.eq.n_type) cycle       ! Excluding heads thermostat for wall-fluid
#                                   endif                    
#                                  if PINNED == 1

              if(i_type.eq.3) cycle       !  Excluding  fixed particles
              if(i_type.eq.4) cycle       ! 
              if(j_type.eq.3) cycle       ! 
              if(j_type.eq.4) cycle       ! 
#                                  endif 

#                              endif
#                       endif 

              if(r_2 < r_cut_dpd_2 ) then 
                  call dpd_forces(inv_sqrt_r_2)
              end if

              !!#              endif DPD_EMBEDDED
#       endif /*THERMOSTAT == 0 */

#endif /* DPD_CUT_OFF */


          end do ! loop over particles neighbors (j_part)

#       if BIN_TYPE == 0 /* binning.f90 */
          force(1,i_part) = force(1,i_part)  + f_ipart(1)
          force(2,i_part) = force(2,i_part)  + f_ipart(2)
          force(3,i_part) = force(3,i_part)  + f_ipart(3)
#       endif

#   if THERMOSTAT == 1 /* Langevin */
          call lgv_forces() 
#   endif

      end do ! loop over particles

#ifdef FLUKT 
    if((n_mon_d.or.n_chain_d.or.n_mon_e.or.n_chain_e).gt.2) then
    write(*,*) "the file fort.555 gonna be too big! undefine flag FLUKT or decrease number of particles4 or -3"
    stop
    endif
    if (mod(i_time, 100).eq.0) then
    write(555,'(1i,18f15.4)') i_time,  r0(part_init_d+1:n_mon_tot, :), & 
                              v(:,part_init_d+1:n_mon_tot), force(part_init_d+1:n_mon_tot, :)
    end if
#endif
#                 if BIN_TYPE == 1                  
                  v_fluid_fluid = 0.5*v_fluid_fluid
#                 endif
#if SYSTEM == 2 || SYSTEM == 3
        v_coul = v_coul + v_coul_self
#endif

! Debugging
if(debug) then
!print '(3f15.5)', force(:,1:n_chain*n_mon:n_mon)
          print *,"v_fluid_fluid=",i_time,v_fluid_fluid/dble(n_mon_tot)
r_dummy =     sum( sum(force(:,:)**2,dim=1) ,dim =1 )
print *,"[fluid_fluid]quad_m_force", i_time,sqrt(r_dummy)/dble(n_mon_tot)
print *, "[fluid_fluid]V=",v_fluid_fluid
!do i_part=1,n_mon_tot
!write(88,'(3i,3f15.4)') i_time, i_part,a_type(i_part),force(:,i_part)
!end do
end if
          !
end subroutine fluid_fluid

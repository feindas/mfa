      program md_pb 
#include "control_simulation.h"
      use commons 
      use util  
      
      implicit none
      real (kind=8) :: t0,t1 !,t2,t3,iter,t_inst ! times to profile the program 
      integer :: tot_time
      logical,parameter :: debug =.false.

! [4/2010] Heavy changes in data structures and improvement of some routines 
! [5/2009]: Additions by Leonid Spirin: pinning, stars, LJ wall and shear
! protocols 
! [10/07/08]: Extended DPD cutoff, for better Temperature control.WARN: It changes
! a lot the friction 
! [20/11/07] Coulomb interactions added. It performs Ewald sum with dipole
! corrections to take into account the walls in Z. SYSTEM = 2 added.
! [18/11/07] Now with SYMMETRY = 1 the program has PBC in three directions for
! BULK simulations 
! [10/10/07] Particle 4 interaction added to the program. It can be also a  chain.
! [09/05/07] Fix an old error when writing confs. f_film discarded
! [2006] Incorporation of C-preprocessor directives
! [18/6/04] Pure velocity verlet  for integration and DPD thermostat implementation of
!           Random and Dissipative forces
! [10/6/04] Shear: constant velocity for the heads of brushes. Force on heads re written.    
! [6/6/04]  Added a 9-3 flat wall potential (instead of explicit wall atoms) 
! [5/5/04]  Dynamic allocation of forces and position. New input file: system_input
! [14/4/04] claudio: using now ifort as the main compiler
! for development
! [then] Torsten Kreer : Brushes and melt with shear
! [originally from] Martin Muesser

        call get_walltime(t0)
        print '(/a/)'," *** MFA_PROG *** MD, DPD, and more. [Compilation Date: "//__DATE__//"]"


! ---- Write out current program compilation settings ----------     
!
           call messages()

! ---- Setup some flags

!      if(ad_flag.eq.0) then
!       f_angle=0
!      end if

! --------------- ENDS writing  out of compilation settings ---------------------------     
      
! ---- Inizialization  routines

           call init_system() ! physical system. Reads system_input
           call init_params()   ! md simulation parameters, box dimensions. Read mfa_input
           call init_config()   ! initial conf or read an old one
           call init_obser()    ! Measurable variables init.

!!!!Initialize Bending Routines
#ifdef BENDING     
          call bending(0) ! writes out chain bending constants parameters
#endif
#ifdef ORIENTATION     
          call orientation(0) ! writes out chain bending orientation parameters
#endif
!!!! End Initialize Bending Routines
  

           tot_time = n_relax + n_obser 

       do i_time = 1 , tot_time !  MAIN TIME LOOP 
           r_time = dble(s_time+i_time-1)*dt

! ---- > > > >  Propagate coordinates < < < < ----

           call verlet_positions()

!deb           if (i_time == 1 ) then
!deb               call write_conf(1,r0,10)
!deb               stop
!deb           end if

           call check_skin  !calculates if it is necessary to update the verlet list. If it is,
#          if BIN_TYPE == 0 
                if (f_skin.eq.1) call binning
#           elif BIN_TYPE == 1 
                if (f_skin.eq.1) call my_binning
#           endif

! ----- Forces to zero
           force(:,:) = 0.0  


! NOTE: fluid-fluid calculates LJ always, DPD and LGV forces and calls ewald in
! real space
           call fluid_fluid() ! Calculates various forces, including thermostat and LJ

#   if SYMMETRY == 0
#       if WALL != 1 
           call fluid_wall(wall_flag) ! 1= wall atoms, 2= 9-3 potenti , 3 and 4 also valid
#       endif

           !try without                   call wall_wall(wall_flag)  ! 1= wall atoms, 2= 9-3 potential
#       if WALL == 1
           !try without                   call intra_wall
#       endif
#   endif

           call intra_molec

!Add Bending Forces and energy
#ifdef BENDING 

          call bending(2)  ! adds chain bending forces and bending energy.
#   if SYMMETRY == 1
                  
          call orientation(3) ! 3 Does nothing
                              ! 2 is with orientation bending 
                              !adds chain bending forces and bending energy for
                              !first and last monomers, with PBC in 3D
#   else if SYMMETRY == 0
          call orientation(1) ! adds chain bending forces and bending energy for
                              !first and last monomers, with PBC in 2D
#   endif
#       if EXT_FORCE == 1 
            call external_force(1)
#       endif
#endif


#   if SYSTEM == 2 || SYSTEM == 3
           call ewald_k(1)  ! coulomb force calculation in K-space for Ewald sum
#         if SYMMETRY == 0
           call dipolar_correction()
#         endif
#   else if SYSTEM == 4
! Claudio, Sept. 15            
#   if RINGS == 1
           call ring_net_force(1) ! Acumulate force over rings

           call fix_force_CM(2) ! set zero force over rings
#          endif
#   endif
#           ifdef POISEUILLE
           call constant_force() ! Poseuille flow generation
#           endif

           ! -----  Update  velocities

           call verlet_velocities()

#ifdef DPD_VV                     

!Note: this recalculates Fd with the new velocities and updates F for the begining og the next cycle 
!       with this new value.  

           call new_dpd_fd()  
#endif


!----  Observe system after equilibration
           if(i_time.gt.n_relax) then 
               call observation 
               if(mod(i_time,n_safe_fftw).eq.0) then
!WARN                   
!DEB test without                  call chain_fftw(2)
               end if
           end if
!----  Make safety copies  to recover from crashes and write out of configurations
        
           if(mod(i_time,n_safe).eq.0) then
               call store_config(2)  ! writes out conf_xmol and conf_new
#           if SYSTEM == 4
!               call chain_fftw(2) ! this initialize the cumulant for the mean fftw
#              if RINGS != 0
                  call ring_net_force(3) ! saving forces in ring_force.dat
#              endif
#           endif

#       if STORE == 0
               call store_config(3)  ! Writes out film_xmol and vel.dat
#       elif STORE == 1
               call store_config(4)  ! Writes out conf_unfold_new (for system == 4), film_xmol and vel.dat UNFOLDED
#       endif
           end if

       end do   ! --------------  ENDS TIME LOOP ----------------

       call obser_out()  
       close(20) ! closing mfa_output

        call get_walltime(t1)
        print *,' WALL TIME (s)= ',t1-t0
        print *,' WALL TIME (min)= ',(t1-t0)/60
        print *,' WALL TIME (hs)= ',(t1-t0)/3600
  end program md_pb

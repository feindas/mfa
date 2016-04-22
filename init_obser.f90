      subroutine init_obser
      use commons 
      use functions 

#include "control_simulation.h"
! OBSOLETE #ifndef PROFILES
! OBSOLETE ! conflicts with vel_prof from functions.f90 
! OBSOLETE       use util , only: velocity_prof
! OBSOLETE #else
#ifdef PROFILES
      use functions
#endif      
      implicit none
      integer :: i_chain_m
      character (len=4) , parameter , dimension(8) :: nam1= (/"Cl  ","O   ","He  ","N   ","S   ","Al   ","H   ","Cu   "/)
      logical :: ufile 
      real (kind=8) :: rgx
!

      count_obs = 0
!
      v_wall_wall_1 = 0.
      v_wall_wall_2 = 0.
      t_wall_1 = 0.
      t_wall_2 = 0.
      e_wall_wall_1 = 0.
      e_wall_wall_2 = 0.
!
      v_fluid_fluid_1 = 0.
      v_fluid_fluid_2 = 0.
      t_fluid_1 = 0.
      t_fluid_2 = 0.
      e_fluid_fluid_1 = 0.
      e_fluid_fluid_2 = 0.
!
      v_fluid_wall_1 = 0.
      v_fluid_wall_2 = 0.
!
      v_total_1 = 0.
      v_total_2 = 0.
      t_total_1 = 0.
      t_total_2 = 0.
      e_total_1 = 0.
      e_total_2 = 0.
!clau: intra
      v_tot_intra_1 =   0.      
      v_tot_intra_2 =   0.
      v_tot_intra_d_1 = 0.
      v_tot_intra_d_2 = 0.

! Total momentum of the system 

      tot_P(:)   = 0.
      tot_P_2(:) = 0.


      time_ave_count = 0
      time_ave_count_2 = 0

!! obso      if(f_logscale.eq.0) then
!! obso          n_time_ave = 200 !if not log scale counts first 200 steps before calculating things in obser
!! obso      end if
!! obso      if(f_logscale.eq.1) then
!! obso       n_time_ave = 0
!! obso      end if

      do i_dim = 1,n_dim
       r0_twall_1(i_dim) = 0.
       f0_spring_1(i_dim) = 0.
      end do
!
      n_fcell_max = 0
      n_wcell_max = 0
!
      n_neigh_fl_max = 0
      n_neigh_fw_max = 0
      n_neigh_ww_max = 0
!
      do i_dim = 1,n_dim
       r0_last_pinned(i_dim) = r0_twall(i_dim)
      end do
!
      z_1 = 0.
      z_2 = 0.
! Radius of Gyration mean variables
        r_g2_mean(:,:) = 0.0
        r_par_per2(:) = 0.0
! End to end radius:
        re_2_mean(:,:) = 0.0
! Force on heads variable 
        f_on_heads(:,:) = 0.0 


! Initialization of store_config routine
 
 call store_config(1)
 
 ! *** For diffusive quantities, drop the initial reference configuration 
 !      of the run (after a first equilibration run)
    
      if( s_time == 1 ) then
print '(/a/)',"  *  Writing the very first generated configuration to conf0.xyz file" 
        open(unit=23,file='conf0.xyz',status='unknown')
          write(23,*) n_part !ORI n_mon_tot
          write(23,*) "    First configuration" 
!
          do i_part = 1 , n_part
          write(23,'(a4,3f15.6)') nam1(a_type(i_part)),r0(:,i_part) 
          end do
          write(23,*) 

        close(23)
      end if

#ifdef PROFILES
      print  '(/a/)',"  *  Calculating profiles during the run" 
!
!  --- Parameters
!
      dim_prof_dens = 200 
      inv_n_mon = 1./dble(n_mon)
      inv_n_chain = 1./dble(n_chain)
      
#if SYSTEM == 0 || SYSTEM == 2 || SYSTEM == 3 /* channel or charges*/
!  --- Allocation       
!
#ifndef PARTICLE_4
        allocate ( v_prof(dim_prof_dens,3,2),v_prof_2(dim_prof_dens,3,2),histo_out(dim_prof_dens,2) )
#else
        allocate ( v_prof(dim_prof_dens,3,3),v_prof_2(dim_prof_dens,3,3),histo_out(dim_prof_dens,3) )
#endif
!
!  --- Init
!
!       call histo(0,dim_prof_dens,n_mon,n_chain,n_mon_d,n_chain_d,surface,boundary(3)) !,histo_out) 

        call histo(0,dim_prof_dens)

       zz_min = 0.

       call  velocity_prof(0,dim_prof_dens,zz_min,boundary(3)) !,v_prof(:,:,:),v_prof_2(:,:,:))

#endif /*channel or charges*/

#if SYSTEM==1 /*droplet*/

! Init droplet measuring variables

        vcm_d(:) = 0.   ! droplet CM velocity
        vcm_d2(:)= 0.
 
! Files for droplet CM's position and velocity 
! common var drop_cm = droplet CM position 

! Find a good value for the droplet CM
!  Also defined the first value for N_CM (which third of the box will contain the CM)

       call rebuild_drop2(0)
       call rebuild_drop2(1)

       guessed_drop_cm(:) = drop_cm(:) 
       print '(a,i/)' ,"        DROP CM should be in the third number = ",n_cm 

!       call write_conf(2,r0,10)
!       stop

! ******************************* IMPORTANT ********************************************

! Profile and initialization from mide_drop_force

! Init. size variables for profiles'

        Lx_o_2(1) = boundary(1) ! For translation of the droplet to comoving frame 
        Lx_o_2(2:3) = 0.

        binning_box_length = 1. ! approx. dimension of the binning box in units of sigma
        
        n_box(2:3) = int (boundary(2:3)/binning_box_length) ! number of boxes in each direction
        n_box(1)=   int (2.*boundary(1)/binning_box_length) ! 
        
        n_box_xz(:) = n_box(:) ; n_box_xz(2) = 1

! Note:  In x dir, we add  boxes needed because of the unfolding. We double the number of boxes  

! Range covered by the binning boxes of the histograms [DEFINED. AVOID CHANGING THEM]. 
! They are function of n_box

        r_box_xz(1) = 2.*boundary(1)/dble(n_box(1)) ! For proper copy-pasting, the range in x is twice the box size
        r_box_xz(2) = boundary(2)
        r_box_xz(3) = boundary(3)/dble(n_box(3))

! Allocate vectors 

      allocate ( dens_xz(n_box(1),1,n_box(3)) , dens_xz_2(n_box(1),1,n_box(3)) )
      allocate (  dens_xz_b(n_box(1),1,n_box(3)) ,  dens_xz_b_2(n_box(1),1,n_box(3)) )
      allocate (  r0_unfold(3,n_part), r_cms(n_chain+n_chain_d,3) )
      allocate ( vel_xz(n_box(1),n_box(3),3), vel_xz_b(n_box(1),n_box(3),3)  )
      allocate ( vel_xz2(n_box(1),n_box(3),3), vel_xz_b2(n_box(1),n_box(3),3)  )

!! NOTE YET #ifdef FORCE_PROFILE      
!! NOTE YET       allocate ( force_prof(n_box(1),3) )
!! NOTE YET #endif

! Init. routines (all in util.f90 )

      call dens_prof(0,dens_xz(:,:,:),dens_xz_b(:,:,:),dens_xz_2(:,:,:),dens_xz_b_2(:,:,:),n_box_xz, &
      r_box_xz) ! end 
        
      call vel_prof_2d(0,vel_xz,vel_xz2,vel_xz_b,vel_xz2,n_box_xz,r_box_xz,dens_xz) ! 


        print '(/a,3f10.5)', "First drop CM = " ,  guessed_drop_cm(:)
        print '(a/)', "  **** CHECK IF IT IS NOT TOO BAD ****** " 

!  Initialize 

      call calc_cms(1) 

      call rebuild_drop(0,rgx) 
      
! Define fix Z coordinate shift     

      
      fix_shift = boundary(3)/2.

      print '(/a,f10.7/)', "  * Constant shift in drop profile in z of ",fix_shift

#endif /*system=1: droplet*/

#endif /*profiles*/

#if STORE == 1 /* I/O refold  */
    print '(/a/)', "  * Storing unfolded coordinates "
!
!  Read unfold conf
!
        inquire(file='conf_unfold',exist=ufile )
        if(ufile) then    ! if conf_unfold exists, read it  
            print '(/a/)' ,'  *  Reading file conf_unfold '
            open (unit=180,file='conf_unfold')
            do i_part = 1 , n_part
            read(180,'(3e16.8)') r0_unfold(:,i_part)
            end do
            close (180)
#           if SYSTEM == 4
                vec_dummy(:) = 0.
                r_dummy = 0.
                do i_part = 1 , n_mon*n_chain
                  vec_dummy(:) = vec_dummy(:) + mass(i_part)*r0_unfold(:,i_part)
                  r_dummy = r_dummy + mass(i_part)
                end do
                  vec_dummy(:) = vec_dummy(:)/r_dummy
                  print '(a,3f16.5)', "[init_obser] S4 from conf_unfold Rcm(chain): ", vec_dummy(:)   
#           endif            
        else ! if conf_unfold does not exist, generate the first  unfolded configuration 
    
        print '(/a/)' ,'  *  conf_unfold not present in dir. Unfolding first conf and starting from here '
    
! Init refold 
    
#if SYMMETRY == 0 /* channel */
        call refold(0,n_mon,1,part_init_d) 
        call refold(1,n_mon,1,part_init_d) 
        call refold(1,n_mon_d,1+part_init_d,part_init_e)

#ifdef PARTICLE_4
        call refold(1,n_mon_e,1+part_init_e,n_mon_tot) 
#endif

#elif SYMMETRY == 1 /*bulk */
#       if SYSTEM == 4
! WARNING! reading this lines if conf_unfold doesn't exists

! if SYSTEM 4 and STORE == 1 
! ************************** !
! WARNING! if you call refold(1,n_mon,1,part_init_d) 
! in r0_unfold you only have the chain without rings
        call refold(0,n_mon,1,part_init_d) 
        call refold(1,n_mon,1,part_init_d) ! this line is for chain
! ************************************************* !
! Now, Creating r0_unfold from r0. then fixing Rcm
#if RINGS == 1
! Now refolding the rings 
        r0_unfold(:,part_init_d+1:n_mon_tot) = r0(:,part_init_d+1:n_mon_tot) ! copy the rings

        ! Now looking for the closest bead of the chain to the ring 1
        j_part = int( abs( (r0(1,1)-r0(1,part_init_d+1))*n_mon*inv_boundary(1) ) )
        ! j_part is aprox. the number of beads between ring 1 and chain's bead 1 

        if (n_ring.eq.1) then
     ! here j_part is the bead closest (in X direction) to the ring 1
     ! Now the rigid shift
          if (r0(1,part_init_d+1) .lt. r0(1,1) ) then
             print '(//a,i3)', "[init_obser] S4 - Bead near to ring 1: ", n_mon-j_part
             do while( abs( r0(1,n_mon-j_part)-r0(1,part_init_d+1) ).lt.abs( r0(1,n_mon-j_part+1)-r0(1,part_init_d+1) ) )
             j_part = j_part + 1
             end do
             j_part = n_mon-j_part+1 ! The key bead
             print '(a,i3)', "[init_obser] S4 - Bead nearest to ring 1: ", j_part

              print*, "[init_obser] S4 - Moving X position of ring 1"
              do i_part =part_init_d+1, n_mon_tot
               r0_unfold(1,i_part)=r0_unfold(1,i_part) + boundary(1) ! to do in x directions
               r0_unfold(2,i_part)=r0_unfold(2,i_part)-r0(2,j_part)+r0_unfold(2,j_part)
               r0_unfold(3,i_part)=r0_unfold(3,i_part)-r0(3,j_part)+r0_unfold(3,j_part)
              end do

          else
              print '(a,i3)', "[init_obser] S4 - Bead near to ring 1: ", j_part
              do while( abs( r0(1,j_part)-r0(1,part_init_d+1) ).gt.abs( r0(1,j_part+1)-r0(1,part_init_d+1) ) )
              j_part = j_part + 1
              end do
              print '(a,i3)', "[init_obser] S4 - Bead nearest to ring 1: ", j_part
           ! The key bead is j_part
              do i_part =part_init_d+1, n_mon_tot
               r0_unfold(2,i_part)=r0_unfold(2,i_part)-r0(2,j_part)+r0_unfold(2,j_part)
               r0_unfold(3,i_part)=r0_unfold(3,i_part)-r0(3,j_part)+r0_unfold(3,j_part)
              end do
          end if

        else ! n_ring == 2
        ! j_part is aprox. the number of beads between ring 1 and chain's bead 1 
        ! Now the rigid shift for ring 1
          if ( r0(1,part_init_d+1).lt.r0(1,1) ) then
             print '(a,i3)', "[init_obser] S4 - Bead near to ring 1: ", n_mon-j_part
             do while( abs( r0(1,n_mon-j_part)-r0(1,part_init_d+1) ).lt.abs( r0(1,n_mon-j_part+1)-r0(1,part_init_d+1) ) )
                  j_part = j_part + 1
             end do
             j_part = n_mon-j_part+1 ! The key bead
             print '(a,i3)', "[init_obsre] S4 - Bead nearest to ring 1: ", j_part
  
             do i_part =part_init_d+1, part_init_e
              r0_unfold(1,i_part)=r0_unfold(1,i_part) + boundary(1) ! X directions
              r0_unfold(2,i_part)=r0_unfold(2,i_part)-r0(2,j_part)+r0_unfold(2,j_part)
              r0_unfold(3,i_part)=r0_unfold(3,i_part)-r0(3,j_part)+r0_unfold(3,j_part)
             end do
          else
             print '(a,i3)', "[init_obser] S4 - Bead near to ring 1: ", j_part
             do while( abs( r0(1,j_part)-r0(1,part_init_d+1) ).gt.abs( r0(1,j_part+1)-r0(1,part_init_d+1) ) )
                  j_part = j_part + 1
             end do
             print '(a,i3)', "[init_obser] S4 - Bead nearest to ring 1: ", j_part
             ! The key bead is j_part 
             do i_part =part_init_d+1,part_init_e 
              r0_unfold(2,i_part)=r0_unfold(2,i_part)-r0(2,j_part)+r0_unfold(2,j_part)
              r0_unfold(3,i_part)=r0_unfold(3,i_part)-r0(3,j_part)+r0_unfold(3,j_part)
             end do
          end if

        ! Now looking for the closest bead of the chain to the ring 2 
        j_part = int( abs( (r0(1,1)-r0(1,part_init_e+1) )*n_mon*inv_boundary(1) ) )
        ! j_part is aprox. the number of beads between ring 2 and chain's bead 1 
          if (r0(1,part_init_e+1) .lt. r0(1,1) ) then
             print '(a,i3)', "[init_obser] S4 - Bead near to ring 2: ", n_mon-j_part
             do while( abs( r0(1,n_mon-j_part)-r0(1,part_init_e+1) ).lt.abs( r0(1,n_mon-j_part+1)-r0(1,part_init_e+1) ) )
                  j_part = j_part + 1
             end do
             j_part = n_mon-j_part+1 ! The key bead
             print '(a,i3)', "[init_obser] S4 - Bead nearest to ring 2: ", j_part
  
             do i_part =part_init_e+1, n_mon_tot
              r0_unfold(1,i_part)=r0_unfold(1,i_part) + boundary(1) ! X directions
              r0_unfold(2,i_part)=r0_unfold(2,i_part)-r0(2,j_part)+r0_unfold(2,j_part)
              r0_unfold(3,i_part)=r0_unfold(3,i_part)-r0(3,j_part)+r0_unfold(3,j_part)
             end do
          else
             print '(a,i3)', "[init_obser] S4 - Bead near to ring 2: ", j_part
             do while( abs( r0(1,j_part)-r0(1,part_init_e+1) ).gt.abs( r0(1,j_part+1)-r0(1,part_init_e+1) ) )
                  j_part = j_part + 1
             end do
  
             ! The key bead is j_part it self
             print '(a,i3/)', "[init_obser] S4 - Bead nearest to ring 2: ", j_part
             do i_part =part_init_e+1,n_mon_tot
              r0_unfold(2,i_part)=r0_unfold(2,i_part)-r0(2,j_part)+r0_unfold(2,j_part)
              r0_unfold(3,i_part)=r0_unfold(3,i_part)-r0(3,j_part)+r0_unfold(3,j_part)
             end do
          end if
     end if
#endif /* RINGS == 1*/

#   if FIXCM == 1
    call fix_CM() ! Rigid shift all the system (chain+rings) in r0_unfold
! Checking the CM of the chain
    vec_dummy(:) = 0.
    r_dummy = 0.
    do i_part = 1 , n_mon*n_chain
      vec_dummy(:) = vec_dummy(:) + mass(i_part)*r0_unfold(:,i_part)
      r_dummy = r_dummy + mass(i_part)
    end do
      vec_dummy(:) = vec_dummy(:)/r_dummy
      print '(a,3f16.5)', "[init_obser,FIXCM] From r0_unfold Corrected Rcm (chain)=", vec_dummy(:)
#   endif

#       else /* system not 4 */
        call refold(0,n_mon_d,1,part_init_e) 
        call refold(1,n_mon_d,1+part_init_d,part_init_e)
#           ifdef PARTICLE_4
        call refold(1,n_mon_e,1+part_init_e,n_mon_tot) 
#           endif
#       endif

#endif /* symmetry */
        
        endif ! close, if conf_unfold does not exist
#endif /* STORE == 1*/

#if STORE == 0 && SYSTEM == 4
    	call refold(0,n_mon,1,part_init_d) 
        call refold(1,n_mon,1,part_init_d) ! this line is for chain + rings refold all 
        print *, " * Refolding for STORE = 0, chain... ok"

! if SYSTEM == 4 and STORE == 0 
! ****************************************************************************** !
! WARNING! refold copy r0 in r0_unfold from 1 to n_mon_tot 
! this means that if you call refold(1,n_mon,1,part_init_d) 
! in r0_unfold you only have the chain without rings
        call refold(0,n_mon,1,part_init_d) 
        call refold(1,n_mon,1,part_init_d) ! this line is for chain

        vec_dummy(:) = 0.0 
        do i_part = 1 , n_mon*n_chain
            vec_dummy(:) = vec_dummy(:) + mass(i_part)*r0_unfold(:,i_part)
            r_dummy = r_dummy + mass(i_part)
            vec_dummy(:) = vec_dummy(:) + mass(i_part)*r0_unfold(:,i_part)
            r_dummy = r_dummy + mass(i_part)
        end do

      vec_dummy(:) = vec_dummy(:)/r_dummy
      print '(a,3f16.5)', "[init_obser,STORE==0] S4 - Rcm (chain)=", vec_dummy(:) 

#if RINGS == 1
        print *, " * Refolding for STORE = 0, now rings"
! Now refolding the rings 
        r0_unfold(:,part_init_d+1:n_mon_tot) = r0(:,part_init_d+1:n_mon_tot) ! copy the rings
        ! first step looking for the closest bead of the chain to the ring 1
        i_part = 1
        j_part = int( abs( (r0(1,i_part)-r0(1,part_init_d+1))*n_mon*inv_boundary(1)  ) )
        do while ( abs(r0(1,j_part)-r0(1,part_init_d+1)).gt.abs(r0(1,j_part+1)-r0(1,part_init_d+1)) )
             j_part = j_part + 1
        end do
        if (n_ring.eq.1) then
        ! here j_part is the bead closest (in X direction) to the ring 1
        ! Now the rigid shift 
        do i_part =part_init_d+1, n_mon_tot
            do i_dim=1,3
            r0_unfold(i_dim,i_part)=r0_unfold(i_dim,i_part)-r0(i_dim,j_part)+r0_unfold(i_dim,j_part)
            end do
        end do

        else
        ! here j_part is the bead closest (in X direction) to the ring 1
        ! Now the rigid shift
        do i_part =1+part_init_d, n_mon_e
          do i_dim=1,3
          r0_unfold(i_dim,i_part)=r0_unfold(i_dim,i_part)-r0(i_dim,j_part)+r0_unfold(i_dim,j_part)
          end do
        end do
          i_part = 1
          j_part = int( abs( (r0(1,i_part)-r0(1,part_init_d+1))*n_mon*inv_boundary(1)  ) )
          do while ( abs(r0(1,j_part)-r0(1,part_init_d+1)).gt.abs(r0(1,j_part+1)-r0(1,part_init_d+1)) )
                j_part = j_part + 1
          end do
          ! here j_part is the bead closest (in X direction) to the ring 2
          ! Now the rigid shift
           do i_part =1+part_init_e, n_mon_tot
               do i_dim=1,3
               r0_unfold(i_dim,i_part)=r0_unfold(i_dim,i_part)-r0(i_dim,j_part)+r0_unfold(i_dim,j_part)
               end do
           end do
        end if
        print *, " * Refolding for STORE = 0, rings... ok"
#endif /* RINGS == 1*/
#   if FIXCM == 1
    call fix_CM() ! Rigid shift all the system (chain+rings) in r0_unfold
! Checking the CM of the chain
    vec_dummy(:) = 0.
    r_dummy = 0.
    do i_part = 1 , n_mon*n_chain
      vec_dummy(:) = vec_dummy(:) + mass(i_part)*r0_unfold(:,i_part)
      r_dummy = r_dummy + mass(i_part)
    end do
      vec_dummy(:) = vec_dummy(:)/r_dummy
      print '(a,3f16.5)', "[init_obser,FIXCM] From r0_unfold Corrected Rcm (chain)=", vec_dummy(:)
#   endif

#endif /* STORE == 0 & SYSTEM == 4 */

#if SYSTEM == 0
#       ifdef PARTICLE_4
!        Mean velocity of particle 4 

               v4_mean(:) = 0.0 
               v4_mean2(:) = 0.0 
                 c4 = 0 

#       endif
#endif
#        ifdef DIFF
!       
!        Init diffusion Routine (only for particle 4 )
!       

             call diff_coef(1,r_time)
#       endif            

#if SYSTEM == 4
        call chain_fftw(1)
        call bond_distance(0)
#       if RINGS != 0
        call ring_net_force(0) ! Openning file ring_force.dat
#       endif /* RINGS != 0 */

!deb
        !fr1(:) = 0.0
        !fr2(:) = 0.0
        fr_CM(:) = 0.0 ! force in CM chain

#endif

end subroutine init_obser

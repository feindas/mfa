  subroutine gen_brush(mode)
      use commons ;use util ; implicit none
      integer, intent(in) :: mode

!---  Specify position of first monomer in chain
#include "control_simulation.h"
  select case (mode)

  case(1)  ! ------ Brush in top and bottom walls: CHANNEL

! NOTE: This works with brushes only. Free polymers must be generated in another
! context       

    do i_chain = 1,n_chain  !loop over all the chains
        i_part = 1 + (i_chain-1)*n_mon
        call r250(mz,random,n_part,n_dim,kptr)
        do i_dim = 1,n_dim
            if(i_dim.eq.n_dim) then             ! z coor of polymer head
                if(i_chain.le.n_chain/2) then !top wall
                    r0(i_dim,i_part) = z_space_wall-z_head     !ori2.**(1./6.)
                end if
                if(i_chain.gt.n_chain/2) then !bottom wall
                    r0(i_dim,i_part) = z_head !ori2.**(1./6.)
                end if
            end if

            if(i_dim.ne.n_dim) then                 ! x y coordinates of heads: at random in the plane
                r0(i_dim,i_part) = random(i_dim)*boundary(i_dim)
            end if
        end do

        !----  generate random walk for next monomers on chain

        do i_mon = 1,n_mon-1
            i_part = i_part + 1
500   continue
            call r250(mz,random,n_part,n_dim,kptr)

            do i_dim = 1,n_dim
                if(i_dim.lt.n_dim) then ! for x and y coordinates, random
                    r0(i_dim,i_part) = r0(i_dim,i_part-1) +    (2*random(i_dim)-1.)*r_chain/sqrt3
                    ! PBC in the plane          
                    if(r0(i_dim,i_part).ge.boundary(i_dim)) then
                        r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                    else if(r0(i_dim,i_part).le.0.) then
                        r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                    end if
                end if !end if x or y 
                if(i_dim.eq.n_dim) then             ! if z coor 

                    if(i_chain.le.n_chain/2) then  ! top wall
#            if SOLVENT == 1 || SOLVENT == 2 
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.4
#           endif                

                    end if  ! end z_coor 
                    if(i_chain.gt.n_chain/2) then !bottom wall

#            if SOLVENT == 1 || SOLVENT == 2 
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.4
#           endif                
     !OBSOLETE                if(solv_flag.eq.0) then
     !OBSOLETE                r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
     !OBSOLETE                else
     !OBSOLETE                r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.4
     !OBSOLETE                end if

                    end if
                end if !end z coor
                !      end if !ad_flag

                !!    if(ad_flag.eq.0) then !if adsorbed
                !!             r0(i_dim,i_part) = r0(i_dim,i_part-1) +                       &
                !!        &    (2*random(i_dim)-1.)*r_chain/sqrt3
                !!             if(r0(i_dim,i_part).gt.boundary(i_dim)) then
                !!              r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                !!             else if(r0(i_dim,i_part).lt.0) then
                !!              r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                !!            end if
                !!        end if !end adsorbed

            end do
            i_dummy = 0
            ! check whether new z coordinates is allright

            if((r0(n_dim,i_part).le.0.).or.(r0(n_dim,i_part).ge.boundary(n_dim)))  i_dummy = 1
            if(i_dummy.eq.1) then
                write(*,*) " New setup attempt for monomer (z-pos)",  i_part,r0(n_dim,i_part)
                goto 500
            end if
        end do

    end do  ! ends loop over chains

!debug        call write_conf(1,r0(:,1:n_chain*n_mon),10) ; stop   

case(2)  ! ------ Brush ONLY in bottom wall: DROPLET

    do i_chain = 1,n_chain  !loop over all the chains
        i_part = 1 + (i_chain-1)*n_mon
        call r250(mz,random,n_part,n_dim,kptr)
        do i_dim = 1,n_dim
            if(i_dim.eq.n_dim) then             ! z coor of polymer head
!
! DROPLET: only the bottom wall is populated 
!
                r0(i_dim,i_part) = z_head !ori2.**(1./6.)
            end if
            if(i_dim.ne.n_dim) then                 ! x y coordinates of heads: at random in the plane
                r0(i_dim,i_part) = random(i_dim)*boundary(i_dim)
            end if
        end do

!----  generate random walk for next monomers on chain

        do i_mon = 1,n_mon-1
            i_part = i_part + 1
            501   continue

        call r250(mz,random,n_part,n_dim,kptr)

        do i_dim = 1,n_dim
! PBC in the plane          
                if(r0(i_dim,i_part).ge.boundary(i_dim)) then
                r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                else if(r0(i_dim,i_part).le.0.) then
                 r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                end if
!assumed            end if !emd if x or y 
            if(i_dim.eq.n_dim) then             ! if z coor 
#            if SOLVENT == 1 || SOLVENT == 2
                   r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                   r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.5
#            endif
            end if
      !not in drop       end if !end z coor
!          end if
!!obs       if(ad_flag.eq.0) then !if adsorbed
!!obs                r0(i_dim,i_part) = r0(i_dim,i_part-1) + (2*random(i_dim)-1.)*r_chain/sqrt3
!!obs                if(r0(i_dim,i_part).gt.boundary(i_dim)) then
!!obs                 r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
!!obs                else if(r0(i_dim,i_part).lt.0) then
!!obs                 r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
!!obs               end if
!!obs           end if !end adsorbed

        end do
         i_dummy = 0
      ! check whether new z coordinates is allright
        if((r0(n_dim,i_part).le.0.).or.(r0(n_dim,i_part).ge.boundary(n_dim)))  i_dummy = 1
         if(i_dummy.eq.1) then
          write(*,*) "new setup attempt for monomer (z-pos)",  i_part,r0(n_dim,i_part)
          goto 501
         end if
        end do
        end do  ! ends loop over chains

    end select
end subroutine gen_brush

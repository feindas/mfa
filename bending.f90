! Routine that allows to add rigidity to the grafted chains
subroutine bending(mode_bend)
#include "control_simulation.h"
use commons
implicit none
integer, intent(in) :: mode_bend
integer ::i,j,l
real(kind=8) ::r_prev(3),r_next(3),dot_pr=0,r_prev_2=0,r_next_2=0, cos_alpha,dir_next(3)=0,dir_prev(3)=0,dir_prev_2=0, dir_next_2=0, F_bend(6),F_mod,delta_alpha
#ifdef BENDING
select case (mode_bend)

case(0)  ! Init  variables 
print*,""
print*," * Simulation with bending stiffness"
print*," * Bending elastic constant ",k_bend
print*," * Equilubrium angle:",alpha_eq
print*,""

case(1)
!SYMETRY == 0 dont correcto for PBC in 3D

case(2)
!SYMETRY == 1 dont correcto for PBC in 3D

!!!!         !  debug
!!!!         print*,"r0",i_time,sum(r0(:,:))/dble(n_part)
!!!!         print*,"v", i_time,sum(v(:,:))/dble(n_part)
!!!!         print*,"a", i_time,sum(a(:,:))/dble(n_part)
!!!!         print*,"f", i_time,sum(force(:,:))/dble(n_part)
!!!!         ! end debug

v_bend=0.  !Reset Bending Potential Energy
Do l=0,n_chain-1 !n_chain , loop over chains
    r_next=r0(:,l*n_mon+2)-r0(:,l*n_mon+1)
    !This lines below correct for perioduc boundry conditions
    r_next(1) = r_next(1) - boundary(1) * int(2.*r_next(1)*inv_boundary(1))
    r_next(2) = r_next(2) - boundary(2) * int(2.*r_next(2)*inv_boundary(2))
    !------ get boundaries right. PBC in x, y, z
    r_next(3) = r_next(3) - boundary(3) * int(2.*r_next(3)*inv_boundary(3))


    !End correction fo periodic boundry conditions
    Do i=2,n_mon-1 !n_mon, loop over particles in chain
        !Reset dummy variables
        dot_pr=0
        r_prev_2=0
        r_next_2=0
        dir_prev_2=0
        dir_next_2=0
        !Calculation of bending forces for the neighbors of i
        r_prev=r_next
        r_next=r0(:,l*n_mon+i+1)-r0(:,l*n_mon+i)
        !This lines below correct for perioduc boundry conditions
        r_next(1) = r_next(1) - boundary(1) * int(2.*r_next(1)*inv_boundary(1))
        r_next(2) = r_next(2) - boundary(2) * int(2.*r_next(2)*inv_boundary(2))
        !------ get boundaries right. PBC in x, y, z
        r_next(3) = r_next(3) - boundary(3) * int(2.*r_next(3)*inv_boundary(3))
        !End correction fo periodic boundry conditions
        Do j=1,3
            dot_pr=dot_pr+r_next(j)*r_prev(j)
            r_prev_2=r_prev_2+r_prev(j)*r_prev(j)
            r_next_2=r_next_2+r_next(j)*r_next(j)
        End do
        cos_alpha=dot_pr/sqrt(r_prev_2*r_next_2)
        dir_next=r_prev*r_next_2-r_next*dot_pr
        dir_prev=-r_next*r_prev_2+r_prev*dot_pr
        Do j=1,3
            dir_prev_2=dir_prev_2+dir_prev(j)*dir_prev(j)
            dir_next_2=dir_next_2+dir_next(j)*dir_next(j)
        End do
        if(dir_next_2.lt.0.0000001) then !if alpha is small make force 0
            dir_next=dir_next*0
           !print*, "ENTRO ACA1"
        else
        !Below I divide by |r_next|, because the bending force is proportional to 
        !k_bend*delta_alpha/|r_next|. This comes from V_vend=1/2*k_bend*delta_alpha_2
            dir_next=dir_next/sqrt(dir_next_2*r_next_2)
        end if
        if(dir_prev_2.lt.0.0000001) then !if alpha is small make force 0
            dir_prev=dir_prev*0
           !print*, "ENTRO ACA2"
        else
        !Below I divide by |r_prev|, because the bending force is proportional to 
        !k_bend*delta_alpha/|r_prev|. This comes from V_vend=1/2*k_bend*delta_alpha_2
            dir_prev=dir_prev/sqrt(dir_prev_2*r_prev_2)
        end if
        !debug
        if(cos_alpha.ge.1.0) then
            print*, "Error cos_alpha =",cos_alpha, " >= 0"
            delta_alpha=0 
        else
            delta_alpha=acos(cos_alpha)-alpha_eq
        end if
        !debug
        F_mod=k_bend*delta_alpha
        v_bend=v_bend+.5*F_mod*delta_alpha
        Do j=1,3
            F_bend(j)=F_mod*dir_prev(j)
            F_bend(3+j)=F_mod*dir_next(j)

            !debug
            !if(F_bend(j).ne.F_bend(j)) then
            !    print*,"Error found, monomer",i,", component",j
            !    print*,"F_mod",F_mod,", dir_prev(j)",dir_prev(j),", delta_alpha", delta_alpha,", cos_alpha ",cos_alpha
            !end if
            !if(F_bend(3+j).ne.F_bend(3+j)) then
            !    print*,"Error found, monomer",i,", component",j+3
            !    print*,"F_mod",F_mod,", dir_next(j)",dir_next(j)
            !end if
            !debug

        End do
        Do j=1,3
            force(j,l*n_mon+i-1) = force(j,l*n_mon+i-1) + F_bend(j)
            force(j,l*n_mon+i)   = force(j,l*n_mon+i) - F_bend(j) - F_bend(3+j)
            force(j,l*n_mon+i+1) = force(j,l*n_mon+i+1) + F_bend(3+j)
        End do
        !!!!!!!!--------------DEBUGGING STUFF--------------------!!!!!!!
        !  TEST TO CHECK IF THE ROUTINE IS CALCULATING FORCES AND ENERGY
        !  CORRECTLY. 
        !print*,"Inside chaind Bending Forces"
        ! print*, F_bend 
        !print*, v_bend, "bending ENERGY"
        ! ERASE AFTER CHECKING
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
    End Do
End do

 !!!!!!!!---------------DEBUGGING STUFF--------------------!!!!!!!!!!!!!!
        !  TEST TO CHECK IF THE ROUTINE IS CALCULATING FORCES AND ENERGY
        !  CORRECTLY. 
        !print*,"f", i_time,sum(force(:,:))/dble(n_part)
        !print*, "positions"
        !print*, r0(1,1)-1.,r0(2,1),r0(3,1)
        !do i=1,n_mon
        !print*,r0(:,i)
        !end do
        !print*, r0(1,n_mon)+1.,r0(2,n_mon),r0(3,n_mon)
        !print*, r0(:,2)
        !print*, r0(:,n_mon-1)
        !print*, r0(:,n_mon)
        !print*, k_bend, "bending constant"
        ! ERASE AFTER CHECKING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end select

#endif /*close BENDING*/
end subroutine bending

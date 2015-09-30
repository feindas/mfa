subroutine chain_fftw(mode)
    use commons

 implicit none
 include 'fftw3.f'
! include 'fftw.f'

integer, intent(in) :: mode
integer ::  n,i,j
real (kind = 8) :: PIii,amp,peri
real (kind = 8) :: lambda_fftw
real (kind = 8) :: n_box_fftw_real
!integer, save (kind=8) :: plan_fftw, j_fftw = 0, i_fftw !con esto ni compila


select case (mode)
!Initialize variables
    case(1) ! this mode is called from init_obser.f90 only. 

! h_tot_fftw(j) will be the total height in the j-th box

! now start the inicialization variables
       call read_mfa() 
       call store_fft_parameters(0)

!  print *, "Este programa calcula bien la altura de las particulas respecto del CM. "

        allocate ( h_fftw(n_mon), h2p_fftw(n_mon), cajitas_fftw(n_box_fftw), h_tot_fftw(n_box_fftw) ) 
        allocate (Ipar_fftw(n_box_fftw/2+1),Imean_fftw(n_mon/2+1)) 
        allocate (Amp_fftw(n_box_fftw/2+1),Im_fftw(n_box_fftw/2+1))
        Imean_fftw(:) = 0.
        Amp_fftw(:) = 0.
        Im_fftw(:) = 0.
        
        call fftw_start(0)
       

!****** end of initialization variables **** 

    case(2) !****** start the FFT process ******
    j_fftw = j_fftw + 1
!    print *, "Times the FFT process is calculated:", j_fftw
    flag_bad_screen = .FALSE.

	    call create_chain_height() 

!	    print *, "calculando la chain_fftw, j: ", j_fftw 
	    call chain_fft() !calculate  mesoscopic chain for the fftw
        ! chain_fft modifie the boolean variable flag_bad_screen, wich indicates if the screen had
        ! a NaN in the h_tot_fftw array.
              
! summing the number of screen where was an empty box. (NaN in h_tot_fftw array)
! if flag_bad_screen is true, then accumulate n_bad_screen
!does not work        if(flag_bad_screen) then
!does not work            n_bad_screen = n_bad_screen + 1 
!does not work            print *, "Number of screen where at least one box was empty (n_bad_screen):", n_bad_screen
!does not work        endif



!	    call fftw_start(j_fftw-1) ! when j == n_screen do fft and destroy plan
	    call fftw_start(2)
            
! now calculate the mean difractogram, and delete Ipar_fftw()

n_box_fftw_real=n_box_fftw

        do i_fftw=1,(n_box_fftw/2+1)
!             if (.not.flag_bad_screen) then
             Imean_fftw(i_fftw) = Imean_fftw(i_fftw) + (REAL(Ipar_fftw(i_fftw))**2 + DIMAG(Ipar_fftw(i_fftw))**2)/(n_box_fftw)

             !integer :: n_box_fftw=1000, n_screen

             Amp_fftw(i_fftw) =   Amp_fftw(i_fftw) + ( REAL(Ipar_fftw(i_fftw)) / SQRT(n_box_fftw_real) )
             Im_fftw(i_fftw) = Im_fftw(i_fftw) + (DIMAG(Ipar_fftw(i_fftw)) / SQRT(n_box_fftw_real)   )

!             end if
           Ipar_fftw(i_fftw) = 0.
        enddo
        
    case(3)
    ! Now store the mean difractogram,
	call save_mean_difrac()
!    print *, " * difrac_promedio is saved"
    print *, "Times the FFT process is calculated:", j_fftw
    case default
    print *, "Ok! something is wrong with or everything is ok! check that 4 &
    >", mode
    call fftw_start(1)
 
end select

CONTAINS

subroutine read_mfa()

       if(n_mon < n_box_fftw) then
          print *, "********** ***************************************************************" 
          print *, "* WARNING: Number of box n_box_fftw is greater than the number of beads! *"
          print *, "********** ***************************************************************" 
          stop
       end if
       delta_x = (boundary(1) )/n_box_fftw
       n_screen = int( (n_relax+n_obser) /n_safe)
!print *, "[read_mfa() 1/3] * chain_fftw: mfa read... ok!"
!print '("[read_mfa() 2/3] * n_screen:" (i4) )', n_screen
!print '("[read_mfa() 3/3] * Lx, Ly, Lz, N: " (4f10.5) )', boundary(1), boundary(2), boundary(3), n_mon
print '("[read_mfa() ] Difracs information:")'
print '("[read_mfa() ] * n_screen:" (i4) )', n_screen
print '("[read_mfa() ] * Lx, Ly, Lz, N: " (3f10.5,x,i3/) )', boundary(1), boundary(2), boundary(3), n_mon

end subroutine read_mfa

subroutine  store_fft_parameters(j)
    integer, intent(in) :: j
    
    open(unit=5, file = 'parametros_fft')!,iostat=ierror )
!	write(5,*) ""
        write(5,"(1a,f10.5)") "[fft]Frec_de_corte: ", n_box_fftw/(2*boundary(1)) 
        write(5,"(1a,f10.5)") "[fft]Delta_de_frecuencia: ", n_box_fftw/(boundary(1)*n_mon)
        write(5,*) "[fft]n_box_fftw: ", n_box_fftw
        write(5,*) "[fft]num_of_beads: ", n_mon
        write(5,"(1a,f10.5)") "[fft]max_x=Lx: ", boundary(1)
        write(5,"(1a,f10.5)") "[fft]Ly: ", boundary(2)
        write(5,*) "[fft]n_screen: ", n_screen
        write(5,"(1a,f10.5)") "[fft]delta_x: ", delta_x

    close(unit=5)
end subroutine store_fft_parameters

subroutine create_chain_height()
     h2pprom = 0
     h2pdesv = 0
	 hprom = 0
     hprom2 = 0
     
    open(unit=14, file = 'h_stat', position = 'append')
! to calculate h_fftw is necessary the r0_unfold variable. [SYSTEM == 4, rings+chain]
    do i_part=1,n_mon 
    ! IMPORTAN_ h_fftw(:) is a global variable needed in subroutine chain_fft()
! originally      h_fftw(i_part)= sqrt( (r0_unfold(2,i_part)-boundary(2)/2)**2 + (r0_unfold(3,i_part)-boundary(3)/2)**2 )
! este h_fftw no es el eje x de difrac promedio, es el eje y
      h_fftw(i_part)= sqrt( (r0_unfold(2,i_part)-boundary(2)/2)**2 + (r0_unfold(3,i_part)-boundary(3)/2)**2 )
	    hprom = hprom + h_fftw(i_part)
	   	hprom2 = hprom2 + h_fftw(i_part)**2
        h2p_fftw(i_part)= (r0_unfold(2,i_part)-boundary(2)/2) + (r0_unfold(3,i_part)-boundary(3)/2)
        h2pprom = h2pprom + h2p_fftw(i_part)
        h2pdesv = h2pdesv + h2p_fftw(i_part)**2

    enddo
    ! j_fftw se actualiza below case(2)
         
        write(14,"(I10, e10.3, e12.3, e12.3, e12.3)") j_fftw, hprom/dble(n_mon), sqrt(hprom2/dble(n_mon) - (hprom/dble(n_mon) )**2), h2pprom/dble(n_mon), sqrt( h2pdesv/dble(n_mon) - (h2pprom/dble(n_mon) )**2 )
    close(unit=14)
end subroutine create_chain_height

subroutine chain_fft()
integer :: p

	cajitas_fftw(:) = 0
	h_tot_fftw(:) = 0
!       
! Here calculate the number of beads per n_box_fftw and the total height in each
! n_box_fftw.
!
       do i_part=1,n_mon    
!******************************************************************
! the height in the box was calculated with r0_unfold       
! while boxes in x direction use r0 (folded configuration)
!******************************************************************
            p = int(r0(1,i_part)/delta_x) + 1 ! calculate the box
            ! increment number of beads in  box p
            cajitas_fftw(p) = cajitas_fftw(p) + 1 
            ! calculate total height  
            h_tot_fftw(p) = h_tot_fftw(p) + h_fftw(i_part)
       enddo

       do i_fftw=1,n_box_fftw
       ! normalize height, be carefull with empty boxes
       if ( cajitas_fftw(i_fftw).eq.0 ) then
                if(i_fftw.eq.1 ) then
!***************************************************************
!the next line is important, maybe forcing modes?
! If found a empty box
                    if(cajitas_fftw(i_fftw+1).eq.0 ) then
                       h_tot_fftw(i_fftw)= h_tot_fftw(n_box_fftw)/cajitas_fftw(n_box_fftw) 
                       else
                       h_tot_fftw(i_fftw)= (h_tot_fftw(n_box_fftw)/cajitas_fftw(n_box_fftw) + h_tot_fftw(i_fftw+1)/cajitas_fftw(i_fftw+1))/2
                    endif
                elseif(i_fftw.eq.n_box_fftw) then
                    !this line is important, maybe forcing modes?
                    h_tot_fftw(i_fftw) = (h_tot_fftw(i_fftw-1) + h_tot_fftw(1))/2
                else     
                    !this line is important, maybe forcing modes?
                    if( cajitas_fftw(i_fftw+1).eq.0 ) then
                    h_tot_fftw(i_fftw)= h_tot_fftw(i_fftw-1) /2
                    else
                    h_tot_fftw(i_fftw)= (h_tot_fftw(i_fftw-1) + h_tot_fftw(i_fftw+1)/cajitas_fftw(i_fftw+1))/2
                    endif
                end if

               ! if h_tot_fftw(i_fftw)  still is a NaN then put it to zero
               ! Is this working?
               ! DOES NOT WORK              if ( h_tot_fftw(i_fftw).NE. h_tot_fftw(i_fftw) ) then
               ! DOES NOT WORK                  h_tot_fftw(i_fftw) = 0.
               ! DOES NOT WORK                  flag_bad_screen = .TRUE.
               ! DOES NOT WORK              end if
            else
               ! normalize height
               h_tot_fftw(i_fftw) = h_tot_fftw(i_fftw)/cajitas_fftw(i_fftw)

       end if

       enddo    
end subroutine chain_fft

subroutine fftw_start(mode)
! this function activate all fftw routines
    integer (kind = 8), intent(in) :: mode
!     if(fin.eq.0) then
!        call dfftw_plan_dft_r2c_1d(plan_fftw, n_box_fftw, h_tot_fftw, Ipar_fftw, FFTW_ESTIMATE)
!       call dfftw_execute(plan_fftw)
!        elseif(fin.eq. (n_screen-1) ) then
!            call dfftw_execute(plan_fftw)
!            call dfftw_destroy_plan(plan_fftw)
!        else
!            call dfftw_execute(plan_fftw)
!     end if
 
    
select case(mode)
  case(0) !init
    call dfftw_plan_dft_r2c_1d(plan_fftw, n_box_fftw, h_tot_fftw, Ipar_fftw, FFTW_ESTIMATE)
  case(1) ! finalize
    call dfftw_destroy_plan(plan_fftw)
  case(2) ! computew
    call dfftw_execute(plan_fftw)
  end select
end subroutine fftw_start

subroutine save_mean_difrac()
!print *, " * Storing difrac_promedio..."
        open(unit=24, file='difrac_promedio')
       
        i_part=0
        q_fftw = 0
        lambda_fftw = 0
        
        
        do i_part = 1, n_box_fftw/2+1
! originally        do while (q_fftw < 1./(delta_x*2) )!i=1,(n_box_fftw/2+1)
! originally            i_part=i_part+1
! originally            q_fftw =real(i_part-1)* (1./boundary(1))

       PIii=4.D0*DATAN(1.D0)

!            q_fftw =real(i_part-1)*((2.*PIii)/2.)* (1./boundary(1))*boundary(1)/(1.*n_mon) 
            q_fftw =real(i_part-1)*((2.*PIii)/2.)* 1./n_mon
            lambda_fftw = (2*PIii)/q_fftw
            ! the true mean intensities, correcting the number of good screen,
            ! (without NaN in h_tot_fftw array)
            ! Imean_fftw(i_part) = Imean_fftw(i_part) / (n_screen - n_bad_screen)
            Imean_fftw(i_part) = Imean_fftw(i_part) / (j_fftw - n_bad_screen)
            
            ! Mean value of Fourier complex coefficients accumulated
            Amp_fftw(i_part) = Amp_fftw(i_part) / (j_fftw - n_bad_screen)
            Im_fftw(i_part) = Im_fftw(i_part) / (j_fftw - n_bad_screen)
            
            write(24,*) q_fftw, Imean_fftw(i_part)
        
                        
        enddo
        
        
        close(unit=26)
        close(unit=24)
!        close(unit=25)

    n = n_mon*n_chain
 
!   do j=1,n
!            print '(/a,3f16.8/)',"chain",r0(:,j)
!            print '(/a,3f16.8/)',"chain unfold",r0_unfold(:,j)
!   end do
 
 !   print *, "  * difrac_promedio was correctly store, n_bad_screen =  ", n_bad_screen
!    print *, "  * n_screen =  ", n_screen 
   
    open(unit=24, file='parametros_fft', position ='append')
        write(24,*) "[fft]n_bad_screen: ", n_bad_screen
    close(24)


end subroutine save_mean_difrac

end subroutine chain_fftw

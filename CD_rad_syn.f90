module CD_rad_syn
  use global
  implicit none
  
  public :: CD_rad_syn_run, CD_rad_syn_data_run, vec_nu
  private :: SSA
  
  contains

  !============================================================================
  ! CD_rad_syn_run (START)
  !============================================================================
  ! This routine calculates the intrinsic (i.e. unabsorbed) synchrotron emission
  ! using as an input the electron distribution in diste(number).dat. The 
  ! Integrated luminosity is stored in rad_syn.dat and emi_syn.dat.
  ! the emission Eps*L_sy(l,Eps_i), which can be corrected for absorption. 
  ! The frequency vector nu(i) is created by the subroutine vec_nu in this module, 
  ! so that it matches for both S1 and S2
  !============================================================================

  subroutine CD_rad_syn_run(number)
    ! ml - CD, mE_e - electrons, mnu - emitted photons and abs_ff table
    character, intent(in) :: number
    logical :: check_SSA = .True.
    integer :: il
    real(dp) :: Eps(mnu),dEps,Eps_int, gamma_e(mE_e), Esy_eV(mnu)
    real(dp) :: nuc,sum,sum_sy
    real(dp) :: ne,nu_RT,x_sy,cte_sy, dummy
    real(dp) :: Lsy_SSA(ml_2,mnu), Lsy_tot_SSA(mnu), sum_sy_ssa
    call validate_shock_number(number,'CD_rad_syn_run')
    800  format(A11,X,Es10.2)
    allocate( N_e(mE_e), N_etot(ml,mE_e) )
    allocate( Lsy_tot(mnu), Lsy(ml_2,mnu) )
    Lsy_tot = 0.d0

    print *, 'Calculate (unabsorbed) synchrotron emission in S',number

    !------------------------------------------------------------------
    ! Read positions (x(l),y(l)) and thermodynamical quantities at them
    
    open(11,file=output_file('thermo'//number//'.dat'))
    do l = 1, ml
      read(11,*) il,rho(l),P(l),v(l),T(l),B(l),n_a(l) ! cgs units
    end do
    close(11)
    B_par = B*par

    open(12,file=output_file('H'//number//'.dat'))
    do l = 1, ml_2
      read(12,*) dummy, H_sh(l), dummy, dummy, dVol(l)
    end do
    close(12)

    ! Electron distribution at position l with energy E_e(i): N_etot(l,i) 
    open (14,file=output_file('dist_e'//number//'.dat'))
    do l = 1, ml
      do i = 1, mE_e
        read(14,*) j, E_e(i),N_etot(l,i)
        ! Energy from eV to erg
        E_e(i) = E_e(i)*eV    
      end do
      read(14,*) 
    end do
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    gamma_e = E_e/mec2

    Eps = E_sy/mec2
    Eps_int = eps(2)/eps(1)
    
    !------------------------------------------------------------------------
    ! Begin synchrotron calculations 

    ! At each position calculate the synchrotron emission and then add everything
    sum_sy = 0.d0
    do l = 1, ml_2
      N_e = N_etot(l,:)      ! N_e = electron distribution at position "l"
      ne = n_a(l)
      nu_RT = 20.d0*ne/B(l)  ! characteristic frequency for R-T
      cte_sy = 4.39d-22*B_par(l)/h
      ! The photon energies vary between Eps_min and max(nu_c(min;max))
      do i = 1, mnu           ! Eps(i), nu(i)
        dEps = Eps(i)*(Eps_int-1.d0)
        sum = 0.d0
        ! Lsy = P(epsilon) = int_Emin^Emax N(E) P(E,eps) dE [erg/s/eps],
        do j = 1, mE_e         ! Ne(j), nuc(j), gamma_e(j)
          nuc = 4.22d6*B_par(l)*gamma_e(j)**2
          x_sy = nu(i)/nuc
          if (x_sy <= 10.d0) then
            P_Eps = cte_sy*x_sy**(1.d0/3.d0)*exp(-x_sy)
            dE_e = E_e(j)*(Eint_e - 1.d0)
            sum = sum + N_e(j)*P_eps*dE_e
          end if
        end do    
        Lsy(l,i) = sum*exp(-nu_RT/nu(i))   ! Lsy corrected for R-T effect
        Lsy_tot(i) = Lsy_tot(i) + Lsy(l,i) ! [erg/s/erg]
        sum_sy = sum_sy + Lsy(l,i)*mec2*dEps
      end do  ! Eps(i),nu(i)
    end do     ! l
    
    write(*,800) 'L_syn(unabs)=',sum_sy

    Lsy_tot_SSA = 0.d0
    Lsy_SSA = 0.d0
    if (check_SSA) then
      call SSA(E_e, N_etot, E_sy, Lsy, B, dVol, Lsy_SSA, number)
      sum_sy_ssa = 0.d0
      do l = 1, ml_2
        do i = 1, mnu           
          dEps = Eps(i)*(Eps_int-1.d0)
          Lsy_tot_SSA(i) = Lsy_tot_SSA(i) + Lsy_SSA(l,i) 
          sum_sy_ssa = sum_sy_ssa + Lsy_SSA(l,i)*mec2*dEps
        end do  ! Eps(i),nu(i)
      end do     ! l
      
      write(*,800) 'L_syn_SSA=',sum_sy_ssa
    endif

    !-----------------------------------------------------------

    Esy_eV = E_sy/eV
    open(17,file=output_file('rad_syn'//number//'.dat'))
    do i = 1, mnu
      write(17,*) Esy_eV(i), Lsy_tot(i), Lsy_tot_SSA(i)
    end do
    close(17)

    open(79,file=output_file('emi_syn'//number//'.dat'))
    do l = 1, ml_2
      do i = 1, mnu
        write(79,*) l, Esy_eV(i), Lsy(l,i), Lsy_SSA(l,i)
      end do
      write(79,*) 
    end do
    do l = ml_2+1, ml
      do i = 1, mnu
        write(79,*) l, Esy_eV(i), Lsy(l-ml_2,i), Lsy_SSA(l-ml_2,i)
      end do
      write(79,*) 
    end do
    close(79)

    deallocate( Lsy_tot, Lsy )
    deallocate( N_e, N_etot )

  END subroutine CD_rad_syn_run

  !============================================================================
  ! CD_rad_syn_data_run (START)
  !============================================================================
  ! Same as CD_rad_syn_run but only for specific frequencies
  !============================================================================

  subroutine CD_rad_syn_data_run(number, nu_data)
    ! ml - CD, mE_e - electrons, mnu - emitted photons and abs_ff table
    character, intent(in) :: number
    real(dp), intent(in) :: nu_data(mnu_data)
    integer :: il
    real(dp) :: gamma_e(mE_e), Esy_eV(mnu_data)
    real(dp) :: nuc,sum
    real(dp) :: ne,nu_RT,x_sy,cte_sy
    call validate_shock_number(number,'CD_rad_syn_data_run')
    allocate( N_e(mE_e), N_etot(ml,mE_e) )
    allocate( Lsy(ml,mnu_data) )

    if (verbose) print *, 'Calculate (unabsorbed) synchrotron emission in S',number

    !------------------------------------------------------------------
    ! Read positions (x(l),y(l)) and thermodynamical quantities at them
    
    open(11,file=output_file('thermo'//number//'.dat'))
    do l = 1, ml
      read(11,*) il,rho(l),P(l),v(l),T(l),B(l),n_a(l) ! cgs units
    end do
    close(11)
    B_par = B*par

    ! Electron distribution at position l with energy E_e(i): N_etot(l,i) 
    open (14,file=output_file('dist_e'//number//'.dat'))
    do l = 1, ml
      do i = 1, mE_e
        read(14,*) j, E_e(i),N_etot(l,i)
        ! Energy from eV to erg
        E_e(i) = E_e(i)*eV    
      end do
      read(14,*) 
    end do
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    gamma_e = E_e/mec2

    !------------------------------------------------------------------------
    ! Begin synchrotron calculations 

    ! At each position calculate the synchrotron emission and then add everything
    do l = 1, ml_2
      N_e = N_etot(l,:)      ! N_e = electron distribution at position "l"
      ne = n_a(l)
      nu_RT = 20.d0*ne/B(l)  ! characteristic frequency for R-T
      cte_sy = 4.39d-22*B_par(l)/h
      ! The photon energies vary between Eps_min and max(nu_c(min;max))
      do i = 1, mnu_data           ! Eps(i), nu(i)
        sum = 0.d0
        ! Lsy = P(epsilon) = int_Emin^Emax N(E) P(E,eps) dE [erg/s/eps],
        do j = 1, mE_e         ! Ne(j), nuc(j), gamma_e(j)
          nuc = 4.22d6*B_par(l)*gamma_e(j)**2
          x_sy = nu_data(i)/nuc
          if (x_sy <= 10.d0) then
            P_Eps = cte_sy*x_sy**(1.d0/3.d0)*exp(-x_sy)
            dE_e = E_e(j)*(Eint_e - 1.d0)
            sum = sum + N_e(j)*P_eps*dE_e
          end if
        end do    
        Lsy(l,i) = sum*exp(-nu_RT/nu_data(i))   ! Lsy corrected for R-T effect
      end do  ! Eps(i),nu(i)
    end do     ! l
    
    !-----------------------------------------------------------
    Esy_eV = nu_data * (h/eV)

    open(79,file=output_file('emi_syn_data'//number//'.dat'))
    do l = 1, ml_2
      do i = 1, mnu_data
        write(79,*) l, Esy_eV(i), Lsy(l,i)
      end do
      write(79,*) 
    end do
    do l = ml_2+1, ml
      do i = 1, mnu_data
        write(79,*) l, Esy_eV(i), Lsy(l-ml_2,i)
      end do
      write(79,*) 
    end do
    close(79)

    deallocate( Lsy )
    deallocate( N_e, N_etot )

  END subroutine CD_rad_syn_data_run

  !============================================================================
  ! vec_nu
  !============================================================================
  
  subroutine vec_nu()!E_nu)
    ! Build the nu(mnu) and E_sy(mnu) vectors
    !real(dp), intent(out) :: E_nu(mnu)
    real(dp) :: dummy, B_aux, B_min, B_max
    real(dp) :: Eaux, Ee_max, gamma_e_min, gamma_e_max
    real(dp) :: nuc_max, nuc_min, nu_int, nu_min_Hz, nu_max_Hz

    !---------------------------------------------------------
    ! Obtain B_min and B_max
    B_min = 1.d10 ! Dummy (big) value
    B_max = 0.d0
    open(11,file=output_file('thermo1.dat'))
    open(12,file=output_file('thermo2.dat'))
    do l = 1, ml_2
      read(11,*) j,dummy,dummy,dummy,dummy,B_aux,dummy
      B_max = max(B_aux,B_max)
      B_min = min(B_aux,B_min)
      read(12,*) j,dummy,dummy,dummy,dummy,B_aux,dummy
      B_max = max(B_aux,B_max)
      B_min = min(B_aux,B_min)
    end do
    close(11)
    close(12)
    B_max = B_max*par 
    B_min = B_min*par 

    !--------------------------------------------
    ! Obtain gamma_e_min and gamma_e_max
    Ee_max = 0.d0 !(this is the cutoff energy, not the max energy)
    open (14,file=output_file('Emax_e1.dat'))
    open (15,file=output_file('Emax_e2.dat'))
    do l=1,ml_2
      read(14,*) dummy, Eaux
      Ee_max = max(Eaux,Ee_max)
      read(15,*) dummy, Eaux
      Ee_max = max(Eaux,Ee_max)
    end do
    close(14) 
    close(15) 
    gamma_e_max = 3.d0 * Ee_max * eV/mec2
    gamma_e_min = Emin_e /mec2

    !--------------------------------------------
    ! Obtain nuc_min and nuc_max, including the hard-coded values in globals
    nuc_min = 4.22d6*B_min*gamma_e_min**2
    nuc_max = 4.22d6*B_max*gamma_e_max**2
    
    nu_min_Hz = max(nuc_min/5.d0, nu_min_syn_Hz)   
    nu_max_Hz = min(10.d0*nuc_max, nu_max_syn_Hz) 

    ! Build nu vector in [Hz]
    call vector_log(nu_min_Hz, nu_max_Hz, nu_int, nu)
    E_sy = h*nu

  END subroutine vec_nu

  !============================================================================
  ! SSA
  !============================================================================
  ! Luminosity corrected for synchrotron self-absorption
  !============================================================================
  subroutine SSA(Ee, Ne_tot, E_sy, L_sy, B, dVol, L_sy_SSA, number)
  real(dp), intent(in) :: Ee(mE_e), Ne_tot(ml,mE_e), E_sy(mnu), L_sy(ml_2, mnu), B(ml), dVol(ml)
  character, intent(in) :: number
  real(dp), intent(out) :: L_sy_SSA(ml_2, mnu)
  real(dp) :: dV(ml), tau_SSA, depth, K_SSA, attenuation(ml,mnu)
  real(dp) :: N_l(mE_e)

  ! Note that dVol is the volume of the whole "ring"
  dV = dVol / (2.d0*mphi)

  attenuation = 1.d0
  !Loop for the segments
  do l = 1, ml_2
    N_l = Ne_tot(l,:)

    !Characteristic length that the photon travel
    depth = dVol(l)**(1./3.)

    !Absorption coefficient (Pacholczyck, 1970)
    do i = 1, mnu
      K_SSA = 2.62d0 * qe**3 * h**2 * B(l)/(4d0*pi*E_sy(i)**2*me*dV(l)) * EIntSSA( B(l), E_sy(i), Ee, N_l )

      !Optical depth for SSA
      tau_SSA = K_SSA * depth

      !Luminosity corrected for synchrotron self-absorption
      if (tau_SSA > 1d-2) attenuation(l,i) = ( 1d0 - exp(-tau_SSA) ) / tau_SSA
      L_sy_SSA(l,i) = L_sy(l,i) * attenuation(l,i)
    enddo

  enddo

  open (88,file=output_file('opac_SSA'//number//'.dat'))
  do l = 1, ml
    do i = 1, mnu
      write(88,*) l, nu(i), attenuation(l,i)
    enddo
  end do
  close(88)

  end subroutine SSA

  !============================================================================
  !Function to calculate the synchrotron-self absorption opacity.
  real(dp) function F_SSA(E, NE, B, Eph)
  real(dp), intent(in) :: E, NE, B, Eph
  real(dp) :: Ec, xc, cte_Ec
  cte_Ec = sqrt(3./2.) * (qe * h) / (2 * pi * mec2**3 / c)
  Ec = cte_Ec * B * E**2        !Critical photon energy and function to integrate for each particle energy
  xc = Eph/Ec
  if (xc > 1d-3 .and. xc < 1d1) then ! Approximation range implemented here.
    F_SSA = NE/E * exp(-xc) * xc**(1./3.) * (xc + 2./3.)
  else 
    F_SSA = 0.0 
  endif 
  end function F_SSA

  !============================================================================
  ! Integral for the particle energy in SSA using the trapezoid method
  real(dp) function EintSSA(B, Eph, E_l, N_l)
  real(dp), intent(in) :: B, Eph, E_l(mE_e), N_l(mE_e)
  real(dp) :: y1, y2, y_int, dE_factor, dE_l
  integer :: i9
  dE_factor = E_l(2)/E_l(1) - 1.d0
  EintSSA = 0.d0
  y1 = F_SSA(E_l(1), N_l(1), B, Eph)
  do i9 = 1, mE_e-1
    y2 = F_SSA(E_l(i9+1), N_l(i9+1), B, Eph)
    y_int = 0.5 * (y1 + y2)
    dE_l = E_l(i9) * dE_factor
    EintSSA = EintSSA + y_int*dE_l
    y1 = y2
  end do
  end function EintSSA

!================================================================

END module CD_rad_syn

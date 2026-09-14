module CD_dist_e
  use global
  implicit none
  
  public :: CD_dist_e_run
  
contains

  !============================================================================
  ! CD_dist_e_run (START)
  !============================================================================
  ! This program calculates, using CD.dat and thermo(SHOCK).dat, the non-thermal
  ! electron energy distribution at each position of a curve at either side of 
  ! the CD. The output is Ne(SHOCK).dat, which stores Ne(l,E(i)) of size ml_2*mE_e.
  ! The maximum and minimum electron energies are set in global.f90
  !============================================================================

  subroutine CD_dist_e_run(number,D)

    character, intent(in) :: number
    real(dp), intent(in) :: D
    integer :: l_min,l2
    real(dp) :: t_cell(ml_2),t_loss(ml_2,mE_e),t_adi(ml_2)
    real(dp) :: D_diff,t_diff,integral,t_coole(ml_2,mE_e),t_esc
    real(dp) :: t_ic1,t_ic2,t_ic,t_sy,t_ac,t_adv(ml_2),t_br,t_br_adi
    real(dp) :: gamm,E_corr,C_ic1_k,C_ic1,C_ic2_k,C_ic2,cte_sy
    real(dp) :: log10Emaxl,R_acel(ml_2),R1,R2,Q
    real(dp) :: E_dot(ml_2,mE_e),Emax_e0(ml_2),Ka,E_NT,sum,U_NT,U_B
    real(dp) :: N_e0(mE_e),n_a(ml_2),L_inj(ml_2),Dl(ml_2)
    real(dp) :: dE_approx,cool,dummy, f1, f2, f
    integer :: index1,index2,count,steps
    logical :: new_count

    interface
      subroutine bisect_e(Be,logEmax,c1,c2,tbr,tad,xi,Rdiff)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp), intent(in) :: Be,c1,c2,tbr,tad,xi,Rdiff
        real(dp), intent(out) :: logEmax
      end subroutine bisect_e
    end interface
    
    interface
      function fQ(E,Estiff,Emin,Emax,Ecorr) ! function Q(E)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp) :: fQ
        real(dp), intent(in) :: E, Estiff, Emin, Emax, Ecorr
      end function fQ
    end interface

    call validate_shock_number(number,'CD_dist_e_run')

    allocate( N_e(mE_e),N_etot(ml_2,mE_e) )
    
    if (verbose) print *, 'Calculate the non-thermal electron distribution in S',number

    !--------------------------------------------------------------------------
    ! Read thermodynamical quantities at each shock position (x(l),y(l))
    open(11,file=output_file('thermo'//number//'.dat'))
    open(12,file=output_file('Linj'//number//'.dat'))
    open(13,file=output_file('eta_acc'//number//'.dat'))
    open(14,file=output_file('H'//number//'.dat'))
    do l = 1, ml_2
      read(11,*) l2,rho(l),P(l),v(l),T(l),B(l),n_a(l) ! CGS
      read(12,*) x(l),y(l),L_inj(l),Dl(l) ! the upper side of the CD
      read(13,*) dummy,eta_acc(l)
      read(14,*) dummy, H_sh(l), dummy, dummy, dVol(l)
    end do
    close(11)
    close(12)
    close(13)
    close(14)
    B_par = B*par        ! only the parallel component of B is involved in sync.
    L_inj = inj_e*L_inj  ! L_inj(l) has the injected energy at position l

    if (number=='1') then 
      z_br = z_br1
    else
      z_br = z_br2
    end if
    
    !--------------------------------------------------------------------
    ! Build an energy vector that is common for all linear-emitters

    call vector_log(Emin_e,Emax_e,Eint_e,E_e)
    

    !========================================================================
    ! First calculate the cooling time at all energies and at each position 
    ! (and also the maximum energy)

    C_ic1 = ( pi*h_bar**3*c**2/(2.d0*re**2*mec2**3) )/(T1mcc**2)
    C_ic2 = ( pi*h_bar**3*c**2/(2.d0*re**2*mec2**3) )/(T2mcc**2)

    open(70,file=output_file('CD_t_cool_e'//number//'.dat'))
    open(78,file=output_file('Emax_e'//number//'.dat'))
    do l = 1, ml_2
      R1 = sqrt(x(l)**2+y(l)**2)     ! Distance to star 1
      R2 = sqrt((D-x(l))**2+y(l)**2) ! Distance to star 2
      R_acel(l) = min(R1,R2)    ! Accelerator size ~ distance to closest star 
       
      t_cell(l) = Dl(l)/v(l)    ! Cell advection time
      t_adv(l) = R_acel(l)/v(l) ! Characteristic advection time
    end do

    do l = 1, ml_2-1
      t_adi(l) = 3.*Dl(l)/(v(l)*log(rho(l)/rho(l+1))) ! Characteristic adiabatic losses time
      if (t_adi(l) <= 0.d0 .or. t_adi(l) > 3.d16) t_adi(l) = 3.d16   ! To avoid numerical errors set to 1 Gyr
    end do
    t_adi(ml_2) = t_adi(ml_2-1)

    do l = 1,ml_2
      R1 = sqrt(x(l)**2+y(l)**2)     ! Distance to star 1
      R2 = sqrt((D-x(l))**2+y(l)**2) ! Distance to star 2
      Kappa1 = (Rst1/(2.*R1))**2
      Kappa2 = (Rst2/(2.*R2))**2
      C_ic1_k = C_ic1/Kappa1
      C_ic2_k = C_ic2/Kappa2  
      cte_sy = 6.*pi*mec2**2/(c*sigma_T*B(l)**2)   ! Valid for isotropic e- distribution

      t_br = 1.d15/n_a(l)/z_br       ! Energy independent
      t_br_adi = 1.d0/(1.d0/t_br + 1.d0/t_adi(l))

      call bisect_e(B(l),log10Emaxl,C_ic1_k,C_ic2_k,t_br_adi,t_adv(l),eta_acc(l),R_acel(l))
      Emax_e0(l) = MIN( 10.d0**log10Emaxl , qe*B(l)*R_acel(l)*sqrt(1.5d0/eta_acc(l)) ) ! Hillas
      if (Emax_e0(l) >= 0.5*Emax_e) then 
        print *, 'Emax_e in global.f90 is too small!'
        stop
      endif

      do i = 1,mE_e
        D_diff = E_e(i)*c/(3.d0*qe*B(l))
        t_diff = R_acel(l)**2/(2.d0*D_diff) ! Diffusion losses
        t_esc = 1./( 1./t_adv(l) + 1.d0/t_diff ) ! escape
          
        t_sy = cte_sy/E_e(i) 
        gamm = E_e(i)/mec2
        t_ic1 = C_ic1_k*gamm/F_iso(4.*gamm*T1mcc)
        t_ic2 = C_ic2_k*gamm/F_iso(4.*gamm*T2mcc)
        t_ic = 1./(1.d0/t_ic1 + 1.d0/t_ic2)
        t_coole(l,i) = 1./(1.d0/t_sy + 1.d0/t_ic + 1.d0/t_br + 1./t_adi(l)) 
        
        t_loss(l,i) = min(t_coole(l,i),t_esc) ! general losses
        
        t_ac = eta_acc(l)*E_e(i)/(B(l)*c*qe)
          
        write(70,*) l,E_e(i)/eV,t_ic1,t_ic2,t_ic,t_sy,t_br,t_diff,t_ac,t_coole(l,i),t_adv(l),t_cell(l),t_esc,t_loss(l,i),t_adi(l)

        E_dot(l,i) = E_e(i)/t_coole(l,i)
      end do

      write(78,*) l, Emax_e0(l)/eV
    end do
    close(70)
    close(78)

    !========================================================================
    
    N_etot = 0.d0
    E_corr = Estiff_e**(alpha2-alpha)
    ! At the starting point of the line we calculate the distribution 
    ! Ne_0(E_0), the cooling rate Edot0, and we sum it to Ne_tot(l,i)

    Do l_min=1,ml_2

      ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. The normalization constant Ka is obtained 
      ! by the condition: integral Q(E')*E'dE' = Linj. [Linj]=erg/s, [Ka]=erg/s, [Q]=1/erg/seg
       
      ! The integral is solved with the trapezoid method; Q(E) is obtained from fQ
      integral = 0.d0
      i = 1
      f1 = E_e(1) * fQ(E_e(1), Estiff_e, Emin_e_inj, Emax_e0(l_min), E_corr)
      do while (i < mE_e .and. E_e(i+1) < 10.d0*Emax_e0(l_min))
        f2 = E_e(i+1) * fQ(E_e(i+1), Estiff_e, Emin_e_inj, Emax_e0(l_min), E_corr)
        f = 0.5 * (f1 + f2)
        dE_e = E_e(i)*(Eint_e-1.d0)
        integral = integral + f*dE_e
        f1 = f2
        i = i + 1
      end do
      Ka = L_inj(l_min)/integral ! Injected at the starting point of the line

      !------------------------------------------------------------------------
       
      N_e0 = 0.d0            ! Initialize to zero at each loop in l_min

      i = 1
      do while (i < mE_e+1 .and. E_e(i) < 10.d0*Emax_e0(l_min)) 
        ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' = Q(E)*dt, with |dE/dt|= E/t_cool 
        ! in [erg/seg], dt = min(t_cell,t_cool), and [N(E)]=1/erg
        Q = Ka * fQ(E_e(i), Estiff_e, Emin_e_inj, Emax_e0(l_min), E_corr)
        N_e0(i) = Q*min(t_cell(l_min),t_coole(l_min,i))
        N_etot(l_min,i) = N_etot(l_min,i)+N_e0(i)  

        ! We calculate the evolved version of the injected electron distribution at each
        ! position "l" as: N_l(l) = N_l-1* (Edot_l-1/Edot_l), where the energy also evolves
        ! from E0 to El given by: t_cell(l) = int_E0^El dE/Edot
          
        index1 = i            ! For the cell "l"
        index2 = i            ! For the previous cell, "l-1"
        N_e = 0.d0            ! Initialize to zero at each loop in l_min
        N_e(index1) = N_e0(i) ! This will be evolved
        l = l_min + 1  
        count = 0
        cool = 0.d0
        new_count = .true.      

        do while( l <= ml_2 .and. index1 > 1) ! Stop at the lowest energy bin considered.
          ! At the position l, the particle has E_e(index1) instead of E_e(i)
          ! We find index1 <= i according to the cooling time:

          if (l == l_min+1) then ! It will cool from l_min to ml_2
            do l2=l,ml_2
              cool = cool + t_cell(l2)/t_coole(l2,index1)
            end do
          else ! cool(l -> ml_2) = cool(1 -> ml_2) - cool(1 -> l-1)
            cool = cool - t_cell(l-1)/t_coole(l-1,index1)
          end if
          dE_approx = (1.-t_cell(l)/t_coole(l,index1))

          if ( 1.-cool > 0.99/Eint_e ) then ! negligible cooling
            continue 

          elseif ( 1.-cool <= 0.99/Eint_e .and. dE_approx > 1./Eint_e .and. index1 > 1 ) then
            if (new_count) then
              sum = 1.d0
              steps = 0 ! number of steps until it changes E(index1) -> E(index-1)
              l2 = l
              do while (sum > 1./Eint_e .and. l2 <= ml_2)
                steps = steps + 1 
                sum = sum*(1.-t_cell(l2)/t_coole(l2,index1))
                l2 = l2 + 1
                if (sum > 1./Eint_e .and. l2 == ml_2) steps = ml_2 ! NOT COOLING
              end do
              new_count = .false.
            end if
            count = count + 1

            if (count >= steps) then
              count = 0
              index1 = index1 - 1
              new_count = .true. 
            end if

          elseif ( dE_approx <= 1./Eint_e ) then
            sum = 0.d0
            do while (sum < t_cell(l) .and. index1 > 1)
              index1 = index1 - 1    
              dE_e = E_e(index1)*(Eint_e-1.d0)
              sum = sum + dE_e/E_dot(l,index1)        
            end do

          end if

          ! The following takes into account energy losses and variable cell size and velocity
          N_e(index1) = N_e(index2) * ( E_dot(l,index2)/E_dot(l,index1) ) * (t_cell(l)/t_cell(l-1))

          N_etot(l,index1) = N_etot(l,index1) + N_e(index1)  ! The energies match by construction

          index2 = index1 ! index2 keeps the value of the previous step

          l = l + 1 
        end do               ! End loop in shock position (l) 
        i = i + 1
      end do                  ! End loop in electron starting energies (i)

    end do                     ! End loop in starting position (l_min = # linear-emitters)

    ! Save the electron distribution at each position 
        open(79,file=output_file('dist_e'//number//'.dat'))
    do l = 1,ml_2
      do i = 1,mE_e
        write(79,*) l,E_e(i)/eV,N_etot(l,i)
      end do
      write(79,*) 
    end do
    do l = ml_2+1,ml
      do i = 1,mE_e
        write(79,*) l,E_e(i)/eV,N_etot(l-ml_2,i)
      end do
      write(79,*) 
    end do
    close(79)

    open(80,file=output_file('ENTe'//number//'.dat'))
    open(81,file=output_file('PNTe'//number//'.dat'))
    do l = 1,ml_2
      E_NT = 0.d0
      U_B = B(l)**2/(8.d0*pi)
      do i = 1,mE_e       
        dE_e = E_e(i)*(Eint_e-1.d0)
        E_NT = E_NT + E_e(i)*N_etot(l,i)*dE_e
      end do
      E_NT = E_NT * 2.d0*mphi ! Total NT energy in WCR "ring"
      U_NT = E_NT/dVol(l)
      write(80,*) l, E_NT, U_NT, U_NT/U_B ! Store the energy density as well
      write(81,*) l, (U_NT/3.d0)/U_B      ! Ratio P_e,rel/P_mag
    end do
    close(80)
    close(81)   

    !------------------------------------------------------------
    ! Group in bins according to distance
    if (calc_N_bin) call Ne_bin(number,E_e,N_etot,Dl,D)

    !------------------------------------------------------------

    deallocate( N_e,N_etot ) 

  END subroutine CD_dist_e_run

  !============================================================================
  ! CD_dist_e_run (END)
  !============================================================================
  
END module CD_dist_e


  !--------------------------------------------------------------------------!
  !                                                                          !
  !                            SUBROUTINES                                   !
  !                                                                          !
  !--------------------------------------------------------------------------!

  subroutine Ne_bin(number,Ee,Ne_tot,dL,D_orb)
    use global
    implicit none    
    character, intent(in) :: number
    real(dp), intent(in) :: Ee(mE_e),Ne_tot(ml_2,mE_e),dL(ml_2),D_orb
    real(dp) :: length, N_total_i
    integer, parameter :: N_bins = 5 ! For the plot diste_sw_bins
    real(dp) :: N_bin_i(N_bins,mE_e)

    call validate_shock_number(number,'Ne_bin')
    
    N_bin_i = 0.d0
    open(81,file=output_file('dist_e'//number//'_total.dat'))
    do i=1,mE_e
      N_total_i = 0.d0
      length = 0.d0
      do l=1,ml_2
        ! Make some distributions binned in intervals of 0.1*D
        length = length + dL(l)
        if (length < 0.1*D_orb) then
          N_bin_i(1,i) = N_bin_i(1,i) + Ne_tot(l,i)
          elseif (length >= 0.1*D_orb .and. length < 0.2*D_orb) then
          N_bin_i(2,i) = N_bin_i(2,i) + Ne_tot(l,i)
          elseif (length >= 0.2*D_orb .and. length < 0.4*D_orb) then
          N_bin_i(3,i) = N_bin_i(3,i) + Ne_tot(l,i)
          elseif (length >= 0.4*D_orb .and. length < 0.8*D_orb) then
          N_bin_i(4,i) = N_bin_i(4,i) + Ne_tot(l,i)
          elseif (length >= 0.8*D_orb .and. length < 1.6*D_orb) then
          N_bin_i(5,i) = N_bin_i(5,i) + Ne_tot(l,i)
        end if  
        N_total_i = N_total_i + Ne_tot(l,i)
      end do

      write(81,*) Ee(i)/eV, N_total_i* 2.d0*mphi ! For the whole WCR   
    end do
    close(81)

    open(82,file=output_file('dist_e_wcr'//number//'_bin.dat'))
    do j = 1,N_bins
      do i = 1, mE_e
        write(82,*) j,Ee(i)/eV, N_bin_i(j,i)* 2.d0*mphi ! For the whole WCR   
      end do
      write(82,*) 
    end do
    close(82)

  end subroutine Ne_bin


  !=========================================================================

   Subroutine bisect_e(Be,x3,c1,c2,t_br,t_adv,xi,Rdiff)
    use global
    implicit none
    integer t_err
    real(dp) :: a_,b_,c_,err
    real(dp), intent(in) :: Be,c1,c2,t_br,t_adv,xi,Rdiff
    real(dp), intent(out) :: x3
    interface
      function diff_e(a,B,ic1,ic2,t_br,t_adv,efi,Rd)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp) :: diff_e
        real(dp), intent(in) :: a,B,ic1,ic2,t_br,t_adv,efi,Rd
      end function diff_e
    end interface

    ! The error in the n-th step is < |a-b|/2**n
    a_ = log10(1.d8*eV) ! log10(Emin [CGS]) ~ 1 GeV
    b_ = log10(Emax_e)  ! log10(Emax [CGS])
    t_err = 2
    err = 10.d0**(-t_err) ! error < 10**-t
    i = 0
    do while ( abs(b_ - a_) > err )
      i = i+1
      c_ = (a_+b_)/2.d0
      if (abs(diff_e(c_,Be,c1,c2,t_br,t_adv,xi,Rdiff)) <= 1.d-6) then
        b_ = c_
        else if (diff_e(a_,Be,c1,c2,t_br,t_adv,xi,Rdiff)*diff_e(c_,Be,c1,c2,t_br,t_adv,xi,Rdiff) < 0.d0) then
        b_ = c_
        else if (diff_e(a_,Be,c1,c2,t_br,t_adv,xi,Rdiff)*diff_e(c_,Be,c1,c2,t_br,t_adv,xi,Rdiff) > 0.d0) then
        a_ = c_
      end if
    end do
    x3=b_

   end Subroutine bisect_e

  !=========================================================================


   Function diff_e(log10E,B_l,c1,c2,t_br,t_adv,xi,Rdiff)
    use global
    implicit none
    real(dp), intent(in) :: log10E,B_l,c1,c2,t_br,t_adv,xi,Rdiff
    real(dp) :: diff_e,E_l,gamm,t_sy,t_ic1,t_ic2,t_ac,t_coole
    real(dp) :: D_diff,t_diff,t_loss,t_esc

    E_l = 10.d0**log10E

    D_diff = E_l*c/(3.d0*qe*B_l) ! Bohm
    t_diff = Rdiff**2/(2.d0*D_diff) ! Diffusion losses APPROX
    t_esc = 1.d0/(1.d0/t_adv + 1.d0/t_diff) ! escape
    
    gamm = E_l/mec2
    t_sy = 1.d0/(1.6d-3 * (B_l*par)**2 * E_l)
    t_ic1 = c1*gamm/F_iso(4.d0*gamm*T1mcc)
    t_ic2 = c2*gamm/F_iso(4.d0*gamm*T2mcc)
    
    t_coole = 1.d0/(1.d0/t_sy + 1.d0/t_ic1 + 1.d0/t_ic2 + 1.d0/t_br) 
    
    t_loss = min(t_coole,t_esc) ! general losses
    
    t_ac = xi*E_l/(B_l*c*qe)  !  t_ac = (efiac*B_l*c*qe/E_l)**(-1)
    diff_e = t_loss - t_ac
    
   end Function diff_e

  !=========================================================================

  Function fQ(E, Estiff, Emin, Emax, Ecorr)
    use global  
    real(dp) :: fQ
    real(dp),intent(in) :: E, Estiff, Emin, Emax, Ecorr
    if (E < Emin) then     ! No injection for E < E_min
      fQ = 0.d0 
    else if (E >= Emin .and. E <= Estiff) then
      fQ = E**(-alpha)*exp(-E/Emax)
    else                   ! Hardening at E > Estiff
      fQ = E**(-alpha2)*exp(-E/Emax)*Ecorr
    endif

  end Function fQ

  !=====================================================================

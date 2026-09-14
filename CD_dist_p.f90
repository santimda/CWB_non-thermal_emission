module CD_dist_p
  use global
  implicit none
  
  public :: CD_dist_p_run

contains

  !============================================================================
  ! CD_dist_p_run (START)
  !============================================================================
  ! This program calculates, using CD.dat and thermo(SHOCK).dat, the non-thermal
  ! proton energy distribution at each position of a curve at either side of 
  ! the CD. The output is Np(SHOCK).dat, which stores Np(l,E(i)) of size ml_2*mE_p.
  ! The maximum and minimum proton energies are set in global.f90
  !============================================================================
  
  subroutine CD_dist_p_run(number,D)

    character, intent(in) :: number
    real(dp), intent(in) :: D
    integer :: l_min, l2
    real(dp) :: t_cell(ml_2), t_adv(ml_2), t_coolp(ml_2,mE_p), t_adi(ml_2)
    real(dp) :: t_ac, t_pp, t_diff, t_esc, t_loss(ml_2,mE_p), t_adi_pp
    real(dp) :: b_pp, L1, D_diff, log10Emaxl, R_acel(ml_2), R1, R2
    real(dp) :: E_corr, integral, Q
    real(dp) :: N_p0(mE_p), n_a(ml_2), L_inj(ml_2), Dl(ml_2) 
    real(dp) :: E_dotp(ml_2,mE_p), Emax_p0(ml_2), Ka, E_NT, U_NT, U_B
    real(dp) :: L_inj_p(ml_2), sum, dummy, f1, f2, f
    real(dp) :: dE_approx, cool
    integer :: index1, index2, count, steps 
    logical :: new_count

    interface
      subroutine bisect_p(Be,logEmax,tpp,tadv,xi,Racel)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp), intent(in) :: Be, tpp, tadv, xi, Racel
        real(dp), intent(out) :: logEmax
      end subroutine bisect_p
    end interface

    interface
      function fQp(E,Estiff,Emin,Emax,Ecorr) ! function Q(E)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp) :: fQp
        real(dp), intent(in) :: E, Estiff, Emin, Emax, Ecorr
      end function fQp
    end interface

    call validate_shock_number(number,'CD_dist_p_run')

    allocate( N_p(mE_p),N_ptot(ml_2,mE_p) )

    if (verbose) print *, 'Calculate the non-thermal proton distribution in S',number

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
    L_inj_p = inj_p*L_inj ! L_inj(l) is the total injected power at the position l

    if (number=='1') then 
      z_pp = z_pp1
     else
      z_pp = z_pp2
     end if
    
    !--------------------------------------------------------------------
    ! Build an energy vector that is common for all linear-emitters

    call vector_log(Emin_p,Emax_p,Eint_p,E_p)
    
    !--------------------------------------------------------------------
    ! First we calculate the cooling time at all energies and at each position 
    ! (and also the maximum energy, Emax(l))

    open(70,file=output_file('CD_t_cool_p'//number//'.dat'))
    open(78,file=output_file('Emax_p'//number//'.dat'))
    do l = 1, ml_2
      R1 = sqrt(x(l)**2 + y(l)**2)     ! Distance to star 1
      R2 = sqrt((D-x(l))**2 + y(l)**2) ! Distance to star 2
      R_acel(l) = min(R1,R2)    ! Accelerator size ~ distance to closest star 

      t_cell(l) = Dl(l)/v(l)    ! Cell advection time
      t_adv(l) = R_acel(l)/v(l) ! Characteristic advection time
    end do

    do l = 1, ml_2-1
      t_adi(l) = 3.*Dl(l)/(v(l)*log(rho(l)/rho(l+1))) ! Characteristic adiabatic losses time
      if (t_adi(l) <= 0.d0 .or. t_adi(l) > 3.d16) t_adi(l) = 3.d16   ! To avoid numerical errors set to 1 Gyr
    end do
    t_adi(ml_2) = t_adi(ml_2-1)

    do l = 1, ml_2
      t_pp = 1.d15/n_a(l)/z_pp      ! p-p losses (very approx, just to calculate E_max)
      t_adi_pp = 1.d0/(1.d0/t_adi(l) + 1.d0/t_pp)   
       
      call bisect_p(B(l),log10Emaxl,t_adi_pp,t_adv(l),eta_acc(l),R_acel(l))
      Emax_p0(l) = MIN( 10.d0**log10Emaxl, qe*B(l)*R_acel(l)*sqrt(1.5d0/eta_acc(l)) ) ! Hillas
      if (Emax_p0(l) >= 0.5*Emax_p) then
        print *, 'Emax_p in global.f90 is too small!'
        stop
      endif
       
      do i = 1, mE_p
        if (E_p(i) <= E_th) then   ! E_th= 1.22 GeV
          b_pp = 1.d-30            ! just not to divide by zero
        else            
          L1 = log(E_p(i)/TeV)     ! L = ln(Ep/1TeV)
          b_pp = c*k_pp*n_a(l)*sigma_pp(L1,E_p(i))*E_p(i)
        endif
        t_pp = E_p(i)/(b_pp*z_pp)

        D_diff = E_p(i)*c/(3.d0*qe*B(l))
        t_diff = R_acel(l)**2/(2.d0*D_diff) ! Diffusion losses
        t_esc = 1.d0/( 1.d0/t_adv(l) + 1.d0/t_diff ) ! escape

        t_coolp(l,i) = 1./(1.d0/t_pp + 1./t_adi(l))
        
        t_loss(l,i) = 1./(1./t_coolp(l,i) + 1./t_esc) ! general losses
        
        t_ac = eta_acc(l)*E_p(i)/(B(l)*c*qe)
        
        write(70,*) l,E_p(i)/eV,t_ac,t_pp,t_adi(l),t_adv(l),t_diff,t_coolp(l,i),t_cell(l),t_esc,t_loss(l,i)

        E_dotp(l,i) = E_p(i)/t_coolp(l,i)
      end do

      write(78,*) l, Emax_p0(l)/eV
    end do
    close(70)
    close(78)

    !========================================================================
    ! Consider different starting points for each linear emitter

    N_ptot = 0.d0 ! Array containing N_p(l,i) with all the contributions at a cell "l"
    E_corr = Estiff_p**(alpha2-alpha) ! For a break in the injected spectrum

    ! At the starting point of the line we calculate the distribution 
    ! Np_0(E_0), the cooling rate Edot0, and we sum it to Np_tot(l,i)

    Do l_min = 1, ml_2

      ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. The normalization constant Ka is obtained 
      ! by the condition: int Q(E')*E'dE' = Linj. [Linj]=erg/s, [Ka]=erg/s, [Q]=1/erg/seg

      ! The integral is solved with the trapezoid method; Q(E) is obtained from fQp
      integral = 0.d0
      i = 1
      f1 = E_p(1) * fQp(E_p(1), Estiff_p, Emin_p_inj, Emax_p0(l_min), E_corr)
      do while (i < mE_p .and. E_p(i+1) < 10.d0*Emax_p0(l_min))
        f2 = E_p(i+1) * fQp(E_p(i+1), Estiff_p, Emin_p_inj, Emax_p0(l_min), E_corr)
        f = (f1 + f2)/2.d0
        dE_p = E_p(i)*(Eint_p-1.d0)
        integral = integral + f*dE_p
        f1 = f2
        i = i + 1
      end do
      Ka = L_inj_p(l_min)/integral ! Injected at the starting point of the line
      
      !------------------------------------------------------------------------

      N_p0 = 0.d0        ! Initialize to zero at each loop in l_min

      i = 1
      do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emax_p0(l_min)) 
        ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' = Q(E)*dt, with |dE/dt|= E/t_cool 
        ! in [erg/seg], dt = min(t_cell,t_cool), and [N(E)]=1/erg
        Q = Ka * fQp(E_p(i), Estiff_p, Emin_p_inj, Emax_p0(l_min), E_corr)
        N_p0(i) = Q * min(t_cell(l_min),t_coolp(l_min,i))
        N_ptot(l_min,i) = N_ptot(l_min,i) + N_p0(i)

        ! We calculate the evolved version of the injected proton distribution at each
        ! position "l" as: N_l(l) = N_l-1* (Edot_l-1/Edot_l), where the energy also evolves
        ! from E0 to El given by: t_cell(l) = int_E0^El dE/Edot
        
        index1 = i            ! For the cell "l"
        index2 = i            ! For the previous cell, "l-1"
        N_p = 0.d0            ! Initialize to zero at each loop in l_min
        N_p(index1) = N_p0(i) ! This will be evolved
        l = l_min + 1  
        count = 0             ! Counter for when to change the index if losses are small
        cool = 0.d0
        new_count = .true.     

        do while( l <= ml_2 .and. index1 > 1) ! Stop at the lowest energy bin considered.
          ! At the position l, the particle has E_p(index1) instead of E_p(i)
          ! We find index1 <= i according to the cooling time:

          ! The cooling parameter "cool" determines whether particles cool efficiently or not
          if (l == l_min+1) then ! The total cooling parameter from l_min to ml_2
            do l2 = l, ml_2
              cool = cool + t_cell(l2)/t_coolp(l2,index1)
            end do
          else ! cool(l -> ml_2) = cool(1 -> ml_2) - cool(1 -> l-1)
            cool = cool - t_cell(l-1)/t_coolp(l-1,index1)
          end if
          dE_approx = (1.-t_cell(l)/t_coolp(l,index1))

          if ( 1.-cool > 0.99/Eint_p ) then ! negligible cooling
            continue

          elseif ( 1.-cool <= 0.99/Eint_p .and. dE_approx > 1./Eint_p .and. index1 > 1 ) then
            if (new_count) then
              sum = 1.d0
              steps = 0 ! number of steps until it changes E(index1) -> E(index-1)
              l2 = l
              do while (sum > 1./Eint_p .and. l2 <= ml_2)
                steps = steps + 1
                sum = sum*(1.-t_cell(l2)/t_coolp(l2,index1))
                l2 = l2 + 1
                if (sum > 1./Eint_p .and. l2 == ml_2) steps = ml_2 ! NOT COOLING
              end do
              new_count = .false.
            end if
            count = count + 1

            if (count >= steps) then
              count = 0
              index1 = index1 - 1
              new_count = .true. 
            end if

          elseif ( dE_approx <= 1./Eint_p ) then
            sum = 0.d0
            do while (sum < t_cell(l) .and. index1 > 1)
              index1 = index1 - 1    
              dE_p = E_p(index1)*(Eint_p-1.d0)
              sum = sum + dE_p/E_dotp(l,index1)        
            end do
          end if

          ! The following takes into account energy losses and variable cell size and velocity
          N_p(index1) = N_p(index2) * ( E_dotp(l,index2)/E_dotp(l,index1) ) * (t_cell(l)/t_cell(l-1))

          N_ptot(l,index1) = N_ptot(l,index1) + N_p(index1)  ! The energies match by construction

          index2 = index1 ! index2 keeps the value of the previous step

          l = l + 1 
        end do               ! End loop in shock position (l) 
        i = i + 1
      end do                  ! End loop in proton starting energies (i)
     
    end do                     ! End loop in starting position (l_min = # linear-emitters)
  
    ! Save the proton distribution at each position 
    open(79,file=output_file('dist_p'//number//'.dat'))
    do l = 1,ml_2
      do i = 1,mE_p
        write(79,*) l,E_p(i)/eV,N_ptot(l,i)
      end do
      write(79,*) 
    end do
    do l = ml_2+1,ml
      do i=1,mE_p
        write(79,*) l,E_p(i)/eV,N_ptot(l-ml_2,i)
      end do
      write(79,*) 
    end do
    close(79)
    
    open(80,file=output_file('ENTp'//number//'.dat'))
    open(81,file=output_file('PNTp'//number//'.dat'))
    do l = 1,ml_2
      E_NT = 0.d0
      U_B = B(l)**2/(8.d0*pi)
      do i = 1,mE_p
        dE_p = E_p(i)*(Eint_p-1.d0)
        E_NT = E_NT + E_p(i)*N_ptot(l,i)*dE_p
      end do
      E_NT = E_NT * 2.d0*mphi ! Total NT energy in WCR "ring"
      U_NT = E_NT/dVol(l)  
      write(80,*) l, E_NT, U_NT, U_NT/U_B! Store the energy density as well
      write(81,*) l, (U_NT/3.d0)/U_B     ! Ratio P_p,rel/P_mag
    end do
    close(80)
    close(81)

    !------------------------------------------------------------
    ! Group in bins according to distance
    if (calc_N_bin) call Np_bin(number,E_p,N_ptot,Dl,D)

    !------------------------------------------------------------

    deallocate( N_p, N_ptot )

  END subroutine CD_dist_p_run

  !============================================================================
  ! CD_dist_p_run (END)
  !============================================================================
  
END module CD_dist_p


  !--------------------------------------------------------------------------!
  !                                                                          !
  !                            SUBROUTINES                                   !
  !                                                                          !
  !--------------------------------------------------------------------------!
  
  subroutine Np_bin(number,Ep,Np_tot,dL,D_orb)
    use global
    implicit none    
    character, intent(in) :: number
    real(dp), intent(in) :: Ep(mE_p),Np_tot(ml_2,mE_p),dL(ml_2),D_orb
    real(dp) :: length, N_total_i
    integer, parameter :: N_bins = 5 ! For the plot diste_sw_bins
    real(dp) :: N_bin_i(N_bins,mE_p)

    call validate_shock_number(number,'Np_bin')
     
    N_bin_i = 0.d0
        open(81,file=output_file('dist_p'//number//'_total.dat'))
    do i = 1,mE_p
      N_total_i = 0.d0
      length = 0.d0
      do l = 1,ml_2
        ! Make some distributions binned in intervals of 0.1*D
        length = length + dL(l)
        if (length < 0.1*D_orb) then
          N_bin_i(1,i) = N_bin_i(1,i) + Np_tot(l,i)
          elseif (length >= 0.1*D_orb .and. length < 0.2*D_orb) then
          N_bin_i(2,i) = N_bin_i(2,i) + Np_tot(l,i)
          elseif (length >= 0.2*D_orb .and. length < 0.4*D_orb) then
          N_bin_i(3,i) = N_bin_i(3,i) + Np_tot(l,i)
          elseif (length >= 0.4*D_orb .and. length < 0.8*D_orb) then
          N_bin_i(4,i) = N_bin_i(4,i) + Np_tot(l,i)
          elseif (length >= 0.8*D_orb .and. length < 1.6*D_orb) then
          N_bin_i(5,i) = N_bin_i(5,i) + Np_tot(l,i)
        end if  
        N_total_i = N_total_i + Np_tot(l,i)
      end do

      write(81,*) Ep(i)/eV, N_total_i* 2.d0*mphi ! For the whole WCR   
    end do
    close(81)

    open(82,file=output_file('dist_p_wcr'//number//'_bin.dat'))
    do j = 1,N_bins
      do i = 1, mE_p
        write(82,*) j,Ep(i)/eV, N_bin_i(j,i)* 2.d0*mphi ! For the whole WCR   
      end do
      write(82,*) 
    end do
    close(82)

  end subroutine Np_bin


  !=========================================================================

   Subroutine bisect_p(Be,x3,t_coolp,t_adv,xi,R_ac)
    use global
    implicit none
    integer :: t_err
    real(dp) :: a_,b_,c_,err
    real(dp), intent(in) :: Be,t_coolp,t_adv,xi,R_ac
    real(dp), intent(out) :: x3
    interface
      function diff_p(a,B,t_pp,t_adv,efi,Rac)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp) :: diff_p
        real(dp), intent(in) :: a,B,t_pp,t_adv,efi,Rac
      end function diff_p
    end interface

    ! The error in the n-th step is < |a-b|/2**n
    a_ = log10(1.d11*eV) ! log10(Emin [CGS]) ~ 100 GeV
    b_ = log10(Emax_p) 
    t_err = 3
    err = 10.d0**(-t_err) ! error < 10**-t
    i = 0
    do while ( abs(b_ - a_) > err )
      i = i+1
      c_ = (a_+b_)/2.
      if (abs(diff_p(c_,Be,t_coolp,t_adv,xi,R_ac)) <= 1.d-6) then
        b_ = c_
      elseif (diff_p(a_,Be,t_coolp,t_adv,xi,R_ac)*diff_p(c_,Be,t_coolp,t_adv,xi,R_ac) < 0.d0) then
        b_ = c_
      elseif (diff_p(a_,Be,t_coolp,t_adv,xi,R_ac)*diff_p(c_,Be,t_coolp,t_adv,xi,R_ac) > 0.d0) then
        a_ = c_
      end if
    end do
    x3=b_

   end Subroutine bisect_p

  !=========================================================================

   Function diff_p(log10E,B_l,t_coolp,t_adv,xi,R_ac)
    use global
    implicit none
    real(dp), intent(in) :: log10E,B_l,t_coolp,t_adv,xi,R_ac
    real(dp) :: diff_p,E_l,t_diff,t_ac
    real(dp) :: D_diff,t_loss,t_esc

    E_l = 10.d0**log10E

    D_diff = E_l*c/(3.d0*qe*B_l) ! Bohm
    t_diff = R_ac**2/(2.d0*D_diff) ! Diffusion losses APPROX
    t_esc = 1.d0/(1.d0/t_adv + 1.d0/t_diff) ! escape

    t_loss = min(t_coolp,t_esc) ! general losses

    t_ac = xi*E_l/(B_l*c*qe)    
    diff_p = t_loss - t_ac
    
   end Function diff_p

  !=========================================================================

  Function fQp(E, Estiff, Emin, Emax, Ecorr)
    use global  
    real(dp) :: fQp
    real(dp),intent(in) :: E, Estiff, Emin, Emax, Ecorr
    if (E < Emin) then     ! No injection for E < E_min
      fQp = 0.d0                   
    else if (E >= Emin .and. E <= Estiff) then
      fQp = E**(-alpha)*exp(-E/Emax)
    else                   ! Hardening at E > Estiff
      fQp = E**(-alpha2)*exp(-E/Emax)*Ecorr
    endif

  end Function fQp

  !=====================================================================

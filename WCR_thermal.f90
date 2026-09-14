module WCR_thermal
  ! Calculate thermal bremsstrahlung emission from the WCR.
  use global

  implicit none
  private
  public :: WCR_thermal_run

contains

  !============================================================================
  ! WCR_thermal_run
  !============================================================================

  subroutine WCR_thermal_run()
    implicit none
    integer, parameter :: mE_ph2 = 400
    real(dp) :: Eph_min, Eph_max, Eph_int, Eph(mE_ph2)
    real(dp) :: L1_X(mE_ph2), L2_X(mE_ph2), L_X(mE_ph2)
    real(dp) :: L1_bol, L2_bol, Lwind1, Lwind2
    integer :: i

    write(*,*) 'Calculate thermal emission from the WCR'

    ! Define a broad energy range from UV to hard X-rays
    Eph_min = 5 * eV
    Eph_max = 3.d5 * eV
    call vector_log(Eph_min,Eph_max,Eph_int,Eph)

    call shock_calc_LX('1',Eph,L1_X)
    call shock_calc_LX('2',Eph,L2_X)

    L1_bol = sum(L1_X)*(Eph_int-1.d0)
    L2_bol = sum(L2_X)*(Eph_int-1.d0)
    Lwind1 = 0.5d0*Mdot1*vinf1**2
    Lwind2 = 0.5d0*Mdot2*vinf2**2
    write(*,'(A,ES10.2)') 'L_X bolometric (S1) = ',L1_bol
    write(*,'(A,ES10.2)') 'L_X bolometric (S2) = ',L2_bol
    if (L1_bol >= Lwind1) write(*,*) 'WARNING: L_X(S1) >= L_wind(S1)'
    if (L2_bol >= Lwind2) write(*,*) 'WARNING: L_X(S2) >= L_wind(S2)'

    L_X = L1_X + L2_X
    open(19,file=output_file('WCR_ff.dat'))
    do i=1,mE_ph2
      write(19,*) Eph(i)/eV,L1_X(i),L2_X(i),L_X(i)
    enddo
    close(19)

  end subroutine WCR_thermal_run

  !============================================================================
  ! Calculate the thermal spectrum for one shock.
  !============================================================================

  subroutine shock_calc_LX(number,Eph,L_X)
    implicit none
    character(len=1), intent(in) :: number
    real(dp), intent(in) :: Eph(:)
    real(dp), intent(out) :: L_X(:)
    real(dp), allocatable :: n_amb(:)
    real(dp) :: ne, ni, emi_ff, dV
    real(dp) :: gamma, mu_i, zq
    integer :: i, l, index
    real(dp) :: dummy

    call validate_shock_number(number,'shock_calc_LX')
    allocate(n_amb(ml_2))
    L_X = 0.d0

    if (number == '1') then
      gamma = gamma_1
      mu_i = mu_i1
      zq = Zq_1
    else
      gamma = gamma_2
      mu_i = mu_i2
      zq = Zq_2
    endif

    open(11,file=output_file('thermo'//number//'.dat'))
    open(12,file=output_file('Linj'//number//'.dat'))
    open(14,file=output_file('H'//number//'.dat'))
    do l=1,ml_2
      read(11,*) index,rho(l),P(l),v(l),T(l),B(l),n_a(l)
      read(12,*) x(l),y(l),dummy,dummy
      read(14,*) index,H_sh(l),dummy,dummy,dVol(l)
      n_amb(l) = rho(l)/(mp*mu_i)
    enddo
    close(11)
    close(12)
    close(14)

    do l=1,ml_2
      ni = n_amb(l)
      ne = gamma*n_amb(l)
      dV = dVol(l)
      do i=1,size(Eph)
        emi_ff = 6.8d-38*zq**2*ne*ni*T(l)**(-0.5d0) &
          *exp(-Eph(i)/(k*T(l)))*g_ff(Eph(i),T(l))
        L_X(i) = L_X(i) + emi_ff*(Eph(i)/h)*dV
      enddo
    enddo

    deallocate(n_amb)

  end subroutine shock_calc_LX

  !============================================================================
  ! WCR_thermal_run (END)
  !============================================================================

  function g_ff(Eph_X,T_X)
    real(dp), intent(in) :: Eph_X,T_X
    real(dp) :: g_ff,zeta

    g_ff = 0.d0
    if (Eph_X/(k*T_X) >= 1.d0 .and. Eph_X/(k*T_X) < 1.d3) then
      g_ff = sqrt((3.d0/pi) * ((k*T_X)/Eph_X))
    else if (Eph_X/(k*T_X) < 1.d0 .and. Eph_X/(k*T_X) >= 1.d-3) then
      zeta = 1.781d0     ! Euler's constant in the astrophysical formula
      g_ff = (sqrt(3.d0)/pi) * log(4.d0*k*T_X/(zeta*Eph_X))
    else if (Eph_X/(k*T_X) > 1.d3 .or. Eph_X/(k*T_X) < 1.d-3) then
      g_ff = 0.d0
    endif
  end function g_ff

END module WCR_THERMAL

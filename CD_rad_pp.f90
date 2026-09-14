module CD_rad_pp
  use global
  implicit none
  
  public :: CD_rad_pp_run,vec_Epp
  

contains

  !============================================================================
  ! CD_rad_pp_run (START)
  !============================================================================
  ! This routine calculates the intrinsic and absorbed p-p emission using as an 
  ! input the proton distribution in distp(number).dat. The integrated 
  ! luminosity (unabs and abs) for the linear emitter is stored in radpp.dat.
  !============================================================================
  
  subroutine CD_rad_pp_run(number,Ef,Lpp_sum,Lpp_sum_abs,D,Xobs)
    character, intent(in) :: number
    real(dp), intent(in) :: Ef(mE_f), D,Xobs(3)
    real(dp), intent(out) :: Lpp_sum(mE_f), Lpp_sum_abs(mE_f)
    real(dp) :: d1, d2, Q_Ef, Q_phi, dx_pp, x_pp, qpi
    real(dp) :: Ef_int, dEf, sum_pp, Epp, L1_low, L1(mE_p), integral, E_p_read
    real(dp) :: Epi_min, Epi_max, Epi_int, dE_pi
    real(dp) :: ang_min, ang_int, ang1, ang2
    real(dp) :: dil1, dil2
    real(dp), parameter :: K_pi=0.17d0, n_pp=0.9d0, kappa_pp=K_pi*n_pp! n_pp=1.1 for alfa=2, n_pp=0.86 for alpha=2.5
    integer :: nl, iang1, iang2, iEf, mE_pi, i_p
    !     ml - CD ; mE_p - protons ; mE_f - photons ; mE_pi(variable) - pions
    real(dp), dimension(:), allocatable :: Epi               
    interface
      function F(a,b)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp) :: F
        real(dp), intent(in) :: a,b
      end function F
    end interface

    call validate_shock_number(number,'CD_rad_pp_run')

    allocate( Lpp(mE_f) )
    allocate( N_ptot(ml,mE_p),N_p(mE_p) )
    allocate( tau1(mE_gg,m_ang_gg),tau2(mE_gg,m_ang_gg) )

80  format(A10,X,Es10.2)

    Lpp_sum_abs = 0.d0
    Lpp_sum = 0.d0

    print *, 'Calculate p-p emission in S',number

    !------------------------------------------------------------------
    ! Read positions (x(l),y(l)) and thermodynamical quantities at them
    
    open(11,file=output_file('thermo'//number//'.dat'))
    do l = 1,ml
      read(11,*) nl,rho(l),P(l),v(l),T(l),B(l),n_a(l) ! cgs units
    end do
    close(11)

    if (number=='1') then 
      z_pp = z_pp1
    else
      z_pp = z_pp2
    end if

    ! Proton distribution at position l with energy E_p(i): N_ptot(l,i) 
    open (14,file=output_file('dist_p'//number//'.dat'))
    do l = 1,ml
      do i = 1,mE_p
        read(14,*) j, E_p_read,N_ptot(l,i)
        if (l==1) E_p(i) = E_p_read*eV    !   convert energy from eV to erg
      end do
      read(14,*)
    end do
    close(14) 
    Eint_p = E_p(2)/E_p(1)
    L1 = log(E_p/TeV)         ! Pre-compute: L = ln(Ep/1TeV), needed later

    ! Matrix with opacity coefficients at a distance of R_star for different angles
    ! between the axis between the emitter and the star, and the LOS; and for 
    ! different energies of the emitted photons: tau(Ef,angle) [one per star]
    
    Ef_int_gg = (Ef_max_gg/Ef_min_gg)**(1.d0/dfloat(mE_gg-1))
    
    open (188,file=output_file('opac_gg1.dat'))
    open (189,file=output_file('opac_gg2.dat'))
    do i = 1,mE_gg
      read(188,*) (tau1(i,j), j=1,m_ang_gg)
      read(189,*) (tau2(i,j), j=1,m_ang_gg)
    end do
    close(188)
    close(189)
    ang_int = pi/float(m_ang_gg)
    ang_min = ang_int/2.d0

    Xst1 = (/ 0.d0, 0.d0, 0.d0 /)
    Xst2 = (/ D, 0.d0, 0.d0 /)


    !=======================================================================
    !  Begin p-p calculations 
    !=======================================================================

    Ef_int = Ef(2)/Ef(1)
    ! Calculate p-p emission at each position and then correct for GGA and sum
    Epi_max = E_p(mE_p)
    sum_pp = 0.d0
    DO l = 1,ml          
      N_p = N_ptot(l,:) ! Work with a vector instead of array for simplicity
      Xemi = (/ x(l),y(l),z(l) /) 

      d1 = sqrt( Xemi(1)**2 + Xemi(2)**2 + Xemi(3)**2 )
      call unit_vector_diff(Xemi,Xst1,3,XEe) ! unit vector from 2nd to 1st
      ang1 = acos(DOT_PRODUCT(Xobs,XEe))! angle between Xobs and XEe
      iang1 = int(m_ang_gg*abs(ang1-ang_min)/pi)+1
      dil1 = Rst1/d1 

      d2 = sqrt( (D-Xemi(1))**2 + Xemi(2)**2 + Xemi(3)**2 )
      call unit_vector_diff(Xemi,Xst2,3,XEe) ! unit vector from 2nd to 1st
      ang2 = acos(DOT_PRODUCT(Xobs,XEe))! angle between Xobs and XEe
      iang2 = int(m_ang_gg*abs(ang2-ang_min)/pi)+1
      dil2 = Rst2/d2

      !------------------------------------------------------------------------

      mE_pi = mE_f
      do i = 1,mE_f           ! Ef(i) 
        dEf = Ef(i)*(Ef_int-1.d0)

        ! Different formalism according to energy range:
        if (Ef(i) < 0.01d0*TeV) then    ! Ep < 100 GeV -> Ef < 0.01 TeV
          allocate(Epi(mE_pi))
          Epi_min = Ef(i) + (m_pic2)**2/(4.d0*Ef(i))
          Epi_min = Epi_min*1.005d0     ! Factor to avoid numerical errors very close to Epi_min
          call vector_log(Epi_min,Epi_max,Epi_int,Epi)
          integral = 0.d0
          do j = 1,mE_pi       ! E_pi
            dE_pi = Epi(j)*(Epi_int-1.d0)
            Epp = Epi(j)/kappa_pp + mpc2           ! Eq.75 Kelner+2006
            i_p = int(log(Epp/Emin_p)/log(Eint_p))+1
            if (i_p <= mE_p .and. Epp <= E_p(mE_p) .and. Epp > E_th) then
              L1_low = log(Epp/TeV)
              qpi = n_pp*N_p(i_p)*sigma_pp(L1_low,Epp)/K_pi ! Eq.77 Kelner+2006 [1/erg/s]
              integral = integral + qpi*dE_pi/sqrt(Epi(j)**2-m_pic2**2)
            endif 
          end do              ! j,E_pi
          Q_Ef = 2.d0*c*n_a(l)*integral*z_pp    !Eq.78 Kelner+2006
          Lpp(i) = Ef(i)*Q_Ef
          deallocate(Epi)
          mE_pi = MAX(mE_pi-1,32) ! Diminish the number of points 

        elseif (Ef(i) >= 0.01d0*TeV) then  ! Ep > 100 GeV -> Ef > 10 GeV
          integral=0.d0
          do j = 1,mE_p
            if (Ef(i) < E_p(j)) then
              dE_p = E_p(j)*(Eint_p-1.d0)
              x_pp = Ef(i)/E_p(j)
              dx_pp = Ef(i)*dE_p/(E_p(j)**2)
              integral = integral + sigma_pp(L1(j),E_p(j))*N_p(j)*F(x_pp,L1(j))*dx_pp/x_pp
            end if
          end do
          Q_phi = c*n_a(l)*integral*z_pp
          Lpp(i) = Ef(i)*Q_phi

        end if
         
        sum_pp = sum_pp + Lpp(i)*dEf
        Lpp_sum(i) = Lpp_sum(i) + Lpp(i)

        ! Correct for GGA
        if (Ef(i) > Ef_min_gg .and. Ef(i) < Ef_max_gg) then
          iEf = int(log(Ef(i)/Ef_min_gg)/log(Ef_int_gg)) + 1
          Lpp_sum_abs(i) = Lpp_sum_abs(i) + Lpp(i)*exp( -(Rst1/d1)*&
            tau1(iEf,iang1) -(Rst2/d2)*tau2(iEf,iang2) )
        else
          Lpp_sum_abs(i) = Lpp_sum_abs(i) + Lpp(i)
        end if

      end do    ! i

    end do      ! l

    write(*,80) 'L_pp=',sum_pp

    !------------------------------------------------------------

    open(16,file=output_file('rad_pp'//number//'.dat'))
    do i=1,mE_f
      write(16,*) Ef(i)/eV,Ef(i)*Lpp_sum(i),Ef(i)*Lpp_sum_abs(i)
    end do
    close(16)

    deallocate( Lpp )
    deallocate( N_ptot,N_p )
    deallocate( tau1,tau2 )


  END subroutine CD_rad_pp_run

  !============================================================================
   ! CD_rad_pp_run (END)
  !============================================================================
  

  !============================================================================
  ! vec_Epp
  !============================================================================
  
  subroutine vec_Epp(Epp)
    use global
    implicit none
    real(dp), intent(out) :: Epp(mE_f)
    real(dp) :: Eaux, Ep_max, Epp_max, Epp_int, dummy

    !--------------------------------------------
    ! Obtain Ep_max (this is the cutoff energy, not the max energy)
    Epp_max = 0.d0
    open (14,file=output_file('Emax_p1.dat'))
    open (15,file=output_file('Emax_p2.dat'))
    do l=1,ml_2
      read(14,*) dummy, Eaux
      Ep_max = max(Eaux,Epp_max)
      read(15,*) dummy, Eaux
      Ep_max = max(Eaux,Epp_max)
    end do
    close(14) 
    close(15) 
    Epp_max = 100.d0 * Ep_max * eV

    !--------------------------------------------------
    ! p-p photons have E > 10 MeV
    call vector_log(1.d7*eV,Epp_max,Epp_int,Epp)

  END subroutine vec_Epp

!================================================================
  

END module CD_rad_pp


!========================================

Function F(x,L)
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: x,L
  real(dp) :: B,k,beta,F,L2
  
  L2 = L**2

  B = 1.30d0 + 0.14d0*L + 0.011d0*L2
  beta = 1.d0/(1.79d0 + 0.11d0*L + 0.008d0*L2)
  k = 1.d0/(0.801d0 + 0.049d0*L + 0.014d0*L2)
  F = B*log(x)/x*((1.d0-x**beta)/(1.d0+k*x**beta*(1.d0-x**beta)))**4 &
    *(1.d0/log(x)-(4.d0*beta*x**beta)/(1.d0-x**beta) &
    -(4.d0*k*beta*x**beta*(1.d0-2.d0*x**beta)/(1.d0+k*x**beta*(1.d0-x**beta))))
  
  return
end Function F

module opac_ff
  use global
  implicit none

  public :: opac_ff_run, opac_ff_data_run

  contains

  !============================================================================
  ! opac_ff_run
  !============================================================================
  ! Calculate f-f emission and absorption coefficients for ionised stellar winds
  ! Assume spherically symmetric stellar winds --> tau_ff(r), emi_ff(r)
  ! Output files:
  ! 1) opac_ff<i>.dat contains a table for tau_ff_<i>(nu,ang)
  ! 2) windemiff<i>.dat contains tables with nu(Hz) and emissivity coefficient 
  ! epsff(l,nu)*dl^3

  subroutine opac_ff_run(D)
    !$ use omp_lib
    real(dp), intent(in) :: D
    real(dp) :: nu_ff(mE_ff),S,S1,S2,dummy

    interface
      subroutine tau_ff_calc(number,nu_ff,D)
        use global
        implicit none
        real(dp), intent(in) :: nu_ff(mE_ff),D
        character, intent(in) :: number
      end subroutine tau_ff_calc
    end interface 

    interface
      subroutine interpol(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
        real(dp), intent(out) :: finterpol
      end subroutine interpol
    end interface

    interface
      subroutine S_ff_map(B_ff_D2,cte_Tau,r_i,r_f,S_ff_cm2,r_map)
        use global
        real(dp), intent(in) :: B_ff_D2,cte_Tau,r_i,r_f
        real(dp), intent(out) :: r_map(mr_map),S_ff_cm2(mr_map)
      end subroutine S_ff_map
    end interface

    interface
      subroutine S_ff_pointmap(S_tot,r_i,r_f,S_map,r_map)
        use global
        real(dp), intent(in) :: S_tot,r_i,r_f
        real(dp), intent(out) :: r_map(mr_map),S_map(mr_map)
      end subroutine S_ff_pointmap
    end interface

    !============================================================

    print *, 'Calculate f-f absorption and emission coefficients of the stellar winds'
    
    ! Frequencies at which tau_ff is calculated
    call vector_log(nu_min_ff,nu_max_ff,nu_int_ff,nu_ff)

    ! Calculate absorption only for nu <~ 200 GHz (scaling with D).
    ! nu_abs_max = 200.d9*((100.d0*UA)/D)**2 ! more absorption when closer
    ! inu_max = int( log(nu_abs_max/nu_ff(1))/log_nu_int ) + 1

    ! Calculate opacity tables for each stellar wind 
    call tau_ff_calc('1',nu_ff,D)     ! stellar wind 1
    call tau_ff_calc('2',nu_ff,D)     ! stellar wind 2
    
    S1 = 0.d0
    S2 = 0.d0
    S = 0.d0
    open(19,file=output_file('windemiff1.dat'))
    open(20,file=output_file('windemiff2.dat'))
    open(21,file=output_file('windemiff.dat'))
    do i = 1, mE_ff
      read(19,*) dummy,S1 ! stored in mJy
      read(20,*) dummy,S2 ! stored in mJy
      S = S1 + S2
      write(21,*) nu_ff(i),S   
    end do
    close(19)
    close(20)
    close(21)

  END subroutine opac_ff_run

  !============================================================================
  ! opac_ff_data_run
  !============================================================================
  ! Calculate f-f emission and absorption coefficients at the observed frequencies
  ! Output files:
  ! 1) opac_ff_data<i>.dat contains a table for tau_ff_<i>(nu,ang)
  ! 2) windemiff_data<i>.dat contains tables with nu(Hz) and emissivity coefficient 
  ! epsff(l,nu)*dl^3

  subroutine opac_ff_data_run(nu_data)
    !$ use omp_lib
    real(dp), intent(in) :: nu_data(mnu_data)
    real(dp) :: S, S1, S2, dummy

    ! Calculate opacity tables for each stellar wind 

    call tau_ff_data_calc('1', nu_data)     ! stellar wind 1
    call tau_ff_data_calc('2', nu_data)     ! stellar wind 2

    S1 = 0.d0
    S2 = 0.d0
    S = 0.d0
    open(19,file=output_file('windemiff_data1.dat'))
    open(20,file=output_file('windemiff_data2.dat'))
    open(21,file=output_file('windemiff_data.dat'))
    do i = 1, mnu_data
      read(19,*) dummy,S1 ! stored in mJy
      read(20,*) dummy,S2 ! stored in mJy
      S = S1 + S2
      write(21,*) nu_data(i),S   
    end do
    close(19)
    close(20)
    close(21)

  END subroutine opac_ff_data_run

!==============================================================================

END module opac_ff


  !============================================================================
  ! gaunt_ff_calc 
  !============================================================================

  Subroutine gaunt_ff_calc(nu_ff,T_w,Zq_w,gaunt_ff)
    use global  
    implicit none
    real(dp), intent(in) :: nu_ff(:),T_w,Zq_w
    real(dp), intent(out) :: gaunt_ff(:)
    integer, parameter :: nu_mu = 146, nu_mga=81 ! data dimensions in gauntff.dat
    real(dp) :: gff(nu_mu,nu_mga),gamma2,dg,dmu,gff_inter
    real(dp), allocatable :: mu_ff(:)
    integer :: j_mu,j_g
    integer :: m_nu_ff

    interface
      subroutine interpol(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
        real(dp), intent(out) :: finterpol
      end subroutine interpol
    end interface
    
    m_nu_ff = SIZE(nu_ff)

    !--------------------------------------------------------------------------
    ! gauntff.dat (van Hoof 2014) allows to compute the Gaunt factor:
    ! <gff(u,gamma2)>, with -8 < log10u < 8, -4 < log10gamma < 4, step 0.2
      
    open(13,file='gauntff.dat')
    do i=1,16
      read(13,*)
    end do
    do i=1,nu_mu
      read(13,*) (gff(i,j), j=1,nu_mga)
    end do
    close(13)

    gamma2 = Zq_w**2*Ry/k/T_w
    j_g = int( (log10(gamma2)+6.d0)/0.2d0 ) + 1
    dg = ( (log10(gamma2)+6.d0)/0.2d0+1.d0 - float(j_g) )*0.2d0
 
    allocate( mu_ff(m_nu_ff) )
    mu_ff = h*nu_ff/(k*T_w)

    do i=1,m_nu_ff 
      if (h*nu_ff(i) > 3*k*T_w) print*, 'WARNING: h nu > kT in opac_ff (gaunt_ff_calc)'
      ! index = (value - table_value(1))/(table_step) + 1, then do linear interpolation
      j_mu = int( (log10(mu_ff(i))+16.d0)/0.2d0 ) + 1
      dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2_dp
      call interpol(dg,dmu,gff(j_mu, j_g),gff(j_mu, j_g+1),gff(j_mu+1, j_g),gff_inter)      
      gaunt_ff(i) = gff_inter
      !gaunt_ff(i) = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w/nu_ff(i))) ! approx. in De Becker+07
    end do  
    deallocate( mu_ff )

  End subroutine gaunt_ff_calc

  !============================================================================
  ! tau_ff_calc<i>  (i=1,2 is the shock number)
  !============================================================================
  ! The outputs of the program are:
  ! 1) the opacity tau_ff<i>(nu,ang)
  ! 2) windemiff<i>.dat, containing: nu [Hz], S<i>_nu [mJy]

  subroutine tau_ff_calc(number,nu_ff,D)
    use global
    implicit none
    character, intent(in) :: number
    real(dp),intent(in) :: nu_ff(mE_ff),D
    integer, parameter :: ir_max = 300   ! steps in integration in r
    real(dp) :: Mdot,Mdot_eff,vinf,Rst,T_w,mu_w,Zq_w,gamma_w
    real(dp) :: g_ff(mE_ff),S(mE_ff)
    real(dp) :: alpha_ff,cte_alpha_ff,cte_nu,dr,A_ff,ne,ni
    real(dp) :: tau,dang,ang,cang,sang,xi,yi,ri
    real(dp) :: f(ir_max)
    integer :: iang,ie,ir,imap
    !Variables for radio maps:
    real(dp) :: log_nu_int
    !real(dp) :: B_ff,B_ff_D2,K_ff,cte_Tau ! Needed for extended stellar winds

    interface
      subroutine gaunt_ff_calc(nu_ff,Tw,Zw,gaunt)
        use global
        implicit none
        real(dp), intent(in) :: nu_ff(:),Tw,Zw
        real(dp), intent(out) :: gaunt(:)
    end subroutine gaunt_ff_calc
    end interface

    call validate_shock_number(number,'tau_ff_calc')

    allocate( tau_ff(mE_ff,m_ang_ff) )

    !----------------------------------------------------------
    ! INITIALIZATION

    if (number=='1') then 
      Mdot = Mdot1
      vinf = vinf1
      Rst = Rst1
      mu_w = mu_w1
      Zq_w = Zq_w1
      gamma_w = gamma_w1
      T_w = T_w1    ! assume T_w = cte for simplicity (the distances are r >> R_*)

    else
      Mdot = Mdot2
      vinf = vinf2
      Rst = Rst2
      mu_w = mu_w2
      Zq_w = Zq_w2
      gamma_w = gamma_w2
      T_w = T_w2    ! assume T_w = cte for simplicity (the distances are r >> R_*)

    end if
    Mdot_eff = Mdot/sqrt(filling) ! effective mass-loss rate (account for clumping)

    !if (ri > Rstar) v_w = 1d5 + vW_inf*(1d0 - Rstar/ri)**0.8d0 !vW = vW_0 + vW_inf*(1-R/r)^0.8
    A_ff = Mdot_eff/(4.d0*pi*vinf*mu_w*mp)  ! Assume r >> Rstr --> v_w = v_infty
    cte_alpha_ff = 3.7d8*Zq_w**2/sqrt(T_w)
    log_nu_int = log(nu_ff(2)/nu_ff(1))

    !----------------------------------------------------------
    ! CALCULATIONS BEGIN

    call gaunt_ff_calc(nu_ff,T_w,Zq_w,g_ff)

    dang = pi/m_ang_ff
    dr = Rst/32.d0    ! step for integration: R_final ~ Rstar + Rstar*ir_max/16*cos(ang)
    do ie=1,mE_ff      ! loop in emitted photon energy: Ef(ie) = h*nu_ff(ie)
      ang = dang/2.d0  ! initial angle
      cte_nu = nu_ff(ie)**(-3)*(1.d0-exp(-h*nu_ff(ie)/k/T_w))*g_ff(ie)

      do iang=1,m_ang_ff ! loop in interaction angle: ang
        cang = cos(ang)
        sang = sin(ang)
        xi = Rst ! initial position (taken > Rstar for simplicity)
        yi = 0.d0
        tau = 0.d0
        do ir=1,ir_max  ! loop in emitted photon trayectory: xi,yi
          xi = dr*cang+xi
          yi = dr*sang+yi
          ri = sqrt(xi**2 + yi**2)
          ni = A_ff/ri**2/(gamma_w+1.d0)  ! mui = mu*(g_e+1)
          ne = gamma_w*ni                      
          alpha_ff = cte_alpha_ff*ne*ni*cte_nu !1/cm
          f(ir) = alpha_ff ! The photon energies are the same
          ! tau = tau + alpha_ff*dr
        end do
        call integrate(ir_max,dr,f,tau)
        tau_ff(ie,iang)=tau
        ang=ang+dang
      end do
    end do

    open (88,file=output_file('opac_ff'//number//'.dat'))
    do i=1,mE_ff
      write(88,*) (tau_ff(i,j), j=1,m_ang_ff)
    end do
    close(88)

    deallocate( tau_ff )

    !------------------------------------------------------------------------------
    !     
    ! The total radio flux from the wind is given by Eq. 8 in Wright 1975:

    open(19,file=output_file('windemiff'//number//'.dat'))
    S = 0.d0

    do i = 1, mE_ff     
      S(i) = (23.2d0/(distance)**2) * ( ((Mdot_eff/M_sun_yr)/&
        (vinf/1.d5)/mu_w)**2 *nu_ff(i)*gamma_w*g_ff(i)*Zq_w**2 )**(2.d0/3.d0) ! Jy = 1d-26 W/m/Hz
      if (nu_ff(i) > 1.d-3*k*T_w/h) then    ! original formula valid for h*nu << k*T
        S(i) = S(i)*exp(-h*nu_ff(i)/k/T_w)  ! else assume exponential correction factor
      endif
      write(19,*) nu_ff(i),S(i)*1.d3 ! Stored in mJy
    end do
    close(19)

    !------------------------------------------------------------------------------
    ! Make the radio maps for all map frequencies
      
    do imap = 1, n_map_freq
      call write_ff_map(number,D,map_freq(imap),nu_ff,S,Rst,log_nu_int)
    end do

  END subroutine tau_ff_calc


  !============================================================================
  ! write_ff_pointmap
  !============================================================================
  ! Write the f-f flux from the stellar winds in a grid to add to the WCR emission maps.
  ! The stellar winds are spherical, so we use polar coordinates. 

  subroutine write_ff_map(number,D,map_nu,nu_ff,S,Rst,log_nu_int)
    use global
    implicit none
    character, intent(in) :: number
    real(dp), intent(in) :: D,map_nu,nu_ff(mE_ff),S(mE_ff),Rst,log_nu_int
    integer :: ii,jj,i_nu
     character(len=16) :: freq_label
    real(dp) :: r_map(mr_map),S_ff_cm2(mr_map)
    real(dp) :: convert,phi_map(mphi),x_map,y_map

    interface
      subroutine S_ff_pointmap(S_tot,r_i,r_f,S_map,r_map)
        use global
        implicit none
        real(dp), intent(in) :: S_tot,r_i,r_f
        real(dp), intent(out) :: S_map(mr_map),r_map(mr_map)
      end subroutine S_ff_pointmap
    end interface

    i_nu = int(log(map_nu/nu_ff(1))/log_nu_int) + 1
    ! Alternative extended-wind map:
    !B_ff = 2.*nu_ff(i_nu)**2*k*T_w/c**2    ! Planck function at the h nu << kT limit
    !B_ff_D2 = B_ff/(distance*kpc)**2*1.d26 ! in mJy
    !K_ff = 3.7d8*h*Zq_w**2*g_ff(i_nu)/(k*T_w**1.5*nu_ff(i_nu)**2)
    !cte_Tau = K_ff*gamma_w*A_ff**2*pi/2.
    !call S_ff_map(B_ff_D2,cte_Tau,Rst,D/(1.d0+sqrt(eta)),S_ff_cm2,r_map)
    call S_ff_pointmap(S(i_nu),Rst,D/(1.d0+sqrt(eta)),S_ff_cm2,r_map)

    freq_label = map_frequency_label(map_nu)
    open(160,file=output_file('S'//number//'_ff_map_'//trim(freq_label)//'G.dat'))

    phi_map(1) = 0.d0
    do jj = 2, mphi
      phi_map(jj) = phi_map(jj-1) + delta_phi
    end do
    convert = 206265.d0/(distance*kpc) * 1.d3

    do ii = 1, mr_map
      !A_cell(i) = r_map(i)*Delta_r*d_phi
      !S_map(i) = S_ff_cm2(i)*A_cell(i)
      do jj = 1, mphi
        y_map = r_map(ii) * sin(phi_map(jj)) * convert
        if (number == '1') then
          x_map = r_map(ii) * cos(phi_map(jj)) * convert
        else
          x_map = (r_map(ii) * cos(phi_map(jj)) + D) * convert
        endif
        write(160,*) x_map,y_map,S_ff_cm2(ii)/mphi*1.d3  ! Stored in mJy
      end do
    end do
    close(160)
  end subroutine write_ff_map

  !===============================================================================

  !============================================================================
  ! tau_ff_data_calc<i> 
  !============================================================================
  ! The outputs of the program are:
  ! 1) 'opac_ff_data<i>.dat' --> contains: opacity tau_ff<i>(nu_data,ang)
  ! 2) 'windemiff_data<i>.dat' --> contains: nu_data [Hz], S<i>_nu [mJy]

  subroutine tau_ff_data_calc(number,nu_ff)
    use global
    implicit none
    character, intent(in) :: number
    real(dp),intent(in) :: nu_ff(mnu_data)
    integer, parameter :: ir_max = 300   ! steps in integration in r
    real(dp) :: Mdot,Mdot_eff,vinf,Rst,T_w,mu_w,Zq_w,gamma_w
    real(dp) :: g_ff(mnu_data),S(mnu_data)
    real(dp) :: alpha_ff,cte_alpha_ff,cte_nu,dr,A_ff,ne,ni
    real(dp) :: tau,dang,ang,cang,sang,xi,yi,ri
    real(dp) :: f(ir_max)
    integer :: iang,ie,ir
    interface
      subroutine gaunt_ff_calc(nu_ff,Tw,Zw,gaunt)
        use global
        implicit none
        real(dp), intent(in) :: nu_ff(:),Tw,Zw
        real(dp), intent(out) :: gaunt(:)
      end subroutine gaunt_ff_calc
    end interface 

    call validate_shock_number(number,'tau_ff_data_calc')

    allocate( tau_ff(mnu_data,m_ang_ff) )

    !----------------------------------------------------------
    ! INITIALIZATION

    if (number=='1') then 
      Mdot = Mdot1
      vinf = vinf1
      Rst = Rst1
      mu_w = mu_w1
      Zq_w = Zq_w1
      gamma_w = gamma_w1
      T_w = T_w1    ! assume T_w = cte for simplicity (the distances are r >> R_*)

    else
      Mdot = Mdot2
      vinf = vinf2
      Rst = Rst2
      mu_w = mu_w2
      Zq_w = Zq_w2
      gamma_w = gamma_w2
      T_w = T_w2    ! assume T_w = cte for simplicity (the distances are r >> R_*)

    end if
    Mdot_eff = Mdot/sqrt(filling) ! effective mass-loss rate (account for clumping)

    !if (ri > Rstar) v_w = 1d5 + vW_inf*(1d0 - Rstar/ri)**0.8d0 !vW = vW_0 + vW_inf*(1-R/r)^0.8
    A_ff = Mdot_eff/(4.d0*pi*vinf*mu_w*mp)  ! Assume r >> Rstr --> v_w = v_infty
    cte_alpha_ff = 3.7d8*Zq_w**2/sqrt(T_w)

    !----------------------------------------------------------
    ! CALCULATIONS BEGIN

    call gaunt_ff_calc(nu_ff,T_w,Zq_w,g_ff)

    dang = pi/m_ang_ff
    dr = Rst/32.d0    ! step for integration: R_final ~ Rstar + Rstar*ir_max/16*cos(ang)
    do ie=1,mnu_data      ! loop in emitted photon energy: Ef(ie) = h*nu_ff(ie)
      ang = dang/2.d0  ! initial angle
      cte_nu = nu_ff(ie)**(-3)*(1.d0-exp(-h*nu_ff(ie)/k/T_w))*g_ff(ie)
      do iang=1,m_ang_ff ! loop in interaction angle: ang
        cang = cos(ang)
        sang = sin(ang)
        xi = Rst ! initial position (taken > Rstar for simplicity)
        yi = 0.d0
        tau = 0.d0
        do ir=1,ir_max  ! loop in emitted photon trayectory: xi,yi
          xi = dr*cang+xi
          yi = dr*sang+yi
          ri = sqrt(xi**2 + yi**2)
          ni = A_ff/ri**2/(gamma_w+1.d0)  ! mui = mu*(g_e+1)
          ne = gamma_w*ni                      
          alpha_ff = cte_alpha_ff*ne*ni*cte_nu !1/cm
          f(ir) = alpha_ff ! The photon energies are the same
          ! tau = tau + alpha_ff*dr
        end do
        call integrate(ir_max,dr,f,tau)
        tau_ff(ie,iang)=tau
        ang=ang+dang
      end do
    end do

    open (88,file=output_file('opac_ff_data'//number//'.dat'))
    do i=1,mnu_data
      write(88,*) (tau_ff(i,j), j=1,m_ang_ff)
    end do
    close(88)

    deallocate( tau_ff )

    !------------------------------------------------------------------------------
    !     
    ! The total radio flux from the wind is given by Eq. 8 in Wright 1975:

    open(19,file=output_file('windemiff_data'//number//'.dat'))
    S = 0.d0
    do i = 1, mnu_data     
      S(i) = (23.2d0/(distance)**2) * ( ((Mdot_eff/M_sun_yr)/&
        (vinf/1.d5)/mu_w)**2 *nu_ff(i)*gamma_w*g_ff(i)*Zq_w**2 )**(2.d0/3.d0) ! Jy = 1d-26 W/m/Hz
      if (nu_ff(i) > 1.d-3*k*T_w/h) then    ! original formula valid for h*nu << k*T
        S(i) = S(i)*exp(-h*nu_ff(i)/k/T_w)  ! else assume exponential correction factor
      endif
      write(19,*) nu_ff(i), S(i)*1.d3 ! Stored in mJy
    end do
    close(19)

  END subroutine tau_ff_data_calc

  !===============================================================================

  Subroutine interpol(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
    implicit none 
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
    real(dp), intent(out) :: finterpol
    real(dp) :: fx,fy

    fx = (fx2y1-fx1y1)/0.2    ! x2-x1 = 0.2
    fy = (fx1y2-fx1y1)/0.2    ! y2-y1 = 0.2
    finterpol = fx1y1 + fx*dx + fy*dy
  end Subroutine interpol
    
  !===============================================================================


  Subroutine S_ff_map(B_ff_D2,cte_Tau,r_i,r_f,S_ff_cm2,r_map)
    ! f-f emission from a spherical wind with a boundary at R_c following Wright+1975
    ! relevant for emission maps in which the stars are not point-like.
    ! WARNING: lacks testing 
    use global
    implicit none
    real(dp), intent(in) :: B_ff_D2,cte_Tau,r_i,r_f
    real(dp), intent(out) :: S_ff_cm2(mr_map),r_map(mr_map)
    real(dp) :: delta_r,Tau_max,R_c,coc_r

    R_c = 2.*r_i 
    delta_r = (r_f - r_i)/float(mr_map-1)
    r_map(1) = r_i
    do i = 2,mr_map
      r_map(i) = r_map(i-1) + delta_r
    end do
    
    do i = 1,mr_map
      ! Eq. 5 from Wright+1975
      if ( r_map(i) >= R_c ) then 
        Tau_max = cte_Tau/r_map(i)**3
      else if ( r_map(i) < R_c ) then
        coc_r = r_map(i)/R_c
        Tau_max = cte_Tau/r_map(i)**3* ( 0.5 - acos(coc_r)/pi -&
          sqrt(1.-coc_r**2)*coc_r/pi )
      endif
      S_ff_cm2(i) = B_ff_D2*(1.-exp(-Tau_max))
    end do
    
  END Subroutine S_ff_map
    
  !===============================================================================


  Subroutine S_ff_pointmap(S_tot,r_i,r_f,S_map,r_map)
    ! set S=S_tot at r=R_st (central point), and S=0 elsewhere. 
    ! Valid for R_ff << theta_beam.
    use global
    implicit none
    real(dp), intent(in) :: S_tot,r_i,r_f
    real(dp), intent(out) :: S_map(mr_map),r_map(mr_map)
    real(dp) :: delta_r

    delta_r = (r_f - r_i)/float(mr_map-1)
    r_map(1) = r_i
    do i = 2,mr_map
      r_map(i) = r_map(i-1) + delta_r
    end do
    
    S_map = 0.d0
    S_map(1) = S_tot
  END Subroutine S_ff_pointmap

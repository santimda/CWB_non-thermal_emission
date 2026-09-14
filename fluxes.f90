MODULE fluxes
  use global
  implicit none

  private
  public :: fluxes_calc

contains

  !========================================================================
  !        fluxes_calc        
  !========================================================================
  !	Computation of the fluxes at different bands
  !************************************************************************

  subroutine fluxes_calc(E_sy,L_sy,E_br,L_br,E_IC,L_IC,E_pp,L_pp,flux)
    real(dp), intent(in) :: E_sy(mnu),L_sy(mnu)
    real(dp), intent(in) :: E_br(mE_f),L_br(mE_f)
    real(dp), intent(in) :: E_IC(mE_f),L_IC(mE_f)
    real(dp), intent(in) :: E_pp(mE_f),L_pp(mE_f)
    real(dp), intent(out) :: Flux(N_bands)
    real(dp) :: Emins(N_bands),Emaxs(N_bands),phFlux,phflux2
    real(dp) :: F_sy(N_bands),F_br(N_bands),F_IC(N_bands),F_pp(N_bands)
    real(dp) :: dilution,Ein_sy,Ein_br,Ein_IC,Ein_pp
    real(dp) :: dE_sy(mnu),dE_br(mE_f),dE_IC(mE_f),dE_pp(mE_f)
    dilution = 4d0*pi*(distance*kpc)**2

    Ein_sy = E_sy(2)/E_sy(1) !Energy increase factor
    do i = 1,mnu
      dE_sy(i) = E_sy(i)*(Ein_sy-1d0)			!Energy interval
    end do

    Ein_br = E_br(2)/E_br(1) !Energy increase factor
    Ein_IC = E_IC(2)/E_IC(1) !Energy increase factor
    Ein_pp = E_pp(2)/E_pp(1) !Energy increase factor
    do i = 1,mE_f
      dE_br(i) = E_br(i)*(Ein_br-1d0)     !Energy interval
      dE_IC(i) = E_IC(i)*(Ein_IC-1d0)     !Energy interval
      dE_pp(i) = E_pp(i)*(Ein_pp-1d0)     !Energy interval
    end do
    
    !Bands:  3-10 keV, 10-30 keV, 0.1-100 MeV, 0.1-10 GeV, 0.1-100 GeV, 0.1-100 TeV
    Emins = (/ 3d3*eV,1.d4*eV,1d5*eV,1d8*eV,1.d8*eV,1d11*eV /)
    Emaxs = (/ 1d4*eV,3.d4*eV,1d8*eV,1d10*eV,1d11*eV,1d13*eV /)
    !Bands:  3-10 keV, 8-18 keV, 0.1-100 MeV, 0.1-10 GeV, 0.1-100 GeV, 0.1-100 TeV
    !Emins = (/ 3d3*eV,8.d3*eV,1d5*eV,1d8*eV,1.d8*eV,1d11*eV /)
    !Emaxs = (/ 1d4*eV,18d3*eV,1d8*eV,1d10*eV,1d11*eV,1d13*eV /)

    !Calculate the energy flux in each band [erg cm^-2 s^-1]
    do j = 1,N_bands
      F_sy(j) = flux_int(mnu,E_sy,dE_sy,L_sy,Emins(j),Emaxs(j))   
      F_br(j) = flux_int(mE_f,E_br,dE_br,L_br,Emins(j),Emaxs(j)) 
      F_IC(j) = flux_int(mE_f,E_IC,dE_IC,L_IC,Emins(j),Emaxs(j)) 
      F_pp(j) = flux_int(mE_f,E_pp,dE_pp,L_pp,Emins(j),Emaxs(j))
    end do
    Flux = (F_sy + F_br + F_IC + F_pp)/dilution

    !In gamma-rays we care about photon fluxes [cm^-2 s^-1]
    j=5
    F_sy(j) = phflux_int(mnu,E_sy,dE_sy,L_sy,Emins(j),Emaxs(j))   
    F_br(j) = phflux_int(mE_f,E_br,dE_br,L_br,Emins(j),Emaxs(j)) 
    F_IC(j) = phflux_int(mE_f,E_IC,dE_IC,L_IC,Emins(j),Emaxs(j)) 
    F_pp(j) = phflux_int(mE_f,E_pp,dE_pp,L_pp,Emins(j),Emaxs(j))
    phFlux = (F_sy(j) + F_br(j) + F_IC(j) + F_pp(j))/dilution

    j=6
    F_sy(j) = phflux_int(mnu,E_sy,dE_sy,L_sy,Emins(j),Emaxs(j))   
    F_br(j) = phflux_int(mE_f,E_br,dE_br,L_br,Emins(j),Emaxs(j)) 
    F_IC(j) = phflux_int(mE_f,E_IC,dE_IC,L_IC,Emins(j),Emaxs(j)) 
    F_pp(j) = phflux_int(mE_f,E_pp,dE_pp,L_pp,Emins(j),Emaxs(j))
    phFlux2 = (F_sy(j) + F_br(j) + F_IC(j) + F_pp(j))/dilution

    99  format (A29,X,ES10.2)
    write(*,99) 'F_3-10keV  [erg/s/cm2] =', Flux(1)
    !write(*,99) 'F_8-18keV [erg/s/cm2] =', Flux(2)
    write(*,99) 'F_10-30keV [erg/s/cm2] =', Flux(2)
    !write(*,99) 'F_GeV  [1/s/cm2] =', phFlux
    !write(*,99) 'F_TeV  [1/s/cm2] =', phFlux2
    write(*,99) 'F_0.1-100GeV [erg/s/cm2] =', Flux(5)
    write(*,99) 'F_0.1-10TeV [erg/s/cm2] =', Flux(6)

  End subroutine fluxes_calc

  !****************************************************************************************
  !	Integrate the energy fluxes at a given band
  !****************************************************************************************

  Function flux_int(N_int,Eph,dEph,Lum,Eph_min,Eph_max)

    use global
    real(dp) :: flux_int
    integer, intent(in) :: N_int
    real(dp), intent(in) :: Eph(N_int),dEph(N_int),Lum(N_int),Eph_min,Eph_max
    integer :: i_int

    !print *, N_int, Eph_min/eV, Eph_max/eV
    flux_int = 0d0
    do i_int = 1,N_int
      if (Eph(i_int) >= Eph_min .and. Eph(i_int) < Eph_max) then
        flux_int = flux_int + Lum(i_int)*dEph(i_int)
        !flux_int = flux_int + 1.d0 * (dEph(i_int)/eV)
        !print*, i_int
      end if
    end do
    !print*, flux_int, (Eph_max - Eph_min)/eV, flux_int/((Eph_max - Eph_min)/eV)

  End Function flux_int

  
  !****************************************************************************************
  ! Integrate the photon fluxes at a given band
  !****************************************************************************************

  Function phflux_int(N_int,Eph,dEph,Lum,Eph_min,Eph_max)

    use global
    real(dp) :: phflux_int
    integer, intent(in) :: N_int
    real(dp), intent(in) :: Eph(N_int),dEph(N_int),Lum(N_int),Eph_min,Eph_max
    integer :: i_int

    phflux_int = 0d0
    do i_int = 1,N_int
      if (Eph(i_int) >= Eph_min .and. Eph(i_int) < Eph_max) then
        phflux_int = phflux_int + Lum(i_int)*dEph(i_int)/Eph(i_int)
      end if
    end do

  End Function phflux_int

  !----------------------------------------------------------------------------------------

END MODULE fluxes

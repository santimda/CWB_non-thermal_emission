module constants

implicit none
public

integer, parameter :: dp = kind(1.d0)
real(dp), parameter :: pi = acos(-1.d0)
real(dp), parameter :: par = sqrt(2.d0/3.d0)
real(dp), parameter :: deg2rad = pi/180.d0

!===================================================
!Physical constants in CGS units

real(dp), parameter :: eV = 1.602176d-12  ! electron-volt
real(dp), parameter :: c = 2.99792458d10  ! speed of light in vacuum
real(dp), parameter :: k = 1.380658d-16   ! boltzmann's constant 
real(dp), parameter :: qe = 4.8032068d-10 ! electron charge
real(dp), parameter :: mec2 = 0.511d6*eV  ! electron rest energy
real(dp), parameter :: me = mec2 / c**2   ! electron mass
real(dp), parameter :: re = 2.8d-13       ! electron classical radius
real(dp), parameter :: sigma = 5.6705d-5  ! stephan-boltzmann's constant
real(dp), parameter :: sigma_T = 6.652459d-25  ! Thomson cross-section for e-
real(dp), parameter :: h = 6.6260755d-27  ! planck's constant
real(dp), parameter :: h_bar = h/(2.*pi)  ! 
real(dp), parameter :: alfa = 7.29735308d-3 ! fine structure constant
real(dp), parameter :: M_sun_yr = 1.989d33/31556926.d0 ! M_sun/yr
real(dp), parameter :: R_sun = 6.96d10    ! solar radius 
real(dp), parameter :: L_sun = 4.d33      ! solar luminosity (CGS)
real(dp), parameter :: AU = 1.496d13      ! Astronomical Unit
real(dp), parameter :: pc = 3.0856d18     ! parsec
real(dp), parameter :: kpc = 3.0856d21    ! kiloparsec
real(dp), parameter :: mp = 1.673d-24     ! proton mass
real(dp), parameter :: g_ad = 5.d0/3.d0   ! adiabatic coefficient 
real(dp), parameter :: m_pic2 = 134.98d6*eV ! pion rest energy
real(dp), parameter :: mpc2 = 1836.d0*mec2  ! proton rest energy
real(dp), parameter :: T_CMB = 2.73d0             
real(dp), parameter :: T_CMB_ad = k*T_CMB/mec2 
real(dp), parameter :: u_CMB = 10.d0**(-236.d0/19.d0) ! CITE?
real(dp), parameter :: Ry = 2.17987d-11   ! Rydberg unit = 2.18e-18 J

!============================================================================
! Some useful constants for p-p
real(dp), parameter :: TeV = 1.d12*eV  
real(dp), parameter :: E_th = 1.22d9*eV   ! threshold energy for pi^0 production
real(dp), parameter :: k_pp = 0.5d0

!============================================================================

END module constants

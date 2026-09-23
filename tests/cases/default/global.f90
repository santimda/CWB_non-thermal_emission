module global
   use constants
   use system_parameters
   use controls
   implicit none
   public

   !====================================================================
   ! This module contains constants (from constants.f90), variables and 
   ! parameters that all other modules use. 
   ! It also has defined some subroutines that many modules call.
   !====================================================================

   integer, parameter :: nth = 4           ! number of threads (parallelization)
   character(len=*), parameter :: output_dir = 'output/'

   !========================================================
   ! Simulation grid sizes

   integer, parameter :: ml = 200      ! spatial resolution for the CD
   integer, parameter :: ml_2 = ml/2   ! used for the upper side of the CD
   integer, parameter :: mE_e = 150    ! energy resolution for the electron distribution
   integer, parameter :: mE_p = 150    ! energy resolution for the proton distribution
   integer, parameter :: mphi = 4     ! number of linear-emitters (without mirroring)
   real(dp) :: delta_phi = (pi/2.d0)/(mphi-1.d0) ! angular separation between linear-emitters
   real(dp) :: d_phi = (pi/4.d0)/mphi  ! angular width of each linear-emitter or wedge
   integer, parameter :: m_orb = 1     ! number of orbital phases
   integer, parameter :: mnu = 800     ! resolution in photon frequency nu (for synchrotron)
   integer, parameter :: mE_f = 440    ! resolution in photon energy (IC, brem, p-p)
   integer, parameter :: N_bands = 6   ! number of flux bands
   integer, parameter :: mnu_data = 4  ! number of radio bands observed

   integer, parameter :: n_map_freq = 3
   real(dp), parameter :: map_freq(n_map_freq) = (/ 2.3d9, 8.6d9, 230.d9 /)


   !========================================================

   integer :: i,j,l

   ! CD
   real(dp), dimension(:), allocatable :: x, y, z

    ! thermo
   real(dp), dimension(:), allocatable :: P, rho, T, v, B, B_par, n_a, H_sh, dVol
    real(dp), dimension(:), allocatable :: eta_acc

   ! non-thermal particle distributions
   real(dp), parameter :: Emin_e = 2.d0*mec2  
   real(dp), parameter :: Emin_e_inj = 1.5*Emin_e   ! Minimum energy for injected spectrum
   real(dp), parameter :: Estiff_e = 1.d18*eV       ! Possible hardening
   real(dp), parameter :: Emax_e = 1.d14*eV         ! 100 TeV
   real(dp) :: Eint_e,dE_e
   real(dp), dimension(:,:), allocatable :: N_etot  ! N_e(l,E_i)
   real(dp), dimension(:), allocatable :: E_e,N_e   ! N_e(E_i) at a fixed l

   real(dp), parameter :: Emin_p = 1.01d0*mpc2      ! ~1 GeV
   real(dp), parameter :: Emin_p_inj = Emin_p       ! Minimum energy for injected spectrum
   real(dp), parameter :: Estiff_p = 1.d6*mpc2      ! ~can be harder close to E_min
   real(dp), parameter :: Emax_p = 1.d16*eV         ! 10 PeV
   real(dp) :: Eint_p,dE_p
   real(dp), dimension(:,:), allocatable :: N_ptot  ! Np(l,E_i)
   real(dp), dimension(:), allocatable :: E_p,N_p   ! N_p(E_i) at a fixed l

   !========================================================
   ! emission

   real(dp) :: P_Eps, dP_Eps
   real(dp) :: Xst1(3),Xst2(3),Xemi(3),XeE(3)    ! Positions for anisotropic processes

   ! syn
   real(dp), dimension(:,:), allocatable :: Lsy  ! L_sy(l,E_i)
   real(dp), dimension(:), allocatable :: nu, E_sy, Lsy_tot, Lsy_tot_abs 

   ! IC
   real(dp), dimension(:), allocatable :: LIC,LIC_tot,LIC_tot_abs
   real(dp) :: Kappa1,Kappa2 ! = (R/2D)**2
   real(dp), parameter :: cte_Ndot = 2.*re**2*mec2**3/(pi*h_bar**3*c**2)
   logical, parameter :: ani = .True., iso = .False.     ! IC interaction regime
   integer, parameter :: n_ic_lookup = 3000
   real(dp), parameter :: ic_lookup_min = 1.d-7
   real(dp), parameter :: ic_lookup_max = 2.d3
   real(dp), parameter :: log_ic_lookup_min = log(ic_lookup_min)
   real(dp), parameter :: log_ic_lookup_max = log(ic_lookup_max)

   ! br
   real(dp), dimension(:), allocatable :: Lbr,Lbr_tot,Lbr_tot_abs
   logical, parameter :: scr = .True., unscr = .False.   ! Type of nuclei

   ! pp
   real(dp), dimension(:), allocatable :: Lpp,Lpp_tot,Lpp_tot_abs

   ! factors for p-p and rel br. emission
   real(dp) :: z_br, z_pp 

   !========================================================
   ! absorption

   ! g-g
   integer, parameter :: mE_gg = 600, m_ang_gg = 400
   real(dp), dimension(:,:), allocatable :: tau1,tau2
   real(dp), parameter :: Ef_min_gg = 5.d7*eV    ! 50 MeV 
   real(dp), parameter :: Ef_max_gg = 5.d13*eV   ! 50 TeV 
   real(dp) :: Ef_gg(mE_gg),Ef_int_gg

   ! f-f
   integer, parameter :: mE_ff = 800, m_ang_ff = 800, mr_map = 1000
   real(dp), dimension(:,:), allocatable :: tau_ff,tau_ff1,tau_ff2
   real(dp), parameter :: nu_min_ff = 1.d7    ! 10 MHz 
   real(dp), parameter :: nu_max_ff = 1.d12   ! 1 THz
   real(dp) :: Lsy1(ml,mnu),Lsy2(ml,mnu),nu_int_ff

   ! syn
   real(dp), parameter :: nu_min_syn_Hz = 5.d7    ! 50 MHz 
   real(dp), parameter :: nu_max_syn_Hz = 1.d19   ! --> modify to restrict the SED if needed

   !========================================================

   real(dp), parameter :: f_Ny = 0.6d0      ! sampling frequency for maps

   !=========================================================


 contains

   function output_file(name) result(path)
      character(len=*), intent(in) :: name
      character(len=256) :: path

      path = output_dir // name
   end function output_file

   function map_frequency_label(nu_hz) result(label)
      real(dp), intent(in) :: nu_hz
      character(len=16) :: label

      write(label, '(I0)') int(nu_hz/1.d9)
   end function map_frequency_label

    !============================================================================
   !  Cooling_IC

   Function F_ani(u)  ! Used for anisotropic IC (Khangulyan 2014)
      integer, parameter :: dp = kind(1.d0)
      real(dp) :: F_ani
      real(dp), intent(in):: u
      real(dp), parameter:: c_ani = 6.13
      F_ani = (c_ani*u*log(1.+2.16*u/c_ani)) / (1.+c_ani*u/0.822)
   end Function F_ani

   Function F_iso(u)  ! Used for isotropic IC (Khangulyan 2014)
      integer, parameter :: dp = kind(1.d0)
      real(dp) :: F_iso,g_u
      real(dp), intent(in):: u
      real(dp), parameter:: c_iso = 5.68
      g_u = 1./(1. + (-0.362*u**0.682)/(1. + 0.826*u**1.281))
       F_iso = (c_iso*u*log(1.+0.722*u/c_iso)) / (1.+c_iso*u/0.822)*g_u ! Corrected fit
   end Function F_iso

   !   N^dot_IC

   Function F_1(x0)  ! Used for anisotropic IC (Khangulyan 2014)
      integer, parameter :: dp = kind(1.d0)
      real(dp) :: F_1,g
      real(dp), intent(in):: x0
      g = 1./(1. + (0.153*x0**0.857)/(1. + 0.254*x0**1.84))
      F_1 = ( (pi**2/6. + x0)*exp(-x0) )*g   
   end Function F_1

   Function F_2(x0)  ! Used for anisotropic IC (Khangulyan 2014)
      integer, parameter :: dp = kind(1.d0)
      real(dp) :: F_2,g
      real(dp), intent(in):: x0
      g = 1./(1. + (1.33*x0**0.691)/(1. + 0.534*x0**1.668))
      F_2 = ( (pi**2/6. + x0)*exp(-x0) )*g   
   end Function F_2

   Function F_3(x0)  ! Used for isotropic IC (Khangulyan 2014)
      integer, parameter :: dp = kind(1.d0)
      real(dp) :: F_3,g
      real(dp), intent(in):: x0
      real(dp), parameter:: c_3 = 0.319
      g = 1./(1. + (0.443*x0**0.606)/(1. + 0.54*x0**1.481))
      F_3 = (pi**2/6.*(1.+c_3*x0)/(1.+pi**2*c_3*x0/6.)*exp(-x0) )*g   
   end Function F_3

   Function F_4(x0)  ! Used for isotropic IC (Khangulyan 2014)
      integer, parameter :: dp = kind(1.d0)
      real(dp) :: F_4,g
      real(dp), intent(in):: x0
      real(dp), parameter:: c_4 = 6.62
      g = 1./(1. + (0.726*x0**0.461)/(1. + 0.382*x0**1.457))
      F_4 = (pi**2/6.*(1.+c_4*x0)/(1.+pi**2*c_4*x0/6.)*exp(-x0) )*g   
   end Function F_4

   !============================================================================
   !  Black-Body spectrum

   Function nBB(Epho,T)
      real(dp) :: nBB
      real(dp), intent(in) :: Epho,T
      nBB = 2.d0*pi/((h*c)**3) * (Epho**2/(exp(Epho/(k*T))-1.d0))
   end Function nBB

   !============================================================================

     ! Validate shock identifiers used by per-shock routines.

    subroutine validate_shock_number(number,routine)
       character(len=*), intent(in) :: number, routine

       if (number /= '1' .and. number /= '2') then
          write(*,*) 'ERROR: invalid shock number in ',trim(routine), &
             ': expected 1 or 2, received ',trim(number)
          error stop
       endif
    end subroutine validate_shock_number

    !============================================================================
    ! p-p cross section

   Function sigma_pp(L1,Ep)
      real(dp), intent(in) :: L1, Ep
      real(dp) :: sigma_pp
      if (Ep < E_th) then !E_th= 1.22 GeV
         sigma_pp = 0.d0
      else           
         sigma_pp = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27 !cm²
         sigma_pp = sigma_pp*(1.d0-(E_th/Ep)**4)**2 ! Add factor for better fit at low energies 
      endif
   end Function sigma_pp


   !============================================================================


END module global

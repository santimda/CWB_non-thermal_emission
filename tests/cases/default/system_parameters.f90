module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   !----------------------------------------------------------------------------
   ! System parameters
   real(dp), parameter :: distance = 2.4d0            ! in kpc (Callingham2018)
    real(dp), parameter :: date_ref_epoch = 2020.0
    real(dp), parameter :: psi_ref_epoch = 85.d0 * deg2rad
    real(dp), parameter :: D_proj_ref_epoch = 47.d0*distance*AU
    real(dp), parameter :: a_orbit = 47.d0*distance*AU ! placeholder for multi-epoch mode
   real(dp), parameter :: e_orbit = 0.7d0             ! eccentricity (Bloot+2021, Han+2021)
   real(dp), parameter :: w_orbit = -8.d0*deg2rad     ! argument of periastron (Bloot+2021, iso)
   real(dp), parameter :: i_orbit = -21.d0*deg2rad    ! inclination (Bloot+2021, iso)
   real(dp), parameter :: Omega_orbit = 88.d0*deg2rad ! position angle of ascending node (Han+2021)
   real(dp), parameter :: v_orbit = -173.d0*deg2rad   ! true anomaly at 2018 epoch (Han+2021)
   real(dp), parameter :: eta = 0.32d0                ! wind momentum rates ratio (Marcote2020->0.44+-0.08; Zhekov+25->0.2)
   ! Eq.3 in Eichler and Usov 1983 gives asymptotic angle theta:
   real(dp), parameter::  theta_inf = 2.1d0*(1.d0-0.25*eta**0.4d0)*eta**(1.d0/3.d0)
   real(dp), parameter :: filling = 0.3_dp            ! filling factor of the stellar winds (Runacres2002->0.2?)


   !----------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: alpha = 2.45d0       ! spectral index at injection (alpha = -2*alpha_radio + 1, alpha_radio = -0.71, Callingham2018)
   real(dp), parameter :: alpha2 = 2.42d0      ! spectral index at injection for high energies (not really used)
   real(dp), parameter :: K_ep = 0.022         ! Fraction L_NT_e/L_NT


   !----------------------------------------------------------------------------
   ! Free parameters   
   real(dp), parameter :: eta_B1 = 5.3d-3, eta_B2 = eta_B1   ! U_mag = eta_B*P_th
   real(dp), parameter :: frac_NT = 2.0d-1     ! Fraction L_NT/L_inj_perp_norm

   !----------------------------------------------------------------------------
   ! Global variables tied to free parameters
   real(dp) :: frac_B1, frac_B2
   real(dp) :: inj_e, inj_p

   !----------------------------------------------------------------------------
   ! Primary star (WN)
   real(dp), parameter :: T1 = 86.08d3            ! Zhekov+25; Crowther2007->65kK
   real(dp), parameter :: Mdot1 = 2.5d-5*M_sun_yr ! intermediate value (Zhekov+25=4.5e-5, dP+22=3e-5)
   real(dp), parameter :: vinf1 = 3350.d5         ! Zhekov+25->3332, Callingham+20->3500
   real(dp), parameter :: T1mcc = k*T1/mec2
   real(dp), parameter :: Rst1 = 6.d0*R_sun       ! 
   real(dp), parameter :: Lst1 = 4.d0*pi*Rst1**2*sigma*T1**4
   !real(dp), parameter :: Lst1 = 10**5.7 * L_sun
   !real(dp), parameter :: Rst1 = sqrt(Lst1/(4.d0*pi*sigma*T1**4))
   real(dp), parameter :: v_rot1 = 0.1d0*vinf1     ! Ekstrom for M > 85 Mo, White & Chen 1995 use generic 250 km/s
   real(dp), parameter :: gamma_w1 = 1.0d0         ! 
   real(dp), parameter :: r_A_st1 = 1.0d0          ! r_A/R_st
   real(dp), parameter :: T_w1 = 0.3d0*T1          ! Drew1990 for r >> Rst
   ! Relative abundances (for solar composition, X=0.73, Y=0.25, Z=0.02; for a WR, X=0.1-0.2, Y=0.4-0.6, Z=0.2-0.8)
   real(dp), parameter :: X_ab1 = 0.1d0, Y_ab1 = 0.5d0 ! No H was detected in Callingham2018, used a value from Crowther2007
   real(dp), parameter :: Z_ab1 = 1.d0 - X_ab1 - Y_ab1 ! X + Y + Z = 1
   ! The mean atomic weight is = 1/(X + Y/4 + Z/15.5) for a neutral gas; the factor 15.5 is for a solar composition
   real(dp), parameter :: mu_w1 = 2.d0             ! Michaël (the wind is not neutral nor totally ionised)
   real(dp), parameter :: Zq_w1 = 1.0d0            ! Z = rms ionic charge 


   !----------------------------------------------------------------------------
   ! Secondary star (WC)
   real(dp), parameter :: T2 = 60.14d3            ! Zhekov+25; Crowther2007->60kK
   real(dp), parameter :: Rst2 = 6.d0*R_sun       ! 
   real(dp), parameter :: vinf2 = 2350.d5         ! Zhekov+25->2374, Callingham+->2200
   real(dp), parameter :: Mdot2 = Mdot1*(vinf1/vinf2)*eta ! using eta 
   !real(dp), parameter :: Mdot2 = 1.3d-5*M_sun_yr ! Zhekov+25->1.3e-5, delPalacio+->2e-5; trying to keep eta~0.44
   real(dp), parameter :: T2mcc = k*T2/mec2
   real(dp), parameter :: Lst2 = 4.d0*pi*Rst2**2*sigma*T2**4
   real(dp), parameter :: v_rot2 = 0.1d0*vinf2     ! rule-of-thumb
   ! Relative abundances (for solar composition, X=0.73, Y=0.25, Z=0.02; for a WR, X=0.1-0.2, Y=0.4-0.6, Z=0.2-0.8)
   real(dp), parameter :: X_ab2 = 0.001d0, Y_ab2 = 0.4d0 ! No H was detected in Callingham2018
   real(dp), parameter :: Z_ab2 = 1.d0 - X_ab2 - Y_ab2 ! X + Y + Z = 1
   ! The mean atomic weight is = 1/(X + Y/4 + Z/15.5) for a neutral gas; the factor 15.5 is for a solar composition
   real(dp), parameter :: mu_w2 = 4.d0             ! Michaël (the wind is not neutral nor totally ionised)
   real(dp), parameter :: Zq_w2 = 1.005d0          ! Z = rms ionic charge, ????
   real(dp), parameter :: gamma_w2 = 1.01d0        ! Leitherer 1995, ????
   real(dp), parameter :: r_A_st2 = 1.0d0          ! r_A/R_st
   real(dp), parameter :: T_w2 = 0.3d0*T2          ! Drew1990 for r >> Rst

   !----------------------------------------------------------------------------
   ! Wind-collision region
   real(dp), parameter :: mu_1 = 4.d0/(3.d0 + 5.d0*X_ab1 - Z_ab1) ! mean atomic weight = 4/(3+5X-Z) for complete ionization
   real(dp), parameter :: mu_e1 = 2.d0/(1.d0 + X_ab1)             ! mu_e = 2/(1+X) only for Z=0; mu_e > mu! 
   real(dp), parameter :: mu_i1 = 1.d0/(1.d0/mu_1 - 1.d0/mu_e1)  
   real(dp), parameter :: Zq_1 = 1.23_dp        
   real(dp), parameter :: gamma_1 = 1.34_dp
   real(dp), parameter :: mu_2 = 4.d0/(3.d0 + 5.d0*X_ab2 - Z_ab2) ! mean atomic weight = 4/(3+5X-Z) for complete ionization
   real(dp), parameter :: mu_e2 = 2.d0/(1.d0 + X_ab2)             ! mu_e = 2/(1+X) only for Z=0; mu_e > mu! 
   real(dp), parameter :: mu_i2 = 1.d0/(1.d0/mu_2 - 1.d0/mu_e2)   
   real(dp), parameter :: Zq_2 = 1.23_dp        
   real(dp), parameter :: gamma_2 = 1.34_dp
   
   real(dp), parameter :: z_br1 = 3.04           ! Zhekov+2025 WN composition (mean_z.py)
   real(dp), parameter :: z_pp1 = 2.89           ! Zhekov+2025 WN composition (mean_z.py)
   real(dp), parameter :: z_br2 = 12.3           ! Zhekov+2025 WC composition (mean_z.py)
   real(dp), parameter :: z_pp2 = 4.84           ! Zhekov+2025 WC composition (mean_z.py)
   !real(dp), parameter :: z_br = 2.24                 ! Padovani+2018 for Willms+2000 abundances
   !real(dp), parameter :: z_pp = 2.17                 ! Padovani+2018 for Willms+2000 abundances
   

   !============================================================================


END module system_parameters

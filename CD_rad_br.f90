module CD_rad_br
   use global
   use numerical_utils, only: safe_acos_dot, unit_vector_diff, vector_log
   implicit none
   
   public :: CD_rad_br_run, vec_Ebr

   contains

   !============================================================================
   ! CD_rad_br_run  (START)
   !============================================================================
   ! This module calculates the intrinsic bremsstrahlung emission. 
   ! If screened = .True. it is screend, if =.False. it is unscreend 
   ! The program reads as an input the electron energy distribution. 
   !============================================================================

   subroutine CD_rad_br_run(number,Ef,Lbr_sum,Lbr_sum_abs,D,Xobs,screened)
      character, intent(in) :: number
      logical, intent(in) :: screened
      Real(dp), intent(in) :: Ef(mE_f), D,Xobs(3)
      Real(dp), intent(out) :: Lbr_sum(mE_f), Lbr_sum_abs(mE_f)
      Real(dp) :: d1, d2, dil1, dil2
      Real(dp) :: Ef_int, dEf, sum, sum_br, dsigma_Br
      Real(dp) :: gamma_e(mE_e)
      Real(dp) :: ang_min, ang_int, ang1, ang2, cte_sigma, u_br
      Integer :: nl, iang1, iang2, iEf
      interface
         function Phi_Br(u_br,gamma_e,screened)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp) :: Phi_Br
            real(dp), intent(in) :: u_br,gamma_e
            logical, intent(in) :: screened
         end function Phi_Br
      end interface

      call validate_shock_number(number,'CD_rad_br_run')

      allocate( Lbr(mE_f) )
      allocate( N_etot(ml,mE_e), N_e(mE_e) )
      allocate( tau1(mE_gg,m_ang_gg), tau2(mE_gg,m_ang_gg) )

      80  format(A10,X,Es10.2)

      print *, 'Calculate rel. br. emission in S',number

      !------------------------------------------------------------------
      ! Read thermodynamical quantities at the positions x(l),y(l)
      open(11,file=output_file('thermo'//number//'.dat'))
      do l=1,ml
         read(11,*) nl,rho(l),P(l),v(l),T(l),B(l),n_a(l) ! cgs units
         B_par(l) = par*B(l)
      end do
      close(11)

      if (number=='1') then 
         z_br = z_br1
      else
         z_br = z_br2
      end if

      cte_sigma = 4.d0*alfa*re**2*z_br

      ! Electron distribution at position l with energy E_e(i): N_etot(l,i) 
      open (14,file=output_file('dist_e'//number//'.dat'))
      do l=1,ml
         do i=1,mE_e
            read(14,*) j, E_e(i), N_etot(l,i)
            E_e(i) = E_e(i)*eV    !   Energy from eV to erg
         end do
         read(14,*)
      end do
      close(14) 
      Eint_e = E_e(2)/E_e(1)

      ! Matrix with opacity coefficients at a distance of R_star for different angles
      ! between the axis between the emitter and the star, and the LOS; and for 
      ! different energies of the emitted photons: tau(Ef,ang) [one per star]

      Ef_int_gg = (Ef_max_gg/Ef_min_gg)**(1.d0/dfloat(mE_gg-1))

      open (188,file=output_file('opac_gg1.dat'))
      do i=1,mE_gg
         read(188,*) (tau1(i,j), j=1,m_ang_gg)
      end do
      close(188)
      open (188,file=output_file('opac_gg2.dat'))
      do i=1,mE_gg
         read(188,*) (tau2(i,j), j=1,m_ang_gg)
      end do
      close(188)
      ang_int = pi/float(m_ang_gg)
      ang_min = ang_int/2.d0

      Xst1 = (/ 0.d0, 0.d0, 0.d0 /)
      Xst2 = (/ D, 0.d0, 0.d0 /)


      !=======================================================================
      !  Begin rel. br. calculations 
      !=======================================================================

      Lbr_sum_abs = 0.d0
      Lbr_sum = 0.d0

      gamma_e = E_e/mec2
      Ef_int = Ef(2)/Ef(1)
      ! Calculate rel. Brem. emission at each position and then correct for GGA and sum
      sum_br=0.d0
      DO l = 1, ml
         N_e = N_etot(l,:) ! Work with a vector instead of array for simplicity
         Xemi = (/ x(l),y(l),z(l) /) 

         d1 = norm2(Xemi-Xst1)
         call unit_vector_diff(Xemi,Xst1,3,XEe) ! unit vector from 2nd to 1st
         ang1 = safe_acos_dot(Xobs,XEe) ! angle between Xobs and XEe
         iang1 = int(m_ang_gg*abs(ang1-ang_min)/pi)+1
         dil1 = Rst1/d1 

         d2 = norm2(Xemi-Xst2)
         call unit_vector_diff(Xemi,Xst2,3,XEe) ! unit vector from 2nd to 1st
         ang2 = safe_acos_dot(Xobs,XEe) ! angle between Xobs and XEe
         iang2 = int(m_ang_gg*abs(ang2-ang_min)/pi)+1
         dil2 = Rst2/d2

         !------------------------------------------------------------------------

         do i = 1,mE_f           ! Ef(i) 
            dEf = Ef(i)*(Ef_int-1.d0)
            sum = 0.d0
            do j = 1,mE_e            ! N_e(j),gamma_e(j)
               u_br = Ef(i)/E_e(j) ! energy ratio
               if (u_br > 1.d0) then   ! Ef <= Ee from energy conservation
                  dsigma_Br = 0.d0  ![cm**2/m_e c2]
               else
                  dsigma_Br= cte_sigma*Phi_Br(u_br,gamma_e(j),screened) !/Ef(i), cancels later
                  if (dsigma_Br > 0.d0) then !Phi_Br can yield negative values in bad E ranges
                     P_eps = n_a(l)*c*dsigma_Br ! *Ef(i), cancels in dsigma_Br; [P_eps] = 1/seg
                     dE_e = E_e(j)*(Eint_e-1.d0)
                     dP_Eps = N_e(j)*P_eps*dE_e
                     sum = dP_Eps + sum ! P = int dP_Eps
                  endif
               endif
            end do ! j
            ! Transform photon emissivity to specific luminosity
            Lbr(i) = Ef(i)*sum  ! [erg/seg/erg]
            Lbr_sum(i) = Lbr_sum(i) + Lbr(i)

            if (Ef(i) > Ef_min_gg .and. Ef(i) < Ef_max_gg) then
               iEf = int(log(Ef(i)/Ef_min_gg)/log(Ef_int_gg)) + 1
               Lbr_sum_abs(i) = Lbr_sum_abs(i) + Lbr(i)&
               *exp( -dil1*tau1(iEf,iang1) - dil2*tau2(iEf,iang2) )
            else
               Lbr_sum_abs(i) = Lbr_sum_abs(i) + Lbr(i)
            end if

            sum_br = sum_br + Lbr(i)*dEf

         end do ! i

      END DO ! l

      write(*,80) 'L_Brem=',sum_br

      !------------------------------------------------------------

      open(16,file=output_file('rad_br'//number//'.dat'))
      do i=1,mE_f
         write(16,*) Ef(i)/eV,Ef(i)*Lbr_sum(i),Ef(i)*Lbr_sum_abs(i)
      end do
      close(16)

      deallocate( N_etot,N_e )
      deallocate( Lbr )
      deallocate( tau1,tau2 )

   END subroutine CD_rad_br_run

   !============================================================================
   ! CD_rad_br_run   (END)
   !============================================================================


   !============================================================================
   ! vec_Ebr
   !============================================================================

   subroutine vec_Ebr(Ebr)
      real(dp), intent(out) :: Ebr(mE_f)
      real(dp) :: Eaux, Ee_max, Ebr_min, Ebr_max, Ebr_int, dummy
 
      !--------------------------------------------
      ! Obtain Ee_max (this is the cutoff energy, not the max energy)
      Ee_max = 0.d0
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
      Ebr_max = 10.d0 * Ee_max * eV ! Ebr_max ~ Ee_max

      Ebr_min = 1.d4*eV ! Rel. br. is for photon energies >> 1 keV
      
      !--------------------------------------------------
      ! Build Ebr vector
      call vector_log(Ebr_min,Ebr_max,Ebr_int,Ebr)

   END subroutine vec_Ebr

   !================================================================

END module CD_rad_br


!========================================

Function Phi_Br(u,g_e,screened)
   ! Calculates Phi function in rel.Br. cross-section assuming a factor z_br given in global
   ! It can give negative values for un-screened nuclei in certain E ranges.
   use global  
   implicit none
   real(dp) :: Phi_Br
   real(dp), intent(in) :: u,g_e
   logical, intent(in) :: screened

   if (screened) then      ! screened nuclei
      Phi_Br = ( 1. + (1.-u)**2 - (2./3.d0)*(1.-u)) *log(191.d0/z_br**(1./3.d0)) + (1.-u)/9.d0
   else                    ! un-screened nuclei 
      Phi_Br = ( 1. + (1.-u)**2 - (2./3.d0)*(1.-u)) *( log(2.*g_e*(1./u-1.)) -0.5d0 )
   endif

   return
End function Phi_Br

!========================================

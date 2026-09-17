module CD_rad_ic
   use global
   implicit none
   
   private
   public :: CD_rad_ic_run, vec_Eic, initialize_ic_lookup, finalize_ic_lookup
   public :: lookup_ic_functions

   real(dp), allocatable :: ic_lookup_f1(:), ic_lookup_f2(:)

   contains

   !============================================================================
   ! CD_rad_ic_run   (START)
   !============================================================================
   ! This module calculates the intrinsic and absorbed IC emission using as an 
   ! input the electron distribution in diste(number).dat. The integrated 
   ! luminosity (unabs and abs) for the linear emitter is stored in radic.dat.
   !============================================================================

   subroutine CD_rad_ic_run(number,Ef,LIC_sum,LIC_sum_abs,D,Xobs)
      character, intent(in) :: number
      real(dp), intent(in) :: Ef(mE_f),D,Xobs(3)
      real(dp), intent(out) :: LIC_sum(mE_f),LIC_sum_abs(mE_f)
      real(dp) :: Ef_int,ang_min,ang_int,ang1,ang2,d1,d2
      real(dp) :: gamma_e(mE_e),Eps(mE_f),NedEe(mE_e),Eps_int,dEe
      real(dp) :: L_unabs,L_unabs_tot,L_abs_tot,dil1,dil2
      integer :: iang1,iang2,iEf,nl
      call validate_shock_number(number,'CD_rad_ic_run')
      allocate( LIC(mE_f) )
      allocate( N_etot(ml,mE_e),N_e(mE_e) )
      allocate( tau1(mE_gg,m_ang_gg),tau2(mE_gg,m_ang_gg) )

      80  format(A12,X,Es10.2)

      LIC_sum_abs = 0.d0
      LIC_sum = 0.d0
      LIC = 0.d0 

      print *, 'Calculate IC emission in S',number

      !------------------------------------------------------------------
      ! Read positions (x(l),y(l)) and thermodynamical quantities at them

      open(11,file=output_file('thermo'//number//'.dat'))
      do l=1,ml
         read(11,*) nl,rho(l),P(l),v(l),T(l),B(l),n_a(l) ! cgs units
      end do
      close(11)

      ! Electron distribution at position l with energy E_e(i): N_etot(l,i) 
      open (14,file=output_file('dist_e'//number//'.dat'))
      do l=1,ml
         do i=1,mE_e
            read(14,*) j, E_e(i),N_etot(l,i)
         end do
         read(14,*)
      end do
      close(14) 
      E_e = E_e*eV         ! Energy from eV to erg
      Eint_e = E_e(2)/E_e(1)
      gamma_e = E_e/mec2   ! Pre-compute to speed-up IC calculation
      Eps = Ef/mec2
      Eps_int = Eps(2)/Eps(1)

      !-----------------------------------------------------------------
      ! Quantities required for posterior IC and gamma-gamma absorption

      Ef_int = Ef(2)/Ef(1)

      call vector_log(Ef_min_gg,Ef_max_gg,Ef_int_gg,Ef_gg)

      ! Matrix with opacity coefficients at a distance of R_star for different angles
      ! between the axis between the emitter and the star, and the LOS; and for 
      ! different energies of the emitted photons: tau(Ef,�ngulo) [one per star]

      open (188,file=output_file('opac_gg1.dat'))
      open (189,file=output_file('opac_gg2.dat'))
      do i=1,mE_gg
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
      !  Begin IC calculations 
      !=======================================================================

      ! Calculate IC emission at each position and then correct for GGA and sum

      DO l=1,ml
         N_e = N_etot(l,:)          ! Work with a vector instead of an array
         do j = 1,mE_e
            dEe = E_e(j)*(Eint_e-1.d0)
            NedEe(j) = N_e(j)*dEe   ! Pre-compute to speed-up the IC calculation 
         end do
         Xemi = (/ x(l),y(l),z(l) /)

         ! Calculate IC-star1, IC-star2, add and correct for GGA

         !------------------------------------------------------------------------
         ! Star 1
         d1 = sqrt(Xemi(1)**2+Xemi(2)**2+Xemi(3)**2)
         call unit_vector_diff(Xemi,Xst1,3,XEe)  ! unit vector from 2nd to 1st
         ang1 = acos(DOT_PRODUCT(Xobs,XEe)) ! angle between Xobs and XEe
         iang1 = int(m_ang_gg*abs(ang1-ang_min)/pi)+1
         Kappa1 = ( Rst1/(2.d0*d1) )**2
         dil1 = Rst1/d1 

         !------------------------------------------------------------------------
         ! Star 2
         d2 = sqrt((D-Xemi(1))**2+Xemi(2)**2+Xemi(3)**2) 
         call unit_vector_diff(Xemi,Xst2,3,XEe)  ! unit vector from 2nd to 1st
         ang2 = acos(DOT_PRODUCT(Xobs,XEe)) ! angle between Xobs and XEe
         iang2 = int(m_ang_gg*abs(ang2-ang_min)/pi)+1
         Kappa2 = ( Rst2/(2.d0*d2) )**2
         dil2 = Rst2/d2
         !------------------------------------------------------------------------

         call IC(gamma_e,Eps,Eps_int,NedEe,Kappa1,T1,ang1,LIC,L_unabs) ! St1
         do i = 1,mE_f
            LIC_sum(i) = LIC_sum(i)+LIC(i)
            if (Ef_min_gg < Ef(i) .and. Ef(i) < Ef_max_gg) then
               iEf = int(log(Ef(i)/Ef_min_gg)/log(Ef_int_gg)) + 1
               LIC_sum_abs(i) = LIC_sum_abs(i) + LIC(i)&
               *exp( -dil1*tau1(iEf,iang1) - dil2*tau2(iEf,iang2) )
            else
               LIC_sum_abs(i) = LIC_sum(i)
            end if
         end do

         call IC(gamma_e,Eps,Eps_int,NedEe,Kappa2,T2,ang2,LIC,L_unabs) ! St2
         do i = 1,mE_f
            LIC_sum(i) = LIC_sum(i)+LIC(i)
            if (Ef_min_gg < Ef(i) .and. Ef(i) < Ef_max_gg) then
               iEf = int(log(Ef(i)/Ef_min_gg)/log(Ef_int_gg)) + 1
               LIC_sum_abs(i) = LIC_sum_abs(i) + LIC(i)&
               *exp( -dil1*tau1(iEf,iang1) - dil2*tau2(iEf,iang2) )
            else
               LIC_sum_abs(i) = LIC_sum(i)
            end if
         end do

         !if (mod(l,100)==0) print *, 'IC, l=',l, '/', ml

      end do

      !------------------------------------------------------------

      L_abs_tot = 0.d0  
      L_unabs_tot = 0.d0
      open(16,file=output_file('rad_ic'//number//'.dat'))
      do i = 1,mE_f
         write(16,*) Ef(i)/eV,Ef(i)*LIC_sum(i),Ef(i)*LIC_sum_abs(i)
         L_unabs_tot = L_unabs_tot + LIC_sum(i)*Ef(i)*(Ef_int-1.d0)
         L_abs_tot = L_abs_tot + LIC_sum_abs(i)*Ef(i)*(Ef_int-1.d0)
      end do
      close(16)
      write(*,80) 'L_IC(un)=', L_unabs_tot
      write(*,80) 'L_IC(abs)=', L_abs_tot

      deallocate( LIC )
      deallocate( tau1,tau2 )
      deallocate( N_etot,N_e )

   END subroutine CD_rad_ic_run

   !============================================================================
   ! CD_rad_ic_run   (END)
   !============================================================================

   !============================================================================
   ! vec_Eic
   !============================================================================

   subroutine vec_Eic(Eic)
      real(dp), intent(out) :: Eic(mE_f)
      real(dp) :: dummy, Eaux, Ee_max, Eic_min, Eic_max, Eic_int

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
      Eic_max = 5.d0 * Ee_max * eV ! Eic_max ~ Ee_max

      Eic_min = 0.1d0 * k * max(T1,T2)

      !--------------------------------------------------
      ! Build Eic vector
      call vector_log(Eic_min, Eic_max, Eic_int, Eic)

   END subroutine vec_Eic

   !============================================================================

   subroutine initialize_ic_lookup()
      integer :: i_lookup
      real(dp), allocatable :: xlookup(:)
      real(dp) :: xlookup_int

      if (allocated(ic_lookup_f1)) return

      allocate(ic_lookup_f1(n_ic_lookup), ic_lookup_f2(n_ic_lookup))
      allocate(xlookup(n_ic_lookup))
      call vector_log(ic_lookup_min,ic_lookup_max,xlookup_int,xlookup)
      do i_lookup = 1,n_ic_lookup
         ic_lookup_f1(i_lookup) = F_1(xlookup(i_lookup))
         ic_lookup_f2(i_lookup) = F_2(xlookup(i_lookup))
      end do
      deallocate(xlookup)
   end subroutine initialize_ic_lookup

   subroutine finalize_ic_lookup()
      if (allocated(ic_lookup_f1)) deallocate(ic_lookup_f1,ic_lookup_f2)
   end subroutine finalize_ic_lookup

   subroutine lookup_ic_functions(x0, f1_value, f2_value)
      real(dp), intent(in) :: x0
      real(dp), intent(out) :: f1_value, f2_value
      real(dp) :: log_x0
      integer :: i_lookup

      if (x0 <= ic_lookup_min) then
         i_lookup = 1
      else if (x0 >= ic_lookup_max) then
         i_lookup = n_ic_lookup
      else
         log_x0 = log(x0)
         i_lookup = nint((n_ic_lookup - 1.d0) * &
                     (log_x0 - log_ic_lookup_min) / &
                     (log_ic_lookup_max - log_ic_lookup_min)) + 1
      end if

      f1_value = ic_lookup_f1(i_lookup)
      f2_value = ic_lookup_f2(i_lookup)
   end subroutine lookup_ic_functions

   !================================================================

END module CD_rad_ic


!============================================================================
! SUBROUTINE IC
!============================================================================

SUBROUTINE IC(gamma_e,Eps,Eps_int,NedEe,Kappa,Tst,ang,L_IC,L_unabs)
   ! Calculate anisotropic IC emission from interactions with BB photon field,
   ! using the formalism of Khangulyan+ 2014 and precomputed energy-grid terms.
   use global
   use CD_rad_ic, only: lookup_ic_functions
   implicit none
   real(dp), intent(in) :: gamma_e(mE_e),Eps(mE_f),Eps_int,NedEe(mE_e)
   real(dp), intent(in) :: Kappa,Tst,ang
   real(dp), intent(out) :: L_IC(mE_f),L_unabs
   real(dp) :: Tmcc,Tmcc_ang,Tmcc_ang2,x0_K,z_ic,t_ic
   real(dp) :: F1_value,F2_value
   real(dp) :: dEps,Ndot,sum,cte_Ndot_ani
   real(dp) :: P_eps_local,dP_eps_local
   integer :: i_ic,j_ic
   Tmcc = k*Tst/mec2
   Tmcc_ang = Tmcc*(1.-cos(ang))
   Tmcc_ang2 = Tmcc_ang*2.d0
   cte_Ndot_ani = cte_Ndot*Tmcc**2*Kappa

   L_unabs = 0.d0
   !$omp parallel do default(none) reduction(+:L_unabs) schedule(static) &
   !$omp& shared(Eps,gamma_e,Eps_int,NedEe,cte_Ndot_ani,Tmcc_ang2,L_IC) &
   !$omp& private(i_ic,j_ic,dEps,sum,z_ic,t_ic,x0_K,Ndot,F1_value,F2_value,P_eps_local,dP_eps_local)
   do i_ic=1,mE_f           ! Ef (emitted photons)
      dEps=Eps(i_ic)*(Eps_int-1.d0)
      sum=0.d0

      do j_ic=1,mE_e        ! Ee
         z_ic = Eps(i_ic)/gamma_e(j_ic)
         if (z_ic < 1.d0) then
            t_ic = gamma_e(j_ic)*Tmcc_ang2
            x0_K = z_ic/( (1.-z_ic)*t_ic )
            call lookup_ic_functions(x0_K,F1_value,F2_value)
            Ndot = (cte_Ndot_ani/gamma_e(j_ic)**2) * &
                   ( (z_ic**2/(2.*(1.-z_ic)))*F1_value + F2_value )
            P_eps_local = Eps(i_ic) * Ndot   ! N^dot = nf*c*dsigmaIC
            dP_eps_local = P_eps_local * NedEe(j_ic)  ! P_eps*Ne(Ee)*dEe
            sum = dP_eps_local + sum
         end if
      end do

      L_IC(i_ic) = sum                 ! Especific luminosity (plot E(j)*LIC(Ej))
      L_unabs = L_unabs + sum*dEps*mec2
    end do
   !$omp end parallel do

end SUBROUTINE IC


!================================================================

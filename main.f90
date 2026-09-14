PROGRAM main
   use global
   use controls
   implicit none
   real(dp) :: nu_data(mnu_data) ! Frequencies at which we have observations

   interface
      subroutine radio(nu_data)
      use global 
         real(dp), intent(in) :: nu_data(mnu_data)
      end subroutine radio
   end interface

   !----------------------------------------------------------------
   ! Program to calculate the broadband SED and emission maps

   if (sed_all_calc) call sed_all

   !----------------------------------------------------------------
   ! Short version to calculate specific radio fluxes
    
   if (radio_calc) then      
      nu_data = (/1.4d9, 2.2d9, 5.5d9, 9.d9/) 
      call radio(nu_data)
   endif
   
END PROGRAM main

!==================================================================

subroutine radio(nu_data)
   use global
   use controls
   use initialize
   use orbit, only: orbit_run
   use CD, only: CD_run
   use CD_rot, only: CD_rot_run
   use CD_thermo, only: CD_thermo_run
   use CD_dist_e, only: CD_dist_e_run
   use CD_rad_syn, only: CD_rad_syn_data_run
   use opac_ff, only: opac_ff_data_run
   use CD_abssyn, only: CD_abssyn_data_run

   implicit none
   real(dp), intent(in) :: nu_data(mnu_data)
   integer :: i_phi,i_orb
   real(dp) :: phi
   real(dp) :: phase(m_orb), D(m_orb), D_proj(m_orb), psi(m_orb)
   real(dp) :: S_phi_unabs(mnu_data), S_phi_abs(mnu_data)
   real(dp) :: S_WCR_unabs(mnu_data), S_WCR_abs(mnu_data)

   print *, '\n PROGRAM STARTED \n'

   ! Initialize vectors
   call initialize_allocate()

   ! Calculate orbit
   call orbit_run(phase,D,D_proj,psi,single_epoch)

   ! Run for different orbital phases
   open(120,file=output_file('S_data_unabs.dat'))
   open(121,file=output_file('S_data_abs.dat'))
   DO i_orb = 1,m_orb
      print *, 'i(orb)=',i_orb,'/',m_orb
      call CD_run(D(i_orb))
      call CD_thermo_run('1',D(i_orb))
      call CD_thermo_run('2',D(i_orb))
      call CD_dist_e_run('1',D(i_orb))
      call CD_dist_e_run('2',D(i_orb))
      call opac_ff_data_run(nu_data)
      call CD_rad_syn_data_run('1', nu_data)
      call CD_rad_syn_data_run('2', nu_data)

      S_WCR_abs = 0.d0
      S_WCR_unabs = 0.d0
      phi = 0.d0     ! phi = Azimuthal angle 
      do i_phi = 1,mphi
         call CD_rot_run(phi)
         call CD_abssyn_data_run(nu_data, S_phi_unabs, S_phi_abs, D(i_orb), psi(i_orb))
         S_WCR_unabs = S_WCR_unabs + S_phi_unabs
         S_WCR_abs = S_WCR_abs + S_phi_abs
         phi = phi + delta_phi
      end do
      write(120,*) phase(i_orb), (S_WCR_unabs(j), j=1,mnu_data)
      write(121,*) phase(i_orb), (S_WCR_abs(j), j=1,mnu_data)

   ENDDO
   close(120)  ! fluxes_unabs
   close(121)  ! fluxes_abs

   if (verbose) print *, '\n Deallocate vectors'    
   call initialize_deallocate

   print *, '\n PROGRAM ENDED \n'
   
end subroutine radio


!=================================================================


subroutine sed_all()
   use global
   use controls
   use initialize
   use BB, only: BB_run
   use WCR_thermal, only: WCR_thermal_run
   use orbit, only: orbit_run
   use CD, only: CD_run
   use CD_rot, only: CD_rot_run
   use CD_thermo, only: CD_thermo_run
   use CD_dist_e, only: CD_dist_e_run
   use CD_rad_syn, only: CD_rad_syn_run, vec_nu
   use opac_ff, only: opac_ff_run
   use opacity_gg, only: opacity_gg_run
   use CD_rad_ic, only: CD_rad_ic_run, vec_Eic
   use CD_rad_br, only: CD_rad_br_run, vec_Ebr
   use CD_abssyn, only: CD_abssyn_run
   use CD_dist_p, only: CD_dist_p_run
   use CD_rad_pp, only: CD_rad_pp_run, vec_Epp
   use WCR_conv, only: WCR_conv_run
   use fluxes, only: fluxes_calc

   implicit none
   integer :: i_phi,i_orb,i_map
   real(dp) :: phase(m_orb),D(m_orb),D_proj(m_orb),psi(m_orb)
   real(dp) :: phi !,dummy
   real(dp) :: E_IC(mE_f),LIC_WCR_unabs(mE_f),LIC_WCR_abs(mE_f)
   real(dp) :: Lsy_WCR_unabs(mnu),Lsy_WCR_abs(mnu)
   real(dp) :: E_br(mE_f),Lbr_WCR_unabs(mE_f),Lbr_WCR_abs(mE_f)
   real(dp) :: E_pp(mE_f),Lpp_WCR_unabs(mE_f),Lpp_WCR_abs(mE_f)
   real(dp) :: u_proj,v_proj,emi,emi_abs
   real(dp) :: flux(N_bands),Xobs(3)
   integer :: map_input_unit(n_map_freq),map_output_unit(n_map_freq)
   character(len=16) :: freq_label

   ! character(len=300) :: dir_data 
   ! character(len=100) :: prefix_dir_data, cmd, str_test


   print *, '\n ============================'
   print *, ' START'
   print *, '============================\n '

   if (absgg) call opacity_gg_run()  ! same table for all orbital phases

   print *, '\n Initialize vectors'    
   call initialize_allocate()

   ! Consistency check of the wind momentum rates ratio:
   print*, '\n eta =', (Mdot2*vinf2)/(Mdot1*vinf1)

   !===========================================================================
   ! Run for different orbital phases:

   call orbit_run(phase,D,D_proj,psi,single_epoch)
   open(120,file=output_file('fluxes_unabs.dat'))
   open(121,file=output_file('fluxes_abs.dat'))

   DO i_orb = 1,m_orb

   LIC_WCR_abs = 0.d0
   LIC_WCR_unabs = 0.d0 
   Lbr_WCR_abs = 0.d0
   Lbr_WCR_unabs = 0.d0
   Lsy_WCR_abs = 0.d0
   Lsy_WCR_unabs = 0.d0
   Lpp_WCR_abs = 0.d0
   Lpp_WCR_unabs = 0.d0

   if (BBrad) call BB_run()
   if (CD_calc) call CD_run(D(i_orb))
   if (thermo) call CD_thermo_run('1',D(i_orb))
   if (thermo) call CD_thermo_run('2',D(i_orb))
   if (Xrad) call WCR_thermal_run()

   Xobs = (/ cos(psi(i_orb)),sin(psi(i_orb)),0.d0 /) ! Used by various rad processes

   !===============================================================!
   !                                                               ! 
   !                      Electrons                                !
   !                                                               !
   !===============================================================!

   !  open(100,file='Gamma'//gammaj//'/SEDsyn_orb_i'//trim(adjustl(incli))//'.dat')

   IF (electrons) then
       print *, '\n ELECTRONS:'

      if (dist1e) call CD_dist_e_run('1',D(i_orb))
      if (dist2e) call CD_dist_e_run('2',D(i_orb))

      call vec_nu()!E_sy) ! frequencies at which sync. emission is calculated
      if (absff) call opac_ff_run(D(i_orb))

      if (rad_syn1) call CD_rad_syn_run('1')
      if (rad_syn2) call CD_rad_syn_run('2')

      call vec_Eic(E_IC) ! energies at which IC emission is calculated
      call vec_Ebr(E_br) ! energies at which br emission is calculated

      if (maps) then
         do i_map = 1, n_map_freq
            freq_label = map_frequency_label(map_freq(i_map))
            map_input_unit(i_map) = 54 + i_map
            map_output_unit(i_map) = 587 + i_map
            open(map_output_unit(i_map), &
               file=output_file('WCR_radiomap_'//trim(freq_label)//'GHz.dat'))
         end do
      end if

      ! Anisotropic IC emission, and FFA and g-g abs depend on the orientation
      phi = 0.d0     ! phi = Azimuthal angle 
      DO i_phi = 1,mphi
         print *, 'i(phi)=',i_phi,'/',mphi

         if (rotate) call CD_rot_run(phi)

         if (rad_ic1) then
            allocate( LIC_tot(mE_f),LIC_tot_abs(mE_f) )
            call CD_rad_ic_run('1',E_IC,LIC_tot,LIC_tot_abs,D(i_orb),Xobs)
            LIC_WCR_unabs = LIC_WCR_unabs + LIC_tot
            LIC_WCR_abs = LIC_WCR_abs + LIC_tot_abs
            deallocate( LIC_tot,LIC_tot_abs )
         end if

         if (rad_ic2) then
            allocate( LIC_tot(mE_f),LIC_tot_abs(mE_f) )
            call CD_rad_ic_run('2',E_IC,LIC_tot,LIC_tot_abs,D(i_orb),Xobs)
            LIC_WCR_unabs = LIC_WCR_unabs + LIC_tot
            LIC_WCR_abs = LIC_WCR_abs + LIC_tot_abs
            deallocate( LIC_tot,LIC_tot_abs )
         end if

         if (rad_br1) then
            allocate( Lbr_tot(mE_f),Lbr_tot_abs(mE_f) )
            call CD_rad_br_run('1',E_br,Lbr_tot,Lbr_tot_abs,D(i_orb),Xobs,unscr)
            Lbr_WCR_unabs = Lbr_WCR_unabs + Lbr_tot
            Lbr_WCR_abs = Lbr_WCR_abs + Lbr_tot_abs
            deallocate( Lbr_tot,Lbr_tot_abs )
         end if

         if (rad_br2) then
            allocate( Lbr_tot(mE_f),Lbr_tot_abs(mE_f) )
            call CD_rad_br_run('2',E_br,Lbr_tot,Lbr_tot_abs,D(i_orb),Xobs,unscr)
            Lbr_WCR_unabs = Lbr_WCR_unabs + Lbr_tot
            Lbr_WCR_abs = Lbr_WCR_abs + Lbr_tot_abs
            deallocate( Lbr_tot,Lbr_tot_abs )
         end if

         if (abssyn) then
            allocate( Lsy_tot(mnu),Lsy_tot_abs(mnu) )
            call CD_abssyn_run(Lsy_tot,Lsy_tot_abs,D(i_orb),psi(i_orb))
            Lsy_WCR_unabs = Lsy_WCR_unabs + Lsy_tot
            Lsy_WCR_abs = Lsy_WCR_abs + Lsy_tot_abs
            deallocate( Lsy_tot,Lsy_tot_abs )
         end if

         if (maps) then           
            print *, 'Save synthetic radio maps'
            do i_map = 1, n_map_freq
               freq_label = map_frequency_label(map_freq(i_map))
               open(map_input_unit(i_map), &
                  file=output_file('CD_radiomap_'//trim(freq_label)//'GHz.dat'))
            end do
            do i=1,ml
               do i_map = 1, n_map_freq
                  read(map_input_unit(i_map),*) u_proj,v_proj,emi,emi_abs
                  write(map_output_unit(i_map),*) u_proj,v_proj,emi,emi_abs
                  if (i_phi.ne.1 .and. i_phi.ne.mphi) then
                     write(map_output_unit(i_map),*) &
                           u_proj,-v_proj,emi,emi_abs ! The symmetric one
                  end if
               end do
            end do
            do i_map = 1, n_map_freq
               close(map_input_unit(i_map))
            end do
         end if

         phi = phi + delta_phi

      END DO

      if (maps) then
         do i_map = 1, n_map_freq
            close(map_output_unit(i_map))
         end do
      end if

      ! Save total luminosities multiplying by a factor 2 due to axial symmetry
      if (write_e) then
         if (rad_ic1.and.rad_ic2) then
            open (50,file=output_file('rad_ic_WCR.dat'))
            do i=1,mE_f
               write(50,*) E_IC(i)/eV,2.d0*LIC_WCR_unabs(i),2.d0*LIC_WCR_abs(i)
            end do
            close(50)
         end if

         if (abssyn) then
            open (51,file=output_file('rad_syn_WCR.dat'))
            do i=1,mnu
               write(51,*) E_sy(i)/eV,2.d0*Lsy_WCR_unabs(i),2.d0*Lsy_WCR_abs(i)
            end do
            close(51)
         end if

          if (rad_br1.and.rad_br2) then
            open (52,file=output_file('rad_br_WCR.dat'))
            do i=1,mE_f
               write(52,*) E_br(i)/eV,2.d0*Lbr_WCR_unabs(i),2.d0*Lbr_WCR_abs(i)
            end do
            close(52)
         end if
      end if

   END IF


   !===============================================================!
   !                                                               ! 
   !                      Protons                                  !
   !                                                               !
   !===============================================================!

   IF (protons) then

      print *, '\n PROTONS:'

      if (dist1p) call CD_dist_p_run('1',D(i_orb))
      if (dist2p) call CD_dist_p_run('2',D(i_orb))

      call vec_Epp(E_pp) ! energies at which pp emission is calculated

      ! The only angle dependence comes from g-g abs
      phi = 0.d0           ! Azimuthal angle phi
      DO i_phi = 1,mphi
         print *, 'i(phi)=',i_phi,'/',mphi

         if (rotate) call CD_rot_run(phi)

         if (rad_pp1) then
            allocate( Lpp_tot(mE_f),Lpp_tot_abs(mE_f) )
            call CD_rad_pp_run('1',E_pp,Lpp_tot,Lpp_tot_abs,D(i_orb),Xobs)
            Lpp_WCR_unabs = Lpp_WCR_unabs + Lpp_tot
            Lpp_WCR_abs = Lpp_WCR_abs + Lpp_tot_abs
            deallocate( Lpp_tot,Lpp_tot_abs )
         end if

         if (rad_pp2) then
            allocate( Lpp_tot(mE_f),Lpp_tot_abs(mE_f) )
            call CD_rad_pp_run('2',E_pp,Lpp_tot,Lpp_tot_abs,D(i_orb),Xobs)
            Lpp_WCR_unabs = Lpp_WCR_unabs + Lpp_tot
            Lpp_WCR_abs = Lpp_WCR_abs + Lpp_tot_abs
            deallocate( Lpp_tot,Lpp_tot_abs )
         end if

      end do

      ! Save total luminosities multiplying by a factor 2 due to axial symmetry
      if (write_p) then
         open (60,file=output_file('rad_pp_WCR.dat'))
         do i=1,mE_f
            write(60,*) E_pp(i)/eV,2.d0*Lpp_WCR_unabs(i),2.d0*Lpp_WCR_abs(i)
         end do
         close(60)
      end if

   END IF


   !Calculate integrated fluxes in different bands
   if (write_e.and.write_p) then
      print *, '\n Unabsorbed fluxes'
      !Add a factor two due to axial symmetry
      call fluxes_calc(E_sy,2.d0*Lsy_WCR_unabs, &
                        E_br,2.d0*Lbr_WCR_unabs, &
                        E_IC,2.d0*LIC_WCR_unabs, &
                        E_pp,2.d0*Lpp_WCR_unabs, flux)
      write(120,*) phase(i_orb),D(i_orb),(flux(j), j=1,N_bands)

      print *, '\n Absorption-corrected fluxes'
      call fluxes_calc(E_sy,2.d0*Lsy_WCR_abs, &
                        E_br,2.d0*Lbr_WCR_abs, &
                        E_IC,2.d0*LIC_WCR_abs, &
                        E_pp,2.d0*Lpp_WCR_abs,flux)   
      write(121,*) phase(i_orb),D(i_orb),(flux(j), j=1,N_bands)
   endif

   !================================================================


   if (convolve) then
      print *, '\n Convolution the emission map with the instrumental response'
      do i_map = 1, n_map_freq
         call WCR_conv_run(D_proj(i_orb),psi(i_orb),map_freq(i_map))
      end do
   end if

   ENDDO  ! Orbital phase

   close(120)  ! fluxes_unabs
   close(121)  ! fluxes_abs

   print *, '\n Deallocate vectors'    
   call initialize_deallocate

   print *, '\n PROGRAM FINISHED \n'
   
end subroutine sed_all

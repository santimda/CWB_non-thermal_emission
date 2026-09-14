module orbit
   use global
   implicit none

   public :: orbit_run

contains

   !============================================================================
   ! orbit_run
   !============================================================================
   ! Select either a single reference epoch or the placeholder multi-epoch
   ! calculation. The single-epoch mode requires one orbital point.
   !============================================================================

   subroutine orbit_run(phase,D,D_proj,psi,single_epoch)
      real(dp), intent(out) :: phase(m_orb), D(m_orb), D_proj(m_orb), psi(m_orb)
      logical, intent(in) :: single_epoch

      if (single_epoch .and. m_orb > 1) then
         error stop 'single_epoch requires m_orb = 1'
      endif

      if (single_epoch) then
         call single_epoch_orbit(phase,D,D_proj,psi)
      else
         call multi_epoch_orbit(phase,D,D_proj,psi)
      endif

   contains

      !----------------------------------------------------------------------------

      subroutine single_epoch_orbit(phase,D,D_proj,psi)
         real(dp), intent(out) :: phase(1), D(1), D_proj(1), psi(1)
         real(dp) :: date(1)

         if (verbose) print *, 'Use orbital reference epoch'

         phase(1) = 0.d0
         psi(1) = pi - psi_ref_epoch
         D_proj(1) = D_proj_ref_epoch/(distance*AU)
         D(1) = D_proj_ref_epoch/abs(sin(psi_ref_epoch))
         date(1) = date_ref_epoch

         open(11,file=output_file('orbit.dat'))
         write(11,*) phase(1), date(1), D(1), D_proj(1), psi(1)
         close(11)
      end subroutine single_epoch_orbit

      !----------------------------------------------------------------------------

      subroutine multi_epoch_orbit(phase,D,D_proj,psi)
         real(dp), intent(out) :: phase(m_orb), D(m_orb), D_proj(m_orb), psi(m_orb)
         real(dp) :: x_orbit, y_orbit, z_orbit, psi_orbit
         real(dp) :: date(m_orb)
         integer :: i

         ! PLACEHOLDER: replace this pseudo-orbit with a self-consistent model.
         ! The current values are not a validated orbital solution.
         if (verbose) print *, 'Calculate placeholder orbital evolution'

         open(11,file=output_file('orbit.dat'))
         do i = 1,m_orb
            date(i) = date_ref_epoch + (i-1) ! Dummy values
            phase(i) = (i-1.d0)/float(m_orb) ! Not independent of true anomaly.
            psi(i) = 85.d0 * deg2rad
            D(i) = a_orbit/abs(sin(psi(i))) * (0.5 + e_orbit * phase(i))
            psi(i) = pi - psi(i) ! If star 1 is in front.
            D_proj(i) = D(i)*abs(sin(psi(i)))/(distance*AU)
            write(11,*) phase(i), date(i), D(i), D_proj(i), psi(i)
         enddo
         close(11)

         ! Normalised 3D position (divided by D).
         x_orbit = cos(w_orbit - v_orbit) * cos(i_orbit)
         y_orbit = sin(w_orbit - v_orbit)
         z_orbit = cos(w_orbit - v_orbit) * sin(i_orbit)

         psi_orbit = atan(sqrt(x_orbit**2 + y_orbit**2)/z_orbit)
         if (verbose) write(*,97) 'psi_orbit(deg) = ', psi_orbit/deg2rad
         97 format (A17,F6.2)
      end subroutine multi_epoch_orbit

   end subroutine orbit_run

end module orbit

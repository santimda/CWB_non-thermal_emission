program ffa_map_driver
   use global
   use initialize
   use orbit, only: orbit_run
   use CD, only: CD_run
   use CD_thermo, only: CD_thermo_run
   use CD_dist_e, only: CD_dist_e_run
   use CD_rad_syn, only: vec_nu
   use opac_ff, only: opac_ff_run
   implicit none

   integer, parameter :: nx_map = 100
   integer, parameter :: ny_map = 100
   real(dp), parameter :: x_min = -1.5d0
   real(dp), parameter :: x_max = 2.5d0
   real(dp), parameter :: y_min = -2.d0
   real(dp), parameter :: y_max = 2.d0
   real(dp), parameter :: view_angles(2) = (/75.d0, 105.d0/)*deg2rad
   real(dp), parameter :: frequencies(3) = (/0.3d9, 1.d9, 5.d9/)

   real(dp) :: phase(m_orb), separation(m_orb), projected_separation(m_orb)
   real(dp) :: orbital_angle(m_orb)
   real(dp) :: tauff(ml,mnu)
   real(dp) :: x_fraction, y_fraction, angle
   integer :: ii, jj, kk, view_index, unit
   character(len=3), parameter :: view_labels(2) = (/'75 ', '105'/)

   interface
      subroutine total_tau_calc(D, incli, tauff)
         use global
         real(dp), intent(in) :: D, incli
         real(dp), intent(out) :: tauff(ml,mnu)
      end subroutine total_tau_calc
   end interface

   call initialize_allocate()
   call orbit_run(phase, separation, projected_separation, orbital_angle, single_epoch)
   call CD_run(separation(1))
   call CD_thermo_run('1', separation(1))
   call CD_thermo_run('2', separation(1))
   call CD_dist_e_run('1', separation(1))
   call CD_dist_e_run('2', separation(1))

   ! The validation copy uses mnu=3 and vec_nu supplies the explicit test frequencies.
   call vec_nu()
   call opac_ff_run(separation(1))

   if (maxval(abs(nu - frequencies)) > 1.d-6*maxval(frequencies)) then
      error stop 'Validation frequency vector was not initialized as requested'
   endif

   do view_index = 1, size(view_angles)
      angle = view_angles(view_index)
      unit = 900 + view_index
      open(unit=unit, file=output_file('ffa_map_psi'//trim(adjustl(view_labels(view_index)))//'.dat'), &
           status='replace', action='write')

      do jj = 1, ny_map
         y_fraction = y_min + (y_max-y_min)*(jj-1.d0)/(ny_map-1.d0)
         do ii = 1, nx_map
            x_fraction = x_min + (x_max-x_min)*(ii-1.d0)/(nx_map-1.d0)
            x(ii) = x_fraction*separation(1)
            y(ii) = y_fraction*separation(1)
            z(ii) = 0.d0
         enddo

         call total_tau_calc(separation(1), angle, tauff)
         do ii = 1, nx_map
            write(unit,*) x(ii)/separation(1), y(ii)/separation(1), &
               (exp(-tauff(ii,kk)), kk=1,mnu)
         enddo
      enddo
      close(unit)
   enddo

   call initialize_deallocate()

end program ffa_map_driver

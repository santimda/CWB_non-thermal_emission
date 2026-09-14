program gaunt_driver
  use global, only: dp, T_w1, T_w2, Zq_w1, Zq_w2
  implicit none

  interface
    subroutine gaunt_ff_calc(nu_ff, T_w, Zq_w, gaunt_ff)
      use global, only: dp
      real(dp), intent(in) :: nu_ff(:), T_w, Zq_w
      real(dp), intent(out) :: gaunt_ff(:)
    end subroutine gaunt_ff_calc
  end interface

  integer, parameter :: n_frequency = 50
  real(dp) :: nu_ff(n_frequency), gaunt(n_frequency)
  integer :: i, j

  do i = 1, n_frequency
    nu_ff(i) = 10.d0**(8.d0 + 4.d0*real(i-1, dp)/real(n_frequency-1, dp))
  end do

  open (10, file='validation/free-free_absorption/results/gaunt_fortran.dat', &
        status='replace')
  write (10, '(a)') '# T_w Zq_w frequency_Hz gaunt_factor'

  do j = 1, 2
    if (j == 1) then
      call gaunt_ff_calc(nu_ff, T_w1, Zq_w1, gaunt)
      do i = 1, size(nu_ff)
        write (10, '(4(es24.16,1x))') T_w1, Zq_w1, nu_ff(i), gaunt(i)
      end do
    else
      call gaunt_ff_calc(nu_ff, T_w2, Zq_w2, gaunt)
      do i = 1, size(nu_ff)
        write (10, '(4(es24.16,1x))') T_w2, Zq_w2, nu_ff(i), gaunt(i)
      end do
    end if
  end do

  close (10)
end program gaunt_driver

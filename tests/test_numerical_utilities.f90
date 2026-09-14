program test_numerical_utilities
   use global, only: dp, vector_log, integrate
   implicit none

   real(dp), parameter :: tolerance = 1.d-12
   real(dp), parameter :: integration_tolerance = 1.d-8
   real(dp) :: values(5), expected_values(5), ratio
   real(dp) :: x_odd(5), f_odd(5), integral_odd
   real(dp) :: x_even(4), f_even(4), integral_even
   real(dp), allocatable :: x_inverse(:), f_inverse(:)
   real(dp) :: xmin, xmax, step, integral_inverse, expected_inverse
   integer :: i, n_inverse

   ! Test the logarithmically spaced vector constructor.
   call vector_log(1.d0, 1.d4, ratio, values)
   expected_values = (/ 1.d0, 10.d0, 100.d0, 1000.d0, 10000.d0 /)
   call assert_close('vector_log values', maxval(abs(values-expected_values)), 0.d0, tolerance)
   call assert_close('vector_log ratio', ratio, 10.d0, tolerance)

   ! Test Simpson integration with an odd number of points.
   x_odd = (/ 0.d0, 1.d0, 2.d0, 3.d0, 4.d0 /)
   f_odd = x_odd**2
   call integrate(5, 1.d0, f_odd, integral_odd)
   call assert_close('integrate odd number of points', integral_odd, 64.d0/3.d0, tolerance)

   ! Test Simpson integration with an even number of points.
   x_even = (/ 0.d0, 1.d0, 2.d0, 3.d0 /)
   f_even = x_even**2
   call integrate(4, 1.d0, f_even, integral_even)
   call assert_close('integrate even number of points', integral_even, 9.d0, tolerance)

   ! Test the integral of 1/x**2 against its analytical result.
   xmin = 1.d0
   xmax = 100.d0
   n_inverse = 10001
   step = (xmax-xmin)/dble(n_inverse-1)
   allocate(x_inverse(n_inverse), f_inverse(n_inverse))
   do i = 1, n_inverse
      x_inverse(i) = xmin + dble(i-1)*step
      f_inverse(i) = 1.d0/x_inverse(i)**2
   end do
   call integrate(n_inverse, step, f_inverse, integral_inverse)
   expected_inverse = 1.d0/xmin - 1.d0/xmax
   call assert_close('integrate inverse square', integral_inverse, expected_inverse, integration_tolerance)
   deallocate(x_inverse, f_inverse)

   print *, 'Numerical utility tests passed.'

contains

   subroutine assert_close(name, actual, expected, allowed_error)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual, expected, allowed_error
      real(dp) :: error

      error = abs(actual-expected)
      if (error > allowed_error) then
         print *, 'FAILED:', trim(name)
         print *, 'actual, expected, error =', actual, expected, error
         stop 1
      end if
   end subroutine assert_close

end program test_numerical_utilities

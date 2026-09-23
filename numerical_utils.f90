module numerical_utils
   use constants, only: dp
   implicit none
   public

contains

   function safe_acos_dot(V1,V2) result(angle)
      ! Protect acos from roundoff taking a dot product just outside [-1, 1].
      real(dp), intent(in) :: V1(:), V2(:)
      real(dp) :: angle, cosine

      cosine = max(-1.d0, min(1.d0, dot_product(V1,V2)))
      angle = acos(cosine)
   end function safe_acos_dot


   subroutine vector_log(Vmin,Vmax,Vint,Vvec)
      ! Build logarithmically-spaced vector.
      real(dp), intent(in) :: Vmin,Vmax
      real(dp), intent(out) :: Vint
      real(dp), intent(out) :: Vvec(:)
      integer :: N_size, i_vec

      N_size = size(Vvec)
      Vint = (Vmax/Vmin)**(1.d0/real(N_size-1,dp))
      Vvec(1) = Vmin
      do i_vec = 2, N_size
         Vvec(i_vec) = Vvec(i_vec-1)*Vint
      end do
   end subroutine vector_log


   subroutine integrate(N,Hs,FI,S)
      ! Integrate f(x) with Simpson's rule.
      integer, intent(in) :: N
      real(dp), intent(in) :: Hs,FI(N)
      real(dp), intent(out) :: S
      real(dp) :: S0,S1,S2
      integer :: i_int

      S = 0.d0
      S0 = 0.d0
      S1 = 0.d0
      S2 = 0.d0
      do i_int = 2, N-1, 2
         S1 = S1 + FI(i_int-1)
         S0 = S0 + FI(i_int)
         S2 = S2 + FI(i_int+1)
      end do
      S = Hs*(S1 + 4.d0*S0 + S2)/3.d0
      ! If N is even, add the last slice separately.
      if (mod(N,2) == 0) S = S + Hs*(5.d0*FI(N) + 8.d0*FI(N-1) - FI(N-2))/12.d0
   end subroutine integrate


   subroutine unit_vector_diff(V1,V2,N,versor)
      ! Build the unit vector (V1-V2)/|V1-V2|.
      integer, intent(in) :: N
      real(dp), intent(in) :: V1(N),V2(N)
      real(dp), intent(out) :: versor(N)
      real(dp) :: difference(N), difference_norm

      difference = V1 - V2
      difference_norm = norm2(difference)
      if (difference_norm > 0.d0) then
         versor = difference/difference_norm
      else
         versor = 0.d0
      end if
   end subroutine unit_vector_diff


end module numerical_utils

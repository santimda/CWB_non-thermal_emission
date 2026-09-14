module CD
   use global
   implicit none

   public :: CD_run

contains

   !=============================================================================
   ! CD_run (START)
   !=============================================================================
   ! This program calculates the position of the contact discontinuity as a 2-D 
   ! curve in the X-Y plane. This curve has azimuthal symmetry. The positions (x,y) 
   ! of the CD are stored in vectors x(n),y(n) such that (x(1),y(1)) is the apex 
   ! and for i>1 points in the positive y direction. Use a RK method to obtain y(i) 
   ! given x(i-1), y(i-1) and dx. The function dy/dx is defined as an external g. 
   ! Star1 is at (0,0) and star2 at (D,0); r1 and r2 are the distances from 
   ! (x,y) to star1 and star2, respectively. The stars are assumed to be far enough 
   ! as for the winds to reach v_infty (so that \eta is constant). Orbital motion 
   ! is neglected. 
   ! The program also provides a hyperbolic fit to the curve, stored in a_CD.dat
   !==============================================================================

   subroutine CD_run(D)
      real(dp), intent(in) :: D
      real(dp) :: x_max
      real(dp) :: x_step,dx,dl,dl_max,L_total
      real(dp), parameter :: dx_factor = 1.1  ! This factor controls some resolution adjustment
      real(dp), parameter :: xmax_factor = 1.8 ! This factor controls x_max (thus size of the WCR)
      
      interface
         Subroutine RK(x1,y1,dx,y2,D)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp), intent(out) :: y2
            real(dp), intent(in) :: x1,y1,dx,D
         end Subroutine RK
      end interface

      interface
         function g(a,b,D)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp) :: g
            real(dp), intent(in) :: a,b,D
         end function g
      end interface

      if (verbose) print *, '\n Calculate the CD position'

      !=======================================================================

      ! Fix x_max in function of eta so that the total length of the WCR is ~2D
      x_max = (1.d0 + xmax_factor*cos(theta_inf)) * (D/AU)

      !Work with x and y in AU units
      y(1) = 1.d-5 ! y(1)=0.d0 brings numerical problems
      x(1) = (D/AU)/(1.d0+sqrt(eta))
      x(1+ml_2) = x(1)
      y(1+ml_2) = -y(1)
      if (x(1) < 10.d0*Rst1/AU .or. (D/AU)-x(1) < 10.d0*Rst2/AU) &
         print *, 'WARNING: STARS TOO CLOSE, WIND TERMINAL VELOCITY NOT REACHED'

      ! Initial guess for a variable step that gives higher resolution near the apex; x_max propto (D_proj/UA) is defined in global.f90
      x_step = (x_max/x(1))**(1.d0/float(ml_2-1))
      L_total = 0.d0
      do i = 2,ml_2
         dx = x(i-1)*(x_step-1.d0)
         dl_max = 2.d0*dx
         ! if dl = sqrt(dx^2+dy^2) is too large, use a smaller step
         dl = dl_max*1.01d0   ! just so it enters in the do while
         dx = dx*dx_factor        ! It is divided by dx_factor inside the loop
         do while (dl > dl_max) 
            dx = dx/dx_factor     ! make the step smaller if more resolution is needed
            x(i)= x(i-1) + dx
            call RK(x(i-1),y(i-1),dx,y(i),D)
            dl = sqrt( dx**2 + (y(i)-y(i-1))**2 )
         end do
         x(i+ml_2) = x(i)
         y(i+ml_2) = -y(i)
         ! Re-adapt the step in case it was reduced inside the loop
         if (i < ml_2) then
            x_step = (x_max/x(i))**(1.d0/float(ml_2-i))
         end if
         L_total = L_total + dl
      end do
      if (verbose) write(*,98) 'D [AU] =', D/AU   
      if (verbose) write(*,98) 'L_total [D] =', L_total/(D/AU)      
      if (verbose) write(*,98) 'th_inf [deg] =', theta_inf * 180 / pi   
      if (verbose) write(*,98) 'max_ang [deg] =', 180 - atan2( y(ml_2), D/AU - x(ml_2) ) * 180 / pi     
      98  format (A15,X,F7.2)

      x = x * AU  ! Convert to cgs
      y = y * AU  ! Convert to cgs

      !----------------------------------------------------------------------------
      ! Save the results
      open(11,file=output_file('CD.dat'))
      do l = 1, ml_2
         write(11,*) x(l), y(l) 
      end do
      close(11)
      z = 0.d0   ! The z coordinate is zero (it is rotated later)

      !----------------------------------------------------------------------------
      ! Obtain a tangent vector to the CD at each location using the function that gives dy/dx
      open(11,file=output_file('n_tan.dat'))
      do l = 1, ml_2
         write(11,*) l, g(x(l)/AU,y(l)/AU,D)
      end do
      close(11)

      call CD_approx(x/AU,y/AU)

   end subroutine CD_run

   !=============================================================================
   ! CD_run (END)
   !=============================================================================


   !=============================================================================
    ! CD_approx (START)
   !=============================================================================
    subroutine CD_approx(x_CD,y_CD)
      ! Analytical hyperboloid approximation used for FFA wind selection.
      real(dp), intent(in) :: x_CD(ml), y_CD(ml)
      real(dp) :: a_CD, b_CD, x_approx(ml)
      integer :: il, l_CD
      
      open(190, file=output_file('CD_approx.dat'))  ! Approximate solution for plotting
      open(191, file=output_file('a_CD.dat'))       ! Coefficient(s) for the approximation

      ! Better fit taking b such that the approximation coincides with the CD
      ! at an intermediate position l_CD
      l_CD = int(ml/2.3)-1
      ! Hyperbola: (x/a)**2 - (y/b)**2 = 1
      ! (a,0) = (x(1),0) ; b/a = tan(theta_inf) = tan(y_max/x_max)
      a_CD = x_CD(1)
      b_CD = y_CD(l_CD) / sqrt( (x_CD(l_CD)/a_CD)**2 - 1.d0 )
      write(191,*) a_CD*AU, b_CD*AU
      do il = 1, ml_2
         x_approx(il) = a_CD * sqrt(1.d0 + (y_CD(il)/b_CD)**2)
         x_approx(il+ml_2) = x_approx(il)
      end do
      do il = 1, ml
         write(190,*) x_approx(il)*AU, y_CD(il)*AU
      end do
      close(190)
      close(191)

    end subroutine CD_approx

   !=============================================================================
   ! CD_approx (END)
   !=============================================================================


end module CD


!=====================================================================

subroutine RK(xi, yi, dx, yf, D)
   implicit none
   integer, parameter :: dp = kind(1.d0)
   real(dp), intent(in) :: xi, yi, dx, D
   real(dp), intent(out) :: yf
   real(dp) :: k1, k2, k3, k4, h
   integer :: i
   integer, parameter :: n=1000
   real(dp), dimension(n) :: y  

   interface
      function g(a, b, D)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp) :: g
         real(dp), intent(in) :: a, b, D
      end function g
   end interface

   y(1) = yi
   h = dx / dble(n)

   do i = 1, n - 1
      k1 = h * g(xi + (i - 1) * h, y(i), D)
      k2 = h * g(xi + (i - 1) * h + h / 2.d0, y(i) + k1 / 2.d0, D)
      k3 = h * g(xi + (i - 1) * h + h / 2.d0, y(i) + k2 / 2.d0, D)
      k4 = h * g(xi + i * h, y(i) + k3, D)
      y(i + 1) = y(i) + (k1 + 2.d0 * k2 + 2.d0 * k3 + k4) / 6.d0
   end do

   yf = y(n)

end subroutine RK

!=====================================================================

Function g(a_i, a_f,D)
   ! g = dy/dx (Antokhin+2004, Eq. 6 gives dx/dy)
   use global  
   real(dp) :: g
   real(dp),intent(in) :: a_i, a_f, D
   real(dp) :: r1, r2
   r1 = sqrt( a_i**2 + a_f**2 )
   r2 = sqrt( (D/AU - a_i)**2 + a_f**2 )
   g = a_f/(a_i-(D/AU*r1**2*sqrt(eta)/(r1**2*sqrt(eta)+r2**2)))

End function g

!=====================================================================

module WCR_conv
  use global
  
  implicit none
  private
  public :: WCR_conv_run
  
contains

  !============================================================================
  ! WCR_conv_run
  !============================================================================
  ! This program makes the convolution of the model emission map at the frequency
  ! map_nu and the instrument (gaussian) response

  subroutine WCR_conv_run(D_proj,incli,map_nu)

    use global
    implicit none
    real(dp), intent(in) :: D_proj,incli
    real(dp), intent(in) :: map_nu
    integer, parameter :: ml_tot = 2*(mphi-1)*ml   ! Total number of data points.
    ! Relation between theta and sigma for the beam
    !real(dp), parameter :: coc = sqrt(8.d0*log(2.d0)) 
    !real(dp), parameter :: theta_x = sigma_x*coc, theta_y = sigma_y*coc 
    ! Use directly sigma 
    real(dp), parameter :: sigma_x = 5.6d0, sigma_y = 11.3d0    ! Marcote+2021
    real(dp) :: x_map(ml_tot),y_map(ml_tot)
    real(dp) :: F_abs(ml_tot), dummy
    real(dp) :: width_x, width_y, D_x, D_y, ci, si
    real(dp) :: x_c0,y_c0,weight_distance_sq,attenuation,sum_wcr,mJy_beam
    real(dp), dimension(:), allocatable :: x_c,y_c
    real(dp), dimension(:,:), allocatable :: F_map_abs
    integer :: m_x,m_y
    integer, parameter:: m_map = mphi*mr_map
    real(dp) :: A_beam,A_cell,sum_s1,x_map2(m_map),y_map2(m_map),S1_map2(m_map)
    real(dp) :: sum1, sum2, sum_s2, x_map3(m_map),y_map3(m_map),S2_map3(m_map)
    real(dp) :: w_x, w_y
    character(len=16) :: freq_label

    freq_label = map_frequency_label(map_nu)

    A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
    A_cell = f_Ny**2*sigma_x*sigma_y  ! cell sample size
    
    ci = abs(cos(incli)) ! incli is the angle psi such that D_proj = D * sin(incli)
    si = abs(sin(incli))

    !-------------------------------------------------------
    ! Initialize the data for the convolution region
    ! The 1.5 factor is arbitrary to give an extra margin
    width_x = ( D_proj * (1.d0 - eta + cos(theta_inf)) + D_proj/si * ci ) * 1.5d0 + 6.d0*sigma_x ! [mas]
    width_y = ( 2.d0 * D_proj/si * sin(theta_inf) )*1.5d0 + 6.d0*sigma_y    ! [mas]

    D_x = f_Ny*sigma_x
    D_y = f_Ny*sigma_y
    m_x = int(width_x/D_x) ! Number of ellipses in the x-direction
    m_y = int(width_y/D_y) ! Number of ellipses in the y-direction
    write(*,*) 'Convolving with m_x=',m_x, 'and m_y=', m_y
    x_c0 = -D_proj/si * sin(theta_inf) * ci * 1.1d0 - 3.d0*sigma_x ! Center of the first ellipse
    y_c0 = -width_y/2.d0   ! Symmetric w.r.t. y = 0
    allocate( x_c(m_x), y_c(m_y), F_map_abs(m_x,m_y)  )

    !------------------------------------------------------
    ! Data reading

    sum_wcr = 0.d0
    open(13,file=output_file('WCR_radiomap_'//trim(freq_label)//'GHz.dat'))
    do l = 1,ml_tot      ! The 3rd column contains the unabsorbed flux, not used
      read(13,*) x_map(l),y_map(l),dummy,F_abs(l)
      sum_wcr = sum_wcr + F_abs(l)
    end do
    if (cos(incli) < 0.d0) x_map = - x_map ! Mirror if needed

    write(*,*) 'Pre-conv:'
    write(*,80) 'S('//trim(freq_label)//'GHz)=',sum_wcr,'mJy'
    close(13)
    sum1 = 0.d0
    sum2 = 0.d0
    open(133,file=output_file('S1_ff_map_'//trim(freq_label)//'G.dat')) ! Star 1
    open(134,file=output_file('S2_ff_map_'//trim(freq_label)//'G.dat')) ! Star 2
    do l = 1,m_map
      read(133,*) x_map2(l),y_map2(l),S1_map2(l)
      read(134,*) x_map3(l),y_map3(l),S2_map3(l)
      sum1 = sum1 + S1_map2(l)
      sum2 = sum2 + S2_map3(l)
    end do
    write(*,80) 'S1_ff('//trim(freq_label)//'GHz)=',sum1,'mJy'
    write(*,80) 'S2_ff('//trim(freq_label)//'GHz)=',sum2,'mJy'
    close(133)
    close(134)

    !------------------------------------------

    w_x = 2.d0*sigma_x**2
    w_y = 2.d0*sigma_y**2
    do i = 1,m_x
      x_c(i) = x_c0 + (i-0.5)*D_x  ! x-coord of the center of the i-th ellipse

      do j = 1,m_y
        y_c(j) = y_c0 + (j-0.5)*D_y  ! y-coord of the center of the (x_i,y_j) ellipse

        sum_wcr = 0.d0
        do l = 1,ml_tot
           weight_distance_sq = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
           if (weight_distance_sq < 6.d0) then ! flux drops < 1e-3 farther from there, no need to sum
            ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
            attenuation = exp( -weight_distance_sq ) 
            sum_wcr = sum_wcr + F_abs(l)*attenuation
          end if
        end do

        sum_s1 = 0.d0
        do l = 1,m_map
          weight_distance_sq = ((x_map2(l)-x_c(i))**2)/w_x + ((y_map2(l)-y_c(j))**2)/w_y
          if (weight_distance_sq < 6.d0) then ! flux drops < 1e-3 farther from there, no need to sum
            ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
            attenuation = exp( -weight_distance_sq ) 
            sum_s1 = sum_s1 + S1_map2(l)*attenuation
          end if
        end do

        sum_s2 = 0.d0
        do l = 1,m_map
          weight_distance_sq = ((x_map3(l)-x_c(i))**2)/w_x + ((y_map3(l)-y_c(j))**2)/w_y
          if (weight_distance_sq < 6.d0) then ! flux drops < 1e-3 farther from there, no need to sum
            ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
            attenuation = exp( -weight_distance_sq ) 
            sum_s2 = sum_s2 + S2_map3(l)*attenuation
          end if
        end do

        F_map_abs(i,j) = sum_wcr + sum_s1 + sum_s2

      end do

    end do

    !///////////////////////////////////////////////////////////////
    ! Write data and do consistency check

    mJy_beam = A_cell/A_beam ! Transform mJy/beam -> mJy (for checks)

    sum_wcr = 0.d0
    open (17,file=output_file('WCR_conv_'//trim(freq_label)//'Gabs.dat'))
    do i = 1,m_x
      do j = 1,m_y
        write(17,61) x_c(i),y_c(j),F_map_abs(i,j)
        sum_wcr = sum_wcr + F_map_abs(i,j)*mJy_beam
        !print*, sum_wcr, F_map_abs(i,j), mJy_beam
      end do
      write(17,*) 
    end do
    write(*,*) 'Post-conv:'
    write(*,80) 'S('//trim(freq_label)//'GHz)=',sum_wcr,'mJy'
    close(17)

    61  format(F11.4,X,F11.4,X,E11.4)
    80  format(A11,X,Es10.2,X,A3)

    deallocate(x_c,y_c,F_map_abs)

  END subroutine WCR_conv_run


END module WCR_conv

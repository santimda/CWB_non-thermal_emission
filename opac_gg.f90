module opacity_gg
  use global
  use numerical_utils, only: integrate, vector_log
  
  implicit none
  private
  public :: opacity_gg_run

  contains

  !============================================================================
  ! opacity_gg_run
  !============================================================================
  ! Calculate a matrix with the opacity coefficients at a distance of one stellar
  ! radii for different angles between the emitter-star direction and the observer, 
  ! and for different gamma-ray photon energies: tau<i>(Ef,angle)

    subroutine opacity_gg_run
      implicit none
      allocate( tau1(mE_gg,m_ang_gg),tau2(mE_gg,m_ang_gg) )
      
      write(*,*) 'Calculate gamma-gamma opacities'

      call opac_gg(tau1,T1,Rst1)
      open (88,file=output_file('opac_gg1.dat'))
      do i=1,mE_gg
         write(88,*) (tau1(i,j), j=1,m_ang_gg)
      enddo
      close(88)
      
      !--------------------------------------------------------------

      call opac_gg(tau2,T2,Rst2)
      open (89,file=output_file('opac_gg2.dat'))
      do i=1,mE_gg
         write(89,*) (tau2(i,j), j=1,m_ang_gg)
      enddo
      close(89)
            
      deallocate( tau1,tau2 )
        
    END subroutine opacity_gg_run
    
  END module opacity_gg
  
  
  !--------------------------------------------------------------------
  !
  !     Subroutines
  !
  !--------------------------------------------------------------------
  
  subroutine opac_gg(tau_gg,zT,rstar)
    use global
    use numerical_utils, only: integrate, vector_log
    implicit none
    integer, parameter :: ir_max=100
    integer, parameter :: ie0_max = 20
    real(dp), intent(in) :: zT,rstar
    real(dp), intent(out) :: tau_gg(mE_gg,m_ang_gg)
    integer :: iang,ie,ir,ie0    
    real(dp) :: Cnumc2,tau,Ef(mE_gg),dang,ang,dr,xi,yi,di
    real(dp) :: mu,nu0,Bb,Be,dEps_int,EpsM,Eps0,Eps0_min_factor,Eps0_min,Eps0_max,Eps_int
    real(dp) :: cte_sigma_gg,Nbmc2,Nst,xe,sigma_gg,sum
    real(dp) :: f(ir_max),cang(m_ang_gg),sang(m_ang_gg)
    real(dp) :: dil(m_ang_gg,ir_max),mu_factor(m_ang_gg,ir_max)

    Cnumc2 = 1.23d20
    
    ! Ef_min_gg and Ef_max_gg are defined in global.f90, and also mE_gg

    EpsM = (3.d0*k*zT)/mec2
    Eps0_max = 20.d0*EpsM
    cte_sigma_gg = (0.5d0*pi*(qe**2/mec2)**2)
    call vector_log(Ef_min_gg,Ef_max_gg,Ef_int_gg,Ef)

    Ef = Ef/mec2 ! Adimensional

    ! Pre-compute angular quantities for optimization
    dang = pi/m_ang_gg
    ang = dang/2.d0
    do iang=1, m_ang_gg
      cang(iang) = cos(ang)
      sang(iang) = sin(ang)
      ang = ang+dang
    enddo

    ! Pre-compute paths and the geometrical factors 1-mu and dilution used in the loops
    dr = rstar/15.d0 ! R_final = Rstar + Rstar*ir_max/15
    do iang=1,m_ang_gg
      xi = Rstar
      yi = 0.d0
      do ir=1,ir_max
        xi = dr*cang(iang) + xi
        yi = dr*sang(iang) + yi
        di = sqrt(xi**2 + yi**2)
        dil(iang,ir) = Rstar**2 / di**2
        mu = (xi*cang(iang) + yi*sang(iang)) / di
        mu_factor(iang,ir) = 1.d0 - mu  
      enddo
    enddo

    !$omp parallel do default(none) schedule(static) &
    !$omp& shared(Ef,Eps0_max,cte_sigma_gg,Cnumc2,dil,mu_factor,dr,tau_gg,zT) &
    !$omp& private(iang,ir,ie0,Eps0_min_factor,Eps0_min,Eps_int,Eps0, &
    !$omp& dEps_int,nu0,Bb,nbmc2,nst,xe,Be,sigma_gg,sum,tau,f)
    do ie=1,mE_gg
      Eps0_min_factor = 1.d0/Ef(ie)/0.5d0

      do iang=1,m_ang_gg   ! loop in interaction angle: ang
        tau = 0.d0

        do ir=1,ir_max     ! loop in gamma-ray photon trayectory: xi,yi
          sum = 0.d0    
          Eps0_min = Eps0_min_factor/mu_factor(iang,ir)
          Eps_int = (Eps0_max/Eps0_min)**(1.d0/float(ie0_max))
          
          if (Eps0_min <= Eps0_max) then  ! else no absorption (sum=0)
            Eps0 = Eps0_min*Eps_int**0.5d0  

            do ie0=1,ie0_max ! loop in stellar photon energies, Eps0
              dEps_int = Eps0*(Eps_int-1.d0)

              ! Density of target photons
              nu0 = Eps0*Cnumc2
              Bb = (2.d0*h*nu0**3/c**2)/(exp(h*nu0/(k*zT))-1.d0)
              nbmc2 = Bb*Cnumc2/(nu0*h)
              nst = (nbmc2*pi/c) * dil(iang,ir)

              ! sigma gamma-gamma (Gould & Schreder 1967)
              xe = Ef(ie)*Eps0*0.5d0*mu_factor(iang,ir)
              if (xe >= 1.d0) then
                Be = sqrt(1.d0-1.d0/xe)
                sigma_gg = cte_sigma_gg*(1.d0-Be**2) &
                  *((3.d0-Be**4)*log((1.d0+Be)/(1.d0-Be)) &
                  -2.d0*Be*(2.d0-Be**2))
                sum = sum + sigma_gg*nst*mu_factor(iang,ir)*dEps_int
              end if

              Eps0 = Eps0*Eps_int
            enddo
          endif        
          f(ir) = sum
        enddo
        call integrate(ir_max,dr,f,tau)
        tau_gg(ie,iang) = tau ! Opacity: tau(Efe,ang)

      enddo

    enddo
    !$omp end parallel do

    return
    
  end subroutine opac_gg
  
  

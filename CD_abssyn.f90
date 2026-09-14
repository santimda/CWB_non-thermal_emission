module CD_abssyn
   use global
   implicit none
   
   public :: CD_abssyn_run, CD_abssyn_data_run

contains

   !========================================================================
   !        CD_abssyn_run        
   !========================================================================
   ! Calculates the FFA (+ R-T) opacity for a radio photon emitted at 
   ! (x(l),y(l),z(l)). Output: radsyn.dat -> E(i),L_sy(E(i), L_sy_abs(E(i)). 
   ! Input: tauff<i>.dat from opac_ff.f90, which has the wind opacities at a 
   ! distance of R_star and for different energies and directions: 
   ! tauff<i>(mE_ff,m_ang_ff). It assumes that psi > 0� so that the star 1
   ! is at front. FFA in the WCR is neglected due to high T plasma.
   !========================================================================

   subroutine CD_abssyn_run(Lsy_sum,Lsy_sum_abs,D,incli)  
      ! $ use omp_lib
      Real(dp), intent(in) :: D,incli
      Real(dp), intent(out) :: Lsy_sum(mnu),Lsy_sum_abs(mnu)
      Real(dp) :: sinc,cinc,tinc
      Real(dp) :: Ef(mnu),log_nu_int,flux
      Real(dp) :: u_proj,v_proj ! coordinates in the plane of the sky
      Real(dp) :: emission(n_map_freq),emission_abs(n_map_freq),sum,convert,dummy
      Integer :: nl,i_map,i_nu(n_map_freq),map_unit(n_map_freq)
      Character(len=16) :: freq_label
      allocate( tau_ff(ml,mnu) )
      tau_ff = 0.d0
      Lsy_sum = 0.d0
      Lsy_sum_abs = 0.d0
      convert = 206265.d3/(distance*kpc) ! convert linear to angular scale in mas
      flux = 1.d26*h/(4.d0*pi*(distance*kpc)**2)! L[erg/s/erg] -> F[mJy]

      log_nu_int = log(nu(2)/nu(1))

      sinc = sin(incli)
      cinc = cos(incli)
      tinc = tan(incli)

      ! set omp threads number
      ! $ call omp_set_num_threads(nth)  

      ! Calculate total opacity at each location for different photon energies

      call total_tau_calc(D,incli,tau_ff)

      open (25,file=output_file('CD_tauff.dat'))
      do l=1,ml    
         do i=1,mnu
            write(25,*) l,nu(i),tau_ff(l,i) 
         end do
      end do
      close(25)

      open (31,file=output_file('emi_syn1.dat'))
      open (32,file=output_file('emi_syn2.dat'))
      do l=1,ml
         do i=1,mnu
            read(31,*) nl,Ef(i),dummy,Lsy1(l,i)
            read(32,*) nl,Ef(i),dummy,Lsy2(l,i)
        end do
      end do
      close(31)
      close(32)

      ! Calculate total (abs and unabs) emission for each linear emitter. Also store
      ! the emission projected in the plane of the sky (u_proj,v_proj) for maps

      do i_map = 1, n_map_freq
         i_nu(i_map) = int(log(map_freq(i_map)/nu(1))/log_nu_int) + 1
         freq_label = map_frequency_label(map_freq(i_map))
         map_unit(i_map) = 69 + i_map
         open(map_unit(i_map),file=output_file('CD_radiomap_'//trim(freq_label)//'GHz.dat'))
      end do

      do l=1,ml
         ! u*cos(psi) = -y_proj, with y_proj the projection of y(l) in the u-v plane
         u_proj = -abs(cinc)*( -tinc*x(l)+y(l) )*convert 
         v_proj = z(l)*convert

         do i=1,mnu
            sum = Lsy1(l,i) + Lsy2(l,i)
            Lsy_sum(i) = Lsy_sum(i) + sum
            if (nu(i) < nu_max_ff) then ! frequencies for which absorption was calculated 
               Lsy_sum_abs(i) = Lsy_sum_abs(i) + sum*exp(-tau_ff(l,i))  
            else
               Lsy_sum_abs(i) = Lsy_sum_abs(i) + sum
            end if         
         end do
         
         do i_map = 1, n_map_freq
            emission(i_map) = (Lsy1(l,i_nu(i_map)) + Lsy2(l,i_nu(i_map)))*flux
            emission_abs(i_map) = emission(i_map)*exp(-tau_ff(l,i_nu(i_map)))
            write(map_unit(i_map),*) u_proj, v_proj, emission(i_map), emission_abs(i_map)
         end do
      end do
      do i_map = 1, n_map_freq
         close(map_unit(i_map))
      end do

      open(33,file=output_file('rad_syn_abs.dat'))
      do i=1,mnu
         write(33,*) Ef(i), Lsy_sum(i), Lsy_sum_abs(i)
      end do
      close(33)

      deallocate( tau_ff )

   END subroutine CD_abssyn_run

   !==========================================================================

   !========================================================================
   !        CD_abssyn_data_run        
   !========================================================================
   ! Similar as CD_abssyn_run but only for the specific frequencies 
   ! for which observations are available
   !========================================================================

   subroutine CD_abssyn_data_run(nu_data, S_sum, S_sum_abs, D, incli)
      Real(dp), intent(in) :: D,incli,nu_data(mnu_data)
      Real(dp), intent(out) :: S_sum(mnu_data), S_sum_abs(mnu_data)
      Real(dp) :: sinc,cinc,tinc
      Real(dp) :: Ef(mnu_data),flux
      Real(dp) :: Lsy_sum(mnu_data), Lsy_sum_abs(mnu_data), sum
      Integer :: nl
      allocate( tau_ff(ml,mnu_data) )
      tau_ff = 0.d0
      Lsy_sum = 0.d0
      Lsy_sum_abs = 0.d0
      flux = 1.d26*h/(4.d0*pi*(distance*kpc)**2)! L[erg/s/erg] -> F[mJy]

      sinc = sin(incli)
      cinc = cos(incli)
      tinc = tan(incli)

      ! Calculate total opacity at each location for different photon energies

      call total_tau_data_calc(D,incli,tau_ff)

      open (25,file=output_file('CD_tauff_data.dat'))
      do l=1,ml    
         do i=1,mnu_data
            write(25,*) l,nu_data(i),tau_ff(l,i) 
         end do
      end do
      close(25)

      open (31,file=output_file('emi_syn_data1.dat'))
      open (32,file=output_file('emi_syn_data2.dat'))
      do l=1,ml
         do i=1,mnu_data
            read(31,*) nl, Ef(i), Lsy1(l,i)
            read(32,*) nl, Ef(i), Lsy2(l,i)
        end do
      end do
      close(31)
      close(32)

      ! Calculate total (abs and unabs) emission for each linear emitter. 
      do l=1,ml
         do i=1,mnu_data
            sum = Lsy1(l,i) + Lsy2(l,i)
            Lsy_sum(i) = Lsy_sum(i) + sum
            Lsy_sum_abs(i) = Lsy_sum_abs(i) + sum*exp(-tau_ff(l,i))         
         end do
      end do
      S_sum = Lsy_sum * flux
      S_sum_abs = Lsy_sum_abs * flux

      open(33,file=output_file('rad_syn_abs.dat'))
      do i=1,mnu_data
         write(33,*) Ef(i), Lsy_sum(i), Lsy_sum_abs(i)
      end do
      close(33)

      deallocate( tau_ff )

   END subroutine CD_abssyn_data_run

   !==========================================================================
END module CD_abssyn

!==========================================================================
! Calculate the total f-f opacity (tau_ff) for photons emitted at each 
! location along the CD. 

subroutine total_tau_calc(D,incli,tauff)
   use global
   implicit none
   real(dp), intent(in) :: D,incli
   real(dp), intent(out) :: tauff(ml,mnu)
   real(dp) :: sinc,cinc,r1,r2,r2_i,r2_f,dr0
   real(dp) :: ang_int,ang_min,ang1,ang2_i,ang2_f
   real(dp) :: a_CD,b_CD,x1,y1,z1,r_cross(3),Xobs(3)
   real(dp), parameter :: sh_width = 0.05 ! 5%
   integer :: iang1,iang2_i,iang2_f,iEf
   logical :: cond 
   allocate( tau_ff1(mE_ff,m_ang_ff), tau_ff2(mE_ff,m_ang_ff) )

   sinc = sin(incli)
   cinc = cos(incli)

   open (288,file=output_file('opac_ff1.dat'))
   open (289,file=output_file('opac_ff2.dat'))
   do i=1,mE_ff 
      read(288,*) (tau_ff1(i,j), j=1,m_ang_ff)
      read(289,*) (tau_ff2(i,j), j=1,m_ang_ff)
   end do
   close(288)
   close(289)

   ang_int = pi/float(m_ang_ff)
   ang_min = ang_int/2.d0

   Xobs = (/ cinc,sinc,0.d0 /)
   Xst1 = (/ 0.d0, 0.d0, 0.d0 /)
   Xst2 = (/ D, 0.d0, 0.d0 /)

   if (verbose) print *, 'Calculate the total absorption for synchrotron radiation'

   ! To know which wind absorbs (1 or 2) we check if the position is within 
   ! the (approx.) surface of the CD

   open(191,file=output_file('a_CD.dat'))
   ! (x/a)**2 = 1 + (y/b)**2 + (z/b)**2
   read(191,*) a_CD, b_CD
   close(191)

   tauff = 0.d0

   do l = 1,ml        ! Position along the CD

      Xemi = (/ x(l),y(l),z(l) /) 

      if( incli < theta_inf ) then  ! only wind 2 absorbs
         call unit_vector_diff(Xemi,Xst2,3,XEe) ! XEe vector from the star to the emitter
         ang2_i = acos(DOT_PRODUCT(Xobs,XEe))! angle between Xobs and XEe
         iang2_i = int(m_ang_ff*abs(ang2_i-ang_min)/pi)+1 
         r2 = sqrt( (D-x(l))**2 + y(l)**2 + z(l)**2 )
         do i = 1, mnu      
            iEf = min(int(log(max(1.d0,nu(i)/nu_min_ff))/log(nu_int_ff)) + 1, mE_ff)
            tauff(l,i) = (Rst2/r2)**3*tau_ff2(iEf,iang2_i)
         end do  
      
      elseif( incli > pi - theta_inf ) then ! only wind 1 absorbs
         call unit_vector_diff(Xemi,Xst1,3,XEe) ! XEe vector from the star to the emitter
         ang1 = acos(DOT_PRODUCT(Xobs,XEe))! angle between Xobs and XEe
         iang1 = int(m_ang_ff*abs(ang1-ang_min)/pi)+1 
         r1 = sqrt( x(l)**2 + y(l)**2+z(l)**2 )
         do i = 1, mnu      
            iEf = min(int(log(max(1.d0,nu(i)/nu_min_ff))/log(nu_int_ff)) + 1,mE_ff)
            tauff(l,i) = (Rst1/r1)**3*tau_ff1(iEf,iang1)
         end do  

      else    ! depending on location, only wind 1 absorbs or both do
         r1 = sqrt( x(l)**2 + y(l)**2+z(l)**2 )
         r2 = sqrt( (D-x(l))**2 + y(l)**2+z(l)**2 )
         dr0 = sh_width*min(r1,r2)
         x1 = x(l)+dr0*cinc
         y1 = y(l)+dr0*sinc
         z1 = z(l)

         ! The physical CD is the positive-x hyperbola branch.
         cond = logical( x1 <= a_CD*sqrt(1.d0 + (y1/b_CD)**2 + (z1/b_CD)**2) )

         if( cond ) then ! moves outside of WCR, only wind 1 absorbs
            call unit_vector_diff(Xemi,Xst1,3,XEe) ! XEe vector from the star to the emitter
            ang1 = acos(DOT_PRODUCT(Xobs,XEe))! angle between Xobs and XEe
            iang1 = int(m_ang_ff*abs(ang1-ang_min)/pi)+1 
            r1 = sqrt( x(l)**2 + y(l)**2 + z(l)**2 )
            do i = 1, mnu      
               iEf = min(int(log(max(1.d0,nu(i)/nu_min_ff))/log(nu_int_ff)) + 1,mE_ff)
               tauff(l,i) = (Rst1/r1)**3*tau_ff1(iEf,iang1)
            end do  

         else ! wind 2 also absorbs

            call crossing_point(sinc,cinc,a_CD,b_CD,x(l),y(l),z(l),r_cross)
            ! from XEe to r_cross wind2 absorbs, whereas wind1 from r_cross to "infinity" 

            call unit_vector_diff(r_cross,Xst1,3,XEe) ! XEe = unit vector from the star1 to the crossing point
            ang1 = acos(DOT_PRODUCT(Xobs,XEe))   ! angle between Xobs and XEe
            iang1 = int(m_ang_ff*abs(ang1-ang_min)/pi)+1
            r1 = sqrt( r_cross(1)**2 + r_cross(2)**2 + r_cross(3)**2 )

            call unit_vector_diff(Xemi,Xst2,3,XEe)   ! XEe = unit vector from the star2 to the emitter
            ang2_i = acos(DOT_PRODUCT(Xobs,XEe))! angle between Xobs and XEe_i
            iang2_i = int(m_ang_ff*abs(ang2_i-ang_min)/pi)+1
            r2_i = sqrt( (D-x(l))**2 + y(l)**2 + z(l)**2 )

            call unit_vector_diff(r_cross,Xst2,3,XEe) ! XEe = unit vector from the star2 to the crossing point
            ang2_f = acos(DOT_PRODUCT(Xobs,XEe)) ! angle between Xobs and XEe_f
            iang2_f = int(m_ang_ff*abs(ang2_f-ang_min)/pi)+1
            r2_f = sqrt( (D-r_cross(1))**2 + r_cross(2)**2 + r_cross(3)**2 )

            do i=1,mnu     
               iEf = min(int(log(max(1.d0,nu(i)/nu_min_ff))/log(nu_int_ff)) + 1,mE_ff)
               tauff(l,i) = ( (Rst2/r2_i)**3*tau_ff2(iEf,iang2_i) &
                  - (Rst2/r2_f)**3*tau_ff2(iEf,iang2_f) ) &
               + (Rst1/r1)**3*tau_ff1(iEf,iang1) 
            end do

         end if

      end if 

   end do

   deallocate( tau_ff1, tau_ff2 )

end subroutine total_tau_calc

!==========================================================================

!==========================================================================
! Same as total_tau_calc but for specific frequencies

subroutine total_tau_data_calc(D,incli,tauff)
   use global
   implicit none
   real(dp), intent(in) :: D,incli
   real(dp), intent(out) :: tauff(ml,mnu_data)
   real(dp) :: sinc,cinc,r1,r2,r2_i,r2_f,dr0
   real(dp) :: ang_int,ang_min,ang1,ang2_i,ang2_f
   real(dp) :: a_CD,b_CD,x1,y1,z1,r_cross(3),Xobs(3)
   real(dp), parameter :: sh_width = 0.05 ! 5%
   integer :: iang1,iang2_i,iang2_f
   logical :: cond 
   allocate( tau_ff1(mnu_data,m_ang_ff), tau_ff2(mnu_data,m_ang_ff) )

   sinc = sin(incli)
   cinc = cos(incli)

   open (288,file=output_file('opac_ff_data1.dat'))
   open (289,file=output_file('opac_ff_data2.dat'))
   do i=1,mnu_data 
      read(288,*) (tau_ff1(i,j), j=1,m_ang_ff)
      read(289,*) (tau_ff2(i,j), j=1,m_ang_ff)
   end do
   close(288)
   close(289)

   ang_int = pi/float(m_ang_ff)
   ang_min = ang_int/2.d0

   Xobs = (/ cinc,sinc,0.d0 /)
   Xst1 = (/ 0.d0, 0.d0, 0.d0 /)
   Xst2 = (/ D, 0.d0, 0.d0 /)

   if (verbose) print *, 'Calculate the total absorption for synchrotron radiation'

   ! To know which wind absorbs (1 or 2) we check if the position is within 
   ! the (approx.) surface of the CD

   open(191,file=output_file('a_CD.dat'))
   ! (x/a)**2 = 1 + (y/b)**2 + (z/b)**2
   read(191,*) a_CD, b_CD
   close(191)

   tauff = 0.d0

   do l = 1,ml        ! Position along the CD

      Xemi = (/ x(l),y(l),z(l) /) 

      if( incli < theta_inf ) then  ! only wind 2 absorbs
         call unit_vector_diff(Xemi,Xst2,3,XEe) ! XEe vector from the star to the emitter
         ang2_i = acos(DOT_PRODUCT(Xobs,XEe))! angle between Xobs and XEe
         iang2_i = int(m_ang_ff*abs(ang2_i-ang_min)/pi)+1 
         r2 = sqrt( (D-x(l))**2 + y(l)**2 + z(l)**2 )
         do i = 1, mnu_data      
            tauff(l,i) = (Rst2/r2)**3*tau_ff2(i,iang2_i)
         end do  
      
      elseif( incli > pi - theta_inf ) then ! only wind 1 absorbs
         call unit_vector_diff(Xemi,Xst1,3,XEe) ! XEe vector from the star to the emitter
         ang1 = acos(DOT_PRODUCT(Xobs,XEe))! angle between Xobs and XEe
         iang1 = int(m_ang_ff*abs(ang1-ang_min)/pi)+1 
         r1 = sqrt( x(l)**2 + y(l)**2+z(l)**2 )
         do i = 1, mnu_data
            tauff(l,i) = (Rst1/r1)**3*tau_ff1(i,iang1)
         end do  

      else    ! depending on location, only wind 1 absorbs or both do
         r1 = sqrt( x(l)**2 + y(l)**2+z(l)**2 )
         r2 = sqrt( (D-x(l))**2 + y(l)**2+z(l)**2 )
         dr0 = sh_width*min(r1,r2)
         x1 = x(l)+dr0*cinc
         y1 = y(l)+dr0*sinc
         z1 = z(l)

         ! The physical CD is the positive-x hyperbola branch.
         cond = logical( x1 <= a_CD*sqrt(1.d0 + (y1/b_CD)**2 + (z1/b_CD)**2) )

         if( cond ) then ! moves outside of WCR, only wind 1 absorbs
            call unit_vector_diff(Xemi,Xst1,3,XEe) ! XEe vector from the star to the emitter
            ang1 = acos(DOT_PRODUCT(Xobs,XEe))! angle between Xobs and XEe
            iang1 = int(m_ang_ff*abs(ang1-ang_min)/pi)+1 
            r1 = sqrt( x(l)**2 + y(l)**2 + z(l)**2 )
            do i = 1, mnu_data    
               tauff(l,i) = (Rst1/r1)**3*tau_ff1(i,iang1)
            end do  

         else ! wind 2 also absorbs

            call crossing_point(sinc,cinc,a_CD,b_CD,x(l),y(l),z(l),r_cross)
            ! from XEe to r_cross wind2 absorbs, whereas wind1 from r_cross to "infinity" 

            call unit_vector_diff(r_cross,Xst1,3,XEe) ! XEe = unit vector from the star1 to the crossing point
            ang1 = acos(DOT_PRODUCT(Xobs,XEe))   ! angle between Xobs and XEe
            iang1 = int(m_ang_ff*abs(ang1-ang_min)/pi)+1
            r1 = sqrt( r_cross(1)**2 + r_cross(2)**2 + r_cross(3)**2 )

            call unit_vector_diff(Xemi,Xst2,3,XEe)   ! XEe = unit vector from the star2 to the emitter
            ang2_i = acos(DOT_PRODUCT(Xobs,XEe))! angle between Xobs and XEe_i
            iang2_i = int(m_ang_ff*abs(ang2_i-ang_min)/pi)+1
            r2_i = sqrt( (D-x(l))**2 + y(l)**2 + z(l)**2 )

            call unit_vector_diff(r_cross,Xst2,3,XEe) ! XEe = unit vector from the star2 to the crossing point
            ang2_f = acos(DOT_PRODUCT(Xobs,XEe)) ! angle between Xobs and XEe_f
            iang2_f = int(m_ang_ff*abs(ang2_f-ang_min)/pi)+1
            r2_f = sqrt( (D-r_cross(1))**2 + r_cross(2)**2 + r_cross(3)**2 )

            do i=1,mnu_data
               tauff(l,i) = ( (Rst2/r2_i)**3*tau_ff2(i,iang2_i) &
                  - (Rst2/r2_f)**3*tau_ff2(i,iang2_f) ) &
               + (Rst1/r1)**3*tau_ff1(i,iang1) 
            end do

         end if

      end if 

   end do

   deallocate( tau_ff1, tau_ff2 )

end subroutine total_tau_data_calc

!==========================================================================
subroutine crossing_point(si,ci,aCD,bCD,xl,yl,zl,r)

   ! Find the intersection between the photon trajectory (line) and the 
   ! WCR (approx. surface) using the parametrizations:
   ! x(t) = xl + ci*t , y(t) = yl + si*t, z(t) = zl   [line]

   ! (x/aCD)**2 = 1 + (y/bCD)**2 + (z/bCD)**2  [hyperboloid]

   ! Solution: find "tc" such that x(tc)=x_sur, y(tc)=y_sur, z(tc)=z_sur

   implicit none
   integer, parameter :: dp = kind(1.d0)
   real(dp), intent(in) :: si,ci,aCD,bCD,xl,yl,zl
   real(dp), intent(out) :: r(3)
   real(dp) :: disc,a,b,c,tc,xc,yc,zc

   a = (ci/aCD)**2 - (si/bCD)**2
   b = 2.d0 * ( xl*ci/(aCD**2) - yl*si/(bCD**2) )
   c = (xl/aCD)**2 - (yl/bCD)**2 - (zl/bCD)**2 - 1.d0
   disc = b**2 - 4.d0*a*c

   if (disc < 0.d0) then
      print *, 'ABSURD! Bad parameters in CD_abssyn.f90'
      print *, aCD, bCD, ci, si, xl, yl, zl
      stop
   end if
   tc = (-b - sqrt(disc))/(2.d0*a) ! The conventions used correspond to the negative solution
   xc = xl + ci*tc 
   yc = yl + si*tc
   zc = zl

   r = (/ xc, yc, zc /)

   RETURN

end subroutine crossing_point

module BB
  use global
  implicit none
  
  public :: BB_run
  
contains
  
  !============================================================================
  ! BB_run
  !============================================================================
  ! Routine to calculate the black body emission from the stars
  !============================================================================
  
  subroutine BB_run  
    use global
    implicit none
    integer, parameter :: mE_ph = 200
    real(dp) :: Eph_0, Eph_min, Eph_max, Eph_int, Eph(mE_ph)
    real(dp) :: L1_BB, L2_BB, L_BB(mE_ph)

    print *, 'Calculate the thermal emission from the stars'

    Eph_0 = 3.d0*k*T1
    Eph_min = Eph_0*1.d-6
    Eph_max = Eph_0*10.d0
    
    call vector_log(Eph_min,Eph_max,Eph_int,Eph)

    open(20,file=output_file('LBB.dat'))
    do i=1,mE_ph
       L1_BB = Rst1**2 * nBB(Eph(i),T1)   
       L2_BB = Rst2**2 * nBB(Eph(i),T2)
       L_BB(i) = L1_BB + L2_BB
       write(20,*) Eph(i)/eV, Eph(i)*L_BB(i), Eph(i)*L1_BB, Eph(i)*L2_BB
    enddo
    close(20)


  END subroutine BB_run

  !============================================================================
  ! BB_run (END)
  !============================================================================
  
 
END module BB
 

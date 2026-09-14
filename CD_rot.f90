module CD_rot
  use global
  implicit none
  
  public :: CD_rot_run

contains

  !============================================================================
  ! CD_rot_run
  !============================================================================
  ! Rotate position vectors of the CD wrt the axis that joins both stars (X-axis)
  !============================================================================

  subroutine CD_rot_run(ang_phi)  
    real(dp), intent(in) :: ang_phi
    real(dp) :: x_no_rot,y_no_rot
   
    if (verbose) print *, 'Rotate the CD coordinates'
    open(13,file=output_file('CD.dat')) ! This contains only the upper half
    do l = 1, ml_2
      read(13, *) x_no_rot, y_no_rot 
      call rotate_vector(x_no_rot,y_no_rot,ang_phi,x(l),y(l),z(l))
      call rotate_vector(x_no_rot,-y_no_rot,ang_phi,x(l+ml_2),y(l+ml_2),z(l+ml_2))
    enddo
    close(13)
  END subroutine CD_rot_run
  
  !===================================================  
  
  Subroutine rotate_vector(Vx,Vy,ang_phi,VRx,VRy,VRz)
    ! Rotate vector V in an angle ang_phiwrt X-axis
    real(dp), intent(in):: Vx,Vy,ang_phi
    real(dp), intent(out) :: VRx,VRy,VRz
    VRx = Vx              ! unaltered
    VRy = cos(ang_phi)*Vy ! Vz=0
    VRz = sin(ang_phi)*Vy ! Vz=0
  END Subroutine rotate_vector

!===================================================  
  
END module CD_rot
  
    

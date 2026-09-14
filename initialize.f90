MODULE initialize
    use global
    !use constants
    !use controls
    implicit none
    
    contains 

    !============================================================================
    ! initialize_allocate (START)
    !============================================================================
    subroutine initialize_allocate() !eta_B_value, frac_NT_value)
    !============================================================================
    ! allocate and initialize vectors 
    !============================================================================
        !real(dp), intent(in) :: eta_B_value, frac_NT_value
        !logical :: fexist

        allocate( x(ml), y(ml), z(ml) )
        allocate( P(ml), rho(ml), v(ml), T(ml), B(ml), B_par(ml), n_a(ml), H_sh(ml), dVol(ml) ) 
        allocate( eta_acc(ml_2) )
        allocate( E_e(mE_e), E_p(mE_p) )
        allocate( nu(mnu), E_sy(mnu) )

        ! Initialize values given by free parameters:
        !eta_B1 = eta_B_value    !Update the value of eta_B (default value defined in system_parameters.f90)
        !eta_B2 = eta_B_value    !Update the value of eta_B2 (default value defined in system_parameters.f90)
        !frac_NT = frac_NT_value !Update the value of Lj (default value defined in system_parameters.f90)

        ! Global variables that are tied to the free parameters
        frac_B1 = sqrt(eta_B1)          ! B = frac*B_equipartition
        frac_B2 = sqrt(eta_B2)          ! B = frac*B_equipartition
        inj_e = frac_NT * K_ep          ! Fraction L_NT_e/L_inj
        inj_p = frac_NT * (1.d0-K_ep)   ! Fraction L_NT_p/L_inj

        !write(*,'(A,Es9.2)') 'eta_B1 = ', eta_B1
        !write(*,'(A,Es9.2)') 'frac_NT = ', frac_NT

    end subroutine initialize_allocate
    !============================================================================
    ! initialize_allocate (END)
    !============================================================================
    

    !============================================================================
    ! initialize_deallocate (START)
    !============================================================================
    subroutine initialize_deallocate
    !============================================================================
    ! deallocate vectors 
    !============================================================================

        deallocate( x, y, z )
        deallocate( P, rho, v, T, B, B_par, n_a, H_sh, dVol )
        deallocate( eta_acc )
        deallocate( E_e, E_p )
        deallocate( nu, E_sy )

    end subroutine initialize_deallocate
    !============================================================================
    ! initialize_deallocate (END)
    !============================================================================


END MODULE initialize

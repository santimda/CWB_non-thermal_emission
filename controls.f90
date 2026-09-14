module controls

implicit none
public

!===================================================
! Select the programs that will be executed

logical, parameter :: off = .FALSE.
logical, parameter :: on = .TRUE.

logical, parameter :: verbose = off ! print info or not

! Meta-task for all SED
logical, parameter :: sed_all_calc = on

! Meta-task for the radio SED
logical, parameter :: radio_calc = off

logical, parameter :: single_epoch = on
logical, parameter :: CD_calc = on
logical, parameter :: thermo = on
logical, parameter :: absgg = on
logical, parameter :: absff = on ! also calculates emi_ff(wind)
logical, parameter :: BBrad = off
logical, parameter :: Xrad = on  ! X-ray emission from the WCR
logical, parameter :: rotate = on

logical, parameter :: calc_N_bin = on

logical, parameter :: electrons = on
logical, parameter :: dist1e = on
logical, parameter :: dist2e = on
logical, parameter :: rad_syn1 = on
logical, parameter :: rad_syn2 = on
logical, parameter :: rad_ic1 = on
logical, parameter :: rad_ic2 = on
logical, parameter :: rad_br1 = on
logical, parameter :: rad_br2 = on
logical, parameter :: abssyn = on
logical, parameter :: maps = on

logical, parameter :: protons = on
logical, parameter :: dist1p = on
logical, parameter :: dist2p = on
logical, parameter :: rad_pp1 = on
logical, parameter :: rad_pp2 = on

logical, parameter :: write_e = on ! flux from leptons in the WCR
logical, parameter :: write_p = on ! flux from protons in the WCR

logical, parameter :: convolve = on

!===================================================

END module controls

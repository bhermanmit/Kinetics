module kinetics_solver 

!-module references

!-module options

  implicit none
  private
  public :: kinetics_execute 

!-module variables

contains

!===============================================================================
! KINETICS_EXECUTE
!===============================================================================

  subroutine kinetics_execute(inner_solver)

!---arguments

    external :: inner_solver

!---begin execution

    call init_data()

  end subroutine kinetics_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

!---external references

    use cmfd_header,        only: compute_core_power
    use constants,          only: ONE
    use global,             only: cmfd, geometry, material, kine
    use kinetics_operator,  only: init_K_operator

!---local variables

    real(8) :: pow

!---begin execution

    ! normalize initial power to unity and set initial power
    pow = compute_core_power(cmfd, size(cmfd%phi), geometry, material)
    cmfd % phi = cmfd % phi * ONE / pow

    ! set up matrices
    call init_K_operator(kine)

  end subroutine init_data

!===============================================================================
! EXECUTE_KINETICS_ITER  in the main kinetics iteration routine 
!                         for the cmfd calculation
!===============================================================================

  subroutine execute_kinetics_iter(inner_solver)

!---external references

!---arguments

    external :: inner_solver

!---local variables

!--begin execution

  end subroutine execute_kinetics_iter 

!==============================================================================
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize data objects

  end subroutine finalize

end module kinetics_solver 

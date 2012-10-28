module cmfd_execute

!-module options

  implicit none
  private
  public execute_cmfd

contains

!==============================================================================
! EXECUTE_CMFD
!==============================================================================

  subroutine execute_cmfd()

!---external references

    use error,       only: fatal_error
    use global,      only: solver_type, message, run_kinetics
    use kinetics_solver,  only: kinetics_execute
    use output,      only: header
    use power_iter,  only: power_execute

!---begin execution

    ! print to screen
    call header('MULTIGROUP DIFFUSION', level=1)

    ! execute solver
    call power_execute() 

    ! call kinetics
    if (run_kinetics) call kinetics_execute()

  end subroutine execute_cmfd

end module cmfd_execute

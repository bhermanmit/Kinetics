module cmfd_execute

!-module options

  implicit none
  private
  public execute_cmfd

contains

!===============================================================================
! EXECUTE_CMFD
!===============================================================================

  subroutine execute_cmfd()

!---external references

    use constants,        only: ONE
    use error,            only: fatal_error
    use global,           only: solver_type, message, mode,                    &
                                adjoint
    use kinetics_solver,  only: kinetics_execute
    use output,           only: header
    use point_kinetics,   only: run_pkinetics 
    use power_iter,       only: power_execute

!---begin execution

    ! run calculation
    select case(trim(mode))

      case('static')
        call header('MULTIGROUP FORWARD DIFFUSION', level=1)
        call power_execute()

      case('kinetics')

        ! we always need a forward solution
        call header('MULTIGROUP KINETICS', level=1)
        call header('Forward Solution', level=2)
        call power_execute()

        ! now run kinetics
        call header('Kinetics Solution', level=2)
        call kinetics_execute()

      case('point_kinetics')

        call run_pkinetics()

      case DEFAULT
        message = 'Calculation Mode not Supported!'
        call fatal_error()

    end select

  end subroutine execute_cmfd

end module cmfd_execute

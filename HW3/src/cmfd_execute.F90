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
    use math,        only: csr_jacobi, csr_gauss_seidel
    use output,      only: header
    use power_iter,  only: power_execute

!---begin execution

    ! print to screen
    call header('MULTIGROUP DIFFUSION', level=1)

    ! execute solver
    select case (trim(solver_type))

      case('jacobi')
        call power_execute(csr_jacobi) 
      case('gauss')
        call power_execute(csr_gauss_seidel)
      case DEFAULT      
        message = "Solver type does not exist!"
        call fatal_error()

    end select

    ! call kinetics
    if (run_kinetics) call kinetics_execute(csr_gauss_seidel)

  end subroutine execute_cmfd

end module cmfd_execute

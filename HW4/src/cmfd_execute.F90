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

    use error,            only: fatal_error
    use global,           only: solver_type, message, run_kinetics, adjoint
    use kinetics_solver,  only: kinetics_execute
    use math,             only: csr_jacobi, csr_gauss_seidel
    use output,           only: header
    use point_kinetics,   only: generate_pkparams, run_pkes
    use power_iter,       only: power_execute

!---begin execution

    ! print to screen
    call header('MULTIGROUP DIFFUSION', level=1)

    ! run adjoint
    if (adjoint == 'math' .or. adjoint == 'physical') then
      call power_execute(csr_gauss_seidel)
      adjoint = ''
    end if

    ! run forward calculation
    call power_execute(csr_gauss_seidel)

    ! generate point kinetics parameters
!   call generate_pkparams()

    ! call kinetics
    if (run_kinetics) call kinetics_execute(csr_gauss_seidel)
    if (run_kinetics) call run_pkes()
  end subroutine execute_cmfd

end module cmfd_execute

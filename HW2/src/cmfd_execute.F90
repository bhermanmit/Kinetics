!==============================================================================!
! MODULE: cmfd_execute
!
!> @author Bryan Herman
!>
!> @brief Routine for running the eigenvalue solve
!==============================================================================!

module cmfd_execute

  use global
  use power_solver, only: cmfd_power_execute
  use slepc_solver, only: cmfd_slepc_execute
  use snes_solver,  only: cmfd_snes_execute
 
  implicit none

contains

!==============================================================================
! EXECUTE_CMFD
!==============================================================================

  subroutine execute_cmfd()

    integer :: ierr  ! petsc error code

    ! execute solver
    select case (trim(solver_type))

      case('power')
        call cmfd_power_execute()
      case('slepc')
        call cmfd_slepc_execute()
      case('snes')
        call cmfd_snes_execute()
      case DEFAULT
        call cmfd_power_execute()

    end select

  end subroutine execute_cmfd

end module cmfd_execute

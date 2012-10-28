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
    use general_pkes,     only: run_gpkes, generate_gpkparams, generate_gpke_shapes
    use global,           only: solver_type, message, mode,                    &
                                adjoint, pke_run, weight, cmfd
    use kinetics_solver,  only: kinetics_execute
    use math,             only: csr_jacobi, csr_gauss_seidel
    use output,           only: header
    use point_kinetics,   only: generate_pkparams, run_pkes, generate_pke_shapes
    use power_iter,       only: power_execute

!---begin execution

    ! run calculation
    select case(trim(mode))

      case('static')
        call header('MULTIGROUP FORWARD DIFFUSION', level=1)
        call power_execute(csr_gauss_seidel,'none')

      case('adjoint')
        call header('MULTIGROUP ADJOINT DIFFUSION', level=1)
        call power_execute(csr_gauss_seidel,trim(adjoint))

      case('kinetics')

        ! we always need a forward solution
        call header('MULTIGROUP KINETICS', level=1)
        call header('Forward Solution', level=2)
        call power_execute(csr_gauss_seidel,'none')

        ! we always generate pke parameters, check for weighting function
        if (trim(weight) == 'adjoint') then
          call header('Adjoint Solution', level=2)
          call power_execute(csr_gauss_seidel,trim(adjoint))
        else
          cmfd % phi_adj = ONE
        end if

        ! now run kinetics
        call header('Kinetics Solution', level=2)
        call kinetics_execute(csr_gauss_seidel)

        ! run point kinetics based on these exact generated values
        if (trim(pke_run) == 'classical') then 
          call header('Point Kinetics w/ Kinetics Parameters',level=2)
          call run_pkes()
        else if (trim(pke_run) == 'general') then
          call header('General Kinetics w/ Kinetics Parameters',level=2)
          call run_gpkes()
        end if

      case('point_kinetics')
        call header('Point Kinetics w/ Static Shape Functions', level=1)
        call header('Running Steady State Forward Solution', level=2)
        call power_execute(csr_gauss_seidel,'none')
        call header('Generating Shape Functions', level=2)
        call generate_pke_shapes() 
        call header('Generating Point Kinetics Parameters', level=2)
        call generate_pkparams()
        call header('Running Point Kinetics', level=2)
        call run_pkes()

      case('general_point_kinetics')
        call header('General Point Kinetics w/ Static Shape Functions', level=1)
        call header('Running Steady State Forward Solution', level=2)
        call power_execute(csr_gauss_seidel,'none')
        call header('Generating Shape Functions', level=2)
        call generate_gpke_shapes()
        call header('Generating General Point Kinetics Parameters', level=2)
        call generate_gpkparams()
        call header('Running General Point Kinetics', level=2)
        call run_gpkes()

      case DEFAULT
        message = 'Calculation Mode not Supported!'
        call fatal_error()

    end select

  end subroutine execute_cmfd

end module cmfd_execute

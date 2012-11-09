module implicit_euler

!-module references

  use global,  only: mpi_err

!-module options

  implicit none
  private
  public :: execute_ie1

!-module external references

# include "finclude/petsc.h90"

!-data for solver

  Mat :: A
  Vec :: yinit
  Vec :: dydtinit
  Vec :: k1
  Vec :: k2
  Vec :: k3
  Vec :: k4
  Vec :: rhs
  Vec :: err

!-solver objects

  KSP :: ksp
  PC  :: pc

contains

!===============================================================================
! EXECUTE_IE1
!===============================================================================

  subroutine execute_ie1(y, A, coefmat, n)

!---references

    use constants,  only: ZERO
    use global,     only: time, dt

!---arguments

    integer :: n 

    Vec :: y
    Mat :: A 

    external coefmat 

!---local variables

    real(8) :: h
    real(8) :: t
    real(8), pointer :: yptr(:)

!---begin execution

    ! initialize time steps
    t = ZERO

    ! initialize solver
    call init_solver(n)

    ! set time step
    h = dt

    ! begin loop
    do while (t <= time)

      ! solve for next time step values
      call solve_ts(y, A, n, t, h, coefmat)

      call VecGetArrayF90(y, yptr, mpi_err)
      print *, 'POWER:', yptr(1), 'TIME:', t
      call VecRestoreArrayF90(y, yptr, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif
    end do

  end subroutine execute_ie1

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver(n)

!---references

    use constants,  only: ZERO
    use global,     only: itol

!---arguments

    integer :: n

!---begin execution

    ! create rhs vector
    call VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, rhs, mpi_err)
    call VecSet(rhs, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create ksp object
    call KSPCreate(PETSC_COMM_WORLD, ksp, mpi_err)
    call KSPSetTolerances(ksp, itol, itol, PETSC_DEFAULT_DOUBLE_PRECISION,&
                          PETSC_DEFAULT_INTEGER, mpi_err)
    call KSPSetType(ksp, KSPGMRES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! set up preconditioner
    call KSPGetPC(ksp, pc, mpi_err)
    call PCSetType(pc, PCILU, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine

!===============================================================================
! SOLVE_TS
!===============================================================================

  subroutine solve_ts(y, A, n, t, h, coefmat)

!---references

    use constants,  only: ZERO, ONE
    use error,      only: fatal_error
    use global,     only: message

!---arguments

    integer, intent(in)     :: n        ! size of problem
    real(8), intent(inout)  :: t        ! t value for step
    real(8), intent(in)     :: h        ! trial time step

    Mat :: A 
    Vec :: y

    external coefmat 

!---local variables

    integer :: istep  ! time step trial iteration
    integer :: irow   ! row iteration counter
    real(8) :: tinit  ! initial time step
    real(8) :: errmax ! maximum error ratio

!---begin execution

    ! calcualate new timestep
    t = t + h

    ! evaluate jacobian at beginning of time step
    call coefmat(t,h,y,A,n)

    ! set the operator and finish setting up
    call KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif DEBUG
    call KSPSetFromOptions(ksp, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif
    call KSPSetUp(ksp, mpi_err)
#   ifdef DEBUG
     CHKERRQ(mpi_err)
#   endif

    ! copy right hand side
    call VecCopy(y, rhs, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! solve for k1
    call KSPSolve(ksp, rhs, y, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine solve_ts

end module implicit_euler

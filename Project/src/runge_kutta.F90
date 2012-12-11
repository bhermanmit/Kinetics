module runge_kutta

!-module references

  use global,  only: mpi_err

!-module options

  implicit none
  private
  public :: execute_rk4

!-module external references

# include "finclude/petsc.h90"

!-runge-kutta time step constants

  integer, parameter :: MAXTRY  = 40
  real(8), parameter :: SAFETY  = 0.9_8
  real(8), parameter :: GROW    = 1.5_8
  real(8), parameter :: PGROW   = -0.25_8
  real(8), parameter :: SHRNK   = 0.5_8
  real(8), parameter :: PSHRNK  = -1.0_8/3.0_8
  real(8), parameter :: ERRCON  = 0.1296_8

!-runge-kutta constants

  real(8), parameter :: GAM = 1.0_8/2.0_8
  real(8), parameter :: A21 = 2.0_8
  real(8), parameter :: A31 = 48.0_8/25.0_8
  real(8), parameter :: A32 = 6.0_8/25.0_8
  real(8), parameter :: C21 = -8.0_8
  real(8), parameter :: C31 = 372.0_8/25.0_8
  real(8), parameter :: C32 = 12.0_8/5.0_8
  real(8), parameter :: C41 = -112.0_8/125.0_8
  real(8), parameter :: C42 = -54.0_8/125.0_8
  real(8), parameter :: C43 = -2.0_8/5.0_8
  real(8), parameter :: M1  = 19.0_8/9.0_8
  real(8), parameter :: M2  = 1.0_8/2.0_8
  real(8), parameter :: M3  = 25.0_8/108.0_8
  real(8), parameter :: M4  = 125.0_8/108.0_8
  real(8), parameter :: B1  = 1.0_8/2.0_8
  real(8), parameter :: B2  = -3.0_8/2.0_8
  real(8), parameter :: B3  = 121.0_8/50.0_8
  real(8), parameter :: B4  = 29.0_8/250.0_8
  real(8), parameter :: D2  = 1.0_8
  real(8), parameter :: D3  = 3.0_8/5.0_8
  real(8), parameter :: E1  = 17.0_8/54.0_8
  real(8), parameter :: E2  = 7.0_8/36.0_8
  real(8), parameter :: E3  = 0.0_8
  real(8), parameter :: E4  = 125.0_8/108.0_8

!-data for solver

  Mat :: A
  Vec :: yinit
  Vec :: yscal
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
  PetscViewer :: viewer
contains

!===============================================================================
! EXECUTE_RK4
!===============================================================================

  subroutine execute_rk4(y, dfdy, dfdt, dydt, derivs, jacobn, n, post_timestep, Cfactopt)

!---references

    use constants,  only: ZERO
    use global,     only: time, dt

!---arguments

    integer :: n 

    real(8), optional :: Cfactopt

    Vec :: y
    Mat :: dfdy
    Vec :: dfdt
    Vec :: dydt 

    external derivs
    external jacobn 
    external post_timestep

!---local variables

    integer :: i
    real(8) :: hdid
    real(8) :: hnext
    real(8) :: htry
    real(8) :: eps
    real(8) :: t
    real(8) :: Cfact
    real(8), pointer :: yptr(:)

!---begin execution

    ! initialize time steps
    htry = dt
    hdid = ZERO
    hnext = ZERO
    t = ZERO
    i = 0

    ! set tolerance
    eps = 1.e-4_8

    ! get constant
    if (present(Cfactopt)) then
      Cfact = Cfactopt
    else
      Cfact = ZERO
    end if

    ! initialize solver
    call init_solver(n, dfdy)

    ! begin loop
    do while (t <= time)

      ! solve for next time step values
      call solve_ts(y, dfdy, dfdt, dydt, n, t, htry, eps, hdid, hnext, derivs, jacobn, Cfact)

      ! set next time step
      htry = hnext

      i = i + 1

      ! post timestep routine
      call post_timestep(t, y, hdid, i)

    end do

    ! clean up objects
    call destroy_objects()

  end subroutine execute_rk4

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver(n, dfdy)

!---references

    use constants,  only: ZERO
    use global,     only: itol

!---arguments

    integer :: n
    Mat     :: dfdy

!---begin execution

    ! create y initial vector
    call VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, yinit, mpi_err)
    call VecSet(yinit, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create dydt initial vector
    call VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, dydtinit, mpi_err)
    call VecSet(dydtinit, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create k vectors
    call VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, k1, mpi_err)
    call VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, k2, mpi_err)
    call VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, k3, mpi_err)
    call VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, k4, mpi_err)
    call VecSet(k1, ZERO, mpi_err)
    call VecSet(k2, ZERO, mpi_err)
    call VecSet(k3, ZERO, mpi_err)
    call VecSet(k4, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create rhs vector
    call VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, rhs, mpi_err)
    call VecSet(rhs, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create error vector
    call VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, err, mpi_err)
    call VecSet(err, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create ksp object
    call KSPCreate(PETSC_COMM_WORLD, ksp, mpi_err)
    call KSPSetTolerances(ksp, 1.e-6_8, 1.e-6_8, PETSC_DEFAULT_DOUBLE_PRECISION,&
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

  subroutine solve_ts(y, dfdy, dfdt, dydt, n, t, htry, eps, hdid, hnext, derivs, jacobn, Cfact)

!---references

    use constants,  only: ZERO, ONE
    use error,      only: fatal_error
    use global,     only: message, var_ts

!---arguments

    integer, intent(in)     :: n        ! size of problem
    real(8), intent(inout)  :: t        ! t value for step
    real(8), intent(in)     :: htry     ! trial time step
    real(8), intent(in)     :: eps      ! tolerance on solution
    real(8), intent(out)    :: hdid     ! the time step it actually took
    real(8), intent(out)    :: hnext    ! the next estimated time step
    real(8), intent(in)     :: Cfact    ! C factor

    Mat :: dfdy
    Vec :: y
    Vec :: dfdt
    Vec :: dydt

    external derivs
    external jacobn

!---local variables

    integer :: istep  ! time step trial iteration
    integer :: irow   ! row iteration counter
    real(8) :: tinit  ! initial time step
    real(8) :: h      ! current time step being evaluated
    real(8) :: errmax ! maximum error ratio

    real(8), pointer :: yptr(:)
    real(8), pointer :: errptr(:)
    real(8), pointer :: yinitptr(:)
    real(8), pointer :: k1ptr(:)
    real(8), pointer :: k2ptr(:)
    real(8), pointer :: k3ptr(:)
    real(8), pointer :: k4ptr(:)
    real(8), pointer :: rhsptr(:)
    real(8), pointer :: dydtptr(:)
    real(8), pointer :: dfdtptr(:)
    real(8), pointer :: yscalptr(:)

!---begin execution

    ! compute scaling vector
    call VecDuplicate(y, yscal, mpi_err)
    call VecCopy(y, yscal, mpi_err)
    call VecGetArrayF90(yscal, yscalptr, mpi_err)
    where (abs(yscalptr) < Cfact)
      yscalptr = Cfact
    end where
    call VecRestoreArrayF90(yscal, yscalptr, mpi_err)

    ! evaluate jacobian at beginning of time step
!   call jacobn(t,y,dfdy,dfdt,n)

    ! call derivs to get dydt vector
    call derivs(t, y, dydt, n)

    ! save all important vectors
    call VecCopy(y, yinit, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif
    call VecCopy(dydt, dydtinit, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif
    tinit = t

    ! set up initial trial time step
    h = htry

    ! begin loop over trial time steps
    do istep = 1, MAXTRY

      ! create coefficient matrix
!     call MatConvert(dfdy, MATSAME, MAT_INITIAL_MATRIX, A, mpi_err)
      ! evaluate jacobian at beginning of time step
      call jacobn(t,y,dfdy,dfdt,n, ONE/(GAM*h))
      ! multiply values by -1
!     call MatScale(A, -1.0_8, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! modify diagonal
!     call MatShift(A, ONE/(GAM*h), mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! set the operator and finish setting up
      call KSPSetOperators(ksp, dfdy, dfdy, SAME_NONZERO_PATTERN, mpi_err)
!  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'jacnew.bin', &
!         FILE_MODE_WRITE, viewer, mpi_err)
!  call MatView(dfdy, viewer, mpi_err)
!  call PetscViewerDestroy(viewer, mpi_err)

#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif
      call KSPSetFromOptions(ksp, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif
      call KSPSetUp(ksp, mpi_err)
#     ifdef DEBUG
       CHKERRQ(mpi_err)
#     endif

      ! set up right hand side for k1
      call VecGetArrayF90(rhs, rhsptr, mpi_err)
      call VecGetArrayF90(dydtinit, dydtptr, mpi_err)
      call VecGetArrayF90(dfdt, dfdtptr, mpi_err)
      rhsptr = dydtptr + h*B1*dfdtptr
      call VecRestoreArrayF90(rhs, rhsptr, mpi_err)
      call VecRestoreArrayF90(dydtinit, dydtptr, mpi_err)
      call VecRestoreArrayF90(dfdt, dfdtptr, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! solve for k1
      call KSPSolve(ksp, rhs, k1, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! update y solution vector and change t location
      call VecGetArrayF90(y, yptr, mpi_err)
      call VecGetArrayF90(yinit, yinitptr, mpi_err)
      call VecGetArrayF90(k1, k1ptr, mpi_err)
      yptr = yinitptr + A21*k1ptr
      t = tinit + D2*h
      call VecRestoreArrayF90(y, yptr, mpi_err)
      call VecRestoreArrayF90(yinit, yinitptr, mpi_err)
      call VecRestoreArrayF90(k1, k1ptr, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! recompute derivative vector with new (t,y)
      call derivs(t, y, dydt, n) 

      ! set up right hand side for k2
      call VecGetArrayF90(rhs, rhsptr, mpi_err)
      call VecGetArrayF90(dydt, dydtptr, mpi_err)
      call VecGetArrayF90(dfdt, dfdtptr, mpi_err)
      call VecGetArrayF90(k1, k1ptr, mpi_err)
      rhsptr = dydtptr + h*B2*dfdtptr + C21*k1ptr/h
      call VecRestoreArrayF90(rhs, rhsptr, mpi_err)
      call VecRestoreArrayF90(dydt, dydtptr, mpi_err)
      call VecRestoreArrayF90(dfdt, dfdtptr, mpi_err)
      call VecRestoreArrayF90(k1, k1ptr, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! solve for k2
      call KSPSolve(ksp, rhs, k2, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! update y solution vector and change t location
      call VecGetArrayF90(y, yptr, mpi_err)
      call VecGetArrayF90(yinit, yinitptr, mpi_err)
      call VecGetArrayF90(k1, k1ptr, mpi_err)
      call VecGetArrayF90(k2, k2ptr, mpi_err)
      yptr = yinitptr + A31*k1ptr + A32*k2ptr
      t = tinit + D3*h
      call VecRestoreArrayF90(y, yptr, mpi_err)
      call VecRestoreArrayF90(yinit, yinitptr, mpi_err)
      call VecRestoreArrayF90(k1, k1ptr, mpi_err)
      call VecRestoreArrayF90(k2, k2ptr, mpi_err)

      ! recompute derivative vector with new (x,y)
      call derivs(t, y, dydt, n)

      ! compute right hand side for k3
      call VecGetArrayF90(rhs, rhsptr, mpi_err)
      call VecGetArrayF90(dydt, dydtptr, mpi_err)
      call VecGetArrayF90(dfdt, dfdtptr, mpi_err)
      call VecGetArrayF90(k1, k1ptr, mpi_err)
      call VecGetArrayF90(k2, k2ptr, mpi_err)
      rhsptr = dydtptr + h*B3*dfdtptr + (C31*k1ptr + C32*k2ptr)/h
      call VecRestoreArrayF90(rhs, rhsptr, mpi_err)
      call VecRestoreArrayF90(dydt, dydtptr, mpi_err)
      call VecRestoreArrayF90(dfdt, dfdtptr, mpi_err)
      call VecRestoreArrayF90(k1, k1ptr, mpi_err)
      call VecRestoreArrayF90(k2, k2ptr, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! solve for k3
      call KSPSolve(ksp, rhs, k3, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! compute right hand side for k4
      call VecGetArrayF90(rhs, rhsptr, mpi_err)
      call VecGetArrayF90(dydt, dydtptr, mpi_err)
      call VecGetArrayF90(dfdt, dfdtptr, mpi_err)
      call VecGetArrayF90(k1, k1ptr, mpi_err)
      call VecGetArrayF90(k2, k2ptr, mpi_err)
      call VecGetArrayF90(k3, k3ptr, mpi_err)
      rhsptr = dydtptr + h*B4*dfdtptr + (C41*k1ptr + C42*k2ptr + C43*k3ptr)/h
      call VecRestoreArrayF90(rhs, rhsptr, mpi_err)
      call VecRestoreArrayF90(dydt, dydtptr, mpi_err)
      call VecRestoreArrayF90(dfdt, dfdtptr, mpi_err)
      call VecRestoreArrayF90(k1, k1ptr, mpi_err)
      call VecRestoreArrayF90(k2, k2ptr, mpi_err)
      call VecRestoreArrayF90(k3, k3ptr, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! solve for k4
      call KSPSolve(ksp, rhs, k4, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! calculate 4th order estimate of solution vector and truncation error
      call VecGetArrayF90(y, yptr, mpi_err)
      call VecGetArrayF90(err, errptr, mpi_err)
      call VecGetArrayF90(yinit, yinitptr, mpi_err)
      call VecGetArrayF90(k1, k1ptr, mpi_err)
      call VecGetArrayF90(k2, k2ptr, mpi_err)
      call VecGetArrayF90(k3, k3ptr, mpi_err)
      call VecGetArrayF90(k4, k4ptr, mpi_err)
      yptr = yinitptr + M1*k1ptr + M2*k2ptr + M3*k3ptr + M4*k4ptr
      errptr = E1*k1ptr + E2*k2ptr + E3*k3ptr + E4*k4ptr
      call VecRestoreArrayF90(y, yptr, mpi_err)
      call VecRestoreArrayF90(err, errptr, mpi_err)
      call VecRestoreArrayF90(yinit, yinitptr, mpi_err)
      call VecRestoreArrayF90(k1, k1ptr, mpi_err)
      call VecRestoreArrayF90(k2, k2ptr, mpi_err)
      call VecRestoreArrayF90(k3, k3ptr, mpi_err)
      call VecRestoreArrayF90(k4, k4ptr, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! set final time
      t = tinit + h

      ! compute maximum error ratio
      call VecGetArrayF90(err, errptr, mpi_err)
      call VecGetArrayF90(yscal, yscalptr, mpi_err)
      where (abs(yscalptr) > ZERO)
        errptr = abs(errptr/yscalptr)
      elsewhere
        errptr = ZERO
      end where
      errmax = maxval(errptr)/eps
      call VecRestoreArrayF90(err, errptr, mpi_err)
      call VecRestoreArrayF90(yscal, yscalptr, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! get rid of A matrix
!     call MatDestroy(A, mpi_err)

      ! check for variable time stepping off
      if (.not. var_ts) then
        hnext = h
        hdid = h
        return
      end if

      ! check for successful time step and let it grow
      if (errmax < ONE) then
        hdid = h
        if (errmax > ERRCON) then
          hnext = SAFETY*h*errmax**PGROW
        else
          hnext = GROW*h
        end if
!       call MatDestroy(A, mpi_err)
        return
     else
        hnext = SAFETY*h*errmax**PSHRNK
        h = sign(max(abs(hnext), SHRNK*abs(h)),h)
     end if
   end do

   ! exceed number of tries
   message = 'exceeded number of time step trials'
   call fatal_error()

  end subroutine solve_ts

!===============================================================================
! DESTROY_OBJECTS
!===============================================================================

  subroutine destroy_objects()

    ! destroy all
    call VecDestroy(yinit, mpi_err)
    call VecDestroy(dydtinit, mpi_err)
    call VecDestroy(k1, mpi_err)
    call VecDestroy(k2, mpi_err)
    call VecDestroy(k3, mpi_err)
    call VecDestroy(k4, mpi_err)
    call VecDestroy(rhs, mpi_err)
    call VecDestroy(err, mpi_err)
    call KSPDestroy(ksp, mpi_err) 

  end subroutine destroy_objects

end module runge_kutta

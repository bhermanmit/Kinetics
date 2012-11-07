module point_kinetics 

!-references

  use global,  only: mpi_err

!-module options

  implicit none
  private
  public :: run_pkinetics

!-module external references

# include "finclude/petsc.h90"

contains

!===============================================================================
! RUN_PKINETICS
!===============================================================================

  subroutine run_pkinetics()

!---external references

    use constants,   only: NUM_PRECS
    use global,      only: pke, time
    use output,      only: header
    use runge_kutta, only: execute_rk4

!---local variables

    integer :: i  ! loop counter
    real(8) :: dt ! local time step

    Mat :: dfdy
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---begin execution

    ! print header for run
    call header("POINT KINETICS SIMULATION", level=1)
print *, pke % rho
print *, pke % time
print *, time
stop
    ! initialize data objects
    call init_data(dfdy, dfdt, dydt, y)

    ! set initial values in y
    call set_init(y)

    ! perform 4th order kaps rentrop
    call execute_rk4(y, dfdy, dfdt, dydt, pk_derivs, pk_jacobn, NUM_PRECS+1)

  end subroutine run_pkinetics

!===============================================================================
! INIT_DATA
!===============================================================================

  subroutine init_data(dfdy,dfdt,dydt,y)

!---references

    use constants,  only: NUM_PRECS

!---arguments

    Mat :: dfdy
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---local variables

    integer :: d_nnz(NUM_PRECS+1)
    integer :: o_nnz(NUM_PRECS+1)

!---begin execution

    ! create dfdy matrix
    d_nnz = 2
    d_nnz(1) = NUM_PRECS + 1
    o_nnz = 0
    call MatCreateAIJ(PETSC_COMM_WORLD, NUM_PRECS+1, NUM_PRECS+1,&
                      PETSC_DETERMINE, PETSC_DETERMINE, PETSC_NULL, d_nnz,&
                      PETSC_NULL, o_nnz, dfdy, mpi_err)
    CHKERRQ(mpi_err)

    ! create dfdt vector
    call VecCreateMPI(PETSC_COMM_WORLD, NUM_PRECS+1, PETSC_DETERMINE, dfdt,&
                      mpi_err)
    CHKERRQ(mpi_err)

    ! create dydt vector
    call VecCreateMPI(PETSC_COMM_WORLD, NUM_PRECS+1, PETSC_DETERMINE, dydt,&
                      mpi_err)
    CHKERRQ(mpi_err)

    ! create y solution vector
    call VecCreateMPI(PETSC_COMM_WORLD, NUM_PRECS+1, PETSC_DETERMINE, y,&
                      mpi_err)
    CHKERRQ(mpi_err)

  end subroutine init_data

!===============================================================================
! SETUP_COEFMAT
!===============================================================================

  subroutine pk_jacobn(t,y,dfdy,dfdt,n)

!---external references

    use constants,  only: beta, NUM_PRECS, lambda, pnl
    use global,     only: pke

!---arguments

    integer :: n
    real(8) :: t
    Mat     :: dfdy
    Vec     :: dfdt
    Vec     :: y

!---local variables

    integer :: i   ! loop counter
    real(8) :: rho ! interpolated reactivity
    real(8) :: val ! temp value for matrix setting

    real(8), pointer :: yptr(:)
    real(8), pointer :: dfdtptr(:)

!---begin execution

    ! get reactivity
    rho = get_reactivity(t,dfdt) 

    ! finish setting dfdt
    call VecGetArrayF90(y, yptr, mpi_err)
    call VecGetArrayF90(dfdt, dfdtptr, mpi_err)
    dfdtptr(1) = dfdtptr(1)*sum(beta)/pnl*yptr(1)
    call VecRestoreArrayF90(y, yptr, mpi_err)
    call VecRestoreArrayF90(dfdt, dfdtptr, mpi_err)
    CHKERRQ(mpi_err)

    ! set up jacobian 
    val = (rho*sum(beta) - sum(beta))/pnl
    call MatSetValue(dfdy, 0, 0, val, INSERT_VALUES, mpi_err)
    CHKERRQ(mpi_err)

    ! begin loop around rest of matrix
    do i = 2, NUM_PRECS + 1

      ! set row 1
      val = lambda(i - 1)
      call MatSetValue(dfdy, 0, i-1, val, INSERT_VALUES, mpi_err)
      CHKERRQ(mpi_err)

      ! set diagonal
      val = -lambda(i - 1)
      call MatSetValue(dfdy, i-1, i-1, val, INSERT_VALUES, mpi_err)
      CHKERRQ(mpi_err)

      ! set column 1
      val = beta(i-1) / pnl
      call MatSetValue(dfdy, i-1, 0, val, INSERT_VALUES, mpi_err)
      CHKERRQ(mpi_err)

    end do 

    ! finalize assembly
    call MatAssemblyBegin(dfdy, MAT_FINAL_ASSEMBLY, mpi_err)
    call MatAssemblyEnd(dfdy, MAT_FINAL_ASSEMBLY, mpi_err)
    CHKERRQ(mpi_err)

  end subroutine pk_jacobn 

!===============================================================================
! PK_DERIVS
!===============================================================================

  subroutine pk_derivs(t,y,dydt,n)

!---external references

    use constants,  only: beta, NUM_PRECS, lambda, pnl
    use global,     only: pke

!---arguments

    integer :: n
    real(8) :: t
    Mat     :: dfdy
    Vec     :: dydt
    Vec     :: y

!---local variables

    integer :: i   ! loop counter
    real(8) :: rho ! interpolated reactivity
    real(8) :: val ! temp value for matrix setting
    real(8), pointer :: yptr(:)

!---begin execution

    ! get reactivity
    rho = get_reactivity(t)

    ! get pointer to solution
    call VecGetArrayF90(y, yptr, mpi_err)
    CHKERRQ(mpi_err)

    ! set first row
    val = (rho*sum(beta) - sum(beta))/pnl*yptr(1) + sum(lambda*yptr(2:NUM_PRECS+1))
    call VecSetValue(dydt, 0, val, mpi_err)
    CHKERRQ(mpi_err)

    ! loop around other common rows
    do i = 2, NUM_PRECS+1

      val = beta(i-1)/pnl*yptr(1) - lambda(i-1)*yptr(i)
      call VecSetValue(dydt, i-1, val, mpi_err)
      CHKERRQ(mpi_err)

    end do

    ! put the pointer back
    call VecRestoreArrayF90(y, yptr, mpi_err)
    CHKERRQ(mpi_err)

  end subroutine pk_derivs

!===============================================================================
! SET_INIT
!===============================================================================

  subroutine set_init(y)

!---external references

    use constants,  only: ONE, beta, lambda, pnl, NUM_PRECS
    use global,     only: pke

!---arguments

    Vec :: y

!---local variables

    integer :: i ! loop counter
    real(8), pointer :: yptr(:)

!---begin execution  

    ! get pointer
    call VecGetArrayF90(y, yptr, mpi_err)
    CHKERRQ(mpi_err)

    ! set power at 1.0
    yptr(1) = ONE

    ! loop through precursors
    do i = 2, NUM_PRECS + 1

      ! set initial value
      yptr(i) = beta(i-1)/(pnl*lambda(i-1)) * yptr(1)

    end do

    ! put pointer back
    call VecRestoreArrayF90(y, yptr, mpi_err)
    CHKERRQ(mpi_err)

  end subroutine set_init

!===============================================================================
! SET_REACTIVITY
!===============================================================================

  function get_reactivity(t,dfdt) result(rho)

!---external references

    use constants,  only: beta, pnl, ZERO
    use global,     only: pke

!---arguments

    Vec, optional :: dfdt
    real(8) :: t
    real(8) :: rho

!---local variables

    integer :: idx ! interpolation index
    real(8) :: m   ! slope

!---begin execution

    ! check if index should be moved in input vectors
    if (t < pke % time(idx+1)) pke % idx = pke % idx + 1
    idx = pke % idx

    ! interpolate on reactivity
    m = ((pke % rho(idx+1) - pke % rho(idx))/(pke % time(idx + 1) - pke % time(idx)))
    rho = pke % rho(idx) + m*(t - pke % time(idx))

    ! set dfdt
    if (present(dfdt)) then
      call VecSet(dfdt, ZERO, mpi_err)
      call VecSetValue(dfdt, 0, m, INSERT_VALUES, mpi_err)
      CHKERRQ(mpi_err)
    end if

  end function get_reactivity

end module point_kinetics 

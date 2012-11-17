module nordheim_fuchs 

!-references

  use global,  only: mpi_err

!-module options

  implicit none
  private
  public :: run_nordheimfuchs

!-module external references

# include "finclude/petsc.h90"

contains

!===============================================================================
! RUN_PKINETICS
!===============================================================================

  subroutine run_nordheimfuchs()

!---external references

    use constants,       only: NUM_PRECS
    use global,          only: pke, time, solver_type
    use implicit_euler,  only: execute_ie1
    use output,          only: header
    use runge_kutta,     only: execute_rk4

!---local variables

    integer :: i  ! loop counter
    real(8) :: dt ! local time step

    Mat :: dfdy
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---begin execution

    ! print header for run
    call header("NORDHEIM FUCHS SIMULATION", level=1)

    ! initialize data objects
    call init_data(dfdy, dfdt, dydt, y) 

    ! set initial values in y
    call set_init(y)

    ! perform 4th order kaps rentrop
    select case(trim(solver_type))

      case('rk4')
        call execute_rk4(y, dfdy, dfdt, dydt, pk_derivs, pk_jacobn, 2, pk_post_timestep)

      case DEFAULT

    end select

    call destroy_objects(dfdy, dfdt, dydt, y)

  end subroutine run_nordheimfuchs

!===============================================================================
! INIT_DATA
!===============================================================================

  subroutine init_data(dfdy,dfdt,dydt,y)

!---references

    use constants,  only: ZERO

!---arguments

    Mat :: dfdy
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---local variables

    integer :: d_nnz(2)
    integer :: o_nnz(2)

!---begin execution

    ! create dfdy matrix
    d_nnz = 2
    o_nnz = 0
    call MatCreateAIJ(PETSC_COMM_WORLD, 2, 2,&
                      PETSC_DETERMINE, PETSC_DETERMINE, PETSC_NULL, d_nnz,&
                      PETSC_NULL, o_nnz, dfdy, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create dfdt vector
    call VecCreateMPI(PETSC_COMM_WORLD, 2, PETSC_DETERMINE, dfdt,&
                      mpi_err)
    call VecSet(dfdt, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create dydt vector
    call VecCreateMPI(PETSC_COMM_WORLD, 2, PETSC_DETERMINE, dydt,&
                      mpi_err)
    call VecSet(dydt, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create y solution vector
    call VecCreateMPI(PETSC_COMM_WORLD, 2, PETSC_DETERMINE, y,&
                      mpi_err)
    call VecSet(y, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine init_data

!===============================================================================
! SETUP_COEFMAT
!===============================================================================

  subroutine pk_jacobn(t,y,dfdy,dfdt,n)

!---external references

    use constants,  only: beta, pnl, ONE, ZERO
    use global,     only: pke, fuel_T, fuel_a, fuel_c, fuel_m

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
    pke % rhot = rho

    ! finish setting dfdt
    call VecGetArrayF90(y, yptr, mpi_err)
    call VecGetArrayF90(dfdt, dfdtptr, mpi_err)
    dfdtptr(1) = dfdtptr(1)*sum(beta)/pnl*yptr(1)
    call VecRestoreArrayF90(dfdt, dfdtptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! set up jacobian 
    ! power equation to power
    val = (rho*sum(beta) - fuel_a*(yptr(2) - fuel_T) - sum(beta))/pnl
    call MatSetValue(dfdy, 0, 0, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! power equation to temperature of fuel
    val = -fuel_a/pnl*yptr(1)
    call MatSetValue(dfdy, 0, 1, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! temperature of fuel equation to power
    val = ONE/(fuel_m*fuel_c)
    call MatSetValue(dfdy, 1, 0, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! temperature of fuel equation to temperature of fuel
    val = ZERO
    call MatSetValue(dfdy, 1, 1, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! finalize assembly
    call MatAssemblyBegin(dfdy, MAT_FINAL_ASSEMBLY, mpi_err)
    call MatAssemblyEnd(dfdy, MAT_FINAL_ASSEMBLY, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! restore yptr to y
    call VecRestoreArrayF90(y, yptr, mpi_err)

  end subroutine pk_jacobn 

!===============================================================================
! PK_DERIVS
!===============================================================================

  subroutine pk_derivs(t,y,dydt,n)

!---external references

    use constants,  only: beta, pnl, ONE
    use global,     only: pke, fuel_a, fuel_c, fuel_m, fuel_T

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
    real(8), pointer :: dydtptr(:)

!---begin execution

    ! get reactivity
    rho = get_reactivity(t)
!   rho = pke % rhot

    ! get pointer to solution
    call VecGetArrayF90(y, yptr, mpi_err)
    call VecGetArrayF90(dydt, dydtptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! set power row
    val = (rho*sum(beta) - fuel_a*(yptr(2) - fuel_T) - sum(beta))/pnl*yptr(1)
    dydtptr(1) = val    

    ! set temperature of fuel row
    val = ONE/(fuel_m*fuel_c)*yptr(1)
    dydtptr(2) = val

    ! put the pointer back
    call VecRestoreArrayF90(y, yptr, mpi_err)
    call VecRestoreArrayF90(dydt, dydtptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine pk_derivs

!===============================================================================
! SET_INIT
!===============================================================================

  subroutine set_init(y)

!---external references

    use constants,  only: ONE, beta, pnl
    use global,     only: pke, power, fuel_T

!---arguments

    Vec :: y

!---local variables

    integer :: i ! loop counter
    real(8), pointer :: yptr(:)

!---begin execution  

    ! get pointer
    call VecGetArrayF90(y, yptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! set initial power
    yptr(1) = power 

    ! set initial fuel temp
    yptr(2) = fuel_T

    ! put pointer back
    call VecRestoreArrayF90(y, yptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

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

    integer :: i   ! iteration counter
    integer :: idx ! interpolation index
    real(8) :: m   ! slope

!---begin execution

    ! search for index
    idx = 1
    do i = 1, size(pke % time)
      if (t < pke % time(i)) then
        idx = i - 1
        exit
      end if
    end do

    ! interpolate on reactivity
    m = ((pke % rho(idx+1) - pke % rho(idx))/(pke % time(idx + 1) - pke % time(idx)))
    rho = pke % rho(idx) + m*(t - pke % time(idx))

    ! set dfdt
    if (present(dfdt)) then
      call VecSetValue(dfdt, 0, m, INSERT_VALUES, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif
    end if

  end function get_reactivity

!===============================================================================
! PK_POST_TIMESTEP
!===============================================================================

  subroutine pk_post_timestep(t, y, h)

!---refrences

    use constants, only: beta
    use global,  only: fuel_T, fuel_a

!---arguments

    real(8) :: t
    real(8) :: h
    Vec :: y

!---local variables

    real(8) :: react
    real(8), pointer :: yptr(:)

!---begin execution

    ! get ptr
    call VecGetArrayF90(y, yptr, mpi_err)

    ! compute reactivity
    react = get_reactivity(t)
    react = react - fuel_a*(yptr(2) - fuel_T) 
    react = sum(beta)*react

    ! write output
    print *, 'TIME:', t, 'REACT:', react, 'POWER:', yptr(1), 'FUEL:', yptr(2), 'STEP:', h

    ! restore ptrs
    call VecRestoreArrayF90(y, yptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine pk_post_timestep

!===============================================================================
! DESTROY_OBJECTS
!===============================================================================

  subroutine destroy_objects(dfdy, dfdt, dydt, y)

!---arguments

    Mat :: dfdy
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---begin execution

    ! destroy all
    call MatDestroy(dfdy, mpi_err)
    call VecDestroy(dfdt, mpi_err)
    call VecDestroy(dydt, mpi_err)
    call VecDestroy(y, mpi_err)

  end subroutine destroy_objects

end module nordheim_fuchs 

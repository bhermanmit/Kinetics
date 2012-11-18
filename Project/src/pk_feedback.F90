module pk_feedback

!-references

  use constants,  only: NUM_PRECS
  use global,     only: mpi_err

!-module options

  implicit none
  private
  public :: run_pkfeedback

!-module external references

# include "finclude/petsc.h90"

!-module variables

  integer, parameter :: FTEMP = NUM_PRECS + 2 ! index for fuel temp
  integer, parameter :: CTEMP = NUM_PRECS + 3 ! index for cool temp

contains

!===============================================================================
! RUN_PKFEEDBACK
!===============================================================================

  subroutine run_pkfeedback()

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
    Mat :: A
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---begin execution

    ! print header for run
    call header("POINT KINETICS SIMULATION w/ FEEDBACK", level=1)

    ! initialize data objects
    call init_data(dfdy, dfdt, dydt, y, A)

    ! set initial values in y
    call set_init(y)

    ! perform 4th order kaps rentrop
    select case(trim(solver_type))

      case('rk4')
        call execute_rk4(y, dfdy, dfdt, dydt, pk_derivs, pk_jacobn, NUM_PRECS+3, pk_post_timestep)

      case('ie1')
        call execute_ie1(y, A, pk_coefmat, NUM_PRECS+3)

      case DEFAULT

    end select

    call destroy_objects(dfdy, A, dfdt, dydt, y)

  end subroutine run_pkfeedback

!===============================================================================
! INIT_DATA
!===============================================================================

  subroutine init_data(dfdy,dfdt,dydt,y,A)

!---references

    use constants,  only: NUM_PRECS, ZERO

!---arguments

    Mat :: dfdy
    Mat :: A
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---local variables

    integer :: d_nnz(NUM_PRECS+3)
    integer :: o_nnz(NUM_PRECS+3)

!---begin execution

    ! create dfdy matrix
    d_nnz = 2
    d_nnz(1) = NUM_PRECS + 3
    d_nnz(FTEMP:CTEMP) = 3
    o_nnz = 0
    call MatCreateAIJ(PETSC_COMM_WORLD, NUM_PRECS+3, NUM_PRECS+3,&
                      PETSC_DETERMINE, PETSC_DETERMINE, PETSC_NULL, d_nnz,&
                      PETSC_NULL, o_nnz, dfdy, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create A matrix
    d_nnz = 2
    d_nnz(1) = NUM_PRECS + 3
    d_nnz(FTEMP:CTEMP) = 3
    o_nnz = 0
    call MatCreateAIJ(PETSC_COMM_WORLD, NUM_PRECS+3, NUM_PRECS+3,&
                      PETSC_DETERMINE, PETSC_DETERMINE, PETSC_NULL, d_nnz,&
                      PETSC_NULL, o_nnz, A, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create dfdt vector
    call VecCreateMPI(PETSC_COMM_WORLD, NUM_PRECS+3, PETSC_DETERMINE, dfdt,&
                      mpi_err)
    call VecSet(dfdt, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create dydt vector
    call VecCreateMPI(PETSC_COMM_WORLD, NUM_PRECS+3, PETSC_DETERMINE, dydt,&
                      mpi_err)
    call VecSet(dydt, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create y solution vector
    call VecCreateMPI(PETSC_COMM_WORLD, NUM_PRECS+3, PETSC_DETERMINE, y,&
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

    use constants,  only: beta, NUM_PRECS, lambda, pnl, ONE
    use global,     only: pke, fuel_a, fuel_c, fuel_m, cool_a, cool_c, cool_m, &
                          P_frac, hA, m_dot

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

    ! extract ptrs
    call VecGetArrayF90(y, yptr, mpi_err)
    call VecGetArrayF90(dfdt, dfdtptr, mpi_err)

    ! finish setting dfdt
    dfdtptr(1) = dfdtptr(1)*sum(beta)/pnl*yptr(1)

    ! restore pointer
    call VecRestoreArrayF90(dfdt, dfdtptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! set up jacobian 
    val = (rho*sum(beta) - sum(beta))/pnl
    call MatSetValue(dfdy, 0, 0, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! begin loop around rest of matrix
    do i = 2, NUM_PRECS + 1

      ! set row 1
      val = lambda(i - 1)
      call MatSetValue(dfdy, 0, i-1, val, INSERT_VALUES, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! set diagonal
      val = -lambda(i - 1)
      call MatSetValue(dfdy, i-1, i-1, val, INSERT_VALUES, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! set column 1
      val = beta(i-1) / pnl
      call MatSetValue(dfdy, i-1, 0, val, INSERT_VALUES, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

    end do 

    ! first row fuel temperature
    val = -fuel_a/pnl*yptr(1)
    call MatSetValue(dfdy, 0, FTEMP-1, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! first row coolant temperature
    val = -cool_a/pnl*yptr(1)
    call MatSetValue(dfdy, 0, CTEMP-1, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! fuel temperature row power
    val = (ONE - P_frac)/(fuel_m*fuel_c)
    call MatSetValue(dfdy, FTEMP-1, 0, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! fuel temperature row fuel temperature
    val = -hA/(fuel_m*fuel_c)
    call MatSetValue(dfdy, FTEMP-1, FTEMP-1, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! fuel temperature row coolant temperature
    val = hA/(fuel_m*fuel_c)
    call MatSetValue(dfdy, FTEMP-1, CTEMP-1, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! coolant temperature row power
    val = P_frac/(cool_m*cool_c)
    call MatSetValue(dfdy, CTEMP-1, 0, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! fuel temperature row fuel temperature
    val = hA/(cool_m*cool_c)
    call MatSetValue(dfdy, CTEMP-1, FTEMP-1, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! fuel temperature row coolant temperature
    val = -(hA + 2*m_dot*cool_c)/(cool_m*cool_c)
    call MatSetValue(dfdy, CTEMP-1, CTEMP-1, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! restore pointer
    call VecRestoreArrayF90(y, yptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! finalize assembly
    call MatAssemblyBegin(dfdy, MAT_FINAL_ASSEMBLY, mpi_err)
    call MatAssemblyEnd(dfdy, MAT_FINAL_ASSEMBLY, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine pk_jacobn 

!===============================================================================
! PK_DERIVS
!===============================================================================

  subroutine pk_derivs(t,y,dydt,n)

!---external references

    use constants,  only: beta, NUM_PRECS, lambda, pnl, TWO, ONE
    use global,     only: pke, fuel_a, fuel_c, fuel_m, cool_a, cool_c, cool_m, &
                          P_frac, hA, m_dot, Tin, fuel_T, cool_T

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

    ! get pointer to solution
    call VecGetArrayF90(y, yptr, mpi_err)
    call VecGetArrayF90(dydt, dydtptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! set first row
    val = (rho*sum(beta) - fuel_a*(yptr(FTEMP) - fuel_T) -&
                           cool_a*(yptr(CTEMP) - cool_T) -&
                           sum(beta))/pnl*yptr(1) +&
                           sum(lambda*yptr(2:NUM_PRECS+1))
    dydtptr(1) = val

    ! loop around other common rows
    do i = 2, NUM_PRECS+1

      val = beta(i-1)/pnl*yptr(1) - lambda(i-1)*yptr(i)
      dydtptr(i) = val

    end do

    ! set fuel temperature row
    val = (ONE - P_frac)/(fuel_m*fuel_c)*yptr(1) - &
          hA/(fuel_m*fuel_c)*(yptr(FTEMP) - yptr(CTEMP))
    dydtptr(FTEMP) = val

    ! set coolant temperature row
    val = P_frac/(cool_m*cool_c)*yptr(1) + &
          hA/(cool_m*cool_c)*(yptr(FTEMP) - &
          yptr(CTEMP)) - TWO*m_dot/cool_m*(yptr(CTEMP) - Tin)
    dydtptr(CTEMP) = val

    ! put the pointer back
    call VecRestoreArrayF90(y, yptr, mpi_err)
    call VecRestoreArrayF90(dydt, dydtptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine pk_derivs

!===============================================================================
! COEFMAT
!===============================================================================

  subroutine pk_coefmat(t,h,y,A,n)

!---external references

    use constants,  only: beta, NUM_PRECS, lambda, pnl, ONE
    use global,     only: pke

!---arguments

    integer :: n
    real(8) :: t
    real(8) :: h
    Mat     :: A 
    Vec     :: y

!---local variables

    integer :: i   ! loop counter
    real(8) :: rho ! interpolated reactivity
    real(8) :: val ! temp value for matrix setting

    real(8), pointer :: yptr(:)
    real(8), pointer :: dfdtptr(:)

!---begin execution

    ! get reactivity
    rho = get_reactivity(t) 

    ! set up jacobian 
    val = (rho*sum(beta) - sum(beta))/pnl
    val = ONE - h*val
    call MatSetValue(A, 0, 0, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! begin loop around rest of matrix
    do i = 2, NUM_PRECS + 1

      ! set row 1
      val = lambda(i - 1)
      val = -h*val
      call MatSetValue(A, 0, i-1, val, INSERT_VALUES, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! set diagonal
      val = -lambda(i - 1)
      val = ONE - h*val
      call MatSetValue(A, i-1, i-1, val, INSERT_VALUES, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! set column 1
      val = beta(i-1) / pnl
      val = -h*val
      call MatSetValue(A, i-1, 0, val, INSERT_VALUES, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

    end do 

    ! finalize assembly
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, mpi_err)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine pk_coefmat 

!===============================================================================
! SET_INIT
!===============================================================================

  subroutine set_init(y)

!---external references

    use constants,  only: ONE, beta, lambda, pnl, NUM_PRECS, TWO
    use global,     only: pke, power, fuel_T, cool_T, m_dot, cool_c, Tin, hA, P_frac

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

    ! set power
    yptr(1) = power

    ! loop through precursors
    do i = 2, NUM_PRECS + 1

      ! set initial value
      yptr(i) = beta(i-1)/(pnl*lambda(i-1)) * yptr(1)

    end do

    ! set coolant temperature
    cool_T = power/(TWO*m_dot*cool_c) + Tin
    yptr(CTEMP) = cool_T

    ! set fuel temperature
    fuel_T = (ONE - P_frac)*power/hA + cool_T
    yptr(FTEMP) = fuel_T

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
      idx = i - 1
      if (t < pke % time(i)) exit 
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

!---references

    use constants,  only: beta
    use global,     only: fuel_a, cool_a, fuel_T, cool_T

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
    react = react*sum(beta) - fuel_a*(yptr(FTEMP) - fuel_T) -&
                              cool_a*(yptr(CTEMP) - cool_T)
    react = react/sum(beta)

    ! write output
    print *, 'TIME:', t, 'REACT:', react, 'POWER:', yptr(1), 'FUEL:', yptr(FTEMP), 'COOL:', yptr(CTEMP), 'STEP:', h

    ! restore ptrs
    call VecRestoreArrayF90(y, yptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine pk_post_timestep

!===============================================================================
! DESTROY_OBJECTS
!===============================================================================

  subroutine destroy_objects(dfdy, A, dfdt, dydt, y)

!---arguments

    Mat :: dfdy
    Mat :: A
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---begin execution

    ! destroy all
    call MatDestroy(A, mpi_err)
    call MatDestroy(dfdy, mpi_err)
    call VecDestroy(dfdt, mpi_err)
    call VecDestroy(dydt, mpi_err)
    call VecDestroy(y, mpi_err)

  end subroutine destroy_objects

end module pk_feedback

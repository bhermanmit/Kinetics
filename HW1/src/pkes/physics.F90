module physics

!-module options

  implicit none
  private
  public :: run_kinetics

contains

!===============================================================================
! RUN_KINETICS
!===============================================================================

  subroutine run_kinetics()

!---external references

    use constants,  only: NUM_PRECS
    use expokit,    only: dense_pade
    use global,     only: pke
    use gnufor2,    only: plot
    use output,     only: header

!---local variables

    integer :: i ! loop counter

!---begin execution

    ! print header for run
    call header("POINT KINETICS SIMULATION", level=1)

    ! set up coefficient matrix
    call setup_coefmat()

    ! set up initial conditions
    call set_init()

    ! begin loop through time steps
    do i = 1, pke % nt

      ! compute exponential matrix
      call set_reactivity(i)

      ! solve matrix exponential
      call dense_pade(pke % coef, NUM_PRECS + 1, pke % dt, pke % expm)

      ! get new vector 
      pke % N(:,i+1) = matmul(pke % expm, pke % N(:,i))

    end do

    call plot(pke % time, pke % N(1,:))

  end subroutine run_kinetics

!===============================================================================
! SETUP_COEFMAT
!===============================================================================

  subroutine setup_coefmat()

!---external references

    use constants,  only: beta, NUM_PRECS, lambda, pnl
    use global,     only: pke

!---local variables

    integer :: i ! loop counter

!---begin execution

    ! set up coeffient matrix manually for time 0
    pke % coef(1,1) = (pke % rho(1)*sum(beta) - sum(beta))/pnl

    ! begin loop around rest of matrix
    do i = 2, NUM_PRECS + 1

      ! set row 1
      pke % coef(1,i) = lambda(i - 1)

      ! set diagonal
      pke % coef(i,i) = -lambda(i - 1)

      ! set column 1
      pke % coef(i,1) = beta(i-1) / pnl

    end do 

  end subroutine setup_coefmat

!===============================================================================
! SET_INIT
!===============================================================================

  subroutine set_init()

!---external references

    use constants,  only: ONE, beta, lambda, pnl, NUM_PRECS
    use global,     only: pke

!---local variables

    integer :: i ! loop counter

!---begin execution  

    ! set power at 1.0
    pke % N(1,1) = ONE

    ! loop through precursors
    do i = 1, NUM_PRECS

      ! set initial value
      pke % N(i+1,1) = beta(i)/(pnl*lambda(i)) * pke % N(1,1)

    end do

  end subroutine set_init

!===============================================================================
! SET_REACTIVITY
!===============================================================================

  subroutine set_reactivity(i)

!---external references

    use constants,  only: beta, pnl
    use global,     only: pke

!---arguments

    integer :: i   ! current time step

!---local variables

    integer :: idx ! interpolation index

!---begin execution

    ! compute current time
    pke % time(i) = float(i - 1) * pke % dt

    ! check if index should be moved in input vectors
    if (pke % time(i) > pke % t(pke % idx + 1)) pke % idx = pke % idx + 1
    idx = pke % idx

    ! interpolate on reactivity
    pke % react(i) = pke % rho(idx) + ((pke % rho(idx+1) - pke % rho(idx))  /  &
                                       (pke % t(idx + 1) - pke % t(idx)))   *  &
                                       (pke % time(i) - pke % t(idx))

    ! set values in coefficient matrix
    pke % coef(1,1) = (pke % react(i)*sum(beta) - sum(beta))/pnl

  end subroutine

end module physics

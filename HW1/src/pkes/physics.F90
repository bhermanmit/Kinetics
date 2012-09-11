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

    use global,  only: pke

!---local variables

    integer :: i ! loop counter

!---begin execution

    ! set up coefficient matrix
    call setup_coefmat()

    ! set up initial conditions

    ! begin loop through time steps
    do i = 2, pke % nt

      ! compute exponential matrix

      ! get new vector 

    end do

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

!---begin execution  

  end subroutine set_init

end module physics

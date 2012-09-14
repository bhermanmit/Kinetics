module physics

!-module options

  implicit none
  private
  public :: run_invkinetics

contains

!===============================================================================
! RUN_INVKINETICS
!===============================================================================

  subroutine run_invkinetics()

!---external references

    use constants,  only: NUM_PRECS, pnl, beta, lambda
    use global,     only: ipke
    use output,     only: header

!---local variables

    integer :: i        ! loop counter
    real(8) :: dt       ! local time step
    real(8) :: avgpower ! avg. power between timesteps

!---begin execution

    ! print header for run
    call header("POINT KINETICS SIMULATION", level=1)

    ! set initial precursors
    ipke % N(1,1) = ipke % power(1)
    ipke % N(2:NUM_PRECS+1,1) = beta/(pnl*lambda) * ipke % N(1,1)

    ! begin loop through time steps
    do i = 1, sum(ipke % nt)

      ! compute exponential matrix
      call set_power(i,dt)

      ! compute average power
      avgpower = (ipke % N(1,i) + ipke % N(1,i+1)) / 2.0_8

      ! solve for precursors 
      ipke %N(2:NUM_PRECS+1,i+1) = ipke % N(2:NUM_PRECS+1,i)*exp(-lambda*dt) + &
      beta/(lambda*pnl)*(avgpower - avgpower*exp(-lambda*dt))

      ! solve for reactivity
      ipke % react(i+1) = pnl/(ipke % N(1,i+1)) *                              &
                          ((ipke % N(1,i+1) - ipke % N(1,i))/dt) +             &
                          sum(beta) - pnl/ipke % N(1,i+1) *                    &
                          sum(lambda*ipke %N(2:NUM_PRECS+1,i+1))

      ! change to dollars
      ipke % react(i+1) = ipke % react(i+1)/sum(beta)

    end do

  end subroutine run_invkinetics

!===============================================================================
! SET_POWER
!===============================================================================

  subroutine set_power(i,dt)

!---external references

    use constants,  only: beta, pnl
    use global,     only: ipke

!---arguments

    integer :: i   ! current time step
    real(8) :: dt  ! current dt

!---local variables

    integer :: idx ! interpolation index

!---begin execution

    ! check if index should be moved in input vectors
    if (i > sum(ipke % nt(1:ipke % idx))) ipke % idx = ipke % idx + 1
    idx = ipke % idx

    ! compute current time
    dt = ipke % dt(idx)
    ipke % time(i+1) = ipke % time(i) + dt 

    ! interpolate on power 
    ipke % N(1,i+1) = ((ipke % power(idx+1) - ipke % power(idx))  /  &
                      (ipke % t(idx + 1) - ipke % t(idx)))   *  &
                      (ipke % time(i+1) - ipke % t(idx)) + &
                      ipke % power(idx)

  end subroutine set_power

end module physics

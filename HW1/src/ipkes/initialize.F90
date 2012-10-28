module initialize

!-module options

  implicit none
  private
  public :: initialize_run

contains

!===============================================================================
! INITIALIZE_RUN
!===============================================================================

  subroutine initialize_run()

!---external references

    use data_header,  only: allocate_ipke_type
    use global,       only: ipke, total_time
    use input_xml,    only: read_input_xml
    use output,       only: title, write_input, header, write_physics
    use timing,       only: timer_start

!---local variables

    integer :: i   ! loop counter
    integer :: nt  ! local number of steps

!---begin execution

    ! start timer
    call timer_start(total_time)

    ! print title
    call title() 
    call header("INITIALIZATION", level=1)

    ! read in input
    call read_input_xml()

    ! compute number of time steps
    allocate(ipke % nt(ipke % npts - 1))
    do i = 1, ipke % npts - 1

      ! compute integer number of time steps
      nt = nint((ipke % t(i+1) - ipke % t(i)) / ipke % dt(i))

      ! recompute time step
      ipke % dt(i) = (ipke % t(i+1) - ipke % t(i)) / dble(nt)

      ! append local # of timesteps to full counter
      ipke % nt(i) = nt 

    end do

    ! echo input
    call header("Input Summary", level=2)
    call write_input()

    ! initialize data matrices
    call allocate_ipke_type(ipke)

    ! echo physics
    call header("Physics Summary", level=2)
    call write_physics()

  end subroutine initialize_run

end module initialize

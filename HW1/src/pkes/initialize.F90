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

    use data_header,  only: allocate_pke_type
    use global,       only: pke, total_time
    use input_xml,    only: read_input_xml, read_binary
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
!  call read_input_xml()
   call read_binary()

    ! compute number of time steps
    allocate(pke % nt(pke % npts - 1))
    do i = 1, pke % npts - 1

      ! compute integer number of time steps
      nt = nint((pke % t(i+1) - pke % t(i)) / pke % dt(i))

      ! recompute time step
      pke % dt(i) = (pke % t(i+1) - pke % t(i)) / dble(nt)

      ! append local # of timesteps to full counter
      pke % nt(i) = nt 

    end do

    ! echo input
    call header("Input Summary", level=2)
    call write_input()

    ! initialize data matrices
    call allocate_pke_type(pke)

    ! echo physics
    call header("Physics Summary", level=2)
    call write_physics()

  end subroutine initialize_run

end module initialize

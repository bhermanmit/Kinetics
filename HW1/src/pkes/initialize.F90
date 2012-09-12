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

    use data_header,  only: allocate_pke_type
    use global,       only: pke, total_time
    use input_xml,    only: read_input_xml
    use output,       only: title, write_input, header, write_physics
    use timing,       only: timer_start

    ! start timer
    call timer_start(total_time)

    ! print title
    call title() 
    call header("INITIALIZATION", level=1)

    ! read in input
    call read_input_xml()

    ! compute number of time steps
    pke % nt = floor(pke % maxt / pke % dt) + 1

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

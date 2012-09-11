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
    use global,       only: pke

    ! print heading

    ! compute number of time steps
    pke % nt = floor(pke % maxt / pke % dt) + 1

    ! initialize data matrices
    call allocate_pke_type(pke)

  end subroutine initialize_run

end module initialize

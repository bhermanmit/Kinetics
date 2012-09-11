module finalize 

!-module options

  implicit none
  private
  public :: finalize_run

contains

!===============================================================================
! FINALIZE_RUN
!===============================================================================

  subroutine finalize_run()

    use data_header,  only: deallocate_pke_type
    use global,       only: pke

    ! finalize data matrices
    call deallocate_pke_type(pke)

  end subroutine finalize_run

end module finalize 


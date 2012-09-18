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

    use output,  only: write_results

    call write_results()

  end subroutine finalize_run

end module finalize 


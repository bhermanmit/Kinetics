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
    use global,       only: pke, total_time
    use output,       only: plot_results, write_results
    use timing,       only: timer_stop

    ! stop timer
    call timer_stop(total_time)

    ! run plots
    call plot_results()

    ! write summary
    call write_results()

    ! finalize data matrices
    call deallocate_pke_type(pke)

  end subroutine finalize_run

end module finalize 


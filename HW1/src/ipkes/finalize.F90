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

    use data_header,  only: deallocate_ipke_type
    use global,       only: ipke, total_time
    use output,       only: plot_results, write_results, write_output
    use timing,       only: timer_stop

    ! run plots
    call plot_results()

    ! stop timer
    call timer_stop(total_time)

    ! write summary
    call write_results()

    ! write output
    call write_output()

    ! finalize data matrices
    call deallocate_ipke_type(ipke)

  end subroutine finalize_run

end module finalize 


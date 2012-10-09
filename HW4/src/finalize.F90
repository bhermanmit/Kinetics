module finalize 

!-module options

  implicit none
  private
  public :: finalize_run

!-external packages

# include <finclude/petsc.h90>
# include <finclude/slepcsys.h>
# include <finclude/slepceps.h>


contains

!===============================================================================
! FINALIZE_RUN
!===============================================================================

  subroutine finalize_run()

!---external references

    use cmfd_header,      only: deallocate_cmfd_type
    use geometry_header,  only: deallocate_geometry_type
    use global,           only: mpi_err, geometry, material, cmfd, n_materials,&
                                time_total
    use hdf5_interface,   only: hdf5_finalize
    use material_header,  only: deallocate_material_type
    use output,           only: write_results, header, write_hdf5
    use timing,           only: timer_stop

!---local variables

    integer :: i

!---begin execution

    ! stop timer
    call timer_stop(time_total)

    ! print out to user
    call header('CALCULATION FINISHED',level=1)

    ! write out results to user
    call write_results()
!   call write_hdf5()

    ! call finalization routines
    call SlepcFinalize(PETSC_COMM_WORLD,mpi_err) 
    call hdf5_finalize()
    call deallocate_cmfd_type(cmfd)
    call deallocate_geometry_type(geometry)
    do i = 1, n_materials
      call deallocate_material_type(material(i))
    end do

  end subroutine finalize_run

end module finalize 

module initialize

!-module options

  implicit none
  private
  public :: initialize_run

!-external packages

# include <finclude/petsc.h90>
# include <finclude/slepcsys.h>
# include <finclude/slepceps.h>

contains

!===============================================================================
! INITIALIZE_RUN
!===============================================================================

  subroutine initialize_run()

!---begin execution

    ! initailize PETSc/SLEPc
    call petsc_init()

  end subroutine initialize_run

!===============================================================================
! PETSC_INIT
!===============================================================================

  subroutine petsc_init()

!---external references

    use global,  only: rank, n_procs, mpi_err, master

!---begin execution

    ! initialize it
    call SlepcInitialize(PETSC_NULL_CHARACTER,mpi_err)

    ! get the mpi info
    call MPI_COMM_RANK(PETSC_COMM_WORLD,rank,mpi_err)
    call MPI_COMM_RANK(PETSC_COMM_WORLD,n_procs,mpi_err)

    ! set master
    if (rank == 0) master = .true.

  end subroutine petsc_init

end module initialize

module initialize

!-module options

  implicit none
  private
  public :: initialize_run

!-external packages

# include <finclude/petsc.h90>

contains

!===============================================================================
! INITIALIZE_RUN
!===============================================================================

  subroutine initialize_run()

!---external references

    use cmfd_header,      only: allocate_cmfd_type
    use input_xml,        only: read_input_xml
    use global,           only: geometry, cmfd, time_total, time_init
    use hdf5_interface,   only: hdf5_initialize
    use output,           only: header
    use timing,           only: timer_start, timer_stop

!---begin execution

    ! start timer
    call timer_start(time_total)
    call timer_start(time_init)

    ! set initialization title
    call header('INITIALIZATION',level=1)

    ! initailize PETSc/SLEPc
    call petsc_init()

    ! initialize hdf5
    call hdf5_initialize()

    ! read in input
    call read_input_xml()

    ! initialize materials
    call materials_init()

    ! initialize geometry
    call geometry_init()

    ! initialize cmfd data
    call allocate_cmfd_type(cmfd,geometry)

    ! stop initialization timer
    call timer_stop(time_init)

  end subroutine initialize_run

!===============================================================================
! PETSC_INIT
!===============================================================================

  subroutine petsc_init()

!---external references

    use global,  only: rank, n_procs, mpi_err, master

!---begin execution

    ! initialize it
    call PetscInitialize(PETSC_NULL_CHARACTER,mpi_err)

    ! get the mpi info
    call MPI_COMM_RANK(PETSC_COMM_WORLD,rank,mpi_err)
    call MPI_COMM_RANK(PETSC_COMM_WORLD,n_procs,mpi_err)

    ! set master
    if (rank == 0) master = .true.

  end subroutine petsc_init

!===============================================================================
! MATERIALS_INIT
!===============================================================================

  subroutine materials_init()

!---external references

    use global,           only: material, n_materials, geometry
    use material_header,  only: arrange_xs

!---local variables

    integer :: i ! loop counter

!---begin execution

    ! arrange material cross sections
    do i = 1, n_materials

      call arrange_xs(material(i),geometry % ncg)

    end do

  end subroutine materials_init

!===============================================================================
! GEOMETRY_INIT
!===============================================================================

  subroutine geometry_init()

!---external references

    use geometry_header,  only: generate_fine_map, compute_widths
    use global,           only: geometry, material

!---begin execution

    ! compute cell widths
    call compute_widths(geometry) 

    ! set up geometry fine map
    call generate_fine_map(geometry, material)

    ! compute geometry product
    geometry % nf = geometry % nfx * geometry % nfy * geometry % nfz * geometry % nfg

  end subroutine geometry_init

end module initialize

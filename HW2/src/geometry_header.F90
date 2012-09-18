module geometry_header

!-module options

  implicit none
  private
  public :: allocate_geometry_type, deallocate_geometry_type, generate_fine_map

!-module variables

  type, public :: geometry_type

    ! coarse indices
    integer :: ncx 
    integer :: ncy
    integer :: ncz
    integer :: ncg 

    ! fine indices
    integer :: nfx
    integer :: nfy
    integer :: nfz
    integer :: nfg

    ! material map
    integer, allocatable :: mat_map(:,:,:)

    ! fine to coarse map based on materials
    integer, allocatable :: fine_map(:,:,:)

    ! size of each coarse mesh
    real(8), allocatable :: xgrid(:)
    real(8), allocatable :: ygrid(:)
    real(8), allocatable :: zgrid(:)

    ! fine mesh per coarse mesh
    integer, allocatable :: nnx(:)
    integer, allocatable :: nny(:)
    integer, allocatable :: nnz(:)

    ! fine mesh spacing
    real(8), allocatable :: dx(:)
    real(8), allocatable :: dy(:)
    real(8), allocatable :: dz(:)

    ! boundary conditions
    real(8) :: bc(6)

  end type geometry_type 

contains

!===============================================================================
! ALLOCATE_GEOMETRY_TYPE
!===============================================================================

  subroutine allocate_geometry_type(this)

!---arguments

    type(geometry_type) :: this

!---begin execution

    ! allocate maps
    allocate(this % mat_map (this % ncx, this % ncy, this % ncz))
    allocate(this % fine_map(this % nfx, this % nfy, this % nfz))

    ! allocate grids
    allocate(this % xgrid(this % ncx))
    allocate(this % ygrid(this % ncy))
    allocate(this % zgrid(this % ncz))

    ! allocate mesh points
    allocate(this % nnx(this % ncx))
    allocate(this % nny(this % ncy))
    allocate(this % nnz(this % ncz))

    ! allocate spacings
    allocate(this % dx(this % ncx))
    allocate(this % dy(this % ncy))
    allocate(this % dz(this % ncz))

  end subroutine allocate_geometry_type

!===============================================================================
! GENERATE_FINE_MAP
!===============================================================================

  subroutine generate_fine_map(this)

!---arguments

    type(geometry_type) :: this

!---begin execution

  end subroutine generate_fine_map

!===============================================================================
! DEALLOCATE_GEOMETRY_TYPE
!===============================================================================

  subroutine deallocate_geometry_type(this)

!---arguments

    type(geometry_type) :: this

!---begin execution

    ! deallocate all
    deallocate(this % mat_map)
    deallocate(this % fine_map)
    deallocate(this % xgrid)
    deallocate(this % ygrid)
    deallocate(this % zgrid)
    deallocate(this % dx)
    deallocate(this % dy)
    deallocate(this % dz)

  end subroutine deallocate_geometry_type

end module geometry_header

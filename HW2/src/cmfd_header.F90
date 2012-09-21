module cmfd_header

!-module options

  implicit none
  private
  public :: allocate_cmfd_type, deallocate_cmfd_type, calc_power

!-module variables

  type, public :: cmfd_type

    ! eigenvector/eigenvalue from cmfd run
    real(8), allocatable :: phi(:)
    real(8) :: keff = 0.0_8

    ! nodal powers
    real(8), allocatable :: power_o(:)
    real(8), allocatable :: power_n(:)

    ! final iteration and norm
    integer :: iter
    real(8) :: norm
    real(8) :: norm_o
    real(8) :: dr

  end type cmfd_type

contains

!===============================================================================
! ALLOCATE_CMFD_TYPE
!===============================================================================

  subroutine allocate_cmfd_type(this,geometry)

!---external references

    use geometry_header,  only: geometry_type

!---arguments

    integer :: n
    type(cmfd_type) :: this
    type(geometry_type) :: geometry

!---begin execution

    ! allocate phi
    allocate(this % phi(geometry % nfx * geometry % nfy * geometry % nfz *     &
                        geometry % nfg))

    ! allocate nodal power maps
    allocate(this % power_o(geometry % n_regs))
    allocate(this % power_n(geometry % n_regs))

  end subroutine allocate_cmfd_type

!===============================================================================
! CALC_POWER
!===============================================================================

  subroutine calc_power(this,fsrc,n,geometry)

!---external references

    use constants,        only: ZERO
    use geometry_header,  only: geometry_type

!---arguments

    integer :: n
    type(cmfd_type) :: this
    type(geometry_type) :: geometry
    real(8) :: fsrc(n)

!---local variables

    integer :: irow
    integer :: i
    integer :: j
    integer :: k
    integer :: g
    integer :: reg

!---begin execution

    ! zero out new power
    this % power_n = ZERO

    ! begin loop around rows
    do irow = 1, n

      ! get coordinates
      call matrix_to_indices(irow-1,g,i,j,k,geometry % nfx, geometry % nfy,    &
                                          geometry % nfz, geometry % nfg) 

      ! get region number
      reg = geometry % fine_map(i,j,k) % reg

      ! sum nodal power
      this % power_n(reg) = this % power_n(reg) + fsrc(irow)

    end do

  end subroutine calc_power

!===============================================================================
! MATRIX_TO_INDICES 
!===============================================================================

  subroutine matrix_to_indices(irow,g,i,j,k,nx,ny,nz,ng)

    integer :: i                    ! iteration counter for x
    integer :: j                    ! iteration counter for y
    integer :: k                    ! iteration counter for z
    integer :: g                    ! iteration counter for groups
    integer :: irow                 ! iteration counter over row (0 reference)
    integer :: nx                   ! max x
    integer :: ny                   ! max y
    integer :: nz                   ! max z
    integer :: ng                   ! max g

    ! compute indices
    g = mod(irow,ng) + 1
    i = mod(irow,ng*nx)/ng + 1
    j = mod(irow,ng*nx*ny)/(ng*nx)+ 1
    k = mod(irow,ng*nx*ny*nz)/(ng*nx*ny) + 1

  end subroutine matrix_to_indices

!===============================================================================
! DEALLOCATE_CMFD_TYPE
!===============================================================================

  subroutine deallocate_cmfd_type(this)

!---arguments

    type(cmfd_type) :: this

!---begin execution

    deallocate(this % phi)
    deallocate(this % power_o)
    deallocate(this % power_n)

  end subroutine deallocate_cmfd_type

end module cmfd_header

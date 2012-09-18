module material_header

!-module options

  implicit none
  private
  public :: allocate_material_type, deallocate_material_type

!-module variables

  type, public ::  material_type

    ! cross sections
    real(8), allocatable :: totalxs(:)
    real(8), allocatable :: absorxs(:)
    real(8), allocatable :: removxs(:)
    real(8), allocatable :: fissvec(:)
    real(8), allocatable :: scattxs(:,:)
    real(8), allocatable :: nfissxs(:,:)

    ! diffusion coefficient
    real(8), allocatable :: diffcof(:)

    ! fission spectrum 
    real(8), allocatable :: chi(:)

    ! logicals for preprocessing
    logical :: abs_based = .false.
    logical :: rem_based = .false.
    logical :: chi_based = .false.

  end type material_type 

contains

!===============================================================================
! ALLOCATE_MATERAL_TYPE
!===============================================================================

  subroutine allocate_material_type(this,ng)

!---external references

    use constants, only: ZERO

!---arguments

    integer             :: ng
    type(material_type) :: this

!---begin execution

    ! allocate xs and diffusion coefficients
    allocate(this % totalxs(ng))
    allocate(this % absorxs(ng))
    allocate(this % removxs(ng))
    allocate(this % fissvec(ng))
    allocate(this % scattxs(ng,ng))
    allocate(this % nfissxs(ng,ng))
    allocate(this % diffcof(ng))
    allocate(this % chi(ng))

    ! zero out everything
    this % totalxs = ZERO 
    this % absorxs = ZERO
    this % removxs = ZERO
    this % fissvec = ZERO
    this % scattxs = ZERO
    this % nfissxs = ZERO
    this % diffcof = ZERO
    this % chi     = ZERO

  end subroutine allocate_material_type

!===============================================================================
! DEALLOCATE_MATERIAL_TYPE
!===============================================================================

  subroutine deallocate_material_type(this)

!---arguments

    type(material_type) :: this

!---begin execution

    ! deallocate all
    deallocate(this % totalxs)
    deallocate(this % absorxs)
    deallocate(this % removxs)
    deallocate(this % fissvec)
    deallocate(this % scattxs)
    deallocate(this % nfissxs)
    deallocate(this % diffcof)
    deallocate(this % chi)

  end subroutine deallocate_material_type

end module material_header

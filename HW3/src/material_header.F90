module material_header

!-module options

  implicit none
  private
  public :: allocate_material_type, deallocate_material_type, arrange_xs

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
    real(8), allocatable :: chip(:)
    real(8), allocatable :: chid(:)

    ! axial buckling
    real(8) :: buckling

    ! kinetics factors
    real(8), allocatable :: kinrem(:)
    real(8), allocatable :: kinfis(:)

    ! logicals for preprocessing
    logical :: abs_based = .false.
    logical :: rem_based = .false.
    logical :: chi_based = .false.
    logical :: use_buckling = .false.

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
    allocate(this % chid(ng))
    allocate(this % chip(ng))

    ! zero out everything
    this % totalxs = ZERO 
    this % absorxs = ZERO
    this % removxs = ZERO
    this % fissvec = ZERO
    this % scattxs = ZERO
    this % nfissxs = ZERO
    this % diffcof = ZERO
    this % chi     = ZERO
    this % chid    = ZERO
    this % chip    = ZERO

  end subroutine allocate_material_type

!===============================================================================
! ARRANGE_XS
!===============================================================================

  subroutine arrange_xs(this,ng)

!---arguments

    integer :: ng
    type(material_type) :: this

!---local variables

    integer :: g
    integer :: h

!---begin execution

    ! check for buckling
    if (this % use_buckling) this % absorxs = this % absorxs +                 &
                             this % diffcof*this % buckling

    ! begin loop of target energy
    GROUP: do g = 1, ng

      ! check if absorption based
      if (this % abs_based) then

        ! set removal based to true
        this % rem_based = .true.

        ! compute removal
        this % removxs(g) = this % absorxs(g) +                  &
                            sum(this % scattxs(:,g)) -           & 
                            this % scattxs(g,g)
      end if

      ! check if removal based
      if (this % rem_based) then
        this % scattxs(g,g) = this % totalxs(g) - this % removxs(g)
      end if

      ! begin loop over outgoing energy
      if (this % chi_based) then
        do h = 1, ng
          this % nfissxs(h,g) = this % chi(h) * this % fissvec(g)
        end do
      end if

    end do GROUP 

  end subroutine arrange_xs

!===============================================================================
! DEALLOCATE_MATERIAL_TYPE
!===============================================================================

  subroutine deallocate_material_type(this)

!---arguments

    type(material_type) :: this

!---begin execution

    ! deallocate all
    if (allocated(this % totalxs)) deallocate(this % totalxs)
    if (allocated(this % absorxs)) deallocate(this % absorxs)
    if (allocated(this % removxs)) deallocate(this % removxs)
    if (allocated(this % fissvec)) deallocate(this % fissvec)
    if (allocated(this % scattxs)) deallocate(this % scattxs)
    if (allocated(this % nfissxs)) deallocate(this % nfissxs)
    if (allocated(this % diffcof)) deallocate(this % diffcof)
    if (allocated(this % chi))     deallocate(this % chi)
    if (allocated(this % chid))    deallocate(this % chid)
    if (allocated(this % chip))    deallocate(this % chip)
    if (allocated(this % kinrem))  deallocate(this % kinrem)
    if (allocated(this % kinfis))  deallocate(this % kinfis)

  end subroutine deallocate_material_type

end module material_header

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

    ! axial buckling
    real(8) :: buckling

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
                             this % diffcof*this % buckling**2

    ! begin loop of target energy
    GROUP: do g = 1, ng

      ! check if absorption based
      if (this % abs_based) then

        ! set removal based to true
        this % rem_based = .true.

        ! compute removal
        this % removxs(g) = this % absorxs(g) +                  &
                            sum(this % scattxs(g,:)) -           & 
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

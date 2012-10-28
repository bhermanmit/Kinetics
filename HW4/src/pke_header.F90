module pke_header

!-external references

  use constants,  only: NUM_PRECS, ZERO

!-module options

  implicit none
  private
  public :: allocate_pke_type, deallocate_pke_type

!-type definitions

  type, public :: pke_type

    real(8), allocatable :: N(:,:)    ! power/prec vector for output
    real(8), allocatable :: coef(:,:) ! coeffcient matrix
    real(8), allocatable :: expm(:,:) ! result after matrix exponential

  end type pke_type

contains

!===============================================================================
! ALLOCATE_PKE_TYPE
!===============================================================================

  subroutine allocate_pke_type(this,ng,nt)

!---arguments

    integer :: ng
    integer :: nt
    type(pke_type) :: this

!---begin execution

    ! allocate
    if (.not.allocated(this % N)) allocate(this % N(NUM_PRECS*ng+ng,nt+1))
    if (.not.allocated(this % coef)) allocate(this % coef(NUM_PRECS*ng+ng,NUM_PRECS*ng+ng))
    if (.not.allocated(this % expm)) allocate(this % expm(NUM_PRECS*ng+ng,NUM_PRECS*ng+ng))

    ! set to 0
    this % N    = ZERO
    this % coef = ZERO
    this % expm = ZERO

  end subroutine allocate_pke_type

!===============================================================================
! DEALLOCATE_PKE_TYPE
!===============================================================================

  subroutine deallocate_pke_type(this)

!---arguments

    type(pke_type) :: this

!---begin execution

    ! deallocate
    if (allocated(this % N)) deallocate(this % N)
    if (allocated(this % coef)) deallocate(this % coef)
    if (allocated(this % expm)) deallocate(this % expm)

  end subroutine deallocate_pke_type

end module pke_header 

module data_header

!-external references

  use constants,  only: NUM_PRECS

!-module options

  implicit none
  private
  public :: allocate_pke_type, deallocate_pke_type

!-type definitions

  type, public :: pke_type

    integer              :: npts         ! number of points in inputvector
    integer              :: idx=1        ! index in input vectors for interp
    integer, allocatable :: nt(:)        ! number of time steps
    real(8), allocatable :: dt(:)        ! delta time
    real(8), allocatable :: t(:)         ! time vector
    real(8), allocatable :: rho(:)       ! reactivity vector
    real(8), allocatable :: time(:)      ! time vector for output
    real(8), allocatable :: N(:,:)       ! power/prec vector for output
    real(8), allocatable :: react(:)     ! reactivity vector for output
    real(8), allocatable :: refpower(:)  ! reference power if restart
    real(8) :: coef(NUM_PRECS+1,NUM_PRECS+1) ! coeffcient matrix
    real(8) :: expm(NUM_PRECS+1,NUM_PRECS+1) ! result after matrix exponential

  end type pke_type

contains

!===============================================================================
! ALLOCATE_PKE_TYPE
!===============================================================================

  subroutine allocate_pke_type(this)

!---arguments

    type(pke_type) :: this

!---begin execution

    ! allocate
    if (.not.allocated(this % time))  allocate(this % time(sum(this % nt)+1))
    if (.not.allocated(this % N))     allocate(this % N(NUM_PRECS+1,sum(this % nt)+1))
    if (.not.allocated(this % react)) allocate(this % react(sum(this % nt)+1))

    ! set to 0
    this % time  = 0.0_8
    this % N     = 0.0_8
    this % react = 0.0_8
    this % coef  = 0.0_8
    this % expm  = 0.0_8

  end subroutine allocate_pke_type

!===============================================================================
! DEALLOCATE_PKE_TYPE
!===============================================================================

  subroutine deallocate_pke_type(this)

!---arguments

    type(pke_type) :: this

!---begin execution

    ! deallocate
    if (allocated(this % t))        deallocate(this % t)
    if (allocated(this % rho))      deallocate(this % rho)
    if (allocated(this % time))     deallocate(this % time)
    if (allocated(this % N))        deallocate(this % N)
    if (allocated(this % react))    deallocate(this % react)
    if (allocated(this % dt))       deallocate(this % dt)
    if (allocated(this % nt))       deallocate(this % nt)
    if (allocated(this % refpower)) deallocate(this % refpower)

  end subroutine deallocate_pke_type

end module data_header

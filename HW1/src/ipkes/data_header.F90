module data_header

!-external references

  use constants,  only: NUM_PRECS

!-module options

  implicit none
  private
  public :: allocate_ipke_type, deallocate_ipke_type

!-type definitions

  type, public :: ipke_type

    integer              :: npts         ! number of points in inputvector
    integer              :: idx=1        ! index in input vectors for interp
    integer, allocatable :: nt(:)        ! number of time steps
    real(8), allocatable :: dt(:)        ! delta time
    real(8), allocatable :: t(:)         ! time vector
    real(8), allocatable :: power(:)     ! power vector 
    real(8), allocatable :: time(:)      ! time vector for output
    real(8), allocatable :: N(:,:)       ! power/prec vector for output
    real(8), allocatable :: react(:)     ! reactivity vector for output

  end type ipke_type

contains

!===============================================================================
! ALLOCATE_PKE_TYPE
!===============================================================================

  subroutine allocate_ipke_type(this)

!---arguments

    type(ipke_type) :: this

!---begin execution

    ! allocate
    if (.not.allocated(this % time))  allocate(this % time(sum(this % nt)+1))
    if (.not.allocated(this % N))     allocate(this % N(NUM_PRECS+1,sum(this % nt)+1))
    if (.not.allocated(this % react)) allocate(this % react(sum(this % nt)+1))

    ! set to 0
    this % time  = 0.0_8
    this % N     = 0.0_8
    this % react = 0.0_8

  end subroutine allocate_ipke_type

!===============================================================================
! DEALLOCATE_PKE_TYPE
!===============================================================================

  subroutine deallocate_ipke_type(this)

!---arguments

    type(ipke_type) :: this

!---begin execution

    ! deallocate
    if (allocated(this % t))     deallocate(this % t)
    if (allocated(this % power)) deallocate(this % power)
    if (allocated(this % time))  deallocate(this % time)
    if (allocated(this % N))     deallocate(this % N)
    if (allocated(this % react)) deallocate(this % react)
    if (allocated(this % dt))    deallocate(this % dt)
    if (allocated(this % nt))    deallocate(this % nt)

  end subroutine deallocate_ipke_type

end module data_header

module pke_header

!-module options

  implicit none
  private
  public :: allocate_pke_type, deallocate_pke_type

!-module variables

  type, public :: pke_type

    integer :: n
    integer :: idx = 1
    real(8) :: rhot
    real(8), allocatable :: rho(:)
    real(8), allocatable :: time(:)

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
    allocate(this % rho(this % n))
    allocate(this % time(this % n))

  end subroutine allocate_pke_type

!===============================================================================
! DEALLOCATE_PKE_TYPE
!===============================================================================

  subroutine deallocate_pke_type(this)

!---arguments

    type(pke_type) :: this

!---begin execution

    ! deallocate all
    deallocate(this % rho)
    deallocate(this % time)

  end subroutine deallocate_pke_type

end module pke_header

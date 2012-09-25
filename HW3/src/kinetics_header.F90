module kinetics_header

!-module options

  implicit none
  private
  public :: allocate_kinetics_type, deallocate_kinetics_type

!-module variables

  type, public :: kinetics_type

    character(len=15) :: xs_id
    integer :: mat_id
    integer :: h, g
    integer :: idxt = 1
    real(8), allocatable :: val(:)
    real(8), allocatable :: time(:)    

  end type kinetics_type

contains

!===============================================================================
! ALLOCATE_KINETICS_TYPE
!===============================================================================

  subroutine allocate_kinetics_type(this,n)

!---arguments

    integer :: n
    type(kinetics_type) :: this

!---begin execution

    allocate(this % val(n))
    allocate(this % time(n))

  end subroutine allocate_kinetics_type

!===============================================================================
! DEALLOCATE_CMFD_TYPE
!===============================================================================

  subroutine deallocate_kinetics_type(this)

!---arguments

    type(kinetics_type) :: this

!---begin execution

    deallocate(this % val)
    deallocate(this % time)

  end subroutine deallocate_kinetics_type

end module kinetics_header

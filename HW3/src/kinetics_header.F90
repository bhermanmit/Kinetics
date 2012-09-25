module kinetics_header

!-module options

  implicit none
  private
  public :: allocate_kinetics_type, deallocate_kinetics_type

!-module variables

  type, public :: kinetics_type

    character(len=15) :: xs_id
    integer :: mat_id
    real(8), allocatable :: val(:)
    real(8), allocatable :: time(:)    

  end type kinetics_type

contains

!===============================================================================
! ALLOCATE_KINETICS_TYPE
!===============================================================================

  subroutine allocate_kinetics_type(this)

!---external references

!---arguments

    type(kinetics_type) :: this

!---begin execution

  end subroutine allocate_kinetics_type

!===============================================================================
! DEALLOCATE_CMFD_TYPE
!===============================================================================

  subroutine deallocate_kinetics_type(this)

!---arguments

    type(kinetics_type) :: this

!---begin execution

  end subroutine deallocate_kinetics_type

end module kinetics_header

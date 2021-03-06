module error

  use, intrinsic :: ISO_FORTRAN_ENV

  use global,  only: message, master

  implicit none

contains

!===============================================================================
! WARNING issues a warning to the user in the log file and the standard output
! stream.
!===============================================================================

  subroutine warning()

    integer :: n_lines ! number of lines
    integer :: i       ! loop index for lines

    ! Only allow master to print to screen
    if (.not. master) return

    write(OUTPUT_UNIT, fmt='(1X,A9)', advance='no') 'WARNING: '

    n_lines = (len_trim(message)-1)/70 + 1
    do i = 1, n_lines
       if (i == 1) then
          write(OUTPUT_UNIT, fmt='(A70)') message(70*(i-1)+1:70*i)
       else
          write(OUTPUT_UNIT, fmt='(10X,A70)') message(70*(i-1)+1:70*i)
       end if
    end do

  end subroutine warning

!===============================================================================
! FATAL_ERROR alerts the user that an error has been encountered and displays a
! message about the particular problem. Errors are considered 'fatal' and hence
! the program is aborted.
!===============================================================================

  subroutine fatal_error(error_code)

    integer, optional :: error_code ! error code

    integer :: code    ! error code
    integer :: n_lines ! number of lines
    integer :: i       ! loop index over lines

    ! set default error code
    if (present(error_code)) then
       code = error_code
    else
       code = -1
    end if

    write(ERROR_UNIT, fmt='(1X,A7)', advance='no') 'ERROR: '

    n_lines = (len_trim(message)-1)/72 + 1
    do i = 1, n_lines
       if (i == 1) then
          write(ERROR_UNIT, fmt='(A72)') message(72*(i-1)+1:72*i)
       else
          write(ERROR_UNIT, fmt='(7X,A72)') message(72*(i-1)+1:72*i)
       end if
    end do
    write(ERROR_UNIT,*)

    ! Release memory from all allocatable arrays
!   call free_memory()

    ! Abort program
    stop

  end subroutine fatal_error

end module error

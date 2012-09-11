module output

  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none

contains

!===============================================================================
! TITLE prints the main title banner as well as information about the program
! developers, version, and date/time which the problem was run.
!===============================================================================

  subroutine title()

    use constants,  only: VERSION_MAJOR, VERSION_MINOR
    use global,     only: n_procs

    write(UNIT=OUTPUT_UNIT, FMT='(/11(A/))') &
         '       .d88888b.                             888b     d888  .d8888b.', &
         '      d88P" "Y88b                            8888b   d8888 d88P  Y88b', &
         '      888     888                            88888b.d88888 888    888', &
         '      888     888 88888b.   .d88b.  88888b.  888Y88888P888 888       ', &
         '      888     888 888 "88b d8P  Y8b 888 "88b 888 Y888P 888 888       ', &
         '      888     888 888  888 88888888 888  888 888  Y8P  888 888    888', &
         '      Y88b. .d88P 888 d88P Y8b.     888  888 888   "   888 Y88b  d88P', &
         '       "Y88888P"  88888P"   "Y8888  888  888 888       888  "Y8888P"', &
         '__________________888______________________________________________________', &
         '                  888', &
         '                  888'

    ! Write version information
    write(UNIT=OUTPUT_UNIT, FMT=*) &
         '     Developed At:  Massachusetts Institute of Technology'
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Version:",7X,I1,".",I1)') &
         VERSION_MAJOR, VERSION_MINOR
#ifdef GIT_SHA1
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Git SHA1:",6X,A)') GIT_SHA1
#endif

    ! Write the date and time
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Date/Time:",5X,A)') &
         time_stamp()

#ifdef MPI
    ! Write number of processors
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"MPI Processes:",1X,A)') &
         trim(to_str(n_procs))
#endif

  end subroutine title

!===============================================================================
! TIME_STAMP returns the current date and time in a formatted string
!===============================================================================

  function time_stamp() result(current_time)

    character(19) :: current_time ! ccyy-mm-dd hh:mm:ss
    character(8)  :: date_        ! ccyymmdd
    character(10) :: time_        ! hhmmss.sss

    call date_and_time(DATE=date_, TIME=time_)
    current_time = date_(1:4) // "-" // date_(5:6) // "-" // date_(7:8) // &
         " " // time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)

  end function time_stamp

!===============================================================================
! HEADER displays a header block according to a specified level. If no level is
! specified, it is assumed to be a minor header block (H3).
!===============================================================================

  subroutine header(msg, unit, level)

    use constants,  only: MAX_LINE_LEN
    use string,     only: upper_case

    character(*), intent(in) :: msg ! header message
    integer, optional :: unit       ! unit to write to
    integer, optional :: level      ! specified header level

    integer :: n            ! number of = signs on left
    integer :: m            ! number of = signs on right
    integer :: unit_        ! unit to write to
    integer :: header_level ! actual header level
    character(MAX_LINE_LEN) :: line

    ! set default level
    if (present(level)) then
       header_level = level
    else
       header_level = 3
    end if

    ! set default unit
    if (present(unit)) then
       unit_ = unit
    else
       unit_ = OUTPUT_UNIT
    end if

    ! determine how many times to repeat '=' character
    n = (63 - len_trim(msg))/2
    m = n
    if (mod(len_trim(msg),2) == 0) m = m + 1

    ! convert line to upper case
    line = msg
    call upper_case(line)

    ! print header based on level
    select case (header_level)
    case (1)
       write(UNIT=unit_, FMT='(/3(1X,A/))') repeat('=', 75), & 
            repeat('=', n) // '>     ' // trim(line) // '     <' // &
            repeat('=', m), repeat('=', 75)
    case (2)
       write(UNIT=unit_, FMT='(/2(1X,A/))') trim(line), repeat('-', 75)
    case (3)
       write(UNIT=unit_, FMT='(/1X,A/)') repeat('=', n) // '>     ' // &
            trim(line) // '     <' // repeat('=', m)
    end select

  end subroutine header

!===============================================================================
! WRITE_MESSAGE displays an informational message to the log file and the 
! standard output stream.
!===============================================================================

  subroutine write_message(level)

    use global,  only: master, message, verbosity

    integer, optional :: level ! verbosity level

    integer :: n_lines ! number of lines needed
    integer :: i       ! index for lines

    ! Only allow master to print to screen
    if (.not. master .and. present(level)) return

    ! TODO: Take care of line wrapping so words don't get cut off
    if (.not. present(level) .or. level <= verbosity) then
       n_lines = (len_trim(message)-1)/79 + 1
       do i = 1, n_lines
          write(OUTPUT_UNIT, fmt='(1X,A)') trim(message(79*(i-1)+1:79*i))
       end do
    end if

  end subroutine write_message

end module output

module output

!-external references

  use, intrinsic :: ISO_FORTRAN_ENV

!-module options

  implicit none
  private
  public :: title, header, write_results, write_message, write_hdf5

contains

!===============================================================================
! TITLE prints the main title banner as well as information about the program
! developers, version, and date/time which the problem was run.
!===============================================================================

  subroutine title()

    use constants,  only: VERSION_MAJOR, VERSION_MINOR
    use global,     only: n_procs
    use string,     only: to_str

    write(UNIT=OUTPUT_UNIT, FMT='(/11(A/))') &
    & '    ,o888888o.           ,8.       ,8.          8 8888888888   8 888888888o.      ', &
    & '   8888     `88.        ,888.     ,888.         8 8888         8 8888    `^888.   ', &
    & ',8 8888       `8.      .`8888.   .`8888.        8 8888         8 8888        `88. ', &
    & '88 8888               ,8.`8888. ,8.`8888.       8 8888         8 8888         `88 ', &
    & '88 8888              ,8"8.`8888,8^8.`8888.      8 888888888888 8 8888          88 ', &
    & '88 8888             ,8" `8.`8888" `8.`8888.     8 8888         8 8888          88 ', &
    & '88 8888            ,8"   `8.`88"   `8.`8888.    8 8888         8 8888         ,88 ', &
    & '`8 8888       .8" ,8"     `8.`"     `8.`8888.   8 8888         8 8888        ,88" ', &
    & '   8888     ,88" ,8"       `8        `8.`8888.  8 8888         8 8888    ,o88P"   ', &
    & '    `8888888P"  ,8"         `         `8.`8888. 8 8888         8 888888888P"      ', &
    & '__________________________________________________________________________________'


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

!===============================================================================
! WRITE_RESULTS
!===============================================================================

  subroutine write_results()

!---external references

    use global,  only: time_total, time_init, time_build, time_power,          &
                       time_inner, cmfd

!---begin execution

    ! write timing
    call header("Timing Statistics", level = 3)
    write(OUTPUT_UNIT, 100) 'Time for initialization', time_init % elapsed
    write(OUTPUT_UNIT, 100) 'Time for building matrices', time_build % elapsed
    write(OUTPUT_UNIT, 100) 'Time for source convergence', time_power % elapsed
    write(OUTPUT_UNIT, 100) 'Time for inner iterations', time_inner % elapsed
    write(OUTPUT_UNIT, 100) 'Total simulation time', time_total % elapsed

    ! write results
    call header("Results",level=3)
    write(OUTPUT_UNIT, 201) 'Final k-effective', cmfd % keff
    write(OUTPUT_UNIT, 202) 'L-2 norm of nodal power', cmfd % norm
    write(OUTPUT_UNIT, 200) 'Number of iterations', cmfd % iter
    write(OUTPUT_UNIT, 202) 'Dominance ratio', cmfd % dr
    write(OUTPUT_UNIT, fmt='(/,A)') 'Simulation Finished.'

    ! fomat for write statements
100 format (1X,A,T35,"= ",ES11.4," seconds")
200 format (1X,A,T35,"=  ",I0)
201 format (1X,A,T35,"=  ",F8.6)
202 format (1X,A,T35,"= ",ES11.4)

  end subroutine write_results

!===============================================================================
! WRITE_HDF5
!===============================================================================

  subroutine write_hdf5()

!---external references

    use global,  only: cmfd, geometry, ktol, stol, itol, hdf5_err
    use hdf5_interface

!---local variables

    integer :: g,i,j,k
    integer :: n

!---begin execution

    ! get full size
    n = geometry % nfx * geometry % nfy * geometry % nfz * geometry % nfg

    ! create a new file
    call hdf5_create_file('output.h5')

    ! write out geometry
    call hdf5_make_integer(hdf5_output_file, "ncx", geometry % ncx)
    call hdf5_make_integer(hdf5_output_file, "ncy", geometry % ncy)
    call hdf5_make_integer(hdf5_output_file, "ncz", geometry % ncz)
    call hdf5_make_integer(hdf5_output_file, "ncg", geometry % ncg)
    call hdf5_make_integer(hdf5_output_file, "nfx", geometry % nfx)
    call hdf5_make_integer(hdf5_output_file, "nfy", geometry % nfy)
    call hdf5_make_integer(hdf5_output_file, "nfz", geometry % nfz)
    call hdf5_make_integer(hdf5_output_file, "nfg", geometry % nfg)

    ! write out the mat
    call hdf5_make_array(hdf5_output_file,"mat",geometry % fmat_map,n)

    ! write out the reg
    call hdf5_make_array(hdf5_output_file,"reg",geometry % freg_map,n)

    ! write out tolerances
    call hdf5_make_double(hdf5_output_file, "ktol", ktol)
    call hdf5_make_double(hdf5_output_file, "stol", stol)
    call hdf5_make_double(hdf5_output_file, "itol", itol)

    ! write out iteration count
    call hdf5_make_integer(hdf5_output_file, "iterations", cmfd % iter)

    ! write out norm
    call hdf5_make_double(hdf5_output_file, "norm", cmfd % norm)

    ! write out dominance ratio
    call hdf5_make_double(hdf5_output_file, "dr", cmfd % dr)

    ! write out keff
    call hdf5_make_double(hdf5_output_file, "keff", cmfd % keff)

    ! write out phi
    call hdf5_make_array(hdf5_output_file, "phi", cmfd % phi, n)

    ! write out power
    call hdf5_make_array(hdf5_output_file, "power", cmfd%power_n,              &
                         size(cmfd%power_n))

    ! close file
    call hdf5_close_file()

  end subroutine write_hdf5

end module output

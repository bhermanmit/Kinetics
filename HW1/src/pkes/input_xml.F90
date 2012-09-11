module input_xml

  implicit none

contains

!===============================================================================
! READ_INPUT_XML
!===============================================================================

  subroutine read_input_xml()

!---external references

    use constants,        only: MAX_FILE_LEN
    use error,            only: fatal_error
    use global,           only: message, pke
    use output,           only: write_message
    use xml_data_input_t

    integer                 :: n
    logical                 :: file_exists
    character(MAX_FILE_LEN) :: filename

    ! Display output message
    message = "Reading settings XML file..."
    call write_message(5)

    ! Check if settings.xml exists
    filename = "input.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
       message = "Settings XML file '" // trim(filename) // "' does not exist!"
       call fatal_error()
    end if

    ! read in input file
    call read_xml_file_input_t(filename)

    ! get size of input vectors
    n = size(reactivity_ % time)
    if (n /= size(reactivity_ % rho)) then
      message = "Time and Rho vectors not of same size!"
      call fatal_error()
    end if

    ! allocate vectors
    allocate(pke % t(n))
    allocate(pke % rho(n))

    ! save in object
    pke % npts = n
    pke % t    = reactivity_ % time
    pke % rho  = reactivity_ % rho
    pke % maxt = maxtime_
    pke % dt   = timestep_ 

  end subroutine read_input_xml

end module input_xml

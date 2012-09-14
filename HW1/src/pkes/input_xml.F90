module input_xml

  implicit none

contains

!===============================================================================
! READ_BINARY
!===============================================================================

  subroutine read_binary()

!---external references

    use global,  only: pke

!---local variables

    integer :: i,sz

!---begin execution

    ! open binary file
    open(10,file='output.bin',form='unformatted')

    ! read size 
    read(10) sz

    ! allocate vectors
    allocate(pke % t(sz))
    allocate(pke % rho(sz))
    allocate(pke % refpower(sz))
    allocate(pke % dt(sz-1))

    ! read vectors
    read(10) pke % t
    read(10) pke % rho
    read(10) pke % dt
    read(10) pke % refpower

    ! set size
    pke % npts = sz

  end subroutine

!===============================================================================
! READ_INPUT_XML
!===============================================================================

  subroutine read_input_xml()

!---external references

    use constants,        only: MAX_FILE_LEN
    use error,            only: fatal_error
    use global,           only: message, pke, restart
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

    ! check for restart
    restart = restart_
    if (restart) then
      call read_binary()
      return
    end if

    ! get size of input vectors
    n = size(time_)
    if (n /= size(rho_)) then
      message = "Time and Rho vectors not of same size!"
      call fatal_error()
    end if

    ! allocate vectors
    allocate(pke % t(n))
    allocate(pke % rho(n))
    allocate(pke % dt(n-1))

    ! save in object
    pke % npts = n
    pke % t    = time_
    pke % rho  = rho_
    pke % dt   = timestep_

  end subroutine read_input_xml

end module input_xml

module constants

!-major, minor version numbers
  integer, parameter :: VERSION_MAJOR   = 0
  integer, parameter :: VERSION_MINOR   = 1

!-Numbers

  real(8), parameter :: ZERO = 0.0_8
  real(8), parameter :: ONE  = 1.0_8

!-Files

  integer, parameter :: MAX_WORDS    = 500
  integer, parameter :: MAX_LINE_LEN = 250
  integer, parameter :: MAX_WORD_LEN = 150
  integer, parameter :: MAX_FILE_LEN = 255

!-Codes for read errors

  integer, parameter :: ERROR_INT  = -huge(0)
  real(8), parameter :: ERROR_REAL = -huge(0.0_8) * 0.917826354_8

end module constants

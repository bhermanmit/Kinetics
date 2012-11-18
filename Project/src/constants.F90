module constants

!-major, minor version numbers
  integer, parameter :: VERSION_MAJOR   = 0
  integer, parameter :: VERSION_MINOR   = 1

!-Numbers

  real(8), parameter :: ZERO = 0.0_8
  real(8), parameter :: ONE  = 1.0_8
  real(8), parameter :: TWO  = 2.0_8

!-Files

  integer, parameter :: MAX_WORDS    = 500
  integer, parameter :: MAX_LINE_LEN = 250
  integer, parameter :: MAX_WORD_LEN = 150
  integer, parameter :: MAX_FILE_LEN = 255

!-Codes for read errors

  integer, parameter :: ERROR_INT  = -huge(0)
  real(8), parameter :: ERROR_REAL = -huge(0.0_8) * 0.917826354_8

!-Kinetics parameters

  ! Number of precursor groups
  integer, parameter :: NUM_PRECS = 8
! integer, parameter :: NUM_PRECS = 6

  ! Decay constant
  real(8), parameter :: lambda(8) = (/0.012467_8, &
                                      0.028292_8, &
                                      0.042524_8, &
                                      0.133042_8, &
                                      0.292467_8, &
                                      0.666488_8, &
                                      1.634781_8, &
                                      3.554600_8/)

  ! Delayed Neutron Fraction
  real(8), parameter :: beta(8)   = (/0.000218_8, &
                                      0.001023_8, &
                                      0.000605_8, &
                                      0.001310_8, &
                                      0.002200_8, &
                                      0.000600_8, &
                                      0.000540_8, &
                                      0.000152_8/)

! real(8), parameter :: beta(6)   = (/2.1450E-4_8, &
!                                     1.4235E-3_8, &
!                                     1.2740E-3_8, &
!                                     2.5675E-3_8, &
!                                     7.4750E-4_8, &
!                                     2.7300E-4_8/)

! real(8), parameter :: lambda(6) = (/0.012400_8, &
!                                     0.030500_8, &
!                                     0.111000_8, &
!                                     0.301000_8, &
!                                     1.140000_8, &
!                                     3.010000_8/)


  ! Prompt Neutron Lifetime
  real(8), parameter :: vel(2) = (/2200._8*100._8*(0.100e4_8/0.0253_8)**0.5_8, & 
                                  2200._8*100._8*(0.100_8/0.0253_8)**0.5_8/)

  ! Prompt Neutron Lifetime
  real(8), parameter :: sigf = 0.05_8
  real(8), parameter :: nu   = 2.45_8
  real(8), parameter :: vel1 = 2200._8*100._8*(0.100_8/0.0253_8)**0.5_8
  real(8)            :: pnl  = ONE/(sigf*nu*vel1)
! real(8) :: pnl = 1.0E-3_8

end module constants

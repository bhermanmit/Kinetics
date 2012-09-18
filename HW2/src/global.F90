module global

  use constants,   only: MAX_LINE_LEN
  use cmfd_header, only: cmfd_obj
  use timing,      only: Timer

  implicit none
  save

!-Main object

  type(cmfd_obj) :: cmfd

!-Timing objects

  type(Timer) :: time_total  ! timer for whole calculation
  type(Timer) :: time_build  ! timer for mat building
  type(Timer) :: time_solve  ! timer for power iteration

!-petsc error code

  integer :: ierr

!-mpi parametesr

  logical :: master = .false. ! am i master
  integer :: rank             ! rank of processor
  integer :: n_procs          ! number of processors
  integer :: mpi_err          ! error code

!-solver type

  character(len=25) :: solver_type

!-Message used in message/warning/fatal_error

  character(MAX_LINE_LEN) :: message

!-The verbosity controls how much information will be printed to the
!-screen and in logs

  integer :: verbosity = 7

end module global

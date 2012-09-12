module global

  use constants,    only: MAX_LINE_LEN
  use data_header,  only: pke_type
  use timing,       only: Timer

  implicit none
  save

!-timing variables

  type(Timer) :: total_time       ! timer for total run

!-Message used in message/warning/fatal_error

  character(MAX_LINE_LEN) :: message

!-The verbosity controls how much information will be printed to the
!-screen and in logs

  integer :: verbosity = 7

!-Parallel processing variables

  integer :: n_procs     = 1       ! number of processes
  integer :: rank        = 0       ! rank of process
  logical :: master      = .true.  ! master process?
  logical :: mpi_enabled = .false. ! is MPI in use and initialized?
  integer :: mpi_err               ! MPI error code

!-Point Kinetics object

  type(pke_type) :: pke

end module global

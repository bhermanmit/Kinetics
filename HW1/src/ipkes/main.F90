program main

!-external references

  use finalize,    only: finalize_run
  use initialize,  only: initialize_run
  use physics,     only: run_invkinetics

!-program options

  implicit none

!-begin program

  ! initialize problem
  call initialize_run()

  ! solve inverse point kinetics equations
  call run_invkinetics()

  ! finalize problem
  call finalize_run()

end program main

program main

!-external references

  use finalize,    only: finalize_run
  use initialize,  only: initialize_run
  use physics,     only: run_kinetics

!-program options

  implicit none

!-begin program

  ! initialize problem
  call initialize_run()

  ! solve point kinetics equations
  call run_kinetics()

  ! finalize problem
  call finalize_run()

end program main

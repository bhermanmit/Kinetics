program main

!-external references

  use finalize,    only: finalize_run
  use initialize,  only: initialize_run
  use input_xml,   only: read_input_xml
  use physics,     only: run_kinetics

!-program options

  implicit none

!-begin program

  ! read in input
  call read_input_xml()

  ! initialize problem
  call initialize_run()

  ! solve point kinetics equations
  call run_kinetics()

  ! plot results

  ! finalize problem
  call finalize_run()

end program main

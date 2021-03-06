program main

!-external references

  use cmfd_execute,  only: execute_cmfd 
  use finalize,      only: finalize_run
  use initialize,    only: initialize_run 
  use output,        only: title

!-program options

  implicit none

!-begin execution

  ! print out title
  call title()

  ! initialize
  call initialize_run()

  ! execute cmfd
  call execute_cmfd()

  ! finalize run
  call finalize_run()

end program main

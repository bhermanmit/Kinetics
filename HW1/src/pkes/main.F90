program main

!-external references

  use expokit
  use gnufor2, only: plot

!-program options

  implicit none

!-begin program

  real(kind=8) :: x(4),y(4)
  real(8) :: A(2,2)

  print *, "Hello World!"

  A(1,1) = 1.0
  A(2,1) = 2.0
  A(1,2) = 3.0
  A(2,2) = 4.0

  call dense_pade(A,2,1.0_8)

  print *,A

  x = (/1.0,2.0,3.0,4.0/)
  y = 2.0*x + 4.0

  call plot(x,y)

end program main

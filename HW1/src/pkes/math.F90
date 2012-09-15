module math 

!-module options

  implicit none
  private
  public :: norm_inf, ilog2

contains

!===============================================================================
! NORM_INF 
!===============================================================================

  function norm_inf(A,n) result(norm)

!---arguments

     integer :: n
     real(8) :: A(n,n)
     real(8) :: norm

!---begin execution

    ! compute max of row sums
    norm = maxval(sum(abs(A),2))

  end function norm_inf

!===============================================================================
! ILOG2
!===============================================================================

  function ilog2(x) result(e)

!---arguments

    integer :: e
    real(8) :: x

!---local variables

    integer :: i

!---begin execution

    ! begin loop to compute exponent
    do i = 1,15

      if ( (x/(2.0_8**i)) < 1.0_8 ) then
        e = i
        return
      end if

    end do 

  end function ilog2 

end module math 

module expokit

contains

!===============================================================================
! DENSE_PADE
!===============================================================================

  subroutine dense_pade(A,m,t)

!---arguments

    integer :: m
    real(8) :: A(m,m)
    real(8) :: t

!---local variables

    integer :: ideg = 6
    integer :: lwsp
    integer :: iexp
    integer :: ns
    integer :: iflag
    integer, allocatable :: iwsp(:)
    real(8), allocatable :: wsp(:)

!---begin execution

    ! compute size of output and allocate
    lwsp = 4*m*m + ideg + 1
    allocate(wsp(lwsp))
    allocate(iwsp(m))

    ! compute matrix exponential
    call DGPADM( ideg, m, t, A, m, wsp, lwsp, iwsp, iexp, ns, iflag ) 

    ! set new matrix
    A = reshape(wsp(iexp:iexp+m*m-1),(/m,m/))

    ! deallocate all memory
    deallocate(wsp)
    deallocate(iwsp)

  end subroutine

end module expokit

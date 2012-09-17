module solvers 

!-module options

  implicit none
  private
  public :: expm_pade

contains

!===============================================================================
! EXPM_PADE 
!===============================================================================

  subroutine expm_pade(A,N,dt,EXPM) 

!---external references

    use constants,  only: ZERO, ONE
    use math,       only: norm_inf, ilog2

!---arguments

    integer :: N          ! dimension of square matrix
    real(8) :: A(N,N)     ! the coefficient matrix
    real(8) :: dt         ! the time step
    real(8) :: EXPM(N,N)  ! output matrix

!---local variables

    integer :: ex   ! this is the integer exponent that goes into th scaling
    integer :: s    ! this is the scaling parameter
    integer :: i    ! loop counter
    integer :: j    ! Pade order counter
    integer :: q    ! Pade approximation order
    integer :: info ! did LAPACK work?
    integer, allocatable :: IPIV(:) ! pivot elements from LU factorization
    real(8) :: norm   ! the norm of the coefficient matrix for scaling
    real(8) :: c      ! Pade coefficient
    real(8) :: shift  ! -1 or +1 and a function of J for D term
    real(8), allocatable :: At(:,:)  ! the exponential being approximated
    real(8), allocatable :: X(:,:)   ! accumulates multiplication of A
    real(8), allocatable :: Xt(:,:)  ! temporary matrix for X
    real(8), allocatable :: E(:,:)   ! numerator of Pade
    real(8), allocatable :: D(:,:)   ! denominator of Pade

!---begin execution

    ! allocate temporary matrices
    allocate(At(N,N))
    allocate(X(N,N))
    allocate(Xt(N,N))
    allocate(E(N,N))
    allocate(D(N,N))
    allocate(IPIV(N))

    ! set matrix
    At = A

    ! get matrix to be approximated in exponential 
    At = At * dt

    ! compute norm for scaling
    norm = norm_inf(At,N)

    ! compute integer exponent of formula norm = f*2**(ex)
    ex = ilog2(norm)

    ! determine scaling parameter
    s = max(0,ex+1)

    ! scale matrix by scaling parameter
    At = At / (2**s)

    ! begin pade approximation for numerator and denominator of rational fun.
    ! Pade order 6, diagonal approximates
    q = 6

    ! can easily set j=0 and j=1 parts right away
    ! begin loop to set the A^0 terms in each E and D, not this is a diag of ones
    E = ZERO
    D = ZERO
    do i=1,N
      E(i,i) = ONE
      D(i,i) = ONE
    end do
    c = ONE
    c = 0.5_8    ! this is j=1 Pade coefficient (need to do this before loop)
    E = E + c*At ! accmulate in sum
    D = D - c*At ! accumulate in sum

    ! begin rest of pade loop
    X = At       ! set A to X, X will accumulated the A^j in the Pade eq.
    shift = ONE ! for order 2 term, the shift for denominator is +1
    do j = 2,q

      ! compute Pade coefficient
      c = c * dble(q-j+1) / dble(j*(2*q-j+1))

      ! multiply A^(j-1) by LAPACK ROUTINE
      call DGEMM('N','N',N,N,N,ONE,At,N,X,N,ZERO,Xt,N)

      ! move temporary to actual
      X = Xt

      ! accumulate into numerator and denominator
      E = E + c*X
      D = D + shift*c*X
      shift = -ONE

    end do

    ! compute rational function (E is R after the next two lines)
    call DGETRF(N,N,D,N,IPIV,info)          ! LU factorization of D, LAPACK
    call DGETRS('N',N,N,D,N,IPIV,E,N,info)  ! solves E = D^-1*E, LAPACK

    ! perform squaring
    do i = 1,s
      call DGEMM('N','N',N,N,N,ONE,E,N,E,N,ZERO,Xt,N) ! does E*E LAPACK
      E = Xt
    end do

    ! set output matrix
    EXPM  = E  ! E is really R, didnt want to make a separate matrix

    ! deallocate temporary matrices
    deallocate(At)
    deallocate(X)
    deallocate(Xt)
    deallocate(E)
    deallocate(D)
    deallocate(IPIV)

  end subroutine expm_pade

end module solvers 

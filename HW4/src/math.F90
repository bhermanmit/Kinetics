module math 

!-module options

  implicit none
  private
  public :: sort_csr, csr_matvec_mult, csr_jacobi, csr_gauss_seidel, expm_pade

contains

!===============================================================================
! SORT
!===============================================================================

  recursive subroutine sort_csr(row, col, val, first, last)

!---arguments

    integer :: row(:)
    integer :: col(:)
    integer :: first
    integer :: last
    real(8) :: val(:)

!---local variables

    integer :: mid

!---begin execution

    if (first < last) then
      call split(row, col, val, first, last, mid)     ! split it
      call sort_csr(row, col, val, first, mid-1)      ! sort left half
      call sort_csr(row, col, val, mid+1, last)       ! sort right half
    end if

  end subroutine sort_csr 

!===============================================================================
! SPLIT
!===============================================================================

  subroutine split(row, col, val, low, high, mid)

!---arguments

    integer :: row(:)
    integer :: col(:)
    integer :: low 
    integer :: high
    integer :: mid
    real(8) :: val(:)

!g---local variables

    integer :: left
    integer :: right 
    integer :: iswap
    integer :: pivot
    integer :: row0
    real(8) :: rswap
    real(8) :: val0

!---begin execution

    left = low
    right = high
    pivot = col(low)
    row0 = row(low)
    val0 = val(low)

    ! repeat the following while left and right havent met
    do while (left < right)

      ! scan right to left to find element < pivot
      do while (left < right .and. col(right) >= pivot)
        right = right - 1
      end do

      ! scan left to right to find element > pivot
      do while (left < right .and. col(left) <= pivot)
        left = left + 1
      end do

      ! if left and right havent met, exchange the items
      if (left < right) then
        iswap = col(left)
        col(left) = col(right)
        col(right) = iswap
        iswap = row(left)
        row(left) = row(right)
        row(right) = iswap
        rswap = val(left)
        val(left) = val(right)
        val(right) = rswap
      end if

    end do

    ! swith the element in split position with pivot
    col(low) = col(right)
    col(right) = pivot
    mid = right
    row(low) = row(right)
    row(right) = row0
    val(low) = val(right)
    val(right) = val0

  end subroutine split

!===============================================================================
! CSR_MATVEC_MULT
!===============================================================================

  function csr_matvec_mult(row,col,val,x,n) result(y)

!---external references

    use constants,  only: ZERO

!---arguments

    integer :: n
    integer :: row(:)
    integer :: col(:)
    real(8) :: val(:)
    real(8) :: x(n)
    real(8) :: y(n)

!---local variables

    integer :: i
    integer :: j

!---begin execution

    ! begin loop around rows
    ROWS: do i = 1, n

      ! initialize target location in vector
      y(i) = ZERO

      ! begin loop around columns
      COLS: do j = row(i), row(i+1) - 1

        y(i) = y(i) + val(j)*x(col(j))

      end do COLS

    end do ROWS

  end function csr_matvec_mult

!===============================================================================
! CSR_JACOBI
!===============================================================================

  subroutine csr_jacobi(row,col,val,diag,x,b,n,nz,tol,iter)

!---external arguments

    use constants,  only: ZERO, ONE
    use global,     only: geometry

!---arguments

    integer, intent(in)     :: n
    integer, intent(in)     :: nz
    integer, intent(inout)  :: iter
    integer, intent(in)     :: row(n+1)
    integer, intent(in)     :: col(nz)
    integer, intent(in)     :: diag(n)
    real(8), intent(in)     :: val(nz) 
    real(8), intent(inout)  :: x(n)
    real(8), intent(in)     :: b(n)
    real(8), intent(in)     :: tol

!---local variables

    integer :: i,j,k,g
    integer :: irow, icol
    real(8) :: sum2
    real(8) :: norm
    real(8) :: vol = ONE
    real(8), allocatable :: tmp(:)

!---begin execution

    ! allocate temp
    allocate(tmp(n))
    vol = ONE/dble(n)

    ! start counter
    iter = 1

    ! loop until converged
    do while(iter <= 1000000) 

      ! init norm sum
      sum2 = ZERO 

      ! begin loop over rows
      do irow = 1, n

        ! initialize y
        tmp(irow) = ZERO

        ! loop over columns in that row but skip diagonal
        do j = row(irow), row(irow+1) - 1

          ! continue if this diagonal element
          if (j == diag(irow)) then
            cycle
          end if

          tmp(irow) = tmp(irow) + val(j)*x(col(j))

        end do

        ! subtract RHS value
        tmp(irow) = b(irow) - tmp(irow)

        ! divide by diagonal
        tmp(irow) = tmp(irow)/val(diag(irow))

        ! get region number
!       vol = geometry % fvol_map(ceiling(real(irow)/real(geometry%nfg)))

        ! sum the difference
        sum2 = sum2 + vol*(tmp(irow) - x(irow))**2

      end do

      ! compute point-wise L2 norm 
      norm = sqrt(sum2)

      ! set all temp x to x
      x = tmp 

      ! check convergence
      if (norm < tol) exit 

      ! increase counter
      iter = iter + 1

    end do

    deallocate(tmp)

  end subroutine csr_jacobi 

!===============================================================================
! CSR_GAUSS_SEIDEL
!===============================================================================

  subroutine csr_gauss_seidel(row,col,val,diag,x,b,n,nz,tol,iter)

!---external arguments

    use constants,  only: ZERO, ONE
    use global,     only: geometry

!---arguments

    integer, intent(in)     :: n
    integer, intent(in)     :: nz
    integer, intent(inout)  :: iter
    integer, intent(in)     :: row(n+1)
    integer, intent(in)     :: col(nz)
    integer, intent(in)     :: diag(n)
    real(8), intent(in)     :: val(nz)
    real(8), intent(inout)  :: x(n)
    real(8), intent(in)     :: b(n)
    real(8), intent(in)     :: tol

!---local variables

    integer :: irow, icol
    integer :: i, j, k, g
    integer :: idx
    real(8) :: sum2 
    real(8) :: norm
    real(8) :: vol=ONE
    real(8), allocatable :: tmp(:)
    real(8), allocatable :: tmp1(:)

!---begin execution

    ! allocate temp
    allocate(tmp(n))

    vol = ONE/dble(n)

    ! start counter
    iter = 1

    ! loop until converged
    do while(iter <= 10000000)

      ! set norm sum to zero
      sum2 = ZERO

      ! begin loop over rows
      do irow = 1, n

        ! initialize y
        tmp(irow) = ZERO

        ! loop over columns in that row but skip diagonal
        do icol = row(irow), row(irow+1) - 1

          ! continue if this diagonal element
          if (icol == diag(irow)) then
            cycle
          end if

          tmp(irow) = tmp(irow) + val(icol)*x(col(icol))

        end do

        ! subtract RHS value
        tmp(irow) = b(irow) - tmp(irow)

        ! divide by diagonal
        tmp(irow) = tmp(irow)/val(diag(irow))

        ! get region number
!       idx = ceiling(real(irow)/real(geometry%nfg))
!       vol = geometry % fvol_map(idx)

        ! sum for norm
        sum2 = sum2 + vol*(tmp(irow) - x(irow))**2

        ! set this value in x
        x(irow) = tmp(irow)

      end do

      ! compute point-wise L2 norm 
      norm = sqrt(sum2)

      ! check convergence
      if (norm < tol) exit

      ! increment counter
      iter = iter + 1

    end do

    deallocate(tmp)

  end subroutine csr_gauss_seidel

!===============================================================================
! EXPM_PADE 
!===============================================================================

  subroutine expm_pade(A,N,dt,EXPM) 

!---external references

    use constants,  only: ZERO, ONE

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

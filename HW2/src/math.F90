module math 

!-module options

  implicit none
  private
  public :: sort_csr, csr_matvec_mult, csr_jacobi, csr_gauss_seidel

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
        vol = geometry % fvol_map(ceiling(real(irow)/real(geometry%nfg)))

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

    ! start counter
    iter = 1

    ! loop until converged
    do while(iter <= 1000000)

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
        idx = ceiling(real(irow)/real(geometry%nfg))
        vol = geometry % fvol_map(idx)

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

end module math 

subroutine expm1 ( n, a, e )

!*****************************************************************************80
!
!! EXPM1 is essentially MATLAB's built-in matrix exponential algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    Cleve Moler, Charles Van Loan
!
!  Reference:
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) E(N,N), the estimate for exp(A).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  real ( kind = 8 ) a_norm
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n,n)
  real ( kind = 8 ) e(n,n)
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) k
  logical p
  integer ( kind = 4 ) , parameter :: q = 6
  real ( kind = 8 ) r8_log_2
  real ( kind = 8 ) r8mat_norm_li
  integer ( kind = 4 ) s
  real ( kind = 8 ) x(n,n)

  a2(1:n,1:n) = a(1:n,1:n)

  a_norm = r8mat_norm_li ( n, n, a2 )

  ee = int ( r8_log_2 ( a_norm ) ) + 1

  s = max ( 0, ee + 1 )

  a2(1:n,1:n) = a2(1:n,1:n) / 2.0D+00**s

  x(1:n,1:n) = a2(1:n,1:n)

  c = 0.5D+00

  call r8mat_identity ( n, e )
  e(1:n,1:n) = e(1:n,1:n) + c * a2(1:n,1:n)

  call r8mat_identity ( n, d )
  d(1:n,1:n) = d(1:n,1:n) - c * a2(1:n,1:n)

  p = .true.

  do k = 2, q

    c = c * real ( q - k + 1, kind = 8 ) &
      / real ( k * ( 2 * q - k + 1 ), kind = 8 )

    x(1:n,1:n) = matmul ( a2(1:n,1:n), x(1:n,1:n) )

    e(1:n,1:n) = e(1:n,1:n) + c * x(1:n,1:n)

    if ( p ) then
      d(1:n,1:n) = d(1:n,1:n) + c * x(1:n,1:n)
    else
      d(1:n,1:n) = d(1:n,1:n) - c * x(1:n,1:n)
    end if

    p = .not. p

  end do
!
!  E -> inverse(D) * E
!
  call r8mat_minvm ( n, n, d, e, e )
!
!  E -> E^(2*S)
!
  do k = 1, s
    e(1:n,1:n) = matmul ( e(1:n,1:n), e(1:n,1:n) )
  end do

  return
end
subroutine expm2 ( n, a, e )

!*****************************************************************************80
!
!! EXPM2 uses the Taylor series for the matrix exponential.
!
!  Discussion:
!
!    Formally,
!
!      exp ( A ) = I + A + 1/2 A^2 + 1/3! A^3 + ...
!
!    This function sums the series until a tolerance is satisfied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    Cleve Moler, Charles Van Loan
!
!  Reference:
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) E(N,N), the estimate for exp(A).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) e(n,n)
  real ( kind = 8 ) f(n,n)
  real ( kind = 8 ) g(n,n)
  integer ( kind = 4 ) k
  logical r8mat_insignificant

  e(1:n,1:n) = 0.0D+00

  call r8mat_identity ( n, f )

  k = 1

  do

    if ( r8mat_insignificant ( n, n, e, f ) ) then
      exit
    end if

    e(1:n,1:n) = e(1:n,1:n) + f(1:n,1:n)

    f(1:n,1:n) = matmul ( a(1:n,1:n), f(1:n,1:n) ) / real ( k, kind = 8 )
    k = k + 1

  end do

  return
end
subroutine expm3 ( n, a, e )

!*****************************************************************************80
!
!! EXPM3 approximates the matrix exponential using an eigenvalue approach.
!
!  Discussion:
!
!    exp(A) = V * D * V
!
!    where V is the matrix of eigenvectors of A, and D is the diagonal matrix
!    whose i-th diagonal entry is exp(lambda(i)), for lambda(i) an eigenvalue
!    of A.
!
!    This function is accurate for matrices which are symmetric, orthogonal,
!    or normal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    Cleve Moler, Charles Van Loan
!
!  Reference:
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) E(N,N), the estimate for exp(A).
!
! [ V, D ] = eig ( A );
! E = V * diag ( exp ( diag ( D ) ) ) / V;
  return
end
function r8_log_2 ( x )

!*****************************************************************************80
!
!! R8_LOG_2 returns the logarithm base 2 of an R8.
!
!  Discussion:
!
!    value = Log ( |X| ) / Log ( 2.0 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose base 2 logarithm is desired.
!    X should not be 0.
!
!    Output, real ( kind = 8 ) R8_LOG_2, the logarithm base 2 of the absolute
!    value of X.  It should be true that |X| = 2**R8_LOG_2.
!
  implicit none

  real ( kind = 8 ) r8_log_2
  real ( kind = 8 ) x

  if ( x == 0.0D+00 ) then
    r8_log_2 = - huge ( x )
  else
    r8_log_2 = log ( abs ( x ) ) / log ( 2.0D+00 )
  end if

  return
end
subroutine r8mat_fss ( n, a, nb, b, info )

!*****************************************************************************80
!
!! R8MAT_FSS factors and solves multiple R8GE systems.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    This routine does not save the LU factors of the matrix, and hence cannot
!    be used to efficiently solve multiple linear systems, or even to
!    factor A at one time, and solve a single linear system at a later time.
!
!    This routine uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input/output, real ( kind = 8 ) B(N,NB).
!    On input, the right hand sides of the linear system.
!    On output, the solutions of the linear systems.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  real ( kind = 8 ) piv
  real ( kind = 8 ) row(n)
  real ( kind = 8 ) t(nb)
  real ( kind = 8 ) temp

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a(i,jcol) ) ) then
        piv = abs ( a(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_FSS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a(jcol,1:n)
      a(jcol,1:n) = a(ipiv,1:n)
      a(ipiv,1:n) = row(1:n)

      t(1:nb)      = b(jcol,1:nb)
      b(jcol,1:nb) = b(ipiv,1:nb)
      b(ipiv,1:nb) = t(1:nb)

    end if
!
!  Scale the pivot row.
!
    a(jcol,jcol+1:n) = a(jcol,jcol+1:n) / a(jcol,jcol)
    b(jcol,1:nb) = b(jcol,1:nb) / a(jcol,jcol)
    a(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a(i,jcol) /= 0.0D+00 ) then
        temp = - a(i,jcol)
        a(i,jcol) = 0.0D+00
        a(i,jcol+1:n) = a(i,jcol+1:n) + temp * a(jcol,jcol+1:n)
        b(i,1:nb) = b(i,1:nb) + temp * b(jcol,1:nb)
      end if
    end do

  end do
!
!  Back solve.
!
  do j = 1, nb
    do jcol = n, 2, -1
      b(1:jcol-1,j) = b(1:jcol-1,j) - a(1:jcol-1,jcol) * b(jcol,j)
    end do
  end do

  return
end
subroutine r8mat_identity ( n, a )

!*****************************************************************************80
!
!! R8MAT_IDENTITY stores the identity matrix in an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) A(N,N), the N by N identity matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i

  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 1.0D+00
  end do

  return
end
function r8mat_insignificant ( m, n, r, s )

!*****************************************************************************80
!
!! R8MAT_INSIGNIFICANT determines if an R8MAT is insignificant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the dimension of the matrices.
!
!    Input, real ( kind = 8 ) R(M,N), the vector to be compared against.
!
!    Input, real ( kind = 8 ) S(M,N), the vector to be compared.
!
!    Output, logical R8MAT_INSIGNIFICANT, is TRUE if S is insignificant
!    compared to R.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(m,n)
  logical r8mat_insignificant
  real ( kind = 8 ) s(m,n)
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  logical value

  value = .true.

  do j = 1, n
    do i = 1, m

      t = r(i,j) + s(i,j)
      tol = epsilon ( r(i,j) ) * abs ( r(i,j) )

      if ( tol < abs ( r(i,j) - t ) ) then 
        value = .false.
        exit
      end if

    end do
  end do
  
  r8mat_insignificant = value

  return
end
subroutine r8mat_minvm ( n1, n2, a, b, c )

!*****************************************************************************80
!
!! R8MAT_MINVM computes inverse(A) * B for R8MAT's.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N1,N1), B(N1,N2), the matrices.
!
!    Output, real ( kind = 8 ) C(N1,N2), the result, C = inverse(A) * B.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a(n1,n1)
  real ( kind = 8 ) alu(n1,n1)
  real ( kind = 8 ) b(n1,n2)
  real ( kind = 8 ) c(n1,n2)
  integer ( kind = 4 ) info

  alu(1:n1,1:n1) = a(1:n1,1:n1)
  c(1:n1,1:n2) = b(1:n1,1:n2)

  call r8mat_fss ( n1, alu, n2, c, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_MINVM - Fatal error!'
    write ( *, '(a)' ) '  The matrix A was numerically singular.'
    stop
  end if

  return
end
function r8mat_norm_l1 ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_L1 returns the matrix L1 norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!    The matrix L1 norm is defined as:
!
!      R8MAT_NORM_L1 = max ( 1 <= J <= N )
!        sum ( 1 <= I <= M ) abs ( A(I,J) ).
!
!    The matrix L1 norm is derived from the vector L1 norm, and
!    satisifies:
!
!      r8vec_norm_l1 ( A * x ) <= r8mat_norm_l1 ( A ) * r8vec_norm_l1 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose L1 norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_L1, the L1 norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_norm_l1

  r8mat_norm_l1 = 0.0D+00

  do j = 1, n
    r8mat_norm_l1 = max ( r8mat_norm_l1, sum ( abs ( a(1:m,j) ) ) )
  end do

  return
end
function r8mat_norm_li ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_LI returns the matrix L-oo norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!    The matrix L-oo norm is defined as:
!
!      R8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
!
!    The matrix L-oo norm is derived from the vector L-oo norm,
!    and satisifies:
!
!      r8vec_norm_li ( A * x ) <= r8mat_norm_li ( A ) * r8vec_norm_li ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose L-oo
!    norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_LI, the L-oo norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8mat_norm_li

  r8mat_norm_li = 0.0D+00

  do i = 1, m
    r8mat_norm_li = max ( r8mat_norm_li, sum ( abs ( a(i,1:n) ) ) )
  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

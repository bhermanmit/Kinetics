module power_iter


  use global,        only: ktol, stol, itol
  use loss_operator, only: loss_operator_type,init_M_operator,                 &
 &                         build_loss_matrix,destroy_M_operator
  use prod_operator, only: prod_operator_type,init_F_operator,                 &
 &                         build_prod_matrix,destroy_F_operator

  implicit none
  private
  public :: power_execute 

  type(loss_operator_type) :: loss
  type(prod_operator_type) :: prod

  integer     :: ierr       ! error flag
  logical     :: iconv      ! did the problem converged
  real(8)     :: k_n        ! new k-eigenvalue
  real(8)     :: k_o        ! old k-eigenvalue
  real(8), allocatable :: phi(:)   ! flux vector
  real(8), allocatable :: S_n(:)   ! new source vector
  real(8), allocatable :: S_o(:)   ! old source vector

contains

!===============================================================================
! CMFD_POWER_EXECUTE
!===============================================================================

  subroutine power_execute(inner_solver)

!---external references

    use global,  only: time_build, time_power
    use timing,  only: timer_start, timer_stop

!---arguments

    external :: inner_solver

    ! initialize matrices and vectors
    call init_data()

    ! start build timer
    call timer_start(time_build)

    ! set up M loss matrix
    call build_loss_matrix(loss) 

    ! set up F production matrix
    call build_prod_matrix(prod)

    ! stop build timer and start power iter timer
    call timer_stop(time_build)
    call timer_start(time_power)

    ! begin power iteration 
    call execute_power_iter(inner_solver)

    ! stop power iter timer
    call timer_stop(time_power)

    ! extract results
    call extract_results()

    ! deallocate petsc objects
    call finalize()

  end subroutine power_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

!---external references

    use constants,  only: ONE

!---local variables

    integer :: n      ! problem size
    integer :: i      ! counter
    real(8) :: guess  ! initial guess

!---begin execution

    ! set up matrices
    call init_M_operator(loss)
    call init_F_operator(prod)

    ! get problem size
    n = loss%n

    ! set up flux vector
    allocate(phi(n))

    ! set up source vectors
    allocate(S_n(n))
    allocate(S_o(n))

    ! set initial guess
    k_n = ONE 
    k_o = ONE
    S_n = ONE
    S_o = ONE 

    ! put random guess to excite all harmonics
    do i=1,n
      phi(i) = rand()
    end do

  end subroutine init_data

!===============================================================================
! EXECUTE_POWER_ITER  in the main power iteration routine 
!                     for the cmfd calculation
!===============================================================================

  subroutine execute_power_iter(inner_solver)

!---external references

    use cmfd_header,  only: calc_power
    use global,       only: cmfd, geometry, time_inner
    use math,         only: csr_matvec_mult, csr_jacobi
    use timing,       only: timer_start, timer_stop


!---local variables

    real(8)     :: num       ! numerator for eigenvalue update
    real(8)     :: den       ! denominator for eigenvalue update
    real(8)     :: norm      ! norm
    external :: inner_solver
    integer     :: i         ! iteration counter
    integer     :: n
    integer     :: nz 

!--begin execution

    ! set sizes
    n = loss % n
    nz = size(loss % col)

    ! reset convergence flag
    iconv = .FALSE.

    ! compute source vector
    S_o =  csr_matvec_mult(prod%row_csr,prod%col,prod%val,phi,prod%n)

    ! compute initail nodal power
    call calc_power(cmfd,S_o,n,geometry)

    ! move power to old
    cmfd % power_o = cmfd % power_n

    ! begin power iteration
    do i = 1,10000

      ! compute source vector
      S_o =  csr_matvec_mult(prod%row_csr,prod%col,prod%val,phi,prod%n)

      ! normalize source vector
      S_o = S_o/k_o

      ! compute new flux vector
      call timer_start(time_inner)
      call inner_solver(loss % row_csr, loss % col, loss % val, loss % diag, phi, S_o, n, nz, itol)
      call timer_stop(time_inner)

      ! compute new source vector
      S_n = csr_matvec_mult(prod%row_csr,prod%col,prod%val,phi,prod%n)

      ! compute new power
      call calc_power(cmfd,S_n,n,geometry)

      ! compute new k-eigenvalue
      num = sum(S_n)
      den = sum(S_o) 
      k_n = num/den

      ! renormalize the old source
      S_o = S_o * k_o 

      ! check convergence
      call convergence(cmfd, i, norm)

      ! to break or not to break
      if (iconv) then
        cmfd % phi = phi
        cmfd % iter = i
        cmfd % keff = k_n
        cmfd % norm = norm
        cmfd % dr = cmfd % norm/cmfd % norm_o
        exit
      end if

      ! record old values
      k_o = k_n
      cmfd % power_o = cmfd % power_n
      cmfd % norm_o = norm

    end do

  end subroutine execute_power_iter 

!===============================================================================
! CONVERGENCE checks the convergence of eigenvalue, eigenvector and source
!===============================================================================

  subroutine convergence(cmfd,iter, norm)

!---external references

    use cmfd_header,  only: cmfd_type

!---arguments

    integer :: iter
    real(8) :: norm
    type(cmfd_type) :: cmfd

!---local variables

    real(8)     :: kerr           ! error in keff

!---begin execution

    ! reset convergence flag
    iconv = .FALSE.

    ! calculate error in keff
    kerr = abs(k_o - k_n)/k_n

    ! calculate max error in source
    norm = sqrt(sum((cmfd % power_n - cmfd % power_o)**2)) 

    ! check for convergence
    if(kerr < ktol .and. norm < stol) iconv = .TRUE.

    ! print out to user (TODO: make formatted)
    write(*,100) iter,k_n,norm 

 100 format(I5,5X,"EIG: ",F7.5,5X,"NORM: ",1PE9.3)

  end subroutine convergence

!==============================================================================
! EXTRACT_RESULTS
!==============================================================================

  subroutine extract_results()


  end subroutine extract_results

!==============================================================================
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize data objects
    call destroy_M_operator(loss)
    call destroy_F_operator(prod)

    deallocate(phi)
    deallocate(S_n)
    deallocate(S_o)

  end subroutine finalize

end module power_iter

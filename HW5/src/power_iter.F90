module power_iter

!-module references

  use global,        only: ktol, stol, itol, loss, prod
  use loss_operator
  use prod_operator

!-module options

  implicit none
  private
  public :: power_execute 

!-module external references

# include "finclude/petsc.h90"

!-module variables

  KSP :: ksp
  Vec :: phip
  Vec :: b
  PC  :: pc

!-module variables

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

  subroutine power_execute()

!---external references

    use global,  only: time_power, geometry
    use timing,  only: timer_start, timer_stop

!---begin execution

    ! initialize matrices and vectors
    call init_data()

    ! initialize solver
    call init_solver()

    ! start power iteration timer
    call timer_start(time_power)

    ! begin power iteration 
    call execute_power_iter()

    ! stop power iter timer
    call timer_stop(time_power)

    ! deallocate petsc objects
    call finalize()

  end subroutine power_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

!---external references

    use constants,  only: ONE
    use global,     only: guess, time_build, mpi_err
    use timing,  only: timer_start, timer_stop

!---local variables

    integer :: n      ! problem size
    integer :: i      ! counter

!---begin execution

    ! set up matrices
    call init_M_operator(loss)
    call init_F_operator(prod)

    ! get problem size
    n = loss%n

    ! start build timer
    call timer_start(time_build)

    ! set up M loss matrix
    call build_loss_matrix(loss) 

    ! set up F production matrix
    call build_prod_matrix(prod)

    ! stop build timer
    call timer_stop(time_build)

    ! set up matrices
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,loss%n,loss%n,loss%row_csr,&
                                   loss%col,loss%val,loss%oper,mpi_err)
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,prod%n,prod%n,prod%row_csr,&
                                   prod%col,prod%val,prod%oper,mpi_err)

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
    if (trim(guess) == 'rand') then
      do i=1,n
        phi(i) = rand()
      end do
    else
      phi = ONE
    end if

  end subroutine init_data

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver()

!---external references

    use global,  only: itol, mpi_err

!---begin execution

    ! set up krylov solver
    call KSPCreate(PETSC_COMM_SELF,ksp,mpi_err)
    call KSPSetTolerances(ksp,itol,PETSC_DEFAULT_DOUBLE_PRECISION,     &
   &                      PETSC_DEFAULT_DOUBLE_PRECISION,                      &
   &                      PETSC_DEFAULT_INTEGER,mpi_err)
    call KSPSetType(ksp,KSPGMRES,mpi_err)
    call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,mpi_err)
    call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,mpi_err)
    call KSPGetPC(ksp,pc,mpi_err)
    call PCSetType(pc,PCILU,mpi_err)
    call PCFactorSetLevels(pc,5,mpi_err)
    call KSPSetFromOptions(ksp,mpi_err)

  end subroutine init_solver

!===============================================================================
! EXECUTE_POWER_ITER  in the main power iteration routine 
!                     for the cmfd calculation
!===============================================================================

  subroutine execute_power_iter()

!---external references

    use cmfd_header,  only: calc_power
    use error,        only: fatal_error
    use global,       only: cmfd, geometry, time_inner, message, mpi_err
    use math,         only: csr_matvec_mult, csr_jacobi
    use timing,       only: timer_start, timer_stop


!---local variables

    real(8)     :: num       ! numerator for eigenvalue update
    real(8)     :: den       ! denominator for eigenvalue update
    real(8)     :: norm      ! norm
    integer     :: i         ! iteration counter
    integer     :: n
    integer     :: nz 
    integer     :: inner

!--begin execution

    ! set sizes
    n = loss % n
    nz = size(loss % col)

    ! reset convergence flag
    iconv = .FALSE.

    ! associate petsc vectors
    call VecCreateSeqWithArray(PETSC_COMM_WORLD,1,loss%n,phi,phip,mpi_err)
    call VecCreateSeqWithArray(PETSC_COMM_WORLD,1,loss%n,S_o,b,mpi_err)

    ! set up krylov info
    call KSPSetOperators(ksp, loss%oper, loss%oper, SAME_NONZERO_PATTERN, mpi_err)

    call KSPSetUp(ksp,mpi_err)

    ! calculate preconditioner (ILU)
    call PCFactorGetMatrix(pc,loss%oper,ierr)

    ! compute source vector
    S_o =  csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val,phi,prod%n)

    ! compute initail nodal power
    call calc_power(cmfd,S_o,n,geometry)

    ! move power to old
    cmfd % power_o = cmfd % power_n

    ! begin power iteration
    do i = 1,1000000

      ! compute source vector
      S_o =  csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val,phi,prod%n)

      ! normalize source vector
      S_o = S_o/k_o

      ! compute new flux vector
      call timer_start(time_inner)
      call KSPSolve(ksp,b,phip,mpi_err)
      call timer_stop(time_inner)

      ! compute new source vector
      S_n = csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val,phi,prod%n)

      ! compute new power
      call calc_power(cmfd,S_n,n,geometry)

      ! compute new k-eigenvalue
      num = sum(S_n)
      den = sum(S_o) 
      k_n = num/den

      ! renormalize the old source
      S_o = S_o * k_o 

      ! check convergence
      call convergence(cmfd, i, norm, inner)

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

  subroutine convergence(cmfd,iter,norm,inner)

!---external references

    use cmfd_header,  only: cmfd_type
    use global, only: time_inner

!---arguments

    integer :: iter
    integer :: inner
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

 100 format(I0,5X,"EIG: ",F7.5,5X,"NORM: ",1PE9.3,5X)

  end subroutine convergence

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

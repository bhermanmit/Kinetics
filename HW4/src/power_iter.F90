module power_iter

!-module references

  use global,        only: ktol, stol, itol, loss, prod, cmfd, adjoint, loss_adj, prod_adj
  use loss_operator, only: init_M_operator,                 &
 &                         build_loss_matrix,destroy_M_operator
  use prod_operator, only: init_F_operator,                 &
 &                         build_prod_matrix,destroy_F_operator
  use operator_header

!-module options

  implicit none
  private
  public :: power_execute 

!-module variables

  integer     :: ierr       ! error flag
  logical     :: iconv      ! did the problem converged
  real(8)     :: k_n        ! new k-eigenvalue
  real(8)     :: k_o        ! old k-eigenvalue
  real(8), allocatable :: phi(:)   ! flux vector
  real(8), allocatable :: S_n(:)   ! new source vector
  real(8), allocatable :: S_o(:)   ! old source vector
  type(operator_type), pointer :: left => null()
  type(operator_type), pointer :: right => null()
  real(8), pointer :: flux(:) => null()

contains

!===============================================================================
! CMFD_POWER_EXECUTE
!===============================================================================

  subroutine power_execute(inner_solver,adjoint)

!---external references

    use global,  only: time_build, time_power, geometry
    use timing,  only: timer_start, timer_stop

!---arguments

    external :: inner_solver
    character(len=*) :: adjoint

    ! set pointers
    if (adjoint == 'math' .or. adjoint == 'physical') then
      left => loss_adj
      right => prod_adj
      flux => cmfd % phi_adj
    else
      left => loss
      right => prod
      flux => cmfd % phi
    end if

    ! initialize matrices and vectors
    call init_data()

    ! start build timer
    call timer_start(time_build)

    ! set up M loss matrix
    call build_loss_matrix(left,adjoint)

    ! set up F production matrix
    call build_prod_matrix(right,adjoint)

    ! stop build timer and start power iter timer
    call timer_stop(time_build)
    call timer_start(time_power)

    ! begin power iteration 
    call execute_power_iter(inner_solver)

    ! stop power iter timer
    call timer_stop(time_power)

    ! deallocate petsc objects
!   call finalize()

    ! compute final fission source
    if (.not.(adjoint == 'math' .or. adjoint == 'physical')) call fission_src()

  end subroutine power_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

!---external references

    use constants,  only: ONE
    use global,     only: guess

!---local variables

    integer :: n      ! problem size
    integer :: i      ! counter

!---begin execution

    ! set up matrices
    call init_M_operator(left)
    call init_F_operator(right)

    ! get problem size
    n = left%n

    ! set up flux vector
    if(.not.allocated(phi)) allocate(phi(n))

    ! set up source vectors
    if(.not.allocated(S_n)) allocate(S_n(n))
    if(.not.allocated(S_o)) allocate(S_o(n))

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
! EXECUTE_POWER_ITER  in the main power iteration routine 
!                     for the cmfd calculation
!===============================================================================

  subroutine execute_power_iter(inner_solver)

!---external references

    use cmfd_header,  only: calc_power
    use error,        only: fatal_error
    use global,       only: geometry, time_inner, message
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
    integer     :: inner

!--begin execution

    ! set sizes
    n = left % n
    nz = size(left % col)

    ! reset convergence flag
    iconv = .FALSE.

    ! compute source vector
    S_o =  csr_matvec_mult(right%row_csr+1,right%col+1,right%val,phi,right%n)

    ! compute initail nodal power
    call calc_power(cmfd,S_o,n,geometry)

    ! move power to old
    cmfd % power_o = cmfd % power_n

    ! begin power iteration
    do i = 1,1000000

      ! compute source vector
      S_o =  csr_matvec_mult(right%row_csr+1,right%col+1,right%val,phi,right%n)

      ! normalize source vector
      S_o = S_o/k_o

      ! compute new flux vector
      call timer_start(time_inner)
      call inner_solver(left % row_csr+1, left % col+1, left % val, left % diag, phi, S_o, n, nz, itol,inner)
      if (inner >= 1000000) then
        message = 'Inner max iteration met'
        call fatal_error()
      end if
      call timer_stop(time_inner)

      ! compute new source vector
      S_n = csr_matvec_mult(right%row_csr+1,right%col+1,right%val,phi,right%n)

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
        flux = phi
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
    write(*,100) iter,k_n,norm,inner,time_inner%elapsed

 100 format(I0,5X,"EIG: ",F7.5,5X,"NORM: ",1PE9.3,5X,"INNER:",I0,2X,1PE9.3)

  end subroutine convergence

!===============================================================================
! CONVERGENCE checks the convergence of eigenvalue, eigenvector and source
!===============================================================================

  subroutine fission_src() 

!---external references

    use constants, only: ZERO, ONE
    use global,  only: cmfd, geometry, material
    use material_header,  only: material_type
    use math,  only: csr_matvec_mult

!---local variables

    integer :: irow, g, idx
    real(8) :: vol, pow
    type(material_type), pointer :: m

!---begin execution

    ! zero out fission src
    cmfd % fsrc = ZERO

    ! normalize flux
!   cmfd % phi = cmfd % phi / sum(cmfd % phi)

    ! begin loop rows
    do irow = 1, right % n

      ! get material
      idx = ceiling(real(irow)/real(geometry % nfg)) 
      m => material(geometry % fmat_map(idx))

      ! get group number at this row
      g = mod(irow-1,geometry%nfg) + 1

      ! get volume
      vol = geometry % fdx_map(idx)*geometry % fdy_map(idx)*geometry % fdz_map(idx)

      ! bank fission rate
      cmfd % fsrc(g) = cmfd % fsrc(g) + m % fissvec(g)*cmfd%phi(irow)*vol

    end do

    ! normalize so that power is one
    pow = sum(csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val/cmfd%keff,        &
              cmfd%phi,prod%n))
    cmfd % phi = cmfd % phi * ONE / pow

end subroutine fission_src

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

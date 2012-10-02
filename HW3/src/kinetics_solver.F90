module kinetics_solver 

!-module references

  use global,  only: mpi_err

!-module options

  implicit none
  private
  public :: kinetics_execute 

!-module external references

# include "finclude/petsc.h90"

!-module variables

  KSP :: ksp
  Vec :: phi
  Vec :: b
  PC  :: pc

contains

!===============================================================================
! KINETICS_EXECUTE
!===============================================================================

  subroutine kinetics_execute(inner_solver)

!---external references

    use global,  only: kine

!---arguments

    external :: inner_solver

!---begin execution

    ! initialize the time 0 data
    call init_data()

    ! initialize solver
    call init_solver()

    ! run through the kinetics iterations
    call execute_kinetics_iter(inner_solver)
 
  end subroutine kinetics_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

!---external references

    use constants,          only: ONE
    use global,             only: cmfd, geometry, material, kine, prod, nt
    use kinetics_operator,  only: init_K_operator
    use math,               only: csr_matvec_mult

!---local variables

    real(8) :: pow

!---begin execution

    ! normalize initial power to unity and set initial power
    pow = sum(csr_matvec_mult(prod%row_csr,prod%col,prod%val/cmfd%keff,        &
              cmfd%phi,prod%n))
    cmfd % phi = cmfd % phi * ONE / pow

    ! compute steady state precursors
    call compute_initial_precursors()

    ! set up matrices
    call init_K_operator(kine)

    ! allocate output arrays
    if (.not.allocated(cmfd%core_power)) allocate(cmfd % core_power(nt+1))
    if (.not.allocated(cmfd%time)) allocate(cmfd % time(nt+1))

  end subroutine init_data

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver()

!---external references

    use global,  only: itol

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
! COMPUTE_INITIAL_PRECURSORS
!===============================================================================

  subroutine compute_initial_precursors()

!---external references

    use constants,        only: NUM_PRECS, beta, lambda
    use global,           only: geometry, material, cmfd
    use material_header,  only: material_type

!---local variables

    integer :: i
    integer :: j
    real(8) :: fiss
    type(material_type), pointer :: m => null()

!---begin execution

    ! allocate precurors
    if (.not.allocated(cmfd % C)) allocate(cmfd%C(geometry%nfx*geometry%nfy *  &
                                                  geometry%nfz,NUM_PRECS))

    ! begin loop around space
    do j = 1,size(cmfd % C,1)

      ! compute fission reaction rate over all groups
      m => material(geometry % fmat_map(j))
      fiss = sum(m % fissvec * cmfd % phi(geometry % nfg * (j - 1) + 1: &
                                                   j * geometry % nfg))
                                                   
      ! begin loop around precursor groups
      do i = 1, NUM_PRECS

        cmfd % C(j,i) = beta(i)/(lambda(i) * cmfd % keff) * fiss 

      end do

    end do

  end subroutine compute_initial_precursors

!===============================================================================
! EXECUTE_KINETICS_ITER  in the main kinetics iteration routine 
!                         for the cmfd calculation
!===============================================================================

  subroutine execute_kinetics_iter(inner_solver)

!---external references

    use constants,          only: ZERO, ONE
    use error,              only: fatal_error
    use global,             only: nt, dt, kine, cmfd, geometry, material,      &
                                  itol, prod, message
    use kinetics_operator,  only: build_kinetics_matrix
    use math,               only: csr_matvec_mult

!---arguments

    external :: inner_solver

!---local variables

    integer :: i
    integer :: iters
    integer :: n
    integer :: nz
    real(8) :: curr_time
    real(8) :: pow
    real(8), allocatable :: rhs(:)

!---begin execution

    ! allocate RHS
    n = size(cmfd % phi)
    nz = size(kine % col)
    allocate(rhs(n))

    ! set initial output
    cmfd % time(1) = ZERO
    cmfd % core_power(1) = ONE

    ! associate petsc vectors
    call VecCreateSeqWithArray(PETSC_COMM_WORLD,1,kine%n,cmfd%phi,phi,mpi_err)
    call VecCreateSeqWithArray(PETSC_COMM_WORLD,1,kine%n,rhs,b,mpi_err)

    ! begin loop around time
    do i = 1, nt

      ! compute current time
      curr_time = dble(i)*dt
    
      ! change the material via kinetics mods
      call change_data(curr_time)

      ! build matrices
      call build_kinetics_matrix(kine)

      ! build right hand side
      call build_rhs(rhs,n)

      ! set up krylov info
      call KSPSetOperators(ksp, kine%oper, kine%oper, SAME_NONZERO_PATTERN, mpi_err)
      call KSPSetUp(ksp,mpi_err)

      ! calculate preconditioner (ILU)
      call PCFactorGetMatrix(pc,kine%oper,mpi_err)

      ! solve matrices
      call KSPSolve(ksp,b,phi,mpi_err)     
!     call inner_solver(kine % row_csr+1, kine % col+1, kine % val, kine % diag,    &
!                       cmfd % phi, rhs, n, nz, itol,iters)
!     if (iters >= 10000000) then 
!       message = "Inner Iteration limit exceed during transient!"
!       call fatal_error()
!     end if
 
     ! compute end of time step precursor concentration
     call compute_final_precursors()

     ! compute power
     pow = sum(csr_matvec_mult(prod%row_csr,prod%col,prod%val/cmfd%keff,       &
              cmfd%phi,prod%n))
     write(*,*) 'Step:', i,' /',nt,' POWER:',pow
     cmfd % time(i+1) = curr_time
     cmfd % core_power(i+1) = pow
   end do

   ! deallocate RHS
   deallocate(rhs)

  end subroutine execute_kinetics_iter 

!===============================================================================
! COMPUTE_FINAL_PRECURSORS
!===============================================================================

  subroutine compute_final_precursors()

!---external references

    use constants,        only: NUM_PRECS, beta, lambda, ONE
    use global,           only: geometry, material, cmfd, dt
    use material_header,  only: material_type

!---local variables

    integer :: i
    integer :: j
    real(8) :: fiss
    type(material_type), pointer :: m => null()

!---begin execution

    ! begin loop around space
    do j = 1,size(cmfd%C,1)

      ! compute fission reaction rate over all groups
      m => material(geometry % fmat_map(j))
      fiss = sum(m % fissvec * cmfd % phi(geometry % nfg * (j - 1) + 1: &
                                                   j * geometry % nfg))

      ! begin loop around precursor groups
      do i = 1, NUM_PRECS

        cmfd % C(j,i) = beta(i)*dt/((ONE + dt*lambda(i)) * cmfd % keff) * fiss +  &
                        cmfd % C(j,i)/(ONE + dt*lambda(i))

      end do

    end do

  end subroutine compute_final_precursors

!===============================================================================
! CHANGE_DATA
!===============================================================================

  subroutine change_data(t)

!---external references

    use constants,        only: ONE, vel, beta, lambda
    use error,            only: fatal_error
    use global,           only: material, kinetics, n_kins, message, dt,       &
                                n_materials, cmfd
    use kinetics_header,  only: kinetics_type
    use material_header,  only: material_type

!---arguments

    real(8) :: t

!---local variables

    integer :: i
    real(8) :: val
    type(kinetics_type), pointer :: k => null()
    type(material_type), pointer :: m => null()

!---begin execution

    ! begin loop around kinetics mods
    do i = 1, n_kins

      ! point to object
      k => kinetics(i)

      ! point to material object
      m => material(k % mat_id)

      ! check to shift index
      if (t > k % time(k % idxt+1)) k % idxt = k % idxt + 1

      ! interpolate value
      val = k % val(k%idxt) + (k % val(k%idxt+1) - k % val(k%idxt)) /          &
            (k % time(k%idxt+1) - k % time(k%idxt)) * (t - k % time(k%idxt))

      ! begin case structure to replace value
      select case (trim(k % xs_id))

        case ('absxs')

          ! remove scattering component from total xs
          m % totalxs(k % g) = m % totalxs(k % g) - sum(m % scattxs(:,k % g))

          ! change absorption
          m % absorxs(k % g) = val
          m % totalxs(k % g) = val

          ! re-add back in scattering
          m % totalxs(k % g) = m % totalxs(k % g) + sum(m % scattxs(:,k % g))

        case DEFAULT

          message = 'Kineics modification not supported!'
          call fatal_error()

      end select

    end do

    ! loop through materials and compute kinetics factors
    do i = 1, n_materials

      ! point to materal
      m => material(i)

      ! compute additive factor that goes on to the removal xs
      m % kinrem = ONE/(vel * dt)

      ! compute multiplicative factor that goes on fission xs
      m % kinfis = (ONE - sum(beta)) * m % chip / cmfd % keff +                &
      m % chid * sum((lambda * beta * dt)/((ONE + lambda*dt) * cmfd % keff))

    end do

  end subroutine change_data

!==============================================================================
! BUILD_RHS
!==============================================================================

  subroutine build_rhs(rhs,n)

!---external references

    use constants,        only: vel, lambda, ONE
    use global,           only: dt, material, geometry, cmfd
    use material_header,  only: material_type

!---arguments

    integer :: n
    real(8) :: rhs(n)

!---local variables

    integer :: irow
    integer :: idx
    integer :: g
    type(material_type), pointer :: m => null()

!---begin execution

    ! begin loop around rows
    do irow = 1, n

      ! get index
      idx = ceiling(real(irow)/real(geometry%nfg))
      g = geometry % nfg - mod(irow,geometry%nfg)

      ! set material pointer
      m => material(geometry % fmat_map(idx))

      ! compute the precursor
      rhs(irow) = cmfd % phi(irow)/(vel(g)*dt) + sum((m%chid(g)*lambda*        &
                                            cmfd % C(idx,:))/(ONE + lambda*dt))
    end do

  end subroutine build_rhs

!==============================================================================
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize data objects

  end subroutine finalize

end module kinetics_solver 

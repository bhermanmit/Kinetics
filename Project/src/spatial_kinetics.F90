module spatial_kinetics 

!-references

  use global,  only: mpi_err

!-module options

  implicit none
  private
  public :: run_spkinetics

!-module external references

# include "finclude/petsc.h90"

contains

!===============================================================================
! RUN_PKINETICS
!===============================================================================

  subroutine run_spkinetics()

!---external references

    use constants,       only: NUM_PRECS
    use global,          only: pke, time, solver_type, geometry
    use implicit_euler,  only: execute_ie1
    use output,          only: header
    use runge_kutta,     only: execute_rk4

!---local variables

    integer :: n  ! size of vectors
    integer :: i  ! loop counter
    real(8) :: dt ! local time step

    Mat :: dfdy
    Mat :: A
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---begin execution

    ! print header for run
    call header("Spatial KINETICS SIMULATION", level=1)

    ! compute size
    n = geometry % nf + NUM_PRECS * geometry % nf

    ! initialize data objects
    call init_data(dfdy, dfdt, dydt, y)

    ! set initial values in y
    call set_init(y)

    ! perform 4th order kaps rentrop
    select case(trim(solver_type))

      case('rk4')
        call execute_rk4(y, dfdy, dfdt, dydt, spk_derivs, spk_jacobn, n, spk_post_timestep)

!     case('ie1')
!       call execute_ie1(y, pk_coefmat, NUM_PRECS+1)

      case DEFAULT

    end select

    call destroy_objects(dfdy, dfdt, dydt, y)

  end subroutine run_spkinetics

!===============================================================================
! INIT_DATA
!===============================================================================

  subroutine init_data(dfdy,dfdt,dydt,y)

!---references

    use constants,      only: NUM_PRECS, ZERO
    use global,         only: loss, prod, geometry
    use loss_operator,  only: init_M_operator, build_loss_matrix
    use prod_operator,  only: init_F_operator, build_prod_matrix

!---arguments

    Mat :: dfdy
    Mat :: A
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---local variables

    integer :: n
    integer :: ng
    integer, allocatable :: d_nnz(:)
    integer, allocatable :: o_nnz(:)

!---begin execution

    ! set up matrices
    call init_M_operator(loss)
    call init_F_operator(prod)

    ! set up M loss matrix
    call build_loss_matrix(loss)

    ! set up F production matrix
    call build_prod_matrix(prod)

    ! set up matrices
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,loss%n,loss%n,loss%row_csr,&
                                   loss%col,loss%val,loss%oper,mpi_err)
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,prod%n,prod%n,prod%row_csr,&
                                   prod%col,prod%val,prod%oper,mpi_err)
    
    ! size of data
    n = geometry % nf 
    ng = geometry % nfg

    ! allocate data sizes
    allocate(d_nnz(n + NUM_PRECS*n))
    allocate(o_nnz(n + NUM_PRECS*n))

    ! create dfdy matrix
    d_nnz(1:n) = loss % d_nnz + NUM_PRECS
    d_nnz(n+1:n+NUM_PRECS*n) = ng + 1
    o_nnz = 0
    call MatCreateAIJ(PETSC_COMM_WORLD,n+NUM_PRECS*n, n+NUM_PRECS*n,&
                      PETSC_DETERMINE, PETSC_DETERMINE, PETSC_NULL, d_nnz,&
                      PETSC_NULL, o_nnz, dfdy, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create dfdt vector
    call VecCreateMPI(PETSC_COMM_WORLD, n+NUM_PRECS*n, PETSC_DETERMINE, dfdt,&
                      mpi_err)
    call VecSet(dfdt, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create dydt vector
    call VecCreateMPI(PETSC_COMM_WORLD, n+NUM_PRECS*n, PETSC_DETERMINE, dydt,&
                      mpi_err)
    call VecSet(dydt, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create y solution vector
    call VecCreateMPI(PETSC_COMM_WORLD, n+NUM_PRECS*n, PETSC_DETERMINE, y,&
                      mpi_err)
    call VecSet(y, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! deallocate data
    deallocate(d_nnz)
    deallocate(o_nnz)

  end subroutine init_data

!===============================================================================
! SETUP_COEFMAT
!===============================================================================

  subroutine spk_jacobn(t,y,dfdy,dfdt,n)

!---external references

    use constants,      only: beta, NUM_PRECS, lambda, pnl, vel, ONE, ZERO
    use global,         only: geometry, cmfd, loss, prod
    use loss_operator,  only: build_loss_matrix
    use prod_operator,  only: build_prod_matrix

!---arguments

    integer :: n
    real(8) :: t
    Mat     :: dfdy
    Vec     :: dfdt
    Vec     :: y

!---local variables

    integer :: nf
    integer :: ng
    integer :: ncols
    integer :: i    ! loop counter
    integer :: irow ! row counter
    integer, allocatable :: cols(:)
    real(8) :: val  ! temp value for matrix setting
    real(8), allocatable :: vals(:)

    real(8), pointer :: yptr(:)
    real(8), pointer :: dfdtptr(:)

    PetscViewer :: viewer

!---begin execution

    ! extract objects
    call VecGetArrayF90(y, yptr, mpi_err)
    call VecGetArrayF90(dfdt, dfdtptr, mpi_err)
    dfdtptr = ZERO

    ! get sub size
    nf = geometry % nf
    ng = geometry % nfg

    ! change the data
    call change_data(t)

    ! rebuild matrices
    call build_loss_matrix(loss)
    call build_prod_matrix(prod)

    ! allocate cols and  initialize to zero
    if (.not. allocated(cols)) allocate(cols(&
         maxval(loss%d_nnz + loss%o_nnz)))
    if (.not. allocated(vals)) allocate(vals(&
         maxval(loss%d_nnz + loss%o_nnz)))
    cols = 0
    vals = ZERO

    ! in flux equation compute effective production operator
    call MatAXPY(loss%oper, -(ONE-sum(beta))/cmfd%keff, prod%oper, SUBSET_NONZERO_PATTERN, mpi_err)
    call MatScale(loss%oper, -ONE, mpi_err)

    ! begin loop around rows
    do irow = 1, nf

      ! set dfdt for flux equation since total xs changes
      dfdtptr(irow) = -yptr(irow)

      ! get row of matrix M
      call MatGetRow(loss%oper, irow-1, ncols, cols, vals, mpi_err)

      ! multiply values by group velocity
      vals = vel(mod(irow-1,ng)+1)*vals

      ! put values in jacobian
      call MatSetValues(dfdy, 1, irow-1, ncols, cols(1:ncols), vals, &
                        INSERT_VALUES, mpi_err)

      ! put the row back
      call MatRestoreRow(loss%oper, irow-1, ncols, cols, vals, mpi_err)

      ! begin loop around precursors
      do i = 1, NUM_PRECS

        ! flux by precursors
        call MatSetValue(dfdy, irow-1, i*nf + (irow - 1), lambda(i), INSERT_VALUES, mpi_err)

        ! precursors by flux
        call MatGetRow(prod%oper, irow-1, ncols, cols, vals, mpi_err)
        vals = vals*beta(i)/cmfd%keff
        call MatSetValues(dfdy, 1, i*nf + (irow-1), ncols, cols(1:ncols), vals, &
                          INSERT_VALUES, mpi_err)
        call MatRestoreRow(prod%oper, irow-1, ncols, cols, vals, mpi_err)

        ! precursors by precursors
        call MatSetValue(dfdy, i*nf + (irow - 1), i*nf + (irow - 1), -lambda(i), INSERT_VALUES, mpi_err)

      end do

    end do

    ! assemble matrix
    call MatAssemblyBegin(dfdy, MAT_FINAL_ASSEMBLY, mpi_err)
    call MatAssemblyEnd(dfdy, MAT_FINAL_ASSEMBLY, mpi_err)

    ! restore vectors
    call VecRestoreArrayF90(y, yptr, mpi_err)
    call VecRestoreArrayF90(dfdt, dfdtptr, mpi_err)

    ! free memory
    deallocate(cols)
    deallocate(vals)

    ! print out matrix
!   call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'jacobian.bin', &
!          FILE_MODE_WRITE, viewer, mpi_err)
!   call MatView(dfdy, viewer, mpi_err)
!   call PetscViewerDestroy(viewer, mpi_err)
stop
  end subroutine spk_jacobn 

!===============================================================================
! PK_DERIVS
!===============================================================================

  subroutine spk_derivs(t,y,dydt,n)

!---external references

    use constants,      only: beta, NUM_PRECS, lambda, ONE, vel
    use global,         only: prod, loss, geometry, cmfd 
    use loss_operator,  only: build_loss_matrix
    use math,           only: csr_matvec_mult
    use prod_operator,  only: build_prod_matrix

!---arguments

    integer :: n
    real(8) :: t
    Mat     :: dfdy
    Vec     :: dydt
    Vec     :: y

!---local variables

    integer :: nf
    integer :: ng
    integer :: i   ! loop counter
    real(8), pointer :: yptr(:)
    real(8), pointer :: dydtptr(:)

!---begin execution

    ! get sub size
    nf = geometry % nf
    ng = geometry % nfg

    ! change the data
    call change_data(t)

    ! rebuild matrices
    call build_loss_matrix(loss)
    call build_prod_matrix(prod)

    ! get pointer to solution
    call VecGetArrayF90(y, yptr, mpi_err)
    call VecGetArrayF90(dydt, dydtptr, mpi_err)

    ! in flux equation compute effective production operator
    call MatAXPY(loss%oper, -(ONE-sum(beta))/cmfd%keff, prod%oper, SUBSET_NONZERO_PATTERN, mpi_err)
    call MatScale(loss%oper, -ONE, mpi_err)

    ! muliply by flux and store
    dydtptr(1:nf) = csr_matvec_mult(loss%row_csr+1,loss%col+1,loss%val,yptr(1:nf),prod%n)

    ! multiply by velocity
    do i = 1,nf
      dydtptr(i) = dydtptr(i)*vel(mod(i-1,ng)+1)
    end do

    ! loop through precursors
    do i = 1,NUM_PRECS
 
      ! from flux equation
      dydtptr(1:nf) = dydtptr(1:nf) + lambda(i)*yptr((i-1)*n + n + 1:(i-1)*n + 2*n)

      ! from precursor equation
      dydtptr((i-1)*n + n + 1:(i-1)*n + 2*n) = beta(i)/cmfd%keff*csr_matvec_mult(loss%row_csr+1, &
                                               loss%col+1,loss%val,yptr(1:nf),prod%n) - lambda(i)*&
                                               yptr((i-1)*n + n + 1:(i-1)*n + 2*n)

    end do

    ! put the pointer back
    call VecRestoreArrayF90(y, yptr, mpi_err)
    call VecRestoreArrayF90(dydt, dydtptr, mpi_err)

  end subroutine spk_derivs

!===============================================================================
! COEFMAT
!===============================================================================

  subroutine spk_coefmat(t,h,y,A,n)

!---external references

    use constants,  only: beta, NUM_PRECS, lambda, pnl, ONE
    use global,     only: pke

!---arguments

    integer :: n
    real(8) :: t
    real(8) :: h
    Mat     :: A 
    Vec     :: y

!---local variables

    integer :: i   ! loop counter
    real(8) :: rho ! interpolated reactivity
    real(8) :: val ! temp value for matrix setting

    real(8), pointer :: yptr(:)
    real(8), pointer :: dfdtptr(:)

!---begin execution

    ! get reactivity
!    rho = get_reactivity(t) 

    ! set up jacobian 
    val = (rho*sum(beta) - sum(beta))/pnl
    val = ONE - h*val
    call MatSetValue(A, 0, 0, val, INSERT_VALUES, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! begin loop around rest of matrix
    do i = 2, NUM_PRECS + 1

      ! set row 1
      val = lambda(i - 1)
      val = -h*val
      call MatSetValue(A, 0, i-1, val, INSERT_VALUES, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! set diagonal
      val = -lambda(i - 1)
      val = ONE - h*val
      call MatSetValue(A, i-1, i-1, val, INSERT_VALUES, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

      ! set column 1
      val = beta(i-1) / pnl
      val = -h*val
      call MatSetValue(A, i-1, 0, val, INSERT_VALUES, mpi_err)
#     ifdef DEBUG
        CHKERRQ(mpi_err)
#     endif

    end do 

    ! finalize assembly
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, mpi_err)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine spk_coefmat 

!===============================================================================
! SET_INIT
!===============================================================================

  subroutine set_init(y)

!---external references

    use constants,  only: ONE, beta, lambda, NUM_PRECS
    use global,     only: cmfd, prod, power, geometry
    use math,       only: csr_matvec_mult

!---arguments

    Vec :: y

!---local variables

    integer :: i
    integer :: n
    real(8) :: pow
    real(8), pointer :: yptr(:)

!---begin execution  

    ! get size
    n = geometry % nf

    ! get y pointer
    call VecGetArrayF90(y, yptr, mpi_err)

    ! compute power
    pow = sum(csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val/cmfd%keff,       &
              cmfd%phi,prod%n))

    ! divide flux by power
    cmfd % phi = cmfd % phi * power / pow

    ! place flux in pointer
    yptr(1:n) = cmfd % phi 

    ! begin loop around precursors
    do i = 1, NUM_PRECS

      yptr((i-1)*n + n + 1:(i-1)*n + 2*n) = beta(i)/lambda(i) * &
       csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val/cmfd%keff,       &
       cmfd%phi,prod%n)

    end do

    ! place back pointer
    call VecRestoreArrayF90(y, yptr, mpi_err)

  end subroutine set_init

!===============================================================================
! PK_POST_TIMESTEP
!===============================================================================

  subroutine spk_post_timestep(t, y, h)

!---arguments

    real(8) :: t
    real(8) :: h
    Vec :: y

!---local variables

    real(8) :: react
    real(8), pointer :: yptr(:)

!---begin execution

    ! get ptr
    call VecGetArrayF90(y, yptr, mpi_err)

    ! compute reactivity
!   react = get_reactivity(t)

    ! write output
    print *, 'TIME:', t, 'REACT:', react, 'POWER:', yptr(1), 'STEP:', h

    ! restore ptrs
    call VecRestoreArrayF90(y, yptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine spk_post_timestep

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

  end subroutine change_data

!===============================================================================
! DESTROY_OBJECTS
!===============================================================================

  subroutine destroy_objects(dfdy, dfdt, dydt, y)

!---references

    use global,         only: loss, prod
    use loss_operator,  only: destroy_M_operator
    use prod_operator,  only: destroy_F_operator

!---arguments

    Mat :: dfdy
    Vec :: dfdt
    Vec :: dydt
    Vec :: y

!---begin execution

    ! destroy all
    call MatDestroy(dfdy, mpi_err)
    call VecDestroy(dfdt, mpi_err)
    call VecDestroy(dydt, mpi_err)
    call VecDestroy(y, mpi_err)

    ! finalize data objects
    call destroy_M_operator(loss)
    call destroy_F_operator(prod)

  end subroutine destroy_objects

end module spatial_kinetics 

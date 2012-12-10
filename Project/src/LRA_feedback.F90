module LRA_feedback

!-references

  use global,  only: mpi_err

!-module options

  implicit none
  private
  public :: run_LRAfeedback

!-module external references

# include "finclude/petsc.h90"

!-constants

  real(8), parameter :: VEL(2) = (/3.0e7_8, 3.0e5_8/)
  real(8), parameter :: BETA(2) = (/0.0054_8, 0.001087_8/)
  real(8), parameter :: LAMBDA(2) = (/0.00654_8, 1.35_8/)
  real(8), parameter :: ALPHA = 3.83e-11_8
  real(8), parameter :: GAM = 2.034e-3_8
  real(8), parameter :: KAPPA = 3.204e-11_8
  real(8), parameter :: NU = 3.24_8
  integer, parameter :: NUM_PRECS = 2

contains

!===============================================================================
! RUN_PKINETICS
!===============================================================================

  subroutine run_LRAfeedback()

!---external references

    use global,          only: pke, time, solver_type, geometry
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

    ! compute size (flux + precursors + temperature)
    n = geometry % nf + NUM_PRECS * geometry % nf + geometry % nf / geometry % nfg

    ! initialize data objects
    call init_data(dfdy, dfdt, dydt, y)

    ! set initial values in y
    call set_init(y)

    ! perform 4th order kaps rentrop
    select case(trim(solver_type))

      case('rk4')
        call execute_rk4(y, dfdy, dfdt, dydt, spk_derivs, spk_jacobn, n, spk_post_timestep,0.01_8)

      case DEFAULT

    end select

    call destroy_objects(dfdy, dfdt, dydt, y)

  end subroutine run_LRAfeedback

!===============================================================================
! INIT_DATA
!===============================================================================

  subroutine init_data(dfdy,dfdt,dydt,y)

!---references

    use constants,      only: ZERO
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
    allocate(d_nnz(n + NUM_PRECS*n + n/ng))
    allocate(o_nnz(n + NUM_PRECS*n + n/ng))

    ! create dfdy matrix
    d_nnz(1:n) = loss % d_nnz + NUM_PRECS + 1
    d_nnz(n+1:n+NUM_PRECS*n) = ng + 1
    d_nnz((1+NUM_PRECS)*n+1:(1+NUM_PRECS)*n+n/ng) = ng + 1 
    o_nnz = 0
    call MatCreateAIJ(PETSC_COMM_WORLD,n+NUM_PRECS*n+n/ng, n+NUM_PRECS*n+n/ng,&
                      PETSC_DETERMINE, PETSC_DETERMINE, PETSC_NULL, d_nnz,&
                      PETSC_NULL, o_nnz, dfdy, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create dfdt vector
    call VecCreateMPI(PETSC_COMM_WORLD, n+NUM_PRECS*n+n/ng, PETSC_DETERMINE, dfdt,&
                      mpi_err)
    call VecSet(dfdt, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create dydt vector
    call VecCreateMPI(PETSC_COMM_WORLD, n+NUM_PRECS*n+n/ng, PETSC_DETERMINE, dydt,&
                      mpi_err)
    call VecSet(dydt, ZERO, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

    ! create y solution vector
    call VecCreateMPI(PETSC_COMM_WORLD, n+NUM_PRECS*n+n/ng, PETSC_DETERMINE, y,&
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

  subroutine spk_jacobn(t,y,dfdy,dfdt,n,const)

!---external references

    use constants,      only: ONE, ZERO
    use global,         only: geometry, cmfd, loss, prod, material
    use loss_operator,  only: build_loss_matrix
    use material_header,  only: material_type
    use prod_operator,  only: build_prod_matrix

!---arguments

    integer :: n
    integer :: ii
    real(8) :: t
    real(8) :: const
    Mat     :: dfdy
    Vec     :: dfdt
    Vec     :: y

!---local variables

    integer :: sidx ! spatial index
    integer :: nf
    integer :: ng
    integer :: ncols
    integer :: i    ! loop counter
    integer :: irow ! row counter
    integer :: g
    integer, allocatable :: cols(:)
    real(8) :: val  ! temp value for matrix setting
    real(8), allocatable :: vals(:)

    real(8), pointer :: yptr(:)
    real(8), pointer :: dfdtptr(:)
    real(8) :: absxs

    type(material_type), pointer :: m => null()

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
    call change_data(t, yptr, dfdtptr)

    ! rebuild matrices
    call build_loss_matrix(loss)
    call build_prod_matrix(prod)

    ! perform temperature feedback
    call perform_feedback(yptr)

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
      dfdtptr(irow) = -yptr(irow)*vel(mod(irow-1,ng)+1)*dfdtptr(irow)

      ! get row of matrix M
      call MatGetRow(loss%oper, irow-1, ncols, cols, vals, mpi_err)

      ! multiply values by group velocity
      vals = vel(mod(irow-1,ng)+1)*vals
      do ii = 1, ncols
        if (cols(ii) == irow - 1) then
          vals(ii) = vals(ii) - const
          exit
        end if
      end do

      ! put values in jacobian
      call MatSetValues(dfdy, 1, irow-1, ncols, cols(1:ncols), -vals(1:ncols), &
                        INSERT_VALUES, mpi_err)

      ! put the row back
      call MatRestoreRow(loss%oper, irow-1, ncols, cols, vals, mpi_err)

      ! begin loop around precursors
      do i = 1, NUM_PRECS

        ! flux by precursors
        call MatSetValue(dfdy, irow-1, i*nf + (irow - 1), -lambda(i)*vel(mod(irow-1,ng)+1), INSERT_VALUES, mpi_err)

        ! precursors by flux
        vals = ZERO
        call MatGetRow(prod%oper, irow-1, ncols, cols, vals, mpi_err)
        vals = vals*beta(i)/cmfd%keff
        call MatSetValues(dfdy, 1, i*nf + (irow-1), ncols, cols(1:ncols), -vals, &
                          INSERT_VALUES, mpi_err)
        call MatRestoreRow(prod%oper, irow-1, ncols, cols, vals, mpi_err)

        ! precursors by precursors
        call MatSetValue(dfdy, i*nf + (irow - 1), i*nf + (irow - 1), lambda(i) + const, INSERT_VALUES, mpi_err)

      end do

      ! put in temperature on diagonal in flux equations
      sidx = (irow-1)/ng + 1
      g = mod(irow-1,ng)+1
      if (g == 1) then
        m => material(geometry % fmat_map(sidx))
        absxs = m % absorxs(g) - m % diffcof(g) * m % buckling ! get base absxs 
        val = -0.5_8*absxs*GAM*yptr((1+NUM_PRECS)*nf+sidx)**(-0.5) * &
              vel(mod(irow-1,ng)+1) * yptr(irow)
        call MatSetValue(dfdy, irow-1, (1+NUM_PRECS)*nf+sidx-1, -val, INSERT_VALUES, mpi_err)
      end if

    end do

    ! loop around temperature part
    do irow = 1,nf/ng
      vals = ZERO
      m => material(geometry % fmat_map(irow))
      vals(1) = m % fissvec(1)/NU*ALPHA
      vals(2) = m % fissvec(2)/NU*ALPHA
      call MatSetValues(dfdy, 1, (NUM_PRECS+1)*nf+irow-1, 2, &
           (/irow-1,irow/), -vals(1:2), INSERT_VALUES, mpi_err)
      call MatSetValue(dfdy, (NUM_PRECS+1)*nf+irow-1, &
           (NUM_PRECS+1)*nf+irow-1, const, INSERT_VALUES, mpi_err)
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
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'jacobian.bin', &
           FILE_MODE_WRITE, viewer, mpi_err)
    call MatView(dfdy, viewer, mpi_err)
    call PetscViewerDestroy(viewer, mpi_err)

  end subroutine spk_jacobn 

!===============================================================================
! PK_DERIVS
!===============================================================================

  subroutine spk_derivs(t,y,dydt,n)

!---external references

    use constants,      only: ONE
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
    PetscViewer :: viewer
!---begin execution

    ! get sub size
    nf = geometry % nf
    ng = geometry % nfg

    ! change the data
    call change_data(t, yptr)

    ! rebuild matrices
    call build_loss_matrix(loss)
    call build_prod_matrix(prod)
!   call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'lossmat.bin', &
!          FILE_MODE_WRITE, viewer, mpi_err)
!   call MatView(loss%oper, viewer, mpi_err)
!   call PetscViewerDestroy(viewer, mpi_err)
!   call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'prodmat.bin', &
!          FILE_MODE_WRITE, viewer, mpi_err)
!   call MatView(prod%oper, viewer, mpi_err)
!   call PetscViewerDestroy(viewer, mpi_err)

    ! get pointer to solution
    call VecGetArrayF90(y, yptr, mpi_err)
    call VecGetArrayF90(dydt, dydtptr, mpi_err)

    ! perform feedback on operators
    call perform_feedback(yptr)

    ! in flux equation compute effective production operator
    call MatAXPY(loss%oper, -(ONE-sum(beta))/cmfd%keff, prod%oper, SUBSET_NONZERO_PATTERN, mpi_err)
    call MatScale(loss%oper, -ONE, mpi_err)

    ! muliply by flux and store
    dydtptr(1:nf) = csr_matvec_mult(loss%row_csr+1,loss%col+1,loss%val,yptr(1:nf),prod%n)

    ! loop through precursors
    do i = 1,NUM_PRECS
 
      ! from flux equation
      dydtptr(1:nf) = dydtptr(1:nf) + lambda(i)*yptr((i-1)*nf + nf + 1:(i-1)*nf + 2*nf)

      ! from precursor equation
      dydtptr((i-1)*nf + nf + 1:(i-1)*nf + 2*nf) = beta(i)/cmfd%keff*csr_matvec_mult(prod%row_csr+1, &
                                               prod%col+1,prod%val,yptr(1:nf),prod%n) - lambda(i)*&
                                               yptr((i-1)*nf + nf + 1:(i-1)*nf + 2*nf)

    end do

    ! multiply flux equation by velocity
    do i = 1,nf
      dydtptr(i) = dydtptr(i)*vel(mod(i-1,ng)+1)
    end do

    ! calculate temperatures
    dydtptr((1+NUM_PRECS)*nf + 1: (1+NUM_PRECS)*nf+nf/ng) = &
                      calc_fiss_rate(yptr(1:nf), ALPHA/NU/cmfd%keff, nf, nf/ng, vol=.false.)

    ! put the pointer back
    call VecRestoreArrayF90(y, yptr, mpi_err)
    call VecRestoreArrayF90(dydt, dydtptr, mpi_err)

  end subroutine spk_derivs

!===============================================================================
! SET_INIT
!===============================================================================

  subroutine set_init(y)

!---external references

    use constants,  only: ONE
    use global,     only: cmfd, prod, power, geometry, fuel_T
    use math,       only: csr_matvec_mult

!---arguments

    Vec :: y

!---local variables

    integer :: i
    integer :: n, ng
    real(8) :: pow
    real(8), pointer :: yptr(:)

!---begin execution  

    ! get size
    n = geometry % nf
    ng = geometry % nfg

    ! get y pointer
    call VecGetArrayF90(y, yptr, mpi_err)

    ! compute power
    pow = sum(calc_fiss_rate(cmfd%phi, KAPPA/NU/cmfd%keff, n, n/ng, vol=.true.))

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

    ! set initial temperature
    yptr((1+NUM_PRECS)*n + 1:(1+NUM_PRECS)*n + n/ng) = fuel_T

    ! place back pointer
    call VecRestoreArrayF90(y, yptr, mpi_err)

  end subroutine set_init

!===============================================================================
! PK_POST_TIMESTEP
!===============================================================================

  subroutine spk_post_timestep(t, y, h)

!---references

    use global,  only: prod, cmfd, geometry, material
    use math,    only: csr_matvec_mult

!---arguments

    real(8) :: t
    real(8) :: h
    Vec :: y

!---local variables

    integer :: n, ng
    real(8) :: pow
    real(8), pointer :: yptr(:)

!---begin execution

    ! get size
    n = geometry % nf
    ng = geometry % nfg

    ! get ptr
    call VecGetArrayF90(y, yptr, mpi_err)

    pow = sum(calc_fiss_rate(yptr(1:n), KAPPA/NU/cmfd%keff, n, n/ng, vol=.true.)) &
         * cmfd % pfactor

    ! write output
    print *, 'TIME:', t, 'POWER:', pow, 'STEP:', h, material(6)%absorxs(2)

    ! restore ptrs
    call VecRestoreArrayF90(y, yptr, mpi_err)
#   ifdef DEBUG
      CHKERRQ(mpi_err)
#   endif

  end subroutine spk_post_timestep

!===============================================================================
! CHANGE_DATA
!===============================================================================

  subroutine change_data(t, yptr, dfdtptr)

!---external references

    use constants,        only: ONE, ZERO
    use global,           only: material, kinetics, geometry
    use kinetics_header,  only: kinetics_type
    use material_header,  only: material_type

!---arguments

    real(8) :: t
    real(8), pointer :: yptr(:)
    real(8), optional, pointer :: dfdtptr(:)

!---local variables

    integer :: i
    integer :: g
    integer :: nf, ng
    integer :: sidx
    integer :: matidx
    type(kinetics_type), pointer :: k => null()
    type(material_type), pointer :: m => null()

!---begin execution

    ! get size
    nf = geometry % nf
    ng = geometry % nfg

    ! change material 6
    if (t <= 2) then

      ! set pointers
      m => material(6)
      k => kinetics(1)

      ! undo buckling for all material 6
      m % absorxs(2) = m % absorxs(2) - m % diffcof(2) * m % buckling

      ! change value of absxs
      m % absorxs(2) = k%val(1)*(ONE - 0.0606184_8*t)

      ! put all buckling back in
      m % absorxs(2) = m % absorxs(2) + m % diffcof(2) * m % buckling

      ! calculate new removal
      m % removxs(2) = m % absorxs(2) + sum(m % scattxs(:,2)) - m % scattxs(2,2)

      ! adjust scattering to account for this removal
      m % scattxs(2,2) = m % totalxs(2) - m % removxs(2)

      ! dfdt 
      if (present(dfdtptr)) then

        ! loop over all flux points
        do i = 1, nf

          ! compute group
          g = mod(i-1,ng)+1

          ! compute spatial idx
          sidx = (i-1)/ng + 1

          ! get material id
          matidx = geometry % fmat_map(sidx)

          ! check for material six and compute slope
          if (matidx == 6 .and. g == 2) dfdtptr(i) = -k%val(1)*0.0606184_8

        end do

      end if

    end if

  end subroutine change_data

!===============================================================================
! PERFORM FEEDBACK
!===============================================================================

  subroutine perform_feedback(yptr)

!---references

    use constants,  only: ZERO, ONE
    use global,  only: loss, geometry, material, fuel_T
    use material_header,  only: material_type

!---arguments

    real(8), pointer :: yptr(:)

!---local variables

    integer :: i, ii
    integer :: g, sidx
    integer :: nf, ng
    real(8) :: val(1)
    real(8) :: absxs 
    type(material_type), pointer :: m => null()
    real(8), allocatable :: val_tmp(:)
    Vec :: D

!---begin execution

    ! get sizes
    nf = geometry % nf
    ng = geometry % nfg

    ! allocate temporary vector
    allocate(val_tmp(nf))
    val_tmp = ZERO

    ! begin loop around space
    do i = 1,nf

      ! get energy group
      g = mod(i-1,ng)+1

      ! check for group 1
      if (g == 1) then

        ! get spatial location
        sidx = (i-1)/ng + 1

        ! get pointer to material
        m => material(geometry % fmat_map(sidx))

        ! take away buckling
        m % absorxs = m % absorxs - m % diffcof * m % buckling

        ! calculate what absorption should be
        absxs = m % absorxs(g) * (ONE + GAM*(sqrt(yptr((1+NUM_PRECS)*nf+sidx)) - sqrt(fuel_T)))

        ! adjust absorption xs
        val_tmp(i) = absxs - m % absorxs(g)

        ! put back in buckling
        m % absorxs = m % absorxs + m % diffcof * m % buckling

      end if
    end do

    ! place values in petsc vector
    call VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nf, val_tmp, D, mpi_err)

    ! put in diagonal
    call MatDiagonalSet(loss%oper, D, ADD_VALUES, mpi_err)

    ! free memory
    deallocate(val_tmp)

  end subroutine perform_feedback

!===============================================================================
! CALC_FISS_RATE
!===============================================================================

  function calc_fiss_rate(flux, const, nin, nout, vol) result(fissrate)

!---references

    use global,  only: geometry, material
    use material_header,  only: material_type

!---arguments

    integer, intent(in)   :: nin, nout
    logical, intent(in)   :: vol
    real(8), intent(in)   :: flux(nin)
    real(8), intent(in)   :: const
    real(8) :: fissrate(nout)

!---local variables

    integer :: i
    integer :: nf, ng
    type(material_type), pointer :: m => null()

!---begin execution

    ! get sizes
    nf = geometry % nf
    ng = geometry % nfg

    ! begin loop over space
    do i = 1,nf/ng

      ! point to material
      m => material(geometry % fmat_map(i)) 
      fissrate(i) = sum(m%fissvec*flux(ng*(i-1)+1:i*ng))*const

      ! check for volume adjust
      if (vol) fissrate(i) = fissrate(i) * geometry % fvol_map(i)

    end do

  end function calc_fiss_rate

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

end module LRA_feedback

module point_kinetics

!-module options

  implicit none
  private
  public :: run_pkes, generate_pkparams, generate_pke_shapes

!-module external references

# include "finclude/petsc.h90"

contains

!===============================================================================
! RUN_PKES
!===============================================================================

  subroutine run_pkes()

!---external references

    use constants,  only: NUM_PRECS
    use global,     only: nt, dt, pke 
    use output,     only: header
    use pke_header, only: allocate_pke_type
    use math,       only: expm_pade

!---local variables

    integer :: i  ! loop counter

!---begin execution

    ! allocate point kinetics
    call allocate_pke_type(pke,nt)

    ! set up initial conditions
    call set_init()

    ! begin loop through time steps
    do i = 1, nt

      ! set up coefficient matrix
      call setup_coefmat(i)

      ! solve matrix exponential
      call expm_pade(pke % coef, NUM_PRECS + 1, dt, pke % expm)  ! my code

      ! get new vector
      pke % N(:,i+1) = matmul(pke % expm, pke % N(:,i))
write(774,*) pke%N(1,i+1)
    end do

  end subroutine run_pkes

!===============================================================================
! SETUP_COEFMAT
!===============================================================================

  subroutine setup_coefmat(iter)

!---external references

    use constants,  only: beta, NUM_PRECS, lambda
    use global,     only: pke, cmfd

!---arguments

    integer :: iter

!---local variables

    integer :: i ! loop counter

!---begin execution

    ! set up coeffient matrix manually for time 0
    pke % coef(1,1) = (cmfd % rho(iter)*sum(beta) - sum(beta))/cmfd % pnl(1)

    ! begin loop around rest of matrix
    do i = 2, NUM_PRECS + 1

      ! set row 1
      pke % coef(1,i) = lambda(i - 1)

      ! set diagonal
      pke % coef(i,i) = -lambda(i - 1)

      ! set column 1
      pke % coef(i,1) = beta(i-1) / cmfd % pnl(1)

    end do
write(432,*) cmfd % rho(iter),sum(beta),cmfd % pnl(iter) 
  end subroutine setup_coefmat

!===============================================================================
! SET_INIT
!===============================================================================

  subroutine set_init()

!---external references

    use constants,  only: ONE, beta, lambda, NUM_PRECS
    use global,     only: pke, cmfd

!---local variables

    integer :: i ! loop counter

!---begin execution  

    ! set power at 1.0
    pke % N(1,1) = ONE

    ! loop through precursors
    do i = 1, NUM_PRECS

      ! set initial value
      pke % N(i+1,1) = beta(i)/(cmfd % pnl(1)*lambda(i)) * pke % N(1,1)

    end do

  end subroutine set_init

!===============================================================================
! GENERATE PKE SHAPES
!===============================================================================

  subroutine generate_pke_shapes()

!---references

    use error,  only: fatal_error
    use global,  only: n_pkes, pke_shape, material, message
    use kinetics_header,  only: kinetics_type
    use material_header,  only: material_type
    use math,     only: csr_gauss_seidel
    use power_iter, only: power_execute

!---local variables

    integer :: i
    type(kinetics_type), pointer :: k
    type(material_type), pointer :: m

!---begin execution

    ! modify the data
    do i = 1, n_pkes

      ! point to kinetics object
      k => pke_shape(i)

      ! point to material object
      m => material(k % mat_id)

      ! begin case structure to replace value
      select case (trim(k % xs_id))

        case ('absxs')

          ! remove scattering component from total xs
          m % totalxs(k % g) = m % totalxs(k % g) - sum(m % scattxs(:,k % g))

          ! change absorption
          m % absorxs(k % g) = pke_shape(i) % val(1)
          m % totalxs(k % g) = pke_shape(i) % val(1)

          ! re-add back in scattering
          m % totalxs(k % g) = m % totalxs(k % g) + sum(m % scattxs(:,k % g))

        case DEFAULT

          message = 'Kinetics modification not supported!'
          call fatal_error()

      end select

    end do

    ! call power iteration
    call power_execute(csr_gauss_seidel,'none')
stop
  end subroutine generate_pke_shapes

!===============================================================================
! GENERATE_PKPARAMS
!===============================================================================

  subroutine generate_pkparams()

!---references

    use constants,        only: beta
    use global,           only: nt, dt, loss, prod, cmfd, mpi_err, kinetics
    use kinetics_solver,  only: change_data
    use loss_operator,    only: init_M_operator, build_loss_matrix,            &
                                destroy_M_operator
    use math,             only: csr_matvec_mult
    use prod_operator,    only: init_F_operator, build_prod_matrix,            &
                                destroy_F_operator

!---local variables

    integer :: i
    real(8) :: curr_time
    real(8) :: rho_num
    real(8) :: rho_den
    real(8), allocatable :: temp1(:)
    real(8), allocatable :: temp2(:)

!---begin execution

    ! initialize operators
    call init_M_operator(loss)
    call init_F_operator(prod)

    ! allocate vectors
    allocate(temp1(loss % n))
    allocate(temp2(prod % n))
    if(.not.allocated(cmfd % rho)) allocate(cmfd % rho(nt))

    ! build petsc matrices
    call build_loss_matrix(loss,'')
    call build_prod_matrix(prod,'')
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,loss%n,loss%n,loss%row_csr,&
                                   loss%col,loss%val,loss%oper,mpi_err)
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,prod%n,prod%n,prod%row_csr,&
                                   prod%col,prod%val,prod%oper,mpi_err)

    ! begin loop around timesteps
    do i = 1, nt

      ! compute current time
      curr_time = dble(i)*dt

      ! change the material via kinetics mods
      call change_data(kinetics,curr_time)

      ! build operators
      call build_loss_matrix(loss,'')
      call build_prod_matrix(prod,'')

      ! we need F - M, store it back in M
      call MatAXPY(loss%oper,-1.0_8/cmfd%keff,prod%oper, SUBSET_NONZERO_PATTERN, mpi_err)
      call MatScale(loss%oper,-1.0_8,mpi_err)

      ! perform reactivity numerator operator multiplication
      temp1 =  csr_matvec_mult(loss%row_csr+1,loss%col+1,loss%val,cmfd%phi,loss%n)
      temp2 =  csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val,cmfd%phi,prod%n)

      ! multiply by volume
      call multiply_volume(temp1,loss%n)
      call multiply_volume(temp2,loss%n)

      ! perform dot product
      rho_num = dot_product(cmfd%phi_adj,temp1)
      rho_den = dot_product(cmfd%phi_adj,temp2)*1.0_8/cmfd%keff

      ! calc reactivity
      cmfd % rho(i) = rho_num/rho_den
write(835,*) cmfd % rho(i) / sum(beta)
    end do

    ! deallocate vectors
    deallocate(temp1)
    deallocate(temp2)

  end subroutine generate_pkparams 

!===============================================================================
! MULTIPLY_VOLUME
!===============================================================================

  subroutine multiply_volume(vec,n)

!---external references

    use global,  only: geometry 

!---arguments

    integer :: n
    real(8) :: vec(n) 

!---local variables

    integer :: irow
    integer :: idx
    real(8) :: vol

!---begin execution

    ! begin loop around rows
    do irow = 1, n

      ! compute index
      idx = ceiling(real(irow)/real(geometry%nfg))

      ! get region number
      vol = geometry % fvol_map(idx)

      ! sum nodal power
      vec(irow) = vec(irow) * vol

    end do

  end subroutine multiply_volume



end module point_kinetics

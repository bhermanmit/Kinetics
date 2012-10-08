module point_kinetics

!-module options

  implicit none
  private
  public :: generate_pkparams

!-module external references

# include "finclude/petsc.h90"

contains

!===============================================================================
! GENERATE_PKPARAMS
!===============================================================================

  subroutine generate_pkparams()

!---references

    use global,           only: nt, dt, loss, prod, cmfd, mpi_err
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
    call build_loss_matrix(loss)
    call build_prod_matrix(prod)
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,loss%n,loss%n,loss%row_csr,&
                                   loss%col,loss%val,loss%oper,mpi_err)
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,prod%n,prod%n,prod%row_csr,&
                                   prod%col,prod%val,prod%oper,mpi_err)

    ! begin loop around timesteps
    do i = 1, nt

      ! compute current time
      curr_time = dble(i)*dt

      ! change the material via kinetics mods
      call change_data(curr_time)

      ! build operators
      call build_loss_matrix(loss)
      call build_prod_matrix(prod)

      ! we need F - M, store it back in M
      call MatAXPY(loss%oper,-1.0_8/cmfd%keff,prod%oper, SUBSET_NONZERO_PATTERN, mpi_err)

      ! perform reactivity numerator operator multiplication
      temp1 =  csr_matvec_mult(loss%row_csr+1,loss%col+1,loss%val,cmfd%phi,loss%n)
      temp2 =  csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val,cmfd%phi,prod%n)

      ! multiply by volume
      call multiply_volume(temp1,loss%n)
      call multiply_volume(temp2,loss%n)

      ! perform dot product
      rho_num = dot_product(cmfd%phi_adj,temp1)
      rho_den = dot_product(cmfd%phi_adj,temp2)

      ! calc reactivity
      cmfd % rho(i) = rho_num/rho_den
write(834,*) cmfd % rho(i)
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

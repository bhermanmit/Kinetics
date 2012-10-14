module general_pkes 

!-module options

  implicit none
  private
  public :: run_gpkes, generate_gpkparams, generate_gpke_shapes

contains

!===============================================================================
! RUN_PKES
!===============================================================================

  subroutine run_gpkes()

!---external references

    use constants,  only: NUM_PRECS
    use global,     only: nt, dt, gpke, pke_grp, cmfd 
    use pke_header, only: allocate_pke_type

!---local variables

    integer :: i,j,k  ! loop counter

!---begin execution

    ! allocate point kinetics
    call allocate_pke_type(gpke,pke_grp,nt)

    ! set up initial conditions
    call set_init()
write(888,*) sum(gpke%N(1:pke_grp,1))
    ! begin loop through time steps
    do i = 1, nt

      ! set up coefficient matrix
      call setup_coefmat(i)
! write(333,*) gpke % coef

      ! solve matrix exponential
      call expm1(pke_grp + NUM_PRECS*pke_grp, gpke%coef*dt, gpke % expm)

      ! get new vector
      gpke % N(:,i+1) = matmul(gpke % expm, gpke % N(:,i))
write(888,*) sum(gpke%N(1:pke_grp,i+1))

    end do

  end subroutine run_gpkes

!===============================================================================
! SETUP_COEFMAT
!===============================================================================

  subroutine setup_coefmat(iter)

!---external references

    use constants,  only: beta, NUM_PRECS, lambda, ZERO
    use global,     only: gpke, cmfd, pke_grp

!---arguments

    integer :: iter

!---local variables

    integer :: g,h,i 

!---begin execution

    ! ZERO out coefficient matrix
    gpke % coef = ZERO

    ! begin loop over group amplitude functions
    GROUPG: do g = 1,pke_grp

      ! loop around energy groups
      do h = 1, pke_grp

        ! do amplitude
        gpke % coef(g,h) = cmfd % vel(g,iter)*cmfd % prompt(g,h,iter)

      end do

      ! loop around precursor groups
      do i = 1, NUM_PRECS 

        ! set value in amplitude rows
        gpke % coef(g,pke_grp+g+pke_grp*(i-1)) = cmfd % vel(g,iter)*lambda(i)

        ! set value in precurosr row
        gpke % coef(pke_grp+g+pke_grp*(i-1),pke_grp+g+pke_grp*(i-1)) = -lambda(i) 

        ! loop around energy group
        do h = 1,pke_grp

          gpke % coef(pke_grp+g+pke_grp*(i-1),h) = cmfd % delay(i,g,h,iter)

        end do

      end do
 
    end do GROUPG

  end subroutine setup_coefmat

!===============================================================================
! SET_INIT
!===============================================================================

  subroutine set_init()

!---external references

    use constants,  only: ONE, beta, lambda, NUM_PRECS, ZERO
    use global,     only: gpke, cmfd, pke_grp

!---local variables

    integer :: i, g, h ! loop counter

!---begin execution  

    ! set power at 1.0
    gpke % N(1:pke_grp,1) = cmfd % fsrc / sum(cmfd % fsrc) 
    gpke % N(pke_grp+1:pke_grp+pke_grp*NUM_PRECS,1) = ZERO
    gpke % N(1:pke_grp,1) = 0.5_8
    ! loop through precursors
    do i = 1, NUM_PRECS

      ! loop around group g
      do g = 1, pke_grp

        ! loop around group h
        do h = 1, pke_grp

          ! sum value
          gpke % N(pke_grp+g+pke_grp*(i-1),1) = gpke % N(pke_grp+g+pke_grp*(i-1),1) + &
          ONE/lambda(i) * cmfd % delay(i,g,h,1) * gpke % N(h,1)

        end do

      end do

    end do
print *,'HERE'
print *, gpke % N(:,1)
  end subroutine set_init

!===============================================================================
! GENERATE PKE SHAPES
!===============================================================================

  subroutine generate_gpke_shapes()

!---references

    use constants,  only: ONE
    use error,  only: fatal_error
    use global,  only: n_pkes_for, pke_shape_for, n_pkes_adj, pke_shape_adj,   &
                       material, message, cmfd, weight, adjoint
    use kinetics_header,  only: kinetics_type
    use material_header,  only: material_type
    use math,     only: csr_gauss_seidel
    use power_iter, only: power_execute

!---local variables

    integer :: i
    real(8) :: temp
    type(kinetics_type), pointer :: k
    type(material_type), pointer :: m

!---begin execution

    ! bank keff in case shape functions are different
    cmfd % kcrit = cmfd % keff

    ! modify the data
    do i = 1, n_pkes_for

      ! point to kinetics object
      k => pke_shape_for(i)

      ! point to material object
      m => material(k % mat_id)

      ! begin case structure to replace value
      select case (trim(k % xs_id))

        case ('absxs')

          ! remove scattering component from total xs
          m % totalxs(k % g) = m % totalxs(k % g) - sum(m % scattxs(:,k % g))

          ! change absorption
          temp = m % absorxs(k % g)
          m % absorxs(k % g) = pke_shape_for(i) % val(1)
          m % totalxs(k % g) = pke_shape_for(i) % val(1)

          ! re-add back in scattering
          m % totalxs(k % g) = m % totalxs(k % g) + sum(m % scattxs(:,k % g))

        case DEFAULT

          message = 'Kinetics modification not supported!'
          call fatal_error()

      end select

    end do

    ! call power iteration
    call power_execute(csr_gauss_seidel,'none')

    ! we always generate pke parameters, check for weighting function
    if (trim(weight) == 'adjoint') then

      ! modify the data
      do i = 1, n_pkes_for

        ! point to kinetics object
        k => pke_shape_for(i)

        ! point to material object
        m => material(k % mat_id)

        ! begin case structure to replace value
        select case (trim(k % xs_id))

          case ('absxs')

            ! remove scattering component from total xs
            m % totalxs(k % g) = m % totalxs(k % g) - sum(m % scattxs(:,k % g))

            ! change absorption
            temp = m % absorxs(k % g)
            m % absorxs(k % g) = pke_shape_adj(i) % val(1)
            m % totalxs(k % g) = pke_shape_adj(i) % val(1)

            ! re-add back in scattering
            m % totalxs(k % g) = m % totalxs(k % g) + sum(m % scattxs(:,k % g))

          case DEFAULT

            message = 'Kinetics modification not supported!'
            call fatal_error()

        end select

      end do

      call power_execute(csr_gauss_seidel,trim(adjoint))
    else
      cmfd % phi_adj = ONE
    end if

  end subroutine generate_gpke_shapes

!===============================================================================
! GENERATE_PKPARAMS
!===============================================================================

  subroutine generate_gpkparams()

!---references

    use constants,        only: beta, NUM_PRECS
    use global,           only: nt, dt, loss, prod, cmfd, mpi_err, kinetics, pke_grp 
    use kinetics_solver,  only: change_data, compute_gpkes
    use loss_operator,    only: init_M_operator, build_loss_matrix,            &
                                destroy_M_operator
    use math,             only: csr_matvec_mult
    use prod_operator,    only: init_F_operator, build_prod_matrix,            &
                                destroy_F_operator

!---local variables

    integer :: i
    real(8) :: curr_time

!---begin execution

    ! initialize operators
    call init_M_operator(loss)
    call init_F_operator(prod)

    ! allocate
    if(.not.allocated(cmfd % prompt)) allocate(cmfd % prompt(pke_grp,pke_grp,nt))
    if(.not.allocated(cmfd % delay)) allocate(cmfd % delay(NUM_PRECS,pke_grp,pke_grp,nt))
    if(.not.allocated(cmfd % vel)) allocate(cmfd % vel(pke_grp,nt))

    ! begin loop around timestep
    do i = 1, nt

      ! compute current time
      curr_time = dble(i)*dt

      ! change the material via kinetics mods
      call change_data(kinetics,curr_time)

      ! compute the pke parameters
      call compute_gpkes(i)

    end do

  end subroutine generate_gpkparams 

end module general_pkes 

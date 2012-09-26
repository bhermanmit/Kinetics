module kinetics_solver 

!-module references

!-module options

  implicit none
  private
  public :: kinetics_execute 

!-module variables

contains

!===============================================================================
! KINETICS_EXECUTE
!===============================================================================

  subroutine kinetics_execute(inner_solver)

!---arguments

    external :: inner_solver

!---begin execution

    ! initialize the time 0 data
    call init_data()

    ! run through the kinetics iterations
    call execute_kinetics_iter(inner_solver)
 
  end subroutine kinetics_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

!---external references

    use cmfd_header,        only: compute_core_power
    use constants,          only: ONE
    use global,             only: cmfd, geometry, material, kine
    use kinetics_operator,  only: init_K_operator

!---local variables

    real(8) :: pow

!---begin execution

    ! normalize initial power to unity and set initial power
    pow = compute_core_power(cmfd, size(cmfd%phi), geometry, material)
    cmfd % phi = cmfd % phi * ONE / pow

    ! compute steady state precursors
    call compute_initial_precursors()

    ! set up matrices
    call init_K_operator(kine)

  end subroutine init_data

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

    use global,             only: nt, dt, kine, cmfd
    use kinetics_operator,  only: build_kinetics_matrix

!---arguments

    external :: inner_solver

!---local variables

    integer :: i
    integer :: n
    real(8) :: curr_time
    real(8), allocatable :: rhs(:)

!--begin execution

   ! allocate RHS
   n = size(cmfd % phi)
   allocate(rhs(n))

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
 
     ! solve matrices

     ! compute end of time step precursor concentration

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

        cmfd % C(j,i) = beta(i)/((ONE + dt*lambda(i)) * cmfd % keff) * fiss +  &
                        cmfd % C(j,i)/(ONE + dt*lambda(i))

      end do

    end do

  end subroutine compute_final_precursors

!===============================================================================
! CHANGE_DATA
!===============================================================================

  subroutine change_data(t)

!---external references

    use error,            only: fatal_error
    use global,           only: material, kinetics, n_kins, message
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
      g = mod(irow,geometry%nfg) + 1

      ! set material pointer
      m => material(geometry % fmat_map(idx))

      ! compute the precursor
      rhs(irow) = cmfd % phi(irow)/(vel(g)*dt) + sum((m%chid(g)*lambda*        &
                                            cmfd % C(irow,:))/(ONE + lambda*dt))
    end do


  end subroutine build_rhs

!==============================================================================
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize data objects

  end subroutine finalize

end module kinetics_solver 

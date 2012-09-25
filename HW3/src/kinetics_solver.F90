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

    call init_data()
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

    ! set up matrices
    call init_K_operator(kine)

  end subroutine init_data

!===============================================================================
! EXECUTE_KINETICS_ITER  in the main kinetics iteration routine 
!                         for the cmfd calculation
!===============================================================================

  subroutine execute_kinetics_iter(inner_solver)

!---external references

    use global,  only: nt, dt, material

!---arguments

    external :: inner_solver

!---local variables

    integer :: i
    real(8) :: curr_time

!--begin execution

   ! begin loop around time
   do i = 1, nt

     ! compute current time
     curr_time = dble(i)*dt
    
     ! change the material via kinetics mods
     call change_data(curr_time)

     write(44,*) material(5) % absorxs(2)

   end do

  end subroutine execute_kinetics_iter 

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
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize data objects

  end subroutine finalize

end module kinetics_solver 

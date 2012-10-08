module cmfd_header

!-module options

  implicit none
  private
  public :: allocate_cmfd_type, deallocate_cmfd_type, calc_power,              &
            compute_core_power

!-module variables

  type, public :: cmfd_type

    ! eigenvector/eigenvalue from cmfd run
    real(8), allocatable :: phi(:)
    real(8), allocatable :: phi_o(:)
    real(8), allocatable :: phi_adj(:)
    real(8) :: keff = 0.0_8

    ! nodal powers
    real(8), allocatable :: power_o(:)
    real(8), allocatable :: power_n(:)

    ! precursors
    real(8), allocatable :: C(:,:)

    ! transient
    real(8), allocatable :: time(:)
    real(8), allocatable :: core_power(:)

    ! final iteration and norm
    integer :: iter
    real(8) :: norm
    real(8) :: norm_o
    real(8) :: dr

    ! reactivity for point kinetics
    real(8), allocatable :: rho(:)

  end type cmfd_type

contains

!===============================================================================
! ALLOCATE_CMFD_TYPE
!===============================================================================

  subroutine allocate_cmfd_type(this,geometry)

!---external references

    use geometry_header,  only: geometry_type

!---arguments

    integer :: n
    type(cmfd_type) :: this
    type(geometry_type) :: geometry

!---begin execution

    ! allocate phi
    allocate(this % phi(geometry % nfx * geometry % nfy * geometry % nfz *     &
                        geometry % nfg))
    allocate(this % phi_adj(geometry % nfx * geometry % nfy * geometry % nfz * &
                        geometry % nfg))

    ! allocate nodal power maps
    allocate(this % power_o(geometry % n_regs))
    allocate(this % power_n(geometry % n_regs))

  end subroutine allocate_cmfd_type

!===============================================================================
! CALC_POWER
!===============================================================================

  subroutine calc_power(this,fsrc,n,geometry)

!---external references

    use constants,        only: ZERO
    use geometry_header,  only: geometry_type

!---arguments

    integer :: n
    type(cmfd_type) :: this
    type(geometry_type) :: geometry
    real(8) :: fsrc(n)

!---local variables

    integer :: irow
    integer :: i
    integer :: j
    integer :: k
    integer :: g
    integer :: idx
    integer :: reg
    real(8) :: vol

!---begin execution

    ! zero out new power
    this % power_n = ZERO

    ! begin loop around rows
    do irow = 1, n

      ! compute index
      idx = ceiling(real(irow)/real(geometry%nfg))

      ! get region number
      reg = geometry % freg_map(idx)
      vol = geometry % fvol_map(idx)

      ! sum nodal power
      this % power_n(reg) = this % power_n(reg) + fsrc(irow)*vol

    end do

    ! normalize power
    this % power_n = this % power_n / sum(this % power_n) *                    &
                     count(this % power_n > 1.e-11_8)

  end subroutine calc_power

!===============================================================================
! COMPUTE_CORE_POWER
!===============================================================================

  function compute_core_power(this, n, geometry, material) result(power)

!---external references

    use constants,        only: ZERO
    use geometry_header,  only: geometry_type
    use material_header,  only: material_type

!---arguments

    integer :: n
    type(cmfd_type) :: this
    type(geometry_type)          :: geometry
    type(material_type), target  :: material(:)
    real(8) :: power

!---local variables

    integer :: irow
    integer :: g
    integer :: idx
    real(8) :: vol
    type(material_type), pointer :: m

!---begin execution

    ! zero out the core power
    power = ZERO

    ! begin loop around rows
    do irow = 1, n

      ! get index
      idx = ceiling(real(irow)/real(geometry%nfg))
      g = mod(irow,geometry%nfg) + 1

      ! get volume
      vol = geometry % fvol_map(idx)

      ! set material pointer
      m => material(geometry % fmat_map(idx))

      ! accumulate power
      power = power + m % fissvec(g) * this % phi(irow) * vol / this % keff

    end do

  end function compute_core_power

!===============================================================================
! DEALLOCATE_CMFD_TYPE
!===============================================================================

  subroutine deallocate_cmfd_type(this)

!---arguments

    type(cmfd_type) :: this

!---begin execution

    if(allocated(this % phi)) deallocate(this % phi)
    if(allocated(this % phi_o)) deallocate(this % phi_o)
    if(allocated(this % phi_adj)) deallocate(this % phi_adj)
    if(allocated(this % power_o)) deallocate(this % power_o)
    if(allocated(this % power_n)) deallocate(this % power_n)

  end subroutine deallocate_cmfd_type

end module cmfd_header

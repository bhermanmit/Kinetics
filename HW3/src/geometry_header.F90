module geometry_header

!-module options

  implicit none
  private
  public :: allocate_geometry_type, deallocate_geometry_type,                  &
            generate_fine_map, compute_widths

!-module variables

  type, public :: geometry_type

    ! coarse indices
    integer :: ncx 
    integer :: ncy
    integer :: ncz
    integer :: ncg 

    ! fine indices
    integer :: nfx
    integer :: nfy
    integer :: nfz
    integer :: nfg

    ! number of regions
    integer :: n_regs

    ! material map
    integer, allocatable :: mat_map(:,:,:)

    ! region map
    integer, allocatable :: reg_map(:,:,:)

    ! fine to coarse map
    integer, allocatable :: fmat_map(:)
    integer, allocatable :: freg_map(:)
    real(8), allocatable :: fdx_map(:)
    real(8), allocatable :: fdy_map(:)
    real(8), allocatable :: fdz_map(:)
    real(8), allocatable :: fvol_map(:)

    ! size of each coarse mesh
    real(8), allocatable :: xgrid(:)
    real(8), allocatable :: ygrid(:)
    real(8), allocatable :: zgrid(:)
    real(8) :: vol

    ! fine mesh per coarse mesh
    integer, allocatable :: nnx(:)
    integer, allocatable :: nny(:)
    integer, allocatable :: nnz(:)

    ! fine mesh spacing
    real(8), allocatable :: dx(:)
    real(8), allocatable :: dy(:)
    real(8), allocatable :: dz(:)

    ! boundary conditions
    real(8) :: bc(6)

  end type geometry_type 


contains

!===============================================================================
! ALLOCATE_GEOMETRY_TYPE
!===============================================================================

  subroutine allocate_geometry_type(this)

!---arguments

    type(geometry_type) :: this

!---begin execution

    ! allocate maps
    allocate(this % mat_map (this % ncx, this % ncy, this % ncz))
    allocate(this % reg_map (this % ncx, this % ncy, this % ncz))
    allocate(this % fmat_map(this % nfx * this % nfy * this % nfz))
    allocate(this % freg_map(this % nfx * this % nfy * this % nfz))
    allocate(this % fdx_map(this % nfx * this % nfy * this % nfz))
    allocate(this % fdy_map(this % nfx * this % nfy * this % nfz))
    allocate(this % fdz_map(this % nfx * this % nfy * this % nfz))
    allocate(this % fvol_map(this % nfx * this % nfy * this % nfz))

    ! allocate grids
    allocate(this % xgrid(this % ncx))
    allocate(this % ygrid(this % ncy))
    allocate(this % zgrid(this % ncz))

    ! allocate mesh points
    allocate(this % nnx(this % ncx))
    allocate(this % nny(this % ncy))
    allocate(this % nnz(this % ncz))

    ! allocate spacings
    allocate(this % dx(this % ncx))
    allocate(this % dy(this % ncy))
    allocate(this % dz(this % ncz))

  end subroutine allocate_geometry_type

!===============================================================================
! GENERATE_FINE_MAP
!===============================================================================

  subroutine generate_fine_map(this)

!---arguments

    type(geometry_type) :: this

!---local variables

    integer :: i  ! loop counter
    integer :: j  ! loop counter
    integer :: k  ! loop counter
    integer :: ii ! loop counter
    integer :: jj ! loop cuonter
    integer :: kk ! loop counter
    integer :: ix ! index counter
    integer :: iy ! index counter
    integer :: iz ! index counter
    integer :: n

!---begin execution

    ! begin loop over coarse map
    iz = 0
    CZLOOP: do k = 1, this % ncz

      iy = 0
      CYLOOP: do j = 1, this % ncy

        ix = 0
        CXLOOP: do i = 1, this % ncx

          ! begin loop over fine mesh
          FZLOOP: do kk = 1, this % nnz(k)

            FYLOOP: do jj = 1, this % nny(j)

              FXLOOP: do ii = 1, this % nnx(i)

                ! row
                n = (ix+ii) + this % nfx * (iy+jj - 1) + this % nfx *          &
                     this % nfy * (iz+kk - 1)

                ! save coarse to fine mesh
                this % fmat_map(n) = this % mat_map(i,j,k)
                this % freg_map(n) = this % reg_map(i,j,k)
                this % fdx_map(n)  = this % dx(i) 
                this % fdy_map(n)  = this % dy(j)
                this % fdz_map(n)  = this % dz(k)
                this % fvol_map(n) = (this % dx(i) *      &
                                   this % dy(j) * this % dz(k))

              end do FXLOOP

            end do FYLOOP

          end do FZLOOP
          ix = ix + this % nnx(i)

        end do CXLOOP
        iy = iy + this % nny(j)

      end do CYLOOP
      iz = iz + this % nnz(k)

    end do CZLOOP

  end subroutine generate_fine_map

!===============================================================================
! COMPUTE_WIDTHS
!===============================================================================

  subroutine compute_widths(this)

!---arguments

    type(geometry_type) :: this

!---begin execution

    ! compute x
    this % dx = this % xgrid / dble(this % nnx)

    ! compute y
    this % dy = this % ygrid / dble(this % nny)

    ! compute z
    this % dz = this % zgrid / dble(this % nnz)

  end subroutine compute_widths

!===============================================================================
! DEALLOCATE_GEOMETRY_TYPE
!===============================================================================

  subroutine deallocate_geometry_type(this)

!---arguments

    type(geometry_type) :: this

!---begin execution

    ! deallocate all
    deallocate(this % mat_map)
    deallocate(this % reg_map)
    deallocate(this % fmat_map)
    deallocate(this % freg_map)
    deallocate(this % fdx_map)
    deallocate(this % fdy_map)
    deallocate(this % fdz_map)
    deallocate(this % fvol_map)
    deallocate(this % xgrid)
    deallocate(this % ygrid)
    deallocate(this % zgrid)
    deallocate(this % dx)
    deallocate(this % dy)
    deallocate(this % dz)

  end subroutine deallocate_geometry_type

end module geometry_header

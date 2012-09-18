module material_header

!-module options

  implicit none
  private

!-module variables

  type, public ::  material_type

    ! cross sections
    real(8), allocatable :: totalxs(:)
    real(8), allocatable :: absorxs(:)
    real(8), allocatable :: removxs(:)
    real(8), allocatable :: scattxs(:,:)
    real(8), allocatable :: nfissxs(:,:)

    ! diffusion coefficient
    real(8), allocatable :: diffcof(:,:,:,:)

  end type material_type 

end module material_header

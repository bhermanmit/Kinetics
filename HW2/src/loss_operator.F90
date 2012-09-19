module loss_operator

  implicit none
  private
  public :: init_M_operator,build_loss_matrix,destroy_M_operator

#include "finclude/petsc.h90"

    integer  :: nx   ! maximum number of x cells
    integer  :: ny   ! maximum number of y cells
    integer  :: nz   ! maximum number of z cells
    integer  :: ng   ! maximum number of groups
    integer  :: ierr ! petsc error code

  type, public :: loss_operator_type

    Mat      :: M      ! petsc matrix for neutronic loss operator
    integer  :: n      ! dimensions of matrix
    integer  :: nnz    ! max number of nonzeros
    integer  :: localn ! local size on proc
    integer, allocatable :: d_nnz(:) ! vector of diagonal preallocation
    integer, allocatable :: o_nnz(:)   ! vector of off-diagonal preallocation
    integer, allocatable :: row(:)     ! vector of row indices
    integer, allocatable :: col(:)     ! vector of column indices
    integer, allocatable :: row_csr(:) ! row indexing for CSR
    real(8), allocatable :: val(:)     ! the value

  end type loss_operator_type

contains

!===============================================================================
! INIT_M_OPERATOR
!===============================================================================

  subroutine init_M_operator(this)

    type(loss_operator_type) :: this

    ! get indices
    call get_M_indices(this)

    ! get preallocation
    call preallocate_loss_matrix(this)

    ! set up M operator
    call MatCreateAIJ(PETSC_COMM_WORLD,this%localn,this%localn,PETSC_DECIDE,&
   & PETSC_DECIDE,PETSC_NULL_INTEGER,this%d_nnz,PETSC_NULL_INTEGER,this%o_nnz, &
   & this%M,ierr)
    call MatSetOption(this%M,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,ierr)
    call MatSetOption(this%M,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)

  end subroutine init_M_operator

!===============================================================================
! GET_M_INDICES
!===============================================================================

  subroutine get_M_indices(this)

    use global, only: geometry 

    type(loss_operator_type) :: this

    ! get maximum number of cells in each direction
    nx = geometry % nfx 
    ny = geometry % nfy
    nz = geometry % nfz
    ng = geometry % nfg

    ! get number of nonzeros
    this%nnz = 7 + ng - 1

    ! calculate dimensions of matrix
    this%n = nx*ny*nz*ng

  end subroutine get_M_indices

!===============================================================================
! PREALLOCATE_LOSS_MATRIX
!===============================================================================

  subroutine preallocate_loss_matrix(this)

    type(loss_operator_type) :: this

    integer :: rank          ! rank of processor
    integer :: sizen         ! number of procs
    integer :: i             ! iteration counter for x
    integer :: j             ! iteration counter for y
    integer :: k             ! iteration counter for z
    integer :: g             ! iteration counter for groups
    integer :: l             ! iteration counter for leakages
    integer :: h             ! energy group when doing scattering
    integer :: n             ! the extent of the matrix
    integer :: irow          ! row counter
    integer :: bound(6)      ! vector for comparing when looking for bound
    integer :: xyz_idx       ! index for determining if x,y or z leakage
    integer :: dir_idx       ! index for determining - or + face of cell
    integer :: neig_idx(3)   ! spatial indices of neighbour
    integer :: nxyz(3,2)     ! single vector containing bound. locations
    integer :: shift_idx     ! parameter to shift index by +1 or -1
    integer :: row_start     ! index of local starting row
    integer :: row_end       ! index of local final row
    integer :: neig_mat_idx  ! matrix index of neighbor cell
    integer :: scatt_mat_idx ! matrix index for h-->g scattering terms

    ! get rank and max rank of procs
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,sizen,ierr)

    ! get local problem size
    n = this%n

    ! determine local size, divide evenly between all other procs
    this%localn = n/(sizen)
        
    ! add 1 more if less proc id is less than mod
    if (rank < mod(n,sizen)) this%localn = this%localn + 1

    ! determine local starting row
    row_start = 0
    if (rank < mod(n,sizen)) then
      row_start = rank*(n/sizen+1)
    else
      row_start = min(mod(n,sizen)*(n/sizen+1)+(rank - mod(n,sizen))*(n/sizen),n) 
    end if

    ! determine local final row
    row_end = row_start + this%localn - 1

    ! allocate counters
    if (.not. allocated(this%d_nnz)) allocate(this%d_nnz(row_start:row_end))
    if (.not. allocated(this%o_nnz)) allocate(this%o_nnz(row_start:row_end))
    this % d_nnz = 0
    this % o_nnz = 0

    ! create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! begin loop around local rows
    ROWS: do irow = row_start,row_end

      ! initialize counters 
      this%d_nnz(irow) = 1 ! already add in matrix diagonal
      this%o_nnz(irow) = 0

      ! get location indices
      call matrix_to_indices(irow,g,i,j,k)

      ! create boundary vector
      bound = (/i,i,j,j,k,k/)

      ! begin loop over leakages
      LEAK: do l = 1,6

        ! define (x,y,z) and (-,+) indices
        xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
        dir_idx = 2 - mod(l,2) ! -=1, +=2

        ! calculate spatial indices of neighbor
        neig_idx = (/i,j,k/)                ! begin with i,j,k
        shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1
        neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

        ! check for global boundary
        if (bound(l) /= nxyz(xyz_idx,dir_idx)) then

          ! get neighbor matrix index
          call indices_to_matrix(g,neig_idx(1),neig_idx(2),neig_idx(3),    &
         &                       neig_mat_idx)

          ! record nonzero
          if (((neig_mat_idx-1) >= row_start) .and.                        &
         &   ((neig_mat_idx-1) <= row_end)) then
            this%d_nnz(irow) = this%d_nnz(irow) + 1
          else
            this%o_nnz(irow) = this%o_nnz(irow) + 1
          end if

        end if

      end do LEAK

      ! begin loop over off diagonal in-scattering
      SCATTR: do h = 1,ng

        ! cycle though if h=g
        if (h == g) then
          cycle
        end if

        ! get neighbor matrix index
        call indices_to_matrix(h,i,j,k,scatt_mat_idx)

        ! record nonzero
        if (((scatt_mat_idx-1) >= row_start) .and.                        &
       &   ((scatt_mat_idx-1) <= row_end)) then
          this%d_nnz(irow) = this%d_nnz(irow) + 1
        else
          this%o_nnz(irow) = this%o_nnz(irow) + 1
        end if

      end do SCATTR

    end do ROWS

    ! allocate CSR objects
    allocate(this % row(sum(this % d_nnz) + sum(this % o_nnz)))
    allocate(this % col(sum(this % d_nnz) + sum(this % o_nnz)))
    allocate(this % val(sum(this % d_nnz) + sum(this % o_nnz)))
    allocate(this % row_csr(2*(row_end - row_start + 1)))

  end subroutine preallocate_loss_matrix

!===============================================================================
! BUILD_LOSS_MATRIX creates the matrix representing loss of neutrons
!===============================================================================

  subroutine build_loss_matrix(this)

    use global,           only: geometry, material, mpi_err
    use material_header,  only: material_type

    type(loss_operator_type) :: this

    integer :: nxyz(3,2)            ! single vector containing bound. locations
    integer :: i                    ! iteration counter for x
    integer :: j                    ! iteration counter for y
    integer :: k                    ! iteration counter for z
    integer :: g                    ! iteration counter for groups
    integer :: l                    ! iteration counter for leakages
    integer :: h                    ! energy group when doing scattering
    integer :: neig_mat_idx         ! matrix index of neighbor cell
    integer :: scatt_mat_idx        ! matrix index for h-->g scattering terms
    integer :: bound(6)             ! vector for comparing when looking for bound
    integer :: xyz_idx              ! index for determining if x,y or z leakage
    integer :: dir_idx              ! index for determining - or + face of cell
    integer :: neig_idx(3)          ! spatial indices of neighbour
    integer :: shift_idx            ! parameter to shift index by +1 or -1
    integer :: kount                ! integer for counting values in vector
    integer :: row_start            ! the first local row on the processor
    integer :: row_finish           ! the last local row on the processor
    integer :: irow                 ! iteration counter over row
    real(8) :: hxyz(3)              ! cell lengths in each direction
    real(8) :: hxyzn(3)             ! cell lengths for neighbor
    real(8) :: dtilde               ! coupling factor
    real(8) :: jn                   ! direction dependent leakage coeff to neig
    real(8) :: jo(6)                ! leakage coeff in front of cell flux
    real(8) :: jnet                 ! net leakage from jo
    real(8) :: val                  ! temporary variable before saving to 
    type(material_type), pointer :: m  ! current cell materials
    type(material_type), pointer :: mn ! neighboring cell materials

    ! create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! get row bounds for this processor
    call MatGetOwnershipRange(this%M,row_start,row_finish,ierr)

    ! initialize counter
    kount = 1

    ! begin iteration loops
    ROWS: do irow = row_start,row_finish-1 

      ! set csr value of row
      this % row_csr(2*irow+1) = kount

      ! get indices for that row
      call matrix_to_indices(irow,g,i,j,k)

      ! set material pointer
      m => material(geometry % fine_map(i,j,k) % mat)

      ! retrieve cell widths 
      hxyz = (/geometry % dx(geometry % fine_map(i,j,k) % x),                  &
               geometry % dy(geometry % fine_map(i,j,k) % y),                  &
               geometry % dz(geometry % fine_map(i,j,k) % z)/)

      ! create boundary vector 
      bound = (/i,i,j,j,k,k/)

      ! begin loop over leakages
      ! 1=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z 
      LEAK: do l = 1,6

        ! define (x,y,z) and (-,+) indices
        xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
        dir_idx = 2 - mod(l,2) ! -=1, +=2

        ! calculate spatial indices of neighbor
        neig_idx = (/i,j,k/)                ! begin with i,j,k
        shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1
        neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

        ! check for global boundary
        if (bound(l) /= nxyz(xyz_idx,dir_idx)) then

          ! set neighbor material pointer
          mn => material(geometry % fine_map(neig_idx(1), neig_idx(2),         &
                                              neig_idx(3)) % mat)

          ! retrive neighbor cell widths
          hxyzn = (/geometry % dx(geometry % fine_map(neig_idx(1),neig_idx(2), &
                                                      neig_idx(3)) % x),       &
                    geometry % dy(geometry % fine_map(neig_idx(1),neig_idx(2), &
                                                      neig_idx(3)) % y),       &
                    geometry % dz(geometry % fine_map(neig_idx(1),neig_idx(2), &
                                                      neig_idx(3)) % z)/)

          ! compute dtilde
          dtilde = (2.0_8*m % diffcof(g)*mn % diffcof(g)) /                    &
                (hxyzn(xyz_idx)*m % diffcof(g) + hxyz(xyz_idx)*mn % diffcof(g))

          ! compute leakage coefficient for neighbor
          jn = -dtilde

          ! get neighbor matrix index
          call indices_to_matrix(g,neig_idx(1),neig_idx(2),neig_idx(3),        &
         &                       neig_mat_idx)

          ! compute value and record to bank
          val = jn/hxyz(xyz_idx)

          ! record value in matrix
          call MatSetValue(this%M,irow,neig_mat_idx-1,val,                     &
          &                 INSERT_VALUES,ierr)

          this % row(kount) = irow + 1
          this % col(kount) = neig_mat_idx
          this % val(kount) = val
          kount = kount + 1

          ! compute leakage coefficient for target to cell
          jo(l) = shift_idx*dtilde

        else

          ! calculate dtilde
          dtilde = (2.0_8*m % diffcof(g)*(1.0_8-geometry % bc(l))) /           &
                   (4.0_8*m % diffcof(g)*(1.0_8+geometry % bc(l)) +                   &
                   (1.0_8-geometry % bc(l))*hxyz(xyz_idx))

          ! compute leakage coefficient for target next to boundary
          jo(l) = shift_idx*dtilde

        end if

      end do LEAK

      ! calate net leakage coefficient for target
      jnet = (jo(2) - jo(1))/hxyz(1) + (jo(4) - jo(3))/hxyz(2) +               &
     &       (jo(6) - jo(5))/hxyz(3)

      ! calculate loss of neutrons
      val = jnet + m % totalxs(g) - m % scattxs(g,g)

      ! record diagonal term
      call MatSetValue(this%M,irow,irow,val,INSERT_VALUES,ierr)
      this % row(kount) = irow + 1
      this % col(kount) = irow + 1 
      this % row_csr(2*irow+2) = kount ! sets diagnal csr index
      this % val(kount) = val
      kount = kount + 1

      ! begin loop over off diagonal in-scattering
      SCATTR: do h = 1,ng

        ! cycle though if h=g
        if (h == g) then
          cycle
        end if

        ! get neighbor matrix index
        call indices_to_matrix(h,i,j,k,scatt_mat_idx)

        ! record value in matrix (negate it)
        val = -m % scattxs(h,g)

        call MatSetValue(this%M,irow,scatt_mat_idx-1,val, INSERT_VALUES,ierr)
        this % row(kount) = irow + 1
        this % col(kount) = scatt_mat_idx
        this % val(kount) = val
        kount = kount + 1

      end do SCATTR

    end do ROWS 

    ! assemble matrix
    call MatAssemblyBegin(this%M,MAT_FLUSH_ASSEMBLY,ierr)
    call MatAssemblyEnd(this%M,MAT_FINAL_ASSEMBLY,ierr)
    call csr_sort_vectors(this)

    ! print out operator to file
    call print_M_operator(this)

  end subroutine build_loss_matrix

!===============================================================================
! INDICES_TO_MATRIX takes (x,y,z,g) indices and computes location in matrix 
!===============================================================================

  subroutine indices_to_matrix(g,i,j,k,matidx)

    use global, only: cmfd

    integer :: matidx         ! the index location in matrix
    integer :: i               ! current x index
    integer :: j               ! current y index
    integer :: k               ! current z index
    integer :: g               ! current group index

    ! compute index
    matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

  end subroutine indices_to_matrix

!===============================================================================
! MATRIX_TO_INDICES 
!===============================================================================

  subroutine matrix_to_indices(irow,g,i,j,k)

    use global, only: cmfd

    integer :: i                    ! iteration counter for x
    integer :: j                    ! iteration counter for y
    integer :: k                    ! iteration counter for z
    integer :: g                    ! iteration counter for groups
    integer :: irow                 ! iteration counter over row (0 reference)

    ! compute indices
    g = mod(irow,ng) + 1 
    i = mod(irow,ng*nx)/ng + 1
    j = mod(irow,ng*nx*ny)/(ng*nx)+ 1
    k = mod(irow,ng*nx*ny*nz)/(ng*nx*ny) + 1 

  end subroutine matrix_to_indices

!===============================================================================
! CSR_SORT_VECTORS
!===============================================================================

  subroutine csr_sort_vectors(this)

!---external references

    use math,  only: sort_csr

!---arguments

    type(loss_operator_type) :: this

!---local variables

    integer :: i
    integer :: j
    integer :: first
    integer :: last
    integer :: sr(4) = (/1,1,1,1/)
    integer :: sc(4) = (/5,2,9,1/)
    real(8) :: sv(4) = (/2.0_8,1.0_8,3.0_8,8.0_8/)

!---begin execution

    ! loop around row csr vector
    do i = 1, size(this % row_csr)/2 - 1

      ! get bounds
      first = this % row_csr(2*(i-1) + 1)
      last =  this %row_csr(2*i+1) - 1

      ! sort a row
      call sort_csr(this % row, this % col, this % val, first, last) 

    end do

  end subroutine csr_sort_vectors

!===============================================================================
! PRINT_M_OPERATOR 
!===============================================================================

  subroutine print_M_operator(this)

    type(loss_operator_type) :: this

    PetscViewer :: viewer

    ! write out matrix in binary file (debugging)
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'lossmat.bin' &
   &     ,FILE_MODE_WRITE,viewer,ierr)
    call MatView(this%M,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine print_M_operator

!==============================================================================
! DESTROY_M_OPERATOR
!==============================================================================

  subroutine destroy_M_operator(this)

    type(loss_operator_type) :: this

    ! deallocate matrix
    call MatDestroy(this%M,ierr)

    ! deallocate CSR objects
    if (allocated(this % row)) deallocate(this % row)
    if (allocated(this % col)) deallocate(this % col)
    if (allocated(this % val)) deallocate(this % val)
    if (allocated(this % row_csr)) deallocate(this % row_csr)

    ! deallocate other parameters
    if (allocated(this%d_nnz)) deallocate(this%d_nnz)
    if (allocated(this%o_nnz)) deallocate(this%o_nnz)

  end subroutine destroy_M_operator

end module loss_operator

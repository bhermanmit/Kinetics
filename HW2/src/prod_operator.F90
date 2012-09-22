module prod_operator

  implicit none
  private
  public :: init_F_operator,build_prod_matrix,destroy_F_operator

#include <finclude/petsc.h90>

    integer  :: nx   ! maximum number of x cells
    integer  :: ny   ! maximum number of y cells
    integer  :: nz   ! maximum number of z cells
    integer  :: ng   ! maximum number of groups
    integer  :: ierr ! petsc error code

  type, public :: prod_operator_type

    Mat      :: F    ! petsc matrix for neutronic prod operator
    integer  :: n    ! dimensions of matrix
    integer  :: nnz  ! max number of nonzeros
    integer  :: localn ! local size on proc
    integer, allocatable :: d_nnz(:) ! vector of diagonal preallocation
    integer, allocatable :: o_nnz(:) ! vector of off-diagonal preallocation
    integer, allocatable :: row(:)     ! vector of row indices
    integer, allocatable :: col(:)     ! vector of column indices
    integer, allocatable :: row_csr(:) ! row indexing for CSR
    real(8), allocatable :: val(:)     ! the value

  end type prod_operator_type

contains

!==============================================================================
! INIT_F_OPERATOR
!==============================================================================

  subroutine init_F_operator(this)

    type(prod_operator_type) :: this

    ! get indices
    call get_F_indices(this)

    ! get preallocation
    call preallocate_prod_matrix(this)

    ! set up M operator
    call MatCreateAIJ(PETSC_COMM_WORLD,this%localn,this%localn,PETSC_DECIDE,&
   & PETSC_DECIDE,PETSC_NULL_INTEGER,this%d_nnz,PETSC_NULL_INTEGER,this%o_nnz, &
   & this%F,ierr)
    call MatSetOption(this%F,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,ierr)
    call MatSetOption(this%F,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)

  end subroutine init_F_operator

!==============================================================================
! GET_F_INDICES
!==============================================================================

  subroutine get_F_indices(this)

    use global, only: geometry 

    type(prod_operator_type) :: this

    ! get maximum number of cells in each direction
    nx = geometry % nfx
    ny = geometry % nfy
    nz = geometry % nfz
    ng = geometry % nfg

    ! get number of nonzeros
    this%nnz = 7 + ng - 1

    ! calculate dimensions of matrix
    this%n = nx*ny*nz*ng

  end subroutine get_F_indices

!===============================================================================
! PREALLOCATE_PROD_MATRIX
!===============================================================================

  subroutine preallocate_prod_matrix(this)

    type(prod_operator_type) :: this

    integer :: rank          ! rank of processor
    integer :: sizen         ! number of procs
    integer :: i             ! iteration counter for x
    integer :: j             ! iteration counter for y
    integer :: k             ! iteration counter for z
    integer :: g             ! iteration counter for groups
    integer :: h             ! energy group when doing scattering
    integer :: n             ! the extent of the matrix
    integer :: irow          ! row counter
    integer :: row_start     ! index of local starting row
    integer :: row_end       ! index of local final row
    integer :: hmat_idx      ! index in matrix for energy group h

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

    ! begin loop around local rows
    ROWS: do irow = row_start,row_end

      ! initialize counters 
      this%d_nnz(irow) = 1 ! already add in matrix diagonal
      this%o_nnz(irow) = 0

      ! get location indices
      call matrix_to_indices(irow,g,i,j,k)

      ! begin loop over off diagonal in-scattering
      NFISS: do h = 1,ng

        ! cycle though if h=g
        if (h == g) then
          cycle
        end if

        ! get neighbor matrix index
        call indices_to_matrix(h,i,j,k,hmat_idx)

        ! record nonzero
        if (((hmat_idx-1) >= row_start) .and.                        &
       &   ((hmat_idx-1) <= row_end)) then
          this%d_nnz(irow) = this%d_nnz(irow) + 1
        else
          this%o_nnz(irow) = this%o_nnz(irow) + 1
        end if

      end do NFISS

    end do ROWS

    ! allocate CSR objects
    allocate(this % row(sum(this % d_nnz) + sum(this % o_nnz)))
    allocate(this % col(sum(this % d_nnz) + sum(this % o_nnz)))
    allocate(this % val(sum(this % d_nnz) + sum(this % o_nnz)))
    allocate(this % row_csr(row_end - row_start + 2))

  end subroutine preallocate_prod_matrix

!===============================================================================
! BUILD_PROD_MATRIX creates the matrix representing loss of neutrons
!===============================================================================

  subroutine build_prod_matrix(this)

    use global,           only: geometry, material 
    use material_header,  only: material_type

    type(prod_operator_type) :: this

    integer :: i                  ! iteration counter for x
    integer :: j                  ! iteration counter for y
    integer :: k                  ! iteration counter for z
    integer :: g                  ! iteration counter for groups
    integer :: l                  ! iteration counter for leakages
    integer :: h                  ! energy group when doing scattering
    integer :: gmat_idx           ! index in matrix for energy group g
    integer :: hmat_idx           ! index in matrix for energy group h
    integer :: ierr               ! Petsc error code
    integer :: row_start          ! the first local row on the processor
    integer :: row_finish         ! the last local row on the processor
    integer :: irow               ! iteration counter over row
    integer :: kount              ! csr counter
    integer :: idx
    real(8) :: nfissxs            ! nufission cross section h-->g
    real(8) :: val                ! temporary variable for nfissxs
    type(material_type), pointer :: m

    ! get row bounds for this processor
    call MatGetOwnershipRange(this%F,row_start,row_finish,ierr)

    ! initialize csr counter
    kount = 1

    ! begin iteration loops
    ROWS: do irow = row_start,row_finish-1

      ! set csr value of row
      this % row_csr(irow+1) = kount

      ! get indices for that row
      call matrix_to_indices(irow,g,i,j,k)

      ! set material pointer
      idx = ceiling(real(irow+1)/real(ng))
      m => material(geometry % fmat_map(idx))

      ! loop around all other groups 
      NFISS: do h = 1,ng

        ! get matrix column location
        call indices_to_matrix(h,i,j,k,hmat_idx)

        ! reocrd value in matrix
        val = m % nfissxs(g,h)
        call MatSetValue(this%F,irow,hmat_idx-1,val,INSERT_VALUES,ierr)
        this % row(kount) = irow + 1
        this % col(kount) = hmat_idx 
        this % val(kount) = val
        kount = kount + 1

      end do NFISS

    end do ROWS 

    ! put in last row index
    this % row_csr(row_finish - row_start) = kount

    ! assemble matrix 
    call MatAssemblyBegin(this%F,MAT_FLUSH_ASSEMBLY,ierr)
    call MatAssemblyEnd(this%F,MAT_FINAL_ASSEMBLY,ierr)
    call csr_sort_vectors(this)

    ! print out operator to file
    call print_F_operator(this)

  end subroutine build_prod_matrix

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

    type(prod_operator_type) :: this

!---local variables

    integer :: i
    integer :: j
    integer :: first
    integer :: last

!---begin execution

    ! loop around row csr vector
    do i = 1, size(this % row_csr) - 1

      ! get bounds
      first = this % row_csr(i)
      last =  this %row_csr(i+1) - 1

      ! sort a row
      call sort_csr(this % row, this % col, this % val, first, last)

    end do

  end subroutine csr_sort_vectors

!===============================================================================
! PRINT_F_OPERATOR 
!===============================================================================

  subroutine print_F_operator(this)

    PetscViewer :: viewer

    type(prod_operator_type) :: this

    ! write out matrix in binary file (debugging)
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'prodmat.bin'&
   &     ,FILE_MODE_WRITE,viewer,ierr)
    call MatView(this%F,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine print_F_operator

!==============================================================================
! DESTROY_F_OPERATOR
!==============================================================================

  subroutine destroy_F_operator(this)

    type(prod_operator_type) :: this

    ! deallocate matrix
    call MatDestroy(this%F,ierr)

    ! deallocate other parameters
    if (allocated(this%d_nnz)) deallocate(this%d_nnz)
    if (allocated(this%o_nnz)) deallocate(this%o_nnz)

  end subroutine destroy_F_operator

end module prod_operator

module operator_header

!-module options

  implicit none

!-module external references

# include "finclude/petsc.h90"

!-defintions

  type, public :: operator_type

    Mat      :: oper   ! petsc matrix
    integer  :: n      ! dimensions of matrix
    integer  :: nnz    ! max number of nonzeros
    integer  :: localn ! local size on proc
    integer, allocatable :: d_nnz(:) ! vector of diagonal preallocation
    integer, allocatable :: o_nnz(:)   ! vector of off-diagonal preallocation
    integer, allocatable :: row(:)     ! vector of row indices
    integer, allocatable :: col(:)     ! vector of column indices
    integer, allocatable :: row_csr(:) ! row indexing for CSR
    integer, allocatable :: diag(:)    ! diagonal index in CSR
    real(8), allocatable :: val(:)     ! the value

  end type operator_type


end module operator_header

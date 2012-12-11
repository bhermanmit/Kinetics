module hdf5_interface

  use global, only: hdf5_err, hdf5_output_file

#ifdef HDF5
  use hdf5
  use h5lt
  use, intrinsic :: ISO_C_BINDING
#endif

  implicit none

#ifdef HDF5

  ! define array interface
  interface hdf5_make_array
    module procedure hdf5_make_array_int4_r1
    module procedure hdf5_make_array_real8_r1
  end interface hdf5_make_array

contains

!===============================================================================
! HDF5_INITIALIZE
!===============================================================================

  subroutine hdf5_initialize()

    call h5open_f(hdf5_err)

  end subroutine hdf5_initialize

!===============================================================================
! HDF5_FINALIZE
!===============================================================================

  subroutine hdf5_finalize()

    call h5close_f(hdf5_err)

  end subroutine hdf5_finalize

!===============================================================================
! HDF5_create_file
!===============================================================================

  subroutine hdf5_create_file(filename)

    character(*) :: filename

    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5_output_file, hdf5_err)

  end subroutine hdf5_create_file

!===============================================================================
! HDF5_close_file
!===============================================================================

  subroutine hdf5_close_file()

    call h5fclose_f(hdf5_output_file, hdf5_err)

  end subroutine hdf5_close_file

!===============================================================================
! HDF5_create_group
!===============================================================================

  subroutine hdf5_create_group(root, group, name)

   integer(HID_T) :: root
   integer(HID_T) :: group
   character(*)   :: name

   call h5gcreate_f(root, name, group, hdf5_err)

  end subroutine hdf5_create_group

!===============================================================================
! HDF5_close_group
!===============================================================================

  subroutine hdf5_close_group(group)

    integer(HID_T) :: group

    call h5gclose_f(group, hdf5_err)

  end subroutine hdf5_close_group

!===============================================================================
! HDF5_MAKE_INTEGER
!===============================================================================

  subroutine hdf5_make_integer(group, name, buffer)

    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    integer,        intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltmake_dataset_int_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_make_integer

!===============================================================================
! HDF5_MAKE_DOUBLE
!===============================================================================

  subroutine hdf5_make_double(group, name, buffer)

    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    real(8),        intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltmake_dataset_double_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_make_double

!===============================================================================
! HDF5_MAKE_ARRAY_INT4_R1
!===============================================================================

  subroutine hdf5_make_array_int4_r1(group, name, buffer, n)

    integer (HID_T), intent(in) :: group
    character(*),    intent(in) :: name
    integer,         intent(in) :: buffer(:)
    integer,         intent(in) :: n

    integer :: rank = 1
    integer(HSIZE_T) :: dims(1)

    dims = (/n/)

    call h5ltmake_dataset_int_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_make_array_int4_r1

!===============================================================================
! HDF5_MAKE_ARRAY_REAL8_R1
!===============================================================================

  subroutine hdf5_make_array_real8_r1(group, name, buffer, n)

    integer (HID_T), intent(in) :: group
    character(*),    intent(in) :: name
    real(8),         intent(in) :: buffer(:)
    integer,         intent(in) :: n

    integer :: rank = 1
    integer(HSIZE_T) :: dims(1)

    dims = (/n/)

    call h5ltmake_dataset_double_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_make_array_real8_r1

#endif

end module hdf5_interface

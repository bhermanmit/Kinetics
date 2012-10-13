module global

!-module external references

  use constants,        only: MAX_LINE_LEN
  use cmfd_header,      only: cmfd_type
  use geometry_header,  only: geometry_type
  use kinetics_header,  only: kinetics_type
  use material_header,  only: material_type
  use operator_header,  only: operator_type
  use pke_header,       only: pke_type
  use timing,           only: Timer

# ifdef HDF5
    use hdf5
# endif

!-module options

  implicit none
  save

!-main objects

  type(cmfd_type), target                  :: cmfd        ! holds result
  type(geometry_type), target              :: geometry    ! holds geometry info
  type(material_type), allocatable, target :: material(:) ! holds material info
  type(kinetics_type), allocatable, target :: kinetics(:) ! holds kinetics info
  type(kinetics_type), allocatable, target :: pke_shape_for(:) ! forward shape function data
  type(kinetics_type), allocatable, target :: pke_shape_adj(:) ! adjoint shape function data

!-material information

  integer :: n_materials

!-timing objects

  type(Timer) :: time_total  ! timer for whole calculation
  type(Timer) :: time_init   ! timer for initialization
  type(Timer) :: time_build  ! timer for mat building
  type(Timer) :: time_power  ! timer for power iteration
  type(Timer) :: time_inner  ! timer for inner iterations
  type(Timer) :: time_kine   ! timer for kinetics

!-petsc error code

  integer :: ierr

!-mpi parameters

  logical :: master = .false. ! am i master
  integer :: rank             ! rank of processor
  integer :: n_procs          ! number of processors
  integer :: mpi_err          ! error code

!-running options

  character(len=50) :: solver_type = "jacobi"
  character(len=50) :: guess = "flat"
  character(len=50) :: adjoint = "none"
  character(len=50) :: mode = "static"
  character(len=50) :: weight = "unity"
  logical           :: pke_run = .false. 
  integer           :: pke_grp

!-solver tolerances

  real(8) :: ktol = 1.e-8_8
  real(8) :: stol = 1.e-6_8
  real(8) :: itol = 1.e-10_8

!-hdf5

  integer :: hdf5_err
  integer(HID_T) :: hdf5_output_file

!-Message used in message/warning/fatal_error

  character(MAX_LINE_LEN) :: message

!-The verbosity controls how much information will be printed to the
!-screen and in logs

  integer :: verbosity = 7

!-operators

  type(operator_type), target :: loss     ! static calculation loss operator
  type(operator_type), target :: loss_adj ! static adjoint calculation loss operator
  type(operator_type), target :: prod     ! static calculation production operator
  type(operator_type), target :: prod_adj ! static adjoint calculation production operator
  type(operator_type), target :: kine     ! kinetics operator (implicit Euler)

!-number of kinetics mods

  integer :: n_kins
  integer :: n_pkes_for
  integer :: n_pkes_adj

!-time step info

  integer :: nt
  real(8) :: time
  real(8) :: dt

!-Point Kinetics object

  type(pke_type) :: pke
  type(pke_type) :: gpke

end module global

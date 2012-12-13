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

  type(cmfd_type)                          :: cmfd        ! holds result
  type(geometry_type), target              :: geometry    ! holds geometry info
  type(material_type), allocatable, target :: material(:) ! holds material info
  type(kinetics_type), allocatable, target :: kinetics(:) ! holds kinetics info

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

!-mpi parametesr

  logical :: master = .false. ! am i master
  integer :: rank             ! rank of processor
  integer :: n_procs          ! number of processors
  integer :: mpi_err          ! error code

!-running options

  character(len=50) :: solver_type = "jacobi"
  character(len=50) :: guess = "flat"
  character(len=50) :: adjoint = ""
  logical           :: run_kinetics = .false.

!-solver tolerances

  real(8) :: ktol = 1.e-8_8
  real(8) :: stol = 1.e-6_8
  real(8) :: itol = 1.e-10_8

!-hdf5

  integer :: hdf5_err
  integer(HID_T) :: hdf5_output_file
  integer(HID_T) :: time_group

!-Message used in message/warning/fatal_error

  character(MAX_LINE_LEN) :: message

!-The verbosity controls how much information will be printed to the
!-screen and in logs

  integer :: verbosity = 7

!-operators

  type(operator_type) :: loss  ! static calculation loss operator
  type(operator_type) :: prod  ! static calculation production operator
  type(operator_type) :: kine  ! kinetics operator (implicit Euler)

!-number of kinetics mods

  integer :: n_kins

!-time step info

  integer :: nt
  real(8) :: time
  real(8) :: dt

!-mode

  character(len=50) :: mode = "static"

!-point kinetics object

  type(pke_type) :: pke

!-variable time_step for rk

  logical :: var_ts = .true.

!-initial power, temperatures

  real(8) :: power = 1.0_8
  real(8) :: fuel_T
  real(8) :: cool_T

!-fuel properties

  real(8) :: fuel_a ! temperature reactivity coefficient
  real(8) :: fuel_c ! specific heat capacity
  real(8) :: fuel_m ! mass of fuel

!-coolant properties

  real(8) :: cool_a ! temperature reactivity coefficient
  real(8) :: cool_c ! specific heat capacity
  real(8) :: cool_m ! mass of coolant

!-other properaties of fuel-coolant system

  real(8) :: hA     ! effective heat transfer from fuel to coolant
  real(8) :: m_dot  ! mass flow rate of coolant
  real(8) :: Tin    ! inlet temperature of coolant
  real(8) :: P_frac ! fraction of power deposited directly to COOLANT

!-truncation error tolerance

  real(8) :: eps

end module global

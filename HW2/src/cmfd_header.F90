module cmfd_header

!-module options

  implicit none

!-module variables

  type cmfd_type

    ! eigenvector/eigenvalue from cmfd run
    real(8), allocatable :: phi(:)
    real(8) :: keff = 0.0_8

  end type cmfd_type

end module cmfd_header

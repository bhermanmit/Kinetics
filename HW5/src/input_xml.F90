module input_xml

!-module options

  implicit none
  private
  public :: read_input_xml

contains

!===============================================================================
! READ_INPUT reads the CMFD input file and organizes it into a data structure
!===============================================================================

  subroutine read_input_xml()

!---external references

    use constants,        only: MAX_FILE_LEN
    use error,            only: fatal_error
    use geometry_header,  only: allocate_geometry_type 
    use global,           only: material, geometry, message, n_materials,      &
                                solver_type, ktol, stol, itol, guess, adjoint, &
                                run_kinetics, nt, time, dt, kinetics, n_kins,  &
                                poi, poi_tol
    use kinetics_header,  only: allocate_kinetics_type
    use material_header,  only: material_type, allocate_material_type
    use output,           only: write_message
    use xml_data_input_t

!---local variables

    type(material_type), pointer :: m => null()
    logical :: file_exists
    character(MAX_FILE_LEN) :: filename
    integer :: i

!---begin execution

    ! display output message
    message = "Reading input XML file..."
    call write_message(5)

    ! check if settings.xml exists
    filename = "input.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
       message = "Input XML file '" // trim(filename) // "' does not exist!"
       call fatal_error()
    end if

    ! read in input file
    call read_xml_file_input_t(filename)

    ! get solver type
    solver_type = trim(solver_)

    ! get starting guess
    guess = trim(guess_)

    ! get adjoint
    adjoint = trim(adjoint_)

    ! get solver tolerances
    ktol = ktol_
    stol = stol_
    itol = itol_

    ! kinetics run
    run_kinetics = run_kinetics_

    ! get time info if running kinetics
    if (run_kinetics) then
      nt = nt_
      time = time_
      dt = time / dble(nt)
    end if

    ! read in geometry indices
    geometry % ncx = geometry_ % nx
    geometry % ncy = geometry_ % ny
    geometry % ncz = geometry_ % nz
    geometry % ncg = geometry_ % ng

    ! set fine grid
    geometry % nfx = sum(geometry_ % nnx)
    geometry % nfy = sum(geometry_ % nny)
    geometry % nfz = sum(geometry_ % nnz)
    geometry % nfg = geometry_ % ng

    ! allocate geometry object
    call allocate_geometry_type(geometry)

    ! read in coarse material map first check dimensions
    if (geometry % ncx*geometry % ncy*geometry % ncz /=                        &
        size(geometry_ % mat)) then
      message = "Map dimensions do not match indices!"
      call fatal_error()
    end if 
    geometry % mat_map = reshape(geometry_ % mat, (/geometry % ncx,            &
                                                    geometry % ncy,            &
                                                    geometry % ncz/))

    ! read in coarse region map first check dimensions
    if (geometry % ncx*geometry % ncy*geometry % ncz /=                        &
        size(geometry_ % reg)) then
      message = "Reg dimensions do not match indices!"
      call fatal_error()
    end if
    geometry % reg_map = reshape(geometry_ % reg, (/geometry % ncx,            &
                                                    geometry % ncy,            &
                                                    geometry % ncz/))

    ! number of regions
    geometry % n_regs = maxval(geometry % reg_map)

    ! read in grid info
    geometry % xgrid = geometry_ % xgrid
    geometry % ygrid = geometry_ % ygrid
    geometry % zgrid = geometry_ % zgrid
    geometry % vol = sum(geometry % xgrid) + sum(geometry % ygrid) +           &
                     sum(geometry % zgrid)
   
    ! read in fine mesh per coarse mesh 
    geometry % nnx = geometry_ % nnx 
    geometry % nny = geometry_ % nny
    geometry % nnz = geometry_ % nnz

    ! read in boundary conditions
    geometry % bc = geometry_ % bc

    ! get size of materials
    n_materials = size(material_)

    ! allocate material
    allocate(material(n_materials))

    ! begin loop and read in material information
    do i = 1, n_materials

      ! set material pointer
      m => material(i)

      ! check that uid matches loop
      if (i /= material_(i) % uid) then
        message = "Materials need to go in ascending order in order from 1" // &
                  " to number of materials"
        call fatal_error()
      end if

      ! allocate material
      call allocate_material_type(m, geometry % ncg)

      ! save materials if they are associated
      if (associated(material_(i) % absxs)) then
        m % abs_based = .true.
        m % absorxs = material_(i) % absxs
      end if

      if (associated(material_(i) % remxs)) then
        m % rem_based = .true.
        m % removxs = material_(i) % remxs
      end if

      if (associated(material_(i) % totalxs)) then
        m % totalxs = material_(i) % totalxs
      else
        if (.not. m % abs_based .and. .not. m % rem_based) then
          message ="Need to specify absorption or total xs!"
          call fatal_error()
        end if
      end if

      if (.not.associated(material_(i) % diffcoef)) then
        message = "Need to specify diffusion coefficients!"
        call fatal_error()
      else
        m % diffcof = material_(i) % diffcoef
      end if

      if (.not.associated(material_(i) % scattxs)) then
        message = "Need to specify scattering xs!"
        call fatal_error()
      else
        m % scattxs = reshape(material_(i) % scattxs,(/geometry % ncg,         &
                                                       geometry % ncg/))
      end if

      if (associated(material_(i) % chi)) then
        m % chi = material_(i) % chi
        m % chi_based = .true.
      end if

      if (run_kinetics) then
        if (.not.associated(material_(i) % chip)) then
          m % chip = m % chi
        else if (m % chi_based) then
          m % chip = material_(i) % chip
        else
          message = "Need to set a chi for kinetics!"
          call fatal_error()
        end if

        if (.not.associated(material_(i) % chid)) then
          m % chid = m % chi
        else if (m % chi_based) then
          m % chid = material_(i) % chid
        else
          message = "Need to set a chi for Kinetics!"
        end if
      end if

      if (associated(material_(i) % nfissxs) .and. m % chi_based) then
        m % fissvec = material_(i) % nfissxs
      else
        m % nfissxs = reshape(material_(i) % nfissxs,(/geometry % ncg,         &
                                                       geometry % ncg/))
      end if

      if (material_(i) % buckling > 1.e-20_8) then
        m % use_buckling = .true.
        m % buckling = material_(i) % buckling
        if (.not. m % abs_based) then
          message = "Cant use buckling unless absxs is specified!"
          call fatal_error()
        end if
      else
        m % buckling = 0.0_8
      end if

      if (run_kinetics) then
        allocate(m % kinrem(geometry % ncg))
        allocate(m % kinfis(geometry % ncg))
      end if

    end do

    ! read in kinetics mods
    n_kins = size(kinetics_)
    allocate(kinetics(n_kins))
    do i = 1, n_kins

      ! allocate kineics 
      call allocate_kinetics_type(kinetics(i),size(kinetics_(i)%time))

      ! set values
      kinetics(i) % mat_id  = kinetics_(i) % mat_id
      kinetics(i) % xs_id   = trim(kinetics_(i) % xs_id)
      kinetics(i) % g       = kinetics_(i) % g
      kinetics(i) % h       = kinetics_(i) % h
      kinetics(i) % time    = kinetics_(i) % time
      kinetics(i) % val     = kinetics_(i) % value

    end do

    ! read in point of interest
    poi = poi_
    poi_tol = poi_tol_

  end subroutine read_input_xml

end module input_xml 

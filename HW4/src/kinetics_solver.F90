module kinetics_solver 

!-module references

  use global,  only: mpi_err

!-module options

  implicit none
  private
  public :: kinetics_execute, change_data, compute_pkes, compute_gpkes 

!-module external references

# include "finclude/petsc.h90"

!-module variables

  KSP :: ksp
  Vec :: phi
  Vec :: b
  PC  :: pc

contains

!===============================================================================
! KINETICS_EXECUTE
!===============================================================================

  subroutine kinetics_execute(inner_solver)

!---external references

    use global,  only: kine

!---arguments

    external :: inner_solver

!---begin execution

    ! initialize the time 0 data
    call init_data()

    ! initialize solver
    call init_solver()

    ! run through the kinetics iterations
    call execute_kinetics_iter(inner_solver)
 
  end subroutine kinetics_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

!---external references

    use constants,          only: ONE
    use global,             only: cmfd, geometry, material, kine, prod, nt
    use kinetics_operator,  only: init_K_operator, build_kinetics_matrix
    use math,               only: csr_matvec_mult

!---local variables

    real(8) :: pow

!---begin execution

    ! normalize initial power to unity and set initial power
!   pow = sum(csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val/cmfd%keff,        &
!             cmfd%phi,prod%n))
!   cmfd % phi = cmfd % phi * ONE / pow

    ! compute steady state precursors
    call compute_initial_precursors()

    ! set up matrices
    call init_K_operator(kine)
    call build_kinetics_matrix(kine)
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,kine%n,kine%n,kine%row_csr,&
                                   kine%col,kine%val,kine%oper,mpi_err)

    ! allocate output arrays
    if (.not.allocated(cmfd%core_power)) allocate(cmfd % core_power(nt+1))
    if (.not.allocated(cmfd%time)) allocate(cmfd % time(nt+1))

    ! bank static keff
    cmfd % kcrit = cmfd % keff

  end subroutine init_data

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver()

!---external references

    use global,  only: itol

!---begin execution

    ! set up krylov solver
    call KSPCreate(PETSC_COMM_SELF,ksp,mpi_err)
    call KSPSetTolerances(ksp,itol,PETSC_DEFAULT_DOUBLE_PRECISION,     &
   &                      PETSC_DEFAULT_DOUBLE_PRECISION,                      &
   &                      PETSC_DEFAULT_INTEGER,mpi_err)
    call KSPSetType(ksp,KSPGMRES,mpi_err)
    call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,mpi_err)
    call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,mpi_err)
    call KSPGetPC(ksp,pc,mpi_err)
    call PCSetType(pc,PCILU,mpi_err)
    call PCFactorSetLevels(pc,5,mpi_err)
    call KSPSetFromOptions(ksp,mpi_err)

  end subroutine init_solver

!===============================================================================
! COMPUTE_INITIAL_PRECURSORS
!===============================================================================

  subroutine compute_initial_precursors()

!---external references

    use constants,        only: NUM_PRECS, beta, lambda
    use global,           only: geometry, material, cmfd
    use material_header,  only: material_type

!---local variables

    integer :: i
    integer :: j
    real(8) :: fiss
    type(material_type), pointer :: m => null()

!---begin execution

    ! allocate precurors
    if (.not.allocated(cmfd % C)) allocate(cmfd%C(geometry%nfx*geometry%nfy *  &
                                                  geometry%nfz,NUM_PRECS))

    ! begin loop around space
    do j = 1,size(cmfd % C,1)

      ! compute fission reaction rate over all groups
      m => material(geometry % fmat_map(j))
      fiss = sum(m % fissvec * cmfd % phi(geometry % nfg * (j - 1) + 1: &
                                                   j * geometry % nfg))
                                                   
      ! begin loop around precursor groups
      do i = 1, NUM_PRECS

        cmfd % C(j,i) = beta(i)/(lambda(i) * cmfd % keff) * fiss 

      end do

    end do

  end subroutine compute_initial_precursors

!===============================================================================
! EXECUTE_KINETICS_ITER  in the main kinetics iteration routine 
!                         for the cmfd calculation
!===============================================================================

  subroutine execute_kinetics_iter(inner_solver)

!---external references

    use constants,          only: ZERO, ONE, beta, NUM_PRECS
    use error,              only: fatal_error
    use global,             only: nt, dt, kine, cmfd, geometry, material,      &
                                  itol, prod, message, time_kine, time_inner,  &
                                  loss, prod, kinetics
    use kinetics_operator,  only: build_kinetics_matrix
    use loss_operator,      only: init_M_operator, build_loss_matrix
    use prod_operator,      only: init_F_operator, build_prod_matrix
    use math,               only: csr_matvec_mult
    use timing,             only: timer_start, timer_stop, timer_reset

!---arguments

    external :: inner_solver

!---local variables

    integer :: i
    integer :: iters
    integer :: n
    integer :: nz
    integer :: ng
    real(8) :: curr_time
    real(8) :: pow
    real(8), allocatable :: rhs(:), temp1(:), temp2(:)
    real(8) :: rho_num, rho_den

!---begin execution

    ! allocate RHS
    n = size(cmfd % phi)
    nz = size(kine % col)
    allocate(rhs(n))

    ! get number of groups
    ng = geometry % nfg

    ! set initial output
    cmfd % time(1) = ZERO
    cmfd % core_power(1) = ONE

    ! associate petsc vectors
    call VecCreateSeqWithArray(PETSC_COMM_WORLD,1,kine%n,cmfd%phi,phi,mpi_err)
    call VecCreateSeqWithArray(PETSC_COMM_WORLD,1,kine%n,rhs,b,mpi_err)

    ! set up krylov info
    call KSPSetOperators(ksp, kine%oper, kine%oper, SAME_NONZERO_PATTERN, mpi_err)
    call KSPSetUp(ksp,mpi_err)

    ! initialize operators
    call init_M_operator(loss)
    call init_F_operator(prod)

    ! allocate vectors
    allocate(temp1(loss % n))
    allocate(temp2(prod % n))
    if(.not.allocated(cmfd % rho)) allocate(cmfd % rho(nt))
    if(.not.allocated(cmfd % pnl)) allocate(cmfd % pnl(nt))
    if(.not.allocated(cmfd % prompt)) allocate(cmfd % prompt(ng,ng,nt))
    if(.not.allocated(cmfd % delay)) allocate(cmfd % delay(NUM_PRECS,ng,ng,nt))
    if(.not.allocated(cmfd % vel)) allocate(cmfd % vel(ng,nt))

    ! build petsc matrices
    call build_loss_matrix(loss,'')
    call build_prod_matrix(prod,'')
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,loss%n,loss%n,loss%row_csr,&
                                   loss%col,loss%val,loss%oper,mpi_err)
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,prod%n,prod%n,prod%row_csr,&
                                   prod%col,prod%val,prod%oper,mpi_err)

    ! begin loop around time
    do i = 1, nt

      ! compute current time
      curr_time = dble(i)*dt
    
      ! change the material via kinetics mods
      call change_data(kinetics,curr_time)
      call  change_kinetics()

      ! build matrices
      call build_kinetics_matrix(kine)

      ! build right hand side
      call build_rhs(rhs,n)

      ! calculate preconditioner (ILU)
      call PCFactorGetMatrix(pc,kine%oper,mpi_err)

      ! solve matrices
      call KSPSolve(ksp,b,phi,mpi_err)

     ! compute end of time step precursor concentration
     call compute_final_precursors()

     ! compute power
     pow = sum(csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val/cmfd%keff,       &
              cmfd%phi,prod%n))
     write(*,*) 'Step:', i,' /',nt,' POWER:',pow
     write(779,*) pow
     cmfd % time(i+1) = curr_time
     cmfd % core_power(i+1) = pow

     call compute_pkes(i,temp1,temp2)
     call compute_gpkes(i)

   end do

   ! deallocate RHS
   deallocate(rhs)
   deallocate(temp1)
   deallocate(temp2)

  end subroutine execute_kinetics_iter 

!===============================================================================
! COMPUTE_FINAL_PRECURSORS
!===============================================================================

  subroutine compute_final_precursors()

!---external references

    use constants,        only: NUM_PRECS, beta, lambda, ONE
    use global,           only: geometry, material, cmfd, dt
    use material_header,  only: material_type

!---local variables

    integer :: i
    integer :: j
    real(8) :: fiss
    type(material_type), pointer :: m => null()

!---begin execution

    ! begin loop around space
    do j = 1,size(cmfd%C,1)

      ! compute fission reaction rate over all groups
      m => material(geometry % fmat_map(j))
      fiss = sum(m % fissvec * cmfd % phi(geometry % nfg * (j - 1) + 1: &
                                                   j * geometry % nfg))

      ! begin loop around precursor groups
      do i = 1, NUM_PRECS

        cmfd % C(j,i) = beta(i)*dt/((ONE + dt*lambda(i)) * cmfd % keff) * fiss +  &
                        cmfd % C(j,i)/(ONE + dt*lambda(i))

      end do

    end do

  end subroutine compute_final_precursors

!===============================================================================
! CHANGE_DATA
!===============================================================================

  subroutine change_data(kinetics,t)

!---external references

    use constants,        only: ONE, vel, beta, lambda
    use error,            only: fatal_error
    use global,           only: material, n_kins, message, dt,       &
                                n_materials, cmfd
    use kinetics_header,  only: kinetics_type
    use material_header,  only: material_type

!---arguments

    real(8) :: t
    type(kinetics_type), target :: kinetics(:)

!---local variables

    integer :: i
    real(8) :: val
    type(kinetics_type), pointer :: k => null()
    type(material_type), pointer :: m => null()

!---begin execution

    ! begin loop around kinetics mods
    do i = 1, n_kins

      ! point to object
      k => kinetics(i)

      ! point to material object
      m => material(k % mat_id)

      ! check to shift index
      if (t > k % time(k % idxt+1)) k % idxt = k % idxt + 1

      ! interpolate value
      val = k % val(k%idxt) + (k % val(k%idxt+1) - k % val(k%idxt)) /          &
            (k % time(k%idxt+1) - k % time(k%idxt)) * (t - k % time(k%idxt))

      ! begin case structure to replace value
      select case (trim(k % xs_id))

        case ('absxs')

          ! remove scattering component from total xs
          m % totalxs(k % g) = m % totalxs(k % g) - sum(m % scattxs(:,k % g))

          ! change absorption
          m % absorxs(k % g) = val
          m % totalxs(k % g) = val

          ! re-add back in scattering
          m % totalxs(k % g) = m % totalxs(k % g) + sum(m % scattxs(:,k % g))

        case DEFAULT

          message = 'Kineics modification not supported!'
          call fatal_error()

      end select

    end do

  end subroutine change_data

!==============================================================================
! CHANGE_KINETICS
!==============================================================================

  subroutine change_kinetics()

!---external references

    use constants,        only: ONE, vel, beta, lambda
    use global,           only: material, n_kins, message, dt,       &
                                n_materials, cmfd
    use material_header,  only: material_type

!---local variables

    integer :: i
    type(material_type), pointer :: m => null()

!---begin execution

    ! loop through materials and compute kinetics factors
    do i = 1, n_materials

      ! point to materal
      m => material(i)

      ! compute additive factor that goes on to the removal xs
      m % kinrem = ONE/(vel * dt)

      ! compute multiplicative factor that goes on fission xs
      m % kinfis = (ONE - sum(beta)) * m % chip / cmfd % keff +                &
      m % chid * sum((lambda * beta * dt)/((ONE + lambda*dt) * cmfd % keff))

    end do

  end subroutine change_kinetics

!==============================================================================
! BUILD_RHS
!==============================================================================

  subroutine build_rhs(rhs,n)

!---external references

    use constants,        only: vel, lambda, ONE
    use global,           only: dt, material, geometry, cmfd
    use material_header,  only: material_type

!---arguments

    integer :: n
    real(8) :: rhs(n)

!---local variables

    integer :: irow
    integer :: idx
    integer :: g
    type(material_type), pointer :: m => null()

!---begin execution

    ! begin loop around rows
    do irow = 1, n

      ! get index
      idx = ceiling(real(irow)/real(geometry%nfg))
      g = geometry % nfg - mod(irow,geometry%nfg)

      ! set material pointer
      m => material(geometry % fmat_map(idx))

      ! compute the precursor
      rhs(irow) = cmfd % phi(irow)/(vel(g)*dt) + sum((m%chid(g)*lambda*        &
                                            cmfd % C(idx,:))/(ONE + lambda*dt))
    end do

  end subroutine build_rhs

!===============================================================================
! COMPUTE PKES
!===============================================================================

  subroutine compute_pkes(i,temp1,temp2)

!---references

    use constants,          only: beta,vel
    use global,             only: cmfd, loss, prod
    use math,               only: csr_matvec_mult
    use loss_operator,      only: init_M_operator, build_loss_matrix
    use prod_operator,      only: init_F_operator, build_prod_matrix

!---arguments

    integer :: i
    real(8) :: temp1(:), temp2(:)

!---local variables

    integer :: irow
    real(8) :: rho_num, rho_den, pnl_num, pnl_den 

!---begin execution

    ! build operators
    call build_loss_matrix(loss,'')
    call build_prod_matrix(prod,'')

    ! we need F - M, store it back in M
    call MatAXPY(loss%oper,-1.0_8/cmfd%kcrit,prod%oper, SUBSET_NONZERO_PATTERN, mpi_err)
    call MatScale(loss%oper,-1.0_8,mpi_err)

    ! perform reactivity numerator operator multiplication
    temp1 =  csr_matvec_mult(loss%row_csr+1,loss%col+1,loss%val,cmfd%phi,loss%n)
    temp2 =  csr_matvec_mult(prod%row_csr+1,prod%col+1,prod%val,cmfd%phi,prod%n)

    ! multiply by volume
    call multiply_volume(temp1,loss%n)
    call multiply_volume(temp2,loss%n)

    ! perform weighting (unity will occur if adjoint flux is 1)
    rho_num = dot_product(cmfd%phi_adj,temp1) 
    rho_den = dot_product(cmfd%phi_adj,temp2)/cmfd % kcrit

    ! calc reactivity
    cmfd % rho(i) = rho_num/rho_den/sum(beta)

    ! loop around forward shape function and multiply by inverse velocity
    do irow = 1,loss%n
      if (mod(irow,2) == 0) then
        temp1(irow) = 1.0_8/vel(2)*cmfd%phi(irow)*cmfd%phi_adj(irow)
      else
        temp1(irow) = 1.0_8/vel(1)*cmfd%phi(irow)*cmfd%phi_adj(irow)
      end if
    end do

    ! get rate by multiplying by volume
    call multiply_volume(temp1,loss%n)

    ! compute pnl
    pnl_num = sum(temp1)
    pnl_den = rho_den 
    cmfd % pnl(i) = pnl_num/pnl_den 
    write(834,*) cmfd % rho(i), cmfd % pnl(i)

  end subroutine compute_pkes

!===============================================================================
! COMPUTE GPKES
!===============================================================================

  subroutine compute_gpkes(t)

!---references

    use constants,          only: beta, vel, ONE, ZERO, NUM_PRECS
    use global,             only: cmfd, loss, prod, geometry, material
    use loss_operator,      only: build_loss_matrix
    use prod_operator,      only: build_prod_matrix
    use material_header,    only: material_type

!---arguments

    integer :: t

!---local variables

    integer :: irow, icol, igrp, iprec
    integer :: h, g, i, j, k, ni, nj, nk
    integer :: first, last
    integer :: idx, idxn, matidx
    real(8) :: vol, voln
    real(8) :: norm
    type(material_type), pointer :: m

!---begin execution

    ! build operators
    call build_loss_matrix(loss,'')
    call build_prod_matrix(prod,'')

    ! zero out values
    cmfd % prompt(:,:,t) = ZERO
    cmfd % delay(:,:,:,t) = ZERO
    cmfd % vel(:,t) = ZERO 

    ! compute normalization factor
    norm = cmfd % factor / sum(cmfd % phi)

    ! begin loop around column in operator
    ROWS: do irow = 1, size(loss % row_csr) - 1

      ! set bounds
      first = loss % row_csr(irow) + 1
      last = loss % row_csr(irow+1)

      ! get indices of this location
      call matrix_to_indices(irow-1,g,i,j,k,geometry%nfg, geometry%nfx,  &
                             geometry%nfy, geometry%nfz)

      ! get index to get fine mesh indices
      idx = ceiling(real(irow)/real(geometry % nfg))

      ! get volume
      vol = geometry % fdx_map(idx)*geometry % fdy_map(idx)*geometry % fdz_map(idx)

      ! get material pointer
      m => material(geometry % fmat_map(idx)) 

      ! loop around columns in row
      COLS: do icol = first, last

        ! get indices of this location
        call matrix_to_indices(loss%col(icol),h,ni,nj,nk,geometry%nfg, geometry%nfx,  &
                               geometry%nfy, geometry%nfz)

        ! get fine map index
        idxn = ceiling(real(loss%col(icol)+1)/real(geometry % nfg))

        ! get volume
        voln = geometry % fdx_map(idxn)*geometry % fdy_map(idxn)*geometry % fdz_map(idxn)

        ! take value multiply volume and add to appropriate location
        cmfd % prompt(g,h,t) = cmfd % prompt(g,h,t) - cmfd%phi_adj(loss%col(icol)+1)* &
                            loss%val(icol)*cmfd%phi(loss%col(icol)+1)*norm*voln

      end do COLS

      ! loop around groups for fission
      GRPS: do igrp = 1, geometry % nfg

        ! get matrix index for this group
        call indices_to_matrix(igrp,i,j,k,matidx,geometry%nfg,geometry%nfx,geometry%nfy,geometry%nfz) 

        ! put in prompt fission value
        cmfd % prompt(g,igrp,t) = cmfd % prompt(g,igrp,t) + cmfd%phi_adj(matidx)*    &
          (ONE - sum(beta))*m%chip(g)/cmfd%kcrit*m%fissvec(igrp)*cmfd%phi(matidx)*norm*vol

        ! put in delayed fission value
        do iprec = 1, NUM_PRECS
          cmfd % delay(iprec,g,igrp,t) = cmfd % delay(iprec,g,igrp,t) + cmfd % phi_adj(matidx)*    &
            beta(iprec)*m%chid(g)/cmfd%kcrit*m%fissvec(igrp)*cmfd%phi(matidx)*norm*vol
        end do

      end do GRPS
        if(g == 1) write(877,*) cmfd % phi(irow)
        if(g == 2) write(878,*) cmfd % phi(irow)


      cmfd % vel(g,t) = cmfd % vel(g,t) + cmfd%phi_adj(irow)*cmfd%phi(irow)*norm*vol

    end do ROWS

    ! divide velocity parameter
    cmfd % vel(:,t) = vel / cmfd % vel(:,t)

  end subroutine compute_gpkes

!===============================================================================
! MULTIPLY_VOLUME
!===============================================================================

  subroutine multiply_volume(vec,n)

!---external references

    use global,  only: geometry

!---arguments

    integer :: n
    real(8) :: vec(n)

!---local variables

    integer :: irow
    integer :: idx
    real(8) :: vol

!---begin execution

    ! begin loop around rows
    do irow = 1, n

      ! compute index
      idx = ceiling(real(irow)/real(geometry%nfg))

      ! get region number
      vol = geometry % fvol_map(idx)

      ! sum nodal power
      vec(irow) = vec(irow) * vol

    end do

  end subroutine multiply_volume

!===============================================================================
! MATRIX_TO_INDICES 
!===============================================================================

  subroutine matrix_to_indices(irow,g,i,j,k,ng,nx,ny,nz)

    integer :: i                    ! iteration counter for x
    integer :: j                    ! iteration counter for y
    integer :: k                    ! iteration counter for z
    integer :: g                    ! iteration counter for groups
    integer :: ng
    integer :: nx
    integer :: ny
    integer :: nz
    integer, intent(in) :: irow     ! iteration counter over row (0 reference)

    ! compute indices
    g = mod(irow,ng) + 1
    i = mod(irow,ng*nx)/ng + 1
    j = mod(irow,ng*nx*ny)/(ng*nx)+ 1
    k = mod(irow,ng*nx*ny*nz)/(ng*nx*ny) + 1

  end subroutine matrix_to_indices

!===============================================================================
! INDICES_TO_MATRIX takes (x,y,z,g) indices and computes location in matrix 
!===============================================================================

  subroutine indices_to_matrix(g,i,j,k,matidx,ng,nx,ny,nz)

    use global, only: cmfd

    integer :: matidx         ! the index location in matrix
    integer :: i               ! current x index
    integer :: j               ! current y index
    integer :: k               ! current z index
    integer :: g               ! current group index
    integer :: ng
    integer :: nx
    integer :: ny
    integer :: nz

    ! compute index
    matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

  end subroutine indices_to_matrix

!==============================================================================
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize data objects

  end subroutine finalize

end module kinetics_solver 

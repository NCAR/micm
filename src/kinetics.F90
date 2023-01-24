! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_kinetics module

!> The kinetics_t type and related functions
!!
!! \todo Consider rewriting micm_kinetics to follow musica style conventions
module micm_kinetics

  use micm_environment,                only : environment_t
  use musica_constants,                only : musica_dk, musica_ik
  use constants,                       only : length, VLEN, STREAM0

  implicit none

  private
  public :: kinetics_t

! rateConst are computed at the beginning of the
!   chemistry_box_solver time step.
!   They are not computed for internal steps of the
!   box-model time-step advancer
! rateConst will be thread-safe memory provided elsewhere.
! rate_constant_store will be an accessor to memory
! For now, it is allocated here. It is not thread safe

type kinetics_t
  private
  integer :: nReact
  integer :: nSpecies
  integer, allocatable         :: Pivot(:,:)
  real(musica_dk), allocatable, public :: rateConst(:,:)
  real(musica_dk), allocatable, public :: MBOdeJac(:,:)      ! ODE solver jacobian
  real(musica_dk), allocatable, public :: chemJac(:,:)       ! chemistry forcing jacobian
  real(musica_dk), allocatable, public :: rates(:,:)         ! rates of reactions
  type(environment_t), allocatable, public :: environment(:)
contains
  procedure, public :: species_names
  procedure, public :: reaction_names
  procedure, public :: photolysis_reaction_names
  procedure, public :: update
  procedure, public :: rateConst_print
  procedure, public :: LinFactor
  procedure, public :: LinSolve
  procedure, public :: calc_force
  procedure, public :: calc_reaction_rates
  procedure, public :: calc_reaction_rate_constants
  procedure, public :: calc_dforce_dy
  procedure, public :: dForcedyxForce
!  procedure, private :: LinFactor
  final             :: DasEnder
end type kinetics_t

  !> Constructor for kinetics_t objects
  interface kinetics_t
    module procedure :: constructor
  end interface kinetics_t

contains

  !---------------------------
  ! Chemical species names
  !---------------------------
  subroutine species_names( this, names )

    use factor_solve_utilities, only : number_of_species
    use kinetics_utilities,     only : get_names => species_names
    use musica_string,          only : string_t

    class(kinetics_t), intent(in) :: this
    type(string_t), allocatable, intent(inout) :: names(:)

    integer :: i_spec
    character(len=128) :: raw_names(number_of_species)

    raw_names = get_names()
    if( allocated( names ) ) deallocate( names )
    allocate( names( number_of_species ) )
    do i_spec = 1, number_of_species
      names(i_spec) = raw_names(i_spec)
    end do

  end subroutine species_names

  !---------------------------
  ! Reaction labels
  !---------------------------
  subroutine reaction_names( this, names )

    use kinetics_utilities, only : number_of_reactions, &
                                   get_names => reaction_names
    use musica_string,      only : string_t

    class(kinetics_t), intent(in) :: this
    type(string_t), allocatable, intent(inout) :: names(:)

    integer :: i_rxn
    character(len=128) :: raw_names(number_of_reactions)

    raw_names = get_names()
    if( allocated( names ) ) deallocate( names )
    allocate( names( number_of_reactions ) )
    do i_rxn = 1, number_of_reactions
      names(i_rxn) = raw_names(i_rxn)
    end do

  end subroutine reaction_names

  !---------------------------
  ! Photolysis reaction names
  !---------------------------
  subroutine photolysis_reaction_names( this, names )

    use kinetics_utilities, only : number_of_photolysis_reactions, &
                                   get_names => photolysis_names
    use musica_string,      only : string_t

    class(kinetics_t), intent(in) :: this
    type(string_t), allocatable, intent(inout) :: names(:)

    integer :: i_rxn
    character(len=128) :: raw_names(number_of_photolysis_reactions)

    raw_names = get_names()
    if( allocated( names ) ) deallocate( names )
    allocate( names( number_of_photolysis_reactions ) )
    do i_rxn = 1, number_of_photolysis_reactions
      names(i_rxn) = raw_names(i_rxn)
    end do

  end subroutine photolysis_reaction_names

  !---------------------------
  ! Compute time rate of change of each molecule (vmr) given reaction rates
  !---------------------------
  subroutine calc_force( this, vmr, force )

    use kinetics_utilities,  only :  p_force
    use factor_solve_utilities, only: number_of_species

    class(kinetics_t) :: this
    real(musica_dk), intent(in)  ::  vmr(length,number_of_species)      ! volume mixing ratios of each component in order
    real(musica_dk), intent(out) ::  force(length,number_of_species)    ! rate of change of each molecule

    ! Local variables
    real(musica_dk)              ::  number_density_air(length)
    integer                      ::  i

    !$acc enter data create(number_density_air) async(STREAM0)

    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector
    do i = 1, length 
       number_density_air(i) = this%environment(i)%number_density_air
    end do
    !$acc end parallel

    !force = p_force( vmr, this%rates, this%number_density, this%rateConst )
    call p_force( this%rateConst, vmr, number_density_air, force)

    !$acc exit data delete(number_density_air) async(STREAM0)

  end subroutine calc_force

  !---------------------------
  ! Calculate the rates for each chemical reaction
  !---------------------------
  subroutine calc_reaction_rates( this, number_density, reaction_rates )

     use kinetics_utilities, only : rxn_rates => calc_reaction_rates, &
                                    nRxn => number_of_reactions
     use factor_solve_utilities, only : number_of_species

     class(kinetics_t), intent(in) :: this
     real(musica_dk), intent(in)   :: number_density(length,number_of_species)  ! number densities of each component (#/cm^3)
     real(musica_dk), intent(out)  :: reaction_rates(length,nRxn)               ! reaction rates

     ! Local variables
     real(musica_dk)               :: number_density_air(length)
     integer                       :: i

     !$acc data copyin  (this,this%rateConst,this%environment, &
     !$acc               number_density) &
     !$acc      copyout (reaction_rates) &
     !$acc      create  (number_density_air)

     !$acc parallel default(present) vector_length(VLEN)
     !$acc loop gang vector
     do i = 1, length 
        number_density_air(i) = this%environment(i)%number_density_air
     end do
     !$acc end parallel

     call rxn_rates( this%rateConst, number_density, number_density_air, reaction_rates )

     !$acc end data 

  end subroutine calc_reaction_rates

  !---------------------------
  ! Get the rate constants for each chemical reaction
  !---------------------------
  subroutine calc_reaction_rate_constants( this, reaction_rate_constants )

    use kinetics_utilities, only : nRxn => number_of_reactions

    class(kinetics_t), intent(in) :: this
    real(musica_dk),  intent(out) :: reaction_rate_constants(length,nRxn) ! reaction rate constants

    integer                       :: i, j

    do j = 1, nRxn
       do i = 1, length
          reaction_rate_constants(i,j) = this%rateConst(i,j)
       end do
    end do

  end subroutine calc_reaction_rate_constants

  subroutine calc_dforce_dy ( this, vmr, dforce_dy )

    use kinetics_utilities, only : p_dforce_dy=>dforce_dy
    use factor_solve_utilities, only: number_sparse_factor_elements, &
                                      number_of_species

    class(kinetics_t) :: this
    real(musica_dk), intent(in)  ::  vmr(length,number_of_species)                    ! volume mixing ratios of each component in order
    real(musica_dk), intent(out) ::  dforce_dy(length,number_sparse_factor_elements)  ! sensitivity of forcing to changes in each vmr

    ! Local variables
    real(musica_dk) :: number_density_air(length)
    integer         :: i

    !$acc enter data create(number_density_air) async(STREAM0)

    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector
    do i = 1, length 
       number_density_air(i) = this%environment(i)%number_density_air
    end do
    !$acc end parallel

    call p_dforce_dy(dforce_dy, this%rateConst, vmr, number_density_air)

    !$acc exit data delete(number_density_air) async(STREAM0)

  end subroutine calc_dforce_dy



  !------------------------------------------------------
  ! allocate rateConst member array
  !------------------------------------------------------
  function constructor( ) result( this )
    use factor_solve_utilities, only: number_sparse_factor_elements, &
                                      number_of_species
    use kinetics_utilities,     only: number_of_reactions
    class(kinetics_t), pointer :: this

    integer :: nRxt    ! total number of reactions
    integer :: nSpecies! total number of reactions

    allocate( this )

    nRxt     = number_of_reactions
    nSpecies = number_of_species

    this%nReact = nRxt
    this%nSpecies = nSpecies

    if( .not. allocated(this%rates)) then
      allocate( this%rates(length,nRxt) )
    else
      write(*,*) 'rates_init: rateConst already allocated'
    endif

    if( .not. allocated(this%rateConst)) then
      allocate( this%rateConst(length,nRxt) )
    else
      write(*,*) 'rateConst_init: rateConst already allocated'
    endif

    if( .not. allocated(this%MBOdeJac)) then
      allocate( this%MBOdeJac(length,number_sparse_factor_elements) )
    else
      write(*,*) 'jacobian_init: MBOdeJac already allocated'
    endif

    if( .not. allocated(this%chemJac)) then
      allocate( this%chemJac(length,number_sparse_factor_elements) )
    else
      write(*,*) 'jacobian_init: chemJac already allocated'
    endif

    if( .not. allocated(this%Pivot)) then
      allocate( this%Pivot(length,nSpecies) )
    else
      write(*,*) 'jacobian_init: Pivot already allocated'
    endif

    !$acc enter data copyin(this) &
    !$acc            create(this%chemJac,this%MBOdeJac,this%rateConst) &
    !$acc            async(STREAM0)

  end function constructor

  !------------------------------------------------------
  ! prepare the rosenbrock ode solver matrix
  !------------------------------------------------------
  subroutine LinFactor( this, H, gam, Y, Singular, istatus )

   use kinetics_utilities, only : factored_alpha_minus_jac
   use factor_solve_utilities, only: number_sparse_factor_elements, &
                                     number_of_species

   class(kinetics_t) :: this
   real(musica_dk), intent(inout) :: H          ! time step (seconds)
   real(musica_dk), intent(in)    :: gam        ! time step factor for specific rosenbrock method
   real(musica_dk), intent(in)    :: Y(length,number_of_species)     ! constituent concentration (molec/cm^3)
   logical, intent(inout)         :: Singular   ! singularity flag (T or F)
   integer, intent(inout)         :: istatus(:) ! rosenbrock status vector

   ! Local variables

   integer, parameter :: Ndec = 6, Nsng = 8
   real(musica_dk), parameter :: ONE  = 1._musica_dk
   real(musica_dk), parameter :: HALF = .5_musica_dk

   INTEGER  :: i, j, k, ising, Nconsecutive
   REAL(musica_dk) :: ghinv
   REAL(musica_dk) :: LU_factored(length,number_sparse_factor_elements)

   !$acc enter data create(LU_factored) async(STREAM0)

! Set the chemical entries for the Ode Jacobian
   call this%calc_dforce_dy( Y, this%MBOdeJac )

   !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
   !$acc loop gang vector collapse(2)
   do k = 1, number_sparse_factor_elements
      do j = 1, length
         this%chemJac(j,k) = this%MBOdeJac(j,k)
      end do
   end do
   !$acc end parallel

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)
     ghinv = ONE/(H*gam)
!    Compute LU decomposition of [ghinv*I - Ghimj]
     call factored_alpha_minus_jac( LU_factored, ghinv, this%MBOdeJac )
     ising = 0;
     istatus(Ndec) = istatus(Ndec) + 1
     IF (ising == 0) THEN
!~~~>    If successful done
       !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
       !$acc loop gang vector collapse(2)
       do k = 1, number_sparse_factor_elements
          do j = 1, length
             this%MBOdeJac(j,k) = LU_factored(j,k)
          end do
       end do
       !$acc end parallel
       Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
       istatus(Nsng) = istatus(Nsng) + 1
       Nconsecutive = Nconsecutive + 1
       Singular = .TRUE.
       write(*,*) 'Warning: LU Decomposition returning ISING = ',ising
       IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
         H = H*HALF
       ELSE  ! More than 5 consecutive failed decompositions
         RETURN
       END IF  ! Nconsecutive
     END IF
   END DO

   !$acc exit data delete(LU_factored) async(STREAM0)

  end subroutine LinFactor

  !------------------------------------------------------
  ! update rateConst
  ! Execute once for the chemistry-time-step advance
  ! Not called from the solver
  !------------------------------------------------------
  subroutine update( this, environment )

    use rate_constants_utility,        only : calculate_rate_constants
    use kinetics_utilities,            only : number_of_reactions

    !> Kinetics calculator
    class(kinetics_t), intent(inout) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment(length)

    ! Local variables
    integer :: i, k

    ! save the environmental conditions
    if( .not. allocated( this%environment ) ) then
      allocate( this%environment(length), source = environment )
      !$acc enter data copyin(this%environment) &
      !$acc            async(STREAM0)
    else
      !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
      !$acc loop gang vector
      do i = 1, length
         this%environment(i) = environment(i)
      end do
      !$acc end parallel
    end if

    ! update the reaction rate constants
    call calculate_rate_constants( this%rateConst, environment )

  end subroutine update

  subroutine rateConst_print( this )

    class(kinetics_t) :: this

    write(*,*) 'rate constants:'
    write(*,'(1p,5(1x,g0))') this%rateConst(1,:)

  end subroutine rateConst_print

  !---------------------------
  !  cleanup when k_rateConst type is removed
  !---------------------------
  subroutine DasEnder( this )

    type(kinetics_t) :: this

    if( allocated( this%rates ) ) then
       deallocate( this%rates )
    endif
    if( allocated( this%rateConst ) ) then
       deallocate( this%rateConst )
    endif
    if( allocated( this%MBOdeJac ) ) then
       deallocate( this%MBOdeJac )
    endif
    if( allocated( this%chemJac ) ) then
       deallocate( this%chemJac )
    endif
    if( allocated( this%Pivot ) ) then
       deallocate( this%Pivot )
    endif

  end subroutine DasEnder


  !---------------------------
  ! Compute dforce/dy = y''
  !---------------------------
  function dForcedyxForce( this, force ) result(d2Fdy2)

   use kinetics_utilities, only: dforce_dy_times_vector
   use factor_solve_utilities, only: number_of_species

   class(kinetics_t)           :: this
   real(musica_dk), intent(in) :: force(length,number_of_species)         ! chem forcing; dy/dt

   ! Local variables
   real(musica_dk) ::  d2Fdy2(length,number_of_species)

   !$acc data copyout (d2Fdy2) &
   !$acc      copyin  (this,this%chemJac,force)

    call dforce_dy_times_vector( this%chemJac, force, d2Fdy2 )

   !$acc end data

  end function dForcedyxForce

  SUBROUTINE LinSolve( this, B)

    use factor_solve_utilities, only : solve, number_of_species

    class(kinetics_t)              :: this
    REAL(musica_dk), INTENT(INOUT) :: B(length,number_of_species)

    ! Local variables
    REAL(musica_dk)                :: x(length,number_of_species)
    integer                        :: i, j

    !$acc enter data create(x) async(STREAM0)

    call solve ( this%MBOdeJac, x, B )

    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector collapse(2)
    do j = 1, number_of_species
       do i = 1, length
          B(i,j) = x(i,j)
       end do
    end do
    !$acc end parallel

    !$acc exit data delete(x) async(STREAM0)

  END SUBROUTINE LinSolve


end module micm_kinetics

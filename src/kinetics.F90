!> \file
!> The micm_kinetics module

!> The kinetics_t type and related functions
!!
!! \todo Consider rewriting micm_kinetics to follow musica style conventions
module micm_kinetics

  use micm_environment,                only : environment_t
  use musica_constants,                only : musica_dk, musica_ik

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
  integer, allocatable         :: Pivot(:)
  real(musica_dk), allocatable :: rateConst(:)
  real(musica_dk), allocatable :: MBOdeJac(:)      ! ODE solver jacobian
  real(musica_dk), allocatable :: chemJac(:)       ! chemistry forcing jacobian
  real(musica_dk), allocatable :: rates(:)         ! rates of reactions
  type(environment_t) :: environment
contains
  procedure, public :: species_names
  procedure, public :: reaction_names
  procedure, public :: update
  procedure, public :: rateConst_print
  procedure, public :: LinFactor
  procedure, public :: LinSolve
  procedure, public :: force
  procedure, public :: reaction_rates
  procedure, public :: reaction_rate_constants
  procedure, public :: dforce_dy
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

    use kinetics_utilities, only : number_of_species, &
                                   get_names => species_names
    use musica_string,      only : string_t

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
  ! Compute time rate of change of each molecule (vmr) given reaction rates
  !---------------------------
  function force( this, vmr )

    use kinetics_utilities,only :  p_force

    class(kinetics_t) :: this
    real(musica_dk), intent(in)  ::  vmr(:)              ! volume mixing ratios of each component in order

    real(musica_dk)              ::  force(size(vmr))    ! rate of change of each molecule

    !force = p_force( vmr, this%rates, this%number_density, this%rateConst )
    call p_force( this%rateConst, vmr, this%environment%number_density_air, force)

  end function force

  !---------------------------
  ! Calculate the rates for each chemical reaction
  !---------------------------
  function reaction_rates( this, number_density )

     use kinetics_utilities, only : rxn_rates => reaction_rates, &
                                    nRxn => number_of_reactions

     class(kinetics_t), intent(in) :: this
     real(musica_dk), intent(in)      ::  number_density(:)    ! number densities of each component (#/cm^3)
     real(musica_dk)                  ::  reaction_rates(nRxn) ! reaction rates

     reaction_rates = rxn_rates( this%rateConst, number_density, this%environment%number_density_air )

  end function reaction_rates

  !---------------------------
  ! Get the rate constants for each chemical reaction
  !---------------------------
  function reaction_rate_constants( this )

    use kinetics_utilities, only : nRxn => number_of_reactions

    class(kinetics_t), intent(in) :: this
    real(musica_dk)                  :: reaction_rate_constants(nRxn) ! reaction rate constants

    reaction_rate_constants(:) = this%rateConst(:)

  end function reaction_rate_constants

  function dforce_dy( this, vmr)

    use kinetics_utilities, only : p_dforce_dy=>dforce_dy
    use factor_solve_utilities, only: number_sparse_factor_elements

    class(kinetics_t) :: this
    real(musica_dk), intent(in)  ::  vmr(:)              ! volume mixing ratios of each component in order

    real(musica_dk) :: dforce_dy(number_sparse_factor_elements)   ! sensitivity of forcing to changes in each vmr

    call p_dforce_dy(dforce_dy, this%rateConst, vmr, this%environment%number_density_air)

  end function dforce_dy



  !------------------------------------------------------
  ! allocate rateConst member array
  !------------------------------------------------------
  function constructor( ) result( this )
    use factor_solve_utilities, only: number_sparse_factor_elements
    use kinetics_utilities,     only: number_of_species, number_of_reactions
    class(kinetics_t), pointer :: this

    integer :: nRxt    ! total number of reactions
    integer :: nSpecies! total number of reactions

    allocate( this )

    nRxt     = number_of_reactions
    nSpecies = number_of_species

    this%nReact = nRxt
    this%nSpecies = nSpecies

    if( .not. allocated(this%rates)) then
      allocate( this%rates(nRxt) )
    else
      write(*,*) 'rates_init: rateConst already allocated'
    endif

    if( .not. allocated(this%rateConst)) then
      allocate( this%rateConst(nRxt) )
    else
      write(*,*) 'rateConst_init: rateConst already allocated'
    endif

    if( .not. allocated(this%MBOdeJac)) then
      allocate( this%MBOdeJac(number_sparse_factor_elements) )
    else
      write(*,*) 'jacobian_init: MBOdeJac already allocated'
    endif

    if( .not. allocated(this%chemJac)) then
      allocate( this%chemJac(number_sparse_factor_elements) )
    else
      write(*,*) 'jacobian_init: chemJac already allocated'
    endif

    if( .not. allocated(this%Pivot)) then
      allocate( this%Pivot(nSpecies) )
    else
      write(*,*) 'jacobian_init: Pivot already allocated'
    endif

  end function constructor

  !------------------------------------------------------
  ! prepare the rosenbrock ode solver matrix
  !------------------------------------------------------
  subroutine LinFactor( this, H, gam, Y, Singular, istatus )

   use kinetics_utilities, only : factored_alpha_minus_jac
   use factor_solve_utilities, only: number_sparse_factor_elements

    class(kinetics_t) :: this
    real(musica_dk), intent(inout) :: H          ! time step (seconds)
    real(musica_dk), intent(in)    :: gam        ! time step factor for specific rosenbrock method
    real(musica_dk), intent(in)    :: Y(:)       ! constituent concentration (molec/cm^3)
    logical, intent(inout)         :: Singular   ! singularity flag (T or F)
    integer, intent(inout)         :: istatus(:) ! rosenbrock status vector

    integer, parameter :: Ndec = 6, Nsng = 8
    real(musica_dk), parameter :: ONE  = 1._musica_dk
    real(musica_dk), parameter :: HALF = .5_musica_dk

    INTEGER  :: i, ising, Nconsecutive
    REAL(musica_dk) :: ghinv
    REAL(musica_dk) :: LU_factored(number_sparse_factor_elements)

   associate( Ghimj => this%MBOdeJac )

! Set the chemical entries for the Ode Jacobian
    Ghimj(:) = this%dforce_dy( Y )
    this%chemJac(:) = Ghimj(:)

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)
     ghinv = ONE/(H*gam)
!    Compute LU decomposition of [ghinv*I - Ghimj]
     call factored_alpha_minus_jac( LU_factored, ghinv, Ghimj )
     ising = 0;
     istatus(Ndec) = istatus(Ndec) + 1
     IF (ising == 0) THEN
!~~~>    If successful done
       Ghimj(:) = LU_factored(:)
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

   end associate

  end subroutine LinFactor

  !------------------------------------------------------
  ! update rateConst
  ! Execute once for the chemistry-time-step advance
  ! Not called from the solver
  !------------------------------------------------------
  subroutine update( this, environment )

    use rate_constants_utility,        only : k_rate_constant

    !> Kinetics calculator
    class(kinetics_t), intent(inout) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment

    integer :: i, size_krateConst, size_jrateConst

    ! save the environmental conditions
    this%environment = environment

    !> \todo accept external photolysis rates into kinetics_t

    ! recalculate the reaction rates
    this%rateConst(:) = 0.0
    call k_rate_constant( this%rateConst, environment )

  end subroutine update

  subroutine rateConst_print( this )

    class(kinetics_t) :: this

    write(*,*) 'rate constants:'
    write(*,'(1p,5(1x,g0))') this%rateConst(:)

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
    class(kinetics_t) :: this
    real(musica_dk), intent(in)::  force(:)         ! chem forcing; dy/dt

    real(musica_dk) ::  d2Fdy2(size(force))

    call dforce_dy_times_vector( this%chemJac, force, d2Fdy2 )

  end function dForcedyxForce

    SUBROUTINE LinSolve( this, B)

      use factor_solve_utilities, only : solve
      class(kinetics_t) :: this
      REAL(musica_dk), INTENT(INOUT) :: B(:)
      REAL(musica_dk)                :: x(size(B))

      call solve ( this%MBOdeJac, x, B )

      B(:) = x(:)


    END SUBROUTINE LinSolve


end module micm_kinetics

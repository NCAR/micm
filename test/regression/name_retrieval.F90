module reaction_photolysis_species_names
  use iso_c_binding

  implicit none

  character(kind=c_char,len=:), allocatable, target, save :: rnames(:)
  character(kind=c_char,len=:), allocatable, target, save :: pnames(:)
  character(kind=c_char,len=:), allocatable, target, save :: snames(:)

contains

  subroutine Finit() bind(C, name="Finit")
    use kinetics_utilities, only : reaction_names, photolysis_names, species_names
    rnames = reaction_names()
    pnames = photolysis_names()
    snames = species_names()
  end subroutine

  subroutine get_reaction_names( r_names ) bind(C, name="get_reaction_names")
    ! Argument list
    character(kind=c_char,len=:), pointer, intent(inout) :: r_names(:)
    r_names => rnames
  end subroutine

  subroutine get_photolysis_names( p_names ) bind(C, name="get_photolysis_names")
    ! Argument list
    character(kind=c_char,len=:), pointer, intent(inout) :: p_names(:)
    p_names => pnames
  end subroutine

  subroutine get_species_names( s_names ) bind(C, name="get_species_names")
    ! Argument list
    character(kind=c_char,len=:), pointer, intent(inout) :: s_names(:)
    s_names => snames
  end subroutine

end module reaction_photolysis_species_names
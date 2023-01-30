module reaction_photolysis_species_names
  use iso_c_binding

  implicit none

  character(kind=c_char,len=:), allocatable, target, save :: names(:)

contains

  subroutine Finit() bind(C, name="Finit")
    use kinetics_utilities, only : reaction_names
    names = reaction_names()
  end subroutine

  subroutine get_names( pnames ) bind(C, name="get_names")
    ! Argument list
    character(kind=c_char,len=:), pointer, intent(inout) :: pnames(:)
    pnames => names
  end subroutine

end module reaction_photolysis_species_names
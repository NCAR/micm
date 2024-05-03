
      module mo_make_grp_vmr

      private
      public :: mak_grp_vmr

      contains

      subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs, chnkpnts )

      use chem_mods, only : grpcnt, gas_pcnst
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------------------------------------
!        ... dummy arguments
!----------------------------------------------------------------------------
      integer, intent(in) :: chnkpnts
      real(r8), intent(in)    :: vmr(chnkpnts,max(1,gas_pcnst))
      real(r8), intent(in)    :: group_ratios(chnkpnts,max(1,grpcnt))
      real(r8), intent(out)   :: group_vmrs(chnkpnts,max(1,grpcnt))

!----------------------------------------------------------------------------
!        ... local variables
!----------------------------------------------------------------------------
      integer ::  k

      end subroutine mak_grp_vmr

      end module mo_make_grp_vmr

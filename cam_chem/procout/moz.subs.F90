
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )
 
      use ppgrid, only : pcols, pver


      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol*pver)
      real(r8), intent(inout) :: rate(ncol*pver,max(1,rxntot))

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      integer   ::  offset
      real(r8)  :: itemp(ncol*pver)
      real(r8)  :: exp_fac(ncol*pver)
      real(r8)  :: ko(ncol*pver)
      real(r8)  :: kinf(ncol*pver)

      rate(:,124) = 1.29e-07_r8
      rate(:,125) = 1.2e-10_r8
      rate(:,129) = 1.2e-10_r8
      rate(:,130) = 1.2e-10_r8
      rate(:,136) = 6.9e-12_r8
      rate(:,137) = 7.2e-11_r8
      rate(:,138) = 1.6e-12_r8
      rate(:,144) = 1.8e-12_r8
      rate(:,148) = 1.8e-12_r8
      rate(:,160) = 3.5e-12_r8
      rate(:,162) = 1.3e-11_r8
      rate(:,163) = 2.2e-11_r8
      rate(:,164) = 5e-11_r8
      rate(:,199) = 1.7e-13_r8
      rate(:,201) = 2.607e-10_r8
      rate(:,202) = 9.75e-11_r8
      rate(:,203) = 2.07e-10_r8
      rate(:,204) = 2.088e-10_r8
      rate(:,205) = 1.17e-10_r8
      rate(:,206) = 4.644e-11_r8
      rate(:,207) = 1.204e-10_r8
      rate(:,208) = 9.9e-11_r8
      rate(:,209) = 3.3e-12_r8
      rate(:,228) = 4.5e-11_r8
      rate(:,229) = 4.62e-10_r8
      rate(:,230) = 1.2e-10_r8
      rate(:,231) = 9e-11_r8
      rate(:,232) = 3e-11_r8
      rate(:,237) = 2.14e-11_r8
      rate(:,238) = 1.9e-10_r8
      rate(:,251) = 2.57e-10_r8
      rate(:,252) = 1.8e-10_r8
      rate(:,253) = 1.794e-10_r8
      rate(:,254) = 1.3e-10_r8
      rate(:,255) = 7.65e-11_r8
      rate(:,268) = 4e-13_r8
      rate(:,272) = 1.31e-10_r8
      rate(:,273) = 3.5e-11_r8
      rate(:,274) = 9e-12_r8
      rate(:,281) = 6.8e-14_r8
      rate(:,282) = 2e-13_r8
      rate(:,297) = 1e-12_r8
      rate(:,301) = 1e-14_r8
      rate(:,302) = 1e-11_r8
      rate(:,303) = 1.15e-11_r8
      rate(:,304) = 4e-14_r8
      rate(:,317) = 3e-12_r8
      rate(:,318) = 6.7e-13_r8
      rate(:,328) = 3.5e-13_r8
      rate(:,329) = 5.4e-11_r8
      rate(:,332) = 2e-12_r8
      rate(:,333) = 1.4e-11_r8
      rate(:,336) = 2.4e-12_r8
      rate(:,347) = 5e-12_r8
      rate(:,357) = 1.6e-12_r8
      rate(:,359) = 6.7e-12_r8
      rate(:,362) = 3.5e-12_r8
      rate(:,365) = 1.3e-11_r8
      rate(:,366) = 1.4e-11_r8
      rate(:,370) = 2.4e-12_r8
      rate(:,371) = 1.4e-11_r8
      rate(:,376) = 2.4e-12_r8
      rate(:,377) = 4e-11_r8
      rate(:,378) = 4e-11_r8
      rate(:,380) = 1.4e-11_r8
      rate(:,384) = 2.4e-12_r8
      rate(:,385) = 4e-11_r8
      rate(:,389) = 7e-11_r8
      rate(:,390) = 1e-10_r8
      rate(:,395) = 2.4e-12_r8
      rate(:,410) = 4.7e-11_r8
      rate(:,423) = 2.1e-12_r8
      rate(:,424) = 2.8e-13_r8
      rate(:,432) = 1.7e-11_r8
      rate(:,438) = 8.4e-11_r8
      rate(:,440) = 1.9e-11_r8
      rate(:,441) = 1.2e-14_r8
      rate(:,442) = 2e-10_r8
      rate(:,449) = 2.4e-12_r8
      rate(:,450) = 2e-11_r8
      rate(:,454) = 2.3e-11_r8
      rate(:,455) = 2e-11_r8
      rate(:,459) = 3.3e-11_r8
      rate(:,460) = 1e-12_r8
      rate(:,461) = 5.7e-11_r8
      rate(:,462) = 3.4e-11_r8
      rate(:,467) = 2.3e-12_r8
      rate(:,469) = 1.2e-11_r8
      rate(:,470) = 5.7e-11_r8
      rate(:,471) = 2.8e-11_r8
      rate(:,472) = 6.6e-11_r8
      rate(:,473) = 1.4e-11_r8
      rate(:,476) = 1.9e-12_r8
      rate(:,488) = 6.34e-08_r8
      rate(:,494) = 1.9e-11_r8
      rate(:,497) = 1.2e-14_r8
      rate(:,498) = 2e-10_r8
      rate(:,509) = 1.34e-11_r8
      rate(:,515) = 1.34e-11_r8
      rate(:,520) = 1.7e-11_r8
      rate(:,540) = 2.31e-07_r8
      rate(:,541) = 2.31e-06_r8
      rate(:,542) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,126) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,127) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:,128) = 3.3e-11_r8 * exp( 55._r8 * itemp(:) )
      rate(:,131) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,134) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,135) = 1.4e-12_r8 * exp_fac(:)
      rate(:,386) = 1.05e-14_r8 * exp_fac(:)
      rate(:,505) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,140) = 3e-11_r8 * exp_fac(:)
      rate(:,226) = 5.5e-12_r8 * exp_fac(:)
      rate(:,265) = 3.8e-12_r8 * exp_fac(:)
      rate(:,286) = 3.8e-12_r8 * exp_fac(:)
      rate(:,313) = 3.8e-12_r8 * exp_fac(:)
      rate(:,321) = 3.8e-12_r8 * exp_fac(:)
      rate(:,325) = 3.8e-12_r8 * exp_fac(:)
      rate(:,341) = 2.3e-11_r8 * exp_fac(:)
      rate(:,351) = 3.8e-12_r8 * exp_fac(:)
      rate(:,361) = 3.8e-12_r8 * exp_fac(:)
      rate(:,388) = 1.52e-11_r8 * exp_fac(:)
      rate(:,396) = 1.52e-12_r8 * exp_fac(:)
      rate(:,402) = 3.8e-12_r8 * exp_fac(:)
      rate(:,405) = 3.8e-12_r8 * exp_fac(:)
      rate(:,409) = 3.8e-12_r8 * exp_fac(:)
      rate(:,425) = 3.8e-12_r8 * exp_fac(:)
      rate(:,429) = 3.8e-12_r8 * exp_fac(:)
      rate(:,435) = 3.8e-12_r8 * exp_fac(:)
      rate(:,439) = 3.8e-12_r8 * exp_fac(:)
      rate(:,141) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,142) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,143) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,145) = 4.8e-11_r8 * exp_fac(:)
      rate(:,224) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,146) = 1.8e-11_r8 * exp_fac(:)
      rate(:,299) = 4.2e-12_r8 * exp_fac(:)
      rate(:,312) = 4.2e-12_r8 * exp_fac(:)
      rate(:,320) = 4.2e-12_r8 * exp_fac(:)
      rate(:,349) = 4.2e-12_r8 * exp_fac(:)
      rate(:,369) = 4.4e-12_r8 * exp_fac(:)
      rate(:,375) = 4.4e-12_r8 * exp_fac(:)
      rate(:,448) = 4.2e-12_r8 * exp_fac(:)
      rate(:,453) = 4.2e-12_r8 * exp_fac(:)
      rate(:,458) = 4.2e-12_r8 * exp_fac(:)
      rate(:,147) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,151) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,152) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,153) = 2.9e-12_r8 * exp_fac(:)
      rate(:,154) = 1.45e-12_r8 * exp_fac(:)
      rate(:,155) = 1.45e-12_r8 * exp_fac(:)
      rate(:,156) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,157) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,158) = 1.2e-13_r8 * exp_fac(:)
      rate(:,184) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,161) = 1.7e-11_r8 * exp_fac(:)
      rate(:,259) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,165) = 3.44e-12_r8 * exp_fac(:)
      rate(:,217) = 2.3e-12_r8 * exp_fac(:)
      rate(:,220) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,166) = 3e-12_r8 * exp_fac(:)
      rate(:,225) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,168) = 7.26e-11_r8 * exp_fac(:)
      rate(:,169) = 4.64e-11_r8 * exp_fac(:)
      rate(:,176) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,177) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,178) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,179) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,180) = 1.4e-11_r8 * exp_fac(:)
      rate(:,194) = 7.4e-12_r8 * exp_fac(:)
      rate(:,295) = 8.1e-12_r8 * exp_fac(:)
      rate(:,181) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,182) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,183) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,185) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,186) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,187) = 2.6e-12_r8 * exp_fac(:)
      rate(:,188) = 6.4e-12_r8 * exp_fac(:)
      rate(:,218) = 4.1e-13_r8 * exp_fac(:)
      rate(:,398) = 7.5e-12_r8 * exp_fac(:)
      rate(:,412) = 7.5e-12_r8 * exp_fac(:)
      rate(:,415) = 7.5e-12_r8 * exp_fac(:)
      rate(:,418) = 7.5e-12_r8 * exp_fac(:)
      rate(:,189) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,191) = 3.6e-12_r8 * exp_fac(:)
      rate(:,240) = 2e-12_r8 * exp_fac(:)
      rate(:,192) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,193) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,195) = 6e-13_r8 * exp_fac(:)
      rate(:,215) = 1.5e-12_r8 * exp_fac(:)
      rate(:,223) = 1.9e-11_r8 * exp_fac(:)
      rate(:,196) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,197) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,198) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,200) = 3e-12_r8 * exp_fac(:)
      rate(:,234) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,212) = 1.7e-11_r8 * exp_fac(:)
      rate(:,239) = 6.3e-12_r8 * exp_fac(:)
      rate(:,213) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,214) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,216) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,219) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,222) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,227) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,233) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,235) = 1.4e-11_r8 * exp_fac(:)
      rate(:,237) = 2.14e-11_r8 * exp_fac(:)
      rate(:,238) = 1.9e-10_r8 * exp_fac(:)
      rate(:,251) = 2.57e-10_r8 * exp_fac(:)
      rate(:,252) = 1.8e-10_r8 * exp_fac(:)
      rate(:,253) = 1.794e-10_r8 * exp_fac(:)
      rate(:,254) = 1.3e-10_r8 * exp_fac(:)
      rate(:,255) = 7.65e-11_r8 * exp_fac(:)
      rate(:,268) = 4e-13_r8 * exp_fac(:)
      rate(:,272) = 1.31e-10_r8 * exp_fac(:)
      rate(:,273) = 3.5e-11_r8 * exp_fac(:)
      rate(:,274) = 9e-12_r8 * exp_fac(:)
      rate(:,281) = 6.8e-14_r8 * exp_fac(:)
      rate(:,282) = 2e-13_r8 * exp_fac(:)
      rate(:,297) = 1e-12_r8 * exp_fac(:)
      rate(:,301) = 1e-14_r8 * exp_fac(:)
      rate(:,302) = 1e-11_r8 * exp_fac(:)
      rate(:,303) = 1.15e-11_r8 * exp_fac(:)
      rate(:,304) = 4e-14_r8 * exp_fac(:)
      rate(:,317) = 3e-12_r8 * exp_fac(:)
      rate(:,318) = 6.7e-13_r8 * exp_fac(:)
      rate(:,328) = 3.5e-13_r8 * exp_fac(:)
      rate(:,329) = 5.4e-11_r8 * exp_fac(:)
      rate(:,332) = 2e-12_r8 * exp_fac(:)
      rate(:,333) = 1.4e-11_r8 * exp_fac(:)
      rate(:,336) = 2.4e-12_r8 * exp_fac(:)
      rate(:,347) = 5e-12_r8 * exp_fac(:)
      rate(:,357) = 1.6e-12_r8 * exp_fac(:)
      rate(:,359) = 6.7e-12_r8 * exp_fac(:)
      rate(:,362) = 3.5e-12_r8 * exp_fac(:)
      rate(:,365) = 1.3e-11_r8 * exp_fac(:)
      rate(:,366) = 1.4e-11_r8 * exp_fac(:)
      rate(:,370) = 2.4e-12_r8 * exp_fac(:)
      rate(:,371) = 1.4e-11_r8 * exp_fac(:)
      rate(:,376) = 2.4e-12_r8 * exp_fac(:)
      rate(:,377) = 4e-11_r8 * exp_fac(:)
      rate(:,378) = 4e-11_r8 * exp_fac(:)
      rate(:,380) = 1.4e-11_r8 * exp_fac(:)
      rate(:,384) = 2.4e-12_r8 * exp_fac(:)
      rate(:,385) = 4e-11_r8 * exp_fac(:)
      rate(:,389) = 7e-11_r8 * exp_fac(:)
      rate(:,390) = 1e-10_r8 * exp_fac(:)
      rate(:,395) = 2.4e-12_r8 * exp_fac(:)
      rate(:,410) = 4.7e-11_r8 * exp_fac(:)
      rate(:,423) = 2.1e-12_r8 * exp_fac(:)
      rate(:,424) = 2.8e-13_r8 * exp_fac(:)
      rate(:,432) = 1.7e-11_r8 * exp_fac(:)
      rate(:,438) = 8.4e-11_r8 * exp_fac(:)
      rate(:,440) = 1.9e-11_r8 * exp_fac(:)
      rate(:,441) = 1.2e-14_r8 * exp_fac(:)
      rate(:,442) = 2e-10_r8 * exp_fac(:)
      rate(:,449) = 2.4e-12_r8 * exp_fac(:)
      rate(:,450) = 2e-11_r8 * exp_fac(:)
      rate(:,454) = 2.3e-11_r8 * exp_fac(:)
      rate(:,455) = 2e-11_r8 * exp_fac(:)
      rate(:,459) = 3.3e-11_r8 * exp_fac(:)
      rate(:,460) = 1e-12_r8 * exp_fac(:)
      rate(:,461) = 5.7e-11_r8 * exp_fac(:)
      rate(:,462) = 3.4e-11_r8 * exp_fac(:)
      rate(:,467) = 2.3e-12_r8 * exp_fac(:)
      rate(:,469) = 1.2e-11_r8 * exp_fac(:)
      rate(:,470) = 5.7e-11_r8 * exp_fac(:)
      rate(:,471) = 2.8e-11_r8 * exp_fac(:)
      rate(:,472) = 6.6e-11_r8 * exp_fac(:)
      rate(:,473) = 1.4e-11_r8 * exp_fac(:)
      rate(:,476) = 1.9e-12_r8 * exp_fac(:)
      rate(:,488) = 6.34e-08_r8 * exp_fac(:)
      rate(:,494) = 1.9e-11_r8 * exp_fac(:)
      rate(:,497) = 1.2e-14_r8 * exp_fac(:)
      rate(:,498) = 2e-10_r8 * exp_fac(:)
      rate(:,509) = 1.34e-11_r8 * exp_fac(:)
      rate(:,515) = 1.34e-11_r8 * exp_fac(:)
      rate(:,520) = 1.7e-11_r8 * exp_fac(:)
      rate(:,540) = 2.31e-07_r8 * exp_fac(:)
      rate(:,541) = 2.31e-06_r8 * exp_fac(:)
      rate(:,542) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,236) = 6e-12_r8 * exp_fac(:)
      rate(:,334) = 5e-13_r8 * exp_fac(:)
      rate(:,367) = 5e-13_r8 * exp_fac(:)
      rate(:,372) = 5e-13_r8 * exp_fac(:)
      rate(:,381) = 5e-13_r8 * exp_fac(:)
      rate(:,392) = 5e-13_r8 * exp_fac(:)
      rate(:,241) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,242) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,243) = 1.64e-12_r8 * exp_fac(:)
      rate(:,353) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,244) = 2.03e-11_r8 * exp_fac(:)
      rate(:,475) = 3.4e-12_r8 * exp_fac(:)
      rate(:,245) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,246) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,247) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,248) = 1.25e-12_r8 * exp_fac(:)
      rate(:,258) = 3.4e-11_r8 * exp_fac(:)
      rate(:,249) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,250) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,256) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,257) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,260) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,261) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,262) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,263) = 2.8e-12_r8 * exp_fac(:)
      rate(:,324) = 2.9e-12_r8 * exp_fac(:)
      rate(:,264) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,266) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,269) = 7.5e-13_r8 * exp_fac(:)
      rate(:,283) = 7.5e-13_r8 * exp_fac(:)
      rate(:,298) = 7.5e-13_r8 * exp_fac(:)
      rate(:,311) = 7.5e-13_r8 * exp_fac(:)
      rate(:,319) = 7.5e-13_r8 * exp_fac(:)
      rate(:,323) = 8.6e-13_r8 * exp_fac(:)
      rate(:,335) = 8e-13_r8 * exp_fac(:)
      rate(:,348) = 7.5e-13_r8 * exp_fac(:)
      rate(:,358) = 7.5e-13_r8 * exp_fac(:)
      rate(:,368) = 8e-13_r8 * exp_fac(:)
      rate(:,373) = 8e-13_r8 * exp_fac(:)
      rate(:,382) = 8e-13_r8 * exp_fac(:)
      rate(:,393) = 8e-13_r8 * exp_fac(:)
      rate(:,400) = 7.5e-13_r8 * exp_fac(:)
      rate(:,404) = 7.5e-13_r8 * exp_fac(:)
      rate(:,407) = 7.5e-13_r8 * exp_fac(:)
      rate(:,420) = 7.5e-13_r8 * exp_fac(:)
      rate(:,427) = 7.5e-13_r8 * exp_fac(:)
      rate(:,433) = 7.5e-13_r8 * exp_fac(:)
      rate(:,436) = 7.5e-13_r8 * exp_fac(:)
      rate(:,447) = 7.5e-13_r8 * exp_fac(:)
      rate(:,452) = 7.5e-13_r8 * exp_fac(:)
      rate(:,457) = 7.5e-13_r8 * exp_fac(:)
      rate(:,500) = 7.5e-13_r8 * exp_fac(:)
      rate(:,507) = 7.5e-13_r8 * exp_fac(:)
      rate(:,517) = 7.5e-13_r8 * exp_fac(:)
      rate(:,521) = 7.5e-13_r8 * exp_fac(:)
      rate(:,270) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,271) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,275) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,280) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,284) = 2.6e-12_r8 * exp_fac(:)
      rate(:,401) = 2.6e-12_r8 * exp_fac(:)
      rate(:,406) = 2.6e-12_r8 * exp_fac(:)
      rate(:,408) = 2.6e-12_r8 * exp_fac(:)
      rate(:,421) = 2.6e-12_r8 * exp_fac(:)
      rate(:,428) = 2.6e-12_r8 * exp_fac(:)
      rate(:,434) = 2.6e-12_r8 * exp_fac(:)
      rate(:,437) = 2.6e-12_r8 * exp_fac(:)
      rate(:,501) = 2.6e-12_r8 * exp_fac(:)
      rate(:,508) = 2.6e-12_r8 * exp_fac(:)
      rate(:,518) = 2.6e-12_r8 * exp_fac(:)
      rate(:,522) = 2.6e-12_r8 * exp_fac(:)
      rate(:,285) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,287) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,288) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,289) = 1.4e-12_r8 * exp_fac(:)
      rate(:,309) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,290) = 4.63e-12_r8 * exp_fac(:)
      rate(:,504) = 2.7e-12_r8 * exp_fac(:)
      rate(:,291) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,292) = 2.9e-12_r8 * exp_fac(:)
      rate(:,293) = 2e-12_r8 * exp_fac(:)
      rate(:,322) = 7.1e-13_r8 * exp_fac(:)
      rate(:,343) = 2e-12_r8 * exp_fac(:)
      rate(:,446) = 2e-12_r8 * exp_fac(:)
      rate(:,451) = 2e-12_r8 * exp_fac(:)
      rate(:,456) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,294) = 4.3e-13_r8 * exp_fac(:)
      rate(:,344) = 4.3e-13_r8 * exp_fac(:)
      rate(:,397) = 4.3e-13_r8 * exp_fac(:)
      rate(:,411) = 4.3e-13_r8 * exp_fac(:)
      rate(:,414) = 4.3e-13_r8 * exp_fac(:)
      rate(:,417) = 4.3e-13_r8 * exp_fac(:)
      rate(:,296) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,300) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,308) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,310) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,314) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,315) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,316) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,330) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,331) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,337) = 2.7e-12_r8 * exp_fac(:)
      rate(:,338) = 1.3e-13_r8 * exp_fac(:)
      rate(:,340) = 9.6e-12_r8 * exp_fac(:)
      rate(:,346) = 5.3e-12_r8 * exp_fac(:)
      rate(:,383) = 2.7e-12_r8 * exp_fac(:)
      rate(:,394) = 2.7e-12_r8 * exp_fac(:)
      rate(:,496) = 2.7e-12_r8 * exp_fac(:)
      rate(:,512) = 2.7e-12_r8 * exp_fac(:)
      rate(:,339) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,342) = 4.6e-12_r8 * exp_fac(:)
      rate(:,345) = 2.3e-12_r8 * exp_fac(:)
      rate(:,350) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,354) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,360) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,363) = 1.86e-11_r8 * exp_fac(:)
      rate(:,364) = 1.86e-11_r8 * exp_fac(:)
      rate(:,374) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,379) = 3.03e-12_r8 * exp_fac(:)
      rate(:,502) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,387) = 2.54e-11_r8 * exp_fac(:)
      rate(:,506) = 2.54e-11_r8 * exp_fac(:)
      rate(:,391) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,399) = 2.3e-12_r8 * exp_fac(:)
      rate(:,499) = 2.3e-12_r8 * exp_fac(:)
      rate(:,403) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,422) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,430) = 1.7e-12_r8 * exp_fac(:)
      rate(:,516) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,443) = 1.2e-12_r8 * exp_fac(:)
      rate(:,510) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,444) = 6.3e-16_r8 * exp_fac(:)
      rate(:,513) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,445) = 1.2e-11_r8 * exp_fac(:)
      rate(:,514) = 1.2e-11_r8 * exp_fac(:)
      rate(:,463) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,464) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,465) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,466) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,474) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,477) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,480) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,495) = 2.75e-13_r8 * exp_fac(:)
      rate(:,503) = 2.12e-13_r8 * exp_fac(:)
      rate(:,511) = 2.6e-13_r8 * exp_fac(:)

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,139), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,149), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,159), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,167), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,170), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,171), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,172), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,190), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,210), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,221), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,267), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,277), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,278), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,279), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,305), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,306), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,326), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,352), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,355), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,413), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,416), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,419), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,426), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,468), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )
 
      use ppgrid, only : pcols, pver


      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol*pver)
      real(r8), intent(inout) :: rate(ncol*pver,max(1,rxntot))

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      integer   ::  offset
      integer   ::  k
      real(r8)  :: itemp(ncol*kbot)
      real(r8)  :: exp_fac(ncol*kbot)
      real(r8)  :: ko(ncol*kbot)
      real(r8)  :: kinf(ncol*kbot)
      real(r8)  :: wrk(ncol*kbot)
 
      n = ncol*kbot

      rate(:n,136) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,127) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:n,131) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,140) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,141) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,142) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,145) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,146) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,147) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,152) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,156) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,157) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,165) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,166) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,139) = wrk(:)
























      end subroutine setrxt_hrates

      end module mo_setrxt

      module mo_adjrxt

      private
      public :: adjrxt

      contains

      subroutine adjrxt( rate, inv, m, ncol, nlev )

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : nfs, rxntot

      implicit none

!--------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol, nlev
      real(r8), intent(in)    :: inv(ncol,nlev,nfs)
      real(r8), intent(in)    :: m(ncol,nlev)
      real(r8), intent(inout) :: rate(ncol,nlev,rxntot)

      rate(:,:,   127) = rate(:,:,   127) * inv(:,:, 3)
      rate(:,:,   128) = rate(:,:,   128) * inv(:,:, 2)
      rate(:,:,   132) = rate(:,:,   132) * inv(:,:, 1)
      rate(:,:,   149) = rate(:,:,   149) * inv(:,:, 1)
      rate(:,:,   156) = rate(:,:,   156) * inv(:,:, 2)
      rate(:,:,   159) = rate(:,:,   159) * inv(:,:, 1)
      rate(:,:,   167) = rate(:,:,   167) * inv(:,:, 1)
      rate(:,:,   170) = rate(:,:,   170) * inv(:,:, 1)
      rate(:,:,   171) = rate(:,:,   171) * inv(:,:, 1)
      rate(:,:,   172) = rate(:,:,   172) * inv(:,:, 1)
      rate(:,:,   174) = rate(:,:,   174) * inv(:,:, 1)
      rate(:,:,   175) = rate(:,:,   175) * inv(:,:, 1)
      rate(:,:,   190) = rate(:,:,   190) * inv(:,:, 1)
      rate(:,:,   210) = rate(:,:,   210) * inv(:,:, 1)
      rate(:,:,   211) = rate(:,:,   211) * inv(:,:, 1)
      rate(:,:,   221) = rate(:,:,   221) * inv(:,:, 1)
      rate(:,:,   267) = rate(:,:,   267) * inv(:,:, 1)
      rate(:,:,   277) = rate(:,:,   277) * inv(:,:, 1)
      rate(:,:,   278) = rate(:,:,   278) * inv(:,:, 1)
      rate(:,:,   279) = rate(:,:,   279) * inv(:,:, 1)
      rate(:,:,   301) = rate(:,:,   301) * inv(:,:, 2)
      rate(:,:,   305) = rate(:,:,   305) * inv(:,:, 1)
      rate(:,:,   306) = rate(:,:,   306) * inv(:,:, 1)
      rate(:,:,   307) = rate(:,:,   307) * inv(:,:, 1)
      rate(:,:,   326) = rate(:,:,   326) * inv(:,:, 1)
      rate(:,:,   352) = rate(:,:,   352) * inv(:,:, 1)
      rate(:,:,   355) = rate(:,:,   355) * inv(:,:, 1)
      rate(:,:,   356) = rate(:,:,   356) * inv(:,:, 1)
      rate(:,:,   413) = rate(:,:,   413) * inv(:,:, 1)
      rate(:,:,   416) = rate(:,:,   416) * inv(:,:, 1)
      rate(:,:,   419) = rate(:,:,   419) * inv(:,:, 1)
      rate(:,:,   426) = rate(:,:,   426) * inv(:,:, 1)
      rate(:,:,   431) = rate(:,:,   431) * inv(:,:, 1)
      rate(:,:,   467) = rate(:,:,   467) * inv(:,:, 2)
      rate(:,:,   468) = rate(:,:,   468) * inv(:,:, 1)
      rate(:,:,   474) = rate(:,:,   474) * inv(:,:, 2)
      rate(:,:,   133) = rate(:,:,   133) * inv(:,:, 2) * inv(:,:, 1)
      rate(:,:,   139) = rate(:,:,   139) * inv(:,:, 2) * inv(:,:, 1)
      rate(:,:,   125) = rate(:,:,   125) * m(:,:)
      rate(:,:,   126) = rate(:,:,   126) * m(:,:)
      rate(:,:,   129) = rate(:,:,   129) * m(:,:)
      rate(:,:,   130) = rate(:,:,   130) * m(:,:)
      rate(:,:,   131) = rate(:,:,   131) * m(:,:)
      rate(:,:,   132) = rate(:,:,   132) * m(:,:)
      rate(:,:,   134) = rate(:,:,   134) * m(:,:)
      rate(:,:,   135) = rate(:,:,   135) * m(:,:)
      rate(:,:,   136) = rate(:,:,   136) * m(:,:)
      rate(:,:,   137) = rate(:,:,   137) * m(:,:)
      rate(:,:,   138) = rate(:,:,   138) * m(:,:)
      rate(:,:,   140) = rate(:,:,   140) * m(:,:)
      rate(:,:,   141) = rate(:,:,   141) * m(:,:)
      rate(:,:,   142) = rate(:,:,   142) * m(:,:)
      rate(:,:,   143) = rate(:,:,   143) * m(:,:)
      rate(:,:,   144) = rate(:,:,   144) * m(:,:)
      rate(:,:,   145) = rate(:,:,   145) * m(:,:)
      rate(:,:,   146) = rate(:,:,   146) * m(:,:)
      rate(:,:,   147) = rate(:,:,   147) * m(:,:)
      rate(:,:,   148) = rate(:,:,   148) * m(:,:)
      rate(:,:,   149) = rate(:,:,   149) * m(:,:)
      rate(:,:,   150) = rate(:,:,   150) * m(:,:)
      rate(:,:,   151) = rate(:,:,   151) * m(:,:)
      rate(:,:,   152) = rate(:,:,   152) * m(:,:)
      rate(:,:,   153) = rate(:,:,   153) * m(:,:)
      rate(:,:,   154) = rate(:,:,   154) * m(:,:)
      rate(:,:,   155) = rate(:,:,   155) * m(:,:)
      rate(:,:,   157) = rate(:,:,   157) * m(:,:)
      rate(:,:,   158) = rate(:,:,   158) * m(:,:)
      rate(:,:,   159) = rate(:,:,   159) * m(:,:)
      rate(:,:,   160) = rate(:,:,   160) * m(:,:)
      rate(:,:,   161) = rate(:,:,   161) * m(:,:)
      rate(:,:,   162) = rate(:,:,   162) * m(:,:)
      rate(:,:,   163) = rate(:,:,   163) * m(:,:)
      rate(:,:,   164) = rate(:,:,   164) * m(:,:)
      rate(:,:,   165) = rate(:,:,   165) * m(:,:)
      rate(:,:,   166) = rate(:,:,   166) * m(:,:)
      rate(:,:,   167) = rate(:,:,   167) * m(:,:)
      rate(:,:,   168) = rate(:,:,   168) * m(:,:)
      rate(:,:,   169) = rate(:,:,   169) * m(:,:)
      rate(:,:,   170) = rate(:,:,   170) * m(:,:)
      rate(:,:,   171) = rate(:,:,   171) * m(:,:)
      rate(:,:,   172) = rate(:,:,   172) * m(:,:)
      rate(:,:,   173) = rate(:,:,   173) * m(:,:)
      rate(:,:,   176) = rate(:,:,   176) * m(:,:)
      rate(:,:,   177) = rate(:,:,   177) * m(:,:)
      rate(:,:,   178) = rate(:,:,   178) * m(:,:)
      rate(:,:,   179) = rate(:,:,   179) * m(:,:)
      rate(:,:,   180) = rate(:,:,   180) * m(:,:)
      rate(:,:,   181) = rate(:,:,   181) * m(:,:)
      rate(:,:,   182) = rate(:,:,   182) * m(:,:)
      rate(:,:,   183) = rate(:,:,   183) * m(:,:)
      rate(:,:,   184) = rate(:,:,   184) * m(:,:)
      rate(:,:,   185) = rate(:,:,   185) * m(:,:)
      rate(:,:,   186) = rate(:,:,   186) * m(:,:)
      rate(:,:,   187) = rate(:,:,   187) * m(:,:)
      rate(:,:,   188) = rate(:,:,   188) * m(:,:)
      rate(:,:,   189) = rate(:,:,   189) * m(:,:)
      rate(:,:,   190) = rate(:,:,   190) * m(:,:)
      rate(:,:,   191) = rate(:,:,   191) * m(:,:)
      rate(:,:,   192) = rate(:,:,   192) * m(:,:)
      rate(:,:,   193) = rate(:,:,   193) * m(:,:)
      rate(:,:,   194) = rate(:,:,   194) * m(:,:)
      rate(:,:,   195) = rate(:,:,   195) * m(:,:)
      rate(:,:,   196) = rate(:,:,   196) * m(:,:)
      rate(:,:,   197) = rate(:,:,   197) * m(:,:)
      rate(:,:,   198) = rate(:,:,   198) * m(:,:)
      rate(:,:,   199) = rate(:,:,   199) * m(:,:)
      rate(:,:,   200) = rate(:,:,   200) * m(:,:)
      rate(:,:,   201) = rate(:,:,   201) * m(:,:)
      rate(:,:,   202) = rate(:,:,   202) * m(:,:)
      rate(:,:,   203) = rate(:,:,   203) * m(:,:)
      rate(:,:,   204) = rate(:,:,   204) * m(:,:)
      rate(:,:,   205) = rate(:,:,   205) * m(:,:)
      rate(:,:,   206) = rate(:,:,   206) * m(:,:)
      rate(:,:,   207) = rate(:,:,   207) * m(:,:)
      rate(:,:,   208) = rate(:,:,   208) * m(:,:)
      rate(:,:,   209) = rate(:,:,   209) * m(:,:)
      rate(:,:,   210) = rate(:,:,   210) * m(:,:)
      rate(:,:,   212) = rate(:,:,   212) * m(:,:)
      rate(:,:,   213) = rate(:,:,   213) * m(:,:)
      rate(:,:,   214) = rate(:,:,   214) * m(:,:)
      rate(:,:,   215) = rate(:,:,   215) * m(:,:)
      rate(:,:,   216) = rate(:,:,   216) * m(:,:)
      rate(:,:,   217) = rate(:,:,   217) * m(:,:)
      rate(:,:,   218) = rate(:,:,   218) * m(:,:)
      rate(:,:,   219) = rate(:,:,   219) * m(:,:)
      rate(:,:,   220) = rate(:,:,   220) * m(:,:)
      rate(:,:,   221) = rate(:,:,   221) * m(:,:)
      rate(:,:,   222) = rate(:,:,   222) * m(:,:)
      rate(:,:,   223) = rate(:,:,   223) * m(:,:)
      rate(:,:,   224) = rate(:,:,   224) * m(:,:)
      rate(:,:,   225) = rate(:,:,   225) * m(:,:)
      rate(:,:,   226) = rate(:,:,   226) * m(:,:)
      rate(:,:,   227) = rate(:,:,   227) * m(:,:)
      rate(:,:,   228) = rate(:,:,   228) * m(:,:)
      rate(:,:,   229) = rate(:,:,   229) * m(:,:)
      rate(:,:,   230) = rate(:,:,   230) * m(:,:)
      rate(:,:,   231) = rate(:,:,   231) * m(:,:)
      rate(:,:,   232) = rate(:,:,   232) * m(:,:)
      rate(:,:,   233) = rate(:,:,   233) * m(:,:)
      rate(:,:,   234) = rate(:,:,   234) * m(:,:)
      rate(:,:,   235) = rate(:,:,   235) * m(:,:)
      rate(:,:,   236) = rate(:,:,   236) * m(:,:)
      rate(:,:,   237) = rate(:,:,   237) * m(:,:)
      rate(:,:,   238) = rate(:,:,   238) * m(:,:)
      rate(:,:,   239) = rate(:,:,   239) * m(:,:)
      rate(:,:,   240) = rate(:,:,   240) * m(:,:)
      rate(:,:,   241) = rate(:,:,   241) * m(:,:)
      rate(:,:,   242) = rate(:,:,   242) * m(:,:)
      rate(:,:,   243) = rate(:,:,   243) * m(:,:)
      rate(:,:,   244) = rate(:,:,   244) * m(:,:)
      rate(:,:,   245) = rate(:,:,   245) * m(:,:)
      rate(:,:,   246) = rate(:,:,   246) * m(:,:)
      rate(:,:,   247) = rate(:,:,   247) * m(:,:)
      rate(:,:,   248) = rate(:,:,   248) * m(:,:)
      rate(:,:,   249) = rate(:,:,   249) * m(:,:)
      rate(:,:,   250) = rate(:,:,   250) * m(:,:)
      rate(:,:,   251) = rate(:,:,   251) * m(:,:)
      rate(:,:,   252) = rate(:,:,   252) * m(:,:)
      rate(:,:,   253) = rate(:,:,   253) * m(:,:)
      rate(:,:,   254) = rate(:,:,   254) * m(:,:)
      rate(:,:,   255) = rate(:,:,   255) * m(:,:)
      rate(:,:,   256) = rate(:,:,   256) * m(:,:)
      rate(:,:,   257) = rate(:,:,   257) * m(:,:)
      rate(:,:,   258) = rate(:,:,   258) * m(:,:)
      rate(:,:,   259) = rate(:,:,   259) * m(:,:)
      rate(:,:,   260) = rate(:,:,   260) * m(:,:)
      rate(:,:,   261) = rate(:,:,   261) * m(:,:)
      rate(:,:,   262) = rate(:,:,   262) * m(:,:)
      rate(:,:,   263) = rate(:,:,   263) * m(:,:)
      rate(:,:,   264) = rate(:,:,   264) * m(:,:)
      rate(:,:,   265) = rate(:,:,   265) * m(:,:)
      rate(:,:,   266) = rate(:,:,   266) * m(:,:)
      rate(:,:,   267) = rate(:,:,   267) * m(:,:)
      rate(:,:,   268) = rate(:,:,   268) * m(:,:)
      rate(:,:,   269) = rate(:,:,   269) * m(:,:)
      rate(:,:,   271) = rate(:,:,   271) * m(:,:)
      rate(:,:,   272) = rate(:,:,   272) * m(:,:)
      rate(:,:,   273) = rate(:,:,   273) * m(:,:)
      rate(:,:,   274) = rate(:,:,   274) * m(:,:)
      rate(:,:,   275) = rate(:,:,   275) * m(:,:)
      rate(:,:,   276) = rate(:,:,   276) * m(:,:)
      rate(:,:,   277) = rate(:,:,   277) * m(:,:)
      rate(:,:,   278) = rate(:,:,   278) * m(:,:)
      rate(:,:,   279) = rate(:,:,   279) * m(:,:)
      rate(:,:,   280) = rate(:,:,   280) * m(:,:)
      rate(:,:,   281) = rate(:,:,   281) * m(:,:)
      rate(:,:,   282) = rate(:,:,   282) * m(:,:)
      rate(:,:,   283) = rate(:,:,   283) * m(:,:)
      rate(:,:,   284) = rate(:,:,   284) * m(:,:)
      rate(:,:,   285) = rate(:,:,   285) * m(:,:)
      rate(:,:,   286) = rate(:,:,   286) * m(:,:)
      rate(:,:,   287) = rate(:,:,   287) * m(:,:)
      rate(:,:,   288) = rate(:,:,   288) * m(:,:)
      rate(:,:,   289) = rate(:,:,   289) * m(:,:)
      rate(:,:,   290) = rate(:,:,   290) * m(:,:)
      rate(:,:,   291) = rate(:,:,   291) * m(:,:)
      rate(:,:,   292) = rate(:,:,   292) * m(:,:)
      rate(:,:,   293) = rate(:,:,   293) * m(:,:)
      rate(:,:,   294) = rate(:,:,   294) * m(:,:)
      rate(:,:,   295) = rate(:,:,   295) * m(:,:)
      rate(:,:,   296) = rate(:,:,   296) * m(:,:)
      rate(:,:,   297) = rate(:,:,   297) * m(:,:)
      rate(:,:,   298) = rate(:,:,   298) * m(:,:)
      rate(:,:,   299) = rate(:,:,   299) * m(:,:)
      rate(:,:,   302) = rate(:,:,   302) * m(:,:)
      rate(:,:,   303) = rate(:,:,   303) * m(:,:)
      rate(:,:,   304) = rate(:,:,   304) * m(:,:)
      rate(:,:,   305) = rate(:,:,   305) * m(:,:)
      rate(:,:,   306) = rate(:,:,   306) * m(:,:)
      rate(:,:,   308) = rate(:,:,   308) * m(:,:)
      rate(:,:,   309) = rate(:,:,   309) * m(:,:)
      rate(:,:,   310) = rate(:,:,   310) * m(:,:)
      rate(:,:,   311) = rate(:,:,   311) * m(:,:)
      rate(:,:,   312) = rate(:,:,   312) * m(:,:)
      rate(:,:,   313) = rate(:,:,   313) * m(:,:)
      rate(:,:,   314) = rate(:,:,   314) * m(:,:)
      rate(:,:,   315) = rate(:,:,   315) * m(:,:)
      rate(:,:,   316) = rate(:,:,   316) * m(:,:)
      rate(:,:,   317) = rate(:,:,   317) * m(:,:)
      rate(:,:,   318) = rate(:,:,   318) * m(:,:)
      rate(:,:,   319) = rate(:,:,   319) * m(:,:)
      rate(:,:,   320) = rate(:,:,   320) * m(:,:)
      rate(:,:,   321) = rate(:,:,   321) * m(:,:)
      rate(:,:,   322) = rate(:,:,   322) * m(:,:)
      rate(:,:,   323) = rate(:,:,   323) * m(:,:)
      rate(:,:,   324) = rate(:,:,   324) * m(:,:)
      rate(:,:,   325) = rate(:,:,   325) * m(:,:)
      rate(:,:,   326) = rate(:,:,   326) * m(:,:)
      rate(:,:,   327) = rate(:,:,   327) * m(:,:)
      rate(:,:,   328) = rate(:,:,   328) * m(:,:)
      rate(:,:,   329) = rate(:,:,   329) * m(:,:)
      rate(:,:,   330) = rate(:,:,   330) * m(:,:)
      rate(:,:,   331) = rate(:,:,   331) * m(:,:)
      rate(:,:,   332) = rate(:,:,   332) * m(:,:)
      rate(:,:,   333) = rate(:,:,   333) * m(:,:)
      rate(:,:,   334) = rate(:,:,   334) * m(:,:)
      rate(:,:,   335) = rate(:,:,   335) * m(:,:)
      rate(:,:,   336) = rate(:,:,   336) * m(:,:)
      rate(:,:,   337) = rate(:,:,   337) * m(:,:)
      rate(:,:,   338) = rate(:,:,   338) * m(:,:)
      rate(:,:,   339) = rate(:,:,   339) * m(:,:)
      rate(:,:,   340) = rate(:,:,   340) * m(:,:)
      rate(:,:,   341) = rate(:,:,   341) * m(:,:)
      rate(:,:,   342) = rate(:,:,   342) * m(:,:)
      rate(:,:,   343) = rate(:,:,   343) * m(:,:)
      rate(:,:,   344) = rate(:,:,   344) * m(:,:)
      rate(:,:,   345) = rate(:,:,   345) * m(:,:)
      rate(:,:,   346) = rate(:,:,   346) * m(:,:)
      rate(:,:,   347) = rate(:,:,   347) * m(:,:)
      rate(:,:,   348) = rate(:,:,   348) * m(:,:)
      rate(:,:,   349) = rate(:,:,   349) * m(:,:)
      rate(:,:,   350) = rate(:,:,   350) * m(:,:)
      rate(:,:,   351) = rate(:,:,   351) * m(:,:)
      rate(:,:,   352) = rate(:,:,   352) * m(:,:)
      rate(:,:,   353) = rate(:,:,   353) * m(:,:)
      rate(:,:,   354) = rate(:,:,   354) * m(:,:)
      rate(:,:,   355) = rate(:,:,   355) * m(:,:)
      rate(:,:,   357) = rate(:,:,   357) * m(:,:)
      rate(:,:,   358) = rate(:,:,   358) * m(:,:)
      rate(:,:,   359) = rate(:,:,   359) * m(:,:)
      rate(:,:,   360) = rate(:,:,   360) * m(:,:)
      rate(:,:,   361) = rate(:,:,   361) * m(:,:)
      rate(:,:,   362) = rate(:,:,   362) * m(:,:)
      rate(:,:,   363) = rate(:,:,   363) * m(:,:)
      rate(:,:,   364) = rate(:,:,   364) * m(:,:)
      rate(:,:,   365) = rate(:,:,   365) * m(:,:)
      rate(:,:,   366) = rate(:,:,   366) * m(:,:)
      rate(:,:,   367) = rate(:,:,   367) * m(:,:)
      rate(:,:,   368) = rate(:,:,   368) * m(:,:)
      rate(:,:,   369) = rate(:,:,   369) * m(:,:)
      rate(:,:,   370) = rate(:,:,   370) * m(:,:)
      rate(:,:,   371) = rate(:,:,   371) * m(:,:)
      rate(:,:,   372) = rate(:,:,   372) * m(:,:)
      rate(:,:,   373) = rate(:,:,   373) * m(:,:)
      rate(:,:,   375) = rate(:,:,   375) * m(:,:)
      rate(:,:,   376) = rate(:,:,   376) * m(:,:)
      rate(:,:,   377) = rate(:,:,   377) * m(:,:)
      rate(:,:,   378) = rate(:,:,   378) * m(:,:)
      rate(:,:,   379) = rate(:,:,   379) * m(:,:)
      rate(:,:,   380) = rate(:,:,   380) * m(:,:)
      rate(:,:,   381) = rate(:,:,   381) * m(:,:)
      rate(:,:,   382) = rate(:,:,   382) * m(:,:)
      rate(:,:,   383) = rate(:,:,   383) * m(:,:)
      rate(:,:,   384) = rate(:,:,   384) * m(:,:)
      rate(:,:,   385) = rate(:,:,   385) * m(:,:)
      rate(:,:,   386) = rate(:,:,   386) * m(:,:)
      rate(:,:,   387) = rate(:,:,   387) * m(:,:)
      rate(:,:,   388) = rate(:,:,   388) * m(:,:)
      rate(:,:,   389) = rate(:,:,   389) * m(:,:)
      rate(:,:,   390) = rate(:,:,   390) * m(:,:)
      rate(:,:,   391) = rate(:,:,   391) * m(:,:)
      rate(:,:,   392) = rate(:,:,   392) * m(:,:)
      rate(:,:,   393) = rate(:,:,   393) * m(:,:)
      rate(:,:,   394) = rate(:,:,   394) * m(:,:)
      rate(:,:,   395) = rate(:,:,   395) * m(:,:)
      rate(:,:,   396) = rate(:,:,   396) * m(:,:)
      rate(:,:,   397) = rate(:,:,   397) * m(:,:)
      rate(:,:,   398) = rate(:,:,   398) * m(:,:)
      rate(:,:,   399) = rate(:,:,   399) * m(:,:)
      rate(:,:,   400) = rate(:,:,   400) * m(:,:)
      rate(:,:,   401) = rate(:,:,   401) * m(:,:)
      rate(:,:,   402) = rate(:,:,   402) * m(:,:)
      rate(:,:,   403) = rate(:,:,   403) * m(:,:)
      rate(:,:,   404) = rate(:,:,   404) * m(:,:)
      rate(:,:,   405) = rate(:,:,   405) * m(:,:)
      rate(:,:,   406) = rate(:,:,   406) * m(:,:)
      rate(:,:,   407) = rate(:,:,   407) * m(:,:)
      rate(:,:,   408) = rate(:,:,   408) * m(:,:)
      rate(:,:,   409) = rate(:,:,   409) * m(:,:)
      rate(:,:,   410) = rate(:,:,   410) * m(:,:)
      rate(:,:,   411) = rate(:,:,   411) * m(:,:)
      rate(:,:,   412) = rate(:,:,   412) * m(:,:)
      rate(:,:,   413) = rate(:,:,   413) * m(:,:)
      rate(:,:,   414) = rate(:,:,   414) * m(:,:)
      rate(:,:,   415) = rate(:,:,   415) * m(:,:)
      rate(:,:,   416) = rate(:,:,   416) * m(:,:)
      rate(:,:,   417) = rate(:,:,   417) * m(:,:)
      rate(:,:,   418) = rate(:,:,   418) * m(:,:)
      rate(:,:,   419) = rate(:,:,   419) * m(:,:)
      rate(:,:,   420) = rate(:,:,   420) * m(:,:)
      rate(:,:,   421) = rate(:,:,   421) * m(:,:)
      rate(:,:,   422) = rate(:,:,   422) * m(:,:)
      rate(:,:,   423) = rate(:,:,   423) * m(:,:)
      rate(:,:,   424) = rate(:,:,   424) * m(:,:)
      rate(:,:,   425) = rate(:,:,   425) * m(:,:)
      rate(:,:,   426) = rate(:,:,   426) * m(:,:)
      rate(:,:,   427) = rate(:,:,   427) * m(:,:)
      rate(:,:,   428) = rate(:,:,   428) * m(:,:)
      rate(:,:,   429) = rate(:,:,   429) * m(:,:)
      rate(:,:,   430) = rate(:,:,   430) * m(:,:)
      rate(:,:,   432) = rate(:,:,   432) * m(:,:)
      rate(:,:,   433) = rate(:,:,   433) * m(:,:)
      rate(:,:,   434) = rate(:,:,   434) * m(:,:)
      rate(:,:,   435) = rate(:,:,   435) * m(:,:)
      rate(:,:,   436) = rate(:,:,   436) * m(:,:)
      rate(:,:,   437) = rate(:,:,   437) * m(:,:)
      rate(:,:,   438) = rate(:,:,   438) * m(:,:)
      rate(:,:,   439) = rate(:,:,   439) * m(:,:)
      rate(:,:,   440) = rate(:,:,   440) * m(:,:)
      rate(:,:,   441) = rate(:,:,   441) * m(:,:)
      rate(:,:,   442) = rate(:,:,   442) * m(:,:)
      rate(:,:,   443) = rate(:,:,   443) * m(:,:)
      rate(:,:,   444) = rate(:,:,   444) * m(:,:)
      rate(:,:,   445) = rate(:,:,   445) * m(:,:)
      rate(:,:,   446) = rate(:,:,   446) * m(:,:)
      rate(:,:,   447) = rate(:,:,   447) * m(:,:)
      rate(:,:,   448) = rate(:,:,   448) * m(:,:)
      rate(:,:,   449) = rate(:,:,   449) * m(:,:)
      rate(:,:,   450) = rate(:,:,   450) * m(:,:)
      rate(:,:,   451) = rate(:,:,   451) * m(:,:)
      rate(:,:,   452) = rate(:,:,   452) * m(:,:)
      rate(:,:,   453) = rate(:,:,   453) * m(:,:)
      rate(:,:,   454) = rate(:,:,   454) * m(:,:)
      rate(:,:,   455) = rate(:,:,   455) * m(:,:)
      rate(:,:,   456) = rate(:,:,   456) * m(:,:)
      rate(:,:,   457) = rate(:,:,   457) * m(:,:)
      rate(:,:,   458) = rate(:,:,   458) * m(:,:)
      rate(:,:,   459) = rate(:,:,   459) * m(:,:)
      rate(:,:,   460) = rate(:,:,   460) * m(:,:)
      rate(:,:,   461) = rate(:,:,   461) * m(:,:)
      rate(:,:,   462) = rate(:,:,   462) * m(:,:)
      rate(:,:,   463) = rate(:,:,   463) * m(:,:)
      rate(:,:,   464) = rate(:,:,   464) * m(:,:)
      rate(:,:,   465) = rate(:,:,   465) * m(:,:)
      rate(:,:,   466) = rate(:,:,   466) * m(:,:)
      rate(:,:,   468) = rate(:,:,   468) * m(:,:)
      rate(:,:,   469) = rate(:,:,   469) * m(:,:)
      rate(:,:,   470) = rate(:,:,   470) * m(:,:)
      rate(:,:,   471) = rate(:,:,   471) * m(:,:)
      rate(:,:,   472) = rate(:,:,   472) * m(:,:)
      rate(:,:,   473) = rate(:,:,   473) * m(:,:)
      rate(:,:,   475) = rate(:,:,   475) * m(:,:)
      rate(:,:,   476) = rate(:,:,   476) * m(:,:)
      rate(:,:,   477) = rate(:,:,   477) * m(:,:)
      rate(:,:,   478) = rate(:,:,   478) * m(:,:)
      rate(:,:,   479) = rate(:,:,   479) * m(:,:)
      rate(:,:,   480) = rate(:,:,   480) * m(:,:)
      rate(:,:,   494) = rate(:,:,   494) * m(:,:)
      rate(:,:,   495) = rate(:,:,   495) * m(:,:)
      rate(:,:,   496) = rate(:,:,   496) * m(:,:)
      rate(:,:,   497) = rate(:,:,   497) * m(:,:)
      rate(:,:,   498) = rate(:,:,   498) * m(:,:)
      rate(:,:,   499) = rate(:,:,   499) * m(:,:)
      rate(:,:,   500) = rate(:,:,   500) * m(:,:)
      rate(:,:,   501) = rate(:,:,   501) * m(:,:)
      rate(:,:,   502) = rate(:,:,   502) * m(:,:)
      rate(:,:,   503) = rate(:,:,   503) * m(:,:)
      rate(:,:,   504) = rate(:,:,   504) * m(:,:)
      rate(:,:,   505) = rate(:,:,   505) * m(:,:)
      rate(:,:,   506) = rate(:,:,   506) * m(:,:)
      rate(:,:,   507) = rate(:,:,   507) * m(:,:)
      rate(:,:,   508) = rate(:,:,   508) * m(:,:)
      rate(:,:,   509) = rate(:,:,   509) * m(:,:)
      rate(:,:,   510) = rate(:,:,   510) * m(:,:)
      rate(:,:,   511) = rate(:,:,   511) * m(:,:)
      rate(:,:,   512) = rate(:,:,   512) * m(:,:)
      rate(:,:,   513) = rate(:,:,   513) * m(:,:)
      rate(:,:,   514) = rate(:,:,   514) * m(:,:)
      rate(:,:,   515) = rate(:,:,   515) * m(:,:)
      rate(:,:,   516) = rate(:,:,   516) * m(:,:)
      rate(:,:,   517) = rate(:,:,   517) * m(:,:)
      rate(:,:,   518) = rate(:,:,   518) * m(:,:)
      rate(:,:,   520) = rate(:,:,   520) * m(:,:)
      rate(:,:,   521) = rate(:,:,   521) * m(:,:)
      rate(:,:,   522) = rate(:,:,   522) * m(:,:)
      rate(:,:,   524) = rate(:,:,   524) * m(:,:)
      rate(:,:,   529) = rate(:,:,   529) * m(:,:)
      rate(:,:,   530) = rate(:,:,   530) * m(:,:)
      rate(:,:,   531) = rate(:,:,   531) * m(:,:)
      rate(:,:,   534) = rate(:,:,   534) * m(:,:)
      rate(:,:,   535) = rate(:,:,   535) * m(:,:)
      rate(:,:,   536) = rate(:,:,   536) * m(:,:)
      rate(:,:,   539) = rate(:,:,   539) * m(:,:)

      end subroutine adjrxt

      end module mo_adjrxt

      module mo_phtadj

      private
      public :: phtadj

      contains

      subroutine phtadj( p_rate, inv, m, ncol, nlev )

      use chem_mods,    only : nfs, phtcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!--------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol, nlev
      real(r8), intent(in)    :: inv(ncol,nlev,max(1,nfs))
      real(r8), intent(in)    :: m(ncol,nlev)
      real(r8), intent(inout) :: p_rate(ncol,nlev,max(1,phtcnt))

!--------------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------------
      integer  ::  k
      real(r8) ::  im(ncol,nlev)

      do k = 1,nlev
         im(:ncol,k) = 1._r8 / m(:ncol,k)
         p_rate(:,k,  5) = p_rate(:,k,  5)  * inv(:,k, 2) * im(:,k)
         p_rate(:,k,  6) = p_rate(:,k,  6)  * inv(:,k, 2) * im(:,k)
      end do

      end subroutine phtadj

      end module mo_phtadj

      module mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat

      use chem_mods,     only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass, crb_mass
      use chem_mods,     only : diag_map
      use chem_mods,     only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods,     only : pht_alias_lst, pht_alias_mult
      use chem_mods,     only : extfrc_lst, inv_lst, slvd_lst
      use chem_mods,     only : enthalpy_cnt, cph_enthalpy, cph_rid, num_rnts, rxntot
      use cam_abortutils,only : endrun
      use mo_tracname,   only : solsym
      use chem_mods,     only : frc_from_dataset
      use chem_mods,     only : is_scalar, is_vector
      use shr_kind_mod,  only : r8 => shr_kind_r8
      use cam_logfile,   only : iulog

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios

      is_scalar = .false.
      is_vector = .true.

      clscnt(:) = (/      2,     0,     0,   227,     0 /)

      cls_rxt_cnt(:,1) = (/      9,     0,     0,     2 /)
      cls_rxt_cnt(:,4) = (/      2,   165,   375,   227 /)

      solsym(:229) = (/ 'ALKNIT          ','ALKOOH          ','AOA_NH          ','bc_a1           ','bc_a4           ', &
                        'BCARY           ','BENZENE         ','BENZOOH         ','BEPOMUC         ','BIGALD          ', &
                        'BIGALD1         ','BIGALD2         ','BIGALD3         ','BIGALD4         ','BIGALK          ', &
                        'BIGENE          ','BR              ','BRCL            ','BRO             ','BRONO2          ', &
                        'BRY             ','BZALD           ','BZOOH           ','C2H2            ','C2H4            ', &
                        'C2H5OH          ','C2H5OOH         ','C2H6            ','C3H6            ','C3H7OOH         ', &
                        'C3H8            ','C6H5OOH         ','CCL4            ','CF2CLBR         ','CF3BR           ', &
                        'CFC11           ','CFC113          ','CFC114          ','CFC115          ','CFC12           ', &
                        'CH2BR2          ','CH2O            ','CH3BR           ','CH3CCL3         ','CH3CHO          ', &
                        'CH3CL           ','CH3CN           ','CH3COCH3        ','CH3COCHO        ','CH3COOH         ', &
                        'CH3COOOH        ','CH3OH           ','CH3OOH          ','CH4             ','CHBR3           ', &
                        'CL              ','CL2             ','CL2O2           ','CLO             ','CLONO2          ', &
                        'CLY             ','CO              ','CO2             ','COF2            ','COFCL           ', &
                        'CRESOL          ','DMS             ','dst_a1          ','dst_a2          ','dst_a3          ', &
                        'E90             ','EOOH            ','F               ','GLYALD          ','GLYOXAL         ', &
                        'H               ','H2              ','H2402           ','H2O2            ','H2SO4           ', &
                        'HBR             ','HCFC141B        ','HCFC142B        ','HCFC22          ','HCL             ', &
                        'HCN             ','HCOOH           ','HF              ','HNO3            ','HO2NO2          ', &
                        'HOBR            ','HOCL            ','HONITR          ','HPALD           ','HYAC            ', &
                        'HYDRALD         ','IEPOX           ','ISOP            ','ISOPNITA        ','ISOPNITB        ', &
                        'ISOPNO3         ','ISOPNOOH        ','ISOPOOH         ','IVOC            ','MACR            ', &
                        'MACROOH         ','MEK             ','MEKOOH          ','MPAN            ','MTERP           ', &
                        'MVK             ','N               ','N2O             ','N2O5            ','NC4CH2OH        ', &
                        'NC4CHO          ','ncl_a1          ','ncl_a2          ','ncl_a3          ','NH3             ', &
                        'NH4             ','NH_5            ','NH_50           ','NO              ','NO2             ', &
                        'NO3             ','NOA             ','NTERPOOH        ','num_a1          ','num_a2          ', &
                        'num_a3          ','num_a4          ','O               ','O3              ','O3S             ', &
                        'OCLO            ','OCS             ','ONITR           ','PAN             ','PBZNIT          ', &
                        'PHENO           ','PHENOL          ','PHENOOH         ','pom_a1          ','pom_a4          ', &
                        'POOH            ','ROOH            ','S               ','SF6             ','SO              ', &
                        'SO2             ','SO3             ','so4_a1          ','so4_a2          ','so4_a3          ', &
                        'soa1_a1         ','soa1_a2         ','soa2_a1         ','soa2_a2         ','soa3_a1         ', &
                        'soa3_a2         ','soa4_a1         ','soa4_a2         ','soa5_a1         ','soa5_a2         ', &
                        'SOAG0           ','SOAG1           ','SOAG2           ','SOAG3           ','SOAG4           ', &
                        'ST80_25         ','SVOC            ','TEPOMUC         ','TERP2OOH        ','TERPNIT         ', &
                        'TERPOOH         ','TERPROD1        ','TERPROD2        ','TOLOOH          ','TOLUENE         ', &
                        'XOOH            ','XYLENES         ','XYLENOOH        ','XYLOL           ','XYLOLOOH        ', &
                        'NHDEP           ','NDEP            ','ACBZO2          ','ALKO2           ','BCARYO2VBS      ', &
                        'BENZO2          ','BENZO2VBS       ','BZOO            ','C2H5O2          ','C3H7O2          ', &
                        'C6H5O2          ','CH3CO3          ','CH3O2           ','DICARBO2        ','ENEO2           ', &
                        'EO              ','EO2             ','HO2             ','HOCH2OO         ','ISOPAO2         ', &
                        'ISOPBO2         ','ISOPO2VBS       ','IVOCO2VBS       ','MACRO2          ','MALO2           ', &
                        'MCO3            ','MDIALO2         ','MEKO2           ','MTERPO2VBS      ','NTERPO2         ', &
                        'O1D             ','OH              ','PHENO2          ','PO2             ','RO2             ', &
                        'TERP2O2         ','TERPO2          ','TOLO2           ','TOLUO2VBS       ','XO2             ', &
                        'XYLENO2         ','XYLEO2VBS       ','XYLOLO2         ','H2O             ' /)

      adv_mass(:229) = (/   133.141340_r8,   104.142600_r8,    28.010400_r8,    12.011000_r8,    12.011000_r8, &
                            204.342600_r8,    78.110400_r8,   160.122200_r8,   126.108600_r8,    98.098200_r8, &
                             84.072400_r8,    98.098200_r8,    98.098200_r8,   112.124000_r8,    72.143800_r8, &
                             56.103200_r8,    79.904000_r8,   115.356700_r8,    95.903400_r8,   141.908940_r8, &
                             99.716850_r8,   106.120800_r8,   124.135000_r8,    26.036800_r8,    28.051600_r8, &
                             46.065800_r8,    62.065200_r8,    30.066400_r8,    42.077400_r8,    76.091000_r8, &
                             44.092200_r8,   110.109200_r8,   153.821800_r8,   165.364506_r8,   148.910210_r8, &
                            137.367503_r8,   187.375310_r8,   170.921013_r8,   154.466716_r8,   120.913206_r8, &
                            173.833800_r8,    30.025200_r8,    94.937200_r8,   133.402300_r8,    44.051000_r8, &
                             50.485900_r8,    41.050940_r8,    58.076800_r8,    72.061400_r8,    60.050400_r8, &
                             76.049800_r8,    32.040000_r8,    48.039400_r8,    16.040600_r8,   252.730400_r8, &
                             35.452700_r8,    70.905400_r8,   102.904200_r8,    51.452100_r8,    97.457640_r8, &
                            100.916850_r8,    28.010400_r8,    44.009800_r8,    66.007206_r8,    82.461503_r8, &
                            108.135600_r8,    62.132400_r8,   135.064039_r8,   135.064039_r8,   135.064039_r8, &
                             28.010400_r8,    78.064600_r8,    18.998403_r8,    60.050400_r8,    58.035600_r8, &
                              1.007400_r8,     2.014800_r8,   259.823613_r8,    34.013600_r8,    98.078400_r8, &
                             80.911400_r8,   116.948003_r8,   100.493706_r8,    86.467906_r8,    36.460100_r8, &
                             27.025140_r8,    46.024600_r8,    20.005803_r8,    63.012340_r8,    79.011740_r8, &
                             96.910800_r8,    52.459500_r8,   135.114940_r8,   116.112400_r8,    74.076200_r8, &
                            100.113000_r8,   118.127200_r8,    68.114200_r8,   147.125940_r8,   147.125940_r8, &
                            162.117940_r8,   163.125340_r8,   118.127200_r8,   184.350200_r8,    70.087800_r8, &
                            120.100800_r8,    72.102600_r8,   104.101400_r8,   147.084740_r8,   136.228400_r8, &
                             70.087800_r8,    14.006740_r8,    44.012880_r8,   108.010480_r8,   147.125940_r8, &
                            145.111140_r8,    58.442468_r8,    58.442468_r8,    58.442468_r8,    17.028940_r8, &
                             18.036340_r8,    28.010400_r8,    28.010400_r8,    30.006140_r8,    46.005540_r8, &
                             62.004940_r8,   119.074340_r8,   231.239540_r8,     1.007400_r8,     1.007400_r8, &
                              1.007400_r8,     1.007400_r8,    15.999400_r8,    47.998200_r8,    47.998200_r8, &
                             67.451500_r8,    60.076400_r8,   133.100140_r8,   121.047940_r8,   183.117740_r8, &
                             93.102400_r8,    94.109800_r8,   176.121600_r8,    12.011000_r8,    12.011000_r8, &
                             92.090400_r8,    90.075600_r8,    32.066000_r8,   146.056419_r8,    48.065400_r8, &
                             64.064800_r8,    80.064200_r8,   115.107340_r8,   115.107340_r8,   115.107340_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                             28.010400_r8,   310.582400_r8,   140.134400_r8,   200.226000_r8,   215.240140_r8, &
                            186.241400_r8,   168.227200_r8,   154.201400_r8,   174.148000_r8,    92.136200_r8, &
                            150.126000_r8,   106.162000_r8,   188.173800_r8,   122.161400_r8,   204.173200_r8, &
                             14.006740_r8,    14.006740_r8,   137.112200_r8,   103.135200_r8,   253.348200_r8, &
                            159.114800_r8,   159.114800_r8,   123.127600_r8,    61.057800_r8,    75.083600_r8, &
                            109.101800_r8,    75.042400_r8,    47.032000_r8,   129.089600_r8,   105.108800_r8, &
                             61.057800_r8,    77.057200_r8,    33.006200_r8,    63.031400_r8,   117.119800_r8, &
                            117.119800_r8,   117.119800_r8,   233.355800_r8,   119.093400_r8,   115.063800_r8, &
                            101.079200_r8,   117.078600_r8,   103.094000_r8,   185.234000_r8,   230.232140_r8, &
                             15.999400_r8,    17.006800_r8,   175.114200_r8,    91.083000_r8,    89.068200_r8, &
                            199.218600_r8,   185.234000_r8,   173.140600_r8,   173.140600_r8,   149.118600_r8, &
                            187.166400_r8,   187.166400_r8,   203.165800_r8,    18.014200_r8 /)

      crb_mass(:229) = (/    60.055000_r8,    60.055000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                            180.165000_r8,    72.066000_r8,    72.066000_r8,    72.066000_r8,    60.055000_r8, &
                             48.044000_r8,    60.055000_r8,    60.055000_r8,    72.066000_r8,    60.055000_r8, &
                             48.044000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    84.077000_r8,    84.077000_r8,    24.022000_r8,    24.022000_r8, &
                             24.022000_r8,    24.022000_r8,    24.022000_r8,    36.033000_r8,    36.033000_r8, &
                             36.033000_r8,    72.066000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                             12.011000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8,    12.011000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    24.022000_r8,    24.022000_r8, &
                             12.011000_r8,    24.022000_r8,    36.033000_r8,    36.033000_r8,    24.022000_r8, &
                             24.022000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                             84.077000_r8,    24.022000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,    24.022000_r8,     0.000000_r8,    24.022000_r8,    24.022000_r8, &
                              0.000000_r8,     0.000000_r8,    24.022000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    24.022000_r8,    24.022000_r8,    12.011000_r8,     0.000000_r8, &
                             12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,    48.044000_r8,    60.055000_r8,    36.033000_r8, &
                             60.055000_r8,    60.055000_r8,    60.055000_r8,    60.055000_r8,    60.055000_r8, &
                             60.055000_r8,    60.055000_r8,    60.055000_r8,   156.143000_r8,    48.044000_r8, &
                             48.044000_r8,    48.044000_r8,    48.044000_r8,    48.044000_r8,   120.110000_r8, &
                             48.044000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,    60.055000_r8, &
                             60.055000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    36.033000_r8,   120.110000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    12.011000_r8,    48.044000_r8,    24.022000_r8,    84.077000_r8, &
                             72.066000_r8,    72.066000_r8,    72.066000_r8,    12.011000_r8,    12.011000_r8, &
                             36.033000_r8,    36.033000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                             12.011000_r8,   264.242000_r8,    84.077000_r8,   120.110000_r8,   120.110000_r8, &
                            120.110000_r8,   120.110000_r8,   108.099000_r8,    84.077000_r8,    84.077000_r8, &
                             60.055000_r8,    96.088000_r8,    96.088000_r8,    96.088000_r8,    96.088000_r8, &
                              0.000000_r8,     0.000000_r8,    84.077000_r8,    60.055000_r8,   180.165000_r8, &
                             72.066000_r8,    72.066000_r8,    84.077000_r8,    24.022000_r8,    36.033000_r8, &
                             72.066000_r8,    24.022000_r8,    12.011000_r8,    60.055000_r8,    48.044000_r8, &
                             24.022000_r8,    24.022000_r8,     0.000000_r8,    12.011000_r8,    60.055000_r8, &
                             60.055000_r8,    60.055000_r8,   156.143000_r8,    48.044000_r8,    48.044000_r8, &
                             48.044000_r8,    48.044000_r8,    48.044000_r8,   120.110000_r8,   120.110000_r8, &
                              0.000000_r8,     0.000000_r8,    72.066000_r8,    36.033000_r8,    36.033000_r8, &
                            120.110000_r8,   120.110000_r8,    84.077000_r8,    84.077000_r8,    60.055000_r8, &
                             96.088000_r8,    96.088000_r8,    96.088000_r8,     0.000000_r8 /)

      fix_mass(:  3) = (/ 0.00000000_r8, 31.9988000_r8, 28.0134800_r8 /)

      clsmap(:  2,1) = (/  186, 187 /)
      clsmap(:227,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                            21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                            31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                            41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                            51,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                            61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                            71,  72,  73,  74,  75,  76,  77,  78,  79,  80, &
                            81,  82,  83,  84,  85,  86,  87,  88,  89,  90, &
                            91,  92,  93,  94,  95,  96,  97,  98,  99, 100, &
                           101, 102, 103, 104, 105, 106, 107, 108, 109, 110, &
                           111, 112, 113, 114, 115, 116, 117, 118, 119, 120, &
                           121, 122, 123, 124, 125, 126, 127, 128, 129, 130, &
                           131, 132, 133, 134, 135, 136, 137, 138, 139, 140, &
                           141, 142, 143, 144, 145, 146, 147, 148, 149, 150, &
                           151, 152, 153, 154, 155, 156, 157, 158, 159, 160, &
                           161, 162, 163, 164, 165, 166, 167, 168, 169, 170, &
                           171, 172, 173, 174, 175, 176, 177, 178, 179, 180, &
                           181, 182, 183, 184, 185, 188, 189, 190, 191, 192, &
                           193, 194, 195, 196, 197, 198, 199, 200, 201, 202, &
                           203, 204, 205, 206, 207, 208, 209, 210, 211, 212, &
                           213, 214, 215, 216, 217, 218, 219, 220, 221, 222, &
                           223, 224, 225, 226, 227, 228, 229 /)

      permute(:227,4) = (/  154, 151,   1,   2,   3, 185,  71, 121,  72, 114, &
                            127,  96, 145, 105,  86, 110, 209,  87, 225, 139, &
                              4,  90, 108,  98, 140,  92, 109,  99, 186, 120, &
                             57,  93,  54,  66,  67,  58,  68,  59,  69,  60, &
                            129, 213, 146,  61, 190, 112,  55, 180, 200, 156, &
                            148, 166, 115, 210, 125, 220,  70,  52, 219, 177, &
                              5, 192, 168,  85,  83,  77, 100,   6,   7,   8, &
                              9,  62, 175, 191, 184, 212, 208,  56, 147,  63, &
                            169,  82,  89, 101, 223,  74, 181,  97, 211, 118, &
                            165, 171, 196,  84, 195, 104,  64, 173, 144, 141, &
                            198, 117, 157,  48, 199, 102, 134, 103, 143, 182, &
                            205, 132,  75,  95, 113, 189,  10,  11,  12,  53, &
                             13,  14,  15, 217, 224, 216, 174, 116,  16,  17, &
                             18,  19, 226, 222,  20, 106, 111,  88, 137,  65, &
                            128,  73, 107,  21,  22, 138, 119, 135,  23, 201, &
                            172,  91,  24,  25,  26,  27,  28,  29,  30,  31, &
                             32,  33,  34,  35,  36,  37,  38,  39,  40,  41, &
                             42,  43,  78, 152, 149, 130, 183, 188, 153,  76, &
                             79,  80, 158,  81, 122, 136, 178,  44, 131,  45, &
                            123, 170, 167, 150, 207, 221, 163, 142,  94, 159, &
                            218, 124, 202, 203,  46,  47, 204, 160, 206, 176, &
                            155,  49, 187, 214, 215, 126, 164, 194, 193, 179, &
                            161,  50, 197, 162,  51, 133, 227 /)

      diag_map(:227) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                            21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                            31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                            41,  42,  48,  54,  60,  66,  72,  74,  80,  86, &
                            92,  93,  96,  99, 102, 105, 109, 113, 117, 121, &
                           125, 129, 132, 135, 138, 141, 146, 151, 156, 161, &
                           164, 170, 174, 179, 183, 187, 194, 199, 203, 208, &
                           216, 221, 226, 230, 235, 238, 241, 244, 248, 253, &
                           257, 261, 265, 269, 273, 279, 282, 285, 291, 297, &
                           304, 310, 315, 320, 323, 328, 333, 339, 344, 349, &
                           357, 365, 373, 379, 385, 391, 397, 403, 409, 415, &
                           421, 427, 435, 441, 447, 454, 460, 463, 467, 474, &
                           483, 491, 499, 506, 511, 518, 524, 532, 540, 548, &
                           556, 564, 572, 581, 590, 594, 603, 610, 617, 625, &
                           632, 642, 655, 666, 677, 684, 690, 701, 714, 721, &
                           732, 748, 759, 768, 778, 786, 791, 801, 804, 815, &
                           824, 832, 841, 857, 864, 873, 884, 899, 912, 922, &
                           929, 948, 969, 979, 999,1024,1045,1059,1073,1085, &
                          1096,1103,1115,1129,1140,1153,1173,1193,1209,1221, &
                          1232,1256,1288,1311,1332,1354,1386,1401,1415,1430, &
                          1447,1463,1485,1526,1691,1749,1842,1950,1977,2017, &
                          2070,2132,2156,2201,2226,2258,2285 /)

      extfrc_lst(: 14) = (/ 'num_a4          ','pom_a4          ','bc_a4           ','SVOC            ','so4_a1          ', &
                            'so4_a2          ','CO              ','SO2             ','NO2             ','num_a1          ', &
                            'num_a2          ','AOA_NH          ','NO              ','N               ' /)

      frc_from_dataset(: 14) = (/ .true., .true., .true., .true., .true., &
                                  .true., .true., .true., .true., .true., &
                                  .true., .false., .false., .false. /)

      inv_lst(:  3) = (/ 'M               ', 'O2              ', 'N2              ' /)

      slvd_lst(: 41) = (/ 'ACBZO2          ', 'ALKO2           ', 'BCARYO2VBS      ', 'BENZO2          ', 'BENZO2VBS       ', &
                          'BZOO            ', 'C2H5O2          ', 'C3H7O2          ', 'C6H5O2          ', 'CH3CO3          ', &
                          'CH3O2           ', 'DICARBO2        ', 'ENEO2           ', 'EO              ', 'EO2             ', &
                          'HO2             ', 'HOCH2OO         ', 'ISOPAO2         ', 'ISOPBO2         ', 'ISOPO2VBS       ', &
                          'IVOCO2VBS       ', 'MACRO2          ', 'MALO2           ', 'MCO3            ', 'MDIALO2         ', &
                          'MEKO2           ', 'MTERPO2VBS      ', 'NTERPO2         ', 'O1D             ', 'OH              ', &
                          'PHENO2          ', 'PO2             ', 'RO2             ', 'TERP2O2         ', 'TERPO2          ', &
                          'TOLO2           ', 'TOLUO2VBS       ', 'XO2             ', 'XYLENO2         ', 'XYLEO2VBS       ', &
                          'XYLOLO2         ' /)

      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(     1:   200) = (/ 'jh2o_b                          ', 'jh2o_a                          ', &
                                      'jh2o_c                          ', 'jh2o2                           ', &
                                      'jo2_a                           ', 'jo2_b                           ', &
                                      'jo3_a                           ', 'jo3_b                           ', &
                                      'jhno3                           ', 'jho2no2_a                       ', &
                                      'jho2no2_b                       ', 'jn2o                            ', &
                                      'jn2o5_a                         ', 'jn2o5_b                         ', &
                                      'jno                             ', 'jno2                            ', &
                                      'jno3_b                          ', 'jno3_a                          ', &
                                      'jalknit                         ', 'jalkooh                         ', &
                                      'jbenzooh                        ', 'jbepomuc                        ', &
                                      'jbigald                         ', 'jbigald1                        ', &
                                      'jbigald2                        ', 'jbigald3                        ', &
                                      'jbigald4                        ', 'jbzooh                          ', &
                                      'jc2h5ooh                        ', 'jc3h7ooh                        ', &
                                      'jc6h5ooh                        ', 'jch2o_b                         ', &
                                      'jch2o_a                         ', 'jch3cho                         ', &
                                      'jacet                           ', 'jmgly                           ', &
                                      'jch3co3h                        ', 'jch3ooh                         ', &
                                      'jch4_b                          ', 'jch4_a                          ', &
                                      'jco2                            ', 'jeooh                           ', &
                                      'jglyald                         ', 'jglyoxal                        ', &
                                      'jhonitr                         ', 'jhpald                          ', &
                                      'jhyac                           ', 'jisopnooh                       ', &
                                      'jisopooh                        ', 'jmacr_a                         ', &
                                      'jmacr_b                         ', 'jmek                            ', &
                                      'jmekooh                         ', 'jmpan                           ', &
                                      'jmvk                            ', 'jnc4cho                         ', &
                                      'jnoa                            ', 'jnterpooh                       ', &
                                      'jonitr                          ', 'jpan                            ', &
                                      'jphenooh                        ', 'jpooh                           ', &
                                      'jrooh                           ', 'jtepomuc                        ', &
                                      'jterp2ooh                       ', 'jterpnit                        ', &
                                      'jterpooh                        ', 'jterprd1                        ', &
                                      'jterprd2                        ', 'jtolooh                         ', &
                                      'jxooh                           ', 'jxylenooh                       ', &
                                      'jxylolooh                       ', 'jbrcl                           ', &
                                      'jbro                            ', 'jbrono2_b                       ', &
                                      'jbrono2_a                       ', 'jccl4                           ', &
                                      'jcf2clbr                        ', 'jcf3br                          ', &
                                      'jcfcl3                          ', 'jcfc113                         ', &
                                      'jcfc114                         ', 'jcfc115                         ', &
                                      'jcf2cl2                         ', 'jch2br2                         ', &
                                      'jch3br                          ', 'jch3ccl3                        ', &
                                      'jch3cl                          ', 'jchbr3                          ', &
                                      'jcl2                            ', 'jcl2o2                          ', &
                                      'jclo                            ', 'jclono2_a                       ', &
                                      'jclono2_b                       ', 'jcof2                           ', &
                                      'jcofcl                          ', 'jh2402                          ', &
                                      'jhbr                            ', 'jhcfc141b                       ', &
                                      'jhcfc142b                       ', 'jhcfc22                         ', &
                                      'jhcl                            ', 'jhf                             ', &
                                      'jhobr                           ', 'jhocl                           ', &
                                      'joclo                           ', 'jsf6                            ', &
                                      'jh2so4                          ', 'jocs                            ', &
                                      'jso                             ', 'jso2                            ', &
                                      'jso3                            ', 'jsoa1_a1                        ', &
                                      'jsoa1_a2                        ', 'jsoa2_a1                        ', &
                                      'jsoa2_a2                        ', 'jsoa3_a1                        ', &
                                      'jsoa3_a2                        ', 'jsoa4_a1                        ', &
                                      'jsoa4_a2                        ', 'jsoa5_a1                        ', &
                                      'jsoa5_a2                        ', 'E90_tau                         ', &
                                      'O1D_H2                          ', 'O1D_H2O                         ', &
                                      'O1D_N2                          ', 'O1D_O2ab                        ', &
                                      'O1D_O3                          ', 'O1D_O3a                         ', &
                                      'O_O3                            ', 'usr_O_O                         ', &
                                      'usr_O_O2                        ', 'H2_O                            ', &
                                      'H2O2_O                          ', 'H_HO2                           ', &
                                      'H_HO2a                          ', 'H_HO2b                          ', &
                                      'H_O2                            ', 'HO2_O                           ', &
                                      'HO2_O3                          ', 'H_O3                            ', &
                                      'OH_H2                           ', 'OH_H2O2                         ', &
                                      'OH_HO2                          ', 'OH_O                            ', &
                                      'OH_O3                           ', 'OH_OH                           ', &
                                      'OH_OH_M                         ', 'usr_HO2_HO2                     ', &
                                      'HO2NO2_OH                       ', 'N_NO                            ', &
                                      'N_NO2a                          ', 'N_NO2b                          ', &
                                      'N_NO2c                          ', 'N_O2                            ', &
                                      'NO2_O                           ', 'NO2_O3                          ', &
                                      'NO2_O_M                         ', 'NO3_HO2                         ', &
                                      'NO3_NO                          ', 'NO3_O                           ', &
                                      'NO3_OH                          ', 'N_OH                            ', &
                                      'NO_HO2                          ', 'NO_O3                           ', &
                                      'NO_O_M                          ', 'O1D_N2Oa                        ', &
                                      'O1D_N2Ob                        ', 'tag_NO2_HO2                     ', &
                                      'tag_NO2_NO3                     ', 'tag_NO2_OH                      ', &
                                      'usr_HNO3_OH                     ', 'usr_HO2NO2_M                    ', &
                                      'usr_N2O5_M                      ', 'CL_CH2O                         ', &
                                      'CL_CH4                          ', 'CL_H2                           ', &
                                      'CL_H2O2                         ', 'CL_HO2a                         ', &
                                      'CL_HO2b                         ', 'CL_O3                           ', &
                                      'CLO_CH3O2                       ', 'CLO_CLOa                        ', &
                                      'CLO_CLOb                        ', 'CLO_CLOc                        ', &
                                      'CLO_HO2                         ', 'CLO_NO                          ', &
                                      'CLONO2_CL                       ', 'CLO_NO2_M                       ', &
                                      'CLONO2_O                        ', 'CLONO2_OH                       ', &
                                      'CLO_O                           ', 'CLO_OHa                         ', &
                                      'CLO_OHb                         ', 'HCL_O                           ', &
                                      'HCL_OH                          ', 'HOCL_CL                         ', &
                                      'HOCL_O                          ', 'HOCL_OH                         ' /)
      rxt_tag_lst(   201:   400) = (/ 'O1D_CCL4                        ', 'O1D_CF2CLBR                     ', &
                                      'O1D_CFC11                       ', 'O1D_CFC113                      ', &
                                      'O1D_CFC114                      ', 'O1D_CFC115                      ', &
                                      'O1D_CFC12                       ', 'O1D_HCLa                        ', &
                                      'O1D_HCLb                        ', 'tag_CLO_CLO_M                   ', &
                                      'usr_CL2O2_M                     ', 'BR_CH2O                         ', &
                                      'BR_HO2                          ', 'BR_O3                           ', &
                                      'BRO_BRO                         ', 'BRO_CLOa                        ', &
                                      'BRO_CLOb                        ', 'BRO_CLOc                        ', &
                                      'BRO_HO2                         ', 'BRO_NO                          ', &
                                      'BRO_NO2_M                       ', 'BRONO2_O                        ', &
                                      'BRO_O                           ', 'BRO_OH                          ', &
                                      'HBR_O                           ', 'HBR_OH                          ', &
                                      'HOBR_O                          ', 'O1D_CF3BR                       ', &
                                      'O1D_CHBR3                       ', 'O1D_H2402                       ', &
                                      'O1D_HBRa                        ', 'O1D_HBRb                        ', &
                                      'F_CH4                           ', 'F_H2                            ', &
                                      'F_H2O                           ', 'F_HNO3                          ', &
                                      'O1D_COF2                        ', 'O1D_COFCL                       ', &
                                      'CH2BR2_CL                       ', 'CH2BR2_OH                       ', &
                                      'CH3BR_CL                        ', 'CH3BR_OH                        ', &
                                      'CH3CCL3_OH                      ', 'CH3CL_CL                        ', &
                                      'CH3CL_OH                        ', 'CHBR3_CL                        ', &
                                      'CHBR3_OH                        ', 'HCFC141B_OH                     ', &
                                      'HCFC142B_OH                     ', 'HCFC22_OH                       ', &
                                      'O1D_CH2BR2                      ', 'O1D_CH3BR                       ', &
                                      'O1D_HCFC141B                    ', 'O1D_HCFC142B                    ', &
                                      'O1D_HCFC22                      ', 'CH2O_HO2                        ', &
                                      'CH2O_NO3                        ', 'CH2O_O                          ', &
                                      'CH2O_OH                         ', 'CH3O2_CH3O2a                    ', &
                                      'CH3O2_CH3O2b                    ', 'CH3O2_HO2                       ', &
                                      'CH3O2_NO                        ', 'CH3OH_OH                        ', &
                                      'CH3OOH_OH                       ', 'CH4_OH                          ', &
                                      'HCN_OH                          ', 'HCOOH_OH                        ', &
                                      'HOCH2OO_HO2                     ', 'HOCH2OO_M                       ', &
                                      'HOCH2OO_NO                      ', 'O1D_CH4a                        ', &
                                      'O1D_CH4b                        ', 'O1D_CH4c                        ', &
                                      'O1D_HCN                         ', 'usr_CO_OH                       ', &
                                      'C2H2_CL_M                       ', 'C2H2_OH_M                       ', &
                                      'C2H4_CL_M                       ', 'C2H4_O3                         ', &
                                      'C2H5O2_C2H5O2                   ', 'C2H5O2_CH3O2                    ', &
                                      'C2H5O2_HO2                      ', 'C2H5O2_NO                       ', &
                                      'C2H5OH_OH                       ', 'C2H5OOH_OH                      ', &
                                      'C2H6_CL                         ', 'C2H6_OH                         ', &
                                      'CH3CHO_NO3                      ', 'CH3CHO_OH                       ', &
                                      'CH3CN_OH                        ', 'CH3CO3_CH3CO3                   ', &
                                      'CH3CO3_CH3O2                    ', 'CH3CO3_HO2                      ', &
                                      'CH3CO3_NO                       ', 'CH3COOH_OH                      ', &
                                      'CH3COOOH_OH                     ', 'EO2_HO2                         ', &
                                      'EO2_NO                          ', 'EO_M                            ', &
                                      'EO_O2                           ', 'GLYALD_OH                       ', &
                                      'GLYOXAL_OH                      ', 'PAN_OH                          ', &
                                      'tag_C2H4_OH                     ', 'tag_CH3CO3_NO2                  ', &
                                      'usr_PAN_M                       ', 'C3H6_NO3                        ', &
                                      'C3H6_O3                         ', 'C3H7O2_CH3O2                    ', &
                                      'C3H7O2_HO2                      ', 'C3H7O2_NO                       ', &
                                      'C3H7OOH_OH                      ', 'C3H8_OH                         ', &
                                      'CH3COCHO_NO3                    ', 'CH3COCHO_OH                     ', &
                                      'HYAC_OH                         ', 'NOA_OH                          ', &
                                      'PO2_HO2                         ', 'PO2_NO                          ', &
                                      'POOH_OH                         ', 'RO2_CH3O2                       ', &
                                      'RO2_HO2                         ', 'RO2_NO                          ', &
                                      'ROOH_OH                         ', 'tag_C3H6_OH                     ', &
                                      'usr_CH3COCH3_OH                 ', 'BIGENE_NO3                      ', &
                                      'BIGENE_OH                       ', 'ENEO2_NO                        ', &
                                      'ENEO2_NOb                       ', 'HONITR_OH                       ', &
                                      'MACRO2_CH3CO3                   ', 'MACRO2_CH3O2                    ', &
                                      'MACRO2_HO2                      ', 'MACRO2_NO3                      ', &
                                      'MACRO2_NOa                      ', 'MACRO2_NOb                      ', &
                                      'MACR_O3                         ', 'MACR_OH                         ', &
                                      'MACROOH_OH                      ', 'MCO3_CH3CO3                     ', &
                                      'MCO3_CH3O2                      ', 'MCO3_HO2                        ', &
                                      'MCO3_MCO3                       ', 'MCO3_NO                         ', &
                                      'MCO3_NO3                        ', 'MEKO2_HO2                       ', &
                                      'MEKO2_NO                        ', 'MEK_OH                          ', &
                                      'MEKOOH_OH                       ', 'MPAN_OH_M                       ', &
                                      'MVK_O3                          ', 'MVK_OH                          ', &
                                      'tag_MCO3_NO2                    ', 'usr_MPAN_M                      ', &
                                      'ALKNIT_OH                       ', 'ALKO2_HO2                       ', &
                                      'ALKO2_NO                        ', 'ALKO2_NOb                       ', &
                                      'ALKOOH_OH                       ', 'BIGALK_OH                       ', &
                                      'HPALD_OH                        ', 'HYDRALD_OH                      ', &
                                      'IEPOX_OH                        ', 'ISOPAO2_CH3CO3                  ', &
                                      'ISOPAO2_CH3O2                   ', 'ISOPAO2_HO2                     ', &
                                      'ISOPAO2_NO                      ', 'ISOPAO2_NO3                     ', &
                                      'ISOPBO2_CH3CO3                  ', 'ISOPBO2_CH3O2                   ', &
                                      'ISOPBO2_HO2                     ', 'ISOPBO2_M                       ', &
                                      'ISOPBO2_NO                      ', 'ISOPBO2_NO3                     ', &
                                      'ISOPNITA_OH                     ', 'ISOPNITB_OH                     ', &
                                      'ISOP_NO3                        ', 'ISOPNO3_CH3CO3                  ', &
                                      'ISOPNO3_CH3O2                   ', 'ISOPNO3_HO2                     ', &
                                      'ISOPNO3_NO                      ', 'ISOPNO3_NO3                     ', &
                                      'ISOPNOOH_OH                     ', 'ISOP_O3                         ', &
                                      'ISOP_OH                         ', 'ISOPOOH_OH                      ', &
                                      'NC4CH2OH_OH                     ', 'NC4CHO_OH                       ', &
                                      'XO2_CH3CO3                      ', 'XO2_CH3O2                       ', &
                                      'XO2_HO2                         ', 'XO2_NO                          ', &
                                      'XO2_NO3                         ', 'XOOH_OH                         ', &
                                      'ACBZO2_HO2                      ', 'ACBZO2_NO                       ', &
                                      'BENZENE_OH                      ', 'BENZO2_HO2                      ' /)
      rxt_tag_lst(   401:   542) = (/ 'BENZO2_NO                       ', 'BENZOOH_OH                      ', &
                                      'BZALD_OH                        ', 'BZOO_HO2                        ', &
                                      'BZOOH_OH                        ', 'BZOO_NO                         ', &
                                      'C6H5O2_HO2                      ', 'C6H5O2_NO                       ', &
                                      'C6H5OOH_OH                      ', 'CRESOL_OH                       ', &
                                      'DICARBO2_HO2                    ', 'DICARBO2_NO                     ', &
                                      'DICARBO2_NO2                    ', 'MALO2_HO2                       ', &
                                      'MALO2_NO                        ', 'MALO2_NO2                       ', &
                                      'MDIALO2_HO2                     ', 'MDIALO2_NO                      ', &
                                      'MDIALO2_NO2                     ', 'PHENO2_HO2                      ', &
                                      'PHENO2_NO                       ', 'PHENOL_OH                       ', &
                                      'PHENO_NO2                       ', 'PHENO_O3                        ', &
                                      'PHENOOH_OH                      ', 'tag_ACBZO2_NO2                  ', &
                                      'TOLO2_HO2                       ', 'TOLO2_NO                        ', &
                                      'TOLOOH_OH                       ', 'TOLUENE_OH                      ', &
                                      'usr_PBZNIT_M                    ', 'XYLENES_OH                      ', &
                                      'XYLENO2_HO2                     ', 'XYLENO2_NO                      ', &
                                      'XYLENOOH_OH                     ', 'XYLOLO2_HO2                     ', &
                                      'XYLOLO2_NO                      ', 'XYLOL_OH                        ', &
                                      'XYLOLOOH_OH                     ', 'BCARY_NO3                       ', &
                                      'BCARY_O3                        ', 'BCARY_OH                        ', &
                                      'MTERP_NO3                       ', 'MTERP_O3                        ', &
                                      'MTERP_OH                        ', 'NTERPO2_CH3O2                   ', &
                                      'NTERPO2_HO2                     ', 'NTERPO2_NO                      ', &
                                      'NTERPO2_NO3                     ', 'NTERPOOH_OH                     ', &
                                      'TERP2O2_CH3O2                   ', 'TERP2O2_HO2                     ', &
                                      'TERP2O2_NO                      ', 'TERP2OOH_OH                     ', &
                                      'TERPNIT_OH                      ', 'TERPO2_CH3O2                    ', &
                                      'TERPO2_HO2                      ', 'TERPO2_NO                       ', &
                                      'TERPOOH_OH                      ', 'TERPROD1_NO3                    ', &
                                      'TERPROD1_OH                     ', 'TERPROD2_OH                     ', &
                                      'DMS_NO3                         ', 'DMS_OHa                         ', &
                                      'OCS_O                           ', 'OCS_OH                          ', &
                                      'S_O2                            ', 'SO2_OH_M                        ', &
                                      'S_O3                            ', 'SO_BRO                          ', &
                                      'SO_CLO                          ', 'S_OH                            ', &
                                      'SO_NO2                          ', 'SO_O2                           ', &
                                      'SO_O3                           ', 'SO_OCLO                         ', &
                                      'SO_OH                           ', 'usr_DMS_OH                      ', &
                                      'usr_SO3_H2O                     ', 'NH3_OH                          ', &
                                      'usr_HO2_aer                     ', 'usr_HONITR_aer                  ', &
                                      'usr_ISOPNITA_aer                ', 'usr_ISOPNITB_aer                ', &
                                      'usr_N2O5_aer                    ', 'usr_NC4CH2OH_aer                ', &
                                      'usr_NC4CHO_aer                  ', 'usr_NH4_strat_tau               ', &
                                      'usr_NO2_aer                     ', 'usr_NO3_aer                     ', &
                                      'usr_NTERPOOH_aer                ', 'usr_ONITR_aer                   ', &
                                      'usr_TERPNIT_aer                 ', 'BCARY_NO3_vbs                   ', &
                                      'BCARYO2_HO2_vbs                 ', 'BCARYO2_NO_vbs                  ', &
                                      'BCARY_O3_vbs                    ', 'BCARY_OH_vbs                    ', &
                                      'BENZENE_OH_vbs                  ', 'BENZO2_HO2_vbs                  ', &
                                      'BENZO2_NO_vbs                   ', 'ISOP_NO3_vbs                    ', &
                                      'ISOPO2_HO2_vbs                  ', 'ISOPO2_NO_vbs                   ', &
                                      'ISOP_O3_vbs                     ', 'ISOP_OH_vbs                     ', &
                                      'IVOCO2_HO2_vbs                  ', 'IVOCO2_NO_vbs                   ', &
                                      'IVOC_OH_vbs                     ', 'MTERP_NO3_vbs                   ', &
                                      'MTERPO2_HO2_vbs                 ', 'MTERPO2_NO_vbs                  ', &
                                      'MTERP_O3_vbs                    ', 'MTERP_OH_vbs                    ', &
                                      'SVOC_OH                         ', 'TOLUENE_OH_vbs                  ', &
                                      'TOLUO2_HO2_vbs                  ', 'TOLUO2_NO_vbs                   ', &
                                      'usr_GLYOXAL_aer                 ', 'XYLENES_OH_vbs                  ', &
                                      'XYLEO2_HO2_vbs                  ', 'XYLEO2_NO_vbs                   ', &
                                      'het1                            ', 'het10                           ', &
                                      'het11                           ', 'het12                           ', &
                                      'het13                           ', 'het14                           ', &
                                      'het15                           ', 'het16                           ', &
                                      'het17                           ', 'het2                            ', &
                                      'het3                            ', 'het4                            ', &
                                      'het5                            ', 'het6                            ', &
                                      'het7                            ', 'het8                            ', &
                                      'het9                            ', 'NH_50_tau                       ', &
                                      'NH_5_tau                        ', 'ST80_25_tau                     ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                                       11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                                       21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                                       31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                                       41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                                       51,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                                       61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                                       71,  72,  73,  74,  75,  76,  77,  78,  79,  80, &
                                       81,  82,  83,  84,  85,  86,  87,  88,  89,  90, &
                                       91,  92,  93,  94,  95,  96,  97,  98,  99, 100, &
                                      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, &
                                      111, 112, 113, 114, 115, 116, 117, 118, 119, 120, &
                                      121, 122, 123, 124, 125, 126, 127, 128, 129, 130, &
                                      131, 132, 133, 134, 135, 136, 137, 138, 139, 140, &
                                      141, 142, 143, 144, 145, 146, 147, 148, 149, 150, &
                                      151, 152, 153, 154, 155, 156, 157, 158, 159, 160, &
                                      161, 162, 163, 164, 165, 166, 167, 168, 169, 170, &
                                      171, 172, 173, 174, 175, 176, 177, 178, 179, 180, &
                                      181, 182, 183, 184, 185, 186, 187, 188, 189, 190, &
                                      191, 192, 193, 194, 195, 196, 197, 198, 199, 200, &
                                      201, 202, 203, 204, 205, 206, 207, 208, 209, 210, &
                                      211, 212, 213, 214, 215, 216, 217, 218, 219, 220, &
                                      221, 222, 223, 224, 225, 226, 227, 228, 229, 230, &
                                      231, 232, 233, 234, 235, 236, 237, 238, 239, 240, &
                                      241, 242, 243, 244, 245, 246, 247, 248, 249, 250, &
                                      251, 252, 253, 254, 255, 256, 257, 258, 259, 260, &
                                      261, 262, 263, 264, 265, 266, 267, 268, 269, 270, &
                                      271, 272, 273, 274, 275, 276, 277, 278, 279, 280, &
                                      281, 282, 283, 284, 285, 286, 287, 288, 289, 290, &
                                      291, 292, 293, 294, 295, 296, 297, 298, 299, 300, &
                                      301, 302, 303, 304, 305, 306, 307, 308, 309, 310, &
                                      311, 312, 313, 314, 315, 316, 317, 318, 319, 320, &
                                      321, 322, 323, 324, 325, 326, 327, 328, 329, 330, &
                                      331, 332, 333, 334, 335, 336, 337, 338, 339, 340, &
                                      341, 342, 343, 344, 345, 346, 347, 348, 349, 350, &
                                      351, 352, 353, 354, 355, 356, 357, 358, 359, 360, &
                                      361, 362, 363, 364, 365, 366, 367, 368, 369, 370, &
                                      371, 372, 373, 374, 375, 376, 377, 378, 379, 380, &
                                      381, 382, 383, 384, 385, 386, 387, 388, 389, 390, &
                                      391, 392, 393, 394, 395, 396, 397, 398, 399, 400, &
                                      401, 402, 403, 404, 405, 406, 407, 408, 409, 410, &
                                      411, 412, 413, 414, 415, 416, 417, 418, 419, 420, &
                                      421, 422, 423, 424, 425, 426, 427, 428, 429, 430, &
                                      431, 432, 433, 434, 435, 436, 437, 438, 439, 440, &
                                      441, 442, 443, 444, 445, 446, 447, 448, 449, 450, &
                                      451, 452, 453, 454, 455, 456, 457, 458, 459, 460, &
                                      461, 462, 463, 464, 465, 466, 467, 468, 469, 470, &
                                      471, 472, 473, 474, 475, 476, 477, 478, 479, 480, &
                                      481, 482, 483, 484, 485, 486, 487, 488, 489, 490, &
                                      491, 492, 493, 494, 495, 496, 497, 498, 499, 500, &
                                      501, 502, 503, 504, 505, 506, 507, 508, 509, 510, &
                                      511, 512, 513, 514, 515, 516, 517, 518, 519, 520, &
                                      521, 522, 523, 524, 525, 526, 527, 528, 529, 530, &
                                      531, 532, 533, 534, 535, 536, 537, 538, 539, 540, &
                                      541, 542 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_lst; error = ',ios
         call endrun
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_mult; error = ',ios
         call endrun
      end if
      pht_alias_lst(:,1) = (/ '                ', '                ', '                ', '                ', &
                              'userdefined     ', 'userdefined     ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', 'userdefined     ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ' /)
      pht_alias_lst(:,2) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', 'jch3ooh         ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jno2            ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              'jh2o2           ', '                ', '                ', '                ', &
                              '                ', 'jch3ooh         ', '                ', 'jmgly           ', &
                              'jch2o_a         ', 'jno2            ', '                ', 'jch3ooh         ', &
                              'jch3ooh         ', '                ', '                ', 'jacet           ', &
                              'jch3ooh         ', 'jpan            ', '                ', 'jch2o_a         ', &
                              'jch2o_a         ', 'jch3ooh         ', 'jch3cho         ', '                ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', 'jno2            ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', 'jch3cho         ', &
                              'jch3cho         ', 'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', &
                              'jch3ooh         ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jno2            ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, .10_r8, 0.2_r8, .14_r8, .20_r8, &
                          .20_r8, .006_r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 0.28_r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          .006_r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, .10_r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, .0004_r8, .0004_r8, &
                          .0004_r8, .0004_r8, .0004_r8, .0004_r8, .0004_r8, &
                          .0004_r8, .0004_r8, .0004_r8 /)
      allocate( cph_enthalpy(enthalpy_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate cph_enthalpy; error = ',ios
         call endrun
      end if
      allocate( cph_rid(enthalpy_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate cph_rid; error = ',ios
         call endrun
      end if
      cph_rid(:)      = (/             127,            131,            132,            133,            136, &
                                       139,            140,            141,            142,            145, &
                                       146,            147,            150,            152,            156, &
                                       157,            165,            166 /)
      cph_enthalpy(:) = (/   189.810000_r8,  392.190000_r8,  493.580000_r8,  101.390000_r8,  232.590000_r8, &
                             203.400000_r8,  226.580000_r8,  120.100000_r8,  194.710000_r8,  293.620000_r8, &
                              67.670000_r8,  165.300000_r8,  165.510000_r8,  313.750000_r8,  133.750000_r8, &
                             193.020000_r8,   34.470000_r8,  199.170000_r8 /)
      allocate( num_rnts(rxntot-phtcnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate num_rnts; error = ',ios
         call endrun
      end if
      num_rnts(:) = (/      1,     2,     2,     2,     2,     2,     2,     2,     3,     3, &
                            2,     2,     2,     2,     2,     3,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     3,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     3,     2,     2,     2,     2, &
                            2,     2,     2,     3,     2,     2,     3,     3,     3,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     3,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     3,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     3,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     3,     2,     2,     1,     2,     2,     2, &
                            2,     2,     2,     3,     3,     3,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     1,     2,     2,     2, &
                            2,     3,     3,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     3,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     3,     2, &
                            2,     3,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            1,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     3, &
                            2,     2,     3,     2,     2,     3,     2,     2,     2,     2, &
                            2,     2,     3,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     3,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     1,     1,     1, &
                            1,     1,     1,     1,     1,     1,     1,     1,     1,     1, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     2,     2,     2,     2,     2, &
                            2,     2,     2,     2,     2,     1,     2,     2,     2,     1, &
                            2,     1,     1,     1,     1,     2,     2,     2,     1,     1, &
                            2,     2,     2,     1,     1,     2,     1,     1,     1 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
module mo_imp_sol
  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : clscnt4, gas_pcnst, clsmap, veclen
  use cam_logfile, only : iulog
  implicit none
  private
  public :: imp_slv_inti, imp_sol
  save
  real(r8), parameter :: rel_err = 1.e-3_r8
  real(r8), parameter :: high_rel_err = 1.e-4_r8
  !-----------------------------------------------------------------------
  ! Newton-Raphson iteration limits
  !-----------------------------------------------------------------------
  integer, parameter :: itermax = 11
  integer, parameter :: cut_limit = 5
  real(r8), parameter :: sol_min = 1.e-20_r8
  real(r8), parameter :: small = 1.e-40_r8
  real(r8) :: epsilon(clscnt4)
  logical :: factor(itermax)
contains
  subroutine imp_slv_inti
    !-----------------------------------------------------------------------
    ! ... Initialize the implict solver
    !-----------------------------------------------------------------------
    use mo_chem_utls, only : get_spc_ndx
    implicit none
    !-----------------------------------------------------------------------
    ! ... Local variables
    !-----------------------------------------------------------------------
    integer :: m, ox_ndx, o3a_ndx
    real(r8) :: eps(gas_pcnst)
    factor(:) = .true.
    eps(:) = rel_err
    ox_ndx = get_spc_ndx( 'OX' )
    if( ox_ndx < 1 ) then
       ox_ndx = get_spc_ndx( 'O3' )
    end if
    if( ox_ndx > 0 ) then
       eps(ox_ndx) = high_rel_err
    end if
    m = get_spc_ndx( 'NO' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HO2NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'N2O5' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'OH' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    o3a_ndx = get_spc_ndx( 'O3A' )
    if( o3a_ndx > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XHNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XHO2NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO2NO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO2XNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    do m = 1,clscnt4
       epsilon(m) = eps(clsmap(m,4))
    end do
  end subroutine imp_slv_inti
  subroutine imp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, &
                      ncol, nlev, lchnk, prod_out, loss_out )
    !-----------------------------------------------------------------------
    ! ... imp_sol advances the volumetric mixing ratio
    ! forward one time step via the fully implicit euler scheme.
    ! this source is meant for vector architectures such as the
    ! nec sx6 and cray x1
    !-----------------------------------------------------------------------
    use chem_mods, only : rxntot, extcnt, nzcnt, permute, cls_rxt_cnt
    use mo_tracname, only : solsym
    use mo_lin_matrix, only : linmat
    use mo_nln_matrix, only : nlnmat
    use mo_lu_factor, only : lu_fac
    use mo_lu_solve, only : lu_slv
    use mo_prod_loss, only : imp_prod_loss
    use mo_indprd, only : indprd
    use time_manager, only : get_nstep
    use perf_mod, only : t_startf, t_stopf
    implicit none
    !-----------------------------------------------------------------------
    ! ... dummy args
    !-----------------------------------------------------------------------
    integer, intent(in) :: ncol ! columns in chunck
    integer, intent(in) :: nlev
    integer, intent(in) :: lchnk ! chunk id
    real(r8), intent(in) :: delt ! time step (s)
    real(r8), intent(in) :: reaction_rates(ncol*nlev,max(1,rxntot)) ! rxt rates (1/cm^3/s)
    real(r8), intent(in) :: extfrc(ncol*nlev,max(1,extcnt)) ! external in-situ forcing (1/cm^3/s)
    real(r8), intent(in) :: het_rates(ncol*nlev,max(1,gas_pcnst)) ! washout rates (1/s)
    real(r8), intent(inout) :: base_sol(ncol*nlev,gas_pcnst) ! species mixing ratios (vmr)
    real(r8), intent(out) :: prod_out(ncol*nlev,max(1,clscnt4))
    real(r8), intent(out) :: loss_out(ncol*nlev,max(1,clscnt4))
    !-----------------------------------------------------------------------
    ! ... local variables
    !-----------------------------------------------------------------------
    integer :: nr_iter
    integer :: ofl
    integer :: ofu
    integer :: avec_len
    integer :: bndx ! base index
    integer :: cndx ! class index
    integer :: pndx ! permuted class index
    integer :: i,m
    integer :: fail_cnt(veclen)
    integer :: cut_cnt(veclen)
    integer :: stp_con_cnt(veclen)
    integer :: nstep
    real(r8) :: interval_done(veclen)
    real(r8) :: dt(veclen)
    real(r8) :: dti(veclen)
    real(r8) :: max_delta(max(1,clscnt4))
    real(r8) :: ind_prd(ncol*nlev,max(1,clscnt4))
    logical :: convergence
    integer :: chnkpnts ! total spatial points in chunk; ncol*ncol
    logical :: diags_out(ncol*nlev,max(1,clscnt4))
    real(r8) :: sys_jac_blk(veclen,max(1,nzcnt))
    real(r8) :: lin_jac_blk(veclen,max(1,nzcnt))
    real(r8) :: solution_blk(veclen,max(1,clscnt4))
    real(r8) :: forcing_blk(veclen,max(1,clscnt4))
    real(r8) :: iter_invariant_blk(veclen,max(1,clscnt4))
    real(r8) :: prod_blk(veclen,max(1,clscnt4))
    real(r8) :: loss_blk(veclen,max(1,clscnt4))
    real(r8) :: ind_prd_blk(veclen,max(1,clscnt4))
    real(r8) :: sbase_sol_blk(veclen,gas_pcnst)
    real(r8) :: wrk_blk(veclen)
    logical :: spc_conv_blk(veclen,max(1,clscnt4))
    logical :: cls_conv_blk(veclen)
    logical :: time_stp_done_blk(veclen)
    real(r8) :: reaction_rates_blk(veclen,max(1,rxntot))
    real(r8) :: extfrc_blk(veclen,max(1,extcnt))
    real(r8) :: het_rates_blk(veclen,max(1,gas_pcnst))
    real(r8) :: base_sol_blk(veclen,gas_pcnst)
    chnkpnts = ncol*nlev
    prod_out = 0._r8
    loss_out = 0._r8
    diags_out = .false.
    !-----------------------------------------------------------------------
    ! ... class independent forcing
    !-----------------------------------------------------------------------
    if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
       call indprd( 4, ind_prd, clscnt4, base_sol, extfrc, &
            reaction_rates, chnkpnts )
    else
       do m = 1,clscnt4
          ind_prd(:,m) = 0._r8
       end do
    end if
    nstep = get_nstep()
    ofl = 1
    chnkpnts_loop : do
       ofu = min( chnkpnts,ofl + veclen - 1 )
       avec_len = (ofu - ofl) + 1
       reaction_rates_blk(1:avec_len,:) = reaction_rates(ofl:ofu,:)
       extfrc_blk(1:avec_len,:) = extfrc(ofl:ofu,:)
       het_rates_blk(1:avec_len,:) = het_rates(ofl:ofu,:)
       ind_prd_blk(1:avec_len,:) = ind_prd(ofl:ofu,:)
       base_sol_blk(1:avec_len,:) = base_sol(ofl:ofu,:)
       cls_conv_blk(1:avec_len) = .false.
       dt(1:avec_len) = delt
       cut_cnt(1:avec_len)  = 0
       fail_cnt(1:avec_len) = 0
       stp_con_cnt(1:avec_len) = 0
       interval_done(1:avec_len) = 0._r8
       time_stp_done_blk(1:avec_len) = .false.
       !-----------------------------------------------------------------------
       ! ... time step loop
       !-----------------------------------------------------------------------
       time_step_loop : do
          dti(1:avec_len) = 1._r8 / dt(1:avec_len)
          !-----------------------------------------------------------------------
          ! ... transfer from base to class array
          !-----------------------------------------------------------------------
          do cndx = 1,clscnt4
             bndx = clsmap(cndx,4)
             pndx = permute(cndx,4)
             do i = 1, avec_len
                solution_blk(i,pndx) = base_sol_blk(i,bndx)
             end do
          end do
          do m = 1,gas_pcnst
            sbase_sol_blk(1:avec_len,m) = base_sol_blk(1:avec_len,m)
          end do
          !-----------------------------------------------------------------------
          ! ... set the iteration invariant part of the function f(y)
          !-----------------------------------------------------------------------
          if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
             do m = 1,clscnt4
                do i = 1, avec_len
                   iter_invariant_blk(i,m) = dti(i) * solution_blk(i,m) + ind_prd_blk(i,m)
                end do
             end do
          else
             do m = 1,clscnt4
                do i = 1, avec_len
                    iter_invariant_blk(i,m) = dti(i) * solution_blk(i,m)
                end do
             end do
          end if
          !-----------------------------------------------------------------------
          ! ... the linear component
          !-----------------------------------------------------------------------
          if( cls_rxt_cnt(2,4) > 0 ) then
             call t_startf( 'lin_mat' )
             call linmat( avec_len, lin_jac_blk, base_sol_blk, &
                  reaction_rates_blk, het_rates_blk )
             call t_stopf( 'lin_mat' )
          end if
          !=======================================================================
          ! the newton-raphson iteration for f(y) = 0
          !=======================================================================
          iter_loop : do nr_iter = 1,itermax
             !-----------------------------------------------------------------------
             ! ... the non-linear component
             !-----------------------------------------------------------------------
             if( factor(nr_iter) ) then
                call t_startf( 'nln_mat' )
                call nlnmat( avec_len, sys_jac_blk, base_sol_blk, &
                     reaction_rates_blk, lin_jac_blk, dti )
                call t_stopf( 'nln_mat' )
                !-----------------------------------------------------------------------
                ! ... factor the "system" matrix
                !-----------------------------------------------------------------------
                call t_startf( 'lu_fac' )
                call lu_fac( avec_len, sys_jac_blk )
                call t_stopf( 'lu_fac' )
             end if
             !-----------------------------------------------------------------------
             ! ... form f(y)
             !-----------------------------------------------------------------------
             call t_startf( 'prod_loss' )
             call imp_prod_loss( avec_len, prod_blk, loss_blk, &
                  base_sol_blk, reaction_rates_blk, het_rates_blk )
             call t_stopf( 'prod_loss' )
             do m = 1,clscnt4
                do i = 1, avec_len
                   forcing_blk(i,m) = solution_blk(i,m)*dti(i) &
                                    - (iter_invariant_blk(i,m) + prod_blk(i,m) - loss_blk(i,m))
                end do
             end do
             !-----------------------------------------------------------------------
             ! ... solve for the mixing ratio at t(n+1)
             !-----------------------------------------------------------------------
             call t_startf( 'lu_slv' )
             call lu_slv( avec_len, sys_jac_blk, forcing_blk )
             call t_stopf( 'lu_slv' )
             do m = 1,clscnt4
                do i = 1, avec_len
                   if( .not. cls_conv_blk(i) )then
                      solution_blk(i,m) = solution_blk(i,m) + forcing_blk(i,m)
                   else
                      forcing_blk(i,m) = 0._r8
                   endif
                end do
             end do
             !-----------------------------------------------------------------------
             ! ... convergence measures and test
             !-----------------------------------------------------------------------
             conv_chk : if( nr_iter > 1 ) then
                !-----------------------------------------------------------------------
                ! ... check for convergence
                !-----------------------------------------------------------------------
                do cndx = 1,clscnt4
                   pndx = permute(cndx,4)
                   bndx = clsmap(cndx,4)
                   do i = 1, avec_len
                     if ( abs( solution_blk(i,pndx) ) > sol_min ) then
                        wrk_blk(i) = abs( forcing_blk(i,pndx)/solution_blk(i,pndx) )
                     else
                      wrk_blk(i) = 0._r8
                     endif
                   enddo
                   max_delta(cndx) = maxval( wrk_blk(1:avec_len) )
                   do i = 1, avec_len
                     solution_blk(i,pndx) = max( 0._r8,solution_blk(i,pndx) )
                     base_sol_blk(i,bndx) = solution_blk(i,pndx)
                     if ( abs( forcing_blk(i,pndx) ) > small ) then
                       spc_conv_blk(i,cndx) = abs(forcing_blk(i,pndx)) <= epsilon(cndx)*abs(solution_blk(i,pndx))
                     else
                       spc_conv_blk(i,cndx) = .true.
                     endif
                   enddo
                   where( spc_conv_blk(1:avec_len,cndx) .and. .not.diags_out(ofl:ofu,cndx) )
                      ! capture output production and loss diagnostics at converged ponits
                      prod_out(ofl:ofu,cndx) = prod_blk(1:avec_len,cndx) + ind_prd_blk(1:avec_len,cndx)
                      loss_out(ofl:ofu,cndx) = loss_blk(1:avec_len,cndx)
                      diags_out(ofl:ofu,cndx) = .true.
                   endwhere
                end do
                do i = 1, avec_len
                  if( .not. cls_conv_blk(i) ) then
                    cls_conv_blk(i) = all( spc_conv_blk(i,:) )
                  end if
                end do
                convergence = all( cls_conv_blk(:) )
                if( convergence ) then
                   exit iter_loop
                end if
             else conv_chk
!-----------------------------------------------------------------------
! ... limit iterate
!-----------------------------------------------------------------------
                do m = 1,clscnt4
                  do i = 1, avec_len
                    solution_blk(i,m) = max( 0._r8,solution_blk(i,m) )
                  end do
                end do
!-----------------------------------------------------------------------
! ... transfer latest solution back to base array
!-----------------------------------------------------------------------
                do cndx = 1,clscnt4
                   pndx = permute(cndx,4)
                   bndx = clsmap(cndx,4)
                   do i = 1, avec_len
                     base_sol_blk(i,bndx) = solution_blk(i,pndx)
                   end do
                end do
             end if conv_chk
          end do iter_loop
          !-----------------------------------------------------------------------
          ! ... check for newton-raphson convergence
          !-----------------------------------------------------------------------
          do i = 1,avec_len
            if( .not. cls_conv_blk(i) ) then
              fail_cnt(i) = fail_cnt(i) + 1
              write(iulog,'('' imp_sol: time step '',1p,g15.7,'' failed to converge @ (lchnk,vctrpos,nstep) = '',3i8)') &
                    dt(i),lchnk,ofl+i-1,nstep
              stp_con_cnt(i) = 0
              if( cut_cnt(i) < cut_limit ) then
                cut_cnt(i) = cut_cnt(i) + 1
                if( cut_cnt(i) < cut_limit ) then
                  dt(i) = .5_r8 * dt(i)
                else
                  dt(i) = .1_r8 * dt(i)
                end if
                base_sol_blk(i,:) = sbase_sol_blk(i,:)
              else
                write(iulog,'('' imp_sol: step failed to converge @ (lchnk,vctrpos,nstep,dt,time) = '',3i8,1p,2g15.7)') &
                      lchnk,ofl+i-1,nstep,dt(i),interval_done+dt(i)
                do m = 1,clscnt4
                   if( .not. spc_conv_blk(i,m) ) then
                      write(iulog,'(1x,a16,1x,1pe10.3)') solsym(clsmap(m,4)), max_delta(m)
                   end if
                end do
                cls_conv_blk(i) = .true.
                if( .not. time_stp_done_blk(i) ) then
                  interval_done(i) = interval_done(i) + dt(i)
                  time_stp_done_blk(i) = abs( delt - interval_done(i) ) <= .0001_r8
                endif
              end if
            elseif( .not. time_stp_done_blk(i) ) then
               interval_done(i) = interval_done(i) + dt(i)
               time_stp_done_blk(i) = abs( delt - interval_done(i) ) <= .0001_r8
               stp_con_cnt(i) = stp_con_cnt(i) + 1
               if( .not. time_stp_done_blk(i) ) then
                 if( stp_con_cnt(i) >= 2 ) then
                   dt(i) = 2._r8*dt(i)
                   stp_con_cnt(i) = 0
                 end if
                 dt(i) = min( dt(i),delt-interval_done(i) )
               else
                 base_sol(ofl+i-1,1:gas_pcnst) = base_sol_blk(i,1:gas_pcnst)
               endif
            endif
          end do
          convergence = all( cls_conv_blk(:) )
          do i = 1,avec_len
            if( cls_conv_blk(i) .and. .not. time_stp_done_blk(i) ) then
              cls_conv_blk(i) = .false.
            endif
          end do
          if( .not. convergence ) then
            cycle time_step_loop
          endif
          !-----------------------------------------------------------------------
          ! ... check for time step done
          !-----------------------------------------------------------------------
          if( all( time_stp_done_blk(1:avec_len) ) ) then
             exit time_step_loop
          end if
       end do time_step_loop

       ofl = ofu + 1
       if( ofl > chnkpnts ) then
          exit chnkpnts_loop
       end if
    end do chnkpnts_loop
  end subroutine imp_sol
end module mo_imp_sol

module mo_exp_sol

  private
  public :: exp_sol
  public :: exp_sol_inti

contains

  subroutine exp_sol_inti

    use mo_tracname, only : solsym
    use chem_mods,   only : clscnt1, clsmap
    use cam_history, only : addfld

    implicit none

    integer :: i,j

    do i = 1,clscnt1

       j = clsmap(i,1)
       call addfld( trim(solsym(j))//'_CHMP', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate' )
       call addfld( trim(solsym(j))//'_CHML', (/ 'lev' /), 'I', '/cm3/s', 'chemical loss rate' )

    enddo
  end subroutine exp_sol_inti


  subroutine exp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, xhnm, ncol, lchnk, ltrop )
    !-----------------------------------------------------------------------
    !      	... Exp_sol advances the volumetric mixing ratio
    !           forward one time step via the fully explicit
    !           Euler scheme
    !-----------------------------------------------------------------------

    use chem_mods,     only : clscnt1, extcnt, gas_pcnst, clsmap, rxntot
    use ppgrid,        only : pcols, pver
    use mo_prod_loss,  only : exp_prod_loss
    use mo_indprd,     only : indprd
    use shr_kind_mod,  only : r8 => shr_kind_r8
    use cam_history,   only : outfld
    use mo_tracname,   only : solsym

    implicit none
    !-----------------------------------------------------------------------
    !     	... Dummy arguments
    !-----------------------------------------------------------------------
    integer,  intent(in)    ::  ncol                                ! columns in chunck
    integer,  intent(in)    ::  lchnk                               ! chunk id
    real(r8), intent(in)    ::  delt                                ! time step (s)
    real(r8), intent(in)    ::  het_rates(ncol,pver,max(1,gas_pcnst))  ! het rates (1/cm^3/s)
    real(r8), intent(in)    ::  reaction_rates(ncol,pver,rxntot)    ! rxt rates (1/cm^3/s)
    real(r8), intent(in)    ::  extfrc(ncol,pver,extcnt)            ! "external insitu forcing" (1/cm^3/s)
    real(r8), intent(in)    ::  xhnm(ncol,pver)
    integer,  intent(in)    ::  ltrop(pcols)                        ! chemistry troposphere boundary (index)
    real(r8), intent(inout) ::  base_sol(ncol,pver,gas_pcnst)       ! working mixing ratios (vmr)

    !-----------------------------------------------------------------------
    !     	... Local variables
    !-----------------------------------------------------------------------
    integer  ::  i, k, l, m
    integer  ::  chnkpnts
    real(r8), dimension(ncol,pver,max(1,clscnt1)) :: &
         prod, &
         loss
    real(r8), dimension(ncol,pver,clscnt1) :: ind_prd

    real(r8), dimension(ncol,pver) :: wrk

    chnkpnts = ncol*pver
    !-----------------------------------------------------------------------      
    !        ... Put "independent" production in the forcing
    !-----------------------------------------------------------------------      
    call indprd( 1, ind_prd, clscnt1, base_sol, extfrc, &
                 reaction_rates, chnkpnts )

    !-----------------------------------------------------------------------      
    !      	... Form F(y)
    !-----------------------------------------------------------------------      
    call exp_prod_loss( 1, chnkpnts, prod, loss, base_sol, reaction_rates, &
                        het_rates, chnkpnts )

    !-----------------------------------------------------------------------      
    !    	... Solve for the mixing ratio at t(n+1)
    !-----------------------------------------------------------------------      
    do m = 1,clscnt1
       l = clsmap(m,1)
       do i = 1,ncol
          do k = ltrop(i)+1,pver
             base_sol(i,k,l)  = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
          end do
       end do

       wrk(:,:) = (prod(:,:,m) + ind_prd(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHMP', wrk(:,:), ncol, lchnk )
       wrk(:,:) = (loss(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHML', wrk(:,:), ncol, lchnk )
       
    end do

  end subroutine exp_sol

end module mo_exp_sol


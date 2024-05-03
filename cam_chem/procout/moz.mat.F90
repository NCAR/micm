      module mo_prod_loss

      use shr_kind_mod, only : r8 => shr_kind_r8

      use chem_mods, only : veclen
      private
      public :: exp_prod_loss
      public :: imp_prod_loss

      contains

      subroutine exp_prod_loss( ofl, ofu, prod, loss, y, &
                                rxt, het_rates, chnkpnts )

      use chem_mods,    only : gas_pcnst,rxntot,clscnt1

      implicit none

!--------------------------------------------------------------------
!     ... dummy args                                                                      
!--------------------------------------------------------------------
      integer, intent(in) :: ofl, ofu, chnkpnts
      real(r8), dimension(chnkpnts,max(1,clscnt1)), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in)    ::  y(chnkpnts,gas_pcnst)
      real(r8), intent(in)    ::  rxt(chnkpnts,rxntot)
      real(r8), intent(in)    ::  het_rates(chnkpnts,gas_pcnst)


!--------------------------------------------------------------------
!     ... local variables                                                                 
!--------------------------------------------------------------------
      integer :: k

!--------------------------------------------------------------------
!       ... loss and production for Explicit method
!--------------------------------------------------------------------

      do k = ofl,ofu
         loss(k,1) = ( + het_rates(k,186))* y(k,186)
         prod(k,1) = 0._r8
         loss(k,2) = ( + het_rates(k,187))* y(k,187)
         prod(k,2) = 0._r8
      end do

      end subroutine exp_prod_loss

      subroutine imp_prod_loss( avec_len, prod, loss, y, &
                                rxt, het_rates )

      use chem_mods,    only : gas_pcnst,rxntot,clscnt4

      implicit none

!--------------------------------------------------------------------
!     ... dummy args                                                                      
!--------------------------------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), dimension(veclen,clscnt4), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in)       ::  y(veclen,gas_pcnst)
      real(r8), intent(in)       ::  rxt(veclen,rxntot)
      real(r8), intent(in)       ::  het_rates(veclen,gas_pcnst)


!--------------------------------------------------------------------
!     ... local variables                                                                 
!--------------------------------------------------------------------
      integer :: k

!--------------------------------------------------------------------
!       ... loss and production for Implicit method
!--------------------------------------------------------------------

      do k = 1,avec_len
         loss(k,154) = (rxt(k,357)* y(k,217) + rxt(k,19) + het_rates(k,1))* y(k,1)
         prod(k,154) =rxt(k,360)*y(k,189)*y(k,124)
         loss(k,151) = (rxt(k,361)* y(k,217) + rxt(k,20) + het_rates(k,2))* y(k,2)
         prod(k,151) =rxt(k,358)*y(k,203)*y(k,189)
         loss(k,1) = ( + het_rates(k,3))* y(k,3)
         prod(k,1) = 0._r8
         loss(k,2) = ( + het_rates(k,4))* y(k,4)
         prod(k,2) = 0._r8
         loss(k,3) = ( + het_rates(k,5))* y(k,5)
         prod(k,3) = 0._r8
         loss(k,185) = (rxt(k,440)* y(k,126) +rxt(k,441)* y(k,134) +rxt(k,442) &
                 * y(k,217) + het_rates(k,6))* y(k,6)
         prod(k,185) = 0._r8
         loss(k,71) = (rxt(k,399)* y(k,217) + het_rates(k,7))* y(k,7)
         prod(k,71) = 0._r8
         loss(k,121) = (rxt(k,402)* y(k,217) + rxt(k,21) + het_rates(k,8))* y(k,8)
         prod(k,121) =rxt(k,400)*y(k,203)*y(k,191)
         loss(k,72) = ( + rxt(k,22) + het_rates(k,9))* y(k,9)
         prod(k,72) =.120_r8*rxt(k,399)*y(k,217)*y(k,7)
         loss(k,114) = ( + rxt(k,23) + het_rates(k,10))* y(k,10)
         prod(k,114) = (.100_r8*rxt(k,441)*y(k,6) +.100_r8*rxt(k,444)*y(k,110)) &
                 *y(k,134)
         loss(k,127) = ( + rxt(k,24) + het_rates(k,11))* y(k,11)
         prod(k,127) = (.500_r8*rxt(k,401)*y(k,191) +.200_r8*rxt(k,428)*y(k,223) + &
                 .060_r8*rxt(k,434)*y(k,226))*y(k,124) +.500_r8*rxt(k,21)*y(k,8) &
                  +rxt(k,22)*y(k,9) +.200_r8*rxt(k,70)*y(k,179) +.060_r8*rxt(k,72) &
                 *y(k,183)
         loss(k,96) = ( + rxt(k,25) + het_rates(k,12))* y(k,12)
         prod(k,96) = (.200_r8*rxt(k,428)*y(k,223) +.200_r8*rxt(k,434)*y(k,226)) &
                 *y(k,124) +.200_r8*rxt(k,70)*y(k,179) +.200_r8*rxt(k,72)*y(k,183)
         loss(k,145) = ( + rxt(k,26) + het_rates(k,13))* y(k,13)
         prod(k,145) = (.200_r8*rxt(k,428)*y(k,223) +.150_r8*rxt(k,434)*y(k,226)) &
                 *y(k,124) +rxt(k,46)*y(k,94) +rxt(k,56)*y(k,116) +.200_r8*rxt(k,70) &
                 *y(k,179) +.150_r8*rxt(k,72)*y(k,183)
         loss(k,105) = ( + rxt(k,27) + het_rates(k,14))* y(k,14)
         prod(k,105) =.210_r8*rxt(k,434)*y(k,226)*y(k,124) +.210_r8*rxt(k,72)*y(k,183)
         loss(k,86) = (rxt(k,362)* y(k,217) + het_rates(k,15))* y(k,15)
         prod(k,86) = (.050_r8*rxt(k,441)*y(k,6) +.050_r8*rxt(k,444)*y(k,110)) &
                 *y(k,134)
         loss(k,110) = (rxt(k,328)* y(k,126) +rxt(k,329)* y(k,217) + het_rates(k,16)) &
                 * y(k,16)
         prod(k,110) = 0._r8
         loss(k,209) = (rxt(k,212)* y(k,42) +rxt(k,214)* y(k,134) +rxt(k,213) &
                 * y(k,203) + het_rates(k,17))* y(k,17)
         prod(k,209) = (rxt(k,75) +2.000_r8*rxt(k,215)*y(k,19) +rxt(k,216)*y(k,59) + &
                 rxt(k,217)*y(k,59) +rxt(k,220)*y(k,124) +rxt(k,223)*y(k,133) + &
                 rxt(k,224)*y(k,217) +rxt(k,470)*y(k,150))*y(k,19) &
                  + (rxt(k,202)*y(k,34) +rxt(k,228)*y(k,35) + &
                 3.000_r8*rxt(k,229)*y(k,55) +2.000_r8*rxt(k,230)*y(k,78) + &
                 rxt(k,231)*y(k,81) +2.000_r8*rxt(k,251)*y(k,41) +rxt(k,252)*y(k,43)) &
                 *y(k,216) + (rxt(k,226)*y(k,81) +2.000_r8*rxt(k,240)*y(k,41) + &
                 rxt(k,242)*y(k,43) +3.000_r8*rxt(k,247)*y(k,55))*y(k,217) &
                  + (2.000_r8*rxt(k,239)*y(k,41) +rxt(k,241)*y(k,43) + &
                 3.000_r8*rxt(k,246)*y(k,55))*y(k,56) + (rxt(k,99) + &
                 rxt(k,225)*y(k,133))*y(k,81) +rxt(k,74)*y(k,18) +rxt(k,77)*y(k,20) &
                  +rxt(k,79)*y(k,34) +rxt(k,80)*y(k,35) +2.000_r8*rxt(k,86)*y(k,41) &
                  +rxt(k,87)*y(k,43) +3.000_r8*rxt(k,90)*y(k,55) +2.000_r8*rxt(k,98) &
                 *y(k,78) +rxt(k,105)*y(k,91)
         loss(k,87) = ( + rxt(k,74) + het_rates(k,18))* y(k,18)
         prod(k,87) = (rxt(k,531)*y(k,91) +rxt(k,536)*y(k,91))*y(k,85) &
                  +rxt(k,218)*y(k,59)*y(k,19)
         loss(k,225) = (2._r8*rxt(k,215)* y(k,19) + (rxt(k,216) +rxt(k,217) + &
                 rxt(k,218))* y(k,59) +rxt(k,220)* y(k,124) +rxt(k,221)* y(k,125) &
                  +rxt(k,223)* y(k,133) +rxt(k,470)* y(k,150) +rxt(k,219)* y(k,203) &
                  +rxt(k,224)* y(k,217) + rxt(k,75) + het_rates(k,19))* y(k,19)
         prod(k,225) = (rxt(k,76) +rxt(k,222)*y(k,133))*y(k,20) +rxt(k,214)*y(k,134) &
                 *y(k,17) +rxt(k,232)*y(k,216)*y(k,81) +rxt(k,227)*y(k,133)*y(k,91)
         loss(k,139) = (rxt(k,222)* y(k,133) + rxt(k,76) + rxt(k,77) + rxt(k,525) &
                  + rxt(k,528) + rxt(k,533) + het_rates(k,20))* y(k,20)
         prod(k,139) =rxt(k,221)*y(k,125)*y(k,19)
         loss(k,4) = ( + het_rates(k,21))* y(k,21)
         prod(k,4) = 0._r8
         loss(k,90) = (rxt(k,403)* y(k,217) + het_rates(k,22))* y(k,22)
         prod(k,90) =rxt(k,28)*y(k,23) +rxt(k,406)*y(k,193)*y(k,124)
         loss(k,108) = (rxt(k,405)* y(k,217) + rxt(k,28) + het_rates(k,23))* y(k,23)
         prod(k,108) =rxt(k,404)*y(k,203)*y(k,193)
         loss(k,98) = (rxt(k,277)* y(k,56) +rxt(k,278)* y(k,217) + het_rates(k,24)) &
                 * y(k,24)
         prod(k,98) = 0._r8
         loss(k,140) = (rxt(k,279)* y(k,56) +rxt(k,280)* y(k,134) +rxt(k,305) &
                 * y(k,217) + het_rates(k,25))* y(k,25)
         prod(k,140) = 0._r8
         loss(k,92) = (rxt(k,285)* y(k,217) + het_rates(k,26))* y(k,26)
         prod(k,92) = (.400_r8*rxt(k,281)*y(k,194) +.200_r8*rxt(k,282)*y(k,198)) &
                 *y(k,194)
         loss(k,109) = (rxt(k,286)* y(k,217) + rxt(k,29) + het_rates(k,27))* y(k,27)
         prod(k,109) =rxt(k,283)*y(k,203)*y(k,194)
         loss(k,99) = (rxt(k,287)* y(k,56) +rxt(k,288)* y(k,217) + het_rates(k,28)) &
                 * y(k,28)
         prod(k,99) = 0._r8
         loss(k,186) = (rxt(k,308)* y(k,126) +rxt(k,309)* y(k,134) +rxt(k,326) &
                 * y(k,217) + het_rates(k,29))* y(k,29)
         prod(k,186) =.130_r8*rxt(k,386)*y(k,134)*y(k,98) +.700_r8*rxt(k,55)*y(k,111)
         loss(k,120) = (rxt(k,313)* y(k,217) + rxt(k,30) + het_rates(k,30))* y(k,30)
         prod(k,120) =rxt(k,311)*y(k,203)*y(k,195)
         loss(k,57) = (rxt(k,314)* y(k,217) + het_rates(k,31))* y(k,31)
         prod(k,57) = 0._r8
         loss(k,93) = (rxt(k,409)* y(k,217) + rxt(k,31) + het_rates(k,32))* y(k,32)
         prod(k,93) =rxt(k,407)*y(k,203)*y(k,196)
         loss(k,54) = (rxt(k,201)* y(k,216) + rxt(k,78) + het_rates(k,33))* y(k,33)
         prod(k,54) = 0._r8
         loss(k,66) = (rxt(k,202)* y(k,216) + rxt(k,79) + het_rates(k,34))* y(k,34)
         prod(k,66) = 0._r8
         loss(k,67) = (rxt(k,228)* y(k,216) + rxt(k,80) + het_rates(k,35))* y(k,35)
         prod(k,67) = 0._r8
         loss(k,58) = (rxt(k,203)* y(k,216) + rxt(k,81) + het_rates(k,36))* y(k,36)
         prod(k,58) = 0._r8
         loss(k,68) = (rxt(k,204)* y(k,216) + rxt(k,82) + het_rates(k,37))* y(k,37)
         prod(k,68) = 0._r8
         loss(k,59) = (rxt(k,205)* y(k,216) + rxt(k,83) + het_rates(k,38))* y(k,38)
         prod(k,59) = 0._r8
         loss(k,69) = (rxt(k,206)* y(k,216) + rxt(k,84) + het_rates(k,39))* y(k,39)
         prod(k,69) = 0._r8
         loss(k,60) = (rxt(k,207)* y(k,216) + rxt(k,85) + het_rates(k,40))* y(k,40)
         prod(k,60) = 0._r8
         loss(k,129) = (rxt(k,239)* y(k,56) +rxt(k,251)* y(k,216) +rxt(k,240) &
                 * y(k,217) + rxt(k,86) + het_rates(k,41))* y(k,41)
         prod(k,129) = 0._r8
         loss(k,213) = (rxt(k,212)* y(k,17) +rxt(k,176)* y(k,56) +rxt(k,257)* y(k,126) &
                  +rxt(k,258)* y(k,133) +rxt(k,256)* y(k,203) +rxt(k,259)* y(k,217) &
                  + rxt(k,32) + rxt(k,33) + het_rates(k,42))* y(k,42)
         prod(k,213) = (rxt(k,183)*y(k,59) +2.000_r8*rxt(k,260)*y(k,198) + &
                 rxt(k,261)*y(k,198) +rxt(k,263)*y(k,124) + &
                 .700_r8*rxt(k,282)*y(k,194) +rxt(k,293)*y(k,197) + &
                 rxt(k,310)*y(k,195) +.800_r8*rxt(k,322)*y(k,220) + &
                 .880_r8*rxt(k,334)*y(k,209) +2.000_r8*rxt(k,343)*y(k,211) + &
                 1.500_r8*rxt(k,367)*y(k,205) +.750_r8*rxt(k,372)*y(k,206) + &
                 .800_r8*rxt(k,381)*y(k,101) +.800_r8*rxt(k,392)*y(k,225) + &
                 .750_r8*rxt(k,446)*y(k,215) +.930_r8*rxt(k,451)*y(k,221) + &
                 .950_r8*rxt(k,456)*y(k,222))*y(k,198) &
                  + (.500_r8*rxt(k,299)*y(k,202) +rxt(k,320)*y(k,219) + &
                 rxt(k,324)*y(k,220) +.500_r8*rxt(k,330)*y(k,200) + &
                 .250_r8*rxt(k,337)*y(k,209) +rxt(k,346)*y(k,211) + &
                 .100_r8*rxt(k,359)*y(k,189) +.920_r8*rxt(k,369)*y(k,205) + &
                 .250_r8*rxt(k,394)*y(k,225) +.340_r8*rxt(k,453)*y(k,221) + &
                 .320_r8*rxt(k,458)*y(k,222))*y(k,124) + (rxt(k,264)*y(k,52) + &
                 .300_r8*rxt(k,265)*y(k,53) +.500_r8*rxt(k,297)*y(k,51) + &
                 .800_r8*rxt(k,302)*y(k,74) +rxt(k,304)*y(k,139) + &
                 .500_r8*rxt(k,352)*y(k,109) +.400_r8*rxt(k,357)*y(k,1) + &
                 .300_r8*rxt(k,377)*y(k,99) +.680_r8*rxt(k,462)*y(k,178))*y(k,217) &
                  + (rxt(k,280)*y(k,25) +.500_r8*rxt(k,309)*y(k,29) + &
                 .120_r8*rxt(k,339)*y(k,105) +.600_r8*rxt(k,353)*y(k,111) + &
                 .910_r8*rxt(k,386)*y(k,98) +.340_r8*rxt(k,441)*y(k,6) + &
                 .340_r8*rxt(k,444)*y(k,110))*y(k,134) + (.500_r8*rxt(k,328)*y(k,16) + &
                 .250_r8*rxt(k,336)*y(k,209) +rxt(k,347)*y(k,211) + &
                 rxt(k,370)*y(k,205))*y(k,126) + (.250_r8*rxt(k,333)*y(k,209) + &
                 rxt(k,342)*y(k,211) +rxt(k,366)*y(k,205) + &
                 .250_r8*rxt(k,391)*y(k,225))*y(k,197) + (.180_r8*rxt(k,39) + &
                 rxt(k,273)*y(k,216) +rxt(k,274)*y(k,216))*y(k,54) &
                  + (.150_r8*rxt(k,323)*y(k,220) +.450_r8*rxt(k,344)*y(k,211)) &
                 *y(k,203) +.100_r8*rxt(k,19)*y(k,1) +.100_r8*rxt(k,20)*y(k,2) &
                  +rxt(k,38)*y(k,53) +rxt(k,43)*y(k,74) +.330_r8*rxt(k,45)*y(k,93) &
                  +rxt(k,47)*y(k,95) +rxt(k,49)*y(k,103) +1.340_r8*rxt(k,50)*y(k,105) &
                  +rxt(k,57)*y(k,127) +rxt(k,62)*y(k,146) +rxt(k,63)*y(k,147) &
                  +.375_r8*rxt(k,65)*y(k,174) +.400_r8*rxt(k,67)*y(k,176) &
                  +.680_r8*rxt(k,69)*y(k,178) +2.000_r8*rxt(k,300)*y(k,201) &
                  +rxt(k,270)*y(k,204) +2.000_r8*rxt(k,345)*y(k,211)*y(k,211)
         loss(k,146) = (rxt(k,241)* y(k,56) +rxt(k,252)* y(k,216) +rxt(k,242) &
                 * y(k,217) + rxt(k,87) + het_rates(k,43))* y(k,43)
         prod(k,146) = 0._r8
         loss(k,61) = (rxt(k,243)* y(k,217) + rxt(k,88) + het_rates(k,44))* y(k,44)
         prod(k,61) = 0._r8
         loss(k,190) = (rxt(k,289)* y(k,126) +rxt(k,290)* y(k,217) + rxt(k,34) &
                  + het_rates(k,45))* y(k,45)
         prod(k,190) = (rxt(k,284)*y(k,194) +.270_r8*rxt(k,312)*y(k,195) + &
                 rxt(k,320)*y(k,219) +rxt(k,330)*y(k,200) +rxt(k,349)*y(k,213) + &
                 .400_r8*rxt(k,359)*y(k,189))*y(k,124) + (rxt(k,285)*y(k,26) + &
                 .500_r8*rxt(k,286)*y(k,27) +.800_r8*rxt(k,357)*y(k,1))*y(k,217) &
                  + (.500_r8*rxt(k,309)*y(k,29) +.100_r8*rxt(k,353)*y(k,111))*y(k,134) &
                  + (1.600_r8*rxt(k,281)*y(k,194) +.800_r8*rxt(k,282)*y(k,198)) &
                 *y(k,194) +.400_r8*rxt(k,19)*y(k,1) +.400_r8*rxt(k,20)*y(k,2) &
                  +rxt(k,328)*y(k,126)*y(k,16) +rxt(k,29)*y(k,27) +.330_r8*rxt(k,45) &
                 *y(k,93) +rxt(k,53)*y(k,108) +rxt(k,62)*y(k,146) &
                  +.200_r8*rxt(k,348)*y(k,213)*y(k,203)
         loss(k,112) = (rxt(k,244)* y(k,56) +rxt(k,245)* y(k,217) + rxt(k,89) &
                  + het_rates(k,46))* y(k,46)
         prod(k,112) = 0._r8
         loss(k,55) = (rxt(k,291)* y(k,217) + het_rates(k,47))* y(k,47)
         prod(k,55) = 0._r8
         loss(k,180) = (rxt(k,327)* y(k,217) + rxt(k,35) + het_rates(k,48))* y(k,48)
         prod(k,180) = (.820_r8*rxt(k,312)*y(k,195) +.500_r8*rxt(k,330)*y(k,200) + &
                 .250_r8*rxt(k,359)*y(k,189) +.270_r8*rxt(k,453)*y(k,221) + &
                 .040_r8*rxt(k,458)*y(k,222))*y(k,124) &
                  + (.820_r8*rxt(k,310)*y(k,195) +.150_r8*rxt(k,451)*y(k,221) + &
                 .025_r8*rxt(k,456)*y(k,222))*y(k,198) + (.250_r8*rxt(k,19) + &
                 .800_r8*rxt(k,357)*y(k,217))*y(k,1) + (.520_r8*rxt(k,441)*y(k,6) + &
                 .520_r8*rxt(k,444)*y(k,110))*y(k,134) + (.500_r8*rxt(k,69) + &
                 .500_r8*rxt(k,462)*y(k,217))*y(k,178) +.250_r8*rxt(k,20)*y(k,2) &
                  +.500_r8*rxt(k,328)*y(k,126)*y(k,16) +.820_r8*rxt(k,30)*y(k,30) &
                  +.170_r8*rxt(k,45)*y(k,93) +.300_r8*rxt(k,65)*y(k,174) &
                  +.050_r8*rxt(k,67)*y(k,176)
         loss(k,200) = (rxt(k,315)* y(k,126) +rxt(k,316)* y(k,217) + rxt(k,36) &
                  + het_rates(k,49))* y(k,49)
         prod(k,200) = (.250_r8*rxt(k,337)*y(k,209) +.050_r8*rxt(k,375)*y(k,206) + &
                 .250_r8*rxt(k,394)*y(k,225) +.170_r8*rxt(k,412)*y(k,199) + &
                 .170_r8*rxt(k,418)*y(k,212) +.400_r8*rxt(k,428)*y(k,223) + &
                 .540_r8*rxt(k,434)*y(k,226) +.510_r8*rxt(k,437)*y(k,228))*y(k,124) &
                  + (.250_r8*rxt(k,336)*y(k,209) +.050_r8*rxt(k,376)*y(k,206) + &
                 .250_r8*rxt(k,395)*y(k,225))*y(k,126) &
                  + (.500_r8*rxt(k,322)*y(k,220) +.240_r8*rxt(k,334)*y(k,209) + &
                 .100_r8*rxt(k,392)*y(k,225))*y(k,198) &
                  + (.880_r8*rxt(k,339)*y(k,105) +.500_r8*rxt(k,353)*y(k,111)) &
                 *y(k,134) + (.250_r8*rxt(k,333)*y(k,209) + &
                 .250_r8*rxt(k,391)*y(k,225))*y(k,197) &
                  + (.070_r8*rxt(k,411)*y(k,199) +.070_r8*rxt(k,417)*y(k,212)) &
                 *y(k,203) + (rxt(k,317)*y(k,95) +rxt(k,318)*y(k,127))*y(k,217) &
                  +.180_r8*rxt(k,23)*y(k,10) +rxt(k,27)*y(k,14) +.400_r8*rxt(k,70) &
                 *y(k,179) +.540_r8*rxt(k,72)*y(k,183) +.510_r8*rxt(k,73)*y(k,185)
         loss(k,156) = (rxt(k,296)* y(k,217) + het_rates(k,50))* y(k,50)
         prod(k,156) = (.100_r8*rxt(k,293)*y(k,198) +.150_r8*rxt(k,294)*y(k,203)) &
                 *y(k,197) +.120_r8*rxt(k,309)*y(k,134)*y(k,29) &
                  +.150_r8*rxt(k,344)*y(k,211)*y(k,203)
         loss(k,148) = (rxt(k,297)* y(k,217) + rxt(k,37) + het_rates(k,51))* y(k,51)
         prod(k,148) = (.400_r8*rxt(k,294)*y(k,197) +.400_r8*rxt(k,344)*y(k,211)) &
                 *y(k,203)
         loss(k,166) = (rxt(k,264)* y(k,217) + het_rates(k,52))* y(k,52)
         prod(k,166) = (rxt(k,261)*y(k,198) +.300_r8*rxt(k,282)*y(k,194) + &
                 .500_r8*rxt(k,322)*y(k,220) +.250_r8*rxt(k,334)*y(k,209) + &
                 .250_r8*rxt(k,367)*y(k,205) +.250_r8*rxt(k,372)*y(k,206) + &
                 .200_r8*rxt(k,381)*y(k,101) +.300_r8*rxt(k,392)*y(k,225) + &
                 .250_r8*rxt(k,446)*y(k,215) +.250_r8*rxt(k,451)*y(k,221) + &
                 .250_r8*rxt(k,456)*y(k,222))*y(k,198)
         loss(k,115) = (rxt(k,265)* y(k,217) + rxt(k,38) + het_rates(k,53))* y(k,53)
         prod(k,115) =rxt(k,262)*y(k,203)*y(k,198)
         loss(k,210) = (rxt(k,177)* y(k,56) +rxt(k,233)* y(k,73) + (rxt(k,272) + &
                 rxt(k,273) +rxt(k,274))* y(k,216) +rxt(k,266)* y(k,217) + rxt(k,39) &
                  + rxt(k,40) + het_rates(k,54))* y(k,54)
         prod(k,210) =.100_r8*rxt(k,309)*y(k,134)*y(k,29)
         loss(k,125) = (rxt(k,246)* y(k,56) +rxt(k,229)* y(k,216) +rxt(k,247) &
                 * y(k,217) + rxt(k,90) + het_rates(k,55))* y(k,55)
         prod(k,125) = 0._r8
         loss(k,220) = (rxt(k,287)* y(k,28) +rxt(k,239)* y(k,41) +rxt(k,176)* y(k,42) &
                  +rxt(k,241)* y(k,43) +rxt(k,244)* y(k,46) +rxt(k,177)* y(k,54) &
                  +rxt(k,246)* y(k,55) +rxt(k,189)* y(k,60) +rxt(k,178)* y(k,77) &
                  +rxt(k,179)* y(k,79) +rxt(k,198)* y(k,92) +rxt(k,182)* y(k,134) &
                  + (rxt(k,180) +rxt(k,181))* y(k,203) + het_rates(k,56))* y(k,56)
         prod(k,220) = (4.000_r8*rxt(k,201)*y(k,33) +rxt(k,202)*y(k,34) + &
                 2.000_r8*rxt(k,203)*y(k,36) +2.000_r8*rxt(k,204)*y(k,37) + &
                 2.000_r8*rxt(k,205)*y(k,38) +rxt(k,206)*y(k,39) + &
                 2.000_r8*rxt(k,207)*y(k,40) +rxt(k,208)*y(k,85) +rxt(k,238)*y(k,65) + &
                 rxt(k,253)*y(k,82) +rxt(k,254)*y(k,83) +rxt(k,255)*y(k,84))*y(k,216) &
                  + (rxt(k,93) +rxt(k,183)*y(k,198) +2.000_r8*rxt(k,184)*y(k,59) + &
                 rxt(k,186)*y(k,59) +rxt(k,188)*y(k,124) +rxt(k,193)*y(k,133) + &
                 rxt(k,194)*y(k,217) +rxt(k,217)*y(k,19) +rxt(k,471)*y(k,150))*y(k,59) &
                  + (rxt(k,197)*y(k,85) +3.000_r8*rxt(k,243)*y(k,44) + &
                 rxt(k,245)*y(k,46) +rxt(k,248)*y(k,82) +rxt(k,249)*y(k,83) + &
                 rxt(k,250)*y(k,84))*y(k,217) + (rxt(k,103) +rxt(k,196)*y(k,133)) &
                 *y(k,85) +rxt(k,74)*y(k,18) +4.000_r8*rxt(k,78)*y(k,33) +rxt(k,79) &
                 *y(k,34) +2.000_r8*rxt(k,81)*y(k,36) +2.000_r8*rxt(k,82)*y(k,37) &
                  +2.000_r8*rxt(k,83)*y(k,38) +rxt(k,84)*y(k,39) +2.000_r8*rxt(k,85) &
                 *y(k,40) +3.000_r8*rxt(k,88)*y(k,44) +rxt(k,89)*y(k,46) &
                  +2.000_r8*rxt(k,91)*y(k,57) +2.000_r8*rxt(k,92)*y(k,58) +rxt(k,94) &
                 *y(k,60) +rxt(k,97)*y(k,65) +rxt(k,100)*y(k,82) +rxt(k,101)*y(k,83) &
                  +rxt(k,102)*y(k,84) +rxt(k,106)*y(k,92)
         loss(k,70) = ( + rxt(k,91) + het_rates(k,57))* y(k,57)
         prod(k,70) = (rxt(k,524)*y(k,92) +rxt(k,529)*y(k,60) +rxt(k,530)*y(k,92) + &
                 rxt(k,534)*y(k,60) +rxt(k,535)*y(k,92) +rxt(k,539)*y(k,60))*y(k,85) &
                  +rxt(k,189)*y(k,60)*y(k,56) +rxt(k,185)*y(k,59)*y(k,59)
         loss(k,52) = ( + rxt(k,92) + rxt(k,211) + het_rates(k,58))* y(k,58)
         prod(k,52) =rxt(k,210)*y(k,59)*y(k,59)
         loss(k,219) = ((rxt(k,216) +rxt(k,217) +rxt(k,218))* y(k,19) &
                  + 2._r8*(rxt(k,184) +rxt(k,185) +rxt(k,186) +rxt(k,210))* y(k,59) &
                  +rxt(k,188)* y(k,124) +rxt(k,190)* y(k,125) +rxt(k,193)* y(k,133) &
                  +rxt(k,471)* y(k,150) +rxt(k,183)* y(k,198) +rxt(k,187)* y(k,203) &
                  + (rxt(k,194) +rxt(k,195))* y(k,217) + rxt(k,93) + het_rates(k,59)) &
                 * y(k,59)
         prod(k,219) = (rxt(k,181)*y(k,203) +rxt(k,182)*y(k,134) +rxt(k,198)*y(k,92)) &
                 *y(k,56) + (rxt(k,95) +rxt(k,191)*y(k,133))*y(k,60) &
                  + (rxt(k,199)*y(k,133) +rxt(k,200)*y(k,217))*y(k,92) + (rxt(k,107) + &
                 rxt(k,476)*y(k,150))*y(k,136) +2.000_r8*rxt(k,211)*y(k,58) &
                  +rxt(k,209)*y(k,216)*y(k,85)
         loss(k,177) = (rxt(k,189)* y(k,56) + (rxt(k,529) +rxt(k,534) +rxt(k,539)) &
                 * y(k,85) +rxt(k,191)* y(k,133) +rxt(k,192)* y(k,217) + rxt(k,94) &
                  + rxt(k,95) + rxt(k,527) + rxt(k,532) + rxt(k,538) &
                  + het_rates(k,60))* y(k,60)
         prod(k,177) =rxt(k,190)*y(k,125)*y(k,59)
         loss(k,5) = ( + het_rates(k,61))* y(k,61)
         prod(k,5) = 0._r8
         loss(k,192) = (rxt(k,276)* y(k,217) + het_rates(k,62))* y(k,62)
         prod(k,192) = (rxt(k,32) +rxt(k,33) +rxt(k,176)*y(k,56) +rxt(k,212)*y(k,17) + &
                 rxt(k,257)*y(k,126) +rxt(k,258)*y(k,133) +rxt(k,259)*y(k,217)) &
                 *y(k,42) + (.630_r8*rxt(k,280)*y(k,25) +.560_r8*rxt(k,309)*y(k,29) + &
                 .650_r8*rxt(k,339)*y(k,105) +.560_r8*rxt(k,353)*y(k,111) + &
                 .620_r8*rxt(k,386)*y(k,98) +.230_r8*rxt(k,441)*y(k,6) + &
                 .230_r8*rxt(k,444)*y(k,110))*y(k,134) &
                  + (.220_r8*rxt(k,337)*y(k,209) +.250_r8*rxt(k,394)*y(k,225) + &
                 .170_r8*rxt(k,412)*y(k,199) +.400_r8*rxt(k,415)*y(k,210) + &
                 .350_r8*rxt(k,418)*y(k,212) +.225_r8*rxt(k,453)*y(k,221))*y(k,124) &
                  + (.350_r8*rxt(k,278)*y(k,24) +rxt(k,303)*y(k,75) + &
                 rxt(k,316)*y(k,49) +.700_r8*rxt(k,462)*y(k,178) +rxt(k,466)*y(k,137)) &
                 *y(k,217) + (rxt(k,315)*y(k,49) +.220_r8*rxt(k,336)*y(k,209) + &
                 .500_r8*rxt(k,395)*y(k,225))*y(k,126) &
                  + (.110_r8*rxt(k,334)*y(k,209) +.200_r8*rxt(k,392)*y(k,225) + &
                 .125_r8*rxt(k,451)*y(k,221))*y(k,198) &
                  + (.070_r8*rxt(k,411)*y(k,199) +.160_r8*rxt(k,414)*y(k,210) + &
                 .140_r8*rxt(k,417)*y(k,212))*y(k,203) + (rxt(k,110) + &
                 rxt(k,465)*y(k,133))*y(k,137) + (.220_r8*rxt(k,333)*y(k,209) + &
                 .250_r8*rxt(k,391)*y(k,225))*y(k,197) +1.500_r8*rxt(k,22)*y(k,9) &
                  +.450_r8*rxt(k,23)*y(k,10) +.600_r8*rxt(k,26)*y(k,13) +rxt(k,27) &
                 *y(k,14) +rxt(k,34)*y(k,45) +rxt(k,244)*y(k,56)*y(k,46) +rxt(k,36) &
                 *y(k,49) +.380_r8*rxt(k,39)*y(k,54) +rxt(k,41)*y(k,63) +rxt(k,43) &
                 *y(k,74) +2.000_r8*rxt(k,44)*y(k,75) +.330_r8*rxt(k,45)*y(k,93) &
                  +1.340_r8*rxt(k,51)*y(k,105) +.700_r8*rxt(k,55)*y(k,111) &
                  +1.500_r8*rxt(k,64)*y(k,173) +.250_r8*rxt(k,65)*y(k,174) +rxt(k,68) &
                 *y(k,177) +1.700_r8*rxt(k,69)*y(k,178)
         loss(k,168) = ( + rxt(k,41) + het_rates(k,63))* y(k,63)
         prod(k,168) = (rxt(k,268)*y(k,87) +rxt(k,276)*y(k,62) +rxt(k,296)*y(k,50) + &
                 .500_r8*rxt(k,297)*y(k,51) +.800_r8*rxt(k,302)*y(k,74) + &
                 rxt(k,303)*y(k,75) +.500_r8*rxt(k,352)*y(k,109) + &
                 1.800_r8*rxt(k,462)*y(k,178))*y(k,217) &
                  + (2.000_r8*rxt(k,292)*y(k,197) +.900_r8*rxt(k,293)*y(k,198) + &
                 rxt(k,295)*y(k,124) +2.000_r8*rxt(k,342)*y(k,211) + &
                 rxt(k,366)*y(k,205) +rxt(k,391)*y(k,225))*y(k,197) &
                  + (.200_r8*rxt(k,309)*y(k,29) +.100_r8*rxt(k,353)*y(k,111) + &
                 .270_r8*rxt(k,441)*y(k,6) +.270_r8*rxt(k,444)*y(k,110))*y(k,134) &
                  + (rxt(k,343)*y(k,198) +.450_r8*rxt(k,344)*y(k,203) + &
                 2.000_r8*rxt(k,345)*y(k,211))*y(k,211) &
                  + (.500_r8*rxt(k,451)*y(k,198) +.900_r8*rxt(k,453)*y(k,124)) &
                 *y(k,221) +rxt(k,37)*y(k,51) +.440_r8*rxt(k,39)*y(k,54) &
                  +.400_r8*rxt(k,60)*y(k,139) +rxt(k,65)*y(k,174) +.800_r8*rxt(k,69) &
                 *y(k,178)
         loss(k,85) = (rxt(k,237)* y(k,216) + rxt(k,96) + het_rates(k,64))* y(k,64)
         prod(k,85) = (rxt(k,202)*y(k,34) +rxt(k,204)*y(k,37) + &
                 2.000_r8*rxt(k,205)*y(k,38) +2.000_r8*rxt(k,206)*y(k,39) + &
                 rxt(k,207)*y(k,40) +rxt(k,228)*y(k,35) +2.000_r8*rxt(k,230)*y(k,78) + &
                 rxt(k,254)*y(k,83) +rxt(k,255)*y(k,84))*y(k,216) + (rxt(k,101) + &
                 rxt(k,249)*y(k,217))*y(k,83) + (rxt(k,102) +rxt(k,250)*y(k,217)) &
                 *y(k,84) +rxt(k,79)*y(k,34) +rxt(k,80)*y(k,35) +rxt(k,82)*y(k,37) &
                  +2.000_r8*rxt(k,83)*y(k,38) +2.000_r8*rxt(k,84)*y(k,39) +rxt(k,85) &
                 *y(k,40) +2.000_r8*rxt(k,98)*y(k,78)
         loss(k,83) = (rxt(k,238)* y(k,216) + rxt(k,97) + het_rates(k,65))* y(k,65)
         prod(k,83) = (rxt(k,100) +rxt(k,248)*y(k,217) +rxt(k,253)*y(k,216))*y(k,82) &
                  + (rxt(k,81) +rxt(k,203)*y(k,216))*y(k,36) + (rxt(k,82) + &
                 rxt(k,204)*y(k,216))*y(k,37)
         loss(k,77) = (rxt(k,410)* y(k,217) + het_rates(k,66))* y(k,66)
         prod(k,77) =.180_r8*rxt(k,430)*y(k,217)*y(k,180)
         loss(k,100) = (rxt(k,463)* y(k,126) + (rxt(k,464) +rxt(k,478))* y(k,217) &
                  + het_rates(k,67))* y(k,67)
         prod(k,100) = 0._r8
         loss(k,6) = ( + het_rates(k,68))* y(k,68)
         prod(k,6) = 0._r8
         loss(k,7) = ( + het_rates(k,69))* y(k,69)
         prod(k,7) = 0._r8
         loss(k,8) = ( + het_rates(k,70))* y(k,70)
         prod(k,8) = 0._r8
         loss(k,9) = ( + rxt(k,124) + het_rates(k,71))* y(k,71)
         prod(k,9) = 0._r8
         loss(k,62) = ( + rxt(k,42) + het_rates(k,72))* y(k,72)
         prod(k,62) =rxt(k,298)*y(k,203)*y(k,202)
         loss(k,175) = (rxt(k,233)* y(k,54) +rxt(k,234)* y(k,77) +rxt(k,236)* y(k,89) &
                  +rxt(k,235)* y(k,229) + het_rates(k,73))* y(k,73)
         prod(k,175) = (rxt(k,206)*y(k,39) +rxt(k,228)*y(k,35) + &
                 2.000_r8*rxt(k,237)*y(k,64) +rxt(k,238)*y(k,65))*y(k,216) +rxt(k,80) &
                 *y(k,35) +rxt(k,84)*y(k,39) +2.000_r8*rxt(k,96)*y(k,64) +rxt(k,97) &
                 *y(k,65) +rxt(k,104)*y(k,88)
         loss(k,191) = (rxt(k,302)* y(k,217) + rxt(k,43) + het_rates(k,74))* y(k,74)
         prod(k,191) = (.530_r8*rxt(k,337)*y(k,209) +.050_r8*rxt(k,375)*y(k,206) + &
                 .250_r8*rxt(k,394)*y(k,225) +.225_r8*rxt(k,453)*y(k,221))*y(k,124) &
                  + (.530_r8*rxt(k,336)*y(k,209) +.050_r8*rxt(k,376)*y(k,206) + &
                 .250_r8*rxt(k,395)*y(k,225))*y(k,126) &
                  + (.260_r8*rxt(k,334)*y(k,209) +.100_r8*rxt(k,392)*y(k,225) + &
                 .125_r8*rxt(k,451)*y(k,221))*y(k,198) + (.700_r8*rxt(k,377)*y(k,99) + &
                 .500_r8*rxt(k,378)*y(k,100) +rxt(k,389)*y(k,115))*y(k,217) &
                  + (.530_r8*rxt(k,333)*y(k,209) +.250_r8*rxt(k,391)*y(k,225)) &
                 *y(k,197) +.330_r8*rxt(k,45)*y(k,93) +.250_r8*rxt(k,65)*y(k,174) &
                  +rxt(k,301)*y(k,201)
         loss(k,184) = (rxt(k,303)* y(k,217) + rxt(k,44) + rxt(k,519) &
                  + het_rates(k,75))* y(k,75)
         prod(k,184) = (.050_r8*rxt(k,375)*y(k,206) +.250_r8*rxt(k,394)*y(k,225) + &
                 rxt(k,401)*y(k,191) +.400_r8*rxt(k,415)*y(k,210) + &
                 .170_r8*rxt(k,418)*y(k,212) +.700_r8*rxt(k,421)*y(k,218) + &
                 .600_r8*rxt(k,428)*y(k,223) +.340_r8*rxt(k,434)*y(k,226) + &
                 .170_r8*rxt(k,437)*y(k,228))*y(k,124) + (.650_r8*rxt(k,278)*y(k,24) + &
                 .200_r8*rxt(k,302)*y(k,74) +rxt(k,390)*y(k,116))*y(k,217) &
                  + (.250_r8*rxt(k,391)*y(k,197) +.100_r8*rxt(k,392)*y(k,198) + &
                 .250_r8*rxt(k,395)*y(k,126))*y(k,225) &
                  + (.160_r8*rxt(k,414)*y(k,210) +.070_r8*rxt(k,417)*y(k,212)) &
                 *y(k,203) +rxt(k,21)*y(k,8) +.130_r8*rxt(k,23)*y(k,10) &
                  +.050_r8*rxt(k,376)*y(k,206)*y(k,126) +.700_r8*rxt(k,61)*y(k,143) &
                  +.600_r8*rxt(k,70)*y(k,179) +.340_r8*rxt(k,72)*y(k,183) &
                  +.170_r8*rxt(k,73)*y(k,185)
         loss(k,212) = (rxt(k,142)* y(k,134) + (rxt(k,136) +rxt(k,137) +rxt(k,138)) &
                 * y(k,203) + rxt(k,139) + het_rates(k,76))* y(k,76)
         prod(k,212) = (rxt(k,143)*y(k,77) +rxt(k,146)*y(k,133) +rxt(k,164)*y(k,112) + &
                 rxt(k,259)*y(k,42) +rxt(k,466)*y(k,137) +rxt(k,472)*y(k,148) + &
                 rxt(k,477)*y(k,150))*y(k,217) + (rxt(k,125)*y(k,216) + &
                 rxt(k,134)*y(k,133) +rxt(k,178)*y(k,56) +rxt(k,234)*y(k,73))*y(k,77) &
                  + (.330_r8*rxt(k,39) +rxt(k,40) +rxt(k,273)*y(k,216))*y(k,54) &
                  + (rxt(k,99) +rxt(k,232)*y(k,216))*y(k,81) + (rxt(k,103) + &
                 rxt(k,209)*y(k,216))*y(k,85) + (rxt(k,2) +2.000_r8*rxt(k,3))*y(k,229) &
                  +2.000_r8*rxt(k,33)*y(k,42) +rxt(k,38)*y(k,53) +rxt(k,104)*y(k,88)
         loss(k,208) = (rxt(k,178)* y(k,56) +rxt(k,234)* y(k,73) +rxt(k,134)* y(k,133) &
                  +rxt(k,125)* y(k,216) +rxt(k,143)* y(k,217) + het_rates(k,77)) &
                 * y(k,77)
         prod(k,208) = (1.440_r8*rxt(k,39) +rxt(k,274)*y(k,216))*y(k,54) +rxt(k,32) &
                 *y(k,42) +rxt(k,136)*y(k,203)*y(k,76) +rxt(k,1)*y(k,229)
         loss(k,56) = (rxt(k,230)* y(k,216) + rxt(k,98) + het_rates(k,78))* y(k,78)
         prod(k,56) = 0._r8
         loss(k,147) = (rxt(k,179)* y(k,56) +rxt(k,135)* y(k,133) +rxt(k,144) &
                 * y(k,217) + rxt(k,4) + het_rates(k,79))* y(k,79)
         prod(k,147) =rxt(k,150)*y(k,203)*y(k,203) +rxt(k,149)*y(k,217)*y(k,217)
         loss(k,63) = ( + rxt(k,109) + het_rates(k,80))* y(k,80)
         prod(k,63) =rxt(k,479)*y(k,229)*y(k,152)
         loss(k,169) = (rxt(k,225)* y(k,133) + (rxt(k,231) +rxt(k,232))* y(k,216) &
                  +rxt(k,226)* y(k,217) + rxt(k,99) + het_rates(k,81))* y(k,81)
         prod(k,169) = (rxt(k,212)*y(k,42) +rxt(k,213)*y(k,203))*y(k,17)
         loss(k,82) = (rxt(k,253)* y(k,216) +rxt(k,248)* y(k,217) + rxt(k,100) &
                  + het_rates(k,82))* y(k,82)
         prod(k,82) = 0._r8
         loss(k,89) = (rxt(k,254)* y(k,216) +rxt(k,249)* y(k,217) + rxt(k,101) &
                  + het_rates(k,83))* y(k,83)
         prod(k,89) = 0._r8
         loss(k,101) = (rxt(k,255)* y(k,216) +rxt(k,250)* y(k,217) + rxt(k,102) &
                  + het_rates(k,84))* y(k,84)
         prod(k,101) = 0._r8
         loss(k,223) = ((rxt(k,529) +rxt(k,534) +rxt(k,539))* y(k,60) + (rxt(k,531) + &
                 rxt(k,536))* y(k,91) + (rxt(k,524) +rxt(k,530) +rxt(k,535))* y(k,92) &
                  +rxt(k,196)* y(k,133) + (rxt(k,208) +rxt(k,209))* y(k,216) &
                  +rxt(k,197)* y(k,217) + rxt(k,103) + het_rates(k,85))* y(k,85)
         prod(k,223) = (rxt(k,176)*y(k,42) +rxt(k,177)*y(k,54) +rxt(k,178)*y(k,77) + &
                 rxt(k,179)*y(k,79) +rxt(k,180)*y(k,203) +rxt(k,198)*y(k,92) + &
                 rxt(k,239)*y(k,41) +rxt(k,241)*y(k,43) +2.000_r8*rxt(k,244)*y(k,46) + &
                 rxt(k,246)*y(k,55) +rxt(k,287)*y(k,28))*y(k,56) +rxt(k,195)*y(k,217) &
                 *y(k,59)
         loss(k,74) = (rxt(k,275)* y(k,216) +rxt(k,267)* y(k,217) + het_rates(k,86)) &
                 * y(k,86)
         prod(k,74) = 0._r8
         loss(k,181) = (rxt(k,268)* y(k,217) + het_rates(k,87))* y(k,87)
         prod(k,181) = (.370_r8*rxt(k,280)*y(k,25) +.120_r8*rxt(k,309)*y(k,29) + &
                 .330_r8*rxt(k,339)*y(k,105) +.120_r8*rxt(k,353)*y(k,111) + &
                 .110_r8*rxt(k,386)*y(k,98) +.050_r8*rxt(k,441)*y(k,6) + &
                 .050_r8*rxt(k,444)*y(k,110))*y(k,134) + (rxt(k,269)*y(k,203) + &
                 rxt(k,271)*y(k,124))*y(k,204) +.350_r8*rxt(k,278)*y(k,217)*y(k,24)
         loss(k,97) = ( + rxt(k,104) + het_rates(k,88))* y(k,88)
         prod(k,97) = (rxt(k,233)*y(k,54) +rxt(k,234)*y(k,77) +rxt(k,235)*y(k,229) + &
                 rxt(k,236)*y(k,89))*y(k,73)
         loss(k,211) = (rxt(k,236)* y(k,73) +rxt(k,173)* y(k,217) + rxt(k,9) &
                  + het_rates(k,89))* y(k,89)
         prod(k,211) = (rxt(k,527) +rxt(k,532) +rxt(k,538) +rxt(k,529)*y(k,85) + &
                 rxt(k,534)*y(k,85) +rxt(k,539)*y(k,85))*y(k,60) + (rxt(k,490) + &
                 rxt(k,257)*y(k,42) +rxt(k,289)*y(k,45) +rxt(k,315)*y(k,49) + &
                 rxt(k,463)*y(k,67))*y(k,126) + (2.000_r8*rxt(k,485) + &
                 2.000_r8*rxt(k,523) +2.000_r8*rxt(k,526) +2.000_r8*rxt(k,537)) &
                 *y(k,114) + (rxt(k,525) +rxt(k,528) +rxt(k,533))*y(k,20) &
                  + (.500_r8*rxt(k,489) +rxt(k,172)*y(k,217))*y(k,125) +rxt(k,482) &
                 *y(k,93) +rxt(k,483)*y(k,99) +rxt(k,484)*y(k,100) +rxt(k,486) &
                 *y(k,115) +rxt(k,487)*y(k,116) +rxt(k,491)*y(k,128) +rxt(k,492) &
                 *y(k,138) +rxt(k,493)*y(k,175)
         loss(k,118) = (rxt(k,151)* y(k,217) + rxt(k,10) + rxt(k,11) + rxt(k,174) &
                  + het_rates(k,90))* y(k,90)
         prod(k,118) =rxt(k,170)*y(k,203)*y(k,125)
         loss(k,165) = ((rxt(k,531) +rxt(k,536))* y(k,85) +rxt(k,227)* y(k,133) &
                  + rxt(k,105) + het_rates(k,91))* y(k,91)
         prod(k,165) = (rxt(k,525) +rxt(k,528) +rxt(k,533))*y(k,20) &
                  +rxt(k,219)*y(k,203)*y(k,19)
         loss(k,171) = (rxt(k,198)* y(k,56) + (rxt(k,524) +rxt(k,530) +rxt(k,535)) &
                 * y(k,85) +rxt(k,199)* y(k,133) +rxt(k,200)* y(k,217) + rxt(k,106) &
                  + het_rates(k,92))* y(k,92)
         prod(k,171) = (rxt(k,527) +rxt(k,532) +rxt(k,538) +rxt(k,192)*y(k,217)) &
                 *y(k,60) +rxt(k,187)*y(k,203)*y(k,59)
         loss(k,196) = (rxt(k,332)* y(k,217) + rxt(k,45) + rxt(k,482) &
                  + het_rates(k,93))* y(k,93)
         prod(k,196) = (rxt(k,331)*y(k,200) +rxt(k,338)*y(k,209))*y(k,124) &
                  + (.300_r8*rxt(k,377)*y(k,99) +.500_r8*rxt(k,378)*y(k,100))*y(k,217)
         loss(k,84) = (rxt(k,363)* y(k,217) + rxt(k,46) + het_rates(k,94))* y(k,94)
         prod(k,84) =rxt(k,374)*y(k,206)
         loss(k,195) = (rxt(k,317)* y(k,217) + rxt(k,47) + het_rates(k,95))* y(k,95)
         prod(k,195) = (.220_r8*rxt(k,333)*y(k,197) +.230_r8*rxt(k,334)*y(k,198) + &
                 .220_r8*rxt(k,336)*y(k,126) +.220_r8*rxt(k,337)*y(k,124))*y(k,209) &
                  + (.500_r8*rxt(k,321)*y(k,146) +.500_r8*rxt(k,352)*y(k,109) + &
                 .700_r8*rxt(k,377)*y(k,99) +.500_r8*rxt(k,378)*y(k,100))*y(k,217) &
                  + (.250_r8*rxt(k,391)*y(k,197) +.100_r8*rxt(k,392)*y(k,198) + &
                 .250_r8*rxt(k,394)*y(k,124) +.250_r8*rxt(k,395)*y(k,126))*y(k,225) &
                  + (.050_r8*rxt(k,375)*y(k,124) +.050_r8*rxt(k,376)*y(k,126)) &
                 *y(k,206) +.170_r8*rxt(k,45)*y(k,93) +.200_r8*rxt(k,322)*y(k,220) &
                 *y(k,198)
         loss(k,104) = (rxt(k,364)* y(k,217) + het_rates(k,96))* y(k,96)
         prod(k,104) = (rxt(k,371)*y(k,197) +.750_r8*rxt(k,372)*y(k,198) + &
                 .870_r8*rxt(k,375)*y(k,124) +.950_r8*rxt(k,376)*y(k,126))*y(k,206)
         loss(k,64) = (rxt(k,365)* y(k,217) + het_rates(k,97))* y(k,97)
         prod(k,64) =.600_r8*rxt(k,388)*y(k,217)*y(k,103)
         loss(k,173) = (rxt(k,379)* y(k,126) +rxt(k,386)* y(k,134) +rxt(k,387) &
                 * y(k,217) + het_rates(k,98))* y(k,98)
         prod(k,173) = 0._r8
         loss(k,144) = (rxt(k,377)* y(k,217) + rxt(k,483) + het_rates(k,99))* y(k,99)
         prod(k,144) =.080_r8*rxt(k,369)*y(k,205)*y(k,124)
         loss(k,141) = (rxt(k,378)* y(k,217) + rxt(k,484) + het_rates(k,100)) &
                 * y(k,100)
         prod(k,141) =.080_r8*rxt(k,375)*y(k,206)*y(k,124)
         loss(k,198) = (rxt(k,383)* y(k,124) +rxt(k,384)* y(k,126) +rxt(k,380) &
                 * y(k,197) +rxt(k,381)* y(k,198) +rxt(k,382)* y(k,203) &
                  + het_rates(k,101))* y(k,101)
         prod(k,198) =rxt(k,379)*y(k,126)*y(k,98)
         loss(k,117) = (rxt(k,385)* y(k,217) + rxt(k,48) + het_rates(k,102))* y(k,102)
         prod(k,117) =rxt(k,382)*y(k,203)*y(k,101)
         loss(k,157) = (rxt(k,388)* y(k,217) + rxt(k,49) + het_rates(k,103))* y(k,103)
         prod(k,157) = (rxt(k,368)*y(k,205) +rxt(k,373)*y(k,206))*y(k,203) +rxt(k,48) &
                 *y(k,102)
         loss(k,48) = (rxt(k,509)* y(k,217) + het_rates(k,104))* y(k,104)
         prod(k,48) = 0._r8
         loss(k,199) = (rxt(k,339)* y(k,134) +rxt(k,340)* y(k,217) + rxt(k,50) &
                  + rxt(k,51) + het_rates(k,105))* y(k,105)
         prod(k,199) = (.390_r8*rxt(k,366)*y(k,197) +.310_r8*rxt(k,367)*y(k,198) + &
                 .360_r8*rxt(k,369)*y(k,124) +.400_r8*rxt(k,370)*y(k,126))*y(k,205) &
                  +.300_r8*rxt(k,386)*y(k,134)*y(k,98) +.300_r8*rxt(k,49)*y(k,103)
         loss(k,102) = (rxt(k,341)* y(k,217) + het_rates(k,106))* y(k,106)
         prod(k,102) =rxt(k,335)*y(k,209)*y(k,203)
         loss(k,134) = (rxt(k,350)* y(k,217) + rxt(k,52) + het_rates(k,107))* y(k,107)
         prod(k,134) =.800_r8*rxt(k,19)*y(k,1) +.800_r8*rxt(k,20)*y(k,2) &
                  +.800_r8*rxt(k,359)*y(k,189)*y(k,124)
         loss(k,103) = (rxt(k,351)* y(k,217) + rxt(k,53) + het_rates(k,108))* y(k,108)
         prod(k,103) =.800_r8*rxt(k,348)*y(k,213)*y(k,203)
         loss(k,143) = (rxt(k,352)* y(k,217) + rxt(k,54) + rxt(k,356) &
                  + het_rates(k,109))* y(k,109)
         prod(k,143) =rxt(k,355)*y(k,211)*y(k,125)
         loss(k,182) = (rxt(k,443)* y(k,126) +rxt(k,444)* y(k,134) +rxt(k,445) &
                 * y(k,217) + het_rates(k,110))* y(k,110)
         prod(k,182) = 0._r8
         loss(k,205) = (rxt(k,353)* y(k,134) +rxt(k,354)* y(k,217) + rxt(k,55) &
                  + het_rates(k,111))* y(k,111)
         prod(k,205) = (.610_r8*rxt(k,366)*y(k,197) +.440_r8*rxt(k,367)*y(k,198) + &
                 .560_r8*rxt(k,369)*y(k,124) +.600_r8*rxt(k,370)*y(k,126))*y(k,205) &
                  +.200_r8*rxt(k,386)*y(k,134)*y(k,98) +.700_r8*rxt(k,49)*y(k,103)
         loss(k,132) = (rxt(k,152)* y(k,124) + (rxt(k,153) +rxt(k,154) +rxt(k,155)) &
                 * y(k,125) +rxt(k,164)* y(k,217) + rxt(k,156) + het_rates(k,112)) &
                 * y(k,112)
         prod(k,132) =rxt(k,15)*y(k,124)
         loss(k,75) = ((rxt(k,168) +rxt(k,169))* y(k,216) + rxt(k,12) &
                  + het_rates(k,113))* y(k,113)
         prod(k,75) =rxt(k,153)*y(k,125)*y(k,112)
         loss(k,95) = ( + rxt(k,13) + rxt(k,14) + rxt(k,175) + rxt(k,485) + rxt(k,523) &
                  + rxt(k,526) + rxt(k,537) + het_rates(k,114))* y(k,114)
         prod(k,95) =rxt(k,171)*y(k,126)*y(k,125)
         loss(k,113) = (rxt(k,389)* y(k,217) + rxt(k,486) + het_rates(k,115)) &
                 * y(k,115)
         prod(k,113) =.200_r8*rxt(k,381)*y(k,198)*y(k,101)
         loss(k,189) = (rxt(k,390)* y(k,217) + rxt(k,56) + rxt(k,487) &
                  + het_rates(k,116))* y(k,116)
         prod(k,189) = (rxt(k,380)*y(k,197) +.800_r8*rxt(k,381)*y(k,198) + &
                 rxt(k,383)*y(k,124) +rxt(k,384)*y(k,126))*y(k,101)
         loss(k,10) = ( + het_rates(k,117))* y(k,117)
         prod(k,10) = 0._r8
         loss(k,11) = ( + het_rates(k,118))* y(k,118)
         prod(k,11) = 0._r8
         loss(k,12) = ( + het_rates(k,119))* y(k,119)
         prod(k,12) = 0._r8
         loss(k,53) = (rxt(k,480)* y(k,217) + het_rates(k,120))* y(k,120)
         prod(k,53) = 0._r8
         loss(k,13) = ( + rxt(k,488) + het_rates(k,121))* y(k,121)
         prod(k,13) = 0._r8
         loss(k,14) = ( + rxt(k,541) + het_rates(k,122))* y(k,122)
         prod(k,14) = 0._r8
         loss(k,15) = ( + rxt(k,540) + het_rates(k,123))* y(k,123)
         prod(k,15) = 0._r8
         loss(k,217) = (rxt(k,220)* y(k,19) +rxt(k,188)* y(k,59) +rxt(k,383)* y(k,101) &
                  +rxt(k,152)* y(k,112) +rxt(k,161)* y(k,126) +rxt(k,167)* y(k,133) &
                  +rxt(k,166)* y(k,134) +rxt(k,398)* y(k,188) + (rxt(k,359) + &
                 rxt(k,360))* y(k,189) +rxt(k,401)* y(k,191) +rxt(k,406)* y(k,193) &
                  +rxt(k,284)* y(k,194) +rxt(k,312)* y(k,195) +rxt(k,408)* y(k,196) &
                  +rxt(k,295)* y(k,197) +rxt(k,263)* y(k,198) +rxt(k,412)* y(k,199) &
                  + (rxt(k,330) +rxt(k,331))* y(k,200) +rxt(k,299)* y(k,202) &
                  +rxt(k,165)* y(k,203) +rxt(k,271)* y(k,204) +rxt(k,369)* y(k,205) &
                  +rxt(k,375)* y(k,206) + (rxt(k,337) +rxt(k,338))* y(k,209) &
                  +rxt(k,415)* y(k,210) +rxt(k,346)* y(k,211) +rxt(k,418)* y(k,212) &
                  +rxt(k,349)* y(k,213) +rxt(k,448)* y(k,215) +rxt(k,421)* y(k,218) &
                  +rxt(k,320)* y(k,219) +rxt(k,324)* y(k,220) +rxt(k,453)* y(k,221) &
                  +rxt(k,458)* y(k,222) +rxt(k,428)* y(k,223) +rxt(k,394)* y(k,225) &
                  +rxt(k,434)* y(k,226) +rxt(k,437)* y(k,228) + rxt(k,15) &
                  + het_rates(k,124))* y(k,124)
         prod(k,217) = (rxt(k,16) +.500_r8*rxt(k,489) +2.000_r8*rxt(k,154)*y(k,112) + &
                 rxt(k,157)*y(k,133) +rxt(k,473)*y(k,150))*y(k,125) + (rxt(k,156) + &
                 rxt(k,164)*y(k,217))*y(k,112) +2.000_r8*rxt(k,168)*y(k,216)*y(k,113) &
                  +rxt(k,14)*y(k,114) +rxt(k,17)*y(k,126)
         loss(k,224) = (rxt(k,221)* y(k,19) +rxt(k,190)* y(k,59) + (rxt(k,153) + &
                 rxt(k,154) +rxt(k,155))* y(k,112) +rxt(k,171)* y(k,126) &
                  + (rxt(k,157) +rxt(k,159))* y(k,133) +rxt(k,158)* y(k,134) &
                  +rxt(k,423)* y(k,141) +rxt(k,473)* y(k,150) +rxt(k,426)* y(k,188) &
                  +rxt(k,306)* y(k,197) +rxt(k,413)* y(k,199) +rxt(k,170)* y(k,203) &
                  +rxt(k,416)* y(k,210) +rxt(k,355)* y(k,211) +rxt(k,419)* y(k,212) &
                  +rxt(k,172)* y(k,217) + rxt(k,16) + rxt(k,489) + het_rates(k,125)) &
                 * y(k,125)
         prod(k,224) = (2.000_r8*rxt(k,161)*y(k,126) +rxt(k,165)*y(k,203) + &
                 rxt(k,166)*y(k,134) +rxt(k,167)*y(k,133) +rxt(k,188)*y(k,59) + &
                 rxt(k,220)*y(k,19) +rxt(k,263)*y(k,198) +rxt(k,271)*y(k,204) + &
                 rxt(k,284)*y(k,194) +rxt(k,295)*y(k,197) +rxt(k,299)*y(k,202) + &
                 rxt(k,312)*y(k,195) +rxt(k,320)*y(k,219) +rxt(k,324)*y(k,220) + &
                 rxt(k,330)*y(k,200) +rxt(k,337)*y(k,209) +rxt(k,346)*y(k,211) + &
                 rxt(k,349)*y(k,213) +rxt(k,359)*y(k,189) + &
                 .920_r8*rxt(k,369)*y(k,205) +.920_r8*rxt(k,375)*y(k,206) + &
                 rxt(k,383)*y(k,101) +rxt(k,394)*y(k,225) +rxt(k,398)*y(k,188) + &
                 rxt(k,401)*y(k,191) +rxt(k,406)*y(k,193) +rxt(k,408)*y(k,196) + &
                 rxt(k,412)*y(k,199) +rxt(k,415)*y(k,210) +rxt(k,418)*y(k,212) + &
                 rxt(k,421)*y(k,218) +rxt(k,428)*y(k,223) +rxt(k,434)*y(k,226) + &
                 rxt(k,437)*y(k,228) +1.600_r8*rxt(k,448)*y(k,215) + &
                 .900_r8*rxt(k,453)*y(k,221) +.800_r8*rxt(k,458)*y(k,222))*y(k,124) &
                  + (rxt(k,18) +rxt(k,160)*y(k,203) +rxt(k,162)*y(k,133) + &
                 rxt(k,163)*y(k,217) +rxt(k,328)*y(k,16) +rxt(k,336)*y(k,209) + &
                 rxt(k,347)*y(k,211) +rxt(k,370)*y(k,205) +rxt(k,376)*y(k,206) + &
                 rxt(k,384)*y(k,101) +rxt(k,395)*y(k,225) + &
                 2.000_r8*rxt(k,449)*y(k,215))*y(k,126) + (rxt(k,151)*y(k,90) + &
                 rxt(k,318)*y(k,127) +rxt(k,357)*y(k,1) +.700_r8*rxt(k,377)*y(k,99) + &
                 rxt(k,455)*y(k,175))*y(k,217) + (rxt(k,11) +rxt(k,174))*y(k,90) &
                  + (rxt(k,54) +rxt(k,356))*y(k,109) + (rxt(k,13) +rxt(k,175)) &
                 *y(k,114) + (.600_r8*rxt(k,60) +rxt(k,307))*y(k,139) +rxt(k,19) &
                 *y(k,1) +rxt(k,76)*y(k,20) +rxt(k,95)*y(k,60) +rxt(k,9)*y(k,89) &
                  +rxt(k,45)*y(k,93) +rxt(k,48)*y(k,102) +rxt(k,56)*y(k,116) &
                  +rxt(k,57)*y(k,127) +rxt(k,58)*y(k,128) +rxt(k,59)*y(k,138) &
                  +rxt(k,431)*y(k,140) +rxt(k,66)*y(k,175) &
                  +.500_r8*rxt(k,446)*y(k,215)*y(k,198)
         loss(k,216) = (rxt(k,440)* y(k,6) +rxt(k,328)* y(k,16) +rxt(k,308)* y(k,29) &
                  +rxt(k,257)* y(k,42) +rxt(k,289)* y(k,45) +rxt(k,315)* y(k,49) &
                  +rxt(k,463)* y(k,67) +rxt(k,379)* y(k,98) +rxt(k,384)* y(k,101) &
                  +rxt(k,443)* y(k,110) +rxt(k,161)* y(k,124) +rxt(k,171)* y(k,125) &
                  +rxt(k,162)* y(k,133) +rxt(k,460)* y(k,177) +rxt(k,160)* y(k,203) &
                  +rxt(k,370)* y(k,205) +rxt(k,376)* y(k,206) +rxt(k,336)* y(k,209) &
                  +rxt(k,347)* y(k,211) +rxt(k,449)* y(k,215) +rxt(k,163)* y(k,217) &
                  +rxt(k,395)* y(k,225) + rxt(k,17) + rxt(k,18) + rxt(k,490) &
                  + het_rates(k,126))* y(k,126)
         prod(k,216) = (rxt(k,94) +rxt(k,189)*y(k,56) +rxt(k,191)*y(k,133) + &
                 rxt(k,192)*y(k,217))*y(k,60) + (rxt(k,13) +rxt(k,14) +rxt(k,175)) &
                 *y(k,114) + (rxt(k,173)*y(k,89) +rxt(k,304)*y(k,139) + &
                 .500_r8*rxt(k,352)*y(k,109))*y(k,217) + (rxt(k,77) + &
                 rxt(k,222)*y(k,133))*y(k,20) + (rxt(k,158)*y(k,134) + &
                 rxt(k,159)*y(k,133))*y(k,125) +rxt(k,236)*y(k,89)*y(k,73) +rxt(k,10) &
                 *y(k,90) +.400_r8*rxt(k,60)*y(k,139)
         loss(k,174) = (rxt(k,318)* y(k,217) + rxt(k,57) + het_rates(k,127))* y(k,127)
         prod(k,174) = (.500_r8*rxt(k,378)*y(k,100) +rxt(k,385)*y(k,102) + &
                 rxt(k,389)*y(k,115) +rxt(k,390)*y(k,116))*y(k,217) &
                  +rxt(k,308)*y(k,126)*y(k,29)
         loss(k,116) = (rxt(k,450)* y(k,217) + rxt(k,58) + rxt(k,491) &
                  + het_rates(k,128))* y(k,128)
         prod(k,116) =rxt(k,447)*y(k,215)*y(k,203)
         loss(k,16) = ( + het_rates(k,129))* y(k,129)
         prod(k,16) = 0._r8
         loss(k,17) = ( + het_rates(k,130))* y(k,130)
         prod(k,17) = 0._r8
         loss(k,18) = ( + het_rates(k,131))* y(k,131)
         prod(k,18) = 0._r8
         loss(k,19) = ( + het_rates(k,132))* y(k,132)
         prod(k,19) = 0._r8
         loss(k,226) = (rxt(k,223)* y(k,19) +rxt(k,222)* y(k,20) +rxt(k,258)* y(k,42) &
                  +rxt(k,193)* y(k,59) +rxt(k,191)* y(k,60) +rxt(k,134)* y(k,77) &
                  +rxt(k,135)* y(k,79) +rxt(k,225)* y(k,81) +rxt(k,196)* y(k,85) &
                  +rxt(k,227)* y(k,91) +rxt(k,199)* y(k,92) +rxt(k,167)* y(k,124) &
                  + (rxt(k,157) +rxt(k,159))* y(k,125) +rxt(k,162)* y(k,126) &
                  + 2._r8*rxt(k,132)* y(k,133) +rxt(k,131)* y(k,134) +rxt(k,465) &
                 * y(k,137) +rxt(k,140)* y(k,203) +rxt(k,146)* y(k,217) + rxt(k,133) &
                  + het_rates(k,133))* y(k,133)
         prod(k,226) = (rxt(k,156) +rxt(k,152)*y(k,124) +rxt(k,153)*y(k,125))*y(k,112) &
                  + (rxt(k,127) +rxt(k,128) +2.000_r8*rxt(k,130)*y(k,134))*y(k,216) &
                  + (rxt(k,111) +rxt(k,474))*y(k,150) +rxt(k,75)*y(k,19) &
                  +.180_r8*rxt(k,39)*y(k,54) +rxt(k,93)*y(k,59) +rxt(k,41)*y(k,63) &
                  +rxt(k,138)*y(k,203)*y(k,76) +rxt(k,14)*y(k,114) +rxt(k,15)*y(k,124) &
                  +rxt(k,16)*y(k,125) +rxt(k,18)*y(k,126) +rxt(k,8)*y(k,134) &
                  +rxt(k,107)*y(k,136) +rxt(k,467)*y(k,148) +rxt(k,112)*y(k,151) &
                  +rxt(k,113)*y(k,152) +rxt(k,148)*y(k,217)*y(k,217) +rxt(k,3) &
                 *y(k,229)
         loss(k,222) = (rxt(k,441)* y(k,6) +rxt(k,214)* y(k,17) +rxt(k,280)* y(k,25) &
                  +rxt(k,309)* y(k,29) +rxt(k,182)* y(k,56) +rxt(k,142)* y(k,76) &
                  +rxt(k,386)* y(k,98) +rxt(k,339)* y(k,105) +rxt(k,444)* y(k,110) &
                  +rxt(k,353)* y(k,111) +rxt(k,166)* y(k,124) +rxt(k,158)* y(k,125) &
                  +rxt(k,131)* y(k,133) +rxt(k,424)* y(k,141) +rxt(k,469)* y(k,148) &
                  +rxt(k,475)* y(k,150) +rxt(k,141)* y(k,203) + (rxt(k,129) + &
                 rxt(k,130))* y(k,216) +rxt(k,147)* y(k,217) + rxt(k,7) + rxt(k,8) &
                  + het_rates(k,134))* y(k,134)
         prod(k,222) = (.150_r8*rxt(k,294)*y(k,197) +.150_r8*rxt(k,344)*y(k,211)) &
                 *y(k,203) +rxt(k,133)*y(k,133)
         loss(k,20) = ( + het_rates(k,135))* y(k,135)
         prod(k,20) = 0._r8
         loss(k,106) = (rxt(k,476)* y(k,150) + rxt(k,107) + het_rates(k,136)) &
                 * y(k,136)
         prod(k,106) = (rxt(k,186)*y(k,59) +rxt(k,216)*y(k,19))*y(k,59)
         loss(k,111) = (rxt(k,465)* y(k,133) +rxt(k,466)* y(k,217) + rxt(k,110) &
                  + het_rates(k,137))* y(k,137)
         prod(k,111) = 0._r8
         loss(k,88) = ( + rxt(k,59) + rxt(k,492) + het_rates(k,138))* y(k,138)
         prod(k,88) =rxt(k,332)*y(k,217)*y(k,93) +.100_r8*rxt(k,453)*y(k,221)*y(k,124)
         loss(k,137) = (rxt(k,304)* y(k,217) + rxt(k,60) + rxt(k,307) &
                  + het_rates(k,139))* y(k,139)
         prod(k,137) =rxt(k,306)*y(k,197)*y(k,125)
         loss(k,65) = ( + rxt(k,431) + het_rates(k,140))* y(k,140)
         prod(k,65) =rxt(k,426)*y(k,188)*y(k,125)
         loss(k,128) = (rxt(k,423)* y(k,125) +rxt(k,424)* y(k,134) + het_rates(k,141)) &
                 * y(k,141)
         prod(k,128) = (.070_r8*rxt(k,410)*y(k,66) +.060_r8*rxt(k,422)*y(k,142) + &
                 .070_r8*rxt(k,438)*y(k,184))*y(k,217) +rxt(k,31)*y(k,32) &
                  +rxt(k,408)*y(k,196)*y(k,124)
         loss(k,73) = (rxt(k,422)* y(k,217) + het_rates(k,142))* y(k,142)
         prod(k,73) =.530_r8*rxt(k,399)*y(k,217)*y(k,7)
         loss(k,107) = (rxt(k,425)* y(k,217) + rxt(k,61) + het_rates(k,143))* y(k,143)
         prod(k,107) =rxt(k,420)*y(k,218)*y(k,203)
         loss(k,21) = ( + het_rates(k,144))* y(k,144)
         prod(k,21) = 0._r8
         loss(k,22) = ( + het_rates(k,145))* y(k,145)
         prod(k,22) = 0._r8
         loss(k,138) = (rxt(k,321)* y(k,217) + rxt(k,62) + het_rates(k,146))* y(k,146)
         prod(k,138) =rxt(k,319)*y(k,219)*y(k,203)
         loss(k,119) = (rxt(k,325)* y(k,217) + rxt(k,63) + het_rates(k,147))* y(k,147)
         prod(k,119) =.850_r8*rxt(k,323)*y(k,220)*y(k,203)
         loss(k,135) = (rxt(k,469)* y(k,134) +rxt(k,472)* y(k,217) + rxt(k,467) &
                  + het_rates(k,148))* y(k,148)
         prod(k,135) =rxt(k,110)*y(k,137) +rxt(k,111)*y(k,150)
         loss(k,23) = ( + rxt(k,108) + het_rates(k,149))* y(k,149)
         prod(k,23) = 0._r8
         loss(k,201) = (rxt(k,470)* y(k,19) +rxt(k,471)* y(k,59) +rxt(k,473)* y(k,125) &
                  +rxt(k,475)* y(k,134) +rxt(k,476)* y(k,136) +rxt(k,477)* y(k,217) &
                  + rxt(k,111) + rxt(k,474) + het_rates(k,150))* y(k,150)
         prod(k,201) = (rxt(k,467) +rxt(k,469)*y(k,134) +rxt(k,472)*y(k,217))*y(k,148) &
                  +rxt(k,465)*y(k,137)*y(k,133) +rxt(k,112)*y(k,151)
         loss(k,172) = (rxt(k,468)* y(k,217) + rxt(k,112) + het_rates(k,151)) &
                 * y(k,151)
         prod(k,172) = (rxt(k,474) +rxt(k,470)*y(k,19) +rxt(k,471)*y(k,59) + &
                 rxt(k,473)*y(k,125) +rxt(k,475)*y(k,134) +rxt(k,476)*y(k,136) + &
                 rxt(k,477)*y(k,217))*y(k,150) + (rxt(k,463)*y(k,126) + &
                 rxt(k,464)*y(k,217) +.500_r8*rxt(k,478)*y(k,217))*y(k,67) &
                  +rxt(k,466)*y(k,217)*y(k,137) +rxt(k,113)*y(k,152)
         loss(k,91) = (rxt(k,479)* y(k,229) + rxt(k,113) + het_rates(k,152))* y(k,152)
         prod(k,91) =rxt(k,109)*y(k,80) +rxt(k,468)*y(k,217)*y(k,151)
         loss(k,24) = ( + het_rates(k,153))* y(k,153)
         prod(k,24) = 0._r8
         loss(k,25) = ( + het_rates(k,154))* y(k,154)
         prod(k,25) = 0._r8
         loss(k,26) = ( + het_rates(k,155))* y(k,155)
         prod(k,26) = 0._r8
         loss(k,27) = ( + rxt(k,114) + het_rates(k,156))* y(k,156)
         prod(k,27) = 0._r8
         loss(k,28) = ( + rxt(k,115) + het_rates(k,157))* y(k,157)
         prod(k,28) = 0._r8
         loss(k,29) = ( + rxt(k,116) + het_rates(k,158))* y(k,158)
         prod(k,29) = 0._r8
         loss(k,30) = ( + rxt(k,117) + het_rates(k,159))* y(k,159)
         prod(k,30) = 0._r8
         loss(k,31) = ( + rxt(k,118) + het_rates(k,160))* y(k,160)
         prod(k,31) = 0._r8
         loss(k,32) = ( + rxt(k,119) + het_rates(k,161))* y(k,161)
         prod(k,32) = 0._r8
         loss(k,33) = ( + rxt(k,120) + het_rates(k,162))* y(k,162)
         prod(k,33) = 0._r8
         loss(k,34) = ( + rxt(k,121) + het_rates(k,163))* y(k,163)
         prod(k,34) = 0._r8
         loss(k,35) = ( + rxt(k,122) + het_rates(k,164))* y(k,164)
         prod(k,35) = 0._r8
         loss(k,36) = ( + rxt(k,123) + het_rates(k,165))* y(k,165)
         prod(k,36) = 0._r8
         loss(k,37) = ( + het_rates(k,166))* y(k,166)
         prod(k,37) = (.1279005_r8*rxt(k,496)*y(k,190) + &
                 .0097005_r8*rxt(k,501)*y(k,192) +.0003005_r8*rxt(k,504)*y(k,207) + &
                 .1056005_r8*rxt(k,508)*y(k,208) +.0245005_r8*rxt(k,512)*y(k,214) + &
                 .0154005_r8*rxt(k,518)*y(k,224) +.0063005_r8*rxt(k,522)*y(k,227)) &
                 *y(k,124) + (.2202005_r8*rxt(k,495)*y(k,190) + &
                 .0023005_r8*rxt(k,500)*y(k,192) +.0031005_r8*rxt(k,503)*y(k,207) + &
                 .2381005_r8*rxt(k,507)*y(k,208) +.0508005_r8*rxt(k,511)*y(k,214) + &
                 .1364005_r8*rxt(k,517)*y(k,224) +.1677005_r8*rxt(k,521)*y(k,227)) &
                 *y(k,203) + (.2202005_r8*rxt(k,497)*y(k,6) + &
                 .0508005_r8*rxt(k,513)*y(k,110))*y(k,134) +rxt(k,519)*y(k,75) &
                  +.5931005_r8*rxt(k,515)*y(k,217)*y(k,172)
         loss(k,38) = ( + het_rates(k,167))* y(k,167)
         prod(k,38) = (.1792005_r8*rxt(k,496)*y(k,190) + &
                 .0034005_r8*rxt(k,501)*y(k,192) +.0003005_r8*rxt(k,504)*y(k,207) + &
                 .1026005_r8*rxt(k,508)*y(k,208) +.0082005_r8*rxt(k,512)*y(k,214) + &
                 .0452005_r8*rxt(k,518)*y(k,224) +.0237005_r8*rxt(k,522)*y(k,227)) &
                 *y(k,124) + (.2067005_r8*rxt(k,495)*y(k,190) + &
                 .0008005_r8*rxt(k,500)*y(k,192) +.0035005_r8*rxt(k,503)*y(k,207) + &
                 .1308005_r8*rxt(k,507)*y(k,208) +.1149005_r8*rxt(k,511)*y(k,214) + &
                 .0101005_r8*rxt(k,517)*y(k,224) +.0174005_r8*rxt(k,521)*y(k,227)) &
                 *y(k,203) + (.2067005_r8*rxt(k,497)*y(k,6) + &
                 .1149005_r8*rxt(k,513)*y(k,110))*y(k,134) &
                  +.1534005_r8*rxt(k,515)*y(k,217)*y(k,172)
         loss(k,39) = ( + het_rates(k,168))* y(k,168)
         prod(k,39) = (.0676005_r8*rxt(k,496)*y(k,190) + &
                 .1579005_r8*rxt(k,501)*y(k,192) +.0073005_r8*rxt(k,504)*y(k,207) + &
                 .0521005_r8*rxt(k,508)*y(k,208) +.0772005_r8*rxt(k,512)*y(k,214) + &
                 .0966005_r8*rxt(k,518)*y(k,224) +.0025005_r8*rxt(k,522)*y(k,227)) &
                 *y(k,124) + (.0653005_r8*rxt(k,495)*y(k,190) + &
                 .0843005_r8*rxt(k,500)*y(k,192) +.0003005_r8*rxt(k,503)*y(k,207) + &
                 .0348005_r8*rxt(k,507)*y(k,208) +.0348005_r8*rxt(k,511)*y(k,214) + &
                 .0763005_r8*rxt(k,517)*y(k,224) +.086_r8*rxt(k,521)*y(k,227)) &
                 *y(k,203) + (.0653005_r8*rxt(k,497)*y(k,6) + &
                 .0348005_r8*rxt(k,513)*y(k,110))*y(k,134) &
                  +.0459005_r8*rxt(k,515)*y(k,217)*y(k,172)
         loss(k,40) = ( + het_rates(k,169))* y(k,169)
         prod(k,40) = (.079_r8*rxt(k,496)*y(k,190) +.0059005_r8*rxt(k,501)*y(k,192) + &
                 .0057005_r8*rxt(k,504)*y(k,207) +.0143005_r8*rxt(k,508)*y(k,208) + &
                 .0332005_r8*rxt(k,512)*y(k,214) +.0073005_r8*rxt(k,518)*y(k,224) + &
                 .011_r8*rxt(k,522)*y(k,227))*y(k,124) &
                  + (.1284005_r8*rxt(k,495)*y(k,190) + &
                 .0443005_r8*rxt(k,500)*y(k,192) +.0271005_r8*rxt(k,503)*y(k,207) + &
                 .0076005_r8*rxt(k,507)*y(k,208) +.0554005_r8*rxt(k,511)*y(k,214) + &
                 .2157005_r8*rxt(k,517)*y(k,224) +.0512005_r8*rxt(k,521)*y(k,227)) &
                 *y(k,203) + (.1749305_r8*rxt(k,494)*y(k,6) + &
                 .0590245_r8*rxt(k,502)*y(k,98) +.1749305_r8*rxt(k,510)*y(k,110)) &
                 *y(k,126) + (.1284005_r8*rxt(k,497)*y(k,6) + &
                 .0033005_r8*rxt(k,505)*y(k,98) +.0554005_r8*rxt(k,513)*y(k,110)) &
                 *y(k,134) +.0085005_r8*rxt(k,515)*y(k,217)*y(k,172)
         loss(k,41) = ( + het_rates(k,170))* y(k,170)
         prod(k,41) = (.1254005_r8*rxt(k,496)*y(k,190) + &
                 .0536005_r8*rxt(k,501)*y(k,192) +.0623005_r8*rxt(k,504)*y(k,207) + &
                 .0166005_r8*rxt(k,508)*y(k,208) +.130_r8*rxt(k,512)*y(k,214) + &
                 .238_r8*rxt(k,518)*y(k,224) +.1185005_r8*rxt(k,522)*y(k,227)) &
                 *y(k,124) + (.114_r8*rxt(k,495)*y(k,190) + &
                 .1621005_r8*rxt(k,500)*y(k,192) +.0474005_r8*rxt(k,503)*y(k,207) + &
                 .0113005_r8*rxt(k,507)*y(k,208) +.1278005_r8*rxt(k,511)*y(k,214) + &
                 .0738005_r8*rxt(k,517)*y(k,224) +.1598005_r8*rxt(k,521)*y(k,227)) &
                 *y(k,203) + (.5901905_r8*rxt(k,494)*y(k,6) + &
                 .0250245_r8*rxt(k,502)*y(k,98) +.5901905_r8*rxt(k,510)*y(k,110)) &
                 *y(k,126) + (.114_r8*rxt(k,497)*y(k,6) + &
                 .1278005_r8*rxt(k,513)*y(k,110))*y(k,134) &
                  +.0128005_r8*rxt(k,515)*y(k,217)*y(k,172)
         loss(k,42) = ( + rxt(k,542) + het_rates(k,171))* y(k,171)
         prod(k,42) = 0._r8
         loss(k,43) = (rxt(k,515)* y(k,217) + het_rates(k,172))* y(k,172)
         prod(k,43) = 0._r8
         loss(k,78) = ( + rxt(k,64) + het_rates(k,173))* y(k,173)
         prod(k,78) = (.100_r8*rxt(k,430)*y(k,180) +.230_r8*rxt(k,432)*y(k,182)) &
                 *y(k,217)
         loss(k,152) = (rxt(k,454)* y(k,217) + rxt(k,65) + het_rates(k,174))* y(k,174)
         prod(k,152) =rxt(k,452)*y(k,221)*y(k,203)
         loss(k,149) = (rxt(k,455)* y(k,217) + rxt(k,66) + rxt(k,493) &
                  + het_rates(k,175))* y(k,175)
         prod(k,149) = (.200_r8*rxt(k,448)*y(k,215) +.200_r8*rxt(k,458)*y(k,222)) &
                 *y(k,124) +.500_r8*rxt(k,446)*y(k,215)*y(k,198)
         loss(k,130) = (rxt(k,459)* y(k,217) + rxt(k,67) + het_rates(k,176))* y(k,176)
         prod(k,130) =rxt(k,457)*y(k,222)*y(k,203)
         loss(k,183) = (rxt(k,460)* y(k,126) +rxt(k,461)* y(k,217) + rxt(k,68) &
                  + het_rates(k,177))* y(k,177)
         prod(k,183) = (.500_r8*rxt(k,446)*y(k,198) +.800_r8*rxt(k,448)*y(k,124) + &
                 rxt(k,449)*y(k,126))*y(k,215) + (.330_r8*rxt(k,441)*y(k,6) + &
                 .330_r8*rxt(k,444)*y(k,110))*y(k,134) + (rxt(k,66) + &
                 rxt(k,455)*y(k,217))*y(k,175) + (rxt(k,456)*y(k,198) + &
                 .800_r8*rxt(k,458)*y(k,124))*y(k,222) +rxt(k,58)*y(k,128) +rxt(k,67) &
                 *y(k,176)
         loss(k,188) = (rxt(k,462)* y(k,217) + rxt(k,69) + het_rates(k,178))* y(k,178)
         prod(k,188) = (.300_r8*rxt(k,441)*y(k,6) +.300_r8*rxt(k,444)*y(k,110)) &
                 *y(k,134) + (rxt(k,451)*y(k,198) +.900_r8*rxt(k,453)*y(k,124)) &
                 *y(k,221) +rxt(k,65)*y(k,174) +rxt(k,68)*y(k,177)
         loss(k,153) = (rxt(k,429)* y(k,217) + rxt(k,70) + het_rates(k,179))* y(k,179)
         prod(k,153) =rxt(k,427)*y(k,223)*y(k,203)
         loss(k,76) = (rxt(k,430)* y(k,217) + het_rates(k,180))* y(k,180)
         prod(k,76) = 0._r8
         loss(k,79) = (rxt(k,396)* y(k,217) + rxt(k,71) + het_rates(k,181))* y(k,181)
         prod(k,79) =rxt(k,393)*y(k,225)*y(k,203)
         loss(k,80) = (rxt(k,432)* y(k,217) + het_rates(k,182))* y(k,182)
         prod(k,80) = 0._r8
         loss(k,158) = (rxt(k,435)* y(k,217) + rxt(k,72) + het_rates(k,183))* y(k,183)
         prod(k,158) =rxt(k,433)*y(k,226)*y(k,203)
         loss(k,81) = (rxt(k,438)* y(k,217) + het_rates(k,184))* y(k,184)
         prod(k,81) =.150_r8*rxt(k,432)*y(k,217)*y(k,182)
         loss(k,122) = (rxt(k,439)* y(k,217) + rxt(k,73) + het_rates(k,185))* y(k,185)
         prod(k,122) =rxt(k,436)*y(k,228)*y(k,203)
         loss(k,136) = (rxt(k,398)* y(k,124) +rxt(k,426)* y(k,125) +rxt(k,397) &
                 * y(k,203) + het_rates(k,188))* y(k,188)
         prod(k,136) =rxt(k,403)*y(k,217)*y(k,22) +rxt(k,431)*y(k,140)
         loss(k,178) = ((rxt(k,359) +rxt(k,360))* y(k,124) +rxt(k,358)* y(k,203) &
                  + het_rates(k,189))* y(k,189)
         prod(k,178) = (rxt(k,361)*y(k,2) +rxt(k,362)*y(k,15))*y(k,217)
         loss(k,44) = (rxt(k,496)* y(k,124) +rxt(k,495)* y(k,203) + het_rates(k,190)) &
                 * y(k,190)
         prod(k,44) =rxt(k,498)*y(k,217)*y(k,6)
         loss(k,131) = (rxt(k,401)* y(k,124) +rxt(k,400)* y(k,203) + het_rates(k,191)) &
                 * y(k,191)
         prod(k,131) = (.350_r8*rxt(k,399)*y(k,7) +rxt(k,402)*y(k,8))*y(k,217)
         loss(k,45) = (rxt(k,501)* y(k,124) +rxt(k,500)* y(k,203) + het_rates(k,192)) &
                 * y(k,192)
         prod(k,45) =rxt(k,499)*y(k,217)*y(k,7)
         loss(k,123) = (rxt(k,406)* y(k,124) +rxt(k,404)* y(k,203) + het_rates(k,193)) &
                 * y(k,193)
         prod(k,123) = (rxt(k,405)*y(k,23) +.070_r8*rxt(k,430)*y(k,180) + &
                 .060_r8*rxt(k,432)*y(k,182))*y(k,217)
         loss(k,170) = (rxt(k,284)* y(k,124) + 2._r8*rxt(k,281)* y(k,194) +rxt(k,282) &
                 * y(k,198) +rxt(k,283)* y(k,203) + het_rates(k,194))* y(k,194)
         prod(k,170) = (rxt(k,287)*y(k,56) +rxt(k,288)*y(k,217))*y(k,28) &
                  +.500_r8*rxt(k,286)*y(k,217)*y(k,27) +rxt(k,52)*y(k,107)
         loss(k,167) = (rxt(k,312)* y(k,124) +rxt(k,310)* y(k,198) +rxt(k,311) &
                 * y(k,203) + het_rates(k,195))* y(k,195)
         prod(k,167) = (rxt(k,313)*y(k,30) +rxt(k,314)*y(k,31))*y(k,217)
         loss(k,150) = (rxt(k,408)* y(k,124) +rxt(k,407)* y(k,203) + het_rates(k,196)) &
                 * y(k,196)
         prod(k,150) = (.400_r8*rxt(k,397)*y(k,203) +rxt(k,398)*y(k,124))*y(k,188) &
                  +rxt(k,409)*y(k,217)*y(k,32) +rxt(k,424)*y(k,141)*y(k,134)
         loss(k,207) = (rxt(k,380)* y(k,101) +rxt(k,295)* y(k,124) +rxt(k,306) &
                 * y(k,125) + 2._r8*rxt(k,292)* y(k,197) +rxt(k,293)* y(k,198) &
                  +rxt(k,294)* y(k,203) +rxt(k,366)* y(k,205) +rxt(k,371)* y(k,206) &
                  +rxt(k,333)* y(k,209) +rxt(k,391)* y(k,225) + het_rates(k,197)) &
                 * y(k,197)
         prod(k,207) = (.100_r8*rxt(k,339)*y(k,105) +.280_r8*rxt(k,353)*y(k,111) + &
                 .080_r8*rxt(k,386)*y(k,98) +.060_r8*rxt(k,441)*y(k,6) + &
                 .060_r8*rxt(k,444)*y(k,110))*y(k,134) + (rxt(k,343)*y(k,198) + &
                 .450_r8*rxt(k,344)*y(k,203) +2.000_r8*rxt(k,345)*y(k,211) + &
                 rxt(k,346)*y(k,124) +rxt(k,347)*y(k,126))*y(k,211) &
                  + (.530_r8*rxt(k,333)*y(k,197) +.260_r8*rxt(k,334)*y(k,198) + &
                 .530_r8*rxt(k,336)*y(k,126) +.530_r8*rxt(k,337)*y(k,124))*y(k,209) &
                  + (rxt(k,290)*y(k,45) +.500_r8*rxt(k,297)*y(k,51) + &
                 rxt(k,316)*y(k,49) +.650_r8*rxt(k,462)*y(k,178))*y(k,217) &
                  + (.300_r8*rxt(k,322)*y(k,198) +.150_r8*rxt(k,323)*y(k,203) + &
                 rxt(k,324)*y(k,124))*y(k,220) + (rxt(k,36) +rxt(k,315)*y(k,126)) &
                 *y(k,49) + (.600_r8*rxt(k,60) +rxt(k,307))*y(k,139) &
                  + (.200_r8*rxt(k,348)*y(k,203) +rxt(k,349)*y(k,124))*y(k,213) &
                  +.130_r8*rxt(k,23)*y(k,10) +rxt(k,27)*y(k,14) +rxt(k,289)*y(k,126) &
                 *y(k,45) +rxt(k,35)*y(k,48) +.330_r8*rxt(k,45)*y(k,93) +rxt(k,47) &
                 *y(k,95) +1.340_r8*rxt(k,50)*y(k,105) +rxt(k,52)*y(k,107) +rxt(k,53) &
                 *y(k,108) +.300_r8*rxt(k,55)*y(k,111) +rxt(k,57)*y(k,127) +rxt(k,63) &
                 *y(k,147) +.500_r8*rxt(k,64)*y(k,173) +.650_r8*rxt(k,69)*y(k,178)
         loss(k,221) = (rxt(k,183)* y(k,59) +rxt(k,381)* y(k,101) +rxt(k,263) &
                 * y(k,124) +rxt(k,282)* y(k,194) +rxt(k,310)* y(k,195) +rxt(k,293) &
                 * y(k,197) + 2._r8*(rxt(k,260) +rxt(k,261))* y(k,198) +rxt(k,262) &
                 * y(k,203) +rxt(k,367)* y(k,205) +rxt(k,372)* y(k,206) +rxt(k,334) &
                 * y(k,209) +rxt(k,343)* y(k,211) +rxt(k,446)* y(k,215) +rxt(k,322) &
                 * y(k,220) +rxt(k,451)* y(k,221) +rxt(k,456)* y(k,222) +rxt(k,392) &
                 * y(k,225) + het_rates(k,198))* y(k,198)
         prod(k,221) = (2.000_r8*rxt(k,292)*y(k,197) +.900_r8*rxt(k,293)*y(k,198) + &
                 .450_r8*rxt(k,294)*y(k,203) +rxt(k,295)*y(k,124) + &
                 rxt(k,333)*y(k,209) +rxt(k,342)*y(k,211) +rxt(k,366)*y(k,205) + &
                 rxt(k,371)*y(k,206) +rxt(k,380)*y(k,101) +rxt(k,391)*y(k,225)) &
                 *y(k,197) + (rxt(k,40) +rxt(k,177)*y(k,56) +rxt(k,233)*y(k,73) + &
                 rxt(k,266)*y(k,217) +rxt(k,272)*y(k,216))*y(k,54) &
                  + (.830_r8*rxt(k,412)*y(k,199) +.170_r8*rxt(k,418)*y(k,212)) &
                 *y(k,124) + (.280_r8*rxt(k,309)*y(k,29) +.050_r8*rxt(k,386)*y(k,98)) &
                 *y(k,134) + (.330_r8*rxt(k,411)*y(k,199) + &
                 .070_r8*rxt(k,417)*y(k,212))*y(k,203) + (.700_r8*rxt(k,265)*y(k,53) + &
                 rxt(k,296)*y(k,50))*y(k,217) +rxt(k,87)*y(k,43) +rxt(k,34)*y(k,45) &
                  +rxt(k,89)*y(k,46) +rxt(k,35)*y(k,48) +rxt(k,37)*y(k,51) &
                  +.300_r8*rxt(k,55)*y(k,111) +.400_r8*rxt(k,60)*y(k,139)
         loss(k,163) = (rxt(k,412)* y(k,124) +rxt(k,413)* y(k,125) +rxt(k,411) &
                 * y(k,203) + het_rates(k,199))* y(k,199)
         prod(k,163) =.600_r8*rxt(k,25)*y(k,12)
         loss(k,142) = ((rxt(k,330) +rxt(k,331))* y(k,124) + het_rates(k,200)) &
                 * y(k,200)
         prod(k,142) =rxt(k,329)*y(k,217)*y(k,16)
         loss(k,94) = ( + rxt(k,300) + rxt(k,301) + het_rates(k,201))* y(k,201)
         prod(k,94) =rxt(k,42)*y(k,72) +.750_r8*rxt(k,299)*y(k,202)*y(k,124)
         loss(k,159) = (rxt(k,299)* y(k,124) +rxt(k,298)* y(k,203) + het_rates(k,202)) &
                 * y(k,202)
         prod(k,159) =rxt(k,305)*y(k,217)*y(k,25)
         loss(k,218) = (rxt(k,213)* y(k,17) +rxt(k,219)* y(k,19) +rxt(k,256)* y(k,42) &
                  + (rxt(k,180) +rxt(k,181))* y(k,56) +rxt(k,187)* y(k,59) &
                  + (rxt(k,136) +rxt(k,137) +rxt(k,138))* y(k,76) +rxt(k,382) &
                 * y(k,101) +rxt(k,165)* y(k,124) +rxt(k,170)* y(k,125) +rxt(k,160) &
                 * y(k,126) +rxt(k,140)* y(k,133) +rxt(k,141)* y(k,134) +rxt(k,397) &
                 * y(k,188) +rxt(k,358)* y(k,189) +rxt(k,400)* y(k,191) +rxt(k,404) &
                 * y(k,193) +rxt(k,283)* y(k,194) +rxt(k,311)* y(k,195) +rxt(k,407) &
                 * y(k,196) +rxt(k,294)* y(k,197) +rxt(k,262)* y(k,198) +rxt(k,411) &
                 * y(k,199) +rxt(k,298)* y(k,202) + 2._r8*rxt(k,150)* y(k,203) &
                  +rxt(k,269)* y(k,204) +rxt(k,368)* y(k,205) +rxt(k,373)* y(k,206) &
                  +rxt(k,335)* y(k,209) +rxt(k,414)* y(k,210) +rxt(k,344)* y(k,211) &
                  +rxt(k,417)* y(k,212) +rxt(k,348)* y(k,213) +rxt(k,447)* y(k,215) &
                  +rxt(k,145)* y(k,217) +rxt(k,420)* y(k,218) +rxt(k,319)* y(k,219) &
                  +rxt(k,323)* y(k,220) +rxt(k,452)* y(k,221) +rxt(k,457)* y(k,222) &
                  +rxt(k,427)* y(k,223) +rxt(k,393)* y(k,225) +rxt(k,433)* y(k,226) &
                  +rxt(k,436)* y(k,228) + rxt(k,481) + het_rates(k,203))* y(k,203)
         prod(k,218) = (rxt(k,144)*y(k,79) +rxt(k,147)*y(k,134) +rxt(k,163)*y(k,126) + &
                 rxt(k,194)*y(k,59) +rxt(k,224)*y(k,19) +rxt(k,242)*y(k,43) + &
                 rxt(k,245)*y(k,46) +rxt(k,264)*y(k,52) +rxt(k,267)*y(k,86) + &
                 rxt(k,268)*y(k,87) +rxt(k,276)*y(k,62) +.350_r8*rxt(k,278)*y(k,24) + &
                 rxt(k,285)*y(k,26) +rxt(k,291)*y(k,47) +rxt(k,302)*y(k,74) + &
                 rxt(k,303)*y(k,75) +rxt(k,317)*y(k,95) +rxt(k,332)*y(k,93) + &
                 .200_r8*rxt(k,341)*y(k,106) +.500_r8*rxt(k,352)*y(k,109) + &
                 .300_r8*rxt(k,377)*y(k,99) +rxt(k,378)*y(k,100) + &
                 rxt(k,385)*y(k,102) +rxt(k,389)*y(k,115) +rxt(k,390)*y(k,116) + &
                 .650_r8*rxt(k,399)*y(k,7) +.730_r8*rxt(k,410)*y(k,66) + &
                 .800_r8*rxt(k,422)*y(k,142) +.280_r8*rxt(k,430)*y(k,180) + &
                 .380_r8*rxt(k,432)*y(k,182) +.630_r8*rxt(k,438)*y(k,184) + &
                 .200_r8*rxt(k,462)*y(k,178) +rxt(k,468)*y(k,151) + &
                 .500_r8*rxt(k,478)*y(k,67))*y(k,217) + (rxt(k,263)*y(k,198) + &
                 rxt(k,271)*y(k,204) +rxt(k,284)*y(k,194) + &
                 .250_r8*rxt(k,299)*y(k,202) +rxt(k,312)*y(k,195) + &
                 rxt(k,320)*y(k,219) +rxt(k,330)*y(k,200) + &
                 .470_r8*rxt(k,337)*y(k,209) +rxt(k,359)*y(k,189) + &
                 .920_r8*rxt(k,369)*y(k,205) +.920_r8*rxt(k,375)*y(k,206) + &
                 rxt(k,383)*y(k,101) +rxt(k,394)*y(k,225) +rxt(k,401)*y(k,191) + &
                 rxt(k,406)*y(k,193) +.170_r8*rxt(k,412)*y(k,199) + &
                 .400_r8*rxt(k,415)*y(k,210) +.830_r8*rxt(k,418)*y(k,212) + &
                 rxt(k,421)*y(k,218) +rxt(k,428)*y(k,223) +rxt(k,434)*y(k,226) + &
                 rxt(k,437)*y(k,228) +.900_r8*rxt(k,453)*y(k,221) + &
                 .800_r8*rxt(k,458)*y(k,222))*y(k,124) + (rxt(k,183)*y(k,59) + &
                 2.000_r8*rxt(k,260)*y(k,198) +rxt(k,282)*y(k,194) + &
                 .900_r8*rxt(k,293)*y(k,197) +rxt(k,310)*y(k,195) + &
                 .300_r8*rxt(k,322)*y(k,220) +.730_r8*rxt(k,334)*y(k,209) + &
                 rxt(k,343)*y(k,211) +rxt(k,367)*y(k,205) +rxt(k,372)*y(k,206) + &
                 1.200_r8*rxt(k,381)*y(k,101) +.800_r8*rxt(k,392)*y(k,225) + &
                 .500_r8*rxt(k,446)*y(k,215) +rxt(k,451)*y(k,221) + &
                 rxt(k,456)*y(k,222))*y(k,198) + (.130_r8*rxt(k,280)*y(k,25) + &
                 .280_r8*rxt(k,309)*y(k,29) +.140_r8*rxt(k,339)*y(k,105) + &
                 .280_r8*rxt(k,353)*y(k,111) +.370_r8*rxt(k,386)*y(k,98) + &
                 .570_r8*rxt(k,441)*y(k,6) +.570_r8*rxt(k,444)*y(k,110))*y(k,134) &
                  + (rxt(k,257)*y(k,42) +.470_r8*rxt(k,336)*y(k,209) + &
                 rxt(k,370)*y(k,205) +rxt(k,376)*y(k,206) +rxt(k,384)*y(k,101) + &
                 rxt(k,395)*y(k,225))*y(k,126) + (.470_r8*rxt(k,333)*y(k,209) + &
                 rxt(k,366)*y(k,205) +rxt(k,371)*y(k,206) +rxt(k,380)*y(k,101) + &
                 rxt(k,391)*y(k,225))*y(k,197) + (rxt(k,176)*y(k,42) + &
                 rxt(k,179)*y(k,79) +rxt(k,241)*y(k,43) +rxt(k,244)*y(k,46))*y(k,56) &
                  + (.070_r8*rxt(k,411)*y(k,199) +.160_r8*rxt(k,414)*y(k,210) + &
                 .330_r8*rxt(k,417)*y(k,212))*y(k,203) + (rxt(k,212)*y(k,17) + &
                 rxt(k,258)*y(k,133))*y(k,42) + (rxt(k,11) +rxt(k,174))*y(k,90) &
                  + (1.340_r8*rxt(k,50) +.660_r8*rxt(k,51))*y(k,105) + (rxt(k,300) + &
                 rxt(k,301))*y(k,201) +rxt(k,19)*y(k,1) +.900_r8*rxt(k,20)*y(k,2) &
                  +rxt(k,21)*y(k,8) +1.500_r8*rxt(k,22)*y(k,9) +.560_r8*rxt(k,23) &
                 *y(k,10) +rxt(k,24)*y(k,11) +.600_r8*rxt(k,25)*y(k,12) &
                  +.600_r8*rxt(k,26)*y(k,13) +rxt(k,27)*y(k,14) +rxt(k,28)*y(k,23) &
                  +rxt(k,29)*y(k,27) +rxt(k,30)*y(k,30) +rxt(k,34)*y(k,45) +rxt(k,36) &
                 *y(k,49) +rxt(k,273)*y(k,216)*y(k,54) +2.000_r8*rxt(k,43)*y(k,74) &
                  +2.000_r8*rxt(k,44)*y(k,75) +rxt(k,139)*y(k,76) +rxt(k,135)*y(k,133) &
                 *y(k,79) +.670_r8*rxt(k,45)*y(k,93) +rxt(k,46)*y(k,94) +rxt(k,47) &
                 *y(k,95) +rxt(k,48)*y(k,102) +rxt(k,49)*y(k,103) +rxt(k,56)*y(k,116) &
                  +rxt(k,61)*y(k,143) +rxt(k,62)*y(k,146) +rxt(k,64)*y(k,173) &
                  +rxt(k,65)*y(k,174) +rxt(k,66)*y(k,175) +rxt(k,67)*y(k,176) &
                  +rxt(k,68)*y(k,177) +1.200_r8*rxt(k,69)*y(k,178) +rxt(k,70)*y(k,179) &
                  +rxt(k,72)*y(k,183) +rxt(k,73)*y(k,185) &
                  +1.200_r8*rxt(k,281)*y(k,194)*y(k,194) +rxt(k,270)*y(k,204) &
                  +rxt(k,374)*y(k,206)
         loss(k,124) = (rxt(k,271)* y(k,124) +rxt(k,269)* y(k,203) + rxt(k,270) &
                  + het_rates(k,204))* y(k,204)
         prod(k,124) =rxt(k,256)*y(k,203)*y(k,42)
         loss(k,202) = (rxt(k,369)* y(k,124) +rxt(k,370)* y(k,126) +rxt(k,366) &
                 * y(k,197) +rxt(k,367)* y(k,198) +rxt(k,368)* y(k,203) &
                  + het_rates(k,205))* y(k,205)
         prod(k,202) =.600_r8*rxt(k,387)*y(k,217)*y(k,98)
         loss(k,203) = (rxt(k,375)* y(k,124) +rxt(k,376)* y(k,126) +rxt(k,371) &
                 * y(k,197) +rxt(k,372)* y(k,198) +rxt(k,373)* y(k,203) + rxt(k,374) &
                  + het_rates(k,206))* y(k,206)
         prod(k,203) =.400_r8*rxt(k,387)*y(k,217)*y(k,98)
         loss(k,46) = (rxt(k,504)* y(k,124) +rxt(k,503)* y(k,203) + het_rates(k,207)) &
                 * y(k,207)
         prod(k,46) =rxt(k,506)*y(k,217)*y(k,98)
         loss(k,47) = (rxt(k,508)* y(k,124) +rxt(k,507)* y(k,203) + het_rates(k,208)) &
                 * y(k,208)
         prod(k,47) =rxt(k,509)*y(k,217)*y(k,104)
         loss(k,204) = ((rxt(k,337) +rxt(k,338))* y(k,124) +rxt(k,336)* y(k,126) &
                  +rxt(k,333)* y(k,197) +rxt(k,334)* y(k,198) +rxt(k,335)* y(k,203) &
                  + het_rates(k,209))* y(k,209)
         prod(k,204) = (.500_r8*rxt(k,340)*y(k,105) +.200_r8*rxt(k,341)*y(k,106) + &
                 rxt(k,354)*y(k,111))*y(k,217)
         loss(k,160) = (rxt(k,415)* y(k,124) +rxt(k,416)* y(k,125) +rxt(k,414) &
                 * y(k,203) + het_rates(k,210))* y(k,210)
         prod(k,160) =.600_r8*rxt(k,24)*y(k,11)
         loss(k,206) = (rxt(k,346)* y(k,124) +rxt(k,355)* y(k,125) +rxt(k,347) &
                 * y(k,126) +rxt(k,342)* y(k,197) +rxt(k,343)* y(k,198) +rxt(k,344) &
                 * y(k,203) + 2._r8*rxt(k,345)* y(k,211) + het_rates(k,211))* y(k,211)
         prod(k,206) = (.660_r8*rxt(k,50) +.500_r8*rxt(k,340)*y(k,217))*y(k,105) &
                  + (rxt(k,54) +rxt(k,356))*y(k,109) +.500_r8*rxt(k,341)*y(k,217) &
                 *y(k,106)
         loss(k,176) = (rxt(k,418)* y(k,124) +rxt(k,419)* y(k,125) +rxt(k,417) &
                 * y(k,203) + het_rates(k,212))* y(k,212)
         prod(k,176) =.600_r8*rxt(k,26)*y(k,13)
         loss(k,155) = (rxt(k,349)* y(k,124) +rxt(k,348)* y(k,203) + het_rates(k,213)) &
                 * y(k,213)
         prod(k,155) = (rxt(k,350)*y(k,107) +rxt(k,351)*y(k,108))*y(k,217)
         loss(k,49) = (rxt(k,512)* y(k,124) +rxt(k,511)* y(k,203) + het_rates(k,214)) &
                 * y(k,214)
         prod(k,49) =rxt(k,514)*y(k,217)*y(k,110)
         loss(k,187) = (rxt(k,448)* y(k,124) +rxt(k,449)* y(k,126) +rxt(k,446) &
                 * y(k,198) +rxt(k,447)* y(k,203) + het_rates(k,215))* y(k,215)
         prod(k,187) = (rxt(k,440)*y(k,6) +rxt(k,443)*y(k,110) + &
                 .500_r8*rxt(k,460)*y(k,177))*y(k,126) +rxt(k,450)*y(k,217)*y(k,128)
         loss(k,214) = (rxt(k,201)* y(k,33) +rxt(k,202)* y(k,34) +rxt(k,228)* y(k,35) &
                  +rxt(k,203)* y(k,36) +rxt(k,204)* y(k,37) +rxt(k,205)* y(k,38) &
                  +rxt(k,206)* y(k,39) +rxt(k,207)* y(k,40) +rxt(k,251)* y(k,41) &
                  +rxt(k,252)* y(k,43) + (rxt(k,272) +rxt(k,273) +rxt(k,274))* y(k,54) &
                  +rxt(k,229)* y(k,55) +rxt(k,237)* y(k,64) +rxt(k,238)* y(k,65) &
                  +rxt(k,125)* y(k,77) +rxt(k,230)* y(k,78) + (rxt(k,231) +rxt(k,232)) &
                 * y(k,81) +rxt(k,253)* y(k,82) +rxt(k,254)* y(k,83) +rxt(k,255) &
                 * y(k,84) + (rxt(k,208) +rxt(k,209))* y(k,85) +rxt(k,275)* y(k,86) &
                  + (rxt(k,168) +rxt(k,169))* y(k,113) + (rxt(k,129) +rxt(k,130)) &
                 * y(k,134) +rxt(k,126)* y(k,229) + rxt(k,127) + rxt(k,128) &
                  + het_rates(k,216))* y(k,216)
         prod(k,214) =rxt(k,12)*y(k,113) +rxt(k,7)*y(k,134) +rxt(k,1)*y(k,229)
         loss(k,215) = (rxt(k,357)* y(k,1) +rxt(k,361)* y(k,2) +rxt(k,442)* y(k,6) &
                  +rxt(k,399)* y(k,7) +rxt(k,402)* y(k,8) +rxt(k,362)* y(k,15) &
                  +rxt(k,329)* y(k,16) +rxt(k,224)* y(k,19) +rxt(k,403)* y(k,22) &
                  +rxt(k,405)* y(k,23) +rxt(k,278)* y(k,24) +rxt(k,305)* y(k,25) &
                  +rxt(k,285)* y(k,26) +rxt(k,286)* y(k,27) +rxt(k,288)* y(k,28) &
                  +rxt(k,326)* y(k,29) +rxt(k,313)* y(k,30) +rxt(k,314)* y(k,31) &
                  +rxt(k,409)* y(k,32) +rxt(k,240)* y(k,41) +rxt(k,259)* y(k,42) &
                  +rxt(k,242)* y(k,43) +rxt(k,243)* y(k,44) +rxt(k,290)* y(k,45) &
                  +rxt(k,245)* y(k,46) +rxt(k,291)* y(k,47) +rxt(k,327)* y(k,48) &
                  +rxt(k,316)* y(k,49) +rxt(k,296)* y(k,50) +rxt(k,297)* y(k,51) &
                  +rxt(k,264)* y(k,52) +rxt(k,265)* y(k,53) +rxt(k,266)* y(k,54) &
                  +rxt(k,247)* y(k,55) + (rxt(k,194) +rxt(k,195))* y(k,59) +rxt(k,192) &
                 * y(k,60) +rxt(k,276)* y(k,62) +rxt(k,410)* y(k,66) + (rxt(k,464) + &
                 rxt(k,478))* y(k,67) +rxt(k,302)* y(k,74) +rxt(k,303)* y(k,75) &
                  +rxt(k,143)* y(k,77) +rxt(k,144)* y(k,79) +rxt(k,226)* y(k,81) &
                  +rxt(k,248)* y(k,82) +rxt(k,249)* y(k,83) +rxt(k,250)* y(k,84) &
                  +rxt(k,197)* y(k,85) +rxt(k,267)* y(k,86) +rxt(k,268)* y(k,87) &
                  +rxt(k,173)* y(k,89) +rxt(k,151)* y(k,90) +rxt(k,200)* y(k,92) &
                  +rxt(k,332)* y(k,93) +rxt(k,363)* y(k,94) +rxt(k,317)* y(k,95) &
                  +rxt(k,364)* y(k,96) +rxt(k,365)* y(k,97) +rxt(k,387)* y(k,98) &
                  +rxt(k,377)* y(k,99) +rxt(k,378)* y(k,100) +rxt(k,385)* y(k,102) &
                  +rxt(k,388)* y(k,103) +rxt(k,340)* y(k,105) +rxt(k,341)* y(k,106) &
                  +rxt(k,350)* y(k,107) +rxt(k,351)* y(k,108) +rxt(k,352)* y(k,109) &
                  +rxt(k,445)* y(k,110) +rxt(k,354)* y(k,111) +rxt(k,164)* y(k,112) &
                  +rxt(k,389)* y(k,115) +rxt(k,390)* y(k,116) +rxt(k,480)* y(k,120) &
                  +rxt(k,172)* y(k,125) +rxt(k,163)* y(k,126) +rxt(k,318)* y(k,127) &
                  +rxt(k,450)* y(k,128) +rxt(k,146)* y(k,133) +rxt(k,147)* y(k,134) &
                  +rxt(k,466)* y(k,137) +rxt(k,304)* y(k,139) +rxt(k,422)* y(k,142) &
                  +rxt(k,425)* y(k,143) +rxt(k,321)* y(k,146) +rxt(k,325)* y(k,147) &
                  +rxt(k,472)* y(k,148) +rxt(k,477)* y(k,150) +rxt(k,468)* y(k,151) &
                  +rxt(k,454)* y(k,174) +rxt(k,455)* y(k,175) +rxt(k,459)* y(k,176) &
                  +rxt(k,461)* y(k,177) +rxt(k,462)* y(k,178) +rxt(k,429)* y(k,179) &
                  +rxt(k,430)* y(k,180) +rxt(k,396)* y(k,181) +rxt(k,432)* y(k,182) &
                  +rxt(k,435)* y(k,183) +rxt(k,438)* y(k,184) +rxt(k,439)* y(k,185) &
                  +rxt(k,145)* y(k,203) + 2._r8*(rxt(k,148) +rxt(k,149))* y(k,217) &
                  + het_rates(k,217))* y(k,217)
         prod(k,215) = (2.000_r8*rxt(k,137)*y(k,76) +rxt(k,140)*y(k,133) + &
                 rxt(k,141)*y(k,134) +rxt(k,160)*y(k,126) +rxt(k,165)*y(k,124) + &
                 rxt(k,181)*y(k,56) +.450_r8*rxt(k,294)*y(k,197) + &
                 .150_r8*rxt(k,323)*y(k,220) +.450_r8*rxt(k,344)*y(k,211) + &
                 .200_r8*rxt(k,348)*y(k,213) +.400_r8*rxt(k,397)*y(k,188) + &
                 .400_r8*rxt(k,411)*y(k,199) +.400_r8*rxt(k,417)*y(k,212))*y(k,203) &
                  + (rxt(k,142)*y(k,76) +.130_r8*rxt(k,280)*y(k,25) + &
                 .360_r8*rxt(k,309)*y(k,29) +.240_r8*rxt(k,339)*y(k,105) + &
                 .360_r8*rxt(k,353)*y(k,111) +.320_r8*rxt(k,386)*y(k,98) + &
                 .630_r8*rxt(k,441)*y(k,6) +.630_r8*rxt(k,444)*y(k,110))*y(k,134) &
                  + (rxt(k,134)*y(k,77) +rxt(k,135)*y(k,79) +rxt(k,196)*y(k,85) + &
                 rxt(k,199)*y(k,92) +rxt(k,225)*y(k,81) +rxt(k,227)*y(k,91) + &
                 rxt(k,258)*y(k,42))*y(k,133) + (.300_r8*rxt(k,265)*y(k,53) + &
                 .650_r8*rxt(k,278)*y(k,24) +.500_r8*rxt(k,286)*y(k,27) + &
                 .500_r8*rxt(k,321)*y(k,146) +.100_r8*rxt(k,341)*y(k,106) + &
                 .600_r8*rxt(k,388)*y(k,103) +.500_r8*rxt(k,396)*y(k,181))*y(k,217) &
                  + (rxt(k,125)*y(k,77) +2.000_r8*rxt(k,126)*y(k,229) + &
                 rxt(k,208)*y(k,85) +rxt(k,231)*y(k,81) +rxt(k,272)*y(k,54) + &
                 rxt(k,275)*y(k,86))*y(k,216) + (rxt(k,2) +rxt(k,235)*y(k,73)) &
                 *y(k,229) +rxt(k,20)*y(k,2) +rxt(k,21)*y(k,8) +rxt(k,28)*y(k,23) &
                  +rxt(k,29)*y(k,27) +rxt(k,30)*y(k,30) +rxt(k,31)*y(k,32) +rxt(k,37) &
                 *y(k,51) +rxt(k,38)*y(k,53) +.330_r8*rxt(k,39)*y(k,54) +rxt(k,42) &
                 *y(k,72) +2.000_r8*rxt(k,4)*y(k,79) +rxt(k,9)*y(k,89) +rxt(k,10) &
                 *y(k,90) +rxt(k,105)*y(k,91) +rxt(k,106)*y(k,92) +rxt(k,46)*y(k,94) &
                  +rxt(k,49)*y(k,103) +rxt(k,53)*y(k,108) +.500_r8*rxt(k,489)*y(k,125) &
                  +rxt(k,58)*y(k,128) +rxt(k,61)*y(k,143) +rxt(k,62)*y(k,146) &
                  +rxt(k,63)*y(k,147) +rxt(k,65)*y(k,174) +rxt(k,67)*y(k,176) &
                  +rxt(k,70)*y(k,179) +rxt(k,71)*y(k,181) +rxt(k,72)*y(k,183) &
                  +rxt(k,73)*y(k,185)
         loss(k,126) = (rxt(k,421)* y(k,124) +rxt(k,420)* y(k,203) + het_rates(k,218)) &
                 * y(k,218)
         prod(k,126) = (.200_r8*rxt(k,410)*y(k,66) +.140_r8*rxt(k,422)*y(k,142) + &
                 rxt(k,425)*y(k,143))*y(k,217)
         loss(k,164) = (rxt(k,320)* y(k,124) +rxt(k,319)* y(k,203) + het_rates(k,219)) &
                 * y(k,219)
         prod(k,164) = (.500_r8*rxt(k,321)*y(k,146) +rxt(k,326)*y(k,29))*y(k,217)
         loss(k,194) = (rxt(k,324)* y(k,124) +rxt(k,322)* y(k,198) +rxt(k,323) &
                 * y(k,203) + het_rates(k,220))* y(k,220)
         prod(k,194) = (rxt(k,325)*y(k,147) +rxt(k,327)*y(k,48) + &
                 .150_r8*rxt(k,462)*y(k,178))*y(k,217) + (.060_r8*rxt(k,441)*y(k,6) + &
                 .060_r8*rxt(k,444)*y(k,110))*y(k,134) +.150_r8*rxt(k,69)*y(k,178)
         loss(k,193) = (rxt(k,453)* y(k,124) +rxt(k,451)* y(k,198) +rxt(k,452) &
                 * y(k,203) + het_rates(k,221))* y(k,221)
         prod(k,193) = (.500_r8*rxt(k,460)*y(k,126) +rxt(k,461)*y(k,217))*y(k,177) &
                  +rxt(k,454)*y(k,217)*y(k,174)
         loss(k,179) = (rxt(k,458)* y(k,124) +rxt(k,456)* y(k,198) +rxt(k,457) &
                 * y(k,203) + het_rates(k,222))* y(k,222)
         prod(k,179) = (rxt(k,442)*y(k,6) +rxt(k,445)*y(k,110) +rxt(k,459)*y(k,176)) &
                 *y(k,217)
         loss(k,161) = (rxt(k,428)* y(k,124) +rxt(k,427)* y(k,203) + het_rates(k,223)) &
                 * y(k,223)
         prod(k,161) = (rxt(k,429)*y(k,179) +.650_r8*rxt(k,430)*y(k,180))*y(k,217)
         loss(k,50) = (rxt(k,518)* y(k,124) +rxt(k,517)* y(k,203) + het_rates(k,224)) &
                 * y(k,224)
         prod(k,50) =rxt(k,516)*y(k,217)*y(k,180)
         loss(k,197) = (rxt(k,394)* y(k,124) +rxt(k,395)* y(k,126) +rxt(k,391) &
                 * y(k,197) +rxt(k,392)* y(k,198) +rxt(k,393)* y(k,203) &
                  + het_rates(k,225))* y(k,225)
         prod(k,197) = (rxt(k,363)*y(k,94) +rxt(k,364)*y(k,96) +rxt(k,365)*y(k,97) + &
                 .400_r8*rxt(k,388)*y(k,103) +.500_r8*rxt(k,396)*y(k,181))*y(k,217)
         loss(k,162) = (rxt(k,434)* y(k,124) +rxt(k,433)* y(k,203) + het_rates(k,226)) &
                 * y(k,226)
         prod(k,162) = (.560_r8*rxt(k,432)*y(k,182) +rxt(k,435)*y(k,183))*y(k,217)
         loss(k,51) = (rxt(k,522)* y(k,124) +rxt(k,521)* y(k,203) + het_rates(k,227)) &
                 * y(k,227)
         prod(k,51) =rxt(k,520)*y(k,217)*y(k,182)
         loss(k,133) = (rxt(k,437)* y(k,124) +rxt(k,436)* y(k,203) + het_rates(k,228)) &
                 * y(k,228)
         prod(k,133) = (.300_r8*rxt(k,438)*y(k,184) +rxt(k,439)*y(k,185))*y(k,217)
         loss(k,227) = (rxt(k,235)* y(k,73) +rxt(k,479)* y(k,152) +rxt(k,126) &
                 * y(k,216) + rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,229)) &
                 * y(k,229)
         prod(k,227) = (rxt(k,143)*y(k,77) +rxt(k,144)*y(k,79) +rxt(k,145)*y(k,203) + &
                 rxt(k,148)*y(k,217) +rxt(k,151)*y(k,90) +rxt(k,173)*y(k,89) + &
                 rxt(k,197)*y(k,85) +rxt(k,200)*y(k,92) +rxt(k,226)*y(k,81) + &
                 rxt(k,240)*y(k,41) +rxt(k,242)*y(k,43) +rxt(k,243)*y(k,44) + &
                 rxt(k,245)*y(k,46) +rxt(k,250)*y(k,84) +rxt(k,259)*y(k,42) + &
                 rxt(k,265)*y(k,53) +rxt(k,266)*y(k,54) +rxt(k,268)*y(k,87) + &
                 rxt(k,288)*y(k,28) +rxt(k,290)*y(k,45) +rxt(k,296)*y(k,50) + &
                 rxt(k,297)*y(k,51) +rxt(k,313)*y(k,30) +rxt(k,314)*y(k,31) + &
                 rxt(k,316)*y(k,49) +rxt(k,321)*y(k,146) +rxt(k,325)*y(k,147) + &
                 rxt(k,327)*y(k,48) +.500_r8*rxt(k,340)*y(k,105) +rxt(k,480)*y(k,120)) &
                 *y(k,217) + (rxt(k,524)*y(k,92) +rxt(k,530)*y(k,92) + &
                 rxt(k,531)*y(k,91) +rxt(k,535)*y(k,92) +rxt(k,536)*y(k,91))*y(k,85) &
                  + (rxt(k,481) +rxt(k,138)*y(k,76))*y(k,203) +.050_r8*rxt(k,39) &
                 *y(k,54) +rxt(k,109)*y(k,80)
      end do

      end subroutine imp_prod_loss

      end module mo_prod_loss

      module mo_indprd

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: indprd

      contains

      subroutine indprd( class, prod, nprod, y, extfrc, rxt, chnkpnts )

      use chem_mods, only : gas_pcnst, extcnt, rxntot

      implicit none

!--------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: chnkpnts
      integer, intent(in) :: nprod
      real(r8), intent(in)    :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in)    :: rxt(chnkpnts,rxntot)
      real(r8), intent(in)    :: extfrc(chnkpnts,extcnt)
      real(r8), intent(inout) :: prod(chnkpnts,nprod)

!--------------------------------------------------------------------
!       ... "independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then
         prod(:,1) =rxt(:,480)*y(:,217)*y(:,120) +rxt(:,488)*y(:,121)
                                                                                          
         prod(:,2) = (rxt(:,413)*y(:,199) +rxt(:,416)*y(:,210) +rxt(:,419)*y(:,212) + &
                 rxt(:,423)*y(:,141))*y(:,125) +.500_r8*rxt(:,352)*y(:,217)*y(:,109) &
                  +.200_r8*rxt(:,448)*y(:,215)*y(:,124) +.500_r8*rxt(:,460)*y(:,177) &
                 *y(:,126)
                                                                                          
!--------------------------------------------------------------------
!       ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,154) = 0._r8
                                                                                          
         prod(:,151) = 0._r8
                                                                                          
         prod(:,1) = + extfrc(:,12)
                                                                                          
         prod(:,2) = 0._r8
                                                                                          
         prod(:,3) = + extfrc(:,3)
                                                                                          
         prod(:,185) = 0._r8
                                                                                          
         prod(:,71) = 0._r8
                                                                                          
         prod(:,121) = 0._r8
                                                                                          
         prod(:,72) = 0._r8
                                                                                          
         prod(:,114) = 0._r8
                                                                                          
         prod(:,127) = 0._r8
                                                                                          
         prod(:,96) = 0._r8
                                                                                          
         prod(:,145) = 0._r8
                                                                                          
         prod(:,105) = 0._r8
                                                                                          
         prod(:,86) = 0._r8
                                                                                          
         prod(:,110) = 0._r8
                                                                                          
         prod(:,209) = 0._r8
                                                                                          
         prod(:,87) = 0._r8
                                                                                          
         prod(:,225) = 0._r8
                                                                                          
         prod(:,139) = 0._r8
                                                                                          
         prod(:,4) = 0._r8
                                                                                          
         prod(:,90) = 0._r8
                                                                                          
         prod(:,108) = 0._r8
                                                                                          
         prod(:,98) = 0._r8
                                                                                          
         prod(:,140) = 0._r8
                                                                                          
         prod(:,92) = 0._r8
                                                                                          
         prod(:,109) = 0._r8
                                                                                          
         prod(:,99) = 0._r8
                                                                                          
         prod(:,186) = 0._r8
                                                                                          
         prod(:,120) = 0._r8
                                                                                          
         prod(:,57) = 0._r8
                                                                                          
         prod(:,93) = 0._r8
                                                                                          
         prod(:,54) = 0._r8
                                                                                          
         prod(:,66) = 0._r8
                                                                                          
         prod(:,67) = 0._r8
                                                                                          
         prod(:,58) = 0._r8
                                                                                          
         prod(:,68) = 0._r8
                                                                                          
         prod(:,59) = 0._r8
                                                                                          
         prod(:,69) = 0._r8
                                                                                          
         prod(:,60) = 0._r8
                                                                                          
         prod(:,129) = 0._r8
                                                                                          
         prod(:,213) = 0._r8
                                                                                          
         prod(:,146) = 0._r8
                                                                                          
         prod(:,61) = 0._r8
                                                                                          
         prod(:,190) = 0._r8
                                                                                          
         prod(:,112) = 0._r8
                                                                                          
         prod(:,55) = 0._r8
                                                                                          
         prod(:,180) = 0._r8
                                                                                          
         prod(:,200) = 0._r8
                                                                                          
         prod(:,156) = 0._r8
                                                                                          
         prod(:,148) = 0._r8
                                                                                          
         prod(:,166) = 0._r8
                                                                                          
         prod(:,115) = 0._r8
                                                                                          
         prod(:,210) = 0._r8
                                                                                          
         prod(:,125) = 0._r8
                                                                                          
         prod(:,220) = 0._r8
                                                                                          
         prod(:,70) = 0._r8
                                                                                          
         prod(:,52) = 0._r8
                                                                                          
         prod(:,219) = 0._r8
                                                                                          
         prod(:,177) = 0._r8
                                                                                          
         prod(:,5) = 0._r8
                                                                                          
         prod(:,192) = + extfrc(:,7)
                                                                                          
         prod(:,168) = 0._r8
                                                                                          
         prod(:,85) = 0._r8
                                                                                          
         prod(:,83) = 0._r8
                                                                                          
         prod(:,77) = 0._r8
                                                                                          
         prod(:,100) = 0._r8
                                                                                          
         prod(:,6) = 0._r8
                                                                                          
         prod(:,7) = 0._r8
                                                                                          
         prod(:,8) = 0._r8
                                                                                          
         prod(:,9) = 0._r8
                                                                                          
         prod(:,62) = 0._r8
                                                                                          
         prod(:,175) = 0._r8
                                                                                          
         prod(:,191) = 0._r8
                                                                                          
         prod(:,184) = 0._r8
                                                                                          
         prod(:,212) = 0._r8
                                                                                          
         prod(:,208) = 0._r8
                                                                                          
         prod(:,56) = 0._r8
                                                                                          
         prod(:,147) = 0._r8
                                                                                          
         prod(:,63) = 0._r8
                                                                                          
         prod(:,169) = 0._r8
                                                                                          
         prod(:,82) = 0._r8
                                                                                          
         prod(:,89) = 0._r8
                                                                                          
         prod(:,101) = 0._r8
                                                                                          
         prod(:,223) = 0._r8
                                                                                          
         prod(:,74) = 0._r8
                                                                                          
         prod(:,181) = 0._r8
                                                                                          
         prod(:,97) = 0._r8
                                                                                          
         prod(:,211) = 0._r8
                                                                                          
         prod(:,118) = 0._r8
                                                                                          
         prod(:,165) = 0._r8
                                                                                          
         prod(:,171) = 0._r8
                                                                                          
         prod(:,196) = 0._r8
                                                                                          
         prod(:,84) = 0._r8
                                                                                          
         prod(:,195) = 0._r8
                                                                                          
         prod(:,104) = 0._r8
                                                                                          
         prod(:,64) = 0._r8
                                                                                          
         prod(:,173) = 0._r8
                                                                                          
         prod(:,144) = 0._r8
                                                                                          
         prod(:,141) = 0._r8
                                                                                          
         prod(:,198) = 0._r8
                                                                                          
         prod(:,117) = 0._r8
                                                                                          
         prod(:,157) = 0._r8
                                                                                          
         prod(:,48) = 0._r8
                                                                                          
         prod(:,199) = 0._r8
                                                                                          
         prod(:,102) = 0._r8
                                                                                          
         prod(:,134) = 0._r8
                                                                                          
         prod(:,103) = 0._r8
                                                                                          
         prod(:,143) = 0._r8
                                                                                          
         prod(:,182) = 0._r8
                                                                                          
         prod(:,205) = 0._r8
                                                                                          
         prod(:,132) = + extfrc(:,14)
                                                                                          
         prod(:,75) = 0._r8
                                                                                          
         prod(:,95) = 0._r8
                                                                                          
         prod(:,113) = 0._r8
                                                                                          
         prod(:,189) = 0._r8
                                                                                          
         prod(:,10) = 0._r8
                                                                                          
         prod(:,11) = 0._r8
                                                                                          
         prod(:,12) = 0._r8
                                                                                          
         prod(:,53) = 0._r8
                                                                                          
         prod(:,13) = 0._r8
                                                                                          
         prod(:,14) = 0._r8
                                                                                          
         prod(:,15) = 0._r8
                                                                                          
         prod(:,217) = + extfrc(:,13)
                                                                                          
         prod(:,224) = + extfrc(:,9)
                                                                                          
         prod(:,216) = 0._r8
                                                                                          
         prod(:,174) = 0._r8
                                                                                          
         prod(:,116) = 0._r8
                                                                                          
         prod(:,16) = + extfrc(:,10)
                                                                                          
         prod(:,17) = + extfrc(:,11)
                                                                                          
         prod(:,18) = 0._r8
                                                                                          
         prod(:,19) = + extfrc(:,1)
                                                                                          
         prod(:,226) = (rxt(:,5) +2.000_r8*rxt(:,6))
                                                                                          
         prod(:,222) = 0._r8
                                                                                          
         prod(:,20) = 0._r8
                                                                                          
         prod(:,106) = 0._r8
                                                                                          
         prod(:,111) = 0._r8
                                                                                          
         prod(:,88) = 0._r8
                                                                                          
         prod(:,137) = 0._r8
                                                                                          
         prod(:,65) = 0._r8
                                                                                          
         prod(:,128) = 0._r8
                                                                                          
         prod(:,73) = 0._r8
                                                                                          
         prod(:,107) = 0._r8
                                                                                          
         prod(:,21) = 0._r8
                                                                                          
         prod(:,22) = + extfrc(:,2)
                                                                                          
         prod(:,138) = 0._r8
                                                                                          
         prod(:,119) = 0._r8
                                                                                          
         prod(:,135) = 0._r8
                                                                                          
         prod(:,23) = 0._r8
                                                                                          
         prod(:,201) = 0._r8
                                                                                          
         prod(:,172) = + extfrc(:,8)
                                                                                          
         prod(:,91) = 0._r8
                                                                                          
         prod(:,24) = + extfrc(:,5)
                                                                                          
         prod(:,25) = + extfrc(:,6)
                                                                                          
         prod(:,26) = 0._r8
                                                                                          
         prod(:,27) = 0._r8
                                                                                          
         prod(:,28) = 0._r8
                                                                                          
         prod(:,29) = 0._r8
                                                                                          
         prod(:,30) = 0._r8
                                                                                          
         prod(:,31) = 0._r8
                                                                                          
         prod(:,32) = 0._r8
                                                                                          
         prod(:,33) = 0._r8
                                                                                          
         prod(:,34) = 0._r8
                                                                                          
         prod(:,35) = 0._r8
                                                                                          
         prod(:,36) = 0._r8
                                                                                          
         prod(:,37) = 0._r8
                                                                                          
         prod(:,38) = 0._r8
                                                                                          
         prod(:,39) = 0._r8
                                                                                          
         prod(:,40) = 0._r8
                                                                                          
         prod(:,41) = 0._r8
                                                                                          
         prod(:,42) = 0._r8
                                                                                          
         prod(:,43) = + extfrc(:,4)
                                                                                          
         prod(:,78) = 0._r8
                                                                                          
         prod(:,152) = 0._r8
                                                                                          
         prod(:,149) = 0._r8
                                                                                          
         prod(:,130) = 0._r8
                                                                                          
         prod(:,183) = 0._r8
                                                                                          
         prod(:,188) = 0._r8
                                                                                          
         prod(:,153) = 0._r8
                                                                                          
         prod(:,76) = 0._r8
                                                                                          
         prod(:,79) = 0._r8
                                                                                          
         prod(:,80) = 0._r8
                                                                                          
         prod(:,158) = 0._r8
                                                                                          
         prod(:,81) = 0._r8
                                                                                          
         prod(:,122) = 0._r8
                                                                                          
         prod(:,136) = 0._r8
                                                                                          
         prod(:,178) = 0._r8
                                                                                          
         prod(:,44) = 0._r8
                                                                                          
         prod(:,131) = 0._r8
                                                                                          
         prod(:,45) = 0._r8
                                                                                          
         prod(:,123) = 0._r8
                                                                                          
         prod(:,170) = 0._r8
                                                                                          
         prod(:,167) = 0._r8
                                                                                          
         prod(:,150) = 0._r8
                                                                                          
         prod(:,207) = 0._r8
                                                                                          
         prod(:,221) = 0._r8
                                                                                          
         prod(:,163) = 0._r8
                                                                                          
         prod(:,142) = 0._r8
                                                                                          
         prod(:,94) = 0._r8
                                                                                          
         prod(:,159) = 0._r8
                                                                                          
         prod(:,218) = 0._r8
                                                                                          
         prod(:,124) = 0._r8
                                                                                          
         prod(:,202) = 0._r8
                                                                                          
         prod(:,203) = 0._r8
                                                                                          
         prod(:,46) = 0._r8
                                                                                          
         prod(:,47) = 0._r8
                                                                                          
         prod(:,204) = 0._r8
                                                                                          
         prod(:,160) = 0._r8
                                                                                          
         prod(:,206) = 0._r8
                                                                                          
         prod(:,176) = 0._r8
                                                                                          
         prod(:,155) = 0._r8
                                                                                          
         prod(:,49) = 0._r8
                                                                                          
         prod(:,187) = 0._r8
                                                                                          
         prod(:,214) =rxt(:,5)
                                                                                          
         prod(:,215) = 0._r8
                                                                                          
         prod(:,126) = 0._r8
                                                                                          
         prod(:,164) = 0._r8
                                                                                          
         prod(:,194) = 0._r8
                                                                                          
         prod(:,193) = 0._r8
                                                                                          
         prod(:,179) = 0._r8
                                                                                          
         prod(:,161) = 0._r8
                                                                                          
         prod(:,50) = 0._r8
                                                                                          
         prod(:,197) = 0._r8
                                                                                          
         prod(:,162) = 0._r8
                                                                                          
         prod(:,51) = 0._r8
                                                                                          
         prod(:,133) = 0._r8
                                                                                          
         prod(:,227) = 0._r8
                                                                                          
      end if                                                                              
                                                                                          
      end subroutine indprd                                                               
                                                                                          
      end module mo_indprd                                                                

      module mo_lin_matrix

      use chem_mods, only: veclen
      private
      public :: linmat

      contains

      subroutine linmat01( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer,  intent(in)    ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(in)    ::  het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)
!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k


      do k = 1,avec_len

         mat(k,666) = -( rxt(k,19) + het_rates(k,1) )

         mat(k,632) = -( rxt(k,20) + het_rates(k,2) )

         mat(k,1) = -( het_rates(k,3) )

         mat(k,2) = -( het_rates(k,4) )

         mat(k,3) = -( het_rates(k,5) )

         mat(k,999) = -( het_rates(k,6) )

         mat(k,164) = -( het_rates(k,7) )

         mat(k,421) = -( rxt(k,21) + het_rates(k,8) )

         mat(k,170) = -( rxt(k,22) + het_rates(k,9) )

         mat(k,379) = -( rxt(k,23) + het_rates(k,10) )

         mat(k,460) = -( rxt(k,24) + het_rates(k,11) )
         mat(k,422) = .500_r8*rxt(k,21)
         mat(k,171) = rxt(k,22)
         mat(k,653) = .200_r8*rxt(k,70)
         mat(k,699) = .060_r8*rxt(k,72)

         mat(k,279) = -( rxt(k,25) + het_rates(k,12) )
         mat(k,652) = .200_r8*rxt(k,70)
         mat(k,697) = .200_r8*rxt(k,72)

         mat(k,590) = -( rxt(k,26) + het_rates(k,13) )
         mat(k,231) = rxt(k,46)
         mat(k,1069) = rxt(k,56)
         mat(k,654) = .200_r8*rxt(k,70)
         mat(k,700) = .150_r8*rxt(k,72)

         mat(k,323) = -( rxt(k,27) + het_rates(k,14) )
         mat(k,698) = .210_r8*rxt(k,72)

         mat(k,238) = -( het_rates(k,15) )

         mat(k,349) = -( het_rates(k,16) )

         mat(k,1415) = -( het_rates(k,17) )
         mat(k,242) = rxt(k,74)
         mat(k,2211) = rxt(k,75)
         mat(k,542) = rxt(k,77)
         mat(k,143) = rxt(k,79)
         mat(k,149) = rxt(k,80)
         mat(k,468) = 2.000_r8*rxt(k,86)
         mat(k,595) = rxt(k,87)
         mat(k,448) = 3.000_r8*rxt(k,90)
         mat(k,107) = 2.000_r8*rxt(k,98)
         mat(k,805) = rxt(k,99)
         mat(k,779) = rxt(k,105)

         mat(k,241) = -( rxt(k,74) + het_rates(k,18) )

         mat(k,2226) = -( rxt(k,75) + het_rates(k,19) )
         mat(k,546) = rxt(k,76)

         mat(k,540) = -( rxt(k,76) + rxt(k,77) + rxt(k,525) + rxt(k,528) + rxt(k,533) &
                 + het_rates(k,20) )

         mat(k,4) = -( het_rates(k,21) )

         mat(k,253) = -( het_rates(k,22) )
         mat(k,338) = rxt(k,28)

         mat(k,339) = -( rxt(k,28) + het_rates(k,23) )

         mat(k,285) = -( het_rates(k,24) )

         mat(k,548) = -( het_rates(k,25) )

         mat(k,261) = -( het_rates(k,26) )

         mat(k,344) = -( rxt(k,29) + het_rates(k,27) )

         mat(k,291) = -( het_rates(k,28) )

         mat(k,1024) = -( het_rates(k,29) )
         mat(k,1326) = .700_r8*rxt(k,55)

         mat(k,415) = -( rxt(k,30) + het_rates(k,30) )

         mat(k,109) = -( het_rates(k,31) )

         mat(k,265) = -( rxt(k,31) + het_rates(k,32) )

         mat(k,99) = -( rxt(k,78) + het_rates(k,33) )

         mat(k,141) = -( rxt(k,79) + het_rates(k,34) )

         mat(k,146) = -( rxt(k,80) + het_rates(k,35) )

         mat(k,113) = -( rxt(k,81) + het_rates(k,36) )

         mat(k,151) = -( rxt(k,82) + het_rates(k,37) )

         mat(k,117) = -( rxt(k,83) + het_rates(k,38) )

         mat(k,156) = -( rxt(k,84) + het_rates(k,39) )

         mat(k,121) = -( rxt(k,85) + het_rates(k,40) )

         mat(k,467) = -( rxt(k,86) + het_rates(k,41) )

         mat(k,1485) = -( rxt(k,32) + rxt(k,33) + het_rates(k,42) )
         mat(k,672) = .100_r8*rxt(k,19)
         mat(k,639) = .100_r8*rxt(k,20)
         mat(k,387) = rxt(k,38)
         mat(k,1433) = .180_r8*rxt(k,39)
         mat(k,1098) = rxt(k,43)
         mat(k,1157) = .330_r8*rxt(k,45)
         mat(k,1143) = rxt(k,47)
         mat(k,694) = rxt(k,49)
         mat(k,1214) = 1.340_r8*rxt(k,50)
         mat(k,860) = rxt(k,57)
         mat(k,536) = rxt(k,62)
         mat(k,412) = rxt(k,63)
         mat(k,649) = .375_r8*rxt(k,65)
         mat(k,478) = .400_r8*rxt(k,67)
         mat(k,1063) = .680_r8*rxt(k,69)
         mat(k,443) = rxt(k,270)
         mat(k,271) = 2.000_r8*rxt(k,300)

         mat(k,594) = -( rxt(k,87) + het_rates(k,43) )

         mat(k,125) = -( rxt(k,88) + het_rates(k,44) )

         mat(k,1085) = -( rxt(k,34) + het_rates(k,45) )
         mat(k,670) = .400_r8*rxt(k,19)
         mat(k,637) = .400_r8*rxt(k,20)
         mat(k,346) = rxt(k,29)
         mat(k,1148) = .330_r8*rxt(k,45)
         mat(k,317) = rxt(k,53)
         mat(k,534) = rxt(k,62)

         mat(k,365) = -( rxt(k,89) + het_rates(k,46) )

         mat(k,102) = -( het_rates(k,47) )

         mat(k,922) = -( rxt(k,35) + het_rates(k,48) )
         mat(k,669) = .250_r8*rxt(k,19)
         mat(k,636) = .250_r8*rxt(k,20)
         mat(k,417) = .820_r8*rxt(k,30)
         mat(k,1147) = .170_r8*rxt(k,45)
         mat(k,644) = .300_r8*rxt(k,65)
         mat(k,476) = .050_r8*rxt(k,67)
         mat(k,1058) = .500_r8*rxt(k,69)

         mat(k,1221) = -( rxt(k,36) + het_rates(k,49) )
         mat(k,382) = .180_r8*rxt(k,23)
         mat(k,325) = rxt(k,27)
         mat(k,662) = .400_r8*rxt(k,70)
         mat(k,708) = .540_r8*rxt(k,72)
         mat(k,430) = .510_r8*rxt(k,73)

         mat(k,684) = -( het_rates(k,50) )

         mat(k,610) = -( rxt(k,37) + het_rates(k,51) )

         mat(k,786) = -( het_rates(k,52) )

         mat(k,385) = -( rxt(k,38) + het_rates(k,53) )

         mat(k,1430) = -( rxt(k,39) + rxt(k,40) + het_rates(k,54) )

         mat(k,447) = -( rxt(k,90) + het_rates(k,55) )

         mat(k,2017) = -( het_rates(k,56) )
         mat(k,243) = rxt(k,74)
         mat(k,101) = 4.000_r8*rxt(k,78)
         mat(k,145) = rxt(k,79)
         mat(k,116) = 2.000_r8*rxt(k,81)
         mat(k,155) = 2.000_r8*rxt(k,82)
         mat(k,120) = 2.000_r8*rxt(k,83)
         mat(k,160) = rxt(k,84)
         mat(k,124) = 2.000_r8*rxt(k,85)
         mat(k,127) = 3.000_r8*rxt(k,88)
         mat(k,369) = rxt(k,89)
         mat(k,162) = 2.000_r8*rxt(k,91)
         mat(k,95) = 2.000_r8*rxt(k,92)
         mat(k,1978) = rxt(k,93)
         mat(k,889) = rxt(k,94)
         mat(k,229) = rxt(k,97)
         mat(k,225) = rxt(k,100)
         mat(k,252) = rxt(k,101)
         mat(k,308) = rxt(k,102)
         mat(k,2153) = rxt(k,103)
         mat(k,827) = rxt(k,106)

         mat(k,161) = -( rxt(k,91) + het_rates(k,57) )

         mat(k,93) = -( rxt(k,92) + rxt(k,211) + het_rates(k,58) )

         mat(k,1977) = -( rxt(k,93) + het_rates(k,59) )
         mat(k,888) = rxt(k,95)
         mat(k,331) = rxt(k,107)
         mat(k,94) = 2.000_r8*rxt(k,211)

         mat(k,884) = -( rxt(k,94) + rxt(k,95) + rxt(k,527) + rxt(k,532) + rxt(k,538) &
                 + het_rates(k,60) )

         mat(k,5) = -( het_rates(k,61) )

         mat(k,1103) = -( het_rates(k,62) )
         mat(k,172) = 1.500_r8*rxt(k,22)
         mat(k,381) = .450_r8*rxt(k,23)
         mat(k,592) = .600_r8*rxt(k,26)
         mat(k,324) = rxt(k,27)
         mat(k,1479) = rxt(k,32) + rxt(k,33)
         mat(k,1086) = rxt(k,34)
         mat(k,1220) = rxt(k,36)
         mat(k,1428) = .380_r8*rxt(k,39)
         mat(k,802) = rxt(k,41)
         mat(k,1097) = rxt(k,43)
         mat(k,980) = 2.000_r8*rxt(k,44)
         mat(k,1150) = .330_r8*rxt(k,45)
         mat(k,1208) = 1.340_r8*rxt(k,51)
         mat(k,1328) = .700_r8*rxt(k,55)
         mat(k,200) = 1.500_r8*rxt(k,64)
         mat(k,647) = .250_r8*rxt(k,65)
         mat(k,972) = rxt(k,68)
         mat(k,1060) = 1.700_r8*rxt(k,69)
         mat(k,360) = rxt(k,110)

         mat(k,801) = -( rxt(k,41) + het_rates(k,63) )
         mat(k,611) = rxt(k,37)
         mat(k,1426) = .440_r8*rxt(k,39)
         mat(k,525) = .400_r8*rxt(k,60)
         mat(k,643) = rxt(k,65)
         mat(k,1057) = .800_r8*rxt(k,69)

         mat(k,235) = -( rxt(k,96) + het_rates(k,64) )
         mat(k,142) = rxt(k,79)
         mat(k,147) = rxt(k,80)
         mat(k,153) = rxt(k,82)
         mat(k,118) = 2.000_r8*rxt(k,83)
         mat(k,157) = 2.000_r8*rxt(k,84)
         mat(k,122) = rxt(k,85)
         mat(k,106) = 2.000_r8*rxt(k,98)
         mat(k,247) = rxt(k,101)
         mat(k,303) = rxt(k,102)

         mat(k,226) = -( rxt(k,97) + het_rates(k,65) )
         mat(k,114) = rxt(k,81)
         mat(k,152) = rxt(k,82)
         mat(k,222) = rxt(k,100)

         mat(k,194) = -( het_rates(k,66) )

         mat(k,297) = -( het_rates(k,67) )

         mat(k,6) = -( het_rates(k,68) )

         mat(k,7) = -( het_rates(k,69) )

         mat(k,8) = -( het_rates(k,70) )

         mat(k,9) = -( rxt(k,124) + het_rates(k,71) )

         mat(k,129) = -( rxt(k,42) + het_rates(k,72) )

         mat(k,864) = -( het_rates(k,73) )
         mat(k,148) = rxt(k,80)
         mat(k,158) = rxt(k,84)
         mat(k,236) = 2.000_r8*rxt(k,96)
         mat(k,227) = rxt(k,97)
         mat(k,283) = rxt(k,104)

         mat(k,1096) = -( rxt(k,43) + het_rates(k,74) )
         mat(k,1149) = .330_r8*rxt(k,45)
         mat(k,646) = .250_r8*rxt(k,65)
         mat(k,270) = rxt(k,301)

      end do

      end subroutine linmat01

      subroutine linmat02( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer,  intent(in)    ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(in)    ::  het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)
!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k


      do k = 1,avec_len

         mat(k,979) = -( rxt(k,44) + rxt(k,519) + het_rates(k,75) )
         mat(k,424) = rxt(k,21)
         mat(k,380) = .130_r8*rxt(k,23)
         mat(k,335) = .700_r8*rxt(k,61)
         mat(k,660) = .600_r8*rxt(k,70)
         mat(k,706) = .340_r8*rxt(k,72)
         mat(k,429) = .170_r8*rxt(k,73)

         mat(k,1463) = -( rxt(k,139) + het_rates(k,76) )
         mat(k,2270) = rxt(k,2) + 2.000_r8*rxt(k,3)
         mat(k,1484) = 2.000_r8*rxt(k,33)
         mat(k,386) = rxt(k,38)
         mat(k,1432) = .330_r8*rxt(k,39) + rxt(k,40)
         mat(k,806) = rxt(k,99)
         mat(k,2145) = rxt(k,103)
         mat(k,284) = rxt(k,104)

         mat(k,1401) = -( het_rates(k,77) )
         mat(k,2266) = rxt(k,1)
         mat(k,1480) = rxt(k,32)
         mat(k,1429) = 1.440_r8*rxt(k,39)

         mat(k,105) = -( rxt(k,98) + het_rates(k,78) )

         mat(k,603) = -( rxt(k,4) + het_rates(k,79) )

         mat(k,132) = -( rxt(k,109) + het_rates(k,80) )

         mat(k,804) = -( rxt(k,99) + het_rates(k,81) )

         mat(k,221) = -( rxt(k,100) + het_rates(k,82) )

         mat(k,248) = -( rxt(k,101) + het_rates(k,83) )

         mat(k,304) = -( rxt(k,102) + het_rates(k,84) )

         mat(k,2156) = -( rxt(k,103) + het_rates(k,85) )

         mat(k,179) = -( het_rates(k,86) )

         mat(k,929) = -( het_rates(k,87) )

         mat(k,282) = -( rxt(k,104) + het_rates(k,88) )

         mat(k,1447) = -( rxt(k,9) + het_rates(k,89) )
         mat(k,1156) = rxt(k,482)
         mat(k,585) = rxt(k,483)
         mat(k,561) = rxt(k,484)
         mat(k,274) = 2.000_r8*rxt(k,485) + 2.000_r8*rxt(k,523) + 2.000_r8*rxt(k,526) &
                      + 2.000_r8*rxt(k,537)
         mat(k,376) = rxt(k,486)
         mat(k,1077) = rxt(k,487)
         mat(k,2188) = .500_r8*rxt(k,489)
         mat(k,1744) = rxt(k,490)
         mat(k,394) = rxt(k,491)
         mat(k,245) = rxt(k,492)
         mat(k,619) = rxt(k,493)
         mat(k,543) = rxt(k,525) + rxt(k,528) + rxt(k,533)
         mat(k,885) = rxt(k,527) + rxt(k,532) + rxt(k,538)

         mat(k,403) = -( rxt(k,10) + rxt(k,11) + rxt(k,174) + het_rates(k,90) )

         mat(k,778) = -( rxt(k,105) + het_rates(k,91) )
         mat(k,541) = rxt(k,525) + rxt(k,528) + rxt(k,533)

         mat(k,824) = -( rxt(k,106) + het_rates(k,92) )
         mat(k,883) = rxt(k,527) + rxt(k,532) + rxt(k,538)

         mat(k,1153) = -( rxt(k,45) + rxt(k,482) + het_rates(k,93) )

         mat(k,230) = -( rxt(k,46) + het_rates(k,94) )
         mat(k,1272) = rxt(k,374)

         mat(k,1140) = -( rxt(k,47) + het_rates(k,95) )
         mat(k,1152) = .170_r8*rxt(k,45)

         mat(k,320) = -( het_rates(k,96) )

         mat(k,135) = -( het_rates(k,97) )

         mat(k,841) = -( het_rates(k,98) )

         mat(k,581) = -( rxt(k,483) + het_rates(k,99) )

         mat(k,556) = -( rxt(k,484) + het_rates(k,100) )

         mat(k,1193) = -( het_rates(k,101) )

         mat(k,397) = -( rxt(k,48) + het_rates(k,102) )

         mat(k,690) = -( rxt(k,49) + het_rates(k,103) )
         mat(k,398) = rxt(k,48)

         mat(k,74) = -( het_rates(k,104) )

         mat(k,1209) = -( rxt(k,50) + rxt(k,51) + het_rates(k,105) )
         mat(k,692) = .300_r8*rxt(k,49)

         mat(k,310) = -( het_rates(k,106) )

         mat(k,506) = -( rxt(k,52) + het_rates(k,107) )
         mat(k,665) = .800_r8*rxt(k,19)
         mat(k,631) = .800_r8*rxt(k,20)

         mat(k,315) = -( rxt(k,53) + het_rates(k,108) )

         mat(k,572) = -( rxt(k,54) + rxt(k,356) + het_rates(k,109) )

         mat(k,948) = -( het_rates(k,110) )

         mat(k,1332) = -( rxt(k,55) + het_rates(k,111) )
         mat(k,693) = .700_r8*rxt(k,49)

         mat(k,491) = -( rxt(k,156) + het_rates(k,112) )
         mat(k,1785) = rxt(k,15)

         mat(k,183) = -( rxt(k,12) + het_rates(k,113) )

         mat(k,273) = -( rxt(k,13) + rxt(k,14) + rxt(k,175) + rxt(k,485) + rxt(k,523) &
                      + rxt(k,526) + rxt(k,537) + het_rates(k,114) )

         mat(k,373) = -( rxt(k,486) + het_rates(k,115) )

         mat(k,1073) = -( rxt(k,56) + rxt(k,487) + het_rates(k,116) )

         mat(k,10) = -( het_rates(k,117) )

         mat(k,11) = -( het_rates(k,118) )

         mat(k,12) = -( het_rates(k,119) )

         mat(k,96) = -( het_rates(k,120) )

         mat(k,13) = -( rxt(k,488) + het_rates(k,121) )

         mat(k,14) = -( rxt(k,541) + het_rates(k,122) )

         mat(k,15) = -( rxt(k,540) + het_rates(k,123) )

         mat(k,1842) = -( rxt(k,15) + het_rates(k,124) )
         mat(k,276) = rxt(k,14)
         mat(k,2194) = rxt(k,16) + .500_r8*rxt(k,489)
         mat(k,1750) = rxt(k,17)
         mat(k,495) = rxt(k,156)

         mat(k,2201) = -( rxt(k,16) + rxt(k,489) + het_rates(k,125) )
         mat(k,1457) = rxt(k,9)
         mat(k,407) = rxt(k,11) + rxt(k,174)
         mat(k,277) = rxt(k,13) + rxt(k,175)
         mat(k,1757) = rxt(k,18)
         mat(k,675) = rxt(k,19)
         mat(k,1163) = rxt(k,45)
         mat(k,402) = rxt(k,48)
         mat(k,580) = rxt(k,54) + rxt(k,356)
         mat(k,1083) = rxt(k,56)
         mat(k,862) = rxt(k,57)
         mat(k,396) = rxt(k,58)
         mat(k,246) = rxt(k,59)
         mat(k,531) = .600_r8*rxt(k,60) + rxt(k,307)
         mat(k,622) = rxt(k,66)
         mat(k,545) = rxt(k,76)
         mat(k,891) = rxt(k,95)
         mat(k,140) = rxt(k,431)

         mat(k,1749) = -( rxt(k,17) + rxt(k,18) + rxt(k,490) + het_rates(k,126) )
         mat(k,405) = rxt(k,10)
         mat(k,275) = rxt(k,13) + rxt(k,14) + rxt(k,175)
         mat(k,529) = .400_r8*rxt(k,60)
         mat(k,544) = rxt(k,77)
         mat(k,887) = rxt(k,94)

         mat(k,857) = -( rxt(k,57) + het_rates(k,127) )

         mat(k,391) = -( rxt(k,58) + rxt(k,491) + het_rates(k,128) )

         mat(k,16) = -( het_rates(k,129) )

         mat(k,17) = -( het_rates(k,130) )

         mat(k,18) = -( het_rates(k,131) )

         mat(k,19) = -( het_rates(k,132) )

         mat(k,2258) = -( rxt(k,133) + het_rates(k,133) )
         mat(k,2284) = rxt(k,3)
         mat(k,2136) = rxt(k,8)
         mat(k,278) = rxt(k,14)
         mat(k,1851) = rxt(k,15)
         mat(k,2203) = rxt(k,16)
         mat(k,1759) = rxt(k,18)
         mat(k,1441) = .180_r8*rxt(k,39)
         mat(k,803) = rxt(k,41)
         mat(k,2227) = rxt(k,75)
         mat(k,1984) = rxt(k,93)
         mat(k,332) = rxt(k,107)
         mat(k,1243) = rxt(k,111) + rxt(k,474)
         mat(k,836) = rxt(k,112)
         mat(k,259) = rxt(k,113)
         mat(k,1538) = rxt(k,127) + rxt(k,128)
         mat(k,497) = rxt(k,156)
         mat(k,516) = rxt(k,467)

         mat(k,2132) = -( rxt(k,7) + rxt(k,8) + het_rates(k,134) )
         mat(k,2254) = rxt(k,133)

         mat(k,20) = -( het_rates(k,135) )

         mat(k,328) = -( rxt(k,107) + het_rates(k,136) )

         mat(k,357) = -( rxt(k,110) + het_rates(k,137) )

         mat(k,244) = -( rxt(k,59) + rxt(k,492) + het_rates(k,138) )

         mat(k,524) = -( rxt(k,60) + rxt(k,307) + het_rates(k,139) )

         mat(k,138) = -( rxt(k,431) + het_rates(k,140) )

         mat(k,463) = -( het_rates(k,141) )
         mat(k,266) = rxt(k,31)

         mat(k,174) = -( het_rates(k,142) )

         mat(k,333) = -( rxt(k,61) + het_rates(k,143) )

         mat(k,21) = -( het_rates(k,144) )

         mat(k,22) = -( het_rates(k,145) )

         mat(k,532) = -( rxt(k,62) + het_rates(k,146) )

         mat(k,409) = -( rxt(k,63) + het_rates(k,147) )

         mat(k,511) = -( rxt(k,467) + het_rates(k,148) )
         mat(k,358) = rxt(k,110)
         mat(k,1230) = rxt(k,111)

         mat(k,23) = -( rxt(k,108) + het_rates(k,149) )

         mat(k,1232) = -( rxt(k,111) + rxt(k,474) + het_rates(k,150) )
         mat(k,833) = rxt(k,112)
         mat(k,512) = rxt(k,467)

         mat(k,832) = -( rxt(k,112) + het_rates(k,151) )
         mat(k,258) = rxt(k,113)
         mat(k,1231) = rxt(k,474)

         mat(k,257) = -( rxt(k,113) + het_rates(k,152) )
         mat(k,133) = rxt(k,109)

         mat(k,24) = -( het_rates(k,153) )

         mat(k,25) = -( het_rates(k,154) )

         mat(k,26) = -( het_rates(k,155) )

         mat(k,27) = -( rxt(k,114) + het_rates(k,156) )

         mat(k,28) = -( rxt(k,115) + het_rates(k,157) )

         mat(k,29) = -( rxt(k,116) + het_rates(k,158) )

         mat(k,30) = -( rxt(k,117) + het_rates(k,159) )

         mat(k,31) = -( rxt(k,118) + het_rates(k,160) )

         mat(k,32) = -( rxt(k,119) + het_rates(k,161) )

         mat(k,33) = -( rxt(k,120) + het_rates(k,162) )

         mat(k,34) = -( rxt(k,121) + het_rates(k,163) )

         mat(k,35) = -( rxt(k,122) + het_rates(k,164) )

         mat(k,36) = -( rxt(k,123) + het_rates(k,165) )

         mat(k,37) = -( het_rates(k,166) )
         mat(k,977) = rxt(k,519)

         mat(k,38) = -( het_rates(k,167) )

         mat(k,39) = -( het_rates(k,168) )

         mat(k,40) = -( het_rates(k,169) )

         mat(k,41) = -( het_rates(k,170) )

         mat(k,42) = -( rxt(k,542) + het_rates(k,171) )

         mat(k,48) = -( het_rates(k,172) )

         mat(k,199) = -( rxt(k,64) + het_rates(k,173) )

         mat(k,642) = -( rxt(k,65) + het_rates(k,174) )

         mat(k,617) = -( rxt(k,66) + rxt(k,493) + het_rates(k,175) )

         mat(k,474) = -( rxt(k,67) + het_rates(k,176) )

         mat(k,969) = -( rxt(k,68) + het_rates(k,177) )
         mat(k,392) = rxt(k,58)
         mat(k,618) = rxt(k,66)
         mat(k,477) = rxt(k,67)

         mat(k,1059) = -( rxt(k,69) + het_rates(k,178) )
         mat(k,645) = rxt(k,65)
         mat(k,971) = rxt(k,68)

      end do

      end subroutine linmat02

      subroutine linmat03( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer,  intent(in)    ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(in)    ::  het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)
!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k


      do k = 1,avec_len

         mat(k,655) = -( rxt(k,70) + het_rates(k,179) )

         mat(k,187) = -( het_rates(k,180) )

         mat(k,203) = -( rxt(k,71) + het_rates(k,181) )

         mat(k,208) = -( het_rates(k,182) )

         mat(k,701) = -( rxt(k,72) + het_rates(k,183) )

         mat(k,216) = -( het_rates(k,184) )

         mat(k,427) = -( rxt(k,73) + het_rates(k,185) )

         mat(k,518) = -( het_rates(k,188) )
         mat(k,139) = rxt(k,431)

         mat(k,899) = -( het_rates(k,189) )

         mat(k,54) = -( het_rates(k,190) )

         mat(k,483) = -( het_rates(k,191) )

         mat(k,60) = -( het_rates(k,192) )

         mat(k,435) = -( het_rates(k,193) )

         mat(k,815) = -( het_rates(k,194) )
         mat(k,508) = rxt(k,52)

         mat(k,791) = -( het_rates(k,195) )

         mat(k,625) = -( het_rates(k,196) )

         mat(k,1386) = -( het_rates(k,197) )
         mat(k,383) = .130_r8*rxt(k,23)
         mat(k,326) = rxt(k,27)
         mat(k,924) = rxt(k,35)
         mat(k,1222) = rxt(k,36)
         mat(k,1155) = .330_r8*rxt(k,45)
         mat(k,1142) = rxt(k,47)
         mat(k,1213) = 1.340_r8*rxt(k,50)
         mat(k,509) = rxt(k,52)
         mat(k,318) = rxt(k,53)
         mat(k,1334) = .300_r8*rxt(k,55)
         mat(k,859) = rxt(k,57)
         mat(k,526) = .600_r8*rxt(k,60) + rxt(k,307)
         mat(k,411) = rxt(k,63)
         mat(k,201) = .500_r8*rxt(k,64)
         mat(k,1062) = .650_r8*rxt(k,69)

         mat(k,2070) = -( het_rates(k,198) )
         mat(k,1092) = rxt(k,34)
         mat(k,926) = rxt(k,35)
         mat(k,615) = rxt(k,37)
         mat(k,1439) = rxt(k,40)
         mat(k,1342) = .300_r8*rxt(k,55)
         mat(k,530) = .400_r8*rxt(k,60)
         mat(k,600) = rxt(k,87)
         mat(k,370) = rxt(k,89)

         mat(k,759) = -( het_rates(k,199) )
         mat(k,280) = .600_r8*rxt(k,25)

         mat(k,564) = -( het_rates(k,200) )

         mat(k,269) = -( rxt(k,300) + rxt(k,301) + het_rates(k,201) )
         mat(k,130) = rxt(k,42)

         mat(k,714) = -( het_rates(k,202) )

         mat(k,1950) = -( rxt(k,481) + het_rates(k,203) )
         mat(k,406) = rxt(k,11) + rxt(k,174)
         mat(k,674) = rxt(k,19)
         mat(k,641) = .900_r8*rxt(k,20)
         mat(k,426) = rxt(k,21)
         mat(k,173) = 1.500_r8*rxt(k,22)
         mat(k,384) = .560_r8*rxt(k,23)
         mat(k,462) = rxt(k,24)
         mat(k,281) = .600_r8*rxt(k,25)
         mat(k,593) = .600_r8*rxt(k,26)
         mat(k,327) = rxt(k,27)
         mat(k,343) = rxt(k,28)
         mat(k,348) = rxt(k,29)
         mat(k,419) = rxt(k,30)
         mat(k,1091) = rxt(k,34)
         mat(k,1226) = rxt(k,36)
         mat(k,1100) = 2.000_r8*rxt(k,43)
         mat(k,982) = 2.000_r8*rxt(k,44)
         mat(k,1161) = .670_r8*rxt(k,45)
         mat(k,234) = rxt(k,46)
         mat(k,1145) = rxt(k,47)
         mat(k,401) = rxt(k,48)
         mat(k,696) = rxt(k,49)
         mat(k,1216) = 1.340_r8*rxt(k,50) + .660_r8*rxt(k,51)
         mat(k,1081) = rxt(k,56)
         mat(k,337) = rxt(k,61)
         mat(k,538) = rxt(k,62)
         mat(k,202) = rxt(k,64)
         mat(k,651) = rxt(k,65)
         mat(k,621) = rxt(k,66)
         mat(k,480) = rxt(k,67)
         mat(k,976) = rxt(k,68)
         mat(k,1065) = 1.200_r8*rxt(k,69)
         mat(k,664) = rxt(k,70)
         mat(k,711) = rxt(k,72)
         mat(k,432) = rxt(k,73)
         mat(k,1468) = rxt(k,139)
         mat(k,445) = rxt(k,270)
         mat(k,272) = rxt(k,300) + rxt(k,301)
         mat(k,1298) = rxt(k,374)

         mat(k,441) = -( rxt(k,270) + het_rates(k,204) )

         mat(k,1256) = -( het_rates(k,205) )

         mat(k,1288) = -( rxt(k,374) + het_rates(k,206) )

         mat(k,66) = -( het_rates(k,207) )

         mat(k,72) = -( het_rates(k,208) )

         mat(k,1311) = -( het_rates(k,209) )

         mat(k,721) = -( het_rates(k,210) )
         mat(k,461) = .600_r8*rxt(k,24)

         mat(k,1354) = -( het_rates(k,211) )
         mat(k,1212) = .660_r8*rxt(k,50)
         mat(k,575) = rxt(k,54) + rxt(k,356)

         mat(k,873) = -( het_rates(k,212) )
         mat(k,591) = .600_r8*rxt(k,26)

         mat(k,677) = -( het_rates(k,213) )

         mat(k,80) = -( het_rates(k,214) )

         mat(k,1045) = -( het_rates(k,215) )

         mat(k,1526) = -( rxt(k,127) + rxt(k,128) + het_rates(k,216) )
         mat(k,2272) = rxt(k,1)
         mat(k,2124) = rxt(k,7)
         mat(k,184) = rxt(k,12)

         mat(k,1691) = -( het_rates(k,217) )
         mat(k,2273) = rxt(k,2)
         mat(k,604) = 2.000_r8*rxt(k,4)
         mat(k,1451) = rxt(k,9)
         mat(k,404) = rxt(k,10)
         mat(k,640) = rxt(k,20)
         mat(k,425) = rxt(k,21)
         mat(k,342) = rxt(k,28)
         mat(k,347) = rxt(k,29)
         mat(k,418) = rxt(k,30)
         mat(k,268) = rxt(k,31)
         mat(k,614) = rxt(k,37)
         mat(k,388) = rxt(k,38)
         mat(k,1435) = .330_r8*rxt(k,39)
         mat(k,131) = rxt(k,42)
         mat(k,233) = rxt(k,46)
         mat(k,695) = rxt(k,49)
         mat(k,319) = rxt(k,53)
         mat(k,395) = rxt(k,58)
         mat(k,336) = rxt(k,61)
         mat(k,537) = rxt(k,62)
         mat(k,413) = rxt(k,63)
         mat(k,650) = rxt(k,65)
         mat(k,479) = rxt(k,67)
         mat(k,663) = rxt(k,70)
         mat(k,205) = rxt(k,71)
         mat(k,710) = rxt(k,72)
         mat(k,431) = rxt(k,73)
         mat(k,780) = rxt(k,105)
         mat(k,825) = rxt(k,106)
         mat(k,2192) = .500_r8*rxt(k,489)

         mat(k,454) = -( het_rates(k,218) )

         mat(k,768) = -( het_rates(k,219) )

         mat(k,1129) = -( het_rates(k,220) )
         mat(k,1061) = .150_r8*rxt(k,69)

         mat(k,1115) = -( het_rates(k,221) )

         mat(k,912) = -( het_rates(k,222) )

         mat(k,732) = -( het_rates(k,223) )

         mat(k,86) = -( het_rates(k,224) )

         mat(k,1173) = -( het_rates(k,225) )

         mat(k,748) = -( het_rates(k,226) )

         mat(k,92) = -( het_rates(k,227) )

         mat(k,499) = -( het_rates(k,228) )

         mat(k,2285) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,229) )
         mat(k,1442) = .050_r8*rxt(k,39)
         mat(k,134) = rxt(k,109)
         mat(k,1959) = rxt(k,481)

      end do

      end subroutine linmat03

      subroutine linmat( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer,  intent(in)    ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(in)    ::  het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)

      call linmat01( avec_len, mat, y, rxt, het_rates )
      call linmat02( avec_len, mat, y, rxt, het_rates )
      call linmat03( avec_len, mat, y, rxt, het_rates )

      end subroutine linmat

      end module mo_lin_matrix

      module mo_nln_matrix

      use shr_kind_mod, only : r8 => shr_kind_r8

      use chem_mods, only: veclen
      private
      public :: nlnmat

      contains

      subroutine     nlnmat01( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,666) = -(rxt(k,357)*y(k,217))
         mat(k,1631) = -rxt(k,357)*y(k,1)

         mat(k,1795) = rxt(k,360)*y(k,189)
         mat(k,896) = rxt(k,360)*y(k,124)

         mat(k,632) = -(rxt(k,361)*y(k,217))
         mat(k,1628) = -rxt(k,361)*y(k,2)

         mat(k,895) = rxt(k,358)*y(k,203)
         mat(k,1895) = rxt(k,358)*y(k,189)




         mat(k,999) = -(rxt(k,440)*y(k,126) + rxt(k,441)*y(k,134) + rxt(k,442) &
                      *y(k,217))
         mat(k,1719) = -rxt(k,440)*y(k,6)
         mat(k,2099) = -rxt(k,441)*y(k,6)
         mat(k,1661) = -rxt(k,442)*y(k,6)

         mat(k,164) = -(rxt(k,399)*y(k,217))
         mat(k,1559) = -rxt(k,399)*y(k,7)

         mat(k,421) = -(rxt(k,402)*y(k,217))
         mat(k,1600) = -rxt(k,402)*y(k,8)

         mat(k,481) = rxt(k,400)*y(k,203)
         mat(k,1880) = rxt(k,400)*y(k,191)


         mat(k,165) = .120_r8*rxt(k,399)*y(k,217)
         mat(k,1560) = .120_r8*rxt(k,399)*y(k,7)


         mat(k,991) = .100_r8*rxt(k,441)*y(k,134)
         mat(k,942) = .100_r8*rxt(k,444)*y(k,134)
         mat(k,2083) = .100_r8*rxt(k,441)*y(k,6) + .100_r8*rxt(k,444)*y(k,110)


         mat(k,1782) = .500_r8*rxt(k,401)*y(k,191) + .200_r8*rxt(k,428)*y(k,223)  &
                      + .060_r8*rxt(k,434)*y(k,226)
         mat(k,482) = .500_r8*rxt(k,401)*y(k,124)
         mat(k,728) = .200_r8*rxt(k,428)*y(k,124)
         mat(k,744) = .060_r8*rxt(k,434)*y(k,124)


         mat(k,1776) = .200_r8*rxt(k,428)*y(k,223) + .200_r8*rxt(k,434)*y(k,226)
         mat(k,727) = .200_r8*rxt(k,428)*y(k,124)
         mat(k,742) = .200_r8*rxt(k,434)*y(k,124)


         mat(k,1792) = .200_r8*rxt(k,428)*y(k,223) + .150_r8*rxt(k,434)*y(k,226)
         mat(k,729) = .200_r8*rxt(k,428)*y(k,124)
         mat(k,745) = .150_r8*rxt(k,434)*y(k,124)


         mat(k,1778) = .210_r8*rxt(k,434)*y(k,226)
         mat(k,743) = .210_r8*rxt(k,434)*y(k,124)

         mat(k,238) = -(rxt(k,362)*y(k,217))
         mat(k,1573) = -rxt(k,362)*y(k,15)

         mat(k,990) = .050_r8*rxt(k,441)*y(k,134)
         mat(k,941) = .050_r8*rxt(k,444)*y(k,134)
         mat(k,2082) = .050_r8*rxt(k,441)*y(k,6) + .050_r8*rxt(k,444)*y(k,110)

         mat(k,349) = -(rxt(k,328)*y(k,126) + rxt(k,329)*y(k,217))
         mat(k,1709) = -rxt(k,328)*y(k,16)
         mat(k,1590) = -rxt(k,329)*y(k,16)

         mat(k,1415) = -(rxt(k,212)*y(k,42) + rxt(k,213)*y(k,203) + rxt(k,214) &
                      *y(k,134))
         mat(k,1481) = -rxt(k,212)*y(k,17)
         mat(k,1941) = -rxt(k,213)*y(k,17)
         mat(k,2119) = -rxt(k,214)*y(k,17)

         mat(k,2211) = 4.000_r8*rxt(k,215)*y(k,19) + (rxt(k,216)+rxt(k,217))*y(k,59)  &
                      + rxt(k,220)*y(k,124) + rxt(k,223)*y(k,133) + rxt(k,470) &
                      *y(k,150) + rxt(k,224)*y(k,217)
         mat(k,143) = rxt(k,202)*y(k,216)
         mat(k,149) = rxt(k,228)*y(k,216)
         mat(k,468) = 2.000_r8*rxt(k,239)*y(k,56) + 2.000_r8*rxt(k,251)*y(k,216)  &
                      + 2.000_r8*rxt(k,240)*y(k,217)
         mat(k,595) = rxt(k,241)*y(k,56) + rxt(k,252)*y(k,216) + rxt(k,242)*y(k,217)
         mat(k,448) = 3.000_r8*rxt(k,246)*y(k,56) + 3.000_r8*rxt(k,229)*y(k,216)  &
                      + 3.000_r8*rxt(k,247)*y(k,217)
         mat(k,2006) = 2.000_r8*rxt(k,239)*y(k,41) + rxt(k,241)*y(k,43)  &
                      + 3.000_r8*rxt(k,246)*y(k,55)
         mat(k,1968) = (rxt(k,216)+rxt(k,217))*y(k,19)
         mat(k,107) = 2.000_r8*rxt(k,230)*y(k,216)
         mat(k,805) = rxt(k,225)*y(k,133) + rxt(k,231)*y(k,216) + rxt(k,226)*y(k,217)
         mat(k,1834) = rxt(k,220)*y(k,19)
         mat(k,2241) = rxt(k,223)*y(k,19) + rxt(k,225)*y(k,81)
         mat(k,1233) = rxt(k,470)*y(k,19)
         mat(k,1521) = rxt(k,202)*y(k,34) + rxt(k,228)*y(k,35) + 2.000_r8*rxt(k,251) &
                      *y(k,41) + rxt(k,252)*y(k,43) + 3.000_r8*rxt(k,229)*y(k,55)  &
                      + 2.000_r8*rxt(k,230)*y(k,78) + rxt(k,231)*y(k,81)
         mat(k,1685) = rxt(k,224)*y(k,19) + 2.000_r8*rxt(k,240)*y(k,41) + rxt(k,242) &
                      *y(k,43) + 3.000_r8*rxt(k,247)*y(k,55) + rxt(k,226)*y(k,81)


         mat(k,2205) = rxt(k,218)*y(k,59)
         mat(k,1962) = rxt(k,218)*y(k,19)
         mat(k,2139) = (rxt(k,531)+rxt(k,536))*y(k,91)
         mat(k,777) = (rxt(k,531)+rxt(k,536))*y(k,85)

         mat(k,2226) = -(4._r8*rxt(k,215)*y(k,19) + (rxt(k,216) + rxt(k,217) + rxt(k,218) &
                      ) * y(k,59) + rxt(k,219)*y(k,203) + rxt(k,220)*y(k,124) &
                      + rxt(k,221)*y(k,125) + rxt(k,223)*y(k,133) + rxt(k,224) &
                      *y(k,217) + rxt(k,470)*y(k,150))
         mat(k,1983) = -(rxt(k,216) + rxt(k,217) + rxt(k,218)) * y(k,19)
         mat(k,1957) = -rxt(k,219)*y(k,19)
         mat(k,1850) = -rxt(k,220)*y(k,19)
         mat(k,2202) = -rxt(k,221)*y(k,19)
         mat(k,2257) = -rxt(k,223)*y(k,19)
         mat(k,1701) = -rxt(k,224)*y(k,19)
         mat(k,1242) = -rxt(k,470)*y(k,19)

         mat(k,1422) = rxt(k,214)*y(k,134)
         mat(k,546) = rxt(k,222)*y(k,133)
         mat(k,809) = rxt(k,232)*y(k,216)
         mat(k,783) = rxt(k,227)*y(k,133)
         mat(k,2257) = mat(k,2257) + rxt(k,222)*y(k,20) + rxt(k,227)*y(k,91)
         mat(k,2135) = rxt(k,214)*y(k,17)
         mat(k,1537) = rxt(k,232)*y(k,81)

         mat(k,540) = -(rxt(k,222)*y(k,133))
         mat(k,2231) = -rxt(k,222)*y(k,20)

         mat(k,2207) = rxt(k,221)*y(k,125)
         mat(k,2169) = rxt(k,221)*y(k,19)


         mat(k,253) = -(rxt(k,403)*y(k,217))
         mat(k,1576) = -rxt(k,403)*y(k,22)

         mat(k,1774) = rxt(k,406)*y(k,193)
         mat(k,433) = rxt(k,406)*y(k,124)

         mat(k,339) = -(rxt(k,405)*y(k,217))
         mat(k,1588) = -rxt(k,405)*y(k,23)

         mat(k,434) = rxt(k,404)*y(k,203)
         mat(k,1872) = rxt(k,404)*y(k,193)

         mat(k,285) = -(rxt(k,277)*y(k,56) + rxt(k,278)*y(k,217))
         mat(k,1987) = -rxt(k,277)*y(k,24)
         mat(k,1580) = -rxt(k,278)*y(k,24)

         mat(k,548) = -(rxt(k,279)*y(k,56) + rxt(k,280)*y(k,134) + rxt(k,305)*y(k,217))
         mat(k,1992) = -rxt(k,279)*y(k,25)
         mat(k,2086) = -rxt(k,280)*y(k,25)
         mat(k,1617) = -rxt(k,305)*y(k,25)

         mat(k,261) = -(rxt(k,285)*y(k,217))
         mat(k,1578) = -rxt(k,285)*y(k,26)

         mat(k,812) = .800_r8*rxt(k,281)*y(k,194) + .200_r8*rxt(k,282)*y(k,198)
         mat(k,2025) = .200_r8*rxt(k,282)*y(k,194)

         mat(k,344) = -(rxt(k,286)*y(k,217))
         mat(k,1589) = -rxt(k,286)*y(k,27)

         mat(k,813) = rxt(k,283)*y(k,203)
         mat(k,1873) = rxt(k,283)*y(k,194)

         mat(k,291) = -(rxt(k,287)*y(k,56) + rxt(k,288)*y(k,217))
         mat(k,1988) = -rxt(k,287)*y(k,28)
         mat(k,1581) = -rxt(k,288)*y(k,28)

         mat(k,1024) = -(rxt(k,308)*y(k,126) + rxt(k,309)*y(k,134) + rxt(k,326) &
                      *y(k,217))
         mat(k,1720) = -rxt(k,308)*y(k,29)
         mat(k,2100) = -rxt(k,309)*y(k,29)
         mat(k,1662) = -rxt(k,326)*y(k,29)

         mat(k,843) = .130_r8*rxt(k,386)*y(k,134)
         mat(k,2100) = mat(k,2100) + .130_r8*rxt(k,386)*y(k,98)

         mat(k,415) = -(rxt(k,313)*y(k,217))
         mat(k,1599) = -rxt(k,313)*y(k,30)

         mat(k,790) = rxt(k,311)*y(k,203)
         mat(k,1879) = rxt(k,311)*y(k,195)

         mat(k,109) = -(rxt(k,314)*y(k,217))
         mat(k,1556) = -rxt(k,314)*y(k,31)

         mat(k,265) = -(rxt(k,409)*y(k,217))
         mat(k,1579) = -rxt(k,409)*y(k,32)

         mat(k,623) = rxt(k,407)*y(k,203)
         mat(k,1867) = rxt(k,407)*y(k,196)

         mat(k,99) = -(rxt(k,201)*y(k,216))
         mat(k,1499) = -rxt(k,201)*y(k,33)

         mat(k,141) = -(rxt(k,202)*y(k,216))
         mat(k,1504) = -rxt(k,202)*y(k,34)

         mat(k,146) = -(rxt(k,228)*y(k,216))
         mat(k,1505) = -rxt(k,228)*y(k,35)

         mat(k,113) = -(rxt(k,203)*y(k,216))
         mat(k,1501) = -rxt(k,203)*y(k,36)

         mat(k,151) = -(rxt(k,204)*y(k,216))
         mat(k,1506) = -rxt(k,204)*y(k,37)

         mat(k,117) = -(rxt(k,205)*y(k,216))
         mat(k,1502) = -rxt(k,205)*y(k,38)

         mat(k,156) = -(rxt(k,206)*y(k,216))
         mat(k,1507) = -rxt(k,206)*y(k,39)

         mat(k,121) = -(rxt(k,207)*y(k,216))
         mat(k,1503) = -rxt(k,207)*y(k,40)

         mat(k,467) = -(rxt(k,239)*y(k,56) + rxt(k,240)*y(k,217) + rxt(k,251)*y(k,216))
         mat(k,1991) = -rxt(k,239)*y(k,41)
         mat(k,1607) = -rxt(k,240)*y(k,41)
         mat(k,1516) = -rxt(k,251)*y(k,41)

         mat(k,1485) = -(rxt(k,176)*y(k,56) + rxt(k,212)*y(k,17) + rxt(k,256)*y(k,203) &
                      + rxt(k,257)*y(k,126) + rxt(k,258)*y(k,133) + rxt(k,259) &
                      *y(k,217))
         mat(k,2010) = -rxt(k,176)*y(k,42)
         mat(k,1417) = -rxt(k,212)*y(k,42)
         mat(k,1945) = -rxt(k,256)*y(k,42)
         mat(k,1746) = -rxt(k,257)*y(k,42)
         mat(k,2245) = -rxt(k,258)*y(k,42)
         mat(k,1689) = -rxt(k,259)*y(k,42)

         mat(k,672) = .400_r8*rxt(k,357)*y(k,217)
         mat(k,1009) = .340_r8*rxt(k,441)*y(k,134)
         mat(k,353) = .500_r8*rxt(k,328)*y(k,126)
         mat(k,552) = rxt(k,280)*y(k,134)
         mat(k,1031) = .500_r8*rxt(k,309)*y(k,134)
         mat(k,613) = .500_r8*rxt(k,297)*y(k,217)
         mat(k,787) = rxt(k,264)*y(k,217)
         mat(k,387) = .300_r8*rxt(k,265)*y(k,217)
         mat(k,1433) = (rxt(k,273)+rxt(k,274))*y(k,216)
         mat(k,1971) = rxt(k,183)*y(k,198)
         mat(k,1098) = .800_r8*rxt(k,302)*y(k,217)
         mat(k,851) = .910_r8*rxt(k,386)*y(k,134)
         mat(k,586) = .300_r8*rxt(k,377)*y(k,217)
         mat(k,1199) = .800_r8*rxt(k,381)*y(k,198)
         mat(k,1214) = .120_r8*rxt(k,339)*y(k,134)
         mat(k,576) = .500_r8*rxt(k,352)*y(k,217)
         mat(k,959) = .340_r8*rxt(k,444)*y(k,134)
         mat(k,1337) = .600_r8*rxt(k,353)*y(k,134)
         mat(k,1838) = .100_r8*rxt(k,359)*y(k,189) + rxt(k,263)*y(k,198)  &
                      + .500_r8*rxt(k,330)*y(k,200) + .500_r8*rxt(k,299)*y(k,202)  &
                      + .920_r8*rxt(k,369)*y(k,205) + .250_r8*rxt(k,337)*y(k,209)  &
                      + rxt(k,346)*y(k,211) + rxt(k,320)*y(k,219) + rxt(k,324) &
                      *y(k,220) + .340_r8*rxt(k,453)*y(k,221) + .320_r8*rxt(k,458) &
                      *y(k,222) + .250_r8*rxt(k,394)*y(k,225)
         mat(k,1746) = mat(k,1746) + .500_r8*rxt(k,328)*y(k,16) + rxt(k,370)*y(k,205)  &
                      + .250_r8*rxt(k,336)*y(k,209) + rxt(k,347)*y(k,211)
         mat(k,2123) = .340_r8*rxt(k,441)*y(k,6) + rxt(k,280)*y(k,25)  &
                      + .500_r8*rxt(k,309)*y(k,29) + .910_r8*rxt(k,386)*y(k,98)  &
                      + .120_r8*rxt(k,339)*y(k,105) + .340_r8*rxt(k,444)*y(k,110)  &
                      + .600_r8*rxt(k,353)*y(k,111)
         mat(k,527) = rxt(k,304)*y(k,217)
         mat(k,1063) = .680_r8*rxt(k,462)*y(k,217)
         mat(k,903) = .100_r8*rxt(k,359)*y(k,124)
         mat(k,817) = .700_r8*rxt(k,282)*y(k,198)
         mat(k,794) = rxt(k,310)*y(k,198)
         mat(k,1389) = rxt(k,293)*y(k,198) + rxt(k,366)*y(k,205) + .250_r8*rxt(k,333) &
                      *y(k,209) + rxt(k,342)*y(k,211) + .250_r8*rxt(k,391)*y(k,225)
         mat(k,2062) = rxt(k,183)*y(k,59) + .800_r8*rxt(k,381)*y(k,101) + rxt(k,263) &
                      *y(k,124) + .700_r8*rxt(k,282)*y(k,194) + rxt(k,310)*y(k,195)  &
                      + rxt(k,293)*y(k,197) + (4.000_r8*rxt(k,260)+2.000_r8*rxt(k,261)) &
                      *y(k,198) + 1.500_r8*rxt(k,367)*y(k,205) + .750_r8*rxt(k,372) &
                      *y(k,206) + .880_r8*rxt(k,334)*y(k,209) + 2.000_r8*rxt(k,343) &
                      *y(k,211) + .750_r8*rxt(k,446)*y(k,215) + .800_r8*rxt(k,322) &
                      *y(k,220) + .930_r8*rxt(k,451)*y(k,221) + .950_r8*rxt(k,456) &
                      *y(k,222) + .800_r8*rxt(k,392)*y(k,225)
         mat(k,568) = .500_r8*rxt(k,330)*y(k,124)
         mat(k,716) = .500_r8*rxt(k,299)*y(k,124)
         mat(k,1945) = mat(k,1945) + .450_r8*rxt(k,344)*y(k,211) + .150_r8*rxt(k,323) &
                      *y(k,220)
         mat(k,1262) = .920_r8*rxt(k,369)*y(k,124) + rxt(k,370)*y(k,126) + rxt(k,366) &
                      *y(k,197) + 1.500_r8*rxt(k,367)*y(k,198)
         mat(k,1294) = .750_r8*rxt(k,372)*y(k,198)
         mat(k,1315) = .250_r8*rxt(k,337)*y(k,124) + .250_r8*rxt(k,336)*y(k,126)  &
                      + .250_r8*rxt(k,333)*y(k,197) + .880_r8*rxt(k,334)*y(k,198)
         mat(k,1357) = rxt(k,346)*y(k,124) + rxt(k,347)*y(k,126) + rxt(k,342)*y(k,197)  &
                      + 2.000_r8*rxt(k,343)*y(k,198) + .450_r8*rxt(k,344)*y(k,203)  &
                      + 4.000_r8*rxt(k,345)*y(k,211)
         mat(k,1050) = .750_r8*rxt(k,446)*y(k,198)
         mat(k,1525) = (rxt(k,273)+rxt(k,274))*y(k,54)
         mat(k,1689) = mat(k,1689) + .400_r8*rxt(k,357)*y(k,1) + .500_r8*rxt(k,297) &
                      *y(k,51) + rxt(k,264)*y(k,52) + .300_r8*rxt(k,265)*y(k,53)  &
                      + .800_r8*rxt(k,302)*y(k,74) + .300_r8*rxt(k,377)*y(k,99)  &
                      + .500_r8*rxt(k,352)*y(k,109) + rxt(k,304)*y(k,139)  &
                      + .680_r8*rxt(k,462)*y(k,178)
         mat(k,771) = rxt(k,320)*y(k,124)
         mat(k,1133) = rxt(k,324)*y(k,124) + .800_r8*rxt(k,322)*y(k,198)  &
                      + .150_r8*rxt(k,323)*y(k,203)
         mat(k,1119) = .340_r8*rxt(k,453)*y(k,124) + .930_r8*rxt(k,451)*y(k,198)
         mat(k,916) = .320_r8*rxt(k,458)*y(k,124) + .950_r8*rxt(k,456)*y(k,198)
         mat(k,1176) = .250_r8*rxt(k,394)*y(k,124) + .250_r8*rxt(k,391)*y(k,197)  &
                      + .800_r8*rxt(k,392)*y(k,198)

      end do

      end subroutine     nlnmat01

      subroutine     nlnmat02( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,594) = -(rxt(k,241)*y(k,56) + rxt(k,242)*y(k,217) + rxt(k,252)*y(k,216))
         mat(k,1993) = -rxt(k,241)*y(k,43)
         mat(k,1623) = -rxt(k,242)*y(k,43)
         mat(k,1517) = -rxt(k,252)*y(k,43)

         mat(k,125) = -(rxt(k,243)*y(k,217))
         mat(k,1557) = -rxt(k,243)*y(k,44)

         mat(k,1085) = -(rxt(k,289)*y(k,126) + rxt(k,290)*y(k,217))
         mat(k,1724) = -rxt(k,289)*y(k,45)
         mat(k,1666) = -rxt(k,290)*y(k,45)

         mat(k,670) = .800_r8*rxt(k,357)*y(k,217)
         mat(k,352) = rxt(k,328)*y(k,126)
         mat(k,262) = rxt(k,285)*y(k,217)
         mat(k,346) = .500_r8*rxt(k,286)*y(k,217)
         mat(k,1025) = .500_r8*rxt(k,309)*y(k,134)
         mat(k,1327) = .100_r8*rxt(k,353)*y(k,134)
         mat(k,1817) = .400_r8*rxt(k,359)*y(k,189) + rxt(k,284)*y(k,194)  &
                      + .270_r8*rxt(k,312)*y(k,195) + rxt(k,330)*y(k,200) + rxt(k,349) &
                      *y(k,213) + rxt(k,320)*y(k,219)
         mat(k,1724) = mat(k,1724) + rxt(k,328)*y(k,16)
         mat(k,2103) = .500_r8*rxt(k,309)*y(k,29) + .100_r8*rxt(k,353)*y(k,111)
         mat(k,901) = .400_r8*rxt(k,359)*y(k,124)
         mat(k,816) = rxt(k,284)*y(k,124) + 3.200_r8*rxt(k,281)*y(k,194)  &
                      + .800_r8*rxt(k,282)*y(k,198)
         mat(k,793) = .270_r8*rxt(k,312)*y(k,124)
         mat(k,2043) = .800_r8*rxt(k,282)*y(k,194)
         mat(k,566) = rxt(k,330)*y(k,124)
         mat(k,1924) = .200_r8*rxt(k,348)*y(k,213)
         mat(k,678) = rxt(k,349)*y(k,124) + .200_r8*rxt(k,348)*y(k,203)
         mat(k,1666) = mat(k,1666) + .800_r8*rxt(k,357)*y(k,1) + rxt(k,285)*y(k,26)  &
                      + .500_r8*rxt(k,286)*y(k,27)
         mat(k,769) = rxt(k,320)*y(k,124)

         mat(k,365) = -(rxt(k,244)*y(k,56) + rxt(k,245)*y(k,217))
         mat(k,1989) = -rxt(k,244)*y(k,46)
         mat(k,1592) = -rxt(k,245)*y(k,46)

         mat(k,102) = -(rxt(k,291)*y(k,217))
         mat(k,1555) = -rxt(k,291)*y(k,47)

         mat(k,922) = -(rxt(k,327)*y(k,217))
         mat(k,1656) = -rxt(k,327)*y(k,48)

         mat(k,669) = .800_r8*rxt(k,357)*y(k,217)
         mat(k,995) = .520_r8*rxt(k,441)*y(k,134)
         mat(k,351) = .500_r8*rxt(k,328)*y(k,126)
         mat(k,946) = .520_r8*rxt(k,444)*y(k,134)
         mat(k,1810) = .250_r8*rxt(k,359)*y(k,189) + .820_r8*rxt(k,312)*y(k,195)  &
                      + .500_r8*rxt(k,330)*y(k,200) + .270_r8*rxt(k,453)*y(k,221)  &
                      + .040_r8*rxt(k,458)*y(k,222)
         mat(k,1714) = .500_r8*rxt(k,328)*y(k,16)
         mat(k,2094) = .520_r8*rxt(k,441)*y(k,6) + .520_r8*rxt(k,444)*y(k,110)
         mat(k,1058) = .500_r8*rxt(k,462)*y(k,217)
         mat(k,900) = .250_r8*rxt(k,359)*y(k,124)
         mat(k,792) = .820_r8*rxt(k,312)*y(k,124) + .820_r8*rxt(k,310)*y(k,198)
         mat(k,2037) = .820_r8*rxt(k,310)*y(k,195) + .150_r8*rxt(k,451)*y(k,221)  &
                      + .025_r8*rxt(k,456)*y(k,222)
         mat(k,565) = .500_r8*rxt(k,330)*y(k,124)
         mat(k,1656) = mat(k,1656) + .800_r8*rxt(k,357)*y(k,1) + .500_r8*rxt(k,462) &
                      *y(k,178)
         mat(k,1111) = .270_r8*rxt(k,453)*y(k,124) + .150_r8*rxt(k,451)*y(k,198)
         mat(k,913) = .040_r8*rxt(k,458)*y(k,124) + .025_r8*rxt(k,456)*y(k,198)

         mat(k,1221) = -(rxt(k,315)*y(k,126) + rxt(k,316)*y(k,217))
         mat(k,1734) = -rxt(k,315)*y(k,49)
         mat(k,1676) = -rxt(k,316)*y(k,49)

         mat(k,1141) = rxt(k,317)*y(k,217)
         mat(k,1210) = .880_r8*rxt(k,339)*y(k,134)
         mat(k,1330) = .500_r8*rxt(k,353)*y(k,134)
         mat(k,1827) = .170_r8*rxt(k,412)*y(k,199) + .050_r8*rxt(k,375)*y(k,206)  &
                      + .250_r8*rxt(k,337)*y(k,209) + .170_r8*rxt(k,418)*y(k,212)  &
                      + .400_r8*rxt(k,428)*y(k,223) + .250_r8*rxt(k,394)*y(k,225)  &
                      + .540_r8*rxt(k,434)*y(k,226) + .510_r8*rxt(k,437)*y(k,228)
         mat(k,1734) = mat(k,1734) + .050_r8*rxt(k,376)*y(k,206) + .250_r8*rxt(k,336) &
                      *y(k,209) + .250_r8*rxt(k,395)*y(k,225)
         mat(k,858) = rxt(k,318)*y(k,217)
         mat(k,2111) = .880_r8*rxt(k,339)*y(k,105) + .500_r8*rxt(k,353)*y(k,111)
         mat(k,1380) = .250_r8*rxt(k,333)*y(k,209) + .250_r8*rxt(k,391)*y(k,225)
         mat(k,2052) = .240_r8*rxt(k,334)*y(k,209) + .500_r8*rxt(k,322)*y(k,220)  &
                      + .100_r8*rxt(k,392)*y(k,225)
         mat(k,761) = .170_r8*rxt(k,412)*y(k,124) + .070_r8*rxt(k,411)*y(k,203)
         mat(k,1933) = .070_r8*rxt(k,411)*y(k,199) + .070_r8*rxt(k,417)*y(k,212)
         mat(k,1287) = .050_r8*rxt(k,375)*y(k,124) + .050_r8*rxt(k,376)*y(k,126)
         mat(k,1310) = .250_r8*rxt(k,337)*y(k,124) + .250_r8*rxt(k,336)*y(k,126)  &
                      + .250_r8*rxt(k,333)*y(k,197) + .240_r8*rxt(k,334)*y(k,198)
         mat(k,876) = .170_r8*rxt(k,418)*y(k,124) + .070_r8*rxt(k,417)*y(k,203)
         mat(k,1676) = mat(k,1676) + rxt(k,317)*y(k,95) + rxt(k,318)*y(k,127)
         mat(k,1131) = .500_r8*rxt(k,322)*y(k,198)
         mat(k,737) = .400_r8*rxt(k,428)*y(k,124)
         mat(k,1174) = .250_r8*rxt(k,394)*y(k,124) + .250_r8*rxt(k,395)*y(k,126)  &
                      + .250_r8*rxt(k,391)*y(k,197) + .100_r8*rxt(k,392)*y(k,198)
         mat(k,753) = .540_r8*rxt(k,434)*y(k,124)
         mat(k,501) = .510_r8*rxt(k,437)*y(k,124)

         mat(k,684) = -(rxt(k,296)*y(k,217))
         mat(k,1633) = -rxt(k,296)*y(k,50)

         mat(k,1019) = .120_r8*rxt(k,309)*y(k,134)
         mat(k,2088) = .120_r8*rxt(k,309)*y(k,29)
         mat(k,1370) = .100_r8*rxt(k,293)*y(k,198) + .150_r8*rxt(k,294)*y(k,203)
         mat(k,2030) = .100_r8*rxt(k,293)*y(k,197)
         mat(k,1899) = .150_r8*rxt(k,294)*y(k,197) + .150_r8*rxt(k,344)*y(k,211)
         mat(k,1349) = .150_r8*rxt(k,344)*y(k,203)

         mat(k,610) = -(rxt(k,297)*y(k,217))
         mat(k,1625) = -rxt(k,297)*y(k,51)

         mat(k,1369) = .400_r8*rxt(k,294)*y(k,203)
         mat(k,1893) = .400_r8*rxt(k,294)*y(k,197) + .400_r8*rxt(k,344)*y(k,211)
         mat(k,1348) = .400_r8*rxt(k,344)*y(k,203)

         mat(k,786) = -(rxt(k,264)*y(k,217))
         mat(k,1642) = -rxt(k,264)*y(k,52)

         mat(k,1187) = .200_r8*rxt(k,381)*y(k,198)
         mat(k,814) = .300_r8*rxt(k,282)*y(k,198)
         mat(k,2031) = .200_r8*rxt(k,381)*y(k,101) + .300_r8*rxt(k,282)*y(k,194)  &
                      + 2.000_r8*rxt(k,261)*y(k,198) + .250_r8*rxt(k,367)*y(k,205)  &
                      + .250_r8*rxt(k,372)*y(k,206) + .250_r8*rxt(k,334)*y(k,209)  &
                      + .250_r8*rxt(k,446)*y(k,215) + .500_r8*rxt(k,322)*y(k,220)  &
                      + .250_r8*rxt(k,451)*y(k,221) + .250_r8*rxt(k,456)*y(k,222)  &
                      + .300_r8*rxt(k,392)*y(k,225)
         mat(k,1247) = .250_r8*rxt(k,367)*y(k,198)
         mat(k,1277) = .250_r8*rxt(k,372)*y(k,198)
         mat(k,1305) = .250_r8*rxt(k,334)*y(k,198)
         mat(k,1043) = .250_r8*rxt(k,446)*y(k,198)
         mat(k,1128) = .500_r8*rxt(k,322)*y(k,198)
         mat(k,1109) = .250_r8*rxt(k,451)*y(k,198)
         mat(k,911) = .250_r8*rxt(k,456)*y(k,198)
         mat(k,1167) = .300_r8*rxt(k,392)*y(k,198)

         mat(k,385) = -(rxt(k,265)*y(k,217))
         mat(k,1594) = -rxt(k,265)*y(k,53)

         mat(k,2028) = rxt(k,262)*y(k,203)
         mat(k,1874) = rxt(k,262)*y(k,198)

         mat(k,1430) = -(rxt(k,177)*y(k,56) + rxt(k,233)*y(k,73) + rxt(k,266)*y(k,217) &
                      + (rxt(k,272) + rxt(k,273) + rxt(k,274)) * y(k,216))
         mat(k,2007) = -rxt(k,177)*y(k,54)
         mat(k,866) = -rxt(k,233)*y(k,54)
         mat(k,1686) = -rxt(k,266)*y(k,54)
         mat(k,1522) = -(rxt(k,272) + rxt(k,273) + rxt(k,274)) * y(k,54)

         mat(k,1030) = .100_r8*rxt(k,309)*y(k,134)
         mat(k,2120) = .100_r8*rxt(k,309)*y(k,29)

         mat(k,447) = -(rxt(k,229)*y(k,216) + rxt(k,246)*y(k,56) + rxt(k,247)*y(k,217))
         mat(k,1515) = -rxt(k,229)*y(k,55)
         mat(k,1990) = -rxt(k,246)*y(k,55)
         mat(k,1603) = -rxt(k,247)*y(k,55)

         mat(k,2017) = -(rxt(k,176)*y(k,42) + rxt(k,177)*y(k,54) + rxt(k,178)*y(k,77) &
                      + rxt(k,179)*y(k,79) + (rxt(k,180) + rxt(k,181)) * y(k,203) &
                      + rxt(k,182)*y(k,134) + rxt(k,189)*y(k,60) + rxt(k,198)*y(k,92) &
                      + rxt(k,239)*y(k,41) + rxt(k,241)*y(k,43) + rxt(k,244)*y(k,46) &
                      + rxt(k,246)*y(k,55) + rxt(k,287)*y(k,28))
         mat(k,1491) = -rxt(k,176)*y(k,56)
         mat(k,1438) = -rxt(k,177)*y(k,56)
         mat(k,1408) = -rxt(k,178)*y(k,56)
         mat(k,606) = -rxt(k,179)*y(k,56)
         mat(k,1952) = -(rxt(k,180) + rxt(k,181)) * y(k,56)
         mat(k,2130) = -rxt(k,182)*y(k,56)
         mat(k,889) = -rxt(k,189)*y(k,56)
         mat(k,827) = -rxt(k,198)*y(k,56)
         mat(k,471) = -rxt(k,239)*y(k,56)
         mat(k,599) = -rxt(k,241)*y(k,56)
         mat(k,369) = -rxt(k,244)*y(k,56)
         mat(k,451) = -rxt(k,246)*y(k,56)
         mat(k,294) = -rxt(k,287)*y(k,56)

         mat(k,2221) = rxt(k,217)*y(k,59)
         mat(k,101) = 4.000_r8*rxt(k,201)*y(k,216)
         mat(k,145) = rxt(k,202)*y(k,216)
         mat(k,116) = 2.000_r8*rxt(k,203)*y(k,216)
         mat(k,155) = 2.000_r8*rxt(k,204)*y(k,216)
         mat(k,120) = 2.000_r8*rxt(k,205)*y(k,216)
         mat(k,160) = rxt(k,206)*y(k,216)
         mat(k,124) = 2.000_r8*rxt(k,207)*y(k,216)
         mat(k,127) = 3.000_r8*rxt(k,243)*y(k,217)
         mat(k,369) = mat(k,369) + rxt(k,245)*y(k,217)
         mat(k,1978) = rxt(k,217)*y(k,19) + (4.000_r8*rxt(k,184)+2.000_r8*rxt(k,186)) &
                      *y(k,59) + rxt(k,188)*y(k,124) + rxt(k,193)*y(k,133)  &
                      + rxt(k,471)*y(k,150) + rxt(k,183)*y(k,198) + rxt(k,194) &
                      *y(k,217)
         mat(k,229) = rxt(k,238)*y(k,216)
         mat(k,225) = rxt(k,253)*y(k,216) + rxt(k,248)*y(k,217)
         mat(k,252) = rxt(k,254)*y(k,216) + rxt(k,249)*y(k,217)
         mat(k,308) = rxt(k,255)*y(k,216) + rxt(k,250)*y(k,217)
         mat(k,2153) = rxt(k,196)*y(k,133) + rxt(k,208)*y(k,216) + rxt(k,197)*y(k,217)
         mat(k,1845) = rxt(k,188)*y(k,59)
         mat(k,2252) = rxt(k,193)*y(k,59) + rxt(k,196)*y(k,85)
         mat(k,1239) = rxt(k,471)*y(k,59)
         mat(k,2069) = rxt(k,183)*y(k,59)
         mat(k,1532) = 4.000_r8*rxt(k,201)*y(k,33) + rxt(k,202)*y(k,34)  &
                      + 2.000_r8*rxt(k,203)*y(k,36) + 2.000_r8*rxt(k,204)*y(k,37)  &
                      + 2.000_r8*rxt(k,205)*y(k,38) + rxt(k,206)*y(k,39)  &
                      + 2.000_r8*rxt(k,207)*y(k,40) + rxt(k,238)*y(k,65) + rxt(k,253) &
                      *y(k,82) + rxt(k,254)*y(k,83) + rxt(k,255)*y(k,84) + rxt(k,208) &
                      *y(k,85)
         mat(k,1696) = 3.000_r8*rxt(k,243)*y(k,44) + rxt(k,245)*y(k,46) + rxt(k,194) &
                      *y(k,59) + rxt(k,248)*y(k,82) + rxt(k,249)*y(k,83) + rxt(k,250) &
                      *y(k,84) + rxt(k,197)*y(k,85)


         mat(k,1986) = rxt(k,189)*y(k,60)
         mat(k,1961) = 2.000_r8*rxt(k,185)*y(k,59)
         mat(k,882) = rxt(k,189)*y(k,56) + (rxt(k,529)+rxt(k,534)+rxt(k,539))*y(k,85)
         mat(k,2138) = (rxt(k,529)+rxt(k,534)+rxt(k,539))*y(k,60) + (rxt(k,524) &
                       +rxt(k,530)+rxt(k,535))*y(k,92)
         mat(k,823) = (rxt(k,524)+rxt(k,530)+rxt(k,535))*y(k,85)


         mat(k,1960) = 2.000_r8*rxt(k,210)*y(k,59)

         mat(k,1977) = -(rxt(k,183)*y(k,198) + (4._r8*rxt(k,184) + 4._r8*rxt(k,185) &
                      + 4._r8*rxt(k,186) + 4._r8*rxt(k,210)) * y(k,59) + rxt(k,187) &
                      *y(k,203) + rxt(k,188)*y(k,124) + rxt(k,190)*y(k,125) + rxt(k,193) &
                      *y(k,133) + (rxt(k,194) + rxt(k,195)) * y(k,217) + (rxt(k,216) &
                      + rxt(k,217) + rxt(k,218)) * y(k,19) + rxt(k,471)*y(k,150))
         mat(k,2068) = -rxt(k,183)*y(k,59)
         mat(k,1951) = -rxt(k,187)*y(k,59)
         mat(k,1844) = -rxt(k,188)*y(k,59)
         mat(k,2196) = -rxt(k,190)*y(k,59)
         mat(k,2251) = -rxt(k,193)*y(k,59)
         mat(k,1695) = -(rxt(k,194) + rxt(k,195)) * y(k,59)
         mat(k,2220) = -(rxt(k,216) + rxt(k,217) + rxt(k,218)) * y(k,59)
         mat(k,1238) = -rxt(k,471)*y(k,59)

         mat(k,2016) = rxt(k,198)*y(k,92) + rxt(k,182)*y(k,134) + rxt(k,181)*y(k,203)
         mat(k,888) = rxt(k,191)*y(k,133)
         mat(k,2152) = rxt(k,209)*y(k,216)
         mat(k,826) = rxt(k,198)*y(k,56) + rxt(k,199)*y(k,133) + rxt(k,200)*y(k,217)
         mat(k,2251) = mat(k,2251) + rxt(k,191)*y(k,60) + rxt(k,199)*y(k,92)
         mat(k,2129) = rxt(k,182)*y(k,56)
         mat(k,331) = rxt(k,476)*y(k,150)
         mat(k,1238) = mat(k,1238) + rxt(k,476)*y(k,136)
         mat(k,1951) = mat(k,1951) + rxt(k,181)*y(k,56)
         mat(k,1531) = rxt(k,209)*y(k,85)
         mat(k,1695) = mat(k,1695) + rxt(k,200)*y(k,92)

      end do

      end subroutine     nlnmat02

      subroutine     nlnmat03( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,884) = -(rxt(k,189)*y(k,56) + rxt(k,191)*y(k,133) + rxt(k,192)*y(k,217) &
                      + (rxt(k,529) + rxt(k,534) + rxt(k,539)) * y(k,85))
         mat(k,1998) = -rxt(k,189)*y(k,60)
         mat(k,2237) = -rxt(k,191)*y(k,60)
         mat(k,1653) = -rxt(k,192)*y(k,60)
         mat(k,2142) = -(rxt(k,529) + rxt(k,534) + rxt(k,539)) * y(k,60)

         mat(k,1966) = rxt(k,190)*y(k,125)
         mat(k,2178) = rxt(k,190)*y(k,59)


         mat(k,1103) = -(rxt(k,276)*y(k,217))
         mat(k,1668) = -rxt(k,276)*y(k,62)

         mat(k,1003) = .230_r8*rxt(k,441)*y(k,134)
         mat(k,1414) = rxt(k,212)*y(k,42)
         mat(k,288) = .350_r8*rxt(k,278)*y(k,217)
         mat(k,551) = .630_r8*rxt(k,280)*y(k,134)
         mat(k,1026) = .560_r8*rxt(k,309)*y(k,134)
         mat(k,1479) = rxt(k,212)*y(k,17) + rxt(k,176)*y(k,56) + rxt(k,257)*y(k,126)  &
                      + rxt(k,258)*y(k,133) + rxt(k,259)*y(k,217)
         mat(k,366) = rxt(k,244)*y(k,56)
         mat(k,1220) = rxt(k,315)*y(k,126) + rxt(k,316)*y(k,217)
         mat(k,2003) = rxt(k,176)*y(k,42) + rxt(k,244)*y(k,46)
         mat(k,980) = rxt(k,303)*y(k,217)
         mat(k,844) = .620_r8*rxt(k,386)*y(k,134)
         mat(k,1208) = .650_r8*rxt(k,339)*y(k,134)
         mat(k,954) = .230_r8*rxt(k,444)*y(k,134)
         mat(k,1328) = .560_r8*rxt(k,353)*y(k,134)
         mat(k,1819) = .170_r8*rxt(k,412)*y(k,199) + .220_r8*rxt(k,337)*y(k,209)  &
                      + .400_r8*rxt(k,415)*y(k,210) + .350_r8*rxt(k,418)*y(k,212)  &
                      + .225_r8*rxt(k,453)*y(k,221) + .250_r8*rxt(k,394)*y(k,225)
         mat(k,1726) = rxt(k,257)*y(k,42) + rxt(k,315)*y(k,49) + .220_r8*rxt(k,336) &
                      *y(k,209) + .500_r8*rxt(k,395)*y(k,225)
         mat(k,2238) = rxt(k,258)*y(k,42) + rxt(k,465)*y(k,137)
         mat(k,2105) = .230_r8*rxt(k,441)*y(k,6) + .630_r8*rxt(k,280)*y(k,25)  &
                      + .560_r8*rxt(k,309)*y(k,29) + .620_r8*rxt(k,386)*y(k,98)  &
                      + .650_r8*rxt(k,339)*y(k,105) + .230_r8*rxt(k,444)*y(k,110)  &
                      + .560_r8*rxt(k,353)*y(k,111)
         mat(k,360) = rxt(k,465)*y(k,133) + rxt(k,466)*y(k,217)
         mat(k,1060) = .700_r8*rxt(k,462)*y(k,217)
         mat(k,1375) = .220_r8*rxt(k,333)*y(k,209) + .250_r8*rxt(k,391)*y(k,225)
         mat(k,2045) = .110_r8*rxt(k,334)*y(k,209) + .125_r8*rxt(k,451)*y(k,221)  &
                      + .200_r8*rxt(k,392)*y(k,225)
         mat(k,760) = .170_r8*rxt(k,412)*y(k,124) + .070_r8*rxt(k,411)*y(k,203)
         mat(k,1926) = .070_r8*rxt(k,411)*y(k,199) + .160_r8*rxt(k,414)*y(k,210)  &
                      + .140_r8*rxt(k,417)*y(k,212)
         mat(k,1307) = .220_r8*rxt(k,337)*y(k,124) + .220_r8*rxt(k,336)*y(k,126)  &
                      + .220_r8*rxt(k,333)*y(k,197) + .110_r8*rxt(k,334)*y(k,198)
         mat(k,723) = .400_r8*rxt(k,415)*y(k,124) + .160_r8*rxt(k,414)*y(k,203)
         mat(k,875) = .350_r8*rxt(k,418)*y(k,124) + .140_r8*rxt(k,417)*y(k,203)
         mat(k,1668) = mat(k,1668) + .350_r8*rxt(k,278)*y(k,24) + rxt(k,259)*y(k,42)  &
                      + rxt(k,316)*y(k,49) + rxt(k,303)*y(k,75) + rxt(k,466)*y(k,137)  &
                      + .700_r8*rxt(k,462)*y(k,178)
         mat(k,1114) = .225_r8*rxt(k,453)*y(k,124) + .125_r8*rxt(k,451)*y(k,198)
         mat(k,1171) = .250_r8*rxt(k,394)*y(k,124) + .500_r8*rxt(k,395)*y(k,126)  &
                      + .250_r8*rxt(k,391)*y(k,197) + .200_r8*rxt(k,392)*y(k,198)


         mat(k,992) = .270_r8*rxt(k,441)*y(k,134)
         mat(k,1021) = .200_r8*rxt(k,309)*y(k,134)
         mat(k,685) = rxt(k,296)*y(k,217)
         mat(k,611) = .500_r8*rxt(k,297)*y(k,217)
         mat(k,1102) = rxt(k,276)*y(k,217)
         mat(k,1094) = .800_r8*rxt(k,302)*y(k,217)
         mat(k,978) = rxt(k,303)*y(k,217)
         mat(k,928) = rxt(k,268)*y(k,217)
         mat(k,573) = .500_r8*rxt(k,352)*y(k,217)
         mat(k,943) = .270_r8*rxt(k,444)*y(k,134)
         mat(k,1324) = .100_r8*rxt(k,353)*y(k,134)
         mat(k,1804) = rxt(k,295)*y(k,197) + .900_r8*rxt(k,453)*y(k,221)
         mat(k,2090) = .270_r8*rxt(k,441)*y(k,6) + .200_r8*rxt(k,309)*y(k,29)  &
                      + .270_r8*rxt(k,444)*y(k,110) + .100_r8*rxt(k,353)*y(k,111)
         mat(k,1057) = 1.800_r8*rxt(k,462)*y(k,217)
         mat(k,1371) = rxt(k,295)*y(k,124) + 4.000_r8*rxt(k,292)*y(k,197)  &
                      + .900_r8*rxt(k,293)*y(k,198) + rxt(k,366)*y(k,205)  &
                      + 2.000_r8*rxt(k,342)*y(k,211) + rxt(k,391)*y(k,225)
         mat(k,2033) = .900_r8*rxt(k,293)*y(k,197) + rxt(k,343)*y(k,211)  &
                      + .500_r8*rxt(k,451)*y(k,221)
         mat(k,1910) = .450_r8*rxt(k,344)*y(k,211)
         mat(k,1248) = rxt(k,366)*y(k,197)
         mat(k,1350) = 2.000_r8*rxt(k,342)*y(k,197) + rxt(k,343)*y(k,198)  &
                      + .450_r8*rxt(k,344)*y(k,203) + 4.000_r8*rxt(k,345)*y(k,211)
         mat(k,1644) = rxt(k,296)*y(k,50) + .500_r8*rxt(k,297)*y(k,51) + rxt(k,276) &
                      *y(k,62) + .800_r8*rxt(k,302)*y(k,74) + rxt(k,303)*y(k,75)  &
                      + rxt(k,268)*y(k,87) + .500_r8*rxt(k,352)*y(k,109)  &
                      + 1.800_r8*rxt(k,462)*y(k,178)
         mat(k,1110) = .900_r8*rxt(k,453)*y(k,124) + .500_r8*rxt(k,451)*y(k,198)
         mat(k,1168) = rxt(k,391)*y(k,197)

         mat(k,235) = -(rxt(k,237)*y(k,216))
         mat(k,1512) = -rxt(k,237)*y(k,64)

         mat(k,142) = rxt(k,202)*y(k,216)
         mat(k,147) = rxt(k,228)*y(k,216)
         mat(k,153) = rxt(k,204)*y(k,216)
         mat(k,118) = 2.000_r8*rxt(k,205)*y(k,216)
         mat(k,157) = 2.000_r8*rxt(k,206)*y(k,216)
         mat(k,122) = rxt(k,207)*y(k,216)
         mat(k,106) = 2.000_r8*rxt(k,230)*y(k,216)
         mat(k,247) = rxt(k,254)*y(k,216) + rxt(k,249)*y(k,217)
         mat(k,303) = rxt(k,255)*y(k,216) + rxt(k,250)*y(k,217)
         mat(k,1512) = mat(k,1512) + rxt(k,202)*y(k,34) + rxt(k,228)*y(k,35)  &
                      + rxt(k,204)*y(k,37) + 2.000_r8*rxt(k,205)*y(k,38)  &
                      + 2.000_r8*rxt(k,206)*y(k,39) + rxt(k,207)*y(k,40)  &
                      + 2.000_r8*rxt(k,230)*y(k,78) + rxt(k,254)*y(k,83) + rxt(k,255) &
                      *y(k,84)
         mat(k,1572) = rxt(k,249)*y(k,83) + rxt(k,250)*y(k,84)

         mat(k,226) = -(rxt(k,238)*y(k,216))
         mat(k,1511) = -rxt(k,238)*y(k,65)

         mat(k,114) = rxt(k,203)*y(k,216)
         mat(k,152) = rxt(k,204)*y(k,216)
         mat(k,222) = rxt(k,253)*y(k,216) + rxt(k,248)*y(k,217)
         mat(k,1511) = mat(k,1511) + rxt(k,203)*y(k,36) + rxt(k,204)*y(k,37)  &
                      + rxt(k,253)*y(k,82)
         mat(k,1570) = rxt(k,248)*y(k,82)

         mat(k,194) = -(rxt(k,410)*y(k,217))
         mat(k,1564) = -rxt(k,410)*y(k,66)

         mat(k,188) = .180_r8*rxt(k,430)*y(k,217)
         mat(k,1564) = mat(k,1564) + .180_r8*rxt(k,430)*y(k,180)

         mat(k,297) = -(rxt(k,463)*y(k,126) + (rxt(k,464) + rxt(k,478)) * y(k,217))
         mat(k,1707) = -rxt(k,463)*y(k,67)
         mat(k,1582) = -(rxt(k,464) + rxt(k,478)) * y(k,67)






         mat(k,712) = rxt(k,298)*y(k,203)
         mat(k,1865) = rxt(k,298)*y(k,202)

         mat(k,864) = -(rxt(k,233)*y(k,54) + rxt(k,234)*y(k,77) + rxt(k,235)*y(k,229) &
                      + rxt(k,236)*y(k,89))
         mat(k,1427) = -rxt(k,233)*y(k,73)
         mat(k,1400) = -rxt(k,234)*y(k,73)
         mat(k,2264) = -rxt(k,235)*y(k,73)
         mat(k,1444) = -rxt(k,236)*y(k,73)

         mat(k,148) = rxt(k,228)*y(k,216)
         mat(k,158) = rxt(k,206)*y(k,216)
         mat(k,236) = 2.000_r8*rxt(k,237)*y(k,216)
         mat(k,227) = rxt(k,238)*y(k,216)
         mat(k,1519) = rxt(k,228)*y(k,35) + rxt(k,206)*y(k,39) + 2.000_r8*rxt(k,237) &
                      *y(k,64) + rxt(k,238)*y(k,65)

         mat(k,1096) = -(rxt(k,302)*y(k,217))
         mat(k,1667) = -rxt(k,302)*y(k,74)

         mat(k,582) = .700_r8*rxt(k,377)*y(k,217)
         mat(k,558) = .500_r8*rxt(k,378)*y(k,217)
         mat(k,375) = rxt(k,389)*y(k,217)
         mat(k,1818) = .050_r8*rxt(k,375)*y(k,206) + .530_r8*rxt(k,337)*y(k,209)  &
                      + .225_r8*rxt(k,453)*y(k,221) + .250_r8*rxt(k,394)*y(k,225)
         mat(k,1725) = .050_r8*rxt(k,376)*y(k,206) + .530_r8*rxt(k,336)*y(k,209)  &
                      + .250_r8*rxt(k,395)*y(k,225)
         mat(k,1374) = .530_r8*rxt(k,333)*y(k,209) + .250_r8*rxt(k,391)*y(k,225)
         mat(k,2044) = .260_r8*rxt(k,334)*y(k,209) + .125_r8*rxt(k,451)*y(k,221)  &
                      + .100_r8*rxt(k,392)*y(k,225)
         mat(k,1281) = .050_r8*rxt(k,375)*y(k,124) + .050_r8*rxt(k,376)*y(k,126)
         mat(k,1306) = .530_r8*rxt(k,337)*y(k,124) + .530_r8*rxt(k,336)*y(k,126)  &
                      + .530_r8*rxt(k,333)*y(k,197) + .260_r8*rxt(k,334)*y(k,198)
         mat(k,1667) = mat(k,1667) + .700_r8*rxt(k,377)*y(k,99) + .500_r8*rxt(k,378) &
                      *y(k,100) + rxt(k,389)*y(k,115)
         mat(k,1113) = .225_r8*rxt(k,453)*y(k,124) + .125_r8*rxt(k,451)*y(k,198)
         mat(k,1170) = .250_r8*rxt(k,394)*y(k,124) + .250_r8*rxt(k,395)*y(k,126)  &
                      + .250_r8*rxt(k,391)*y(k,197) + .100_r8*rxt(k,392)*y(k,198)

         mat(k,979) = -(rxt(k,303)*y(k,217))
         mat(k,1660) = -rxt(k,303)*y(k,75)

         mat(k,287) = .650_r8*rxt(k,278)*y(k,217)
         mat(k,1095) = .200_r8*rxt(k,302)*y(k,217)
         mat(k,1072) = rxt(k,390)*y(k,217)
         mat(k,1813) = rxt(k,401)*y(k,191) + .050_r8*rxt(k,375)*y(k,206)  &
                      + .400_r8*rxt(k,415)*y(k,210) + .170_r8*rxt(k,418)*y(k,212)  &
                      + .700_r8*rxt(k,421)*y(k,218) + .600_r8*rxt(k,428)*y(k,223)  &
                      + .250_r8*rxt(k,394)*y(k,225) + .340_r8*rxt(k,434)*y(k,226)  &
                      + .170_r8*rxt(k,437)*y(k,228)
         mat(k,1718) = .050_r8*rxt(k,376)*y(k,206) + .250_r8*rxt(k,395)*y(k,225)
         mat(k,485) = rxt(k,401)*y(k,124)
         mat(k,1372) = .250_r8*rxt(k,391)*y(k,225)
         mat(k,2039) = .100_r8*rxt(k,392)*y(k,225)
         mat(k,1921) = .160_r8*rxt(k,414)*y(k,210) + .070_r8*rxt(k,417)*y(k,212)
         mat(k,1280) = .050_r8*rxt(k,375)*y(k,124) + .050_r8*rxt(k,376)*y(k,126)
         mat(k,722) = .400_r8*rxt(k,415)*y(k,124) + .160_r8*rxt(k,414)*y(k,203)
         mat(k,874) = .170_r8*rxt(k,418)*y(k,124) + .070_r8*rxt(k,417)*y(k,203)
         mat(k,1660) = mat(k,1660) + .650_r8*rxt(k,278)*y(k,24) + .200_r8*rxt(k,302) &
                      *y(k,74) + rxt(k,390)*y(k,116)
         mat(k,455) = .700_r8*rxt(k,421)*y(k,124)
         mat(k,735) = .600_r8*rxt(k,428)*y(k,124)
         mat(k,1169) = .250_r8*rxt(k,394)*y(k,124) + .250_r8*rxt(k,395)*y(k,126)  &
                      + .250_r8*rxt(k,391)*y(k,197) + .100_r8*rxt(k,392)*y(k,198)
         mat(k,751) = .340_r8*rxt(k,434)*y(k,124)
         mat(k,500) = .170_r8*rxt(k,437)*y(k,124)

         mat(k,1463) = -((rxt(k,136) + rxt(k,137) + rxt(k,138)) * y(k,203) + rxt(k,142) &
                      *y(k,134))
         mat(k,1944) = -(rxt(k,136) + rxt(k,137) + rxt(k,138)) * y(k,76)
         mat(k,2122) = -rxt(k,142)*y(k,76)

         mat(k,1484) = rxt(k,259)*y(k,217)
         mat(k,1432) = rxt(k,273)*y(k,216)
         mat(k,2009) = rxt(k,178)*y(k,77)
         mat(k,868) = rxt(k,234)*y(k,77)
         mat(k,1404) = rxt(k,178)*y(k,56) + rxt(k,234)*y(k,73) + rxt(k,134)*y(k,133)  &
                      + rxt(k,125)*y(k,216) + rxt(k,143)*y(k,217)
         mat(k,806) = rxt(k,232)*y(k,216)
         mat(k,2145) = rxt(k,209)*y(k,216)
         mat(k,492) = rxt(k,164)*y(k,217)
         mat(k,2244) = rxt(k,134)*y(k,77) + rxt(k,146)*y(k,217)
         mat(k,362) = rxt(k,466)*y(k,217)
         mat(k,513) = rxt(k,472)*y(k,217)
         mat(k,1234) = rxt(k,477)*y(k,217)
         mat(k,1524) = rxt(k,273)*y(k,54) + rxt(k,125)*y(k,77) + rxt(k,232)*y(k,81)  &
                      + rxt(k,209)*y(k,85)
         mat(k,1688) = rxt(k,259)*y(k,42) + rxt(k,143)*y(k,77) + rxt(k,164)*y(k,112)  &
                      + rxt(k,146)*y(k,133) + rxt(k,466)*y(k,137) + rxt(k,472) &
                      *y(k,148) + rxt(k,477)*y(k,150)

         mat(k,1401) = -(rxt(k,125)*y(k,216) + rxt(k,134)*y(k,133) + rxt(k,143) &
                      *y(k,217) + rxt(k,178)*y(k,56) + rxt(k,234)*y(k,73))
         mat(k,1520) = -rxt(k,125)*y(k,77)
         mat(k,2240) = -rxt(k,134)*y(k,77)
         mat(k,1684) = -rxt(k,143)*y(k,77)
         mat(k,2005) = -rxt(k,178)*y(k,77)
         mat(k,865) = -rxt(k,234)*y(k,77)

         mat(k,1429) = rxt(k,274)*y(k,216)
         mat(k,1460) = rxt(k,136)*y(k,203)
         mat(k,1940) = rxt(k,136)*y(k,76)
         mat(k,1520) = mat(k,1520) + rxt(k,274)*y(k,54)

      end do

      end subroutine     nlnmat03

      subroutine     nlnmat04( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,105) = -(rxt(k,230)*y(k,216))
         mat(k,1500) = -rxt(k,230)*y(k,78)

         mat(k,603) = -(rxt(k,135)*y(k,133) + rxt(k,144)*y(k,217) + rxt(k,179)*y(k,56))
         mat(k,2232) = -rxt(k,135)*y(k,79)
         mat(k,1624) = -rxt(k,144)*y(k,79)
         mat(k,1994) = -rxt(k,179)*y(k,79)

         mat(k,1892) = 2.000_r8*rxt(k,150)*y(k,203)
         mat(k,1624) = mat(k,1624) + 2.000_r8*rxt(k,149)*y(k,217)


         mat(k,256) = rxt(k,479)*y(k,229)
         mat(k,2260) = rxt(k,479)*y(k,152)

         mat(k,804) = -(rxt(k,225)*y(k,133) + rxt(k,226)*y(k,217) + (rxt(k,231) &
                      + rxt(k,232)) * y(k,216))
         mat(k,2234) = -rxt(k,225)*y(k,81)
         mat(k,1645) = -rxt(k,226)*y(k,81)
         mat(k,1518) = -(rxt(k,231) + rxt(k,232)) * y(k,81)

         mat(k,1413) = rxt(k,212)*y(k,42) + rxt(k,213)*y(k,203)
         mat(k,1477) = rxt(k,212)*y(k,17)
         mat(k,1911) = rxt(k,213)*y(k,17)

         mat(k,221) = -(rxt(k,248)*y(k,217) + rxt(k,253)*y(k,216))
         mat(k,1569) = -rxt(k,248)*y(k,82)
         mat(k,1510) = -rxt(k,253)*y(k,82)

         mat(k,248) = -(rxt(k,249)*y(k,217) + rxt(k,254)*y(k,216))
         mat(k,1575) = -rxt(k,249)*y(k,83)
         mat(k,1513) = -rxt(k,254)*y(k,83)

         mat(k,304) = -(rxt(k,250)*y(k,217) + rxt(k,255)*y(k,216))
         mat(k,1583) = -rxt(k,250)*y(k,84)
         mat(k,1514) = -rxt(k,255)*y(k,84)

         mat(k,2156) = -(rxt(k,196)*y(k,133) + rxt(k,197)*y(k,217) + (rxt(k,208) &
                      + rxt(k,209)) * y(k,216) + (rxt(k,524) + rxt(k,530) + rxt(k,535) &
                      ) * y(k,92) + (rxt(k,529) + rxt(k,534) + rxt(k,539)) * y(k,60) &
                      + (rxt(k,531) + rxt(k,536)) * y(k,91))
         mat(k,2255) = -rxt(k,196)*y(k,85)
         mat(k,1699) = -rxt(k,197)*y(k,85)
         mat(k,1535) = -(rxt(k,208) + rxt(k,209)) * y(k,85)
         mat(k,828) = -(rxt(k,524) + rxt(k,530) + rxt(k,535)) * y(k,85)
         mat(k,890) = -(rxt(k,529) + rxt(k,534) + rxt(k,539)) * y(k,85)
         mat(k,782) = -(rxt(k,531) + rxt(k,536)) * y(k,85)

         mat(k,295) = rxt(k,287)*y(k,56)
         mat(k,472) = rxt(k,239)*y(k,56)
         mat(k,1494) = rxt(k,176)*y(k,56)
         mat(k,601) = rxt(k,241)*y(k,56)
         mat(k,371) = 2.000_r8*rxt(k,244)*y(k,56)
         mat(k,1440) = rxt(k,177)*y(k,56)
         mat(k,452) = rxt(k,246)*y(k,56)
         mat(k,2020) = rxt(k,287)*y(k,28) + rxt(k,239)*y(k,41) + rxt(k,176)*y(k,42)  &
                      + rxt(k,241)*y(k,43) + 2.000_r8*rxt(k,244)*y(k,46) + rxt(k,177) &
                      *y(k,54) + rxt(k,246)*y(k,55) + rxt(k,178)*y(k,77) + rxt(k,179) &
                      *y(k,79) + rxt(k,198)*y(k,92) + rxt(k,180)*y(k,203)
         mat(k,1981) = rxt(k,195)*y(k,217)
         mat(k,1410) = rxt(k,178)*y(k,56)
         mat(k,607) = rxt(k,179)*y(k,56)
         mat(k,828) = mat(k,828) + rxt(k,198)*y(k,56)
         mat(k,1955) = rxt(k,180)*y(k,56)
         mat(k,1699) = mat(k,1699) + rxt(k,195)*y(k,59)

         mat(k,179) = -(rxt(k,267)*y(k,217) + rxt(k,275)*y(k,216))
         mat(k,1562) = -rxt(k,267)*y(k,86)
         mat(k,1508) = -rxt(k,275)*y(k,86)

         mat(k,929) = -(rxt(k,268)*y(k,217))
         mat(k,1657) = -rxt(k,268)*y(k,87)

         mat(k,996) = .050_r8*rxt(k,441)*y(k,134)
         mat(k,286) = .350_r8*rxt(k,278)*y(k,217)
         mat(k,550) = .370_r8*rxt(k,280)*y(k,134)
         mat(k,1023) = .120_r8*rxt(k,309)*y(k,134)
         mat(k,842) = .110_r8*rxt(k,386)*y(k,134)
         mat(k,1207) = .330_r8*rxt(k,339)*y(k,134)
         mat(k,947) = .050_r8*rxt(k,444)*y(k,134)
         mat(k,1325) = .120_r8*rxt(k,353)*y(k,134)
         mat(k,1811) = rxt(k,271)*y(k,204)
         mat(k,2095) = .050_r8*rxt(k,441)*y(k,6) + .370_r8*rxt(k,280)*y(k,25)  &
                      + .120_r8*rxt(k,309)*y(k,29) + .110_r8*rxt(k,386)*y(k,98)  &
                      + .330_r8*rxt(k,339)*y(k,105) + .050_r8*rxt(k,444)*y(k,110)  &
                      + .120_r8*rxt(k,353)*y(k,111)
         mat(k,1919) = rxt(k,269)*y(k,204)
         mat(k,442) = rxt(k,271)*y(k,124) + rxt(k,269)*y(k,203)
         mat(k,1657) = mat(k,1657) + .350_r8*rxt(k,278)*y(k,24)


         mat(k,1425) = rxt(k,233)*y(k,73)
         mat(k,863) = rxt(k,233)*y(k,54) + rxt(k,234)*y(k,77) + rxt(k,236)*y(k,89)  &
                      + rxt(k,235)*y(k,229)
         mat(k,1399) = rxt(k,234)*y(k,73)
         mat(k,1443) = rxt(k,236)*y(k,73)
         mat(k,2262) = rxt(k,235)*y(k,73)

         mat(k,1447) = -(rxt(k,173)*y(k,217) + rxt(k,236)*y(k,73))
         mat(k,1687) = -rxt(k,173)*y(k,89)
         mat(k,867) = -rxt(k,236)*y(k,89)

         mat(k,1483) = rxt(k,257)*y(k,126)
         mat(k,1088) = rxt(k,289)*y(k,126)
         mat(k,1223) = rxt(k,315)*y(k,126)
         mat(k,885) = (rxt(k,529)+rxt(k,534)+rxt(k,539))*y(k,85)
         mat(k,299) = rxt(k,463)*y(k,126)
         mat(k,2144) = (rxt(k,529)+rxt(k,534)+rxt(k,539))*y(k,60)
         mat(k,2188) = rxt(k,172)*y(k,217)
         mat(k,1744) = rxt(k,257)*y(k,42) + rxt(k,289)*y(k,45) + rxt(k,315)*y(k,49)  &
                      + rxt(k,463)*y(k,67)
         mat(k,1687) = mat(k,1687) + rxt(k,172)*y(k,125)

         mat(k,403) = -(rxt(k,151)*y(k,217))
         mat(k,1597) = -rxt(k,151)*y(k,90)

         mat(k,2164) = rxt(k,170)*y(k,203)
         mat(k,1877) = rxt(k,170)*y(k,125)

         mat(k,778) = -(rxt(k,227)*y(k,133) + (rxt(k,531) + rxt(k,536)) * y(k,85))
         mat(k,2233) = -rxt(k,227)*y(k,91)
         mat(k,2140) = -(rxt(k,531) + rxt(k,536)) * y(k,91)

         mat(k,2208) = rxt(k,219)*y(k,203)
         mat(k,1908) = rxt(k,219)*y(k,19)

         mat(k,824) = -(rxt(k,198)*y(k,56) + rxt(k,199)*y(k,133) + rxt(k,200)*y(k,217) &
                      + (rxt(k,524) + rxt(k,530) + rxt(k,535)) * y(k,85))
         mat(k,1997) = -rxt(k,198)*y(k,92)
         mat(k,2235) = -rxt(k,199)*y(k,92)
         mat(k,1647) = -rxt(k,200)*y(k,92)
         mat(k,2141) = -(rxt(k,524) + rxt(k,530) + rxt(k,535)) * y(k,92)

         mat(k,1964) = rxt(k,187)*y(k,203)
         mat(k,883) = rxt(k,192)*y(k,217)
         mat(k,1913) = rxt(k,187)*y(k,59)
         mat(k,1647) = mat(k,1647) + rxt(k,192)*y(k,60)

         mat(k,1153) = -(rxt(k,332)*y(k,217))
         mat(k,1672) = -rxt(k,332)*y(k,93)

         mat(k,584) = .300_r8*rxt(k,377)*y(k,217)
         mat(k,560) = .500_r8*rxt(k,378)*y(k,217)
         mat(k,1823) = rxt(k,331)*y(k,200) + rxt(k,338)*y(k,209)
         mat(k,567) = rxt(k,331)*y(k,124)
         mat(k,1309) = rxt(k,338)*y(k,124)
         mat(k,1672) = mat(k,1672) + .300_r8*rxt(k,377)*y(k,99) + .500_r8*rxt(k,378) &
                      *y(k,100)

         mat(k,230) = -(rxt(k,363)*y(k,217))
         mat(k,1571) = -rxt(k,363)*y(k,94)

         mat(k,1140) = -(rxt(k,317)*y(k,217))
         mat(k,1671) = -rxt(k,317)*y(k,95)

         mat(k,583) = .700_r8*rxt(k,377)*y(k,217)
         mat(k,559) = .500_r8*rxt(k,378)*y(k,217)
         mat(k,574) = .500_r8*rxt(k,352)*y(k,217)
         mat(k,1822) = .050_r8*rxt(k,375)*y(k,206) + .220_r8*rxt(k,337)*y(k,209)  &
                      + .250_r8*rxt(k,394)*y(k,225)
         mat(k,1729) = .050_r8*rxt(k,376)*y(k,206) + .220_r8*rxt(k,336)*y(k,209)  &
                      + .250_r8*rxt(k,395)*y(k,225)
         mat(k,535) = .500_r8*rxt(k,321)*y(k,217)
         mat(k,1376) = .220_r8*rxt(k,333)*y(k,209) + .250_r8*rxt(k,391)*y(k,225)
         mat(k,2048) = .230_r8*rxt(k,334)*y(k,209) + .200_r8*rxt(k,322)*y(k,220)  &
                      + .100_r8*rxt(k,392)*y(k,225)
         mat(k,1283) = .050_r8*rxt(k,375)*y(k,124) + .050_r8*rxt(k,376)*y(k,126)
         mat(k,1308) = .220_r8*rxt(k,337)*y(k,124) + .220_r8*rxt(k,336)*y(k,126)  &
                      + .220_r8*rxt(k,333)*y(k,197) + .230_r8*rxt(k,334)*y(k,198)
         mat(k,1671) = mat(k,1671) + .700_r8*rxt(k,377)*y(k,99) + .500_r8*rxt(k,378) &
                      *y(k,100) + .500_r8*rxt(k,352)*y(k,109) + .500_r8*rxt(k,321) &
                      *y(k,146)
         mat(k,1130) = .200_r8*rxt(k,322)*y(k,198)
         mat(k,1172) = .250_r8*rxt(k,394)*y(k,124) + .250_r8*rxt(k,395)*y(k,126)  &
                      + .250_r8*rxt(k,391)*y(k,197) + .100_r8*rxt(k,392)*y(k,198)

         mat(k,320) = -(rxt(k,364)*y(k,217))
         mat(k,1586) = -rxt(k,364)*y(k,96)

         mat(k,1777) = .870_r8*rxt(k,375)*y(k,206)
         mat(k,1708) = .950_r8*rxt(k,376)*y(k,206)
         mat(k,1367) = rxt(k,371)*y(k,206)
         mat(k,2026) = .750_r8*rxt(k,372)*y(k,206)
         mat(k,1273) = .870_r8*rxt(k,375)*y(k,124) + .950_r8*rxt(k,376)*y(k,126)  &
                      + rxt(k,371)*y(k,197) + .750_r8*rxt(k,372)*y(k,198)

         mat(k,135) = -(rxt(k,365)*y(k,217))
         mat(k,1558) = -rxt(k,365)*y(k,97)

         mat(k,689) = .600_r8*rxt(k,388)*y(k,217)
         mat(k,1558) = mat(k,1558) + .600_r8*rxt(k,388)*y(k,103)

         mat(k,841) = -(rxt(k,379)*y(k,126) + rxt(k,386)*y(k,134) + rxt(k,387) &
                      *y(k,217))
         mat(k,1712) = -rxt(k,379)*y(k,98)
         mat(k,2092) = -rxt(k,386)*y(k,98)
         mat(k,1649) = -rxt(k,387)*y(k,98)

         mat(k,581) = -(rxt(k,377)*y(k,217))
         mat(k,1621) = -rxt(k,377)*y(k,99)

         mat(k,1791) = .080_r8*rxt(k,369)*y(k,205)
         mat(k,1245) = .080_r8*rxt(k,369)*y(k,124)

         mat(k,556) = -(rxt(k,378)*y(k,217))
         mat(k,1618) = -rxt(k,378)*y(k,100)

         mat(k,1789) = .080_r8*rxt(k,375)*y(k,206)
         mat(k,1274) = .080_r8*rxt(k,375)*y(k,124)

         mat(k,1193) = -(rxt(k,380)*y(k,197) + rxt(k,381)*y(k,198) + rxt(k,382) &
                      *y(k,203) + rxt(k,383)*y(k,124) + rxt(k,384)*y(k,126))
         mat(k,1378) = -rxt(k,380)*y(k,101)
         mat(k,2050) = -rxt(k,381)*y(k,101)
         mat(k,1931) = -rxt(k,382)*y(k,101)
         mat(k,1825) = -rxt(k,383)*y(k,101)
         mat(k,1732) = -rxt(k,384)*y(k,101)

         mat(k,845) = rxt(k,379)*y(k,126)
         mat(k,1732) = mat(k,1732) + rxt(k,379)*y(k,98)

         mat(k,397) = -(rxt(k,385)*y(k,217))
         mat(k,1596) = -rxt(k,385)*y(k,102)

         mat(k,1185) = rxt(k,382)*y(k,203)
         mat(k,1876) = rxt(k,382)*y(k,101)

         mat(k,690) = -(rxt(k,388)*y(k,217))
         mat(k,1634) = -rxt(k,388)*y(k,103)

         mat(k,1900) = rxt(k,368)*y(k,205) + rxt(k,373)*y(k,206)
         mat(k,1246) = rxt(k,368)*y(k,203)
         mat(k,1276) = rxt(k,373)*y(k,203)

         mat(k,74) = -(rxt(k,509)*y(k,217))
         mat(k,1550) = -rxt(k,509)*y(k,104)

         mat(k,1209) = -(rxt(k,339)*y(k,134) + rxt(k,340)*y(k,217))
         mat(k,2110) = -rxt(k,339)*y(k,105)
         mat(k,1675) = -rxt(k,340)*y(k,105)

         mat(k,846) = .300_r8*rxt(k,386)*y(k,134)
         mat(k,1826) = .360_r8*rxt(k,369)*y(k,205)
         mat(k,1733) = .400_r8*rxt(k,370)*y(k,205)
         mat(k,2110) = mat(k,2110) + .300_r8*rxt(k,386)*y(k,98)
         mat(k,1379) = .390_r8*rxt(k,366)*y(k,205)
         mat(k,2051) = .310_r8*rxt(k,367)*y(k,205)
         mat(k,1254) = .360_r8*rxt(k,369)*y(k,124) + .400_r8*rxt(k,370)*y(k,126)  &
                      + .390_r8*rxt(k,366)*y(k,197) + .310_r8*rxt(k,367)*y(k,198)

      end do

      end subroutine     nlnmat04

      subroutine     nlnmat05( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,310) = -(rxt(k,341)*y(k,217))
         mat(k,1584) = -rxt(k,341)*y(k,106)

         mat(k,1869) = rxt(k,335)*y(k,209)
         mat(k,1304) = rxt(k,335)*y(k,203)

         mat(k,506) = -(rxt(k,350)*y(k,217))
         mat(k,1612) = -rxt(k,350)*y(k,107)

         mat(k,1787) = .800_r8*rxt(k,359)*y(k,189)
         mat(k,894) = .800_r8*rxt(k,359)*y(k,124)

         mat(k,315) = -(rxt(k,351)*y(k,217))
         mat(k,1585) = -rxt(k,351)*y(k,108)

         mat(k,1870) = .800_r8*rxt(k,348)*y(k,213)
         mat(k,676) = .800_r8*rxt(k,348)*y(k,203)

         mat(k,572) = -(rxt(k,352)*y(k,217))
         mat(k,1620) = -rxt(k,352)*y(k,109)

         mat(k,2170) = rxt(k,355)*y(k,211)
         mat(k,1347) = rxt(k,355)*y(k,125)

         mat(k,948) = -(rxt(k,443)*y(k,126) + rxt(k,444)*y(k,134) + rxt(k,445) &
                      *y(k,217))
         mat(k,1716) = -rxt(k,443)*y(k,110)
         mat(k,2096) = -rxt(k,444)*y(k,110)
         mat(k,1658) = -rxt(k,445)*y(k,110)

         mat(k,1332) = -(rxt(k,353)*y(k,134) + rxt(k,354)*y(k,217))
         mat(k,2116) = -rxt(k,353)*y(k,111)
         mat(k,1681) = -rxt(k,354)*y(k,111)

         mat(k,849) = .200_r8*rxt(k,386)*y(k,134)
         mat(k,1831) = .560_r8*rxt(k,369)*y(k,205)
         mat(k,1739) = .600_r8*rxt(k,370)*y(k,205)
         mat(k,2116) = mat(k,2116) + .200_r8*rxt(k,386)*y(k,98)
         mat(k,1384) = .610_r8*rxt(k,366)*y(k,205)
         mat(k,2056) = .440_r8*rxt(k,367)*y(k,205)
         mat(k,1258) = .560_r8*rxt(k,369)*y(k,124) + .600_r8*rxt(k,370)*y(k,126)  &
                      + .610_r8*rxt(k,366)*y(k,197) + .440_r8*rxt(k,367)*y(k,198)

         mat(k,491) = -(rxt(k,152)*y(k,124) + (rxt(k,153) + rxt(k,154) + rxt(k,155) &
                      ) * y(k,125) + rxt(k,164)*y(k,217))
         mat(k,1785) = -rxt(k,152)*y(k,112)
         mat(k,2166) = -(rxt(k,153) + rxt(k,154) + rxt(k,155)) * y(k,112)
         mat(k,1610) = -rxt(k,164)*y(k,112)

         mat(k,183) = -((rxt(k,168) + rxt(k,169)) * y(k,216))
         mat(k,1509) = -(rxt(k,168) + rxt(k,169)) * y(k,113)

         mat(k,490) = rxt(k,153)*y(k,125)
         mat(k,2162) = rxt(k,153)*y(k,112)


         mat(k,2163) = rxt(k,171)*y(k,126)
         mat(k,1706) = rxt(k,171)*y(k,125)

         mat(k,373) = -(rxt(k,389)*y(k,217))
         mat(k,1593) = -rxt(k,389)*y(k,115)

         mat(k,1184) = .200_r8*rxt(k,381)*y(k,198)
         mat(k,2027) = .200_r8*rxt(k,381)*y(k,101)

         mat(k,1073) = -(rxt(k,390)*y(k,217))
         mat(k,1665) = -rxt(k,390)*y(k,116)

         mat(k,1189) = rxt(k,383)*y(k,124) + rxt(k,384)*y(k,126) + rxt(k,380)*y(k,197)  &
                      + .800_r8*rxt(k,381)*y(k,198)
         mat(k,1816) = rxt(k,383)*y(k,101)
         mat(k,1723) = rxt(k,384)*y(k,101)
         mat(k,1373) = rxt(k,380)*y(k,101)
         mat(k,2042) = .800_r8*rxt(k,381)*y(k,101)




         mat(k,96) = -(rxt(k,480)*y(k,217))
         mat(k,1554) = -rxt(k,480)*y(k,120)




         mat(k,1842) = -(rxt(k,152)*y(k,112) + rxt(k,161)*y(k,126) + rxt(k,165) &
                      *y(k,203) + rxt(k,166)*y(k,134) + rxt(k,167)*y(k,133) + rxt(k,188) &
                      *y(k,59) + rxt(k,220)*y(k,19) + rxt(k,263)*y(k,198) + rxt(k,271) &
                      *y(k,204) + rxt(k,284)*y(k,194) + rxt(k,295)*y(k,197) + rxt(k,299) &
                      *y(k,202) + rxt(k,312)*y(k,195) + rxt(k,320)*y(k,219) + rxt(k,324) &
                      *y(k,220) + (rxt(k,330) + rxt(k,331)) * y(k,200) + (rxt(k,337) &
                      + rxt(k,338)) * y(k,209) + rxt(k,346)*y(k,211) + rxt(k,349) &
                      *y(k,213) + (rxt(k,359) + rxt(k,360)) * y(k,189) + rxt(k,369) &
                      *y(k,205) + rxt(k,375)*y(k,206) + rxt(k,383)*y(k,101) + rxt(k,394) &
                      *y(k,225) + rxt(k,398)*y(k,188) + rxt(k,401)*y(k,191) + rxt(k,406) &
                      *y(k,193) + rxt(k,408)*y(k,196) + rxt(k,412)*y(k,199) + rxt(k,415) &
                      *y(k,210) + rxt(k,418)*y(k,212) + rxt(k,421)*y(k,218) + rxt(k,428) &
                      *y(k,223) + rxt(k,434)*y(k,226) + rxt(k,437)*y(k,228) + rxt(k,448) &
                      *y(k,215) + rxt(k,453)*y(k,221) + rxt(k,458)*y(k,222))
         mat(k,495) = -rxt(k,152)*y(k,124)
         mat(k,1750) = -rxt(k,161)*y(k,124)
         mat(k,1949) = -rxt(k,165)*y(k,124)
         mat(k,2127) = -rxt(k,166)*y(k,124)
         mat(k,2249) = -rxt(k,167)*y(k,124)
         mat(k,1975) = -rxt(k,188)*y(k,124)
         mat(k,2218) = -rxt(k,220)*y(k,124)
         mat(k,2066) = -rxt(k,263)*y(k,124)
         mat(k,444) = -rxt(k,271)*y(k,124)
         mat(k,819) = -rxt(k,284)*y(k,124)
         mat(k,1392) = -rxt(k,295)*y(k,124)
         mat(k,718) = -rxt(k,299)*y(k,124)
         mat(k,796) = -rxt(k,312)*y(k,124)
         mat(k,773) = -rxt(k,320)*y(k,124)
         mat(k,1135) = -rxt(k,324)*y(k,124)
         mat(k,569) = -(rxt(k,330) + rxt(k,331)) * y(k,124)
         mat(k,1318) = -(rxt(k,337) + rxt(k,338)) * y(k,124)
         mat(k,1360) = -rxt(k,346)*y(k,124)
         mat(k,681) = -rxt(k,349)*y(k,124)
         mat(k,905) = -(rxt(k,359) + rxt(k,360)) * y(k,124)
         mat(k,1265) = -rxt(k,369)*y(k,124)
         mat(k,1297) = -rxt(k,375)*y(k,124)
         mat(k,1202) = -rxt(k,383)*y(k,124)
         mat(k,1179) = -rxt(k,394)*y(k,124)
         mat(k,521) = -rxt(k,398)*y(k,124)
         mat(k,487) = -rxt(k,401)*y(k,124)
         mat(k,438) = -rxt(k,406)*y(k,124)
         mat(k,627) = -rxt(k,408)*y(k,124)
         mat(k,763) = -rxt(k,412)*y(k,124)
         mat(k,724) = -rxt(k,415)*y(k,124)
         mat(k,878) = -rxt(k,418)*y(k,124)
         mat(k,457) = -rxt(k,421)*y(k,124)
         mat(k,739) = -rxt(k,428)*y(k,124)
         mat(k,756) = -rxt(k,434)*y(k,124)
         mat(k,503) = -rxt(k,437)*y(k,124)
         mat(k,1053) = -rxt(k,448)*y(k,124)
         mat(k,1121) = -rxt(k,453)*y(k,124)
         mat(k,918) = -rxt(k,458)*y(k,124)

         mat(k,495) = mat(k,495) + 2.000_r8*rxt(k,154)*y(k,125) + rxt(k,164)*y(k,217)
         mat(k,185) = 2.000_r8*rxt(k,168)*y(k,216)
         mat(k,2194) = 2.000_r8*rxt(k,154)*y(k,112) + rxt(k,157)*y(k,133) + rxt(k,473) &
                      *y(k,150)
         mat(k,2249) = mat(k,2249) + rxt(k,157)*y(k,125)
         mat(k,1236) = rxt(k,473)*y(k,125)
         mat(k,1529) = 2.000_r8*rxt(k,168)*y(k,113)
         mat(k,1693) = rxt(k,164)*y(k,112)

         mat(k,2201) = -((rxt(k,153) + rxt(k,154) + rxt(k,155)) * y(k,112) + (rxt(k,157) &
                      + rxt(k,159)) * y(k,133) + rxt(k,158)*y(k,134) + rxt(k,170) &
                      *y(k,203) + rxt(k,171)*y(k,126) + rxt(k,172)*y(k,217) + rxt(k,190) &
                      *y(k,59) + rxt(k,221)*y(k,19) + rxt(k,306)*y(k,197) + rxt(k,355) &
                      *y(k,211) + rxt(k,413)*y(k,199) + rxt(k,416)*y(k,210) + rxt(k,419) &
                      *y(k,212) + rxt(k,423)*y(k,141) + rxt(k,426)*y(k,188) + rxt(k,473) &
                      *y(k,150))
         mat(k,496) = -(rxt(k,153) + rxt(k,154) + rxt(k,155)) * y(k,125)
         mat(k,2256) = -(rxt(k,157) + rxt(k,159)) * y(k,125)
         mat(k,2134) = -rxt(k,158)*y(k,125)
         mat(k,1956) = -rxt(k,170)*y(k,125)
         mat(k,1757) = -rxt(k,171)*y(k,125)
         mat(k,1700) = -rxt(k,172)*y(k,125)
         mat(k,1982) = -rxt(k,190)*y(k,125)
         mat(k,2225) = -rxt(k,221)*y(k,125)
         mat(k,1396) = -rxt(k,306)*y(k,125)
         mat(k,1364) = -rxt(k,355)*y(k,125)
         mat(k,766) = -rxt(k,413)*y(k,125)
         mat(k,726) = -rxt(k,416)*y(k,125)
         mat(k,881) = -rxt(k,419)*y(k,125)
         mat(k,466) = -rxt(k,423)*y(k,125)
         mat(k,523) = -rxt(k,426)*y(k,125)
         mat(k,1241) = -rxt(k,473)*y(k,125)

         mat(k,675) = rxt(k,357)*y(k,217)
         mat(k,356) = rxt(k,328)*y(k,126)
         mat(k,2225) = mat(k,2225) + rxt(k,220)*y(k,124)
         mat(k,1982) = mat(k,1982) + rxt(k,188)*y(k,124)
         mat(k,407) = rxt(k,151)*y(k,217)
         mat(k,589) = .700_r8*rxt(k,377)*y(k,217)
         mat(k,1205) = rxt(k,383)*y(k,124) + rxt(k,384)*y(k,126)
         mat(k,1849) = rxt(k,220)*y(k,19) + rxt(k,188)*y(k,59) + rxt(k,383)*y(k,101)  &
                      + 2.000_r8*rxt(k,161)*y(k,126) + rxt(k,167)*y(k,133)  &
                      + rxt(k,166)*y(k,134) + rxt(k,398)*y(k,188) + rxt(k,359) &
                      *y(k,189) + rxt(k,401)*y(k,191) + rxt(k,406)*y(k,193)  &
                      + rxt(k,284)*y(k,194) + rxt(k,312)*y(k,195) + rxt(k,408) &
                      *y(k,196) + rxt(k,295)*y(k,197) + rxt(k,263)*y(k,198)  &
                      + rxt(k,412)*y(k,199) + rxt(k,330)*y(k,200) + rxt(k,299) &
                      *y(k,202) + rxt(k,165)*y(k,203) + rxt(k,271)*y(k,204)  &
                      + .920_r8*rxt(k,369)*y(k,205) + .920_r8*rxt(k,375)*y(k,206)  &
                      + rxt(k,337)*y(k,209) + rxt(k,415)*y(k,210) + rxt(k,346) &
                      *y(k,211) + rxt(k,418)*y(k,212) + rxt(k,349)*y(k,213)  &
                      + 1.600_r8*rxt(k,448)*y(k,215) + rxt(k,421)*y(k,218)  &
                      + rxt(k,320)*y(k,219) + rxt(k,324)*y(k,220) + .900_r8*rxt(k,453) &
                      *y(k,221) + .800_r8*rxt(k,458)*y(k,222) + rxt(k,428)*y(k,223)  &
                      + rxt(k,394)*y(k,225) + rxt(k,434)*y(k,226) + rxt(k,437) &
                      *y(k,228)
         mat(k,1757) = mat(k,1757) + rxt(k,328)*y(k,16) + rxt(k,384)*y(k,101)  &
                      + 2.000_r8*rxt(k,161)*y(k,124) + rxt(k,162)*y(k,133)  &
                      + rxt(k,160)*y(k,203) + rxt(k,370)*y(k,205) + rxt(k,376) &
                      *y(k,206) + rxt(k,336)*y(k,209) + rxt(k,347)*y(k,211)  &
                      + 2.000_r8*rxt(k,449)*y(k,215) + rxt(k,163)*y(k,217)  &
                      + rxt(k,395)*y(k,225)
         mat(k,862) = rxt(k,318)*y(k,217)
         mat(k,2256) = mat(k,2256) + rxt(k,167)*y(k,124) + rxt(k,162)*y(k,126)
         mat(k,2134) = mat(k,2134) + rxt(k,166)*y(k,124)
         mat(k,622) = rxt(k,455)*y(k,217)
         mat(k,523) = mat(k,523) + rxt(k,398)*y(k,124)
         mat(k,908) = rxt(k,359)*y(k,124)
         mat(k,489) = rxt(k,401)*y(k,124)
         mat(k,440) = rxt(k,406)*y(k,124)
         mat(k,822) = rxt(k,284)*y(k,124)
         mat(k,799) = rxt(k,312)*y(k,124)
         mat(k,630) = rxt(k,408)*y(k,124)
         mat(k,1396) = mat(k,1396) + rxt(k,295)*y(k,124)
         mat(k,2073) = rxt(k,263)*y(k,124) + .500_r8*rxt(k,446)*y(k,215)
         mat(k,766) = mat(k,766) + rxt(k,412)*y(k,124)
         mat(k,571) = rxt(k,330)*y(k,124)
         mat(k,720) = rxt(k,299)*y(k,124)
         mat(k,1956) = mat(k,1956) + rxt(k,165)*y(k,124) + rxt(k,160)*y(k,126)
         mat(k,446) = rxt(k,271)*y(k,124)
         mat(k,1269) = .920_r8*rxt(k,369)*y(k,124) + rxt(k,370)*y(k,126)
         mat(k,1301) = .920_r8*rxt(k,375)*y(k,124) + rxt(k,376)*y(k,126)
         mat(k,1321) = rxt(k,337)*y(k,124) + rxt(k,336)*y(k,126)
         mat(k,726) = mat(k,726) + rxt(k,415)*y(k,124)
         mat(k,1364) = mat(k,1364) + rxt(k,346)*y(k,124) + rxt(k,347)*y(k,126)
         mat(k,881) = mat(k,881) + rxt(k,418)*y(k,124)
         mat(k,683) = rxt(k,349)*y(k,124)
         mat(k,1056) = 1.600_r8*rxt(k,448)*y(k,124) + 2.000_r8*rxt(k,449)*y(k,126)  &
                      + .500_r8*rxt(k,446)*y(k,198)
         mat(k,1700) = mat(k,1700) + rxt(k,357)*y(k,1) + rxt(k,151)*y(k,90)  &
                      + .700_r8*rxt(k,377)*y(k,99) + rxt(k,163)*y(k,126) + rxt(k,318) &
                      *y(k,127) + rxt(k,455)*y(k,175)
         mat(k,459) = rxt(k,421)*y(k,124)
         mat(k,775) = rxt(k,320)*y(k,124)
         mat(k,1138) = rxt(k,324)*y(k,124)
         mat(k,1124) = .900_r8*rxt(k,453)*y(k,124)
         mat(k,921) = .800_r8*rxt(k,458)*y(k,124)
         mat(k,741) = rxt(k,428)*y(k,124)
         mat(k,1182) = rxt(k,394)*y(k,124) + rxt(k,395)*y(k,126)
         mat(k,758) = rxt(k,434)*y(k,124)
         mat(k,505) = rxt(k,437)*y(k,124)

      end do

      end subroutine     nlnmat05

      subroutine     nlnmat06( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,1749) = -(rxt(k,160)*y(k,203) + rxt(k,161)*y(k,124) + rxt(k,162) &
                      *y(k,133) + rxt(k,163)*y(k,217) + rxt(k,171)*y(k,125) + rxt(k,257) &
                      *y(k,42) + rxt(k,289)*y(k,45) + rxt(k,308)*y(k,29) + rxt(k,315) &
                      *y(k,49) + rxt(k,328)*y(k,16) + rxt(k,336)*y(k,209) + rxt(k,347) &
                      *y(k,211) + rxt(k,370)*y(k,205) + rxt(k,376)*y(k,206) + rxt(k,379) &
                      *y(k,98) + rxt(k,384)*y(k,101) + rxt(k,395)*y(k,225) + rxt(k,440) &
                      *y(k,6) + rxt(k,443)*y(k,110) + rxt(k,449)*y(k,215) + rxt(k,460) &
                      *y(k,177) + rxt(k,463)*y(k,67))
         mat(k,1948) = -rxt(k,160)*y(k,126)
         mat(k,1841) = -rxt(k,161)*y(k,126)
         mat(k,2248) = -rxt(k,162)*y(k,126)
         mat(k,1692) = -rxt(k,163)*y(k,126)
         mat(k,2193) = -rxt(k,171)*y(k,126)
         mat(k,1488) = -rxt(k,257)*y(k,126)
         mat(k,1090) = -rxt(k,289)*y(k,126)
         mat(k,1033) = -rxt(k,308)*y(k,126)
         mat(k,1225) = -rxt(k,315)*y(k,126)
         mat(k,355) = -rxt(k,328)*y(k,126)
         mat(k,1317) = -rxt(k,336)*y(k,126)
         mat(k,1359) = -rxt(k,347)*y(k,126)
         mat(k,1264) = -rxt(k,370)*y(k,126)
         mat(k,1296) = -rxt(k,376)*y(k,126)
         mat(k,853) = -rxt(k,379)*y(k,126)
         mat(k,1201) = -rxt(k,384)*y(k,126)
         mat(k,1178) = -rxt(k,395)*y(k,126)
         mat(k,1011) = -rxt(k,440)*y(k,126)
         mat(k,961) = -rxt(k,443)*y(k,126)
         mat(k,1052) = -rxt(k,449)*y(k,126)
         mat(k,975) = -rxt(k,460)*y(k,126)
         mat(k,301) = -rxt(k,463)*y(k,126)

         mat(k,544) = rxt(k,222)*y(k,133)
         mat(k,2013) = rxt(k,189)*y(k,60)
         mat(k,887) = rxt(k,189)*y(k,56) + rxt(k,191)*y(k,133) + rxt(k,192)*y(k,217)
         mat(k,870) = rxt(k,236)*y(k,89)
         mat(k,1452) = rxt(k,236)*y(k,73) + rxt(k,173)*y(k,217)
         mat(k,578) = .500_r8*rxt(k,352)*y(k,217)
         mat(k,2193) = mat(k,2193) + rxt(k,159)*y(k,133) + rxt(k,158)*y(k,134)
         mat(k,2248) = mat(k,2248) + rxt(k,222)*y(k,20) + rxt(k,191)*y(k,60)  &
                      + rxt(k,159)*y(k,125)
         mat(k,2126) = rxt(k,158)*y(k,125)
         mat(k,529) = rxt(k,304)*y(k,217)
         mat(k,1692) = mat(k,1692) + rxt(k,192)*y(k,60) + rxt(k,173)*y(k,89)  &
                      + .500_r8*rxt(k,352)*y(k,109) + rxt(k,304)*y(k,139)

         mat(k,857) = -(rxt(k,318)*y(k,217))
         mat(k,1650) = -rxt(k,318)*y(k,127)

         mat(k,1022) = rxt(k,308)*y(k,126)
         mat(k,557) = .500_r8*rxt(k,378)*y(k,217)
         mat(k,399) = rxt(k,385)*y(k,217)
         mat(k,374) = rxt(k,389)*y(k,217)
         mat(k,1070) = rxt(k,390)*y(k,217)
         mat(k,1713) = rxt(k,308)*y(k,29)
         mat(k,1650) = mat(k,1650) + .500_r8*rxt(k,378)*y(k,100) + rxt(k,385)*y(k,102)  &
                      + rxt(k,389)*y(k,115) + rxt(k,390)*y(k,116)

         mat(k,391) = -(rxt(k,450)*y(k,217))
         mat(k,1595) = -rxt(k,450)*y(k,128)

         mat(k,1875) = rxt(k,447)*y(k,215)
         mat(k,1041) = rxt(k,447)*y(k,203)





         mat(k,2258) = -(rxt(k,131)*y(k,134) + 4._r8*rxt(k,132)*y(k,133) + rxt(k,134) &
                      *y(k,77) + rxt(k,135)*y(k,79) + rxt(k,140)*y(k,203) + rxt(k,146) &
                      *y(k,217) + (rxt(k,157) + rxt(k,159)) * y(k,125) + rxt(k,162) &
                      *y(k,126) + rxt(k,167)*y(k,124) + rxt(k,191)*y(k,60) + rxt(k,193) &
                      *y(k,59) + rxt(k,196)*y(k,85) + rxt(k,199)*y(k,92) + rxt(k,222) &
                      *y(k,20) + rxt(k,223)*y(k,19) + rxt(k,225)*y(k,81) + rxt(k,227) &
                      *y(k,91) + rxt(k,258)*y(k,42) + rxt(k,465)*y(k,137))
         mat(k,2136) = -rxt(k,131)*y(k,133)
         mat(k,1411) = -rxt(k,134)*y(k,133)
         mat(k,608) = -rxt(k,135)*y(k,133)
         mat(k,1958) = -rxt(k,140)*y(k,133)
         mat(k,1702) = -rxt(k,146)*y(k,133)
         mat(k,2203) = -(rxt(k,157) + rxt(k,159)) * y(k,133)
         mat(k,1759) = -rxt(k,162)*y(k,133)
         mat(k,1851) = -rxt(k,167)*y(k,133)
         mat(k,892) = -rxt(k,191)*y(k,133)
         mat(k,1984) = -rxt(k,193)*y(k,133)
         mat(k,2159) = -rxt(k,196)*y(k,133)
         mat(k,829) = -rxt(k,199)*y(k,133)
         mat(k,547) = -rxt(k,222)*y(k,133)
         mat(k,2227) = -rxt(k,223)*y(k,133)
         mat(k,810) = -rxt(k,225)*y(k,133)
         mat(k,784) = -rxt(k,227)*y(k,133)
         mat(k,1497) = -rxt(k,258)*y(k,133)
         mat(k,364) = -rxt(k,465)*y(k,133)

         mat(k,1474) = rxt(k,138)*y(k,203)
         mat(k,497) = rxt(k,152)*y(k,124) + rxt(k,153)*y(k,125)
         mat(k,1851) = mat(k,1851) + rxt(k,152)*y(k,112)
         mat(k,2203) = mat(k,2203) + rxt(k,153)*y(k,112)
         mat(k,2136) = mat(k,2136) + 2.000_r8*rxt(k,130)*y(k,216)
         mat(k,1958) = mat(k,1958) + rxt(k,138)*y(k,76)
         mat(k,1538) = 2.000_r8*rxt(k,130)*y(k,134)
         mat(k,1702) = mat(k,1702) + 2.000_r8*rxt(k,148)*y(k,217)

         mat(k,2132) = -((rxt(k,129) + rxt(k,130)) * y(k,216) + rxt(k,131)*y(k,133) &
                      + rxt(k,141)*y(k,203) + rxt(k,142)*y(k,76) + rxt(k,147)*y(k,217) &
                      + rxt(k,158)*y(k,125) + rxt(k,166)*y(k,124) + rxt(k,182)*y(k,56) &
                      + rxt(k,214)*y(k,17) + rxt(k,280)*y(k,25) + rxt(k,309)*y(k,29) &
                      + rxt(k,339)*y(k,105) + rxt(k,353)*y(k,111) + rxt(k,386)*y(k,98) &
                      + rxt(k,424)*y(k,141) + rxt(k,441)*y(k,6) + rxt(k,444)*y(k,110) &
                      + rxt(k,469)*y(k,148) + rxt(k,475)*y(k,150))
         mat(k,1534) = -(rxt(k,129) + rxt(k,130)) * y(k,134)
         mat(k,2254) = -rxt(k,131)*y(k,134)
         mat(k,1954) = -rxt(k,141)*y(k,134)
         mat(k,1471) = -rxt(k,142)*y(k,134)
         mat(k,1698) = -rxt(k,147)*y(k,134)
         mat(k,2199) = -rxt(k,158)*y(k,134)
         mat(k,1847) = -rxt(k,166)*y(k,134)
         mat(k,2019) = -rxt(k,182)*y(k,134)
         mat(k,1421) = -rxt(k,214)*y(k,134)
         mat(k,555) = -rxt(k,280)*y(k,134)
         mat(k,1037) = -rxt(k,309)*y(k,134)
         mat(k,1217) = -rxt(k,339)*y(k,134)
         mat(k,1343) = -rxt(k,353)*y(k,134)
         mat(k,856) = -rxt(k,386)*y(k,134)
         mat(k,465) = -rxt(k,424)*y(k,134)
         mat(k,1015) = -rxt(k,441)*y(k,134)
         mat(k,965) = -rxt(k,444)*y(k,134)
         mat(k,515) = -rxt(k,469)*y(k,134)
         mat(k,1240) = -rxt(k,475)*y(k,134)

         mat(k,1395) = .150_r8*rxt(k,294)*y(k,203)
         mat(k,1954) = mat(k,1954) + .150_r8*rxt(k,294)*y(k,197) + .150_r8*rxt(k,344) &
                      *y(k,211)
         mat(k,1363) = .150_r8*rxt(k,344)*y(k,203)


         mat(k,328) = -(rxt(k,476)*y(k,150))
         mat(k,1229) = -rxt(k,476)*y(k,136)

         mat(k,2206) = rxt(k,216)*y(k,59)
         mat(k,1963) = rxt(k,216)*y(k,19) + 2.000_r8*rxt(k,186)*y(k,59)

         mat(k,357) = -(rxt(k,465)*y(k,133) + rxt(k,466)*y(k,217))
         mat(k,2229) = -rxt(k,465)*y(k,137)
         mat(k,1591) = -rxt(k,466)*y(k,137)


         mat(k,1146) = rxt(k,332)*y(k,217)
         mat(k,1773) = .100_r8*rxt(k,453)*y(k,221)
         mat(k,1574) = rxt(k,332)*y(k,93)
         mat(k,1107) = .100_r8*rxt(k,453)*y(k,124)

         mat(k,524) = -(rxt(k,304)*y(k,217))
         mat(k,1615) = -rxt(k,304)*y(k,139)

         mat(k,2168) = rxt(k,306)*y(k,197)
         mat(k,1368) = rxt(k,306)*y(k,125)


         mat(k,2161) = rxt(k,426)*y(k,188)
         mat(k,517) = rxt(k,426)*y(k,125)

         mat(k,463) = -(rxt(k,423)*y(k,125) + rxt(k,424)*y(k,134))
         mat(k,2165) = -rxt(k,423)*y(k,141)
         mat(k,2084) = -rxt(k,424)*y(k,141)

         mat(k,196) = .070_r8*rxt(k,410)*y(k,217)
         mat(k,1783) = rxt(k,408)*y(k,196)
         mat(k,176) = .060_r8*rxt(k,422)*y(k,217)
         mat(k,217) = .070_r8*rxt(k,438)*y(k,217)
         mat(k,624) = rxt(k,408)*y(k,124)
         mat(k,1606) = .070_r8*rxt(k,410)*y(k,66) + .060_r8*rxt(k,422)*y(k,142)  &
                      + .070_r8*rxt(k,438)*y(k,184)

         mat(k,174) = -(rxt(k,422)*y(k,217))
         mat(k,1561) = -rxt(k,422)*y(k,142)

         mat(k,166) = .530_r8*rxt(k,399)*y(k,217)
         mat(k,1561) = mat(k,1561) + .530_r8*rxt(k,399)*y(k,7)

         mat(k,333) = -(rxt(k,425)*y(k,217))
         mat(k,1587) = -rxt(k,425)*y(k,143)

         mat(k,1871) = rxt(k,420)*y(k,218)
         mat(k,453) = rxt(k,420)*y(k,203)



         mat(k,532) = -(rxt(k,321)*y(k,217))
         mat(k,1616) = -rxt(k,321)*y(k,146)

         mat(k,1891) = rxt(k,319)*y(k,219)
         mat(k,767) = rxt(k,319)*y(k,203)

         mat(k,409) = -(rxt(k,325)*y(k,217))
         mat(k,1598) = -rxt(k,325)*y(k,147)

         mat(k,1878) = .850_r8*rxt(k,323)*y(k,220)
         mat(k,1127) = .850_r8*rxt(k,323)*y(k,203)

         mat(k,511) = -(rxt(k,469)*y(k,134) + rxt(k,472)*y(k,217))
         mat(k,2085) = -rxt(k,469)*y(k,148)
         mat(k,1613) = -rxt(k,472)*y(k,148)


         mat(k,1232) = -(rxt(k,470)*y(k,19) + rxt(k,471)*y(k,59) + rxt(k,473)*y(k,125) &
                      + rxt(k,475)*y(k,134) + rxt(k,476)*y(k,136) + rxt(k,477) &
                      *y(k,217))
         mat(k,2210) = -rxt(k,470)*y(k,150)
         mat(k,1967) = -rxt(k,471)*y(k,150)
         mat(k,2183) = -rxt(k,473)*y(k,150)
         mat(k,2112) = -rxt(k,475)*y(k,150)
         mat(k,330) = -rxt(k,476)*y(k,150)
         mat(k,1677) = -rxt(k,477)*y(k,150)

         mat(k,2239) = rxt(k,465)*y(k,137)
         mat(k,2112) = mat(k,2112) + rxt(k,469)*y(k,148)
         mat(k,361) = rxt(k,465)*y(k,133)
         mat(k,512) = rxt(k,469)*y(k,134) + rxt(k,472)*y(k,217)
         mat(k,1677) = mat(k,1677) + rxt(k,472)*y(k,148)

         mat(k,832) = -(rxt(k,468)*y(k,217))
         mat(k,1648) = -rxt(k,468)*y(k,151)

         mat(k,2209) = rxt(k,470)*y(k,150)
         mat(k,1965) = rxt(k,471)*y(k,150)
         mat(k,298) = rxt(k,463)*y(k,126) + (rxt(k,464)+.500_r8*rxt(k,478))*y(k,217)
         mat(k,2176) = rxt(k,473)*y(k,150)
         mat(k,1711) = rxt(k,463)*y(k,67)
         mat(k,2091) = rxt(k,475)*y(k,150)
         mat(k,329) = rxt(k,476)*y(k,150)
         mat(k,359) = rxt(k,466)*y(k,217)
         mat(k,1231) = rxt(k,470)*y(k,19) + rxt(k,471)*y(k,59) + rxt(k,473)*y(k,125)  &
                      + rxt(k,475)*y(k,134) + rxt(k,476)*y(k,136) + rxt(k,477) &
                      *y(k,217)
         mat(k,1648) = mat(k,1648) + (rxt(k,464)+.500_r8*rxt(k,478))*y(k,67)  &
                      + rxt(k,466)*y(k,137) + rxt(k,477)*y(k,150)

         mat(k,257) = -(rxt(k,479)*y(k,229))
         mat(k,2261) = -rxt(k,479)*y(k,152)

         mat(k,831) = rxt(k,468)*y(k,217)
         mat(k,1577) = rxt(k,468)*y(k,151)















         mat(k,984) = .2202005_r8*rxt(k,497)*y(k,134)
         mat(k,935) = .0508005_r8*rxt(k,513)*y(k,134)
         mat(k,1761) = .1279005_r8*rxt(k,496)*y(k,190) + .0097005_r8*rxt(k,501) &
                      *y(k,192) + .0003005_r8*rxt(k,504)*y(k,207)  &
                      + .1056005_r8*rxt(k,508)*y(k,208) + .0245005_r8*rxt(k,512) &
                      *y(k,214) + .0154005_r8*rxt(k,518)*y(k,224)  &
                      + .0063005_r8*rxt(k,522)*y(k,227)
         mat(k,2077) = .2202005_r8*rxt(k,497)*y(k,6) + .0508005_r8*rxt(k,513)*y(k,110)
         mat(k,43) = .5931005_r8*rxt(k,515)*y(k,217)
         mat(k,49) = .1279005_r8*rxt(k,496)*y(k,124) + .2202005_r8*rxt(k,495)*y(k,203)
         mat(k,55) = .0097005_r8*rxt(k,501)*y(k,124) + .0023005_r8*rxt(k,500)*y(k,203)
         mat(k,1853) = .2202005_r8*rxt(k,495)*y(k,190) + .0023005_r8*rxt(k,500) &
                      *y(k,192) + .0031005_r8*rxt(k,503)*y(k,207)  &
                      + .2381005_r8*rxt(k,507)*y(k,208) + .0508005_r8*rxt(k,511) &
                      *y(k,214) + .1364005_r8*rxt(k,517)*y(k,224)  &
                      + .1677005_r8*rxt(k,521)*y(k,227)
         mat(k,61) = .0003005_r8*rxt(k,504)*y(k,124) + .0031005_r8*rxt(k,503)*y(k,203)
         mat(k,67) = .1056005_r8*rxt(k,508)*y(k,124) + .2381005_r8*rxt(k,507)*y(k,203)
         mat(k,75) = .0245005_r8*rxt(k,512)*y(k,124) + .0508005_r8*rxt(k,511)*y(k,203)
         mat(k,1540) = .5931005_r8*rxt(k,515)*y(k,172)
         mat(k,81) = .0154005_r8*rxt(k,518)*y(k,124) + .1364005_r8*rxt(k,517)*y(k,203)
         mat(k,87) = .0063005_r8*rxt(k,522)*y(k,124) + .1677005_r8*rxt(k,521)*y(k,203)

      end do

      end subroutine     nlnmat06

      subroutine     nlnmat07( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len

         mat(k,985) = .2067005_r8*rxt(k,497)*y(k,134)
         mat(k,936) = .1149005_r8*rxt(k,513)*y(k,134)
         mat(k,1762) = .1792005_r8*rxt(k,496)*y(k,190) + .0034005_r8*rxt(k,501) &
                      *y(k,192) + .0003005_r8*rxt(k,504)*y(k,207)  &
                      + .1026005_r8*rxt(k,508)*y(k,208) + .0082005_r8*rxt(k,512) &
                      *y(k,214) + .0452005_r8*rxt(k,518)*y(k,224)  &
                      + .0237005_r8*rxt(k,522)*y(k,227)
         mat(k,2078) = .2067005_r8*rxt(k,497)*y(k,6) + .1149005_r8*rxt(k,513)*y(k,110)
         mat(k,44) = .1534005_r8*rxt(k,515)*y(k,217)
         mat(k,50) = .1792005_r8*rxt(k,496)*y(k,124) + .2067005_r8*rxt(k,495)*y(k,203)
         mat(k,56) = .0034005_r8*rxt(k,501)*y(k,124) + .0008005_r8*rxt(k,500)*y(k,203)
         mat(k,1854) = .2067005_r8*rxt(k,495)*y(k,190) + .0008005_r8*rxt(k,500) &
                      *y(k,192) + .0035005_r8*rxt(k,503)*y(k,207)  &
                      + .1308005_r8*rxt(k,507)*y(k,208) + .1149005_r8*rxt(k,511) &
                      *y(k,214) + .0101005_r8*rxt(k,517)*y(k,224)  &
                      + .0174005_r8*rxt(k,521)*y(k,227)
         mat(k,62) = .0003005_r8*rxt(k,504)*y(k,124) + .0035005_r8*rxt(k,503)*y(k,203)
         mat(k,68) = .1026005_r8*rxt(k,508)*y(k,124) + .1308005_r8*rxt(k,507)*y(k,203)
         mat(k,76) = .0082005_r8*rxt(k,512)*y(k,124) + .1149005_r8*rxt(k,511)*y(k,203)
         mat(k,1541) = .1534005_r8*rxt(k,515)*y(k,172)
         mat(k,82) = .0452005_r8*rxt(k,518)*y(k,124) + .0101005_r8*rxt(k,517)*y(k,203)
         mat(k,88) = .0237005_r8*rxt(k,522)*y(k,124) + .0174005_r8*rxt(k,521)*y(k,203)


         mat(k,986) = .0653005_r8*rxt(k,497)*y(k,134)
         mat(k,937) = .0348005_r8*rxt(k,513)*y(k,134)
         mat(k,1763) = .0676005_r8*rxt(k,496)*y(k,190) + .1579005_r8*rxt(k,501) &
                      *y(k,192) + .0073005_r8*rxt(k,504)*y(k,207)  &
                      + .0521005_r8*rxt(k,508)*y(k,208) + .0772005_r8*rxt(k,512) &
                      *y(k,214) + .0966005_r8*rxt(k,518)*y(k,224)  &
                      + .0025005_r8*rxt(k,522)*y(k,227)
         mat(k,2079) = .0653005_r8*rxt(k,497)*y(k,6) + .0348005_r8*rxt(k,513)*y(k,110)
         mat(k,45) = .0459005_r8*rxt(k,515)*y(k,217)
         mat(k,51) = .0676005_r8*rxt(k,496)*y(k,124) + .0653005_r8*rxt(k,495)*y(k,203)
         mat(k,57) = .1579005_r8*rxt(k,501)*y(k,124) + .0843005_r8*rxt(k,500)*y(k,203)
         mat(k,1855) = .0653005_r8*rxt(k,495)*y(k,190) + .0843005_r8*rxt(k,500) &
                      *y(k,192) + .0003005_r8*rxt(k,503)*y(k,207)  &
                      + .0348005_r8*rxt(k,507)*y(k,208) + .0348005_r8*rxt(k,511) &
                      *y(k,214) + .0763005_r8*rxt(k,517)*y(k,224) + .086_r8*rxt(k,521) &
                      *y(k,227)
         mat(k,63) = .0073005_r8*rxt(k,504)*y(k,124) + .0003005_r8*rxt(k,503)*y(k,203)
         mat(k,69) = .0521005_r8*rxt(k,508)*y(k,124) + .0348005_r8*rxt(k,507)*y(k,203)
         mat(k,77) = .0772005_r8*rxt(k,512)*y(k,124) + .0348005_r8*rxt(k,511)*y(k,203)
         mat(k,1542) = .0459005_r8*rxt(k,515)*y(k,172)
         mat(k,83) = .0966005_r8*rxt(k,518)*y(k,124) + .0763005_r8*rxt(k,517)*y(k,203)
         mat(k,89) = .0025005_r8*rxt(k,522)*y(k,124) + .086_r8*rxt(k,521)*y(k,203)


         mat(k,987) = .1749305_r8*rxt(k,494)*y(k,126) + .1284005_r8*rxt(k,497) &
                      *y(k,134)
         mat(k,838) = .0590245_r8*rxt(k,502)*y(k,126) + .0033005_r8*rxt(k,505) &
                      *y(k,134)
         mat(k,938) = .1749305_r8*rxt(k,510)*y(k,126) + .0554005_r8*rxt(k,513) &
                      *y(k,134)
         mat(k,1764) = .079_r8*rxt(k,496)*y(k,190) + .0059005_r8*rxt(k,501)*y(k,192)  &
                      + .0057005_r8*rxt(k,504)*y(k,207) + .0143005_r8*rxt(k,508) &
                      *y(k,208) + .0332005_r8*rxt(k,512)*y(k,214)  &
                      + .0073005_r8*rxt(k,518)*y(k,224) + .011_r8*rxt(k,522)*y(k,227)
         mat(k,1704) = .1749305_r8*rxt(k,494)*y(k,6) + .0590245_r8*rxt(k,502)*y(k,98)  &
                      + .1749305_r8*rxt(k,510)*y(k,110)
         mat(k,2080) = .1284005_r8*rxt(k,497)*y(k,6) + .0033005_r8*rxt(k,505)*y(k,98)  &
                      + .0554005_r8*rxt(k,513)*y(k,110)
         mat(k,46) = .0085005_r8*rxt(k,515)*y(k,217)
         mat(k,52) = .079_r8*rxt(k,496)*y(k,124) + .1284005_r8*rxt(k,495)*y(k,203)
         mat(k,58) = .0059005_r8*rxt(k,501)*y(k,124) + .0443005_r8*rxt(k,500)*y(k,203)
         mat(k,1856) = .1284005_r8*rxt(k,495)*y(k,190) + .0443005_r8*rxt(k,500) &
                      *y(k,192) + .0271005_r8*rxt(k,503)*y(k,207)  &
                      + .0076005_r8*rxt(k,507)*y(k,208) + .0554005_r8*rxt(k,511) &
                      *y(k,214) + .2157005_r8*rxt(k,517)*y(k,224)  &
                      + .0512005_r8*rxt(k,521)*y(k,227)
         mat(k,64) = .0057005_r8*rxt(k,504)*y(k,124) + .0271005_r8*rxt(k,503)*y(k,203)
         mat(k,70) = .0143005_r8*rxt(k,508)*y(k,124) + .0076005_r8*rxt(k,507)*y(k,203)
         mat(k,78) = .0332005_r8*rxt(k,512)*y(k,124) + .0554005_r8*rxt(k,511)*y(k,203)
         mat(k,1543) = .0085005_r8*rxt(k,515)*y(k,172)
         mat(k,84) = .0073005_r8*rxt(k,518)*y(k,124) + .2157005_r8*rxt(k,517)*y(k,203)
         mat(k,90) = .011_r8*rxt(k,522)*y(k,124) + .0512005_r8*rxt(k,521)*y(k,203)


         mat(k,988) = .5901905_r8*rxt(k,494)*y(k,126) + .114_r8*rxt(k,497)*y(k,134)
         mat(k,839) = .0250245_r8*rxt(k,502)*y(k,126)
         mat(k,939) = .5901905_r8*rxt(k,510)*y(k,126) + .1278005_r8*rxt(k,513) &
                      *y(k,134)
         mat(k,1765) = .1254005_r8*rxt(k,496)*y(k,190) + .0536005_r8*rxt(k,501) &
                      *y(k,192) + .0623005_r8*rxt(k,504)*y(k,207)  &
                      + .0166005_r8*rxt(k,508)*y(k,208) + .130_r8*rxt(k,512)*y(k,214)  &
                      + .238_r8*rxt(k,518)*y(k,224) + .1185005_r8*rxt(k,522)*y(k,227)
         mat(k,1705) = .5901905_r8*rxt(k,494)*y(k,6) + .0250245_r8*rxt(k,502)*y(k,98)  &
                      + .5901905_r8*rxt(k,510)*y(k,110)
         mat(k,2081) = .114_r8*rxt(k,497)*y(k,6) + .1278005_r8*rxt(k,513)*y(k,110)
         mat(k,47) = .0128005_r8*rxt(k,515)*y(k,217)
         mat(k,53) = .1254005_r8*rxt(k,496)*y(k,124) + .114_r8*rxt(k,495)*y(k,203)
         mat(k,59) = .0536005_r8*rxt(k,501)*y(k,124) + .1621005_r8*rxt(k,500)*y(k,203)
         mat(k,1857) = .114_r8*rxt(k,495)*y(k,190) + .1621005_r8*rxt(k,500)*y(k,192)  &
                      + .0474005_r8*rxt(k,503)*y(k,207) + .0113005_r8*rxt(k,507) &
                      *y(k,208) + .1278005_r8*rxt(k,511)*y(k,214)  &
                      + .0738005_r8*rxt(k,517)*y(k,224) + .1598005_r8*rxt(k,521) &
                      *y(k,227)
         mat(k,65) = .0623005_r8*rxt(k,504)*y(k,124) + .0474005_r8*rxt(k,503)*y(k,203)
         mat(k,71) = .0166005_r8*rxt(k,508)*y(k,124) + .0113005_r8*rxt(k,507)*y(k,203)
         mat(k,79) = .130_r8*rxt(k,512)*y(k,124) + .1278005_r8*rxt(k,511)*y(k,203)
         mat(k,1544) = .0128005_r8*rxt(k,515)*y(k,172)
         mat(k,85) = .238_r8*rxt(k,518)*y(k,124) + .0738005_r8*rxt(k,517)*y(k,203)
         mat(k,91) = .1185005_r8*rxt(k,522)*y(k,124) + .1598005_r8*rxt(k,521)*y(k,203)


         mat(k,48) = -(rxt(k,515)*y(k,217))
         mat(k,1545) = -rxt(k,515)*y(k,172)


         mat(k,189) = .100_r8*rxt(k,430)*y(k,217)
         mat(k,207) = .230_r8*rxt(k,432)*y(k,217)
         mat(k,1565) = .100_r8*rxt(k,430)*y(k,180) + .230_r8*rxt(k,432)*y(k,182)

         mat(k,642) = -(rxt(k,454)*y(k,217))
         mat(k,1629) = -rxt(k,454)*y(k,174)

         mat(k,1896) = rxt(k,452)*y(k,221)
         mat(k,1108) = rxt(k,452)*y(k,203)

         mat(k,617) = -(rxt(k,455)*y(k,217))
         mat(k,1626) = -rxt(k,455)*y(k,175)

         mat(k,1793) = .200_r8*rxt(k,448)*y(k,215) + .200_r8*rxt(k,458)*y(k,222)
         mat(k,2029) = .500_r8*rxt(k,446)*y(k,215)
         mat(k,1042) = .200_r8*rxt(k,448)*y(k,124) + .500_r8*rxt(k,446)*y(k,198)
         mat(k,910) = .200_r8*rxt(k,458)*y(k,124)

         mat(k,474) = -(rxt(k,459)*y(k,217))
         mat(k,1608) = -rxt(k,459)*y(k,176)

         mat(k,1887) = rxt(k,457)*y(k,222)
         mat(k,909) = rxt(k,457)*y(k,203)

         mat(k,969) = -(rxt(k,460)*y(k,126) + rxt(k,461)*y(k,217))
         mat(k,1717) = -rxt(k,460)*y(k,177)
         mat(k,1659) = -rxt(k,461)*y(k,177)

         mat(k,997) = .330_r8*rxt(k,441)*y(k,134)
         mat(k,949) = .330_r8*rxt(k,444)*y(k,134)
         mat(k,1812) = .800_r8*rxt(k,448)*y(k,215) + .800_r8*rxt(k,458)*y(k,222)
         mat(k,1717) = mat(k,1717) + rxt(k,449)*y(k,215)
         mat(k,2097) = .330_r8*rxt(k,441)*y(k,6) + .330_r8*rxt(k,444)*y(k,110)
         mat(k,618) = rxt(k,455)*y(k,217)
         mat(k,2038) = .500_r8*rxt(k,446)*y(k,215) + rxt(k,456)*y(k,222)
         mat(k,1044) = .800_r8*rxt(k,448)*y(k,124) + rxt(k,449)*y(k,126)  &
                      + .500_r8*rxt(k,446)*y(k,198)
         mat(k,1659) = mat(k,1659) + rxt(k,455)*y(k,175)
         mat(k,914) = .800_r8*rxt(k,458)*y(k,124) + rxt(k,456)*y(k,198)

         mat(k,1059) = -(rxt(k,462)*y(k,217))
         mat(k,1664) = -rxt(k,462)*y(k,178)

         mat(k,1001) = .300_r8*rxt(k,441)*y(k,134)
         mat(k,952) = .300_r8*rxt(k,444)*y(k,134)
         mat(k,1815) = .900_r8*rxt(k,453)*y(k,221)
         mat(k,2102) = .300_r8*rxt(k,441)*y(k,6) + .300_r8*rxt(k,444)*y(k,110)
         mat(k,2041) = rxt(k,451)*y(k,221)
         mat(k,1112) = .900_r8*rxt(k,453)*y(k,124) + rxt(k,451)*y(k,198)

         mat(k,655) = -(rxt(k,429)*y(k,217))
         mat(k,1630) = -rxt(k,429)*y(k,179)

         mat(k,1897) = rxt(k,427)*y(k,223)
         mat(k,730) = rxt(k,427)*y(k,203)

         mat(k,187) = -(rxt(k,430)*y(k,217))
         mat(k,1563) = -rxt(k,430)*y(k,180)

         mat(k,203) = -(rxt(k,396)*y(k,217))
         mat(k,1566) = -rxt(k,396)*y(k,181)

         mat(k,1866) = rxt(k,393)*y(k,225)
         mat(k,1166) = rxt(k,393)*y(k,203)

         mat(k,208) = -(rxt(k,432)*y(k,217))
         mat(k,1567) = -rxt(k,432)*y(k,182)

         mat(k,701) = -(rxt(k,435)*y(k,217))
         mat(k,1635) = -rxt(k,435)*y(k,183)

         mat(k,1901) = rxt(k,433)*y(k,226)
         mat(k,746) = rxt(k,433)*y(k,203)

         mat(k,216) = -(rxt(k,438)*y(k,217))
         mat(k,1568) = -rxt(k,438)*y(k,184)

         mat(k,209) = .150_r8*rxt(k,432)*y(k,217)
         mat(k,1568) = mat(k,1568) + .150_r8*rxt(k,432)*y(k,182)

         mat(k,427) = -(rxt(k,439)*y(k,217))
         mat(k,1601) = -rxt(k,439)*y(k,185)

         mat(k,1881) = rxt(k,436)*y(k,228)
         mat(k,498) = rxt(k,436)*y(k,203)

         mat(k,518) = -(rxt(k,397)*y(k,203) + rxt(k,398)*y(k,124) + rxt(k,426) &
                      *y(k,125))
         mat(k,1890) = -rxt(k,397)*y(k,188)
         mat(k,1788) = -rxt(k,398)*y(k,188)
         mat(k,2167) = -rxt(k,426)*y(k,188)

         mat(k,254) = rxt(k,403)*y(k,217)
         mat(k,1614) = rxt(k,403)*y(k,22)

         mat(k,899) = -(rxt(k,358)*y(k,203) + (rxt(k,359) + rxt(k,360)) * y(k,124))
         mat(k,1916) = -rxt(k,358)*y(k,189)
         mat(k,1808) = -(rxt(k,359) + rxt(k,360)) * y(k,189)

         mat(k,635) = rxt(k,361)*y(k,217)
         mat(k,239) = rxt(k,362)*y(k,217)
         mat(k,1654) = rxt(k,361)*y(k,2) + rxt(k,362)*y(k,15)

         mat(k,54) = -(rxt(k,495)*y(k,203) + rxt(k,496)*y(k,124))
         mat(k,1858) = -rxt(k,495)*y(k,190)
         mat(k,1766) = -rxt(k,496)*y(k,190)

         mat(k,989) = rxt(k,498)*y(k,217)
         mat(k,1546) = rxt(k,498)*y(k,6)

         mat(k,483) = -(rxt(k,400)*y(k,203) + rxt(k,401)*y(k,124))
         mat(k,1888) = -rxt(k,400)*y(k,191)
         mat(k,1784) = -rxt(k,401)*y(k,191)

         mat(k,167) = .350_r8*rxt(k,399)*y(k,217)
         mat(k,423) = rxt(k,402)*y(k,217)
         mat(k,1609) = .350_r8*rxt(k,399)*y(k,7) + rxt(k,402)*y(k,8)

         mat(k,60) = -(rxt(k,500)*y(k,203) + rxt(k,501)*y(k,124))
         mat(k,1859) = -rxt(k,500)*y(k,192)
         mat(k,1767) = -rxt(k,501)*y(k,192)

         mat(k,163) = rxt(k,499)*y(k,217)
         mat(k,1547) = rxt(k,499)*y(k,7)

         mat(k,435) = -(rxt(k,404)*y(k,203) + rxt(k,406)*y(k,124))
         mat(k,1882) = -rxt(k,404)*y(k,193)
         mat(k,1779) = -rxt(k,406)*y(k,193)

         mat(k,340) = rxt(k,405)*y(k,217)
         mat(k,190) = .070_r8*rxt(k,430)*y(k,217)
         mat(k,210) = .060_r8*rxt(k,432)*y(k,217)
         mat(k,1602) = rxt(k,405)*y(k,23) + .070_r8*rxt(k,430)*y(k,180)  &
                      + .060_r8*rxt(k,432)*y(k,182)

         mat(k,815) = -(4._r8*rxt(k,281)*y(k,194) + rxt(k,282)*y(k,198) + rxt(k,283) &
                      *y(k,203) + rxt(k,284)*y(k,124))
         mat(k,2034) = -rxt(k,282)*y(k,194)
         mat(k,1912) = -rxt(k,283)*y(k,194)
         mat(k,1805) = -rxt(k,284)*y(k,194)

         mat(k,345) = .500_r8*rxt(k,286)*y(k,217)
         mat(k,292) = rxt(k,287)*y(k,56) + rxt(k,288)*y(k,217)
         mat(k,1996) = rxt(k,287)*y(k,28)
         mat(k,1646) = .500_r8*rxt(k,286)*y(k,27) + rxt(k,288)*y(k,28)

      end do

      end subroutine     nlnmat07

      subroutine     nlnmat08( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,791) = -(rxt(k,310)*y(k,198) + rxt(k,311)*y(k,203) + rxt(k,312) &
                      *y(k,124))
         mat(k,2032) = -rxt(k,310)*y(k,195)
         mat(k,1909) = -rxt(k,311)*y(k,195)
         mat(k,1803) = -rxt(k,312)*y(k,195)

         mat(k,416) = rxt(k,313)*y(k,217)
         mat(k,110) = rxt(k,314)*y(k,217)
         mat(k,1643) = rxt(k,313)*y(k,30) + rxt(k,314)*y(k,31)

         mat(k,625) = -(rxt(k,407)*y(k,203) + rxt(k,408)*y(k,124))
         mat(k,1894) = -rxt(k,407)*y(k,196)
         mat(k,1794) = -rxt(k,408)*y(k,196)

         mat(k,267) = rxt(k,409)*y(k,217)
         mat(k,1794) = mat(k,1794) + rxt(k,398)*y(k,188)
         mat(k,2087) = rxt(k,424)*y(k,141)
         mat(k,464) = rxt(k,424)*y(k,134)
         mat(k,519) = rxt(k,398)*y(k,124) + .400_r8*rxt(k,397)*y(k,203)
         mat(k,1894) = mat(k,1894) + .400_r8*rxt(k,397)*y(k,188)
         mat(k,1627) = rxt(k,409)*y(k,32)

         mat(k,1386) = -(4._r8*rxt(k,292)*y(k,197) + rxt(k,293)*y(k,198) + rxt(k,294) &
                      *y(k,203) + rxt(k,295)*y(k,124) + rxt(k,306)*y(k,125) + rxt(k,333) &
                      *y(k,209) + rxt(k,366)*y(k,205) + rxt(k,371)*y(k,206) + rxt(k,380) &
                      *y(k,101) + rxt(k,391)*y(k,225))
         mat(k,2058) = -rxt(k,293)*y(k,197)
         mat(k,1939) = -rxt(k,294)*y(k,197)
         mat(k,1833) = -rxt(k,295)*y(k,197)
         mat(k,2185) = -rxt(k,306)*y(k,197)
         mat(k,1313) = -rxt(k,333)*y(k,197)
         mat(k,1260) = -rxt(k,366)*y(k,197)
         mat(k,1292) = -rxt(k,371)*y(k,197)
         mat(k,1197) = -rxt(k,380)*y(k,197)
         mat(k,1175) = -rxt(k,391)*y(k,197)

         mat(k,1007) = .060_r8*rxt(k,441)*y(k,134)
         mat(k,1087) = rxt(k,289)*y(k,126) + rxt(k,290)*y(k,217)
         mat(k,1222) = rxt(k,315)*y(k,126) + rxt(k,316)*y(k,217)
         mat(k,612) = .500_r8*rxt(k,297)*y(k,217)
         mat(k,850) = .080_r8*rxt(k,386)*y(k,134)
         mat(k,1213) = .100_r8*rxt(k,339)*y(k,134)
         mat(k,957) = .060_r8*rxt(k,444)*y(k,134)
         mat(k,1334) = .280_r8*rxt(k,353)*y(k,134)
         mat(k,1833) = mat(k,1833) + .530_r8*rxt(k,337)*y(k,209) + rxt(k,346)*y(k,211)  &
                      + rxt(k,349)*y(k,213) + rxt(k,324)*y(k,220)
         mat(k,1741) = rxt(k,289)*y(k,45) + rxt(k,315)*y(k,49) + .530_r8*rxt(k,336) &
                      *y(k,209) + rxt(k,347)*y(k,211)
         mat(k,2118) = .060_r8*rxt(k,441)*y(k,6) + .080_r8*rxt(k,386)*y(k,98)  &
                      + .100_r8*rxt(k,339)*y(k,105) + .060_r8*rxt(k,444)*y(k,110)  &
                      + .280_r8*rxt(k,353)*y(k,111)
         mat(k,1062) = .650_r8*rxt(k,462)*y(k,217)
         mat(k,1386) = mat(k,1386) + .530_r8*rxt(k,333)*y(k,209)
         mat(k,2058) = mat(k,2058) + .260_r8*rxt(k,334)*y(k,209) + rxt(k,343)*y(k,211)  &
                      + .300_r8*rxt(k,322)*y(k,220)
         mat(k,1939) = mat(k,1939) + .450_r8*rxt(k,344)*y(k,211) + .200_r8*rxt(k,348) &
                      *y(k,213) + .150_r8*rxt(k,323)*y(k,220)
         mat(k,1313) = mat(k,1313) + .530_r8*rxt(k,337)*y(k,124) + .530_r8*rxt(k,336) &
                      *y(k,126) + .530_r8*rxt(k,333)*y(k,197) + .260_r8*rxt(k,334) &
                      *y(k,198)
         mat(k,1355) = rxt(k,346)*y(k,124) + rxt(k,347)*y(k,126) + rxt(k,343)*y(k,198)  &
                      + .450_r8*rxt(k,344)*y(k,203) + 4.000_r8*rxt(k,345)*y(k,211)
         mat(k,679) = rxt(k,349)*y(k,124) + .200_r8*rxt(k,348)*y(k,203)
         mat(k,1683) = rxt(k,290)*y(k,45) + rxt(k,316)*y(k,49) + .500_r8*rxt(k,297) &
                      *y(k,51) + .650_r8*rxt(k,462)*y(k,178)
         mat(k,1132) = rxt(k,324)*y(k,124) + .300_r8*rxt(k,322)*y(k,198)  &
                      + .150_r8*rxt(k,323)*y(k,203)

         mat(k,2070) = -(rxt(k,183)*y(k,59) + (4._r8*rxt(k,260) + 4._r8*rxt(k,261) &
                      ) * y(k,198) + rxt(k,262)*y(k,203) + rxt(k,263)*y(k,124) &
                      + rxt(k,282)*y(k,194) + rxt(k,293)*y(k,197) + rxt(k,310) &
                      *y(k,195) + rxt(k,322)*y(k,220) + rxt(k,334)*y(k,209) + rxt(k,343) &
                      *y(k,211) + rxt(k,367)*y(k,205) + rxt(k,372)*y(k,206) + rxt(k,381) &
                      *y(k,101) + rxt(k,392)*y(k,225) + rxt(k,446)*y(k,215) + rxt(k,451) &
                      *y(k,221) + rxt(k,456)*y(k,222))
         mat(k,1979) = -rxt(k,183)*y(k,198)
         mat(k,1953) = -rxt(k,262)*y(k,198)
         mat(k,1846) = -rxt(k,263)*y(k,198)
         mat(k,821) = -rxt(k,282)*y(k,198)
         mat(k,1394) = -rxt(k,293)*y(k,198)
         mat(k,798) = -rxt(k,310)*y(k,198)
         mat(k,1137) = -rxt(k,322)*y(k,198)
         mat(k,1320) = -rxt(k,334)*y(k,198)
         mat(k,1362) = -rxt(k,343)*y(k,198)
         mat(k,1267) = -rxt(k,367)*y(k,198)
         mat(k,1299) = -rxt(k,372)*y(k,198)
         mat(k,1204) = -rxt(k,381)*y(k,198)
         mat(k,1181) = -rxt(k,392)*y(k,198)
         mat(k,1055) = -rxt(k,446)*y(k,198)
         mat(k,1123) = -rxt(k,451)*y(k,198)
         mat(k,920) = -rxt(k,456)*y(k,198)

         mat(k,1036) = .280_r8*rxt(k,309)*y(k,134)
         mat(k,687) = rxt(k,296)*y(k,217)
         mat(k,389) = .700_r8*rxt(k,265)*y(k,217)
         mat(k,1439) = rxt(k,177)*y(k,56) + rxt(k,233)*y(k,73) + rxt(k,272)*y(k,216)  &
                      + rxt(k,266)*y(k,217)
         mat(k,2018) = rxt(k,177)*y(k,54)
         mat(k,871) = rxt(k,233)*y(k,54)
         mat(k,855) = .050_r8*rxt(k,386)*y(k,134)
         mat(k,1204) = mat(k,1204) + rxt(k,380)*y(k,197)
         mat(k,1846) = mat(k,1846) + rxt(k,295)*y(k,197) + .830_r8*rxt(k,412)*y(k,199)  &
                      + .170_r8*rxt(k,418)*y(k,212)
         mat(k,2131) = .280_r8*rxt(k,309)*y(k,29) + .050_r8*rxt(k,386)*y(k,98)
         mat(k,1394) = mat(k,1394) + rxt(k,380)*y(k,101) + rxt(k,295)*y(k,124)  &
                      + 4.000_r8*rxt(k,292)*y(k,197) + .900_r8*rxt(k,293)*y(k,198)  &
                      + .450_r8*rxt(k,294)*y(k,203) + rxt(k,366)*y(k,205) + rxt(k,371) &
                      *y(k,206) + rxt(k,333)*y(k,209) + rxt(k,342)*y(k,211)  &
                      + rxt(k,391)*y(k,225)
         mat(k,2070) = mat(k,2070) + .900_r8*rxt(k,293)*y(k,197)
         mat(k,765) = .830_r8*rxt(k,412)*y(k,124) + .330_r8*rxt(k,411)*y(k,203)
         mat(k,1953) = mat(k,1953) + .450_r8*rxt(k,294)*y(k,197) + .330_r8*rxt(k,411) &
                      *y(k,199) + .070_r8*rxt(k,417)*y(k,212)
         mat(k,1267) = mat(k,1267) + rxt(k,366)*y(k,197)
         mat(k,1299) = mat(k,1299) + rxt(k,371)*y(k,197)
         mat(k,1320) = mat(k,1320) + rxt(k,333)*y(k,197)
         mat(k,1362) = mat(k,1362) + rxt(k,342)*y(k,197)
         mat(k,880) = .170_r8*rxt(k,418)*y(k,124) + .070_r8*rxt(k,417)*y(k,203)
         mat(k,1533) = rxt(k,272)*y(k,54)
         mat(k,1697) = rxt(k,296)*y(k,50) + .700_r8*rxt(k,265)*y(k,53) + rxt(k,266) &
                      *y(k,54)
         mat(k,1181) = mat(k,1181) + rxt(k,391)*y(k,197)

         mat(k,759) = -(rxt(k,411)*y(k,203) + rxt(k,412)*y(k,124) + rxt(k,413) &
                      *y(k,125))
         mat(k,1906) = -rxt(k,411)*y(k,199)
         mat(k,1801) = -rxt(k,412)*y(k,199)
         mat(k,2173) = -rxt(k,413)*y(k,199)

         mat(k,564) = -((rxt(k,330) + rxt(k,331)) * y(k,124))
         mat(k,1790) = -(rxt(k,330) + rxt(k,331)) * y(k,200)

         mat(k,350) = rxt(k,329)*y(k,217)
         mat(k,1619) = rxt(k,329)*y(k,16)


         mat(k,1775) = .750_r8*rxt(k,299)*y(k,202)
         mat(k,713) = .750_r8*rxt(k,299)*y(k,124)

         mat(k,714) = -(rxt(k,298)*y(k,203) + rxt(k,299)*y(k,124))
         mat(k,1902) = -rxt(k,298)*y(k,202)
         mat(k,1797) = -rxt(k,299)*y(k,202)

         mat(k,549) = rxt(k,305)*y(k,217)
         mat(k,1636) = rxt(k,305)*y(k,25)

         mat(k,1950) = -((rxt(k,136) + rxt(k,137) + rxt(k,138)) * y(k,76) + rxt(k,140) &
                      *y(k,133) + rxt(k,141)*y(k,134) + rxt(k,145)*y(k,217) &
                      + 4._r8*rxt(k,150)*y(k,203) + rxt(k,160)*y(k,126) + rxt(k,165) &
                      *y(k,124) + rxt(k,170)*y(k,125) + (rxt(k,180) + rxt(k,181) &
                      ) * y(k,56) + rxt(k,187)*y(k,59) + rxt(k,213)*y(k,17) + rxt(k,219) &
                      *y(k,19) + rxt(k,256)*y(k,42) + rxt(k,262)*y(k,198) + rxt(k,269) &
                      *y(k,204) + rxt(k,283)*y(k,194) + rxt(k,294)*y(k,197) + rxt(k,298) &
                      *y(k,202) + rxt(k,311)*y(k,195) + rxt(k,319)*y(k,219) + rxt(k,323) &
                      *y(k,220) + rxt(k,335)*y(k,209) + rxt(k,344)*y(k,211) + rxt(k,348) &
                      *y(k,213) + rxt(k,358)*y(k,189) + rxt(k,368)*y(k,205) + rxt(k,373) &
                      *y(k,206) + rxt(k,382)*y(k,101) + rxt(k,393)*y(k,225) + rxt(k,397) &
                      *y(k,188) + rxt(k,400)*y(k,191) + rxt(k,404)*y(k,193) + rxt(k,407) &
                      *y(k,196) + rxt(k,411)*y(k,199) + rxt(k,414)*y(k,210) + rxt(k,417) &
                      *y(k,212) + rxt(k,420)*y(k,218) + rxt(k,427)*y(k,223) + rxt(k,433) &
                      *y(k,226) + rxt(k,436)*y(k,228) + rxt(k,447)*y(k,215) + rxt(k,452) &
                      *y(k,221) + rxt(k,457)*y(k,222))
         mat(k,1468) = -(rxt(k,136) + rxt(k,137) + rxt(k,138)) * y(k,203)
         mat(k,2250) = -rxt(k,140)*y(k,203)
         mat(k,2128) = -rxt(k,141)*y(k,203)
         mat(k,1694) = -rxt(k,145)*y(k,203)
         mat(k,1751) = -rxt(k,160)*y(k,203)
         mat(k,1843) = -rxt(k,165)*y(k,203)
         mat(k,2195) = -rxt(k,170)*y(k,203)
         mat(k,2015) = -(rxt(k,180) + rxt(k,181)) * y(k,203)
         mat(k,1976) = -rxt(k,187)*y(k,203)
         mat(k,1420) = -rxt(k,213)*y(k,203)
         mat(k,2219) = -rxt(k,219)*y(k,203)
         mat(k,1490) = -rxt(k,256)*y(k,203)
         mat(k,2067) = -rxt(k,262)*y(k,203)
         mat(k,445) = -rxt(k,269)*y(k,203)
         mat(k,820) = -rxt(k,283)*y(k,203)
         mat(k,1393) = -rxt(k,294)*y(k,203)
         mat(k,719) = -rxt(k,298)*y(k,203)
         mat(k,797) = -rxt(k,311)*y(k,203)
         mat(k,774) = -rxt(k,319)*y(k,203)
         mat(k,1136) = -rxt(k,323)*y(k,203)
         mat(k,1319) = -rxt(k,335)*y(k,203)
         mat(k,1361) = -rxt(k,344)*y(k,203)
         mat(k,682) = -rxt(k,348)*y(k,203)
         mat(k,906) = -rxt(k,358)*y(k,203)
         mat(k,1266) = -rxt(k,368)*y(k,203)
         mat(k,1298) = -rxt(k,373)*y(k,203)
         mat(k,1203) = -rxt(k,382)*y(k,203)
         mat(k,1180) = -rxt(k,393)*y(k,203)
         mat(k,522) = -rxt(k,397)*y(k,203)
         mat(k,488) = -rxt(k,400)*y(k,203)
         mat(k,439) = -rxt(k,404)*y(k,203)
         mat(k,628) = -rxt(k,407)*y(k,203)
         mat(k,764) = -rxt(k,411)*y(k,203)
         mat(k,725) = -rxt(k,414)*y(k,203)
         mat(k,879) = -rxt(k,417)*y(k,203)
         mat(k,458) = -rxt(k,420)*y(k,203)
         mat(k,740) = -rxt(k,427)*y(k,203)
         mat(k,757) = -rxt(k,433)*y(k,203)
         mat(k,504) = -rxt(k,436)*y(k,203)
         mat(k,1054) = -rxt(k,447)*y(k,203)
         mat(k,1122) = -rxt(k,452)*y(k,203)
         mat(k,919) = -rxt(k,457)*y(k,203)

         mat(k,1013) = .570_r8*rxt(k,441)*y(k,134)
         mat(k,169) = .650_r8*rxt(k,399)*y(k,217)
         mat(k,1420) = mat(k,1420) + rxt(k,212)*y(k,42)
         mat(k,2219) = mat(k,2219) + rxt(k,224)*y(k,217)
         mat(k,290) = .350_r8*rxt(k,278)*y(k,217)
         mat(k,554) = .130_r8*rxt(k,280)*y(k,134)
         mat(k,264) = rxt(k,285)*y(k,217)
         mat(k,1035) = .280_r8*rxt(k,309)*y(k,134)
         mat(k,1490) = mat(k,1490) + rxt(k,212)*y(k,17) + rxt(k,176)*y(k,56)  &
                      + rxt(k,257)*y(k,126) + rxt(k,258)*y(k,133)
         mat(k,598) = rxt(k,241)*y(k,56) + rxt(k,242)*y(k,217)
         mat(k,368) = rxt(k,244)*y(k,56) + rxt(k,245)*y(k,217)
         mat(k,104) = rxt(k,291)*y(k,217)
         mat(k,789) = rxt(k,264)*y(k,217)
         mat(k,1437) = rxt(k,273)*y(k,216)
         mat(k,2015) = mat(k,2015) + rxt(k,176)*y(k,42) + rxt(k,241)*y(k,43)  &
                      + rxt(k,244)*y(k,46) + rxt(k,179)*y(k,79)
         mat(k,1976) = mat(k,1976) + rxt(k,183)*y(k,198) + rxt(k,194)*y(k,217)
         mat(k,1105) = rxt(k,276)*y(k,217)
         mat(k,198) = .730_r8*rxt(k,410)*y(k,217)
         mat(k,302) = .500_r8*rxt(k,478)*y(k,217)
         mat(k,1100) = rxt(k,302)*y(k,217)
         mat(k,982) = rxt(k,303)*y(k,217)
         mat(k,605) = rxt(k,179)*y(k,56) + rxt(k,135)*y(k,133) + rxt(k,144)*y(k,217)
         mat(k,182) = rxt(k,267)*y(k,217)
         mat(k,932) = rxt(k,268)*y(k,217)
         mat(k,1161) = rxt(k,332)*y(k,217)
         mat(k,1145) = rxt(k,317)*y(k,217)
         mat(k,854) = .370_r8*rxt(k,386)*y(k,134)
         mat(k,588) = .300_r8*rxt(k,377)*y(k,217)
         mat(k,563) = rxt(k,378)*y(k,217)
         mat(k,1203) = mat(k,1203) + rxt(k,383)*y(k,124) + rxt(k,384)*y(k,126)  &
                      + rxt(k,380)*y(k,197) + 1.200_r8*rxt(k,381)*y(k,198)
         mat(k,401) = rxt(k,385)*y(k,217)
         mat(k,1216) = .140_r8*rxt(k,339)*y(k,134)
         mat(k,314) = .200_r8*rxt(k,341)*y(k,217)
         mat(k,579) = .500_r8*rxt(k,352)*y(k,217)
         mat(k,963) = .570_r8*rxt(k,444)*y(k,134)
         mat(k,1341) = .280_r8*rxt(k,353)*y(k,134)
         mat(k,378) = rxt(k,389)*y(k,217)
         mat(k,1081) = rxt(k,390)*y(k,217)
         mat(k,1843) = mat(k,1843) + rxt(k,383)*y(k,101) + rxt(k,359)*y(k,189)  &
                      + rxt(k,401)*y(k,191) + rxt(k,406)*y(k,193) + rxt(k,284) &
                      *y(k,194) + rxt(k,312)*y(k,195) + rxt(k,263)*y(k,198)  &
                      + .170_r8*rxt(k,412)*y(k,199) + rxt(k,330)*y(k,200)  &
                      + .250_r8*rxt(k,299)*y(k,202) + rxt(k,271)*y(k,204)  &
                      + .920_r8*rxt(k,369)*y(k,205) + .920_r8*rxt(k,375)*y(k,206)  &
                      + .470_r8*rxt(k,337)*y(k,209) + .400_r8*rxt(k,415)*y(k,210)  &
                      + .830_r8*rxt(k,418)*y(k,212) + rxt(k,421)*y(k,218) + rxt(k,320) &
                      *y(k,219) + .900_r8*rxt(k,453)*y(k,221) + .800_r8*rxt(k,458) &
                      *y(k,222) + rxt(k,428)*y(k,223) + rxt(k,394)*y(k,225)  &
                      + rxt(k,434)*y(k,226) + rxt(k,437)*y(k,228)
         mat(k,1751) = mat(k,1751) + rxt(k,257)*y(k,42) + rxt(k,384)*y(k,101)  &
                      + rxt(k,370)*y(k,205) + rxt(k,376)*y(k,206) + .470_r8*rxt(k,336) &
                      *y(k,209) + rxt(k,163)*y(k,217) + rxt(k,395)*y(k,225)
         mat(k,2250) = mat(k,2250) + rxt(k,258)*y(k,42) + rxt(k,135)*y(k,79)
         mat(k,2128) = mat(k,2128) + .570_r8*rxt(k,441)*y(k,6) + .130_r8*rxt(k,280) &
                      *y(k,25) + .280_r8*rxt(k,309)*y(k,29) + .370_r8*rxt(k,386) &
                      *y(k,98) + .140_r8*rxt(k,339)*y(k,105) + .570_r8*rxt(k,444) &
                      *y(k,110) + .280_r8*rxt(k,353)*y(k,111) + rxt(k,147)*y(k,217)
         mat(k,178) = .800_r8*rxt(k,422)*y(k,217)
         mat(k,835) = rxt(k,468)*y(k,217)
         mat(k,1065) = .200_r8*rxt(k,462)*y(k,217)
         mat(k,193) = .280_r8*rxt(k,430)*y(k,217)
         mat(k,215) = .380_r8*rxt(k,432)*y(k,217)
         mat(k,220) = .630_r8*rxt(k,438)*y(k,217)
         mat(k,906) = mat(k,906) + rxt(k,359)*y(k,124)
         mat(k,488) = mat(k,488) + rxt(k,401)*y(k,124)
         mat(k,439) = mat(k,439) + rxt(k,406)*y(k,124)
         mat(k,820) = mat(k,820) + rxt(k,284)*y(k,124) + 2.400_r8*rxt(k,281)*y(k,194)  &
                      + rxt(k,282)*y(k,198)
         mat(k,797) = mat(k,797) + rxt(k,312)*y(k,124) + rxt(k,310)*y(k,198)
         mat(k,1393) = mat(k,1393) + rxt(k,380)*y(k,101) + .900_r8*rxt(k,293)*y(k,198)  &
                      + rxt(k,366)*y(k,205) + rxt(k,371)*y(k,206) + .470_r8*rxt(k,333) &
                      *y(k,209) + rxt(k,391)*y(k,225)
         mat(k,2067) = mat(k,2067) + rxt(k,183)*y(k,59) + 1.200_r8*rxt(k,381)*y(k,101)  &
                      + rxt(k,263)*y(k,124) + rxt(k,282)*y(k,194) + rxt(k,310) &
                      *y(k,195) + .900_r8*rxt(k,293)*y(k,197) + 4.000_r8*rxt(k,260) &
                      *y(k,198) + rxt(k,367)*y(k,205) + rxt(k,372)*y(k,206)  &
                      + .730_r8*rxt(k,334)*y(k,209) + rxt(k,343)*y(k,211)  &
                      + .500_r8*rxt(k,446)*y(k,215) + .300_r8*rxt(k,322)*y(k,220)  &
                      + rxt(k,451)*y(k,221) + rxt(k,456)*y(k,222) + .800_r8*rxt(k,392) &
                      *y(k,225)
         mat(k,764) = mat(k,764) + .170_r8*rxt(k,412)*y(k,124) + .070_r8*rxt(k,411) &
                      *y(k,203)
         mat(k,570) = rxt(k,330)*y(k,124)
         mat(k,719) = mat(k,719) + .250_r8*rxt(k,299)*y(k,124)
         mat(k,1950) = mat(k,1950) + .070_r8*rxt(k,411)*y(k,199) + .160_r8*rxt(k,414) &
                      *y(k,210) + .330_r8*rxt(k,417)*y(k,212)
         mat(k,445) = mat(k,445) + rxt(k,271)*y(k,124)
         mat(k,1266) = mat(k,1266) + .920_r8*rxt(k,369)*y(k,124) + rxt(k,370)*y(k,126)  &
                      + rxt(k,366)*y(k,197) + rxt(k,367)*y(k,198)
         mat(k,1298) = mat(k,1298) + .920_r8*rxt(k,375)*y(k,124) + rxt(k,376)*y(k,126)  &
                      + rxt(k,371)*y(k,197) + rxt(k,372)*y(k,198)
         mat(k,1319) = mat(k,1319) + .470_r8*rxt(k,337)*y(k,124) + .470_r8*rxt(k,336) &
                      *y(k,126) + .470_r8*rxt(k,333)*y(k,197) + .730_r8*rxt(k,334) &
                      *y(k,198)
         mat(k,725) = mat(k,725) + .400_r8*rxt(k,415)*y(k,124) + .160_r8*rxt(k,414) &
                      *y(k,203)
         mat(k,1361) = mat(k,1361) + rxt(k,343)*y(k,198)
         mat(k,879) = mat(k,879) + .830_r8*rxt(k,418)*y(k,124) + .330_r8*rxt(k,417) &
                      *y(k,203)
         mat(k,1054) = mat(k,1054) + .500_r8*rxt(k,446)*y(k,198)
         mat(k,1530) = rxt(k,273)*y(k,54)
         mat(k,1694) = mat(k,1694) + .650_r8*rxt(k,399)*y(k,7) + rxt(k,224)*y(k,19)  &
                      + .350_r8*rxt(k,278)*y(k,24) + rxt(k,285)*y(k,26) + rxt(k,242) &
                      *y(k,43) + rxt(k,245)*y(k,46) + rxt(k,291)*y(k,47) + rxt(k,264) &
                      *y(k,52) + rxt(k,194)*y(k,59) + rxt(k,276)*y(k,62)  &
                      + .730_r8*rxt(k,410)*y(k,66) + .500_r8*rxt(k,478)*y(k,67)  &
                      + rxt(k,302)*y(k,74) + rxt(k,303)*y(k,75) + rxt(k,144)*y(k,79)  &
                      + rxt(k,267)*y(k,86) + rxt(k,268)*y(k,87) + rxt(k,332)*y(k,93)  &
                      + rxt(k,317)*y(k,95) + .300_r8*rxt(k,377)*y(k,99) + rxt(k,378) &
                      *y(k,100) + rxt(k,385)*y(k,102) + .200_r8*rxt(k,341)*y(k,106)  &
                      + .500_r8*rxt(k,352)*y(k,109) + rxt(k,389)*y(k,115) + rxt(k,390) &
                      *y(k,116) + rxt(k,163)*y(k,126) + rxt(k,147)*y(k,134)  &
                      + .800_r8*rxt(k,422)*y(k,142) + rxt(k,468)*y(k,151)  &
                      + .200_r8*rxt(k,462)*y(k,178) + .280_r8*rxt(k,430)*y(k,180)  &
                      + .380_r8*rxt(k,432)*y(k,182) + .630_r8*rxt(k,438)*y(k,184)
         mat(k,458) = mat(k,458) + rxt(k,421)*y(k,124)
         mat(k,774) = mat(k,774) + rxt(k,320)*y(k,124)
         mat(k,1136) = mat(k,1136) + .300_r8*rxt(k,322)*y(k,198)
         mat(k,1122) = mat(k,1122) + .900_r8*rxt(k,453)*y(k,124) + rxt(k,451)*y(k,198)
         mat(k,919) = mat(k,919) + .800_r8*rxt(k,458)*y(k,124) + rxt(k,456)*y(k,198)
         mat(k,740) = mat(k,740) + rxt(k,428)*y(k,124)
         mat(k,1180) = mat(k,1180) + rxt(k,394)*y(k,124) + rxt(k,395)*y(k,126)  &
                      + rxt(k,391)*y(k,197) + .800_r8*rxt(k,392)*y(k,198)
         mat(k,757) = mat(k,757) + rxt(k,434)*y(k,124)
         mat(k,504) = mat(k,504) + rxt(k,437)*y(k,124)

      end do

      end subroutine     nlnmat08

      subroutine     nlnmat09( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,441) = -(rxt(k,269)*y(k,203) + rxt(k,271)*y(k,124))
         mat(k,1883) = -rxt(k,269)*y(k,204)
         mat(k,1780) = -rxt(k,271)*y(k,204)

         mat(k,1476) = rxt(k,256)*y(k,203)
         mat(k,1883) = mat(k,1883) + rxt(k,256)*y(k,42)

         mat(k,1256) = -(rxt(k,366)*y(k,197) + rxt(k,367)*y(k,198) + rxt(k,368) &
                      *y(k,203) + rxt(k,369)*y(k,124) + rxt(k,370)*y(k,126))
         mat(k,1381) = -rxt(k,366)*y(k,205)
         mat(k,2053) = -rxt(k,367)*y(k,205)
         mat(k,1934) = -rxt(k,368)*y(k,205)
         mat(k,1828) = -rxt(k,369)*y(k,205)
         mat(k,1736) = -rxt(k,370)*y(k,205)

         mat(k,847) = .600_r8*rxt(k,387)*y(k,217)
         mat(k,1678) = .600_r8*rxt(k,387)*y(k,98)

         mat(k,1288) = -(rxt(k,371)*y(k,197) + rxt(k,372)*y(k,198) + rxt(k,373) &
                      *y(k,203) + rxt(k,375)*y(k,124) + rxt(k,376)*y(k,126))
         mat(k,1382) = -rxt(k,371)*y(k,206)
         mat(k,2054) = -rxt(k,372)*y(k,206)
         mat(k,1935) = -rxt(k,373)*y(k,206)
         mat(k,1829) = -rxt(k,375)*y(k,206)
         mat(k,1737) = -rxt(k,376)*y(k,206)

         mat(k,848) = .400_r8*rxt(k,387)*y(k,217)
         mat(k,1679) = .400_r8*rxt(k,387)*y(k,98)

         mat(k,66) = -(rxt(k,503)*y(k,203) + rxt(k,504)*y(k,124))
         mat(k,1860) = -rxt(k,503)*y(k,207)
         mat(k,1768) = -rxt(k,504)*y(k,207)

         mat(k,840) = rxt(k,506)*y(k,217)
         mat(k,1548) = rxt(k,506)*y(k,98)

         mat(k,72) = -(rxt(k,507)*y(k,203) + rxt(k,508)*y(k,124))
         mat(k,1861) = -rxt(k,507)*y(k,208)
         mat(k,1769) = -rxt(k,508)*y(k,208)

         mat(k,73) = rxt(k,509)*y(k,217)
         mat(k,1549) = rxt(k,509)*y(k,104)

         mat(k,1311) = -(rxt(k,333)*y(k,197) + rxt(k,334)*y(k,198) + rxt(k,335) &
                      *y(k,203) + rxt(k,336)*y(k,126) + (rxt(k,337) + rxt(k,338) &
                      ) * y(k,124))
         mat(k,1383) = -rxt(k,333)*y(k,209)
         mat(k,2055) = -rxt(k,334)*y(k,209)
         mat(k,1936) = -rxt(k,335)*y(k,209)
         mat(k,1738) = -rxt(k,336)*y(k,209)
         mat(k,1830) = -(rxt(k,337) + rxt(k,338)) * y(k,209)

         mat(k,1211) = .500_r8*rxt(k,340)*y(k,217)
         mat(k,311) = .200_r8*rxt(k,341)*y(k,217)
         mat(k,1331) = rxt(k,354)*y(k,217)
         mat(k,1680) = .500_r8*rxt(k,340)*y(k,105) + .200_r8*rxt(k,341)*y(k,106)  &
                      + rxt(k,354)*y(k,111)

         mat(k,721) = -(rxt(k,414)*y(k,203) + rxt(k,415)*y(k,124) + rxt(k,416) &
                      *y(k,125))
         mat(k,1903) = -rxt(k,414)*y(k,210)
         mat(k,1798) = -rxt(k,415)*y(k,210)
         mat(k,2172) = -rxt(k,416)*y(k,210)

         mat(k,1354) = -(rxt(k,342)*y(k,197) + rxt(k,343)*y(k,198) + rxt(k,344) &
                      *y(k,203) + 4._r8*rxt(k,345)*y(k,211) + rxt(k,346)*y(k,124) &
                      + rxt(k,347)*y(k,126) + rxt(k,355)*y(k,125))
         mat(k,1385) = -rxt(k,342)*y(k,211)
         mat(k,2057) = -rxt(k,343)*y(k,211)
         mat(k,1938) = -rxt(k,344)*y(k,211)
         mat(k,1832) = -rxt(k,346)*y(k,211)
         mat(k,1740) = -rxt(k,347)*y(k,211)
         mat(k,2184) = -rxt(k,355)*y(k,211)

         mat(k,1212) = .500_r8*rxt(k,340)*y(k,217)
         mat(k,312) = .500_r8*rxt(k,341)*y(k,217)
         mat(k,1682) = .500_r8*rxt(k,340)*y(k,105) + .500_r8*rxt(k,341)*y(k,106)

         mat(k,873) = -(rxt(k,417)*y(k,203) + rxt(k,418)*y(k,124) + rxt(k,419) &
                      *y(k,125))
         mat(k,1915) = -rxt(k,417)*y(k,212)
         mat(k,1807) = -rxt(k,418)*y(k,212)
         mat(k,2177) = -rxt(k,419)*y(k,212)

         mat(k,677) = -(rxt(k,348)*y(k,203) + rxt(k,349)*y(k,124))
         mat(k,1898) = -rxt(k,348)*y(k,213)
         mat(k,1796) = -rxt(k,349)*y(k,213)

         mat(k,507) = rxt(k,350)*y(k,217)
         mat(k,316) = rxt(k,351)*y(k,217)
         mat(k,1632) = rxt(k,350)*y(k,107) + rxt(k,351)*y(k,108)

         mat(k,80) = -(rxt(k,511)*y(k,203) + rxt(k,512)*y(k,124))
         mat(k,1862) = -rxt(k,511)*y(k,214)
         mat(k,1770) = -rxt(k,512)*y(k,214)

         mat(k,940) = rxt(k,514)*y(k,217)
         mat(k,1551) = rxt(k,514)*y(k,110)

         mat(k,1045) = -(rxt(k,446)*y(k,198) + rxt(k,447)*y(k,203) + rxt(k,448) &
                      *y(k,124) + rxt(k,449)*y(k,126))
         mat(k,2040) = -rxt(k,446)*y(k,215)
         mat(k,1922) = -rxt(k,447)*y(k,215)
         mat(k,1814) = -rxt(k,448)*y(k,215)
         mat(k,1721) = -rxt(k,449)*y(k,215)

         mat(k,1000) = rxt(k,440)*y(k,126)
         mat(k,951) = rxt(k,443)*y(k,126)
         mat(k,1721) = mat(k,1721) + rxt(k,440)*y(k,6) + rxt(k,443)*y(k,110)  &
                      + .500_r8*rxt(k,460)*y(k,177)
         mat(k,393) = rxt(k,450)*y(k,217)
         mat(k,970) = .500_r8*rxt(k,460)*y(k,126)
         mat(k,1663) = rxt(k,450)*y(k,128)

         mat(k,1526) = -(rxt(k,125)*y(k,77) + rxt(k,126)*y(k,229) + (rxt(k,129) &
                      + rxt(k,130)) * y(k,134) + (rxt(k,168) + rxt(k,169)) * y(k,113) &
                      + rxt(k,201)*y(k,33) + rxt(k,202)*y(k,34) + rxt(k,203)*y(k,36) &
                      + rxt(k,204)*y(k,37) + rxt(k,205)*y(k,38) + rxt(k,206)*y(k,39) &
                      + rxt(k,207)*y(k,40) + (rxt(k,208) + rxt(k,209)) * y(k,85) &
                      + rxt(k,228)*y(k,35) + rxt(k,229)*y(k,55) + rxt(k,230)*y(k,78) &
                      + (rxt(k,231) + rxt(k,232)) * y(k,81) + rxt(k,237)*y(k,64) &
                      + rxt(k,238)*y(k,65) + rxt(k,251)*y(k,41) + rxt(k,252)*y(k,43) &
                      + rxt(k,253)*y(k,82) + rxt(k,254)*y(k,83) + rxt(k,255)*y(k,84) &
                      + (rxt(k,272) + rxt(k,273) + rxt(k,274)) * y(k,54) + rxt(k,275) &
                      *y(k,86))
         mat(k,1405) = -rxt(k,125)*y(k,216)
         mat(k,2272) = -rxt(k,126)*y(k,216)
         mat(k,2124) = -(rxt(k,129) + rxt(k,130)) * y(k,216)
         mat(k,184) = -(rxt(k,168) + rxt(k,169)) * y(k,216)
         mat(k,100) = -rxt(k,201)*y(k,216)
         mat(k,144) = -rxt(k,202)*y(k,216)
         mat(k,115) = -rxt(k,203)*y(k,216)
         mat(k,154) = -rxt(k,204)*y(k,216)
         mat(k,119) = -rxt(k,205)*y(k,216)
         mat(k,159) = -rxt(k,206)*y(k,216)
         mat(k,123) = -rxt(k,207)*y(k,216)
         mat(k,2147) = -(rxt(k,208) + rxt(k,209)) * y(k,216)
         mat(k,150) = -rxt(k,228)*y(k,216)
         mat(k,449) = -rxt(k,229)*y(k,216)
         mat(k,108) = -rxt(k,230)*y(k,216)
         mat(k,807) = -(rxt(k,231) + rxt(k,232)) * y(k,216)
         mat(k,237) = -rxt(k,237)*y(k,216)
         mat(k,228) = -rxt(k,238)*y(k,216)
         mat(k,469) = -rxt(k,251)*y(k,216)
         mat(k,596) = -rxt(k,252)*y(k,216)
         mat(k,223) = -rxt(k,253)*y(k,216)
         mat(k,250) = -rxt(k,254)*y(k,216)
         mat(k,306) = -rxt(k,255)*y(k,216)
         mat(k,1434) = -(rxt(k,272) + rxt(k,273) + rxt(k,274)) * y(k,216)
         mat(k,180) = -rxt(k,275)*y(k,216)

         mat(k,1691) = -(rxt(k,143)*y(k,77) + rxt(k,144)*y(k,79) + rxt(k,145)*y(k,203) &
                      + rxt(k,146)*y(k,133) + rxt(k,147)*y(k,134) + (4._r8*rxt(k,148) &
                      + 4._r8*rxt(k,149)) * y(k,217) + rxt(k,151)*y(k,90) + rxt(k,163) &
                      *y(k,126) + rxt(k,164)*y(k,112) + rxt(k,172)*y(k,125) + rxt(k,173) &
                      *y(k,89) + rxt(k,192)*y(k,60) + (rxt(k,194) + rxt(k,195) &
                      ) * y(k,59) + rxt(k,197)*y(k,85) + rxt(k,200)*y(k,92) + rxt(k,224) &
                      *y(k,19) + rxt(k,226)*y(k,81) + rxt(k,240)*y(k,41) + rxt(k,242) &
                      *y(k,43) + rxt(k,243)*y(k,44) + rxt(k,245)*y(k,46) + rxt(k,247) &
                      *y(k,55) + rxt(k,248)*y(k,82) + rxt(k,249)*y(k,83) + rxt(k,250) &
                      *y(k,84) + rxt(k,259)*y(k,42) + rxt(k,264)*y(k,52) + rxt(k,265) &
                      *y(k,53) + rxt(k,266)*y(k,54) + rxt(k,267)*y(k,86) + rxt(k,268) &
                      *y(k,87) + rxt(k,276)*y(k,62) + rxt(k,278)*y(k,24) + rxt(k,285) &
                      *y(k,26) + rxt(k,286)*y(k,27) + rxt(k,288)*y(k,28) + rxt(k,290) &
                      *y(k,45) + rxt(k,291)*y(k,47) + rxt(k,296)*y(k,50) + rxt(k,297) &
                      *y(k,51) + rxt(k,302)*y(k,74) + rxt(k,303)*y(k,75) + rxt(k,304) &
                      *y(k,139) + rxt(k,305)*y(k,25) + rxt(k,313)*y(k,30) + rxt(k,314) &
                      *y(k,31) + rxt(k,316)*y(k,49) + rxt(k,317)*y(k,95) + rxt(k,318) &
                      *y(k,127) + rxt(k,321)*y(k,146) + rxt(k,325)*y(k,147) + rxt(k,326) &
                      *y(k,29) + rxt(k,327)*y(k,48) + rxt(k,329)*y(k,16) + rxt(k,332) &
                      *y(k,93) + rxt(k,340)*y(k,105) + rxt(k,341)*y(k,106) + rxt(k,350) &
                      *y(k,107) + rxt(k,351)*y(k,108) + rxt(k,352)*y(k,109) + rxt(k,354) &
                      *y(k,111) + rxt(k,357)*y(k,1) + rxt(k,361)*y(k,2) + rxt(k,362) &
                      *y(k,15) + rxt(k,363)*y(k,94) + rxt(k,364)*y(k,96) + rxt(k,365) &
                      *y(k,97) + rxt(k,377)*y(k,99) + rxt(k,378)*y(k,100) + rxt(k,385) &
                      *y(k,102) + rxt(k,387)*y(k,98) + rxt(k,388)*y(k,103) + rxt(k,389) &
                      *y(k,115) + rxt(k,390)*y(k,116) + rxt(k,396)*y(k,181) + rxt(k,399) &
                      *y(k,7) + rxt(k,402)*y(k,8) + rxt(k,403)*y(k,22) + rxt(k,405) &
                      *y(k,23) + rxt(k,409)*y(k,32) + rxt(k,410)*y(k,66) + rxt(k,422) &
                      *y(k,142) + rxt(k,425)*y(k,143) + rxt(k,429)*y(k,179) + rxt(k,430) &
                      *y(k,180) + rxt(k,432)*y(k,182) + rxt(k,435)*y(k,183) + rxt(k,438) &
                      *y(k,184) + rxt(k,439)*y(k,185) + rxt(k,442)*y(k,6) + rxt(k,445) &
                      *y(k,110) + rxt(k,450)*y(k,128) + rxt(k,454)*y(k,174) + rxt(k,455) &
                      *y(k,175) + rxt(k,459)*y(k,176) + rxt(k,461)*y(k,177) + rxt(k,462) &
                      *y(k,178) + (rxt(k,464) + rxt(k,478)) * y(k,67) + rxt(k,466) &
                      *y(k,137) + rxt(k,468)*y(k,151) + rxt(k,472)*y(k,148) + rxt(k,477) &
                      *y(k,150) + rxt(k,480)*y(k,120))
         mat(k,1406) = -rxt(k,143)*y(k,217)
         mat(k,604) = -rxt(k,144)*y(k,217)
         mat(k,1947) = -rxt(k,145)*y(k,217)
         mat(k,2247) = -rxt(k,146)*y(k,217)
         mat(k,2125) = -rxt(k,147)*y(k,217)
         mat(k,404) = -rxt(k,151)*y(k,217)
         mat(k,1748) = -rxt(k,163)*y(k,217)
         mat(k,494) = -rxt(k,164)*y(k,217)
         mat(k,2192) = -rxt(k,172)*y(k,217)
         mat(k,1451) = -rxt(k,173)*y(k,217)
         mat(k,886) = -rxt(k,192)*y(k,217)
         mat(k,1973) = -(rxt(k,194) + rxt(k,195)) * y(k,217)
         mat(k,2148) = -rxt(k,197)*y(k,217)
         mat(k,825) = -rxt(k,200)*y(k,217)
         mat(k,2216) = -rxt(k,224)*y(k,217)
         mat(k,808) = -rxt(k,226)*y(k,217)
         mat(k,470) = -rxt(k,240)*y(k,217)
         mat(k,597) = -rxt(k,242)*y(k,217)
         mat(k,126) = -rxt(k,243)*y(k,217)
         mat(k,367) = -rxt(k,245)*y(k,217)
         mat(k,450) = -rxt(k,247)*y(k,217)
         mat(k,224) = -rxt(k,248)*y(k,217)
         mat(k,251) = -rxt(k,249)*y(k,217)
         mat(k,307) = -rxt(k,250)*y(k,217)
         mat(k,1487) = -rxt(k,259)*y(k,217)
         mat(k,788) = -rxt(k,264)*y(k,217)
         mat(k,388) = -rxt(k,265)*y(k,217)
         mat(k,1435) = -rxt(k,266)*y(k,217)
         mat(k,181) = -rxt(k,267)*y(k,217)
         mat(k,931) = -rxt(k,268)*y(k,217)
         mat(k,1104) = -rxt(k,276)*y(k,217)
         mat(k,289) = -rxt(k,278)*y(k,217)
         mat(k,263) = -rxt(k,285)*y(k,217)
         mat(k,347) = -rxt(k,286)*y(k,217)
         mat(k,293) = -rxt(k,288)*y(k,217)
         mat(k,1089) = -rxt(k,290)*y(k,217)
         mat(k,103) = -rxt(k,291)*y(k,217)
         mat(k,686) = -rxt(k,296)*y(k,217)
         mat(k,614) = -rxt(k,297)*y(k,217)
         mat(k,1099) = -rxt(k,302)*y(k,217)
         mat(k,981) = -rxt(k,303)*y(k,217)
         mat(k,528) = -rxt(k,304)*y(k,217)
         mat(k,553) = -rxt(k,305)*y(k,217)
         mat(k,418) = -rxt(k,313)*y(k,217)
         mat(k,111) = -rxt(k,314)*y(k,217)
         mat(k,1224) = -rxt(k,316)*y(k,217)
         mat(k,1144) = -rxt(k,317)*y(k,217)
         mat(k,861) = -rxt(k,318)*y(k,217)
         mat(k,537) = -rxt(k,321)*y(k,217)
         mat(k,413) = -rxt(k,325)*y(k,217)
         mat(k,1032) = -rxt(k,326)*y(k,217)
         mat(k,925) = -rxt(k,327)*y(k,217)
         mat(k,354) = -rxt(k,329)*y(k,217)
         mat(k,1158) = -rxt(k,332)*y(k,217)
         mat(k,1215) = -rxt(k,340)*y(k,217)
         mat(k,313) = -rxt(k,341)*y(k,217)
         mat(k,510) = -rxt(k,350)*y(k,217)
         mat(k,319) = -rxt(k,351)*y(k,217)
         mat(k,577) = -rxt(k,352)*y(k,217)
         mat(k,1338) = -rxt(k,354)*y(k,217)
         mat(k,673) = -rxt(k,357)*y(k,217)
         mat(k,640) = -rxt(k,361)*y(k,217)
         mat(k,240) = -rxt(k,362)*y(k,217)
         mat(k,233) = -rxt(k,363)*y(k,217)
         mat(k,322) = -rxt(k,364)*y(k,217)
         mat(k,137) = -rxt(k,365)*y(k,217)
         mat(k,587) = -rxt(k,377)*y(k,217)
         mat(k,562) = -rxt(k,378)*y(k,217)
         mat(k,400) = -rxt(k,385)*y(k,217)
         mat(k,852) = -rxt(k,387)*y(k,217)
         mat(k,695) = -rxt(k,388)*y(k,217)
         mat(k,377) = -rxt(k,389)*y(k,217)
         mat(k,1079) = -rxt(k,390)*y(k,217)
         mat(k,205) = -rxt(k,396)*y(k,217)
         mat(k,168) = -rxt(k,399)*y(k,217)
         mat(k,425) = -rxt(k,402)*y(k,217)
         mat(k,255) = -rxt(k,403)*y(k,217)
         mat(k,342) = -rxt(k,405)*y(k,217)
         mat(k,268) = -rxt(k,409)*y(k,217)
         mat(k,197) = -rxt(k,410)*y(k,217)
         mat(k,177) = -rxt(k,422)*y(k,217)
         mat(k,336) = -rxt(k,425)*y(k,217)
         mat(k,663) = -rxt(k,429)*y(k,217)
         mat(k,192) = -rxt(k,430)*y(k,217)
         mat(k,214) = -rxt(k,432)*y(k,217)
         mat(k,710) = -rxt(k,435)*y(k,217)
         mat(k,219) = -rxt(k,438)*y(k,217)
         mat(k,431) = -rxt(k,439)*y(k,217)
         mat(k,1010) = -rxt(k,442)*y(k,217)
         mat(k,960) = -rxt(k,445)*y(k,217)
         mat(k,395) = -rxt(k,450)*y(k,217)
         mat(k,650) = -rxt(k,454)*y(k,217)
         mat(k,620) = -rxt(k,455)*y(k,217)
         mat(k,479) = -rxt(k,459)*y(k,217)
         mat(k,974) = -rxt(k,461)*y(k,217)
         mat(k,1064) = -rxt(k,462)*y(k,217)
         mat(k,300) = -(rxt(k,464) + rxt(k,478)) * y(k,217)
         mat(k,363) = -rxt(k,466)*y(k,217)
         mat(k,834) = -rxt(k,468)*y(k,217)
         mat(k,514) = -rxt(k,472)*y(k,217)
         mat(k,1235) = -rxt(k,477)*y(k,217)
         mat(k,97) = -rxt(k,480)*y(k,217)

         mat(k,1010) = mat(k,1010) + .630_r8*rxt(k,441)*y(k,134)
         mat(k,289) = mat(k,289) + .650_r8*rxt(k,278)*y(k,217)
         mat(k,553) = mat(k,553) + .130_r8*rxt(k,280)*y(k,134)
         mat(k,347) = mat(k,347) + .500_r8*rxt(k,286)*y(k,217)
         mat(k,1032) = mat(k,1032) + .360_r8*rxt(k,309)*y(k,134)
         mat(k,1487) = mat(k,1487) + rxt(k,258)*y(k,133)
         mat(k,388) = mat(k,388) + .300_r8*rxt(k,265)*y(k,217)
         mat(k,1435) = mat(k,1435) + rxt(k,272)*y(k,216)
         mat(k,2012) = rxt(k,181)*y(k,203)
         mat(k,869) = rxt(k,235)*y(k,229)
         mat(k,1466) = rxt(k,142)*y(k,134) + 2.000_r8*rxt(k,137)*y(k,203)
         mat(k,1406) = mat(k,1406) + rxt(k,134)*y(k,133) + rxt(k,125)*y(k,216)
         mat(k,604) = mat(k,604) + rxt(k,135)*y(k,133)
         mat(k,808) = mat(k,808) + rxt(k,225)*y(k,133) + rxt(k,231)*y(k,216)
         mat(k,2148) = mat(k,2148) + rxt(k,196)*y(k,133) + rxt(k,208)*y(k,216)
         mat(k,181) = mat(k,181) + rxt(k,275)*y(k,216)
         mat(k,780) = rxt(k,227)*y(k,133)
         mat(k,825) = mat(k,825) + rxt(k,199)*y(k,133)
         mat(k,852) = mat(k,852) + .320_r8*rxt(k,386)*y(k,134)
         mat(k,695) = mat(k,695) + .600_r8*rxt(k,388)*y(k,217)
         mat(k,1215) = mat(k,1215) + .240_r8*rxt(k,339)*y(k,134)
         mat(k,313) = mat(k,313) + .100_r8*rxt(k,341)*y(k,217)
         mat(k,960) = mat(k,960) + .630_r8*rxt(k,444)*y(k,134)
         mat(k,1338) = mat(k,1338) + .360_r8*rxt(k,353)*y(k,134)
         mat(k,1840) = rxt(k,165)*y(k,203)
         mat(k,1748) = mat(k,1748) + rxt(k,160)*y(k,203)
         mat(k,2247) = mat(k,2247) + rxt(k,258)*y(k,42) + rxt(k,134)*y(k,77)  &
                      + rxt(k,135)*y(k,79) + rxt(k,225)*y(k,81) + rxt(k,196)*y(k,85)  &
                      + rxt(k,227)*y(k,91) + rxt(k,199)*y(k,92) + rxt(k,140)*y(k,203)
         mat(k,2125) = mat(k,2125) + .630_r8*rxt(k,441)*y(k,6) + .130_r8*rxt(k,280) &
                      *y(k,25) + .360_r8*rxt(k,309)*y(k,29) + rxt(k,142)*y(k,76)  &
                      + .320_r8*rxt(k,386)*y(k,98) + .240_r8*rxt(k,339)*y(k,105)  &
                      + .630_r8*rxt(k,444)*y(k,110) + .360_r8*rxt(k,353)*y(k,111)  &
                      + rxt(k,141)*y(k,203)
         mat(k,537) = mat(k,537) + .500_r8*rxt(k,321)*y(k,217)
         mat(k,205) = mat(k,205) + .500_r8*rxt(k,396)*y(k,217)
         mat(k,520) = .400_r8*rxt(k,397)*y(k,203)
         mat(k,1390) = .450_r8*rxt(k,294)*y(k,203)
         mat(k,762) = .400_r8*rxt(k,411)*y(k,203)
         mat(k,1947) = mat(k,1947) + rxt(k,181)*y(k,56) + 2.000_r8*rxt(k,137)*y(k,76)  &
                      + rxt(k,165)*y(k,124) + rxt(k,160)*y(k,126) + rxt(k,140) &
                      *y(k,133) + rxt(k,141)*y(k,134) + .400_r8*rxt(k,397)*y(k,188)  &
                      + .450_r8*rxt(k,294)*y(k,197) + .400_r8*rxt(k,411)*y(k,199)  &
                      + .450_r8*rxt(k,344)*y(k,211) + .400_r8*rxt(k,417)*y(k,212)  &
                      + .200_r8*rxt(k,348)*y(k,213) + .150_r8*rxt(k,323)*y(k,220)
         mat(k,1358) = .450_r8*rxt(k,344)*y(k,203)
         mat(k,877) = .400_r8*rxt(k,417)*y(k,203)
         mat(k,680) = .200_r8*rxt(k,348)*y(k,203)
         mat(k,1527) = rxt(k,272)*y(k,54) + rxt(k,125)*y(k,77) + rxt(k,231)*y(k,81)  &
                      + rxt(k,208)*y(k,85) + rxt(k,275)*y(k,86) + 2.000_r8*rxt(k,126) &
                      *y(k,229)
         mat(k,1691) = mat(k,1691) + .650_r8*rxt(k,278)*y(k,24) + .500_r8*rxt(k,286) &
                      *y(k,27) + .300_r8*rxt(k,265)*y(k,53) + .600_r8*rxt(k,388) &
                      *y(k,103) + .100_r8*rxt(k,341)*y(k,106) + .500_r8*rxt(k,321) &
                      *y(k,146) + .500_r8*rxt(k,396)*y(k,181)
         mat(k,1134) = .150_r8*rxt(k,323)*y(k,203)
         mat(k,2273) = rxt(k,235)*y(k,73) + 2.000_r8*rxt(k,126)*y(k,216)

      end do

      end subroutine     nlnmat09

      subroutine     nlnmat10( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,454) = -(rxt(k,420)*y(k,203) + rxt(k,421)*y(k,124))
         mat(k,1884) = -rxt(k,420)*y(k,218)
         mat(k,1781) = -rxt(k,421)*y(k,218)

         mat(k,195) = .200_r8*rxt(k,410)*y(k,217)
         mat(k,175) = .140_r8*rxt(k,422)*y(k,217)
         mat(k,334) = rxt(k,425)*y(k,217)
         mat(k,1604) = .200_r8*rxt(k,410)*y(k,66) + .140_r8*rxt(k,422)*y(k,142)  &
                      + rxt(k,425)*y(k,143)

         mat(k,768) = -(rxt(k,319)*y(k,203) + rxt(k,320)*y(k,124))
         mat(k,1907) = -rxt(k,319)*y(k,219)
         mat(k,1802) = -rxt(k,320)*y(k,219)

         mat(k,1020) = rxt(k,326)*y(k,217)
         mat(k,533) = .500_r8*rxt(k,321)*y(k,217)
         mat(k,1641) = rxt(k,326)*y(k,29) + .500_r8*rxt(k,321)*y(k,146)

         mat(k,1129) = -(rxt(k,322)*y(k,198) + rxt(k,323)*y(k,203) + rxt(k,324) &
                      *y(k,124))
         mat(k,2047) = -rxt(k,322)*y(k,220)
         mat(k,1928) = -rxt(k,323)*y(k,220)
         mat(k,1821) = -rxt(k,324)*y(k,220)

         mat(k,1005) = .060_r8*rxt(k,441)*y(k,134)
         mat(k,923) = rxt(k,327)*y(k,217)
         mat(k,955) = .060_r8*rxt(k,444)*y(k,134)
         mat(k,2107) = .060_r8*rxt(k,441)*y(k,6) + .060_r8*rxt(k,444)*y(k,110)
         mat(k,410) = rxt(k,325)*y(k,217)
         mat(k,1061) = .150_r8*rxt(k,462)*y(k,217)
         mat(k,1670) = rxt(k,327)*y(k,48) + rxt(k,325)*y(k,147) + .150_r8*rxt(k,462) &
                      *y(k,178)

         mat(k,1115) = -(rxt(k,451)*y(k,198) + rxt(k,452)*y(k,203) + rxt(k,453) &
                      *y(k,124))
         mat(k,2046) = -rxt(k,451)*y(k,221)
         mat(k,1927) = -rxt(k,452)*y(k,221)
         mat(k,1820) = -rxt(k,453)*y(k,221)

         mat(k,1727) = .500_r8*rxt(k,460)*y(k,177)
         mat(k,648) = rxt(k,454)*y(k,217)
         mat(k,973) = .500_r8*rxt(k,460)*y(k,126) + rxt(k,461)*y(k,217)
         mat(k,1669) = rxt(k,454)*y(k,174) + rxt(k,461)*y(k,177)

         mat(k,912) = -(rxt(k,456)*y(k,198) + rxt(k,457)*y(k,203) + rxt(k,458) &
                      *y(k,124))
         mat(k,2036) = -rxt(k,456)*y(k,222)
         mat(k,1917) = -rxt(k,457)*y(k,222)
         mat(k,1809) = -rxt(k,458)*y(k,222)

         mat(k,994) = rxt(k,442)*y(k,217)
         mat(k,945) = rxt(k,445)*y(k,217)
         mat(k,475) = rxt(k,459)*y(k,217)
         mat(k,1655) = rxt(k,442)*y(k,6) + rxt(k,445)*y(k,110) + rxt(k,459)*y(k,176)

         mat(k,732) = -(rxt(k,427)*y(k,203) + rxt(k,428)*y(k,124))
         mat(k,1904) = -rxt(k,427)*y(k,223)
         mat(k,1799) = -rxt(k,428)*y(k,223)

         mat(k,657) = rxt(k,429)*y(k,217)
         mat(k,191) = .650_r8*rxt(k,430)*y(k,217)
         mat(k,1638) = rxt(k,429)*y(k,179) + .650_r8*rxt(k,430)*y(k,180)

         mat(k,86) = -(rxt(k,517)*y(k,203) + rxt(k,518)*y(k,124))
         mat(k,1863) = -rxt(k,517)*y(k,224)
         mat(k,1771) = -rxt(k,518)*y(k,224)

         mat(k,186) = rxt(k,516)*y(k,217)
         mat(k,1552) = rxt(k,516)*y(k,180)

         mat(k,1173) = -(rxt(k,391)*y(k,197) + rxt(k,392)*y(k,198) + rxt(k,393) &
                      *y(k,203) + rxt(k,394)*y(k,124) + rxt(k,395)*y(k,126))
         mat(k,1377) = -rxt(k,391)*y(k,225)
         mat(k,2049) = -rxt(k,392)*y(k,225)
         mat(k,1930) = -rxt(k,393)*y(k,225)
         mat(k,1824) = -rxt(k,394)*y(k,225)
         mat(k,1731) = -rxt(k,395)*y(k,225)

         mat(k,232) = rxt(k,363)*y(k,217)
         mat(k,321) = rxt(k,364)*y(k,217)
         mat(k,136) = rxt(k,365)*y(k,217)
         mat(k,691) = .400_r8*rxt(k,388)*y(k,217)
         mat(k,204) = .500_r8*rxt(k,396)*y(k,217)
         mat(k,1673) = rxt(k,363)*y(k,94) + rxt(k,364)*y(k,96) + rxt(k,365)*y(k,97)  &
                      + .400_r8*rxt(k,388)*y(k,103) + .500_r8*rxt(k,396)*y(k,181)

         mat(k,748) = -(rxt(k,433)*y(k,203) + rxt(k,434)*y(k,124))
         mat(k,1905) = -rxt(k,433)*y(k,226)
         mat(k,1800) = -rxt(k,434)*y(k,226)

         mat(k,211) = .560_r8*rxt(k,432)*y(k,217)
         mat(k,703) = rxt(k,435)*y(k,217)
         mat(k,1639) = .560_r8*rxt(k,432)*y(k,182) + rxt(k,435)*y(k,183)

         mat(k,92) = -(rxt(k,521)*y(k,203) + rxt(k,522)*y(k,124))
         mat(k,1864) = -rxt(k,521)*y(k,227)
         mat(k,1772) = -rxt(k,522)*y(k,227)

         mat(k,206) = rxt(k,520)*y(k,217)
         mat(k,1553) = rxt(k,520)*y(k,182)

         mat(k,499) = -(rxt(k,436)*y(k,203) + rxt(k,437)*y(k,124))
         mat(k,1889) = -rxt(k,436)*y(k,228)
         mat(k,1786) = -rxt(k,437)*y(k,228)

         mat(k,218) = .300_r8*rxt(k,438)*y(k,217)
         mat(k,428) = rxt(k,439)*y(k,217)
         mat(k,1611) = .300_r8*rxt(k,438)*y(k,184) + rxt(k,439)*y(k,185)

         mat(k,2285) = -(rxt(k,126)*y(k,216) + rxt(k,235)*y(k,73) + rxt(k,479) &
                      *y(k,152))
         mat(k,1539) = -rxt(k,126)*y(k,229)
         mat(k,872) = -rxt(k,235)*y(k,229)
         mat(k,260) = -rxt(k,479)*y(k,229)

         mat(k,296) = rxt(k,288)*y(k,217)
         mat(k,420) = rxt(k,313)*y(k,217)
         mat(k,112) = rxt(k,314)*y(k,217)
         mat(k,473) = rxt(k,240)*y(k,217)
         mat(k,1498) = rxt(k,259)*y(k,217)
         mat(k,602) = rxt(k,242)*y(k,217)
         mat(k,128) = rxt(k,243)*y(k,217)
         mat(k,1093) = rxt(k,290)*y(k,217)
         mat(k,372) = rxt(k,245)*y(k,217)
         mat(k,927) = rxt(k,327)*y(k,217)
         mat(k,1228) = rxt(k,316)*y(k,217)
         mat(k,688) = rxt(k,296)*y(k,217)
         mat(k,616) = rxt(k,297)*y(k,217)
         mat(k,390) = rxt(k,265)*y(k,217)
         mat(k,1442) = rxt(k,266)*y(k,217)
         mat(k,1475) = rxt(k,138)*y(k,203)
         mat(k,1412) = rxt(k,143)*y(k,217)
         mat(k,609) = rxt(k,144)*y(k,217)
         mat(k,811) = rxt(k,226)*y(k,217)
         mat(k,309) = rxt(k,250)*y(k,217)
         mat(k,2160) = (rxt(k,531)+rxt(k,536))*y(k,91) + (rxt(k,524)+rxt(k,530) &
                       +rxt(k,535))*y(k,92) + rxt(k,197)*y(k,217)
         mat(k,934) = rxt(k,268)*y(k,217)
         mat(k,1459) = rxt(k,173)*y(k,217)
         mat(k,408) = rxt(k,151)*y(k,217)
         mat(k,785) = (rxt(k,531)+rxt(k,536))*y(k,85)
         mat(k,830) = (rxt(k,524)+rxt(k,530)+rxt(k,535))*y(k,85) + rxt(k,200)*y(k,217)
         mat(k,1219) = .500_r8*rxt(k,340)*y(k,217)
         mat(k,98) = rxt(k,480)*y(k,217)
         mat(k,539) = rxt(k,321)*y(k,217)
         mat(k,414) = rxt(k,325)*y(k,217)
         mat(k,1959) = rxt(k,138)*y(k,76) + rxt(k,145)*y(k,217)
         mat(k,1703) = rxt(k,288)*y(k,28) + rxt(k,313)*y(k,30) + rxt(k,314)*y(k,31)  &
                      + rxt(k,240)*y(k,41) + rxt(k,259)*y(k,42) + rxt(k,242)*y(k,43)  &
                      + rxt(k,243)*y(k,44) + rxt(k,290)*y(k,45) + rxt(k,245)*y(k,46)  &
                      + rxt(k,327)*y(k,48) + rxt(k,316)*y(k,49) + rxt(k,296)*y(k,50)  &
                      + rxt(k,297)*y(k,51) + rxt(k,265)*y(k,53) + rxt(k,266)*y(k,54)  &
                      + rxt(k,143)*y(k,77) + rxt(k,144)*y(k,79) + rxt(k,226)*y(k,81)  &
                      + rxt(k,250)*y(k,84) + rxt(k,197)*y(k,85) + rxt(k,268)*y(k,87)  &
                      + rxt(k,173)*y(k,89) + rxt(k,151)*y(k,90) + rxt(k,200)*y(k,92)  &
                      + .500_r8*rxt(k,340)*y(k,105) + rxt(k,480)*y(k,120) + rxt(k,321) &
                      *y(k,146) + rxt(k,325)*y(k,147) + rxt(k,145)*y(k,203)  &
                      + 2.000_r8*rxt(k,148)*y(k,217)

      end do

      end subroutine     nlnmat10

      subroutine     nlnmat_finit( avec_len, mat, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  dti(veclen)
      real(r8), intent(in)    ::  lmat(veclen,nzcnt)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,   1) = lmat(k,   1)
         mat(k,   2) = lmat(k,   2)
         mat(k,   3) = lmat(k,   3)
         mat(k,   4) = lmat(k,   4)
         mat(k,   5) = lmat(k,   5)
         mat(k,   6) = lmat(k,   6)
         mat(k,   7) = lmat(k,   7)
         mat(k,   8) = lmat(k,   8)
         mat(k,   9) = lmat(k,   9)
         mat(k,  10) = lmat(k,  10)
         mat(k,  11) = lmat(k,  11)
         mat(k,  12) = lmat(k,  12)
         mat(k,  13) = lmat(k,  13)
         mat(k,  14) = lmat(k,  14)
         mat(k,  15) = lmat(k,  15)
         mat(k,  16) = lmat(k,  16)
         mat(k,  17) = lmat(k,  17)
         mat(k,  18) = lmat(k,  18)
         mat(k,  19) = lmat(k,  19)
         mat(k,  20) = lmat(k,  20)
         mat(k,  21) = lmat(k,  21)
         mat(k,  22) = lmat(k,  22)
         mat(k,  23) = lmat(k,  23)
         mat(k,  24) = lmat(k,  24)
         mat(k,  25) = lmat(k,  25)
         mat(k,  26) = lmat(k,  26)
         mat(k,  27) = lmat(k,  27)
         mat(k,  28) = lmat(k,  28)
         mat(k,  29) = lmat(k,  29)
         mat(k,  30) = lmat(k,  30)
         mat(k,  31) = lmat(k,  31)
         mat(k,  32) = lmat(k,  32)
         mat(k,  33) = lmat(k,  33)
         mat(k,  34) = lmat(k,  34)
         mat(k,  35) = lmat(k,  35)
         mat(k,  36) = lmat(k,  36)
         mat(k,  37) = lmat(k,  37)
         mat(k,  38) = lmat(k,  38)
         mat(k,  39) = lmat(k,  39)
         mat(k,  40) = lmat(k,  40)
         mat(k,  41) = lmat(k,  41)
         mat(k,  42) = lmat(k,  42)
         mat(k,  48) = mat(k,  48) + lmat(k,  48)
         mat(k,  54) = mat(k,  54) + lmat(k,  54)
         mat(k,  60) = mat(k,  60) + lmat(k,  60)
         mat(k,  66) = mat(k,  66) + lmat(k,  66)
         mat(k,  72) = mat(k,  72) + lmat(k,  72)
         mat(k,  74) = mat(k,  74) + lmat(k,  74)
         mat(k,  80) = mat(k,  80) + lmat(k,  80)
         mat(k,  86) = mat(k,  86) + lmat(k,  86)
         mat(k,  92) = mat(k,  92) + lmat(k,  92)
         mat(k,  93) = lmat(k,  93)
         mat(k,  94) = lmat(k,  94)
         mat(k,  95) = lmat(k,  95)
         mat(k,  96) = mat(k,  96) + lmat(k,  96)
         mat(k,  99) = mat(k,  99) + lmat(k,  99)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 106) = mat(k, 106) + lmat(k, 106)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 109) = mat(k, 109) + lmat(k, 109)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 117) = mat(k, 117) + lmat(k, 117)
         mat(k, 118) = mat(k, 118) + lmat(k, 118)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 121) = mat(k, 121) + lmat(k, 121)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 127) = mat(k, 127) + lmat(k, 127)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = lmat(k, 131)
         mat(k, 132) = lmat(k, 132)
         mat(k, 133) = lmat(k, 133)
         mat(k, 134) = lmat(k, 134)
         mat(k, 135) = mat(k, 135) + lmat(k, 135)
         mat(k, 138) = lmat(k, 138)
         mat(k, 139) = lmat(k, 139)
         mat(k, 140) = lmat(k, 140)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 145) = mat(k, 145) + lmat(k, 145)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 161) = lmat(k, 161)
         mat(k, 162) = lmat(k, 162)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 170) = lmat(k, 170)
         mat(k, 171) = lmat(k, 171)
         mat(k, 172) = lmat(k, 172)
         mat(k, 173) = lmat(k, 173)
         mat(k, 174) = mat(k, 174) + lmat(k, 174)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 187) = mat(k, 187) + lmat(k, 187)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 199) = lmat(k, 199)
         mat(k, 200) = lmat(k, 200)
         mat(k, 201) = lmat(k, 201)
         mat(k, 202) = lmat(k, 202)
         mat(k, 203) = mat(k, 203) + lmat(k, 203)
         mat(k, 205) = mat(k, 205) + lmat(k, 205)
         mat(k, 208) = mat(k, 208) + lmat(k, 208)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 226) = mat(k, 226) + lmat(k, 226)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 229) = mat(k, 229) + lmat(k, 229)
         mat(k, 230) = mat(k, 230) + lmat(k, 230)
         mat(k, 231) = lmat(k, 231)
         mat(k, 233) = mat(k, 233) + lmat(k, 233)
         mat(k, 234) = lmat(k, 234)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 236) = mat(k, 236) + lmat(k, 236)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 241) = lmat(k, 241)
         mat(k, 242) = lmat(k, 242)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = lmat(k, 246)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 258) = lmat(k, 258)
         mat(k, 259) = lmat(k, 259)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 266) = lmat(k, 266)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 269) = lmat(k, 269)
         mat(k, 270) = lmat(k, 270)
         mat(k, 271) = lmat(k, 271)
         mat(k, 272) = lmat(k, 272)
         mat(k, 273) = lmat(k, 273)
         mat(k, 274) = lmat(k, 274)
         mat(k, 275) = lmat(k, 275)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 278) = lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = lmat(k, 280)
         mat(k, 281) = lmat(k, 281)
         mat(k, 282) = lmat(k, 282)
         mat(k, 283) = lmat(k, 283)
         mat(k, 284) = lmat(k, 284)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 297) = mat(k, 297) + lmat(k, 297)
         mat(k, 303) = mat(k, 303) + lmat(k, 303)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 315) = mat(k, 315) + lmat(k, 315)
         mat(k, 317) = lmat(k, 317)
         mat(k, 318) = lmat(k, 318)
         mat(k, 319) = mat(k, 319) + lmat(k, 319)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = lmat(k, 325)
         mat(k, 326) = lmat(k, 326)
         mat(k, 327) = lmat(k, 327)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 331) = mat(k, 331) + lmat(k, 331)
         mat(k, 332) = lmat(k, 332)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 335) = lmat(k, 335)
         mat(k, 336) = mat(k, 336) + lmat(k, 336)
         mat(k, 337) = lmat(k, 337)
         mat(k, 338) = lmat(k, 338)
         mat(k, 339) = mat(k, 339) + lmat(k, 339)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 343) = lmat(k, 343)
         mat(k, 344) = mat(k, 344) + lmat(k, 344)
         mat(k, 346) = mat(k, 346) + lmat(k, 346)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 348) = lmat(k, 348)
         mat(k, 349) = mat(k, 349) + lmat(k, 349)
         mat(k, 357) = mat(k, 357) + lmat(k, 357)
         mat(k, 358) = lmat(k, 358)
         mat(k, 360) = mat(k, 360) + lmat(k, 360)
         mat(k, 365) = mat(k, 365) + lmat(k, 365)
         mat(k, 369) = mat(k, 369) + lmat(k, 369)
         mat(k, 370) = lmat(k, 370)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 376) = lmat(k, 376)
         mat(k, 379) = lmat(k, 379)
         mat(k, 380) = lmat(k, 380)
         mat(k, 381) = lmat(k, 381)
         mat(k, 382) = lmat(k, 382)
         mat(k, 383) = lmat(k, 383)
         mat(k, 384) = lmat(k, 384)
         mat(k, 385) = mat(k, 385) + lmat(k, 385)
         mat(k, 386) = lmat(k, 386)
         mat(k, 387) = mat(k, 387) + lmat(k, 387)
         mat(k, 388) = mat(k, 388) + lmat(k, 388)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 392) = lmat(k, 392)
         mat(k, 394) = lmat(k, 394)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 396) = lmat(k, 396)
         mat(k, 397) = mat(k, 397) + lmat(k, 397)
         mat(k, 398) = lmat(k, 398)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 402) = lmat(k, 402)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 405) = lmat(k, 405)
         mat(k, 406) = lmat(k, 406)
         mat(k, 407) = mat(k, 407) + lmat(k, 407)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 411) = lmat(k, 411)
         mat(k, 412) = lmat(k, 412)
         mat(k, 413) = mat(k, 413) + lmat(k, 413)
         mat(k, 415) = mat(k, 415) + lmat(k, 415)
         mat(k, 417) = lmat(k, 417)
         mat(k, 418) = mat(k, 418) + lmat(k, 418)
         mat(k, 419) = lmat(k, 419)
         mat(k, 421) = mat(k, 421) + lmat(k, 421)
         mat(k, 422) = lmat(k, 422)
         mat(k, 424) = lmat(k, 424)
         mat(k, 425) = mat(k, 425) + lmat(k, 425)
         mat(k, 426) = lmat(k, 426)
         mat(k, 427) = mat(k, 427) + lmat(k, 427)
         mat(k, 429) = lmat(k, 429)
         mat(k, 430) = lmat(k, 430)
         mat(k, 431) = mat(k, 431) + lmat(k, 431)
         mat(k, 432) = lmat(k, 432)
         mat(k, 435) = mat(k, 435) + lmat(k, 435)
         mat(k, 441) = mat(k, 441) + lmat(k, 441)
         mat(k, 443) = lmat(k, 443)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 460) = lmat(k, 460)
         mat(k, 461) = lmat(k, 461)
         mat(k, 462) = lmat(k, 462)
         mat(k, 463) = mat(k, 463) + lmat(k, 463)
         mat(k, 467) = mat(k, 467) + lmat(k, 467)
         mat(k, 468) = mat(k, 468) + lmat(k, 468)
         mat(k, 474) = mat(k, 474) + lmat(k, 474)
         mat(k, 476) = lmat(k, 476)
         mat(k, 477) = lmat(k, 477)
         mat(k, 478) = lmat(k, 478)
         mat(k, 479) = mat(k, 479) + lmat(k, 479)
         mat(k, 480) = lmat(k, 480)
         mat(k, 483) = mat(k, 483) + lmat(k, 483)
         mat(k, 491) = mat(k, 491) + lmat(k, 491)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 497) = mat(k, 497) + lmat(k, 497)
         mat(k, 499) = mat(k, 499) + lmat(k, 499)
         mat(k, 506) = mat(k, 506) + lmat(k, 506)
         mat(k, 508) = lmat(k, 508)
         mat(k, 509) = lmat(k, 509)
         mat(k, 511) = mat(k, 511) + lmat(k, 511)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 516) = lmat(k, 516)
         mat(k, 518) = mat(k, 518) + lmat(k, 518)
         mat(k, 524) = mat(k, 524) + lmat(k, 524)
         mat(k, 525) = lmat(k, 525)
         mat(k, 526) = lmat(k, 526)
         mat(k, 529) = mat(k, 529) + lmat(k, 529)
         mat(k, 530) = lmat(k, 530)
         mat(k, 531) = lmat(k, 531)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 534) = lmat(k, 534)
         mat(k, 536) = lmat(k, 536)
         mat(k, 537) = mat(k, 537) + lmat(k, 537)
         mat(k, 538) = lmat(k, 538)
         mat(k, 540) = mat(k, 540) + lmat(k, 540)
         mat(k, 541) = lmat(k, 541)
         mat(k, 542) = lmat(k, 542)
         mat(k, 543) = lmat(k, 543)
         mat(k, 544) = mat(k, 544) + lmat(k, 544)
         mat(k, 545) = lmat(k, 545)
         mat(k, 546) = mat(k, 546) + lmat(k, 546)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 561) = lmat(k, 561)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 572) = mat(k, 572) + lmat(k, 572)
         mat(k, 575) = lmat(k, 575)
         mat(k, 580) = lmat(k, 580)
         mat(k, 581) = mat(k, 581) + lmat(k, 581)
         mat(k, 585) = lmat(k, 585)
         mat(k, 590) = lmat(k, 590)
         mat(k, 591) = lmat(k, 591)
         mat(k, 592) = lmat(k, 592)
         mat(k, 593) = lmat(k, 593)
         mat(k, 594) = mat(k, 594) + lmat(k, 594)
         mat(k, 595) = mat(k, 595) + lmat(k, 595)
         mat(k, 600) = lmat(k, 600)
         mat(k, 603) = mat(k, 603) + lmat(k, 603)
         mat(k, 604) = mat(k, 604) + lmat(k, 604)
         mat(k, 610) = mat(k, 610) + lmat(k, 610)
         mat(k, 611) = mat(k, 611) + lmat(k, 611)
         mat(k, 614) = mat(k, 614) + lmat(k, 614)
         mat(k, 615) = lmat(k, 615)
         mat(k, 617) = mat(k, 617) + lmat(k, 617)
         mat(k, 618) = mat(k, 618) + lmat(k, 618)
         mat(k, 619) = lmat(k, 619)
         mat(k, 621) = lmat(k, 621)
         mat(k, 622) = mat(k, 622) + lmat(k, 622)
         mat(k, 625) = mat(k, 625) + lmat(k, 625)
         mat(k, 631) = lmat(k, 631)
         mat(k, 632) = mat(k, 632) + lmat(k, 632)
         mat(k, 636) = lmat(k, 636)
         mat(k, 637) = lmat(k, 637)
         mat(k, 639) = lmat(k, 639)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 641) = lmat(k, 641)
         mat(k, 642) = mat(k, 642) + lmat(k, 642)
         mat(k, 643) = lmat(k, 643)
         mat(k, 644) = lmat(k, 644)
         mat(k, 645) = lmat(k, 645)
         mat(k, 646) = lmat(k, 646)
         mat(k, 647) = lmat(k, 647)
         mat(k, 649) = lmat(k, 649)
         mat(k, 650) = mat(k, 650) + lmat(k, 650)
         mat(k, 651) = lmat(k, 651)
         mat(k, 652) = lmat(k, 652)
         mat(k, 653) = lmat(k, 653)
         mat(k, 654) = lmat(k, 654)
         mat(k, 655) = mat(k, 655) + lmat(k, 655)
         mat(k, 660) = lmat(k, 660)
         mat(k, 662) = lmat(k, 662)
         mat(k, 663) = mat(k, 663) + lmat(k, 663)
         mat(k, 664) = lmat(k, 664)
         mat(k, 665) = lmat(k, 665)
         mat(k, 666) = mat(k, 666) + lmat(k, 666)
         mat(k, 669) = mat(k, 669) + lmat(k, 669)
         mat(k, 670) = mat(k, 670) + lmat(k, 670)
         mat(k, 672) = mat(k, 672) + lmat(k, 672)
         mat(k, 674) = lmat(k, 674)
         mat(k, 675) = mat(k, 675) + lmat(k, 675)
         mat(k, 677) = mat(k, 677) + lmat(k, 677)
         mat(k, 684) = mat(k, 684) + lmat(k, 684)
         mat(k, 690) = mat(k, 690) + lmat(k, 690)
         mat(k, 692) = lmat(k, 692)
         mat(k, 693) = lmat(k, 693)
         mat(k, 694) = lmat(k, 694)
         mat(k, 695) = mat(k, 695) + lmat(k, 695)
         mat(k, 696) = lmat(k, 696)
         mat(k, 697) = lmat(k, 697)
         mat(k, 698) = lmat(k, 698)
         mat(k, 699) = lmat(k, 699)
         mat(k, 700) = lmat(k, 700)
         mat(k, 701) = mat(k, 701) + lmat(k, 701)
         mat(k, 706) = lmat(k, 706)
         mat(k, 708) = lmat(k, 708)
         mat(k, 710) = mat(k, 710) + lmat(k, 710)
         mat(k, 711) = lmat(k, 711)
         mat(k, 714) = mat(k, 714) + lmat(k, 714)
         mat(k, 721) = mat(k, 721) + lmat(k, 721)
         mat(k, 732) = mat(k, 732) + lmat(k, 732)
         mat(k, 748) = mat(k, 748) + lmat(k, 748)
         mat(k, 759) = mat(k, 759) + lmat(k, 759)
         mat(k, 768) = mat(k, 768) + lmat(k, 768)
         mat(k, 778) = mat(k, 778) + lmat(k, 778)
         mat(k, 779) = lmat(k, 779)
         mat(k, 780) = mat(k, 780) + lmat(k, 780)
         mat(k, 786) = mat(k, 786) + lmat(k, 786)
         mat(k, 791) = mat(k, 791) + lmat(k, 791)
         mat(k, 801) = lmat(k, 801)
         mat(k, 802) = lmat(k, 802)
         mat(k, 803) = lmat(k, 803)
         mat(k, 804) = mat(k, 804) + lmat(k, 804)
         mat(k, 805) = mat(k, 805) + lmat(k, 805)
         mat(k, 806) = mat(k, 806) + lmat(k, 806)
         mat(k, 815) = mat(k, 815) + lmat(k, 815)
         mat(k, 824) = mat(k, 824) + lmat(k, 824)
         mat(k, 825) = mat(k, 825) + lmat(k, 825)
         mat(k, 827) = mat(k, 827) + lmat(k, 827)
         mat(k, 832) = mat(k, 832) + lmat(k, 832)
         mat(k, 833) = lmat(k, 833)
         mat(k, 836) = lmat(k, 836)
         mat(k, 841) = mat(k, 841) + lmat(k, 841)
         mat(k, 857) = mat(k, 857) + lmat(k, 857)
         mat(k, 859) = lmat(k, 859)
         mat(k, 860) = lmat(k, 860)
         mat(k, 862) = mat(k, 862) + lmat(k, 862)
         mat(k, 864) = mat(k, 864) + lmat(k, 864)
         mat(k, 873) = mat(k, 873) + lmat(k, 873)
         mat(k, 883) = mat(k, 883) + lmat(k, 883)
         mat(k, 884) = mat(k, 884) + lmat(k, 884)
         mat(k, 885) = mat(k, 885) + lmat(k, 885)
         mat(k, 887) = mat(k, 887) + lmat(k, 887)
         mat(k, 888) = mat(k, 888) + lmat(k, 888)
         mat(k, 889) = mat(k, 889) + lmat(k, 889)
         mat(k, 891) = lmat(k, 891)
         mat(k, 899) = mat(k, 899) + lmat(k, 899)
         mat(k, 912) = mat(k, 912) + lmat(k, 912)
         mat(k, 922) = mat(k, 922) + lmat(k, 922)
         mat(k, 924) = lmat(k, 924)
         mat(k, 926) = lmat(k, 926)
         mat(k, 929) = mat(k, 929) + lmat(k, 929)
         mat(k, 948) = mat(k, 948) + lmat(k, 948)
         mat(k, 969) = mat(k, 969) + lmat(k, 969)
         mat(k, 971) = lmat(k, 971)
         mat(k, 972) = lmat(k, 972)
         mat(k, 976) = lmat(k, 976)
         mat(k, 977) = lmat(k, 977)
         mat(k, 979) = mat(k, 979) + lmat(k, 979)
         mat(k, 980) = mat(k, 980) + lmat(k, 980)
         mat(k, 982) = mat(k, 982) + lmat(k, 982)
         mat(k, 999) = mat(k, 999) + lmat(k, 999)
         mat(k,1024) = mat(k,1024) + lmat(k,1024)
         mat(k,1045) = mat(k,1045) + lmat(k,1045)
         mat(k,1057) = mat(k,1057) + lmat(k,1057)
         mat(k,1058) = mat(k,1058) + lmat(k,1058)
         mat(k,1059) = mat(k,1059) + lmat(k,1059)
         mat(k,1060) = mat(k,1060) + lmat(k,1060)
         mat(k,1061) = mat(k,1061) + lmat(k,1061)
         mat(k,1062) = mat(k,1062) + lmat(k,1062)
         mat(k,1063) = mat(k,1063) + lmat(k,1063)
         mat(k,1065) = mat(k,1065) + lmat(k,1065)
         mat(k,1069) = lmat(k,1069)
         mat(k,1073) = mat(k,1073) + lmat(k,1073)
         mat(k,1077) = lmat(k,1077)
         mat(k,1081) = mat(k,1081) + lmat(k,1081)
         mat(k,1083) = lmat(k,1083)
         mat(k,1085) = mat(k,1085) + lmat(k,1085)
         mat(k,1086) = lmat(k,1086)
         mat(k,1091) = lmat(k,1091)
         mat(k,1092) = lmat(k,1092)
         mat(k,1096) = mat(k,1096) + lmat(k,1096)
         mat(k,1097) = lmat(k,1097)
         mat(k,1098) = mat(k,1098) + lmat(k,1098)
         mat(k,1100) = mat(k,1100) + lmat(k,1100)
         mat(k,1103) = mat(k,1103) + lmat(k,1103)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1129) = mat(k,1129) + lmat(k,1129)
         mat(k,1140) = mat(k,1140) + lmat(k,1140)
         mat(k,1142) = lmat(k,1142)
         mat(k,1143) = lmat(k,1143)
         mat(k,1145) = mat(k,1145) + lmat(k,1145)
         mat(k,1147) = lmat(k,1147)
         mat(k,1148) = lmat(k,1148)
         mat(k,1149) = lmat(k,1149)
         mat(k,1150) = lmat(k,1150)
         mat(k,1152) = lmat(k,1152)
         mat(k,1153) = mat(k,1153) + lmat(k,1153)
         mat(k,1155) = lmat(k,1155)
         mat(k,1156) = lmat(k,1156)
         mat(k,1157) = lmat(k,1157)
         mat(k,1161) = mat(k,1161) + lmat(k,1161)
         mat(k,1163) = lmat(k,1163)
         mat(k,1173) = mat(k,1173) + lmat(k,1173)
         mat(k,1193) = mat(k,1193) + lmat(k,1193)
         mat(k,1208) = mat(k,1208) + lmat(k,1208)
         mat(k,1209) = mat(k,1209) + lmat(k,1209)
         mat(k,1212) = mat(k,1212) + lmat(k,1212)
         mat(k,1213) = mat(k,1213) + lmat(k,1213)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1216) = mat(k,1216) + lmat(k,1216)
         mat(k,1220) = mat(k,1220) + lmat(k,1220)
         mat(k,1221) = mat(k,1221) + lmat(k,1221)
         mat(k,1222) = mat(k,1222) + lmat(k,1222)
         mat(k,1226) = lmat(k,1226)
         mat(k,1230) = lmat(k,1230)
         mat(k,1231) = mat(k,1231) + lmat(k,1231)
         mat(k,1232) = mat(k,1232) + lmat(k,1232)
         mat(k,1243) = lmat(k,1243)
         mat(k,1256) = mat(k,1256) + lmat(k,1256)
         mat(k,1272) = lmat(k,1272)
         mat(k,1288) = mat(k,1288) + lmat(k,1288)
         mat(k,1298) = mat(k,1298) + lmat(k,1298)
         mat(k,1311) = mat(k,1311) + lmat(k,1311)
         mat(k,1326) = lmat(k,1326)
         mat(k,1328) = mat(k,1328) + lmat(k,1328)
         mat(k,1332) = mat(k,1332) + lmat(k,1332)
         mat(k,1334) = mat(k,1334) + lmat(k,1334)
         mat(k,1342) = lmat(k,1342)
         mat(k,1354) = mat(k,1354) + lmat(k,1354)
         mat(k,1386) = mat(k,1386) + lmat(k,1386)
         mat(k,1401) = mat(k,1401) + lmat(k,1401)
         mat(k,1415) = mat(k,1415) + lmat(k,1415)
         mat(k,1426) = lmat(k,1426)
         mat(k,1428) = lmat(k,1428)
         mat(k,1429) = mat(k,1429) + lmat(k,1429)
         mat(k,1430) = mat(k,1430) + lmat(k,1430)
         mat(k,1432) = mat(k,1432) + lmat(k,1432)
         mat(k,1433) = mat(k,1433) + lmat(k,1433)
         mat(k,1435) = mat(k,1435) + lmat(k,1435)
         mat(k,1439) = mat(k,1439) + lmat(k,1439)
         mat(k,1441) = lmat(k,1441)
         mat(k,1442) = mat(k,1442) + lmat(k,1442)
         mat(k,1447) = mat(k,1447) + lmat(k,1447)
         mat(k,1451) = mat(k,1451) + lmat(k,1451)
         mat(k,1457) = lmat(k,1457)
         mat(k,1463) = mat(k,1463) + lmat(k,1463)
         mat(k,1468) = mat(k,1468) + lmat(k,1468)
         mat(k,1479) = mat(k,1479) + lmat(k,1479)
         mat(k,1480) = lmat(k,1480)
         mat(k,1484) = mat(k,1484) + lmat(k,1484)
         mat(k,1485) = mat(k,1485) + lmat(k,1485)
         mat(k,1526) = mat(k,1526) + lmat(k,1526)
         mat(k,1538) = mat(k,1538) + lmat(k,1538)
         mat(k,1691) = mat(k,1691) + lmat(k,1691)
         mat(k,1744) = mat(k,1744) + lmat(k,1744)
         mat(k,1749) = mat(k,1749) + lmat(k,1749)
         mat(k,1750) = mat(k,1750) + lmat(k,1750)
         mat(k,1757) = mat(k,1757) + lmat(k,1757)
         mat(k,1759) = mat(k,1759) + lmat(k,1759)
         mat(k,1785) = mat(k,1785) + lmat(k,1785)
         mat(k,1842) = mat(k,1842) + lmat(k,1842)
         mat(k,1851) = mat(k,1851) + lmat(k,1851)
         mat(k,1950) = mat(k,1950) + lmat(k,1950)
         mat(k,1959) = mat(k,1959) + lmat(k,1959)
         mat(k,1977) = mat(k,1977) + lmat(k,1977)
         mat(k,1978) = mat(k,1978) + lmat(k,1978)
         mat(k,1984) = mat(k,1984) + lmat(k,1984)
         mat(k,2017) = mat(k,2017) + lmat(k,2017)
         mat(k,2070) = mat(k,2070) + lmat(k,2070)
         mat(k,2124) = mat(k,2124) + lmat(k,2124)
         mat(k,2132) = mat(k,2132) + lmat(k,2132)
         mat(k,2136) = mat(k,2136) + lmat(k,2136)
         mat(k,2145) = mat(k,2145) + lmat(k,2145)
         mat(k,2153) = mat(k,2153) + lmat(k,2153)
         mat(k,2156) = mat(k,2156) + lmat(k,2156)
         mat(k,2188) = mat(k,2188) + lmat(k,2188)
         mat(k,2192) = mat(k,2192) + lmat(k,2192)
         mat(k,2194) = mat(k,2194) + lmat(k,2194)
         mat(k,2201) = mat(k,2201) + lmat(k,2201)
         mat(k,2203) = mat(k,2203) + lmat(k,2203)
         mat(k,2211) = mat(k,2211) + lmat(k,2211)
         mat(k,2226) = mat(k,2226) + lmat(k,2226)
         mat(k,2227) = mat(k,2227) + lmat(k,2227)
         mat(k,2254) = mat(k,2254) + lmat(k,2254)
         mat(k,2258) = mat(k,2258) + lmat(k,2258)
         mat(k,2266) = lmat(k,2266)
         mat(k,2270) = lmat(k,2270)
         mat(k,2272) = mat(k,2272) + lmat(k,2272)
         mat(k,2273) = mat(k,2273) + lmat(k,2273)
         mat(k,2284) = lmat(k,2284)
         mat(k,2285) = mat(k,2285) + lmat(k,2285)
         mat(k, 212) = 0._r8
         mat(k, 213) = 0._r8
         mat(k, 249) = 0._r8
         mat(k, 305) = 0._r8
         mat(k, 341) = 0._r8
         mat(k, 436) = 0._r8
         mat(k, 437) = 0._r8
         mat(k, 456) = 0._r8
         mat(k, 484) = 0._r8
         mat(k, 486) = 0._r8
         mat(k, 493) = 0._r8
         mat(k, 502) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 633) = 0._r8
         mat(k, 634) = 0._r8
         mat(k, 638) = 0._r8
         mat(k, 656) = 0._r8
         mat(k, 658) = 0._r8
         mat(k, 659) = 0._r8
         mat(k, 661) = 0._r8
         mat(k, 667) = 0._r8
         mat(k, 668) = 0._r8
         mat(k, 671) = 0._r8
         mat(k, 702) = 0._r8
         mat(k, 704) = 0._r8
         mat(k, 705) = 0._r8
         mat(k, 707) = 0._r8
         mat(k, 709) = 0._r8
         mat(k, 715) = 0._r8
         mat(k, 717) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 733) = 0._r8
         mat(k, 734) = 0._r8
         mat(k, 736) = 0._r8
         mat(k, 738) = 0._r8
         mat(k, 747) = 0._r8
         mat(k, 749) = 0._r8
         mat(k, 750) = 0._r8
         mat(k, 752) = 0._r8
         mat(k, 754) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 776) = 0._r8
         mat(k, 781) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 800) = 0._r8
         mat(k, 818) = 0._r8
         mat(k, 837) = 0._r8
         mat(k, 893) = 0._r8
         mat(k, 897) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 902) = 0._r8
         mat(k, 904) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 915) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 930) = 0._r8
         mat(k, 933) = 0._r8
         mat(k, 944) = 0._r8
         mat(k, 950) = 0._r8
         mat(k, 953) = 0._r8
         mat(k, 956) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 962) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 966) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 968) = 0._r8
         mat(k, 983) = 0._r8
         mat(k, 993) = 0._r8
         mat(k, 998) = 0._r8
         mat(k,1002) = 0._r8
         mat(k,1004) = 0._r8
         mat(k,1006) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1012) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1016) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1028) = 0._r8
         mat(k,1029) = 0._r8
         mat(k,1034) = 0._r8
         mat(k,1038) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1046) = 0._r8
         mat(k,1047) = 0._r8
         mat(k,1048) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1066) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1068) = 0._r8
         mat(k,1071) = 0._r8
         mat(k,1074) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1076) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1080) = 0._r8
         mat(k,1082) = 0._r8
         mat(k,1084) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1116) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1118) = 0._r8
         mat(k,1120) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1126) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1162) = 0._r8
         mat(k,1164) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1192) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1206) = 0._r8
         mat(k,1218) = 0._r8
         mat(k,1227) = 0._r8
         mat(k,1237) = 0._r8
         mat(k,1244) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1250) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1255) = 0._r8
         mat(k,1257) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1263) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1270) = 0._r8
         mat(k,1271) = 0._r8
         mat(k,1275) = 0._r8
         mat(k,1278) = 0._r8
         mat(k,1279) = 0._r8
         mat(k,1282) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1285) = 0._r8
         mat(k,1286) = 0._r8
         mat(k,1289) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1312) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1316) = 0._r8
         mat(k,1322) = 0._r8
         mat(k,1323) = 0._r8
         mat(k,1329) = 0._r8
         mat(k,1333) = 0._r8
         mat(k,1335) = 0._r8
         mat(k,1336) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1345) = 0._r8
         mat(k,1346) = 0._r8
         mat(k,1351) = 0._r8
         mat(k,1352) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1387) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1391) = 0._r8
         mat(k,1397) = 0._r8
         mat(k,1398) = 0._r8
         mat(k,1402) = 0._r8
         mat(k,1403) = 0._r8
         mat(k,1407) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1416) = 0._r8
         mat(k,1418) = 0._r8
         mat(k,1419) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1431) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1446) = 0._r8
         mat(k,1448) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1454) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1456) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1461) = 0._r8
         mat(k,1462) = 0._r8
         mat(k,1464) = 0._r8
         mat(k,1465) = 0._r8
         mat(k,1467) = 0._r8
         mat(k,1469) = 0._r8
         mat(k,1470) = 0._r8
         mat(k,1472) = 0._r8
         mat(k,1473) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1486) = 0._r8
         mat(k,1489) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1493) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1523) = 0._r8
         mat(k,1528) = 0._r8
         mat(k,1536) = 0._r8
         mat(k,1605) = 0._r8
         mat(k,1622) = 0._r8
         mat(k,1637) = 0._r8
         mat(k,1640) = 0._r8
         mat(k,1651) = 0._r8
         mat(k,1652) = 0._r8
         mat(k,1674) = 0._r8
         mat(k,1690) = 0._r8
         mat(k,1710) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1722) = 0._r8
         mat(k,1728) = 0._r8
         mat(k,1730) = 0._r8
         mat(k,1735) = 0._r8
         mat(k,1742) = 0._r8
         mat(k,1743) = 0._r8
         mat(k,1745) = 0._r8
         mat(k,1747) = 0._r8
         mat(k,1752) = 0._r8
         mat(k,1753) = 0._r8
         mat(k,1754) = 0._r8
         mat(k,1755) = 0._r8
         mat(k,1756) = 0._r8
         mat(k,1758) = 0._r8
         mat(k,1760) = 0._r8
         mat(k,1806) = 0._r8
         mat(k,1835) = 0._r8
         mat(k,1836) = 0._r8
         mat(k,1837) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1848) = 0._r8
         mat(k,1852) = 0._r8
         mat(k,1868) = 0._r8
         mat(k,1885) = 0._r8
         mat(k,1886) = 0._r8
         mat(k,1914) = 0._r8
         mat(k,1918) = 0._r8
         mat(k,1920) = 0._r8
         mat(k,1923) = 0._r8
         mat(k,1925) = 0._r8
         mat(k,1929) = 0._r8
         mat(k,1932) = 0._r8
         mat(k,1937) = 0._r8
         mat(k,1942) = 0._r8
         mat(k,1943) = 0._r8
         mat(k,1946) = 0._r8
         mat(k,1969) = 0._r8
         mat(k,1970) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1980) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,1999) = 0._r8
         mat(k,2000) = 0._r8
         mat(k,2001) = 0._r8
         mat(k,2002) = 0._r8
         mat(k,2004) = 0._r8
         mat(k,2008) = 0._r8
         mat(k,2011) = 0._r8
         mat(k,2014) = 0._r8
         mat(k,2021) = 0._r8
         mat(k,2022) = 0._r8
         mat(k,2023) = 0._r8
         mat(k,2024) = 0._r8
         mat(k,2035) = 0._r8
         mat(k,2059) = 0._r8
         mat(k,2060) = 0._r8
         mat(k,2061) = 0._r8
         mat(k,2063) = 0._r8
         mat(k,2064) = 0._r8
         mat(k,2065) = 0._r8
         mat(k,2071) = 0._r8
         mat(k,2072) = 0._r8
         mat(k,2074) = 0._r8
         mat(k,2075) = 0._r8
         mat(k,2076) = 0._r8
         mat(k,2089) = 0._r8
         mat(k,2093) = 0._r8
         mat(k,2098) = 0._r8
         mat(k,2101) = 0._r8
         mat(k,2104) = 0._r8
         mat(k,2106) = 0._r8
         mat(k,2108) = 0._r8
         mat(k,2109) = 0._r8
         mat(k,2113) = 0._r8
         mat(k,2114) = 0._r8
         mat(k,2115) = 0._r8
         mat(k,2117) = 0._r8
         mat(k,2121) = 0._r8
         mat(k,2133) = 0._r8
         mat(k,2137) = 0._r8
         mat(k,2143) = 0._r8
         mat(k,2146) = 0._r8
         mat(k,2149) = 0._r8
         mat(k,2150) = 0._r8
         mat(k,2151) = 0._r8
         mat(k,2154) = 0._r8
         mat(k,2155) = 0._r8
         mat(k,2157) = 0._r8
         mat(k,2158) = 0._r8
         mat(k,2171) = 0._r8
         mat(k,2174) = 0._r8
         mat(k,2175) = 0._r8
         mat(k,2179) = 0._r8
         mat(k,2180) = 0._r8
         mat(k,2181) = 0._r8
         mat(k,2182) = 0._r8
         mat(k,2186) = 0._r8
         mat(k,2187) = 0._r8
         mat(k,2189) = 0._r8
         mat(k,2190) = 0._r8
         mat(k,2191) = 0._r8
         mat(k,2197) = 0._r8
         mat(k,2198) = 0._r8
         mat(k,2200) = 0._r8
         mat(k,2204) = 0._r8
         mat(k,2212) = 0._r8
         mat(k,2213) = 0._r8
         mat(k,2214) = 0._r8
         mat(k,2215) = 0._r8
         mat(k,2217) = 0._r8
         mat(k,2222) = 0._r8
         mat(k,2223) = 0._r8
         mat(k,2224) = 0._r8
         mat(k,2228) = 0._r8
         mat(k,2230) = 0._r8
         mat(k,2236) = 0._r8
         mat(k,2242) = 0._r8
         mat(k,2243) = 0._r8
         mat(k,2246) = 0._r8
         mat(k,2253) = 0._r8
         mat(k,2259) = 0._r8
         mat(k,2263) = 0._r8
         mat(k,2265) = 0._r8
         mat(k,2267) = 0._r8
         mat(k,2268) = 0._r8
         mat(k,2269) = 0._r8
         mat(k,2271) = 0._r8
         mat(k,2274) = 0._r8
         mat(k,2275) = 0._r8
         mat(k,2276) = 0._r8
         mat(k,2277) = 0._r8
         mat(k,2278) = 0._r8
         mat(k,2279) = 0._r8
         mat(k,2280) = 0._r8
         mat(k,2281) = 0._r8
         mat(k,2282) = 0._r8
         mat(k,2283) = 0._r8
         mat(k,   1) = mat(k,   1) - dti(k)
         mat(k,   2) = mat(k,   2) - dti(k)
         mat(k,   3) = mat(k,   3) - dti(k)
         mat(k,   4) = mat(k,   4) - dti(k)
         mat(k,   5) = mat(k,   5) - dti(k)
         mat(k,   6) = mat(k,   6) - dti(k)
         mat(k,   7) = mat(k,   7) - dti(k)
         mat(k,   8) = mat(k,   8) - dti(k)
         mat(k,   9) = mat(k,   9) - dti(k)
         mat(k,  10) = mat(k,  10) - dti(k)
         mat(k,  11) = mat(k,  11) - dti(k)
         mat(k,  12) = mat(k,  12) - dti(k)
         mat(k,  13) = mat(k,  13) - dti(k)
         mat(k,  14) = mat(k,  14) - dti(k)
         mat(k,  15) = mat(k,  15) - dti(k)
         mat(k,  16) = mat(k,  16) - dti(k)
         mat(k,  17) = mat(k,  17) - dti(k)
         mat(k,  18) = mat(k,  18) - dti(k)
         mat(k,  19) = mat(k,  19) - dti(k)
         mat(k,  20) = mat(k,  20) - dti(k)
         mat(k,  21) = mat(k,  21) - dti(k)
         mat(k,  22) = mat(k,  22) - dti(k)
         mat(k,  23) = mat(k,  23) - dti(k)
         mat(k,  24) = mat(k,  24) - dti(k)
         mat(k,  25) = mat(k,  25) - dti(k)
         mat(k,  26) = mat(k,  26) - dti(k)
         mat(k,  27) = mat(k,  27) - dti(k)
         mat(k,  28) = mat(k,  28) - dti(k)
         mat(k,  29) = mat(k,  29) - dti(k)
         mat(k,  30) = mat(k,  30) - dti(k)
         mat(k,  31) = mat(k,  31) - dti(k)
         mat(k,  32) = mat(k,  32) - dti(k)
         mat(k,  33) = mat(k,  33) - dti(k)
         mat(k,  34) = mat(k,  34) - dti(k)
         mat(k,  35) = mat(k,  35) - dti(k)
         mat(k,  36) = mat(k,  36) - dti(k)
         mat(k,  37) = mat(k,  37) - dti(k)
         mat(k,  38) = mat(k,  38) - dti(k)
         mat(k,  39) = mat(k,  39) - dti(k)
         mat(k,  40) = mat(k,  40) - dti(k)
         mat(k,  41) = mat(k,  41) - dti(k)
         mat(k,  42) = mat(k,  42) - dti(k)
         mat(k,  48) = mat(k,  48) - dti(k)
         mat(k,  54) = mat(k,  54) - dti(k)
         mat(k,  60) = mat(k,  60) - dti(k)
         mat(k,  66) = mat(k,  66) - dti(k)
         mat(k,  72) = mat(k,  72) - dti(k)
         mat(k,  74) = mat(k,  74) - dti(k)
         mat(k,  80) = mat(k,  80) - dti(k)
         mat(k,  86) = mat(k,  86) - dti(k)
         mat(k,  92) = mat(k,  92) - dti(k)
         mat(k,  93) = mat(k,  93) - dti(k)
         mat(k,  96) = mat(k,  96) - dti(k)
         mat(k,  99) = mat(k,  99) - dti(k)
         mat(k, 102) = mat(k, 102) - dti(k)
         mat(k, 105) = mat(k, 105) - dti(k)
         mat(k, 109) = mat(k, 109) - dti(k)
         mat(k, 113) = mat(k, 113) - dti(k)
         mat(k, 117) = mat(k, 117) - dti(k)
         mat(k, 121) = mat(k, 121) - dti(k)
         mat(k, 125) = mat(k, 125) - dti(k)
         mat(k, 129) = mat(k, 129) - dti(k)
         mat(k, 132) = mat(k, 132) - dti(k)
         mat(k, 135) = mat(k, 135) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 141) = mat(k, 141) - dti(k)
         mat(k, 146) = mat(k, 146) - dti(k)
         mat(k, 151) = mat(k, 151) - dti(k)
         mat(k, 156) = mat(k, 156) - dti(k)
         mat(k, 161) = mat(k, 161) - dti(k)
         mat(k, 164) = mat(k, 164) - dti(k)
         mat(k, 170) = mat(k, 170) - dti(k)
         mat(k, 174) = mat(k, 174) - dti(k)
         mat(k, 179) = mat(k, 179) - dti(k)
         mat(k, 183) = mat(k, 183) - dti(k)
         mat(k, 187) = mat(k, 187) - dti(k)
         mat(k, 194) = mat(k, 194) - dti(k)
         mat(k, 199) = mat(k, 199) - dti(k)
         mat(k, 203) = mat(k, 203) - dti(k)
         mat(k, 208) = mat(k, 208) - dti(k)
         mat(k, 216) = mat(k, 216) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 226) = mat(k, 226) - dti(k)
         mat(k, 230) = mat(k, 230) - dti(k)
         mat(k, 235) = mat(k, 235) - dti(k)
         mat(k, 238) = mat(k, 238) - dti(k)
         mat(k, 241) = mat(k, 241) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 248) = mat(k, 248) - dti(k)
         mat(k, 253) = mat(k, 253) - dti(k)
         mat(k, 257) = mat(k, 257) - dti(k)
         mat(k, 261) = mat(k, 261) - dti(k)
         mat(k, 265) = mat(k, 265) - dti(k)
         mat(k, 269) = mat(k, 269) - dti(k)
         mat(k, 273) = mat(k, 273) - dti(k)
         mat(k, 279) = mat(k, 279) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 291) = mat(k, 291) - dti(k)
         mat(k, 297) = mat(k, 297) - dti(k)
         mat(k, 304) = mat(k, 304) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 315) = mat(k, 315) - dti(k)
         mat(k, 320) = mat(k, 320) - dti(k)
         mat(k, 323) = mat(k, 323) - dti(k)
         mat(k, 328) = mat(k, 328) - dti(k)
         mat(k, 333) = mat(k, 333) - dti(k)
         mat(k, 339) = mat(k, 339) - dti(k)
         mat(k, 344) = mat(k, 344) - dti(k)
         mat(k, 349) = mat(k, 349) - dti(k)
         mat(k, 357) = mat(k, 357) - dti(k)
         mat(k, 365) = mat(k, 365) - dti(k)
         mat(k, 373) = mat(k, 373) - dti(k)
         mat(k, 379) = mat(k, 379) - dti(k)
         mat(k, 385) = mat(k, 385) - dti(k)
         mat(k, 391) = mat(k, 391) - dti(k)
         mat(k, 397) = mat(k, 397) - dti(k)
         mat(k, 403) = mat(k, 403) - dti(k)
         mat(k, 409) = mat(k, 409) - dti(k)
         mat(k, 415) = mat(k, 415) - dti(k)
         mat(k, 421) = mat(k, 421) - dti(k)
         mat(k, 427) = mat(k, 427) - dti(k)
         mat(k, 435) = mat(k, 435) - dti(k)
         mat(k, 441) = mat(k, 441) - dti(k)
         mat(k, 447) = mat(k, 447) - dti(k)
         mat(k, 454) = mat(k, 454) - dti(k)
         mat(k, 460) = mat(k, 460) - dti(k)
         mat(k, 463) = mat(k, 463) - dti(k)
         mat(k, 467) = mat(k, 467) - dti(k)
         mat(k, 474) = mat(k, 474) - dti(k)
         mat(k, 483) = mat(k, 483) - dti(k)
         mat(k, 491) = mat(k, 491) - dti(k)
         mat(k, 499) = mat(k, 499) - dti(k)
         mat(k, 506) = mat(k, 506) - dti(k)
         mat(k, 511) = mat(k, 511) - dti(k)
         mat(k, 518) = mat(k, 518) - dti(k)
         mat(k, 524) = mat(k, 524) - dti(k)
         mat(k, 532) = mat(k, 532) - dti(k)
         mat(k, 540) = mat(k, 540) - dti(k)
         mat(k, 548) = mat(k, 548) - dti(k)
         mat(k, 556) = mat(k, 556) - dti(k)
         mat(k, 564) = mat(k, 564) - dti(k)
         mat(k, 572) = mat(k, 572) - dti(k)
         mat(k, 581) = mat(k, 581) - dti(k)
         mat(k, 590) = mat(k, 590) - dti(k)
         mat(k, 594) = mat(k, 594) - dti(k)
         mat(k, 603) = mat(k, 603) - dti(k)
         mat(k, 610) = mat(k, 610) - dti(k)
         mat(k, 617) = mat(k, 617) - dti(k)
         mat(k, 625) = mat(k, 625) - dti(k)
         mat(k, 632) = mat(k, 632) - dti(k)
         mat(k, 642) = mat(k, 642) - dti(k)
         mat(k, 655) = mat(k, 655) - dti(k)
         mat(k, 666) = mat(k, 666) - dti(k)
         mat(k, 677) = mat(k, 677) - dti(k)
         mat(k, 684) = mat(k, 684) - dti(k)
         mat(k, 690) = mat(k, 690) - dti(k)
         mat(k, 701) = mat(k, 701) - dti(k)
         mat(k, 714) = mat(k, 714) - dti(k)
         mat(k, 721) = mat(k, 721) - dti(k)
         mat(k, 732) = mat(k, 732) - dti(k)
         mat(k, 748) = mat(k, 748) - dti(k)
         mat(k, 759) = mat(k, 759) - dti(k)
         mat(k, 768) = mat(k, 768) - dti(k)
         mat(k, 778) = mat(k, 778) - dti(k)
         mat(k, 786) = mat(k, 786) - dti(k)
         mat(k, 791) = mat(k, 791) - dti(k)
         mat(k, 801) = mat(k, 801) - dti(k)
         mat(k, 804) = mat(k, 804) - dti(k)
         mat(k, 815) = mat(k, 815) - dti(k)
         mat(k, 824) = mat(k, 824) - dti(k)
         mat(k, 832) = mat(k, 832) - dti(k)
         mat(k, 841) = mat(k, 841) - dti(k)
         mat(k, 857) = mat(k, 857) - dti(k)
         mat(k, 864) = mat(k, 864) - dti(k)
         mat(k, 873) = mat(k, 873) - dti(k)
         mat(k, 884) = mat(k, 884) - dti(k)
         mat(k, 899) = mat(k, 899) - dti(k)
         mat(k, 912) = mat(k, 912) - dti(k)
         mat(k, 922) = mat(k, 922) - dti(k)
         mat(k, 929) = mat(k, 929) - dti(k)
         mat(k, 948) = mat(k, 948) - dti(k)
         mat(k, 969) = mat(k, 969) - dti(k)
         mat(k, 979) = mat(k, 979) - dti(k)
         mat(k, 999) = mat(k, 999) - dti(k)
         mat(k,1024) = mat(k,1024) - dti(k)
         mat(k,1045) = mat(k,1045) - dti(k)
         mat(k,1059) = mat(k,1059) - dti(k)
         mat(k,1073) = mat(k,1073) - dti(k)
         mat(k,1085) = mat(k,1085) - dti(k)
         mat(k,1096) = mat(k,1096) - dti(k)
         mat(k,1103) = mat(k,1103) - dti(k)
         mat(k,1115) = mat(k,1115) - dti(k)
         mat(k,1129) = mat(k,1129) - dti(k)
         mat(k,1140) = mat(k,1140) - dti(k)
         mat(k,1153) = mat(k,1153) - dti(k)
         mat(k,1173) = mat(k,1173) - dti(k)
         mat(k,1193) = mat(k,1193) - dti(k)
         mat(k,1209) = mat(k,1209) - dti(k)
         mat(k,1221) = mat(k,1221) - dti(k)
         mat(k,1232) = mat(k,1232) - dti(k)
         mat(k,1256) = mat(k,1256) - dti(k)
         mat(k,1288) = mat(k,1288) - dti(k)
         mat(k,1311) = mat(k,1311) - dti(k)
         mat(k,1332) = mat(k,1332) - dti(k)
         mat(k,1354) = mat(k,1354) - dti(k)
         mat(k,1386) = mat(k,1386) - dti(k)
         mat(k,1401) = mat(k,1401) - dti(k)
         mat(k,1415) = mat(k,1415) - dti(k)
         mat(k,1430) = mat(k,1430) - dti(k)
         mat(k,1447) = mat(k,1447) - dti(k)
         mat(k,1463) = mat(k,1463) - dti(k)
         mat(k,1485) = mat(k,1485) - dti(k)
         mat(k,1526) = mat(k,1526) - dti(k)
         mat(k,1691) = mat(k,1691) - dti(k)
         mat(k,1749) = mat(k,1749) - dti(k)
         mat(k,1842) = mat(k,1842) - dti(k)
         mat(k,1950) = mat(k,1950) - dti(k)
         mat(k,1977) = mat(k,1977) - dti(k)
         mat(k,2017) = mat(k,2017) - dti(k)
         mat(k,2070) = mat(k,2070) - dti(k)
         mat(k,2132) = mat(k,2132) - dti(k)
         mat(k,2156) = mat(k,2156) - dti(k)
         mat(k,2201) = mat(k,2201) - dti(k)
         mat(k,2226) = mat(k,2226) - dti(k)
         mat(k,2258) = mat(k,2258) - dti(k)
         mat(k,2285) = mat(k,2285) - dti(k)
      end do

      end subroutine nlnmat_finit

      subroutine     nlnmat( avec_len, mat, y, rxt, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  dti(veclen)
      real(r8), intent(in)    ::  lmat(veclen,nzcnt)
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)

      call     nlnmat01( avec_len, mat, y, rxt )
      call     nlnmat02( avec_len, mat, y, rxt )
      call     nlnmat03( avec_len, mat, y, rxt )
      call     nlnmat04( avec_len, mat, y, rxt )
      call     nlnmat05( avec_len, mat, y, rxt )
      call     nlnmat06( avec_len, mat, y, rxt )
      call     nlnmat07( avec_len, mat, y, rxt )
      call     nlnmat08( avec_len, mat, y, rxt )
      call     nlnmat09( avec_len, mat, y, rxt )
      call     nlnmat10( avec_len, mat, y, rxt )
      call     nlnmat_finit( avec_len, mat, lmat, dti )

      end subroutine nlnmat

      end module mo_nln_matrix

      module mo_lu_factor

      use chem_mods, only: veclen
      private
      public :: lu_fac

      contains
                                                                        
      subroutine lu_fac01( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,1) = 1._r8 / lu(k,1)
                                                                        
         lu(k,2) = 1._r8 / lu(k,2)
                                                                        
         lu(k,3) = 1._r8 / lu(k,3)
                                                                        
         lu(k,4) = 1._r8 / lu(k,4)
                                                                        
         lu(k,5) = 1._r8 / lu(k,5)
                                                                        
         lu(k,6) = 1._r8 / lu(k,6)
                                                                        
         lu(k,7) = 1._r8 / lu(k,7)
                                                                        
         lu(k,8) = 1._r8 / lu(k,8)
                                                                        
         lu(k,9) = 1._r8 / lu(k,9)
                                                                        
         lu(k,10) = 1._r8 / lu(k,10)
                                                                        
         lu(k,11) = 1._r8 / lu(k,11)
                                                                        
         lu(k,12) = 1._r8 / lu(k,12)
                                                                        
         lu(k,13) = 1._r8 / lu(k,13)
                                                                        
         lu(k,14) = 1._r8 / lu(k,14)
                                                                        
         lu(k,15) = 1._r8 / lu(k,15)
                                                                        
         lu(k,16) = 1._r8 / lu(k,16)
                                                                        
         lu(k,17) = 1._r8 / lu(k,17)
                                                                        
         lu(k,18) = 1._r8 / lu(k,18)
                                                                        
         lu(k,19) = 1._r8 / lu(k,19)
                                                                        
         lu(k,20) = 1._r8 / lu(k,20)
                                                                        
         lu(k,21) = 1._r8 / lu(k,21)
                                                                        
         lu(k,22) = 1._r8 / lu(k,22)
                                                                        
         lu(k,23) = 1._r8 / lu(k,23)
                                                                        
         lu(k,24) = 1._r8 / lu(k,24)
                                                                        
         lu(k,25) = 1._r8 / lu(k,25)
                                                                        
         lu(k,26) = 1._r8 / lu(k,26)
                                                                        
         lu(k,27) = 1._r8 / lu(k,27)
                                                                        
         lu(k,28) = 1._r8 / lu(k,28)
                                                                        
         lu(k,29) = 1._r8 / lu(k,29)
                                                                        
         lu(k,30) = 1._r8 / lu(k,30)
                                                                        
         lu(k,31) = 1._r8 / lu(k,31)
                                                                        
         lu(k,32) = 1._r8 / lu(k,32)
                                                                        
         lu(k,33) = 1._r8 / lu(k,33)
                                                                        
         lu(k,34) = 1._r8 / lu(k,34)
                                                                        
         lu(k,35) = 1._r8 / lu(k,35)
                                                                        
         lu(k,36) = 1._r8 / lu(k,36)
                                                                        
         lu(k,37) = 1._r8 / lu(k,37)
                                                                        
         lu(k,38) = 1._r8 / lu(k,38)
                                                                        
         lu(k,39) = 1._r8 / lu(k,39)
                                                                        
         lu(k,40) = 1._r8 / lu(k,40)
                                                                        
         lu(k,41) = 1._r8 / lu(k,41)
                                                                        
         lu(k,42) = 1._r8 / lu(k,42)
                                                                        
         lu(k,48) = 1._r8 / lu(k,48)
                                                                        
         lu(k,54) = 1._r8 / lu(k,54)
                                                                        
         lu(k,60) = 1._r8 / lu(k,60)
                                                                        
         lu(k,66) = 1._r8 / lu(k,66)
                                                                        
         lu(k,72) = 1._r8 / lu(k,72)
                                                                        
         lu(k,74) = 1._r8 / lu(k,74)
                                                                        
         lu(k,80) = 1._r8 / lu(k,80)
                                                                        
         lu(k,86) = 1._r8 / lu(k,86)
                                                                        
         lu(k,92) = 1._r8 / lu(k,92)
                                                                        
      end do
                                                                        
      end subroutine lu_fac01
                                                                        
      subroutine lu_fac02( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,93) = 1._r8 / lu(k,93)
         lu(k,94) = lu(k,94) * lu(k,93)
         lu(k,95) = lu(k,95) * lu(k,93)
         lu(k,1977) = lu(k,1977) - lu(k,94) * lu(k,1960)
         lu(k,1978) = lu(k,1978) - lu(k,95) * lu(k,1960)
                                                                        
         lu(k,96) = 1._r8 / lu(k,96)
         lu(k,97) = lu(k,97) * lu(k,96)
         lu(k,98) = lu(k,98) * lu(k,96)
         lu(k,1691) = lu(k,1691) - lu(k,97) * lu(k,1554)
         lu(k,1703) = lu(k,1703) - lu(k,98) * lu(k,1554)
                                                                        
         lu(k,99) = 1._r8 / lu(k,99)
         lu(k,100) = lu(k,100) * lu(k,99)
         lu(k,101) = lu(k,101) * lu(k,99)
         lu(k,1526) = lu(k,1526) - lu(k,100) * lu(k,1499)
         lu(k,1532) = lu(k,1532) - lu(k,101) * lu(k,1499)
                                                                        
         lu(k,102) = 1._r8 / lu(k,102)
         lu(k,103) = lu(k,103) * lu(k,102)
         lu(k,104) = lu(k,104) * lu(k,102)
         lu(k,1691) = lu(k,1691) - lu(k,103) * lu(k,1555)
         lu(k,1694) = lu(k,1694) - lu(k,104) * lu(k,1555)
                                                                        
         lu(k,105) = 1._r8 / lu(k,105)
         lu(k,106) = lu(k,106) * lu(k,105)
         lu(k,107) = lu(k,107) * lu(k,105)
         lu(k,108) = lu(k,108) * lu(k,105)
         lu(k,1512) = lu(k,1512) - lu(k,106) * lu(k,1500)
         lu(k,1521) = lu(k,1521) - lu(k,107) * lu(k,1500)
         lu(k,1526) = lu(k,1526) - lu(k,108) * lu(k,1500)
                                                                        
         lu(k,109) = 1._r8 / lu(k,109)
         lu(k,110) = lu(k,110) * lu(k,109)
         lu(k,111) = lu(k,111) * lu(k,109)
         lu(k,112) = lu(k,112) * lu(k,109)
         lu(k,1643) = lu(k,1643) - lu(k,110) * lu(k,1556)
         lu(k,1691) = lu(k,1691) - lu(k,111) * lu(k,1556)
         lu(k,1703) = lu(k,1703) - lu(k,112) * lu(k,1556)
                                                                        
         lu(k,113) = 1._r8 / lu(k,113)
         lu(k,114) = lu(k,114) * lu(k,113)
         lu(k,115) = lu(k,115) * lu(k,113)
         lu(k,116) = lu(k,116) * lu(k,113)
         lu(k,1511) = lu(k,1511) - lu(k,114) * lu(k,1501)
         lu(k,1526) = lu(k,1526) - lu(k,115) * lu(k,1501)
         lu(k,1532) = lu(k,1532) - lu(k,116) * lu(k,1501)
                                                                        
         lu(k,117) = 1._r8 / lu(k,117)
         lu(k,118) = lu(k,118) * lu(k,117)
         lu(k,119) = lu(k,119) * lu(k,117)
         lu(k,120) = lu(k,120) * lu(k,117)
         lu(k,1512) = lu(k,1512) - lu(k,118) * lu(k,1502)
         lu(k,1526) = lu(k,1526) - lu(k,119) * lu(k,1502)
         lu(k,1532) = lu(k,1532) - lu(k,120) * lu(k,1502)
                                                                        
         lu(k,121) = 1._r8 / lu(k,121)
         lu(k,122) = lu(k,122) * lu(k,121)
         lu(k,123) = lu(k,123) * lu(k,121)
         lu(k,124) = lu(k,124) * lu(k,121)
         lu(k,1512) = lu(k,1512) - lu(k,122) * lu(k,1503)
         lu(k,1526) = lu(k,1526) - lu(k,123) * lu(k,1503)
         lu(k,1532) = lu(k,1532) - lu(k,124) * lu(k,1503)
                                                                        
         lu(k,125) = 1._r8 / lu(k,125)
         lu(k,126) = lu(k,126) * lu(k,125)
         lu(k,127) = lu(k,127) * lu(k,125)
         lu(k,128) = lu(k,128) * lu(k,125)
         lu(k,1691) = lu(k,1691) - lu(k,126) * lu(k,1557)
         lu(k,1696) = lu(k,1696) - lu(k,127) * lu(k,1557)
         lu(k,1703) = lu(k,1703) - lu(k,128) * lu(k,1557)
                                                                        
         lu(k,129) = 1._r8 / lu(k,129)
         lu(k,130) = lu(k,130) * lu(k,129)
         lu(k,131) = lu(k,131) * lu(k,129)
         lu(k,713) = lu(k,713) - lu(k,130) * lu(k,712)
         lu(k,717) = - lu(k,131) * lu(k,712)
         lu(k,1868) = - lu(k,130) * lu(k,1865)
         lu(k,1947) = lu(k,1947) - lu(k,131) * lu(k,1865)
                                                                        
         lu(k,132) = 1._r8 / lu(k,132)
         lu(k,133) = lu(k,133) * lu(k,132)
         lu(k,134) = lu(k,134) * lu(k,132)
         lu(k,257) = lu(k,257) - lu(k,133) * lu(k,256)
         lu(k,260) = lu(k,260) - lu(k,134) * lu(k,256)
         lu(k,2261) = lu(k,2261) - lu(k,133) * lu(k,2260)
         lu(k,2285) = lu(k,2285) - lu(k,134) * lu(k,2260)
                                                                        
         lu(k,135) = 1._r8 / lu(k,135)
         lu(k,136) = lu(k,136) * lu(k,135)
         lu(k,137) = lu(k,137) * lu(k,135)
         lu(k,691) = lu(k,691) - lu(k,136) * lu(k,689)
         lu(k,695) = lu(k,695) - lu(k,137) * lu(k,689)
         lu(k,1673) = lu(k,1673) - lu(k,136) * lu(k,1558)
         lu(k,1691) = lu(k,1691) - lu(k,137) * lu(k,1558)
                                                                        
         lu(k,138) = 1._r8 / lu(k,138)
         lu(k,139) = lu(k,139) * lu(k,138)
         lu(k,140) = lu(k,140) * lu(k,138)
         lu(k,518) = lu(k,518) - lu(k,139) * lu(k,517)
         lu(k,523) = lu(k,523) - lu(k,140) * lu(k,517)
         lu(k,2167) = lu(k,2167) - lu(k,139) * lu(k,2161)
         lu(k,2201) = lu(k,2201) - lu(k,140) * lu(k,2161)
                                                                        
         lu(k,141) = 1._r8 / lu(k,141)
         lu(k,142) = lu(k,142) * lu(k,141)
         lu(k,143) = lu(k,143) * lu(k,141)
         lu(k,144) = lu(k,144) * lu(k,141)
         lu(k,145) = lu(k,145) * lu(k,141)
         lu(k,1512) = lu(k,1512) - lu(k,142) * lu(k,1504)
         lu(k,1521) = lu(k,1521) - lu(k,143) * lu(k,1504)
         lu(k,1526) = lu(k,1526) - lu(k,144) * lu(k,1504)
         lu(k,1532) = lu(k,1532) - lu(k,145) * lu(k,1504)
                                                                        
      end do
                                                                        
      end subroutine lu_fac02
                                                                        
      subroutine lu_fac03( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,146) = 1._r8 / lu(k,146)
         lu(k,147) = lu(k,147) * lu(k,146)
         lu(k,148) = lu(k,148) * lu(k,146)
         lu(k,149) = lu(k,149) * lu(k,146)
         lu(k,150) = lu(k,150) * lu(k,146)
         lu(k,1512) = lu(k,1512) - lu(k,147) * lu(k,1505)
         lu(k,1519) = lu(k,1519) - lu(k,148) * lu(k,1505)
         lu(k,1521) = lu(k,1521) - lu(k,149) * lu(k,1505)
         lu(k,1526) = lu(k,1526) - lu(k,150) * lu(k,1505)
                                                                        
         lu(k,151) = 1._r8 / lu(k,151)
         lu(k,152) = lu(k,152) * lu(k,151)
         lu(k,153) = lu(k,153) * lu(k,151)
         lu(k,154) = lu(k,154) * lu(k,151)
         lu(k,155) = lu(k,155) * lu(k,151)
         lu(k,1511) = lu(k,1511) - lu(k,152) * lu(k,1506)
         lu(k,1512) = lu(k,1512) - lu(k,153) * lu(k,1506)
         lu(k,1526) = lu(k,1526) - lu(k,154) * lu(k,1506)
         lu(k,1532) = lu(k,1532) - lu(k,155) * lu(k,1506)
                                                                        
         lu(k,156) = 1._r8 / lu(k,156)
         lu(k,157) = lu(k,157) * lu(k,156)
         lu(k,158) = lu(k,158) * lu(k,156)
         lu(k,159) = lu(k,159) * lu(k,156)
         lu(k,160) = lu(k,160) * lu(k,156)
         lu(k,1512) = lu(k,1512) - lu(k,157) * lu(k,1507)
         lu(k,1519) = lu(k,1519) - lu(k,158) * lu(k,1507)
         lu(k,1526) = lu(k,1526) - lu(k,159) * lu(k,1507)
         lu(k,1532) = lu(k,1532) - lu(k,160) * lu(k,1507)
                                                                        
         lu(k,161) = 1._r8 / lu(k,161)
         lu(k,162) = lu(k,162) * lu(k,161)
         lu(k,827) = lu(k,827) - lu(k,162) * lu(k,823)
         lu(k,889) = lu(k,889) - lu(k,162) * lu(k,882)
         lu(k,1978) = lu(k,1978) - lu(k,162) * lu(k,1961)
         lu(k,2017) = lu(k,2017) - lu(k,162) * lu(k,1986)
         lu(k,2153) = lu(k,2153) - lu(k,162) * lu(k,2138)
                                                                        
         lu(k,164) = 1._r8 / lu(k,164)
         lu(k,165) = lu(k,165) * lu(k,164)
         lu(k,166) = lu(k,166) * lu(k,164)
         lu(k,167) = lu(k,167) * lu(k,164)
         lu(k,168) = lu(k,168) * lu(k,164)
         lu(k,169) = lu(k,169) * lu(k,164)
         lu(k,1560) = lu(k,1560) - lu(k,165) * lu(k,1559)
         lu(k,1561) = lu(k,1561) - lu(k,166) * lu(k,1559)
         lu(k,1609) = lu(k,1609) - lu(k,167) * lu(k,1559)
         lu(k,1691) = lu(k,1691) - lu(k,168) * lu(k,1559)
         lu(k,1694) = lu(k,1694) - lu(k,169) * lu(k,1559)
                                                                        
         lu(k,170) = 1._r8 / lu(k,170)
         lu(k,171) = lu(k,171) * lu(k,170)
         lu(k,172) = lu(k,172) * lu(k,170)
         lu(k,173) = lu(k,173) * lu(k,170)
         lu(k,1605) = - lu(k,171) * lu(k,1560)
         lu(k,1668) = lu(k,1668) - lu(k,172) * lu(k,1560)
         lu(k,1694) = lu(k,1694) - lu(k,173) * lu(k,1560)
                                                                        
         lu(k,174) = 1._r8 / lu(k,174)
         lu(k,175) = lu(k,175) * lu(k,174)
         lu(k,176) = lu(k,176) * lu(k,174)
         lu(k,177) = lu(k,177) * lu(k,174)
         lu(k,178) = lu(k,178) * lu(k,174)
         lu(k,1604) = lu(k,1604) - lu(k,175) * lu(k,1561)
         lu(k,1606) = lu(k,1606) - lu(k,176) * lu(k,1561)
         lu(k,1691) = lu(k,1691) - lu(k,177) * lu(k,1561)
         lu(k,1694) = lu(k,1694) - lu(k,178) * lu(k,1561)
                                                                        
         lu(k,179) = 1._r8 / lu(k,179)
         lu(k,180) = lu(k,180) * lu(k,179)
         lu(k,181) = lu(k,181) * lu(k,179)
         lu(k,182) = lu(k,182) * lu(k,179)
         lu(k,1526) = lu(k,1526) - lu(k,180) * lu(k,1508)
         lu(k,1527) = lu(k,1527) - lu(k,181) * lu(k,1508)
         lu(k,1530) = lu(k,1530) - lu(k,182) * lu(k,1508)
         lu(k,1690) = - lu(k,180) * lu(k,1562)
         lu(k,1691) = lu(k,1691) - lu(k,181) * lu(k,1562)
         lu(k,1694) = lu(k,1694) - lu(k,182) * lu(k,1562)
                                                                        
         lu(k,183) = 1._r8 / lu(k,183)
         lu(k,184) = lu(k,184) * lu(k,183)
         lu(k,185) = lu(k,185) * lu(k,183)
         lu(k,493) = - lu(k,184) * lu(k,490)
         lu(k,495) = lu(k,495) - lu(k,185) * lu(k,490)
         lu(k,1526) = lu(k,1526) - lu(k,184) * lu(k,1509)
         lu(k,1529) = lu(k,1529) - lu(k,185) * lu(k,1509)
         lu(k,2191) = - lu(k,184) * lu(k,2162)
         lu(k,2194) = lu(k,2194) - lu(k,185) * lu(k,2162)
                                                                        
         lu(k,187) = 1._r8 / lu(k,187)
         lu(k,188) = lu(k,188) * lu(k,187)
         lu(k,189) = lu(k,189) * lu(k,187)
         lu(k,190) = lu(k,190) * lu(k,187)
         lu(k,191) = lu(k,191) * lu(k,187)
         lu(k,192) = lu(k,192) * lu(k,187)
         lu(k,193) = lu(k,193) * lu(k,187)
         lu(k,1564) = lu(k,1564) - lu(k,188) * lu(k,1563)
         lu(k,1565) = lu(k,1565) - lu(k,189) * lu(k,1563)
         lu(k,1602) = lu(k,1602) - lu(k,190) * lu(k,1563)
         lu(k,1638) = lu(k,1638) - lu(k,191) * lu(k,1563)
         lu(k,1691) = lu(k,1691) - lu(k,192) * lu(k,1563)
         lu(k,1694) = lu(k,1694) - lu(k,193) * lu(k,1563)
                                                                        
         lu(k,194) = 1._r8 / lu(k,194)
         lu(k,195) = lu(k,195) * lu(k,194)
         lu(k,196) = lu(k,196) * lu(k,194)
         lu(k,197) = lu(k,197) * lu(k,194)
         lu(k,198) = lu(k,198) * lu(k,194)
         lu(k,1604) = lu(k,1604) - lu(k,195) * lu(k,1564)
         lu(k,1606) = lu(k,1606) - lu(k,196) * lu(k,1564)
         lu(k,1691) = lu(k,1691) - lu(k,197) * lu(k,1564)
         lu(k,1694) = lu(k,1694) - lu(k,198) * lu(k,1564)
                                                                        
      end do
                                                                        
      end subroutine lu_fac03
                                                                        
      subroutine lu_fac04( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,199) = 1._r8 / lu(k,199)
         lu(k,200) = lu(k,200) * lu(k,199)
         lu(k,201) = lu(k,201) * lu(k,199)
         lu(k,202) = lu(k,202) * lu(k,199)
         lu(k,212) = - lu(k,200) * lu(k,207)
         lu(k,213) = - lu(k,201) * lu(k,207)
         lu(k,215) = lu(k,215) - lu(k,202) * lu(k,207)
         lu(k,1668) = lu(k,1668) - lu(k,200) * lu(k,1565)
         lu(k,1683) = lu(k,1683) - lu(k,201) * lu(k,1565)
         lu(k,1694) = lu(k,1694) - lu(k,202) * lu(k,1565)
                                                                        
         lu(k,203) = 1._r8 / lu(k,203)
         lu(k,204) = lu(k,204) * lu(k,203)
         lu(k,205) = lu(k,205) * lu(k,203)
         lu(k,1173) = lu(k,1173) - lu(k,204) * lu(k,1166)
         lu(k,1177) = - lu(k,205) * lu(k,1166)
         lu(k,1673) = lu(k,1673) - lu(k,204) * lu(k,1566)
         lu(k,1691) = lu(k,1691) - lu(k,205) * lu(k,1566)
         lu(k,1930) = lu(k,1930) - lu(k,204) * lu(k,1866)
         lu(k,1947) = lu(k,1947) - lu(k,205) * lu(k,1866)
                                                                        
         lu(k,208) = 1._r8 / lu(k,208)
         lu(k,209) = lu(k,209) * lu(k,208)
         lu(k,210) = lu(k,210) * lu(k,208)
         lu(k,211) = lu(k,211) * lu(k,208)
         lu(k,212) = lu(k,212) * lu(k,208)
         lu(k,213) = lu(k,213) * lu(k,208)
         lu(k,214) = lu(k,214) * lu(k,208)
         lu(k,215) = lu(k,215) * lu(k,208)
         lu(k,1568) = lu(k,1568) - lu(k,209) * lu(k,1567)
         lu(k,1602) = lu(k,1602) - lu(k,210) * lu(k,1567)
         lu(k,1639) = lu(k,1639) - lu(k,211) * lu(k,1567)
         lu(k,1668) = lu(k,1668) - lu(k,212) * lu(k,1567)
         lu(k,1683) = lu(k,1683) - lu(k,213) * lu(k,1567)
         lu(k,1691) = lu(k,1691) - lu(k,214) * lu(k,1567)
         lu(k,1694) = lu(k,1694) - lu(k,215) * lu(k,1567)
                                                                        
         lu(k,216) = 1._r8 / lu(k,216)
         lu(k,217) = lu(k,217) * lu(k,216)
         lu(k,218) = lu(k,218) * lu(k,216)
         lu(k,219) = lu(k,219) * lu(k,216)
         lu(k,220) = lu(k,220) * lu(k,216)
         lu(k,1606) = lu(k,1606) - lu(k,217) * lu(k,1568)
         lu(k,1611) = lu(k,1611) - lu(k,218) * lu(k,1568)
         lu(k,1691) = lu(k,1691) - lu(k,219) * lu(k,1568)
         lu(k,1694) = lu(k,1694) - lu(k,220) * lu(k,1568)
                                                                        
         lu(k,221) = 1._r8 / lu(k,221)
         lu(k,222) = lu(k,222) * lu(k,221)
         lu(k,223) = lu(k,223) * lu(k,221)
         lu(k,224) = lu(k,224) * lu(k,221)
         lu(k,225) = lu(k,225) * lu(k,221)
         lu(k,1511) = lu(k,1511) - lu(k,222) * lu(k,1510)
         lu(k,1526) = lu(k,1526) - lu(k,223) * lu(k,1510)
         lu(k,1527) = lu(k,1527) - lu(k,224) * lu(k,1510)
         lu(k,1532) = lu(k,1532) - lu(k,225) * lu(k,1510)
         lu(k,1570) = lu(k,1570) - lu(k,222) * lu(k,1569)
         lu(k,1690) = lu(k,1690) - lu(k,223) * lu(k,1569)
         lu(k,1691) = lu(k,1691) - lu(k,224) * lu(k,1569)
         lu(k,1696) = lu(k,1696) - lu(k,225) * lu(k,1569)
                                                                        
         lu(k,226) = 1._r8 / lu(k,226)
         lu(k,227) = lu(k,227) * lu(k,226)
         lu(k,228) = lu(k,228) * lu(k,226)
         lu(k,229) = lu(k,229) * lu(k,226)
         lu(k,1519) = lu(k,1519) - lu(k,227) * lu(k,1511)
         lu(k,1526) = lu(k,1526) - lu(k,228) * lu(k,1511)
         lu(k,1532) = lu(k,1532) - lu(k,229) * lu(k,1511)
         lu(k,1651) = - lu(k,227) * lu(k,1570)
         lu(k,1690) = lu(k,1690) - lu(k,228) * lu(k,1570)
         lu(k,1696) = lu(k,1696) - lu(k,229) * lu(k,1570)
                                                                        
         lu(k,230) = 1._r8 / lu(k,230)
         lu(k,231) = lu(k,231) * lu(k,230)
         lu(k,232) = lu(k,232) * lu(k,230)
         lu(k,233) = lu(k,233) * lu(k,230)
         lu(k,234) = lu(k,234) * lu(k,230)
         lu(k,1275) = - lu(k,231) * lu(k,1272)
         lu(k,1285) = - lu(k,232) * lu(k,1272)
         lu(k,1295) = - lu(k,233) * lu(k,1272)
         lu(k,1298) = lu(k,1298) - lu(k,234) * lu(k,1272)
         lu(k,1622) = - lu(k,231) * lu(k,1571)
         lu(k,1673) = lu(k,1673) - lu(k,232) * lu(k,1571)
         lu(k,1691) = lu(k,1691) - lu(k,233) * lu(k,1571)
         lu(k,1694) = lu(k,1694) - lu(k,234) * lu(k,1571)
                                                                        
         lu(k,235) = 1._r8 / lu(k,235)
         lu(k,236) = lu(k,236) * lu(k,235)
         lu(k,237) = lu(k,237) * lu(k,235)
         lu(k,249) = - lu(k,236) * lu(k,247)
         lu(k,250) = lu(k,250) - lu(k,237) * lu(k,247)
         lu(k,305) = - lu(k,236) * lu(k,303)
         lu(k,306) = lu(k,306) - lu(k,237) * lu(k,303)
         lu(k,1519) = lu(k,1519) - lu(k,236) * lu(k,1512)
         lu(k,1526) = lu(k,1526) - lu(k,237) * lu(k,1512)
         lu(k,1651) = lu(k,1651) - lu(k,236) * lu(k,1572)
         lu(k,1690) = lu(k,1690) - lu(k,237) * lu(k,1572)
                                                                        
         lu(k,238) = 1._r8 / lu(k,238)
         lu(k,239) = lu(k,239) * lu(k,238)
         lu(k,240) = lu(k,240) * lu(k,238)
         lu(k,944) = - lu(k,239) * lu(k,941)
         lu(k,960) = lu(k,960) - lu(k,240) * lu(k,941)
         lu(k,993) = - lu(k,239) * lu(k,990)
         lu(k,1010) = lu(k,1010) - lu(k,240) * lu(k,990)
         lu(k,1654) = lu(k,1654) - lu(k,239) * lu(k,1573)
         lu(k,1691) = lu(k,1691) - lu(k,240) * lu(k,1573)
         lu(k,2093) = - lu(k,239) * lu(k,2082)
         lu(k,2125) = lu(k,2125) - lu(k,240) * lu(k,2082)
                                                                        
      end do
                                                                        
      end subroutine lu_fac04
                                                                        
      subroutine lu_fac05( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,241) = 1._r8 / lu(k,241)
         lu(k,242) = lu(k,242) * lu(k,241)
         lu(k,243) = lu(k,243) * lu(k,241)
         lu(k,779) = lu(k,779) - lu(k,242) * lu(k,777)
         lu(k,781) = - lu(k,243) * lu(k,777)
         lu(k,1968) = lu(k,1968) - lu(k,242) * lu(k,1962)
         lu(k,1978) = lu(k,1978) - lu(k,243) * lu(k,1962)
         lu(k,2143) = - lu(k,242) * lu(k,2139)
         lu(k,2153) = lu(k,2153) - lu(k,243) * lu(k,2139)
         lu(k,2211) = lu(k,2211) - lu(k,242) * lu(k,2205)
         lu(k,2221) = lu(k,2221) - lu(k,243) * lu(k,2205)
                                                                        
         lu(k,244) = 1._r8 / lu(k,244)
         lu(k,245) = lu(k,245) * lu(k,244)
         lu(k,246) = lu(k,246) * lu(k,244)
         lu(k,1118) = - lu(k,245) * lu(k,1107)
         lu(k,1124) = lu(k,1124) - lu(k,246) * lu(k,1107)
         lu(k,1156) = lu(k,1156) - lu(k,245) * lu(k,1146)
         lu(k,1163) = lu(k,1163) - lu(k,246) * lu(k,1146)
         lu(k,1687) = lu(k,1687) - lu(k,245) * lu(k,1574)
         lu(k,1700) = lu(k,1700) - lu(k,246) * lu(k,1574)
         lu(k,1836) = - lu(k,245) * lu(k,1773)
         lu(k,1849) = lu(k,1849) - lu(k,246) * lu(k,1773)
                                                                        
         lu(k,248) = 1._r8 / lu(k,248)
         lu(k,249) = lu(k,249) * lu(k,248)
         lu(k,250) = lu(k,250) * lu(k,248)
         lu(k,251) = lu(k,251) * lu(k,248)
         lu(k,252) = lu(k,252) * lu(k,248)
         lu(k,1519) = lu(k,1519) - lu(k,249) * lu(k,1513)
         lu(k,1526) = lu(k,1526) - lu(k,250) * lu(k,1513)
         lu(k,1527) = lu(k,1527) - lu(k,251) * lu(k,1513)
         lu(k,1532) = lu(k,1532) - lu(k,252) * lu(k,1513)
         lu(k,1651) = lu(k,1651) - lu(k,249) * lu(k,1575)
         lu(k,1690) = lu(k,1690) - lu(k,250) * lu(k,1575)
         lu(k,1691) = lu(k,1691) - lu(k,251) * lu(k,1575)
         lu(k,1696) = lu(k,1696) - lu(k,252) * lu(k,1575)
                                                                        
         lu(k,253) = 1._r8 / lu(k,253)
         lu(k,254) = lu(k,254) * lu(k,253)
         lu(k,255) = lu(k,255) * lu(k,253)
         lu(k,341) = - lu(k,254) * lu(k,338)
         lu(k,342) = lu(k,342) - lu(k,255) * lu(k,338)
         lu(k,436) = - lu(k,254) * lu(k,433)
         lu(k,437) = - lu(k,255) * lu(k,433)
         lu(k,1614) = lu(k,1614) - lu(k,254) * lu(k,1576)
         lu(k,1691) = lu(k,1691) - lu(k,255) * lu(k,1576)
         lu(k,1788) = lu(k,1788) - lu(k,254) * lu(k,1774)
         lu(k,1840) = lu(k,1840) - lu(k,255) * lu(k,1774)
                                                                        
         lu(k,257) = 1._r8 / lu(k,257)
         lu(k,258) = lu(k,258) * lu(k,257)
         lu(k,259) = lu(k,259) * lu(k,257)
         lu(k,260) = lu(k,260) * lu(k,257)
         lu(k,832) = lu(k,832) - lu(k,258) * lu(k,831)
         lu(k,836) = lu(k,836) - lu(k,259) * lu(k,831)
         lu(k,837) = - lu(k,260) * lu(k,831)
         lu(k,1648) = lu(k,1648) - lu(k,258) * lu(k,1577)
         lu(k,1702) = lu(k,1702) - lu(k,259) * lu(k,1577)
         lu(k,1703) = lu(k,1703) - lu(k,260) * lu(k,1577)
         lu(k,2263) = - lu(k,258) * lu(k,2261)
         lu(k,2284) = lu(k,2284) - lu(k,259) * lu(k,2261)
         lu(k,2285) = lu(k,2285) - lu(k,260) * lu(k,2261)
                                                                        
         lu(k,261) = 1._r8 / lu(k,261)
         lu(k,262) = lu(k,262) * lu(k,261)
         lu(k,263) = lu(k,263) * lu(k,261)
         lu(k,264) = lu(k,264) * lu(k,261)
         lu(k,816) = lu(k,816) - lu(k,262) * lu(k,812)
         lu(k,818) = - lu(k,263) * lu(k,812)
         lu(k,820) = lu(k,820) - lu(k,264) * lu(k,812)
         lu(k,1666) = lu(k,1666) - lu(k,262) * lu(k,1578)
         lu(k,1691) = lu(k,1691) - lu(k,263) * lu(k,1578)
         lu(k,1694) = lu(k,1694) - lu(k,264) * lu(k,1578)
         lu(k,2043) = lu(k,2043) - lu(k,262) * lu(k,2025)
         lu(k,2064) = - lu(k,263) * lu(k,2025)
         lu(k,2067) = lu(k,2067) - lu(k,264) * lu(k,2025)
                                                                        
         lu(k,265) = 1._r8 / lu(k,265)
         lu(k,266) = lu(k,266) * lu(k,265)
         lu(k,267) = lu(k,267) * lu(k,265)
         lu(k,268) = lu(k,268) * lu(k,265)
         lu(k,624) = lu(k,624) - lu(k,266) * lu(k,623)
         lu(k,625) = lu(k,625) - lu(k,267) * lu(k,623)
         lu(k,626) = - lu(k,268) * lu(k,623)
         lu(k,1606) = lu(k,1606) - lu(k,266) * lu(k,1579)
         lu(k,1627) = lu(k,1627) - lu(k,267) * lu(k,1579)
         lu(k,1691) = lu(k,1691) - lu(k,268) * lu(k,1579)
         lu(k,1886) = - lu(k,266) * lu(k,1867)
         lu(k,1894) = lu(k,1894) - lu(k,267) * lu(k,1867)
         lu(k,1947) = lu(k,1947) - lu(k,268) * lu(k,1867)
                                                                        
         lu(k,269) = 1._r8 / lu(k,269)
         lu(k,270) = lu(k,270) * lu(k,269)
         lu(k,271) = lu(k,271) * lu(k,269)
         lu(k,272) = lu(k,272) * lu(k,269)
         lu(k,715) = - lu(k,270) * lu(k,713)
         lu(k,716) = lu(k,716) - lu(k,271) * lu(k,713)
         lu(k,719) = lu(k,719) - lu(k,272) * lu(k,713)
         lu(k,1818) = lu(k,1818) - lu(k,270) * lu(k,1775)
         lu(k,1838) = lu(k,1838) - lu(k,271) * lu(k,1775)
         lu(k,1843) = lu(k,1843) - lu(k,272) * lu(k,1775)
         lu(k,1925) = - lu(k,270) * lu(k,1868)
         lu(k,1945) = lu(k,1945) - lu(k,271) * lu(k,1868)
         lu(k,1950) = lu(k,1950) - lu(k,272) * lu(k,1868)
                                                                        
         lu(k,273) = 1._r8 / lu(k,273)
         lu(k,274) = lu(k,274) * lu(k,273)
         lu(k,275) = lu(k,275) * lu(k,273)
         lu(k,276) = lu(k,276) * lu(k,273)
         lu(k,277) = lu(k,277) * lu(k,273)
         lu(k,278) = lu(k,278) * lu(k,273)
         lu(k,1744) = lu(k,1744) - lu(k,274) * lu(k,1706)
         lu(k,1749) = lu(k,1749) - lu(k,275) * lu(k,1706)
         lu(k,1750) = lu(k,1750) - lu(k,276) * lu(k,1706)
         lu(k,1757) = lu(k,1757) - lu(k,277) * lu(k,1706)
         lu(k,1759) = lu(k,1759) - lu(k,278) * lu(k,1706)
         lu(k,2188) = lu(k,2188) - lu(k,274) * lu(k,2163)
         lu(k,2193) = lu(k,2193) - lu(k,275) * lu(k,2163)
         lu(k,2194) = lu(k,2194) - lu(k,276) * lu(k,2163)
         lu(k,2201) = lu(k,2201) - lu(k,277) * lu(k,2163)
         lu(k,2203) = lu(k,2203) - lu(k,278) * lu(k,2163)
                                                                        
      end do
                                                                        
      end subroutine lu_fac05
                                                                        
      subroutine lu_fac06( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,279) = 1._r8 / lu(k,279)
         lu(k,280) = lu(k,280) * lu(k,279)
         lu(k,281) = lu(k,281) * lu(k,279)
         lu(k,658) = - lu(k,280) * lu(k,652)
         lu(k,664) = lu(k,664) - lu(k,281) * lu(k,652)
         lu(k,704) = - lu(k,280) * lu(k,697)
         lu(k,711) = lu(k,711) - lu(k,281) * lu(k,697)
         lu(k,733) = - lu(k,280) * lu(k,727)
         lu(k,740) = lu(k,740) - lu(k,281) * lu(k,727)
         lu(k,749) = - lu(k,280) * lu(k,742)
         lu(k,757) = lu(k,757) - lu(k,281) * lu(k,742)
         lu(k,1801) = lu(k,1801) - lu(k,280) * lu(k,1776)
         lu(k,1843) = lu(k,1843) - lu(k,281) * lu(k,1776)
                                                                        
         lu(k,282) = 1._r8 / lu(k,282)
         lu(k,283) = lu(k,283) * lu(k,282)
         lu(k,284) = lu(k,284) * lu(k,282)
         lu(k,864) = lu(k,864) - lu(k,283) * lu(k,863)
         lu(k,868) = lu(k,868) - lu(k,284) * lu(k,863)
         lu(k,1400) = lu(k,1400) - lu(k,283) * lu(k,1399)
         lu(k,1404) = lu(k,1404) - lu(k,284) * lu(k,1399)
         lu(k,1427) = lu(k,1427) - lu(k,283) * lu(k,1425)
         lu(k,1432) = lu(k,1432) - lu(k,284) * lu(k,1425)
         lu(k,1444) = lu(k,1444) - lu(k,283) * lu(k,1443)
         lu(k,1448) = - lu(k,284) * lu(k,1443)
         lu(k,2264) = lu(k,2264) - lu(k,283) * lu(k,2262)
         lu(k,2270) = lu(k,2270) - lu(k,284) * lu(k,2262)
                                                                        
         lu(k,285) = 1._r8 / lu(k,285)
         lu(k,286) = lu(k,286) * lu(k,285)
         lu(k,287) = lu(k,287) * lu(k,285)
         lu(k,288) = lu(k,288) * lu(k,285)
         lu(k,289) = lu(k,289) * lu(k,285)
         lu(k,290) = lu(k,290) * lu(k,285)
         lu(k,1657) = lu(k,1657) - lu(k,286) * lu(k,1580)
         lu(k,1660) = lu(k,1660) - lu(k,287) * lu(k,1580)
         lu(k,1668) = lu(k,1668) - lu(k,288) * lu(k,1580)
         lu(k,1691) = lu(k,1691) - lu(k,289) * lu(k,1580)
         lu(k,1694) = lu(k,1694) - lu(k,290) * lu(k,1580)
         lu(k,1999) = - lu(k,286) * lu(k,1987)
         lu(k,2000) = - lu(k,287) * lu(k,1987)
         lu(k,2003) = lu(k,2003) - lu(k,288) * lu(k,1987)
         lu(k,2012) = lu(k,2012) - lu(k,289) * lu(k,1987)
         lu(k,2015) = lu(k,2015) - lu(k,290) * lu(k,1987)
                                                                        
         lu(k,291) = 1._r8 / lu(k,291)
         lu(k,292) = lu(k,292) * lu(k,291)
         lu(k,293) = lu(k,293) * lu(k,291)
         lu(k,294) = lu(k,294) * lu(k,291)
         lu(k,295) = lu(k,295) * lu(k,291)
         lu(k,296) = lu(k,296) * lu(k,291)
         lu(k,1646) = lu(k,1646) - lu(k,292) * lu(k,1581)
         lu(k,1691) = lu(k,1691) - lu(k,293) * lu(k,1581)
         lu(k,1696) = lu(k,1696) - lu(k,294) * lu(k,1581)
         lu(k,1699) = lu(k,1699) - lu(k,295) * lu(k,1581)
         lu(k,1703) = lu(k,1703) - lu(k,296) * lu(k,1581)
         lu(k,1996) = lu(k,1996) - lu(k,292) * lu(k,1988)
         lu(k,2012) = lu(k,2012) - lu(k,293) * lu(k,1988)
         lu(k,2017) = lu(k,2017) - lu(k,294) * lu(k,1988)
         lu(k,2020) = lu(k,2020) - lu(k,295) * lu(k,1988)
         lu(k,2024) = - lu(k,296) * lu(k,1988)
                                                                        
         lu(k,297) = 1._r8 / lu(k,297)
         lu(k,298) = lu(k,298) * lu(k,297)
         lu(k,299) = lu(k,299) * lu(k,297)
         lu(k,300) = lu(k,300) * lu(k,297)
         lu(k,301) = lu(k,301) * lu(k,297)
         lu(k,302) = lu(k,302) * lu(k,297)
         lu(k,1648) = lu(k,1648) - lu(k,298) * lu(k,1582)
         lu(k,1687) = lu(k,1687) - lu(k,299) * lu(k,1582)
         lu(k,1691) = lu(k,1691) - lu(k,300) * lu(k,1582)
         lu(k,1692) = lu(k,1692) - lu(k,301) * lu(k,1582)
         lu(k,1694) = lu(k,1694) - lu(k,302) * lu(k,1582)
         lu(k,1711) = lu(k,1711) - lu(k,298) * lu(k,1707)
         lu(k,1744) = lu(k,1744) - lu(k,299) * lu(k,1707)
         lu(k,1748) = lu(k,1748) - lu(k,300) * lu(k,1707)
         lu(k,1749) = lu(k,1749) - lu(k,301) * lu(k,1707)
         lu(k,1751) = lu(k,1751) - lu(k,302) * lu(k,1707)
                                                                        
         lu(k,304) = 1._r8 / lu(k,304)
         lu(k,305) = lu(k,305) * lu(k,304)
         lu(k,306) = lu(k,306) * lu(k,304)
         lu(k,307) = lu(k,307) * lu(k,304)
         lu(k,308) = lu(k,308) * lu(k,304)
         lu(k,309) = lu(k,309) * lu(k,304)
         lu(k,1519) = lu(k,1519) - lu(k,305) * lu(k,1514)
         lu(k,1526) = lu(k,1526) - lu(k,306) * lu(k,1514)
         lu(k,1527) = lu(k,1527) - lu(k,307) * lu(k,1514)
         lu(k,1532) = lu(k,1532) - lu(k,308) * lu(k,1514)
         lu(k,1539) = lu(k,1539) - lu(k,309) * lu(k,1514)
         lu(k,1651) = lu(k,1651) - lu(k,305) * lu(k,1583)
         lu(k,1690) = lu(k,1690) - lu(k,306) * lu(k,1583)
         lu(k,1691) = lu(k,1691) - lu(k,307) * lu(k,1583)
         lu(k,1696) = lu(k,1696) - lu(k,308) * lu(k,1583)
         lu(k,1703) = lu(k,1703) - lu(k,309) * lu(k,1583)
                                                                        
         lu(k,310) = 1._r8 / lu(k,310)
         lu(k,311) = lu(k,311) * lu(k,310)
         lu(k,312) = lu(k,312) * lu(k,310)
         lu(k,313) = lu(k,313) * lu(k,310)
         lu(k,314) = lu(k,314) * lu(k,310)
         lu(k,1311) = lu(k,1311) - lu(k,311) * lu(k,1304)
         lu(k,1312) = - lu(k,312) * lu(k,1304)
         lu(k,1316) = - lu(k,313) * lu(k,1304)
         lu(k,1319) = lu(k,1319) - lu(k,314) * lu(k,1304)
         lu(k,1680) = lu(k,1680) - lu(k,311) * lu(k,1584)
         lu(k,1682) = lu(k,1682) - lu(k,312) * lu(k,1584)
         lu(k,1691) = lu(k,1691) - lu(k,313) * lu(k,1584)
         lu(k,1694) = lu(k,1694) - lu(k,314) * lu(k,1584)
         lu(k,1936) = lu(k,1936) - lu(k,311) * lu(k,1869)
         lu(k,1938) = lu(k,1938) - lu(k,312) * lu(k,1869)
         lu(k,1947) = lu(k,1947) - lu(k,313) * lu(k,1869)
         lu(k,1950) = lu(k,1950) - lu(k,314) * lu(k,1869)
                                                                        
         lu(k,315) = 1._r8 / lu(k,315)
         lu(k,316) = lu(k,316) * lu(k,315)
         lu(k,317) = lu(k,317) * lu(k,315)
         lu(k,318) = lu(k,318) * lu(k,315)
         lu(k,319) = lu(k,319) * lu(k,315)
         lu(k,677) = lu(k,677) - lu(k,316) * lu(k,676)
         lu(k,678) = lu(k,678) - lu(k,317) * lu(k,676)
         lu(k,679) = lu(k,679) - lu(k,318) * lu(k,676)
         lu(k,680) = lu(k,680) - lu(k,319) * lu(k,676)
         lu(k,1632) = lu(k,1632) - lu(k,316) * lu(k,1585)
         lu(k,1666) = lu(k,1666) - lu(k,317) * lu(k,1585)
         lu(k,1683) = lu(k,1683) - lu(k,318) * lu(k,1585)
         lu(k,1691) = lu(k,1691) - lu(k,319) * lu(k,1585)
         lu(k,1898) = lu(k,1898) - lu(k,316) * lu(k,1870)
         lu(k,1924) = lu(k,1924) - lu(k,317) * lu(k,1870)
         lu(k,1939) = lu(k,1939) - lu(k,318) * lu(k,1870)
         lu(k,1947) = lu(k,1947) - lu(k,319) * lu(k,1870)
                                                                        
         lu(k,320) = 1._r8 / lu(k,320)
         lu(k,321) = lu(k,321) * lu(k,320)
         lu(k,322) = lu(k,322) * lu(k,320)
         lu(k,1285) = lu(k,1285) - lu(k,321) * lu(k,1273)
         lu(k,1295) = lu(k,1295) - lu(k,322) * lu(k,1273)
         lu(k,1377) = lu(k,1377) - lu(k,321) * lu(k,1367)
         lu(k,1390) = lu(k,1390) - lu(k,322) * lu(k,1367)
         lu(k,1673) = lu(k,1673) - lu(k,321) * lu(k,1586)
         lu(k,1691) = lu(k,1691) - lu(k,322) * lu(k,1586)
         lu(k,1731) = lu(k,1731) - lu(k,321) * lu(k,1708)
         lu(k,1748) = lu(k,1748) - lu(k,322) * lu(k,1708)
         lu(k,1824) = lu(k,1824) - lu(k,321) * lu(k,1777)
         lu(k,1840) = lu(k,1840) - lu(k,322) * lu(k,1777)
         lu(k,2049) = lu(k,2049) - lu(k,321) * lu(k,2026)
         lu(k,2064) = lu(k,2064) - lu(k,322) * lu(k,2026)
                                                                        
      end do
                                                                        
      end subroutine lu_fac06
                                                                        
      subroutine lu_fac07( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,323) = 1._r8 / lu(k,323)
         lu(k,324) = lu(k,324) * lu(k,323)
         lu(k,325) = lu(k,325) * lu(k,323)
         lu(k,326) = lu(k,326) * lu(k,323)
         lu(k,327) = lu(k,327) * lu(k,323)
         lu(k,707) = - lu(k,324) * lu(k,698)
         lu(k,708) = lu(k,708) - lu(k,325) * lu(k,698)
         lu(k,709) = - lu(k,326) * lu(k,698)
         lu(k,711) = lu(k,711) - lu(k,327) * lu(k,698)
         lu(k,752) = - lu(k,324) * lu(k,743)
         lu(k,753) = lu(k,753) - lu(k,325) * lu(k,743)
         lu(k,754) = - lu(k,326) * lu(k,743)
         lu(k,757) = lu(k,757) - lu(k,327) * lu(k,743)
         lu(k,1819) = lu(k,1819) - lu(k,324) * lu(k,1778)
         lu(k,1827) = lu(k,1827) - lu(k,325) * lu(k,1778)
         lu(k,1833) = lu(k,1833) - lu(k,326) * lu(k,1778)
         lu(k,1843) = lu(k,1843) - lu(k,327) * lu(k,1778)
                                                                        
         lu(k,328) = 1._r8 / lu(k,328)
         lu(k,329) = lu(k,329) * lu(k,328)
         lu(k,330) = lu(k,330) * lu(k,328)
         lu(k,331) = lu(k,331) * lu(k,328)
         lu(k,332) = lu(k,332) * lu(k,328)
         lu(k,1231) = lu(k,1231) - lu(k,329) * lu(k,1229)
         lu(k,1232) = lu(k,1232) - lu(k,330) * lu(k,1229)
         lu(k,1238) = lu(k,1238) - lu(k,331) * lu(k,1229)
         lu(k,1243) = lu(k,1243) - lu(k,332) * lu(k,1229)
         lu(k,1965) = lu(k,1965) - lu(k,329) * lu(k,1963)
         lu(k,1967) = lu(k,1967) - lu(k,330) * lu(k,1963)
         lu(k,1977) = lu(k,1977) - lu(k,331) * lu(k,1963)
         lu(k,1984) = lu(k,1984) - lu(k,332) * lu(k,1963)
         lu(k,2209) = lu(k,2209) - lu(k,329) * lu(k,2206)
         lu(k,2210) = lu(k,2210) - lu(k,330) * lu(k,2206)
         lu(k,2220) = lu(k,2220) - lu(k,331) * lu(k,2206)
         lu(k,2227) = lu(k,2227) - lu(k,332) * lu(k,2206)
                                                                        
         lu(k,333) = 1._r8 / lu(k,333)
         lu(k,334) = lu(k,334) * lu(k,333)
         lu(k,335) = lu(k,335) * lu(k,333)
         lu(k,336) = lu(k,336) * lu(k,333)
         lu(k,337) = lu(k,337) * lu(k,333)
         lu(k,454) = lu(k,454) - lu(k,334) * lu(k,453)
         lu(k,455) = lu(k,455) - lu(k,335) * lu(k,453)
         lu(k,456) = - lu(k,336) * lu(k,453)
         lu(k,458) = lu(k,458) - lu(k,337) * lu(k,453)
         lu(k,1604) = lu(k,1604) - lu(k,334) * lu(k,1587)
         lu(k,1660) = lu(k,1660) - lu(k,335) * lu(k,1587)
         lu(k,1691) = lu(k,1691) - lu(k,336) * lu(k,1587)
         lu(k,1694) = lu(k,1694) - lu(k,337) * lu(k,1587)
         lu(k,1884) = lu(k,1884) - lu(k,334) * lu(k,1871)
         lu(k,1921) = lu(k,1921) - lu(k,335) * lu(k,1871)
         lu(k,1947) = lu(k,1947) - lu(k,336) * lu(k,1871)
         lu(k,1950) = lu(k,1950) - lu(k,337) * lu(k,1871)
                                                                        
         lu(k,339) = 1._r8 / lu(k,339)
         lu(k,340) = lu(k,340) * lu(k,339)
         lu(k,341) = lu(k,341) * lu(k,339)
         lu(k,342) = lu(k,342) * lu(k,339)
         lu(k,343) = lu(k,343) * lu(k,339)
         lu(k,435) = lu(k,435) - lu(k,340) * lu(k,434)
         lu(k,436) = lu(k,436) - lu(k,341) * lu(k,434)
         lu(k,437) = lu(k,437) - lu(k,342) * lu(k,434)
         lu(k,439) = lu(k,439) - lu(k,343) * lu(k,434)
         lu(k,1602) = lu(k,1602) - lu(k,340) * lu(k,1588)
         lu(k,1614) = lu(k,1614) - lu(k,341) * lu(k,1588)
         lu(k,1691) = lu(k,1691) - lu(k,342) * lu(k,1588)
         lu(k,1694) = lu(k,1694) - lu(k,343) * lu(k,1588)
         lu(k,1882) = lu(k,1882) - lu(k,340) * lu(k,1872)
         lu(k,1890) = lu(k,1890) - lu(k,341) * lu(k,1872)
         lu(k,1947) = lu(k,1947) - lu(k,342) * lu(k,1872)
         lu(k,1950) = lu(k,1950) - lu(k,343) * lu(k,1872)
                                                                        
         lu(k,344) = 1._r8 / lu(k,344)
         lu(k,345) = lu(k,345) * lu(k,344)
         lu(k,346) = lu(k,346) * lu(k,344)
         lu(k,347) = lu(k,347) * lu(k,344)
         lu(k,348) = lu(k,348) * lu(k,344)
         lu(k,815) = lu(k,815) - lu(k,345) * lu(k,813)
         lu(k,816) = lu(k,816) - lu(k,346) * lu(k,813)
         lu(k,818) = lu(k,818) - lu(k,347) * lu(k,813)
         lu(k,820) = lu(k,820) - lu(k,348) * lu(k,813)
         lu(k,1646) = lu(k,1646) - lu(k,345) * lu(k,1589)
         lu(k,1666) = lu(k,1666) - lu(k,346) * lu(k,1589)
         lu(k,1691) = lu(k,1691) - lu(k,347) * lu(k,1589)
         lu(k,1694) = lu(k,1694) - lu(k,348) * lu(k,1589)
         lu(k,1912) = lu(k,1912) - lu(k,345) * lu(k,1873)
         lu(k,1924) = lu(k,1924) - lu(k,346) * lu(k,1873)
         lu(k,1947) = lu(k,1947) - lu(k,347) * lu(k,1873)
         lu(k,1950) = lu(k,1950) - lu(k,348) * lu(k,1873)
                                                                        
         lu(k,349) = 1._r8 / lu(k,349)
         lu(k,350) = lu(k,350) * lu(k,349)
         lu(k,351) = lu(k,351) * lu(k,349)
         lu(k,352) = lu(k,352) * lu(k,349)
         lu(k,353) = lu(k,353) * lu(k,349)
         lu(k,354) = lu(k,354) * lu(k,349)
         lu(k,355) = lu(k,355) * lu(k,349)
         lu(k,356) = lu(k,356) * lu(k,349)
         lu(k,1619) = lu(k,1619) - lu(k,350) * lu(k,1590)
         lu(k,1656) = lu(k,1656) - lu(k,351) * lu(k,1590)
         lu(k,1666) = lu(k,1666) - lu(k,352) * lu(k,1590)
         lu(k,1689) = lu(k,1689) - lu(k,353) * lu(k,1590)
         lu(k,1691) = lu(k,1691) - lu(k,354) * lu(k,1590)
         lu(k,1692) = lu(k,1692) - lu(k,355) * lu(k,1590)
         lu(k,1700) = lu(k,1700) - lu(k,356) * lu(k,1590)
         lu(k,1710) = - lu(k,350) * lu(k,1709)
         lu(k,1714) = lu(k,1714) - lu(k,351) * lu(k,1709)
         lu(k,1724) = lu(k,1724) - lu(k,352) * lu(k,1709)
         lu(k,1746) = lu(k,1746) - lu(k,353) * lu(k,1709)
         lu(k,1748) = lu(k,1748) - lu(k,354) * lu(k,1709)
         lu(k,1749) = lu(k,1749) - lu(k,355) * lu(k,1709)
         lu(k,1757) = lu(k,1757) - lu(k,356) * lu(k,1709)
                                                                        
         lu(k,357) = 1._r8 / lu(k,357)
         lu(k,358) = lu(k,358) * lu(k,357)
         lu(k,359) = lu(k,359) * lu(k,357)
         lu(k,360) = lu(k,360) * lu(k,357)
         lu(k,361) = lu(k,361) * lu(k,357)
         lu(k,362) = lu(k,362) * lu(k,357)
         lu(k,363) = lu(k,363) * lu(k,357)
         lu(k,364) = lu(k,364) * lu(k,357)
         lu(k,1613) = lu(k,1613) - lu(k,358) * lu(k,1591)
         lu(k,1648) = lu(k,1648) - lu(k,359) * lu(k,1591)
         lu(k,1668) = lu(k,1668) - lu(k,360) * lu(k,1591)
         lu(k,1677) = lu(k,1677) - lu(k,361) * lu(k,1591)
         lu(k,1688) = lu(k,1688) - lu(k,362) * lu(k,1591)
         lu(k,1691) = lu(k,1691) - lu(k,363) * lu(k,1591)
         lu(k,1702) = lu(k,1702) - lu(k,364) * lu(k,1591)
         lu(k,2230) = - lu(k,358) * lu(k,2229)
         lu(k,2236) = - lu(k,359) * lu(k,2229)
         lu(k,2238) = lu(k,2238) - lu(k,360) * lu(k,2229)
         lu(k,2239) = lu(k,2239) - lu(k,361) * lu(k,2229)
         lu(k,2244) = lu(k,2244) - lu(k,362) * lu(k,2229)
         lu(k,2247) = lu(k,2247) - lu(k,363) * lu(k,2229)
         lu(k,2258) = lu(k,2258) - lu(k,364) * lu(k,2229)
                                                                        
         lu(k,365) = 1._r8 / lu(k,365)
         lu(k,366) = lu(k,366) * lu(k,365)
         lu(k,367) = lu(k,367) * lu(k,365)
         lu(k,368) = lu(k,368) * lu(k,365)
         lu(k,369) = lu(k,369) * lu(k,365)
         lu(k,370) = lu(k,370) * lu(k,365)
         lu(k,371) = lu(k,371) * lu(k,365)
         lu(k,372) = lu(k,372) * lu(k,365)
         lu(k,1668) = lu(k,1668) - lu(k,366) * lu(k,1592)
         lu(k,1691) = lu(k,1691) - lu(k,367) * lu(k,1592)
         lu(k,1694) = lu(k,1694) - lu(k,368) * lu(k,1592)
         lu(k,1696) = lu(k,1696) - lu(k,369) * lu(k,1592)
         lu(k,1697) = lu(k,1697) - lu(k,370) * lu(k,1592)
         lu(k,1699) = lu(k,1699) - lu(k,371) * lu(k,1592)
         lu(k,1703) = lu(k,1703) - lu(k,372) * lu(k,1592)
         lu(k,2003) = lu(k,2003) - lu(k,366) * lu(k,1989)
         lu(k,2012) = lu(k,2012) - lu(k,367) * lu(k,1989)
         lu(k,2015) = lu(k,2015) - lu(k,368) * lu(k,1989)
         lu(k,2017) = lu(k,2017) - lu(k,369) * lu(k,1989)
         lu(k,2018) = lu(k,2018) - lu(k,370) * lu(k,1989)
         lu(k,2020) = lu(k,2020) - lu(k,371) * lu(k,1989)
         lu(k,2024) = lu(k,2024) - lu(k,372) * lu(k,1989)
                                                                        
      end do
                                                                        
      end subroutine lu_fac07
                                                                        
      subroutine lu_fac08( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,373) = 1._r8 / lu(k,373)
         lu(k,374) = lu(k,374) * lu(k,373)
         lu(k,375) = lu(k,375) * lu(k,373)
         lu(k,376) = lu(k,376) * lu(k,373)
         lu(k,377) = lu(k,377) * lu(k,373)
         lu(k,378) = lu(k,378) * lu(k,373)
         lu(k,1188) = - lu(k,374) * lu(k,1184)
         lu(k,1190) = - lu(k,375) * lu(k,1184)
         lu(k,1198) = - lu(k,376) * lu(k,1184)
         lu(k,1200) = - lu(k,377) * lu(k,1184)
         lu(k,1203) = lu(k,1203) - lu(k,378) * lu(k,1184)
         lu(k,1650) = lu(k,1650) - lu(k,374) * lu(k,1593)
         lu(k,1667) = lu(k,1667) - lu(k,375) * lu(k,1593)
         lu(k,1687) = lu(k,1687) - lu(k,376) * lu(k,1593)
         lu(k,1691) = lu(k,1691) - lu(k,377) * lu(k,1593)
         lu(k,1694) = lu(k,1694) - lu(k,378) * lu(k,1593)
         lu(k,2035) = - lu(k,374) * lu(k,2027)
         lu(k,2044) = lu(k,2044) - lu(k,375) * lu(k,2027)
         lu(k,2060) = - lu(k,376) * lu(k,2027)
         lu(k,2064) = lu(k,2064) - lu(k,377) * lu(k,2027)
         lu(k,2067) = lu(k,2067) - lu(k,378) * lu(k,2027)
                                                                        
         lu(k,379) = 1._r8 / lu(k,379)
         lu(k,380) = lu(k,380) * lu(k,379)
         lu(k,381) = lu(k,381) * lu(k,379)
         lu(k,382) = lu(k,382) * lu(k,379)
         lu(k,383) = lu(k,383) * lu(k,379)
         lu(k,384) = lu(k,384) * lu(k,379)
         lu(k,950) = - lu(k,380) * lu(k,942)
         lu(k,954) = lu(k,954) - lu(k,381) * lu(k,942)
         lu(k,956) = - lu(k,382) * lu(k,942)
         lu(k,957) = lu(k,957) - lu(k,383) * lu(k,942)
         lu(k,963) = lu(k,963) - lu(k,384) * lu(k,942)
         lu(k,998) = - lu(k,380) * lu(k,991)
         lu(k,1003) = lu(k,1003) - lu(k,381) * lu(k,991)
         lu(k,1006) = - lu(k,382) * lu(k,991)
         lu(k,1007) = lu(k,1007) - lu(k,383) * lu(k,991)
         lu(k,1013) = lu(k,1013) - lu(k,384) * lu(k,991)
         lu(k,2098) = - lu(k,380) * lu(k,2083)
         lu(k,2105) = lu(k,2105) - lu(k,381) * lu(k,2083)
         lu(k,2111) = lu(k,2111) - lu(k,382) * lu(k,2083)
         lu(k,2118) = lu(k,2118) - lu(k,383) * lu(k,2083)
         lu(k,2128) = lu(k,2128) - lu(k,384) * lu(k,2083)
                                                                        
         lu(k,385) = 1._r8 / lu(k,385)
         lu(k,386) = lu(k,386) * lu(k,385)
         lu(k,387) = lu(k,387) * lu(k,385)
         lu(k,388) = lu(k,388) * lu(k,385)
         lu(k,389) = lu(k,389) * lu(k,385)
         lu(k,390) = lu(k,390) * lu(k,385)
         lu(k,1688) = lu(k,1688) - lu(k,386) * lu(k,1594)
         lu(k,1689) = lu(k,1689) - lu(k,387) * lu(k,1594)
         lu(k,1691) = lu(k,1691) - lu(k,388) * lu(k,1594)
         lu(k,1697) = lu(k,1697) - lu(k,389) * lu(k,1594)
         lu(k,1703) = lu(k,1703) - lu(k,390) * lu(k,1594)
         lu(k,1944) = lu(k,1944) - lu(k,386) * lu(k,1874)
         lu(k,1945) = lu(k,1945) - lu(k,387) * lu(k,1874)
         lu(k,1947) = lu(k,1947) - lu(k,388) * lu(k,1874)
         lu(k,1953) = lu(k,1953) - lu(k,389) * lu(k,1874)
         lu(k,1959) = lu(k,1959) - lu(k,390) * lu(k,1874)
         lu(k,2061) = - lu(k,386) * lu(k,2028)
         lu(k,2062) = lu(k,2062) - lu(k,387) * lu(k,2028)
         lu(k,2064) = lu(k,2064) - lu(k,388) * lu(k,2028)
         lu(k,2070) = lu(k,2070) - lu(k,389) * lu(k,2028)
         lu(k,2076) = - lu(k,390) * lu(k,2028)
                                                                        
         lu(k,391) = 1._r8 / lu(k,391)
         lu(k,392) = lu(k,392) * lu(k,391)
         lu(k,393) = lu(k,393) * lu(k,391)
         lu(k,394) = lu(k,394) * lu(k,391)
         lu(k,395) = lu(k,395) * lu(k,391)
         lu(k,396) = lu(k,396) * lu(k,391)
         lu(k,1044) = lu(k,1044) - lu(k,392) * lu(k,1041)
         lu(k,1045) = lu(k,1045) - lu(k,393) * lu(k,1041)
         lu(k,1049) = - lu(k,394) * lu(k,1041)
         lu(k,1051) = - lu(k,395) * lu(k,1041)
         lu(k,1056) = lu(k,1056) - lu(k,396) * lu(k,1041)
         lu(k,1659) = lu(k,1659) - lu(k,392) * lu(k,1595)
         lu(k,1663) = lu(k,1663) - lu(k,393) * lu(k,1595)
         lu(k,1687) = lu(k,1687) - lu(k,394) * lu(k,1595)
         lu(k,1691) = lu(k,1691) - lu(k,395) * lu(k,1595)
         lu(k,1700) = lu(k,1700) - lu(k,396) * lu(k,1595)
         lu(k,1920) = - lu(k,392) * lu(k,1875)
         lu(k,1922) = lu(k,1922) - lu(k,393) * lu(k,1875)
         lu(k,1943) = - lu(k,394) * lu(k,1875)
         lu(k,1947) = lu(k,1947) - lu(k,395) * lu(k,1875)
         lu(k,1956) = lu(k,1956) - lu(k,396) * lu(k,1875)
                                                                        
         lu(k,397) = 1._r8 / lu(k,397)
         lu(k,398) = lu(k,398) * lu(k,397)
         lu(k,399) = lu(k,399) * lu(k,397)
         lu(k,400) = lu(k,400) * lu(k,397)
         lu(k,401) = lu(k,401) * lu(k,397)
         lu(k,402) = lu(k,402) * lu(k,397)
         lu(k,1186) = - lu(k,398) * lu(k,1185)
         lu(k,1188) = lu(k,1188) - lu(k,399) * lu(k,1185)
         lu(k,1200) = lu(k,1200) - lu(k,400) * lu(k,1185)
         lu(k,1203) = lu(k,1203) - lu(k,401) * lu(k,1185)
         lu(k,1205) = lu(k,1205) - lu(k,402) * lu(k,1185)
         lu(k,1634) = lu(k,1634) - lu(k,398) * lu(k,1596)
         lu(k,1650) = lu(k,1650) - lu(k,399) * lu(k,1596)
         lu(k,1691) = lu(k,1691) - lu(k,400) * lu(k,1596)
         lu(k,1694) = lu(k,1694) - lu(k,401) * lu(k,1596)
         lu(k,1700) = lu(k,1700) - lu(k,402) * lu(k,1596)
         lu(k,1900) = lu(k,1900) - lu(k,398) * lu(k,1876)
         lu(k,1914) = - lu(k,399) * lu(k,1876)
         lu(k,1947) = lu(k,1947) - lu(k,400) * lu(k,1876)
         lu(k,1950) = lu(k,1950) - lu(k,401) * lu(k,1876)
         lu(k,1956) = lu(k,1956) - lu(k,402) * lu(k,1876)
                                                                        
         lu(k,403) = 1._r8 / lu(k,403)
         lu(k,404) = lu(k,404) * lu(k,403)
         lu(k,405) = lu(k,405) * lu(k,403)
         lu(k,406) = lu(k,406) * lu(k,403)
         lu(k,407) = lu(k,407) * lu(k,403)
         lu(k,408) = lu(k,408) * lu(k,403)
         lu(k,1691) = lu(k,1691) - lu(k,404) * lu(k,1597)
         lu(k,1692) = lu(k,1692) - lu(k,405) * lu(k,1597)
         lu(k,1694) = lu(k,1694) - lu(k,406) * lu(k,1597)
         lu(k,1700) = lu(k,1700) - lu(k,407) * lu(k,1597)
         lu(k,1703) = lu(k,1703) - lu(k,408) * lu(k,1597)
         lu(k,1947) = lu(k,1947) - lu(k,404) * lu(k,1877)
         lu(k,1948) = lu(k,1948) - lu(k,405) * lu(k,1877)
         lu(k,1950) = lu(k,1950) - lu(k,406) * lu(k,1877)
         lu(k,1956) = lu(k,1956) - lu(k,407) * lu(k,1877)
         lu(k,1959) = lu(k,1959) - lu(k,408) * lu(k,1877)
         lu(k,2192) = lu(k,2192) - lu(k,404) * lu(k,2164)
         lu(k,2193) = lu(k,2193) - lu(k,405) * lu(k,2164)
         lu(k,2195) = lu(k,2195) - lu(k,406) * lu(k,2164)
         lu(k,2201) = lu(k,2201) - lu(k,407) * lu(k,2164)
         lu(k,2204) = - lu(k,408) * lu(k,2164)
                                                                        
      end do
                                                                        
      end subroutine lu_fac08
                                                                        
      subroutine lu_fac09( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,409) = 1._r8 / lu(k,409)
         lu(k,410) = lu(k,410) * lu(k,409)
         lu(k,411) = lu(k,411) * lu(k,409)
         lu(k,412) = lu(k,412) * lu(k,409)
         lu(k,413) = lu(k,413) * lu(k,409)
         lu(k,414) = lu(k,414) * lu(k,409)
         lu(k,1129) = lu(k,1129) - lu(k,410) * lu(k,1127)
         lu(k,1132) = lu(k,1132) - lu(k,411) * lu(k,1127)
         lu(k,1133) = lu(k,1133) - lu(k,412) * lu(k,1127)
         lu(k,1134) = lu(k,1134) - lu(k,413) * lu(k,1127)
         lu(k,1139) = - lu(k,414) * lu(k,1127)
         lu(k,1670) = lu(k,1670) - lu(k,410) * lu(k,1598)
         lu(k,1683) = lu(k,1683) - lu(k,411) * lu(k,1598)
         lu(k,1689) = lu(k,1689) - lu(k,412) * lu(k,1598)
         lu(k,1691) = lu(k,1691) - lu(k,413) * lu(k,1598)
         lu(k,1703) = lu(k,1703) - lu(k,414) * lu(k,1598)
         lu(k,1928) = lu(k,1928) - lu(k,410) * lu(k,1878)
         lu(k,1939) = lu(k,1939) - lu(k,411) * lu(k,1878)
         lu(k,1945) = lu(k,1945) - lu(k,412) * lu(k,1878)
         lu(k,1947) = lu(k,1947) - lu(k,413) * lu(k,1878)
         lu(k,1959) = lu(k,1959) - lu(k,414) * lu(k,1878)
                                                                        
         lu(k,415) = 1._r8 / lu(k,415)
         lu(k,416) = lu(k,416) * lu(k,415)
         lu(k,417) = lu(k,417) * lu(k,415)
         lu(k,418) = lu(k,418) * lu(k,415)
         lu(k,419) = lu(k,419) * lu(k,415)
         lu(k,420) = lu(k,420) * lu(k,415)
         lu(k,791) = lu(k,791) - lu(k,416) * lu(k,790)
         lu(k,792) = lu(k,792) - lu(k,417) * lu(k,790)
         lu(k,795) = - lu(k,418) * lu(k,790)
         lu(k,797) = lu(k,797) - lu(k,419) * lu(k,790)
         lu(k,800) = - lu(k,420) * lu(k,790)
         lu(k,1643) = lu(k,1643) - lu(k,416) * lu(k,1599)
         lu(k,1656) = lu(k,1656) - lu(k,417) * lu(k,1599)
         lu(k,1691) = lu(k,1691) - lu(k,418) * lu(k,1599)
         lu(k,1694) = lu(k,1694) - lu(k,419) * lu(k,1599)
         lu(k,1703) = lu(k,1703) - lu(k,420) * lu(k,1599)
         lu(k,1909) = lu(k,1909) - lu(k,416) * lu(k,1879)
         lu(k,1918) = - lu(k,417) * lu(k,1879)
         lu(k,1947) = lu(k,1947) - lu(k,418) * lu(k,1879)
         lu(k,1950) = lu(k,1950) - lu(k,419) * lu(k,1879)
         lu(k,1959) = lu(k,1959) - lu(k,420) * lu(k,1879)
                                                                        
         lu(k,421) = 1._r8 / lu(k,421)
         lu(k,422) = lu(k,422) * lu(k,421)
         lu(k,423) = lu(k,423) * lu(k,421)
         lu(k,424) = lu(k,424) * lu(k,421)
         lu(k,425) = lu(k,425) * lu(k,421)
         lu(k,426) = lu(k,426) * lu(k,421)
         lu(k,482) = lu(k,482) - lu(k,422) * lu(k,481)
         lu(k,483) = lu(k,483) - lu(k,423) * lu(k,481)
         lu(k,485) = lu(k,485) - lu(k,424) * lu(k,481)
         lu(k,486) = - lu(k,425) * lu(k,481)
         lu(k,488) = lu(k,488) - lu(k,426) * lu(k,481)
         lu(k,1605) = lu(k,1605) - lu(k,422) * lu(k,1600)
         lu(k,1609) = lu(k,1609) - lu(k,423) * lu(k,1600)
         lu(k,1660) = lu(k,1660) - lu(k,424) * lu(k,1600)
         lu(k,1691) = lu(k,1691) - lu(k,425) * lu(k,1600)
         lu(k,1694) = lu(k,1694) - lu(k,426) * lu(k,1600)
         lu(k,1885) = - lu(k,422) * lu(k,1880)
         lu(k,1888) = lu(k,1888) - lu(k,423) * lu(k,1880)
         lu(k,1921) = lu(k,1921) - lu(k,424) * lu(k,1880)
         lu(k,1947) = lu(k,1947) - lu(k,425) * lu(k,1880)
         lu(k,1950) = lu(k,1950) - lu(k,426) * lu(k,1880)
                                                                        
         lu(k,427) = 1._r8 / lu(k,427)
         lu(k,428) = lu(k,428) * lu(k,427)
         lu(k,429) = lu(k,429) * lu(k,427)
         lu(k,430) = lu(k,430) * lu(k,427)
         lu(k,431) = lu(k,431) * lu(k,427)
         lu(k,432) = lu(k,432) * lu(k,427)
         lu(k,499) = lu(k,499) - lu(k,428) * lu(k,498)
         lu(k,500) = lu(k,500) - lu(k,429) * lu(k,498)
         lu(k,501) = lu(k,501) - lu(k,430) * lu(k,498)
         lu(k,502) = - lu(k,431) * lu(k,498)
         lu(k,504) = lu(k,504) - lu(k,432) * lu(k,498)
         lu(k,1611) = lu(k,1611) - lu(k,428) * lu(k,1601)
         lu(k,1660) = lu(k,1660) - lu(k,429) * lu(k,1601)
         lu(k,1676) = lu(k,1676) - lu(k,430) * lu(k,1601)
         lu(k,1691) = lu(k,1691) - lu(k,431) * lu(k,1601)
         lu(k,1694) = lu(k,1694) - lu(k,432) * lu(k,1601)
         lu(k,1889) = lu(k,1889) - lu(k,428) * lu(k,1881)
         lu(k,1921) = lu(k,1921) - lu(k,429) * lu(k,1881)
         lu(k,1933) = lu(k,1933) - lu(k,430) * lu(k,1881)
         lu(k,1947) = lu(k,1947) - lu(k,431) * lu(k,1881)
         lu(k,1950) = lu(k,1950) - lu(k,432) * lu(k,1881)
                                                                        
         lu(k,435) = 1._r8 / lu(k,435)
         lu(k,436) = lu(k,436) * lu(k,435)
         lu(k,437) = lu(k,437) * lu(k,435)
         lu(k,438) = lu(k,438) * lu(k,435)
         lu(k,439) = lu(k,439) * lu(k,435)
         lu(k,440) = lu(k,440) * lu(k,435)
         lu(k,1614) = lu(k,1614) - lu(k,436) * lu(k,1602)
         lu(k,1691) = lu(k,1691) - lu(k,437) * lu(k,1602)
         lu(k,1693) = lu(k,1693) - lu(k,438) * lu(k,1602)
         lu(k,1694) = lu(k,1694) - lu(k,439) * lu(k,1602)
         lu(k,1700) = lu(k,1700) - lu(k,440) * lu(k,1602)
         lu(k,1788) = lu(k,1788) - lu(k,436) * lu(k,1779)
         lu(k,1840) = lu(k,1840) - lu(k,437) * lu(k,1779)
         lu(k,1842) = lu(k,1842) - lu(k,438) * lu(k,1779)
         lu(k,1843) = lu(k,1843) - lu(k,439) * lu(k,1779)
         lu(k,1849) = lu(k,1849) - lu(k,440) * lu(k,1779)
         lu(k,1890) = lu(k,1890) - lu(k,436) * lu(k,1882)
         lu(k,1947) = lu(k,1947) - lu(k,437) * lu(k,1882)
         lu(k,1949) = lu(k,1949) - lu(k,438) * lu(k,1882)
         lu(k,1950) = lu(k,1950) - lu(k,439) * lu(k,1882)
         lu(k,1956) = lu(k,1956) - lu(k,440) * lu(k,1882)
                                                                        
         lu(k,441) = 1._r8 / lu(k,441)
         lu(k,442) = lu(k,442) * lu(k,441)
         lu(k,443) = lu(k,443) * lu(k,441)
         lu(k,444) = lu(k,444) * lu(k,441)
         lu(k,445) = lu(k,445) * lu(k,441)
         lu(k,446) = lu(k,446) * lu(k,441)
         lu(k,1478) = - lu(k,442) * lu(k,1476)
         lu(k,1485) = lu(k,1485) - lu(k,443) * lu(k,1476)
         lu(k,1489) = - lu(k,444) * lu(k,1476)
         lu(k,1490) = lu(k,1490) - lu(k,445) * lu(k,1476)
         lu(k,1495) = - lu(k,446) * lu(k,1476)
         lu(k,1811) = lu(k,1811) - lu(k,442) * lu(k,1780)
         lu(k,1838) = lu(k,1838) - lu(k,443) * lu(k,1780)
         lu(k,1842) = lu(k,1842) - lu(k,444) * lu(k,1780)
         lu(k,1843) = lu(k,1843) - lu(k,445) * lu(k,1780)
         lu(k,1849) = lu(k,1849) - lu(k,446) * lu(k,1780)
         lu(k,1919) = lu(k,1919) - lu(k,442) * lu(k,1883)
         lu(k,1945) = lu(k,1945) - lu(k,443) * lu(k,1883)
         lu(k,1949) = lu(k,1949) - lu(k,444) * lu(k,1883)
         lu(k,1950) = lu(k,1950) - lu(k,445) * lu(k,1883)
         lu(k,1956) = lu(k,1956) - lu(k,446) * lu(k,1883)
                                                                        
         lu(k,447) = 1._r8 / lu(k,447)
         lu(k,448) = lu(k,448) * lu(k,447)
         lu(k,449) = lu(k,449) * lu(k,447)
         lu(k,450) = lu(k,450) * lu(k,447)
         lu(k,451) = lu(k,451) * lu(k,447)
         lu(k,452) = lu(k,452) * lu(k,447)
         lu(k,1521) = lu(k,1521) - lu(k,448) * lu(k,1515)
         lu(k,1526) = lu(k,1526) - lu(k,449) * lu(k,1515)
         lu(k,1527) = lu(k,1527) - lu(k,450) * lu(k,1515)
         lu(k,1532) = lu(k,1532) - lu(k,451) * lu(k,1515)
         lu(k,1535) = lu(k,1535) - lu(k,452) * lu(k,1515)
         lu(k,1685) = lu(k,1685) - lu(k,448) * lu(k,1603)
         lu(k,1690) = lu(k,1690) - lu(k,449) * lu(k,1603)
         lu(k,1691) = lu(k,1691) - lu(k,450) * lu(k,1603)
         lu(k,1696) = lu(k,1696) - lu(k,451) * lu(k,1603)
         lu(k,1699) = lu(k,1699) - lu(k,452) * lu(k,1603)
         lu(k,2006) = lu(k,2006) - lu(k,448) * lu(k,1990)
         lu(k,2011) = - lu(k,449) * lu(k,1990)
         lu(k,2012) = lu(k,2012) - lu(k,450) * lu(k,1990)
         lu(k,2017) = lu(k,2017) - lu(k,451) * lu(k,1990)
         lu(k,2020) = lu(k,2020) - lu(k,452) * lu(k,1990)
                                                                        
      end do
                                                                        
      end subroutine lu_fac09
                                                                        
      subroutine lu_fac10( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,454) = 1._r8 / lu(k,454)
         lu(k,455) = lu(k,455) * lu(k,454)
         lu(k,456) = lu(k,456) * lu(k,454)
         lu(k,457) = lu(k,457) * lu(k,454)
         lu(k,458) = lu(k,458) * lu(k,454)
         lu(k,459) = lu(k,459) * lu(k,454)
         lu(k,1660) = lu(k,1660) - lu(k,455) * lu(k,1604)
         lu(k,1691) = lu(k,1691) - lu(k,456) * lu(k,1604)
         lu(k,1693) = lu(k,1693) - lu(k,457) * lu(k,1604)
         lu(k,1694) = lu(k,1694) - lu(k,458) * lu(k,1604)
         lu(k,1700) = lu(k,1700) - lu(k,459) * lu(k,1604)
         lu(k,1813) = lu(k,1813) - lu(k,455) * lu(k,1781)
         lu(k,1840) = lu(k,1840) - lu(k,456) * lu(k,1781)
         lu(k,1842) = lu(k,1842) - lu(k,457) * lu(k,1781)
         lu(k,1843) = lu(k,1843) - lu(k,458) * lu(k,1781)
         lu(k,1849) = lu(k,1849) - lu(k,459) * lu(k,1781)
         lu(k,1921) = lu(k,1921) - lu(k,455) * lu(k,1884)
         lu(k,1947) = lu(k,1947) - lu(k,456) * lu(k,1884)
         lu(k,1949) = lu(k,1949) - lu(k,457) * lu(k,1884)
         lu(k,1950) = lu(k,1950) - lu(k,458) * lu(k,1884)
         lu(k,1956) = lu(k,1956) - lu(k,459) * lu(k,1884)
                                                                        
         lu(k,460) = 1._r8 / lu(k,460)
         lu(k,461) = lu(k,461) * lu(k,460)
         lu(k,462) = lu(k,462) * lu(k,460)
         lu(k,484) = - lu(k,461) * lu(k,482)
         lu(k,488) = lu(k,488) - lu(k,462) * lu(k,482)
         lu(k,656) = - lu(k,461) * lu(k,653)
         lu(k,664) = lu(k,664) - lu(k,462) * lu(k,653)
         lu(k,702) = - lu(k,461) * lu(k,699)
         lu(k,711) = lu(k,711) - lu(k,462) * lu(k,699)
         lu(k,731) = - lu(k,461) * lu(k,728)
         lu(k,740) = lu(k,740) - lu(k,462) * lu(k,728)
         lu(k,747) = - lu(k,461) * lu(k,744)
         lu(k,757) = lu(k,757) - lu(k,462) * lu(k,744)
         lu(k,1637) = - lu(k,461) * lu(k,1605)
         lu(k,1694) = lu(k,1694) - lu(k,462) * lu(k,1605)
         lu(k,1798) = lu(k,1798) - lu(k,461) * lu(k,1782)
         lu(k,1843) = lu(k,1843) - lu(k,462) * lu(k,1782)
         lu(k,1903) = lu(k,1903) - lu(k,461) * lu(k,1885)
         lu(k,1950) = lu(k,1950) - lu(k,462) * lu(k,1885)
                                                                        
         lu(k,463) = 1._r8 / lu(k,463)
         lu(k,464) = lu(k,464) * lu(k,463)
         lu(k,465) = lu(k,465) * lu(k,463)
         lu(k,466) = lu(k,466) * lu(k,463)
         lu(k,625) = lu(k,625) - lu(k,464) * lu(k,624)
         lu(k,629) = - lu(k,465) * lu(k,624)
         lu(k,630) = lu(k,630) - lu(k,466) * lu(k,624)
         lu(k,1627) = lu(k,1627) - lu(k,464) * lu(k,1606)
         lu(k,1698) = lu(k,1698) - lu(k,465) * lu(k,1606)
         lu(k,1700) = lu(k,1700) - lu(k,466) * lu(k,1606)
         lu(k,1794) = lu(k,1794) - lu(k,464) * lu(k,1783)
         lu(k,1847) = lu(k,1847) - lu(k,465) * lu(k,1783)
         lu(k,1849) = lu(k,1849) - lu(k,466) * lu(k,1783)
         lu(k,1894) = lu(k,1894) - lu(k,464) * lu(k,1886)
         lu(k,1954) = lu(k,1954) - lu(k,465) * lu(k,1886)
         lu(k,1956) = lu(k,1956) - lu(k,466) * lu(k,1886)
         lu(k,2087) = lu(k,2087) - lu(k,464) * lu(k,2084)
         lu(k,2132) = lu(k,2132) - lu(k,465) * lu(k,2084)
         lu(k,2134) = lu(k,2134) - lu(k,466) * lu(k,2084)
         lu(k,2171) = - lu(k,464) * lu(k,2165)
         lu(k,2199) = lu(k,2199) - lu(k,465) * lu(k,2165)
         lu(k,2201) = lu(k,2201) - lu(k,466) * lu(k,2165)
                                                                        
         lu(k,467) = 1._r8 / lu(k,467)
         lu(k,468) = lu(k,468) * lu(k,467)
         lu(k,469) = lu(k,469) * lu(k,467)
         lu(k,470) = lu(k,470) * lu(k,467)
         lu(k,471) = lu(k,471) * lu(k,467)
         lu(k,472) = lu(k,472) * lu(k,467)
         lu(k,473) = lu(k,473) * lu(k,467)
         lu(k,1521) = lu(k,1521) - lu(k,468) * lu(k,1516)
         lu(k,1526) = lu(k,1526) - lu(k,469) * lu(k,1516)
         lu(k,1527) = lu(k,1527) - lu(k,470) * lu(k,1516)
         lu(k,1532) = lu(k,1532) - lu(k,471) * lu(k,1516)
         lu(k,1535) = lu(k,1535) - lu(k,472) * lu(k,1516)
         lu(k,1539) = lu(k,1539) - lu(k,473) * lu(k,1516)
         lu(k,1685) = lu(k,1685) - lu(k,468) * lu(k,1607)
         lu(k,1690) = lu(k,1690) - lu(k,469) * lu(k,1607)
         lu(k,1691) = lu(k,1691) - lu(k,470) * lu(k,1607)
         lu(k,1696) = lu(k,1696) - lu(k,471) * lu(k,1607)
         lu(k,1699) = lu(k,1699) - lu(k,472) * lu(k,1607)
         lu(k,1703) = lu(k,1703) - lu(k,473) * lu(k,1607)
         lu(k,2006) = lu(k,2006) - lu(k,468) * lu(k,1991)
         lu(k,2011) = lu(k,2011) - lu(k,469) * lu(k,1991)
         lu(k,2012) = lu(k,2012) - lu(k,470) * lu(k,1991)
         lu(k,2017) = lu(k,2017) - lu(k,471) * lu(k,1991)
         lu(k,2020) = lu(k,2020) - lu(k,472) * lu(k,1991)
         lu(k,2024) = lu(k,2024) - lu(k,473) * lu(k,1991)
                                                                        
         lu(k,474) = 1._r8 / lu(k,474)
         lu(k,475) = lu(k,475) * lu(k,474)
         lu(k,476) = lu(k,476) * lu(k,474)
         lu(k,477) = lu(k,477) * lu(k,474)
         lu(k,478) = lu(k,478) * lu(k,474)
         lu(k,479) = lu(k,479) * lu(k,474)
         lu(k,480) = lu(k,480) * lu(k,474)
         lu(k,912) = lu(k,912) - lu(k,475) * lu(k,909)
         lu(k,913) = lu(k,913) - lu(k,476) * lu(k,909)
         lu(k,914) = lu(k,914) - lu(k,477) * lu(k,909)
         lu(k,916) = lu(k,916) - lu(k,478) * lu(k,909)
         lu(k,917) = - lu(k,479) * lu(k,909)
         lu(k,919) = lu(k,919) - lu(k,480) * lu(k,909)
         lu(k,1655) = lu(k,1655) - lu(k,475) * lu(k,1608)
         lu(k,1656) = lu(k,1656) - lu(k,476) * lu(k,1608)
         lu(k,1659) = lu(k,1659) - lu(k,477) * lu(k,1608)
         lu(k,1689) = lu(k,1689) - lu(k,478) * lu(k,1608)
         lu(k,1691) = lu(k,1691) - lu(k,479) * lu(k,1608)
         lu(k,1694) = lu(k,1694) - lu(k,480) * lu(k,1608)
         lu(k,1917) = lu(k,1917) - lu(k,475) * lu(k,1887)
         lu(k,1918) = lu(k,1918) - lu(k,476) * lu(k,1887)
         lu(k,1920) = lu(k,1920) - lu(k,477) * lu(k,1887)
         lu(k,1945) = lu(k,1945) - lu(k,478) * lu(k,1887)
         lu(k,1947) = lu(k,1947) - lu(k,479) * lu(k,1887)
         lu(k,1950) = lu(k,1950) - lu(k,480) * lu(k,1887)
                                                                        
         lu(k,483) = 1._r8 / lu(k,483)
         lu(k,484) = lu(k,484) * lu(k,483)
         lu(k,485) = lu(k,485) * lu(k,483)
         lu(k,486) = lu(k,486) * lu(k,483)
         lu(k,487) = lu(k,487) * lu(k,483)
         lu(k,488) = lu(k,488) * lu(k,483)
         lu(k,489) = lu(k,489) * lu(k,483)
         lu(k,1637) = lu(k,1637) - lu(k,484) * lu(k,1609)
         lu(k,1660) = lu(k,1660) - lu(k,485) * lu(k,1609)
         lu(k,1691) = lu(k,1691) - lu(k,486) * lu(k,1609)
         lu(k,1693) = lu(k,1693) - lu(k,487) * lu(k,1609)
         lu(k,1694) = lu(k,1694) - lu(k,488) * lu(k,1609)
         lu(k,1700) = lu(k,1700) - lu(k,489) * lu(k,1609)
         lu(k,1798) = lu(k,1798) - lu(k,484) * lu(k,1784)
         lu(k,1813) = lu(k,1813) - lu(k,485) * lu(k,1784)
         lu(k,1840) = lu(k,1840) - lu(k,486) * lu(k,1784)
         lu(k,1842) = lu(k,1842) - lu(k,487) * lu(k,1784)
         lu(k,1843) = lu(k,1843) - lu(k,488) * lu(k,1784)
         lu(k,1849) = lu(k,1849) - lu(k,489) * lu(k,1784)
         lu(k,1903) = lu(k,1903) - lu(k,484) * lu(k,1888)
         lu(k,1921) = lu(k,1921) - lu(k,485) * lu(k,1888)
         lu(k,1947) = lu(k,1947) - lu(k,486) * lu(k,1888)
         lu(k,1949) = lu(k,1949) - lu(k,487) * lu(k,1888)
         lu(k,1950) = lu(k,1950) - lu(k,488) * lu(k,1888)
         lu(k,1956) = lu(k,1956) - lu(k,489) * lu(k,1888)
                                                                        
         lu(k,491) = 1._r8 / lu(k,491)
         lu(k,492) = lu(k,492) * lu(k,491)
         lu(k,493) = lu(k,493) * lu(k,491)
         lu(k,494) = lu(k,494) * lu(k,491)
         lu(k,495) = lu(k,495) * lu(k,491)
         lu(k,496) = lu(k,496) * lu(k,491)
         lu(k,497) = lu(k,497) * lu(k,491)
         lu(k,1688) = lu(k,1688) - lu(k,492) * lu(k,1610)
         lu(k,1690) = lu(k,1690) - lu(k,493) * lu(k,1610)
         lu(k,1691) = lu(k,1691) - lu(k,494) * lu(k,1610)
         lu(k,1693) = lu(k,1693) - lu(k,495) * lu(k,1610)
         lu(k,1700) = lu(k,1700) - lu(k,496) * lu(k,1610)
         lu(k,1702) = lu(k,1702) - lu(k,497) * lu(k,1610)
         lu(k,1837) = - lu(k,492) * lu(k,1785)
         lu(k,1839) = - lu(k,493) * lu(k,1785)
         lu(k,1840) = lu(k,1840) - lu(k,494) * lu(k,1785)
         lu(k,1842) = lu(k,1842) - lu(k,495) * lu(k,1785)
         lu(k,1849) = lu(k,1849) - lu(k,496) * lu(k,1785)
         lu(k,1851) = lu(k,1851) - lu(k,497) * lu(k,1785)
         lu(k,2189) = - lu(k,492) * lu(k,2166)
         lu(k,2191) = lu(k,2191) - lu(k,493) * lu(k,2166)
         lu(k,2192) = lu(k,2192) - lu(k,494) * lu(k,2166)
         lu(k,2194) = lu(k,2194) - lu(k,495) * lu(k,2166)
         lu(k,2201) = lu(k,2201) - lu(k,496) * lu(k,2166)
         lu(k,2203) = lu(k,2203) - lu(k,497) * lu(k,2166)
                                                                        
      end do
                                                                        
      end subroutine lu_fac10
                                                                        
      subroutine lu_fac11( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,499) = 1._r8 / lu(k,499)
         lu(k,500) = lu(k,500) * lu(k,499)
         lu(k,501) = lu(k,501) * lu(k,499)
         lu(k,502) = lu(k,502) * lu(k,499)
         lu(k,503) = lu(k,503) * lu(k,499)
         lu(k,504) = lu(k,504) * lu(k,499)
         lu(k,505) = lu(k,505) * lu(k,499)
         lu(k,1660) = lu(k,1660) - lu(k,500) * lu(k,1611)
         lu(k,1676) = lu(k,1676) - lu(k,501) * lu(k,1611)
         lu(k,1691) = lu(k,1691) - lu(k,502) * lu(k,1611)
         lu(k,1693) = lu(k,1693) - lu(k,503) * lu(k,1611)
         lu(k,1694) = lu(k,1694) - lu(k,504) * lu(k,1611)
         lu(k,1700) = lu(k,1700) - lu(k,505) * lu(k,1611)
         lu(k,1813) = lu(k,1813) - lu(k,500) * lu(k,1786)
         lu(k,1827) = lu(k,1827) - lu(k,501) * lu(k,1786)
         lu(k,1840) = lu(k,1840) - lu(k,502) * lu(k,1786)
         lu(k,1842) = lu(k,1842) - lu(k,503) * lu(k,1786)
         lu(k,1843) = lu(k,1843) - lu(k,504) * lu(k,1786)
         lu(k,1849) = lu(k,1849) - lu(k,505) * lu(k,1786)
         lu(k,1921) = lu(k,1921) - lu(k,500) * lu(k,1889)
         lu(k,1933) = lu(k,1933) - lu(k,501) * lu(k,1889)
         lu(k,1947) = lu(k,1947) - lu(k,502) * lu(k,1889)
         lu(k,1949) = lu(k,1949) - lu(k,503) * lu(k,1889)
         lu(k,1950) = lu(k,1950) - lu(k,504) * lu(k,1889)
         lu(k,1956) = lu(k,1956) - lu(k,505) * lu(k,1889)
                                                                        
         lu(k,506) = 1._r8 / lu(k,506)
         lu(k,507) = lu(k,507) * lu(k,506)
         lu(k,508) = lu(k,508) * lu(k,506)
         lu(k,509) = lu(k,509) * lu(k,506)
         lu(k,510) = lu(k,510) * lu(k,506)
         lu(k,633) = - lu(k,507) * lu(k,631)
         lu(k,634) = - lu(k,508) * lu(k,631)
         lu(k,638) = - lu(k,509) * lu(k,631)
         lu(k,640) = lu(k,640) - lu(k,510) * lu(k,631)
         lu(k,667) = - lu(k,507) * lu(k,665)
         lu(k,668) = - lu(k,508) * lu(k,665)
         lu(k,671) = - lu(k,509) * lu(k,665)
         lu(k,673) = lu(k,673) - lu(k,510) * lu(k,665)
         lu(k,897) = - lu(k,507) * lu(k,894)
         lu(k,898) = - lu(k,508) * lu(k,894)
         lu(k,902) = - lu(k,509) * lu(k,894)
         lu(k,904) = - lu(k,510) * lu(k,894)
         lu(k,1632) = lu(k,1632) - lu(k,507) * lu(k,1612)
         lu(k,1646) = lu(k,1646) - lu(k,508) * lu(k,1612)
         lu(k,1683) = lu(k,1683) - lu(k,509) * lu(k,1612)
         lu(k,1691) = lu(k,1691) - lu(k,510) * lu(k,1612)
         lu(k,1796) = lu(k,1796) - lu(k,507) * lu(k,1787)
         lu(k,1805) = lu(k,1805) - lu(k,508) * lu(k,1787)
         lu(k,1833) = lu(k,1833) - lu(k,509) * lu(k,1787)
         lu(k,1840) = lu(k,1840) - lu(k,510) * lu(k,1787)
                                                                        
         lu(k,511) = 1._r8 / lu(k,511)
         lu(k,512) = lu(k,512) * lu(k,511)
         lu(k,513) = lu(k,513) * lu(k,511)
         lu(k,514) = lu(k,514) * lu(k,511)
         lu(k,515) = lu(k,515) * lu(k,511)
         lu(k,516) = lu(k,516) * lu(k,511)
         lu(k,1232) = lu(k,1232) - lu(k,512) * lu(k,1230)
         lu(k,1234) = lu(k,1234) - lu(k,513) * lu(k,1230)
         lu(k,1235) = lu(k,1235) - lu(k,514) * lu(k,1230)
         lu(k,1240) = lu(k,1240) - lu(k,515) * lu(k,1230)
         lu(k,1243) = lu(k,1243) - lu(k,516) * lu(k,1230)
         lu(k,1677) = lu(k,1677) - lu(k,512) * lu(k,1613)
         lu(k,1688) = lu(k,1688) - lu(k,513) * lu(k,1613)
         lu(k,1691) = lu(k,1691) - lu(k,514) * lu(k,1613)
         lu(k,1698) = lu(k,1698) - lu(k,515) * lu(k,1613)
         lu(k,1702) = lu(k,1702) - lu(k,516) * lu(k,1613)
         lu(k,2112) = lu(k,2112) - lu(k,512) * lu(k,2085)
         lu(k,2122) = lu(k,2122) - lu(k,513) * lu(k,2085)
         lu(k,2125) = lu(k,2125) - lu(k,514) * lu(k,2085)
         lu(k,2132) = lu(k,2132) - lu(k,515) * lu(k,2085)
         lu(k,2136) = lu(k,2136) - lu(k,516) * lu(k,2085)
         lu(k,2239) = lu(k,2239) - lu(k,512) * lu(k,2230)
         lu(k,2244) = lu(k,2244) - lu(k,513) * lu(k,2230)
         lu(k,2247) = lu(k,2247) - lu(k,514) * lu(k,2230)
         lu(k,2254) = lu(k,2254) - lu(k,515) * lu(k,2230)
         lu(k,2258) = lu(k,2258) - lu(k,516) * lu(k,2230)
                                                                        
         lu(k,518) = 1._r8 / lu(k,518)
         lu(k,519) = lu(k,519) * lu(k,518)
         lu(k,520) = lu(k,520) * lu(k,518)
         lu(k,521) = lu(k,521) * lu(k,518)
         lu(k,522) = lu(k,522) * lu(k,518)
         lu(k,523) = lu(k,523) * lu(k,518)
         lu(k,1627) = lu(k,1627) - lu(k,519) * lu(k,1614)
         lu(k,1691) = lu(k,1691) - lu(k,520) * lu(k,1614)
         lu(k,1693) = lu(k,1693) - lu(k,521) * lu(k,1614)
         lu(k,1694) = lu(k,1694) - lu(k,522) * lu(k,1614)
         lu(k,1700) = lu(k,1700) - lu(k,523) * lu(k,1614)
         lu(k,1794) = lu(k,1794) - lu(k,519) * lu(k,1788)
         lu(k,1840) = lu(k,1840) - lu(k,520) * lu(k,1788)
         lu(k,1842) = lu(k,1842) - lu(k,521) * lu(k,1788)
         lu(k,1843) = lu(k,1843) - lu(k,522) * lu(k,1788)
         lu(k,1849) = lu(k,1849) - lu(k,523) * lu(k,1788)
         lu(k,1894) = lu(k,1894) - lu(k,519) * lu(k,1890)
         lu(k,1947) = lu(k,1947) - lu(k,520) * lu(k,1890)
         lu(k,1949) = lu(k,1949) - lu(k,521) * lu(k,1890)
         lu(k,1950) = lu(k,1950) - lu(k,522) * lu(k,1890)
         lu(k,1956) = lu(k,1956) - lu(k,523) * lu(k,1890)
         lu(k,2171) = lu(k,2171) - lu(k,519) * lu(k,2167)
         lu(k,2192) = lu(k,2192) - lu(k,520) * lu(k,2167)
         lu(k,2194) = lu(k,2194) - lu(k,521) * lu(k,2167)
         lu(k,2195) = lu(k,2195) - lu(k,522) * lu(k,2167)
         lu(k,2201) = lu(k,2201) - lu(k,523) * lu(k,2167)
                                                                        
         lu(k,524) = 1._r8 / lu(k,524)
         lu(k,525) = lu(k,525) * lu(k,524)
         lu(k,526) = lu(k,526) * lu(k,524)
         lu(k,527) = lu(k,527) * lu(k,524)
         lu(k,528) = lu(k,528) * lu(k,524)
         lu(k,529) = lu(k,529) * lu(k,524)
         lu(k,530) = lu(k,530) * lu(k,524)
         lu(k,531) = lu(k,531) * lu(k,524)
         lu(k,1371) = lu(k,1371) - lu(k,525) * lu(k,1368)
         lu(k,1386) = lu(k,1386) - lu(k,526) * lu(k,1368)
         lu(k,1389) = lu(k,1389) - lu(k,527) * lu(k,1368)
         lu(k,1390) = lu(k,1390) - lu(k,528) * lu(k,1368)
         lu(k,1391) = - lu(k,529) * lu(k,1368)
         lu(k,1394) = lu(k,1394) - lu(k,530) * lu(k,1368)
         lu(k,1396) = lu(k,1396) - lu(k,531) * lu(k,1368)
         lu(k,1644) = lu(k,1644) - lu(k,525) * lu(k,1615)
         lu(k,1683) = lu(k,1683) - lu(k,526) * lu(k,1615)
         lu(k,1689) = lu(k,1689) - lu(k,527) * lu(k,1615)
         lu(k,1691) = lu(k,1691) - lu(k,528) * lu(k,1615)
         lu(k,1692) = lu(k,1692) - lu(k,529) * lu(k,1615)
         lu(k,1697) = lu(k,1697) - lu(k,530) * lu(k,1615)
         lu(k,1700) = lu(k,1700) - lu(k,531) * lu(k,1615)
         lu(k,2175) = - lu(k,525) * lu(k,2168)
         lu(k,2185) = lu(k,2185) - lu(k,526) * lu(k,2168)
         lu(k,2190) = - lu(k,527) * lu(k,2168)
         lu(k,2192) = lu(k,2192) - lu(k,528) * lu(k,2168)
         lu(k,2193) = lu(k,2193) - lu(k,529) * lu(k,2168)
         lu(k,2198) = - lu(k,530) * lu(k,2168)
         lu(k,2201) = lu(k,2201) - lu(k,531) * lu(k,2168)
                                                                        
         lu(k,532) = 1._r8 / lu(k,532)
         lu(k,533) = lu(k,533) * lu(k,532)
         lu(k,534) = lu(k,534) * lu(k,532)
         lu(k,535) = lu(k,535) * lu(k,532)
         lu(k,536) = lu(k,536) * lu(k,532)
         lu(k,537) = lu(k,537) * lu(k,532)
         lu(k,538) = lu(k,538) * lu(k,532)
         lu(k,539) = lu(k,539) * lu(k,532)
         lu(k,768) = lu(k,768) - lu(k,533) * lu(k,767)
         lu(k,769) = lu(k,769) - lu(k,534) * lu(k,767)
         lu(k,770) = - lu(k,535) * lu(k,767)
         lu(k,771) = lu(k,771) - lu(k,536) * lu(k,767)
         lu(k,772) = - lu(k,537) * lu(k,767)
         lu(k,774) = lu(k,774) - lu(k,538) * lu(k,767)
         lu(k,776) = - lu(k,539) * lu(k,767)
         lu(k,1641) = lu(k,1641) - lu(k,533) * lu(k,1616)
         lu(k,1666) = lu(k,1666) - lu(k,534) * lu(k,1616)
         lu(k,1671) = lu(k,1671) - lu(k,535) * lu(k,1616)
         lu(k,1689) = lu(k,1689) - lu(k,536) * lu(k,1616)
         lu(k,1691) = lu(k,1691) - lu(k,537) * lu(k,1616)
         lu(k,1694) = lu(k,1694) - lu(k,538) * lu(k,1616)
         lu(k,1703) = lu(k,1703) - lu(k,539) * lu(k,1616)
         lu(k,1907) = lu(k,1907) - lu(k,533) * lu(k,1891)
         lu(k,1924) = lu(k,1924) - lu(k,534) * lu(k,1891)
         lu(k,1929) = - lu(k,535) * lu(k,1891)
         lu(k,1945) = lu(k,1945) - lu(k,536) * lu(k,1891)
         lu(k,1947) = lu(k,1947) - lu(k,537) * lu(k,1891)
         lu(k,1950) = lu(k,1950) - lu(k,538) * lu(k,1891)
         lu(k,1959) = lu(k,1959) - lu(k,539) * lu(k,1891)
                                                                        
      end do
                                                                        
      end subroutine lu_fac11
                                                                        
      subroutine lu_fac12( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,540) = 1._r8 / lu(k,540)
         lu(k,541) = lu(k,541) * lu(k,540)
         lu(k,542) = lu(k,542) * lu(k,540)
         lu(k,543) = lu(k,543) * lu(k,540)
         lu(k,544) = lu(k,544) * lu(k,540)
         lu(k,545) = lu(k,545) * lu(k,540)
         lu(k,546) = lu(k,546) * lu(k,540)
         lu(k,547) = lu(k,547) * lu(k,540)
         lu(k,2174) = - lu(k,541) * lu(k,2169)
         lu(k,2186) = - lu(k,542) * lu(k,2169)
         lu(k,2188) = lu(k,2188) - lu(k,543) * lu(k,2169)
         lu(k,2193) = lu(k,2193) - lu(k,544) * lu(k,2169)
         lu(k,2201) = lu(k,2201) - lu(k,545) * lu(k,2169)
         lu(k,2202) = lu(k,2202) - lu(k,546) * lu(k,2169)
         lu(k,2203) = lu(k,2203) - lu(k,547) * lu(k,2169)
         lu(k,2208) = lu(k,2208) - lu(k,541) * lu(k,2207)
         lu(k,2211) = lu(k,2211) - lu(k,542) * lu(k,2207)
         lu(k,2212) = - lu(k,543) * lu(k,2207)
         lu(k,2217) = - lu(k,544) * lu(k,2207)
         lu(k,2225) = lu(k,2225) - lu(k,545) * lu(k,2207)
         lu(k,2226) = lu(k,2226) - lu(k,546) * lu(k,2207)
         lu(k,2227) = lu(k,2227) - lu(k,547) * lu(k,2207)
         lu(k,2233) = lu(k,2233) - lu(k,541) * lu(k,2231)
         lu(k,2241) = lu(k,2241) - lu(k,542) * lu(k,2231)
         lu(k,2243) = - lu(k,543) * lu(k,2231)
         lu(k,2248) = lu(k,2248) - lu(k,544) * lu(k,2231)
         lu(k,2256) = lu(k,2256) - lu(k,545) * lu(k,2231)
         lu(k,2257) = lu(k,2257) - lu(k,546) * lu(k,2231)
         lu(k,2258) = lu(k,2258) - lu(k,547) * lu(k,2231)
                                                                        
         lu(k,548) = 1._r8 / lu(k,548)
         lu(k,549) = lu(k,549) * lu(k,548)
         lu(k,550) = lu(k,550) * lu(k,548)
         lu(k,551) = lu(k,551) * lu(k,548)
         lu(k,552) = lu(k,552) * lu(k,548)
         lu(k,553) = lu(k,553) * lu(k,548)
         lu(k,554) = lu(k,554) * lu(k,548)
         lu(k,555) = lu(k,555) * lu(k,548)
         lu(k,1636) = lu(k,1636) - lu(k,549) * lu(k,1617)
         lu(k,1657) = lu(k,1657) - lu(k,550) * lu(k,1617)
         lu(k,1668) = lu(k,1668) - lu(k,551) * lu(k,1617)
         lu(k,1689) = lu(k,1689) - lu(k,552) * lu(k,1617)
         lu(k,1691) = lu(k,1691) - lu(k,553) * lu(k,1617)
         lu(k,1694) = lu(k,1694) - lu(k,554) * lu(k,1617)
         lu(k,1698) = lu(k,1698) - lu(k,555) * lu(k,1617)
         lu(k,1995) = - lu(k,549) * lu(k,1992)
         lu(k,1999) = lu(k,1999) - lu(k,550) * lu(k,1992)
         lu(k,2003) = lu(k,2003) - lu(k,551) * lu(k,1992)
         lu(k,2010) = lu(k,2010) - lu(k,552) * lu(k,1992)
         lu(k,2012) = lu(k,2012) - lu(k,553) * lu(k,1992)
         lu(k,2015) = lu(k,2015) - lu(k,554) * lu(k,1992)
         lu(k,2019) = lu(k,2019) - lu(k,555) * lu(k,1992)
         lu(k,2089) = - lu(k,549) * lu(k,2086)
         lu(k,2095) = lu(k,2095) - lu(k,550) * lu(k,2086)
         lu(k,2105) = lu(k,2105) - lu(k,551) * lu(k,2086)
         lu(k,2123) = lu(k,2123) - lu(k,552) * lu(k,2086)
         lu(k,2125) = lu(k,2125) - lu(k,553) * lu(k,2086)
         lu(k,2128) = lu(k,2128) - lu(k,554) * lu(k,2086)
         lu(k,2132) = lu(k,2132) - lu(k,555) * lu(k,2086)
                                                                        
         lu(k,556) = 1._r8 / lu(k,556)
         lu(k,557) = lu(k,557) * lu(k,556)
         lu(k,558) = lu(k,558) * lu(k,556)
         lu(k,559) = lu(k,559) * lu(k,556)
         lu(k,560) = lu(k,560) * lu(k,556)
         lu(k,561) = lu(k,561) * lu(k,556)
         lu(k,562) = lu(k,562) * lu(k,556)
         lu(k,563) = lu(k,563) * lu(k,556)
         lu(k,1278) = - lu(k,557) * lu(k,1274)
         lu(k,1281) = lu(k,1281) - lu(k,558) * lu(k,1274)
         lu(k,1283) = lu(k,1283) - lu(k,559) * lu(k,1274)
         lu(k,1284) = - lu(k,560) * lu(k,1274)
         lu(k,1293) = - lu(k,561) * lu(k,1274)
         lu(k,1295) = lu(k,1295) - lu(k,562) * lu(k,1274)
         lu(k,1298) = lu(k,1298) - lu(k,563) * lu(k,1274)
         lu(k,1650) = lu(k,1650) - lu(k,557) * lu(k,1618)
         lu(k,1667) = lu(k,1667) - lu(k,558) * lu(k,1618)
         lu(k,1671) = lu(k,1671) - lu(k,559) * lu(k,1618)
         lu(k,1672) = lu(k,1672) - lu(k,560) * lu(k,1618)
         lu(k,1687) = lu(k,1687) - lu(k,561) * lu(k,1618)
         lu(k,1691) = lu(k,1691) - lu(k,562) * lu(k,1618)
         lu(k,1694) = lu(k,1694) - lu(k,563) * lu(k,1618)
         lu(k,1806) = - lu(k,557) * lu(k,1789)
         lu(k,1818) = lu(k,1818) - lu(k,558) * lu(k,1789)
         lu(k,1822) = lu(k,1822) - lu(k,559) * lu(k,1789)
         lu(k,1823) = lu(k,1823) - lu(k,560) * lu(k,1789)
         lu(k,1836) = lu(k,1836) - lu(k,561) * lu(k,1789)
         lu(k,1840) = lu(k,1840) - lu(k,562) * lu(k,1789)
         lu(k,1843) = lu(k,1843) - lu(k,563) * lu(k,1789)
                                                                        
         lu(k,564) = 1._r8 / lu(k,564)
         lu(k,565) = lu(k,565) * lu(k,564)
         lu(k,566) = lu(k,566) * lu(k,564)
         lu(k,567) = lu(k,567) * lu(k,564)
         lu(k,568) = lu(k,568) * lu(k,564)
         lu(k,569) = lu(k,569) * lu(k,564)
         lu(k,570) = lu(k,570) * lu(k,564)
         lu(k,571) = lu(k,571) * lu(k,564)
         lu(k,1656) = lu(k,1656) - lu(k,565) * lu(k,1619)
         lu(k,1666) = lu(k,1666) - lu(k,566) * lu(k,1619)
         lu(k,1672) = lu(k,1672) - lu(k,567) * lu(k,1619)
         lu(k,1689) = lu(k,1689) - lu(k,568) * lu(k,1619)
         lu(k,1693) = lu(k,1693) - lu(k,569) * lu(k,1619)
         lu(k,1694) = lu(k,1694) - lu(k,570) * lu(k,1619)
         lu(k,1700) = lu(k,1700) - lu(k,571) * lu(k,1619)
         lu(k,1714) = lu(k,1714) - lu(k,565) * lu(k,1710)
         lu(k,1724) = lu(k,1724) - lu(k,566) * lu(k,1710)
         lu(k,1730) = - lu(k,567) * lu(k,1710)
         lu(k,1746) = lu(k,1746) - lu(k,568) * lu(k,1710)
         lu(k,1750) = lu(k,1750) - lu(k,569) * lu(k,1710)
         lu(k,1751) = lu(k,1751) - lu(k,570) * lu(k,1710)
         lu(k,1757) = lu(k,1757) - lu(k,571) * lu(k,1710)
         lu(k,1810) = lu(k,1810) - lu(k,565) * lu(k,1790)
         lu(k,1817) = lu(k,1817) - lu(k,566) * lu(k,1790)
         lu(k,1823) = lu(k,1823) - lu(k,567) * lu(k,1790)
         lu(k,1838) = lu(k,1838) - lu(k,568) * lu(k,1790)
         lu(k,1842) = lu(k,1842) - lu(k,569) * lu(k,1790)
         lu(k,1843) = lu(k,1843) - lu(k,570) * lu(k,1790)
         lu(k,1849) = lu(k,1849) - lu(k,571) * lu(k,1790)
                                                                        
         lu(k,572) = 1._r8 / lu(k,572)
         lu(k,573) = lu(k,573) * lu(k,572)
         lu(k,574) = lu(k,574) * lu(k,572)
         lu(k,575) = lu(k,575) * lu(k,572)
         lu(k,576) = lu(k,576) * lu(k,572)
         lu(k,577) = lu(k,577) * lu(k,572)
         lu(k,578) = lu(k,578) * lu(k,572)
         lu(k,579) = lu(k,579) * lu(k,572)
         lu(k,580) = lu(k,580) * lu(k,572)
         lu(k,1350) = lu(k,1350) - lu(k,573) * lu(k,1347)
         lu(k,1352) = - lu(k,574) * lu(k,1347)
         lu(k,1354) = lu(k,1354) - lu(k,575) * lu(k,1347)
         lu(k,1357) = lu(k,1357) - lu(k,576) * lu(k,1347)
         lu(k,1358) = lu(k,1358) - lu(k,577) * lu(k,1347)
         lu(k,1359) = lu(k,1359) - lu(k,578) * lu(k,1347)
         lu(k,1361) = lu(k,1361) - lu(k,579) * lu(k,1347)
         lu(k,1364) = lu(k,1364) - lu(k,580) * lu(k,1347)
         lu(k,1644) = lu(k,1644) - lu(k,573) * lu(k,1620)
         lu(k,1671) = lu(k,1671) - lu(k,574) * lu(k,1620)
         lu(k,1682) = lu(k,1682) - lu(k,575) * lu(k,1620)
         lu(k,1689) = lu(k,1689) - lu(k,576) * lu(k,1620)
         lu(k,1691) = lu(k,1691) - lu(k,577) * lu(k,1620)
         lu(k,1692) = lu(k,1692) - lu(k,578) * lu(k,1620)
         lu(k,1694) = lu(k,1694) - lu(k,579) * lu(k,1620)
         lu(k,1700) = lu(k,1700) - lu(k,580) * lu(k,1620)
         lu(k,2175) = lu(k,2175) - lu(k,573) * lu(k,2170)
         lu(k,2181) = - lu(k,574) * lu(k,2170)
         lu(k,2184) = lu(k,2184) - lu(k,575) * lu(k,2170)
         lu(k,2190) = lu(k,2190) - lu(k,576) * lu(k,2170)
         lu(k,2192) = lu(k,2192) - lu(k,577) * lu(k,2170)
         lu(k,2193) = lu(k,2193) - lu(k,578) * lu(k,2170)
         lu(k,2195) = lu(k,2195) - lu(k,579) * lu(k,2170)
         lu(k,2201) = lu(k,2201) - lu(k,580) * lu(k,2170)
                                                                        
      end do
                                                                        
      end subroutine lu_fac12
                                                                        
      subroutine lu_fac13( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,581) = 1._r8 / lu(k,581)
         lu(k,582) = lu(k,582) * lu(k,581)
         lu(k,583) = lu(k,583) * lu(k,581)
         lu(k,584) = lu(k,584) * lu(k,581)
         lu(k,585) = lu(k,585) * lu(k,581)
         lu(k,586) = lu(k,586) * lu(k,581)
         lu(k,587) = lu(k,587) * lu(k,581)
         lu(k,588) = lu(k,588) * lu(k,581)
         lu(k,589) = lu(k,589) * lu(k,581)
         lu(k,1249) = - lu(k,582) * lu(k,1245)
         lu(k,1251) = - lu(k,583) * lu(k,1245)
         lu(k,1252) = - lu(k,584) * lu(k,1245)
         lu(k,1261) = - lu(k,585) * lu(k,1245)
         lu(k,1262) = lu(k,1262) - lu(k,586) * lu(k,1245)
         lu(k,1263) = - lu(k,587) * lu(k,1245)
         lu(k,1266) = lu(k,1266) - lu(k,588) * lu(k,1245)
         lu(k,1269) = lu(k,1269) - lu(k,589) * lu(k,1245)
         lu(k,1667) = lu(k,1667) - lu(k,582) * lu(k,1621)
         lu(k,1671) = lu(k,1671) - lu(k,583) * lu(k,1621)
         lu(k,1672) = lu(k,1672) - lu(k,584) * lu(k,1621)
         lu(k,1687) = lu(k,1687) - lu(k,585) * lu(k,1621)
         lu(k,1689) = lu(k,1689) - lu(k,586) * lu(k,1621)
         lu(k,1691) = lu(k,1691) - lu(k,587) * lu(k,1621)
         lu(k,1694) = lu(k,1694) - lu(k,588) * lu(k,1621)
         lu(k,1700) = lu(k,1700) - lu(k,589) * lu(k,1621)
         lu(k,1818) = lu(k,1818) - lu(k,582) * lu(k,1791)
         lu(k,1822) = lu(k,1822) - lu(k,583) * lu(k,1791)
         lu(k,1823) = lu(k,1823) - lu(k,584) * lu(k,1791)
         lu(k,1836) = lu(k,1836) - lu(k,585) * lu(k,1791)
         lu(k,1838) = lu(k,1838) - lu(k,586) * lu(k,1791)
         lu(k,1840) = lu(k,1840) - lu(k,587) * lu(k,1791)
         lu(k,1843) = lu(k,1843) - lu(k,588) * lu(k,1791)
         lu(k,1849) = lu(k,1849) - lu(k,589) * lu(k,1791)
                                                                        
         lu(k,590) = 1._r8 / lu(k,590)
         lu(k,591) = lu(k,591) * lu(k,590)
         lu(k,592) = lu(k,592) * lu(k,590)
         lu(k,593) = lu(k,593) * lu(k,590)
         lu(k,659) = - lu(k,591) * lu(k,654)
         lu(k,661) = - lu(k,592) * lu(k,654)
         lu(k,664) = lu(k,664) - lu(k,593) * lu(k,654)
         lu(k,705) = - lu(k,591) * lu(k,700)
         lu(k,707) = lu(k,707) - lu(k,592) * lu(k,700)
         lu(k,711) = lu(k,711) - lu(k,593) * lu(k,700)
         lu(k,734) = - lu(k,591) * lu(k,729)
         lu(k,736) = - lu(k,592) * lu(k,729)
         lu(k,740) = lu(k,740) - lu(k,593) * lu(k,729)
         lu(k,750) = - lu(k,591) * lu(k,745)
         lu(k,752) = lu(k,752) - lu(k,592) * lu(k,745)
         lu(k,757) = lu(k,757) - lu(k,593) * lu(k,745)
         lu(k,1071) = - lu(k,591) * lu(k,1069)
         lu(k,1074) = - lu(k,592) * lu(k,1069)
         lu(k,1081) = lu(k,1081) - lu(k,593) * lu(k,1069)
         lu(k,1279) = - lu(k,591) * lu(k,1275)
         lu(k,1282) = - lu(k,592) * lu(k,1275)
         lu(k,1298) = lu(k,1298) - lu(k,593) * lu(k,1275)
         lu(k,1652) = - lu(k,591) * lu(k,1622)
         lu(k,1668) = lu(k,1668) - lu(k,592) * lu(k,1622)
         lu(k,1694) = lu(k,1694) - lu(k,593) * lu(k,1622)
         lu(k,1807) = lu(k,1807) - lu(k,591) * lu(k,1792)
         lu(k,1819) = lu(k,1819) - lu(k,592) * lu(k,1792)
         lu(k,1843) = lu(k,1843) - lu(k,593) * lu(k,1792)
                                                                        
         lu(k,594) = 1._r8 / lu(k,594)
         lu(k,595) = lu(k,595) * lu(k,594)
         lu(k,596) = lu(k,596) * lu(k,594)
         lu(k,597) = lu(k,597) * lu(k,594)
         lu(k,598) = lu(k,598) * lu(k,594)
         lu(k,599) = lu(k,599) * lu(k,594)
         lu(k,600) = lu(k,600) * lu(k,594)
         lu(k,601) = lu(k,601) * lu(k,594)
         lu(k,602) = lu(k,602) * lu(k,594)
         lu(k,1521) = lu(k,1521) - lu(k,595) * lu(k,1517)
         lu(k,1526) = lu(k,1526) - lu(k,596) * lu(k,1517)
         lu(k,1527) = lu(k,1527) - lu(k,597) * lu(k,1517)
         lu(k,1530) = lu(k,1530) - lu(k,598) * lu(k,1517)
         lu(k,1532) = lu(k,1532) - lu(k,599) * lu(k,1517)
         lu(k,1533) = lu(k,1533) - lu(k,600) * lu(k,1517)
         lu(k,1535) = lu(k,1535) - lu(k,601) * lu(k,1517)
         lu(k,1539) = lu(k,1539) - lu(k,602) * lu(k,1517)
         lu(k,1685) = lu(k,1685) - lu(k,595) * lu(k,1623)
         lu(k,1690) = lu(k,1690) - lu(k,596) * lu(k,1623)
         lu(k,1691) = lu(k,1691) - lu(k,597) * lu(k,1623)
         lu(k,1694) = lu(k,1694) - lu(k,598) * lu(k,1623)
         lu(k,1696) = lu(k,1696) - lu(k,599) * lu(k,1623)
         lu(k,1697) = lu(k,1697) - lu(k,600) * lu(k,1623)
         lu(k,1699) = lu(k,1699) - lu(k,601) * lu(k,1623)
         lu(k,1703) = lu(k,1703) - lu(k,602) * lu(k,1623)
         lu(k,2006) = lu(k,2006) - lu(k,595) * lu(k,1993)
         lu(k,2011) = lu(k,2011) - lu(k,596) * lu(k,1993)
         lu(k,2012) = lu(k,2012) - lu(k,597) * lu(k,1993)
         lu(k,2015) = lu(k,2015) - lu(k,598) * lu(k,1993)
         lu(k,2017) = lu(k,2017) - lu(k,599) * lu(k,1993)
         lu(k,2018) = lu(k,2018) - lu(k,600) * lu(k,1993)
         lu(k,2020) = lu(k,2020) - lu(k,601) * lu(k,1993)
         lu(k,2024) = lu(k,2024) - lu(k,602) * lu(k,1993)
                                                                        
         lu(k,603) = 1._r8 / lu(k,603)
         lu(k,604) = lu(k,604) * lu(k,603)
         lu(k,605) = lu(k,605) * lu(k,603)
         lu(k,606) = lu(k,606) * lu(k,603)
         lu(k,607) = lu(k,607) * lu(k,603)
         lu(k,608) = lu(k,608) * lu(k,603)
         lu(k,609) = lu(k,609) * lu(k,603)
         lu(k,1691) = lu(k,1691) - lu(k,604) * lu(k,1624)
         lu(k,1694) = lu(k,1694) - lu(k,605) * lu(k,1624)
         lu(k,1696) = lu(k,1696) - lu(k,606) * lu(k,1624)
         lu(k,1699) = lu(k,1699) - lu(k,607) * lu(k,1624)
         lu(k,1702) = lu(k,1702) - lu(k,608) * lu(k,1624)
         lu(k,1703) = lu(k,1703) - lu(k,609) * lu(k,1624)
         lu(k,1947) = lu(k,1947) - lu(k,604) * lu(k,1892)
         lu(k,1950) = lu(k,1950) - lu(k,605) * lu(k,1892)
         lu(k,1952) = lu(k,1952) - lu(k,606) * lu(k,1892)
         lu(k,1955) = lu(k,1955) - lu(k,607) * lu(k,1892)
         lu(k,1958) = lu(k,1958) - lu(k,608) * lu(k,1892)
         lu(k,1959) = lu(k,1959) - lu(k,609) * lu(k,1892)
         lu(k,2012) = lu(k,2012) - lu(k,604) * lu(k,1994)
         lu(k,2015) = lu(k,2015) - lu(k,605) * lu(k,1994)
         lu(k,2017) = lu(k,2017) - lu(k,606) * lu(k,1994)
         lu(k,2020) = lu(k,2020) - lu(k,607) * lu(k,1994)
         lu(k,2023) = - lu(k,608) * lu(k,1994)
         lu(k,2024) = lu(k,2024) - lu(k,609) * lu(k,1994)
         lu(k,2247) = lu(k,2247) - lu(k,604) * lu(k,2232)
         lu(k,2250) = lu(k,2250) - lu(k,605) * lu(k,2232)
         lu(k,2252) = lu(k,2252) - lu(k,606) * lu(k,2232)
         lu(k,2255) = lu(k,2255) - lu(k,607) * lu(k,2232)
         lu(k,2258) = lu(k,2258) - lu(k,608) * lu(k,2232)
         lu(k,2259) = - lu(k,609) * lu(k,2232)
                                                                        
         lu(k,610) = 1._r8 / lu(k,610)
         lu(k,611) = lu(k,611) * lu(k,610)
         lu(k,612) = lu(k,612) * lu(k,610)
         lu(k,613) = lu(k,613) * lu(k,610)
         lu(k,614) = lu(k,614) * lu(k,610)
         lu(k,615) = lu(k,615) * lu(k,610)
         lu(k,616) = lu(k,616) * lu(k,610)
         lu(k,1350) = lu(k,1350) - lu(k,611) * lu(k,1348)
         lu(k,1355) = lu(k,1355) - lu(k,612) * lu(k,1348)
         lu(k,1357) = lu(k,1357) - lu(k,613) * lu(k,1348)
         lu(k,1358) = lu(k,1358) - lu(k,614) * lu(k,1348)
         lu(k,1362) = lu(k,1362) - lu(k,615) * lu(k,1348)
         lu(k,1366) = - lu(k,616) * lu(k,1348)
         lu(k,1371) = lu(k,1371) - lu(k,611) * lu(k,1369)
         lu(k,1386) = lu(k,1386) - lu(k,612) * lu(k,1369)
         lu(k,1389) = lu(k,1389) - lu(k,613) * lu(k,1369)
         lu(k,1390) = lu(k,1390) - lu(k,614) * lu(k,1369)
         lu(k,1394) = lu(k,1394) - lu(k,615) * lu(k,1369)
         lu(k,1398) = - lu(k,616) * lu(k,1369)
         lu(k,1644) = lu(k,1644) - lu(k,611) * lu(k,1625)
         lu(k,1683) = lu(k,1683) - lu(k,612) * lu(k,1625)
         lu(k,1689) = lu(k,1689) - lu(k,613) * lu(k,1625)
         lu(k,1691) = lu(k,1691) - lu(k,614) * lu(k,1625)
         lu(k,1697) = lu(k,1697) - lu(k,615) * lu(k,1625)
         lu(k,1703) = lu(k,1703) - lu(k,616) * lu(k,1625)
         lu(k,1910) = lu(k,1910) - lu(k,611) * lu(k,1893)
         lu(k,1939) = lu(k,1939) - lu(k,612) * lu(k,1893)
         lu(k,1945) = lu(k,1945) - lu(k,613) * lu(k,1893)
         lu(k,1947) = lu(k,1947) - lu(k,614) * lu(k,1893)
         lu(k,1953) = lu(k,1953) - lu(k,615) * lu(k,1893)
         lu(k,1959) = lu(k,1959) - lu(k,616) * lu(k,1893)
                                                                        
      end do
                                                                        
      end subroutine lu_fac13
                                                                        
      subroutine lu_fac14( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,617) = 1._r8 / lu(k,617)
         lu(k,618) = lu(k,618) * lu(k,617)
         lu(k,619) = lu(k,619) * lu(k,617)
         lu(k,620) = lu(k,620) * lu(k,617)
         lu(k,621) = lu(k,621) * lu(k,617)
         lu(k,622) = lu(k,622) * lu(k,617)
         lu(k,914) = lu(k,914) - lu(k,618) * lu(k,910)
         lu(k,915) = - lu(k,619) * lu(k,910)
         lu(k,917) = lu(k,917) - lu(k,620) * lu(k,910)
         lu(k,919) = lu(k,919) - lu(k,621) * lu(k,910)
         lu(k,921) = lu(k,921) - lu(k,622) * lu(k,910)
         lu(k,1044) = lu(k,1044) - lu(k,618) * lu(k,1042)
         lu(k,1049) = lu(k,1049) - lu(k,619) * lu(k,1042)
         lu(k,1051) = lu(k,1051) - lu(k,620) * lu(k,1042)
         lu(k,1054) = lu(k,1054) - lu(k,621) * lu(k,1042)
         lu(k,1056) = lu(k,1056) - lu(k,622) * lu(k,1042)
         lu(k,1659) = lu(k,1659) - lu(k,618) * lu(k,1626)
         lu(k,1687) = lu(k,1687) - lu(k,619) * lu(k,1626)
         lu(k,1691) = lu(k,1691) - lu(k,620) * lu(k,1626)
         lu(k,1694) = lu(k,1694) - lu(k,621) * lu(k,1626)
         lu(k,1700) = lu(k,1700) - lu(k,622) * lu(k,1626)
         lu(k,1812) = lu(k,1812) - lu(k,618) * lu(k,1793)
         lu(k,1836) = lu(k,1836) - lu(k,619) * lu(k,1793)
         lu(k,1840) = lu(k,1840) - lu(k,620) * lu(k,1793)
         lu(k,1843) = lu(k,1843) - lu(k,621) * lu(k,1793)
         lu(k,1849) = lu(k,1849) - lu(k,622) * lu(k,1793)
         lu(k,2038) = lu(k,2038) - lu(k,618) * lu(k,2029)
         lu(k,2060) = lu(k,2060) - lu(k,619) * lu(k,2029)
         lu(k,2064) = lu(k,2064) - lu(k,620) * lu(k,2029)
         lu(k,2067) = lu(k,2067) - lu(k,621) * lu(k,2029)
         lu(k,2073) = lu(k,2073) - lu(k,622) * lu(k,2029)
                                                                        
         lu(k,625) = 1._r8 / lu(k,625)
         lu(k,626) = lu(k,626) * lu(k,625)
         lu(k,627) = lu(k,627) * lu(k,625)
         lu(k,628) = lu(k,628) * lu(k,625)
         lu(k,629) = lu(k,629) * lu(k,625)
         lu(k,630) = lu(k,630) * lu(k,625)
         lu(k,1691) = lu(k,1691) - lu(k,626) * lu(k,1627)
         lu(k,1693) = lu(k,1693) - lu(k,627) * lu(k,1627)
         lu(k,1694) = lu(k,1694) - lu(k,628) * lu(k,1627)
         lu(k,1698) = lu(k,1698) - lu(k,629) * lu(k,1627)
         lu(k,1700) = lu(k,1700) - lu(k,630) * lu(k,1627)
         lu(k,1840) = lu(k,1840) - lu(k,626) * lu(k,1794)
         lu(k,1842) = lu(k,1842) - lu(k,627) * lu(k,1794)
         lu(k,1843) = lu(k,1843) - lu(k,628) * lu(k,1794)
         lu(k,1847) = lu(k,1847) - lu(k,629) * lu(k,1794)
         lu(k,1849) = lu(k,1849) - lu(k,630) * lu(k,1794)
         lu(k,1947) = lu(k,1947) - lu(k,626) * lu(k,1894)
         lu(k,1949) = lu(k,1949) - lu(k,627) * lu(k,1894)
         lu(k,1950) = lu(k,1950) - lu(k,628) * lu(k,1894)
         lu(k,1954) = lu(k,1954) - lu(k,629) * lu(k,1894)
         lu(k,1956) = lu(k,1956) - lu(k,630) * lu(k,1894)
         lu(k,2125) = lu(k,2125) - lu(k,626) * lu(k,2087)
         lu(k,2127) = lu(k,2127) - lu(k,627) * lu(k,2087)
         lu(k,2128) = lu(k,2128) - lu(k,628) * lu(k,2087)
         lu(k,2132) = lu(k,2132) - lu(k,629) * lu(k,2087)
         lu(k,2134) = lu(k,2134) - lu(k,630) * lu(k,2087)
         lu(k,2192) = lu(k,2192) - lu(k,626) * lu(k,2171)
         lu(k,2194) = lu(k,2194) - lu(k,627) * lu(k,2171)
         lu(k,2195) = lu(k,2195) - lu(k,628) * lu(k,2171)
         lu(k,2199) = lu(k,2199) - lu(k,629) * lu(k,2171)
         lu(k,2201) = lu(k,2201) - lu(k,630) * lu(k,2171)
                                                                        
         lu(k,632) = 1._r8 / lu(k,632)
         lu(k,633) = lu(k,633) * lu(k,632)
         lu(k,634) = lu(k,634) * lu(k,632)
         lu(k,635) = lu(k,635) * lu(k,632)
         lu(k,636) = lu(k,636) * lu(k,632)
         lu(k,637) = lu(k,637) * lu(k,632)
         lu(k,638) = lu(k,638) * lu(k,632)
         lu(k,639) = lu(k,639) * lu(k,632)
         lu(k,640) = lu(k,640) * lu(k,632)
         lu(k,641) = lu(k,641) * lu(k,632)
         lu(k,897) = lu(k,897) - lu(k,633) * lu(k,895)
         lu(k,898) = lu(k,898) - lu(k,634) * lu(k,895)
         lu(k,899) = lu(k,899) - lu(k,635) * lu(k,895)
         lu(k,900) = lu(k,900) - lu(k,636) * lu(k,895)
         lu(k,901) = lu(k,901) - lu(k,637) * lu(k,895)
         lu(k,902) = lu(k,902) - lu(k,638) * lu(k,895)
         lu(k,903) = lu(k,903) - lu(k,639) * lu(k,895)
         lu(k,904) = lu(k,904) - lu(k,640) * lu(k,895)
         lu(k,906) = lu(k,906) - lu(k,641) * lu(k,895)
         lu(k,1632) = lu(k,1632) - lu(k,633) * lu(k,1628)
         lu(k,1646) = lu(k,1646) - lu(k,634) * lu(k,1628)
         lu(k,1654) = lu(k,1654) - lu(k,635) * lu(k,1628)
         lu(k,1656) = lu(k,1656) - lu(k,636) * lu(k,1628)
         lu(k,1666) = lu(k,1666) - lu(k,637) * lu(k,1628)
         lu(k,1683) = lu(k,1683) - lu(k,638) * lu(k,1628)
         lu(k,1689) = lu(k,1689) - lu(k,639) * lu(k,1628)
         lu(k,1691) = lu(k,1691) - lu(k,640) * lu(k,1628)
         lu(k,1694) = lu(k,1694) - lu(k,641) * lu(k,1628)
         lu(k,1898) = lu(k,1898) - lu(k,633) * lu(k,1895)
         lu(k,1912) = lu(k,1912) - lu(k,634) * lu(k,1895)
         lu(k,1916) = lu(k,1916) - lu(k,635) * lu(k,1895)
         lu(k,1918) = lu(k,1918) - lu(k,636) * lu(k,1895)
         lu(k,1924) = lu(k,1924) - lu(k,637) * lu(k,1895)
         lu(k,1939) = lu(k,1939) - lu(k,638) * lu(k,1895)
         lu(k,1945) = lu(k,1945) - lu(k,639) * lu(k,1895)
         lu(k,1947) = lu(k,1947) - lu(k,640) * lu(k,1895)
         lu(k,1950) = lu(k,1950) - lu(k,641) * lu(k,1895)
                                                                        
         lu(k,642) = 1._r8 / lu(k,642)
         lu(k,643) = lu(k,643) * lu(k,642)
         lu(k,644) = lu(k,644) * lu(k,642)
         lu(k,645) = lu(k,645) * lu(k,642)
         lu(k,646) = lu(k,646) * lu(k,642)
         lu(k,647) = lu(k,647) * lu(k,642)
         lu(k,648) = lu(k,648) * lu(k,642)
         lu(k,649) = lu(k,649) * lu(k,642)
         lu(k,650) = lu(k,650) * lu(k,642)
         lu(k,651) = lu(k,651) * lu(k,642)
         lu(k,1110) = lu(k,1110) - lu(k,643) * lu(k,1108)
         lu(k,1111) = lu(k,1111) - lu(k,644) * lu(k,1108)
         lu(k,1112) = lu(k,1112) - lu(k,645) * lu(k,1108)
         lu(k,1113) = lu(k,1113) - lu(k,646) * lu(k,1108)
         lu(k,1114) = lu(k,1114) - lu(k,647) * lu(k,1108)
         lu(k,1115) = lu(k,1115) - lu(k,648) * lu(k,1108)
         lu(k,1119) = lu(k,1119) - lu(k,649) * lu(k,1108)
         lu(k,1120) = - lu(k,650) * lu(k,1108)
         lu(k,1122) = lu(k,1122) - lu(k,651) * lu(k,1108)
         lu(k,1644) = lu(k,1644) - lu(k,643) * lu(k,1629)
         lu(k,1656) = lu(k,1656) - lu(k,644) * lu(k,1629)
         lu(k,1664) = lu(k,1664) - lu(k,645) * lu(k,1629)
         lu(k,1667) = lu(k,1667) - lu(k,646) * lu(k,1629)
         lu(k,1668) = lu(k,1668) - lu(k,647) * lu(k,1629)
         lu(k,1669) = lu(k,1669) - lu(k,648) * lu(k,1629)
         lu(k,1689) = lu(k,1689) - lu(k,649) * lu(k,1629)
         lu(k,1691) = lu(k,1691) - lu(k,650) * lu(k,1629)
         lu(k,1694) = lu(k,1694) - lu(k,651) * lu(k,1629)
         lu(k,1910) = lu(k,1910) - lu(k,643) * lu(k,1896)
         lu(k,1918) = lu(k,1918) - lu(k,644) * lu(k,1896)
         lu(k,1923) = - lu(k,645) * lu(k,1896)
         lu(k,1925) = lu(k,1925) - lu(k,646) * lu(k,1896)
         lu(k,1926) = lu(k,1926) - lu(k,647) * lu(k,1896)
         lu(k,1927) = lu(k,1927) - lu(k,648) * lu(k,1896)
         lu(k,1945) = lu(k,1945) - lu(k,649) * lu(k,1896)
         lu(k,1947) = lu(k,1947) - lu(k,650) * lu(k,1896)
         lu(k,1950) = lu(k,1950) - lu(k,651) * lu(k,1896)
                                                                        
         lu(k,655) = 1._r8 / lu(k,655)
         lu(k,656) = lu(k,656) * lu(k,655)
         lu(k,657) = lu(k,657) * lu(k,655)
         lu(k,658) = lu(k,658) * lu(k,655)
         lu(k,659) = lu(k,659) * lu(k,655)
         lu(k,660) = lu(k,660) * lu(k,655)
         lu(k,661) = lu(k,661) * lu(k,655)
         lu(k,662) = lu(k,662) * lu(k,655)
         lu(k,663) = lu(k,663) * lu(k,655)
         lu(k,664) = lu(k,664) * lu(k,655)
         lu(k,731) = lu(k,731) - lu(k,656) * lu(k,730)
         lu(k,732) = lu(k,732) - lu(k,657) * lu(k,730)
         lu(k,733) = lu(k,733) - lu(k,658) * lu(k,730)
         lu(k,734) = lu(k,734) - lu(k,659) * lu(k,730)
         lu(k,735) = lu(k,735) - lu(k,660) * lu(k,730)
         lu(k,736) = lu(k,736) - lu(k,661) * lu(k,730)
         lu(k,737) = lu(k,737) - lu(k,662) * lu(k,730)
         lu(k,738) = - lu(k,663) * lu(k,730)
         lu(k,740) = lu(k,740) - lu(k,664) * lu(k,730)
         lu(k,1637) = lu(k,1637) - lu(k,656) * lu(k,1630)
         lu(k,1638) = lu(k,1638) - lu(k,657) * lu(k,1630)
         lu(k,1640) = - lu(k,658) * lu(k,1630)
         lu(k,1652) = lu(k,1652) - lu(k,659) * lu(k,1630)
         lu(k,1660) = lu(k,1660) - lu(k,660) * lu(k,1630)
         lu(k,1668) = lu(k,1668) - lu(k,661) * lu(k,1630)
         lu(k,1676) = lu(k,1676) - lu(k,662) * lu(k,1630)
         lu(k,1691) = lu(k,1691) - lu(k,663) * lu(k,1630)
         lu(k,1694) = lu(k,1694) - lu(k,664) * lu(k,1630)
         lu(k,1903) = lu(k,1903) - lu(k,656) * lu(k,1897)
         lu(k,1904) = lu(k,1904) - lu(k,657) * lu(k,1897)
         lu(k,1906) = lu(k,1906) - lu(k,658) * lu(k,1897)
         lu(k,1915) = lu(k,1915) - lu(k,659) * lu(k,1897)
         lu(k,1921) = lu(k,1921) - lu(k,660) * lu(k,1897)
         lu(k,1926) = lu(k,1926) - lu(k,661) * lu(k,1897)
         lu(k,1933) = lu(k,1933) - lu(k,662) * lu(k,1897)
         lu(k,1947) = lu(k,1947) - lu(k,663) * lu(k,1897)
         lu(k,1950) = lu(k,1950) - lu(k,664) * lu(k,1897)
                                                                        
         lu(k,666) = 1._r8 / lu(k,666)
         lu(k,667) = lu(k,667) * lu(k,666)
         lu(k,668) = lu(k,668) * lu(k,666)
         lu(k,669) = lu(k,669) * lu(k,666)
         lu(k,670) = lu(k,670) * lu(k,666)
         lu(k,671) = lu(k,671) * lu(k,666)
         lu(k,672) = lu(k,672) * lu(k,666)
         lu(k,673) = lu(k,673) * lu(k,666)
         lu(k,674) = lu(k,674) * lu(k,666)
         lu(k,675) = lu(k,675) * lu(k,666)
         lu(k,897) = lu(k,897) - lu(k,667) * lu(k,896)
         lu(k,898) = lu(k,898) - lu(k,668) * lu(k,896)
         lu(k,900) = lu(k,900) - lu(k,669) * lu(k,896)
         lu(k,901) = lu(k,901) - lu(k,670) * lu(k,896)
         lu(k,902) = lu(k,902) - lu(k,671) * lu(k,896)
         lu(k,903) = lu(k,903) - lu(k,672) * lu(k,896)
         lu(k,904) = lu(k,904) - lu(k,673) * lu(k,896)
         lu(k,906) = lu(k,906) - lu(k,674) * lu(k,896)
         lu(k,908) = lu(k,908) - lu(k,675) * lu(k,896)
         lu(k,1632) = lu(k,1632) - lu(k,667) * lu(k,1631)
         lu(k,1646) = lu(k,1646) - lu(k,668) * lu(k,1631)
         lu(k,1656) = lu(k,1656) - lu(k,669) * lu(k,1631)
         lu(k,1666) = lu(k,1666) - lu(k,670) * lu(k,1631)
         lu(k,1683) = lu(k,1683) - lu(k,671) * lu(k,1631)
         lu(k,1689) = lu(k,1689) - lu(k,672) * lu(k,1631)
         lu(k,1691) = lu(k,1691) - lu(k,673) * lu(k,1631)
         lu(k,1694) = lu(k,1694) - lu(k,674) * lu(k,1631)
         lu(k,1700) = lu(k,1700) - lu(k,675) * lu(k,1631)
         lu(k,1796) = lu(k,1796) - lu(k,667) * lu(k,1795)
         lu(k,1805) = lu(k,1805) - lu(k,668) * lu(k,1795)
         lu(k,1810) = lu(k,1810) - lu(k,669) * lu(k,1795)
         lu(k,1817) = lu(k,1817) - lu(k,670) * lu(k,1795)
         lu(k,1833) = lu(k,1833) - lu(k,671) * lu(k,1795)
         lu(k,1838) = lu(k,1838) - lu(k,672) * lu(k,1795)
         lu(k,1840) = lu(k,1840) - lu(k,673) * lu(k,1795)
         lu(k,1843) = lu(k,1843) - lu(k,674) * lu(k,1795)
         lu(k,1849) = lu(k,1849) - lu(k,675) * lu(k,1795)
                                                                        
      end do
                                                                        
      end subroutine lu_fac14
                                                                        
      subroutine lu_fac15( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,677) = 1._r8 / lu(k,677)
         lu(k,678) = lu(k,678) * lu(k,677)
         lu(k,679) = lu(k,679) * lu(k,677)
         lu(k,680) = lu(k,680) * lu(k,677)
         lu(k,681) = lu(k,681) * lu(k,677)
         lu(k,682) = lu(k,682) * lu(k,677)
         lu(k,683) = lu(k,683) * lu(k,677)
         lu(k,901) = lu(k,901) - lu(k,678) * lu(k,897)
         lu(k,902) = lu(k,902) - lu(k,679) * lu(k,897)
         lu(k,904) = lu(k,904) - lu(k,680) * lu(k,897)
         lu(k,905) = lu(k,905) - lu(k,681) * lu(k,897)
         lu(k,906) = lu(k,906) - lu(k,682) * lu(k,897)
         lu(k,908) = lu(k,908) - lu(k,683) * lu(k,897)
         lu(k,1666) = lu(k,1666) - lu(k,678) * lu(k,1632)
         lu(k,1683) = lu(k,1683) - lu(k,679) * lu(k,1632)
         lu(k,1691) = lu(k,1691) - lu(k,680) * lu(k,1632)
         lu(k,1693) = lu(k,1693) - lu(k,681) * lu(k,1632)
         lu(k,1694) = lu(k,1694) - lu(k,682) * lu(k,1632)
         lu(k,1700) = lu(k,1700) - lu(k,683) * lu(k,1632)
         lu(k,1817) = lu(k,1817) - lu(k,678) * lu(k,1796)
         lu(k,1833) = lu(k,1833) - lu(k,679) * lu(k,1796)
         lu(k,1840) = lu(k,1840) - lu(k,680) * lu(k,1796)
         lu(k,1842) = lu(k,1842) - lu(k,681) * lu(k,1796)
         lu(k,1843) = lu(k,1843) - lu(k,682) * lu(k,1796)
         lu(k,1849) = lu(k,1849) - lu(k,683) * lu(k,1796)
         lu(k,1924) = lu(k,1924) - lu(k,678) * lu(k,1898)
         lu(k,1939) = lu(k,1939) - lu(k,679) * lu(k,1898)
         lu(k,1947) = lu(k,1947) - lu(k,680) * lu(k,1898)
         lu(k,1949) = lu(k,1949) - lu(k,681) * lu(k,1898)
         lu(k,1950) = lu(k,1950) - lu(k,682) * lu(k,1898)
         lu(k,1956) = lu(k,1956) - lu(k,683) * lu(k,1898)
                                                                        
         lu(k,684) = 1._r8 / lu(k,684)
         lu(k,685) = lu(k,685) * lu(k,684)
         lu(k,686) = lu(k,686) * lu(k,684)
         lu(k,687) = lu(k,687) * lu(k,684)
         lu(k,688) = lu(k,688) * lu(k,684)
         lu(k,1021) = lu(k,1021) - lu(k,685) * lu(k,1019)
         lu(k,1032) = lu(k,1032) - lu(k,686) * lu(k,1019)
         lu(k,1036) = lu(k,1036) - lu(k,687) * lu(k,1019)
         lu(k,1040) = - lu(k,688) * lu(k,1019)
         lu(k,1350) = lu(k,1350) - lu(k,685) * lu(k,1349)
         lu(k,1358) = lu(k,1358) - lu(k,686) * lu(k,1349)
         lu(k,1362) = lu(k,1362) - lu(k,687) * lu(k,1349)
         lu(k,1366) = lu(k,1366) - lu(k,688) * lu(k,1349)
         lu(k,1371) = lu(k,1371) - lu(k,685) * lu(k,1370)
         lu(k,1390) = lu(k,1390) - lu(k,686) * lu(k,1370)
         lu(k,1394) = lu(k,1394) - lu(k,687) * lu(k,1370)
         lu(k,1398) = lu(k,1398) - lu(k,688) * lu(k,1370)
         lu(k,1644) = lu(k,1644) - lu(k,685) * lu(k,1633)
         lu(k,1691) = lu(k,1691) - lu(k,686) * lu(k,1633)
         lu(k,1697) = lu(k,1697) - lu(k,687) * lu(k,1633)
         lu(k,1703) = lu(k,1703) - lu(k,688) * lu(k,1633)
         lu(k,1910) = lu(k,1910) - lu(k,685) * lu(k,1899)
         lu(k,1947) = lu(k,1947) - lu(k,686) * lu(k,1899)
         lu(k,1953) = lu(k,1953) - lu(k,687) * lu(k,1899)
         lu(k,1959) = lu(k,1959) - lu(k,688) * lu(k,1899)
         lu(k,2033) = lu(k,2033) - lu(k,685) * lu(k,2030)
         lu(k,2064) = lu(k,2064) - lu(k,686) * lu(k,2030)
         lu(k,2070) = lu(k,2070) - lu(k,687) * lu(k,2030)
         lu(k,2076) = lu(k,2076) - lu(k,688) * lu(k,2030)
         lu(k,2090) = lu(k,2090) - lu(k,685) * lu(k,2088)
         lu(k,2125) = lu(k,2125) - lu(k,686) * lu(k,2088)
         lu(k,2131) = lu(k,2131) - lu(k,687) * lu(k,2088)
         lu(k,2137) = - lu(k,688) * lu(k,2088)
                                                                        
         lu(k,690) = 1._r8 / lu(k,690)
         lu(k,691) = lu(k,691) * lu(k,690)
         lu(k,692) = lu(k,692) * lu(k,690)
         lu(k,693) = lu(k,693) * lu(k,690)
         lu(k,694) = lu(k,694) * lu(k,690)
         lu(k,695) = lu(k,695) * lu(k,690)
         lu(k,696) = lu(k,696) * lu(k,690)
         lu(k,1192) = - lu(k,691) * lu(k,1186)
         lu(k,1194) = - lu(k,692) * lu(k,1186)
         lu(k,1196) = - lu(k,693) * lu(k,1186)
         lu(k,1199) = lu(k,1199) - lu(k,694) * lu(k,1186)
         lu(k,1200) = lu(k,1200) - lu(k,695) * lu(k,1186)
         lu(k,1203) = lu(k,1203) - lu(k,696) * lu(k,1186)
         lu(k,1253) = - lu(k,691) * lu(k,1246)
         lu(k,1254) = lu(k,1254) - lu(k,692) * lu(k,1246)
         lu(k,1258) = lu(k,1258) - lu(k,693) * lu(k,1246)
         lu(k,1262) = lu(k,1262) - lu(k,694) * lu(k,1246)
         lu(k,1263) = lu(k,1263) - lu(k,695) * lu(k,1246)
         lu(k,1266) = lu(k,1266) - lu(k,696) * lu(k,1246)
         lu(k,1285) = lu(k,1285) - lu(k,691) * lu(k,1276)
         lu(k,1286) = - lu(k,692) * lu(k,1276)
         lu(k,1290) = - lu(k,693) * lu(k,1276)
         lu(k,1294) = lu(k,1294) - lu(k,694) * lu(k,1276)
         lu(k,1295) = lu(k,1295) - lu(k,695) * lu(k,1276)
         lu(k,1298) = lu(k,1298) - lu(k,696) * lu(k,1276)
         lu(k,1673) = lu(k,1673) - lu(k,691) * lu(k,1634)
         lu(k,1675) = lu(k,1675) - lu(k,692) * lu(k,1634)
         lu(k,1681) = lu(k,1681) - lu(k,693) * lu(k,1634)
         lu(k,1689) = lu(k,1689) - lu(k,694) * lu(k,1634)
         lu(k,1691) = lu(k,1691) - lu(k,695) * lu(k,1634)
         lu(k,1694) = lu(k,1694) - lu(k,696) * lu(k,1634)
         lu(k,1930) = lu(k,1930) - lu(k,691) * lu(k,1900)
         lu(k,1932) = - lu(k,692) * lu(k,1900)
         lu(k,1937) = - lu(k,693) * lu(k,1900)
         lu(k,1945) = lu(k,1945) - lu(k,694) * lu(k,1900)
         lu(k,1947) = lu(k,1947) - lu(k,695) * lu(k,1900)
         lu(k,1950) = lu(k,1950) - lu(k,696) * lu(k,1900)
                                                                        
         lu(k,701) = 1._r8 / lu(k,701)
         lu(k,702) = lu(k,702) * lu(k,701)
         lu(k,703) = lu(k,703) * lu(k,701)
         lu(k,704) = lu(k,704) * lu(k,701)
         lu(k,705) = lu(k,705) * lu(k,701)
         lu(k,706) = lu(k,706) * lu(k,701)
         lu(k,707) = lu(k,707) * lu(k,701)
         lu(k,708) = lu(k,708) * lu(k,701)
         lu(k,709) = lu(k,709) * lu(k,701)
         lu(k,710) = lu(k,710) * lu(k,701)
         lu(k,711) = lu(k,711) * lu(k,701)
         lu(k,747) = lu(k,747) - lu(k,702) * lu(k,746)
         lu(k,748) = lu(k,748) - lu(k,703) * lu(k,746)
         lu(k,749) = lu(k,749) - lu(k,704) * lu(k,746)
         lu(k,750) = lu(k,750) - lu(k,705) * lu(k,746)
         lu(k,751) = lu(k,751) - lu(k,706) * lu(k,746)
         lu(k,752) = lu(k,752) - lu(k,707) * lu(k,746)
         lu(k,753) = lu(k,753) - lu(k,708) * lu(k,746)
         lu(k,754) = lu(k,754) - lu(k,709) * lu(k,746)
         lu(k,755) = - lu(k,710) * lu(k,746)
         lu(k,757) = lu(k,757) - lu(k,711) * lu(k,746)
         lu(k,1637) = lu(k,1637) - lu(k,702) * lu(k,1635)
         lu(k,1639) = lu(k,1639) - lu(k,703) * lu(k,1635)
         lu(k,1640) = lu(k,1640) - lu(k,704) * lu(k,1635)
         lu(k,1652) = lu(k,1652) - lu(k,705) * lu(k,1635)
         lu(k,1660) = lu(k,1660) - lu(k,706) * lu(k,1635)
         lu(k,1668) = lu(k,1668) - lu(k,707) * lu(k,1635)
         lu(k,1676) = lu(k,1676) - lu(k,708) * lu(k,1635)
         lu(k,1683) = lu(k,1683) - lu(k,709) * lu(k,1635)
         lu(k,1691) = lu(k,1691) - lu(k,710) * lu(k,1635)
         lu(k,1694) = lu(k,1694) - lu(k,711) * lu(k,1635)
         lu(k,1903) = lu(k,1903) - lu(k,702) * lu(k,1901)
         lu(k,1905) = lu(k,1905) - lu(k,703) * lu(k,1901)
         lu(k,1906) = lu(k,1906) - lu(k,704) * lu(k,1901)
         lu(k,1915) = lu(k,1915) - lu(k,705) * lu(k,1901)
         lu(k,1921) = lu(k,1921) - lu(k,706) * lu(k,1901)
         lu(k,1926) = lu(k,1926) - lu(k,707) * lu(k,1901)
         lu(k,1933) = lu(k,1933) - lu(k,708) * lu(k,1901)
         lu(k,1939) = lu(k,1939) - lu(k,709) * lu(k,1901)
         lu(k,1947) = lu(k,1947) - lu(k,710) * lu(k,1901)
         lu(k,1950) = lu(k,1950) - lu(k,711) * lu(k,1901)
                                                                        
         lu(k,714) = 1._r8 / lu(k,714)
         lu(k,715) = lu(k,715) * lu(k,714)
         lu(k,716) = lu(k,716) * lu(k,714)
         lu(k,717) = lu(k,717) * lu(k,714)
         lu(k,718) = lu(k,718) * lu(k,714)
         lu(k,719) = lu(k,719) * lu(k,714)
         lu(k,720) = lu(k,720) * lu(k,714)
         lu(k,1667) = lu(k,1667) - lu(k,715) * lu(k,1636)
         lu(k,1689) = lu(k,1689) - lu(k,716) * lu(k,1636)
         lu(k,1691) = lu(k,1691) - lu(k,717) * lu(k,1636)
         lu(k,1693) = lu(k,1693) - lu(k,718) * lu(k,1636)
         lu(k,1694) = lu(k,1694) - lu(k,719) * lu(k,1636)
         lu(k,1700) = lu(k,1700) - lu(k,720) * lu(k,1636)
         lu(k,1818) = lu(k,1818) - lu(k,715) * lu(k,1797)
         lu(k,1838) = lu(k,1838) - lu(k,716) * lu(k,1797)
         lu(k,1840) = lu(k,1840) - lu(k,717) * lu(k,1797)
         lu(k,1842) = lu(k,1842) - lu(k,718) * lu(k,1797)
         lu(k,1843) = lu(k,1843) - lu(k,719) * lu(k,1797)
         lu(k,1849) = lu(k,1849) - lu(k,720) * lu(k,1797)
         lu(k,1925) = lu(k,1925) - lu(k,715) * lu(k,1902)
         lu(k,1945) = lu(k,1945) - lu(k,716) * lu(k,1902)
         lu(k,1947) = lu(k,1947) - lu(k,717) * lu(k,1902)
         lu(k,1949) = lu(k,1949) - lu(k,718) * lu(k,1902)
         lu(k,1950) = lu(k,1950) - lu(k,719) * lu(k,1902)
         lu(k,1956) = lu(k,1956) - lu(k,720) * lu(k,1902)
         lu(k,2002) = - lu(k,715) * lu(k,1995)
         lu(k,2010) = lu(k,2010) - lu(k,716) * lu(k,1995)
         lu(k,2012) = lu(k,2012) - lu(k,717) * lu(k,1995)
         lu(k,2014) = - lu(k,718) * lu(k,1995)
         lu(k,2015) = lu(k,2015) - lu(k,719) * lu(k,1995)
         lu(k,2021) = - lu(k,720) * lu(k,1995)
         lu(k,2104) = - lu(k,715) * lu(k,2089)
         lu(k,2123) = lu(k,2123) - lu(k,716) * lu(k,2089)
         lu(k,2125) = lu(k,2125) - lu(k,717) * lu(k,2089)
         lu(k,2127) = lu(k,2127) - lu(k,718) * lu(k,2089)
         lu(k,2128) = lu(k,2128) - lu(k,719) * lu(k,2089)
         lu(k,2134) = lu(k,2134) - lu(k,720) * lu(k,2089)
                                                                        
      end do
                                                                        
      end subroutine lu_fac15
                                                                        
      subroutine lu_fac16( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,721) = 1._r8 / lu(k,721)
         lu(k,722) = lu(k,722) * lu(k,721)
         lu(k,723) = lu(k,723) * lu(k,721)
         lu(k,724) = lu(k,724) * lu(k,721)
         lu(k,725) = lu(k,725) * lu(k,721)
         lu(k,726) = lu(k,726) * lu(k,721)
         lu(k,735) = lu(k,735) - lu(k,722) * lu(k,731)
         lu(k,736) = lu(k,736) - lu(k,723) * lu(k,731)
         lu(k,739) = lu(k,739) - lu(k,724) * lu(k,731)
         lu(k,740) = lu(k,740) - lu(k,725) * lu(k,731)
         lu(k,741) = lu(k,741) - lu(k,726) * lu(k,731)
         lu(k,751) = lu(k,751) - lu(k,722) * lu(k,747)
         lu(k,752) = lu(k,752) - lu(k,723) * lu(k,747)
         lu(k,756) = lu(k,756) - lu(k,724) * lu(k,747)
         lu(k,757) = lu(k,757) - lu(k,725) * lu(k,747)
         lu(k,758) = lu(k,758) - lu(k,726) * lu(k,747)
         lu(k,1660) = lu(k,1660) - lu(k,722) * lu(k,1637)
         lu(k,1668) = lu(k,1668) - lu(k,723) * lu(k,1637)
         lu(k,1693) = lu(k,1693) - lu(k,724) * lu(k,1637)
         lu(k,1694) = lu(k,1694) - lu(k,725) * lu(k,1637)
         lu(k,1700) = lu(k,1700) - lu(k,726) * lu(k,1637)
         lu(k,1813) = lu(k,1813) - lu(k,722) * lu(k,1798)
         lu(k,1819) = lu(k,1819) - lu(k,723) * lu(k,1798)
         lu(k,1842) = lu(k,1842) - lu(k,724) * lu(k,1798)
         lu(k,1843) = lu(k,1843) - lu(k,725) * lu(k,1798)
         lu(k,1849) = lu(k,1849) - lu(k,726) * lu(k,1798)
         lu(k,1921) = lu(k,1921) - lu(k,722) * lu(k,1903)
         lu(k,1926) = lu(k,1926) - lu(k,723) * lu(k,1903)
         lu(k,1949) = lu(k,1949) - lu(k,724) * lu(k,1903)
         lu(k,1950) = lu(k,1950) - lu(k,725) * lu(k,1903)
         lu(k,1956) = lu(k,1956) - lu(k,726) * lu(k,1903)
         lu(k,2179) = - lu(k,722) * lu(k,2172)
         lu(k,2180) = - lu(k,723) * lu(k,2172)
         lu(k,2194) = lu(k,2194) - lu(k,724) * lu(k,2172)
         lu(k,2195) = lu(k,2195) - lu(k,725) * lu(k,2172)
         lu(k,2201) = lu(k,2201) - lu(k,726) * lu(k,2172)
                                                                        
         lu(k,732) = 1._r8 / lu(k,732)
         lu(k,733) = lu(k,733) * lu(k,732)
         lu(k,734) = lu(k,734) * lu(k,732)
         lu(k,735) = lu(k,735) * lu(k,732)
         lu(k,736) = lu(k,736) * lu(k,732)
         lu(k,737) = lu(k,737) * lu(k,732)
         lu(k,738) = lu(k,738) * lu(k,732)
         lu(k,739) = lu(k,739) * lu(k,732)
         lu(k,740) = lu(k,740) * lu(k,732)
         lu(k,741) = lu(k,741) * lu(k,732)
         lu(k,1640) = lu(k,1640) - lu(k,733) * lu(k,1638)
         lu(k,1652) = lu(k,1652) - lu(k,734) * lu(k,1638)
         lu(k,1660) = lu(k,1660) - lu(k,735) * lu(k,1638)
         lu(k,1668) = lu(k,1668) - lu(k,736) * lu(k,1638)
         lu(k,1676) = lu(k,1676) - lu(k,737) * lu(k,1638)
         lu(k,1691) = lu(k,1691) - lu(k,738) * lu(k,1638)
         lu(k,1693) = lu(k,1693) - lu(k,739) * lu(k,1638)
         lu(k,1694) = lu(k,1694) - lu(k,740) * lu(k,1638)
         lu(k,1700) = lu(k,1700) - lu(k,741) * lu(k,1638)
         lu(k,1801) = lu(k,1801) - lu(k,733) * lu(k,1799)
         lu(k,1807) = lu(k,1807) - lu(k,734) * lu(k,1799)
         lu(k,1813) = lu(k,1813) - lu(k,735) * lu(k,1799)
         lu(k,1819) = lu(k,1819) - lu(k,736) * lu(k,1799)
         lu(k,1827) = lu(k,1827) - lu(k,737) * lu(k,1799)
         lu(k,1840) = lu(k,1840) - lu(k,738) * lu(k,1799)
         lu(k,1842) = lu(k,1842) - lu(k,739) * lu(k,1799)
         lu(k,1843) = lu(k,1843) - lu(k,740) * lu(k,1799)
         lu(k,1849) = lu(k,1849) - lu(k,741) * lu(k,1799)
         lu(k,1906) = lu(k,1906) - lu(k,733) * lu(k,1904)
         lu(k,1915) = lu(k,1915) - lu(k,734) * lu(k,1904)
         lu(k,1921) = lu(k,1921) - lu(k,735) * lu(k,1904)
         lu(k,1926) = lu(k,1926) - lu(k,736) * lu(k,1904)
         lu(k,1933) = lu(k,1933) - lu(k,737) * lu(k,1904)
         lu(k,1947) = lu(k,1947) - lu(k,738) * lu(k,1904)
         lu(k,1949) = lu(k,1949) - lu(k,739) * lu(k,1904)
         lu(k,1950) = lu(k,1950) - lu(k,740) * lu(k,1904)
         lu(k,1956) = lu(k,1956) - lu(k,741) * lu(k,1904)
                                                                        
         lu(k,748) = 1._r8 / lu(k,748)
         lu(k,749) = lu(k,749) * lu(k,748)
         lu(k,750) = lu(k,750) * lu(k,748)
         lu(k,751) = lu(k,751) * lu(k,748)
         lu(k,752) = lu(k,752) * lu(k,748)
         lu(k,753) = lu(k,753) * lu(k,748)
         lu(k,754) = lu(k,754) * lu(k,748)
         lu(k,755) = lu(k,755) * lu(k,748)
         lu(k,756) = lu(k,756) * lu(k,748)
         lu(k,757) = lu(k,757) * lu(k,748)
         lu(k,758) = lu(k,758) * lu(k,748)
         lu(k,1640) = lu(k,1640) - lu(k,749) * lu(k,1639)
         lu(k,1652) = lu(k,1652) - lu(k,750) * lu(k,1639)
         lu(k,1660) = lu(k,1660) - lu(k,751) * lu(k,1639)
         lu(k,1668) = lu(k,1668) - lu(k,752) * lu(k,1639)
         lu(k,1676) = lu(k,1676) - lu(k,753) * lu(k,1639)
         lu(k,1683) = lu(k,1683) - lu(k,754) * lu(k,1639)
         lu(k,1691) = lu(k,1691) - lu(k,755) * lu(k,1639)
         lu(k,1693) = lu(k,1693) - lu(k,756) * lu(k,1639)
         lu(k,1694) = lu(k,1694) - lu(k,757) * lu(k,1639)
         lu(k,1700) = lu(k,1700) - lu(k,758) * lu(k,1639)
         lu(k,1801) = lu(k,1801) - lu(k,749) * lu(k,1800)
         lu(k,1807) = lu(k,1807) - lu(k,750) * lu(k,1800)
         lu(k,1813) = lu(k,1813) - lu(k,751) * lu(k,1800)
         lu(k,1819) = lu(k,1819) - lu(k,752) * lu(k,1800)
         lu(k,1827) = lu(k,1827) - lu(k,753) * lu(k,1800)
         lu(k,1833) = lu(k,1833) - lu(k,754) * lu(k,1800)
         lu(k,1840) = lu(k,1840) - lu(k,755) * lu(k,1800)
         lu(k,1842) = lu(k,1842) - lu(k,756) * lu(k,1800)
         lu(k,1843) = lu(k,1843) - lu(k,757) * lu(k,1800)
         lu(k,1849) = lu(k,1849) - lu(k,758) * lu(k,1800)
         lu(k,1906) = lu(k,1906) - lu(k,749) * lu(k,1905)
         lu(k,1915) = lu(k,1915) - lu(k,750) * lu(k,1905)
         lu(k,1921) = lu(k,1921) - lu(k,751) * lu(k,1905)
         lu(k,1926) = lu(k,1926) - lu(k,752) * lu(k,1905)
         lu(k,1933) = lu(k,1933) - lu(k,753) * lu(k,1905)
         lu(k,1939) = lu(k,1939) - lu(k,754) * lu(k,1905)
         lu(k,1947) = lu(k,1947) - lu(k,755) * lu(k,1905)
         lu(k,1949) = lu(k,1949) - lu(k,756) * lu(k,1905)
         lu(k,1950) = lu(k,1950) - lu(k,757) * lu(k,1905)
         lu(k,1956) = lu(k,1956) - lu(k,758) * lu(k,1905)
                                                                        
         lu(k,759) = 1._r8 / lu(k,759)
         lu(k,760) = lu(k,760) * lu(k,759)
         lu(k,761) = lu(k,761) * lu(k,759)
         lu(k,762) = lu(k,762) * lu(k,759)
         lu(k,763) = lu(k,763) * lu(k,759)
         lu(k,764) = lu(k,764) * lu(k,759)
         lu(k,765) = lu(k,765) * lu(k,759)
         lu(k,766) = lu(k,766) * lu(k,759)
         lu(k,1668) = lu(k,1668) - lu(k,760) * lu(k,1640)
         lu(k,1676) = lu(k,1676) - lu(k,761) * lu(k,1640)
         lu(k,1691) = lu(k,1691) - lu(k,762) * lu(k,1640)
         lu(k,1693) = lu(k,1693) - lu(k,763) * lu(k,1640)
         lu(k,1694) = lu(k,1694) - lu(k,764) * lu(k,1640)
         lu(k,1697) = lu(k,1697) - lu(k,765) * lu(k,1640)
         lu(k,1700) = lu(k,1700) - lu(k,766) * lu(k,1640)
         lu(k,1819) = lu(k,1819) - lu(k,760) * lu(k,1801)
         lu(k,1827) = lu(k,1827) - lu(k,761) * lu(k,1801)
         lu(k,1840) = lu(k,1840) - lu(k,762) * lu(k,1801)
         lu(k,1842) = lu(k,1842) - lu(k,763) * lu(k,1801)
         lu(k,1843) = lu(k,1843) - lu(k,764) * lu(k,1801)
         lu(k,1846) = lu(k,1846) - lu(k,765) * lu(k,1801)
         lu(k,1849) = lu(k,1849) - lu(k,766) * lu(k,1801)
         lu(k,1926) = lu(k,1926) - lu(k,760) * lu(k,1906)
         lu(k,1933) = lu(k,1933) - lu(k,761) * lu(k,1906)
         lu(k,1947) = lu(k,1947) - lu(k,762) * lu(k,1906)
         lu(k,1949) = lu(k,1949) - lu(k,763) * lu(k,1906)
         lu(k,1950) = lu(k,1950) - lu(k,764) * lu(k,1906)
         lu(k,1953) = lu(k,1953) - lu(k,765) * lu(k,1906)
         lu(k,1956) = lu(k,1956) - lu(k,766) * lu(k,1906)
         lu(k,2180) = lu(k,2180) - lu(k,760) * lu(k,2173)
         lu(k,2182) = - lu(k,761) * lu(k,2173)
         lu(k,2192) = lu(k,2192) - lu(k,762) * lu(k,2173)
         lu(k,2194) = lu(k,2194) - lu(k,763) * lu(k,2173)
         lu(k,2195) = lu(k,2195) - lu(k,764) * lu(k,2173)
         lu(k,2198) = lu(k,2198) - lu(k,765) * lu(k,2173)
         lu(k,2201) = lu(k,2201) - lu(k,766) * lu(k,2173)
                                                                        
         lu(k,768) = 1._r8 / lu(k,768)
         lu(k,769) = lu(k,769) * lu(k,768)
         lu(k,770) = lu(k,770) * lu(k,768)
         lu(k,771) = lu(k,771) * lu(k,768)
         lu(k,772) = lu(k,772) * lu(k,768)
         lu(k,773) = lu(k,773) * lu(k,768)
         lu(k,774) = lu(k,774) * lu(k,768)
         lu(k,775) = lu(k,775) * lu(k,768)
         lu(k,776) = lu(k,776) * lu(k,768)
         lu(k,1025) = lu(k,1025) - lu(k,769) * lu(k,1020)
         lu(k,1027) = - lu(k,770) * lu(k,1020)
         lu(k,1031) = lu(k,1031) - lu(k,771) * lu(k,1020)
         lu(k,1032) = lu(k,1032) - lu(k,772) * lu(k,1020)
         lu(k,1034) = - lu(k,773) * lu(k,1020)
         lu(k,1035) = lu(k,1035) - lu(k,774) * lu(k,1020)
         lu(k,1038) = - lu(k,775) * lu(k,1020)
         lu(k,1040) = lu(k,1040) - lu(k,776) * lu(k,1020)
         lu(k,1666) = lu(k,1666) - lu(k,769) * lu(k,1641)
         lu(k,1671) = lu(k,1671) - lu(k,770) * lu(k,1641)
         lu(k,1689) = lu(k,1689) - lu(k,771) * lu(k,1641)
         lu(k,1691) = lu(k,1691) - lu(k,772) * lu(k,1641)
         lu(k,1693) = lu(k,1693) - lu(k,773) * lu(k,1641)
         lu(k,1694) = lu(k,1694) - lu(k,774) * lu(k,1641)
         lu(k,1700) = lu(k,1700) - lu(k,775) * lu(k,1641)
         lu(k,1703) = lu(k,1703) - lu(k,776) * lu(k,1641)
         lu(k,1817) = lu(k,1817) - lu(k,769) * lu(k,1802)
         lu(k,1822) = lu(k,1822) - lu(k,770) * lu(k,1802)
         lu(k,1838) = lu(k,1838) - lu(k,771) * lu(k,1802)
         lu(k,1840) = lu(k,1840) - lu(k,772) * lu(k,1802)
         lu(k,1842) = lu(k,1842) - lu(k,773) * lu(k,1802)
         lu(k,1843) = lu(k,1843) - lu(k,774) * lu(k,1802)
         lu(k,1849) = lu(k,1849) - lu(k,775) * lu(k,1802)
         lu(k,1852) = - lu(k,776) * lu(k,1802)
         lu(k,1924) = lu(k,1924) - lu(k,769) * lu(k,1907)
         lu(k,1929) = lu(k,1929) - lu(k,770) * lu(k,1907)
         lu(k,1945) = lu(k,1945) - lu(k,771) * lu(k,1907)
         lu(k,1947) = lu(k,1947) - lu(k,772) * lu(k,1907)
         lu(k,1949) = lu(k,1949) - lu(k,773) * lu(k,1907)
         lu(k,1950) = lu(k,1950) - lu(k,774) * lu(k,1907)
         lu(k,1956) = lu(k,1956) - lu(k,775) * lu(k,1907)
         lu(k,1959) = lu(k,1959) - lu(k,776) * lu(k,1907)
                                                                        
      end do
                                                                        
      end subroutine lu_fac16
                                                                        
      subroutine lu_fac17( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,778) = 1._r8 / lu(k,778)
         lu(k,779) = lu(k,779) * lu(k,778)
         lu(k,780) = lu(k,780) * lu(k,778)
         lu(k,781) = lu(k,781) * lu(k,778)
         lu(k,782) = lu(k,782) * lu(k,778)
         lu(k,783) = lu(k,783) * lu(k,778)
         lu(k,784) = lu(k,784) * lu(k,778)
         lu(k,785) = lu(k,785) * lu(k,778)
         lu(k,1941) = lu(k,1941) - lu(k,779) * lu(k,1908)
         lu(k,1947) = lu(k,1947) - lu(k,780) * lu(k,1908)
         lu(k,1952) = lu(k,1952) - lu(k,781) * lu(k,1908)
         lu(k,1955) = lu(k,1955) - lu(k,782) * lu(k,1908)
         lu(k,1957) = lu(k,1957) - lu(k,783) * lu(k,1908)
         lu(k,1958) = lu(k,1958) - lu(k,784) * lu(k,1908)
         lu(k,1959) = lu(k,1959) - lu(k,785) * lu(k,1908)
         lu(k,2143) = lu(k,2143) - lu(k,779) * lu(k,2140)
         lu(k,2148) = lu(k,2148) - lu(k,780) * lu(k,2140)
         lu(k,2153) = lu(k,2153) - lu(k,781) * lu(k,2140)
         lu(k,2156) = lu(k,2156) - lu(k,782) * lu(k,2140)
         lu(k,2158) = - lu(k,783) * lu(k,2140)
         lu(k,2159) = lu(k,2159) - lu(k,784) * lu(k,2140)
         lu(k,2160) = lu(k,2160) - lu(k,785) * lu(k,2140)
         lu(k,2186) = lu(k,2186) - lu(k,779) * lu(k,2174)
         lu(k,2192) = lu(k,2192) - lu(k,780) * lu(k,2174)
         lu(k,2197) = - lu(k,781) * lu(k,2174)
         lu(k,2200) = - lu(k,782) * lu(k,2174)
         lu(k,2202) = lu(k,2202) - lu(k,783) * lu(k,2174)
         lu(k,2203) = lu(k,2203) - lu(k,784) * lu(k,2174)
         lu(k,2204) = lu(k,2204) - lu(k,785) * lu(k,2174)
         lu(k,2211) = lu(k,2211) - lu(k,779) * lu(k,2208)
         lu(k,2216) = lu(k,2216) - lu(k,780) * lu(k,2208)
         lu(k,2221) = lu(k,2221) - lu(k,781) * lu(k,2208)
         lu(k,2224) = - lu(k,782) * lu(k,2208)
         lu(k,2226) = lu(k,2226) - lu(k,783) * lu(k,2208)
         lu(k,2227) = lu(k,2227) - lu(k,784) * lu(k,2208)
         lu(k,2228) = - lu(k,785) * lu(k,2208)
         lu(k,2241) = lu(k,2241) - lu(k,779) * lu(k,2233)
         lu(k,2247) = lu(k,2247) - lu(k,780) * lu(k,2233)
         lu(k,2252) = lu(k,2252) - lu(k,781) * lu(k,2233)
         lu(k,2255) = lu(k,2255) - lu(k,782) * lu(k,2233)
         lu(k,2257) = lu(k,2257) - lu(k,783) * lu(k,2233)
         lu(k,2258) = lu(k,2258) - lu(k,784) * lu(k,2233)
         lu(k,2259) = lu(k,2259) - lu(k,785) * lu(k,2233)
                                                                        
         lu(k,786) = 1._r8 / lu(k,786)
         lu(k,787) = lu(k,787) * lu(k,786)
         lu(k,788) = lu(k,788) * lu(k,786)
         lu(k,789) = lu(k,789) * lu(k,786)
         lu(k,817) = lu(k,817) - lu(k,787) * lu(k,814)
         lu(k,818) = lu(k,818) - lu(k,788) * lu(k,814)
         lu(k,820) = lu(k,820) - lu(k,789) * lu(k,814)
         lu(k,916) = lu(k,916) - lu(k,787) * lu(k,911)
         lu(k,917) = lu(k,917) - lu(k,788) * lu(k,911)
         lu(k,919) = lu(k,919) - lu(k,789) * lu(k,911)
         lu(k,1050) = lu(k,1050) - lu(k,787) * lu(k,1043)
         lu(k,1051) = lu(k,1051) - lu(k,788) * lu(k,1043)
         lu(k,1054) = lu(k,1054) - lu(k,789) * lu(k,1043)
         lu(k,1119) = lu(k,1119) - lu(k,787) * lu(k,1109)
         lu(k,1120) = lu(k,1120) - lu(k,788) * lu(k,1109)
         lu(k,1122) = lu(k,1122) - lu(k,789) * lu(k,1109)
         lu(k,1133) = lu(k,1133) - lu(k,787) * lu(k,1128)
         lu(k,1134) = lu(k,1134) - lu(k,788) * lu(k,1128)
         lu(k,1136) = lu(k,1136) - lu(k,789) * lu(k,1128)
         lu(k,1176) = lu(k,1176) - lu(k,787) * lu(k,1167)
         lu(k,1177) = lu(k,1177) - lu(k,788) * lu(k,1167)
         lu(k,1180) = lu(k,1180) - lu(k,789) * lu(k,1167)
         lu(k,1199) = lu(k,1199) - lu(k,787) * lu(k,1187)
         lu(k,1200) = lu(k,1200) - lu(k,788) * lu(k,1187)
         lu(k,1203) = lu(k,1203) - lu(k,789) * lu(k,1187)
         lu(k,1262) = lu(k,1262) - lu(k,787) * lu(k,1247)
         lu(k,1263) = lu(k,1263) - lu(k,788) * lu(k,1247)
         lu(k,1266) = lu(k,1266) - lu(k,789) * lu(k,1247)
         lu(k,1294) = lu(k,1294) - lu(k,787) * lu(k,1277)
         lu(k,1295) = lu(k,1295) - lu(k,788) * lu(k,1277)
         lu(k,1298) = lu(k,1298) - lu(k,789) * lu(k,1277)
         lu(k,1315) = lu(k,1315) - lu(k,787) * lu(k,1305)
         lu(k,1316) = lu(k,1316) - lu(k,788) * lu(k,1305)
         lu(k,1319) = lu(k,1319) - lu(k,789) * lu(k,1305)
         lu(k,1689) = lu(k,1689) - lu(k,787) * lu(k,1642)
         lu(k,1691) = lu(k,1691) - lu(k,788) * lu(k,1642)
         lu(k,1694) = lu(k,1694) - lu(k,789) * lu(k,1642)
         lu(k,2062) = lu(k,2062) - lu(k,787) * lu(k,2031)
         lu(k,2064) = lu(k,2064) - lu(k,788) * lu(k,2031)
         lu(k,2067) = lu(k,2067) - lu(k,789) * lu(k,2031)
                                                                        
         lu(k,791) = 1._r8 / lu(k,791)
         lu(k,792) = lu(k,792) * lu(k,791)
         lu(k,793) = lu(k,793) * lu(k,791)
         lu(k,794) = lu(k,794) * lu(k,791)
         lu(k,795) = lu(k,795) * lu(k,791)
         lu(k,796) = lu(k,796) * lu(k,791)
         lu(k,797) = lu(k,797) * lu(k,791)
         lu(k,798) = lu(k,798) * lu(k,791)
         lu(k,799) = lu(k,799) * lu(k,791)
         lu(k,800) = lu(k,800) * lu(k,791)
         lu(k,1656) = lu(k,1656) - lu(k,792) * lu(k,1643)
         lu(k,1666) = lu(k,1666) - lu(k,793) * lu(k,1643)
         lu(k,1689) = lu(k,1689) - lu(k,794) * lu(k,1643)
         lu(k,1691) = lu(k,1691) - lu(k,795) * lu(k,1643)
         lu(k,1693) = lu(k,1693) - lu(k,796) * lu(k,1643)
         lu(k,1694) = lu(k,1694) - lu(k,797) * lu(k,1643)
         lu(k,1697) = lu(k,1697) - lu(k,798) * lu(k,1643)
         lu(k,1700) = lu(k,1700) - lu(k,799) * lu(k,1643)
         lu(k,1703) = lu(k,1703) - lu(k,800) * lu(k,1643)
         lu(k,1810) = lu(k,1810) - lu(k,792) * lu(k,1803)
         lu(k,1817) = lu(k,1817) - lu(k,793) * lu(k,1803)
         lu(k,1838) = lu(k,1838) - lu(k,794) * lu(k,1803)
         lu(k,1840) = lu(k,1840) - lu(k,795) * lu(k,1803)
         lu(k,1842) = lu(k,1842) - lu(k,796) * lu(k,1803)
         lu(k,1843) = lu(k,1843) - lu(k,797) * lu(k,1803)
         lu(k,1846) = lu(k,1846) - lu(k,798) * lu(k,1803)
         lu(k,1849) = lu(k,1849) - lu(k,799) * lu(k,1803)
         lu(k,1852) = lu(k,1852) - lu(k,800) * lu(k,1803)
         lu(k,1918) = lu(k,1918) - lu(k,792) * lu(k,1909)
         lu(k,1924) = lu(k,1924) - lu(k,793) * lu(k,1909)
         lu(k,1945) = lu(k,1945) - lu(k,794) * lu(k,1909)
         lu(k,1947) = lu(k,1947) - lu(k,795) * lu(k,1909)
         lu(k,1949) = lu(k,1949) - lu(k,796) * lu(k,1909)
         lu(k,1950) = lu(k,1950) - lu(k,797) * lu(k,1909)
         lu(k,1953) = lu(k,1953) - lu(k,798) * lu(k,1909)
         lu(k,1956) = lu(k,1956) - lu(k,799) * lu(k,1909)
         lu(k,1959) = lu(k,1959) - lu(k,800) * lu(k,1909)
         lu(k,2037) = lu(k,2037) - lu(k,792) * lu(k,2032)
         lu(k,2043) = lu(k,2043) - lu(k,793) * lu(k,2032)
         lu(k,2062) = lu(k,2062) - lu(k,794) * lu(k,2032)
         lu(k,2064) = lu(k,2064) - lu(k,795) * lu(k,2032)
         lu(k,2066) = lu(k,2066) - lu(k,796) * lu(k,2032)
         lu(k,2067) = lu(k,2067) - lu(k,797) * lu(k,2032)
         lu(k,2070) = lu(k,2070) - lu(k,798) * lu(k,2032)
         lu(k,2073) = lu(k,2073) - lu(k,799) * lu(k,2032)
         lu(k,2076) = lu(k,2076) - lu(k,800) * lu(k,2032)
                                                                        
         lu(k,801) = 1._r8 / lu(k,801)
         lu(k,802) = lu(k,802) * lu(k,801)
         lu(k,803) = lu(k,803) * lu(k,801)
         lu(k,930) = - lu(k,802) * lu(k,928)
         lu(k,933) = - lu(k,803) * lu(k,928)
         lu(k,954) = lu(k,954) - lu(k,802) * lu(k,943)
         lu(k,967) = - lu(k,803) * lu(k,943)
         lu(k,980) = lu(k,980) - lu(k,802) * lu(k,978)
         lu(k,983) = - lu(k,803) * lu(k,978)
         lu(k,1003) = lu(k,1003) - lu(k,802) * lu(k,992)
         lu(k,1017) = - lu(k,803) * lu(k,992)
         lu(k,1026) = lu(k,1026) - lu(k,802) * lu(k,1021)
         lu(k,1039) = - lu(k,803) * lu(k,1021)
         lu(k,1060) = lu(k,1060) - lu(k,802) * lu(k,1057)
         lu(k,1067) = - lu(k,803) * lu(k,1057)
         lu(k,1097) = lu(k,1097) - lu(k,802) * lu(k,1094)
         lu(k,1101) = - lu(k,803) * lu(k,1094)
         lu(k,1103) = lu(k,1103) - lu(k,802) * lu(k,1102)
         lu(k,1106) = - lu(k,803) * lu(k,1102)
         lu(k,1114) = lu(k,1114) - lu(k,802) * lu(k,1110)
         lu(k,1125) = - lu(k,803) * lu(k,1110)
         lu(k,1171) = lu(k,1171) - lu(k,802) * lu(k,1168)
         lu(k,1183) = - lu(k,803) * lu(k,1168)
         lu(k,1250) = - lu(k,802) * lu(k,1248)
         lu(k,1270) = - lu(k,803) * lu(k,1248)
         lu(k,1328) = lu(k,1328) - lu(k,802) * lu(k,1324)
         lu(k,1345) = - lu(k,803) * lu(k,1324)
         lu(k,1351) = - lu(k,802) * lu(k,1350)
         lu(k,1365) = - lu(k,803) * lu(k,1350)
         lu(k,1375) = lu(k,1375) - lu(k,802) * lu(k,1371)
         lu(k,1397) = - lu(k,803) * lu(k,1371)
         lu(k,1428) = lu(k,1428) - lu(k,802) * lu(k,1426)
         lu(k,1441) = lu(k,1441) - lu(k,803) * lu(k,1426)
         lu(k,1668) = lu(k,1668) - lu(k,802) * lu(k,1644)
         lu(k,1702) = lu(k,1702) - lu(k,803) * lu(k,1644)
         lu(k,1819) = lu(k,1819) - lu(k,802) * lu(k,1804)
         lu(k,1851) = lu(k,1851) - lu(k,803) * lu(k,1804)
         lu(k,1926) = lu(k,1926) - lu(k,802) * lu(k,1910)
         lu(k,1958) = lu(k,1958) - lu(k,803) * lu(k,1910)
         lu(k,2045) = lu(k,2045) - lu(k,802) * lu(k,2033)
         lu(k,2075) = - lu(k,803) * lu(k,2033)
         lu(k,2105) = lu(k,2105) - lu(k,802) * lu(k,2090)
         lu(k,2136) = lu(k,2136) - lu(k,803) * lu(k,2090)
         lu(k,2180) = lu(k,2180) - lu(k,802) * lu(k,2175)
         lu(k,2203) = lu(k,2203) - lu(k,803) * lu(k,2175)
                                                                        
         lu(k,804) = 1._r8 / lu(k,804)
         lu(k,805) = lu(k,805) * lu(k,804)
         lu(k,806) = lu(k,806) * lu(k,804)
         lu(k,807) = lu(k,807) * lu(k,804)
         lu(k,808) = lu(k,808) * lu(k,804)
         lu(k,809) = lu(k,809) * lu(k,804)
         lu(k,810) = lu(k,810) * lu(k,804)
         lu(k,811) = lu(k,811) * lu(k,804)
         lu(k,1415) = lu(k,1415) - lu(k,805) * lu(k,1413)
         lu(k,1416) = - lu(k,806) * lu(k,1413)
         lu(k,1418) = - lu(k,807) * lu(k,1413)
         lu(k,1419) = - lu(k,808) * lu(k,1413)
         lu(k,1422) = lu(k,1422) - lu(k,809) * lu(k,1413)
         lu(k,1423) = - lu(k,810) * lu(k,1413)
         lu(k,1424) = - lu(k,811) * lu(k,1413)
         lu(k,1481) = lu(k,1481) - lu(k,805) * lu(k,1477)
         lu(k,1484) = lu(k,1484) - lu(k,806) * lu(k,1477)
         lu(k,1486) = - lu(k,807) * lu(k,1477)
         lu(k,1487) = lu(k,1487) - lu(k,808) * lu(k,1477)
         lu(k,1496) = - lu(k,809) * lu(k,1477)
         lu(k,1497) = lu(k,1497) - lu(k,810) * lu(k,1477)
         lu(k,1498) = lu(k,1498) - lu(k,811) * lu(k,1477)
         lu(k,1521) = lu(k,1521) - lu(k,805) * lu(k,1518)
         lu(k,1524) = lu(k,1524) - lu(k,806) * lu(k,1518)
         lu(k,1526) = lu(k,1526) - lu(k,807) * lu(k,1518)
         lu(k,1527) = lu(k,1527) - lu(k,808) * lu(k,1518)
         lu(k,1537) = lu(k,1537) - lu(k,809) * lu(k,1518)
         lu(k,1538) = lu(k,1538) - lu(k,810) * lu(k,1518)
         lu(k,1539) = lu(k,1539) - lu(k,811) * lu(k,1518)
         lu(k,1685) = lu(k,1685) - lu(k,805) * lu(k,1645)
         lu(k,1688) = lu(k,1688) - lu(k,806) * lu(k,1645)
         lu(k,1690) = lu(k,1690) - lu(k,807) * lu(k,1645)
         lu(k,1691) = lu(k,1691) - lu(k,808) * lu(k,1645)
         lu(k,1701) = lu(k,1701) - lu(k,809) * lu(k,1645)
         lu(k,1702) = lu(k,1702) - lu(k,810) * lu(k,1645)
         lu(k,1703) = lu(k,1703) - lu(k,811) * lu(k,1645)
         lu(k,1941) = lu(k,1941) - lu(k,805) * lu(k,1911)
         lu(k,1944) = lu(k,1944) - lu(k,806) * lu(k,1911)
         lu(k,1946) = - lu(k,807) * lu(k,1911)
         lu(k,1947) = lu(k,1947) - lu(k,808) * lu(k,1911)
         lu(k,1957) = lu(k,1957) - lu(k,809) * lu(k,1911)
         lu(k,1958) = lu(k,1958) - lu(k,810) * lu(k,1911)
         lu(k,1959) = lu(k,1959) - lu(k,811) * lu(k,1911)
         lu(k,2241) = lu(k,2241) - lu(k,805) * lu(k,2234)
         lu(k,2244) = lu(k,2244) - lu(k,806) * lu(k,2234)
         lu(k,2246) = - lu(k,807) * lu(k,2234)
         lu(k,2247) = lu(k,2247) - lu(k,808) * lu(k,2234)
         lu(k,2257) = lu(k,2257) - lu(k,809) * lu(k,2234)
         lu(k,2258) = lu(k,2258) - lu(k,810) * lu(k,2234)
         lu(k,2259) = lu(k,2259) - lu(k,811) * lu(k,2234)
                                                                        
      end do
                                                                        
      end subroutine lu_fac17
                                                                        
      subroutine lu_fac18( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,815) = 1._r8 / lu(k,815)
         lu(k,816) = lu(k,816) * lu(k,815)
         lu(k,817) = lu(k,817) * lu(k,815)
         lu(k,818) = lu(k,818) * lu(k,815)
         lu(k,819) = lu(k,819) * lu(k,815)
         lu(k,820) = lu(k,820) * lu(k,815)
         lu(k,821) = lu(k,821) * lu(k,815)
         lu(k,822) = lu(k,822) * lu(k,815)
         lu(k,901) = lu(k,901) - lu(k,816) * lu(k,898)
         lu(k,903) = lu(k,903) - lu(k,817) * lu(k,898)
         lu(k,904) = lu(k,904) - lu(k,818) * lu(k,898)
         lu(k,905) = lu(k,905) - lu(k,819) * lu(k,898)
         lu(k,906) = lu(k,906) - lu(k,820) * lu(k,898)
         lu(k,907) = - lu(k,821) * lu(k,898)
         lu(k,908) = lu(k,908) - lu(k,822) * lu(k,898)
         lu(k,1666) = lu(k,1666) - lu(k,816) * lu(k,1646)
         lu(k,1689) = lu(k,1689) - lu(k,817) * lu(k,1646)
         lu(k,1691) = lu(k,1691) - lu(k,818) * lu(k,1646)
         lu(k,1693) = lu(k,1693) - lu(k,819) * lu(k,1646)
         lu(k,1694) = lu(k,1694) - lu(k,820) * lu(k,1646)
         lu(k,1697) = lu(k,1697) - lu(k,821) * lu(k,1646)
         lu(k,1700) = lu(k,1700) - lu(k,822) * lu(k,1646)
         lu(k,1817) = lu(k,1817) - lu(k,816) * lu(k,1805)
         lu(k,1838) = lu(k,1838) - lu(k,817) * lu(k,1805)
         lu(k,1840) = lu(k,1840) - lu(k,818) * lu(k,1805)
         lu(k,1842) = lu(k,1842) - lu(k,819) * lu(k,1805)
         lu(k,1843) = lu(k,1843) - lu(k,820) * lu(k,1805)
         lu(k,1846) = lu(k,1846) - lu(k,821) * lu(k,1805)
         lu(k,1849) = lu(k,1849) - lu(k,822) * lu(k,1805)
         lu(k,1924) = lu(k,1924) - lu(k,816) * lu(k,1912)
         lu(k,1945) = lu(k,1945) - lu(k,817) * lu(k,1912)
         lu(k,1947) = lu(k,1947) - lu(k,818) * lu(k,1912)
         lu(k,1949) = lu(k,1949) - lu(k,819) * lu(k,1912)
         lu(k,1950) = lu(k,1950) - lu(k,820) * lu(k,1912)
         lu(k,1953) = lu(k,1953) - lu(k,821) * lu(k,1912)
         lu(k,1956) = lu(k,1956) - lu(k,822) * lu(k,1912)
         lu(k,2001) = - lu(k,816) * lu(k,1996)
         lu(k,2010) = lu(k,2010) - lu(k,817) * lu(k,1996)
         lu(k,2012) = lu(k,2012) - lu(k,818) * lu(k,1996)
         lu(k,2014) = lu(k,2014) - lu(k,819) * lu(k,1996)
         lu(k,2015) = lu(k,2015) - lu(k,820) * lu(k,1996)
         lu(k,2018) = lu(k,2018) - lu(k,821) * lu(k,1996)
         lu(k,2021) = lu(k,2021) - lu(k,822) * lu(k,1996)
         lu(k,2043) = lu(k,2043) - lu(k,816) * lu(k,2034)
         lu(k,2062) = lu(k,2062) - lu(k,817) * lu(k,2034)
         lu(k,2064) = lu(k,2064) - lu(k,818) * lu(k,2034)
         lu(k,2066) = lu(k,2066) - lu(k,819) * lu(k,2034)
         lu(k,2067) = lu(k,2067) - lu(k,820) * lu(k,2034)
         lu(k,2070) = lu(k,2070) - lu(k,821) * lu(k,2034)
         lu(k,2073) = lu(k,2073) - lu(k,822) * lu(k,2034)
                                                                        
         lu(k,824) = 1._r8 / lu(k,824)
         lu(k,825) = lu(k,825) * lu(k,824)
         lu(k,826) = lu(k,826) * lu(k,824)
         lu(k,827) = lu(k,827) * lu(k,824)
         lu(k,828) = lu(k,828) * lu(k,824)
         lu(k,829) = lu(k,829) * lu(k,824)
         lu(k,830) = lu(k,830) * lu(k,824)
         lu(k,886) = lu(k,886) - lu(k,825) * lu(k,883)
         lu(k,888) = lu(k,888) - lu(k,826) * lu(k,883)
         lu(k,889) = lu(k,889) - lu(k,827) * lu(k,883)
         lu(k,890) = lu(k,890) - lu(k,828) * lu(k,883)
         lu(k,892) = lu(k,892) - lu(k,829) * lu(k,883)
         lu(k,893) = - lu(k,830) * lu(k,883)
         lu(k,1691) = lu(k,1691) - lu(k,825) * lu(k,1647)
         lu(k,1695) = lu(k,1695) - lu(k,826) * lu(k,1647)
         lu(k,1696) = lu(k,1696) - lu(k,827) * lu(k,1647)
         lu(k,1699) = lu(k,1699) - lu(k,828) * lu(k,1647)
         lu(k,1702) = lu(k,1702) - lu(k,829) * lu(k,1647)
         lu(k,1703) = lu(k,1703) - lu(k,830) * lu(k,1647)
         lu(k,1947) = lu(k,1947) - lu(k,825) * lu(k,1913)
         lu(k,1951) = lu(k,1951) - lu(k,826) * lu(k,1913)
         lu(k,1952) = lu(k,1952) - lu(k,827) * lu(k,1913)
         lu(k,1955) = lu(k,1955) - lu(k,828) * lu(k,1913)
         lu(k,1958) = lu(k,1958) - lu(k,829) * lu(k,1913)
         lu(k,1959) = lu(k,1959) - lu(k,830) * lu(k,1913)
         lu(k,1973) = lu(k,1973) - lu(k,825) * lu(k,1964)
         lu(k,1977) = lu(k,1977) - lu(k,826) * lu(k,1964)
         lu(k,1978) = lu(k,1978) - lu(k,827) * lu(k,1964)
         lu(k,1981) = lu(k,1981) - lu(k,828) * lu(k,1964)
         lu(k,1984) = lu(k,1984) - lu(k,829) * lu(k,1964)
         lu(k,1985) = - lu(k,830) * lu(k,1964)
         lu(k,2012) = lu(k,2012) - lu(k,825) * lu(k,1997)
         lu(k,2016) = lu(k,2016) - lu(k,826) * lu(k,1997)
         lu(k,2017) = lu(k,2017) - lu(k,827) * lu(k,1997)
         lu(k,2020) = lu(k,2020) - lu(k,828) * lu(k,1997)
         lu(k,2023) = lu(k,2023) - lu(k,829) * lu(k,1997)
         lu(k,2024) = lu(k,2024) - lu(k,830) * lu(k,1997)
         lu(k,2148) = lu(k,2148) - lu(k,825) * lu(k,2141)
         lu(k,2152) = lu(k,2152) - lu(k,826) * lu(k,2141)
         lu(k,2153) = lu(k,2153) - lu(k,827) * lu(k,2141)
         lu(k,2156) = lu(k,2156) - lu(k,828) * lu(k,2141)
         lu(k,2159) = lu(k,2159) - lu(k,829) * lu(k,2141)
         lu(k,2160) = lu(k,2160) - lu(k,830) * lu(k,2141)
         lu(k,2247) = lu(k,2247) - lu(k,825) * lu(k,2235)
         lu(k,2251) = lu(k,2251) - lu(k,826) * lu(k,2235)
         lu(k,2252) = lu(k,2252) - lu(k,827) * lu(k,2235)
         lu(k,2255) = lu(k,2255) - lu(k,828) * lu(k,2235)
         lu(k,2258) = lu(k,2258) - lu(k,829) * lu(k,2235)
         lu(k,2259) = lu(k,2259) - lu(k,830) * lu(k,2235)
                                                                        
         lu(k,832) = 1._r8 / lu(k,832)
         lu(k,833) = lu(k,833) * lu(k,832)
         lu(k,834) = lu(k,834) * lu(k,832)
         lu(k,835) = lu(k,835) * lu(k,832)
         lu(k,836) = lu(k,836) * lu(k,832)
         lu(k,837) = lu(k,837) * lu(k,832)
         lu(k,1232) = lu(k,1232) - lu(k,833) * lu(k,1231)
         lu(k,1235) = lu(k,1235) - lu(k,834) * lu(k,1231)
         lu(k,1237) = - lu(k,835) * lu(k,1231)
         lu(k,1243) = lu(k,1243) - lu(k,836) * lu(k,1231)
         lu(k,1244) = - lu(k,837) * lu(k,1231)
         lu(k,1677) = lu(k,1677) - lu(k,833) * lu(k,1648)
         lu(k,1691) = lu(k,1691) - lu(k,834) * lu(k,1648)
         lu(k,1694) = lu(k,1694) - lu(k,835) * lu(k,1648)
         lu(k,1702) = lu(k,1702) - lu(k,836) * lu(k,1648)
         lu(k,1703) = lu(k,1703) - lu(k,837) * lu(k,1648)
         lu(k,1735) = - lu(k,833) * lu(k,1711)
         lu(k,1748) = lu(k,1748) - lu(k,834) * lu(k,1711)
         lu(k,1751) = lu(k,1751) - lu(k,835) * lu(k,1711)
         lu(k,1759) = lu(k,1759) - lu(k,836) * lu(k,1711)
         lu(k,1760) = - lu(k,837) * lu(k,1711)
         lu(k,1967) = lu(k,1967) - lu(k,833) * lu(k,1965)
         lu(k,1973) = lu(k,1973) - lu(k,834) * lu(k,1965)
         lu(k,1976) = lu(k,1976) - lu(k,835) * lu(k,1965)
         lu(k,1984) = lu(k,1984) - lu(k,836) * lu(k,1965)
         lu(k,1985) = lu(k,1985) - lu(k,837) * lu(k,1965)
         lu(k,2112) = lu(k,2112) - lu(k,833) * lu(k,2091)
         lu(k,2125) = lu(k,2125) - lu(k,834) * lu(k,2091)
         lu(k,2128) = lu(k,2128) - lu(k,835) * lu(k,2091)
         lu(k,2136) = lu(k,2136) - lu(k,836) * lu(k,2091)
         lu(k,2137) = lu(k,2137) - lu(k,837) * lu(k,2091)
         lu(k,2183) = lu(k,2183) - lu(k,833) * lu(k,2176)
         lu(k,2192) = lu(k,2192) - lu(k,834) * lu(k,2176)
         lu(k,2195) = lu(k,2195) - lu(k,835) * lu(k,2176)
         lu(k,2203) = lu(k,2203) - lu(k,836) * lu(k,2176)
         lu(k,2204) = lu(k,2204) - lu(k,837) * lu(k,2176)
         lu(k,2210) = lu(k,2210) - lu(k,833) * lu(k,2209)
         lu(k,2216) = lu(k,2216) - lu(k,834) * lu(k,2209)
         lu(k,2219) = lu(k,2219) - lu(k,835) * lu(k,2209)
         lu(k,2227) = lu(k,2227) - lu(k,836) * lu(k,2209)
         lu(k,2228) = lu(k,2228) - lu(k,837) * lu(k,2209)
         lu(k,2239) = lu(k,2239) - lu(k,833) * lu(k,2236)
         lu(k,2247) = lu(k,2247) - lu(k,834) * lu(k,2236)
         lu(k,2250) = lu(k,2250) - lu(k,835) * lu(k,2236)
         lu(k,2258) = lu(k,2258) - lu(k,836) * lu(k,2236)
         lu(k,2259) = lu(k,2259) - lu(k,837) * lu(k,2236)
         lu(k,2265) = - lu(k,833) * lu(k,2263)
         lu(k,2273) = lu(k,2273) - lu(k,834) * lu(k,2263)
         lu(k,2276) = - lu(k,835) * lu(k,2263)
         lu(k,2284) = lu(k,2284) - lu(k,836) * lu(k,2263)
         lu(k,2285) = lu(k,2285) - lu(k,837) * lu(k,2263)
                                                                        
         lu(k,841) = 1._r8 / lu(k,841)
         lu(k,842) = lu(k,842) * lu(k,841)
         lu(k,843) = lu(k,843) * lu(k,841)
         lu(k,844) = lu(k,844) * lu(k,841)
         lu(k,845) = lu(k,845) * lu(k,841)
         lu(k,846) = lu(k,846) * lu(k,841)
         lu(k,847) = lu(k,847) * lu(k,841)
         lu(k,848) = lu(k,848) * lu(k,841)
         lu(k,849) = lu(k,849) * lu(k,841)
         lu(k,850) = lu(k,850) * lu(k,841)
         lu(k,851) = lu(k,851) * lu(k,841)
         lu(k,852) = lu(k,852) * lu(k,841)
         lu(k,853) = lu(k,853) * lu(k,841)
         lu(k,854) = lu(k,854) * lu(k,841)
         lu(k,855) = lu(k,855) * lu(k,841)
         lu(k,856) = lu(k,856) * lu(k,841)
         lu(k,1657) = lu(k,1657) - lu(k,842) * lu(k,1649)
         lu(k,1662) = lu(k,1662) - lu(k,843) * lu(k,1649)
         lu(k,1668) = lu(k,1668) - lu(k,844) * lu(k,1649)
         lu(k,1674) = - lu(k,845) * lu(k,1649)
         lu(k,1675) = lu(k,1675) - lu(k,846) * lu(k,1649)
         lu(k,1678) = lu(k,1678) - lu(k,847) * lu(k,1649)
         lu(k,1679) = lu(k,1679) - lu(k,848) * lu(k,1649)
         lu(k,1681) = lu(k,1681) - lu(k,849) * lu(k,1649)
         lu(k,1683) = lu(k,1683) - lu(k,850) * lu(k,1649)
         lu(k,1689) = lu(k,1689) - lu(k,851) * lu(k,1649)
         lu(k,1691) = lu(k,1691) - lu(k,852) * lu(k,1649)
         lu(k,1692) = lu(k,1692) - lu(k,853) * lu(k,1649)
         lu(k,1694) = lu(k,1694) - lu(k,854) * lu(k,1649)
         lu(k,1697) = lu(k,1697) - lu(k,855) * lu(k,1649)
         lu(k,1698) = lu(k,1698) - lu(k,856) * lu(k,1649)
         lu(k,1715) = - lu(k,842) * lu(k,1712)
         lu(k,1720) = lu(k,1720) - lu(k,843) * lu(k,1712)
         lu(k,1726) = lu(k,1726) - lu(k,844) * lu(k,1712)
         lu(k,1732) = lu(k,1732) - lu(k,845) * lu(k,1712)
         lu(k,1733) = lu(k,1733) - lu(k,846) * lu(k,1712)
         lu(k,1736) = lu(k,1736) - lu(k,847) * lu(k,1712)
         lu(k,1737) = lu(k,1737) - lu(k,848) * lu(k,1712)
         lu(k,1739) = lu(k,1739) - lu(k,849) * lu(k,1712)
         lu(k,1741) = lu(k,1741) - lu(k,850) * lu(k,1712)
         lu(k,1746) = lu(k,1746) - lu(k,851) * lu(k,1712)
         lu(k,1748) = lu(k,1748) - lu(k,852) * lu(k,1712)
         lu(k,1749) = lu(k,1749) - lu(k,853) * lu(k,1712)
         lu(k,1751) = lu(k,1751) - lu(k,854) * lu(k,1712)
         lu(k,1754) = - lu(k,855) * lu(k,1712)
         lu(k,1755) = - lu(k,856) * lu(k,1712)
         lu(k,2095) = lu(k,2095) - lu(k,842) * lu(k,2092)
         lu(k,2100) = lu(k,2100) - lu(k,843) * lu(k,2092)
         lu(k,2105) = lu(k,2105) - lu(k,844) * lu(k,2092)
         lu(k,2109) = - lu(k,845) * lu(k,2092)
         lu(k,2110) = lu(k,2110) - lu(k,846) * lu(k,2092)
         lu(k,2113) = - lu(k,847) * lu(k,2092)
         lu(k,2114) = - lu(k,848) * lu(k,2092)
         lu(k,2116) = lu(k,2116) - lu(k,849) * lu(k,2092)
         lu(k,2118) = lu(k,2118) - lu(k,850) * lu(k,2092)
         lu(k,2123) = lu(k,2123) - lu(k,851) * lu(k,2092)
         lu(k,2125) = lu(k,2125) - lu(k,852) * lu(k,2092)
         lu(k,2126) = lu(k,2126) - lu(k,853) * lu(k,2092)
         lu(k,2128) = lu(k,2128) - lu(k,854) * lu(k,2092)
         lu(k,2131) = lu(k,2131) - lu(k,855) * lu(k,2092)
         lu(k,2132) = lu(k,2132) - lu(k,856) * lu(k,2092)
                                                                        
      end do
                                                                        
      end subroutine lu_fac18
                                                                        
      subroutine lu_fac19( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,857) = 1._r8 / lu(k,857)
         lu(k,858) = lu(k,858) * lu(k,857)
         lu(k,859) = lu(k,859) * lu(k,857)
         lu(k,860) = lu(k,860) * lu(k,857)
         lu(k,861) = lu(k,861) * lu(k,857)
         lu(k,862) = lu(k,862) * lu(k,857)
         lu(k,1028) = - lu(k,858) * lu(k,1022)
         lu(k,1029) = - lu(k,859) * lu(k,1022)
         lu(k,1031) = lu(k,1031) - lu(k,860) * lu(k,1022)
         lu(k,1032) = lu(k,1032) - lu(k,861) * lu(k,1022)
         lu(k,1038) = lu(k,1038) - lu(k,862) * lu(k,1022)
         lu(k,1075) = - lu(k,858) * lu(k,1070)
         lu(k,1076) = - lu(k,859) * lu(k,1070)
         lu(k,1078) = - lu(k,860) * lu(k,1070)
         lu(k,1079) = lu(k,1079) - lu(k,861) * lu(k,1070)
         lu(k,1083) = lu(k,1083) - lu(k,862) * lu(k,1070)
         lu(k,1195) = - lu(k,858) * lu(k,1188)
         lu(k,1197) = lu(k,1197) - lu(k,859) * lu(k,1188)
         lu(k,1199) = lu(k,1199) - lu(k,860) * lu(k,1188)
         lu(k,1200) = lu(k,1200) - lu(k,861) * lu(k,1188)
         lu(k,1205) = lu(k,1205) - lu(k,862) * lu(k,1188)
         lu(k,1287) = lu(k,1287) - lu(k,858) * lu(k,1278)
         lu(k,1292) = lu(k,1292) - lu(k,859) * lu(k,1278)
         lu(k,1294) = lu(k,1294) - lu(k,860) * lu(k,1278)
         lu(k,1295) = lu(k,1295) - lu(k,861) * lu(k,1278)
         lu(k,1301) = lu(k,1301) - lu(k,862) * lu(k,1278)
         lu(k,1676) = lu(k,1676) - lu(k,858) * lu(k,1650)
         lu(k,1683) = lu(k,1683) - lu(k,859) * lu(k,1650)
         lu(k,1689) = lu(k,1689) - lu(k,860) * lu(k,1650)
         lu(k,1691) = lu(k,1691) - lu(k,861) * lu(k,1650)
         lu(k,1700) = lu(k,1700) - lu(k,862) * lu(k,1650)
         lu(k,1734) = lu(k,1734) - lu(k,858) * lu(k,1713)
         lu(k,1741) = lu(k,1741) - lu(k,859) * lu(k,1713)
         lu(k,1746) = lu(k,1746) - lu(k,860) * lu(k,1713)
         lu(k,1748) = lu(k,1748) - lu(k,861) * lu(k,1713)
         lu(k,1757) = lu(k,1757) - lu(k,862) * lu(k,1713)
         lu(k,1827) = lu(k,1827) - lu(k,858) * lu(k,1806)
         lu(k,1833) = lu(k,1833) - lu(k,859) * lu(k,1806)
         lu(k,1838) = lu(k,1838) - lu(k,860) * lu(k,1806)
         lu(k,1840) = lu(k,1840) - lu(k,861) * lu(k,1806)
         lu(k,1849) = lu(k,1849) - lu(k,862) * lu(k,1806)
         lu(k,1933) = lu(k,1933) - lu(k,858) * lu(k,1914)
         lu(k,1939) = lu(k,1939) - lu(k,859) * lu(k,1914)
         lu(k,1945) = lu(k,1945) - lu(k,860) * lu(k,1914)
         lu(k,1947) = lu(k,1947) - lu(k,861) * lu(k,1914)
         lu(k,1956) = lu(k,1956) - lu(k,862) * lu(k,1914)
         lu(k,2052) = lu(k,2052) - lu(k,858) * lu(k,2035)
         lu(k,2058) = lu(k,2058) - lu(k,859) * lu(k,2035)
         lu(k,2062) = lu(k,2062) - lu(k,860) * lu(k,2035)
         lu(k,2064) = lu(k,2064) - lu(k,861) * lu(k,2035)
         lu(k,2073) = lu(k,2073) - lu(k,862) * lu(k,2035)
                                                                        
         lu(k,864) = 1._r8 / lu(k,864)
         lu(k,865) = lu(k,865) * lu(k,864)
         lu(k,866) = lu(k,866) * lu(k,864)
         lu(k,867) = lu(k,867) * lu(k,864)
         lu(k,868) = lu(k,868) * lu(k,864)
         lu(k,869) = lu(k,869) * lu(k,864)
         lu(k,870) = lu(k,870) * lu(k,864)
         lu(k,871) = lu(k,871) * lu(k,864)
         lu(k,872) = lu(k,872) * lu(k,864)
         lu(k,1401) = lu(k,1401) - lu(k,865) * lu(k,1400)
         lu(k,1402) = - lu(k,866) * lu(k,1400)
         lu(k,1403) = - lu(k,867) * lu(k,1400)
         lu(k,1404) = lu(k,1404) - lu(k,868) * lu(k,1400)
         lu(k,1406) = lu(k,1406) - lu(k,869) * lu(k,1400)
         lu(k,1407) = - lu(k,870) * lu(k,1400)
         lu(k,1409) = - lu(k,871) * lu(k,1400)
         lu(k,1412) = lu(k,1412) - lu(k,872) * lu(k,1400)
         lu(k,1429) = lu(k,1429) - lu(k,865) * lu(k,1427)
         lu(k,1430) = lu(k,1430) - lu(k,866) * lu(k,1427)
         lu(k,1431) = - lu(k,867) * lu(k,1427)
         lu(k,1432) = lu(k,1432) - lu(k,868) * lu(k,1427)
         lu(k,1435) = lu(k,1435) - lu(k,869) * lu(k,1427)
         lu(k,1436) = - lu(k,870) * lu(k,1427)
         lu(k,1439) = lu(k,1439) - lu(k,871) * lu(k,1427)
         lu(k,1442) = lu(k,1442) - lu(k,872) * lu(k,1427)
         lu(k,1445) = - lu(k,865) * lu(k,1444)
         lu(k,1446) = - lu(k,866) * lu(k,1444)
         lu(k,1447) = lu(k,1447) - lu(k,867) * lu(k,1444)
         lu(k,1448) = lu(k,1448) - lu(k,868) * lu(k,1444)
         lu(k,1451) = lu(k,1451) - lu(k,869) * lu(k,1444)
         lu(k,1452) = lu(k,1452) - lu(k,870) * lu(k,1444)
         lu(k,1455) = - lu(k,871) * lu(k,1444)
         lu(k,1459) = lu(k,1459) - lu(k,872) * lu(k,1444)
         lu(k,1520) = lu(k,1520) - lu(k,865) * lu(k,1519)
         lu(k,1522) = lu(k,1522) - lu(k,866) * lu(k,1519)
         lu(k,1523) = - lu(k,867) * lu(k,1519)
         lu(k,1524) = lu(k,1524) - lu(k,868) * lu(k,1519)
         lu(k,1527) = lu(k,1527) - lu(k,869) * lu(k,1519)
         lu(k,1528) = - lu(k,870) * lu(k,1519)
         lu(k,1533) = lu(k,1533) - lu(k,871) * lu(k,1519)
         lu(k,1539) = lu(k,1539) - lu(k,872) * lu(k,1519)
         lu(k,1684) = lu(k,1684) - lu(k,865) * lu(k,1651)
         lu(k,1686) = lu(k,1686) - lu(k,866) * lu(k,1651)
         lu(k,1687) = lu(k,1687) - lu(k,867) * lu(k,1651)
         lu(k,1688) = lu(k,1688) - lu(k,868) * lu(k,1651)
         lu(k,1691) = lu(k,1691) - lu(k,869) * lu(k,1651)
         lu(k,1692) = lu(k,1692) - lu(k,870) * lu(k,1651)
         lu(k,1697) = lu(k,1697) - lu(k,871) * lu(k,1651)
         lu(k,1703) = lu(k,1703) - lu(k,872) * lu(k,1651)
         lu(k,2266) = lu(k,2266) - lu(k,865) * lu(k,2264)
         lu(k,2268) = - lu(k,866) * lu(k,2264)
         lu(k,2269) = - lu(k,867) * lu(k,2264)
         lu(k,2270) = lu(k,2270) - lu(k,868) * lu(k,2264)
         lu(k,2273) = lu(k,2273) - lu(k,869) * lu(k,2264)
         lu(k,2274) = - lu(k,870) * lu(k,2264)
         lu(k,2279) = - lu(k,871) * lu(k,2264)
         lu(k,2285) = lu(k,2285) - lu(k,872) * lu(k,2264)
                                                                        
         lu(k,873) = 1._r8 / lu(k,873)
         lu(k,874) = lu(k,874) * lu(k,873)
         lu(k,875) = lu(k,875) * lu(k,873)
         lu(k,876) = lu(k,876) * lu(k,873)
         lu(k,877) = lu(k,877) * lu(k,873)
         lu(k,878) = lu(k,878) * lu(k,873)
         lu(k,879) = lu(k,879) * lu(k,873)
         lu(k,880) = lu(k,880) * lu(k,873)
         lu(k,881) = lu(k,881) * lu(k,873)
         lu(k,1072) = lu(k,1072) - lu(k,874) * lu(k,1071)
         lu(k,1074) = lu(k,1074) - lu(k,875) * lu(k,1071)
         lu(k,1075) = lu(k,1075) - lu(k,876) * lu(k,1071)
         lu(k,1079) = lu(k,1079) - lu(k,877) * lu(k,1071)
         lu(k,1080) = - lu(k,878) * lu(k,1071)
         lu(k,1081) = lu(k,1081) - lu(k,879) * lu(k,1071)
         lu(k,1082) = - lu(k,880) * lu(k,1071)
         lu(k,1083) = lu(k,1083) - lu(k,881) * lu(k,1071)
         lu(k,1280) = lu(k,1280) - lu(k,874) * lu(k,1279)
         lu(k,1282) = lu(k,1282) - lu(k,875) * lu(k,1279)
         lu(k,1287) = lu(k,1287) - lu(k,876) * lu(k,1279)
         lu(k,1295) = lu(k,1295) - lu(k,877) * lu(k,1279)
         lu(k,1297) = lu(k,1297) - lu(k,878) * lu(k,1279)
         lu(k,1298) = lu(k,1298) - lu(k,879) * lu(k,1279)
         lu(k,1299) = lu(k,1299) - lu(k,880) * lu(k,1279)
         lu(k,1301) = lu(k,1301) - lu(k,881) * lu(k,1279)
         lu(k,1660) = lu(k,1660) - lu(k,874) * lu(k,1652)
         lu(k,1668) = lu(k,1668) - lu(k,875) * lu(k,1652)
         lu(k,1676) = lu(k,1676) - lu(k,876) * lu(k,1652)
         lu(k,1691) = lu(k,1691) - lu(k,877) * lu(k,1652)
         lu(k,1693) = lu(k,1693) - lu(k,878) * lu(k,1652)
         lu(k,1694) = lu(k,1694) - lu(k,879) * lu(k,1652)
         lu(k,1697) = lu(k,1697) - lu(k,880) * lu(k,1652)
         lu(k,1700) = lu(k,1700) - lu(k,881) * lu(k,1652)
         lu(k,1813) = lu(k,1813) - lu(k,874) * lu(k,1807)
         lu(k,1819) = lu(k,1819) - lu(k,875) * lu(k,1807)
         lu(k,1827) = lu(k,1827) - lu(k,876) * lu(k,1807)
         lu(k,1840) = lu(k,1840) - lu(k,877) * lu(k,1807)
         lu(k,1842) = lu(k,1842) - lu(k,878) * lu(k,1807)
         lu(k,1843) = lu(k,1843) - lu(k,879) * lu(k,1807)
         lu(k,1846) = lu(k,1846) - lu(k,880) * lu(k,1807)
         lu(k,1849) = lu(k,1849) - lu(k,881) * lu(k,1807)
         lu(k,1921) = lu(k,1921) - lu(k,874) * lu(k,1915)
         lu(k,1926) = lu(k,1926) - lu(k,875) * lu(k,1915)
         lu(k,1933) = lu(k,1933) - lu(k,876) * lu(k,1915)
         lu(k,1947) = lu(k,1947) - lu(k,877) * lu(k,1915)
         lu(k,1949) = lu(k,1949) - lu(k,878) * lu(k,1915)
         lu(k,1950) = lu(k,1950) - lu(k,879) * lu(k,1915)
         lu(k,1953) = lu(k,1953) - lu(k,880) * lu(k,1915)
         lu(k,1956) = lu(k,1956) - lu(k,881) * lu(k,1915)
         lu(k,2179) = lu(k,2179) - lu(k,874) * lu(k,2177)
         lu(k,2180) = lu(k,2180) - lu(k,875) * lu(k,2177)
         lu(k,2182) = lu(k,2182) - lu(k,876) * lu(k,2177)
         lu(k,2192) = lu(k,2192) - lu(k,877) * lu(k,2177)
         lu(k,2194) = lu(k,2194) - lu(k,878) * lu(k,2177)
         lu(k,2195) = lu(k,2195) - lu(k,879) * lu(k,2177)
         lu(k,2198) = lu(k,2198) - lu(k,880) * lu(k,2177)
         lu(k,2201) = lu(k,2201) - lu(k,881) * lu(k,2177)
                                                                        
         lu(k,884) = 1._r8 / lu(k,884)
         lu(k,885) = lu(k,885) * lu(k,884)
         lu(k,886) = lu(k,886) * lu(k,884)
         lu(k,887) = lu(k,887) * lu(k,884)
         lu(k,888) = lu(k,888) * lu(k,884)
         lu(k,889) = lu(k,889) * lu(k,884)
         lu(k,890) = lu(k,890) * lu(k,884)
         lu(k,891) = lu(k,891) * lu(k,884)
         lu(k,892) = lu(k,892) * lu(k,884)
         lu(k,893) = lu(k,893) * lu(k,884)
         lu(k,1687) = lu(k,1687) - lu(k,885) * lu(k,1653)
         lu(k,1691) = lu(k,1691) - lu(k,886) * lu(k,1653)
         lu(k,1692) = lu(k,1692) - lu(k,887) * lu(k,1653)
         lu(k,1695) = lu(k,1695) - lu(k,888) * lu(k,1653)
         lu(k,1696) = lu(k,1696) - lu(k,889) * lu(k,1653)
         lu(k,1699) = lu(k,1699) - lu(k,890) * lu(k,1653)
         lu(k,1700) = lu(k,1700) - lu(k,891) * lu(k,1653)
         lu(k,1702) = lu(k,1702) - lu(k,892) * lu(k,1653)
         lu(k,1703) = lu(k,1703) - lu(k,893) * lu(k,1653)
         lu(k,1969) = - lu(k,885) * lu(k,1966)
         lu(k,1973) = lu(k,1973) - lu(k,886) * lu(k,1966)
         lu(k,1974) = - lu(k,887) * lu(k,1966)
         lu(k,1977) = lu(k,1977) - lu(k,888) * lu(k,1966)
         lu(k,1978) = lu(k,1978) - lu(k,889) * lu(k,1966)
         lu(k,1981) = lu(k,1981) - lu(k,890) * lu(k,1966)
         lu(k,1982) = lu(k,1982) - lu(k,891) * lu(k,1966)
         lu(k,1984) = lu(k,1984) - lu(k,892) * lu(k,1966)
         lu(k,1985) = lu(k,1985) - lu(k,893) * lu(k,1966)
         lu(k,2008) = - lu(k,885) * lu(k,1998)
         lu(k,2012) = lu(k,2012) - lu(k,886) * lu(k,1998)
         lu(k,2013) = lu(k,2013) - lu(k,887) * lu(k,1998)
         lu(k,2016) = lu(k,2016) - lu(k,888) * lu(k,1998)
         lu(k,2017) = lu(k,2017) - lu(k,889) * lu(k,1998)
         lu(k,2020) = lu(k,2020) - lu(k,890) * lu(k,1998)
         lu(k,2021) = lu(k,2021) - lu(k,891) * lu(k,1998)
         lu(k,2023) = lu(k,2023) - lu(k,892) * lu(k,1998)
         lu(k,2024) = lu(k,2024) - lu(k,893) * lu(k,1998)
         lu(k,2144) = lu(k,2144) - lu(k,885) * lu(k,2142)
         lu(k,2148) = lu(k,2148) - lu(k,886) * lu(k,2142)
         lu(k,2149) = - lu(k,887) * lu(k,2142)
         lu(k,2152) = lu(k,2152) - lu(k,888) * lu(k,2142)
         lu(k,2153) = lu(k,2153) - lu(k,889) * lu(k,2142)
         lu(k,2156) = lu(k,2156) - lu(k,890) * lu(k,2142)
         lu(k,2157) = - lu(k,891) * lu(k,2142)
         lu(k,2159) = lu(k,2159) - lu(k,892) * lu(k,2142)
         lu(k,2160) = lu(k,2160) - lu(k,893) * lu(k,2142)
         lu(k,2188) = lu(k,2188) - lu(k,885) * lu(k,2178)
         lu(k,2192) = lu(k,2192) - lu(k,886) * lu(k,2178)
         lu(k,2193) = lu(k,2193) - lu(k,887) * lu(k,2178)
         lu(k,2196) = lu(k,2196) - lu(k,888) * lu(k,2178)
         lu(k,2197) = lu(k,2197) - lu(k,889) * lu(k,2178)
         lu(k,2200) = lu(k,2200) - lu(k,890) * lu(k,2178)
         lu(k,2201) = lu(k,2201) - lu(k,891) * lu(k,2178)
         lu(k,2203) = lu(k,2203) - lu(k,892) * lu(k,2178)
         lu(k,2204) = lu(k,2204) - lu(k,893) * lu(k,2178)
         lu(k,2243) = lu(k,2243) - lu(k,885) * lu(k,2237)
         lu(k,2247) = lu(k,2247) - lu(k,886) * lu(k,2237)
         lu(k,2248) = lu(k,2248) - lu(k,887) * lu(k,2237)
         lu(k,2251) = lu(k,2251) - lu(k,888) * lu(k,2237)
         lu(k,2252) = lu(k,2252) - lu(k,889) * lu(k,2237)
         lu(k,2255) = lu(k,2255) - lu(k,890) * lu(k,2237)
         lu(k,2256) = lu(k,2256) - lu(k,891) * lu(k,2237)
         lu(k,2258) = lu(k,2258) - lu(k,892) * lu(k,2237)
         lu(k,2259) = lu(k,2259) - lu(k,893) * lu(k,2237)
                                                                        
      end do
                                                                        
      end subroutine lu_fac19
                                                                        
      subroutine lu_fac20( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,899) = 1._r8 / lu(k,899)
         lu(k,900) = lu(k,900) * lu(k,899)
         lu(k,901) = lu(k,901) * lu(k,899)
         lu(k,902) = lu(k,902) * lu(k,899)
         lu(k,903) = lu(k,903) * lu(k,899)
         lu(k,904) = lu(k,904) * lu(k,899)
         lu(k,905) = lu(k,905) * lu(k,899)
         lu(k,906) = lu(k,906) * lu(k,899)
         lu(k,907) = lu(k,907) * lu(k,899)
         lu(k,908) = lu(k,908) * lu(k,899)
         lu(k,946) = lu(k,946) - lu(k,900) * lu(k,944)
         lu(k,953) = - lu(k,901) * lu(k,944)
         lu(k,957) = lu(k,957) - lu(k,902) * lu(k,944)
         lu(k,959) = lu(k,959) - lu(k,903) * lu(k,944)
         lu(k,960) = lu(k,960) - lu(k,904) * lu(k,944)
         lu(k,962) = - lu(k,905) * lu(k,944)
         lu(k,963) = lu(k,963) - lu(k,906) * lu(k,944)
         lu(k,964) = - lu(k,907) * lu(k,944)
         lu(k,966) = - lu(k,908) * lu(k,944)
         lu(k,995) = lu(k,995) - lu(k,900) * lu(k,993)
         lu(k,1002) = - lu(k,901) * lu(k,993)
         lu(k,1007) = lu(k,1007) - lu(k,902) * lu(k,993)
         lu(k,1009) = lu(k,1009) - lu(k,903) * lu(k,993)
         lu(k,1010) = lu(k,1010) - lu(k,904) * lu(k,993)
         lu(k,1012) = - lu(k,905) * lu(k,993)
         lu(k,1013) = lu(k,1013) - lu(k,906) * lu(k,993)
         lu(k,1014) = - lu(k,907) * lu(k,993)
         lu(k,1016) = - lu(k,908) * lu(k,993)
         lu(k,1656) = lu(k,1656) - lu(k,900) * lu(k,1654)
         lu(k,1666) = lu(k,1666) - lu(k,901) * lu(k,1654)
         lu(k,1683) = lu(k,1683) - lu(k,902) * lu(k,1654)
         lu(k,1689) = lu(k,1689) - lu(k,903) * lu(k,1654)
         lu(k,1691) = lu(k,1691) - lu(k,904) * lu(k,1654)
         lu(k,1693) = lu(k,1693) - lu(k,905) * lu(k,1654)
         lu(k,1694) = lu(k,1694) - lu(k,906) * lu(k,1654)
         lu(k,1697) = lu(k,1697) - lu(k,907) * lu(k,1654)
         lu(k,1700) = lu(k,1700) - lu(k,908) * lu(k,1654)
         lu(k,1810) = lu(k,1810) - lu(k,900) * lu(k,1808)
         lu(k,1817) = lu(k,1817) - lu(k,901) * lu(k,1808)
         lu(k,1833) = lu(k,1833) - lu(k,902) * lu(k,1808)
         lu(k,1838) = lu(k,1838) - lu(k,903) * lu(k,1808)
         lu(k,1840) = lu(k,1840) - lu(k,904) * lu(k,1808)
         lu(k,1842) = lu(k,1842) - lu(k,905) * lu(k,1808)
         lu(k,1843) = lu(k,1843) - lu(k,906) * lu(k,1808)
         lu(k,1846) = lu(k,1846) - lu(k,907) * lu(k,1808)
         lu(k,1849) = lu(k,1849) - lu(k,908) * lu(k,1808)
         lu(k,1918) = lu(k,1918) - lu(k,900) * lu(k,1916)
         lu(k,1924) = lu(k,1924) - lu(k,901) * lu(k,1916)
         lu(k,1939) = lu(k,1939) - lu(k,902) * lu(k,1916)
         lu(k,1945) = lu(k,1945) - lu(k,903) * lu(k,1916)
         lu(k,1947) = lu(k,1947) - lu(k,904) * lu(k,1916)
         lu(k,1949) = lu(k,1949) - lu(k,905) * lu(k,1916)
         lu(k,1950) = lu(k,1950) - lu(k,906) * lu(k,1916)
         lu(k,1953) = lu(k,1953) - lu(k,907) * lu(k,1916)
         lu(k,1956) = lu(k,1956) - lu(k,908) * lu(k,1916)
         lu(k,2094) = lu(k,2094) - lu(k,900) * lu(k,2093)
         lu(k,2103) = lu(k,2103) - lu(k,901) * lu(k,2093)
         lu(k,2118) = lu(k,2118) - lu(k,902) * lu(k,2093)
         lu(k,2123) = lu(k,2123) - lu(k,903) * lu(k,2093)
         lu(k,2125) = lu(k,2125) - lu(k,904) * lu(k,2093)
         lu(k,2127) = lu(k,2127) - lu(k,905) * lu(k,2093)
         lu(k,2128) = lu(k,2128) - lu(k,906) * lu(k,2093)
         lu(k,2131) = lu(k,2131) - lu(k,907) * lu(k,2093)
         lu(k,2134) = lu(k,2134) - lu(k,908) * lu(k,2093)
                                                                        
         lu(k,912) = 1._r8 / lu(k,912)
         lu(k,913) = lu(k,913) * lu(k,912)
         lu(k,914) = lu(k,914) * lu(k,912)
         lu(k,915) = lu(k,915) * lu(k,912)
         lu(k,916) = lu(k,916) * lu(k,912)
         lu(k,917) = lu(k,917) * lu(k,912)
         lu(k,918) = lu(k,918) * lu(k,912)
         lu(k,919) = lu(k,919) * lu(k,912)
         lu(k,920) = lu(k,920) * lu(k,912)
         lu(k,921) = lu(k,921) * lu(k,912)
         lu(k,946) = lu(k,946) - lu(k,913) * lu(k,945)
         lu(k,949) = lu(k,949) - lu(k,914) * lu(k,945)
         lu(k,958) = - lu(k,915) * lu(k,945)
         lu(k,959) = lu(k,959) - lu(k,916) * lu(k,945)
         lu(k,960) = lu(k,960) - lu(k,917) * lu(k,945)
         lu(k,962) = lu(k,962) - lu(k,918) * lu(k,945)
         lu(k,963) = lu(k,963) - lu(k,919) * lu(k,945)
         lu(k,964) = lu(k,964) - lu(k,920) * lu(k,945)
         lu(k,966) = lu(k,966) - lu(k,921) * lu(k,945)
         lu(k,995) = lu(k,995) - lu(k,913) * lu(k,994)
         lu(k,997) = lu(k,997) - lu(k,914) * lu(k,994)
         lu(k,1008) = - lu(k,915) * lu(k,994)
         lu(k,1009) = lu(k,1009) - lu(k,916) * lu(k,994)
         lu(k,1010) = lu(k,1010) - lu(k,917) * lu(k,994)
         lu(k,1012) = lu(k,1012) - lu(k,918) * lu(k,994)
         lu(k,1013) = lu(k,1013) - lu(k,919) * lu(k,994)
         lu(k,1014) = lu(k,1014) - lu(k,920) * lu(k,994)
         lu(k,1016) = lu(k,1016) - lu(k,921) * lu(k,994)
         lu(k,1656) = lu(k,1656) - lu(k,913) * lu(k,1655)
         lu(k,1659) = lu(k,1659) - lu(k,914) * lu(k,1655)
         lu(k,1687) = lu(k,1687) - lu(k,915) * lu(k,1655)
         lu(k,1689) = lu(k,1689) - lu(k,916) * lu(k,1655)
         lu(k,1691) = lu(k,1691) - lu(k,917) * lu(k,1655)
         lu(k,1693) = lu(k,1693) - lu(k,918) * lu(k,1655)
         lu(k,1694) = lu(k,1694) - lu(k,919) * lu(k,1655)
         lu(k,1697) = lu(k,1697) - lu(k,920) * lu(k,1655)
         lu(k,1700) = lu(k,1700) - lu(k,921) * lu(k,1655)
         lu(k,1810) = lu(k,1810) - lu(k,913) * lu(k,1809)
         lu(k,1812) = lu(k,1812) - lu(k,914) * lu(k,1809)
         lu(k,1836) = lu(k,1836) - lu(k,915) * lu(k,1809)
         lu(k,1838) = lu(k,1838) - lu(k,916) * lu(k,1809)
         lu(k,1840) = lu(k,1840) - lu(k,917) * lu(k,1809)
         lu(k,1842) = lu(k,1842) - lu(k,918) * lu(k,1809)
         lu(k,1843) = lu(k,1843) - lu(k,919) * lu(k,1809)
         lu(k,1846) = lu(k,1846) - lu(k,920) * lu(k,1809)
         lu(k,1849) = lu(k,1849) - lu(k,921) * lu(k,1809)
         lu(k,1918) = lu(k,1918) - lu(k,913) * lu(k,1917)
         lu(k,1920) = lu(k,1920) - lu(k,914) * lu(k,1917)
         lu(k,1943) = lu(k,1943) - lu(k,915) * lu(k,1917)
         lu(k,1945) = lu(k,1945) - lu(k,916) * lu(k,1917)
         lu(k,1947) = lu(k,1947) - lu(k,917) * lu(k,1917)
         lu(k,1949) = lu(k,1949) - lu(k,918) * lu(k,1917)
         lu(k,1950) = lu(k,1950) - lu(k,919) * lu(k,1917)
         lu(k,1953) = lu(k,1953) - lu(k,920) * lu(k,1917)
         lu(k,1956) = lu(k,1956) - lu(k,921) * lu(k,1917)
         lu(k,2037) = lu(k,2037) - lu(k,913) * lu(k,2036)
         lu(k,2038) = lu(k,2038) - lu(k,914) * lu(k,2036)
         lu(k,2060) = lu(k,2060) - lu(k,915) * lu(k,2036)
         lu(k,2062) = lu(k,2062) - lu(k,916) * lu(k,2036)
         lu(k,2064) = lu(k,2064) - lu(k,917) * lu(k,2036)
         lu(k,2066) = lu(k,2066) - lu(k,918) * lu(k,2036)
         lu(k,2067) = lu(k,2067) - lu(k,919) * lu(k,2036)
         lu(k,2070) = lu(k,2070) - lu(k,920) * lu(k,2036)
         lu(k,2073) = lu(k,2073) - lu(k,921) * lu(k,2036)
                                                                        
         lu(k,922) = 1._r8 / lu(k,922)
         lu(k,923) = lu(k,923) * lu(k,922)
         lu(k,924) = lu(k,924) * lu(k,922)
         lu(k,925) = lu(k,925) * lu(k,922)
         lu(k,926) = lu(k,926) * lu(k,922)
         lu(k,927) = lu(k,927) * lu(k,922)
         lu(k,955) = lu(k,955) - lu(k,923) * lu(k,946)
         lu(k,957) = lu(k,957) - lu(k,924) * lu(k,946)
         lu(k,960) = lu(k,960) - lu(k,925) * lu(k,946)
         lu(k,964) = lu(k,964) - lu(k,926) * lu(k,946)
         lu(k,968) = - lu(k,927) * lu(k,946)
         lu(k,1005) = lu(k,1005) - lu(k,923) * lu(k,995)
         lu(k,1007) = lu(k,1007) - lu(k,924) * lu(k,995)
         lu(k,1010) = lu(k,1010) - lu(k,925) * lu(k,995)
         lu(k,1014) = lu(k,1014) - lu(k,926) * lu(k,995)
         lu(k,1018) = - lu(k,927) * lu(k,995)
         lu(k,1061) = lu(k,1061) - lu(k,923) * lu(k,1058)
         lu(k,1062) = lu(k,1062) - lu(k,924) * lu(k,1058)
         lu(k,1064) = lu(k,1064) - lu(k,925) * lu(k,1058)
         lu(k,1066) = - lu(k,926) * lu(k,1058)
         lu(k,1068) = - lu(k,927) * lu(k,1058)
         lu(k,1116) = - lu(k,923) * lu(k,1111)
         lu(k,1117) = - lu(k,924) * lu(k,1111)
         lu(k,1120) = lu(k,1120) - lu(k,925) * lu(k,1111)
         lu(k,1123) = lu(k,1123) - lu(k,926) * lu(k,1111)
         lu(k,1126) = - lu(k,927) * lu(k,1111)
         lu(k,1151) = - lu(k,923) * lu(k,1147)
         lu(k,1155) = lu(k,1155) - lu(k,924) * lu(k,1147)
         lu(k,1158) = lu(k,1158) - lu(k,925) * lu(k,1147)
         lu(k,1162) = - lu(k,926) * lu(k,1147)
         lu(k,1165) = - lu(k,927) * lu(k,1147)
         lu(k,1670) = lu(k,1670) - lu(k,923) * lu(k,1656)
         lu(k,1683) = lu(k,1683) - lu(k,924) * lu(k,1656)
         lu(k,1691) = lu(k,1691) - lu(k,925) * lu(k,1656)
         lu(k,1697) = lu(k,1697) - lu(k,926) * lu(k,1656)
         lu(k,1703) = lu(k,1703) - lu(k,927) * lu(k,1656)
         lu(k,1728) = - lu(k,923) * lu(k,1714)
         lu(k,1741) = lu(k,1741) - lu(k,924) * lu(k,1714)
         lu(k,1748) = lu(k,1748) - lu(k,925) * lu(k,1714)
         lu(k,1754) = lu(k,1754) - lu(k,926) * lu(k,1714)
         lu(k,1760) = lu(k,1760) - lu(k,927) * lu(k,1714)
         lu(k,1821) = lu(k,1821) - lu(k,923) * lu(k,1810)
         lu(k,1833) = lu(k,1833) - lu(k,924) * lu(k,1810)
         lu(k,1840) = lu(k,1840) - lu(k,925) * lu(k,1810)
         lu(k,1846) = lu(k,1846) - lu(k,926) * lu(k,1810)
         lu(k,1852) = lu(k,1852) - lu(k,927) * lu(k,1810)
         lu(k,1928) = lu(k,1928) - lu(k,923) * lu(k,1918)
         lu(k,1939) = lu(k,1939) - lu(k,924) * lu(k,1918)
         lu(k,1947) = lu(k,1947) - lu(k,925) * lu(k,1918)
         lu(k,1953) = lu(k,1953) - lu(k,926) * lu(k,1918)
         lu(k,1959) = lu(k,1959) - lu(k,927) * lu(k,1918)
         lu(k,2047) = lu(k,2047) - lu(k,923) * lu(k,2037)
         lu(k,2058) = lu(k,2058) - lu(k,924) * lu(k,2037)
         lu(k,2064) = lu(k,2064) - lu(k,925) * lu(k,2037)
         lu(k,2070) = lu(k,2070) - lu(k,926) * lu(k,2037)
         lu(k,2076) = lu(k,2076) - lu(k,927) * lu(k,2037)
         lu(k,2107) = lu(k,2107) - lu(k,923) * lu(k,2094)
         lu(k,2118) = lu(k,2118) - lu(k,924) * lu(k,2094)
         lu(k,2125) = lu(k,2125) - lu(k,925) * lu(k,2094)
         lu(k,2131) = lu(k,2131) - lu(k,926) * lu(k,2094)
         lu(k,2137) = lu(k,2137) - lu(k,927) * lu(k,2094)
                                                                        
         lu(k,929) = 1._r8 / lu(k,929)
         lu(k,930) = lu(k,930) * lu(k,929)
         lu(k,931) = lu(k,931) * lu(k,929)
         lu(k,932) = lu(k,932) * lu(k,929)
         lu(k,933) = lu(k,933) * lu(k,929)
         lu(k,934) = lu(k,934) * lu(k,929)
         lu(k,954) = lu(k,954) - lu(k,930) * lu(k,947)
         lu(k,960) = lu(k,960) - lu(k,931) * lu(k,947)
         lu(k,963) = lu(k,963) - lu(k,932) * lu(k,947)
         lu(k,967) = lu(k,967) - lu(k,933) * lu(k,947)
         lu(k,968) = lu(k,968) - lu(k,934) * lu(k,947)
         lu(k,1003) = lu(k,1003) - lu(k,930) * lu(k,996)
         lu(k,1010) = lu(k,1010) - lu(k,931) * lu(k,996)
         lu(k,1013) = lu(k,1013) - lu(k,932) * lu(k,996)
         lu(k,1017) = lu(k,1017) - lu(k,933) * lu(k,996)
         lu(k,1018) = lu(k,1018) - lu(k,934) * lu(k,996)
         lu(k,1026) = lu(k,1026) - lu(k,930) * lu(k,1023)
         lu(k,1032) = lu(k,1032) - lu(k,931) * lu(k,1023)
         lu(k,1035) = lu(k,1035) - lu(k,932) * lu(k,1023)
         lu(k,1039) = lu(k,1039) - lu(k,933) * lu(k,1023)
         lu(k,1040) = lu(k,1040) - lu(k,934) * lu(k,1023)
         lu(k,1208) = lu(k,1208) - lu(k,930) * lu(k,1207)
         lu(k,1215) = lu(k,1215) - lu(k,931) * lu(k,1207)
         lu(k,1216) = lu(k,1216) - lu(k,932) * lu(k,1207)
         lu(k,1218) = - lu(k,933) * lu(k,1207)
         lu(k,1219) = lu(k,1219) - lu(k,934) * lu(k,1207)
         lu(k,1328) = lu(k,1328) - lu(k,930) * lu(k,1325)
         lu(k,1338) = lu(k,1338) - lu(k,931) * lu(k,1325)
         lu(k,1341) = lu(k,1341) - lu(k,932) * lu(k,1325)
         lu(k,1345) = lu(k,1345) - lu(k,933) * lu(k,1325)
         lu(k,1346) = - lu(k,934) * lu(k,1325)
         lu(k,1479) = lu(k,1479) - lu(k,930) * lu(k,1478)
         lu(k,1487) = lu(k,1487) - lu(k,931) * lu(k,1478)
         lu(k,1490) = lu(k,1490) - lu(k,932) * lu(k,1478)
         lu(k,1497) = lu(k,1497) - lu(k,933) * lu(k,1478)
         lu(k,1498) = lu(k,1498) - lu(k,934) * lu(k,1478)
         lu(k,1668) = lu(k,1668) - lu(k,930) * lu(k,1657)
         lu(k,1691) = lu(k,1691) - lu(k,931) * lu(k,1657)
         lu(k,1694) = lu(k,1694) - lu(k,932) * lu(k,1657)
         lu(k,1702) = lu(k,1702) - lu(k,933) * lu(k,1657)
         lu(k,1703) = lu(k,1703) - lu(k,934) * lu(k,1657)
         lu(k,1726) = lu(k,1726) - lu(k,930) * lu(k,1715)
         lu(k,1748) = lu(k,1748) - lu(k,931) * lu(k,1715)
         lu(k,1751) = lu(k,1751) - lu(k,932) * lu(k,1715)
         lu(k,1759) = lu(k,1759) - lu(k,933) * lu(k,1715)
         lu(k,1760) = lu(k,1760) - lu(k,934) * lu(k,1715)
         lu(k,1819) = lu(k,1819) - lu(k,930) * lu(k,1811)
         lu(k,1840) = lu(k,1840) - lu(k,931) * lu(k,1811)
         lu(k,1843) = lu(k,1843) - lu(k,932) * lu(k,1811)
         lu(k,1851) = lu(k,1851) - lu(k,933) * lu(k,1811)
         lu(k,1852) = lu(k,1852) - lu(k,934) * lu(k,1811)
         lu(k,1926) = lu(k,1926) - lu(k,930) * lu(k,1919)
         lu(k,1947) = lu(k,1947) - lu(k,931) * lu(k,1919)
         lu(k,1950) = lu(k,1950) - lu(k,932) * lu(k,1919)
         lu(k,1958) = lu(k,1958) - lu(k,933) * lu(k,1919)
         lu(k,1959) = lu(k,1959) - lu(k,934) * lu(k,1919)
         lu(k,2003) = lu(k,2003) - lu(k,930) * lu(k,1999)
         lu(k,2012) = lu(k,2012) - lu(k,931) * lu(k,1999)
         lu(k,2015) = lu(k,2015) - lu(k,932) * lu(k,1999)
         lu(k,2023) = lu(k,2023) - lu(k,933) * lu(k,1999)
         lu(k,2024) = lu(k,2024) - lu(k,934) * lu(k,1999)
         lu(k,2105) = lu(k,2105) - lu(k,930) * lu(k,2095)
         lu(k,2125) = lu(k,2125) - lu(k,931) * lu(k,2095)
         lu(k,2128) = lu(k,2128) - lu(k,932) * lu(k,2095)
         lu(k,2136) = lu(k,2136) - lu(k,933) * lu(k,2095)
         lu(k,2137) = lu(k,2137) - lu(k,934) * lu(k,2095)
                                                                        
      end do
                                                                        
      end subroutine lu_fac20
                                                                        
      subroutine lu_fac21( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,948) = 1._r8 / lu(k,948)
         lu(k,949) = lu(k,949) * lu(k,948)
         lu(k,950) = lu(k,950) * lu(k,948)
         lu(k,951) = lu(k,951) * lu(k,948)
         lu(k,952) = lu(k,952) * lu(k,948)
         lu(k,953) = lu(k,953) * lu(k,948)
         lu(k,954) = lu(k,954) * lu(k,948)
         lu(k,955) = lu(k,955) * lu(k,948)
         lu(k,956) = lu(k,956) * lu(k,948)
         lu(k,957) = lu(k,957) * lu(k,948)
         lu(k,958) = lu(k,958) * lu(k,948)
         lu(k,959) = lu(k,959) * lu(k,948)
         lu(k,960) = lu(k,960) * lu(k,948)
         lu(k,961) = lu(k,961) * lu(k,948)
         lu(k,962) = lu(k,962) * lu(k,948)
         lu(k,963) = lu(k,963) * lu(k,948)
         lu(k,964) = lu(k,964) * lu(k,948)
         lu(k,965) = lu(k,965) * lu(k,948)
         lu(k,966) = lu(k,966) * lu(k,948)
         lu(k,967) = lu(k,967) * lu(k,948)
         lu(k,968) = lu(k,968) * lu(k,948)
         lu(k,1659) = lu(k,1659) - lu(k,949) * lu(k,1658)
         lu(k,1660) = lu(k,1660) - lu(k,950) * lu(k,1658)
         lu(k,1663) = lu(k,1663) - lu(k,951) * lu(k,1658)
         lu(k,1664) = lu(k,1664) - lu(k,952) * lu(k,1658)
         lu(k,1666) = lu(k,1666) - lu(k,953) * lu(k,1658)
         lu(k,1668) = lu(k,1668) - lu(k,954) * lu(k,1658)
         lu(k,1670) = lu(k,1670) - lu(k,955) * lu(k,1658)
         lu(k,1676) = lu(k,1676) - lu(k,956) * lu(k,1658)
         lu(k,1683) = lu(k,1683) - lu(k,957) * lu(k,1658)
         lu(k,1687) = lu(k,1687) - lu(k,958) * lu(k,1658)
         lu(k,1689) = lu(k,1689) - lu(k,959) * lu(k,1658)
         lu(k,1691) = lu(k,1691) - lu(k,960) * lu(k,1658)
         lu(k,1692) = lu(k,1692) - lu(k,961) * lu(k,1658)
         lu(k,1693) = lu(k,1693) - lu(k,962) * lu(k,1658)
         lu(k,1694) = lu(k,1694) - lu(k,963) * lu(k,1658)
         lu(k,1697) = lu(k,1697) - lu(k,964) * lu(k,1658)
         lu(k,1698) = lu(k,1698) - lu(k,965) * lu(k,1658)
         lu(k,1700) = lu(k,1700) - lu(k,966) * lu(k,1658)
         lu(k,1702) = lu(k,1702) - lu(k,967) * lu(k,1658)
         lu(k,1703) = lu(k,1703) - lu(k,968) * lu(k,1658)
         lu(k,1717) = lu(k,1717) - lu(k,949) * lu(k,1716)
         lu(k,1718) = lu(k,1718) - lu(k,950) * lu(k,1716)
         lu(k,1721) = lu(k,1721) - lu(k,951) * lu(k,1716)
         lu(k,1722) = - lu(k,952) * lu(k,1716)
         lu(k,1724) = lu(k,1724) - lu(k,953) * lu(k,1716)
         lu(k,1726) = lu(k,1726) - lu(k,954) * lu(k,1716)
         lu(k,1728) = lu(k,1728) - lu(k,955) * lu(k,1716)
         lu(k,1734) = lu(k,1734) - lu(k,956) * lu(k,1716)
         lu(k,1741) = lu(k,1741) - lu(k,957) * lu(k,1716)
         lu(k,1744) = lu(k,1744) - lu(k,958) * lu(k,1716)
         lu(k,1746) = lu(k,1746) - lu(k,959) * lu(k,1716)
         lu(k,1748) = lu(k,1748) - lu(k,960) * lu(k,1716)
         lu(k,1749) = lu(k,1749) - lu(k,961) * lu(k,1716)
         lu(k,1750) = lu(k,1750) - lu(k,962) * lu(k,1716)
         lu(k,1751) = lu(k,1751) - lu(k,963) * lu(k,1716)
         lu(k,1754) = lu(k,1754) - lu(k,964) * lu(k,1716)
         lu(k,1755) = lu(k,1755) - lu(k,965) * lu(k,1716)
         lu(k,1757) = lu(k,1757) - lu(k,966) * lu(k,1716)
         lu(k,1759) = lu(k,1759) - lu(k,967) * lu(k,1716)
         lu(k,1760) = lu(k,1760) - lu(k,968) * lu(k,1716)
         lu(k,2097) = lu(k,2097) - lu(k,949) * lu(k,2096)
         lu(k,2098) = lu(k,2098) - lu(k,950) * lu(k,2096)
         lu(k,2101) = - lu(k,951) * lu(k,2096)
         lu(k,2102) = lu(k,2102) - lu(k,952) * lu(k,2096)
         lu(k,2103) = lu(k,2103) - lu(k,953) * lu(k,2096)
         lu(k,2105) = lu(k,2105) - lu(k,954) * lu(k,2096)
         lu(k,2107) = lu(k,2107) - lu(k,955) * lu(k,2096)
         lu(k,2111) = lu(k,2111) - lu(k,956) * lu(k,2096)
         lu(k,2118) = lu(k,2118) - lu(k,957) * lu(k,2096)
         lu(k,2121) = - lu(k,958) * lu(k,2096)
         lu(k,2123) = lu(k,2123) - lu(k,959) * lu(k,2096)
         lu(k,2125) = lu(k,2125) - lu(k,960) * lu(k,2096)
         lu(k,2126) = lu(k,2126) - lu(k,961) * lu(k,2096)
         lu(k,2127) = lu(k,2127) - lu(k,962) * lu(k,2096)
         lu(k,2128) = lu(k,2128) - lu(k,963) * lu(k,2096)
         lu(k,2131) = lu(k,2131) - lu(k,964) * lu(k,2096)
         lu(k,2132) = lu(k,2132) - lu(k,965) * lu(k,2096)
         lu(k,2134) = lu(k,2134) - lu(k,966) * lu(k,2096)
         lu(k,2136) = lu(k,2136) - lu(k,967) * lu(k,2096)
         lu(k,2137) = lu(k,2137) - lu(k,968) * lu(k,2096)
                                                                        
         lu(k,969) = 1._r8 / lu(k,969)
         lu(k,970) = lu(k,970) * lu(k,969)
         lu(k,971) = lu(k,971) * lu(k,969)
         lu(k,972) = lu(k,972) * lu(k,969)
         lu(k,973) = lu(k,973) * lu(k,969)
         lu(k,974) = lu(k,974) * lu(k,969)
         lu(k,975) = lu(k,975) * lu(k,969)
         lu(k,976) = lu(k,976) * lu(k,969)
         lu(k,1000) = lu(k,1000) - lu(k,970) * lu(k,997)
         lu(k,1001) = lu(k,1001) - lu(k,971) * lu(k,997)
         lu(k,1003) = lu(k,1003) - lu(k,972) * lu(k,997)
         lu(k,1004) = - lu(k,973) * lu(k,997)
         lu(k,1010) = lu(k,1010) - lu(k,974) * lu(k,997)
         lu(k,1011) = lu(k,1011) - lu(k,975) * lu(k,997)
         lu(k,1013) = lu(k,1013) - lu(k,976) * lu(k,997)
         lu(k,1045) = lu(k,1045) - lu(k,970) * lu(k,1044)
         lu(k,1046) = - lu(k,971) * lu(k,1044)
         lu(k,1047) = - lu(k,972) * lu(k,1044)
         lu(k,1048) = - lu(k,973) * lu(k,1044)
         lu(k,1051) = lu(k,1051) - lu(k,974) * lu(k,1044)
         lu(k,1052) = lu(k,1052) - lu(k,975) * lu(k,1044)
         lu(k,1054) = lu(k,1054) - lu(k,976) * lu(k,1044)
         lu(k,1663) = lu(k,1663) - lu(k,970) * lu(k,1659)
         lu(k,1664) = lu(k,1664) - lu(k,971) * lu(k,1659)
         lu(k,1668) = lu(k,1668) - lu(k,972) * lu(k,1659)
         lu(k,1669) = lu(k,1669) - lu(k,973) * lu(k,1659)
         lu(k,1691) = lu(k,1691) - lu(k,974) * lu(k,1659)
         lu(k,1692) = lu(k,1692) - lu(k,975) * lu(k,1659)
         lu(k,1694) = lu(k,1694) - lu(k,976) * lu(k,1659)
         lu(k,1721) = lu(k,1721) - lu(k,970) * lu(k,1717)
         lu(k,1722) = lu(k,1722) - lu(k,971) * lu(k,1717)
         lu(k,1726) = lu(k,1726) - lu(k,972) * lu(k,1717)
         lu(k,1727) = lu(k,1727) - lu(k,973) * lu(k,1717)
         lu(k,1748) = lu(k,1748) - lu(k,974) * lu(k,1717)
         lu(k,1749) = lu(k,1749) - lu(k,975) * lu(k,1717)
         lu(k,1751) = lu(k,1751) - lu(k,976) * lu(k,1717)
         lu(k,1814) = lu(k,1814) - lu(k,970) * lu(k,1812)
         lu(k,1815) = lu(k,1815) - lu(k,971) * lu(k,1812)
         lu(k,1819) = lu(k,1819) - lu(k,972) * lu(k,1812)
         lu(k,1820) = lu(k,1820) - lu(k,973) * lu(k,1812)
         lu(k,1840) = lu(k,1840) - lu(k,974) * lu(k,1812)
         lu(k,1841) = lu(k,1841) - lu(k,975) * lu(k,1812)
         lu(k,1843) = lu(k,1843) - lu(k,976) * lu(k,1812)
         lu(k,1922) = lu(k,1922) - lu(k,970) * lu(k,1920)
         lu(k,1923) = lu(k,1923) - lu(k,971) * lu(k,1920)
         lu(k,1926) = lu(k,1926) - lu(k,972) * lu(k,1920)
         lu(k,1927) = lu(k,1927) - lu(k,973) * lu(k,1920)
         lu(k,1947) = lu(k,1947) - lu(k,974) * lu(k,1920)
         lu(k,1948) = lu(k,1948) - lu(k,975) * lu(k,1920)
         lu(k,1950) = lu(k,1950) - lu(k,976) * lu(k,1920)
         lu(k,2040) = lu(k,2040) - lu(k,970) * lu(k,2038)
         lu(k,2041) = lu(k,2041) - lu(k,971) * lu(k,2038)
         lu(k,2045) = lu(k,2045) - lu(k,972) * lu(k,2038)
         lu(k,2046) = lu(k,2046) - lu(k,973) * lu(k,2038)
         lu(k,2064) = lu(k,2064) - lu(k,974) * lu(k,2038)
         lu(k,2065) = - lu(k,975) * lu(k,2038)
         lu(k,2067) = lu(k,2067) - lu(k,976) * lu(k,2038)
         lu(k,2101) = lu(k,2101) - lu(k,970) * lu(k,2097)
         lu(k,2102) = lu(k,2102) - lu(k,971) * lu(k,2097)
         lu(k,2105) = lu(k,2105) - lu(k,972) * lu(k,2097)
         lu(k,2106) = - lu(k,973) * lu(k,2097)
         lu(k,2125) = lu(k,2125) - lu(k,974) * lu(k,2097)
         lu(k,2126) = lu(k,2126) - lu(k,975) * lu(k,2097)
         lu(k,2128) = lu(k,2128) - lu(k,976) * lu(k,2097)
                                                                        
         lu(k,979) = 1._r8 / lu(k,979)
         lu(k,980) = lu(k,980) * lu(k,979)
         lu(k,981) = lu(k,981) * lu(k,979)
         lu(k,982) = lu(k,982) * lu(k,979)
         lu(k,983) = lu(k,983) * lu(k,979)
         lu(k,1003) = lu(k,1003) - lu(k,980) * lu(k,998)
         lu(k,1010) = lu(k,1010) - lu(k,981) * lu(k,998)
         lu(k,1013) = lu(k,1013) - lu(k,982) * lu(k,998)
         lu(k,1017) = lu(k,1017) - lu(k,983) * lu(k,998)
         lu(k,1074) = lu(k,1074) - lu(k,980) * lu(k,1072)
         lu(k,1079) = lu(k,1079) - lu(k,981) * lu(k,1072)
         lu(k,1081) = lu(k,1081) - lu(k,982) * lu(k,1072)
         lu(k,1084) = - lu(k,983) * lu(k,1072)
         lu(k,1097) = lu(k,1097) - lu(k,980) * lu(k,1095)
         lu(k,1099) = lu(k,1099) - lu(k,981) * lu(k,1095)
         lu(k,1100) = lu(k,1100) - lu(k,982) * lu(k,1095)
         lu(k,1101) = lu(k,1101) - lu(k,983) * lu(k,1095)
         lu(k,1171) = lu(k,1171) - lu(k,980) * lu(k,1169)
         lu(k,1177) = lu(k,1177) - lu(k,981) * lu(k,1169)
         lu(k,1180) = lu(k,1180) - lu(k,982) * lu(k,1169)
         lu(k,1183) = lu(k,1183) - lu(k,983) * lu(k,1169)
         lu(k,1282) = lu(k,1282) - lu(k,980) * lu(k,1280)
         lu(k,1295) = lu(k,1295) - lu(k,981) * lu(k,1280)
         lu(k,1298) = lu(k,1298) - lu(k,982) * lu(k,1280)
         lu(k,1302) = - lu(k,983) * lu(k,1280)
         lu(k,1375) = lu(k,1375) - lu(k,980) * lu(k,1372)
         lu(k,1390) = lu(k,1390) - lu(k,981) * lu(k,1372)
         lu(k,1393) = lu(k,1393) - lu(k,982) * lu(k,1372)
         lu(k,1397) = lu(k,1397) - lu(k,983) * lu(k,1372)
         lu(k,1668) = lu(k,1668) - lu(k,980) * lu(k,1660)
         lu(k,1691) = lu(k,1691) - lu(k,981) * lu(k,1660)
         lu(k,1694) = lu(k,1694) - lu(k,982) * lu(k,1660)
         lu(k,1702) = lu(k,1702) - lu(k,983) * lu(k,1660)
         lu(k,1726) = lu(k,1726) - lu(k,980) * lu(k,1718)
         lu(k,1748) = lu(k,1748) - lu(k,981) * lu(k,1718)
         lu(k,1751) = lu(k,1751) - lu(k,982) * lu(k,1718)
         lu(k,1759) = lu(k,1759) - lu(k,983) * lu(k,1718)
         lu(k,1819) = lu(k,1819) - lu(k,980) * lu(k,1813)
         lu(k,1840) = lu(k,1840) - lu(k,981) * lu(k,1813)
         lu(k,1843) = lu(k,1843) - lu(k,982) * lu(k,1813)
         lu(k,1851) = lu(k,1851) - lu(k,983) * lu(k,1813)
         lu(k,1926) = lu(k,1926) - lu(k,980) * lu(k,1921)
         lu(k,1947) = lu(k,1947) - lu(k,981) * lu(k,1921)
         lu(k,1950) = lu(k,1950) - lu(k,982) * lu(k,1921)
         lu(k,1958) = lu(k,1958) - lu(k,983) * lu(k,1921)
         lu(k,2003) = lu(k,2003) - lu(k,980) * lu(k,2000)
         lu(k,2012) = lu(k,2012) - lu(k,981) * lu(k,2000)
         lu(k,2015) = lu(k,2015) - lu(k,982) * lu(k,2000)
         lu(k,2023) = lu(k,2023) - lu(k,983) * lu(k,2000)
         lu(k,2045) = lu(k,2045) - lu(k,980) * lu(k,2039)
         lu(k,2064) = lu(k,2064) - lu(k,981) * lu(k,2039)
         lu(k,2067) = lu(k,2067) - lu(k,982) * lu(k,2039)
         lu(k,2075) = lu(k,2075) - lu(k,983) * lu(k,2039)
         lu(k,2105) = lu(k,2105) - lu(k,980) * lu(k,2098)
         lu(k,2125) = lu(k,2125) - lu(k,981) * lu(k,2098)
         lu(k,2128) = lu(k,2128) - lu(k,982) * lu(k,2098)
         lu(k,2136) = lu(k,2136) - lu(k,983) * lu(k,2098)
         lu(k,2180) = lu(k,2180) - lu(k,980) * lu(k,2179)
         lu(k,2192) = lu(k,2192) - lu(k,981) * lu(k,2179)
         lu(k,2195) = lu(k,2195) - lu(k,982) * lu(k,2179)
         lu(k,2203) = lu(k,2203) - lu(k,983) * lu(k,2179)
                                                                        
         lu(k,999) = 1._r8 / lu(k,999)
         lu(k,1000) = lu(k,1000) * lu(k,999)
         lu(k,1001) = lu(k,1001) * lu(k,999)
         lu(k,1002) = lu(k,1002) * lu(k,999)
         lu(k,1003) = lu(k,1003) * lu(k,999)
         lu(k,1004) = lu(k,1004) * lu(k,999)
         lu(k,1005) = lu(k,1005) * lu(k,999)
         lu(k,1006) = lu(k,1006) * lu(k,999)
         lu(k,1007) = lu(k,1007) * lu(k,999)
         lu(k,1008) = lu(k,1008) * lu(k,999)
         lu(k,1009) = lu(k,1009) * lu(k,999)
         lu(k,1010) = lu(k,1010) * lu(k,999)
         lu(k,1011) = lu(k,1011) * lu(k,999)
         lu(k,1012) = lu(k,1012) * lu(k,999)
         lu(k,1013) = lu(k,1013) * lu(k,999)
         lu(k,1014) = lu(k,1014) * lu(k,999)
         lu(k,1015) = lu(k,1015) * lu(k,999)
         lu(k,1016) = lu(k,1016) * lu(k,999)
         lu(k,1017) = lu(k,1017) * lu(k,999)
         lu(k,1018) = lu(k,1018) * lu(k,999)
         lu(k,1663) = lu(k,1663) - lu(k,1000) * lu(k,1661)
         lu(k,1664) = lu(k,1664) - lu(k,1001) * lu(k,1661)
         lu(k,1666) = lu(k,1666) - lu(k,1002) * lu(k,1661)
         lu(k,1668) = lu(k,1668) - lu(k,1003) * lu(k,1661)
         lu(k,1669) = lu(k,1669) - lu(k,1004) * lu(k,1661)
         lu(k,1670) = lu(k,1670) - lu(k,1005) * lu(k,1661)
         lu(k,1676) = lu(k,1676) - lu(k,1006) * lu(k,1661)
         lu(k,1683) = lu(k,1683) - lu(k,1007) * lu(k,1661)
         lu(k,1687) = lu(k,1687) - lu(k,1008) * lu(k,1661)
         lu(k,1689) = lu(k,1689) - lu(k,1009) * lu(k,1661)
         lu(k,1691) = lu(k,1691) - lu(k,1010) * lu(k,1661)
         lu(k,1692) = lu(k,1692) - lu(k,1011) * lu(k,1661)
         lu(k,1693) = lu(k,1693) - lu(k,1012) * lu(k,1661)
         lu(k,1694) = lu(k,1694) - lu(k,1013) * lu(k,1661)
         lu(k,1697) = lu(k,1697) - lu(k,1014) * lu(k,1661)
         lu(k,1698) = lu(k,1698) - lu(k,1015) * lu(k,1661)
         lu(k,1700) = lu(k,1700) - lu(k,1016) * lu(k,1661)
         lu(k,1702) = lu(k,1702) - lu(k,1017) * lu(k,1661)
         lu(k,1703) = lu(k,1703) - lu(k,1018) * lu(k,1661)
         lu(k,1721) = lu(k,1721) - lu(k,1000) * lu(k,1719)
         lu(k,1722) = lu(k,1722) - lu(k,1001) * lu(k,1719)
         lu(k,1724) = lu(k,1724) - lu(k,1002) * lu(k,1719)
         lu(k,1726) = lu(k,1726) - lu(k,1003) * lu(k,1719)
         lu(k,1727) = lu(k,1727) - lu(k,1004) * lu(k,1719)
         lu(k,1728) = lu(k,1728) - lu(k,1005) * lu(k,1719)
         lu(k,1734) = lu(k,1734) - lu(k,1006) * lu(k,1719)
         lu(k,1741) = lu(k,1741) - lu(k,1007) * lu(k,1719)
         lu(k,1744) = lu(k,1744) - lu(k,1008) * lu(k,1719)
         lu(k,1746) = lu(k,1746) - lu(k,1009) * lu(k,1719)
         lu(k,1748) = lu(k,1748) - lu(k,1010) * lu(k,1719)
         lu(k,1749) = lu(k,1749) - lu(k,1011) * lu(k,1719)
         lu(k,1750) = lu(k,1750) - lu(k,1012) * lu(k,1719)
         lu(k,1751) = lu(k,1751) - lu(k,1013) * lu(k,1719)
         lu(k,1754) = lu(k,1754) - lu(k,1014) * lu(k,1719)
         lu(k,1755) = lu(k,1755) - lu(k,1015) * lu(k,1719)
         lu(k,1757) = lu(k,1757) - lu(k,1016) * lu(k,1719)
         lu(k,1759) = lu(k,1759) - lu(k,1017) * lu(k,1719)
         lu(k,1760) = lu(k,1760) - lu(k,1018) * lu(k,1719)
         lu(k,2101) = lu(k,2101) - lu(k,1000) * lu(k,2099)
         lu(k,2102) = lu(k,2102) - lu(k,1001) * lu(k,2099)
         lu(k,2103) = lu(k,2103) - lu(k,1002) * lu(k,2099)
         lu(k,2105) = lu(k,2105) - lu(k,1003) * lu(k,2099)
         lu(k,2106) = lu(k,2106) - lu(k,1004) * lu(k,2099)
         lu(k,2107) = lu(k,2107) - lu(k,1005) * lu(k,2099)
         lu(k,2111) = lu(k,2111) - lu(k,1006) * lu(k,2099)
         lu(k,2118) = lu(k,2118) - lu(k,1007) * lu(k,2099)
         lu(k,2121) = lu(k,2121) - lu(k,1008) * lu(k,2099)
         lu(k,2123) = lu(k,2123) - lu(k,1009) * lu(k,2099)
         lu(k,2125) = lu(k,2125) - lu(k,1010) * lu(k,2099)
         lu(k,2126) = lu(k,2126) - lu(k,1011) * lu(k,2099)
         lu(k,2127) = lu(k,2127) - lu(k,1012) * lu(k,2099)
         lu(k,2128) = lu(k,2128) - lu(k,1013) * lu(k,2099)
         lu(k,2131) = lu(k,2131) - lu(k,1014) * lu(k,2099)
         lu(k,2132) = lu(k,2132) - lu(k,1015) * lu(k,2099)
         lu(k,2134) = lu(k,2134) - lu(k,1016) * lu(k,2099)
         lu(k,2136) = lu(k,2136) - lu(k,1017) * lu(k,2099)
         lu(k,2137) = lu(k,2137) - lu(k,1018) * lu(k,2099)
                                                                        
      end do
                                                                        
      end subroutine lu_fac21
                                                                        
      subroutine lu_fac22( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,1024) = 1._r8 / lu(k,1024)
         lu(k,1025) = lu(k,1025) * lu(k,1024)
         lu(k,1026) = lu(k,1026) * lu(k,1024)
         lu(k,1027) = lu(k,1027) * lu(k,1024)
         lu(k,1028) = lu(k,1028) * lu(k,1024)
         lu(k,1029) = lu(k,1029) * lu(k,1024)
         lu(k,1030) = lu(k,1030) * lu(k,1024)
         lu(k,1031) = lu(k,1031) * lu(k,1024)
         lu(k,1032) = lu(k,1032) * lu(k,1024)
         lu(k,1033) = lu(k,1033) * lu(k,1024)
         lu(k,1034) = lu(k,1034) * lu(k,1024)
         lu(k,1035) = lu(k,1035) * lu(k,1024)
         lu(k,1036) = lu(k,1036) * lu(k,1024)
         lu(k,1037) = lu(k,1037) * lu(k,1024)
         lu(k,1038) = lu(k,1038) * lu(k,1024)
         lu(k,1039) = lu(k,1039) * lu(k,1024)
         lu(k,1040) = lu(k,1040) * lu(k,1024)
         lu(k,1327) = lu(k,1327) - lu(k,1025) * lu(k,1326)
         lu(k,1328) = lu(k,1328) - lu(k,1026) * lu(k,1326)
         lu(k,1329) = - lu(k,1027) * lu(k,1326)
         lu(k,1330) = lu(k,1330) - lu(k,1028) * lu(k,1326)
         lu(k,1334) = lu(k,1334) - lu(k,1029) * lu(k,1326)
         lu(k,1335) = - lu(k,1030) * lu(k,1326)
         lu(k,1337) = lu(k,1337) - lu(k,1031) * lu(k,1326)
         lu(k,1338) = lu(k,1338) - lu(k,1032) * lu(k,1326)
         lu(k,1339) = - lu(k,1033) * lu(k,1326)
         lu(k,1340) = - lu(k,1034) * lu(k,1326)
         lu(k,1341) = lu(k,1341) - lu(k,1035) * lu(k,1326)
         lu(k,1342) = lu(k,1342) - lu(k,1036) * lu(k,1326)
         lu(k,1343) = lu(k,1343) - lu(k,1037) * lu(k,1326)
         lu(k,1344) = - lu(k,1038) * lu(k,1326)
         lu(k,1345) = lu(k,1345) - lu(k,1039) * lu(k,1326)
         lu(k,1346) = lu(k,1346) - lu(k,1040) * lu(k,1326)
         lu(k,1666) = lu(k,1666) - lu(k,1025) * lu(k,1662)
         lu(k,1668) = lu(k,1668) - lu(k,1026) * lu(k,1662)
         lu(k,1671) = lu(k,1671) - lu(k,1027) * lu(k,1662)
         lu(k,1676) = lu(k,1676) - lu(k,1028) * lu(k,1662)
         lu(k,1683) = lu(k,1683) - lu(k,1029) * lu(k,1662)
         lu(k,1686) = lu(k,1686) - lu(k,1030) * lu(k,1662)
         lu(k,1689) = lu(k,1689) - lu(k,1031) * lu(k,1662)
         lu(k,1691) = lu(k,1691) - lu(k,1032) * lu(k,1662)
         lu(k,1692) = lu(k,1692) - lu(k,1033) * lu(k,1662)
         lu(k,1693) = lu(k,1693) - lu(k,1034) * lu(k,1662)
         lu(k,1694) = lu(k,1694) - lu(k,1035) * lu(k,1662)
         lu(k,1697) = lu(k,1697) - lu(k,1036) * lu(k,1662)
         lu(k,1698) = lu(k,1698) - lu(k,1037) * lu(k,1662)
         lu(k,1700) = lu(k,1700) - lu(k,1038) * lu(k,1662)
         lu(k,1702) = lu(k,1702) - lu(k,1039) * lu(k,1662)
         lu(k,1703) = lu(k,1703) - lu(k,1040) * lu(k,1662)
         lu(k,1724) = lu(k,1724) - lu(k,1025) * lu(k,1720)
         lu(k,1726) = lu(k,1726) - lu(k,1026) * lu(k,1720)
         lu(k,1729) = lu(k,1729) - lu(k,1027) * lu(k,1720)
         lu(k,1734) = lu(k,1734) - lu(k,1028) * lu(k,1720)
         lu(k,1741) = lu(k,1741) - lu(k,1029) * lu(k,1720)
         lu(k,1743) = - lu(k,1030) * lu(k,1720)
         lu(k,1746) = lu(k,1746) - lu(k,1031) * lu(k,1720)
         lu(k,1748) = lu(k,1748) - lu(k,1032) * lu(k,1720)
         lu(k,1749) = lu(k,1749) - lu(k,1033) * lu(k,1720)
         lu(k,1750) = lu(k,1750) - lu(k,1034) * lu(k,1720)
         lu(k,1751) = lu(k,1751) - lu(k,1035) * lu(k,1720)
         lu(k,1754) = lu(k,1754) - lu(k,1036) * lu(k,1720)
         lu(k,1755) = lu(k,1755) - lu(k,1037) * lu(k,1720)
         lu(k,1757) = lu(k,1757) - lu(k,1038) * lu(k,1720)
         lu(k,1759) = lu(k,1759) - lu(k,1039) * lu(k,1720)
         lu(k,1760) = lu(k,1760) - lu(k,1040) * lu(k,1720)
         lu(k,2103) = lu(k,2103) - lu(k,1025) * lu(k,2100)
         lu(k,2105) = lu(k,2105) - lu(k,1026) * lu(k,2100)
         lu(k,2108) = - lu(k,1027) * lu(k,2100)
         lu(k,2111) = lu(k,2111) - lu(k,1028) * lu(k,2100)
         lu(k,2118) = lu(k,2118) - lu(k,1029) * lu(k,2100)
         lu(k,2120) = lu(k,2120) - lu(k,1030) * lu(k,2100)
         lu(k,2123) = lu(k,2123) - lu(k,1031) * lu(k,2100)
         lu(k,2125) = lu(k,2125) - lu(k,1032) * lu(k,2100)
         lu(k,2126) = lu(k,2126) - lu(k,1033) * lu(k,2100)
         lu(k,2127) = lu(k,2127) - lu(k,1034) * lu(k,2100)
         lu(k,2128) = lu(k,2128) - lu(k,1035) * lu(k,2100)
         lu(k,2131) = lu(k,2131) - lu(k,1036) * lu(k,2100)
         lu(k,2132) = lu(k,2132) - lu(k,1037) * lu(k,2100)
         lu(k,2134) = lu(k,2134) - lu(k,1038) * lu(k,2100)
         lu(k,2136) = lu(k,2136) - lu(k,1039) * lu(k,2100)
         lu(k,2137) = lu(k,2137) - lu(k,1040) * lu(k,2100)
                                                                        
         lu(k,1045) = 1._r8 / lu(k,1045)
         lu(k,1046) = lu(k,1046) * lu(k,1045)
         lu(k,1047) = lu(k,1047) * lu(k,1045)
         lu(k,1048) = lu(k,1048) * lu(k,1045)
         lu(k,1049) = lu(k,1049) * lu(k,1045)
         lu(k,1050) = lu(k,1050) * lu(k,1045)
         lu(k,1051) = lu(k,1051) * lu(k,1045)
         lu(k,1052) = lu(k,1052) * lu(k,1045)
         lu(k,1053) = lu(k,1053) * lu(k,1045)
         lu(k,1054) = lu(k,1054) * lu(k,1045)
         lu(k,1055) = lu(k,1055) * lu(k,1045)
         lu(k,1056) = lu(k,1056) * lu(k,1045)
         lu(k,1664) = lu(k,1664) - lu(k,1046) * lu(k,1663)
         lu(k,1668) = lu(k,1668) - lu(k,1047) * lu(k,1663)
         lu(k,1669) = lu(k,1669) - lu(k,1048) * lu(k,1663)
         lu(k,1687) = lu(k,1687) - lu(k,1049) * lu(k,1663)
         lu(k,1689) = lu(k,1689) - lu(k,1050) * lu(k,1663)
         lu(k,1691) = lu(k,1691) - lu(k,1051) * lu(k,1663)
         lu(k,1692) = lu(k,1692) - lu(k,1052) * lu(k,1663)
         lu(k,1693) = lu(k,1693) - lu(k,1053) * lu(k,1663)
         lu(k,1694) = lu(k,1694) - lu(k,1054) * lu(k,1663)
         lu(k,1697) = lu(k,1697) - lu(k,1055) * lu(k,1663)
         lu(k,1700) = lu(k,1700) - lu(k,1056) * lu(k,1663)
         lu(k,1722) = lu(k,1722) - lu(k,1046) * lu(k,1721)
         lu(k,1726) = lu(k,1726) - lu(k,1047) * lu(k,1721)
         lu(k,1727) = lu(k,1727) - lu(k,1048) * lu(k,1721)
         lu(k,1744) = lu(k,1744) - lu(k,1049) * lu(k,1721)
         lu(k,1746) = lu(k,1746) - lu(k,1050) * lu(k,1721)
         lu(k,1748) = lu(k,1748) - lu(k,1051) * lu(k,1721)
         lu(k,1749) = lu(k,1749) - lu(k,1052) * lu(k,1721)
         lu(k,1750) = lu(k,1750) - lu(k,1053) * lu(k,1721)
         lu(k,1751) = lu(k,1751) - lu(k,1054) * lu(k,1721)
         lu(k,1754) = lu(k,1754) - lu(k,1055) * lu(k,1721)
         lu(k,1757) = lu(k,1757) - lu(k,1056) * lu(k,1721)
         lu(k,1815) = lu(k,1815) - lu(k,1046) * lu(k,1814)
         lu(k,1819) = lu(k,1819) - lu(k,1047) * lu(k,1814)
         lu(k,1820) = lu(k,1820) - lu(k,1048) * lu(k,1814)
         lu(k,1836) = lu(k,1836) - lu(k,1049) * lu(k,1814)
         lu(k,1838) = lu(k,1838) - lu(k,1050) * lu(k,1814)
         lu(k,1840) = lu(k,1840) - lu(k,1051) * lu(k,1814)
         lu(k,1841) = lu(k,1841) - lu(k,1052) * lu(k,1814)
         lu(k,1842) = lu(k,1842) - lu(k,1053) * lu(k,1814)
         lu(k,1843) = lu(k,1843) - lu(k,1054) * lu(k,1814)
         lu(k,1846) = lu(k,1846) - lu(k,1055) * lu(k,1814)
         lu(k,1849) = lu(k,1849) - lu(k,1056) * lu(k,1814)
         lu(k,1923) = lu(k,1923) - lu(k,1046) * lu(k,1922)
         lu(k,1926) = lu(k,1926) - lu(k,1047) * lu(k,1922)
         lu(k,1927) = lu(k,1927) - lu(k,1048) * lu(k,1922)
         lu(k,1943) = lu(k,1943) - lu(k,1049) * lu(k,1922)
         lu(k,1945) = lu(k,1945) - lu(k,1050) * lu(k,1922)
         lu(k,1947) = lu(k,1947) - lu(k,1051) * lu(k,1922)
         lu(k,1948) = lu(k,1948) - lu(k,1052) * lu(k,1922)
         lu(k,1949) = lu(k,1949) - lu(k,1053) * lu(k,1922)
         lu(k,1950) = lu(k,1950) - lu(k,1054) * lu(k,1922)
         lu(k,1953) = lu(k,1953) - lu(k,1055) * lu(k,1922)
         lu(k,1956) = lu(k,1956) - lu(k,1056) * lu(k,1922)
         lu(k,2041) = lu(k,2041) - lu(k,1046) * lu(k,2040)
         lu(k,2045) = lu(k,2045) - lu(k,1047) * lu(k,2040)
         lu(k,2046) = lu(k,2046) - lu(k,1048) * lu(k,2040)
         lu(k,2060) = lu(k,2060) - lu(k,1049) * lu(k,2040)
         lu(k,2062) = lu(k,2062) - lu(k,1050) * lu(k,2040)
         lu(k,2064) = lu(k,2064) - lu(k,1051) * lu(k,2040)
         lu(k,2065) = lu(k,2065) - lu(k,1052) * lu(k,2040)
         lu(k,2066) = lu(k,2066) - lu(k,1053) * lu(k,2040)
         lu(k,2067) = lu(k,2067) - lu(k,1054) * lu(k,2040)
         lu(k,2070) = lu(k,2070) - lu(k,1055) * lu(k,2040)
         lu(k,2073) = lu(k,2073) - lu(k,1056) * lu(k,2040)
         lu(k,2102) = lu(k,2102) - lu(k,1046) * lu(k,2101)
         lu(k,2105) = lu(k,2105) - lu(k,1047) * lu(k,2101)
         lu(k,2106) = lu(k,2106) - lu(k,1048) * lu(k,2101)
         lu(k,2121) = lu(k,2121) - lu(k,1049) * lu(k,2101)
         lu(k,2123) = lu(k,2123) - lu(k,1050) * lu(k,2101)
         lu(k,2125) = lu(k,2125) - lu(k,1051) * lu(k,2101)
         lu(k,2126) = lu(k,2126) - lu(k,1052) * lu(k,2101)
         lu(k,2127) = lu(k,2127) - lu(k,1053) * lu(k,2101)
         lu(k,2128) = lu(k,2128) - lu(k,1054) * lu(k,2101)
         lu(k,2131) = lu(k,2131) - lu(k,1055) * lu(k,2101)
         lu(k,2134) = lu(k,2134) - lu(k,1056) * lu(k,2101)
                                                                        
         lu(k,1059) = 1._r8 / lu(k,1059)
         lu(k,1060) = lu(k,1060) * lu(k,1059)
         lu(k,1061) = lu(k,1061) * lu(k,1059)
         lu(k,1062) = lu(k,1062) * lu(k,1059)
         lu(k,1063) = lu(k,1063) * lu(k,1059)
         lu(k,1064) = lu(k,1064) * lu(k,1059)
         lu(k,1065) = lu(k,1065) * lu(k,1059)
         lu(k,1066) = lu(k,1066) * lu(k,1059)
         lu(k,1067) = lu(k,1067) * lu(k,1059)
         lu(k,1068) = lu(k,1068) * lu(k,1059)
         lu(k,1114) = lu(k,1114) - lu(k,1060) * lu(k,1112)
         lu(k,1116) = lu(k,1116) - lu(k,1061) * lu(k,1112)
         lu(k,1117) = lu(k,1117) - lu(k,1062) * lu(k,1112)
         lu(k,1119) = lu(k,1119) - lu(k,1063) * lu(k,1112)
         lu(k,1120) = lu(k,1120) - lu(k,1064) * lu(k,1112)
         lu(k,1122) = lu(k,1122) - lu(k,1065) * lu(k,1112)
         lu(k,1123) = lu(k,1123) - lu(k,1066) * lu(k,1112)
         lu(k,1125) = lu(k,1125) - lu(k,1067) * lu(k,1112)
         lu(k,1126) = lu(k,1126) - lu(k,1068) * lu(k,1112)
         lu(k,1668) = lu(k,1668) - lu(k,1060) * lu(k,1664)
         lu(k,1670) = lu(k,1670) - lu(k,1061) * lu(k,1664)
         lu(k,1683) = lu(k,1683) - lu(k,1062) * lu(k,1664)
         lu(k,1689) = lu(k,1689) - lu(k,1063) * lu(k,1664)
         lu(k,1691) = lu(k,1691) - lu(k,1064) * lu(k,1664)
         lu(k,1694) = lu(k,1694) - lu(k,1065) * lu(k,1664)
         lu(k,1697) = lu(k,1697) - lu(k,1066) * lu(k,1664)
         lu(k,1702) = lu(k,1702) - lu(k,1067) * lu(k,1664)
         lu(k,1703) = lu(k,1703) - lu(k,1068) * lu(k,1664)
         lu(k,1726) = lu(k,1726) - lu(k,1060) * lu(k,1722)
         lu(k,1728) = lu(k,1728) - lu(k,1061) * lu(k,1722)
         lu(k,1741) = lu(k,1741) - lu(k,1062) * lu(k,1722)
         lu(k,1746) = lu(k,1746) - lu(k,1063) * lu(k,1722)
         lu(k,1748) = lu(k,1748) - lu(k,1064) * lu(k,1722)
         lu(k,1751) = lu(k,1751) - lu(k,1065) * lu(k,1722)
         lu(k,1754) = lu(k,1754) - lu(k,1066) * lu(k,1722)
         lu(k,1759) = lu(k,1759) - lu(k,1067) * lu(k,1722)
         lu(k,1760) = lu(k,1760) - lu(k,1068) * lu(k,1722)
         lu(k,1819) = lu(k,1819) - lu(k,1060) * lu(k,1815)
         lu(k,1821) = lu(k,1821) - lu(k,1061) * lu(k,1815)
         lu(k,1833) = lu(k,1833) - lu(k,1062) * lu(k,1815)
         lu(k,1838) = lu(k,1838) - lu(k,1063) * lu(k,1815)
         lu(k,1840) = lu(k,1840) - lu(k,1064) * lu(k,1815)
         lu(k,1843) = lu(k,1843) - lu(k,1065) * lu(k,1815)
         lu(k,1846) = lu(k,1846) - lu(k,1066) * lu(k,1815)
         lu(k,1851) = lu(k,1851) - lu(k,1067) * lu(k,1815)
         lu(k,1852) = lu(k,1852) - lu(k,1068) * lu(k,1815)
         lu(k,1926) = lu(k,1926) - lu(k,1060) * lu(k,1923)
         lu(k,1928) = lu(k,1928) - lu(k,1061) * lu(k,1923)
         lu(k,1939) = lu(k,1939) - lu(k,1062) * lu(k,1923)
         lu(k,1945) = lu(k,1945) - lu(k,1063) * lu(k,1923)
         lu(k,1947) = lu(k,1947) - lu(k,1064) * lu(k,1923)
         lu(k,1950) = lu(k,1950) - lu(k,1065) * lu(k,1923)
         lu(k,1953) = lu(k,1953) - lu(k,1066) * lu(k,1923)
         lu(k,1958) = lu(k,1958) - lu(k,1067) * lu(k,1923)
         lu(k,1959) = lu(k,1959) - lu(k,1068) * lu(k,1923)
         lu(k,2045) = lu(k,2045) - lu(k,1060) * lu(k,2041)
         lu(k,2047) = lu(k,2047) - lu(k,1061) * lu(k,2041)
         lu(k,2058) = lu(k,2058) - lu(k,1062) * lu(k,2041)
         lu(k,2062) = lu(k,2062) - lu(k,1063) * lu(k,2041)
         lu(k,2064) = lu(k,2064) - lu(k,1064) * lu(k,2041)
         lu(k,2067) = lu(k,2067) - lu(k,1065) * lu(k,2041)
         lu(k,2070) = lu(k,2070) - lu(k,1066) * lu(k,2041)
         lu(k,2075) = lu(k,2075) - lu(k,1067) * lu(k,2041)
         lu(k,2076) = lu(k,2076) - lu(k,1068) * lu(k,2041)
         lu(k,2105) = lu(k,2105) - lu(k,1060) * lu(k,2102)
         lu(k,2107) = lu(k,2107) - lu(k,1061) * lu(k,2102)
         lu(k,2118) = lu(k,2118) - lu(k,1062) * lu(k,2102)
         lu(k,2123) = lu(k,2123) - lu(k,1063) * lu(k,2102)
         lu(k,2125) = lu(k,2125) - lu(k,1064) * lu(k,2102)
         lu(k,2128) = lu(k,2128) - lu(k,1065) * lu(k,2102)
         lu(k,2131) = lu(k,2131) - lu(k,1066) * lu(k,2102)
         lu(k,2136) = lu(k,2136) - lu(k,1067) * lu(k,2102)
         lu(k,2137) = lu(k,2137) - lu(k,1068) * lu(k,2102)
                                                                        
         lu(k,1073) = 1._r8 / lu(k,1073)
         lu(k,1074) = lu(k,1074) * lu(k,1073)
         lu(k,1075) = lu(k,1075) * lu(k,1073)
         lu(k,1076) = lu(k,1076) * lu(k,1073)
         lu(k,1077) = lu(k,1077) * lu(k,1073)
         lu(k,1078) = lu(k,1078) * lu(k,1073)
         lu(k,1079) = lu(k,1079) * lu(k,1073)
         lu(k,1080) = lu(k,1080) * lu(k,1073)
         lu(k,1081) = lu(k,1081) * lu(k,1073)
         lu(k,1082) = lu(k,1082) * lu(k,1073)
         lu(k,1083) = lu(k,1083) * lu(k,1073)
         lu(k,1084) = lu(k,1084) * lu(k,1073)
         lu(k,1191) = - lu(k,1074) * lu(k,1189)
         lu(k,1195) = lu(k,1195) - lu(k,1075) * lu(k,1189)
         lu(k,1197) = lu(k,1197) - lu(k,1076) * lu(k,1189)
         lu(k,1198) = lu(k,1198) - lu(k,1077) * lu(k,1189)
         lu(k,1199) = lu(k,1199) - lu(k,1078) * lu(k,1189)
         lu(k,1200) = lu(k,1200) - lu(k,1079) * lu(k,1189)
         lu(k,1202) = lu(k,1202) - lu(k,1080) * lu(k,1189)
         lu(k,1203) = lu(k,1203) - lu(k,1081) * lu(k,1189)
         lu(k,1204) = lu(k,1204) - lu(k,1082) * lu(k,1189)
         lu(k,1205) = lu(k,1205) - lu(k,1083) * lu(k,1189)
         lu(k,1206) = - lu(k,1084) * lu(k,1189)
         lu(k,1375) = lu(k,1375) - lu(k,1074) * lu(k,1373)
         lu(k,1380) = lu(k,1380) - lu(k,1075) * lu(k,1373)
         lu(k,1386) = lu(k,1386) - lu(k,1076) * lu(k,1373)
         lu(k,1388) = - lu(k,1077) * lu(k,1373)
         lu(k,1389) = lu(k,1389) - lu(k,1078) * lu(k,1373)
         lu(k,1390) = lu(k,1390) - lu(k,1079) * lu(k,1373)
         lu(k,1392) = lu(k,1392) - lu(k,1080) * lu(k,1373)
         lu(k,1393) = lu(k,1393) - lu(k,1081) * lu(k,1373)
         lu(k,1394) = lu(k,1394) - lu(k,1082) * lu(k,1373)
         lu(k,1396) = lu(k,1396) - lu(k,1083) * lu(k,1373)
         lu(k,1397) = lu(k,1397) - lu(k,1084) * lu(k,1373)
         lu(k,1668) = lu(k,1668) - lu(k,1074) * lu(k,1665)
         lu(k,1676) = lu(k,1676) - lu(k,1075) * lu(k,1665)
         lu(k,1683) = lu(k,1683) - lu(k,1076) * lu(k,1665)
         lu(k,1687) = lu(k,1687) - lu(k,1077) * lu(k,1665)
         lu(k,1689) = lu(k,1689) - lu(k,1078) * lu(k,1665)
         lu(k,1691) = lu(k,1691) - lu(k,1079) * lu(k,1665)
         lu(k,1693) = lu(k,1693) - lu(k,1080) * lu(k,1665)
         lu(k,1694) = lu(k,1694) - lu(k,1081) * lu(k,1665)
         lu(k,1697) = lu(k,1697) - lu(k,1082) * lu(k,1665)
         lu(k,1700) = lu(k,1700) - lu(k,1083) * lu(k,1665)
         lu(k,1702) = lu(k,1702) - lu(k,1084) * lu(k,1665)
         lu(k,1726) = lu(k,1726) - lu(k,1074) * lu(k,1723)
         lu(k,1734) = lu(k,1734) - lu(k,1075) * lu(k,1723)
         lu(k,1741) = lu(k,1741) - lu(k,1076) * lu(k,1723)
         lu(k,1744) = lu(k,1744) - lu(k,1077) * lu(k,1723)
         lu(k,1746) = lu(k,1746) - lu(k,1078) * lu(k,1723)
         lu(k,1748) = lu(k,1748) - lu(k,1079) * lu(k,1723)
         lu(k,1750) = lu(k,1750) - lu(k,1080) * lu(k,1723)
         lu(k,1751) = lu(k,1751) - lu(k,1081) * lu(k,1723)
         lu(k,1754) = lu(k,1754) - lu(k,1082) * lu(k,1723)
         lu(k,1757) = lu(k,1757) - lu(k,1083) * lu(k,1723)
         lu(k,1759) = lu(k,1759) - lu(k,1084) * lu(k,1723)
         lu(k,1819) = lu(k,1819) - lu(k,1074) * lu(k,1816)
         lu(k,1827) = lu(k,1827) - lu(k,1075) * lu(k,1816)
         lu(k,1833) = lu(k,1833) - lu(k,1076) * lu(k,1816)
         lu(k,1836) = lu(k,1836) - lu(k,1077) * lu(k,1816)
         lu(k,1838) = lu(k,1838) - lu(k,1078) * lu(k,1816)
         lu(k,1840) = lu(k,1840) - lu(k,1079) * lu(k,1816)
         lu(k,1842) = lu(k,1842) - lu(k,1080) * lu(k,1816)
         lu(k,1843) = lu(k,1843) - lu(k,1081) * lu(k,1816)
         lu(k,1846) = lu(k,1846) - lu(k,1082) * lu(k,1816)
         lu(k,1849) = lu(k,1849) - lu(k,1083) * lu(k,1816)
         lu(k,1851) = lu(k,1851) - lu(k,1084) * lu(k,1816)
         lu(k,2045) = lu(k,2045) - lu(k,1074) * lu(k,2042)
         lu(k,2052) = lu(k,2052) - lu(k,1075) * lu(k,2042)
         lu(k,2058) = lu(k,2058) - lu(k,1076) * lu(k,2042)
         lu(k,2060) = lu(k,2060) - lu(k,1077) * lu(k,2042)
         lu(k,2062) = lu(k,2062) - lu(k,1078) * lu(k,2042)
         lu(k,2064) = lu(k,2064) - lu(k,1079) * lu(k,2042)
         lu(k,2066) = lu(k,2066) - lu(k,1080) * lu(k,2042)
         lu(k,2067) = lu(k,2067) - lu(k,1081) * lu(k,2042)
         lu(k,2070) = lu(k,2070) - lu(k,1082) * lu(k,2042)
         lu(k,2073) = lu(k,2073) - lu(k,1083) * lu(k,2042)
         lu(k,2075) = lu(k,2075) - lu(k,1084) * lu(k,2042)
                                                                        
      end do
                                                                        
      end subroutine lu_fac22
                                                                        
      subroutine lu_fac23( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,1085) = 1._r8 / lu(k,1085)
         lu(k,1086) = lu(k,1086) * lu(k,1085)
         lu(k,1087) = lu(k,1087) * lu(k,1085)
         lu(k,1088) = lu(k,1088) * lu(k,1085)
         lu(k,1089) = lu(k,1089) * lu(k,1085)
         lu(k,1090) = lu(k,1090) * lu(k,1085)
         lu(k,1091) = lu(k,1091) * lu(k,1085)
         lu(k,1092) = lu(k,1092) * lu(k,1085)
         lu(k,1093) = lu(k,1093) * lu(k,1085)
         lu(k,1150) = lu(k,1150) - lu(k,1086) * lu(k,1148)
         lu(k,1155) = lu(k,1155) - lu(k,1087) * lu(k,1148)
         lu(k,1156) = lu(k,1156) - lu(k,1088) * lu(k,1148)
         lu(k,1158) = lu(k,1158) - lu(k,1089) * lu(k,1148)
         lu(k,1159) = - lu(k,1090) * lu(k,1148)
         lu(k,1161) = lu(k,1161) - lu(k,1091) * lu(k,1148)
         lu(k,1162) = lu(k,1162) - lu(k,1092) * lu(k,1148)
         lu(k,1165) = lu(k,1165) - lu(k,1093) * lu(k,1148)
         lu(k,1328) = lu(k,1328) - lu(k,1086) * lu(k,1327)
         lu(k,1334) = lu(k,1334) - lu(k,1087) * lu(k,1327)
         lu(k,1336) = - lu(k,1088) * lu(k,1327)
         lu(k,1338) = lu(k,1338) - lu(k,1089) * lu(k,1327)
         lu(k,1339) = lu(k,1339) - lu(k,1090) * lu(k,1327)
         lu(k,1341) = lu(k,1341) - lu(k,1091) * lu(k,1327)
         lu(k,1342) = lu(k,1342) - lu(k,1092) * lu(k,1327)
         lu(k,1346) = lu(k,1346) - lu(k,1093) * lu(k,1327)
         lu(k,1668) = lu(k,1668) - lu(k,1086) * lu(k,1666)
         lu(k,1683) = lu(k,1683) - lu(k,1087) * lu(k,1666)
         lu(k,1687) = lu(k,1687) - lu(k,1088) * lu(k,1666)
         lu(k,1691) = lu(k,1691) - lu(k,1089) * lu(k,1666)
         lu(k,1692) = lu(k,1692) - lu(k,1090) * lu(k,1666)
         lu(k,1694) = lu(k,1694) - lu(k,1091) * lu(k,1666)
         lu(k,1697) = lu(k,1697) - lu(k,1092) * lu(k,1666)
         lu(k,1703) = lu(k,1703) - lu(k,1093) * lu(k,1666)
         lu(k,1726) = lu(k,1726) - lu(k,1086) * lu(k,1724)
         lu(k,1741) = lu(k,1741) - lu(k,1087) * lu(k,1724)
         lu(k,1744) = lu(k,1744) - lu(k,1088) * lu(k,1724)
         lu(k,1748) = lu(k,1748) - lu(k,1089) * lu(k,1724)
         lu(k,1749) = lu(k,1749) - lu(k,1090) * lu(k,1724)
         lu(k,1751) = lu(k,1751) - lu(k,1091) * lu(k,1724)
         lu(k,1754) = lu(k,1754) - lu(k,1092) * lu(k,1724)
         lu(k,1760) = lu(k,1760) - lu(k,1093) * lu(k,1724)
         lu(k,1819) = lu(k,1819) - lu(k,1086) * lu(k,1817)
         lu(k,1833) = lu(k,1833) - lu(k,1087) * lu(k,1817)
         lu(k,1836) = lu(k,1836) - lu(k,1088) * lu(k,1817)
         lu(k,1840) = lu(k,1840) - lu(k,1089) * lu(k,1817)
         lu(k,1841) = lu(k,1841) - lu(k,1090) * lu(k,1817)
         lu(k,1843) = lu(k,1843) - lu(k,1091) * lu(k,1817)
         lu(k,1846) = lu(k,1846) - lu(k,1092) * lu(k,1817)
         lu(k,1852) = lu(k,1852) - lu(k,1093) * lu(k,1817)
         lu(k,1926) = lu(k,1926) - lu(k,1086) * lu(k,1924)
         lu(k,1939) = lu(k,1939) - lu(k,1087) * lu(k,1924)
         lu(k,1943) = lu(k,1943) - lu(k,1088) * lu(k,1924)
         lu(k,1947) = lu(k,1947) - lu(k,1089) * lu(k,1924)
         lu(k,1948) = lu(k,1948) - lu(k,1090) * lu(k,1924)
         lu(k,1950) = lu(k,1950) - lu(k,1091) * lu(k,1924)
         lu(k,1953) = lu(k,1953) - lu(k,1092) * lu(k,1924)
         lu(k,1959) = lu(k,1959) - lu(k,1093) * lu(k,1924)
         lu(k,2003) = lu(k,2003) - lu(k,1086) * lu(k,2001)
         lu(k,2004) = - lu(k,1087) * lu(k,2001)
         lu(k,2008) = lu(k,2008) - lu(k,1088) * lu(k,2001)
         lu(k,2012) = lu(k,2012) - lu(k,1089) * lu(k,2001)
         lu(k,2013) = lu(k,2013) - lu(k,1090) * lu(k,2001)
         lu(k,2015) = lu(k,2015) - lu(k,1091) * lu(k,2001)
         lu(k,2018) = lu(k,2018) - lu(k,1092) * lu(k,2001)
         lu(k,2024) = lu(k,2024) - lu(k,1093) * lu(k,2001)
         lu(k,2045) = lu(k,2045) - lu(k,1086) * lu(k,2043)
         lu(k,2058) = lu(k,2058) - lu(k,1087) * lu(k,2043)
         lu(k,2060) = lu(k,2060) - lu(k,1088) * lu(k,2043)
         lu(k,2064) = lu(k,2064) - lu(k,1089) * lu(k,2043)
         lu(k,2065) = lu(k,2065) - lu(k,1090) * lu(k,2043)
         lu(k,2067) = lu(k,2067) - lu(k,1091) * lu(k,2043)
         lu(k,2070) = lu(k,2070) - lu(k,1092) * lu(k,2043)
         lu(k,2076) = lu(k,2076) - lu(k,1093) * lu(k,2043)
         lu(k,2105) = lu(k,2105) - lu(k,1086) * lu(k,2103)
         lu(k,2118) = lu(k,2118) - lu(k,1087) * lu(k,2103)
         lu(k,2121) = lu(k,2121) - lu(k,1088) * lu(k,2103)
         lu(k,2125) = lu(k,2125) - lu(k,1089) * lu(k,2103)
         lu(k,2126) = lu(k,2126) - lu(k,1090) * lu(k,2103)
         lu(k,2128) = lu(k,2128) - lu(k,1091) * lu(k,2103)
         lu(k,2131) = lu(k,2131) - lu(k,1092) * lu(k,2103)
         lu(k,2137) = lu(k,2137) - lu(k,1093) * lu(k,2103)
                                                                        
         lu(k,1096) = 1._r8 / lu(k,1096)
         lu(k,1097) = lu(k,1097) * lu(k,1096)
         lu(k,1098) = lu(k,1098) * lu(k,1096)
         lu(k,1099) = lu(k,1099) * lu(k,1096)
         lu(k,1100) = lu(k,1100) * lu(k,1096)
         lu(k,1101) = lu(k,1101) * lu(k,1096)
         lu(k,1114) = lu(k,1114) - lu(k,1097) * lu(k,1113)
         lu(k,1119) = lu(k,1119) - lu(k,1098) * lu(k,1113)
         lu(k,1120) = lu(k,1120) - lu(k,1099) * lu(k,1113)
         lu(k,1122) = lu(k,1122) - lu(k,1100) * lu(k,1113)
         lu(k,1125) = lu(k,1125) - lu(k,1101) * lu(k,1113)
         lu(k,1150) = lu(k,1150) - lu(k,1097) * lu(k,1149)
         lu(k,1157) = lu(k,1157) - lu(k,1098) * lu(k,1149)
         lu(k,1158) = lu(k,1158) - lu(k,1099) * lu(k,1149)
         lu(k,1161) = lu(k,1161) - lu(k,1100) * lu(k,1149)
         lu(k,1164) = - lu(k,1101) * lu(k,1149)
         lu(k,1171) = lu(k,1171) - lu(k,1097) * lu(k,1170)
         lu(k,1176) = lu(k,1176) - lu(k,1098) * lu(k,1170)
         lu(k,1177) = lu(k,1177) - lu(k,1099) * lu(k,1170)
         lu(k,1180) = lu(k,1180) - lu(k,1100) * lu(k,1170)
         lu(k,1183) = lu(k,1183) - lu(k,1101) * lu(k,1170)
         lu(k,1191) = lu(k,1191) - lu(k,1097) * lu(k,1190)
         lu(k,1199) = lu(k,1199) - lu(k,1098) * lu(k,1190)
         lu(k,1200) = lu(k,1200) - lu(k,1099) * lu(k,1190)
         lu(k,1203) = lu(k,1203) - lu(k,1100) * lu(k,1190)
         lu(k,1206) = lu(k,1206) - lu(k,1101) * lu(k,1190)
         lu(k,1250) = lu(k,1250) - lu(k,1097) * lu(k,1249)
         lu(k,1262) = lu(k,1262) - lu(k,1098) * lu(k,1249)
         lu(k,1263) = lu(k,1263) - lu(k,1099) * lu(k,1249)
         lu(k,1266) = lu(k,1266) - lu(k,1100) * lu(k,1249)
         lu(k,1270) = lu(k,1270) - lu(k,1101) * lu(k,1249)
         lu(k,1282) = lu(k,1282) - lu(k,1097) * lu(k,1281)
         lu(k,1294) = lu(k,1294) - lu(k,1098) * lu(k,1281)
         lu(k,1295) = lu(k,1295) - lu(k,1099) * lu(k,1281)
         lu(k,1298) = lu(k,1298) - lu(k,1100) * lu(k,1281)
         lu(k,1302) = lu(k,1302) - lu(k,1101) * lu(k,1281)
         lu(k,1307) = lu(k,1307) - lu(k,1097) * lu(k,1306)
         lu(k,1315) = lu(k,1315) - lu(k,1098) * lu(k,1306)
         lu(k,1316) = lu(k,1316) - lu(k,1099) * lu(k,1306)
         lu(k,1319) = lu(k,1319) - lu(k,1100) * lu(k,1306)
         lu(k,1322) = - lu(k,1101) * lu(k,1306)
         lu(k,1375) = lu(k,1375) - lu(k,1097) * lu(k,1374)
         lu(k,1389) = lu(k,1389) - lu(k,1098) * lu(k,1374)
         lu(k,1390) = lu(k,1390) - lu(k,1099) * lu(k,1374)
         lu(k,1393) = lu(k,1393) - lu(k,1100) * lu(k,1374)
         lu(k,1397) = lu(k,1397) - lu(k,1101) * lu(k,1374)
         lu(k,1668) = lu(k,1668) - lu(k,1097) * lu(k,1667)
         lu(k,1689) = lu(k,1689) - lu(k,1098) * lu(k,1667)
         lu(k,1691) = lu(k,1691) - lu(k,1099) * lu(k,1667)
         lu(k,1694) = lu(k,1694) - lu(k,1100) * lu(k,1667)
         lu(k,1702) = lu(k,1702) - lu(k,1101) * lu(k,1667)
         lu(k,1726) = lu(k,1726) - lu(k,1097) * lu(k,1725)
         lu(k,1746) = lu(k,1746) - lu(k,1098) * lu(k,1725)
         lu(k,1748) = lu(k,1748) - lu(k,1099) * lu(k,1725)
         lu(k,1751) = lu(k,1751) - lu(k,1100) * lu(k,1725)
         lu(k,1759) = lu(k,1759) - lu(k,1101) * lu(k,1725)
         lu(k,1819) = lu(k,1819) - lu(k,1097) * lu(k,1818)
         lu(k,1838) = lu(k,1838) - lu(k,1098) * lu(k,1818)
         lu(k,1840) = lu(k,1840) - lu(k,1099) * lu(k,1818)
         lu(k,1843) = lu(k,1843) - lu(k,1100) * lu(k,1818)
         lu(k,1851) = lu(k,1851) - lu(k,1101) * lu(k,1818)
         lu(k,1926) = lu(k,1926) - lu(k,1097) * lu(k,1925)
         lu(k,1945) = lu(k,1945) - lu(k,1098) * lu(k,1925)
         lu(k,1947) = lu(k,1947) - lu(k,1099) * lu(k,1925)
         lu(k,1950) = lu(k,1950) - lu(k,1100) * lu(k,1925)
         lu(k,1958) = lu(k,1958) - lu(k,1101) * lu(k,1925)
         lu(k,2003) = lu(k,2003) - lu(k,1097) * lu(k,2002)
         lu(k,2010) = lu(k,2010) - lu(k,1098) * lu(k,2002)
         lu(k,2012) = lu(k,2012) - lu(k,1099) * lu(k,2002)
         lu(k,2015) = lu(k,2015) - lu(k,1100) * lu(k,2002)
         lu(k,2023) = lu(k,2023) - lu(k,1101) * lu(k,2002)
         lu(k,2045) = lu(k,2045) - lu(k,1097) * lu(k,2044)
         lu(k,2062) = lu(k,2062) - lu(k,1098) * lu(k,2044)
         lu(k,2064) = lu(k,2064) - lu(k,1099) * lu(k,2044)
         lu(k,2067) = lu(k,2067) - lu(k,1100) * lu(k,2044)
         lu(k,2075) = lu(k,2075) - lu(k,1101) * lu(k,2044)
         lu(k,2105) = lu(k,2105) - lu(k,1097) * lu(k,2104)
         lu(k,2123) = lu(k,2123) - lu(k,1098) * lu(k,2104)
         lu(k,2125) = lu(k,2125) - lu(k,1099) * lu(k,2104)
         lu(k,2128) = lu(k,2128) - lu(k,1100) * lu(k,2104)
         lu(k,2136) = lu(k,2136) - lu(k,1101) * lu(k,2104)
                                                                        
         lu(k,1103) = 1._r8 / lu(k,1103)
         lu(k,1104) = lu(k,1104) * lu(k,1103)
         lu(k,1105) = lu(k,1105) * lu(k,1103)
         lu(k,1106) = lu(k,1106) * lu(k,1103)
         lu(k,1120) = lu(k,1120) - lu(k,1104) * lu(k,1114)
         lu(k,1122) = lu(k,1122) - lu(k,1105) * lu(k,1114)
         lu(k,1125) = lu(k,1125) - lu(k,1106) * lu(k,1114)
         lu(k,1158) = lu(k,1158) - lu(k,1104) * lu(k,1150)
         lu(k,1161) = lu(k,1161) - lu(k,1105) * lu(k,1150)
         lu(k,1164) = lu(k,1164) - lu(k,1106) * lu(k,1150)
         lu(k,1177) = lu(k,1177) - lu(k,1104) * lu(k,1171)
         lu(k,1180) = lu(k,1180) - lu(k,1105) * lu(k,1171)
         lu(k,1183) = lu(k,1183) - lu(k,1106) * lu(k,1171)
         lu(k,1200) = lu(k,1200) - lu(k,1104) * lu(k,1191)
         lu(k,1203) = lu(k,1203) - lu(k,1105) * lu(k,1191)
         lu(k,1206) = lu(k,1206) - lu(k,1106) * lu(k,1191)
         lu(k,1215) = lu(k,1215) - lu(k,1104) * lu(k,1208)
         lu(k,1216) = lu(k,1216) - lu(k,1105) * lu(k,1208)
         lu(k,1218) = lu(k,1218) - lu(k,1106) * lu(k,1208)
         lu(k,1224) = lu(k,1224) - lu(k,1104) * lu(k,1220)
         lu(k,1226) = lu(k,1226) - lu(k,1105) * lu(k,1220)
         lu(k,1227) = - lu(k,1106) * lu(k,1220)
         lu(k,1263) = lu(k,1263) - lu(k,1104) * lu(k,1250)
         lu(k,1266) = lu(k,1266) - lu(k,1105) * lu(k,1250)
         lu(k,1270) = lu(k,1270) - lu(k,1106) * lu(k,1250)
         lu(k,1295) = lu(k,1295) - lu(k,1104) * lu(k,1282)
         lu(k,1298) = lu(k,1298) - lu(k,1105) * lu(k,1282)
         lu(k,1302) = lu(k,1302) - lu(k,1106) * lu(k,1282)
         lu(k,1316) = lu(k,1316) - lu(k,1104) * lu(k,1307)
         lu(k,1319) = lu(k,1319) - lu(k,1105) * lu(k,1307)
         lu(k,1322) = lu(k,1322) - lu(k,1106) * lu(k,1307)
         lu(k,1338) = lu(k,1338) - lu(k,1104) * lu(k,1328)
         lu(k,1341) = lu(k,1341) - lu(k,1105) * lu(k,1328)
         lu(k,1345) = lu(k,1345) - lu(k,1106) * lu(k,1328)
         lu(k,1358) = lu(k,1358) - lu(k,1104) * lu(k,1351)
         lu(k,1361) = lu(k,1361) - lu(k,1105) * lu(k,1351)
         lu(k,1365) = lu(k,1365) - lu(k,1106) * lu(k,1351)
         lu(k,1390) = lu(k,1390) - lu(k,1104) * lu(k,1375)
         lu(k,1393) = lu(k,1393) - lu(k,1105) * lu(k,1375)
         lu(k,1397) = lu(k,1397) - lu(k,1106) * lu(k,1375)
         lu(k,1419) = lu(k,1419) - lu(k,1104) * lu(k,1414)
         lu(k,1420) = lu(k,1420) - lu(k,1105) * lu(k,1414)
         lu(k,1423) = lu(k,1423) - lu(k,1106) * lu(k,1414)
         lu(k,1435) = lu(k,1435) - lu(k,1104) * lu(k,1428)
         lu(k,1437) = lu(k,1437) - lu(k,1105) * lu(k,1428)
         lu(k,1441) = lu(k,1441) - lu(k,1106) * lu(k,1428)
         lu(k,1487) = lu(k,1487) - lu(k,1104) * lu(k,1479)
         lu(k,1490) = lu(k,1490) - lu(k,1105) * lu(k,1479)
         lu(k,1497) = lu(k,1497) - lu(k,1106) * lu(k,1479)
         lu(k,1691) = lu(k,1691) - lu(k,1104) * lu(k,1668)
         lu(k,1694) = lu(k,1694) - lu(k,1105) * lu(k,1668)
         lu(k,1702) = lu(k,1702) - lu(k,1106) * lu(k,1668)
         lu(k,1748) = lu(k,1748) - lu(k,1104) * lu(k,1726)
         lu(k,1751) = lu(k,1751) - lu(k,1105) * lu(k,1726)
         lu(k,1759) = lu(k,1759) - lu(k,1106) * lu(k,1726)
         lu(k,1840) = lu(k,1840) - lu(k,1104) * lu(k,1819)
         lu(k,1843) = lu(k,1843) - lu(k,1105) * lu(k,1819)
         lu(k,1851) = lu(k,1851) - lu(k,1106) * lu(k,1819)
         lu(k,1947) = lu(k,1947) - lu(k,1104) * lu(k,1926)
         lu(k,1950) = lu(k,1950) - lu(k,1105) * lu(k,1926)
         lu(k,1958) = lu(k,1958) - lu(k,1106) * lu(k,1926)
         lu(k,2012) = lu(k,2012) - lu(k,1104) * lu(k,2003)
         lu(k,2015) = lu(k,2015) - lu(k,1105) * lu(k,2003)
         lu(k,2023) = lu(k,2023) - lu(k,1106) * lu(k,2003)
         lu(k,2064) = lu(k,2064) - lu(k,1104) * lu(k,2045)
         lu(k,2067) = lu(k,2067) - lu(k,1105) * lu(k,2045)
         lu(k,2075) = lu(k,2075) - lu(k,1106) * lu(k,2045)
         lu(k,2125) = lu(k,2125) - lu(k,1104) * lu(k,2105)
         lu(k,2128) = lu(k,2128) - lu(k,1105) * lu(k,2105)
         lu(k,2136) = lu(k,2136) - lu(k,1106) * lu(k,2105)
         lu(k,2192) = lu(k,2192) - lu(k,1104) * lu(k,2180)
         lu(k,2195) = lu(k,2195) - lu(k,1105) * lu(k,2180)
         lu(k,2203) = lu(k,2203) - lu(k,1106) * lu(k,2180)
         lu(k,2247) = lu(k,2247) - lu(k,1104) * lu(k,2238)
         lu(k,2250) = lu(k,2250) - lu(k,1105) * lu(k,2238)
         lu(k,2258) = lu(k,2258) - lu(k,1106) * lu(k,2238)
                                                                        
         lu(k,1115) = 1._r8 / lu(k,1115)
         lu(k,1116) = lu(k,1116) * lu(k,1115)
         lu(k,1117) = lu(k,1117) * lu(k,1115)
         lu(k,1118) = lu(k,1118) * lu(k,1115)
         lu(k,1119) = lu(k,1119) * lu(k,1115)
         lu(k,1120) = lu(k,1120) * lu(k,1115)
         lu(k,1121) = lu(k,1121) * lu(k,1115)
         lu(k,1122) = lu(k,1122) * lu(k,1115)
         lu(k,1123) = lu(k,1123) * lu(k,1115)
         lu(k,1124) = lu(k,1124) * lu(k,1115)
         lu(k,1125) = lu(k,1125) * lu(k,1115)
         lu(k,1126) = lu(k,1126) * lu(k,1115)
         lu(k,1670) = lu(k,1670) - lu(k,1116) * lu(k,1669)
         lu(k,1683) = lu(k,1683) - lu(k,1117) * lu(k,1669)
         lu(k,1687) = lu(k,1687) - lu(k,1118) * lu(k,1669)
         lu(k,1689) = lu(k,1689) - lu(k,1119) * lu(k,1669)
         lu(k,1691) = lu(k,1691) - lu(k,1120) * lu(k,1669)
         lu(k,1693) = lu(k,1693) - lu(k,1121) * lu(k,1669)
         lu(k,1694) = lu(k,1694) - lu(k,1122) * lu(k,1669)
         lu(k,1697) = lu(k,1697) - lu(k,1123) * lu(k,1669)
         lu(k,1700) = lu(k,1700) - lu(k,1124) * lu(k,1669)
         lu(k,1702) = lu(k,1702) - lu(k,1125) * lu(k,1669)
         lu(k,1703) = lu(k,1703) - lu(k,1126) * lu(k,1669)
         lu(k,1728) = lu(k,1728) - lu(k,1116) * lu(k,1727)
         lu(k,1741) = lu(k,1741) - lu(k,1117) * lu(k,1727)
         lu(k,1744) = lu(k,1744) - lu(k,1118) * lu(k,1727)
         lu(k,1746) = lu(k,1746) - lu(k,1119) * lu(k,1727)
         lu(k,1748) = lu(k,1748) - lu(k,1120) * lu(k,1727)
         lu(k,1750) = lu(k,1750) - lu(k,1121) * lu(k,1727)
         lu(k,1751) = lu(k,1751) - lu(k,1122) * lu(k,1727)
         lu(k,1754) = lu(k,1754) - lu(k,1123) * lu(k,1727)
         lu(k,1757) = lu(k,1757) - lu(k,1124) * lu(k,1727)
         lu(k,1759) = lu(k,1759) - lu(k,1125) * lu(k,1727)
         lu(k,1760) = lu(k,1760) - lu(k,1126) * lu(k,1727)
         lu(k,1821) = lu(k,1821) - lu(k,1116) * lu(k,1820)
         lu(k,1833) = lu(k,1833) - lu(k,1117) * lu(k,1820)
         lu(k,1836) = lu(k,1836) - lu(k,1118) * lu(k,1820)
         lu(k,1838) = lu(k,1838) - lu(k,1119) * lu(k,1820)
         lu(k,1840) = lu(k,1840) - lu(k,1120) * lu(k,1820)
         lu(k,1842) = lu(k,1842) - lu(k,1121) * lu(k,1820)
         lu(k,1843) = lu(k,1843) - lu(k,1122) * lu(k,1820)
         lu(k,1846) = lu(k,1846) - lu(k,1123) * lu(k,1820)
         lu(k,1849) = lu(k,1849) - lu(k,1124) * lu(k,1820)
         lu(k,1851) = lu(k,1851) - lu(k,1125) * lu(k,1820)
         lu(k,1852) = lu(k,1852) - lu(k,1126) * lu(k,1820)
         lu(k,1928) = lu(k,1928) - lu(k,1116) * lu(k,1927)
         lu(k,1939) = lu(k,1939) - lu(k,1117) * lu(k,1927)
         lu(k,1943) = lu(k,1943) - lu(k,1118) * lu(k,1927)
         lu(k,1945) = lu(k,1945) - lu(k,1119) * lu(k,1927)
         lu(k,1947) = lu(k,1947) - lu(k,1120) * lu(k,1927)
         lu(k,1949) = lu(k,1949) - lu(k,1121) * lu(k,1927)
         lu(k,1950) = lu(k,1950) - lu(k,1122) * lu(k,1927)
         lu(k,1953) = lu(k,1953) - lu(k,1123) * lu(k,1927)
         lu(k,1956) = lu(k,1956) - lu(k,1124) * lu(k,1927)
         lu(k,1958) = lu(k,1958) - lu(k,1125) * lu(k,1927)
         lu(k,1959) = lu(k,1959) - lu(k,1126) * lu(k,1927)
         lu(k,2047) = lu(k,2047) - lu(k,1116) * lu(k,2046)
         lu(k,2058) = lu(k,2058) - lu(k,1117) * lu(k,2046)
         lu(k,2060) = lu(k,2060) - lu(k,1118) * lu(k,2046)
         lu(k,2062) = lu(k,2062) - lu(k,1119) * lu(k,2046)
         lu(k,2064) = lu(k,2064) - lu(k,1120) * lu(k,2046)
         lu(k,2066) = lu(k,2066) - lu(k,1121) * lu(k,2046)
         lu(k,2067) = lu(k,2067) - lu(k,1122) * lu(k,2046)
         lu(k,2070) = lu(k,2070) - lu(k,1123) * lu(k,2046)
         lu(k,2073) = lu(k,2073) - lu(k,1124) * lu(k,2046)
         lu(k,2075) = lu(k,2075) - lu(k,1125) * lu(k,2046)
         lu(k,2076) = lu(k,2076) - lu(k,1126) * lu(k,2046)
         lu(k,2107) = lu(k,2107) - lu(k,1116) * lu(k,2106)
         lu(k,2118) = lu(k,2118) - lu(k,1117) * lu(k,2106)
         lu(k,2121) = lu(k,2121) - lu(k,1118) * lu(k,2106)
         lu(k,2123) = lu(k,2123) - lu(k,1119) * lu(k,2106)
         lu(k,2125) = lu(k,2125) - lu(k,1120) * lu(k,2106)
         lu(k,2127) = lu(k,2127) - lu(k,1121) * lu(k,2106)
         lu(k,2128) = lu(k,2128) - lu(k,1122) * lu(k,2106)
         lu(k,2131) = lu(k,2131) - lu(k,1123) * lu(k,2106)
         lu(k,2134) = lu(k,2134) - lu(k,1124) * lu(k,2106)
         lu(k,2136) = lu(k,2136) - lu(k,1125) * lu(k,2106)
         lu(k,2137) = lu(k,2137) - lu(k,1126) * lu(k,2106)
                                                                        
         lu(k,1129) = 1._r8 / lu(k,1129)
         lu(k,1130) = lu(k,1130) * lu(k,1129)
         lu(k,1131) = lu(k,1131) * lu(k,1129)
         lu(k,1132) = lu(k,1132) * lu(k,1129)
         lu(k,1133) = lu(k,1133) * lu(k,1129)
         lu(k,1134) = lu(k,1134) * lu(k,1129)
         lu(k,1135) = lu(k,1135) * lu(k,1129)
         lu(k,1136) = lu(k,1136) * lu(k,1129)
         lu(k,1137) = lu(k,1137) * lu(k,1129)
         lu(k,1138) = lu(k,1138) * lu(k,1129)
         lu(k,1139) = lu(k,1139) * lu(k,1129)
         lu(k,1152) = lu(k,1152) - lu(k,1130) * lu(k,1151)
         lu(k,1154) = - lu(k,1131) * lu(k,1151)
         lu(k,1155) = lu(k,1155) - lu(k,1132) * lu(k,1151)
         lu(k,1157) = lu(k,1157) - lu(k,1133) * lu(k,1151)
         lu(k,1158) = lu(k,1158) - lu(k,1134) * lu(k,1151)
         lu(k,1160) = - lu(k,1135) * lu(k,1151)
         lu(k,1161) = lu(k,1161) - lu(k,1136) * lu(k,1151)
         lu(k,1162) = lu(k,1162) - lu(k,1137) * lu(k,1151)
         lu(k,1163) = lu(k,1163) - lu(k,1138) * lu(k,1151)
         lu(k,1165) = lu(k,1165) - lu(k,1139) * lu(k,1151)
         lu(k,1671) = lu(k,1671) - lu(k,1130) * lu(k,1670)
         lu(k,1676) = lu(k,1676) - lu(k,1131) * lu(k,1670)
         lu(k,1683) = lu(k,1683) - lu(k,1132) * lu(k,1670)
         lu(k,1689) = lu(k,1689) - lu(k,1133) * lu(k,1670)
         lu(k,1691) = lu(k,1691) - lu(k,1134) * lu(k,1670)
         lu(k,1693) = lu(k,1693) - lu(k,1135) * lu(k,1670)
         lu(k,1694) = lu(k,1694) - lu(k,1136) * lu(k,1670)
         lu(k,1697) = lu(k,1697) - lu(k,1137) * lu(k,1670)
         lu(k,1700) = lu(k,1700) - lu(k,1138) * lu(k,1670)
         lu(k,1703) = lu(k,1703) - lu(k,1139) * lu(k,1670)
         lu(k,1729) = lu(k,1729) - lu(k,1130) * lu(k,1728)
         lu(k,1734) = lu(k,1734) - lu(k,1131) * lu(k,1728)
         lu(k,1741) = lu(k,1741) - lu(k,1132) * lu(k,1728)
         lu(k,1746) = lu(k,1746) - lu(k,1133) * lu(k,1728)
         lu(k,1748) = lu(k,1748) - lu(k,1134) * lu(k,1728)
         lu(k,1750) = lu(k,1750) - lu(k,1135) * lu(k,1728)
         lu(k,1751) = lu(k,1751) - lu(k,1136) * lu(k,1728)
         lu(k,1754) = lu(k,1754) - lu(k,1137) * lu(k,1728)
         lu(k,1757) = lu(k,1757) - lu(k,1138) * lu(k,1728)
         lu(k,1760) = lu(k,1760) - lu(k,1139) * lu(k,1728)
         lu(k,1822) = lu(k,1822) - lu(k,1130) * lu(k,1821)
         lu(k,1827) = lu(k,1827) - lu(k,1131) * lu(k,1821)
         lu(k,1833) = lu(k,1833) - lu(k,1132) * lu(k,1821)
         lu(k,1838) = lu(k,1838) - lu(k,1133) * lu(k,1821)
         lu(k,1840) = lu(k,1840) - lu(k,1134) * lu(k,1821)
         lu(k,1842) = lu(k,1842) - lu(k,1135) * lu(k,1821)
         lu(k,1843) = lu(k,1843) - lu(k,1136) * lu(k,1821)
         lu(k,1846) = lu(k,1846) - lu(k,1137) * lu(k,1821)
         lu(k,1849) = lu(k,1849) - lu(k,1138) * lu(k,1821)
         lu(k,1852) = lu(k,1852) - lu(k,1139) * lu(k,1821)
         lu(k,1929) = lu(k,1929) - lu(k,1130) * lu(k,1928)
         lu(k,1933) = lu(k,1933) - lu(k,1131) * lu(k,1928)
         lu(k,1939) = lu(k,1939) - lu(k,1132) * lu(k,1928)
         lu(k,1945) = lu(k,1945) - lu(k,1133) * lu(k,1928)
         lu(k,1947) = lu(k,1947) - lu(k,1134) * lu(k,1928)
         lu(k,1949) = lu(k,1949) - lu(k,1135) * lu(k,1928)
         lu(k,1950) = lu(k,1950) - lu(k,1136) * lu(k,1928)
         lu(k,1953) = lu(k,1953) - lu(k,1137) * lu(k,1928)
         lu(k,1956) = lu(k,1956) - lu(k,1138) * lu(k,1928)
         lu(k,1959) = lu(k,1959) - lu(k,1139) * lu(k,1928)
         lu(k,2048) = lu(k,2048) - lu(k,1130) * lu(k,2047)
         lu(k,2052) = lu(k,2052) - lu(k,1131) * lu(k,2047)
         lu(k,2058) = lu(k,2058) - lu(k,1132) * lu(k,2047)
         lu(k,2062) = lu(k,2062) - lu(k,1133) * lu(k,2047)
         lu(k,2064) = lu(k,2064) - lu(k,1134) * lu(k,2047)
         lu(k,2066) = lu(k,2066) - lu(k,1135) * lu(k,2047)
         lu(k,2067) = lu(k,2067) - lu(k,1136) * lu(k,2047)
         lu(k,2070) = lu(k,2070) - lu(k,1137) * lu(k,2047)
         lu(k,2073) = lu(k,2073) - lu(k,1138) * lu(k,2047)
         lu(k,2076) = lu(k,2076) - lu(k,1139) * lu(k,2047)
         lu(k,2108) = lu(k,2108) - lu(k,1130) * lu(k,2107)
         lu(k,2111) = lu(k,2111) - lu(k,1131) * lu(k,2107)
         lu(k,2118) = lu(k,2118) - lu(k,1132) * lu(k,2107)
         lu(k,2123) = lu(k,2123) - lu(k,1133) * lu(k,2107)
         lu(k,2125) = lu(k,2125) - lu(k,1134) * lu(k,2107)
         lu(k,2127) = lu(k,2127) - lu(k,1135) * lu(k,2107)
         lu(k,2128) = lu(k,2128) - lu(k,1136) * lu(k,2107)
         lu(k,2131) = lu(k,2131) - lu(k,1137) * lu(k,2107)
         lu(k,2134) = lu(k,2134) - lu(k,1138) * lu(k,2107)
         lu(k,2137) = lu(k,2137) - lu(k,1139) * lu(k,2107)
                                                                        
         lu(k,1140) = 1._r8 / lu(k,1140)
         lu(k,1141) = lu(k,1141) * lu(k,1140)
         lu(k,1142) = lu(k,1142) * lu(k,1140)
         lu(k,1143) = lu(k,1143) * lu(k,1140)
         lu(k,1144) = lu(k,1144) * lu(k,1140)
         lu(k,1145) = lu(k,1145) * lu(k,1140)
         lu(k,1154) = lu(k,1154) - lu(k,1141) * lu(k,1152)
         lu(k,1155) = lu(k,1155) - lu(k,1142) * lu(k,1152)
         lu(k,1157) = lu(k,1157) - lu(k,1143) * lu(k,1152)
         lu(k,1158) = lu(k,1158) - lu(k,1144) * lu(k,1152)
         lu(k,1161) = lu(k,1161) - lu(k,1145) * lu(k,1152)
         lu(k,1174) = lu(k,1174) - lu(k,1141) * lu(k,1172)
         lu(k,1175) = lu(k,1175) - lu(k,1142) * lu(k,1172)
         lu(k,1176) = lu(k,1176) - lu(k,1143) * lu(k,1172)
         lu(k,1177) = lu(k,1177) - lu(k,1144) * lu(k,1172)
         lu(k,1180) = lu(k,1180) - lu(k,1145) * lu(k,1172)
         lu(k,1255) = - lu(k,1141) * lu(k,1251)
         lu(k,1260) = lu(k,1260) - lu(k,1142) * lu(k,1251)
         lu(k,1262) = lu(k,1262) - lu(k,1143) * lu(k,1251)
         lu(k,1263) = lu(k,1263) - lu(k,1144) * lu(k,1251)
         lu(k,1266) = lu(k,1266) - lu(k,1145) * lu(k,1251)
         lu(k,1287) = lu(k,1287) - lu(k,1141) * lu(k,1283)
         lu(k,1292) = lu(k,1292) - lu(k,1142) * lu(k,1283)
         lu(k,1294) = lu(k,1294) - lu(k,1143) * lu(k,1283)
         lu(k,1295) = lu(k,1295) - lu(k,1144) * lu(k,1283)
         lu(k,1298) = lu(k,1298) - lu(k,1145) * lu(k,1283)
         lu(k,1310) = lu(k,1310) - lu(k,1141) * lu(k,1308)
         lu(k,1313) = lu(k,1313) - lu(k,1142) * lu(k,1308)
         lu(k,1315) = lu(k,1315) - lu(k,1143) * lu(k,1308)
         lu(k,1316) = lu(k,1316) - lu(k,1144) * lu(k,1308)
         lu(k,1319) = lu(k,1319) - lu(k,1145) * lu(k,1308)
         lu(k,1330) = lu(k,1330) - lu(k,1141) * lu(k,1329)
         lu(k,1334) = lu(k,1334) - lu(k,1142) * lu(k,1329)
         lu(k,1337) = lu(k,1337) - lu(k,1143) * lu(k,1329)
         lu(k,1338) = lu(k,1338) - lu(k,1144) * lu(k,1329)
         lu(k,1341) = lu(k,1341) - lu(k,1145) * lu(k,1329)
         lu(k,1353) = - lu(k,1141) * lu(k,1352)
         lu(k,1355) = lu(k,1355) - lu(k,1142) * lu(k,1352)
         lu(k,1357) = lu(k,1357) - lu(k,1143) * lu(k,1352)
         lu(k,1358) = lu(k,1358) - lu(k,1144) * lu(k,1352)
         lu(k,1361) = lu(k,1361) - lu(k,1145) * lu(k,1352)
         lu(k,1380) = lu(k,1380) - lu(k,1141) * lu(k,1376)
         lu(k,1386) = lu(k,1386) - lu(k,1142) * lu(k,1376)
         lu(k,1389) = lu(k,1389) - lu(k,1143) * lu(k,1376)
         lu(k,1390) = lu(k,1390) - lu(k,1144) * lu(k,1376)
         lu(k,1393) = lu(k,1393) - lu(k,1145) * lu(k,1376)
         lu(k,1676) = lu(k,1676) - lu(k,1141) * lu(k,1671)
         lu(k,1683) = lu(k,1683) - lu(k,1142) * lu(k,1671)
         lu(k,1689) = lu(k,1689) - lu(k,1143) * lu(k,1671)
         lu(k,1691) = lu(k,1691) - lu(k,1144) * lu(k,1671)
         lu(k,1694) = lu(k,1694) - lu(k,1145) * lu(k,1671)
         lu(k,1734) = lu(k,1734) - lu(k,1141) * lu(k,1729)
         lu(k,1741) = lu(k,1741) - lu(k,1142) * lu(k,1729)
         lu(k,1746) = lu(k,1746) - lu(k,1143) * lu(k,1729)
         lu(k,1748) = lu(k,1748) - lu(k,1144) * lu(k,1729)
         lu(k,1751) = lu(k,1751) - lu(k,1145) * lu(k,1729)
         lu(k,1827) = lu(k,1827) - lu(k,1141) * lu(k,1822)
         lu(k,1833) = lu(k,1833) - lu(k,1142) * lu(k,1822)
         lu(k,1838) = lu(k,1838) - lu(k,1143) * lu(k,1822)
         lu(k,1840) = lu(k,1840) - lu(k,1144) * lu(k,1822)
         lu(k,1843) = lu(k,1843) - lu(k,1145) * lu(k,1822)
         lu(k,1933) = lu(k,1933) - lu(k,1141) * lu(k,1929)
         lu(k,1939) = lu(k,1939) - lu(k,1142) * lu(k,1929)
         lu(k,1945) = lu(k,1945) - lu(k,1143) * lu(k,1929)
         lu(k,1947) = lu(k,1947) - lu(k,1144) * lu(k,1929)
         lu(k,1950) = lu(k,1950) - lu(k,1145) * lu(k,1929)
         lu(k,2052) = lu(k,2052) - lu(k,1141) * lu(k,2048)
         lu(k,2058) = lu(k,2058) - lu(k,1142) * lu(k,2048)
         lu(k,2062) = lu(k,2062) - lu(k,1143) * lu(k,2048)
         lu(k,2064) = lu(k,2064) - lu(k,1144) * lu(k,2048)
         lu(k,2067) = lu(k,2067) - lu(k,1145) * lu(k,2048)
         lu(k,2111) = lu(k,2111) - lu(k,1141) * lu(k,2108)
         lu(k,2118) = lu(k,2118) - lu(k,1142) * lu(k,2108)
         lu(k,2123) = lu(k,2123) - lu(k,1143) * lu(k,2108)
         lu(k,2125) = lu(k,2125) - lu(k,1144) * lu(k,2108)
         lu(k,2128) = lu(k,2128) - lu(k,1145) * lu(k,2108)
         lu(k,2182) = lu(k,2182) - lu(k,1141) * lu(k,2181)
         lu(k,2185) = lu(k,2185) - lu(k,1142) * lu(k,2181)
         lu(k,2190) = lu(k,2190) - lu(k,1143) * lu(k,2181)
         lu(k,2192) = lu(k,2192) - lu(k,1144) * lu(k,2181)
         lu(k,2195) = lu(k,2195) - lu(k,1145) * lu(k,2181)
                                                                        
      end do
                                                                        
      end subroutine lu_fac23
                                                                        
      subroutine lu_fac24( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,1153) = 1._r8 / lu(k,1153)
         lu(k,1154) = lu(k,1154) * lu(k,1153)
         lu(k,1155) = lu(k,1155) * lu(k,1153)
         lu(k,1156) = lu(k,1156) * lu(k,1153)
         lu(k,1157) = lu(k,1157) * lu(k,1153)
         lu(k,1158) = lu(k,1158) * lu(k,1153)
         lu(k,1159) = lu(k,1159) * lu(k,1153)
         lu(k,1160) = lu(k,1160) * lu(k,1153)
         lu(k,1161) = lu(k,1161) * lu(k,1153)
         lu(k,1162) = lu(k,1162) * lu(k,1153)
         lu(k,1163) = lu(k,1163) * lu(k,1153)
         lu(k,1164) = lu(k,1164) * lu(k,1153)
         lu(k,1165) = lu(k,1165) * lu(k,1153)
         lu(k,1255) = lu(k,1255) - lu(k,1154) * lu(k,1252)
         lu(k,1260) = lu(k,1260) - lu(k,1155) * lu(k,1252)
         lu(k,1261) = lu(k,1261) - lu(k,1156) * lu(k,1252)
         lu(k,1262) = lu(k,1262) - lu(k,1157) * lu(k,1252)
         lu(k,1263) = lu(k,1263) - lu(k,1158) * lu(k,1252)
         lu(k,1264) = lu(k,1264) - lu(k,1159) * lu(k,1252)
         lu(k,1265) = lu(k,1265) - lu(k,1160) * lu(k,1252)
         lu(k,1266) = lu(k,1266) - lu(k,1161) * lu(k,1252)
         lu(k,1267) = lu(k,1267) - lu(k,1162) * lu(k,1252)
         lu(k,1269) = lu(k,1269) - lu(k,1163) * lu(k,1252)
         lu(k,1270) = lu(k,1270) - lu(k,1164) * lu(k,1252)
         lu(k,1271) = - lu(k,1165) * lu(k,1252)
         lu(k,1287) = lu(k,1287) - lu(k,1154) * lu(k,1284)
         lu(k,1292) = lu(k,1292) - lu(k,1155) * lu(k,1284)
         lu(k,1293) = lu(k,1293) - lu(k,1156) * lu(k,1284)
         lu(k,1294) = lu(k,1294) - lu(k,1157) * lu(k,1284)
         lu(k,1295) = lu(k,1295) - lu(k,1158) * lu(k,1284)
         lu(k,1296) = lu(k,1296) - lu(k,1159) * lu(k,1284)
         lu(k,1297) = lu(k,1297) - lu(k,1160) * lu(k,1284)
         lu(k,1298) = lu(k,1298) - lu(k,1161) * lu(k,1284)
         lu(k,1299) = lu(k,1299) - lu(k,1162) * lu(k,1284)
         lu(k,1301) = lu(k,1301) - lu(k,1163) * lu(k,1284)
         lu(k,1302) = lu(k,1302) - lu(k,1164) * lu(k,1284)
         lu(k,1303) = - lu(k,1165) * lu(k,1284)
         lu(k,1310) = lu(k,1310) - lu(k,1154) * lu(k,1309)
         lu(k,1313) = lu(k,1313) - lu(k,1155) * lu(k,1309)
         lu(k,1314) = - lu(k,1156) * lu(k,1309)
         lu(k,1315) = lu(k,1315) - lu(k,1157) * lu(k,1309)
         lu(k,1316) = lu(k,1316) - lu(k,1158) * lu(k,1309)
         lu(k,1317) = lu(k,1317) - lu(k,1159) * lu(k,1309)
         lu(k,1318) = lu(k,1318) - lu(k,1160) * lu(k,1309)
         lu(k,1319) = lu(k,1319) - lu(k,1161) * lu(k,1309)
         lu(k,1320) = lu(k,1320) - lu(k,1162) * lu(k,1309)
         lu(k,1321) = lu(k,1321) - lu(k,1163) * lu(k,1309)
         lu(k,1322) = lu(k,1322) - lu(k,1164) * lu(k,1309)
         lu(k,1323) = - lu(k,1165) * lu(k,1309)
         lu(k,1676) = lu(k,1676) - lu(k,1154) * lu(k,1672)
         lu(k,1683) = lu(k,1683) - lu(k,1155) * lu(k,1672)
         lu(k,1687) = lu(k,1687) - lu(k,1156) * lu(k,1672)
         lu(k,1689) = lu(k,1689) - lu(k,1157) * lu(k,1672)
         lu(k,1691) = lu(k,1691) - lu(k,1158) * lu(k,1672)
         lu(k,1692) = lu(k,1692) - lu(k,1159) * lu(k,1672)
         lu(k,1693) = lu(k,1693) - lu(k,1160) * lu(k,1672)
         lu(k,1694) = lu(k,1694) - lu(k,1161) * lu(k,1672)
         lu(k,1697) = lu(k,1697) - lu(k,1162) * lu(k,1672)
         lu(k,1700) = lu(k,1700) - lu(k,1163) * lu(k,1672)
         lu(k,1702) = lu(k,1702) - lu(k,1164) * lu(k,1672)
         lu(k,1703) = lu(k,1703) - lu(k,1165) * lu(k,1672)
         lu(k,1734) = lu(k,1734) - lu(k,1154) * lu(k,1730)
         lu(k,1741) = lu(k,1741) - lu(k,1155) * lu(k,1730)
         lu(k,1744) = lu(k,1744) - lu(k,1156) * lu(k,1730)
         lu(k,1746) = lu(k,1746) - lu(k,1157) * lu(k,1730)
         lu(k,1748) = lu(k,1748) - lu(k,1158) * lu(k,1730)
         lu(k,1749) = lu(k,1749) - lu(k,1159) * lu(k,1730)
         lu(k,1750) = lu(k,1750) - lu(k,1160) * lu(k,1730)
         lu(k,1751) = lu(k,1751) - lu(k,1161) * lu(k,1730)
         lu(k,1754) = lu(k,1754) - lu(k,1162) * lu(k,1730)
         lu(k,1757) = lu(k,1757) - lu(k,1163) * lu(k,1730)
         lu(k,1759) = lu(k,1759) - lu(k,1164) * lu(k,1730)
         lu(k,1760) = lu(k,1760) - lu(k,1165) * lu(k,1730)
         lu(k,1827) = lu(k,1827) - lu(k,1154) * lu(k,1823)
         lu(k,1833) = lu(k,1833) - lu(k,1155) * lu(k,1823)
         lu(k,1836) = lu(k,1836) - lu(k,1156) * lu(k,1823)
         lu(k,1838) = lu(k,1838) - lu(k,1157) * lu(k,1823)
         lu(k,1840) = lu(k,1840) - lu(k,1158) * lu(k,1823)
         lu(k,1841) = lu(k,1841) - lu(k,1159) * lu(k,1823)
         lu(k,1842) = lu(k,1842) - lu(k,1160) * lu(k,1823)
         lu(k,1843) = lu(k,1843) - lu(k,1161) * lu(k,1823)
         lu(k,1846) = lu(k,1846) - lu(k,1162) * lu(k,1823)
         lu(k,1849) = lu(k,1849) - lu(k,1163) * lu(k,1823)
         lu(k,1851) = lu(k,1851) - lu(k,1164) * lu(k,1823)
         lu(k,1852) = lu(k,1852) - lu(k,1165) * lu(k,1823)
                                                                        
         lu(k,1173) = 1._r8 / lu(k,1173)
         lu(k,1174) = lu(k,1174) * lu(k,1173)
         lu(k,1175) = lu(k,1175) * lu(k,1173)
         lu(k,1176) = lu(k,1176) * lu(k,1173)
         lu(k,1177) = lu(k,1177) * lu(k,1173)
         lu(k,1178) = lu(k,1178) * lu(k,1173)
         lu(k,1179) = lu(k,1179) * lu(k,1173)
         lu(k,1180) = lu(k,1180) * lu(k,1173)
         lu(k,1181) = lu(k,1181) * lu(k,1173)
         lu(k,1182) = lu(k,1182) * lu(k,1173)
         lu(k,1183) = lu(k,1183) * lu(k,1173)
         lu(k,1195) = lu(k,1195) - lu(k,1174) * lu(k,1192)
         lu(k,1197) = lu(k,1197) - lu(k,1175) * lu(k,1192)
         lu(k,1199) = lu(k,1199) - lu(k,1176) * lu(k,1192)
         lu(k,1200) = lu(k,1200) - lu(k,1177) * lu(k,1192)
         lu(k,1201) = lu(k,1201) - lu(k,1178) * lu(k,1192)
         lu(k,1202) = lu(k,1202) - lu(k,1179) * lu(k,1192)
         lu(k,1203) = lu(k,1203) - lu(k,1180) * lu(k,1192)
         lu(k,1204) = lu(k,1204) - lu(k,1181) * lu(k,1192)
         lu(k,1205) = lu(k,1205) - lu(k,1182) * lu(k,1192)
         lu(k,1206) = lu(k,1206) - lu(k,1183) * lu(k,1192)
         lu(k,1255) = lu(k,1255) - lu(k,1174) * lu(k,1253)
         lu(k,1260) = lu(k,1260) - lu(k,1175) * lu(k,1253)
         lu(k,1262) = lu(k,1262) - lu(k,1176) * lu(k,1253)
         lu(k,1263) = lu(k,1263) - lu(k,1177) * lu(k,1253)
         lu(k,1264) = lu(k,1264) - lu(k,1178) * lu(k,1253)
         lu(k,1265) = lu(k,1265) - lu(k,1179) * lu(k,1253)
         lu(k,1266) = lu(k,1266) - lu(k,1180) * lu(k,1253)
         lu(k,1267) = lu(k,1267) - lu(k,1181) * lu(k,1253)
         lu(k,1269) = lu(k,1269) - lu(k,1182) * lu(k,1253)
         lu(k,1270) = lu(k,1270) - lu(k,1183) * lu(k,1253)
         lu(k,1287) = lu(k,1287) - lu(k,1174) * lu(k,1285)
         lu(k,1292) = lu(k,1292) - lu(k,1175) * lu(k,1285)
         lu(k,1294) = lu(k,1294) - lu(k,1176) * lu(k,1285)
         lu(k,1295) = lu(k,1295) - lu(k,1177) * lu(k,1285)
         lu(k,1296) = lu(k,1296) - lu(k,1178) * lu(k,1285)
         lu(k,1297) = lu(k,1297) - lu(k,1179) * lu(k,1285)
         lu(k,1298) = lu(k,1298) - lu(k,1180) * lu(k,1285)
         lu(k,1299) = lu(k,1299) - lu(k,1181) * lu(k,1285)
         lu(k,1301) = lu(k,1301) - lu(k,1182) * lu(k,1285)
         lu(k,1302) = lu(k,1302) - lu(k,1183) * lu(k,1285)
         lu(k,1380) = lu(k,1380) - lu(k,1174) * lu(k,1377)
         lu(k,1386) = lu(k,1386) - lu(k,1175) * lu(k,1377)
         lu(k,1389) = lu(k,1389) - lu(k,1176) * lu(k,1377)
         lu(k,1390) = lu(k,1390) - lu(k,1177) * lu(k,1377)
         lu(k,1391) = lu(k,1391) - lu(k,1178) * lu(k,1377)
         lu(k,1392) = lu(k,1392) - lu(k,1179) * lu(k,1377)
         lu(k,1393) = lu(k,1393) - lu(k,1180) * lu(k,1377)
         lu(k,1394) = lu(k,1394) - lu(k,1181) * lu(k,1377)
         lu(k,1396) = lu(k,1396) - lu(k,1182) * lu(k,1377)
         lu(k,1397) = lu(k,1397) - lu(k,1183) * lu(k,1377)
         lu(k,1676) = lu(k,1676) - lu(k,1174) * lu(k,1673)
         lu(k,1683) = lu(k,1683) - lu(k,1175) * lu(k,1673)
         lu(k,1689) = lu(k,1689) - lu(k,1176) * lu(k,1673)
         lu(k,1691) = lu(k,1691) - lu(k,1177) * lu(k,1673)
         lu(k,1692) = lu(k,1692) - lu(k,1178) * lu(k,1673)
         lu(k,1693) = lu(k,1693) - lu(k,1179) * lu(k,1673)
         lu(k,1694) = lu(k,1694) - lu(k,1180) * lu(k,1673)
         lu(k,1697) = lu(k,1697) - lu(k,1181) * lu(k,1673)
         lu(k,1700) = lu(k,1700) - lu(k,1182) * lu(k,1673)
         lu(k,1702) = lu(k,1702) - lu(k,1183) * lu(k,1673)
         lu(k,1734) = lu(k,1734) - lu(k,1174) * lu(k,1731)
         lu(k,1741) = lu(k,1741) - lu(k,1175) * lu(k,1731)
         lu(k,1746) = lu(k,1746) - lu(k,1176) * lu(k,1731)
         lu(k,1748) = lu(k,1748) - lu(k,1177) * lu(k,1731)
         lu(k,1749) = lu(k,1749) - lu(k,1178) * lu(k,1731)
         lu(k,1750) = lu(k,1750) - lu(k,1179) * lu(k,1731)
         lu(k,1751) = lu(k,1751) - lu(k,1180) * lu(k,1731)
         lu(k,1754) = lu(k,1754) - lu(k,1181) * lu(k,1731)
         lu(k,1757) = lu(k,1757) - lu(k,1182) * lu(k,1731)
         lu(k,1759) = lu(k,1759) - lu(k,1183) * lu(k,1731)
         lu(k,1827) = lu(k,1827) - lu(k,1174) * lu(k,1824)
         lu(k,1833) = lu(k,1833) - lu(k,1175) * lu(k,1824)
         lu(k,1838) = lu(k,1838) - lu(k,1176) * lu(k,1824)
         lu(k,1840) = lu(k,1840) - lu(k,1177) * lu(k,1824)
         lu(k,1841) = lu(k,1841) - lu(k,1178) * lu(k,1824)
         lu(k,1842) = lu(k,1842) - lu(k,1179) * lu(k,1824)
         lu(k,1843) = lu(k,1843) - lu(k,1180) * lu(k,1824)
         lu(k,1846) = lu(k,1846) - lu(k,1181) * lu(k,1824)
         lu(k,1849) = lu(k,1849) - lu(k,1182) * lu(k,1824)
         lu(k,1851) = lu(k,1851) - lu(k,1183) * lu(k,1824)
         lu(k,1933) = lu(k,1933) - lu(k,1174) * lu(k,1930)
         lu(k,1939) = lu(k,1939) - lu(k,1175) * lu(k,1930)
         lu(k,1945) = lu(k,1945) - lu(k,1176) * lu(k,1930)
         lu(k,1947) = lu(k,1947) - lu(k,1177) * lu(k,1930)
         lu(k,1948) = lu(k,1948) - lu(k,1178) * lu(k,1930)
         lu(k,1949) = lu(k,1949) - lu(k,1179) * lu(k,1930)
         lu(k,1950) = lu(k,1950) - lu(k,1180) * lu(k,1930)
         lu(k,1953) = lu(k,1953) - lu(k,1181) * lu(k,1930)
         lu(k,1956) = lu(k,1956) - lu(k,1182) * lu(k,1930)
         lu(k,1958) = lu(k,1958) - lu(k,1183) * lu(k,1930)
         lu(k,2052) = lu(k,2052) - lu(k,1174) * lu(k,2049)
         lu(k,2058) = lu(k,2058) - lu(k,1175) * lu(k,2049)
         lu(k,2062) = lu(k,2062) - lu(k,1176) * lu(k,2049)
         lu(k,2064) = lu(k,2064) - lu(k,1177) * lu(k,2049)
         lu(k,2065) = lu(k,2065) - lu(k,1178) * lu(k,2049)
         lu(k,2066) = lu(k,2066) - lu(k,1179) * lu(k,2049)
         lu(k,2067) = lu(k,2067) - lu(k,1180) * lu(k,2049)
         lu(k,2070) = lu(k,2070) - lu(k,1181) * lu(k,2049)
         lu(k,2073) = lu(k,2073) - lu(k,1182) * lu(k,2049)
         lu(k,2075) = lu(k,2075) - lu(k,1183) * lu(k,2049)
                                                                        
         lu(k,1193) = 1._r8 / lu(k,1193)
         lu(k,1194) = lu(k,1194) * lu(k,1193)
         lu(k,1195) = lu(k,1195) * lu(k,1193)
         lu(k,1196) = lu(k,1196) * lu(k,1193)
         lu(k,1197) = lu(k,1197) * lu(k,1193)
         lu(k,1198) = lu(k,1198) * lu(k,1193)
         lu(k,1199) = lu(k,1199) * lu(k,1193)
         lu(k,1200) = lu(k,1200) * lu(k,1193)
         lu(k,1201) = lu(k,1201) * lu(k,1193)
         lu(k,1202) = lu(k,1202) * lu(k,1193)
         lu(k,1203) = lu(k,1203) * lu(k,1193)
         lu(k,1204) = lu(k,1204) * lu(k,1193)
         lu(k,1205) = lu(k,1205) * lu(k,1193)
         lu(k,1206) = lu(k,1206) * lu(k,1193)
         lu(k,1379) = lu(k,1379) - lu(k,1194) * lu(k,1378)
         lu(k,1380) = lu(k,1380) - lu(k,1195) * lu(k,1378)
         lu(k,1384) = lu(k,1384) - lu(k,1196) * lu(k,1378)
         lu(k,1386) = lu(k,1386) - lu(k,1197) * lu(k,1378)
         lu(k,1388) = lu(k,1388) - lu(k,1198) * lu(k,1378)
         lu(k,1389) = lu(k,1389) - lu(k,1199) * lu(k,1378)
         lu(k,1390) = lu(k,1390) - lu(k,1200) * lu(k,1378)
         lu(k,1391) = lu(k,1391) - lu(k,1201) * lu(k,1378)
         lu(k,1392) = lu(k,1392) - lu(k,1202) * lu(k,1378)
         lu(k,1393) = lu(k,1393) - lu(k,1203) * lu(k,1378)
         lu(k,1394) = lu(k,1394) - lu(k,1204) * lu(k,1378)
         lu(k,1396) = lu(k,1396) - lu(k,1205) * lu(k,1378)
         lu(k,1397) = lu(k,1397) - lu(k,1206) * lu(k,1378)
         lu(k,1675) = lu(k,1675) - lu(k,1194) * lu(k,1674)
         lu(k,1676) = lu(k,1676) - lu(k,1195) * lu(k,1674)
         lu(k,1681) = lu(k,1681) - lu(k,1196) * lu(k,1674)
         lu(k,1683) = lu(k,1683) - lu(k,1197) * lu(k,1674)
         lu(k,1687) = lu(k,1687) - lu(k,1198) * lu(k,1674)
         lu(k,1689) = lu(k,1689) - lu(k,1199) * lu(k,1674)
         lu(k,1691) = lu(k,1691) - lu(k,1200) * lu(k,1674)
         lu(k,1692) = lu(k,1692) - lu(k,1201) * lu(k,1674)
         lu(k,1693) = lu(k,1693) - lu(k,1202) * lu(k,1674)
         lu(k,1694) = lu(k,1694) - lu(k,1203) * lu(k,1674)
         lu(k,1697) = lu(k,1697) - lu(k,1204) * lu(k,1674)
         lu(k,1700) = lu(k,1700) - lu(k,1205) * lu(k,1674)
         lu(k,1702) = lu(k,1702) - lu(k,1206) * lu(k,1674)
         lu(k,1733) = lu(k,1733) - lu(k,1194) * lu(k,1732)
         lu(k,1734) = lu(k,1734) - lu(k,1195) * lu(k,1732)
         lu(k,1739) = lu(k,1739) - lu(k,1196) * lu(k,1732)
         lu(k,1741) = lu(k,1741) - lu(k,1197) * lu(k,1732)
         lu(k,1744) = lu(k,1744) - lu(k,1198) * lu(k,1732)
         lu(k,1746) = lu(k,1746) - lu(k,1199) * lu(k,1732)
         lu(k,1748) = lu(k,1748) - lu(k,1200) * lu(k,1732)
         lu(k,1749) = lu(k,1749) - lu(k,1201) * lu(k,1732)
         lu(k,1750) = lu(k,1750) - lu(k,1202) * lu(k,1732)
         lu(k,1751) = lu(k,1751) - lu(k,1203) * lu(k,1732)
         lu(k,1754) = lu(k,1754) - lu(k,1204) * lu(k,1732)
         lu(k,1757) = lu(k,1757) - lu(k,1205) * lu(k,1732)
         lu(k,1759) = lu(k,1759) - lu(k,1206) * lu(k,1732)
         lu(k,1826) = lu(k,1826) - lu(k,1194) * lu(k,1825)
         lu(k,1827) = lu(k,1827) - lu(k,1195) * lu(k,1825)
         lu(k,1831) = lu(k,1831) - lu(k,1196) * lu(k,1825)
         lu(k,1833) = lu(k,1833) - lu(k,1197) * lu(k,1825)
         lu(k,1836) = lu(k,1836) - lu(k,1198) * lu(k,1825)
         lu(k,1838) = lu(k,1838) - lu(k,1199) * lu(k,1825)
         lu(k,1840) = lu(k,1840) - lu(k,1200) * lu(k,1825)
         lu(k,1841) = lu(k,1841) - lu(k,1201) * lu(k,1825)
         lu(k,1842) = lu(k,1842) - lu(k,1202) * lu(k,1825)
         lu(k,1843) = lu(k,1843) - lu(k,1203) * lu(k,1825)
         lu(k,1846) = lu(k,1846) - lu(k,1204) * lu(k,1825)
         lu(k,1849) = lu(k,1849) - lu(k,1205) * lu(k,1825)
         lu(k,1851) = lu(k,1851) - lu(k,1206) * lu(k,1825)
         lu(k,1932) = lu(k,1932) - lu(k,1194) * lu(k,1931)
         lu(k,1933) = lu(k,1933) - lu(k,1195) * lu(k,1931)
         lu(k,1937) = lu(k,1937) - lu(k,1196) * lu(k,1931)
         lu(k,1939) = lu(k,1939) - lu(k,1197) * lu(k,1931)
         lu(k,1943) = lu(k,1943) - lu(k,1198) * lu(k,1931)
         lu(k,1945) = lu(k,1945) - lu(k,1199) * lu(k,1931)
         lu(k,1947) = lu(k,1947) - lu(k,1200) * lu(k,1931)
         lu(k,1948) = lu(k,1948) - lu(k,1201) * lu(k,1931)
         lu(k,1949) = lu(k,1949) - lu(k,1202) * lu(k,1931)
         lu(k,1950) = lu(k,1950) - lu(k,1203) * lu(k,1931)
         lu(k,1953) = lu(k,1953) - lu(k,1204) * lu(k,1931)
         lu(k,1956) = lu(k,1956) - lu(k,1205) * lu(k,1931)
         lu(k,1958) = lu(k,1958) - lu(k,1206) * lu(k,1931)
         lu(k,2051) = lu(k,2051) - lu(k,1194) * lu(k,2050)
         lu(k,2052) = lu(k,2052) - lu(k,1195) * lu(k,2050)
         lu(k,2056) = lu(k,2056) - lu(k,1196) * lu(k,2050)
         lu(k,2058) = lu(k,2058) - lu(k,1197) * lu(k,2050)
         lu(k,2060) = lu(k,2060) - lu(k,1198) * lu(k,2050)
         lu(k,2062) = lu(k,2062) - lu(k,1199) * lu(k,2050)
         lu(k,2064) = lu(k,2064) - lu(k,1200) * lu(k,2050)
         lu(k,2065) = lu(k,2065) - lu(k,1201) * lu(k,2050)
         lu(k,2066) = lu(k,2066) - lu(k,1202) * lu(k,2050)
         lu(k,2067) = lu(k,2067) - lu(k,1203) * lu(k,2050)
         lu(k,2070) = lu(k,2070) - lu(k,1204) * lu(k,2050)
         lu(k,2073) = lu(k,2073) - lu(k,1205) * lu(k,2050)
         lu(k,2075) = lu(k,2075) - lu(k,1206) * lu(k,2050)
         lu(k,2110) = lu(k,2110) - lu(k,1194) * lu(k,2109)
         lu(k,2111) = lu(k,2111) - lu(k,1195) * lu(k,2109)
         lu(k,2116) = lu(k,2116) - lu(k,1196) * lu(k,2109)
         lu(k,2118) = lu(k,2118) - lu(k,1197) * lu(k,2109)
         lu(k,2121) = lu(k,2121) - lu(k,1198) * lu(k,2109)
         lu(k,2123) = lu(k,2123) - lu(k,1199) * lu(k,2109)
         lu(k,2125) = lu(k,2125) - lu(k,1200) * lu(k,2109)
         lu(k,2126) = lu(k,2126) - lu(k,1201) * lu(k,2109)
         lu(k,2127) = lu(k,2127) - lu(k,1202) * lu(k,2109)
         lu(k,2128) = lu(k,2128) - lu(k,1203) * lu(k,2109)
         lu(k,2131) = lu(k,2131) - lu(k,1204) * lu(k,2109)
         lu(k,2134) = lu(k,2134) - lu(k,1205) * lu(k,2109)
         lu(k,2136) = lu(k,2136) - lu(k,1206) * lu(k,2109)
                                                                        
         lu(k,1209) = 1._r8 / lu(k,1209)
         lu(k,1210) = lu(k,1210) * lu(k,1209)
         lu(k,1211) = lu(k,1211) * lu(k,1209)
         lu(k,1212) = lu(k,1212) * lu(k,1209)
         lu(k,1213) = lu(k,1213) * lu(k,1209)
         lu(k,1214) = lu(k,1214) * lu(k,1209)
         lu(k,1215) = lu(k,1215) * lu(k,1209)
         lu(k,1216) = lu(k,1216) * lu(k,1209)
         lu(k,1217) = lu(k,1217) * lu(k,1209)
         lu(k,1218) = lu(k,1218) * lu(k,1209)
         lu(k,1219) = lu(k,1219) * lu(k,1209)
         lu(k,1255) = lu(k,1255) - lu(k,1210) * lu(k,1254)
         lu(k,1257) = - lu(k,1211) * lu(k,1254)
         lu(k,1259) = - lu(k,1212) * lu(k,1254)
         lu(k,1260) = lu(k,1260) - lu(k,1213) * lu(k,1254)
         lu(k,1262) = lu(k,1262) - lu(k,1214) * lu(k,1254)
         lu(k,1263) = lu(k,1263) - lu(k,1215) * lu(k,1254)
         lu(k,1266) = lu(k,1266) - lu(k,1216) * lu(k,1254)
         lu(k,1268) = - lu(k,1217) * lu(k,1254)
         lu(k,1270) = lu(k,1270) - lu(k,1218) * lu(k,1254)
         lu(k,1271) = lu(k,1271) - lu(k,1219) * lu(k,1254)
         lu(k,1287) = lu(k,1287) - lu(k,1210) * lu(k,1286)
         lu(k,1289) = - lu(k,1211) * lu(k,1286)
         lu(k,1291) = - lu(k,1212) * lu(k,1286)
         lu(k,1292) = lu(k,1292) - lu(k,1213) * lu(k,1286)
         lu(k,1294) = lu(k,1294) - lu(k,1214) * lu(k,1286)
         lu(k,1295) = lu(k,1295) - lu(k,1215) * lu(k,1286)
         lu(k,1298) = lu(k,1298) - lu(k,1216) * lu(k,1286)
         lu(k,1300) = - lu(k,1217) * lu(k,1286)
         lu(k,1302) = lu(k,1302) - lu(k,1218) * lu(k,1286)
         lu(k,1303) = lu(k,1303) - lu(k,1219) * lu(k,1286)
         lu(k,1380) = lu(k,1380) - lu(k,1210) * lu(k,1379)
         lu(k,1383) = lu(k,1383) - lu(k,1211) * lu(k,1379)
         lu(k,1385) = lu(k,1385) - lu(k,1212) * lu(k,1379)
         lu(k,1386) = lu(k,1386) - lu(k,1213) * lu(k,1379)
         lu(k,1389) = lu(k,1389) - lu(k,1214) * lu(k,1379)
         lu(k,1390) = lu(k,1390) - lu(k,1215) * lu(k,1379)
         lu(k,1393) = lu(k,1393) - lu(k,1216) * lu(k,1379)
         lu(k,1395) = lu(k,1395) - lu(k,1217) * lu(k,1379)
         lu(k,1397) = lu(k,1397) - lu(k,1218) * lu(k,1379)
         lu(k,1398) = lu(k,1398) - lu(k,1219) * lu(k,1379)
         lu(k,1676) = lu(k,1676) - lu(k,1210) * lu(k,1675)
         lu(k,1680) = lu(k,1680) - lu(k,1211) * lu(k,1675)
         lu(k,1682) = lu(k,1682) - lu(k,1212) * lu(k,1675)
         lu(k,1683) = lu(k,1683) - lu(k,1213) * lu(k,1675)
         lu(k,1689) = lu(k,1689) - lu(k,1214) * lu(k,1675)
         lu(k,1691) = lu(k,1691) - lu(k,1215) * lu(k,1675)
         lu(k,1694) = lu(k,1694) - lu(k,1216) * lu(k,1675)
         lu(k,1698) = lu(k,1698) - lu(k,1217) * lu(k,1675)
         lu(k,1702) = lu(k,1702) - lu(k,1218) * lu(k,1675)
         lu(k,1703) = lu(k,1703) - lu(k,1219) * lu(k,1675)
         lu(k,1734) = lu(k,1734) - lu(k,1210) * lu(k,1733)
         lu(k,1738) = lu(k,1738) - lu(k,1211) * lu(k,1733)
         lu(k,1740) = lu(k,1740) - lu(k,1212) * lu(k,1733)
         lu(k,1741) = lu(k,1741) - lu(k,1213) * lu(k,1733)
         lu(k,1746) = lu(k,1746) - lu(k,1214) * lu(k,1733)
         lu(k,1748) = lu(k,1748) - lu(k,1215) * lu(k,1733)
         lu(k,1751) = lu(k,1751) - lu(k,1216) * lu(k,1733)
         lu(k,1755) = lu(k,1755) - lu(k,1217) * lu(k,1733)
         lu(k,1759) = lu(k,1759) - lu(k,1218) * lu(k,1733)
         lu(k,1760) = lu(k,1760) - lu(k,1219) * lu(k,1733)
         lu(k,1827) = lu(k,1827) - lu(k,1210) * lu(k,1826)
         lu(k,1830) = lu(k,1830) - lu(k,1211) * lu(k,1826)
         lu(k,1832) = lu(k,1832) - lu(k,1212) * lu(k,1826)
         lu(k,1833) = lu(k,1833) - lu(k,1213) * lu(k,1826)
         lu(k,1838) = lu(k,1838) - lu(k,1214) * lu(k,1826)
         lu(k,1840) = lu(k,1840) - lu(k,1215) * lu(k,1826)
         lu(k,1843) = lu(k,1843) - lu(k,1216) * lu(k,1826)
         lu(k,1847) = lu(k,1847) - lu(k,1217) * lu(k,1826)
         lu(k,1851) = lu(k,1851) - lu(k,1218) * lu(k,1826)
         lu(k,1852) = lu(k,1852) - lu(k,1219) * lu(k,1826)
         lu(k,1933) = lu(k,1933) - lu(k,1210) * lu(k,1932)
         lu(k,1936) = lu(k,1936) - lu(k,1211) * lu(k,1932)
         lu(k,1938) = lu(k,1938) - lu(k,1212) * lu(k,1932)
         lu(k,1939) = lu(k,1939) - lu(k,1213) * lu(k,1932)
         lu(k,1945) = lu(k,1945) - lu(k,1214) * lu(k,1932)
         lu(k,1947) = lu(k,1947) - lu(k,1215) * lu(k,1932)
         lu(k,1950) = lu(k,1950) - lu(k,1216) * lu(k,1932)
         lu(k,1954) = lu(k,1954) - lu(k,1217) * lu(k,1932)
         lu(k,1958) = lu(k,1958) - lu(k,1218) * lu(k,1932)
         lu(k,1959) = lu(k,1959) - lu(k,1219) * lu(k,1932)
         lu(k,2052) = lu(k,2052) - lu(k,1210) * lu(k,2051)
         lu(k,2055) = lu(k,2055) - lu(k,1211) * lu(k,2051)
         lu(k,2057) = lu(k,2057) - lu(k,1212) * lu(k,2051)
         lu(k,2058) = lu(k,2058) - lu(k,1213) * lu(k,2051)
         lu(k,2062) = lu(k,2062) - lu(k,1214) * lu(k,2051)
         lu(k,2064) = lu(k,2064) - lu(k,1215) * lu(k,2051)
         lu(k,2067) = lu(k,2067) - lu(k,1216) * lu(k,2051)
         lu(k,2071) = - lu(k,1217) * lu(k,2051)
         lu(k,2075) = lu(k,2075) - lu(k,1218) * lu(k,2051)
         lu(k,2076) = lu(k,2076) - lu(k,1219) * lu(k,2051)
         lu(k,2111) = lu(k,2111) - lu(k,1210) * lu(k,2110)
         lu(k,2115) = - lu(k,1211) * lu(k,2110)
         lu(k,2117) = - lu(k,1212) * lu(k,2110)
         lu(k,2118) = lu(k,2118) - lu(k,1213) * lu(k,2110)
         lu(k,2123) = lu(k,2123) - lu(k,1214) * lu(k,2110)
         lu(k,2125) = lu(k,2125) - lu(k,1215) * lu(k,2110)
         lu(k,2128) = lu(k,2128) - lu(k,1216) * lu(k,2110)
         lu(k,2132) = lu(k,2132) - lu(k,1217) * lu(k,2110)
         lu(k,2136) = lu(k,2136) - lu(k,1218) * lu(k,2110)
         lu(k,2137) = lu(k,2137) - lu(k,1219) * lu(k,2110)
                                                                        
      end do
                                                                        
      end subroutine lu_fac24
                                                                        
      subroutine lu_fac25( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,1221) = 1._r8 / lu(k,1221)
         lu(k,1222) = lu(k,1222) * lu(k,1221)
         lu(k,1223) = lu(k,1223) * lu(k,1221)
         lu(k,1224) = lu(k,1224) * lu(k,1221)
         lu(k,1225) = lu(k,1225) * lu(k,1221)
         lu(k,1226) = lu(k,1226) * lu(k,1221)
         lu(k,1227) = lu(k,1227) * lu(k,1221)
         lu(k,1228) = lu(k,1228) * lu(k,1221)
         lu(k,1260) = lu(k,1260) - lu(k,1222) * lu(k,1255)
         lu(k,1261) = lu(k,1261) - lu(k,1223) * lu(k,1255)
         lu(k,1263) = lu(k,1263) - lu(k,1224) * lu(k,1255)
         lu(k,1264) = lu(k,1264) - lu(k,1225) * lu(k,1255)
         lu(k,1266) = lu(k,1266) - lu(k,1226) * lu(k,1255)
         lu(k,1270) = lu(k,1270) - lu(k,1227) * lu(k,1255)
         lu(k,1271) = lu(k,1271) - lu(k,1228) * lu(k,1255)
         lu(k,1292) = lu(k,1292) - lu(k,1222) * lu(k,1287)
         lu(k,1293) = lu(k,1293) - lu(k,1223) * lu(k,1287)
         lu(k,1295) = lu(k,1295) - lu(k,1224) * lu(k,1287)
         lu(k,1296) = lu(k,1296) - lu(k,1225) * lu(k,1287)
         lu(k,1298) = lu(k,1298) - lu(k,1226) * lu(k,1287)
         lu(k,1302) = lu(k,1302) - lu(k,1227) * lu(k,1287)
         lu(k,1303) = lu(k,1303) - lu(k,1228) * lu(k,1287)
         lu(k,1313) = lu(k,1313) - lu(k,1222) * lu(k,1310)
         lu(k,1314) = lu(k,1314) - lu(k,1223) * lu(k,1310)
         lu(k,1316) = lu(k,1316) - lu(k,1224) * lu(k,1310)
         lu(k,1317) = lu(k,1317) - lu(k,1225) * lu(k,1310)
         lu(k,1319) = lu(k,1319) - lu(k,1226) * lu(k,1310)
         lu(k,1322) = lu(k,1322) - lu(k,1227) * lu(k,1310)
         lu(k,1323) = lu(k,1323) - lu(k,1228) * lu(k,1310)
         lu(k,1334) = lu(k,1334) - lu(k,1222) * lu(k,1330)
         lu(k,1336) = lu(k,1336) - lu(k,1223) * lu(k,1330)
         lu(k,1338) = lu(k,1338) - lu(k,1224) * lu(k,1330)
         lu(k,1339) = lu(k,1339) - lu(k,1225) * lu(k,1330)
         lu(k,1341) = lu(k,1341) - lu(k,1226) * lu(k,1330)
         lu(k,1345) = lu(k,1345) - lu(k,1227) * lu(k,1330)
         lu(k,1346) = lu(k,1346) - lu(k,1228) * lu(k,1330)
         lu(k,1355) = lu(k,1355) - lu(k,1222) * lu(k,1353)
         lu(k,1356) = - lu(k,1223) * lu(k,1353)
         lu(k,1358) = lu(k,1358) - lu(k,1224) * lu(k,1353)
         lu(k,1359) = lu(k,1359) - lu(k,1225) * lu(k,1353)
         lu(k,1361) = lu(k,1361) - lu(k,1226) * lu(k,1353)
         lu(k,1365) = lu(k,1365) - lu(k,1227) * lu(k,1353)
         lu(k,1366) = lu(k,1366) - lu(k,1228) * lu(k,1353)
         lu(k,1386) = lu(k,1386) - lu(k,1222) * lu(k,1380)
         lu(k,1388) = lu(k,1388) - lu(k,1223) * lu(k,1380)
         lu(k,1390) = lu(k,1390) - lu(k,1224) * lu(k,1380)
         lu(k,1391) = lu(k,1391) - lu(k,1225) * lu(k,1380)
         lu(k,1393) = lu(k,1393) - lu(k,1226) * lu(k,1380)
         lu(k,1397) = lu(k,1397) - lu(k,1227) * lu(k,1380)
         lu(k,1398) = lu(k,1398) - lu(k,1228) * lu(k,1380)
         lu(k,1683) = lu(k,1683) - lu(k,1222) * lu(k,1676)
         lu(k,1687) = lu(k,1687) - lu(k,1223) * lu(k,1676)
         lu(k,1691) = lu(k,1691) - lu(k,1224) * lu(k,1676)
         lu(k,1692) = lu(k,1692) - lu(k,1225) * lu(k,1676)
         lu(k,1694) = lu(k,1694) - lu(k,1226) * lu(k,1676)
         lu(k,1702) = lu(k,1702) - lu(k,1227) * lu(k,1676)
         lu(k,1703) = lu(k,1703) - lu(k,1228) * lu(k,1676)
         lu(k,1741) = lu(k,1741) - lu(k,1222) * lu(k,1734)
         lu(k,1744) = lu(k,1744) - lu(k,1223) * lu(k,1734)
         lu(k,1748) = lu(k,1748) - lu(k,1224) * lu(k,1734)
         lu(k,1749) = lu(k,1749) - lu(k,1225) * lu(k,1734)
         lu(k,1751) = lu(k,1751) - lu(k,1226) * lu(k,1734)
         lu(k,1759) = lu(k,1759) - lu(k,1227) * lu(k,1734)
         lu(k,1760) = lu(k,1760) - lu(k,1228) * lu(k,1734)
         lu(k,1833) = lu(k,1833) - lu(k,1222) * lu(k,1827)
         lu(k,1836) = lu(k,1836) - lu(k,1223) * lu(k,1827)
         lu(k,1840) = lu(k,1840) - lu(k,1224) * lu(k,1827)
         lu(k,1841) = lu(k,1841) - lu(k,1225) * lu(k,1827)
         lu(k,1843) = lu(k,1843) - lu(k,1226) * lu(k,1827)
         lu(k,1851) = lu(k,1851) - lu(k,1227) * lu(k,1827)
         lu(k,1852) = lu(k,1852) - lu(k,1228) * lu(k,1827)
         lu(k,1939) = lu(k,1939) - lu(k,1222) * lu(k,1933)
         lu(k,1943) = lu(k,1943) - lu(k,1223) * lu(k,1933)
         lu(k,1947) = lu(k,1947) - lu(k,1224) * lu(k,1933)
         lu(k,1948) = lu(k,1948) - lu(k,1225) * lu(k,1933)
         lu(k,1950) = lu(k,1950) - lu(k,1226) * lu(k,1933)
         lu(k,1958) = lu(k,1958) - lu(k,1227) * lu(k,1933)
         lu(k,1959) = lu(k,1959) - lu(k,1228) * lu(k,1933)
         lu(k,2058) = lu(k,2058) - lu(k,1222) * lu(k,2052)
         lu(k,2060) = lu(k,2060) - lu(k,1223) * lu(k,2052)
         lu(k,2064) = lu(k,2064) - lu(k,1224) * lu(k,2052)
         lu(k,2065) = lu(k,2065) - lu(k,1225) * lu(k,2052)
         lu(k,2067) = lu(k,2067) - lu(k,1226) * lu(k,2052)
         lu(k,2075) = lu(k,2075) - lu(k,1227) * lu(k,2052)
         lu(k,2076) = lu(k,2076) - lu(k,1228) * lu(k,2052)
         lu(k,2118) = lu(k,2118) - lu(k,1222) * lu(k,2111)
         lu(k,2121) = lu(k,2121) - lu(k,1223) * lu(k,2111)
         lu(k,2125) = lu(k,2125) - lu(k,1224) * lu(k,2111)
         lu(k,2126) = lu(k,2126) - lu(k,1225) * lu(k,2111)
         lu(k,2128) = lu(k,2128) - lu(k,1226) * lu(k,2111)
         lu(k,2136) = lu(k,2136) - lu(k,1227) * lu(k,2111)
         lu(k,2137) = lu(k,2137) - lu(k,1228) * lu(k,2111)
         lu(k,2185) = lu(k,2185) - lu(k,1222) * lu(k,2182)
         lu(k,2188) = lu(k,2188) - lu(k,1223) * lu(k,2182)
         lu(k,2192) = lu(k,2192) - lu(k,1224) * lu(k,2182)
         lu(k,2193) = lu(k,2193) - lu(k,1225) * lu(k,2182)
         lu(k,2195) = lu(k,2195) - lu(k,1226) * lu(k,2182)
         lu(k,2203) = lu(k,2203) - lu(k,1227) * lu(k,2182)
         lu(k,2204) = lu(k,2204) - lu(k,1228) * lu(k,2182)
                                                                        
         lu(k,1232) = 1._r8 / lu(k,1232)
         lu(k,1233) = lu(k,1233) * lu(k,1232)
         lu(k,1234) = lu(k,1234) * lu(k,1232)
         lu(k,1235) = lu(k,1235) * lu(k,1232)
         lu(k,1236) = lu(k,1236) * lu(k,1232)
         lu(k,1237) = lu(k,1237) * lu(k,1232)
         lu(k,1238) = lu(k,1238) * lu(k,1232)
         lu(k,1239) = lu(k,1239) * lu(k,1232)
         lu(k,1240) = lu(k,1240) * lu(k,1232)
         lu(k,1241) = lu(k,1241) * lu(k,1232)
         lu(k,1242) = lu(k,1242) * lu(k,1232)
         lu(k,1243) = lu(k,1243) * lu(k,1232)
         lu(k,1244) = lu(k,1244) * lu(k,1232)
         lu(k,1685) = lu(k,1685) - lu(k,1233) * lu(k,1677)
         lu(k,1688) = lu(k,1688) - lu(k,1234) * lu(k,1677)
         lu(k,1691) = lu(k,1691) - lu(k,1235) * lu(k,1677)
         lu(k,1693) = lu(k,1693) - lu(k,1236) * lu(k,1677)
         lu(k,1694) = lu(k,1694) - lu(k,1237) * lu(k,1677)
         lu(k,1695) = lu(k,1695) - lu(k,1238) * lu(k,1677)
         lu(k,1696) = lu(k,1696) - lu(k,1239) * lu(k,1677)
         lu(k,1698) = lu(k,1698) - lu(k,1240) * lu(k,1677)
         lu(k,1700) = lu(k,1700) - lu(k,1241) * lu(k,1677)
         lu(k,1701) = lu(k,1701) - lu(k,1242) * lu(k,1677)
         lu(k,1702) = lu(k,1702) - lu(k,1243) * lu(k,1677)
         lu(k,1703) = lu(k,1703) - lu(k,1244) * lu(k,1677)
         lu(k,1742) = - lu(k,1233) * lu(k,1735)
         lu(k,1745) = - lu(k,1234) * lu(k,1735)
         lu(k,1748) = lu(k,1748) - lu(k,1235) * lu(k,1735)
         lu(k,1750) = lu(k,1750) - lu(k,1236) * lu(k,1735)
         lu(k,1751) = lu(k,1751) - lu(k,1237) * lu(k,1735)
         lu(k,1752) = - lu(k,1238) * lu(k,1735)
         lu(k,1753) = - lu(k,1239) * lu(k,1735)
         lu(k,1755) = lu(k,1755) - lu(k,1240) * lu(k,1735)
         lu(k,1757) = lu(k,1757) - lu(k,1241) * lu(k,1735)
         lu(k,1758) = - lu(k,1242) * lu(k,1735)
         lu(k,1759) = lu(k,1759) - lu(k,1243) * lu(k,1735)
         lu(k,1760) = lu(k,1760) - lu(k,1244) * lu(k,1735)
         lu(k,1968) = lu(k,1968) - lu(k,1233) * lu(k,1967)
         lu(k,1970) = - lu(k,1234) * lu(k,1967)
         lu(k,1973) = lu(k,1973) - lu(k,1235) * lu(k,1967)
         lu(k,1975) = lu(k,1975) - lu(k,1236) * lu(k,1967)
         lu(k,1976) = lu(k,1976) - lu(k,1237) * lu(k,1967)
         lu(k,1977) = lu(k,1977) - lu(k,1238) * lu(k,1967)
         lu(k,1978) = lu(k,1978) - lu(k,1239) * lu(k,1967)
         lu(k,1980) = - lu(k,1240) * lu(k,1967)
         lu(k,1982) = lu(k,1982) - lu(k,1241) * lu(k,1967)
         lu(k,1983) = lu(k,1983) - lu(k,1242) * lu(k,1967)
         lu(k,1984) = lu(k,1984) - lu(k,1243) * lu(k,1967)
         lu(k,1985) = lu(k,1985) - lu(k,1244) * lu(k,1967)
         lu(k,2119) = lu(k,2119) - lu(k,1233) * lu(k,2112)
         lu(k,2122) = lu(k,2122) - lu(k,1234) * lu(k,2112)
         lu(k,2125) = lu(k,2125) - lu(k,1235) * lu(k,2112)
         lu(k,2127) = lu(k,2127) - lu(k,1236) * lu(k,2112)
         lu(k,2128) = lu(k,2128) - lu(k,1237) * lu(k,2112)
         lu(k,2129) = lu(k,2129) - lu(k,1238) * lu(k,2112)
         lu(k,2130) = lu(k,2130) - lu(k,1239) * lu(k,2112)
         lu(k,2132) = lu(k,2132) - lu(k,1240) * lu(k,2112)
         lu(k,2134) = lu(k,2134) - lu(k,1241) * lu(k,2112)
         lu(k,2135) = lu(k,2135) - lu(k,1242) * lu(k,2112)
         lu(k,2136) = lu(k,2136) - lu(k,1243) * lu(k,2112)
         lu(k,2137) = lu(k,2137) - lu(k,1244) * lu(k,2112)
         lu(k,2186) = lu(k,2186) - lu(k,1233) * lu(k,2183)
         lu(k,2189) = lu(k,2189) - lu(k,1234) * lu(k,2183)
         lu(k,2192) = lu(k,2192) - lu(k,1235) * lu(k,2183)
         lu(k,2194) = lu(k,2194) - lu(k,1236) * lu(k,2183)
         lu(k,2195) = lu(k,2195) - lu(k,1237) * lu(k,2183)
         lu(k,2196) = lu(k,2196) - lu(k,1238) * lu(k,2183)
         lu(k,2197) = lu(k,2197) - lu(k,1239) * lu(k,2183)
         lu(k,2199) = lu(k,2199) - lu(k,1240) * lu(k,2183)
         lu(k,2201) = lu(k,2201) - lu(k,1241) * lu(k,2183)
         lu(k,2202) = lu(k,2202) - lu(k,1242) * lu(k,2183)
         lu(k,2203) = lu(k,2203) - lu(k,1243) * lu(k,2183)
         lu(k,2204) = lu(k,2204) - lu(k,1244) * lu(k,2183)
         lu(k,2211) = lu(k,2211) - lu(k,1233) * lu(k,2210)
         lu(k,2213) = - lu(k,1234) * lu(k,2210)
         lu(k,2216) = lu(k,2216) - lu(k,1235) * lu(k,2210)
         lu(k,2218) = lu(k,2218) - lu(k,1236) * lu(k,2210)
         lu(k,2219) = lu(k,2219) - lu(k,1237) * lu(k,2210)
         lu(k,2220) = lu(k,2220) - lu(k,1238) * lu(k,2210)
         lu(k,2221) = lu(k,2221) - lu(k,1239) * lu(k,2210)
         lu(k,2223) = - lu(k,1240) * lu(k,2210)
         lu(k,2225) = lu(k,2225) - lu(k,1241) * lu(k,2210)
         lu(k,2226) = lu(k,2226) - lu(k,1242) * lu(k,2210)
         lu(k,2227) = lu(k,2227) - lu(k,1243) * lu(k,2210)
         lu(k,2228) = lu(k,2228) - lu(k,1244) * lu(k,2210)
         lu(k,2241) = lu(k,2241) - lu(k,1233) * lu(k,2239)
         lu(k,2244) = lu(k,2244) - lu(k,1234) * lu(k,2239)
         lu(k,2247) = lu(k,2247) - lu(k,1235) * lu(k,2239)
         lu(k,2249) = lu(k,2249) - lu(k,1236) * lu(k,2239)
         lu(k,2250) = lu(k,2250) - lu(k,1237) * lu(k,2239)
         lu(k,2251) = lu(k,2251) - lu(k,1238) * lu(k,2239)
         lu(k,2252) = lu(k,2252) - lu(k,1239) * lu(k,2239)
         lu(k,2254) = lu(k,2254) - lu(k,1240) * lu(k,2239)
         lu(k,2256) = lu(k,2256) - lu(k,1241) * lu(k,2239)
         lu(k,2257) = lu(k,2257) - lu(k,1242) * lu(k,2239)
         lu(k,2258) = lu(k,2258) - lu(k,1243) * lu(k,2239)
         lu(k,2259) = lu(k,2259) - lu(k,1244) * lu(k,2239)
         lu(k,2267) = - lu(k,1233) * lu(k,2265)
         lu(k,2270) = lu(k,2270) - lu(k,1234) * lu(k,2265)
         lu(k,2273) = lu(k,2273) - lu(k,1235) * lu(k,2265)
         lu(k,2275) = - lu(k,1236) * lu(k,2265)
         lu(k,2276) = lu(k,2276) - lu(k,1237) * lu(k,2265)
         lu(k,2277) = - lu(k,1238) * lu(k,2265)
         lu(k,2278) = - lu(k,1239) * lu(k,2265)
         lu(k,2280) = - lu(k,1240) * lu(k,2265)
         lu(k,2282) = - lu(k,1241) * lu(k,2265)
         lu(k,2283) = - lu(k,1242) * lu(k,2265)
         lu(k,2284) = lu(k,2284) - lu(k,1243) * lu(k,2265)
         lu(k,2285) = lu(k,2285) - lu(k,1244) * lu(k,2265)
                                                                        
         lu(k,1256) = 1._r8 / lu(k,1256)
         lu(k,1257) = lu(k,1257) * lu(k,1256)
         lu(k,1258) = lu(k,1258) * lu(k,1256)
         lu(k,1259) = lu(k,1259) * lu(k,1256)
         lu(k,1260) = lu(k,1260) * lu(k,1256)
         lu(k,1261) = lu(k,1261) * lu(k,1256)
         lu(k,1262) = lu(k,1262) * lu(k,1256)
         lu(k,1263) = lu(k,1263) * lu(k,1256)
         lu(k,1264) = lu(k,1264) * lu(k,1256)
         lu(k,1265) = lu(k,1265) * lu(k,1256)
         lu(k,1266) = lu(k,1266) * lu(k,1256)
         lu(k,1267) = lu(k,1267) * lu(k,1256)
         lu(k,1268) = lu(k,1268) * lu(k,1256)
         lu(k,1269) = lu(k,1269) * lu(k,1256)
         lu(k,1270) = lu(k,1270) * lu(k,1256)
         lu(k,1271) = lu(k,1271) * lu(k,1256)
         lu(k,1383) = lu(k,1383) - lu(k,1257) * lu(k,1381)
         lu(k,1384) = lu(k,1384) - lu(k,1258) * lu(k,1381)
         lu(k,1385) = lu(k,1385) - lu(k,1259) * lu(k,1381)
         lu(k,1386) = lu(k,1386) - lu(k,1260) * lu(k,1381)
         lu(k,1388) = lu(k,1388) - lu(k,1261) * lu(k,1381)
         lu(k,1389) = lu(k,1389) - lu(k,1262) * lu(k,1381)
         lu(k,1390) = lu(k,1390) - lu(k,1263) * lu(k,1381)
         lu(k,1391) = lu(k,1391) - lu(k,1264) * lu(k,1381)
         lu(k,1392) = lu(k,1392) - lu(k,1265) * lu(k,1381)
         lu(k,1393) = lu(k,1393) - lu(k,1266) * lu(k,1381)
         lu(k,1394) = lu(k,1394) - lu(k,1267) * lu(k,1381)
         lu(k,1395) = lu(k,1395) - lu(k,1268) * lu(k,1381)
         lu(k,1396) = lu(k,1396) - lu(k,1269) * lu(k,1381)
         lu(k,1397) = lu(k,1397) - lu(k,1270) * lu(k,1381)
         lu(k,1398) = lu(k,1398) - lu(k,1271) * lu(k,1381)
         lu(k,1680) = lu(k,1680) - lu(k,1257) * lu(k,1678)
         lu(k,1681) = lu(k,1681) - lu(k,1258) * lu(k,1678)
         lu(k,1682) = lu(k,1682) - lu(k,1259) * lu(k,1678)
         lu(k,1683) = lu(k,1683) - lu(k,1260) * lu(k,1678)
         lu(k,1687) = lu(k,1687) - lu(k,1261) * lu(k,1678)
         lu(k,1689) = lu(k,1689) - lu(k,1262) * lu(k,1678)
         lu(k,1691) = lu(k,1691) - lu(k,1263) * lu(k,1678)
         lu(k,1692) = lu(k,1692) - lu(k,1264) * lu(k,1678)
         lu(k,1693) = lu(k,1693) - lu(k,1265) * lu(k,1678)
         lu(k,1694) = lu(k,1694) - lu(k,1266) * lu(k,1678)
         lu(k,1697) = lu(k,1697) - lu(k,1267) * lu(k,1678)
         lu(k,1698) = lu(k,1698) - lu(k,1268) * lu(k,1678)
         lu(k,1700) = lu(k,1700) - lu(k,1269) * lu(k,1678)
         lu(k,1702) = lu(k,1702) - lu(k,1270) * lu(k,1678)
         lu(k,1703) = lu(k,1703) - lu(k,1271) * lu(k,1678)
         lu(k,1738) = lu(k,1738) - lu(k,1257) * lu(k,1736)
         lu(k,1739) = lu(k,1739) - lu(k,1258) * lu(k,1736)
         lu(k,1740) = lu(k,1740) - lu(k,1259) * lu(k,1736)
         lu(k,1741) = lu(k,1741) - lu(k,1260) * lu(k,1736)
         lu(k,1744) = lu(k,1744) - lu(k,1261) * lu(k,1736)
         lu(k,1746) = lu(k,1746) - lu(k,1262) * lu(k,1736)
         lu(k,1748) = lu(k,1748) - lu(k,1263) * lu(k,1736)
         lu(k,1749) = lu(k,1749) - lu(k,1264) * lu(k,1736)
         lu(k,1750) = lu(k,1750) - lu(k,1265) * lu(k,1736)
         lu(k,1751) = lu(k,1751) - lu(k,1266) * lu(k,1736)
         lu(k,1754) = lu(k,1754) - lu(k,1267) * lu(k,1736)
         lu(k,1755) = lu(k,1755) - lu(k,1268) * lu(k,1736)
         lu(k,1757) = lu(k,1757) - lu(k,1269) * lu(k,1736)
         lu(k,1759) = lu(k,1759) - lu(k,1270) * lu(k,1736)
         lu(k,1760) = lu(k,1760) - lu(k,1271) * lu(k,1736)
         lu(k,1830) = lu(k,1830) - lu(k,1257) * lu(k,1828)
         lu(k,1831) = lu(k,1831) - lu(k,1258) * lu(k,1828)
         lu(k,1832) = lu(k,1832) - lu(k,1259) * lu(k,1828)
         lu(k,1833) = lu(k,1833) - lu(k,1260) * lu(k,1828)
         lu(k,1836) = lu(k,1836) - lu(k,1261) * lu(k,1828)
         lu(k,1838) = lu(k,1838) - lu(k,1262) * lu(k,1828)
         lu(k,1840) = lu(k,1840) - lu(k,1263) * lu(k,1828)
         lu(k,1841) = lu(k,1841) - lu(k,1264) * lu(k,1828)
         lu(k,1842) = lu(k,1842) - lu(k,1265) * lu(k,1828)
         lu(k,1843) = lu(k,1843) - lu(k,1266) * lu(k,1828)
         lu(k,1846) = lu(k,1846) - lu(k,1267) * lu(k,1828)
         lu(k,1847) = lu(k,1847) - lu(k,1268) * lu(k,1828)
         lu(k,1849) = lu(k,1849) - lu(k,1269) * lu(k,1828)
         lu(k,1851) = lu(k,1851) - lu(k,1270) * lu(k,1828)
         lu(k,1852) = lu(k,1852) - lu(k,1271) * lu(k,1828)
         lu(k,1936) = lu(k,1936) - lu(k,1257) * lu(k,1934)
         lu(k,1937) = lu(k,1937) - lu(k,1258) * lu(k,1934)
         lu(k,1938) = lu(k,1938) - lu(k,1259) * lu(k,1934)
         lu(k,1939) = lu(k,1939) - lu(k,1260) * lu(k,1934)
         lu(k,1943) = lu(k,1943) - lu(k,1261) * lu(k,1934)
         lu(k,1945) = lu(k,1945) - lu(k,1262) * lu(k,1934)
         lu(k,1947) = lu(k,1947) - lu(k,1263) * lu(k,1934)
         lu(k,1948) = lu(k,1948) - lu(k,1264) * lu(k,1934)
         lu(k,1949) = lu(k,1949) - lu(k,1265) * lu(k,1934)
         lu(k,1950) = lu(k,1950) - lu(k,1266) * lu(k,1934)
         lu(k,1953) = lu(k,1953) - lu(k,1267) * lu(k,1934)
         lu(k,1954) = lu(k,1954) - lu(k,1268) * lu(k,1934)
         lu(k,1956) = lu(k,1956) - lu(k,1269) * lu(k,1934)
         lu(k,1958) = lu(k,1958) - lu(k,1270) * lu(k,1934)
         lu(k,1959) = lu(k,1959) - lu(k,1271) * lu(k,1934)
         lu(k,2055) = lu(k,2055) - lu(k,1257) * lu(k,2053)
         lu(k,2056) = lu(k,2056) - lu(k,1258) * lu(k,2053)
         lu(k,2057) = lu(k,2057) - lu(k,1259) * lu(k,2053)
         lu(k,2058) = lu(k,2058) - lu(k,1260) * lu(k,2053)
         lu(k,2060) = lu(k,2060) - lu(k,1261) * lu(k,2053)
         lu(k,2062) = lu(k,2062) - lu(k,1262) * lu(k,2053)
         lu(k,2064) = lu(k,2064) - lu(k,1263) * lu(k,2053)
         lu(k,2065) = lu(k,2065) - lu(k,1264) * lu(k,2053)
         lu(k,2066) = lu(k,2066) - lu(k,1265) * lu(k,2053)
         lu(k,2067) = lu(k,2067) - lu(k,1266) * lu(k,2053)
         lu(k,2070) = lu(k,2070) - lu(k,1267) * lu(k,2053)
         lu(k,2071) = lu(k,2071) - lu(k,1268) * lu(k,2053)
         lu(k,2073) = lu(k,2073) - lu(k,1269) * lu(k,2053)
         lu(k,2075) = lu(k,2075) - lu(k,1270) * lu(k,2053)
         lu(k,2076) = lu(k,2076) - lu(k,1271) * lu(k,2053)
         lu(k,2115) = lu(k,2115) - lu(k,1257) * lu(k,2113)
         lu(k,2116) = lu(k,2116) - lu(k,1258) * lu(k,2113)
         lu(k,2117) = lu(k,2117) - lu(k,1259) * lu(k,2113)
         lu(k,2118) = lu(k,2118) - lu(k,1260) * lu(k,2113)
         lu(k,2121) = lu(k,2121) - lu(k,1261) * lu(k,2113)
         lu(k,2123) = lu(k,2123) - lu(k,1262) * lu(k,2113)
         lu(k,2125) = lu(k,2125) - lu(k,1263) * lu(k,2113)
         lu(k,2126) = lu(k,2126) - lu(k,1264) * lu(k,2113)
         lu(k,2127) = lu(k,2127) - lu(k,1265) * lu(k,2113)
         lu(k,2128) = lu(k,2128) - lu(k,1266) * lu(k,2113)
         lu(k,2131) = lu(k,2131) - lu(k,1267) * lu(k,2113)
         lu(k,2132) = lu(k,2132) - lu(k,1268) * lu(k,2113)
         lu(k,2134) = lu(k,2134) - lu(k,1269) * lu(k,2113)
         lu(k,2136) = lu(k,2136) - lu(k,1270) * lu(k,2113)
         lu(k,2137) = lu(k,2137) - lu(k,1271) * lu(k,2113)
                                                                        
      end do
                                                                        
      end subroutine lu_fac25
                                                                        
      subroutine lu_fac26( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,1288) = 1._r8 / lu(k,1288)
         lu(k,1289) = lu(k,1289) * lu(k,1288)
         lu(k,1290) = lu(k,1290) * lu(k,1288)
         lu(k,1291) = lu(k,1291) * lu(k,1288)
         lu(k,1292) = lu(k,1292) * lu(k,1288)
         lu(k,1293) = lu(k,1293) * lu(k,1288)
         lu(k,1294) = lu(k,1294) * lu(k,1288)
         lu(k,1295) = lu(k,1295) * lu(k,1288)
         lu(k,1296) = lu(k,1296) * lu(k,1288)
         lu(k,1297) = lu(k,1297) * lu(k,1288)
         lu(k,1298) = lu(k,1298) * lu(k,1288)
         lu(k,1299) = lu(k,1299) * lu(k,1288)
         lu(k,1300) = lu(k,1300) * lu(k,1288)
         lu(k,1301) = lu(k,1301) * lu(k,1288)
         lu(k,1302) = lu(k,1302) * lu(k,1288)
         lu(k,1303) = lu(k,1303) * lu(k,1288)
         lu(k,1383) = lu(k,1383) - lu(k,1289) * lu(k,1382)
         lu(k,1384) = lu(k,1384) - lu(k,1290) * lu(k,1382)
         lu(k,1385) = lu(k,1385) - lu(k,1291) * lu(k,1382)
         lu(k,1386) = lu(k,1386) - lu(k,1292) * lu(k,1382)
         lu(k,1388) = lu(k,1388) - lu(k,1293) * lu(k,1382)
         lu(k,1389) = lu(k,1389) - lu(k,1294) * lu(k,1382)
         lu(k,1390) = lu(k,1390) - lu(k,1295) * lu(k,1382)
         lu(k,1391) = lu(k,1391) - lu(k,1296) * lu(k,1382)
         lu(k,1392) = lu(k,1392) - lu(k,1297) * lu(k,1382)
         lu(k,1393) = lu(k,1393) - lu(k,1298) * lu(k,1382)
         lu(k,1394) = lu(k,1394) - lu(k,1299) * lu(k,1382)
         lu(k,1395) = lu(k,1395) - lu(k,1300) * lu(k,1382)
         lu(k,1396) = lu(k,1396) - lu(k,1301) * lu(k,1382)
         lu(k,1397) = lu(k,1397) - lu(k,1302) * lu(k,1382)
         lu(k,1398) = lu(k,1398) - lu(k,1303) * lu(k,1382)
         lu(k,1680) = lu(k,1680) - lu(k,1289) * lu(k,1679)
         lu(k,1681) = lu(k,1681) - lu(k,1290) * lu(k,1679)
         lu(k,1682) = lu(k,1682) - lu(k,1291) * lu(k,1679)
         lu(k,1683) = lu(k,1683) - lu(k,1292) * lu(k,1679)
         lu(k,1687) = lu(k,1687) - lu(k,1293) * lu(k,1679)
         lu(k,1689) = lu(k,1689) - lu(k,1294) * lu(k,1679)
         lu(k,1691) = lu(k,1691) - lu(k,1295) * lu(k,1679)
         lu(k,1692) = lu(k,1692) - lu(k,1296) * lu(k,1679)
         lu(k,1693) = lu(k,1693) - lu(k,1297) * lu(k,1679)
         lu(k,1694) = lu(k,1694) - lu(k,1298) * lu(k,1679)
         lu(k,1697) = lu(k,1697) - lu(k,1299) * lu(k,1679)
         lu(k,1698) = lu(k,1698) - lu(k,1300) * lu(k,1679)
         lu(k,1700) = lu(k,1700) - lu(k,1301) * lu(k,1679)
         lu(k,1702) = lu(k,1702) - lu(k,1302) * lu(k,1679)
         lu(k,1703) = lu(k,1703) - lu(k,1303) * lu(k,1679)
         lu(k,1738) = lu(k,1738) - lu(k,1289) * lu(k,1737)
         lu(k,1739) = lu(k,1739) - lu(k,1290) * lu(k,1737)
         lu(k,1740) = lu(k,1740) - lu(k,1291) * lu(k,1737)
         lu(k,1741) = lu(k,1741) - lu(k,1292) * lu(k,1737)
         lu(k,1744) = lu(k,1744) - lu(k,1293) * lu(k,1737)
         lu(k,1746) = lu(k,1746) - lu(k,1294) * lu(k,1737)
         lu(k,1748) = lu(k,1748) - lu(k,1295) * lu(k,1737)
         lu(k,1749) = lu(k,1749) - lu(k,1296) * lu(k,1737)
         lu(k,1750) = lu(k,1750) - lu(k,1297) * lu(k,1737)
         lu(k,1751) = lu(k,1751) - lu(k,1298) * lu(k,1737)
         lu(k,1754) = lu(k,1754) - lu(k,1299) * lu(k,1737)
         lu(k,1755) = lu(k,1755) - lu(k,1300) * lu(k,1737)
         lu(k,1757) = lu(k,1757) - lu(k,1301) * lu(k,1737)
         lu(k,1759) = lu(k,1759) - lu(k,1302) * lu(k,1737)
         lu(k,1760) = lu(k,1760) - lu(k,1303) * lu(k,1737)
         lu(k,1830) = lu(k,1830) - lu(k,1289) * lu(k,1829)
         lu(k,1831) = lu(k,1831) - lu(k,1290) * lu(k,1829)
         lu(k,1832) = lu(k,1832) - lu(k,1291) * lu(k,1829)
         lu(k,1833) = lu(k,1833) - lu(k,1292) * lu(k,1829)
         lu(k,1836) = lu(k,1836) - lu(k,1293) * lu(k,1829)
         lu(k,1838) = lu(k,1838) - lu(k,1294) * lu(k,1829)
         lu(k,1840) = lu(k,1840) - lu(k,1295) * lu(k,1829)
         lu(k,1841) = lu(k,1841) - lu(k,1296) * lu(k,1829)
         lu(k,1842) = lu(k,1842) - lu(k,1297) * lu(k,1829)
         lu(k,1843) = lu(k,1843) - lu(k,1298) * lu(k,1829)
         lu(k,1846) = lu(k,1846) - lu(k,1299) * lu(k,1829)
         lu(k,1847) = lu(k,1847) - lu(k,1300) * lu(k,1829)
         lu(k,1849) = lu(k,1849) - lu(k,1301) * lu(k,1829)
         lu(k,1851) = lu(k,1851) - lu(k,1302) * lu(k,1829)
         lu(k,1852) = lu(k,1852) - lu(k,1303) * lu(k,1829)
         lu(k,1936) = lu(k,1936) - lu(k,1289) * lu(k,1935)
         lu(k,1937) = lu(k,1937) - lu(k,1290) * lu(k,1935)
         lu(k,1938) = lu(k,1938) - lu(k,1291) * lu(k,1935)
         lu(k,1939) = lu(k,1939) - lu(k,1292) * lu(k,1935)
         lu(k,1943) = lu(k,1943) - lu(k,1293) * lu(k,1935)
         lu(k,1945) = lu(k,1945) - lu(k,1294) * lu(k,1935)
         lu(k,1947) = lu(k,1947) - lu(k,1295) * lu(k,1935)
         lu(k,1948) = lu(k,1948) - lu(k,1296) * lu(k,1935)
         lu(k,1949) = lu(k,1949) - lu(k,1297) * lu(k,1935)
         lu(k,1950) = lu(k,1950) - lu(k,1298) * lu(k,1935)
         lu(k,1953) = lu(k,1953) - lu(k,1299) * lu(k,1935)
         lu(k,1954) = lu(k,1954) - lu(k,1300) * lu(k,1935)
         lu(k,1956) = lu(k,1956) - lu(k,1301) * lu(k,1935)
         lu(k,1958) = lu(k,1958) - lu(k,1302) * lu(k,1935)
         lu(k,1959) = lu(k,1959) - lu(k,1303) * lu(k,1935)
         lu(k,2055) = lu(k,2055) - lu(k,1289) * lu(k,2054)
         lu(k,2056) = lu(k,2056) - lu(k,1290) * lu(k,2054)
         lu(k,2057) = lu(k,2057) - lu(k,1291) * lu(k,2054)
         lu(k,2058) = lu(k,2058) - lu(k,1292) * lu(k,2054)
         lu(k,2060) = lu(k,2060) - lu(k,1293) * lu(k,2054)
         lu(k,2062) = lu(k,2062) - lu(k,1294) * lu(k,2054)
         lu(k,2064) = lu(k,2064) - lu(k,1295) * lu(k,2054)
         lu(k,2065) = lu(k,2065) - lu(k,1296) * lu(k,2054)
         lu(k,2066) = lu(k,2066) - lu(k,1297) * lu(k,2054)
         lu(k,2067) = lu(k,2067) - lu(k,1298) * lu(k,2054)
         lu(k,2070) = lu(k,2070) - lu(k,1299) * lu(k,2054)
         lu(k,2071) = lu(k,2071) - lu(k,1300) * lu(k,2054)
         lu(k,2073) = lu(k,2073) - lu(k,1301) * lu(k,2054)
         lu(k,2075) = lu(k,2075) - lu(k,1302) * lu(k,2054)
         lu(k,2076) = lu(k,2076) - lu(k,1303) * lu(k,2054)
         lu(k,2115) = lu(k,2115) - lu(k,1289) * lu(k,2114)
         lu(k,2116) = lu(k,2116) - lu(k,1290) * lu(k,2114)
         lu(k,2117) = lu(k,2117) - lu(k,1291) * lu(k,2114)
         lu(k,2118) = lu(k,2118) - lu(k,1292) * lu(k,2114)
         lu(k,2121) = lu(k,2121) - lu(k,1293) * lu(k,2114)
         lu(k,2123) = lu(k,2123) - lu(k,1294) * lu(k,2114)
         lu(k,2125) = lu(k,2125) - lu(k,1295) * lu(k,2114)
         lu(k,2126) = lu(k,2126) - lu(k,1296) * lu(k,2114)
         lu(k,2127) = lu(k,2127) - lu(k,1297) * lu(k,2114)
         lu(k,2128) = lu(k,2128) - lu(k,1298) * lu(k,2114)
         lu(k,2131) = lu(k,2131) - lu(k,1299) * lu(k,2114)
         lu(k,2132) = lu(k,2132) - lu(k,1300) * lu(k,2114)
         lu(k,2134) = lu(k,2134) - lu(k,1301) * lu(k,2114)
         lu(k,2136) = lu(k,2136) - lu(k,1302) * lu(k,2114)
         lu(k,2137) = lu(k,2137) - lu(k,1303) * lu(k,2114)
                                                                        
         lu(k,1311) = 1._r8 / lu(k,1311)
         lu(k,1312) = lu(k,1312) * lu(k,1311)
         lu(k,1313) = lu(k,1313) * lu(k,1311)
         lu(k,1314) = lu(k,1314) * lu(k,1311)
         lu(k,1315) = lu(k,1315) * lu(k,1311)
         lu(k,1316) = lu(k,1316) * lu(k,1311)
         lu(k,1317) = lu(k,1317) * lu(k,1311)
         lu(k,1318) = lu(k,1318) * lu(k,1311)
         lu(k,1319) = lu(k,1319) * lu(k,1311)
         lu(k,1320) = lu(k,1320) * lu(k,1311)
         lu(k,1321) = lu(k,1321) * lu(k,1311)
         lu(k,1322) = lu(k,1322) * lu(k,1311)
         lu(k,1323) = lu(k,1323) * lu(k,1311)
         lu(k,1333) = - lu(k,1312) * lu(k,1331)
         lu(k,1334) = lu(k,1334) - lu(k,1313) * lu(k,1331)
         lu(k,1336) = lu(k,1336) - lu(k,1314) * lu(k,1331)
         lu(k,1337) = lu(k,1337) - lu(k,1315) * lu(k,1331)
         lu(k,1338) = lu(k,1338) - lu(k,1316) * lu(k,1331)
         lu(k,1339) = lu(k,1339) - lu(k,1317) * lu(k,1331)
         lu(k,1340) = lu(k,1340) - lu(k,1318) * lu(k,1331)
         lu(k,1341) = lu(k,1341) - lu(k,1319) * lu(k,1331)
         lu(k,1342) = lu(k,1342) - lu(k,1320) * lu(k,1331)
         lu(k,1344) = lu(k,1344) - lu(k,1321) * lu(k,1331)
         lu(k,1345) = lu(k,1345) - lu(k,1322) * lu(k,1331)
         lu(k,1346) = lu(k,1346) - lu(k,1323) * lu(k,1331)
         lu(k,1385) = lu(k,1385) - lu(k,1312) * lu(k,1383)
         lu(k,1386) = lu(k,1386) - lu(k,1313) * lu(k,1383)
         lu(k,1388) = lu(k,1388) - lu(k,1314) * lu(k,1383)
         lu(k,1389) = lu(k,1389) - lu(k,1315) * lu(k,1383)
         lu(k,1390) = lu(k,1390) - lu(k,1316) * lu(k,1383)
         lu(k,1391) = lu(k,1391) - lu(k,1317) * lu(k,1383)
         lu(k,1392) = lu(k,1392) - lu(k,1318) * lu(k,1383)
         lu(k,1393) = lu(k,1393) - lu(k,1319) * lu(k,1383)
         lu(k,1394) = lu(k,1394) - lu(k,1320) * lu(k,1383)
         lu(k,1396) = lu(k,1396) - lu(k,1321) * lu(k,1383)
         lu(k,1397) = lu(k,1397) - lu(k,1322) * lu(k,1383)
         lu(k,1398) = lu(k,1398) - lu(k,1323) * lu(k,1383)
         lu(k,1682) = lu(k,1682) - lu(k,1312) * lu(k,1680)
         lu(k,1683) = lu(k,1683) - lu(k,1313) * lu(k,1680)
         lu(k,1687) = lu(k,1687) - lu(k,1314) * lu(k,1680)
         lu(k,1689) = lu(k,1689) - lu(k,1315) * lu(k,1680)
         lu(k,1691) = lu(k,1691) - lu(k,1316) * lu(k,1680)
         lu(k,1692) = lu(k,1692) - lu(k,1317) * lu(k,1680)
         lu(k,1693) = lu(k,1693) - lu(k,1318) * lu(k,1680)
         lu(k,1694) = lu(k,1694) - lu(k,1319) * lu(k,1680)
         lu(k,1697) = lu(k,1697) - lu(k,1320) * lu(k,1680)
         lu(k,1700) = lu(k,1700) - lu(k,1321) * lu(k,1680)
         lu(k,1702) = lu(k,1702) - lu(k,1322) * lu(k,1680)
         lu(k,1703) = lu(k,1703) - lu(k,1323) * lu(k,1680)
         lu(k,1740) = lu(k,1740) - lu(k,1312) * lu(k,1738)
         lu(k,1741) = lu(k,1741) - lu(k,1313) * lu(k,1738)
         lu(k,1744) = lu(k,1744) - lu(k,1314) * lu(k,1738)
         lu(k,1746) = lu(k,1746) - lu(k,1315) * lu(k,1738)
         lu(k,1748) = lu(k,1748) - lu(k,1316) * lu(k,1738)
         lu(k,1749) = lu(k,1749) - lu(k,1317) * lu(k,1738)
         lu(k,1750) = lu(k,1750) - lu(k,1318) * lu(k,1738)
         lu(k,1751) = lu(k,1751) - lu(k,1319) * lu(k,1738)
         lu(k,1754) = lu(k,1754) - lu(k,1320) * lu(k,1738)
         lu(k,1757) = lu(k,1757) - lu(k,1321) * lu(k,1738)
         lu(k,1759) = lu(k,1759) - lu(k,1322) * lu(k,1738)
         lu(k,1760) = lu(k,1760) - lu(k,1323) * lu(k,1738)
         lu(k,1832) = lu(k,1832) - lu(k,1312) * lu(k,1830)
         lu(k,1833) = lu(k,1833) - lu(k,1313) * lu(k,1830)
         lu(k,1836) = lu(k,1836) - lu(k,1314) * lu(k,1830)
         lu(k,1838) = lu(k,1838) - lu(k,1315) * lu(k,1830)
         lu(k,1840) = lu(k,1840) - lu(k,1316) * lu(k,1830)
         lu(k,1841) = lu(k,1841) - lu(k,1317) * lu(k,1830)
         lu(k,1842) = lu(k,1842) - lu(k,1318) * lu(k,1830)
         lu(k,1843) = lu(k,1843) - lu(k,1319) * lu(k,1830)
         lu(k,1846) = lu(k,1846) - lu(k,1320) * lu(k,1830)
         lu(k,1849) = lu(k,1849) - lu(k,1321) * lu(k,1830)
         lu(k,1851) = lu(k,1851) - lu(k,1322) * lu(k,1830)
         lu(k,1852) = lu(k,1852) - lu(k,1323) * lu(k,1830)
         lu(k,1938) = lu(k,1938) - lu(k,1312) * lu(k,1936)
         lu(k,1939) = lu(k,1939) - lu(k,1313) * lu(k,1936)
         lu(k,1943) = lu(k,1943) - lu(k,1314) * lu(k,1936)
         lu(k,1945) = lu(k,1945) - lu(k,1315) * lu(k,1936)
         lu(k,1947) = lu(k,1947) - lu(k,1316) * lu(k,1936)
         lu(k,1948) = lu(k,1948) - lu(k,1317) * lu(k,1936)
         lu(k,1949) = lu(k,1949) - lu(k,1318) * lu(k,1936)
         lu(k,1950) = lu(k,1950) - lu(k,1319) * lu(k,1936)
         lu(k,1953) = lu(k,1953) - lu(k,1320) * lu(k,1936)
         lu(k,1956) = lu(k,1956) - lu(k,1321) * lu(k,1936)
         lu(k,1958) = lu(k,1958) - lu(k,1322) * lu(k,1936)
         lu(k,1959) = lu(k,1959) - lu(k,1323) * lu(k,1936)
         lu(k,2057) = lu(k,2057) - lu(k,1312) * lu(k,2055)
         lu(k,2058) = lu(k,2058) - lu(k,1313) * lu(k,2055)
         lu(k,2060) = lu(k,2060) - lu(k,1314) * lu(k,2055)
         lu(k,2062) = lu(k,2062) - lu(k,1315) * lu(k,2055)
         lu(k,2064) = lu(k,2064) - lu(k,1316) * lu(k,2055)
         lu(k,2065) = lu(k,2065) - lu(k,1317) * lu(k,2055)
         lu(k,2066) = lu(k,2066) - lu(k,1318) * lu(k,2055)
         lu(k,2067) = lu(k,2067) - lu(k,1319) * lu(k,2055)
         lu(k,2070) = lu(k,2070) - lu(k,1320) * lu(k,2055)
         lu(k,2073) = lu(k,2073) - lu(k,1321) * lu(k,2055)
         lu(k,2075) = lu(k,2075) - lu(k,1322) * lu(k,2055)
         lu(k,2076) = lu(k,2076) - lu(k,1323) * lu(k,2055)
         lu(k,2117) = lu(k,2117) - lu(k,1312) * lu(k,2115)
         lu(k,2118) = lu(k,2118) - lu(k,1313) * lu(k,2115)
         lu(k,2121) = lu(k,2121) - lu(k,1314) * lu(k,2115)
         lu(k,2123) = lu(k,2123) - lu(k,1315) * lu(k,2115)
         lu(k,2125) = lu(k,2125) - lu(k,1316) * lu(k,2115)
         lu(k,2126) = lu(k,2126) - lu(k,1317) * lu(k,2115)
         lu(k,2127) = lu(k,2127) - lu(k,1318) * lu(k,2115)
         lu(k,2128) = lu(k,2128) - lu(k,1319) * lu(k,2115)
         lu(k,2131) = lu(k,2131) - lu(k,1320) * lu(k,2115)
         lu(k,2134) = lu(k,2134) - lu(k,1321) * lu(k,2115)
         lu(k,2136) = lu(k,2136) - lu(k,1322) * lu(k,2115)
         lu(k,2137) = lu(k,2137) - lu(k,1323) * lu(k,2115)
                                                                        
         lu(k,1332) = 1._r8 / lu(k,1332)
         lu(k,1333) = lu(k,1333) * lu(k,1332)
         lu(k,1334) = lu(k,1334) * lu(k,1332)
         lu(k,1335) = lu(k,1335) * lu(k,1332)
         lu(k,1336) = lu(k,1336) * lu(k,1332)
         lu(k,1337) = lu(k,1337) * lu(k,1332)
         lu(k,1338) = lu(k,1338) * lu(k,1332)
         lu(k,1339) = lu(k,1339) * lu(k,1332)
         lu(k,1340) = lu(k,1340) * lu(k,1332)
         lu(k,1341) = lu(k,1341) * lu(k,1332)
         lu(k,1342) = lu(k,1342) * lu(k,1332)
         lu(k,1343) = lu(k,1343) * lu(k,1332)
         lu(k,1344) = lu(k,1344) * lu(k,1332)
         lu(k,1345) = lu(k,1345) * lu(k,1332)
         lu(k,1346) = lu(k,1346) * lu(k,1332)
         lu(k,1385) = lu(k,1385) - lu(k,1333) * lu(k,1384)
         lu(k,1386) = lu(k,1386) - lu(k,1334) * lu(k,1384)
         lu(k,1387) = - lu(k,1335) * lu(k,1384)
         lu(k,1388) = lu(k,1388) - lu(k,1336) * lu(k,1384)
         lu(k,1389) = lu(k,1389) - lu(k,1337) * lu(k,1384)
         lu(k,1390) = lu(k,1390) - lu(k,1338) * lu(k,1384)
         lu(k,1391) = lu(k,1391) - lu(k,1339) * lu(k,1384)
         lu(k,1392) = lu(k,1392) - lu(k,1340) * lu(k,1384)
         lu(k,1393) = lu(k,1393) - lu(k,1341) * lu(k,1384)
         lu(k,1394) = lu(k,1394) - lu(k,1342) * lu(k,1384)
         lu(k,1395) = lu(k,1395) - lu(k,1343) * lu(k,1384)
         lu(k,1396) = lu(k,1396) - lu(k,1344) * lu(k,1384)
         lu(k,1397) = lu(k,1397) - lu(k,1345) * lu(k,1384)
         lu(k,1398) = lu(k,1398) - lu(k,1346) * lu(k,1384)
         lu(k,1682) = lu(k,1682) - lu(k,1333) * lu(k,1681)
         lu(k,1683) = lu(k,1683) - lu(k,1334) * lu(k,1681)
         lu(k,1686) = lu(k,1686) - lu(k,1335) * lu(k,1681)
         lu(k,1687) = lu(k,1687) - lu(k,1336) * lu(k,1681)
         lu(k,1689) = lu(k,1689) - lu(k,1337) * lu(k,1681)
         lu(k,1691) = lu(k,1691) - lu(k,1338) * lu(k,1681)
         lu(k,1692) = lu(k,1692) - lu(k,1339) * lu(k,1681)
         lu(k,1693) = lu(k,1693) - lu(k,1340) * lu(k,1681)
         lu(k,1694) = lu(k,1694) - lu(k,1341) * lu(k,1681)
         lu(k,1697) = lu(k,1697) - lu(k,1342) * lu(k,1681)
         lu(k,1698) = lu(k,1698) - lu(k,1343) * lu(k,1681)
         lu(k,1700) = lu(k,1700) - lu(k,1344) * lu(k,1681)
         lu(k,1702) = lu(k,1702) - lu(k,1345) * lu(k,1681)
         lu(k,1703) = lu(k,1703) - lu(k,1346) * lu(k,1681)
         lu(k,1740) = lu(k,1740) - lu(k,1333) * lu(k,1739)
         lu(k,1741) = lu(k,1741) - lu(k,1334) * lu(k,1739)
         lu(k,1743) = lu(k,1743) - lu(k,1335) * lu(k,1739)
         lu(k,1744) = lu(k,1744) - lu(k,1336) * lu(k,1739)
         lu(k,1746) = lu(k,1746) - lu(k,1337) * lu(k,1739)
         lu(k,1748) = lu(k,1748) - lu(k,1338) * lu(k,1739)
         lu(k,1749) = lu(k,1749) - lu(k,1339) * lu(k,1739)
         lu(k,1750) = lu(k,1750) - lu(k,1340) * lu(k,1739)
         lu(k,1751) = lu(k,1751) - lu(k,1341) * lu(k,1739)
         lu(k,1754) = lu(k,1754) - lu(k,1342) * lu(k,1739)
         lu(k,1755) = lu(k,1755) - lu(k,1343) * lu(k,1739)
         lu(k,1757) = lu(k,1757) - lu(k,1344) * lu(k,1739)
         lu(k,1759) = lu(k,1759) - lu(k,1345) * lu(k,1739)
         lu(k,1760) = lu(k,1760) - lu(k,1346) * lu(k,1739)
         lu(k,1832) = lu(k,1832) - lu(k,1333) * lu(k,1831)
         lu(k,1833) = lu(k,1833) - lu(k,1334) * lu(k,1831)
         lu(k,1835) = - lu(k,1335) * lu(k,1831)
         lu(k,1836) = lu(k,1836) - lu(k,1336) * lu(k,1831)
         lu(k,1838) = lu(k,1838) - lu(k,1337) * lu(k,1831)
         lu(k,1840) = lu(k,1840) - lu(k,1338) * lu(k,1831)
         lu(k,1841) = lu(k,1841) - lu(k,1339) * lu(k,1831)
         lu(k,1842) = lu(k,1842) - lu(k,1340) * lu(k,1831)
         lu(k,1843) = lu(k,1843) - lu(k,1341) * lu(k,1831)
         lu(k,1846) = lu(k,1846) - lu(k,1342) * lu(k,1831)
         lu(k,1847) = lu(k,1847) - lu(k,1343) * lu(k,1831)
         lu(k,1849) = lu(k,1849) - lu(k,1344) * lu(k,1831)
         lu(k,1851) = lu(k,1851) - lu(k,1345) * lu(k,1831)
         lu(k,1852) = lu(k,1852) - lu(k,1346) * lu(k,1831)
         lu(k,1938) = lu(k,1938) - lu(k,1333) * lu(k,1937)
         lu(k,1939) = lu(k,1939) - lu(k,1334) * lu(k,1937)
         lu(k,1942) = - lu(k,1335) * lu(k,1937)
         lu(k,1943) = lu(k,1943) - lu(k,1336) * lu(k,1937)
         lu(k,1945) = lu(k,1945) - lu(k,1337) * lu(k,1937)
         lu(k,1947) = lu(k,1947) - lu(k,1338) * lu(k,1937)
         lu(k,1948) = lu(k,1948) - lu(k,1339) * lu(k,1937)
         lu(k,1949) = lu(k,1949) - lu(k,1340) * lu(k,1937)
         lu(k,1950) = lu(k,1950) - lu(k,1341) * lu(k,1937)
         lu(k,1953) = lu(k,1953) - lu(k,1342) * lu(k,1937)
         lu(k,1954) = lu(k,1954) - lu(k,1343) * lu(k,1937)
         lu(k,1956) = lu(k,1956) - lu(k,1344) * lu(k,1937)
         lu(k,1958) = lu(k,1958) - lu(k,1345) * lu(k,1937)
         lu(k,1959) = lu(k,1959) - lu(k,1346) * lu(k,1937)
         lu(k,2057) = lu(k,2057) - lu(k,1333) * lu(k,2056)
         lu(k,2058) = lu(k,2058) - lu(k,1334) * lu(k,2056)
         lu(k,2059) = - lu(k,1335) * lu(k,2056)
         lu(k,2060) = lu(k,2060) - lu(k,1336) * lu(k,2056)
         lu(k,2062) = lu(k,2062) - lu(k,1337) * lu(k,2056)
         lu(k,2064) = lu(k,2064) - lu(k,1338) * lu(k,2056)
         lu(k,2065) = lu(k,2065) - lu(k,1339) * lu(k,2056)
         lu(k,2066) = lu(k,2066) - lu(k,1340) * lu(k,2056)
         lu(k,2067) = lu(k,2067) - lu(k,1341) * lu(k,2056)
         lu(k,2070) = lu(k,2070) - lu(k,1342) * lu(k,2056)
         lu(k,2071) = lu(k,2071) - lu(k,1343) * lu(k,2056)
         lu(k,2073) = lu(k,2073) - lu(k,1344) * lu(k,2056)
         lu(k,2075) = lu(k,2075) - lu(k,1345) * lu(k,2056)
         lu(k,2076) = lu(k,2076) - lu(k,1346) * lu(k,2056)
         lu(k,2117) = lu(k,2117) - lu(k,1333) * lu(k,2116)
         lu(k,2118) = lu(k,2118) - lu(k,1334) * lu(k,2116)
         lu(k,2120) = lu(k,2120) - lu(k,1335) * lu(k,2116)
         lu(k,2121) = lu(k,2121) - lu(k,1336) * lu(k,2116)
         lu(k,2123) = lu(k,2123) - lu(k,1337) * lu(k,2116)
         lu(k,2125) = lu(k,2125) - lu(k,1338) * lu(k,2116)
         lu(k,2126) = lu(k,2126) - lu(k,1339) * lu(k,2116)
         lu(k,2127) = lu(k,2127) - lu(k,1340) * lu(k,2116)
         lu(k,2128) = lu(k,2128) - lu(k,1341) * lu(k,2116)
         lu(k,2131) = lu(k,2131) - lu(k,1342) * lu(k,2116)
         lu(k,2132) = lu(k,2132) - lu(k,1343) * lu(k,2116)
         lu(k,2134) = lu(k,2134) - lu(k,1344) * lu(k,2116)
         lu(k,2136) = lu(k,2136) - lu(k,1345) * lu(k,2116)
         lu(k,2137) = lu(k,2137) - lu(k,1346) * lu(k,2116)
                                                                        
         lu(k,1354) = 1._r8 / lu(k,1354)
         lu(k,1355) = lu(k,1355) * lu(k,1354)
         lu(k,1356) = lu(k,1356) * lu(k,1354)
         lu(k,1357) = lu(k,1357) * lu(k,1354)
         lu(k,1358) = lu(k,1358) * lu(k,1354)
         lu(k,1359) = lu(k,1359) * lu(k,1354)
         lu(k,1360) = lu(k,1360) * lu(k,1354)
         lu(k,1361) = lu(k,1361) * lu(k,1354)
         lu(k,1362) = lu(k,1362) * lu(k,1354)
         lu(k,1363) = lu(k,1363) * lu(k,1354)
         lu(k,1364) = lu(k,1364) * lu(k,1354)
         lu(k,1365) = lu(k,1365) * lu(k,1354)
         lu(k,1366) = lu(k,1366) * lu(k,1354)
         lu(k,1386) = lu(k,1386) - lu(k,1355) * lu(k,1385)
         lu(k,1388) = lu(k,1388) - lu(k,1356) * lu(k,1385)
         lu(k,1389) = lu(k,1389) - lu(k,1357) * lu(k,1385)
         lu(k,1390) = lu(k,1390) - lu(k,1358) * lu(k,1385)
         lu(k,1391) = lu(k,1391) - lu(k,1359) * lu(k,1385)
         lu(k,1392) = lu(k,1392) - lu(k,1360) * lu(k,1385)
         lu(k,1393) = lu(k,1393) - lu(k,1361) * lu(k,1385)
         lu(k,1394) = lu(k,1394) - lu(k,1362) * lu(k,1385)
         lu(k,1395) = lu(k,1395) - lu(k,1363) * lu(k,1385)
         lu(k,1396) = lu(k,1396) - lu(k,1364) * lu(k,1385)
         lu(k,1397) = lu(k,1397) - lu(k,1365) * lu(k,1385)
         lu(k,1398) = lu(k,1398) - lu(k,1366) * lu(k,1385)
         lu(k,1683) = lu(k,1683) - lu(k,1355) * lu(k,1682)
         lu(k,1687) = lu(k,1687) - lu(k,1356) * lu(k,1682)
         lu(k,1689) = lu(k,1689) - lu(k,1357) * lu(k,1682)
         lu(k,1691) = lu(k,1691) - lu(k,1358) * lu(k,1682)
         lu(k,1692) = lu(k,1692) - lu(k,1359) * lu(k,1682)
         lu(k,1693) = lu(k,1693) - lu(k,1360) * lu(k,1682)
         lu(k,1694) = lu(k,1694) - lu(k,1361) * lu(k,1682)
         lu(k,1697) = lu(k,1697) - lu(k,1362) * lu(k,1682)
         lu(k,1698) = lu(k,1698) - lu(k,1363) * lu(k,1682)
         lu(k,1700) = lu(k,1700) - lu(k,1364) * lu(k,1682)
         lu(k,1702) = lu(k,1702) - lu(k,1365) * lu(k,1682)
         lu(k,1703) = lu(k,1703) - lu(k,1366) * lu(k,1682)
         lu(k,1741) = lu(k,1741) - lu(k,1355) * lu(k,1740)
         lu(k,1744) = lu(k,1744) - lu(k,1356) * lu(k,1740)
         lu(k,1746) = lu(k,1746) - lu(k,1357) * lu(k,1740)
         lu(k,1748) = lu(k,1748) - lu(k,1358) * lu(k,1740)
         lu(k,1749) = lu(k,1749) - lu(k,1359) * lu(k,1740)
         lu(k,1750) = lu(k,1750) - lu(k,1360) * lu(k,1740)
         lu(k,1751) = lu(k,1751) - lu(k,1361) * lu(k,1740)
         lu(k,1754) = lu(k,1754) - lu(k,1362) * lu(k,1740)
         lu(k,1755) = lu(k,1755) - lu(k,1363) * lu(k,1740)
         lu(k,1757) = lu(k,1757) - lu(k,1364) * lu(k,1740)
         lu(k,1759) = lu(k,1759) - lu(k,1365) * lu(k,1740)
         lu(k,1760) = lu(k,1760) - lu(k,1366) * lu(k,1740)
         lu(k,1833) = lu(k,1833) - lu(k,1355) * lu(k,1832)
         lu(k,1836) = lu(k,1836) - lu(k,1356) * lu(k,1832)
         lu(k,1838) = lu(k,1838) - lu(k,1357) * lu(k,1832)
         lu(k,1840) = lu(k,1840) - lu(k,1358) * lu(k,1832)
         lu(k,1841) = lu(k,1841) - lu(k,1359) * lu(k,1832)
         lu(k,1842) = lu(k,1842) - lu(k,1360) * lu(k,1832)
         lu(k,1843) = lu(k,1843) - lu(k,1361) * lu(k,1832)
         lu(k,1846) = lu(k,1846) - lu(k,1362) * lu(k,1832)
         lu(k,1847) = lu(k,1847) - lu(k,1363) * lu(k,1832)
         lu(k,1849) = lu(k,1849) - lu(k,1364) * lu(k,1832)
         lu(k,1851) = lu(k,1851) - lu(k,1365) * lu(k,1832)
         lu(k,1852) = lu(k,1852) - lu(k,1366) * lu(k,1832)
         lu(k,1939) = lu(k,1939) - lu(k,1355) * lu(k,1938)
         lu(k,1943) = lu(k,1943) - lu(k,1356) * lu(k,1938)
         lu(k,1945) = lu(k,1945) - lu(k,1357) * lu(k,1938)
         lu(k,1947) = lu(k,1947) - lu(k,1358) * lu(k,1938)
         lu(k,1948) = lu(k,1948) - lu(k,1359) * lu(k,1938)
         lu(k,1949) = lu(k,1949) - lu(k,1360) * lu(k,1938)
         lu(k,1950) = lu(k,1950) - lu(k,1361) * lu(k,1938)
         lu(k,1953) = lu(k,1953) - lu(k,1362) * lu(k,1938)
         lu(k,1954) = lu(k,1954) - lu(k,1363) * lu(k,1938)
         lu(k,1956) = lu(k,1956) - lu(k,1364) * lu(k,1938)
         lu(k,1958) = lu(k,1958) - lu(k,1365) * lu(k,1938)
         lu(k,1959) = lu(k,1959) - lu(k,1366) * lu(k,1938)
         lu(k,2058) = lu(k,2058) - lu(k,1355) * lu(k,2057)
         lu(k,2060) = lu(k,2060) - lu(k,1356) * lu(k,2057)
         lu(k,2062) = lu(k,2062) - lu(k,1357) * lu(k,2057)
         lu(k,2064) = lu(k,2064) - lu(k,1358) * lu(k,2057)
         lu(k,2065) = lu(k,2065) - lu(k,1359) * lu(k,2057)
         lu(k,2066) = lu(k,2066) - lu(k,1360) * lu(k,2057)
         lu(k,2067) = lu(k,2067) - lu(k,1361) * lu(k,2057)
         lu(k,2070) = lu(k,2070) - lu(k,1362) * lu(k,2057)
         lu(k,2071) = lu(k,2071) - lu(k,1363) * lu(k,2057)
         lu(k,2073) = lu(k,2073) - lu(k,1364) * lu(k,2057)
         lu(k,2075) = lu(k,2075) - lu(k,1365) * lu(k,2057)
         lu(k,2076) = lu(k,2076) - lu(k,1366) * lu(k,2057)
         lu(k,2118) = lu(k,2118) - lu(k,1355) * lu(k,2117)
         lu(k,2121) = lu(k,2121) - lu(k,1356) * lu(k,2117)
         lu(k,2123) = lu(k,2123) - lu(k,1357) * lu(k,2117)
         lu(k,2125) = lu(k,2125) - lu(k,1358) * lu(k,2117)
         lu(k,2126) = lu(k,2126) - lu(k,1359) * lu(k,2117)
         lu(k,2127) = lu(k,2127) - lu(k,1360) * lu(k,2117)
         lu(k,2128) = lu(k,2128) - lu(k,1361) * lu(k,2117)
         lu(k,2131) = lu(k,2131) - lu(k,1362) * lu(k,2117)
         lu(k,2132) = lu(k,2132) - lu(k,1363) * lu(k,2117)
         lu(k,2134) = lu(k,2134) - lu(k,1364) * lu(k,2117)
         lu(k,2136) = lu(k,2136) - lu(k,1365) * lu(k,2117)
         lu(k,2137) = lu(k,2137) - lu(k,1366) * lu(k,2117)
         lu(k,2185) = lu(k,2185) - lu(k,1355) * lu(k,2184)
         lu(k,2188) = lu(k,2188) - lu(k,1356) * lu(k,2184)
         lu(k,2190) = lu(k,2190) - lu(k,1357) * lu(k,2184)
         lu(k,2192) = lu(k,2192) - lu(k,1358) * lu(k,2184)
         lu(k,2193) = lu(k,2193) - lu(k,1359) * lu(k,2184)
         lu(k,2194) = lu(k,2194) - lu(k,1360) * lu(k,2184)
         lu(k,2195) = lu(k,2195) - lu(k,1361) * lu(k,2184)
         lu(k,2198) = lu(k,2198) - lu(k,1362) * lu(k,2184)
         lu(k,2199) = lu(k,2199) - lu(k,1363) * lu(k,2184)
         lu(k,2201) = lu(k,2201) - lu(k,1364) * lu(k,2184)
         lu(k,2203) = lu(k,2203) - lu(k,1365) * lu(k,2184)
         lu(k,2204) = lu(k,2204) - lu(k,1366) * lu(k,2184)
                                                                        
      end do
                                                                        
      end subroutine lu_fac26
                                                                        
      subroutine lu_fac27( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,1386) = 1._r8 / lu(k,1386)
         lu(k,1387) = lu(k,1387) * lu(k,1386)
         lu(k,1388) = lu(k,1388) * lu(k,1386)
         lu(k,1389) = lu(k,1389) * lu(k,1386)
         lu(k,1390) = lu(k,1390) * lu(k,1386)
         lu(k,1391) = lu(k,1391) * lu(k,1386)
         lu(k,1392) = lu(k,1392) * lu(k,1386)
         lu(k,1393) = lu(k,1393) * lu(k,1386)
         lu(k,1394) = lu(k,1394) * lu(k,1386)
         lu(k,1395) = lu(k,1395) * lu(k,1386)
         lu(k,1396) = lu(k,1396) * lu(k,1386)
         lu(k,1397) = lu(k,1397) * lu(k,1386)
         lu(k,1398) = lu(k,1398) * lu(k,1386)
         lu(k,1686) = lu(k,1686) - lu(k,1387) * lu(k,1683)
         lu(k,1687) = lu(k,1687) - lu(k,1388) * lu(k,1683)
         lu(k,1689) = lu(k,1689) - lu(k,1389) * lu(k,1683)
         lu(k,1691) = lu(k,1691) - lu(k,1390) * lu(k,1683)
         lu(k,1692) = lu(k,1692) - lu(k,1391) * lu(k,1683)
         lu(k,1693) = lu(k,1693) - lu(k,1392) * lu(k,1683)
         lu(k,1694) = lu(k,1694) - lu(k,1393) * lu(k,1683)
         lu(k,1697) = lu(k,1697) - lu(k,1394) * lu(k,1683)
         lu(k,1698) = lu(k,1698) - lu(k,1395) * lu(k,1683)
         lu(k,1700) = lu(k,1700) - lu(k,1396) * lu(k,1683)
         lu(k,1702) = lu(k,1702) - lu(k,1397) * lu(k,1683)
         lu(k,1703) = lu(k,1703) - lu(k,1398) * lu(k,1683)
         lu(k,1743) = lu(k,1743) - lu(k,1387) * lu(k,1741)
         lu(k,1744) = lu(k,1744) - lu(k,1388) * lu(k,1741)
         lu(k,1746) = lu(k,1746) - lu(k,1389) * lu(k,1741)
         lu(k,1748) = lu(k,1748) - lu(k,1390) * lu(k,1741)
         lu(k,1749) = lu(k,1749) - lu(k,1391) * lu(k,1741)
         lu(k,1750) = lu(k,1750) - lu(k,1392) * lu(k,1741)
         lu(k,1751) = lu(k,1751) - lu(k,1393) * lu(k,1741)
         lu(k,1754) = lu(k,1754) - lu(k,1394) * lu(k,1741)
         lu(k,1755) = lu(k,1755) - lu(k,1395) * lu(k,1741)
         lu(k,1757) = lu(k,1757) - lu(k,1396) * lu(k,1741)
         lu(k,1759) = lu(k,1759) - lu(k,1397) * lu(k,1741)
         lu(k,1760) = lu(k,1760) - lu(k,1398) * lu(k,1741)
         lu(k,1835) = lu(k,1835) - lu(k,1387) * lu(k,1833)
         lu(k,1836) = lu(k,1836) - lu(k,1388) * lu(k,1833)
         lu(k,1838) = lu(k,1838) - lu(k,1389) * lu(k,1833)
         lu(k,1840) = lu(k,1840) - lu(k,1390) * lu(k,1833)
         lu(k,1841) = lu(k,1841) - lu(k,1391) * lu(k,1833)
         lu(k,1842) = lu(k,1842) - lu(k,1392) * lu(k,1833)
         lu(k,1843) = lu(k,1843) - lu(k,1393) * lu(k,1833)
         lu(k,1846) = lu(k,1846) - lu(k,1394) * lu(k,1833)
         lu(k,1847) = lu(k,1847) - lu(k,1395) * lu(k,1833)
         lu(k,1849) = lu(k,1849) - lu(k,1396) * lu(k,1833)
         lu(k,1851) = lu(k,1851) - lu(k,1397) * lu(k,1833)
         lu(k,1852) = lu(k,1852) - lu(k,1398) * lu(k,1833)
         lu(k,1942) = lu(k,1942) - lu(k,1387) * lu(k,1939)
         lu(k,1943) = lu(k,1943) - lu(k,1388) * lu(k,1939)
         lu(k,1945) = lu(k,1945) - lu(k,1389) * lu(k,1939)
         lu(k,1947) = lu(k,1947) - lu(k,1390) * lu(k,1939)
         lu(k,1948) = lu(k,1948) - lu(k,1391) * lu(k,1939)
         lu(k,1949) = lu(k,1949) - lu(k,1392) * lu(k,1939)
         lu(k,1950) = lu(k,1950) - lu(k,1393) * lu(k,1939)
         lu(k,1953) = lu(k,1953) - lu(k,1394) * lu(k,1939)
         lu(k,1954) = lu(k,1954) - lu(k,1395) * lu(k,1939)
         lu(k,1956) = lu(k,1956) - lu(k,1396) * lu(k,1939)
         lu(k,1958) = lu(k,1958) - lu(k,1397) * lu(k,1939)
         lu(k,1959) = lu(k,1959) - lu(k,1398) * lu(k,1939)
         lu(k,2007) = lu(k,2007) - lu(k,1387) * lu(k,2004)
         lu(k,2008) = lu(k,2008) - lu(k,1388) * lu(k,2004)
         lu(k,2010) = lu(k,2010) - lu(k,1389) * lu(k,2004)
         lu(k,2012) = lu(k,2012) - lu(k,1390) * lu(k,2004)
         lu(k,2013) = lu(k,2013) - lu(k,1391) * lu(k,2004)
         lu(k,2014) = lu(k,2014) - lu(k,1392) * lu(k,2004)
         lu(k,2015) = lu(k,2015) - lu(k,1393) * lu(k,2004)
         lu(k,2018) = lu(k,2018) - lu(k,1394) * lu(k,2004)
         lu(k,2019) = lu(k,2019) - lu(k,1395) * lu(k,2004)
         lu(k,2021) = lu(k,2021) - lu(k,1396) * lu(k,2004)
         lu(k,2023) = lu(k,2023) - lu(k,1397) * lu(k,2004)
         lu(k,2024) = lu(k,2024) - lu(k,1398) * lu(k,2004)
         lu(k,2059) = lu(k,2059) - lu(k,1387) * lu(k,2058)
         lu(k,2060) = lu(k,2060) - lu(k,1388) * lu(k,2058)
         lu(k,2062) = lu(k,2062) - lu(k,1389) * lu(k,2058)
         lu(k,2064) = lu(k,2064) - lu(k,1390) * lu(k,2058)
         lu(k,2065) = lu(k,2065) - lu(k,1391) * lu(k,2058)
         lu(k,2066) = lu(k,2066) - lu(k,1392) * lu(k,2058)
         lu(k,2067) = lu(k,2067) - lu(k,1393) * lu(k,2058)
         lu(k,2070) = lu(k,2070) - lu(k,1394) * lu(k,2058)
         lu(k,2071) = lu(k,2071) - lu(k,1395) * lu(k,2058)
         lu(k,2073) = lu(k,2073) - lu(k,1396) * lu(k,2058)
         lu(k,2075) = lu(k,2075) - lu(k,1397) * lu(k,2058)
         lu(k,2076) = lu(k,2076) - lu(k,1398) * lu(k,2058)
         lu(k,2120) = lu(k,2120) - lu(k,1387) * lu(k,2118)
         lu(k,2121) = lu(k,2121) - lu(k,1388) * lu(k,2118)
         lu(k,2123) = lu(k,2123) - lu(k,1389) * lu(k,2118)
         lu(k,2125) = lu(k,2125) - lu(k,1390) * lu(k,2118)
         lu(k,2126) = lu(k,2126) - lu(k,1391) * lu(k,2118)
         lu(k,2127) = lu(k,2127) - lu(k,1392) * lu(k,2118)
         lu(k,2128) = lu(k,2128) - lu(k,1393) * lu(k,2118)
         lu(k,2131) = lu(k,2131) - lu(k,1394) * lu(k,2118)
         lu(k,2132) = lu(k,2132) - lu(k,1395) * lu(k,2118)
         lu(k,2134) = lu(k,2134) - lu(k,1396) * lu(k,2118)
         lu(k,2136) = lu(k,2136) - lu(k,1397) * lu(k,2118)
         lu(k,2137) = lu(k,2137) - lu(k,1398) * lu(k,2118)
         lu(k,2187) = - lu(k,1387) * lu(k,2185)
         lu(k,2188) = lu(k,2188) - lu(k,1388) * lu(k,2185)
         lu(k,2190) = lu(k,2190) - lu(k,1389) * lu(k,2185)
         lu(k,2192) = lu(k,2192) - lu(k,1390) * lu(k,2185)
         lu(k,2193) = lu(k,2193) - lu(k,1391) * lu(k,2185)
         lu(k,2194) = lu(k,2194) - lu(k,1392) * lu(k,2185)
         lu(k,2195) = lu(k,2195) - lu(k,1393) * lu(k,2185)
         lu(k,2198) = lu(k,2198) - lu(k,1394) * lu(k,2185)
         lu(k,2199) = lu(k,2199) - lu(k,1395) * lu(k,2185)
         lu(k,2201) = lu(k,2201) - lu(k,1396) * lu(k,2185)
         lu(k,2203) = lu(k,2203) - lu(k,1397) * lu(k,2185)
         lu(k,2204) = lu(k,2204) - lu(k,1398) * lu(k,2185)
                                                                        
         lu(k,1401) = 1._r8 / lu(k,1401)
         lu(k,1402) = lu(k,1402) * lu(k,1401)
         lu(k,1403) = lu(k,1403) * lu(k,1401)
         lu(k,1404) = lu(k,1404) * lu(k,1401)
         lu(k,1405) = lu(k,1405) * lu(k,1401)
         lu(k,1406) = lu(k,1406) * lu(k,1401)
         lu(k,1407) = lu(k,1407) * lu(k,1401)
         lu(k,1408) = lu(k,1408) * lu(k,1401)
         lu(k,1409) = lu(k,1409) * lu(k,1401)
         lu(k,1410) = lu(k,1410) * lu(k,1401)
         lu(k,1411) = lu(k,1411) * lu(k,1401)
         lu(k,1412) = lu(k,1412) * lu(k,1401)
         lu(k,1430) = lu(k,1430) - lu(k,1402) * lu(k,1429)
         lu(k,1431) = lu(k,1431) - lu(k,1403) * lu(k,1429)
         lu(k,1432) = lu(k,1432) - lu(k,1404) * lu(k,1429)
         lu(k,1434) = lu(k,1434) - lu(k,1405) * lu(k,1429)
         lu(k,1435) = lu(k,1435) - lu(k,1406) * lu(k,1429)
         lu(k,1436) = lu(k,1436) - lu(k,1407) * lu(k,1429)
         lu(k,1438) = lu(k,1438) - lu(k,1408) * lu(k,1429)
         lu(k,1439) = lu(k,1439) - lu(k,1409) * lu(k,1429)
         lu(k,1440) = lu(k,1440) - lu(k,1410) * lu(k,1429)
         lu(k,1441) = lu(k,1441) - lu(k,1411) * lu(k,1429)
         lu(k,1442) = lu(k,1442) - lu(k,1412) * lu(k,1429)
         lu(k,1446) = lu(k,1446) - lu(k,1402) * lu(k,1445)
         lu(k,1447) = lu(k,1447) - lu(k,1403) * lu(k,1445)
         lu(k,1448) = lu(k,1448) - lu(k,1404) * lu(k,1445)
         lu(k,1450) = - lu(k,1405) * lu(k,1445)
         lu(k,1451) = lu(k,1451) - lu(k,1406) * lu(k,1445)
         lu(k,1452) = lu(k,1452) - lu(k,1407) * lu(k,1445)
         lu(k,1454) = - lu(k,1408) * lu(k,1445)
         lu(k,1455) = lu(k,1455) - lu(k,1409) * lu(k,1445)
         lu(k,1456) = - lu(k,1410) * lu(k,1445)
         lu(k,1458) = - lu(k,1411) * lu(k,1445)
         lu(k,1459) = lu(k,1459) - lu(k,1412) * lu(k,1445)
         lu(k,1461) = - lu(k,1402) * lu(k,1460)
         lu(k,1462) = - lu(k,1403) * lu(k,1460)
         lu(k,1463) = lu(k,1463) - lu(k,1404) * lu(k,1460)
         lu(k,1465) = - lu(k,1405) * lu(k,1460)
         lu(k,1466) = lu(k,1466) - lu(k,1406) * lu(k,1460)
         lu(k,1467) = - lu(k,1407) * lu(k,1460)
         lu(k,1469) = - lu(k,1408) * lu(k,1460)
         lu(k,1470) = - lu(k,1409) * lu(k,1460)
         lu(k,1472) = - lu(k,1410) * lu(k,1460)
         lu(k,1474) = lu(k,1474) - lu(k,1411) * lu(k,1460)
         lu(k,1475) = lu(k,1475) - lu(k,1412) * lu(k,1460)
         lu(k,1482) = - lu(k,1402) * lu(k,1480)
         lu(k,1483) = lu(k,1483) - lu(k,1403) * lu(k,1480)
         lu(k,1484) = lu(k,1484) - lu(k,1404) * lu(k,1480)
         lu(k,1486) = lu(k,1486) - lu(k,1405) * lu(k,1480)
         lu(k,1487) = lu(k,1487) - lu(k,1406) * lu(k,1480)
         lu(k,1488) = lu(k,1488) - lu(k,1407) * lu(k,1480)
         lu(k,1491) = lu(k,1491) - lu(k,1408) * lu(k,1480)
         lu(k,1492) = - lu(k,1409) * lu(k,1480)
         lu(k,1494) = lu(k,1494) - lu(k,1410) * lu(k,1480)
         lu(k,1497) = lu(k,1497) - lu(k,1411) * lu(k,1480)
         lu(k,1498) = lu(k,1498) - lu(k,1412) * lu(k,1480)
         lu(k,1522) = lu(k,1522) - lu(k,1402) * lu(k,1520)
         lu(k,1523) = lu(k,1523) - lu(k,1403) * lu(k,1520)
         lu(k,1524) = lu(k,1524) - lu(k,1404) * lu(k,1520)
         lu(k,1526) = lu(k,1526) - lu(k,1405) * lu(k,1520)
         lu(k,1527) = lu(k,1527) - lu(k,1406) * lu(k,1520)
         lu(k,1528) = lu(k,1528) - lu(k,1407) * lu(k,1520)
         lu(k,1532) = lu(k,1532) - lu(k,1408) * lu(k,1520)
         lu(k,1533) = lu(k,1533) - lu(k,1409) * lu(k,1520)
         lu(k,1535) = lu(k,1535) - lu(k,1410) * lu(k,1520)
         lu(k,1538) = lu(k,1538) - lu(k,1411) * lu(k,1520)
         lu(k,1539) = lu(k,1539) - lu(k,1412) * lu(k,1520)
         lu(k,1686) = lu(k,1686) - lu(k,1402) * lu(k,1684)
         lu(k,1687) = lu(k,1687) - lu(k,1403) * lu(k,1684)
         lu(k,1688) = lu(k,1688) - lu(k,1404) * lu(k,1684)
         lu(k,1690) = lu(k,1690) - lu(k,1405) * lu(k,1684)
         lu(k,1691) = lu(k,1691) - lu(k,1406) * lu(k,1684)
         lu(k,1692) = lu(k,1692) - lu(k,1407) * lu(k,1684)
         lu(k,1696) = lu(k,1696) - lu(k,1408) * lu(k,1684)
         lu(k,1697) = lu(k,1697) - lu(k,1409) * lu(k,1684)
         lu(k,1699) = lu(k,1699) - lu(k,1410) * lu(k,1684)
         lu(k,1702) = lu(k,1702) - lu(k,1411) * lu(k,1684)
         lu(k,1703) = lu(k,1703) - lu(k,1412) * lu(k,1684)
         lu(k,1942) = lu(k,1942) - lu(k,1402) * lu(k,1940)
         lu(k,1943) = lu(k,1943) - lu(k,1403) * lu(k,1940)
         lu(k,1944) = lu(k,1944) - lu(k,1404) * lu(k,1940)
         lu(k,1946) = lu(k,1946) - lu(k,1405) * lu(k,1940)
         lu(k,1947) = lu(k,1947) - lu(k,1406) * lu(k,1940)
         lu(k,1948) = lu(k,1948) - lu(k,1407) * lu(k,1940)
         lu(k,1952) = lu(k,1952) - lu(k,1408) * lu(k,1940)
         lu(k,1953) = lu(k,1953) - lu(k,1409) * lu(k,1940)
         lu(k,1955) = lu(k,1955) - lu(k,1410) * lu(k,1940)
         lu(k,1958) = lu(k,1958) - lu(k,1411) * lu(k,1940)
         lu(k,1959) = lu(k,1959) - lu(k,1412) * lu(k,1940)
         lu(k,2007) = lu(k,2007) - lu(k,1402) * lu(k,2005)
         lu(k,2008) = lu(k,2008) - lu(k,1403) * lu(k,2005)
         lu(k,2009) = lu(k,2009) - lu(k,1404) * lu(k,2005)
         lu(k,2011) = lu(k,2011) - lu(k,1405) * lu(k,2005)
         lu(k,2012) = lu(k,2012) - lu(k,1406) * lu(k,2005)
         lu(k,2013) = lu(k,2013) - lu(k,1407) * lu(k,2005)
         lu(k,2017) = lu(k,2017) - lu(k,1408) * lu(k,2005)
         lu(k,2018) = lu(k,2018) - lu(k,1409) * lu(k,2005)
         lu(k,2020) = lu(k,2020) - lu(k,1410) * lu(k,2005)
         lu(k,2023) = lu(k,2023) - lu(k,1411) * lu(k,2005)
         lu(k,2024) = lu(k,2024) - lu(k,1412) * lu(k,2005)
         lu(k,2242) = - lu(k,1402) * lu(k,2240)
         lu(k,2243) = lu(k,2243) - lu(k,1403) * lu(k,2240)
         lu(k,2244) = lu(k,2244) - lu(k,1404) * lu(k,2240)
         lu(k,2246) = lu(k,2246) - lu(k,1405) * lu(k,2240)
         lu(k,2247) = lu(k,2247) - lu(k,1406) * lu(k,2240)
         lu(k,2248) = lu(k,2248) - lu(k,1407) * lu(k,2240)
         lu(k,2252) = lu(k,2252) - lu(k,1408) * lu(k,2240)
         lu(k,2253) = - lu(k,1409) * lu(k,2240)
         lu(k,2255) = lu(k,2255) - lu(k,1410) * lu(k,2240)
         lu(k,2258) = lu(k,2258) - lu(k,1411) * lu(k,2240)
         lu(k,2259) = lu(k,2259) - lu(k,1412) * lu(k,2240)
         lu(k,2268) = lu(k,2268) - lu(k,1402) * lu(k,2266)
         lu(k,2269) = lu(k,2269) - lu(k,1403) * lu(k,2266)
         lu(k,2270) = lu(k,2270) - lu(k,1404) * lu(k,2266)
         lu(k,2272) = lu(k,2272) - lu(k,1405) * lu(k,2266)
         lu(k,2273) = lu(k,2273) - lu(k,1406) * lu(k,2266)
         lu(k,2274) = lu(k,2274) - lu(k,1407) * lu(k,2266)
         lu(k,2278) = lu(k,2278) - lu(k,1408) * lu(k,2266)
         lu(k,2279) = lu(k,2279) - lu(k,1409) * lu(k,2266)
         lu(k,2281) = - lu(k,1410) * lu(k,2266)
         lu(k,2284) = lu(k,2284) - lu(k,1411) * lu(k,2266)
         lu(k,2285) = lu(k,2285) - lu(k,1412) * lu(k,2266)
                                                                        
         lu(k,1415) = 1._r8 / lu(k,1415)
         lu(k,1416) = lu(k,1416) * lu(k,1415)
         lu(k,1417) = lu(k,1417) * lu(k,1415)
         lu(k,1418) = lu(k,1418) * lu(k,1415)
         lu(k,1419) = lu(k,1419) * lu(k,1415)
         lu(k,1420) = lu(k,1420) * lu(k,1415)
         lu(k,1421) = lu(k,1421) * lu(k,1415)
         lu(k,1422) = lu(k,1422) * lu(k,1415)
         lu(k,1423) = lu(k,1423) * lu(k,1415)
         lu(k,1424) = lu(k,1424) * lu(k,1415)
         lu(k,1484) = lu(k,1484) - lu(k,1416) * lu(k,1481)
         lu(k,1485) = lu(k,1485) - lu(k,1417) * lu(k,1481)
         lu(k,1486) = lu(k,1486) - lu(k,1418) * lu(k,1481)
         lu(k,1487) = lu(k,1487) - lu(k,1419) * lu(k,1481)
         lu(k,1490) = lu(k,1490) - lu(k,1420) * lu(k,1481)
         lu(k,1493) = - lu(k,1421) * lu(k,1481)
         lu(k,1496) = lu(k,1496) - lu(k,1422) * lu(k,1481)
         lu(k,1497) = lu(k,1497) - lu(k,1423) * lu(k,1481)
         lu(k,1498) = lu(k,1498) - lu(k,1424) * lu(k,1481)
         lu(k,1524) = lu(k,1524) - lu(k,1416) * lu(k,1521)
         lu(k,1525) = lu(k,1525) - lu(k,1417) * lu(k,1521)
         lu(k,1526) = lu(k,1526) - lu(k,1418) * lu(k,1521)
         lu(k,1527) = lu(k,1527) - lu(k,1419) * lu(k,1521)
         lu(k,1530) = lu(k,1530) - lu(k,1420) * lu(k,1521)
         lu(k,1534) = lu(k,1534) - lu(k,1421) * lu(k,1521)
         lu(k,1537) = lu(k,1537) - lu(k,1422) * lu(k,1521)
         lu(k,1538) = lu(k,1538) - lu(k,1423) * lu(k,1521)
         lu(k,1539) = lu(k,1539) - lu(k,1424) * lu(k,1521)
         lu(k,1688) = lu(k,1688) - lu(k,1416) * lu(k,1685)
         lu(k,1689) = lu(k,1689) - lu(k,1417) * lu(k,1685)
         lu(k,1690) = lu(k,1690) - lu(k,1418) * lu(k,1685)
         lu(k,1691) = lu(k,1691) - lu(k,1419) * lu(k,1685)
         lu(k,1694) = lu(k,1694) - lu(k,1420) * lu(k,1685)
         lu(k,1698) = lu(k,1698) - lu(k,1421) * lu(k,1685)
         lu(k,1701) = lu(k,1701) - lu(k,1422) * lu(k,1685)
         lu(k,1702) = lu(k,1702) - lu(k,1423) * lu(k,1685)
         lu(k,1703) = lu(k,1703) - lu(k,1424) * lu(k,1685)
         lu(k,1745) = lu(k,1745) - lu(k,1416) * lu(k,1742)
         lu(k,1746) = lu(k,1746) - lu(k,1417) * lu(k,1742)
         lu(k,1747) = - lu(k,1418) * lu(k,1742)
         lu(k,1748) = lu(k,1748) - lu(k,1419) * lu(k,1742)
         lu(k,1751) = lu(k,1751) - lu(k,1420) * lu(k,1742)
         lu(k,1755) = lu(k,1755) - lu(k,1421) * lu(k,1742)
         lu(k,1758) = lu(k,1758) - lu(k,1422) * lu(k,1742)
         lu(k,1759) = lu(k,1759) - lu(k,1423) * lu(k,1742)
         lu(k,1760) = lu(k,1760) - lu(k,1424) * lu(k,1742)
         lu(k,1837) = lu(k,1837) - lu(k,1416) * lu(k,1834)
         lu(k,1838) = lu(k,1838) - lu(k,1417) * lu(k,1834)
         lu(k,1839) = lu(k,1839) - lu(k,1418) * lu(k,1834)
         lu(k,1840) = lu(k,1840) - lu(k,1419) * lu(k,1834)
         lu(k,1843) = lu(k,1843) - lu(k,1420) * lu(k,1834)
         lu(k,1847) = lu(k,1847) - lu(k,1421) * lu(k,1834)
         lu(k,1850) = lu(k,1850) - lu(k,1422) * lu(k,1834)
         lu(k,1851) = lu(k,1851) - lu(k,1423) * lu(k,1834)
         lu(k,1852) = lu(k,1852) - lu(k,1424) * lu(k,1834)
         lu(k,1944) = lu(k,1944) - lu(k,1416) * lu(k,1941)
         lu(k,1945) = lu(k,1945) - lu(k,1417) * lu(k,1941)
         lu(k,1946) = lu(k,1946) - lu(k,1418) * lu(k,1941)
         lu(k,1947) = lu(k,1947) - lu(k,1419) * lu(k,1941)
         lu(k,1950) = lu(k,1950) - lu(k,1420) * lu(k,1941)
         lu(k,1954) = lu(k,1954) - lu(k,1421) * lu(k,1941)
         lu(k,1957) = lu(k,1957) - lu(k,1422) * lu(k,1941)
         lu(k,1958) = lu(k,1958) - lu(k,1423) * lu(k,1941)
         lu(k,1959) = lu(k,1959) - lu(k,1424) * lu(k,1941)
         lu(k,1970) = lu(k,1970) - lu(k,1416) * lu(k,1968)
         lu(k,1971) = lu(k,1971) - lu(k,1417) * lu(k,1968)
         lu(k,1972) = - lu(k,1418) * lu(k,1968)
         lu(k,1973) = lu(k,1973) - lu(k,1419) * lu(k,1968)
         lu(k,1976) = lu(k,1976) - lu(k,1420) * lu(k,1968)
         lu(k,1980) = lu(k,1980) - lu(k,1421) * lu(k,1968)
         lu(k,1983) = lu(k,1983) - lu(k,1422) * lu(k,1968)
         lu(k,1984) = lu(k,1984) - lu(k,1423) * lu(k,1968)
         lu(k,1985) = lu(k,1985) - lu(k,1424) * lu(k,1968)
         lu(k,2009) = lu(k,2009) - lu(k,1416) * lu(k,2006)
         lu(k,2010) = lu(k,2010) - lu(k,1417) * lu(k,2006)
         lu(k,2011) = lu(k,2011) - lu(k,1418) * lu(k,2006)
         lu(k,2012) = lu(k,2012) - lu(k,1419) * lu(k,2006)
         lu(k,2015) = lu(k,2015) - lu(k,1420) * lu(k,2006)
         lu(k,2019) = lu(k,2019) - lu(k,1421) * lu(k,2006)
         lu(k,2022) = - lu(k,1422) * lu(k,2006)
         lu(k,2023) = lu(k,2023) - lu(k,1423) * lu(k,2006)
         lu(k,2024) = lu(k,2024) - lu(k,1424) * lu(k,2006)
         lu(k,2122) = lu(k,2122) - lu(k,1416) * lu(k,2119)
         lu(k,2123) = lu(k,2123) - lu(k,1417) * lu(k,2119)
         lu(k,2124) = lu(k,2124) - lu(k,1418) * lu(k,2119)
         lu(k,2125) = lu(k,2125) - lu(k,1419) * lu(k,2119)
         lu(k,2128) = lu(k,2128) - lu(k,1420) * lu(k,2119)
         lu(k,2132) = lu(k,2132) - lu(k,1421) * lu(k,2119)
         lu(k,2135) = lu(k,2135) - lu(k,1422) * lu(k,2119)
         lu(k,2136) = lu(k,2136) - lu(k,1423) * lu(k,2119)
         lu(k,2137) = lu(k,2137) - lu(k,1424) * lu(k,2119)
         lu(k,2145) = lu(k,2145) - lu(k,1416) * lu(k,2143)
         lu(k,2146) = - lu(k,1417) * lu(k,2143)
         lu(k,2147) = lu(k,2147) - lu(k,1418) * lu(k,2143)
         lu(k,2148) = lu(k,2148) - lu(k,1419) * lu(k,2143)
         lu(k,2151) = - lu(k,1420) * lu(k,2143)
         lu(k,2155) = - lu(k,1421) * lu(k,2143)
         lu(k,2158) = lu(k,2158) - lu(k,1422) * lu(k,2143)
         lu(k,2159) = lu(k,2159) - lu(k,1423) * lu(k,2143)
         lu(k,2160) = lu(k,2160) - lu(k,1424) * lu(k,2143)
         lu(k,2189) = lu(k,2189) - lu(k,1416) * lu(k,2186)
         lu(k,2190) = lu(k,2190) - lu(k,1417) * lu(k,2186)
         lu(k,2191) = lu(k,2191) - lu(k,1418) * lu(k,2186)
         lu(k,2192) = lu(k,2192) - lu(k,1419) * lu(k,2186)
         lu(k,2195) = lu(k,2195) - lu(k,1420) * lu(k,2186)
         lu(k,2199) = lu(k,2199) - lu(k,1421) * lu(k,2186)
         lu(k,2202) = lu(k,2202) - lu(k,1422) * lu(k,2186)
         lu(k,2203) = lu(k,2203) - lu(k,1423) * lu(k,2186)
         lu(k,2204) = lu(k,2204) - lu(k,1424) * lu(k,2186)
         lu(k,2213) = lu(k,2213) - lu(k,1416) * lu(k,2211)
         lu(k,2214) = - lu(k,1417) * lu(k,2211)
         lu(k,2215) = - lu(k,1418) * lu(k,2211)
         lu(k,2216) = lu(k,2216) - lu(k,1419) * lu(k,2211)
         lu(k,2219) = lu(k,2219) - lu(k,1420) * lu(k,2211)
         lu(k,2223) = lu(k,2223) - lu(k,1421) * lu(k,2211)
         lu(k,2226) = lu(k,2226) - lu(k,1422) * lu(k,2211)
         lu(k,2227) = lu(k,2227) - lu(k,1423) * lu(k,2211)
         lu(k,2228) = lu(k,2228) - lu(k,1424) * lu(k,2211)
         lu(k,2244) = lu(k,2244) - lu(k,1416) * lu(k,2241)
         lu(k,2245) = lu(k,2245) - lu(k,1417) * lu(k,2241)
         lu(k,2246) = lu(k,2246) - lu(k,1418) * lu(k,2241)
         lu(k,2247) = lu(k,2247) - lu(k,1419) * lu(k,2241)
         lu(k,2250) = lu(k,2250) - lu(k,1420) * lu(k,2241)
         lu(k,2254) = lu(k,2254) - lu(k,1421) * lu(k,2241)
         lu(k,2257) = lu(k,2257) - lu(k,1422) * lu(k,2241)
         lu(k,2258) = lu(k,2258) - lu(k,1423) * lu(k,2241)
         lu(k,2259) = lu(k,2259) - lu(k,1424) * lu(k,2241)
         lu(k,2270) = lu(k,2270) - lu(k,1416) * lu(k,2267)
         lu(k,2271) = - lu(k,1417) * lu(k,2267)
         lu(k,2272) = lu(k,2272) - lu(k,1418) * lu(k,2267)
         lu(k,2273) = lu(k,2273) - lu(k,1419) * lu(k,2267)
         lu(k,2276) = lu(k,2276) - lu(k,1420) * lu(k,2267)
         lu(k,2280) = lu(k,2280) - lu(k,1421) * lu(k,2267)
         lu(k,2283) = lu(k,2283) - lu(k,1422) * lu(k,2267)
         lu(k,2284) = lu(k,2284) - lu(k,1423) * lu(k,2267)
         lu(k,2285) = lu(k,2285) - lu(k,1424) * lu(k,2267)
                                                                        
      end do
                                                                        
      end subroutine lu_fac27
                                                                        
      subroutine lu_fac28( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,1430) = 1._r8 / lu(k,1430)
         lu(k,1431) = lu(k,1431) * lu(k,1430)
         lu(k,1432) = lu(k,1432) * lu(k,1430)
         lu(k,1433) = lu(k,1433) * lu(k,1430)
         lu(k,1434) = lu(k,1434) * lu(k,1430)
         lu(k,1435) = lu(k,1435) * lu(k,1430)
         lu(k,1436) = lu(k,1436) * lu(k,1430)
         lu(k,1437) = lu(k,1437) * lu(k,1430)
         lu(k,1438) = lu(k,1438) * lu(k,1430)
         lu(k,1439) = lu(k,1439) * lu(k,1430)
         lu(k,1440) = lu(k,1440) * lu(k,1430)
         lu(k,1441) = lu(k,1441) * lu(k,1430)
         lu(k,1442) = lu(k,1442) * lu(k,1430)
         lu(k,1447) = lu(k,1447) - lu(k,1431) * lu(k,1446)
         lu(k,1448) = lu(k,1448) - lu(k,1432) * lu(k,1446)
         lu(k,1449) = - lu(k,1433) * lu(k,1446)
         lu(k,1450) = lu(k,1450) - lu(k,1434) * lu(k,1446)
         lu(k,1451) = lu(k,1451) - lu(k,1435) * lu(k,1446)
         lu(k,1452) = lu(k,1452) - lu(k,1436) * lu(k,1446)
         lu(k,1453) = - lu(k,1437) * lu(k,1446)
         lu(k,1454) = lu(k,1454) - lu(k,1438) * lu(k,1446)
         lu(k,1455) = lu(k,1455) - lu(k,1439) * lu(k,1446)
         lu(k,1456) = lu(k,1456) - lu(k,1440) * lu(k,1446)
         lu(k,1458) = lu(k,1458) - lu(k,1441) * lu(k,1446)
         lu(k,1459) = lu(k,1459) - lu(k,1442) * lu(k,1446)
         lu(k,1462) = lu(k,1462) - lu(k,1431) * lu(k,1461)
         lu(k,1463) = lu(k,1463) - lu(k,1432) * lu(k,1461)
         lu(k,1464) = - lu(k,1433) * lu(k,1461)
         lu(k,1465) = lu(k,1465) - lu(k,1434) * lu(k,1461)
         lu(k,1466) = lu(k,1466) - lu(k,1435) * lu(k,1461)
         lu(k,1467) = lu(k,1467) - lu(k,1436) * lu(k,1461)
         lu(k,1468) = lu(k,1468) - lu(k,1437) * lu(k,1461)
         lu(k,1469) = lu(k,1469) - lu(k,1438) * lu(k,1461)
         lu(k,1470) = lu(k,1470) - lu(k,1439) * lu(k,1461)
         lu(k,1472) = lu(k,1472) - lu(k,1440) * lu(k,1461)
         lu(k,1474) = lu(k,1474) - lu(k,1441) * lu(k,1461)
         lu(k,1475) = lu(k,1475) - lu(k,1442) * lu(k,1461)
         lu(k,1483) = lu(k,1483) - lu(k,1431) * lu(k,1482)
         lu(k,1484) = lu(k,1484) - lu(k,1432) * lu(k,1482)
         lu(k,1485) = lu(k,1485) - lu(k,1433) * lu(k,1482)
         lu(k,1486) = lu(k,1486) - lu(k,1434) * lu(k,1482)
         lu(k,1487) = lu(k,1487) - lu(k,1435) * lu(k,1482)
         lu(k,1488) = lu(k,1488) - lu(k,1436) * lu(k,1482)
         lu(k,1490) = lu(k,1490) - lu(k,1437) * lu(k,1482)
         lu(k,1491) = lu(k,1491) - lu(k,1438) * lu(k,1482)
         lu(k,1492) = lu(k,1492) - lu(k,1439) * lu(k,1482)
         lu(k,1494) = lu(k,1494) - lu(k,1440) * lu(k,1482)
         lu(k,1497) = lu(k,1497) - lu(k,1441) * lu(k,1482)
         lu(k,1498) = lu(k,1498) - lu(k,1442) * lu(k,1482)
         lu(k,1523) = lu(k,1523) - lu(k,1431) * lu(k,1522)
         lu(k,1524) = lu(k,1524) - lu(k,1432) * lu(k,1522)
         lu(k,1525) = lu(k,1525) - lu(k,1433) * lu(k,1522)
         lu(k,1526) = lu(k,1526) - lu(k,1434) * lu(k,1522)
         lu(k,1527) = lu(k,1527) - lu(k,1435) * lu(k,1522)
         lu(k,1528) = lu(k,1528) - lu(k,1436) * lu(k,1522)
         lu(k,1530) = lu(k,1530) - lu(k,1437) * lu(k,1522)
         lu(k,1532) = lu(k,1532) - lu(k,1438) * lu(k,1522)
         lu(k,1533) = lu(k,1533) - lu(k,1439) * lu(k,1522)
         lu(k,1535) = lu(k,1535) - lu(k,1440) * lu(k,1522)
         lu(k,1538) = lu(k,1538) - lu(k,1441) * lu(k,1522)
         lu(k,1539) = lu(k,1539) - lu(k,1442) * lu(k,1522)
         lu(k,1687) = lu(k,1687) - lu(k,1431) * lu(k,1686)
         lu(k,1688) = lu(k,1688) - lu(k,1432) * lu(k,1686)
         lu(k,1689) = lu(k,1689) - lu(k,1433) * lu(k,1686)
         lu(k,1690) = lu(k,1690) - lu(k,1434) * lu(k,1686)
         lu(k,1691) = lu(k,1691) - lu(k,1435) * lu(k,1686)
         lu(k,1692) = lu(k,1692) - lu(k,1436) * lu(k,1686)
         lu(k,1694) = lu(k,1694) - lu(k,1437) * lu(k,1686)
         lu(k,1696) = lu(k,1696) - lu(k,1438) * lu(k,1686)
         lu(k,1697) = lu(k,1697) - lu(k,1439) * lu(k,1686)
         lu(k,1699) = lu(k,1699) - lu(k,1440) * lu(k,1686)
         lu(k,1702) = lu(k,1702) - lu(k,1441) * lu(k,1686)
         lu(k,1703) = lu(k,1703) - lu(k,1442) * lu(k,1686)
         lu(k,1744) = lu(k,1744) - lu(k,1431) * lu(k,1743)
         lu(k,1745) = lu(k,1745) - lu(k,1432) * lu(k,1743)
         lu(k,1746) = lu(k,1746) - lu(k,1433) * lu(k,1743)
         lu(k,1747) = lu(k,1747) - lu(k,1434) * lu(k,1743)
         lu(k,1748) = lu(k,1748) - lu(k,1435) * lu(k,1743)
         lu(k,1749) = lu(k,1749) - lu(k,1436) * lu(k,1743)
         lu(k,1751) = lu(k,1751) - lu(k,1437) * lu(k,1743)
         lu(k,1753) = lu(k,1753) - lu(k,1438) * lu(k,1743)
         lu(k,1754) = lu(k,1754) - lu(k,1439) * lu(k,1743)
         lu(k,1756) = - lu(k,1440) * lu(k,1743)
         lu(k,1759) = lu(k,1759) - lu(k,1441) * lu(k,1743)
         lu(k,1760) = lu(k,1760) - lu(k,1442) * lu(k,1743)
         lu(k,1836) = lu(k,1836) - lu(k,1431) * lu(k,1835)
         lu(k,1837) = lu(k,1837) - lu(k,1432) * lu(k,1835)
         lu(k,1838) = lu(k,1838) - lu(k,1433) * lu(k,1835)
         lu(k,1839) = lu(k,1839) - lu(k,1434) * lu(k,1835)
         lu(k,1840) = lu(k,1840) - lu(k,1435) * lu(k,1835)
         lu(k,1841) = lu(k,1841) - lu(k,1436) * lu(k,1835)
         lu(k,1843) = lu(k,1843) - lu(k,1437) * lu(k,1835)
         lu(k,1845) = lu(k,1845) - lu(k,1438) * lu(k,1835)
         lu(k,1846) = lu(k,1846) - lu(k,1439) * lu(k,1835)
         lu(k,1848) = - lu(k,1440) * lu(k,1835)
         lu(k,1851) = lu(k,1851) - lu(k,1441) * lu(k,1835)
         lu(k,1852) = lu(k,1852) - lu(k,1442) * lu(k,1835)
         lu(k,1943) = lu(k,1943) - lu(k,1431) * lu(k,1942)
         lu(k,1944) = lu(k,1944) - lu(k,1432) * lu(k,1942)
         lu(k,1945) = lu(k,1945) - lu(k,1433) * lu(k,1942)
         lu(k,1946) = lu(k,1946) - lu(k,1434) * lu(k,1942)
         lu(k,1947) = lu(k,1947) - lu(k,1435) * lu(k,1942)
         lu(k,1948) = lu(k,1948) - lu(k,1436) * lu(k,1942)
         lu(k,1950) = lu(k,1950) - lu(k,1437) * lu(k,1942)
         lu(k,1952) = lu(k,1952) - lu(k,1438) * lu(k,1942)
         lu(k,1953) = lu(k,1953) - lu(k,1439) * lu(k,1942)
         lu(k,1955) = lu(k,1955) - lu(k,1440) * lu(k,1942)
         lu(k,1958) = lu(k,1958) - lu(k,1441) * lu(k,1942)
         lu(k,1959) = lu(k,1959) - lu(k,1442) * lu(k,1942)
         lu(k,2008) = lu(k,2008) - lu(k,1431) * lu(k,2007)
         lu(k,2009) = lu(k,2009) - lu(k,1432) * lu(k,2007)
         lu(k,2010) = lu(k,2010) - lu(k,1433) * lu(k,2007)
         lu(k,2011) = lu(k,2011) - lu(k,1434) * lu(k,2007)
         lu(k,2012) = lu(k,2012) - lu(k,1435) * lu(k,2007)
         lu(k,2013) = lu(k,2013) - lu(k,1436) * lu(k,2007)
         lu(k,2015) = lu(k,2015) - lu(k,1437) * lu(k,2007)
         lu(k,2017) = lu(k,2017) - lu(k,1438) * lu(k,2007)
         lu(k,2018) = lu(k,2018) - lu(k,1439) * lu(k,2007)
         lu(k,2020) = lu(k,2020) - lu(k,1440) * lu(k,2007)
         lu(k,2023) = lu(k,2023) - lu(k,1441) * lu(k,2007)
         lu(k,2024) = lu(k,2024) - lu(k,1442) * lu(k,2007)
         lu(k,2060) = lu(k,2060) - lu(k,1431) * lu(k,2059)
         lu(k,2061) = lu(k,2061) - lu(k,1432) * lu(k,2059)
         lu(k,2062) = lu(k,2062) - lu(k,1433) * lu(k,2059)
         lu(k,2063) = - lu(k,1434) * lu(k,2059)
         lu(k,2064) = lu(k,2064) - lu(k,1435) * lu(k,2059)
         lu(k,2065) = lu(k,2065) - lu(k,1436) * lu(k,2059)
         lu(k,2067) = lu(k,2067) - lu(k,1437) * lu(k,2059)
         lu(k,2069) = lu(k,2069) - lu(k,1438) * lu(k,2059)
         lu(k,2070) = lu(k,2070) - lu(k,1439) * lu(k,2059)
         lu(k,2072) = - lu(k,1440) * lu(k,2059)
         lu(k,2075) = lu(k,2075) - lu(k,1441) * lu(k,2059)
         lu(k,2076) = lu(k,2076) - lu(k,1442) * lu(k,2059)
         lu(k,2121) = lu(k,2121) - lu(k,1431) * lu(k,2120)
         lu(k,2122) = lu(k,2122) - lu(k,1432) * lu(k,2120)
         lu(k,2123) = lu(k,2123) - lu(k,1433) * lu(k,2120)
         lu(k,2124) = lu(k,2124) - lu(k,1434) * lu(k,2120)
         lu(k,2125) = lu(k,2125) - lu(k,1435) * lu(k,2120)
         lu(k,2126) = lu(k,2126) - lu(k,1436) * lu(k,2120)
         lu(k,2128) = lu(k,2128) - lu(k,1437) * lu(k,2120)
         lu(k,2130) = lu(k,2130) - lu(k,1438) * lu(k,2120)
         lu(k,2131) = lu(k,2131) - lu(k,1439) * lu(k,2120)
         lu(k,2133) = - lu(k,1440) * lu(k,2120)
         lu(k,2136) = lu(k,2136) - lu(k,1441) * lu(k,2120)
         lu(k,2137) = lu(k,2137) - lu(k,1442) * lu(k,2120)
         lu(k,2188) = lu(k,2188) - lu(k,1431) * lu(k,2187)
         lu(k,2189) = lu(k,2189) - lu(k,1432) * lu(k,2187)
         lu(k,2190) = lu(k,2190) - lu(k,1433) * lu(k,2187)
         lu(k,2191) = lu(k,2191) - lu(k,1434) * lu(k,2187)
         lu(k,2192) = lu(k,2192) - lu(k,1435) * lu(k,2187)
         lu(k,2193) = lu(k,2193) - lu(k,1436) * lu(k,2187)
         lu(k,2195) = lu(k,2195) - lu(k,1437) * lu(k,2187)
         lu(k,2197) = lu(k,2197) - lu(k,1438) * lu(k,2187)
         lu(k,2198) = lu(k,2198) - lu(k,1439) * lu(k,2187)
         lu(k,2200) = lu(k,2200) - lu(k,1440) * lu(k,2187)
         lu(k,2203) = lu(k,2203) - lu(k,1441) * lu(k,2187)
         lu(k,2204) = lu(k,2204) - lu(k,1442) * lu(k,2187)
         lu(k,2243) = lu(k,2243) - lu(k,1431) * lu(k,2242)
         lu(k,2244) = lu(k,2244) - lu(k,1432) * lu(k,2242)
         lu(k,2245) = lu(k,2245) - lu(k,1433) * lu(k,2242)
         lu(k,2246) = lu(k,2246) - lu(k,1434) * lu(k,2242)
         lu(k,2247) = lu(k,2247) - lu(k,1435) * lu(k,2242)
         lu(k,2248) = lu(k,2248) - lu(k,1436) * lu(k,2242)
         lu(k,2250) = lu(k,2250) - lu(k,1437) * lu(k,2242)
         lu(k,2252) = lu(k,2252) - lu(k,1438) * lu(k,2242)
         lu(k,2253) = lu(k,2253) - lu(k,1439) * lu(k,2242)
         lu(k,2255) = lu(k,2255) - lu(k,1440) * lu(k,2242)
         lu(k,2258) = lu(k,2258) - lu(k,1441) * lu(k,2242)
         lu(k,2259) = lu(k,2259) - lu(k,1442) * lu(k,2242)
         lu(k,2269) = lu(k,2269) - lu(k,1431) * lu(k,2268)
         lu(k,2270) = lu(k,2270) - lu(k,1432) * lu(k,2268)
         lu(k,2271) = lu(k,2271) - lu(k,1433) * lu(k,2268)
         lu(k,2272) = lu(k,2272) - lu(k,1434) * lu(k,2268)
         lu(k,2273) = lu(k,2273) - lu(k,1435) * lu(k,2268)
         lu(k,2274) = lu(k,2274) - lu(k,1436) * lu(k,2268)
         lu(k,2276) = lu(k,2276) - lu(k,1437) * lu(k,2268)
         lu(k,2278) = lu(k,2278) - lu(k,1438) * lu(k,2268)
         lu(k,2279) = lu(k,2279) - lu(k,1439) * lu(k,2268)
         lu(k,2281) = lu(k,2281) - lu(k,1440) * lu(k,2268)
         lu(k,2284) = lu(k,2284) - lu(k,1441) * lu(k,2268)
         lu(k,2285) = lu(k,2285) - lu(k,1442) * lu(k,2268)
                                                                        
         lu(k,1447) = 1._r8 / lu(k,1447)
         lu(k,1448) = lu(k,1448) * lu(k,1447)
         lu(k,1449) = lu(k,1449) * lu(k,1447)
         lu(k,1450) = lu(k,1450) * lu(k,1447)
         lu(k,1451) = lu(k,1451) * lu(k,1447)
         lu(k,1452) = lu(k,1452) * lu(k,1447)
         lu(k,1453) = lu(k,1453) * lu(k,1447)
         lu(k,1454) = lu(k,1454) * lu(k,1447)
         lu(k,1455) = lu(k,1455) * lu(k,1447)
         lu(k,1456) = lu(k,1456) * lu(k,1447)
         lu(k,1457) = lu(k,1457) * lu(k,1447)
         lu(k,1458) = lu(k,1458) * lu(k,1447)
         lu(k,1459) = lu(k,1459) * lu(k,1447)
         lu(k,1463) = lu(k,1463) - lu(k,1448) * lu(k,1462)
         lu(k,1464) = lu(k,1464) - lu(k,1449) * lu(k,1462)
         lu(k,1465) = lu(k,1465) - lu(k,1450) * lu(k,1462)
         lu(k,1466) = lu(k,1466) - lu(k,1451) * lu(k,1462)
         lu(k,1467) = lu(k,1467) - lu(k,1452) * lu(k,1462)
         lu(k,1468) = lu(k,1468) - lu(k,1453) * lu(k,1462)
         lu(k,1469) = lu(k,1469) - lu(k,1454) * lu(k,1462)
         lu(k,1470) = lu(k,1470) - lu(k,1455) * lu(k,1462)
         lu(k,1472) = lu(k,1472) - lu(k,1456) * lu(k,1462)
         lu(k,1473) = - lu(k,1457) * lu(k,1462)
         lu(k,1474) = lu(k,1474) - lu(k,1458) * lu(k,1462)
         lu(k,1475) = lu(k,1475) - lu(k,1459) * lu(k,1462)
         lu(k,1484) = lu(k,1484) - lu(k,1448) * lu(k,1483)
         lu(k,1485) = lu(k,1485) - lu(k,1449) * lu(k,1483)
         lu(k,1486) = lu(k,1486) - lu(k,1450) * lu(k,1483)
         lu(k,1487) = lu(k,1487) - lu(k,1451) * lu(k,1483)
         lu(k,1488) = lu(k,1488) - lu(k,1452) * lu(k,1483)
         lu(k,1490) = lu(k,1490) - lu(k,1453) * lu(k,1483)
         lu(k,1491) = lu(k,1491) - lu(k,1454) * lu(k,1483)
         lu(k,1492) = lu(k,1492) - lu(k,1455) * lu(k,1483)
         lu(k,1494) = lu(k,1494) - lu(k,1456) * lu(k,1483)
         lu(k,1495) = lu(k,1495) - lu(k,1457) * lu(k,1483)
         lu(k,1497) = lu(k,1497) - lu(k,1458) * lu(k,1483)
         lu(k,1498) = lu(k,1498) - lu(k,1459) * lu(k,1483)
         lu(k,1524) = lu(k,1524) - lu(k,1448) * lu(k,1523)
         lu(k,1525) = lu(k,1525) - lu(k,1449) * lu(k,1523)
         lu(k,1526) = lu(k,1526) - lu(k,1450) * lu(k,1523)
         lu(k,1527) = lu(k,1527) - lu(k,1451) * lu(k,1523)
         lu(k,1528) = lu(k,1528) - lu(k,1452) * lu(k,1523)
         lu(k,1530) = lu(k,1530) - lu(k,1453) * lu(k,1523)
         lu(k,1532) = lu(k,1532) - lu(k,1454) * lu(k,1523)
         lu(k,1533) = lu(k,1533) - lu(k,1455) * lu(k,1523)
         lu(k,1535) = lu(k,1535) - lu(k,1456) * lu(k,1523)
         lu(k,1536) = - lu(k,1457) * lu(k,1523)
         lu(k,1538) = lu(k,1538) - lu(k,1458) * lu(k,1523)
         lu(k,1539) = lu(k,1539) - lu(k,1459) * lu(k,1523)
         lu(k,1688) = lu(k,1688) - lu(k,1448) * lu(k,1687)
         lu(k,1689) = lu(k,1689) - lu(k,1449) * lu(k,1687)
         lu(k,1690) = lu(k,1690) - lu(k,1450) * lu(k,1687)
         lu(k,1691) = lu(k,1691) - lu(k,1451) * lu(k,1687)
         lu(k,1692) = lu(k,1692) - lu(k,1452) * lu(k,1687)
         lu(k,1694) = lu(k,1694) - lu(k,1453) * lu(k,1687)
         lu(k,1696) = lu(k,1696) - lu(k,1454) * lu(k,1687)
         lu(k,1697) = lu(k,1697) - lu(k,1455) * lu(k,1687)
         lu(k,1699) = lu(k,1699) - lu(k,1456) * lu(k,1687)
         lu(k,1700) = lu(k,1700) - lu(k,1457) * lu(k,1687)
         lu(k,1702) = lu(k,1702) - lu(k,1458) * lu(k,1687)
         lu(k,1703) = lu(k,1703) - lu(k,1459) * lu(k,1687)
         lu(k,1745) = lu(k,1745) - lu(k,1448) * lu(k,1744)
         lu(k,1746) = lu(k,1746) - lu(k,1449) * lu(k,1744)
         lu(k,1747) = lu(k,1747) - lu(k,1450) * lu(k,1744)
         lu(k,1748) = lu(k,1748) - lu(k,1451) * lu(k,1744)
         lu(k,1749) = lu(k,1749) - lu(k,1452) * lu(k,1744)
         lu(k,1751) = lu(k,1751) - lu(k,1453) * lu(k,1744)
         lu(k,1753) = lu(k,1753) - lu(k,1454) * lu(k,1744)
         lu(k,1754) = lu(k,1754) - lu(k,1455) * lu(k,1744)
         lu(k,1756) = lu(k,1756) - lu(k,1456) * lu(k,1744)
         lu(k,1757) = lu(k,1757) - lu(k,1457) * lu(k,1744)
         lu(k,1759) = lu(k,1759) - lu(k,1458) * lu(k,1744)
         lu(k,1760) = lu(k,1760) - lu(k,1459) * lu(k,1744)
         lu(k,1837) = lu(k,1837) - lu(k,1448) * lu(k,1836)
         lu(k,1838) = lu(k,1838) - lu(k,1449) * lu(k,1836)
         lu(k,1839) = lu(k,1839) - lu(k,1450) * lu(k,1836)
         lu(k,1840) = lu(k,1840) - lu(k,1451) * lu(k,1836)
         lu(k,1841) = lu(k,1841) - lu(k,1452) * lu(k,1836)
         lu(k,1843) = lu(k,1843) - lu(k,1453) * lu(k,1836)
         lu(k,1845) = lu(k,1845) - lu(k,1454) * lu(k,1836)
         lu(k,1846) = lu(k,1846) - lu(k,1455) * lu(k,1836)
         lu(k,1848) = lu(k,1848) - lu(k,1456) * lu(k,1836)
         lu(k,1849) = lu(k,1849) - lu(k,1457) * lu(k,1836)
         lu(k,1851) = lu(k,1851) - lu(k,1458) * lu(k,1836)
         lu(k,1852) = lu(k,1852) - lu(k,1459) * lu(k,1836)
         lu(k,1944) = lu(k,1944) - lu(k,1448) * lu(k,1943)
         lu(k,1945) = lu(k,1945) - lu(k,1449) * lu(k,1943)
         lu(k,1946) = lu(k,1946) - lu(k,1450) * lu(k,1943)
         lu(k,1947) = lu(k,1947) - lu(k,1451) * lu(k,1943)
         lu(k,1948) = lu(k,1948) - lu(k,1452) * lu(k,1943)
         lu(k,1950) = lu(k,1950) - lu(k,1453) * lu(k,1943)
         lu(k,1952) = lu(k,1952) - lu(k,1454) * lu(k,1943)
         lu(k,1953) = lu(k,1953) - lu(k,1455) * lu(k,1943)
         lu(k,1955) = lu(k,1955) - lu(k,1456) * lu(k,1943)
         lu(k,1956) = lu(k,1956) - lu(k,1457) * lu(k,1943)
         lu(k,1958) = lu(k,1958) - lu(k,1458) * lu(k,1943)
         lu(k,1959) = lu(k,1959) - lu(k,1459) * lu(k,1943)
         lu(k,1970) = lu(k,1970) - lu(k,1448) * lu(k,1969)
         lu(k,1971) = lu(k,1971) - lu(k,1449) * lu(k,1969)
         lu(k,1972) = lu(k,1972) - lu(k,1450) * lu(k,1969)
         lu(k,1973) = lu(k,1973) - lu(k,1451) * lu(k,1969)
         lu(k,1974) = lu(k,1974) - lu(k,1452) * lu(k,1969)
         lu(k,1976) = lu(k,1976) - lu(k,1453) * lu(k,1969)
         lu(k,1978) = lu(k,1978) - lu(k,1454) * lu(k,1969)
         lu(k,1979) = lu(k,1979) - lu(k,1455) * lu(k,1969)
         lu(k,1981) = lu(k,1981) - lu(k,1456) * lu(k,1969)
         lu(k,1982) = lu(k,1982) - lu(k,1457) * lu(k,1969)
         lu(k,1984) = lu(k,1984) - lu(k,1458) * lu(k,1969)
         lu(k,1985) = lu(k,1985) - lu(k,1459) * lu(k,1969)
         lu(k,2009) = lu(k,2009) - lu(k,1448) * lu(k,2008)
         lu(k,2010) = lu(k,2010) - lu(k,1449) * lu(k,2008)
         lu(k,2011) = lu(k,2011) - lu(k,1450) * lu(k,2008)
         lu(k,2012) = lu(k,2012) - lu(k,1451) * lu(k,2008)
         lu(k,2013) = lu(k,2013) - lu(k,1452) * lu(k,2008)
         lu(k,2015) = lu(k,2015) - lu(k,1453) * lu(k,2008)
         lu(k,2017) = lu(k,2017) - lu(k,1454) * lu(k,2008)
         lu(k,2018) = lu(k,2018) - lu(k,1455) * lu(k,2008)
         lu(k,2020) = lu(k,2020) - lu(k,1456) * lu(k,2008)
         lu(k,2021) = lu(k,2021) - lu(k,1457) * lu(k,2008)
         lu(k,2023) = lu(k,2023) - lu(k,1458) * lu(k,2008)
         lu(k,2024) = lu(k,2024) - lu(k,1459) * lu(k,2008)
         lu(k,2061) = lu(k,2061) - lu(k,1448) * lu(k,2060)
         lu(k,2062) = lu(k,2062) - lu(k,1449) * lu(k,2060)
         lu(k,2063) = lu(k,2063) - lu(k,1450) * lu(k,2060)
         lu(k,2064) = lu(k,2064) - lu(k,1451) * lu(k,2060)
         lu(k,2065) = lu(k,2065) - lu(k,1452) * lu(k,2060)
         lu(k,2067) = lu(k,2067) - lu(k,1453) * lu(k,2060)
         lu(k,2069) = lu(k,2069) - lu(k,1454) * lu(k,2060)
         lu(k,2070) = lu(k,2070) - lu(k,1455) * lu(k,2060)
         lu(k,2072) = lu(k,2072) - lu(k,1456) * lu(k,2060)
         lu(k,2073) = lu(k,2073) - lu(k,1457) * lu(k,2060)
         lu(k,2075) = lu(k,2075) - lu(k,1458) * lu(k,2060)
         lu(k,2076) = lu(k,2076) - lu(k,1459) * lu(k,2060)
         lu(k,2122) = lu(k,2122) - lu(k,1448) * lu(k,2121)
         lu(k,2123) = lu(k,2123) - lu(k,1449) * lu(k,2121)
         lu(k,2124) = lu(k,2124) - lu(k,1450) * lu(k,2121)
         lu(k,2125) = lu(k,2125) - lu(k,1451) * lu(k,2121)
         lu(k,2126) = lu(k,2126) - lu(k,1452) * lu(k,2121)
         lu(k,2128) = lu(k,2128) - lu(k,1453) * lu(k,2121)
         lu(k,2130) = lu(k,2130) - lu(k,1454) * lu(k,2121)
         lu(k,2131) = lu(k,2131) - lu(k,1455) * lu(k,2121)
         lu(k,2133) = lu(k,2133) - lu(k,1456) * lu(k,2121)
         lu(k,2134) = lu(k,2134) - lu(k,1457) * lu(k,2121)
         lu(k,2136) = lu(k,2136) - lu(k,1458) * lu(k,2121)
         lu(k,2137) = lu(k,2137) - lu(k,1459) * lu(k,2121)
         lu(k,2145) = lu(k,2145) - lu(k,1448) * lu(k,2144)
         lu(k,2146) = lu(k,2146) - lu(k,1449) * lu(k,2144)
         lu(k,2147) = lu(k,2147) - lu(k,1450) * lu(k,2144)
         lu(k,2148) = lu(k,2148) - lu(k,1451) * lu(k,2144)
         lu(k,2149) = lu(k,2149) - lu(k,1452) * lu(k,2144)
         lu(k,2151) = lu(k,2151) - lu(k,1453) * lu(k,2144)
         lu(k,2153) = lu(k,2153) - lu(k,1454) * lu(k,2144)
         lu(k,2154) = - lu(k,1455) * lu(k,2144)
         lu(k,2156) = lu(k,2156) - lu(k,1456) * lu(k,2144)
         lu(k,2157) = lu(k,2157) - lu(k,1457) * lu(k,2144)
         lu(k,2159) = lu(k,2159) - lu(k,1458) * lu(k,2144)
         lu(k,2160) = lu(k,2160) - lu(k,1459) * lu(k,2144)
         lu(k,2189) = lu(k,2189) - lu(k,1448) * lu(k,2188)
         lu(k,2190) = lu(k,2190) - lu(k,1449) * lu(k,2188)
         lu(k,2191) = lu(k,2191) - lu(k,1450) * lu(k,2188)
         lu(k,2192) = lu(k,2192) - lu(k,1451) * lu(k,2188)
         lu(k,2193) = lu(k,2193) - lu(k,1452) * lu(k,2188)
         lu(k,2195) = lu(k,2195) - lu(k,1453) * lu(k,2188)
         lu(k,2197) = lu(k,2197) - lu(k,1454) * lu(k,2188)
         lu(k,2198) = lu(k,2198) - lu(k,1455) * lu(k,2188)
         lu(k,2200) = lu(k,2200) - lu(k,1456) * lu(k,2188)
         lu(k,2201) = lu(k,2201) - lu(k,1457) * lu(k,2188)
         lu(k,2203) = lu(k,2203) - lu(k,1458) * lu(k,2188)
         lu(k,2204) = lu(k,2204) - lu(k,1459) * lu(k,2188)
         lu(k,2213) = lu(k,2213) - lu(k,1448) * lu(k,2212)
         lu(k,2214) = lu(k,2214) - lu(k,1449) * lu(k,2212)
         lu(k,2215) = lu(k,2215) - lu(k,1450) * lu(k,2212)
         lu(k,2216) = lu(k,2216) - lu(k,1451) * lu(k,2212)
         lu(k,2217) = lu(k,2217) - lu(k,1452) * lu(k,2212)
         lu(k,2219) = lu(k,2219) - lu(k,1453) * lu(k,2212)
         lu(k,2221) = lu(k,2221) - lu(k,1454) * lu(k,2212)
         lu(k,2222) = - lu(k,1455) * lu(k,2212)
         lu(k,2224) = lu(k,2224) - lu(k,1456) * lu(k,2212)
         lu(k,2225) = lu(k,2225) - lu(k,1457) * lu(k,2212)
         lu(k,2227) = lu(k,2227) - lu(k,1458) * lu(k,2212)
         lu(k,2228) = lu(k,2228) - lu(k,1459) * lu(k,2212)
         lu(k,2244) = lu(k,2244) - lu(k,1448) * lu(k,2243)
         lu(k,2245) = lu(k,2245) - lu(k,1449) * lu(k,2243)
         lu(k,2246) = lu(k,2246) - lu(k,1450) * lu(k,2243)
         lu(k,2247) = lu(k,2247) - lu(k,1451) * lu(k,2243)
         lu(k,2248) = lu(k,2248) - lu(k,1452) * lu(k,2243)
         lu(k,2250) = lu(k,2250) - lu(k,1453) * lu(k,2243)
         lu(k,2252) = lu(k,2252) - lu(k,1454) * lu(k,2243)
         lu(k,2253) = lu(k,2253) - lu(k,1455) * lu(k,2243)
         lu(k,2255) = lu(k,2255) - lu(k,1456) * lu(k,2243)
         lu(k,2256) = lu(k,2256) - lu(k,1457) * lu(k,2243)
         lu(k,2258) = lu(k,2258) - lu(k,1458) * lu(k,2243)
         lu(k,2259) = lu(k,2259) - lu(k,1459) * lu(k,2243)
         lu(k,2270) = lu(k,2270) - lu(k,1448) * lu(k,2269)
         lu(k,2271) = lu(k,2271) - lu(k,1449) * lu(k,2269)
         lu(k,2272) = lu(k,2272) - lu(k,1450) * lu(k,2269)
         lu(k,2273) = lu(k,2273) - lu(k,1451) * lu(k,2269)
         lu(k,2274) = lu(k,2274) - lu(k,1452) * lu(k,2269)
         lu(k,2276) = lu(k,2276) - lu(k,1453) * lu(k,2269)
         lu(k,2278) = lu(k,2278) - lu(k,1454) * lu(k,2269)
         lu(k,2279) = lu(k,2279) - lu(k,1455) * lu(k,2269)
         lu(k,2281) = lu(k,2281) - lu(k,1456) * lu(k,2269)
         lu(k,2282) = lu(k,2282) - lu(k,1457) * lu(k,2269)
         lu(k,2284) = lu(k,2284) - lu(k,1458) * lu(k,2269)
         lu(k,2285) = lu(k,2285) - lu(k,1459) * lu(k,2269)
                                                                        
         lu(k,1463) = 1._r8 / lu(k,1463)
         lu(k,1464) = lu(k,1464) * lu(k,1463)
         lu(k,1465) = lu(k,1465) * lu(k,1463)
         lu(k,1466) = lu(k,1466) * lu(k,1463)
         lu(k,1467) = lu(k,1467) * lu(k,1463)
         lu(k,1468) = lu(k,1468) * lu(k,1463)
         lu(k,1469) = lu(k,1469) * lu(k,1463)
         lu(k,1470) = lu(k,1470) * lu(k,1463)
         lu(k,1471) = lu(k,1471) * lu(k,1463)
         lu(k,1472) = lu(k,1472) * lu(k,1463)
         lu(k,1473) = lu(k,1473) * lu(k,1463)
         lu(k,1474) = lu(k,1474) * lu(k,1463)
         lu(k,1475) = lu(k,1475) * lu(k,1463)
         lu(k,1485) = lu(k,1485) - lu(k,1464) * lu(k,1484)
         lu(k,1486) = lu(k,1486) - lu(k,1465) * lu(k,1484)
         lu(k,1487) = lu(k,1487) - lu(k,1466) * lu(k,1484)
         lu(k,1488) = lu(k,1488) - lu(k,1467) * lu(k,1484)
         lu(k,1490) = lu(k,1490) - lu(k,1468) * lu(k,1484)
         lu(k,1491) = lu(k,1491) - lu(k,1469) * lu(k,1484)
         lu(k,1492) = lu(k,1492) - lu(k,1470) * lu(k,1484)
         lu(k,1493) = lu(k,1493) - lu(k,1471) * lu(k,1484)
         lu(k,1494) = lu(k,1494) - lu(k,1472) * lu(k,1484)
         lu(k,1495) = lu(k,1495) - lu(k,1473) * lu(k,1484)
         lu(k,1497) = lu(k,1497) - lu(k,1474) * lu(k,1484)
         lu(k,1498) = lu(k,1498) - lu(k,1475) * lu(k,1484)
         lu(k,1525) = lu(k,1525) - lu(k,1464) * lu(k,1524)
         lu(k,1526) = lu(k,1526) - lu(k,1465) * lu(k,1524)
         lu(k,1527) = lu(k,1527) - lu(k,1466) * lu(k,1524)
         lu(k,1528) = lu(k,1528) - lu(k,1467) * lu(k,1524)
         lu(k,1530) = lu(k,1530) - lu(k,1468) * lu(k,1524)
         lu(k,1532) = lu(k,1532) - lu(k,1469) * lu(k,1524)
         lu(k,1533) = lu(k,1533) - lu(k,1470) * lu(k,1524)
         lu(k,1534) = lu(k,1534) - lu(k,1471) * lu(k,1524)
         lu(k,1535) = lu(k,1535) - lu(k,1472) * lu(k,1524)
         lu(k,1536) = lu(k,1536) - lu(k,1473) * lu(k,1524)
         lu(k,1538) = lu(k,1538) - lu(k,1474) * lu(k,1524)
         lu(k,1539) = lu(k,1539) - lu(k,1475) * lu(k,1524)
         lu(k,1689) = lu(k,1689) - lu(k,1464) * lu(k,1688)
         lu(k,1690) = lu(k,1690) - lu(k,1465) * lu(k,1688)
         lu(k,1691) = lu(k,1691) - lu(k,1466) * lu(k,1688)
         lu(k,1692) = lu(k,1692) - lu(k,1467) * lu(k,1688)
         lu(k,1694) = lu(k,1694) - lu(k,1468) * lu(k,1688)
         lu(k,1696) = lu(k,1696) - lu(k,1469) * lu(k,1688)
         lu(k,1697) = lu(k,1697) - lu(k,1470) * lu(k,1688)
         lu(k,1698) = lu(k,1698) - lu(k,1471) * lu(k,1688)
         lu(k,1699) = lu(k,1699) - lu(k,1472) * lu(k,1688)
         lu(k,1700) = lu(k,1700) - lu(k,1473) * lu(k,1688)
         lu(k,1702) = lu(k,1702) - lu(k,1474) * lu(k,1688)
         lu(k,1703) = lu(k,1703) - lu(k,1475) * lu(k,1688)
         lu(k,1746) = lu(k,1746) - lu(k,1464) * lu(k,1745)
         lu(k,1747) = lu(k,1747) - lu(k,1465) * lu(k,1745)
         lu(k,1748) = lu(k,1748) - lu(k,1466) * lu(k,1745)
         lu(k,1749) = lu(k,1749) - lu(k,1467) * lu(k,1745)
         lu(k,1751) = lu(k,1751) - lu(k,1468) * lu(k,1745)
         lu(k,1753) = lu(k,1753) - lu(k,1469) * lu(k,1745)
         lu(k,1754) = lu(k,1754) - lu(k,1470) * lu(k,1745)
         lu(k,1755) = lu(k,1755) - lu(k,1471) * lu(k,1745)
         lu(k,1756) = lu(k,1756) - lu(k,1472) * lu(k,1745)
         lu(k,1757) = lu(k,1757) - lu(k,1473) * lu(k,1745)
         lu(k,1759) = lu(k,1759) - lu(k,1474) * lu(k,1745)
         lu(k,1760) = lu(k,1760) - lu(k,1475) * lu(k,1745)
         lu(k,1838) = lu(k,1838) - lu(k,1464) * lu(k,1837)
         lu(k,1839) = lu(k,1839) - lu(k,1465) * lu(k,1837)
         lu(k,1840) = lu(k,1840) - lu(k,1466) * lu(k,1837)
         lu(k,1841) = lu(k,1841) - lu(k,1467) * lu(k,1837)
         lu(k,1843) = lu(k,1843) - lu(k,1468) * lu(k,1837)
         lu(k,1845) = lu(k,1845) - lu(k,1469) * lu(k,1837)
         lu(k,1846) = lu(k,1846) - lu(k,1470) * lu(k,1837)
         lu(k,1847) = lu(k,1847) - lu(k,1471) * lu(k,1837)
         lu(k,1848) = lu(k,1848) - lu(k,1472) * lu(k,1837)
         lu(k,1849) = lu(k,1849) - lu(k,1473) * lu(k,1837)
         lu(k,1851) = lu(k,1851) - lu(k,1474) * lu(k,1837)
         lu(k,1852) = lu(k,1852) - lu(k,1475) * lu(k,1837)
         lu(k,1945) = lu(k,1945) - lu(k,1464) * lu(k,1944)
         lu(k,1946) = lu(k,1946) - lu(k,1465) * lu(k,1944)
         lu(k,1947) = lu(k,1947) - lu(k,1466) * lu(k,1944)
         lu(k,1948) = lu(k,1948) - lu(k,1467) * lu(k,1944)
         lu(k,1950) = lu(k,1950) - lu(k,1468) * lu(k,1944)
         lu(k,1952) = lu(k,1952) - lu(k,1469) * lu(k,1944)
         lu(k,1953) = lu(k,1953) - lu(k,1470) * lu(k,1944)
         lu(k,1954) = lu(k,1954) - lu(k,1471) * lu(k,1944)
         lu(k,1955) = lu(k,1955) - lu(k,1472) * lu(k,1944)
         lu(k,1956) = lu(k,1956) - lu(k,1473) * lu(k,1944)
         lu(k,1958) = lu(k,1958) - lu(k,1474) * lu(k,1944)
         lu(k,1959) = lu(k,1959) - lu(k,1475) * lu(k,1944)
         lu(k,1971) = lu(k,1971) - lu(k,1464) * lu(k,1970)
         lu(k,1972) = lu(k,1972) - lu(k,1465) * lu(k,1970)
         lu(k,1973) = lu(k,1973) - lu(k,1466) * lu(k,1970)
         lu(k,1974) = lu(k,1974) - lu(k,1467) * lu(k,1970)
         lu(k,1976) = lu(k,1976) - lu(k,1468) * lu(k,1970)
         lu(k,1978) = lu(k,1978) - lu(k,1469) * lu(k,1970)
         lu(k,1979) = lu(k,1979) - lu(k,1470) * lu(k,1970)
         lu(k,1980) = lu(k,1980) - lu(k,1471) * lu(k,1970)
         lu(k,1981) = lu(k,1981) - lu(k,1472) * lu(k,1970)
         lu(k,1982) = lu(k,1982) - lu(k,1473) * lu(k,1970)
         lu(k,1984) = lu(k,1984) - lu(k,1474) * lu(k,1970)
         lu(k,1985) = lu(k,1985) - lu(k,1475) * lu(k,1970)
         lu(k,2010) = lu(k,2010) - lu(k,1464) * lu(k,2009)
         lu(k,2011) = lu(k,2011) - lu(k,1465) * lu(k,2009)
         lu(k,2012) = lu(k,2012) - lu(k,1466) * lu(k,2009)
         lu(k,2013) = lu(k,2013) - lu(k,1467) * lu(k,2009)
         lu(k,2015) = lu(k,2015) - lu(k,1468) * lu(k,2009)
         lu(k,2017) = lu(k,2017) - lu(k,1469) * lu(k,2009)
         lu(k,2018) = lu(k,2018) - lu(k,1470) * lu(k,2009)
         lu(k,2019) = lu(k,2019) - lu(k,1471) * lu(k,2009)
         lu(k,2020) = lu(k,2020) - lu(k,1472) * lu(k,2009)
         lu(k,2021) = lu(k,2021) - lu(k,1473) * lu(k,2009)
         lu(k,2023) = lu(k,2023) - lu(k,1474) * lu(k,2009)
         lu(k,2024) = lu(k,2024) - lu(k,1475) * lu(k,2009)
         lu(k,2062) = lu(k,2062) - lu(k,1464) * lu(k,2061)
         lu(k,2063) = lu(k,2063) - lu(k,1465) * lu(k,2061)
         lu(k,2064) = lu(k,2064) - lu(k,1466) * lu(k,2061)
         lu(k,2065) = lu(k,2065) - lu(k,1467) * lu(k,2061)
         lu(k,2067) = lu(k,2067) - lu(k,1468) * lu(k,2061)
         lu(k,2069) = lu(k,2069) - lu(k,1469) * lu(k,2061)
         lu(k,2070) = lu(k,2070) - lu(k,1470) * lu(k,2061)
         lu(k,2071) = lu(k,2071) - lu(k,1471) * lu(k,2061)
         lu(k,2072) = lu(k,2072) - lu(k,1472) * lu(k,2061)
         lu(k,2073) = lu(k,2073) - lu(k,1473) * lu(k,2061)
         lu(k,2075) = lu(k,2075) - lu(k,1474) * lu(k,2061)
         lu(k,2076) = lu(k,2076) - lu(k,1475) * lu(k,2061)
         lu(k,2123) = lu(k,2123) - lu(k,1464) * lu(k,2122)
         lu(k,2124) = lu(k,2124) - lu(k,1465) * lu(k,2122)
         lu(k,2125) = lu(k,2125) - lu(k,1466) * lu(k,2122)
         lu(k,2126) = lu(k,2126) - lu(k,1467) * lu(k,2122)
         lu(k,2128) = lu(k,2128) - lu(k,1468) * lu(k,2122)
         lu(k,2130) = lu(k,2130) - lu(k,1469) * lu(k,2122)
         lu(k,2131) = lu(k,2131) - lu(k,1470) * lu(k,2122)
         lu(k,2132) = lu(k,2132) - lu(k,1471) * lu(k,2122)
         lu(k,2133) = lu(k,2133) - lu(k,1472) * lu(k,2122)
         lu(k,2134) = lu(k,2134) - lu(k,1473) * lu(k,2122)
         lu(k,2136) = lu(k,2136) - lu(k,1474) * lu(k,2122)
         lu(k,2137) = lu(k,2137) - lu(k,1475) * lu(k,2122)
         lu(k,2146) = lu(k,2146) - lu(k,1464) * lu(k,2145)
         lu(k,2147) = lu(k,2147) - lu(k,1465) * lu(k,2145)
         lu(k,2148) = lu(k,2148) - lu(k,1466) * lu(k,2145)
         lu(k,2149) = lu(k,2149) - lu(k,1467) * lu(k,2145)
         lu(k,2151) = lu(k,2151) - lu(k,1468) * lu(k,2145)
         lu(k,2153) = lu(k,2153) - lu(k,1469) * lu(k,2145)
         lu(k,2154) = lu(k,2154) - lu(k,1470) * lu(k,2145)
         lu(k,2155) = lu(k,2155) - lu(k,1471) * lu(k,2145)
         lu(k,2156) = lu(k,2156) - lu(k,1472) * lu(k,2145)
         lu(k,2157) = lu(k,2157) - lu(k,1473) * lu(k,2145)
         lu(k,2159) = lu(k,2159) - lu(k,1474) * lu(k,2145)
         lu(k,2160) = lu(k,2160) - lu(k,1475) * lu(k,2145)
         lu(k,2190) = lu(k,2190) - lu(k,1464) * lu(k,2189)
         lu(k,2191) = lu(k,2191) - lu(k,1465) * lu(k,2189)
         lu(k,2192) = lu(k,2192) - lu(k,1466) * lu(k,2189)
         lu(k,2193) = lu(k,2193) - lu(k,1467) * lu(k,2189)
         lu(k,2195) = lu(k,2195) - lu(k,1468) * lu(k,2189)
         lu(k,2197) = lu(k,2197) - lu(k,1469) * lu(k,2189)
         lu(k,2198) = lu(k,2198) - lu(k,1470) * lu(k,2189)
         lu(k,2199) = lu(k,2199) - lu(k,1471) * lu(k,2189)
         lu(k,2200) = lu(k,2200) - lu(k,1472) * lu(k,2189)
         lu(k,2201) = lu(k,2201) - lu(k,1473) * lu(k,2189)
         lu(k,2203) = lu(k,2203) - lu(k,1474) * lu(k,2189)
         lu(k,2204) = lu(k,2204) - lu(k,1475) * lu(k,2189)
         lu(k,2214) = lu(k,2214) - lu(k,1464) * lu(k,2213)
         lu(k,2215) = lu(k,2215) - lu(k,1465) * lu(k,2213)
         lu(k,2216) = lu(k,2216) - lu(k,1466) * lu(k,2213)
         lu(k,2217) = lu(k,2217) - lu(k,1467) * lu(k,2213)
         lu(k,2219) = lu(k,2219) - lu(k,1468) * lu(k,2213)
         lu(k,2221) = lu(k,2221) - lu(k,1469) * lu(k,2213)
         lu(k,2222) = lu(k,2222) - lu(k,1470) * lu(k,2213)
         lu(k,2223) = lu(k,2223) - lu(k,1471) * lu(k,2213)
         lu(k,2224) = lu(k,2224) - lu(k,1472) * lu(k,2213)
         lu(k,2225) = lu(k,2225) - lu(k,1473) * lu(k,2213)
         lu(k,2227) = lu(k,2227) - lu(k,1474) * lu(k,2213)
         lu(k,2228) = lu(k,2228) - lu(k,1475) * lu(k,2213)
         lu(k,2245) = lu(k,2245) - lu(k,1464) * lu(k,2244)
         lu(k,2246) = lu(k,2246) - lu(k,1465) * lu(k,2244)
         lu(k,2247) = lu(k,2247) - lu(k,1466) * lu(k,2244)
         lu(k,2248) = lu(k,2248) - lu(k,1467) * lu(k,2244)
         lu(k,2250) = lu(k,2250) - lu(k,1468) * lu(k,2244)
         lu(k,2252) = lu(k,2252) - lu(k,1469) * lu(k,2244)
         lu(k,2253) = lu(k,2253) - lu(k,1470) * lu(k,2244)
         lu(k,2254) = lu(k,2254) - lu(k,1471) * lu(k,2244)
         lu(k,2255) = lu(k,2255) - lu(k,1472) * lu(k,2244)
         lu(k,2256) = lu(k,2256) - lu(k,1473) * lu(k,2244)
         lu(k,2258) = lu(k,2258) - lu(k,1474) * lu(k,2244)
         lu(k,2259) = lu(k,2259) - lu(k,1475) * lu(k,2244)
         lu(k,2271) = lu(k,2271) - lu(k,1464) * lu(k,2270)
         lu(k,2272) = lu(k,2272) - lu(k,1465) * lu(k,2270)
         lu(k,2273) = lu(k,2273) - lu(k,1466) * lu(k,2270)
         lu(k,2274) = lu(k,2274) - lu(k,1467) * lu(k,2270)
         lu(k,2276) = lu(k,2276) - lu(k,1468) * lu(k,2270)
         lu(k,2278) = lu(k,2278) - lu(k,1469) * lu(k,2270)
         lu(k,2279) = lu(k,2279) - lu(k,1470) * lu(k,2270)
         lu(k,2280) = lu(k,2280) - lu(k,1471) * lu(k,2270)
         lu(k,2281) = lu(k,2281) - lu(k,1472) * lu(k,2270)
         lu(k,2282) = lu(k,2282) - lu(k,1473) * lu(k,2270)
         lu(k,2284) = lu(k,2284) - lu(k,1474) * lu(k,2270)
         lu(k,2285) = lu(k,2285) - lu(k,1475) * lu(k,2270)
                                                                        
      end do
                                                                        
      end subroutine lu_fac28
                                                                        
      subroutine lu_fac29( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,1485) = 1._r8 / lu(k,1485)
         lu(k,1486) = lu(k,1486) * lu(k,1485)
         lu(k,1487) = lu(k,1487) * lu(k,1485)
         lu(k,1488) = lu(k,1488) * lu(k,1485)
         lu(k,1489) = lu(k,1489) * lu(k,1485)
         lu(k,1490) = lu(k,1490) * lu(k,1485)
         lu(k,1491) = lu(k,1491) * lu(k,1485)
         lu(k,1492) = lu(k,1492) * lu(k,1485)
         lu(k,1493) = lu(k,1493) * lu(k,1485)
         lu(k,1494) = lu(k,1494) * lu(k,1485)
         lu(k,1495) = lu(k,1495) * lu(k,1485)
         lu(k,1496) = lu(k,1496) * lu(k,1485)
         lu(k,1497) = lu(k,1497) * lu(k,1485)
         lu(k,1498) = lu(k,1498) * lu(k,1485)
         lu(k,1526) = lu(k,1526) - lu(k,1486) * lu(k,1525)
         lu(k,1527) = lu(k,1527) - lu(k,1487) * lu(k,1525)
         lu(k,1528) = lu(k,1528) - lu(k,1488) * lu(k,1525)
         lu(k,1529) = lu(k,1529) - lu(k,1489) * lu(k,1525)
         lu(k,1530) = lu(k,1530) - lu(k,1490) * lu(k,1525)
         lu(k,1532) = lu(k,1532) - lu(k,1491) * lu(k,1525)
         lu(k,1533) = lu(k,1533) - lu(k,1492) * lu(k,1525)
         lu(k,1534) = lu(k,1534) - lu(k,1493) * lu(k,1525)
         lu(k,1535) = lu(k,1535) - lu(k,1494) * lu(k,1525)
         lu(k,1536) = lu(k,1536) - lu(k,1495) * lu(k,1525)
         lu(k,1537) = lu(k,1537) - lu(k,1496) * lu(k,1525)
         lu(k,1538) = lu(k,1538) - lu(k,1497) * lu(k,1525)
         lu(k,1539) = lu(k,1539) - lu(k,1498) * lu(k,1525)
         lu(k,1690) = lu(k,1690) - lu(k,1486) * lu(k,1689)
         lu(k,1691) = lu(k,1691) - lu(k,1487) * lu(k,1689)
         lu(k,1692) = lu(k,1692) - lu(k,1488) * lu(k,1689)
         lu(k,1693) = lu(k,1693) - lu(k,1489) * lu(k,1689)
         lu(k,1694) = lu(k,1694) - lu(k,1490) * lu(k,1689)
         lu(k,1696) = lu(k,1696) - lu(k,1491) * lu(k,1689)
         lu(k,1697) = lu(k,1697) - lu(k,1492) * lu(k,1689)
         lu(k,1698) = lu(k,1698) - lu(k,1493) * lu(k,1689)
         lu(k,1699) = lu(k,1699) - lu(k,1494) * lu(k,1689)
         lu(k,1700) = lu(k,1700) - lu(k,1495) * lu(k,1689)
         lu(k,1701) = lu(k,1701) - lu(k,1496) * lu(k,1689)
         lu(k,1702) = lu(k,1702) - lu(k,1497) * lu(k,1689)
         lu(k,1703) = lu(k,1703) - lu(k,1498) * lu(k,1689)
         lu(k,1747) = lu(k,1747) - lu(k,1486) * lu(k,1746)
         lu(k,1748) = lu(k,1748) - lu(k,1487) * lu(k,1746)
         lu(k,1749) = lu(k,1749) - lu(k,1488) * lu(k,1746)
         lu(k,1750) = lu(k,1750) - lu(k,1489) * lu(k,1746)
         lu(k,1751) = lu(k,1751) - lu(k,1490) * lu(k,1746)
         lu(k,1753) = lu(k,1753) - lu(k,1491) * lu(k,1746)
         lu(k,1754) = lu(k,1754) - lu(k,1492) * lu(k,1746)
         lu(k,1755) = lu(k,1755) - lu(k,1493) * lu(k,1746)
         lu(k,1756) = lu(k,1756) - lu(k,1494) * lu(k,1746)
         lu(k,1757) = lu(k,1757) - lu(k,1495) * lu(k,1746)
         lu(k,1758) = lu(k,1758) - lu(k,1496) * lu(k,1746)
         lu(k,1759) = lu(k,1759) - lu(k,1497) * lu(k,1746)
         lu(k,1760) = lu(k,1760) - lu(k,1498) * lu(k,1746)
         lu(k,1839) = lu(k,1839) - lu(k,1486) * lu(k,1838)
         lu(k,1840) = lu(k,1840) - lu(k,1487) * lu(k,1838)
         lu(k,1841) = lu(k,1841) - lu(k,1488) * lu(k,1838)
         lu(k,1842) = lu(k,1842) - lu(k,1489) * lu(k,1838)
         lu(k,1843) = lu(k,1843) - lu(k,1490) * lu(k,1838)
         lu(k,1845) = lu(k,1845) - lu(k,1491) * lu(k,1838)
         lu(k,1846) = lu(k,1846) - lu(k,1492) * lu(k,1838)
         lu(k,1847) = lu(k,1847) - lu(k,1493) * lu(k,1838)
         lu(k,1848) = lu(k,1848) - lu(k,1494) * lu(k,1838)
         lu(k,1849) = lu(k,1849) - lu(k,1495) * lu(k,1838)
         lu(k,1850) = lu(k,1850) - lu(k,1496) * lu(k,1838)
         lu(k,1851) = lu(k,1851) - lu(k,1497) * lu(k,1838)
         lu(k,1852) = lu(k,1852) - lu(k,1498) * lu(k,1838)
         lu(k,1946) = lu(k,1946) - lu(k,1486) * lu(k,1945)
         lu(k,1947) = lu(k,1947) - lu(k,1487) * lu(k,1945)
         lu(k,1948) = lu(k,1948) - lu(k,1488) * lu(k,1945)
         lu(k,1949) = lu(k,1949) - lu(k,1489) * lu(k,1945)
         lu(k,1950) = lu(k,1950) - lu(k,1490) * lu(k,1945)
         lu(k,1952) = lu(k,1952) - lu(k,1491) * lu(k,1945)
         lu(k,1953) = lu(k,1953) - lu(k,1492) * lu(k,1945)
         lu(k,1954) = lu(k,1954) - lu(k,1493) * lu(k,1945)
         lu(k,1955) = lu(k,1955) - lu(k,1494) * lu(k,1945)
         lu(k,1956) = lu(k,1956) - lu(k,1495) * lu(k,1945)
         lu(k,1957) = lu(k,1957) - lu(k,1496) * lu(k,1945)
         lu(k,1958) = lu(k,1958) - lu(k,1497) * lu(k,1945)
         lu(k,1959) = lu(k,1959) - lu(k,1498) * lu(k,1945)
         lu(k,1972) = lu(k,1972) - lu(k,1486) * lu(k,1971)
         lu(k,1973) = lu(k,1973) - lu(k,1487) * lu(k,1971)
         lu(k,1974) = lu(k,1974) - lu(k,1488) * lu(k,1971)
         lu(k,1975) = lu(k,1975) - lu(k,1489) * lu(k,1971)
         lu(k,1976) = lu(k,1976) - lu(k,1490) * lu(k,1971)
         lu(k,1978) = lu(k,1978) - lu(k,1491) * lu(k,1971)
         lu(k,1979) = lu(k,1979) - lu(k,1492) * lu(k,1971)
         lu(k,1980) = lu(k,1980) - lu(k,1493) * lu(k,1971)
         lu(k,1981) = lu(k,1981) - lu(k,1494) * lu(k,1971)
         lu(k,1982) = lu(k,1982) - lu(k,1495) * lu(k,1971)
         lu(k,1983) = lu(k,1983) - lu(k,1496) * lu(k,1971)
         lu(k,1984) = lu(k,1984) - lu(k,1497) * lu(k,1971)
         lu(k,1985) = lu(k,1985) - lu(k,1498) * lu(k,1971)
         lu(k,2011) = lu(k,2011) - lu(k,1486) * lu(k,2010)
         lu(k,2012) = lu(k,2012) - lu(k,1487) * lu(k,2010)
         lu(k,2013) = lu(k,2013) - lu(k,1488) * lu(k,2010)
         lu(k,2014) = lu(k,2014) - lu(k,1489) * lu(k,2010)
         lu(k,2015) = lu(k,2015) - lu(k,1490) * lu(k,2010)
         lu(k,2017) = lu(k,2017) - lu(k,1491) * lu(k,2010)
         lu(k,2018) = lu(k,2018) - lu(k,1492) * lu(k,2010)
         lu(k,2019) = lu(k,2019) - lu(k,1493) * lu(k,2010)
         lu(k,2020) = lu(k,2020) - lu(k,1494) * lu(k,2010)
         lu(k,2021) = lu(k,2021) - lu(k,1495) * lu(k,2010)
         lu(k,2022) = lu(k,2022) - lu(k,1496) * lu(k,2010)
         lu(k,2023) = lu(k,2023) - lu(k,1497) * lu(k,2010)
         lu(k,2024) = lu(k,2024) - lu(k,1498) * lu(k,2010)
         lu(k,2063) = lu(k,2063) - lu(k,1486) * lu(k,2062)
         lu(k,2064) = lu(k,2064) - lu(k,1487) * lu(k,2062)
         lu(k,2065) = lu(k,2065) - lu(k,1488) * lu(k,2062)
         lu(k,2066) = lu(k,2066) - lu(k,1489) * lu(k,2062)
         lu(k,2067) = lu(k,2067) - lu(k,1490) * lu(k,2062)
         lu(k,2069) = lu(k,2069) - lu(k,1491) * lu(k,2062)
         lu(k,2070) = lu(k,2070) - lu(k,1492) * lu(k,2062)
         lu(k,2071) = lu(k,2071) - lu(k,1493) * lu(k,2062)
         lu(k,2072) = lu(k,2072) - lu(k,1494) * lu(k,2062)
         lu(k,2073) = lu(k,2073) - lu(k,1495) * lu(k,2062)
         lu(k,2074) = - lu(k,1496) * lu(k,2062)
         lu(k,2075) = lu(k,2075) - lu(k,1497) * lu(k,2062)
         lu(k,2076) = lu(k,2076) - lu(k,1498) * lu(k,2062)
         lu(k,2124) = lu(k,2124) - lu(k,1486) * lu(k,2123)
         lu(k,2125) = lu(k,2125) - lu(k,1487) * lu(k,2123)
         lu(k,2126) = lu(k,2126) - lu(k,1488) * lu(k,2123)
         lu(k,2127) = lu(k,2127) - lu(k,1489) * lu(k,2123)
         lu(k,2128) = lu(k,2128) - lu(k,1490) * lu(k,2123)
         lu(k,2130) = lu(k,2130) - lu(k,1491) * lu(k,2123)
         lu(k,2131) = lu(k,2131) - lu(k,1492) * lu(k,2123)
         lu(k,2132) = lu(k,2132) - lu(k,1493) * lu(k,2123)
         lu(k,2133) = lu(k,2133) - lu(k,1494) * lu(k,2123)
         lu(k,2134) = lu(k,2134) - lu(k,1495) * lu(k,2123)
         lu(k,2135) = lu(k,2135) - lu(k,1496) * lu(k,2123)
         lu(k,2136) = lu(k,2136) - lu(k,1497) * lu(k,2123)
         lu(k,2137) = lu(k,2137) - lu(k,1498) * lu(k,2123)
         lu(k,2147) = lu(k,2147) - lu(k,1486) * lu(k,2146)
         lu(k,2148) = lu(k,2148) - lu(k,1487) * lu(k,2146)
         lu(k,2149) = lu(k,2149) - lu(k,1488) * lu(k,2146)
         lu(k,2150) = - lu(k,1489) * lu(k,2146)
         lu(k,2151) = lu(k,2151) - lu(k,1490) * lu(k,2146)
         lu(k,2153) = lu(k,2153) - lu(k,1491) * lu(k,2146)
         lu(k,2154) = lu(k,2154) - lu(k,1492) * lu(k,2146)
         lu(k,2155) = lu(k,2155) - lu(k,1493) * lu(k,2146)
         lu(k,2156) = lu(k,2156) - lu(k,1494) * lu(k,2146)
         lu(k,2157) = lu(k,2157) - lu(k,1495) * lu(k,2146)
         lu(k,2158) = lu(k,2158) - lu(k,1496) * lu(k,2146)
         lu(k,2159) = lu(k,2159) - lu(k,1497) * lu(k,2146)
         lu(k,2160) = lu(k,2160) - lu(k,1498) * lu(k,2146)
         lu(k,2191) = lu(k,2191) - lu(k,1486) * lu(k,2190)
         lu(k,2192) = lu(k,2192) - lu(k,1487) * lu(k,2190)
         lu(k,2193) = lu(k,2193) - lu(k,1488) * lu(k,2190)
         lu(k,2194) = lu(k,2194) - lu(k,1489) * lu(k,2190)
         lu(k,2195) = lu(k,2195) - lu(k,1490) * lu(k,2190)
         lu(k,2197) = lu(k,2197) - lu(k,1491) * lu(k,2190)
         lu(k,2198) = lu(k,2198) - lu(k,1492) * lu(k,2190)
         lu(k,2199) = lu(k,2199) - lu(k,1493) * lu(k,2190)
         lu(k,2200) = lu(k,2200) - lu(k,1494) * lu(k,2190)
         lu(k,2201) = lu(k,2201) - lu(k,1495) * lu(k,2190)
         lu(k,2202) = lu(k,2202) - lu(k,1496) * lu(k,2190)
         lu(k,2203) = lu(k,2203) - lu(k,1497) * lu(k,2190)
         lu(k,2204) = lu(k,2204) - lu(k,1498) * lu(k,2190)
         lu(k,2215) = lu(k,2215) - lu(k,1486) * lu(k,2214)
         lu(k,2216) = lu(k,2216) - lu(k,1487) * lu(k,2214)
         lu(k,2217) = lu(k,2217) - lu(k,1488) * lu(k,2214)
         lu(k,2218) = lu(k,2218) - lu(k,1489) * lu(k,2214)
         lu(k,2219) = lu(k,2219) - lu(k,1490) * lu(k,2214)
         lu(k,2221) = lu(k,2221) - lu(k,1491) * lu(k,2214)
         lu(k,2222) = lu(k,2222) - lu(k,1492) * lu(k,2214)
         lu(k,2223) = lu(k,2223) - lu(k,1493) * lu(k,2214)
         lu(k,2224) = lu(k,2224) - lu(k,1494) * lu(k,2214)
         lu(k,2225) = lu(k,2225) - lu(k,1495) * lu(k,2214)
         lu(k,2226) = lu(k,2226) - lu(k,1496) * lu(k,2214)
         lu(k,2227) = lu(k,2227) - lu(k,1497) * lu(k,2214)
         lu(k,2228) = lu(k,2228) - lu(k,1498) * lu(k,2214)
         lu(k,2246) = lu(k,2246) - lu(k,1486) * lu(k,2245)
         lu(k,2247) = lu(k,2247) - lu(k,1487) * lu(k,2245)
         lu(k,2248) = lu(k,2248) - lu(k,1488) * lu(k,2245)
         lu(k,2249) = lu(k,2249) - lu(k,1489) * lu(k,2245)
         lu(k,2250) = lu(k,2250) - lu(k,1490) * lu(k,2245)
         lu(k,2252) = lu(k,2252) - lu(k,1491) * lu(k,2245)
         lu(k,2253) = lu(k,2253) - lu(k,1492) * lu(k,2245)
         lu(k,2254) = lu(k,2254) - lu(k,1493) * lu(k,2245)
         lu(k,2255) = lu(k,2255) - lu(k,1494) * lu(k,2245)
         lu(k,2256) = lu(k,2256) - lu(k,1495) * lu(k,2245)
         lu(k,2257) = lu(k,2257) - lu(k,1496) * lu(k,2245)
         lu(k,2258) = lu(k,2258) - lu(k,1497) * lu(k,2245)
         lu(k,2259) = lu(k,2259) - lu(k,1498) * lu(k,2245)
         lu(k,2272) = lu(k,2272) - lu(k,1486) * lu(k,2271)
         lu(k,2273) = lu(k,2273) - lu(k,1487) * lu(k,2271)
         lu(k,2274) = lu(k,2274) - lu(k,1488) * lu(k,2271)
         lu(k,2275) = lu(k,2275) - lu(k,1489) * lu(k,2271)
         lu(k,2276) = lu(k,2276) - lu(k,1490) * lu(k,2271)
         lu(k,2278) = lu(k,2278) - lu(k,1491) * lu(k,2271)
         lu(k,2279) = lu(k,2279) - lu(k,1492) * lu(k,2271)
         lu(k,2280) = lu(k,2280) - lu(k,1493) * lu(k,2271)
         lu(k,2281) = lu(k,2281) - lu(k,1494) * lu(k,2271)
         lu(k,2282) = lu(k,2282) - lu(k,1495) * lu(k,2271)
         lu(k,2283) = lu(k,2283) - lu(k,1496) * lu(k,2271)
         lu(k,2284) = lu(k,2284) - lu(k,1497) * lu(k,2271)
         lu(k,2285) = lu(k,2285) - lu(k,1498) * lu(k,2271)
                                                                        
         lu(k,1526) = 1._r8 / lu(k,1526)
         lu(k,1527) = lu(k,1527) * lu(k,1526)
         lu(k,1528) = lu(k,1528) * lu(k,1526)
         lu(k,1529) = lu(k,1529) * lu(k,1526)
         lu(k,1530) = lu(k,1530) * lu(k,1526)
         lu(k,1531) = lu(k,1531) * lu(k,1526)
         lu(k,1532) = lu(k,1532) * lu(k,1526)
         lu(k,1533) = lu(k,1533) * lu(k,1526)
         lu(k,1534) = lu(k,1534) * lu(k,1526)
         lu(k,1535) = lu(k,1535) * lu(k,1526)
         lu(k,1536) = lu(k,1536) * lu(k,1526)
         lu(k,1537) = lu(k,1537) * lu(k,1526)
         lu(k,1538) = lu(k,1538) * lu(k,1526)
         lu(k,1539) = lu(k,1539) * lu(k,1526)
         lu(k,1691) = lu(k,1691) - lu(k,1527) * lu(k,1690)
         lu(k,1692) = lu(k,1692) - lu(k,1528) * lu(k,1690)
         lu(k,1693) = lu(k,1693) - lu(k,1529) * lu(k,1690)
         lu(k,1694) = lu(k,1694) - lu(k,1530) * lu(k,1690)
         lu(k,1695) = lu(k,1695) - lu(k,1531) * lu(k,1690)
         lu(k,1696) = lu(k,1696) - lu(k,1532) * lu(k,1690)
         lu(k,1697) = lu(k,1697) - lu(k,1533) * lu(k,1690)
         lu(k,1698) = lu(k,1698) - lu(k,1534) * lu(k,1690)
         lu(k,1699) = lu(k,1699) - lu(k,1535) * lu(k,1690)
         lu(k,1700) = lu(k,1700) - lu(k,1536) * lu(k,1690)
         lu(k,1701) = lu(k,1701) - lu(k,1537) * lu(k,1690)
         lu(k,1702) = lu(k,1702) - lu(k,1538) * lu(k,1690)
         lu(k,1703) = lu(k,1703) - lu(k,1539) * lu(k,1690)
         lu(k,1748) = lu(k,1748) - lu(k,1527) * lu(k,1747)
         lu(k,1749) = lu(k,1749) - lu(k,1528) * lu(k,1747)
         lu(k,1750) = lu(k,1750) - lu(k,1529) * lu(k,1747)
         lu(k,1751) = lu(k,1751) - lu(k,1530) * lu(k,1747)
         lu(k,1752) = lu(k,1752) - lu(k,1531) * lu(k,1747)
         lu(k,1753) = lu(k,1753) - lu(k,1532) * lu(k,1747)
         lu(k,1754) = lu(k,1754) - lu(k,1533) * lu(k,1747)
         lu(k,1755) = lu(k,1755) - lu(k,1534) * lu(k,1747)
         lu(k,1756) = lu(k,1756) - lu(k,1535) * lu(k,1747)
         lu(k,1757) = lu(k,1757) - lu(k,1536) * lu(k,1747)
         lu(k,1758) = lu(k,1758) - lu(k,1537) * lu(k,1747)
         lu(k,1759) = lu(k,1759) - lu(k,1538) * lu(k,1747)
         lu(k,1760) = lu(k,1760) - lu(k,1539) * lu(k,1747)
         lu(k,1840) = lu(k,1840) - lu(k,1527) * lu(k,1839)
         lu(k,1841) = lu(k,1841) - lu(k,1528) * lu(k,1839)
         lu(k,1842) = lu(k,1842) - lu(k,1529) * lu(k,1839)
         lu(k,1843) = lu(k,1843) - lu(k,1530) * lu(k,1839)
         lu(k,1844) = lu(k,1844) - lu(k,1531) * lu(k,1839)
         lu(k,1845) = lu(k,1845) - lu(k,1532) * lu(k,1839)
         lu(k,1846) = lu(k,1846) - lu(k,1533) * lu(k,1839)
         lu(k,1847) = lu(k,1847) - lu(k,1534) * lu(k,1839)
         lu(k,1848) = lu(k,1848) - lu(k,1535) * lu(k,1839)
         lu(k,1849) = lu(k,1849) - lu(k,1536) * lu(k,1839)
         lu(k,1850) = lu(k,1850) - lu(k,1537) * lu(k,1839)
         lu(k,1851) = lu(k,1851) - lu(k,1538) * lu(k,1839)
         lu(k,1852) = lu(k,1852) - lu(k,1539) * lu(k,1839)
         lu(k,1947) = lu(k,1947) - lu(k,1527) * lu(k,1946)
         lu(k,1948) = lu(k,1948) - lu(k,1528) * lu(k,1946)
         lu(k,1949) = lu(k,1949) - lu(k,1529) * lu(k,1946)
         lu(k,1950) = lu(k,1950) - lu(k,1530) * lu(k,1946)
         lu(k,1951) = lu(k,1951) - lu(k,1531) * lu(k,1946)
         lu(k,1952) = lu(k,1952) - lu(k,1532) * lu(k,1946)
         lu(k,1953) = lu(k,1953) - lu(k,1533) * lu(k,1946)
         lu(k,1954) = lu(k,1954) - lu(k,1534) * lu(k,1946)
         lu(k,1955) = lu(k,1955) - lu(k,1535) * lu(k,1946)
         lu(k,1956) = lu(k,1956) - lu(k,1536) * lu(k,1946)
         lu(k,1957) = lu(k,1957) - lu(k,1537) * lu(k,1946)
         lu(k,1958) = lu(k,1958) - lu(k,1538) * lu(k,1946)
         lu(k,1959) = lu(k,1959) - lu(k,1539) * lu(k,1946)
         lu(k,1973) = lu(k,1973) - lu(k,1527) * lu(k,1972)
         lu(k,1974) = lu(k,1974) - lu(k,1528) * lu(k,1972)
         lu(k,1975) = lu(k,1975) - lu(k,1529) * lu(k,1972)
         lu(k,1976) = lu(k,1976) - lu(k,1530) * lu(k,1972)
         lu(k,1977) = lu(k,1977) - lu(k,1531) * lu(k,1972)
         lu(k,1978) = lu(k,1978) - lu(k,1532) * lu(k,1972)
         lu(k,1979) = lu(k,1979) - lu(k,1533) * lu(k,1972)
         lu(k,1980) = lu(k,1980) - lu(k,1534) * lu(k,1972)
         lu(k,1981) = lu(k,1981) - lu(k,1535) * lu(k,1972)
         lu(k,1982) = lu(k,1982) - lu(k,1536) * lu(k,1972)
         lu(k,1983) = lu(k,1983) - lu(k,1537) * lu(k,1972)
         lu(k,1984) = lu(k,1984) - lu(k,1538) * lu(k,1972)
         lu(k,1985) = lu(k,1985) - lu(k,1539) * lu(k,1972)
         lu(k,2012) = lu(k,2012) - lu(k,1527) * lu(k,2011)
         lu(k,2013) = lu(k,2013) - lu(k,1528) * lu(k,2011)
         lu(k,2014) = lu(k,2014) - lu(k,1529) * lu(k,2011)
         lu(k,2015) = lu(k,2015) - lu(k,1530) * lu(k,2011)
         lu(k,2016) = lu(k,2016) - lu(k,1531) * lu(k,2011)
         lu(k,2017) = lu(k,2017) - lu(k,1532) * lu(k,2011)
         lu(k,2018) = lu(k,2018) - lu(k,1533) * lu(k,2011)
         lu(k,2019) = lu(k,2019) - lu(k,1534) * lu(k,2011)
         lu(k,2020) = lu(k,2020) - lu(k,1535) * lu(k,2011)
         lu(k,2021) = lu(k,2021) - lu(k,1536) * lu(k,2011)
         lu(k,2022) = lu(k,2022) - lu(k,1537) * lu(k,2011)
         lu(k,2023) = lu(k,2023) - lu(k,1538) * lu(k,2011)
         lu(k,2024) = lu(k,2024) - lu(k,1539) * lu(k,2011)
         lu(k,2064) = lu(k,2064) - lu(k,1527) * lu(k,2063)
         lu(k,2065) = lu(k,2065) - lu(k,1528) * lu(k,2063)
         lu(k,2066) = lu(k,2066) - lu(k,1529) * lu(k,2063)
         lu(k,2067) = lu(k,2067) - lu(k,1530) * lu(k,2063)
         lu(k,2068) = lu(k,2068) - lu(k,1531) * lu(k,2063)
         lu(k,2069) = lu(k,2069) - lu(k,1532) * lu(k,2063)
         lu(k,2070) = lu(k,2070) - lu(k,1533) * lu(k,2063)
         lu(k,2071) = lu(k,2071) - lu(k,1534) * lu(k,2063)
         lu(k,2072) = lu(k,2072) - lu(k,1535) * lu(k,2063)
         lu(k,2073) = lu(k,2073) - lu(k,1536) * lu(k,2063)
         lu(k,2074) = lu(k,2074) - lu(k,1537) * lu(k,2063)
         lu(k,2075) = lu(k,2075) - lu(k,1538) * lu(k,2063)
         lu(k,2076) = lu(k,2076) - lu(k,1539) * lu(k,2063)
         lu(k,2125) = lu(k,2125) - lu(k,1527) * lu(k,2124)
         lu(k,2126) = lu(k,2126) - lu(k,1528) * lu(k,2124)
         lu(k,2127) = lu(k,2127) - lu(k,1529) * lu(k,2124)
         lu(k,2128) = lu(k,2128) - lu(k,1530) * lu(k,2124)
         lu(k,2129) = lu(k,2129) - lu(k,1531) * lu(k,2124)
         lu(k,2130) = lu(k,2130) - lu(k,1532) * lu(k,2124)
         lu(k,2131) = lu(k,2131) - lu(k,1533) * lu(k,2124)
         lu(k,2132) = lu(k,2132) - lu(k,1534) * lu(k,2124)
         lu(k,2133) = lu(k,2133) - lu(k,1535) * lu(k,2124)
         lu(k,2134) = lu(k,2134) - lu(k,1536) * lu(k,2124)
         lu(k,2135) = lu(k,2135) - lu(k,1537) * lu(k,2124)
         lu(k,2136) = lu(k,2136) - lu(k,1538) * lu(k,2124)
         lu(k,2137) = lu(k,2137) - lu(k,1539) * lu(k,2124)
         lu(k,2148) = lu(k,2148) - lu(k,1527) * lu(k,2147)
         lu(k,2149) = lu(k,2149) - lu(k,1528) * lu(k,2147)
         lu(k,2150) = lu(k,2150) - lu(k,1529) * lu(k,2147)
         lu(k,2151) = lu(k,2151) - lu(k,1530) * lu(k,2147)
         lu(k,2152) = lu(k,2152) - lu(k,1531) * lu(k,2147)
         lu(k,2153) = lu(k,2153) - lu(k,1532) * lu(k,2147)
         lu(k,2154) = lu(k,2154) - lu(k,1533) * lu(k,2147)
         lu(k,2155) = lu(k,2155) - lu(k,1534) * lu(k,2147)
         lu(k,2156) = lu(k,2156) - lu(k,1535) * lu(k,2147)
         lu(k,2157) = lu(k,2157) - lu(k,1536) * lu(k,2147)
         lu(k,2158) = lu(k,2158) - lu(k,1537) * lu(k,2147)
         lu(k,2159) = lu(k,2159) - lu(k,1538) * lu(k,2147)
         lu(k,2160) = lu(k,2160) - lu(k,1539) * lu(k,2147)
         lu(k,2192) = lu(k,2192) - lu(k,1527) * lu(k,2191)
         lu(k,2193) = lu(k,2193) - lu(k,1528) * lu(k,2191)
         lu(k,2194) = lu(k,2194) - lu(k,1529) * lu(k,2191)
         lu(k,2195) = lu(k,2195) - lu(k,1530) * lu(k,2191)
         lu(k,2196) = lu(k,2196) - lu(k,1531) * lu(k,2191)
         lu(k,2197) = lu(k,2197) - lu(k,1532) * lu(k,2191)
         lu(k,2198) = lu(k,2198) - lu(k,1533) * lu(k,2191)
         lu(k,2199) = lu(k,2199) - lu(k,1534) * lu(k,2191)
         lu(k,2200) = lu(k,2200) - lu(k,1535) * lu(k,2191)
         lu(k,2201) = lu(k,2201) - lu(k,1536) * lu(k,2191)
         lu(k,2202) = lu(k,2202) - lu(k,1537) * lu(k,2191)
         lu(k,2203) = lu(k,2203) - lu(k,1538) * lu(k,2191)
         lu(k,2204) = lu(k,2204) - lu(k,1539) * lu(k,2191)
         lu(k,2216) = lu(k,2216) - lu(k,1527) * lu(k,2215)
         lu(k,2217) = lu(k,2217) - lu(k,1528) * lu(k,2215)
         lu(k,2218) = lu(k,2218) - lu(k,1529) * lu(k,2215)
         lu(k,2219) = lu(k,2219) - lu(k,1530) * lu(k,2215)
         lu(k,2220) = lu(k,2220) - lu(k,1531) * lu(k,2215)
         lu(k,2221) = lu(k,2221) - lu(k,1532) * lu(k,2215)
         lu(k,2222) = lu(k,2222) - lu(k,1533) * lu(k,2215)
         lu(k,2223) = lu(k,2223) - lu(k,1534) * lu(k,2215)
         lu(k,2224) = lu(k,2224) - lu(k,1535) * lu(k,2215)
         lu(k,2225) = lu(k,2225) - lu(k,1536) * lu(k,2215)
         lu(k,2226) = lu(k,2226) - lu(k,1537) * lu(k,2215)
         lu(k,2227) = lu(k,2227) - lu(k,1538) * lu(k,2215)
         lu(k,2228) = lu(k,2228) - lu(k,1539) * lu(k,2215)
         lu(k,2247) = lu(k,2247) - lu(k,1527) * lu(k,2246)
         lu(k,2248) = lu(k,2248) - lu(k,1528) * lu(k,2246)
         lu(k,2249) = lu(k,2249) - lu(k,1529) * lu(k,2246)
         lu(k,2250) = lu(k,2250) - lu(k,1530) * lu(k,2246)
         lu(k,2251) = lu(k,2251) - lu(k,1531) * lu(k,2246)
         lu(k,2252) = lu(k,2252) - lu(k,1532) * lu(k,2246)
         lu(k,2253) = lu(k,2253) - lu(k,1533) * lu(k,2246)
         lu(k,2254) = lu(k,2254) - lu(k,1534) * lu(k,2246)
         lu(k,2255) = lu(k,2255) - lu(k,1535) * lu(k,2246)
         lu(k,2256) = lu(k,2256) - lu(k,1536) * lu(k,2246)
         lu(k,2257) = lu(k,2257) - lu(k,1537) * lu(k,2246)
         lu(k,2258) = lu(k,2258) - lu(k,1538) * lu(k,2246)
         lu(k,2259) = lu(k,2259) - lu(k,1539) * lu(k,2246)
         lu(k,2273) = lu(k,2273) - lu(k,1527) * lu(k,2272)
         lu(k,2274) = lu(k,2274) - lu(k,1528) * lu(k,2272)
         lu(k,2275) = lu(k,2275) - lu(k,1529) * lu(k,2272)
         lu(k,2276) = lu(k,2276) - lu(k,1530) * lu(k,2272)
         lu(k,2277) = lu(k,2277) - lu(k,1531) * lu(k,2272)
         lu(k,2278) = lu(k,2278) - lu(k,1532) * lu(k,2272)
         lu(k,2279) = lu(k,2279) - lu(k,1533) * lu(k,2272)
         lu(k,2280) = lu(k,2280) - lu(k,1534) * lu(k,2272)
         lu(k,2281) = lu(k,2281) - lu(k,1535) * lu(k,2272)
         lu(k,2282) = lu(k,2282) - lu(k,1536) * lu(k,2272)
         lu(k,2283) = lu(k,2283) - lu(k,1537) * lu(k,2272)
         lu(k,2284) = lu(k,2284) - lu(k,1538) * lu(k,2272)
         lu(k,2285) = lu(k,2285) - lu(k,1539) * lu(k,2272)
                                                                        
         lu(k,1691) = 1._r8 / lu(k,1691)
         lu(k,1692) = lu(k,1692) * lu(k,1691)
         lu(k,1693) = lu(k,1693) * lu(k,1691)
         lu(k,1694) = lu(k,1694) * lu(k,1691)
         lu(k,1695) = lu(k,1695) * lu(k,1691)
         lu(k,1696) = lu(k,1696) * lu(k,1691)
         lu(k,1697) = lu(k,1697) * lu(k,1691)
         lu(k,1698) = lu(k,1698) * lu(k,1691)
         lu(k,1699) = lu(k,1699) * lu(k,1691)
         lu(k,1700) = lu(k,1700) * lu(k,1691)
         lu(k,1701) = lu(k,1701) * lu(k,1691)
         lu(k,1702) = lu(k,1702) * lu(k,1691)
         lu(k,1703) = lu(k,1703) * lu(k,1691)
         lu(k,1749) = lu(k,1749) - lu(k,1692) * lu(k,1748)
         lu(k,1750) = lu(k,1750) - lu(k,1693) * lu(k,1748)
         lu(k,1751) = lu(k,1751) - lu(k,1694) * lu(k,1748)
         lu(k,1752) = lu(k,1752) - lu(k,1695) * lu(k,1748)
         lu(k,1753) = lu(k,1753) - lu(k,1696) * lu(k,1748)
         lu(k,1754) = lu(k,1754) - lu(k,1697) * lu(k,1748)
         lu(k,1755) = lu(k,1755) - lu(k,1698) * lu(k,1748)
         lu(k,1756) = lu(k,1756) - lu(k,1699) * lu(k,1748)
         lu(k,1757) = lu(k,1757) - lu(k,1700) * lu(k,1748)
         lu(k,1758) = lu(k,1758) - lu(k,1701) * lu(k,1748)
         lu(k,1759) = lu(k,1759) - lu(k,1702) * lu(k,1748)
         lu(k,1760) = lu(k,1760) - lu(k,1703) * lu(k,1748)
         lu(k,1841) = lu(k,1841) - lu(k,1692) * lu(k,1840)
         lu(k,1842) = lu(k,1842) - lu(k,1693) * lu(k,1840)
         lu(k,1843) = lu(k,1843) - lu(k,1694) * lu(k,1840)
         lu(k,1844) = lu(k,1844) - lu(k,1695) * lu(k,1840)
         lu(k,1845) = lu(k,1845) - lu(k,1696) * lu(k,1840)
         lu(k,1846) = lu(k,1846) - lu(k,1697) * lu(k,1840)
         lu(k,1847) = lu(k,1847) - lu(k,1698) * lu(k,1840)
         lu(k,1848) = lu(k,1848) - lu(k,1699) * lu(k,1840)
         lu(k,1849) = lu(k,1849) - lu(k,1700) * lu(k,1840)
         lu(k,1850) = lu(k,1850) - lu(k,1701) * lu(k,1840)
         lu(k,1851) = lu(k,1851) - lu(k,1702) * lu(k,1840)
         lu(k,1852) = lu(k,1852) - lu(k,1703) * lu(k,1840)
         lu(k,1948) = lu(k,1948) - lu(k,1692) * lu(k,1947)
         lu(k,1949) = lu(k,1949) - lu(k,1693) * lu(k,1947)
         lu(k,1950) = lu(k,1950) - lu(k,1694) * lu(k,1947)
         lu(k,1951) = lu(k,1951) - lu(k,1695) * lu(k,1947)
         lu(k,1952) = lu(k,1952) - lu(k,1696) * lu(k,1947)
         lu(k,1953) = lu(k,1953) - lu(k,1697) * lu(k,1947)
         lu(k,1954) = lu(k,1954) - lu(k,1698) * lu(k,1947)
         lu(k,1955) = lu(k,1955) - lu(k,1699) * lu(k,1947)
         lu(k,1956) = lu(k,1956) - lu(k,1700) * lu(k,1947)
         lu(k,1957) = lu(k,1957) - lu(k,1701) * lu(k,1947)
         lu(k,1958) = lu(k,1958) - lu(k,1702) * lu(k,1947)
         lu(k,1959) = lu(k,1959) - lu(k,1703) * lu(k,1947)
         lu(k,1974) = lu(k,1974) - lu(k,1692) * lu(k,1973)
         lu(k,1975) = lu(k,1975) - lu(k,1693) * lu(k,1973)
         lu(k,1976) = lu(k,1976) - lu(k,1694) * lu(k,1973)
         lu(k,1977) = lu(k,1977) - lu(k,1695) * lu(k,1973)
         lu(k,1978) = lu(k,1978) - lu(k,1696) * lu(k,1973)
         lu(k,1979) = lu(k,1979) - lu(k,1697) * lu(k,1973)
         lu(k,1980) = lu(k,1980) - lu(k,1698) * lu(k,1973)
         lu(k,1981) = lu(k,1981) - lu(k,1699) * lu(k,1973)
         lu(k,1982) = lu(k,1982) - lu(k,1700) * lu(k,1973)
         lu(k,1983) = lu(k,1983) - lu(k,1701) * lu(k,1973)
         lu(k,1984) = lu(k,1984) - lu(k,1702) * lu(k,1973)
         lu(k,1985) = lu(k,1985) - lu(k,1703) * lu(k,1973)
         lu(k,2013) = lu(k,2013) - lu(k,1692) * lu(k,2012)
         lu(k,2014) = lu(k,2014) - lu(k,1693) * lu(k,2012)
         lu(k,2015) = lu(k,2015) - lu(k,1694) * lu(k,2012)
         lu(k,2016) = lu(k,2016) - lu(k,1695) * lu(k,2012)
         lu(k,2017) = lu(k,2017) - lu(k,1696) * lu(k,2012)
         lu(k,2018) = lu(k,2018) - lu(k,1697) * lu(k,2012)
         lu(k,2019) = lu(k,2019) - lu(k,1698) * lu(k,2012)
         lu(k,2020) = lu(k,2020) - lu(k,1699) * lu(k,2012)
         lu(k,2021) = lu(k,2021) - lu(k,1700) * lu(k,2012)
         lu(k,2022) = lu(k,2022) - lu(k,1701) * lu(k,2012)
         lu(k,2023) = lu(k,2023) - lu(k,1702) * lu(k,2012)
         lu(k,2024) = lu(k,2024) - lu(k,1703) * lu(k,2012)
         lu(k,2065) = lu(k,2065) - lu(k,1692) * lu(k,2064)
         lu(k,2066) = lu(k,2066) - lu(k,1693) * lu(k,2064)
         lu(k,2067) = lu(k,2067) - lu(k,1694) * lu(k,2064)
         lu(k,2068) = lu(k,2068) - lu(k,1695) * lu(k,2064)
         lu(k,2069) = lu(k,2069) - lu(k,1696) * lu(k,2064)
         lu(k,2070) = lu(k,2070) - lu(k,1697) * lu(k,2064)
         lu(k,2071) = lu(k,2071) - lu(k,1698) * lu(k,2064)
         lu(k,2072) = lu(k,2072) - lu(k,1699) * lu(k,2064)
         lu(k,2073) = lu(k,2073) - lu(k,1700) * lu(k,2064)
         lu(k,2074) = lu(k,2074) - lu(k,1701) * lu(k,2064)
         lu(k,2075) = lu(k,2075) - lu(k,1702) * lu(k,2064)
         lu(k,2076) = lu(k,2076) - lu(k,1703) * lu(k,2064)
         lu(k,2126) = lu(k,2126) - lu(k,1692) * lu(k,2125)
         lu(k,2127) = lu(k,2127) - lu(k,1693) * lu(k,2125)
         lu(k,2128) = lu(k,2128) - lu(k,1694) * lu(k,2125)
         lu(k,2129) = lu(k,2129) - lu(k,1695) * lu(k,2125)
         lu(k,2130) = lu(k,2130) - lu(k,1696) * lu(k,2125)
         lu(k,2131) = lu(k,2131) - lu(k,1697) * lu(k,2125)
         lu(k,2132) = lu(k,2132) - lu(k,1698) * lu(k,2125)
         lu(k,2133) = lu(k,2133) - lu(k,1699) * lu(k,2125)
         lu(k,2134) = lu(k,2134) - lu(k,1700) * lu(k,2125)
         lu(k,2135) = lu(k,2135) - lu(k,1701) * lu(k,2125)
         lu(k,2136) = lu(k,2136) - lu(k,1702) * lu(k,2125)
         lu(k,2137) = lu(k,2137) - lu(k,1703) * lu(k,2125)
         lu(k,2149) = lu(k,2149) - lu(k,1692) * lu(k,2148)
         lu(k,2150) = lu(k,2150) - lu(k,1693) * lu(k,2148)
         lu(k,2151) = lu(k,2151) - lu(k,1694) * lu(k,2148)
         lu(k,2152) = lu(k,2152) - lu(k,1695) * lu(k,2148)
         lu(k,2153) = lu(k,2153) - lu(k,1696) * lu(k,2148)
         lu(k,2154) = lu(k,2154) - lu(k,1697) * lu(k,2148)
         lu(k,2155) = lu(k,2155) - lu(k,1698) * lu(k,2148)
         lu(k,2156) = lu(k,2156) - lu(k,1699) * lu(k,2148)
         lu(k,2157) = lu(k,2157) - lu(k,1700) * lu(k,2148)
         lu(k,2158) = lu(k,2158) - lu(k,1701) * lu(k,2148)
         lu(k,2159) = lu(k,2159) - lu(k,1702) * lu(k,2148)
         lu(k,2160) = lu(k,2160) - lu(k,1703) * lu(k,2148)
         lu(k,2193) = lu(k,2193) - lu(k,1692) * lu(k,2192)
         lu(k,2194) = lu(k,2194) - lu(k,1693) * lu(k,2192)
         lu(k,2195) = lu(k,2195) - lu(k,1694) * lu(k,2192)
         lu(k,2196) = lu(k,2196) - lu(k,1695) * lu(k,2192)
         lu(k,2197) = lu(k,2197) - lu(k,1696) * lu(k,2192)
         lu(k,2198) = lu(k,2198) - lu(k,1697) * lu(k,2192)
         lu(k,2199) = lu(k,2199) - lu(k,1698) * lu(k,2192)
         lu(k,2200) = lu(k,2200) - lu(k,1699) * lu(k,2192)
         lu(k,2201) = lu(k,2201) - lu(k,1700) * lu(k,2192)
         lu(k,2202) = lu(k,2202) - lu(k,1701) * lu(k,2192)
         lu(k,2203) = lu(k,2203) - lu(k,1702) * lu(k,2192)
         lu(k,2204) = lu(k,2204) - lu(k,1703) * lu(k,2192)
         lu(k,2217) = lu(k,2217) - lu(k,1692) * lu(k,2216)
         lu(k,2218) = lu(k,2218) - lu(k,1693) * lu(k,2216)
         lu(k,2219) = lu(k,2219) - lu(k,1694) * lu(k,2216)
         lu(k,2220) = lu(k,2220) - lu(k,1695) * lu(k,2216)
         lu(k,2221) = lu(k,2221) - lu(k,1696) * lu(k,2216)
         lu(k,2222) = lu(k,2222) - lu(k,1697) * lu(k,2216)
         lu(k,2223) = lu(k,2223) - lu(k,1698) * lu(k,2216)
         lu(k,2224) = lu(k,2224) - lu(k,1699) * lu(k,2216)
         lu(k,2225) = lu(k,2225) - lu(k,1700) * lu(k,2216)
         lu(k,2226) = lu(k,2226) - lu(k,1701) * lu(k,2216)
         lu(k,2227) = lu(k,2227) - lu(k,1702) * lu(k,2216)
         lu(k,2228) = lu(k,2228) - lu(k,1703) * lu(k,2216)
         lu(k,2248) = lu(k,2248) - lu(k,1692) * lu(k,2247)
         lu(k,2249) = lu(k,2249) - lu(k,1693) * lu(k,2247)
         lu(k,2250) = lu(k,2250) - lu(k,1694) * lu(k,2247)
         lu(k,2251) = lu(k,2251) - lu(k,1695) * lu(k,2247)
         lu(k,2252) = lu(k,2252) - lu(k,1696) * lu(k,2247)
         lu(k,2253) = lu(k,2253) - lu(k,1697) * lu(k,2247)
         lu(k,2254) = lu(k,2254) - lu(k,1698) * lu(k,2247)
         lu(k,2255) = lu(k,2255) - lu(k,1699) * lu(k,2247)
         lu(k,2256) = lu(k,2256) - lu(k,1700) * lu(k,2247)
         lu(k,2257) = lu(k,2257) - lu(k,1701) * lu(k,2247)
         lu(k,2258) = lu(k,2258) - lu(k,1702) * lu(k,2247)
         lu(k,2259) = lu(k,2259) - lu(k,1703) * lu(k,2247)
         lu(k,2274) = lu(k,2274) - lu(k,1692) * lu(k,2273)
         lu(k,2275) = lu(k,2275) - lu(k,1693) * lu(k,2273)
         lu(k,2276) = lu(k,2276) - lu(k,1694) * lu(k,2273)
         lu(k,2277) = lu(k,2277) - lu(k,1695) * lu(k,2273)
         lu(k,2278) = lu(k,2278) - lu(k,1696) * lu(k,2273)
         lu(k,2279) = lu(k,2279) - lu(k,1697) * lu(k,2273)
         lu(k,2280) = lu(k,2280) - lu(k,1698) * lu(k,2273)
         lu(k,2281) = lu(k,2281) - lu(k,1699) * lu(k,2273)
         lu(k,2282) = lu(k,2282) - lu(k,1700) * lu(k,2273)
         lu(k,2283) = lu(k,2283) - lu(k,1701) * lu(k,2273)
         lu(k,2284) = lu(k,2284) - lu(k,1702) * lu(k,2273)
         lu(k,2285) = lu(k,2285) - lu(k,1703) * lu(k,2273)
                                                                        
         lu(k,1749) = 1._r8 / lu(k,1749)
         lu(k,1750) = lu(k,1750) * lu(k,1749)
         lu(k,1751) = lu(k,1751) * lu(k,1749)
         lu(k,1752) = lu(k,1752) * lu(k,1749)
         lu(k,1753) = lu(k,1753) * lu(k,1749)
         lu(k,1754) = lu(k,1754) * lu(k,1749)
         lu(k,1755) = lu(k,1755) * lu(k,1749)
         lu(k,1756) = lu(k,1756) * lu(k,1749)
         lu(k,1757) = lu(k,1757) * lu(k,1749)
         lu(k,1758) = lu(k,1758) * lu(k,1749)
         lu(k,1759) = lu(k,1759) * lu(k,1749)
         lu(k,1760) = lu(k,1760) * lu(k,1749)
         lu(k,1842) = lu(k,1842) - lu(k,1750) * lu(k,1841)
         lu(k,1843) = lu(k,1843) - lu(k,1751) * lu(k,1841)
         lu(k,1844) = lu(k,1844) - lu(k,1752) * lu(k,1841)
         lu(k,1845) = lu(k,1845) - lu(k,1753) * lu(k,1841)
         lu(k,1846) = lu(k,1846) - lu(k,1754) * lu(k,1841)
         lu(k,1847) = lu(k,1847) - lu(k,1755) * lu(k,1841)
         lu(k,1848) = lu(k,1848) - lu(k,1756) * lu(k,1841)
         lu(k,1849) = lu(k,1849) - lu(k,1757) * lu(k,1841)
         lu(k,1850) = lu(k,1850) - lu(k,1758) * lu(k,1841)
         lu(k,1851) = lu(k,1851) - lu(k,1759) * lu(k,1841)
         lu(k,1852) = lu(k,1852) - lu(k,1760) * lu(k,1841)
         lu(k,1949) = lu(k,1949) - lu(k,1750) * lu(k,1948)
         lu(k,1950) = lu(k,1950) - lu(k,1751) * lu(k,1948)
         lu(k,1951) = lu(k,1951) - lu(k,1752) * lu(k,1948)
         lu(k,1952) = lu(k,1952) - lu(k,1753) * lu(k,1948)
         lu(k,1953) = lu(k,1953) - lu(k,1754) * lu(k,1948)
         lu(k,1954) = lu(k,1954) - lu(k,1755) * lu(k,1948)
         lu(k,1955) = lu(k,1955) - lu(k,1756) * lu(k,1948)
         lu(k,1956) = lu(k,1956) - lu(k,1757) * lu(k,1948)
         lu(k,1957) = lu(k,1957) - lu(k,1758) * lu(k,1948)
         lu(k,1958) = lu(k,1958) - lu(k,1759) * lu(k,1948)
         lu(k,1959) = lu(k,1959) - lu(k,1760) * lu(k,1948)
         lu(k,1975) = lu(k,1975) - lu(k,1750) * lu(k,1974)
         lu(k,1976) = lu(k,1976) - lu(k,1751) * lu(k,1974)
         lu(k,1977) = lu(k,1977) - lu(k,1752) * lu(k,1974)
         lu(k,1978) = lu(k,1978) - lu(k,1753) * lu(k,1974)
         lu(k,1979) = lu(k,1979) - lu(k,1754) * lu(k,1974)
         lu(k,1980) = lu(k,1980) - lu(k,1755) * lu(k,1974)
         lu(k,1981) = lu(k,1981) - lu(k,1756) * lu(k,1974)
         lu(k,1982) = lu(k,1982) - lu(k,1757) * lu(k,1974)
         lu(k,1983) = lu(k,1983) - lu(k,1758) * lu(k,1974)
         lu(k,1984) = lu(k,1984) - lu(k,1759) * lu(k,1974)
         lu(k,1985) = lu(k,1985) - lu(k,1760) * lu(k,1974)
         lu(k,2014) = lu(k,2014) - lu(k,1750) * lu(k,2013)
         lu(k,2015) = lu(k,2015) - lu(k,1751) * lu(k,2013)
         lu(k,2016) = lu(k,2016) - lu(k,1752) * lu(k,2013)
         lu(k,2017) = lu(k,2017) - lu(k,1753) * lu(k,2013)
         lu(k,2018) = lu(k,2018) - lu(k,1754) * lu(k,2013)
         lu(k,2019) = lu(k,2019) - lu(k,1755) * lu(k,2013)
         lu(k,2020) = lu(k,2020) - lu(k,1756) * lu(k,2013)
         lu(k,2021) = lu(k,2021) - lu(k,1757) * lu(k,2013)
         lu(k,2022) = lu(k,2022) - lu(k,1758) * lu(k,2013)
         lu(k,2023) = lu(k,2023) - lu(k,1759) * lu(k,2013)
         lu(k,2024) = lu(k,2024) - lu(k,1760) * lu(k,2013)
         lu(k,2066) = lu(k,2066) - lu(k,1750) * lu(k,2065)
         lu(k,2067) = lu(k,2067) - lu(k,1751) * lu(k,2065)
         lu(k,2068) = lu(k,2068) - lu(k,1752) * lu(k,2065)
         lu(k,2069) = lu(k,2069) - lu(k,1753) * lu(k,2065)
         lu(k,2070) = lu(k,2070) - lu(k,1754) * lu(k,2065)
         lu(k,2071) = lu(k,2071) - lu(k,1755) * lu(k,2065)
         lu(k,2072) = lu(k,2072) - lu(k,1756) * lu(k,2065)
         lu(k,2073) = lu(k,2073) - lu(k,1757) * lu(k,2065)
         lu(k,2074) = lu(k,2074) - lu(k,1758) * lu(k,2065)
         lu(k,2075) = lu(k,2075) - lu(k,1759) * lu(k,2065)
         lu(k,2076) = lu(k,2076) - lu(k,1760) * lu(k,2065)
         lu(k,2127) = lu(k,2127) - lu(k,1750) * lu(k,2126)
         lu(k,2128) = lu(k,2128) - lu(k,1751) * lu(k,2126)
         lu(k,2129) = lu(k,2129) - lu(k,1752) * lu(k,2126)
         lu(k,2130) = lu(k,2130) - lu(k,1753) * lu(k,2126)
         lu(k,2131) = lu(k,2131) - lu(k,1754) * lu(k,2126)
         lu(k,2132) = lu(k,2132) - lu(k,1755) * lu(k,2126)
         lu(k,2133) = lu(k,2133) - lu(k,1756) * lu(k,2126)
         lu(k,2134) = lu(k,2134) - lu(k,1757) * lu(k,2126)
         lu(k,2135) = lu(k,2135) - lu(k,1758) * lu(k,2126)
         lu(k,2136) = lu(k,2136) - lu(k,1759) * lu(k,2126)
         lu(k,2137) = lu(k,2137) - lu(k,1760) * lu(k,2126)
         lu(k,2150) = lu(k,2150) - lu(k,1750) * lu(k,2149)
         lu(k,2151) = lu(k,2151) - lu(k,1751) * lu(k,2149)
         lu(k,2152) = lu(k,2152) - lu(k,1752) * lu(k,2149)
         lu(k,2153) = lu(k,2153) - lu(k,1753) * lu(k,2149)
         lu(k,2154) = lu(k,2154) - lu(k,1754) * lu(k,2149)
         lu(k,2155) = lu(k,2155) - lu(k,1755) * lu(k,2149)
         lu(k,2156) = lu(k,2156) - lu(k,1756) * lu(k,2149)
         lu(k,2157) = lu(k,2157) - lu(k,1757) * lu(k,2149)
         lu(k,2158) = lu(k,2158) - lu(k,1758) * lu(k,2149)
         lu(k,2159) = lu(k,2159) - lu(k,1759) * lu(k,2149)
         lu(k,2160) = lu(k,2160) - lu(k,1760) * lu(k,2149)
         lu(k,2194) = lu(k,2194) - lu(k,1750) * lu(k,2193)
         lu(k,2195) = lu(k,2195) - lu(k,1751) * lu(k,2193)
         lu(k,2196) = lu(k,2196) - lu(k,1752) * lu(k,2193)
         lu(k,2197) = lu(k,2197) - lu(k,1753) * lu(k,2193)
         lu(k,2198) = lu(k,2198) - lu(k,1754) * lu(k,2193)
         lu(k,2199) = lu(k,2199) - lu(k,1755) * lu(k,2193)
         lu(k,2200) = lu(k,2200) - lu(k,1756) * lu(k,2193)
         lu(k,2201) = lu(k,2201) - lu(k,1757) * lu(k,2193)
         lu(k,2202) = lu(k,2202) - lu(k,1758) * lu(k,2193)
         lu(k,2203) = lu(k,2203) - lu(k,1759) * lu(k,2193)
         lu(k,2204) = lu(k,2204) - lu(k,1760) * lu(k,2193)
         lu(k,2218) = lu(k,2218) - lu(k,1750) * lu(k,2217)
         lu(k,2219) = lu(k,2219) - lu(k,1751) * lu(k,2217)
         lu(k,2220) = lu(k,2220) - lu(k,1752) * lu(k,2217)
         lu(k,2221) = lu(k,2221) - lu(k,1753) * lu(k,2217)
         lu(k,2222) = lu(k,2222) - lu(k,1754) * lu(k,2217)
         lu(k,2223) = lu(k,2223) - lu(k,1755) * lu(k,2217)
         lu(k,2224) = lu(k,2224) - lu(k,1756) * lu(k,2217)
         lu(k,2225) = lu(k,2225) - lu(k,1757) * lu(k,2217)
         lu(k,2226) = lu(k,2226) - lu(k,1758) * lu(k,2217)
         lu(k,2227) = lu(k,2227) - lu(k,1759) * lu(k,2217)
         lu(k,2228) = lu(k,2228) - lu(k,1760) * lu(k,2217)
         lu(k,2249) = lu(k,2249) - lu(k,1750) * lu(k,2248)
         lu(k,2250) = lu(k,2250) - lu(k,1751) * lu(k,2248)
         lu(k,2251) = lu(k,2251) - lu(k,1752) * lu(k,2248)
         lu(k,2252) = lu(k,2252) - lu(k,1753) * lu(k,2248)
         lu(k,2253) = lu(k,2253) - lu(k,1754) * lu(k,2248)
         lu(k,2254) = lu(k,2254) - lu(k,1755) * lu(k,2248)
         lu(k,2255) = lu(k,2255) - lu(k,1756) * lu(k,2248)
         lu(k,2256) = lu(k,2256) - lu(k,1757) * lu(k,2248)
         lu(k,2257) = lu(k,2257) - lu(k,1758) * lu(k,2248)
         lu(k,2258) = lu(k,2258) - lu(k,1759) * lu(k,2248)
         lu(k,2259) = lu(k,2259) - lu(k,1760) * lu(k,2248)
         lu(k,2275) = lu(k,2275) - lu(k,1750) * lu(k,2274)
         lu(k,2276) = lu(k,2276) - lu(k,1751) * lu(k,2274)
         lu(k,2277) = lu(k,2277) - lu(k,1752) * lu(k,2274)
         lu(k,2278) = lu(k,2278) - lu(k,1753) * lu(k,2274)
         lu(k,2279) = lu(k,2279) - lu(k,1754) * lu(k,2274)
         lu(k,2280) = lu(k,2280) - lu(k,1755) * lu(k,2274)
         lu(k,2281) = lu(k,2281) - lu(k,1756) * lu(k,2274)
         lu(k,2282) = lu(k,2282) - lu(k,1757) * lu(k,2274)
         lu(k,2283) = lu(k,2283) - lu(k,1758) * lu(k,2274)
         lu(k,2284) = lu(k,2284) - lu(k,1759) * lu(k,2274)
         lu(k,2285) = lu(k,2285) - lu(k,1760) * lu(k,2274)
                                                                        
      end do
                                                                        
      end subroutine lu_fac29
                                                                        
      subroutine lu_fac30( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,1842) = 1._r8 / lu(k,1842)
         lu(k,1843) = lu(k,1843) * lu(k,1842)
         lu(k,1844) = lu(k,1844) * lu(k,1842)
         lu(k,1845) = lu(k,1845) * lu(k,1842)
         lu(k,1846) = lu(k,1846) * lu(k,1842)
         lu(k,1847) = lu(k,1847) * lu(k,1842)
         lu(k,1848) = lu(k,1848) * lu(k,1842)
         lu(k,1849) = lu(k,1849) * lu(k,1842)
         lu(k,1850) = lu(k,1850) * lu(k,1842)
         lu(k,1851) = lu(k,1851) * lu(k,1842)
         lu(k,1852) = lu(k,1852) * lu(k,1842)
         lu(k,1950) = lu(k,1950) - lu(k,1843) * lu(k,1949)
         lu(k,1951) = lu(k,1951) - lu(k,1844) * lu(k,1949)
         lu(k,1952) = lu(k,1952) - lu(k,1845) * lu(k,1949)
         lu(k,1953) = lu(k,1953) - lu(k,1846) * lu(k,1949)
         lu(k,1954) = lu(k,1954) - lu(k,1847) * lu(k,1949)
         lu(k,1955) = lu(k,1955) - lu(k,1848) * lu(k,1949)
         lu(k,1956) = lu(k,1956) - lu(k,1849) * lu(k,1949)
         lu(k,1957) = lu(k,1957) - lu(k,1850) * lu(k,1949)
         lu(k,1958) = lu(k,1958) - lu(k,1851) * lu(k,1949)
         lu(k,1959) = lu(k,1959) - lu(k,1852) * lu(k,1949)
         lu(k,1976) = lu(k,1976) - lu(k,1843) * lu(k,1975)
         lu(k,1977) = lu(k,1977) - lu(k,1844) * lu(k,1975)
         lu(k,1978) = lu(k,1978) - lu(k,1845) * lu(k,1975)
         lu(k,1979) = lu(k,1979) - lu(k,1846) * lu(k,1975)
         lu(k,1980) = lu(k,1980) - lu(k,1847) * lu(k,1975)
         lu(k,1981) = lu(k,1981) - lu(k,1848) * lu(k,1975)
         lu(k,1982) = lu(k,1982) - lu(k,1849) * lu(k,1975)
         lu(k,1983) = lu(k,1983) - lu(k,1850) * lu(k,1975)
         lu(k,1984) = lu(k,1984) - lu(k,1851) * lu(k,1975)
         lu(k,1985) = lu(k,1985) - lu(k,1852) * lu(k,1975)
         lu(k,2015) = lu(k,2015) - lu(k,1843) * lu(k,2014)
         lu(k,2016) = lu(k,2016) - lu(k,1844) * lu(k,2014)
         lu(k,2017) = lu(k,2017) - lu(k,1845) * lu(k,2014)
         lu(k,2018) = lu(k,2018) - lu(k,1846) * lu(k,2014)
         lu(k,2019) = lu(k,2019) - lu(k,1847) * lu(k,2014)
         lu(k,2020) = lu(k,2020) - lu(k,1848) * lu(k,2014)
         lu(k,2021) = lu(k,2021) - lu(k,1849) * lu(k,2014)
         lu(k,2022) = lu(k,2022) - lu(k,1850) * lu(k,2014)
         lu(k,2023) = lu(k,2023) - lu(k,1851) * lu(k,2014)
         lu(k,2024) = lu(k,2024) - lu(k,1852) * lu(k,2014)
         lu(k,2067) = lu(k,2067) - lu(k,1843) * lu(k,2066)
         lu(k,2068) = lu(k,2068) - lu(k,1844) * lu(k,2066)
         lu(k,2069) = lu(k,2069) - lu(k,1845) * lu(k,2066)
         lu(k,2070) = lu(k,2070) - lu(k,1846) * lu(k,2066)
         lu(k,2071) = lu(k,2071) - lu(k,1847) * lu(k,2066)
         lu(k,2072) = lu(k,2072) - lu(k,1848) * lu(k,2066)
         lu(k,2073) = lu(k,2073) - lu(k,1849) * lu(k,2066)
         lu(k,2074) = lu(k,2074) - lu(k,1850) * lu(k,2066)
         lu(k,2075) = lu(k,2075) - lu(k,1851) * lu(k,2066)
         lu(k,2076) = lu(k,2076) - lu(k,1852) * lu(k,2066)
         lu(k,2128) = lu(k,2128) - lu(k,1843) * lu(k,2127)
         lu(k,2129) = lu(k,2129) - lu(k,1844) * lu(k,2127)
         lu(k,2130) = lu(k,2130) - lu(k,1845) * lu(k,2127)
         lu(k,2131) = lu(k,2131) - lu(k,1846) * lu(k,2127)
         lu(k,2132) = lu(k,2132) - lu(k,1847) * lu(k,2127)
         lu(k,2133) = lu(k,2133) - lu(k,1848) * lu(k,2127)
         lu(k,2134) = lu(k,2134) - lu(k,1849) * lu(k,2127)
         lu(k,2135) = lu(k,2135) - lu(k,1850) * lu(k,2127)
         lu(k,2136) = lu(k,2136) - lu(k,1851) * lu(k,2127)
         lu(k,2137) = lu(k,2137) - lu(k,1852) * lu(k,2127)
         lu(k,2151) = lu(k,2151) - lu(k,1843) * lu(k,2150)
         lu(k,2152) = lu(k,2152) - lu(k,1844) * lu(k,2150)
         lu(k,2153) = lu(k,2153) - lu(k,1845) * lu(k,2150)
         lu(k,2154) = lu(k,2154) - lu(k,1846) * lu(k,2150)
         lu(k,2155) = lu(k,2155) - lu(k,1847) * lu(k,2150)
         lu(k,2156) = lu(k,2156) - lu(k,1848) * lu(k,2150)
         lu(k,2157) = lu(k,2157) - lu(k,1849) * lu(k,2150)
         lu(k,2158) = lu(k,2158) - lu(k,1850) * lu(k,2150)
         lu(k,2159) = lu(k,2159) - lu(k,1851) * lu(k,2150)
         lu(k,2160) = lu(k,2160) - lu(k,1852) * lu(k,2150)
         lu(k,2195) = lu(k,2195) - lu(k,1843) * lu(k,2194)
         lu(k,2196) = lu(k,2196) - lu(k,1844) * lu(k,2194)
         lu(k,2197) = lu(k,2197) - lu(k,1845) * lu(k,2194)
         lu(k,2198) = lu(k,2198) - lu(k,1846) * lu(k,2194)
         lu(k,2199) = lu(k,2199) - lu(k,1847) * lu(k,2194)
         lu(k,2200) = lu(k,2200) - lu(k,1848) * lu(k,2194)
         lu(k,2201) = lu(k,2201) - lu(k,1849) * lu(k,2194)
         lu(k,2202) = lu(k,2202) - lu(k,1850) * lu(k,2194)
         lu(k,2203) = lu(k,2203) - lu(k,1851) * lu(k,2194)
         lu(k,2204) = lu(k,2204) - lu(k,1852) * lu(k,2194)
         lu(k,2219) = lu(k,2219) - lu(k,1843) * lu(k,2218)
         lu(k,2220) = lu(k,2220) - lu(k,1844) * lu(k,2218)
         lu(k,2221) = lu(k,2221) - lu(k,1845) * lu(k,2218)
         lu(k,2222) = lu(k,2222) - lu(k,1846) * lu(k,2218)
         lu(k,2223) = lu(k,2223) - lu(k,1847) * lu(k,2218)
         lu(k,2224) = lu(k,2224) - lu(k,1848) * lu(k,2218)
         lu(k,2225) = lu(k,2225) - lu(k,1849) * lu(k,2218)
         lu(k,2226) = lu(k,2226) - lu(k,1850) * lu(k,2218)
         lu(k,2227) = lu(k,2227) - lu(k,1851) * lu(k,2218)
         lu(k,2228) = lu(k,2228) - lu(k,1852) * lu(k,2218)
         lu(k,2250) = lu(k,2250) - lu(k,1843) * lu(k,2249)
         lu(k,2251) = lu(k,2251) - lu(k,1844) * lu(k,2249)
         lu(k,2252) = lu(k,2252) - lu(k,1845) * lu(k,2249)
         lu(k,2253) = lu(k,2253) - lu(k,1846) * lu(k,2249)
         lu(k,2254) = lu(k,2254) - lu(k,1847) * lu(k,2249)
         lu(k,2255) = lu(k,2255) - lu(k,1848) * lu(k,2249)
         lu(k,2256) = lu(k,2256) - lu(k,1849) * lu(k,2249)
         lu(k,2257) = lu(k,2257) - lu(k,1850) * lu(k,2249)
         lu(k,2258) = lu(k,2258) - lu(k,1851) * lu(k,2249)
         lu(k,2259) = lu(k,2259) - lu(k,1852) * lu(k,2249)
         lu(k,2276) = lu(k,2276) - lu(k,1843) * lu(k,2275)
         lu(k,2277) = lu(k,2277) - lu(k,1844) * lu(k,2275)
         lu(k,2278) = lu(k,2278) - lu(k,1845) * lu(k,2275)
         lu(k,2279) = lu(k,2279) - lu(k,1846) * lu(k,2275)
         lu(k,2280) = lu(k,2280) - lu(k,1847) * lu(k,2275)
         lu(k,2281) = lu(k,2281) - lu(k,1848) * lu(k,2275)
         lu(k,2282) = lu(k,2282) - lu(k,1849) * lu(k,2275)
         lu(k,2283) = lu(k,2283) - lu(k,1850) * lu(k,2275)
         lu(k,2284) = lu(k,2284) - lu(k,1851) * lu(k,2275)
         lu(k,2285) = lu(k,2285) - lu(k,1852) * lu(k,2275)
                                                                        
         lu(k,1950) = 1._r8 / lu(k,1950)
         lu(k,1951) = lu(k,1951) * lu(k,1950)
         lu(k,1952) = lu(k,1952) * lu(k,1950)
         lu(k,1953) = lu(k,1953) * lu(k,1950)
         lu(k,1954) = lu(k,1954) * lu(k,1950)
         lu(k,1955) = lu(k,1955) * lu(k,1950)
         lu(k,1956) = lu(k,1956) * lu(k,1950)
         lu(k,1957) = lu(k,1957) * lu(k,1950)
         lu(k,1958) = lu(k,1958) * lu(k,1950)
         lu(k,1959) = lu(k,1959) * lu(k,1950)
         lu(k,1977) = lu(k,1977) - lu(k,1951) * lu(k,1976)
         lu(k,1978) = lu(k,1978) - lu(k,1952) * lu(k,1976)
         lu(k,1979) = lu(k,1979) - lu(k,1953) * lu(k,1976)
         lu(k,1980) = lu(k,1980) - lu(k,1954) * lu(k,1976)
         lu(k,1981) = lu(k,1981) - lu(k,1955) * lu(k,1976)
         lu(k,1982) = lu(k,1982) - lu(k,1956) * lu(k,1976)
         lu(k,1983) = lu(k,1983) - lu(k,1957) * lu(k,1976)
         lu(k,1984) = lu(k,1984) - lu(k,1958) * lu(k,1976)
         lu(k,1985) = lu(k,1985) - lu(k,1959) * lu(k,1976)
         lu(k,2016) = lu(k,2016) - lu(k,1951) * lu(k,2015)
         lu(k,2017) = lu(k,2017) - lu(k,1952) * lu(k,2015)
         lu(k,2018) = lu(k,2018) - lu(k,1953) * lu(k,2015)
         lu(k,2019) = lu(k,2019) - lu(k,1954) * lu(k,2015)
         lu(k,2020) = lu(k,2020) - lu(k,1955) * lu(k,2015)
         lu(k,2021) = lu(k,2021) - lu(k,1956) * lu(k,2015)
         lu(k,2022) = lu(k,2022) - lu(k,1957) * lu(k,2015)
         lu(k,2023) = lu(k,2023) - lu(k,1958) * lu(k,2015)
         lu(k,2024) = lu(k,2024) - lu(k,1959) * lu(k,2015)
         lu(k,2068) = lu(k,2068) - lu(k,1951) * lu(k,2067)
         lu(k,2069) = lu(k,2069) - lu(k,1952) * lu(k,2067)
         lu(k,2070) = lu(k,2070) - lu(k,1953) * lu(k,2067)
         lu(k,2071) = lu(k,2071) - lu(k,1954) * lu(k,2067)
         lu(k,2072) = lu(k,2072) - lu(k,1955) * lu(k,2067)
         lu(k,2073) = lu(k,2073) - lu(k,1956) * lu(k,2067)
         lu(k,2074) = lu(k,2074) - lu(k,1957) * lu(k,2067)
         lu(k,2075) = lu(k,2075) - lu(k,1958) * lu(k,2067)
         lu(k,2076) = lu(k,2076) - lu(k,1959) * lu(k,2067)
         lu(k,2129) = lu(k,2129) - lu(k,1951) * lu(k,2128)
         lu(k,2130) = lu(k,2130) - lu(k,1952) * lu(k,2128)
         lu(k,2131) = lu(k,2131) - lu(k,1953) * lu(k,2128)
         lu(k,2132) = lu(k,2132) - lu(k,1954) * lu(k,2128)
         lu(k,2133) = lu(k,2133) - lu(k,1955) * lu(k,2128)
         lu(k,2134) = lu(k,2134) - lu(k,1956) * lu(k,2128)
         lu(k,2135) = lu(k,2135) - lu(k,1957) * lu(k,2128)
         lu(k,2136) = lu(k,2136) - lu(k,1958) * lu(k,2128)
         lu(k,2137) = lu(k,2137) - lu(k,1959) * lu(k,2128)
         lu(k,2152) = lu(k,2152) - lu(k,1951) * lu(k,2151)
         lu(k,2153) = lu(k,2153) - lu(k,1952) * lu(k,2151)
         lu(k,2154) = lu(k,2154) - lu(k,1953) * lu(k,2151)
         lu(k,2155) = lu(k,2155) - lu(k,1954) * lu(k,2151)
         lu(k,2156) = lu(k,2156) - lu(k,1955) * lu(k,2151)
         lu(k,2157) = lu(k,2157) - lu(k,1956) * lu(k,2151)
         lu(k,2158) = lu(k,2158) - lu(k,1957) * lu(k,2151)
         lu(k,2159) = lu(k,2159) - lu(k,1958) * lu(k,2151)
         lu(k,2160) = lu(k,2160) - lu(k,1959) * lu(k,2151)
         lu(k,2196) = lu(k,2196) - lu(k,1951) * lu(k,2195)
         lu(k,2197) = lu(k,2197) - lu(k,1952) * lu(k,2195)
         lu(k,2198) = lu(k,2198) - lu(k,1953) * lu(k,2195)
         lu(k,2199) = lu(k,2199) - lu(k,1954) * lu(k,2195)
         lu(k,2200) = lu(k,2200) - lu(k,1955) * lu(k,2195)
         lu(k,2201) = lu(k,2201) - lu(k,1956) * lu(k,2195)
         lu(k,2202) = lu(k,2202) - lu(k,1957) * lu(k,2195)
         lu(k,2203) = lu(k,2203) - lu(k,1958) * lu(k,2195)
         lu(k,2204) = lu(k,2204) - lu(k,1959) * lu(k,2195)
         lu(k,2220) = lu(k,2220) - lu(k,1951) * lu(k,2219)
         lu(k,2221) = lu(k,2221) - lu(k,1952) * lu(k,2219)
         lu(k,2222) = lu(k,2222) - lu(k,1953) * lu(k,2219)
         lu(k,2223) = lu(k,2223) - lu(k,1954) * lu(k,2219)
         lu(k,2224) = lu(k,2224) - lu(k,1955) * lu(k,2219)
         lu(k,2225) = lu(k,2225) - lu(k,1956) * lu(k,2219)
         lu(k,2226) = lu(k,2226) - lu(k,1957) * lu(k,2219)
         lu(k,2227) = lu(k,2227) - lu(k,1958) * lu(k,2219)
         lu(k,2228) = lu(k,2228) - lu(k,1959) * lu(k,2219)
         lu(k,2251) = lu(k,2251) - lu(k,1951) * lu(k,2250)
         lu(k,2252) = lu(k,2252) - lu(k,1952) * lu(k,2250)
         lu(k,2253) = lu(k,2253) - lu(k,1953) * lu(k,2250)
         lu(k,2254) = lu(k,2254) - lu(k,1954) * lu(k,2250)
         lu(k,2255) = lu(k,2255) - lu(k,1955) * lu(k,2250)
         lu(k,2256) = lu(k,2256) - lu(k,1956) * lu(k,2250)
         lu(k,2257) = lu(k,2257) - lu(k,1957) * lu(k,2250)
         lu(k,2258) = lu(k,2258) - lu(k,1958) * lu(k,2250)
         lu(k,2259) = lu(k,2259) - lu(k,1959) * lu(k,2250)
         lu(k,2277) = lu(k,2277) - lu(k,1951) * lu(k,2276)
         lu(k,2278) = lu(k,2278) - lu(k,1952) * lu(k,2276)
         lu(k,2279) = lu(k,2279) - lu(k,1953) * lu(k,2276)
         lu(k,2280) = lu(k,2280) - lu(k,1954) * lu(k,2276)
         lu(k,2281) = lu(k,2281) - lu(k,1955) * lu(k,2276)
         lu(k,2282) = lu(k,2282) - lu(k,1956) * lu(k,2276)
         lu(k,2283) = lu(k,2283) - lu(k,1957) * lu(k,2276)
         lu(k,2284) = lu(k,2284) - lu(k,1958) * lu(k,2276)
         lu(k,2285) = lu(k,2285) - lu(k,1959) * lu(k,2276)
                                                                        
         lu(k,1977) = 1._r8 / lu(k,1977)
         lu(k,1978) = lu(k,1978) * lu(k,1977)
         lu(k,1979) = lu(k,1979) * lu(k,1977)
         lu(k,1980) = lu(k,1980) * lu(k,1977)
         lu(k,1981) = lu(k,1981) * lu(k,1977)
         lu(k,1982) = lu(k,1982) * lu(k,1977)
         lu(k,1983) = lu(k,1983) * lu(k,1977)
         lu(k,1984) = lu(k,1984) * lu(k,1977)
         lu(k,1985) = lu(k,1985) * lu(k,1977)
         lu(k,2017) = lu(k,2017) - lu(k,1978) * lu(k,2016)
         lu(k,2018) = lu(k,2018) - lu(k,1979) * lu(k,2016)
         lu(k,2019) = lu(k,2019) - lu(k,1980) * lu(k,2016)
         lu(k,2020) = lu(k,2020) - lu(k,1981) * lu(k,2016)
         lu(k,2021) = lu(k,2021) - lu(k,1982) * lu(k,2016)
         lu(k,2022) = lu(k,2022) - lu(k,1983) * lu(k,2016)
         lu(k,2023) = lu(k,2023) - lu(k,1984) * lu(k,2016)
         lu(k,2024) = lu(k,2024) - lu(k,1985) * lu(k,2016)
         lu(k,2069) = lu(k,2069) - lu(k,1978) * lu(k,2068)
         lu(k,2070) = lu(k,2070) - lu(k,1979) * lu(k,2068)
         lu(k,2071) = lu(k,2071) - lu(k,1980) * lu(k,2068)
         lu(k,2072) = lu(k,2072) - lu(k,1981) * lu(k,2068)
         lu(k,2073) = lu(k,2073) - lu(k,1982) * lu(k,2068)
         lu(k,2074) = lu(k,2074) - lu(k,1983) * lu(k,2068)
         lu(k,2075) = lu(k,2075) - lu(k,1984) * lu(k,2068)
         lu(k,2076) = lu(k,2076) - lu(k,1985) * lu(k,2068)
         lu(k,2130) = lu(k,2130) - lu(k,1978) * lu(k,2129)
         lu(k,2131) = lu(k,2131) - lu(k,1979) * lu(k,2129)
         lu(k,2132) = lu(k,2132) - lu(k,1980) * lu(k,2129)
         lu(k,2133) = lu(k,2133) - lu(k,1981) * lu(k,2129)
         lu(k,2134) = lu(k,2134) - lu(k,1982) * lu(k,2129)
         lu(k,2135) = lu(k,2135) - lu(k,1983) * lu(k,2129)
         lu(k,2136) = lu(k,2136) - lu(k,1984) * lu(k,2129)
         lu(k,2137) = lu(k,2137) - lu(k,1985) * lu(k,2129)
         lu(k,2153) = lu(k,2153) - lu(k,1978) * lu(k,2152)
         lu(k,2154) = lu(k,2154) - lu(k,1979) * lu(k,2152)
         lu(k,2155) = lu(k,2155) - lu(k,1980) * lu(k,2152)
         lu(k,2156) = lu(k,2156) - lu(k,1981) * lu(k,2152)
         lu(k,2157) = lu(k,2157) - lu(k,1982) * lu(k,2152)
         lu(k,2158) = lu(k,2158) - lu(k,1983) * lu(k,2152)
         lu(k,2159) = lu(k,2159) - lu(k,1984) * lu(k,2152)
         lu(k,2160) = lu(k,2160) - lu(k,1985) * lu(k,2152)
         lu(k,2197) = lu(k,2197) - lu(k,1978) * lu(k,2196)
         lu(k,2198) = lu(k,2198) - lu(k,1979) * lu(k,2196)
         lu(k,2199) = lu(k,2199) - lu(k,1980) * lu(k,2196)
         lu(k,2200) = lu(k,2200) - lu(k,1981) * lu(k,2196)
         lu(k,2201) = lu(k,2201) - lu(k,1982) * lu(k,2196)
         lu(k,2202) = lu(k,2202) - lu(k,1983) * lu(k,2196)
         lu(k,2203) = lu(k,2203) - lu(k,1984) * lu(k,2196)
         lu(k,2204) = lu(k,2204) - lu(k,1985) * lu(k,2196)
         lu(k,2221) = lu(k,2221) - lu(k,1978) * lu(k,2220)
         lu(k,2222) = lu(k,2222) - lu(k,1979) * lu(k,2220)
         lu(k,2223) = lu(k,2223) - lu(k,1980) * lu(k,2220)
         lu(k,2224) = lu(k,2224) - lu(k,1981) * lu(k,2220)
         lu(k,2225) = lu(k,2225) - lu(k,1982) * lu(k,2220)
         lu(k,2226) = lu(k,2226) - lu(k,1983) * lu(k,2220)
         lu(k,2227) = lu(k,2227) - lu(k,1984) * lu(k,2220)
         lu(k,2228) = lu(k,2228) - lu(k,1985) * lu(k,2220)
         lu(k,2252) = lu(k,2252) - lu(k,1978) * lu(k,2251)
         lu(k,2253) = lu(k,2253) - lu(k,1979) * lu(k,2251)
         lu(k,2254) = lu(k,2254) - lu(k,1980) * lu(k,2251)
         lu(k,2255) = lu(k,2255) - lu(k,1981) * lu(k,2251)
         lu(k,2256) = lu(k,2256) - lu(k,1982) * lu(k,2251)
         lu(k,2257) = lu(k,2257) - lu(k,1983) * lu(k,2251)
         lu(k,2258) = lu(k,2258) - lu(k,1984) * lu(k,2251)
         lu(k,2259) = lu(k,2259) - lu(k,1985) * lu(k,2251)
         lu(k,2278) = lu(k,2278) - lu(k,1978) * lu(k,2277)
         lu(k,2279) = lu(k,2279) - lu(k,1979) * lu(k,2277)
         lu(k,2280) = lu(k,2280) - lu(k,1980) * lu(k,2277)
         lu(k,2281) = lu(k,2281) - lu(k,1981) * lu(k,2277)
         lu(k,2282) = lu(k,2282) - lu(k,1982) * lu(k,2277)
         lu(k,2283) = lu(k,2283) - lu(k,1983) * lu(k,2277)
         lu(k,2284) = lu(k,2284) - lu(k,1984) * lu(k,2277)
         lu(k,2285) = lu(k,2285) - lu(k,1985) * lu(k,2277)
                                                                        
         lu(k,2017) = 1._r8 / lu(k,2017)
         lu(k,2018) = lu(k,2018) * lu(k,2017)
         lu(k,2019) = lu(k,2019) * lu(k,2017)
         lu(k,2020) = lu(k,2020) * lu(k,2017)
         lu(k,2021) = lu(k,2021) * lu(k,2017)
         lu(k,2022) = lu(k,2022) * lu(k,2017)
         lu(k,2023) = lu(k,2023) * lu(k,2017)
         lu(k,2024) = lu(k,2024) * lu(k,2017)
         lu(k,2070) = lu(k,2070) - lu(k,2018) * lu(k,2069)
         lu(k,2071) = lu(k,2071) - lu(k,2019) * lu(k,2069)
         lu(k,2072) = lu(k,2072) - lu(k,2020) * lu(k,2069)
         lu(k,2073) = lu(k,2073) - lu(k,2021) * lu(k,2069)
         lu(k,2074) = lu(k,2074) - lu(k,2022) * lu(k,2069)
         lu(k,2075) = lu(k,2075) - lu(k,2023) * lu(k,2069)
         lu(k,2076) = lu(k,2076) - lu(k,2024) * lu(k,2069)
         lu(k,2131) = lu(k,2131) - lu(k,2018) * lu(k,2130)
         lu(k,2132) = lu(k,2132) - lu(k,2019) * lu(k,2130)
         lu(k,2133) = lu(k,2133) - lu(k,2020) * lu(k,2130)
         lu(k,2134) = lu(k,2134) - lu(k,2021) * lu(k,2130)
         lu(k,2135) = lu(k,2135) - lu(k,2022) * lu(k,2130)
         lu(k,2136) = lu(k,2136) - lu(k,2023) * lu(k,2130)
         lu(k,2137) = lu(k,2137) - lu(k,2024) * lu(k,2130)
         lu(k,2154) = lu(k,2154) - lu(k,2018) * lu(k,2153)
         lu(k,2155) = lu(k,2155) - lu(k,2019) * lu(k,2153)
         lu(k,2156) = lu(k,2156) - lu(k,2020) * lu(k,2153)
         lu(k,2157) = lu(k,2157) - lu(k,2021) * lu(k,2153)
         lu(k,2158) = lu(k,2158) - lu(k,2022) * lu(k,2153)
         lu(k,2159) = lu(k,2159) - lu(k,2023) * lu(k,2153)
         lu(k,2160) = lu(k,2160) - lu(k,2024) * lu(k,2153)
         lu(k,2198) = lu(k,2198) - lu(k,2018) * lu(k,2197)
         lu(k,2199) = lu(k,2199) - lu(k,2019) * lu(k,2197)
         lu(k,2200) = lu(k,2200) - lu(k,2020) * lu(k,2197)
         lu(k,2201) = lu(k,2201) - lu(k,2021) * lu(k,2197)
         lu(k,2202) = lu(k,2202) - lu(k,2022) * lu(k,2197)
         lu(k,2203) = lu(k,2203) - lu(k,2023) * lu(k,2197)
         lu(k,2204) = lu(k,2204) - lu(k,2024) * lu(k,2197)
         lu(k,2222) = lu(k,2222) - lu(k,2018) * lu(k,2221)
         lu(k,2223) = lu(k,2223) - lu(k,2019) * lu(k,2221)
         lu(k,2224) = lu(k,2224) - lu(k,2020) * lu(k,2221)
         lu(k,2225) = lu(k,2225) - lu(k,2021) * lu(k,2221)
         lu(k,2226) = lu(k,2226) - lu(k,2022) * lu(k,2221)
         lu(k,2227) = lu(k,2227) - lu(k,2023) * lu(k,2221)
         lu(k,2228) = lu(k,2228) - lu(k,2024) * lu(k,2221)
         lu(k,2253) = lu(k,2253) - lu(k,2018) * lu(k,2252)
         lu(k,2254) = lu(k,2254) - lu(k,2019) * lu(k,2252)
         lu(k,2255) = lu(k,2255) - lu(k,2020) * lu(k,2252)
         lu(k,2256) = lu(k,2256) - lu(k,2021) * lu(k,2252)
         lu(k,2257) = lu(k,2257) - lu(k,2022) * lu(k,2252)
         lu(k,2258) = lu(k,2258) - lu(k,2023) * lu(k,2252)
         lu(k,2259) = lu(k,2259) - lu(k,2024) * lu(k,2252)
         lu(k,2279) = lu(k,2279) - lu(k,2018) * lu(k,2278)
         lu(k,2280) = lu(k,2280) - lu(k,2019) * lu(k,2278)
         lu(k,2281) = lu(k,2281) - lu(k,2020) * lu(k,2278)
         lu(k,2282) = lu(k,2282) - lu(k,2021) * lu(k,2278)
         lu(k,2283) = lu(k,2283) - lu(k,2022) * lu(k,2278)
         lu(k,2284) = lu(k,2284) - lu(k,2023) * lu(k,2278)
         lu(k,2285) = lu(k,2285) - lu(k,2024) * lu(k,2278)
                                                                        
         lu(k,2070) = 1._r8 / lu(k,2070)
         lu(k,2071) = lu(k,2071) * lu(k,2070)
         lu(k,2072) = lu(k,2072) * lu(k,2070)
         lu(k,2073) = lu(k,2073) * lu(k,2070)
         lu(k,2074) = lu(k,2074) * lu(k,2070)
         lu(k,2075) = lu(k,2075) * lu(k,2070)
         lu(k,2076) = lu(k,2076) * lu(k,2070)
         lu(k,2132) = lu(k,2132) - lu(k,2071) * lu(k,2131)
         lu(k,2133) = lu(k,2133) - lu(k,2072) * lu(k,2131)
         lu(k,2134) = lu(k,2134) - lu(k,2073) * lu(k,2131)
         lu(k,2135) = lu(k,2135) - lu(k,2074) * lu(k,2131)
         lu(k,2136) = lu(k,2136) - lu(k,2075) * lu(k,2131)
         lu(k,2137) = lu(k,2137) - lu(k,2076) * lu(k,2131)
         lu(k,2155) = lu(k,2155) - lu(k,2071) * lu(k,2154)
         lu(k,2156) = lu(k,2156) - lu(k,2072) * lu(k,2154)
         lu(k,2157) = lu(k,2157) - lu(k,2073) * lu(k,2154)
         lu(k,2158) = lu(k,2158) - lu(k,2074) * lu(k,2154)
         lu(k,2159) = lu(k,2159) - lu(k,2075) * lu(k,2154)
         lu(k,2160) = lu(k,2160) - lu(k,2076) * lu(k,2154)
         lu(k,2199) = lu(k,2199) - lu(k,2071) * lu(k,2198)
         lu(k,2200) = lu(k,2200) - lu(k,2072) * lu(k,2198)
         lu(k,2201) = lu(k,2201) - lu(k,2073) * lu(k,2198)
         lu(k,2202) = lu(k,2202) - lu(k,2074) * lu(k,2198)
         lu(k,2203) = lu(k,2203) - lu(k,2075) * lu(k,2198)
         lu(k,2204) = lu(k,2204) - lu(k,2076) * lu(k,2198)
         lu(k,2223) = lu(k,2223) - lu(k,2071) * lu(k,2222)
         lu(k,2224) = lu(k,2224) - lu(k,2072) * lu(k,2222)
         lu(k,2225) = lu(k,2225) - lu(k,2073) * lu(k,2222)
         lu(k,2226) = lu(k,2226) - lu(k,2074) * lu(k,2222)
         lu(k,2227) = lu(k,2227) - lu(k,2075) * lu(k,2222)
         lu(k,2228) = lu(k,2228) - lu(k,2076) * lu(k,2222)
         lu(k,2254) = lu(k,2254) - lu(k,2071) * lu(k,2253)
         lu(k,2255) = lu(k,2255) - lu(k,2072) * lu(k,2253)
         lu(k,2256) = lu(k,2256) - lu(k,2073) * lu(k,2253)
         lu(k,2257) = lu(k,2257) - lu(k,2074) * lu(k,2253)
         lu(k,2258) = lu(k,2258) - lu(k,2075) * lu(k,2253)
         lu(k,2259) = lu(k,2259) - lu(k,2076) * lu(k,2253)
         lu(k,2280) = lu(k,2280) - lu(k,2071) * lu(k,2279)
         lu(k,2281) = lu(k,2281) - lu(k,2072) * lu(k,2279)
         lu(k,2282) = lu(k,2282) - lu(k,2073) * lu(k,2279)
         lu(k,2283) = lu(k,2283) - lu(k,2074) * lu(k,2279)
         lu(k,2284) = lu(k,2284) - lu(k,2075) * lu(k,2279)
         lu(k,2285) = lu(k,2285) - lu(k,2076) * lu(k,2279)
                                                                        
         lu(k,2132) = 1._r8 / lu(k,2132)
         lu(k,2133) = lu(k,2133) * lu(k,2132)
         lu(k,2134) = lu(k,2134) * lu(k,2132)
         lu(k,2135) = lu(k,2135) * lu(k,2132)
         lu(k,2136) = lu(k,2136) * lu(k,2132)
         lu(k,2137) = lu(k,2137) * lu(k,2132)
         lu(k,2156) = lu(k,2156) - lu(k,2133) * lu(k,2155)
         lu(k,2157) = lu(k,2157) - lu(k,2134) * lu(k,2155)
         lu(k,2158) = lu(k,2158) - lu(k,2135) * lu(k,2155)
         lu(k,2159) = lu(k,2159) - lu(k,2136) * lu(k,2155)
         lu(k,2160) = lu(k,2160) - lu(k,2137) * lu(k,2155)
         lu(k,2200) = lu(k,2200) - lu(k,2133) * lu(k,2199)
         lu(k,2201) = lu(k,2201) - lu(k,2134) * lu(k,2199)
         lu(k,2202) = lu(k,2202) - lu(k,2135) * lu(k,2199)
         lu(k,2203) = lu(k,2203) - lu(k,2136) * lu(k,2199)
         lu(k,2204) = lu(k,2204) - lu(k,2137) * lu(k,2199)
         lu(k,2224) = lu(k,2224) - lu(k,2133) * lu(k,2223)
         lu(k,2225) = lu(k,2225) - lu(k,2134) * lu(k,2223)
         lu(k,2226) = lu(k,2226) - lu(k,2135) * lu(k,2223)
         lu(k,2227) = lu(k,2227) - lu(k,2136) * lu(k,2223)
         lu(k,2228) = lu(k,2228) - lu(k,2137) * lu(k,2223)
         lu(k,2255) = lu(k,2255) - lu(k,2133) * lu(k,2254)
         lu(k,2256) = lu(k,2256) - lu(k,2134) * lu(k,2254)
         lu(k,2257) = lu(k,2257) - lu(k,2135) * lu(k,2254)
         lu(k,2258) = lu(k,2258) - lu(k,2136) * lu(k,2254)
         lu(k,2259) = lu(k,2259) - lu(k,2137) * lu(k,2254)
         lu(k,2281) = lu(k,2281) - lu(k,2133) * lu(k,2280)
         lu(k,2282) = lu(k,2282) - lu(k,2134) * lu(k,2280)
         lu(k,2283) = lu(k,2283) - lu(k,2135) * lu(k,2280)
         lu(k,2284) = lu(k,2284) - lu(k,2136) * lu(k,2280)
         lu(k,2285) = lu(k,2285) - lu(k,2137) * lu(k,2280)
                                                                        
      end do
                                                                        
      end subroutine lu_fac30
                                                                        
      subroutine lu_fac31( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,avec_len
         lu(k,2156) = 1._r8 / lu(k,2156)
         lu(k,2157) = lu(k,2157) * lu(k,2156)
         lu(k,2158) = lu(k,2158) * lu(k,2156)
         lu(k,2159) = lu(k,2159) * lu(k,2156)
         lu(k,2160) = lu(k,2160) * lu(k,2156)
         lu(k,2201) = lu(k,2201) - lu(k,2157) * lu(k,2200)
         lu(k,2202) = lu(k,2202) - lu(k,2158) * lu(k,2200)
         lu(k,2203) = lu(k,2203) - lu(k,2159) * lu(k,2200)
         lu(k,2204) = lu(k,2204) - lu(k,2160) * lu(k,2200)
         lu(k,2225) = lu(k,2225) - lu(k,2157) * lu(k,2224)
         lu(k,2226) = lu(k,2226) - lu(k,2158) * lu(k,2224)
         lu(k,2227) = lu(k,2227) - lu(k,2159) * lu(k,2224)
         lu(k,2228) = lu(k,2228) - lu(k,2160) * lu(k,2224)
         lu(k,2256) = lu(k,2256) - lu(k,2157) * lu(k,2255)
         lu(k,2257) = lu(k,2257) - lu(k,2158) * lu(k,2255)
         lu(k,2258) = lu(k,2258) - lu(k,2159) * lu(k,2255)
         lu(k,2259) = lu(k,2259) - lu(k,2160) * lu(k,2255)
         lu(k,2282) = lu(k,2282) - lu(k,2157) * lu(k,2281)
         lu(k,2283) = lu(k,2283) - lu(k,2158) * lu(k,2281)
         lu(k,2284) = lu(k,2284) - lu(k,2159) * lu(k,2281)
         lu(k,2285) = lu(k,2285) - lu(k,2160) * lu(k,2281)
                                                                        
         lu(k,2201) = 1._r8 / lu(k,2201)
         lu(k,2202) = lu(k,2202) * lu(k,2201)
         lu(k,2203) = lu(k,2203) * lu(k,2201)
         lu(k,2204) = lu(k,2204) * lu(k,2201)
         lu(k,2226) = lu(k,2226) - lu(k,2202) * lu(k,2225)
         lu(k,2227) = lu(k,2227) - lu(k,2203) * lu(k,2225)
         lu(k,2228) = lu(k,2228) - lu(k,2204) * lu(k,2225)
         lu(k,2257) = lu(k,2257) - lu(k,2202) * lu(k,2256)
         lu(k,2258) = lu(k,2258) - lu(k,2203) * lu(k,2256)
         lu(k,2259) = lu(k,2259) - lu(k,2204) * lu(k,2256)
         lu(k,2283) = lu(k,2283) - lu(k,2202) * lu(k,2282)
         lu(k,2284) = lu(k,2284) - lu(k,2203) * lu(k,2282)
         lu(k,2285) = lu(k,2285) - lu(k,2204) * lu(k,2282)
                                                                        
         lu(k,2226) = 1._r8 / lu(k,2226)
         lu(k,2227) = lu(k,2227) * lu(k,2226)
         lu(k,2228) = lu(k,2228) * lu(k,2226)
         lu(k,2258) = lu(k,2258) - lu(k,2227) * lu(k,2257)
         lu(k,2259) = lu(k,2259) - lu(k,2228) * lu(k,2257)
         lu(k,2284) = lu(k,2284) - lu(k,2227) * lu(k,2283)
         lu(k,2285) = lu(k,2285) - lu(k,2228) * lu(k,2283)
                                                                        
         lu(k,2258) = 1._r8 / lu(k,2258)
         lu(k,2259) = lu(k,2259) * lu(k,2258)
         lu(k,2285) = lu(k,2285) - lu(k,2259) * lu(k,2284)
                                                                        
         lu(k,2285) = 1._r8 / lu(k,2285)
                                                                        
      end do
                                                                        
      end subroutine lu_fac31
                                                                        
      subroutine lu_fac( avec_len, lu )
                                                                        
      use chem_mods, only : nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(inout) ::   lu(veclen,nzcnt)
                                                                        
      call lu_fac01( avec_len, lu )
      call lu_fac02( avec_len, lu )
      call lu_fac03( avec_len, lu )
      call lu_fac04( avec_len, lu )
      call lu_fac05( avec_len, lu )
      call lu_fac06( avec_len, lu )
      call lu_fac07( avec_len, lu )
      call lu_fac08( avec_len, lu )
      call lu_fac09( avec_len, lu )
      call lu_fac10( avec_len, lu )
      call lu_fac11( avec_len, lu )
      call lu_fac12( avec_len, lu )
      call lu_fac13( avec_len, lu )
      call lu_fac14( avec_len, lu )
      call lu_fac15( avec_len, lu )
      call lu_fac16( avec_len, lu )
      call lu_fac17( avec_len, lu )
      call lu_fac18( avec_len, lu )
      call lu_fac19( avec_len, lu )
      call lu_fac20( avec_len, lu )
      call lu_fac21( avec_len, lu )
      call lu_fac22( avec_len, lu )
      call lu_fac23( avec_len, lu )
      call lu_fac24( avec_len, lu )
      call lu_fac25( avec_len, lu )
      call lu_fac26( avec_len, lu )
      call lu_fac27( avec_len, lu )
      call lu_fac28( avec_len, lu )
      call lu_fac29( avec_len, lu )
      call lu_fac30( avec_len, lu )
      call lu_fac31( avec_len, lu )
                                                                        
      end subroutine lu_fac
                                                                        
      end module mo_lu_factor

      module mo_lu_solve

      use chem_mods, only: veclen
      private
      public :: lu_slv

      contains
                                                                        
      subroutine lu_slv01( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
         b(k,219) = b(k,219) - lu(k,94) * b(k,52)
         b(k,220) = b(k,220) - lu(k,95) * b(k,52)
                                                                        
         b(k,215) = b(k,215) - lu(k,97) * b(k,53)
         b(k,227) = b(k,227) - lu(k,98) * b(k,53)
                                                                        
         b(k,214) = b(k,214) - lu(k,100) * b(k,54)
         b(k,220) = b(k,220) - lu(k,101) * b(k,54)
                                                                        
         b(k,215) = b(k,215) - lu(k,103) * b(k,55)
         b(k,218) = b(k,218) - lu(k,104) * b(k,55)
                                                                        
         b(k,85) = b(k,85) - lu(k,106) * b(k,56)
         b(k,209) = b(k,209) - lu(k,107) * b(k,56)
         b(k,214) = b(k,214) - lu(k,108) * b(k,56)
                                                                        
         b(k,167) = b(k,167) - lu(k,110) * b(k,57)
         b(k,215) = b(k,215) - lu(k,111) * b(k,57)
         b(k,227) = b(k,227) - lu(k,112) * b(k,57)
                                                                        
         b(k,83) = b(k,83) - lu(k,114) * b(k,58)
         b(k,214) = b(k,214) - lu(k,115) * b(k,58)
         b(k,220) = b(k,220) - lu(k,116) * b(k,58)
                                                                        
         b(k,85) = b(k,85) - lu(k,118) * b(k,59)
         b(k,214) = b(k,214) - lu(k,119) * b(k,59)
         b(k,220) = b(k,220) - lu(k,120) * b(k,59)
                                                                        
         b(k,85) = b(k,85) - lu(k,122) * b(k,60)
         b(k,214) = b(k,214) - lu(k,123) * b(k,60)
         b(k,220) = b(k,220) - lu(k,124) * b(k,60)
                                                                        
         b(k,215) = b(k,215) - lu(k,126) * b(k,61)
         b(k,220) = b(k,220) - lu(k,127) * b(k,61)
         b(k,227) = b(k,227) - lu(k,128) * b(k,61)
                                                                        
         b(k,94) = b(k,94) - lu(k,130) * b(k,62)
         b(k,215) = b(k,215) - lu(k,131) * b(k,62)
                                                                        
         b(k,91) = b(k,91) - lu(k,133) * b(k,63)
         b(k,227) = b(k,227) - lu(k,134) * b(k,63)
                                                                        
         b(k,197) = b(k,197) - lu(k,136) * b(k,64)
         b(k,215) = b(k,215) - lu(k,137) * b(k,64)
                                                                        
         b(k,136) = b(k,136) - lu(k,139) * b(k,65)
         b(k,224) = b(k,224) - lu(k,140) * b(k,65)
                                                                        
         b(k,85) = b(k,85) - lu(k,142) * b(k,66)
         b(k,209) = b(k,209) - lu(k,143) * b(k,66)
         b(k,214) = b(k,214) - lu(k,144) * b(k,66)
         b(k,220) = b(k,220) - lu(k,145) * b(k,66)
                                                                        
         b(k,85) = b(k,85) - lu(k,147) * b(k,67)
         b(k,175) = b(k,175) - lu(k,148) * b(k,67)
         b(k,209) = b(k,209) - lu(k,149) * b(k,67)
         b(k,214) = b(k,214) - lu(k,150) * b(k,67)
                                                                        
         b(k,83) = b(k,83) - lu(k,152) * b(k,68)
         b(k,85) = b(k,85) - lu(k,153) * b(k,68)
         b(k,214) = b(k,214) - lu(k,154) * b(k,68)
         b(k,220) = b(k,220) - lu(k,155) * b(k,68)
                                                                        
         b(k,85) = b(k,85) - lu(k,157) * b(k,69)
         b(k,175) = b(k,175) - lu(k,158) * b(k,69)
         b(k,214) = b(k,214) - lu(k,159) * b(k,69)
         b(k,220) = b(k,220) - lu(k,160) * b(k,69)
                                                                        
         b(k,220) = b(k,220) - lu(k,162) * b(k,70)
                                                                        
         b(k,72) = b(k,72) - lu(k,165) * b(k,71)
         b(k,73) = b(k,73) - lu(k,166) * b(k,71)
         b(k,131) = b(k,131) - lu(k,167) * b(k,71)
         b(k,215) = b(k,215) - lu(k,168) * b(k,71)
         b(k,218) = b(k,218) - lu(k,169) * b(k,71)
                                                                        
         b(k,127) = b(k,127) - lu(k,171) * b(k,72)
         b(k,192) = b(k,192) - lu(k,172) * b(k,72)
         b(k,218) = b(k,218) - lu(k,173) * b(k,72)
                                                                        
         b(k,126) = b(k,126) - lu(k,175) * b(k,73)
         b(k,128) = b(k,128) - lu(k,176) * b(k,73)
         b(k,215) = b(k,215) - lu(k,177) * b(k,73)
         b(k,218) = b(k,218) - lu(k,178) * b(k,73)
                                                                        
         b(k,214) = b(k,214) - lu(k,180) * b(k,74)
         b(k,215) = b(k,215) - lu(k,181) * b(k,74)
         b(k,218) = b(k,218) - lu(k,182) * b(k,74)
                                                                        
         b(k,214) = b(k,214) - lu(k,184) * b(k,75)
         b(k,217) = b(k,217) - lu(k,185) * b(k,75)
                                                                        
         b(k,77) = b(k,77) - lu(k,188) * b(k,76)
         b(k,78) = b(k,78) - lu(k,189) * b(k,76)
         b(k,123) = b(k,123) - lu(k,190) * b(k,76)
         b(k,161) = b(k,161) - lu(k,191) * b(k,76)
         b(k,215) = b(k,215) - lu(k,192) * b(k,76)
         b(k,218) = b(k,218) - lu(k,193) * b(k,76)
                                                                        
         b(k,126) = b(k,126) - lu(k,195) * b(k,77)
         b(k,128) = b(k,128) - lu(k,196) * b(k,77)
         b(k,215) = b(k,215) - lu(k,197) * b(k,77)
         b(k,218) = b(k,218) - lu(k,198) * b(k,77)
                                                                        
         b(k,192) = b(k,192) - lu(k,200) * b(k,78)
         b(k,207) = b(k,207) - lu(k,201) * b(k,78)
         b(k,218) = b(k,218) - lu(k,202) * b(k,78)
                                                                        
         b(k,197) = b(k,197) - lu(k,204) * b(k,79)
         b(k,215) = b(k,215) - lu(k,205) * b(k,79)
                                                                        
         b(k,81) = b(k,81) - lu(k,209) * b(k,80)
         b(k,123) = b(k,123) - lu(k,210) * b(k,80)
         b(k,162) = b(k,162) - lu(k,211) * b(k,80)
         b(k,192) = b(k,192) - lu(k,212) * b(k,80)
         b(k,207) = b(k,207) - lu(k,213) * b(k,80)
         b(k,215) = b(k,215) - lu(k,214) * b(k,80)
         b(k,218) = b(k,218) - lu(k,215) * b(k,80)
                                                                        
         b(k,128) = b(k,128) - lu(k,217) * b(k,81)
         b(k,133) = b(k,133) - lu(k,218) * b(k,81)
         b(k,215) = b(k,215) - lu(k,219) * b(k,81)
         b(k,218) = b(k,218) - lu(k,220) * b(k,81)
                                                                        
         b(k,83) = b(k,83) - lu(k,222) * b(k,82)
         b(k,214) = b(k,214) - lu(k,223) * b(k,82)
         b(k,215) = b(k,215) - lu(k,224) * b(k,82)
         b(k,220) = b(k,220) - lu(k,225) * b(k,82)
                                                                        
         b(k,175) = b(k,175) - lu(k,227) * b(k,83)
         b(k,214) = b(k,214) - lu(k,228) * b(k,83)
         b(k,220) = b(k,220) - lu(k,229) * b(k,83)
                                                                        
         b(k,145) = b(k,145) - lu(k,231) * b(k,84)
         b(k,197) = b(k,197) - lu(k,232) * b(k,84)
         b(k,215) = b(k,215) - lu(k,233) * b(k,84)
         b(k,218) = b(k,218) - lu(k,234) * b(k,84)
                                                                        
         b(k,175) = b(k,175) - lu(k,236) * b(k,85)
         b(k,214) = b(k,214) - lu(k,237) * b(k,85)
                                                                        
         b(k,178) = b(k,178) - lu(k,239) * b(k,86)
         b(k,215) = b(k,215) - lu(k,240) * b(k,86)
                                                                        
         b(k,209) = b(k,209) - lu(k,242) * b(k,87)
         b(k,220) = b(k,220) - lu(k,243) * b(k,87)
                                                                        
         b(k,211) = b(k,211) - lu(k,245) * b(k,88)
         b(k,224) = b(k,224) - lu(k,246) * b(k,88)
                                                                        
         b(k,175) = b(k,175) - lu(k,249) * b(k,89)
         b(k,214) = b(k,214) - lu(k,250) * b(k,89)
         b(k,215) = b(k,215) - lu(k,251) * b(k,89)
         b(k,220) = b(k,220) - lu(k,252) * b(k,89)
                                                                        
         b(k,136) = b(k,136) - lu(k,254) * b(k,90)
         b(k,215) = b(k,215) - lu(k,255) * b(k,90)
                                                                        
         b(k,172) = b(k,172) - lu(k,258) * b(k,91)
         b(k,226) = b(k,226) - lu(k,259) * b(k,91)
         b(k,227) = b(k,227) - lu(k,260) * b(k,91)
                                                                        
         b(k,190) = b(k,190) - lu(k,262) * b(k,92)
         b(k,215) = b(k,215) - lu(k,263) * b(k,92)
         b(k,218) = b(k,218) - lu(k,264) * b(k,92)
                                                                        
         b(k,128) = b(k,128) - lu(k,266) * b(k,93)
         b(k,150) = b(k,150) - lu(k,267) * b(k,93)
         b(k,215) = b(k,215) - lu(k,268) * b(k,93)
                                                                        
         b(k,191) = b(k,191) - lu(k,270) * b(k,94)
         b(k,213) = b(k,213) - lu(k,271) * b(k,94)
         b(k,218) = b(k,218) - lu(k,272) * b(k,94)
                                                                        
         b(k,211) = b(k,211) - lu(k,274) * b(k,95)
         b(k,216) = b(k,216) - lu(k,275) * b(k,95)
         b(k,217) = b(k,217) - lu(k,276) * b(k,95)
         b(k,224) = b(k,224) - lu(k,277) * b(k,95)
         b(k,226) = b(k,226) - lu(k,278) * b(k,95)
                                                                        
         b(k,163) = b(k,163) - lu(k,280) * b(k,96)
         b(k,218) = b(k,218) - lu(k,281) * b(k,96)
                                                                        
         b(k,175) = b(k,175) - lu(k,283) * b(k,97)
         b(k,212) = b(k,212) - lu(k,284) * b(k,97)
                                                                        
         b(k,181) = b(k,181) - lu(k,286) * b(k,98)
         b(k,184) = b(k,184) - lu(k,287) * b(k,98)
         b(k,192) = b(k,192) - lu(k,288) * b(k,98)
         b(k,215) = b(k,215) - lu(k,289) * b(k,98)
         b(k,218) = b(k,218) - lu(k,290) * b(k,98)
                                                                        
         b(k,170) = b(k,170) - lu(k,292) * b(k,99)
         b(k,215) = b(k,215) - lu(k,293) * b(k,99)
         b(k,220) = b(k,220) - lu(k,294) * b(k,99)
         b(k,223) = b(k,223) - lu(k,295) * b(k,99)
         b(k,227) = b(k,227) - lu(k,296) * b(k,99)
                                                                        
         b(k,172) = b(k,172) - lu(k,298) * b(k,100)
         b(k,211) = b(k,211) - lu(k,299) * b(k,100)
         b(k,215) = b(k,215) - lu(k,300) * b(k,100)
         b(k,216) = b(k,216) - lu(k,301) * b(k,100)
         b(k,218) = b(k,218) - lu(k,302) * b(k,100)
                                                                        
         b(k,175) = b(k,175) - lu(k,305) * b(k,101)
         b(k,214) = b(k,214) - lu(k,306) * b(k,101)
         b(k,215) = b(k,215) - lu(k,307) * b(k,101)
         b(k,220) = b(k,220) - lu(k,308) * b(k,101)
         b(k,227) = b(k,227) - lu(k,309) * b(k,101)
                                                                        
         b(k,204) = b(k,204) - lu(k,311) * b(k,102)
         b(k,206) = b(k,206) - lu(k,312) * b(k,102)
         b(k,215) = b(k,215) - lu(k,313) * b(k,102)
         b(k,218) = b(k,218) - lu(k,314) * b(k,102)
                                                                        
         b(k,155) = b(k,155) - lu(k,316) * b(k,103)
         b(k,190) = b(k,190) - lu(k,317) * b(k,103)
         b(k,207) = b(k,207) - lu(k,318) * b(k,103)
         b(k,215) = b(k,215) - lu(k,319) * b(k,103)
                                                                        
         b(k,197) = b(k,197) - lu(k,321) * b(k,104)
         b(k,215) = b(k,215) - lu(k,322) * b(k,104)
                                                                        
         b(k,192) = b(k,192) - lu(k,324) * b(k,105)
         b(k,200) = b(k,200) - lu(k,325) * b(k,105)
         b(k,207) = b(k,207) - lu(k,326) * b(k,105)
         b(k,218) = b(k,218) - lu(k,327) * b(k,105)
                                                                        
         b(k,172) = b(k,172) - lu(k,329) * b(k,106)
         b(k,201) = b(k,201) - lu(k,330) * b(k,106)
         b(k,219) = b(k,219) - lu(k,331) * b(k,106)
         b(k,226) = b(k,226) - lu(k,332) * b(k,106)
                                                                        
         b(k,126) = b(k,126) - lu(k,334) * b(k,107)
         b(k,184) = b(k,184) - lu(k,335) * b(k,107)
         b(k,215) = b(k,215) - lu(k,336) * b(k,107)
         b(k,218) = b(k,218) - lu(k,337) * b(k,107)
                                                                        
         b(k,123) = b(k,123) - lu(k,340) * b(k,108)
         b(k,136) = b(k,136) - lu(k,341) * b(k,108)
         b(k,215) = b(k,215) - lu(k,342) * b(k,108)
         b(k,218) = b(k,218) - lu(k,343) * b(k,108)
                                                                        
         b(k,170) = b(k,170) - lu(k,345) * b(k,109)
         b(k,190) = b(k,190) - lu(k,346) * b(k,109)
         b(k,215) = b(k,215) - lu(k,347) * b(k,109)
         b(k,218) = b(k,218) - lu(k,348) * b(k,109)
                                                                        
         b(k,142) = b(k,142) - lu(k,350) * b(k,110)
         b(k,180) = b(k,180) - lu(k,351) * b(k,110)
         b(k,190) = b(k,190) - lu(k,352) * b(k,110)
         b(k,213) = b(k,213) - lu(k,353) * b(k,110)
         b(k,215) = b(k,215) - lu(k,354) * b(k,110)
         b(k,216) = b(k,216) - lu(k,355) * b(k,110)
         b(k,224) = b(k,224) - lu(k,356) * b(k,110)
                                                                        
         b(k,135) = b(k,135) - lu(k,358) * b(k,111)
         b(k,172) = b(k,172) - lu(k,359) * b(k,111)
         b(k,192) = b(k,192) - lu(k,360) * b(k,111)
         b(k,201) = b(k,201) - lu(k,361) * b(k,111)
         b(k,212) = b(k,212) - lu(k,362) * b(k,111)
         b(k,215) = b(k,215) - lu(k,363) * b(k,111)
         b(k,226) = b(k,226) - lu(k,364) * b(k,111)
                                                                        
      end do
                                                                        
      end subroutine lu_slv01
                                                                        
      subroutine lu_slv02( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
         b(k,192) = b(k,192) - lu(k,366) * b(k,112)
         b(k,215) = b(k,215) - lu(k,367) * b(k,112)
         b(k,218) = b(k,218) - lu(k,368) * b(k,112)
         b(k,220) = b(k,220) - lu(k,369) * b(k,112)
         b(k,221) = b(k,221) - lu(k,370) * b(k,112)
         b(k,223) = b(k,223) - lu(k,371) * b(k,112)
         b(k,227) = b(k,227) - lu(k,372) * b(k,112)
                                                                        
         b(k,174) = b(k,174) - lu(k,374) * b(k,113)
         b(k,191) = b(k,191) - lu(k,375) * b(k,113)
         b(k,211) = b(k,211) - lu(k,376) * b(k,113)
         b(k,215) = b(k,215) - lu(k,377) * b(k,113)
         b(k,218) = b(k,218) - lu(k,378) * b(k,113)
                                                                        
         b(k,184) = b(k,184) - lu(k,380) * b(k,114)
         b(k,192) = b(k,192) - lu(k,381) * b(k,114)
         b(k,200) = b(k,200) - lu(k,382) * b(k,114)
         b(k,207) = b(k,207) - lu(k,383) * b(k,114)
         b(k,218) = b(k,218) - lu(k,384) * b(k,114)
                                                                        
         b(k,212) = b(k,212) - lu(k,386) * b(k,115)
         b(k,213) = b(k,213) - lu(k,387) * b(k,115)
         b(k,215) = b(k,215) - lu(k,388) * b(k,115)
         b(k,221) = b(k,221) - lu(k,389) * b(k,115)
         b(k,227) = b(k,227) - lu(k,390) * b(k,115)
                                                                        
         b(k,183) = b(k,183) - lu(k,392) * b(k,116)
         b(k,187) = b(k,187) - lu(k,393) * b(k,116)
         b(k,211) = b(k,211) - lu(k,394) * b(k,116)
         b(k,215) = b(k,215) - lu(k,395) * b(k,116)
         b(k,224) = b(k,224) - lu(k,396) * b(k,116)
                                                                        
         b(k,157) = b(k,157) - lu(k,398) * b(k,117)
         b(k,174) = b(k,174) - lu(k,399) * b(k,117)
         b(k,215) = b(k,215) - lu(k,400) * b(k,117)
         b(k,218) = b(k,218) - lu(k,401) * b(k,117)
         b(k,224) = b(k,224) - lu(k,402) * b(k,117)
                                                                        
         b(k,215) = b(k,215) - lu(k,404) * b(k,118)
         b(k,216) = b(k,216) - lu(k,405) * b(k,118)
         b(k,218) = b(k,218) - lu(k,406) * b(k,118)
         b(k,224) = b(k,224) - lu(k,407) * b(k,118)
         b(k,227) = b(k,227) - lu(k,408) * b(k,118)
                                                                        
         b(k,194) = b(k,194) - lu(k,410) * b(k,119)
         b(k,207) = b(k,207) - lu(k,411) * b(k,119)
         b(k,213) = b(k,213) - lu(k,412) * b(k,119)
         b(k,215) = b(k,215) - lu(k,413) * b(k,119)
         b(k,227) = b(k,227) - lu(k,414) * b(k,119)
                                                                        
         b(k,167) = b(k,167) - lu(k,416) * b(k,120)
         b(k,180) = b(k,180) - lu(k,417) * b(k,120)
         b(k,215) = b(k,215) - lu(k,418) * b(k,120)
         b(k,218) = b(k,218) - lu(k,419) * b(k,120)
         b(k,227) = b(k,227) - lu(k,420) * b(k,120)
                                                                        
         b(k,127) = b(k,127) - lu(k,422) * b(k,121)
         b(k,131) = b(k,131) - lu(k,423) * b(k,121)
         b(k,184) = b(k,184) - lu(k,424) * b(k,121)
         b(k,215) = b(k,215) - lu(k,425) * b(k,121)
         b(k,218) = b(k,218) - lu(k,426) * b(k,121)
                                                                        
         b(k,133) = b(k,133) - lu(k,428) * b(k,122)
         b(k,184) = b(k,184) - lu(k,429) * b(k,122)
         b(k,200) = b(k,200) - lu(k,430) * b(k,122)
         b(k,215) = b(k,215) - lu(k,431) * b(k,122)
         b(k,218) = b(k,218) - lu(k,432) * b(k,122)
                                                                        
         b(k,136) = b(k,136) - lu(k,436) * b(k,123)
         b(k,215) = b(k,215) - lu(k,437) * b(k,123)
         b(k,217) = b(k,217) - lu(k,438) * b(k,123)
         b(k,218) = b(k,218) - lu(k,439) * b(k,123)
         b(k,224) = b(k,224) - lu(k,440) * b(k,123)
                                                                        
         b(k,181) = b(k,181) - lu(k,442) * b(k,124)
         b(k,213) = b(k,213) - lu(k,443) * b(k,124)
         b(k,217) = b(k,217) - lu(k,444) * b(k,124)
         b(k,218) = b(k,218) - lu(k,445) * b(k,124)
         b(k,224) = b(k,224) - lu(k,446) * b(k,124)
                                                                        
         b(k,209) = b(k,209) - lu(k,448) * b(k,125)
         b(k,214) = b(k,214) - lu(k,449) * b(k,125)
         b(k,215) = b(k,215) - lu(k,450) * b(k,125)
         b(k,220) = b(k,220) - lu(k,451) * b(k,125)
         b(k,223) = b(k,223) - lu(k,452) * b(k,125)
                                                                        
         b(k,184) = b(k,184) - lu(k,455) * b(k,126)
         b(k,215) = b(k,215) - lu(k,456) * b(k,126)
         b(k,217) = b(k,217) - lu(k,457) * b(k,126)
         b(k,218) = b(k,218) - lu(k,458) * b(k,126)
         b(k,224) = b(k,224) - lu(k,459) * b(k,126)
                                                                        
         b(k,160) = b(k,160) - lu(k,461) * b(k,127)
         b(k,218) = b(k,218) - lu(k,462) * b(k,127)
                                                                        
         b(k,150) = b(k,150) - lu(k,464) * b(k,128)
         b(k,222) = b(k,222) - lu(k,465) * b(k,128)
         b(k,224) = b(k,224) - lu(k,466) * b(k,128)
                                                                        
         b(k,209) = b(k,209) - lu(k,468) * b(k,129)
         b(k,214) = b(k,214) - lu(k,469) * b(k,129)
         b(k,215) = b(k,215) - lu(k,470) * b(k,129)
         b(k,220) = b(k,220) - lu(k,471) * b(k,129)
         b(k,223) = b(k,223) - lu(k,472) * b(k,129)
         b(k,227) = b(k,227) - lu(k,473) * b(k,129)
                                                                        
         b(k,179) = b(k,179) - lu(k,475) * b(k,130)
         b(k,180) = b(k,180) - lu(k,476) * b(k,130)
         b(k,183) = b(k,183) - lu(k,477) * b(k,130)
         b(k,213) = b(k,213) - lu(k,478) * b(k,130)
         b(k,215) = b(k,215) - lu(k,479) * b(k,130)
         b(k,218) = b(k,218) - lu(k,480) * b(k,130)
                                                                        
         b(k,160) = b(k,160) - lu(k,484) * b(k,131)
         b(k,184) = b(k,184) - lu(k,485) * b(k,131)
         b(k,215) = b(k,215) - lu(k,486) * b(k,131)
         b(k,217) = b(k,217) - lu(k,487) * b(k,131)
         b(k,218) = b(k,218) - lu(k,488) * b(k,131)
         b(k,224) = b(k,224) - lu(k,489) * b(k,131)
                                                                        
         b(k,212) = b(k,212) - lu(k,492) * b(k,132)
         b(k,214) = b(k,214) - lu(k,493) * b(k,132)
         b(k,215) = b(k,215) - lu(k,494) * b(k,132)
         b(k,217) = b(k,217) - lu(k,495) * b(k,132)
         b(k,224) = b(k,224) - lu(k,496) * b(k,132)
         b(k,226) = b(k,226) - lu(k,497) * b(k,132)
                                                                        
         b(k,184) = b(k,184) - lu(k,500) * b(k,133)
         b(k,200) = b(k,200) - lu(k,501) * b(k,133)
         b(k,215) = b(k,215) - lu(k,502) * b(k,133)
         b(k,217) = b(k,217) - lu(k,503) * b(k,133)
         b(k,218) = b(k,218) - lu(k,504) * b(k,133)
         b(k,224) = b(k,224) - lu(k,505) * b(k,133)
                                                                        
         b(k,155) = b(k,155) - lu(k,507) * b(k,134)
         b(k,170) = b(k,170) - lu(k,508) * b(k,134)
         b(k,207) = b(k,207) - lu(k,509) * b(k,134)
         b(k,215) = b(k,215) - lu(k,510) * b(k,134)
                                                                        
         b(k,201) = b(k,201) - lu(k,512) * b(k,135)
         b(k,212) = b(k,212) - lu(k,513) * b(k,135)
         b(k,215) = b(k,215) - lu(k,514) * b(k,135)
         b(k,222) = b(k,222) - lu(k,515) * b(k,135)
         b(k,226) = b(k,226) - lu(k,516) * b(k,135)
                                                                        
         b(k,150) = b(k,150) - lu(k,519) * b(k,136)
         b(k,215) = b(k,215) - lu(k,520) * b(k,136)
         b(k,217) = b(k,217) - lu(k,521) * b(k,136)
         b(k,218) = b(k,218) - lu(k,522) * b(k,136)
         b(k,224) = b(k,224) - lu(k,523) * b(k,136)
                                                                        
         b(k,168) = b(k,168) - lu(k,525) * b(k,137)
         b(k,207) = b(k,207) - lu(k,526) * b(k,137)
         b(k,213) = b(k,213) - lu(k,527) * b(k,137)
         b(k,215) = b(k,215) - lu(k,528) * b(k,137)
         b(k,216) = b(k,216) - lu(k,529) * b(k,137)
         b(k,221) = b(k,221) - lu(k,530) * b(k,137)
         b(k,224) = b(k,224) - lu(k,531) * b(k,137)
                                                                        
         b(k,164) = b(k,164) - lu(k,533) * b(k,138)
         b(k,190) = b(k,190) - lu(k,534) * b(k,138)
         b(k,195) = b(k,195) - lu(k,535) * b(k,138)
         b(k,213) = b(k,213) - lu(k,536) * b(k,138)
         b(k,215) = b(k,215) - lu(k,537) * b(k,138)
         b(k,218) = b(k,218) - lu(k,538) * b(k,138)
         b(k,227) = b(k,227) - lu(k,539) * b(k,138)
                                                                        
         b(k,165) = b(k,165) - lu(k,541) * b(k,139)
         b(k,209) = b(k,209) - lu(k,542) * b(k,139)
         b(k,211) = b(k,211) - lu(k,543) * b(k,139)
         b(k,216) = b(k,216) - lu(k,544) * b(k,139)
         b(k,224) = b(k,224) - lu(k,545) * b(k,139)
         b(k,225) = b(k,225) - lu(k,546) * b(k,139)
         b(k,226) = b(k,226) - lu(k,547) * b(k,139)
                                                                        
         b(k,159) = b(k,159) - lu(k,549) * b(k,140)
         b(k,181) = b(k,181) - lu(k,550) * b(k,140)
         b(k,192) = b(k,192) - lu(k,551) * b(k,140)
         b(k,213) = b(k,213) - lu(k,552) * b(k,140)
         b(k,215) = b(k,215) - lu(k,553) * b(k,140)
         b(k,218) = b(k,218) - lu(k,554) * b(k,140)
         b(k,222) = b(k,222) - lu(k,555) * b(k,140)
                                                                        
         b(k,174) = b(k,174) - lu(k,557) * b(k,141)
         b(k,191) = b(k,191) - lu(k,558) * b(k,141)
         b(k,195) = b(k,195) - lu(k,559) * b(k,141)
         b(k,196) = b(k,196) - lu(k,560) * b(k,141)
         b(k,211) = b(k,211) - lu(k,561) * b(k,141)
         b(k,215) = b(k,215) - lu(k,562) * b(k,141)
         b(k,218) = b(k,218) - lu(k,563) * b(k,141)
                                                                        
         b(k,180) = b(k,180) - lu(k,565) * b(k,142)
         b(k,190) = b(k,190) - lu(k,566) * b(k,142)
         b(k,196) = b(k,196) - lu(k,567) * b(k,142)
         b(k,213) = b(k,213) - lu(k,568) * b(k,142)
         b(k,217) = b(k,217) - lu(k,569) * b(k,142)
         b(k,218) = b(k,218) - lu(k,570) * b(k,142)
         b(k,224) = b(k,224) - lu(k,571) * b(k,142)
                                                                        
         b(k,168) = b(k,168) - lu(k,573) * b(k,143)
         b(k,195) = b(k,195) - lu(k,574) * b(k,143)
         b(k,206) = b(k,206) - lu(k,575) * b(k,143)
         b(k,213) = b(k,213) - lu(k,576) * b(k,143)
         b(k,215) = b(k,215) - lu(k,577) * b(k,143)
         b(k,216) = b(k,216) - lu(k,578) * b(k,143)
         b(k,218) = b(k,218) - lu(k,579) * b(k,143)
         b(k,224) = b(k,224) - lu(k,580) * b(k,143)
                                                                        
         b(k,191) = b(k,191) - lu(k,582) * b(k,144)
         b(k,195) = b(k,195) - lu(k,583) * b(k,144)
         b(k,196) = b(k,196) - lu(k,584) * b(k,144)
         b(k,211) = b(k,211) - lu(k,585) * b(k,144)
         b(k,213) = b(k,213) - lu(k,586) * b(k,144)
         b(k,215) = b(k,215) - lu(k,587) * b(k,144)
         b(k,218) = b(k,218) - lu(k,588) * b(k,144)
         b(k,224) = b(k,224) - lu(k,589) * b(k,144)
                                                                        
         b(k,176) = b(k,176) - lu(k,591) * b(k,145)
         b(k,192) = b(k,192) - lu(k,592) * b(k,145)
         b(k,218) = b(k,218) - lu(k,593) * b(k,145)
                                                                        
         b(k,209) = b(k,209) - lu(k,595) * b(k,146)
         b(k,214) = b(k,214) - lu(k,596) * b(k,146)
         b(k,215) = b(k,215) - lu(k,597) * b(k,146)
         b(k,218) = b(k,218) - lu(k,598) * b(k,146)
         b(k,220) = b(k,220) - lu(k,599) * b(k,146)
         b(k,221) = b(k,221) - lu(k,600) * b(k,146)
         b(k,223) = b(k,223) - lu(k,601) * b(k,146)
         b(k,227) = b(k,227) - lu(k,602) * b(k,146)
                                                                        
         b(k,215) = b(k,215) - lu(k,604) * b(k,147)
         b(k,218) = b(k,218) - lu(k,605) * b(k,147)
         b(k,220) = b(k,220) - lu(k,606) * b(k,147)
         b(k,223) = b(k,223) - lu(k,607) * b(k,147)
         b(k,226) = b(k,226) - lu(k,608) * b(k,147)
         b(k,227) = b(k,227) - lu(k,609) * b(k,147)
                                                                        
      end do
                                                                        
      end subroutine lu_slv02
                                                                        
      subroutine lu_slv03( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
         b(k,168) = b(k,168) - lu(k,611) * b(k,148)
         b(k,207) = b(k,207) - lu(k,612) * b(k,148)
         b(k,213) = b(k,213) - lu(k,613) * b(k,148)
         b(k,215) = b(k,215) - lu(k,614) * b(k,148)
         b(k,221) = b(k,221) - lu(k,615) * b(k,148)
         b(k,227) = b(k,227) - lu(k,616) * b(k,148)
                                                                        
         b(k,183) = b(k,183) - lu(k,618) * b(k,149)
         b(k,211) = b(k,211) - lu(k,619) * b(k,149)
         b(k,215) = b(k,215) - lu(k,620) * b(k,149)
         b(k,218) = b(k,218) - lu(k,621) * b(k,149)
         b(k,224) = b(k,224) - lu(k,622) * b(k,149)
                                                                        
         b(k,215) = b(k,215) - lu(k,626) * b(k,150)
         b(k,217) = b(k,217) - lu(k,627) * b(k,150)
         b(k,218) = b(k,218) - lu(k,628) * b(k,150)
         b(k,222) = b(k,222) - lu(k,629) * b(k,150)
         b(k,224) = b(k,224) - lu(k,630) * b(k,150)
                                                                        
         b(k,155) = b(k,155) - lu(k,633) * b(k,151)
         b(k,170) = b(k,170) - lu(k,634) * b(k,151)
         b(k,178) = b(k,178) - lu(k,635) * b(k,151)
         b(k,180) = b(k,180) - lu(k,636) * b(k,151)
         b(k,190) = b(k,190) - lu(k,637) * b(k,151)
         b(k,207) = b(k,207) - lu(k,638) * b(k,151)
         b(k,213) = b(k,213) - lu(k,639) * b(k,151)
         b(k,215) = b(k,215) - lu(k,640) * b(k,151)
         b(k,218) = b(k,218) - lu(k,641) * b(k,151)
                                                                        
         b(k,168) = b(k,168) - lu(k,643) * b(k,152)
         b(k,180) = b(k,180) - lu(k,644) * b(k,152)
         b(k,188) = b(k,188) - lu(k,645) * b(k,152)
         b(k,191) = b(k,191) - lu(k,646) * b(k,152)
         b(k,192) = b(k,192) - lu(k,647) * b(k,152)
         b(k,193) = b(k,193) - lu(k,648) * b(k,152)
         b(k,213) = b(k,213) - lu(k,649) * b(k,152)
         b(k,215) = b(k,215) - lu(k,650) * b(k,152)
         b(k,218) = b(k,218) - lu(k,651) * b(k,152)
                                                                        
         b(k,160) = b(k,160) - lu(k,656) * b(k,153)
         b(k,161) = b(k,161) - lu(k,657) * b(k,153)
         b(k,163) = b(k,163) - lu(k,658) * b(k,153)
         b(k,176) = b(k,176) - lu(k,659) * b(k,153)
         b(k,184) = b(k,184) - lu(k,660) * b(k,153)
         b(k,192) = b(k,192) - lu(k,661) * b(k,153)
         b(k,200) = b(k,200) - lu(k,662) * b(k,153)
         b(k,215) = b(k,215) - lu(k,663) * b(k,153)
         b(k,218) = b(k,218) - lu(k,664) * b(k,153)
                                                                        
         b(k,155) = b(k,155) - lu(k,667) * b(k,154)
         b(k,170) = b(k,170) - lu(k,668) * b(k,154)
         b(k,180) = b(k,180) - lu(k,669) * b(k,154)
         b(k,190) = b(k,190) - lu(k,670) * b(k,154)
         b(k,207) = b(k,207) - lu(k,671) * b(k,154)
         b(k,213) = b(k,213) - lu(k,672) * b(k,154)
         b(k,215) = b(k,215) - lu(k,673) * b(k,154)
         b(k,218) = b(k,218) - lu(k,674) * b(k,154)
         b(k,224) = b(k,224) - lu(k,675) * b(k,154)
                                                                        
         b(k,190) = b(k,190) - lu(k,678) * b(k,155)
         b(k,207) = b(k,207) - lu(k,679) * b(k,155)
         b(k,215) = b(k,215) - lu(k,680) * b(k,155)
         b(k,217) = b(k,217) - lu(k,681) * b(k,155)
         b(k,218) = b(k,218) - lu(k,682) * b(k,155)
         b(k,224) = b(k,224) - lu(k,683) * b(k,155)
                                                                        
         b(k,168) = b(k,168) - lu(k,685) * b(k,156)
         b(k,215) = b(k,215) - lu(k,686) * b(k,156)
         b(k,221) = b(k,221) - lu(k,687) * b(k,156)
         b(k,227) = b(k,227) - lu(k,688) * b(k,156)
                                                                        
         b(k,197) = b(k,197) - lu(k,691) * b(k,157)
         b(k,199) = b(k,199) - lu(k,692) * b(k,157)
         b(k,205) = b(k,205) - lu(k,693) * b(k,157)
         b(k,213) = b(k,213) - lu(k,694) * b(k,157)
         b(k,215) = b(k,215) - lu(k,695) * b(k,157)
         b(k,218) = b(k,218) - lu(k,696) * b(k,157)
                                                                        
         b(k,160) = b(k,160) - lu(k,702) * b(k,158)
         b(k,162) = b(k,162) - lu(k,703) * b(k,158)
         b(k,163) = b(k,163) - lu(k,704) * b(k,158)
         b(k,176) = b(k,176) - lu(k,705) * b(k,158)
         b(k,184) = b(k,184) - lu(k,706) * b(k,158)
         b(k,192) = b(k,192) - lu(k,707) * b(k,158)
         b(k,200) = b(k,200) - lu(k,708) * b(k,158)
         b(k,207) = b(k,207) - lu(k,709) * b(k,158)
         b(k,215) = b(k,215) - lu(k,710) * b(k,158)
         b(k,218) = b(k,218) - lu(k,711) * b(k,158)
                                                                        
         b(k,191) = b(k,191) - lu(k,715) * b(k,159)
         b(k,213) = b(k,213) - lu(k,716) * b(k,159)
         b(k,215) = b(k,215) - lu(k,717) * b(k,159)
         b(k,217) = b(k,217) - lu(k,718) * b(k,159)
         b(k,218) = b(k,218) - lu(k,719) * b(k,159)
         b(k,224) = b(k,224) - lu(k,720) * b(k,159)
                                                                        
         b(k,184) = b(k,184) - lu(k,722) * b(k,160)
         b(k,192) = b(k,192) - lu(k,723) * b(k,160)
         b(k,217) = b(k,217) - lu(k,724) * b(k,160)
         b(k,218) = b(k,218) - lu(k,725) * b(k,160)
         b(k,224) = b(k,224) - lu(k,726) * b(k,160)
                                                                        
         b(k,163) = b(k,163) - lu(k,733) * b(k,161)
         b(k,176) = b(k,176) - lu(k,734) * b(k,161)
         b(k,184) = b(k,184) - lu(k,735) * b(k,161)
         b(k,192) = b(k,192) - lu(k,736) * b(k,161)
         b(k,200) = b(k,200) - lu(k,737) * b(k,161)
         b(k,215) = b(k,215) - lu(k,738) * b(k,161)
         b(k,217) = b(k,217) - lu(k,739) * b(k,161)
         b(k,218) = b(k,218) - lu(k,740) * b(k,161)
         b(k,224) = b(k,224) - lu(k,741) * b(k,161)
                                                                        
         b(k,163) = b(k,163) - lu(k,749) * b(k,162)
         b(k,176) = b(k,176) - lu(k,750) * b(k,162)
         b(k,184) = b(k,184) - lu(k,751) * b(k,162)
         b(k,192) = b(k,192) - lu(k,752) * b(k,162)
         b(k,200) = b(k,200) - lu(k,753) * b(k,162)
         b(k,207) = b(k,207) - lu(k,754) * b(k,162)
         b(k,215) = b(k,215) - lu(k,755) * b(k,162)
         b(k,217) = b(k,217) - lu(k,756) * b(k,162)
         b(k,218) = b(k,218) - lu(k,757) * b(k,162)
         b(k,224) = b(k,224) - lu(k,758) * b(k,162)
                                                                        
         b(k,192) = b(k,192) - lu(k,760) * b(k,163)
         b(k,200) = b(k,200) - lu(k,761) * b(k,163)
         b(k,215) = b(k,215) - lu(k,762) * b(k,163)
         b(k,217) = b(k,217) - lu(k,763) * b(k,163)
         b(k,218) = b(k,218) - lu(k,764) * b(k,163)
         b(k,221) = b(k,221) - lu(k,765) * b(k,163)
         b(k,224) = b(k,224) - lu(k,766) * b(k,163)
                                                                        
         b(k,190) = b(k,190) - lu(k,769) * b(k,164)
         b(k,195) = b(k,195) - lu(k,770) * b(k,164)
         b(k,213) = b(k,213) - lu(k,771) * b(k,164)
         b(k,215) = b(k,215) - lu(k,772) * b(k,164)
         b(k,217) = b(k,217) - lu(k,773) * b(k,164)
         b(k,218) = b(k,218) - lu(k,774) * b(k,164)
         b(k,224) = b(k,224) - lu(k,775) * b(k,164)
         b(k,227) = b(k,227) - lu(k,776) * b(k,164)
                                                                        
         b(k,209) = b(k,209) - lu(k,779) * b(k,165)
         b(k,215) = b(k,215) - lu(k,780) * b(k,165)
         b(k,220) = b(k,220) - lu(k,781) * b(k,165)
         b(k,223) = b(k,223) - lu(k,782) * b(k,165)
         b(k,225) = b(k,225) - lu(k,783) * b(k,165)
         b(k,226) = b(k,226) - lu(k,784) * b(k,165)
         b(k,227) = b(k,227) - lu(k,785) * b(k,165)
                                                                        
         b(k,213) = b(k,213) - lu(k,787) * b(k,166)
         b(k,215) = b(k,215) - lu(k,788) * b(k,166)
         b(k,218) = b(k,218) - lu(k,789) * b(k,166)
                                                                        
         b(k,180) = b(k,180) - lu(k,792) * b(k,167)
         b(k,190) = b(k,190) - lu(k,793) * b(k,167)
         b(k,213) = b(k,213) - lu(k,794) * b(k,167)
         b(k,215) = b(k,215) - lu(k,795) * b(k,167)
         b(k,217) = b(k,217) - lu(k,796) * b(k,167)
         b(k,218) = b(k,218) - lu(k,797) * b(k,167)
         b(k,221) = b(k,221) - lu(k,798) * b(k,167)
         b(k,224) = b(k,224) - lu(k,799) * b(k,167)
         b(k,227) = b(k,227) - lu(k,800) * b(k,167)
                                                                        
         b(k,192) = b(k,192) - lu(k,802) * b(k,168)
         b(k,226) = b(k,226) - lu(k,803) * b(k,168)
                                                                        
         b(k,209) = b(k,209) - lu(k,805) * b(k,169)
         b(k,212) = b(k,212) - lu(k,806) * b(k,169)
         b(k,214) = b(k,214) - lu(k,807) * b(k,169)
         b(k,215) = b(k,215) - lu(k,808) * b(k,169)
         b(k,225) = b(k,225) - lu(k,809) * b(k,169)
         b(k,226) = b(k,226) - lu(k,810) * b(k,169)
         b(k,227) = b(k,227) - lu(k,811) * b(k,169)
                                                                        
         b(k,190) = b(k,190) - lu(k,816) * b(k,170)
         b(k,213) = b(k,213) - lu(k,817) * b(k,170)
         b(k,215) = b(k,215) - lu(k,818) * b(k,170)
         b(k,217) = b(k,217) - lu(k,819) * b(k,170)
         b(k,218) = b(k,218) - lu(k,820) * b(k,170)
         b(k,221) = b(k,221) - lu(k,821) * b(k,170)
         b(k,224) = b(k,224) - lu(k,822) * b(k,170)
                                                                        
         b(k,215) = b(k,215) - lu(k,825) * b(k,171)
         b(k,219) = b(k,219) - lu(k,826) * b(k,171)
         b(k,220) = b(k,220) - lu(k,827) * b(k,171)
         b(k,223) = b(k,223) - lu(k,828) * b(k,171)
         b(k,226) = b(k,226) - lu(k,829) * b(k,171)
         b(k,227) = b(k,227) - lu(k,830) * b(k,171)
                                                                        
         b(k,201) = b(k,201) - lu(k,833) * b(k,172)
         b(k,215) = b(k,215) - lu(k,834) * b(k,172)
         b(k,218) = b(k,218) - lu(k,835) * b(k,172)
         b(k,226) = b(k,226) - lu(k,836) * b(k,172)
         b(k,227) = b(k,227) - lu(k,837) * b(k,172)
                                                                        
         b(k,181) = b(k,181) - lu(k,842) * b(k,173)
         b(k,186) = b(k,186) - lu(k,843) * b(k,173)
         b(k,192) = b(k,192) - lu(k,844) * b(k,173)
         b(k,198) = b(k,198) - lu(k,845) * b(k,173)
         b(k,199) = b(k,199) - lu(k,846) * b(k,173)
         b(k,202) = b(k,202) - lu(k,847) * b(k,173)
         b(k,203) = b(k,203) - lu(k,848) * b(k,173)
         b(k,205) = b(k,205) - lu(k,849) * b(k,173)
         b(k,207) = b(k,207) - lu(k,850) * b(k,173)
         b(k,213) = b(k,213) - lu(k,851) * b(k,173)
         b(k,215) = b(k,215) - lu(k,852) * b(k,173)
         b(k,216) = b(k,216) - lu(k,853) * b(k,173)
         b(k,218) = b(k,218) - lu(k,854) * b(k,173)
         b(k,221) = b(k,221) - lu(k,855) * b(k,173)
         b(k,222) = b(k,222) - lu(k,856) * b(k,173)
                                                                        
         b(k,200) = b(k,200) - lu(k,858) * b(k,174)
         b(k,207) = b(k,207) - lu(k,859) * b(k,174)
         b(k,213) = b(k,213) - lu(k,860) * b(k,174)
         b(k,215) = b(k,215) - lu(k,861) * b(k,174)
         b(k,224) = b(k,224) - lu(k,862) * b(k,174)
                                                                        
         b(k,208) = b(k,208) - lu(k,865) * b(k,175)
         b(k,210) = b(k,210) - lu(k,866) * b(k,175)
         b(k,211) = b(k,211) - lu(k,867) * b(k,175)
         b(k,212) = b(k,212) - lu(k,868) * b(k,175)
         b(k,215) = b(k,215) - lu(k,869) * b(k,175)
         b(k,216) = b(k,216) - lu(k,870) * b(k,175)
         b(k,221) = b(k,221) - lu(k,871) * b(k,175)
         b(k,227) = b(k,227) - lu(k,872) * b(k,175)
                                                                        
         b(k,184) = b(k,184) - lu(k,874) * b(k,176)
         b(k,192) = b(k,192) - lu(k,875) * b(k,176)
         b(k,200) = b(k,200) - lu(k,876) * b(k,176)
         b(k,215) = b(k,215) - lu(k,877) * b(k,176)
         b(k,217) = b(k,217) - lu(k,878) * b(k,176)
         b(k,218) = b(k,218) - lu(k,879) * b(k,176)
         b(k,221) = b(k,221) - lu(k,880) * b(k,176)
         b(k,224) = b(k,224) - lu(k,881) * b(k,176)
                                                                        
      end do
                                                                        
      end subroutine lu_slv03
                                                                        
      subroutine lu_slv04( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
         b(k,211) = b(k,211) - lu(k,885) * b(k,177)
         b(k,215) = b(k,215) - lu(k,886) * b(k,177)
         b(k,216) = b(k,216) - lu(k,887) * b(k,177)
         b(k,219) = b(k,219) - lu(k,888) * b(k,177)
         b(k,220) = b(k,220) - lu(k,889) * b(k,177)
         b(k,223) = b(k,223) - lu(k,890) * b(k,177)
         b(k,224) = b(k,224) - lu(k,891) * b(k,177)
         b(k,226) = b(k,226) - lu(k,892) * b(k,177)
         b(k,227) = b(k,227) - lu(k,893) * b(k,177)
                                                                        
         b(k,180) = b(k,180) - lu(k,900) * b(k,178)
         b(k,190) = b(k,190) - lu(k,901) * b(k,178)
         b(k,207) = b(k,207) - lu(k,902) * b(k,178)
         b(k,213) = b(k,213) - lu(k,903) * b(k,178)
         b(k,215) = b(k,215) - lu(k,904) * b(k,178)
         b(k,217) = b(k,217) - lu(k,905) * b(k,178)
         b(k,218) = b(k,218) - lu(k,906) * b(k,178)
         b(k,221) = b(k,221) - lu(k,907) * b(k,178)
         b(k,224) = b(k,224) - lu(k,908) * b(k,178)
                                                                        
         b(k,180) = b(k,180) - lu(k,913) * b(k,179)
         b(k,183) = b(k,183) - lu(k,914) * b(k,179)
         b(k,211) = b(k,211) - lu(k,915) * b(k,179)
         b(k,213) = b(k,213) - lu(k,916) * b(k,179)
         b(k,215) = b(k,215) - lu(k,917) * b(k,179)
         b(k,217) = b(k,217) - lu(k,918) * b(k,179)
         b(k,218) = b(k,218) - lu(k,919) * b(k,179)
         b(k,221) = b(k,221) - lu(k,920) * b(k,179)
         b(k,224) = b(k,224) - lu(k,921) * b(k,179)
                                                                        
         b(k,194) = b(k,194) - lu(k,923) * b(k,180)
         b(k,207) = b(k,207) - lu(k,924) * b(k,180)
         b(k,215) = b(k,215) - lu(k,925) * b(k,180)
         b(k,221) = b(k,221) - lu(k,926) * b(k,180)
         b(k,227) = b(k,227) - lu(k,927) * b(k,180)
                                                                        
         b(k,192) = b(k,192) - lu(k,930) * b(k,181)
         b(k,215) = b(k,215) - lu(k,931) * b(k,181)
         b(k,218) = b(k,218) - lu(k,932) * b(k,181)
         b(k,226) = b(k,226) - lu(k,933) * b(k,181)
         b(k,227) = b(k,227) - lu(k,934) * b(k,181)
                                                                        
         b(k,183) = b(k,183) - lu(k,949) * b(k,182)
         b(k,184) = b(k,184) - lu(k,950) * b(k,182)
         b(k,187) = b(k,187) - lu(k,951) * b(k,182)
         b(k,188) = b(k,188) - lu(k,952) * b(k,182)
         b(k,190) = b(k,190) - lu(k,953) * b(k,182)
         b(k,192) = b(k,192) - lu(k,954) * b(k,182)
         b(k,194) = b(k,194) - lu(k,955) * b(k,182)
         b(k,200) = b(k,200) - lu(k,956) * b(k,182)
         b(k,207) = b(k,207) - lu(k,957) * b(k,182)
         b(k,211) = b(k,211) - lu(k,958) * b(k,182)
         b(k,213) = b(k,213) - lu(k,959) * b(k,182)
         b(k,215) = b(k,215) - lu(k,960) * b(k,182)
         b(k,216) = b(k,216) - lu(k,961) * b(k,182)
         b(k,217) = b(k,217) - lu(k,962) * b(k,182)
         b(k,218) = b(k,218) - lu(k,963) * b(k,182)
         b(k,221) = b(k,221) - lu(k,964) * b(k,182)
         b(k,222) = b(k,222) - lu(k,965) * b(k,182)
         b(k,224) = b(k,224) - lu(k,966) * b(k,182)
         b(k,226) = b(k,226) - lu(k,967) * b(k,182)
         b(k,227) = b(k,227) - lu(k,968) * b(k,182)
                                                                        
         b(k,187) = b(k,187) - lu(k,970) * b(k,183)
         b(k,188) = b(k,188) - lu(k,971) * b(k,183)
         b(k,192) = b(k,192) - lu(k,972) * b(k,183)
         b(k,193) = b(k,193) - lu(k,973) * b(k,183)
         b(k,215) = b(k,215) - lu(k,974) * b(k,183)
         b(k,216) = b(k,216) - lu(k,975) * b(k,183)
         b(k,218) = b(k,218) - lu(k,976) * b(k,183)
                                                                        
         b(k,192) = b(k,192) - lu(k,980) * b(k,184)
         b(k,215) = b(k,215) - lu(k,981) * b(k,184)
         b(k,218) = b(k,218) - lu(k,982) * b(k,184)
         b(k,226) = b(k,226) - lu(k,983) * b(k,184)
                                                                        
         b(k,187) = b(k,187) - lu(k,1000) * b(k,185)
         b(k,188) = b(k,188) - lu(k,1001) * b(k,185)
         b(k,190) = b(k,190) - lu(k,1002) * b(k,185)
         b(k,192) = b(k,192) - lu(k,1003) * b(k,185)
         b(k,193) = b(k,193) - lu(k,1004) * b(k,185)
         b(k,194) = b(k,194) - lu(k,1005) * b(k,185)
         b(k,200) = b(k,200) - lu(k,1006) * b(k,185)
         b(k,207) = b(k,207) - lu(k,1007) * b(k,185)
         b(k,211) = b(k,211) - lu(k,1008) * b(k,185)
         b(k,213) = b(k,213) - lu(k,1009) * b(k,185)
         b(k,215) = b(k,215) - lu(k,1010) * b(k,185)
         b(k,216) = b(k,216) - lu(k,1011) * b(k,185)
         b(k,217) = b(k,217) - lu(k,1012) * b(k,185)
         b(k,218) = b(k,218) - lu(k,1013) * b(k,185)
         b(k,221) = b(k,221) - lu(k,1014) * b(k,185)
         b(k,222) = b(k,222) - lu(k,1015) * b(k,185)
         b(k,224) = b(k,224) - lu(k,1016) * b(k,185)
         b(k,226) = b(k,226) - lu(k,1017) * b(k,185)
         b(k,227) = b(k,227) - lu(k,1018) * b(k,185)
                                                                        
         b(k,190) = b(k,190) - lu(k,1025) * b(k,186)
         b(k,192) = b(k,192) - lu(k,1026) * b(k,186)
         b(k,195) = b(k,195) - lu(k,1027) * b(k,186)
         b(k,200) = b(k,200) - lu(k,1028) * b(k,186)
         b(k,207) = b(k,207) - lu(k,1029) * b(k,186)
         b(k,210) = b(k,210) - lu(k,1030) * b(k,186)
         b(k,213) = b(k,213) - lu(k,1031) * b(k,186)
         b(k,215) = b(k,215) - lu(k,1032) * b(k,186)
         b(k,216) = b(k,216) - lu(k,1033) * b(k,186)
         b(k,217) = b(k,217) - lu(k,1034) * b(k,186)
         b(k,218) = b(k,218) - lu(k,1035) * b(k,186)
         b(k,221) = b(k,221) - lu(k,1036) * b(k,186)
         b(k,222) = b(k,222) - lu(k,1037) * b(k,186)
         b(k,224) = b(k,224) - lu(k,1038) * b(k,186)
         b(k,226) = b(k,226) - lu(k,1039) * b(k,186)
         b(k,227) = b(k,227) - lu(k,1040) * b(k,186)
                                                                        
         b(k,188) = b(k,188) - lu(k,1046) * b(k,187)
         b(k,192) = b(k,192) - lu(k,1047) * b(k,187)
         b(k,193) = b(k,193) - lu(k,1048) * b(k,187)
         b(k,211) = b(k,211) - lu(k,1049) * b(k,187)
         b(k,213) = b(k,213) - lu(k,1050) * b(k,187)
         b(k,215) = b(k,215) - lu(k,1051) * b(k,187)
         b(k,216) = b(k,216) - lu(k,1052) * b(k,187)
         b(k,217) = b(k,217) - lu(k,1053) * b(k,187)
         b(k,218) = b(k,218) - lu(k,1054) * b(k,187)
         b(k,221) = b(k,221) - lu(k,1055) * b(k,187)
         b(k,224) = b(k,224) - lu(k,1056) * b(k,187)
                                                                        
         b(k,192) = b(k,192) - lu(k,1060) * b(k,188)
         b(k,194) = b(k,194) - lu(k,1061) * b(k,188)
         b(k,207) = b(k,207) - lu(k,1062) * b(k,188)
         b(k,213) = b(k,213) - lu(k,1063) * b(k,188)
         b(k,215) = b(k,215) - lu(k,1064) * b(k,188)
         b(k,218) = b(k,218) - lu(k,1065) * b(k,188)
         b(k,221) = b(k,221) - lu(k,1066) * b(k,188)
         b(k,226) = b(k,226) - lu(k,1067) * b(k,188)
         b(k,227) = b(k,227) - lu(k,1068) * b(k,188)
                                                                        
         b(k,192) = b(k,192) - lu(k,1074) * b(k,189)
         b(k,200) = b(k,200) - lu(k,1075) * b(k,189)
         b(k,207) = b(k,207) - lu(k,1076) * b(k,189)
         b(k,211) = b(k,211) - lu(k,1077) * b(k,189)
         b(k,213) = b(k,213) - lu(k,1078) * b(k,189)
         b(k,215) = b(k,215) - lu(k,1079) * b(k,189)
         b(k,217) = b(k,217) - lu(k,1080) * b(k,189)
         b(k,218) = b(k,218) - lu(k,1081) * b(k,189)
         b(k,221) = b(k,221) - lu(k,1082) * b(k,189)
         b(k,224) = b(k,224) - lu(k,1083) * b(k,189)
         b(k,226) = b(k,226) - lu(k,1084) * b(k,189)
                                                                        
         b(k,192) = b(k,192) - lu(k,1086) * b(k,190)
         b(k,207) = b(k,207) - lu(k,1087) * b(k,190)
         b(k,211) = b(k,211) - lu(k,1088) * b(k,190)
         b(k,215) = b(k,215) - lu(k,1089) * b(k,190)
         b(k,216) = b(k,216) - lu(k,1090) * b(k,190)
         b(k,218) = b(k,218) - lu(k,1091) * b(k,190)
         b(k,221) = b(k,221) - lu(k,1092) * b(k,190)
         b(k,227) = b(k,227) - lu(k,1093) * b(k,190)
                                                                        
         b(k,192) = b(k,192) - lu(k,1097) * b(k,191)
         b(k,213) = b(k,213) - lu(k,1098) * b(k,191)
         b(k,215) = b(k,215) - lu(k,1099) * b(k,191)
         b(k,218) = b(k,218) - lu(k,1100) * b(k,191)
         b(k,226) = b(k,226) - lu(k,1101) * b(k,191)
                                                                        
         b(k,215) = b(k,215) - lu(k,1104) * b(k,192)
         b(k,218) = b(k,218) - lu(k,1105) * b(k,192)
         b(k,226) = b(k,226) - lu(k,1106) * b(k,192)
                                                                        
         b(k,194) = b(k,194) - lu(k,1116) * b(k,193)
         b(k,207) = b(k,207) - lu(k,1117) * b(k,193)
         b(k,211) = b(k,211) - lu(k,1118) * b(k,193)
         b(k,213) = b(k,213) - lu(k,1119) * b(k,193)
         b(k,215) = b(k,215) - lu(k,1120) * b(k,193)
         b(k,217) = b(k,217) - lu(k,1121) * b(k,193)
         b(k,218) = b(k,218) - lu(k,1122) * b(k,193)
         b(k,221) = b(k,221) - lu(k,1123) * b(k,193)
         b(k,224) = b(k,224) - lu(k,1124) * b(k,193)
         b(k,226) = b(k,226) - lu(k,1125) * b(k,193)
         b(k,227) = b(k,227) - lu(k,1126) * b(k,193)
                                                                        
         b(k,195) = b(k,195) - lu(k,1130) * b(k,194)
         b(k,200) = b(k,200) - lu(k,1131) * b(k,194)
         b(k,207) = b(k,207) - lu(k,1132) * b(k,194)
         b(k,213) = b(k,213) - lu(k,1133) * b(k,194)
         b(k,215) = b(k,215) - lu(k,1134) * b(k,194)
         b(k,217) = b(k,217) - lu(k,1135) * b(k,194)
         b(k,218) = b(k,218) - lu(k,1136) * b(k,194)
         b(k,221) = b(k,221) - lu(k,1137) * b(k,194)
         b(k,224) = b(k,224) - lu(k,1138) * b(k,194)
         b(k,227) = b(k,227) - lu(k,1139) * b(k,194)
                                                                        
         b(k,200) = b(k,200) - lu(k,1141) * b(k,195)
         b(k,207) = b(k,207) - lu(k,1142) * b(k,195)
         b(k,213) = b(k,213) - lu(k,1143) * b(k,195)
         b(k,215) = b(k,215) - lu(k,1144) * b(k,195)
         b(k,218) = b(k,218) - lu(k,1145) * b(k,195)
                                                                        
         b(k,200) = b(k,200) - lu(k,1154) * b(k,196)
         b(k,207) = b(k,207) - lu(k,1155) * b(k,196)
         b(k,211) = b(k,211) - lu(k,1156) * b(k,196)
         b(k,213) = b(k,213) - lu(k,1157) * b(k,196)
         b(k,215) = b(k,215) - lu(k,1158) * b(k,196)
         b(k,216) = b(k,216) - lu(k,1159) * b(k,196)
         b(k,217) = b(k,217) - lu(k,1160) * b(k,196)
         b(k,218) = b(k,218) - lu(k,1161) * b(k,196)
         b(k,221) = b(k,221) - lu(k,1162) * b(k,196)
         b(k,224) = b(k,224) - lu(k,1163) * b(k,196)
         b(k,226) = b(k,226) - lu(k,1164) * b(k,196)
         b(k,227) = b(k,227) - lu(k,1165) * b(k,196)
                                                                        
         b(k,200) = b(k,200) - lu(k,1174) * b(k,197)
         b(k,207) = b(k,207) - lu(k,1175) * b(k,197)
         b(k,213) = b(k,213) - lu(k,1176) * b(k,197)
         b(k,215) = b(k,215) - lu(k,1177) * b(k,197)
         b(k,216) = b(k,216) - lu(k,1178) * b(k,197)
         b(k,217) = b(k,217) - lu(k,1179) * b(k,197)
         b(k,218) = b(k,218) - lu(k,1180) * b(k,197)
         b(k,221) = b(k,221) - lu(k,1181) * b(k,197)
         b(k,224) = b(k,224) - lu(k,1182) * b(k,197)
         b(k,226) = b(k,226) - lu(k,1183) * b(k,197)
                                                                        
         b(k,199) = b(k,199) - lu(k,1194) * b(k,198)
         b(k,200) = b(k,200) - lu(k,1195) * b(k,198)
         b(k,205) = b(k,205) - lu(k,1196) * b(k,198)
         b(k,207) = b(k,207) - lu(k,1197) * b(k,198)
         b(k,211) = b(k,211) - lu(k,1198) * b(k,198)
         b(k,213) = b(k,213) - lu(k,1199) * b(k,198)
         b(k,215) = b(k,215) - lu(k,1200) * b(k,198)
         b(k,216) = b(k,216) - lu(k,1201) * b(k,198)
         b(k,217) = b(k,217) - lu(k,1202) * b(k,198)
         b(k,218) = b(k,218) - lu(k,1203) * b(k,198)
         b(k,221) = b(k,221) - lu(k,1204) * b(k,198)
         b(k,224) = b(k,224) - lu(k,1205) * b(k,198)
         b(k,226) = b(k,226) - lu(k,1206) * b(k,198)
                                                                        
      end do
                                                                        
      end subroutine lu_slv04
                                                                        
      subroutine lu_slv05( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
         b(k,200) = b(k,200) - lu(k,1210) * b(k,199)
         b(k,204) = b(k,204) - lu(k,1211) * b(k,199)
         b(k,206) = b(k,206) - lu(k,1212) * b(k,199)
         b(k,207) = b(k,207) - lu(k,1213) * b(k,199)
         b(k,213) = b(k,213) - lu(k,1214) * b(k,199)
         b(k,215) = b(k,215) - lu(k,1215) * b(k,199)
         b(k,218) = b(k,218) - lu(k,1216) * b(k,199)
         b(k,222) = b(k,222) - lu(k,1217) * b(k,199)
         b(k,226) = b(k,226) - lu(k,1218) * b(k,199)
         b(k,227) = b(k,227) - lu(k,1219) * b(k,199)
                                                                        
         b(k,207) = b(k,207) - lu(k,1222) * b(k,200)
         b(k,211) = b(k,211) - lu(k,1223) * b(k,200)
         b(k,215) = b(k,215) - lu(k,1224) * b(k,200)
         b(k,216) = b(k,216) - lu(k,1225) * b(k,200)
         b(k,218) = b(k,218) - lu(k,1226) * b(k,200)
         b(k,226) = b(k,226) - lu(k,1227) * b(k,200)
         b(k,227) = b(k,227) - lu(k,1228) * b(k,200)
                                                                        
         b(k,209) = b(k,209) - lu(k,1233) * b(k,201)
         b(k,212) = b(k,212) - lu(k,1234) * b(k,201)
         b(k,215) = b(k,215) - lu(k,1235) * b(k,201)
         b(k,217) = b(k,217) - lu(k,1236) * b(k,201)
         b(k,218) = b(k,218) - lu(k,1237) * b(k,201)
         b(k,219) = b(k,219) - lu(k,1238) * b(k,201)
         b(k,220) = b(k,220) - lu(k,1239) * b(k,201)
         b(k,222) = b(k,222) - lu(k,1240) * b(k,201)
         b(k,224) = b(k,224) - lu(k,1241) * b(k,201)
         b(k,225) = b(k,225) - lu(k,1242) * b(k,201)
         b(k,226) = b(k,226) - lu(k,1243) * b(k,201)
         b(k,227) = b(k,227) - lu(k,1244) * b(k,201)
                                                                        
         b(k,204) = b(k,204) - lu(k,1257) * b(k,202)
         b(k,205) = b(k,205) - lu(k,1258) * b(k,202)
         b(k,206) = b(k,206) - lu(k,1259) * b(k,202)
         b(k,207) = b(k,207) - lu(k,1260) * b(k,202)
         b(k,211) = b(k,211) - lu(k,1261) * b(k,202)
         b(k,213) = b(k,213) - lu(k,1262) * b(k,202)
         b(k,215) = b(k,215) - lu(k,1263) * b(k,202)
         b(k,216) = b(k,216) - lu(k,1264) * b(k,202)
         b(k,217) = b(k,217) - lu(k,1265) * b(k,202)
         b(k,218) = b(k,218) - lu(k,1266) * b(k,202)
         b(k,221) = b(k,221) - lu(k,1267) * b(k,202)
         b(k,222) = b(k,222) - lu(k,1268) * b(k,202)
         b(k,224) = b(k,224) - lu(k,1269) * b(k,202)
         b(k,226) = b(k,226) - lu(k,1270) * b(k,202)
         b(k,227) = b(k,227) - lu(k,1271) * b(k,202)
                                                                        
         b(k,204) = b(k,204) - lu(k,1289) * b(k,203)
         b(k,205) = b(k,205) - lu(k,1290) * b(k,203)
         b(k,206) = b(k,206) - lu(k,1291) * b(k,203)
         b(k,207) = b(k,207) - lu(k,1292) * b(k,203)
         b(k,211) = b(k,211) - lu(k,1293) * b(k,203)
         b(k,213) = b(k,213) - lu(k,1294) * b(k,203)
         b(k,215) = b(k,215) - lu(k,1295) * b(k,203)
         b(k,216) = b(k,216) - lu(k,1296) * b(k,203)
         b(k,217) = b(k,217) - lu(k,1297) * b(k,203)
         b(k,218) = b(k,218) - lu(k,1298) * b(k,203)
         b(k,221) = b(k,221) - lu(k,1299) * b(k,203)
         b(k,222) = b(k,222) - lu(k,1300) * b(k,203)
         b(k,224) = b(k,224) - lu(k,1301) * b(k,203)
         b(k,226) = b(k,226) - lu(k,1302) * b(k,203)
         b(k,227) = b(k,227) - lu(k,1303) * b(k,203)
                                                                        
         b(k,206) = b(k,206) - lu(k,1312) * b(k,204)
         b(k,207) = b(k,207) - lu(k,1313) * b(k,204)
         b(k,211) = b(k,211) - lu(k,1314) * b(k,204)
         b(k,213) = b(k,213) - lu(k,1315) * b(k,204)
         b(k,215) = b(k,215) - lu(k,1316) * b(k,204)
         b(k,216) = b(k,216) - lu(k,1317) * b(k,204)
         b(k,217) = b(k,217) - lu(k,1318) * b(k,204)
         b(k,218) = b(k,218) - lu(k,1319) * b(k,204)
         b(k,221) = b(k,221) - lu(k,1320) * b(k,204)
         b(k,224) = b(k,224) - lu(k,1321) * b(k,204)
         b(k,226) = b(k,226) - lu(k,1322) * b(k,204)
         b(k,227) = b(k,227) - lu(k,1323) * b(k,204)
                                                                        
         b(k,206) = b(k,206) - lu(k,1333) * b(k,205)
         b(k,207) = b(k,207) - lu(k,1334) * b(k,205)
         b(k,210) = b(k,210) - lu(k,1335) * b(k,205)
         b(k,211) = b(k,211) - lu(k,1336) * b(k,205)
         b(k,213) = b(k,213) - lu(k,1337) * b(k,205)
         b(k,215) = b(k,215) - lu(k,1338) * b(k,205)
         b(k,216) = b(k,216) - lu(k,1339) * b(k,205)
         b(k,217) = b(k,217) - lu(k,1340) * b(k,205)
         b(k,218) = b(k,218) - lu(k,1341) * b(k,205)
         b(k,221) = b(k,221) - lu(k,1342) * b(k,205)
         b(k,222) = b(k,222) - lu(k,1343) * b(k,205)
         b(k,224) = b(k,224) - lu(k,1344) * b(k,205)
         b(k,226) = b(k,226) - lu(k,1345) * b(k,205)
         b(k,227) = b(k,227) - lu(k,1346) * b(k,205)
                                                                        
         b(k,207) = b(k,207) - lu(k,1355) * b(k,206)
         b(k,211) = b(k,211) - lu(k,1356) * b(k,206)
         b(k,213) = b(k,213) - lu(k,1357) * b(k,206)
         b(k,215) = b(k,215) - lu(k,1358) * b(k,206)
         b(k,216) = b(k,216) - lu(k,1359) * b(k,206)
         b(k,217) = b(k,217) - lu(k,1360) * b(k,206)
         b(k,218) = b(k,218) - lu(k,1361) * b(k,206)
         b(k,221) = b(k,221) - lu(k,1362) * b(k,206)
         b(k,222) = b(k,222) - lu(k,1363) * b(k,206)
         b(k,224) = b(k,224) - lu(k,1364) * b(k,206)
         b(k,226) = b(k,226) - lu(k,1365) * b(k,206)
         b(k,227) = b(k,227) - lu(k,1366) * b(k,206)
                                                                        
         b(k,210) = b(k,210) - lu(k,1387) * b(k,207)
         b(k,211) = b(k,211) - lu(k,1388) * b(k,207)
         b(k,213) = b(k,213) - lu(k,1389) * b(k,207)
         b(k,215) = b(k,215) - lu(k,1390) * b(k,207)
         b(k,216) = b(k,216) - lu(k,1391) * b(k,207)
         b(k,217) = b(k,217) - lu(k,1392) * b(k,207)
         b(k,218) = b(k,218) - lu(k,1393) * b(k,207)
         b(k,221) = b(k,221) - lu(k,1394) * b(k,207)
         b(k,222) = b(k,222) - lu(k,1395) * b(k,207)
         b(k,224) = b(k,224) - lu(k,1396) * b(k,207)
         b(k,226) = b(k,226) - lu(k,1397) * b(k,207)
         b(k,227) = b(k,227) - lu(k,1398) * b(k,207)
                                                                        
         b(k,210) = b(k,210) - lu(k,1402) * b(k,208)
         b(k,211) = b(k,211) - lu(k,1403) * b(k,208)
         b(k,212) = b(k,212) - lu(k,1404) * b(k,208)
         b(k,214) = b(k,214) - lu(k,1405) * b(k,208)
         b(k,215) = b(k,215) - lu(k,1406) * b(k,208)
         b(k,216) = b(k,216) - lu(k,1407) * b(k,208)
         b(k,220) = b(k,220) - lu(k,1408) * b(k,208)
         b(k,221) = b(k,221) - lu(k,1409) * b(k,208)
         b(k,223) = b(k,223) - lu(k,1410) * b(k,208)
         b(k,226) = b(k,226) - lu(k,1411) * b(k,208)
         b(k,227) = b(k,227) - lu(k,1412) * b(k,208)
                                                                        
         b(k,212) = b(k,212) - lu(k,1416) * b(k,209)
         b(k,213) = b(k,213) - lu(k,1417) * b(k,209)
         b(k,214) = b(k,214) - lu(k,1418) * b(k,209)
         b(k,215) = b(k,215) - lu(k,1419) * b(k,209)
         b(k,218) = b(k,218) - lu(k,1420) * b(k,209)
         b(k,222) = b(k,222) - lu(k,1421) * b(k,209)
         b(k,225) = b(k,225) - lu(k,1422) * b(k,209)
         b(k,226) = b(k,226) - lu(k,1423) * b(k,209)
         b(k,227) = b(k,227) - lu(k,1424) * b(k,209)
                                                                        
         b(k,211) = b(k,211) - lu(k,1431) * b(k,210)
         b(k,212) = b(k,212) - lu(k,1432) * b(k,210)
         b(k,213) = b(k,213) - lu(k,1433) * b(k,210)
         b(k,214) = b(k,214) - lu(k,1434) * b(k,210)
         b(k,215) = b(k,215) - lu(k,1435) * b(k,210)
         b(k,216) = b(k,216) - lu(k,1436) * b(k,210)
         b(k,218) = b(k,218) - lu(k,1437) * b(k,210)
         b(k,220) = b(k,220) - lu(k,1438) * b(k,210)
         b(k,221) = b(k,221) - lu(k,1439) * b(k,210)
         b(k,223) = b(k,223) - lu(k,1440) * b(k,210)
         b(k,226) = b(k,226) - lu(k,1441) * b(k,210)
         b(k,227) = b(k,227) - lu(k,1442) * b(k,210)
                                                                        
         b(k,212) = b(k,212) - lu(k,1448) * b(k,211)
         b(k,213) = b(k,213) - lu(k,1449) * b(k,211)
         b(k,214) = b(k,214) - lu(k,1450) * b(k,211)
         b(k,215) = b(k,215) - lu(k,1451) * b(k,211)
         b(k,216) = b(k,216) - lu(k,1452) * b(k,211)
         b(k,218) = b(k,218) - lu(k,1453) * b(k,211)
         b(k,220) = b(k,220) - lu(k,1454) * b(k,211)
         b(k,221) = b(k,221) - lu(k,1455) * b(k,211)
         b(k,223) = b(k,223) - lu(k,1456) * b(k,211)
         b(k,224) = b(k,224) - lu(k,1457) * b(k,211)
         b(k,226) = b(k,226) - lu(k,1458) * b(k,211)
         b(k,227) = b(k,227) - lu(k,1459) * b(k,211)
                                                                        
         b(k,213) = b(k,213) - lu(k,1464) * b(k,212)
         b(k,214) = b(k,214) - lu(k,1465) * b(k,212)
         b(k,215) = b(k,215) - lu(k,1466) * b(k,212)
         b(k,216) = b(k,216) - lu(k,1467) * b(k,212)
         b(k,218) = b(k,218) - lu(k,1468) * b(k,212)
         b(k,220) = b(k,220) - lu(k,1469) * b(k,212)
         b(k,221) = b(k,221) - lu(k,1470) * b(k,212)
         b(k,222) = b(k,222) - lu(k,1471) * b(k,212)
         b(k,223) = b(k,223) - lu(k,1472) * b(k,212)
         b(k,224) = b(k,224) - lu(k,1473) * b(k,212)
         b(k,226) = b(k,226) - lu(k,1474) * b(k,212)
         b(k,227) = b(k,227) - lu(k,1475) * b(k,212)
                                                                        
         b(k,214) = b(k,214) - lu(k,1486) * b(k,213)
         b(k,215) = b(k,215) - lu(k,1487) * b(k,213)
         b(k,216) = b(k,216) - lu(k,1488) * b(k,213)
         b(k,217) = b(k,217) - lu(k,1489) * b(k,213)
         b(k,218) = b(k,218) - lu(k,1490) * b(k,213)
         b(k,220) = b(k,220) - lu(k,1491) * b(k,213)
         b(k,221) = b(k,221) - lu(k,1492) * b(k,213)
         b(k,222) = b(k,222) - lu(k,1493) * b(k,213)
         b(k,223) = b(k,223) - lu(k,1494) * b(k,213)
         b(k,224) = b(k,224) - lu(k,1495) * b(k,213)
         b(k,225) = b(k,225) - lu(k,1496) * b(k,213)
         b(k,226) = b(k,226) - lu(k,1497) * b(k,213)
         b(k,227) = b(k,227) - lu(k,1498) * b(k,213)
                                                                        
         b(k,215) = b(k,215) - lu(k,1527) * b(k,214)
         b(k,216) = b(k,216) - lu(k,1528) * b(k,214)
         b(k,217) = b(k,217) - lu(k,1529) * b(k,214)
         b(k,218) = b(k,218) - lu(k,1530) * b(k,214)
         b(k,219) = b(k,219) - lu(k,1531) * b(k,214)
         b(k,220) = b(k,220) - lu(k,1532) * b(k,214)
         b(k,221) = b(k,221) - lu(k,1533) * b(k,214)
         b(k,222) = b(k,222) - lu(k,1534) * b(k,214)
         b(k,223) = b(k,223) - lu(k,1535) * b(k,214)
         b(k,224) = b(k,224) - lu(k,1536) * b(k,214)
         b(k,225) = b(k,225) - lu(k,1537) * b(k,214)
         b(k,226) = b(k,226) - lu(k,1538) * b(k,214)
         b(k,227) = b(k,227) - lu(k,1539) * b(k,214)
                                                                        
         b(k,216) = b(k,216) - lu(k,1692) * b(k,215)
         b(k,217) = b(k,217) - lu(k,1693) * b(k,215)
         b(k,218) = b(k,218) - lu(k,1694) * b(k,215)
         b(k,219) = b(k,219) - lu(k,1695) * b(k,215)
         b(k,220) = b(k,220) - lu(k,1696) * b(k,215)
         b(k,221) = b(k,221) - lu(k,1697) * b(k,215)
         b(k,222) = b(k,222) - lu(k,1698) * b(k,215)
         b(k,223) = b(k,223) - lu(k,1699) * b(k,215)
         b(k,224) = b(k,224) - lu(k,1700) * b(k,215)
         b(k,225) = b(k,225) - lu(k,1701) * b(k,215)
         b(k,226) = b(k,226) - lu(k,1702) * b(k,215)
         b(k,227) = b(k,227) - lu(k,1703) * b(k,215)
                                                                        
      end do
                                                                        
      end subroutine lu_slv05
                                                                        
      subroutine lu_slv06( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
         b(k,217) = b(k,217) - lu(k,1750) * b(k,216)
         b(k,218) = b(k,218) - lu(k,1751) * b(k,216)
         b(k,219) = b(k,219) - lu(k,1752) * b(k,216)
         b(k,220) = b(k,220) - lu(k,1753) * b(k,216)
         b(k,221) = b(k,221) - lu(k,1754) * b(k,216)
         b(k,222) = b(k,222) - lu(k,1755) * b(k,216)
         b(k,223) = b(k,223) - lu(k,1756) * b(k,216)
         b(k,224) = b(k,224) - lu(k,1757) * b(k,216)
         b(k,225) = b(k,225) - lu(k,1758) * b(k,216)
         b(k,226) = b(k,226) - lu(k,1759) * b(k,216)
         b(k,227) = b(k,227) - lu(k,1760) * b(k,216)
                                                                        
         b(k,218) = b(k,218) - lu(k,1843) * b(k,217)
         b(k,219) = b(k,219) - lu(k,1844) * b(k,217)
         b(k,220) = b(k,220) - lu(k,1845) * b(k,217)
         b(k,221) = b(k,221) - lu(k,1846) * b(k,217)
         b(k,222) = b(k,222) - lu(k,1847) * b(k,217)
         b(k,223) = b(k,223) - lu(k,1848) * b(k,217)
         b(k,224) = b(k,224) - lu(k,1849) * b(k,217)
         b(k,225) = b(k,225) - lu(k,1850) * b(k,217)
         b(k,226) = b(k,226) - lu(k,1851) * b(k,217)
         b(k,227) = b(k,227) - lu(k,1852) * b(k,217)
                                                                        
         b(k,219) = b(k,219) - lu(k,1951) * b(k,218)
         b(k,220) = b(k,220) - lu(k,1952) * b(k,218)
         b(k,221) = b(k,221) - lu(k,1953) * b(k,218)
         b(k,222) = b(k,222) - lu(k,1954) * b(k,218)
         b(k,223) = b(k,223) - lu(k,1955) * b(k,218)
         b(k,224) = b(k,224) - lu(k,1956) * b(k,218)
         b(k,225) = b(k,225) - lu(k,1957) * b(k,218)
         b(k,226) = b(k,226) - lu(k,1958) * b(k,218)
         b(k,227) = b(k,227) - lu(k,1959) * b(k,218)
                                                                        
         b(k,220) = b(k,220) - lu(k,1978) * b(k,219)
         b(k,221) = b(k,221) - lu(k,1979) * b(k,219)
         b(k,222) = b(k,222) - lu(k,1980) * b(k,219)
         b(k,223) = b(k,223) - lu(k,1981) * b(k,219)
         b(k,224) = b(k,224) - lu(k,1982) * b(k,219)
         b(k,225) = b(k,225) - lu(k,1983) * b(k,219)
         b(k,226) = b(k,226) - lu(k,1984) * b(k,219)
         b(k,227) = b(k,227) - lu(k,1985) * b(k,219)
                                                                        
         b(k,221) = b(k,221) - lu(k,2018) * b(k,220)
         b(k,222) = b(k,222) - lu(k,2019) * b(k,220)
         b(k,223) = b(k,223) - lu(k,2020) * b(k,220)
         b(k,224) = b(k,224) - lu(k,2021) * b(k,220)
         b(k,225) = b(k,225) - lu(k,2022) * b(k,220)
         b(k,226) = b(k,226) - lu(k,2023) * b(k,220)
         b(k,227) = b(k,227) - lu(k,2024) * b(k,220)
                                                                        
         b(k,222) = b(k,222) - lu(k,2071) * b(k,221)
         b(k,223) = b(k,223) - lu(k,2072) * b(k,221)
         b(k,224) = b(k,224) - lu(k,2073) * b(k,221)
         b(k,225) = b(k,225) - lu(k,2074) * b(k,221)
         b(k,226) = b(k,226) - lu(k,2075) * b(k,221)
         b(k,227) = b(k,227) - lu(k,2076) * b(k,221)
                                                                        
         b(k,223) = b(k,223) - lu(k,2133) * b(k,222)
         b(k,224) = b(k,224) - lu(k,2134) * b(k,222)
         b(k,225) = b(k,225) - lu(k,2135) * b(k,222)
         b(k,226) = b(k,226) - lu(k,2136) * b(k,222)
         b(k,227) = b(k,227) - lu(k,2137) * b(k,222)
                                                                        
         b(k,224) = b(k,224) - lu(k,2157) * b(k,223)
         b(k,225) = b(k,225) - lu(k,2158) * b(k,223)
         b(k,226) = b(k,226) - lu(k,2159) * b(k,223)
         b(k,227) = b(k,227) - lu(k,2160) * b(k,223)
                                                                        
         b(k,225) = b(k,225) - lu(k,2202) * b(k,224)
         b(k,226) = b(k,226) - lu(k,2203) * b(k,224)
         b(k,227) = b(k,227) - lu(k,2204) * b(k,224)
                                                                        
         b(k,226) = b(k,226) - lu(k,2227) * b(k,225)
         b(k,227) = b(k,227) - lu(k,2228) * b(k,225)
                                                                        
         b(k,227) = b(k,227) - lu(k,2259) * b(k,226)
                                                                        
      end do
                                                                        
      end subroutine lu_slv06
                                                                        
      subroutine lu_slv07( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
                                                                        
!-----------------------------------------------------------------------
!       ... Solve U * x = y
!-----------------------------------------------------------------------
         b(k,227) = b(k,227) * lu(k,2285)
         b(k,226) = b(k,226) - lu(k,2284) * b(k,227)
         b(k,225) = b(k,225) - lu(k,2283) * b(k,227)
         b(k,224) = b(k,224) - lu(k,2282) * b(k,227)
         b(k,223) = b(k,223) - lu(k,2281) * b(k,227)
         b(k,222) = b(k,222) - lu(k,2280) * b(k,227)
         b(k,221) = b(k,221) - lu(k,2279) * b(k,227)
         b(k,220) = b(k,220) - lu(k,2278) * b(k,227)
         b(k,219) = b(k,219) - lu(k,2277) * b(k,227)
         b(k,218) = b(k,218) - lu(k,2276) * b(k,227)
         b(k,217) = b(k,217) - lu(k,2275) * b(k,227)
         b(k,216) = b(k,216) - lu(k,2274) * b(k,227)
         b(k,215) = b(k,215) - lu(k,2273) * b(k,227)
         b(k,214) = b(k,214) - lu(k,2272) * b(k,227)
         b(k,213) = b(k,213) - lu(k,2271) * b(k,227)
         b(k,212) = b(k,212) - lu(k,2270) * b(k,227)
         b(k,211) = b(k,211) - lu(k,2269) * b(k,227)
         b(k,210) = b(k,210) - lu(k,2268) * b(k,227)
         b(k,209) = b(k,209) - lu(k,2267) * b(k,227)
         b(k,208) = b(k,208) - lu(k,2266) * b(k,227)
         b(k,201) = b(k,201) - lu(k,2265) * b(k,227)
         b(k,175) = b(k,175) - lu(k,2264) * b(k,227)
         b(k,172) = b(k,172) - lu(k,2263) * b(k,227)
         b(k,97) = b(k,97) - lu(k,2262) * b(k,227)
         b(k,91) = b(k,91) - lu(k,2261) * b(k,227)
         b(k,63) = b(k,63) - lu(k,2260) * b(k,227)
                                                                        
         b(k,226) = b(k,226) * lu(k,2258)
         b(k,225) = b(k,225) - lu(k,2257) * b(k,226)
         b(k,224) = b(k,224) - lu(k,2256) * b(k,226)
         b(k,223) = b(k,223) - lu(k,2255) * b(k,226)
         b(k,222) = b(k,222) - lu(k,2254) * b(k,226)
         b(k,221) = b(k,221) - lu(k,2253) * b(k,226)
         b(k,220) = b(k,220) - lu(k,2252) * b(k,226)
         b(k,219) = b(k,219) - lu(k,2251) * b(k,226)
         b(k,218) = b(k,218) - lu(k,2250) * b(k,226)
         b(k,217) = b(k,217) - lu(k,2249) * b(k,226)
         b(k,216) = b(k,216) - lu(k,2248) * b(k,226)
         b(k,215) = b(k,215) - lu(k,2247) * b(k,226)
         b(k,214) = b(k,214) - lu(k,2246) * b(k,226)
         b(k,213) = b(k,213) - lu(k,2245) * b(k,226)
         b(k,212) = b(k,212) - lu(k,2244) * b(k,226)
         b(k,211) = b(k,211) - lu(k,2243) * b(k,226)
         b(k,210) = b(k,210) - lu(k,2242) * b(k,226)
         b(k,209) = b(k,209) - lu(k,2241) * b(k,226)
         b(k,208) = b(k,208) - lu(k,2240) * b(k,226)
         b(k,201) = b(k,201) - lu(k,2239) * b(k,226)
         b(k,192) = b(k,192) - lu(k,2238) * b(k,226)
         b(k,177) = b(k,177) - lu(k,2237) * b(k,226)
         b(k,172) = b(k,172) - lu(k,2236) * b(k,226)
         b(k,171) = b(k,171) - lu(k,2235) * b(k,226)
         b(k,169) = b(k,169) - lu(k,2234) * b(k,226)
         b(k,165) = b(k,165) - lu(k,2233) * b(k,226)
         b(k,147) = b(k,147) - lu(k,2232) * b(k,226)
         b(k,139) = b(k,139) - lu(k,2231) * b(k,226)
         b(k,135) = b(k,135) - lu(k,2230) * b(k,226)
         b(k,111) = b(k,111) - lu(k,2229) * b(k,226)
                                                                        
         b(k,225) = b(k,225) * lu(k,2226)
         b(k,224) = b(k,224) - lu(k,2225) * b(k,225)
         b(k,223) = b(k,223) - lu(k,2224) * b(k,225)
         b(k,222) = b(k,222) - lu(k,2223) * b(k,225)
         b(k,221) = b(k,221) - lu(k,2222) * b(k,225)
         b(k,220) = b(k,220) - lu(k,2221) * b(k,225)
         b(k,219) = b(k,219) - lu(k,2220) * b(k,225)
         b(k,218) = b(k,218) - lu(k,2219) * b(k,225)
         b(k,217) = b(k,217) - lu(k,2218) * b(k,225)
         b(k,216) = b(k,216) - lu(k,2217) * b(k,225)
         b(k,215) = b(k,215) - lu(k,2216) * b(k,225)
         b(k,214) = b(k,214) - lu(k,2215) * b(k,225)
         b(k,213) = b(k,213) - lu(k,2214) * b(k,225)
         b(k,212) = b(k,212) - lu(k,2213) * b(k,225)
         b(k,211) = b(k,211) - lu(k,2212) * b(k,225)
         b(k,209) = b(k,209) - lu(k,2211) * b(k,225)
         b(k,201) = b(k,201) - lu(k,2210) * b(k,225)
         b(k,172) = b(k,172) - lu(k,2209) * b(k,225)
         b(k,165) = b(k,165) - lu(k,2208) * b(k,225)
         b(k,139) = b(k,139) - lu(k,2207) * b(k,225)
         b(k,106) = b(k,106) - lu(k,2206) * b(k,225)
         b(k,87) = b(k,87) - lu(k,2205) * b(k,225)
                                                                        
         b(k,224) = b(k,224) * lu(k,2201)
         b(k,223) = b(k,223) - lu(k,2200) * b(k,224)
         b(k,222) = b(k,222) - lu(k,2199) * b(k,224)
         b(k,221) = b(k,221) - lu(k,2198) * b(k,224)
         b(k,220) = b(k,220) - lu(k,2197) * b(k,224)
         b(k,219) = b(k,219) - lu(k,2196) * b(k,224)
         b(k,218) = b(k,218) - lu(k,2195) * b(k,224)
         b(k,217) = b(k,217) - lu(k,2194) * b(k,224)
         b(k,216) = b(k,216) - lu(k,2193) * b(k,224)
         b(k,215) = b(k,215) - lu(k,2192) * b(k,224)
         b(k,214) = b(k,214) - lu(k,2191) * b(k,224)
         b(k,213) = b(k,213) - lu(k,2190) * b(k,224)
         b(k,212) = b(k,212) - lu(k,2189) * b(k,224)
         b(k,211) = b(k,211) - lu(k,2188) * b(k,224)
         b(k,210) = b(k,210) - lu(k,2187) * b(k,224)
         b(k,209) = b(k,209) - lu(k,2186) * b(k,224)
         b(k,207) = b(k,207) - lu(k,2185) * b(k,224)
         b(k,206) = b(k,206) - lu(k,2184) * b(k,224)
         b(k,201) = b(k,201) - lu(k,2183) * b(k,224)
         b(k,200) = b(k,200) - lu(k,2182) * b(k,224)
         b(k,195) = b(k,195) - lu(k,2181) * b(k,224)
         b(k,192) = b(k,192) - lu(k,2180) * b(k,224)
         b(k,184) = b(k,184) - lu(k,2179) * b(k,224)
         b(k,177) = b(k,177) - lu(k,2178) * b(k,224)
         b(k,176) = b(k,176) - lu(k,2177) * b(k,224)
         b(k,172) = b(k,172) - lu(k,2176) * b(k,224)
         b(k,168) = b(k,168) - lu(k,2175) * b(k,224)
         b(k,165) = b(k,165) - lu(k,2174) * b(k,224)
         b(k,163) = b(k,163) - lu(k,2173) * b(k,224)
         b(k,160) = b(k,160) - lu(k,2172) * b(k,224)
         b(k,150) = b(k,150) - lu(k,2171) * b(k,224)
         b(k,143) = b(k,143) - lu(k,2170) * b(k,224)
         b(k,139) = b(k,139) - lu(k,2169) * b(k,224)
         b(k,137) = b(k,137) - lu(k,2168) * b(k,224)
         b(k,136) = b(k,136) - lu(k,2167) * b(k,224)
         b(k,132) = b(k,132) - lu(k,2166) * b(k,224)
         b(k,128) = b(k,128) - lu(k,2165) * b(k,224)
         b(k,118) = b(k,118) - lu(k,2164) * b(k,224)
         b(k,95) = b(k,95) - lu(k,2163) * b(k,224)
         b(k,75) = b(k,75) - lu(k,2162) * b(k,224)
         b(k,65) = b(k,65) - lu(k,2161) * b(k,224)
                                                                        
         b(k,223) = b(k,223) * lu(k,2156)
         b(k,222) = b(k,222) - lu(k,2155) * b(k,223)
         b(k,221) = b(k,221) - lu(k,2154) * b(k,223)
         b(k,220) = b(k,220) - lu(k,2153) * b(k,223)
         b(k,219) = b(k,219) - lu(k,2152) * b(k,223)
         b(k,218) = b(k,218) - lu(k,2151) * b(k,223)
         b(k,217) = b(k,217) - lu(k,2150) * b(k,223)
         b(k,216) = b(k,216) - lu(k,2149) * b(k,223)
         b(k,215) = b(k,215) - lu(k,2148) * b(k,223)
         b(k,214) = b(k,214) - lu(k,2147) * b(k,223)
         b(k,213) = b(k,213) - lu(k,2146) * b(k,223)
         b(k,212) = b(k,212) - lu(k,2145) * b(k,223)
         b(k,211) = b(k,211) - lu(k,2144) * b(k,223)
         b(k,209) = b(k,209) - lu(k,2143) * b(k,223)
         b(k,177) = b(k,177) - lu(k,2142) * b(k,223)
         b(k,171) = b(k,171) - lu(k,2141) * b(k,223)
         b(k,165) = b(k,165) - lu(k,2140) * b(k,223)
         b(k,87) = b(k,87) - lu(k,2139) * b(k,223)
         b(k,70) = b(k,70) - lu(k,2138) * b(k,223)
                                                                        
         b(k,222) = b(k,222) * lu(k,2132)
         b(k,221) = b(k,221) - lu(k,2131) * b(k,222)
         b(k,220) = b(k,220) - lu(k,2130) * b(k,222)
         b(k,219) = b(k,219) - lu(k,2129) * b(k,222)
         b(k,218) = b(k,218) - lu(k,2128) * b(k,222)
         b(k,217) = b(k,217) - lu(k,2127) * b(k,222)
         b(k,216) = b(k,216) - lu(k,2126) * b(k,222)
         b(k,215) = b(k,215) - lu(k,2125) * b(k,222)
         b(k,214) = b(k,214) - lu(k,2124) * b(k,222)
         b(k,213) = b(k,213) - lu(k,2123) * b(k,222)
         b(k,212) = b(k,212) - lu(k,2122) * b(k,222)
         b(k,211) = b(k,211) - lu(k,2121) * b(k,222)
         b(k,210) = b(k,210) - lu(k,2120) * b(k,222)
         b(k,209) = b(k,209) - lu(k,2119) * b(k,222)
         b(k,207) = b(k,207) - lu(k,2118) * b(k,222)
         b(k,206) = b(k,206) - lu(k,2117) * b(k,222)
         b(k,205) = b(k,205) - lu(k,2116) * b(k,222)
         b(k,204) = b(k,204) - lu(k,2115) * b(k,222)
         b(k,203) = b(k,203) - lu(k,2114) * b(k,222)
         b(k,202) = b(k,202) - lu(k,2113) * b(k,222)
         b(k,201) = b(k,201) - lu(k,2112) * b(k,222)
         b(k,200) = b(k,200) - lu(k,2111) * b(k,222)
         b(k,199) = b(k,199) - lu(k,2110) * b(k,222)
         b(k,198) = b(k,198) - lu(k,2109) * b(k,222)
         b(k,195) = b(k,195) - lu(k,2108) * b(k,222)
         b(k,194) = b(k,194) - lu(k,2107) * b(k,222)
         b(k,193) = b(k,193) - lu(k,2106) * b(k,222)
         b(k,192) = b(k,192) - lu(k,2105) * b(k,222)
         b(k,191) = b(k,191) - lu(k,2104) * b(k,222)
         b(k,190) = b(k,190) - lu(k,2103) * b(k,222)
         b(k,188) = b(k,188) - lu(k,2102) * b(k,222)
         b(k,187) = b(k,187) - lu(k,2101) * b(k,222)
         b(k,186) = b(k,186) - lu(k,2100) * b(k,222)
         b(k,185) = b(k,185) - lu(k,2099) * b(k,222)
         b(k,184) = b(k,184) - lu(k,2098) * b(k,222)
         b(k,183) = b(k,183) - lu(k,2097) * b(k,222)
         b(k,182) = b(k,182) - lu(k,2096) * b(k,222)
         b(k,181) = b(k,181) - lu(k,2095) * b(k,222)
         b(k,180) = b(k,180) - lu(k,2094) * b(k,222)
         b(k,178) = b(k,178) - lu(k,2093) * b(k,222)
         b(k,173) = b(k,173) - lu(k,2092) * b(k,222)
         b(k,172) = b(k,172) - lu(k,2091) * b(k,222)
         b(k,168) = b(k,168) - lu(k,2090) * b(k,222)
         b(k,159) = b(k,159) - lu(k,2089) * b(k,222)
         b(k,156) = b(k,156) - lu(k,2088) * b(k,222)
         b(k,150) = b(k,150) - lu(k,2087) * b(k,222)
         b(k,140) = b(k,140) - lu(k,2086) * b(k,222)
         b(k,135) = b(k,135) - lu(k,2085) * b(k,222)
         b(k,128) = b(k,128) - lu(k,2084) * b(k,222)
         b(k,114) = b(k,114) - lu(k,2083) * b(k,222)
         b(k,86) = b(k,86) - lu(k,2082) * b(k,222)
         b(k,41) = b(k,41) - lu(k,2081) * b(k,222)
         b(k,40) = b(k,40) - lu(k,2080) * b(k,222)
         b(k,39) = b(k,39) - lu(k,2079) * b(k,222)
         b(k,38) = b(k,38) - lu(k,2078) * b(k,222)
         b(k,37) = b(k,37) - lu(k,2077) * b(k,222)
                                                                        
         b(k,221) = b(k,221) * lu(k,2070)
         b(k,220) = b(k,220) - lu(k,2069) * b(k,221)
         b(k,219) = b(k,219) - lu(k,2068) * b(k,221)
         b(k,218) = b(k,218) - lu(k,2067) * b(k,221)
         b(k,217) = b(k,217) - lu(k,2066) * b(k,221)
         b(k,216) = b(k,216) - lu(k,2065) * b(k,221)
         b(k,215) = b(k,215) - lu(k,2064) * b(k,221)
         b(k,214) = b(k,214) - lu(k,2063) * b(k,221)
         b(k,213) = b(k,213) - lu(k,2062) * b(k,221)
         b(k,212) = b(k,212) - lu(k,2061) * b(k,221)
         b(k,211) = b(k,211) - lu(k,2060) * b(k,221)
         b(k,210) = b(k,210) - lu(k,2059) * b(k,221)
         b(k,207) = b(k,207) - lu(k,2058) * b(k,221)
         b(k,206) = b(k,206) - lu(k,2057) * b(k,221)
         b(k,205) = b(k,205) - lu(k,2056) * b(k,221)
         b(k,204) = b(k,204) - lu(k,2055) * b(k,221)
         b(k,203) = b(k,203) - lu(k,2054) * b(k,221)
         b(k,202) = b(k,202) - lu(k,2053) * b(k,221)
         b(k,200) = b(k,200) - lu(k,2052) * b(k,221)
         b(k,199) = b(k,199) - lu(k,2051) * b(k,221)
         b(k,198) = b(k,198) - lu(k,2050) * b(k,221)
         b(k,197) = b(k,197) - lu(k,2049) * b(k,221)
         b(k,195) = b(k,195) - lu(k,2048) * b(k,221)
         b(k,194) = b(k,194) - lu(k,2047) * b(k,221)
         b(k,193) = b(k,193) - lu(k,2046) * b(k,221)
         b(k,192) = b(k,192) - lu(k,2045) * b(k,221)
         b(k,191) = b(k,191) - lu(k,2044) * b(k,221)
         b(k,190) = b(k,190) - lu(k,2043) * b(k,221)
         b(k,189) = b(k,189) - lu(k,2042) * b(k,221)
         b(k,188) = b(k,188) - lu(k,2041) * b(k,221)
         b(k,187) = b(k,187) - lu(k,2040) * b(k,221)
         b(k,184) = b(k,184) - lu(k,2039) * b(k,221)
         b(k,183) = b(k,183) - lu(k,2038) * b(k,221)
         b(k,180) = b(k,180) - lu(k,2037) * b(k,221)
         b(k,179) = b(k,179) - lu(k,2036) * b(k,221)
         b(k,174) = b(k,174) - lu(k,2035) * b(k,221)
         b(k,170) = b(k,170) - lu(k,2034) * b(k,221)
         b(k,168) = b(k,168) - lu(k,2033) * b(k,221)
         b(k,167) = b(k,167) - lu(k,2032) * b(k,221)
         b(k,166) = b(k,166) - lu(k,2031) * b(k,221)
         b(k,156) = b(k,156) - lu(k,2030) * b(k,221)
         b(k,149) = b(k,149) - lu(k,2029) * b(k,221)
         b(k,115) = b(k,115) - lu(k,2028) * b(k,221)
         b(k,113) = b(k,113) - lu(k,2027) * b(k,221)
         b(k,104) = b(k,104) - lu(k,2026) * b(k,221)
         b(k,92) = b(k,92) - lu(k,2025) * b(k,221)
                                                                        
      end do
                                                                        
      end subroutine lu_slv07
                                                                        
      subroutine lu_slv08( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
         b(k,220) = b(k,220) * lu(k,2017)
         b(k,219) = b(k,219) - lu(k,2016) * b(k,220)
         b(k,218) = b(k,218) - lu(k,2015) * b(k,220)
         b(k,217) = b(k,217) - lu(k,2014) * b(k,220)
         b(k,216) = b(k,216) - lu(k,2013) * b(k,220)
         b(k,215) = b(k,215) - lu(k,2012) * b(k,220)
         b(k,214) = b(k,214) - lu(k,2011) * b(k,220)
         b(k,213) = b(k,213) - lu(k,2010) * b(k,220)
         b(k,212) = b(k,212) - lu(k,2009) * b(k,220)
         b(k,211) = b(k,211) - lu(k,2008) * b(k,220)
         b(k,210) = b(k,210) - lu(k,2007) * b(k,220)
         b(k,209) = b(k,209) - lu(k,2006) * b(k,220)
         b(k,208) = b(k,208) - lu(k,2005) * b(k,220)
         b(k,207) = b(k,207) - lu(k,2004) * b(k,220)
         b(k,192) = b(k,192) - lu(k,2003) * b(k,220)
         b(k,191) = b(k,191) - lu(k,2002) * b(k,220)
         b(k,190) = b(k,190) - lu(k,2001) * b(k,220)
         b(k,184) = b(k,184) - lu(k,2000) * b(k,220)
         b(k,181) = b(k,181) - lu(k,1999) * b(k,220)
         b(k,177) = b(k,177) - lu(k,1998) * b(k,220)
         b(k,171) = b(k,171) - lu(k,1997) * b(k,220)
         b(k,170) = b(k,170) - lu(k,1996) * b(k,220)
         b(k,159) = b(k,159) - lu(k,1995) * b(k,220)
         b(k,147) = b(k,147) - lu(k,1994) * b(k,220)
         b(k,146) = b(k,146) - lu(k,1993) * b(k,220)
         b(k,140) = b(k,140) - lu(k,1992) * b(k,220)
         b(k,129) = b(k,129) - lu(k,1991) * b(k,220)
         b(k,125) = b(k,125) - lu(k,1990) * b(k,220)
         b(k,112) = b(k,112) - lu(k,1989) * b(k,220)
         b(k,99) = b(k,99) - lu(k,1988) * b(k,220)
         b(k,98) = b(k,98) - lu(k,1987) * b(k,220)
         b(k,70) = b(k,70) - lu(k,1986) * b(k,220)
                                                                        
         b(k,219) = b(k,219) * lu(k,1977)
         b(k,218) = b(k,218) - lu(k,1976) * b(k,219)
         b(k,217) = b(k,217) - lu(k,1975) * b(k,219)
         b(k,216) = b(k,216) - lu(k,1974) * b(k,219)
         b(k,215) = b(k,215) - lu(k,1973) * b(k,219)
         b(k,214) = b(k,214) - lu(k,1972) * b(k,219)
         b(k,213) = b(k,213) - lu(k,1971) * b(k,219)
         b(k,212) = b(k,212) - lu(k,1970) * b(k,219)
         b(k,211) = b(k,211) - lu(k,1969) * b(k,219)
         b(k,209) = b(k,209) - lu(k,1968) * b(k,219)
         b(k,201) = b(k,201) - lu(k,1967) * b(k,219)
         b(k,177) = b(k,177) - lu(k,1966) * b(k,219)
         b(k,172) = b(k,172) - lu(k,1965) * b(k,219)
         b(k,171) = b(k,171) - lu(k,1964) * b(k,219)
         b(k,106) = b(k,106) - lu(k,1963) * b(k,219)
         b(k,87) = b(k,87) - lu(k,1962) * b(k,219)
         b(k,70) = b(k,70) - lu(k,1961) * b(k,219)
         b(k,52) = b(k,52) - lu(k,1960) * b(k,219)
                                                                        
         b(k,218) = b(k,218) * lu(k,1950)
         b(k,217) = b(k,217) - lu(k,1949) * b(k,218)
         b(k,216) = b(k,216) - lu(k,1948) * b(k,218)
         b(k,215) = b(k,215) - lu(k,1947) * b(k,218)
         b(k,214) = b(k,214) - lu(k,1946) * b(k,218)
         b(k,213) = b(k,213) - lu(k,1945) * b(k,218)
         b(k,212) = b(k,212) - lu(k,1944) * b(k,218)
         b(k,211) = b(k,211) - lu(k,1943) * b(k,218)
         b(k,210) = b(k,210) - lu(k,1942) * b(k,218)
         b(k,209) = b(k,209) - lu(k,1941) * b(k,218)
         b(k,208) = b(k,208) - lu(k,1940) * b(k,218)
         b(k,207) = b(k,207) - lu(k,1939) * b(k,218)
         b(k,206) = b(k,206) - lu(k,1938) * b(k,218)
         b(k,205) = b(k,205) - lu(k,1937) * b(k,218)
         b(k,204) = b(k,204) - lu(k,1936) * b(k,218)
         b(k,203) = b(k,203) - lu(k,1935) * b(k,218)
         b(k,202) = b(k,202) - lu(k,1934) * b(k,218)
         b(k,200) = b(k,200) - lu(k,1933) * b(k,218)
         b(k,199) = b(k,199) - lu(k,1932) * b(k,218)
         b(k,198) = b(k,198) - lu(k,1931) * b(k,218)
         b(k,197) = b(k,197) - lu(k,1930) * b(k,218)
         b(k,195) = b(k,195) - lu(k,1929) * b(k,218)
         b(k,194) = b(k,194) - lu(k,1928) * b(k,218)
         b(k,193) = b(k,193) - lu(k,1927) * b(k,218)
         b(k,192) = b(k,192) - lu(k,1926) * b(k,218)
         b(k,191) = b(k,191) - lu(k,1925) * b(k,218)
         b(k,190) = b(k,190) - lu(k,1924) * b(k,218)
         b(k,188) = b(k,188) - lu(k,1923) * b(k,218)
         b(k,187) = b(k,187) - lu(k,1922) * b(k,218)
         b(k,184) = b(k,184) - lu(k,1921) * b(k,218)
         b(k,183) = b(k,183) - lu(k,1920) * b(k,218)
         b(k,181) = b(k,181) - lu(k,1919) * b(k,218)
         b(k,180) = b(k,180) - lu(k,1918) * b(k,218)
         b(k,179) = b(k,179) - lu(k,1917) * b(k,218)
         b(k,178) = b(k,178) - lu(k,1916) * b(k,218)
         b(k,176) = b(k,176) - lu(k,1915) * b(k,218)
         b(k,174) = b(k,174) - lu(k,1914) * b(k,218)
         b(k,171) = b(k,171) - lu(k,1913) * b(k,218)
         b(k,170) = b(k,170) - lu(k,1912) * b(k,218)
         b(k,169) = b(k,169) - lu(k,1911) * b(k,218)
         b(k,168) = b(k,168) - lu(k,1910) * b(k,218)
         b(k,167) = b(k,167) - lu(k,1909) * b(k,218)
         b(k,165) = b(k,165) - lu(k,1908) * b(k,218)
         b(k,164) = b(k,164) - lu(k,1907) * b(k,218)
         b(k,163) = b(k,163) - lu(k,1906) * b(k,218)
         b(k,162) = b(k,162) - lu(k,1905) * b(k,218)
         b(k,161) = b(k,161) - lu(k,1904) * b(k,218)
         b(k,160) = b(k,160) - lu(k,1903) * b(k,218)
         b(k,159) = b(k,159) - lu(k,1902) * b(k,218)
         b(k,158) = b(k,158) - lu(k,1901) * b(k,218)
         b(k,157) = b(k,157) - lu(k,1900) * b(k,218)
         b(k,156) = b(k,156) - lu(k,1899) * b(k,218)
         b(k,155) = b(k,155) - lu(k,1898) * b(k,218)
         b(k,153) = b(k,153) - lu(k,1897) * b(k,218)
         b(k,152) = b(k,152) - lu(k,1896) * b(k,218)
         b(k,151) = b(k,151) - lu(k,1895) * b(k,218)
         b(k,150) = b(k,150) - lu(k,1894) * b(k,218)
         b(k,148) = b(k,148) - lu(k,1893) * b(k,218)
         b(k,147) = b(k,147) - lu(k,1892) * b(k,218)
         b(k,138) = b(k,138) - lu(k,1891) * b(k,218)
         b(k,136) = b(k,136) - lu(k,1890) * b(k,218)
         b(k,133) = b(k,133) - lu(k,1889) * b(k,218)
         b(k,131) = b(k,131) - lu(k,1888) * b(k,218)
         b(k,130) = b(k,130) - lu(k,1887) * b(k,218)
         b(k,128) = b(k,128) - lu(k,1886) * b(k,218)
         b(k,127) = b(k,127) - lu(k,1885) * b(k,218)
         b(k,126) = b(k,126) - lu(k,1884) * b(k,218)
         b(k,124) = b(k,124) - lu(k,1883) * b(k,218)
         b(k,123) = b(k,123) - lu(k,1882) * b(k,218)
         b(k,122) = b(k,122) - lu(k,1881) * b(k,218)
         b(k,121) = b(k,121) - lu(k,1880) * b(k,218)
         b(k,120) = b(k,120) - lu(k,1879) * b(k,218)
         b(k,119) = b(k,119) - lu(k,1878) * b(k,218)
         b(k,118) = b(k,118) - lu(k,1877) * b(k,218)
         b(k,117) = b(k,117) - lu(k,1876) * b(k,218)
         b(k,116) = b(k,116) - lu(k,1875) * b(k,218)
         b(k,115) = b(k,115) - lu(k,1874) * b(k,218)
         b(k,109) = b(k,109) - lu(k,1873) * b(k,218)
         b(k,108) = b(k,108) - lu(k,1872) * b(k,218)
         b(k,107) = b(k,107) - lu(k,1871) * b(k,218)
         b(k,103) = b(k,103) - lu(k,1870) * b(k,218)
         b(k,102) = b(k,102) - lu(k,1869) * b(k,218)
         b(k,94) = b(k,94) - lu(k,1868) * b(k,218)
         b(k,93) = b(k,93) - lu(k,1867) * b(k,218)
         b(k,79) = b(k,79) - lu(k,1866) * b(k,218)
         b(k,62) = b(k,62) - lu(k,1865) * b(k,218)
         b(k,51) = b(k,51) - lu(k,1864) * b(k,218)
         b(k,50) = b(k,50) - lu(k,1863) * b(k,218)
         b(k,49) = b(k,49) - lu(k,1862) * b(k,218)
         b(k,47) = b(k,47) - lu(k,1861) * b(k,218)
         b(k,46) = b(k,46) - lu(k,1860) * b(k,218)
         b(k,45) = b(k,45) - lu(k,1859) * b(k,218)
         b(k,44) = b(k,44) - lu(k,1858) * b(k,218)
         b(k,41) = b(k,41) - lu(k,1857) * b(k,218)
         b(k,40) = b(k,40) - lu(k,1856) * b(k,218)
         b(k,39) = b(k,39) - lu(k,1855) * b(k,218)
         b(k,38) = b(k,38) - lu(k,1854) * b(k,218)
         b(k,37) = b(k,37) - lu(k,1853) * b(k,218)
                                                                        
         b(k,217) = b(k,217) * lu(k,1842)
         b(k,216) = b(k,216) - lu(k,1841) * b(k,217)
         b(k,215) = b(k,215) - lu(k,1840) * b(k,217)
         b(k,214) = b(k,214) - lu(k,1839) * b(k,217)
         b(k,213) = b(k,213) - lu(k,1838) * b(k,217)
         b(k,212) = b(k,212) - lu(k,1837) * b(k,217)
         b(k,211) = b(k,211) - lu(k,1836) * b(k,217)
         b(k,210) = b(k,210) - lu(k,1835) * b(k,217)
         b(k,209) = b(k,209) - lu(k,1834) * b(k,217)
         b(k,207) = b(k,207) - lu(k,1833) * b(k,217)
         b(k,206) = b(k,206) - lu(k,1832) * b(k,217)
         b(k,205) = b(k,205) - lu(k,1831) * b(k,217)
         b(k,204) = b(k,204) - lu(k,1830) * b(k,217)
         b(k,203) = b(k,203) - lu(k,1829) * b(k,217)
         b(k,202) = b(k,202) - lu(k,1828) * b(k,217)
         b(k,200) = b(k,200) - lu(k,1827) * b(k,217)
         b(k,199) = b(k,199) - lu(k,1826) * b(k,217)
         b(k,198) = b(k,198) - lu(k,1825) * b(k,217)
         b(k,197) = b(k,197) - lu(k,1824) * b(k,217)
         b(k,196) = b(k,196) - lu(k,1823) * b(k,217)
         b(k,195) = b(k,195) - lu(k,1822) * b(k,217)
         b(k,194) = b(k,194) - lu(k,1821) * b(k,217)
         b(k,193) = b(k,193) - lu(k,1820) * b(k,217)
         b(k,192) = b(k,192) - lu(k,1819) * b(k,217)
         b(k,191) = b(k,191) - lu(k,1818) * b(k,217)
         b(k,190) = b(k,190) - lu(k,1817) * b(k,217)
         b(k,189) = b(k,189) - lu(k,1816) * b(k,217)
         b(k,188) = b(k,188) - lu(k,1815) * b(k,217)
         b(k,187) = b(k,187) - lu(k,1814) * b(k,217)
         b(k,184) = b(k,184) - lu(k,1813) * b(k,217)
         b(k,183) = b(k,183) - lu(k,1812) * b(k,217)
         b(k,181) = b(k,181) - lu(k,1811) * b(k,217)
         b(k,180) = b(k,180) - lu(k,1810) * b(k,217)
         b(k,179) = b(k,179) - lu(k,1809) * b(k,217)
         b(k,178) = b(k,178) - lu(k,1808) * b(k,217)
         b(k,176) = b(k,176) - lu(k,1807) * b(k,217)
         b(k,174) = b(k,174) - lu(k,1806) * b(k,217)
         b(k,170) = b(k,170) - lu(k,1805) * b(k,217)
         b(k,168) = b(k,168) - lu(k,1804) * b(k,217)
         b(k,167) = b(k,167) - lu(k,1803) * b(k,217)
         b(k,164) = b(k,164) - lu(k,1802) * b(k,217)
         b(k,163) = b(k,163) - lu(k,1801) * b(k,217)
         b(k,162) = b(k,162) - lu(k,1800) * b(k,217)
         b(k,161) = b(k,161) - lu(k,1799) * b(k,217)
         b(k,160) = b(k,160) - lu(k,1798) * b(k,217)
         b(k,159) = b(k,159) - lu(k,1797) * b(k,217)
         b(k,155) = b(k,155) - lu(k,1796) * b(k,217)
         b(k,154) = b(k,154) - lu(k,1795) * b(k,217)
         b(k,150) = b(k,150) - lu(k,1794) * b(k,217)
         b(k,149) = b(k,149) - lu(k,1793) * b(k,217)
         b(k,145) = b(k,145) - lu(k,1792) * b(k,217)
         b(k,144) = b(k,144) - lu(k,1791) * b(k,217)
         b(k,142) = b(k,142) - lu(k,1790) * b(k,217)
         b(k,141) = b(k,141) - lu(k,1789) * b(k,217)
         b(k,136) = b(k,136) - lu(k,1788) * b(k,217)
         b(k,134) = b(k,134) - lu(k,1787) * b(k,217)
         b(k,133) = b(k,133) - lu(k,1786) * b(k,217)
         b(k,132) = b(k,132) - lu(k,1785) * b(k,217)
         b(k,131) = b(k,131) - lu(k,1784) * b(k,217)
         b(k,128) = b(k,128) - lu(k,1783) * b(k,217)
         b(k,127) = b(k,127) - lu(k,1782) * b(k,217)
         b(k,126) = b(k,126) - lu(k,1781) * b(k,217)
         b(k,124) = b(k,124) - lu(k,1780) * b(k,217)
         b(k,123) = b(k,123) - lu(k,1779) * b(k,217)
         b(k,105) = b(k,105) - lu(k,1778) * b(k,217)
         b(k,104) = b(k,104) - lu(k,1777) * b(k,217)
         b(k,96) = b(k,96) - lu(k,1776) * b(k,217)
         b(k,94) = b(k,94) - lu(k,1775) * b(k,217)
         b(k,90) = b(k,90) - lu(k,1774) * b(k,217)
         b(k,88) = b(k,88) - lu(k,1773) * b(k,217)
         b(k,51) = b(k,51) - lu(k,1772) * b(k,217)
         b(k,50) = b(k,50) - lu(k,1771) * b(k,217)
         b(k,49) = b(k,49) - lu(k,1770) * b(k,217)
         b(k,47) = b(k,47) - lu(k,1769) * b(k,217)
         b(k,46) = b(k,46) - lu(k,1768) * b(k,217)
         b(k,45) = b(k,45) - lu(k,1767) * b(k,217)
         b(k,44) = b(k,44) - lu(k,1766) * b(k,217)
         b(k,41) = b(k,41) - lu(k,1765) * b(k,217)
         b(k,40) = b(k,40) - lu(k,1764) * b(k,217)
         b(k,39) = b(k,39) - lu(k,1763) * b(k,217)
         b(k,38) = b(k,38) - lu(k,1762) * b(k,217)
         b(k,37) = b(k,37) - lu(k,1761) * b(k,217)
                                                                        
      end do
                                                                        
      end subroutine lu_slv08
                                                                        
      subroutine lu_slv09( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
         b(k,216) = b(k,216) * lu(k,1749)
         b(k,215) = b(k,215) - lu(k,1748) * b(k,216)
         b(k,214) = b(k,214) - lu(k,1747) * b(k,216)
         b(k,213) = b(k,213) - lu(k,1746) * b(k,216)
         b(k,212) = b(k,212) - lu(k,1745) * b(k,216)
         b(k,211) = b(k,211) - lu(k,1744) * b(k,216)
         b(k,210) = b(k,210) - lu(k,1743) * b(k,216)
         b(k,209) = b(k,209) - lu(k,1742) * b(k,216)
         b(k,207) = b(k,207) - lu(k,1741) * b(k,216)
         b(k,206) = b(k,206) - lu(k,1740) * b(k,216)
         b(k,205) = b(k,205) - lu(k,1739) * b(k,216)
         b(k,204) = b(k,204) - lu(k,1738) * b(k,216)
         b(k,203) = b(k,203) - lu(k,1737) * b(k,216)
         b(k,202) = b(k,202) - lu(k,1736) * b(k,216)
         b(k,201) = b(k,201) - lu(k,1735) * b(k,216)
         b(k,200) = b(k,200) - lu(k,1734) * b(k,216)
         b(k,199) = b(k,199) - lu(k,1733) * b(k,216)
         b(k,198) = b(k,198) - lu(k,1732) * b(k,216)
         b(k,197) = b(k,197) - lu(k,1731) * b(k,216)
         b(k,196) = b(k,196) - lu(k,1730) * b(k,216)
         b(k,195) = b(k,195) - lu(k,1729) * b(k,216)
         b(k,194) = b(k,194) - lu(k,1728) * b(k,216)
         b(k,193) = b(k,193) - lu(k,1727) * b(k,216)
         b(k,192) = b(k,192) - lu(k,1726) * b(k,216)
         b(k,191) = b(k,191) - lu(k,1725) * b(k,216)
         b(k,190) = b(k,190) - lu(k,1724) * b(k,216)
         b(k,189) = b(k,189) - lu(k,1723) * b(k,216)
         b(k,188) = b(k,188) - lu(k,1722) * b(k,216)
         b(k,187) = b(k,187) - lu(k,1721) * b(k,216)
         b(k,186) = b(k,186) - lu(k,1720) * b(k,216)
         b(k,185) = b(k,185) - lu(k,1719) * b(k,216)
         b(k,184) = b(k,184) - lu(k,1718) * b(k,216)
         b(k,183) = b(k,183) - lu(k,1717) * b(k,216)
         b(k,182) = b(k,182) - lu(k,1716) * b(k,216)
         b(k,181) = b(k,181) - lu(k,1715) * b(k,216)
         b(k,180) = b(k,180) - lu(k,1714) * b(k,216)
         b(k,174) = b(k,174) - lu(k,1713) * b(k,216)
         b(k,173) = b(k,173) - lu(k,1712) * b(k,216)
         b(k,172) = b(k,172) - lu(k,1711) * b(k,216)
         b(k,142) = b(k,142) - lu(k,1710) * b(k,216)
         b(k,110) = b(k,110) - lu(k,1709) * b(k,216)
         b(k,104) = b(k,104) - lu(k,1708) * b(k,216)
         b(k,100) = b(k,100) - lu(k,1707) * b(k,216)
         b(k,95) = b(k,95) - lu(k,1706) * b(k,216)
         b(k,41) = b(k,41) - lu(k,1705) * b(k,216)
         b(k,40) = b(k,40) - lu(k,1704) * b(k,216)
                                                                        
         b(k,215) = b(k,215) * lu(k,1691)
         b(k,214) = b(k,214) - lu(k,1690) * b(k,215)
         b(k,213) = b(k,213) - lu(k,1689) * b(k,215)
         b(k,212) = b(k,212) - lu(k,1688) * b(k,215)
         b(k,211) = b(k,211) - lu(k,1687) * b(k,215)
         b(k,210) = b(k,210) - lu(k,1686) * b(k,215)
         b(k,209) = b(k,209) - lu(k,1685) * b(k,215)
         b(k,208) = b(k,208) - lu(k,1684) * b(k,215)
         b(k,207) = b(k,207) - lu(k,1683) * b(k,215)
         b(k,206) = b(k,206) - lu(k,1682) * b(k,215)
         b(k,205) = b(k,205) - lu(k,1681) * b(k,215)
         b(k,204) = b(k,204) - lu(k,1680) * b(k,215)
         b(k,203) = b(k,203) - lu(k,1679) * b(k,215)
         b(k,202) = b(k,202) - lu(k,1678) * b(k,215)
         b(k,201) = b(k,201) - lu(k,1677) * b(k,215)
         b(k,200) = b(k,200) - lu(k,1676) * b(k,215)
         b(k,199) = b(k,199) - lu(k,1675) * b(k,215)
         b(k,198) = b(k,198) - lu(k,1674) * b(k,215)
         b(k,197) = b(k,197) - lu(k,1673) * b(k,215)
         b(k,196) = b(k,196) - lu(k,1672) * b(k,215)
         b(k,195) = b(k,195) - lu(k,1671) * b(k,215)
         b(k,194) = b(k,194) - lu(k,1670) * b(k,215)
         b(k,193) = b(k,193) - lu(k,1669) * b(k,215)
         b(k,192) = b(k,192) - lu(k,1668) * b(k,215)
         b(k,191) = b(k,191) - lu(k,1667) * b(k,215)
         b(k,190) = b(k,190) - lu(k,1666) * b(k,215)
         b(k,189) = b(k,189) - lu(k,1665) * b(k,215)
         b(k,188) = b(k,188) - lu(k,1664) * b(k,215)
         b(k,187) = b(k,187) - lu(k,1663) * b(k,215)
         b(k,186) = b(k,186) - lu(k,1662) * b(k,215)
         b(k,185) = b(k,185) - lu(k,1661) * b(k,215)
         b(k,184) = b(k,184) - lu(k,1660) * b(k,215)
         b(k,183) = b(k,183) - lu(k,1659) * b(k,215)
         b(k,182) = b(k,182) - lu(k,1658) * b(k,215)
         b(k,181) = b(k,181) - lu(k,1657) * b(k,215)
         b(k,180) = b(k,180) - lu(k,1656) * b(k,215)
         b(k,179) = b(k,179) - lu(k,1655) * b(k,215)
         b(k,178) = b(k,178) - lu(k,1654) * b(k,215)
         b(k,177) = b(k,177) - lu(k,1653) * b(k,215)
         b(k,176) = b(k,176) - lu(k,1652) * b(k,215)
         b(k,175) = b(k,175) - lu(k,1651) * b(k,215)
         b(k,174) = b(k,174) - lu(k,1650) * b(k,215)
         b(k,173) = b(k,173) - lu(k,1649) * b(k,215)
         b(k,172) = b(k,172) - lu(k,1648) * b(k,215)
         b(k,171) = b(k,171) - lu(k,1647) * b(k,215)
         b(k,170) = b(k,170) - lu(k,1646) * b(k,215)
         b(k,169) = b(k,169) - lu(k,1645) * b(k,215)
         b(k,168) = b(k,168) - lu(k,1644) * b(k,215)
         b(k,167) = b(k,167) - lu(k,1643) * b(k,215)
         b(k,166) = b(k,166) - lu(k,1642) * b(k,215)
         b(k,164) = b(k,164) - lu(k,1641) * b(k,215)
         b(k,163) = b(k,163) - lu(k,1640) * b(k,215)
         b(k,162) = b(k,162) - lu(k,1639) * b(k,215)
         b(k,161) = b(k,161) - lu(k,1638) * b(k,215)
         b(k,160) = b(k,160) - lu(k,1637) * b(k,215)
         b(k,159) = b(k,159) - lu(k,1636) * b(k,215)
         b(k,158) = b(k,158) - lu(k,1635) * b(k,215)
         b(k,157) = b(k,157) - lu(k,1634) * b(k,215)
         b(k,156) = b(k,156) - lu(k,1633) * b(k,215)
         b(k,155) = b(k,155) - lu(k,1632) * b(k,215)
         b(k,154) = b(k,154) - lu(k,1631) * b(k,215)
         b(k,153) = b(k,153) - lu(k,1630) * b(k,215)
         b(k,152) = b(k,152) - lu(k,1629) * b(k,215)
         b(k,151) = b(k,151) - lu(k,1628) * b(k,215)
         b(k,150) = b(k,150) - lu(k,1627) * b(k,215)
         b(k,149) = b(k,149) - lu(k,1626) * b(k,215)
         b(k,148) = b(k,148) - lu(k,1625) * b(k,215)
         b(k,147) = b(k,147) - lu(k,1624) * b(k,215)
         b(k,146) = b(k,146) - lu(k,1623) * b(k,215)
         b(k,145) = b(k,145) - lu(k,1622) * b(k,215)
         b(k,144) = b(k,144) - lu(k,1621) * b(k,215)
         b(k,143) = b(k,143) - lu(k,1620) * b(k,215)
         b(k,142) = b(k,142) - lu(k,1619) * b(k,215)
         b(k,141) = b(k,141) - lu(k,1618) * b(k,215)
         b(k,140) = b(k,140) - lu(k,1617) * b(k,215)
         b(k,138) = b(k,138) - lu(k,1616) * b(k,215)
         b(k,137) = b(k,137) - lu(k,1615) * b(k,215)
         b(k,136) = b(k,136) - lu(k,1614) * b(k,215)
         b(k,135) = b(k,135) - lu(k,1613) * b(k,215)
         b(k,134) = b(k,134) - lu(k,1612) * b(k,215)
         b(k,133) = b(k,133) - lu(k,1611) * b(k,215)
         b(k,132) = b(k,132) - lu(k,1610) * b(k,215)
         b(k,131) = b(k,131) - lu(k,1609) * b(k,215)
         b(k,130) = b(k,130) - lu(k,1608) * b(k,215)
         b(k,129) = b(k,129) - lu(k,1607) * b(k,215)
         b(k,128) = b(k,128) - lu(k,1606) * b(k,215)
         b(k,127) = b(k,127) - lu(k,1605) * b(k,215)
         b(k,126) = b(k,126) - lu(k,1604) * b(k,215)
         b(k,125) = b(k,125) - lu(k,1603) * b(k,215)
         b(k,123) = b(k,123) - lu(k,1602) * b(k,215)
         b(k,122) = b(k,122) - lu(k,1601) * b(k,215)
         b(k,121) = b(k,121) - lu(k,1600) * b(k,215)
         b(k,120) = b(k,120) - lu(k,1599) * b(k,215)
         b(k,119) = b(k,119) - lu(k,1598) * b(k,215)
         b(k,118) = b(k,118) - lu(k,1597) * b(k,215)
         b(k,117) = b(k,117) - lu(k,1596) * b(k,215)
         b(k,116) = b(k,116) - lu(k,1595) * b(k,215)
         b(k,115) = b(k,115) - lu(k,1594) * b(k,215)
         b(k,113) = b(k,113) - lu(k,1593) * b(k,215)
         b(k,112) = b(k,112) - lu(k,1592) * b(k,215)
         b(k,111) = b(k,111) - lu(k,1591) * b(k,215)
         b(k,110) = b(k,110) - lu(k,1590) * b(k,215)
         b(k,109) = b(k,109) - lu(k,1589) * b(k,215)
         b(k,108) = b(k,108) - lu(k,1588) * b(k,215)
         b(k,107) = b(k,107) - lu(k,1587) * b(k,215)
         b(k,104) = b(k,104) - lu(k,1586) * b(k,215)
         b(k,103) = b(k,103) - lu(k,1585) * b(k,215)
         b(k,102) = b(k,102) - lu(k,1584) * b(k,215)
         b(k,101) = b(k,101) - lu(k,1583) * b(k,215)
         b(k,100) = b(k,100) - lu(k,1582) * b(k,215)
         b(k,99) = b(k,99) - lu(k,1581) * b(k,215)
         b(k,98) = b(k,98) - lu(k,1580) * b(k,215)
         b(k,93) = b(k,93) - lu(k,1579) * b(k,215)
         b(k,92) = b(k,92) - lu(k,1578) * b(k,215)
         b(k,91) = b(k,91) - lu(k,1577) * b(k,215)
         b(k,90) = b(k,90) - lu(k,1576) * b(k,215)
         b(k,89) = b(k,89) - lu(k,1575) * b(k,215)
         b(k,88) = b(k,88) - lu(k,1574) * b(k,215)
         b(k,86) = b(k,86) - lu(k,1573) * b(k,215)
         b(k,85) = b(k,85) - lu(k,1572) * b(k,215)
         b(k,84) = b(k,84) - lu(k,1571) * b(k,215)
         b(k,83) = b(k,83) - lu(k,1570) * b(k,215)
         b(k,82) = b(k,82) - lu(k,1569) * b(k,215)
         b(k,81) = b(k,81) - lu(k,1568) * b(k,215)
         b(k,80) = b(k,80) - lu(k,1567) * b(k,215)
         b(k,79) = b(k,79) - lu(k,1566) * b(k,215)
         b(k,78) = b(k,78) - lu(k,1565) * b(k,215)
         b(k,77) = b(k,77) - lu(k,1564) * b(k,215)
         b(k,76) = b(k,76) - lu(k,1563) * b(k,215)
         b(k,74) = b(k,74) - lu(k,1562) * b(k,215)
         b(k,73) = b(k,73) - lu(k,1561) * b(k,215)
         b(k,72) = b(k,72) - lu(k,1560) * b(k,215)
         b(k,71) = b(k,71) - lu(k,1559) * b(k,215)
         b(k,64) = b(k,64) - lu(k,1558) * b(k,215)
         b(k,61) = b(k,61) - lu(k,1557) * b(k,215)
         b(k,57) = b(k,57) - lu(k,1556) * b(k,215)
         b(k,55) = b(k,55) - lu(k,1555) * b(k,215)
         b(k,53) = b(k,53) - lu(k,1554) * b(k,215)
         b(k,51) = b(k,51) - lu(k,1553) * b(k,215)
         b(k,50) = b(k,50) - lu(k,1552) * b(k,215)
         b(k,49) = b(k,49) - lu(k,1551) * b(k,215)
         b(k,48) = b(k,48) - lu(k,1550) * b(k,215)
         b(k,47) = b(k,47) - lu(k,1549) * b(k,215)
         b(k,46) = b(k,46) - lu(k,1548) * b(k,215)
         b(k,45) = b(k,45) - lu(k,1547) * b(k,215)
         b(k,44) = b(k,44) - lu(k,1546) * b(k,215)
         b(k,43) = b(k,43) - lu(k,1545) * b(k,215)
         b(k,41) = b(k,41) - lu(k,1544) * b(k,215)
         b(k,40) = b(k,40) - lu(k,1543) * b(k,215)
         b(k,39) = b(k,39) - lu(k,1542) * b(k,215)
         b(k,38) = b(k,38) - lu(k,1541) * b(k,215)
         b(k,37) = b(k,37) - lu(k,1540) * b(k,215)
                                                                        
         b(k,214) = b(k,214) * lu(k,1526)
         b(k,213) = b(k,213) - lu(k,1525) * b(k,214)
         b(k,212) = b(k,212) - lu(k,1524) * b(k,214)
         b(k,211) = b(k,211) - lu(k,1523) * b(k,214)
         b(k,210) = b(k,210) - lu(k,1522) * b(k,214)
         b(k,209) = b(k,209) - lu(k,1521) * b(k,214)
         b(k,208) = b(k,208) - lu(k,1520) * b(k,214)
         b(k,175) = b(k,175) - lu(k,1519) * b(k,214)
         b(k,169) = b(k,169) - lu(k,1518) * b(k,214)
         b(k,146) = b(k,146) - lu(k,1517) * b(k,214)
         b(k,129) = b(k,129) - lu(k,1516) * b(k,214)
         b(k,125) = b(k,125) - lu(k,1515) * b(k,214)
         b(k,101) = b(k,101) - lu(k,1514) * b(k,214)
         b(k,89) = b(k,89) - lu(k,1513) * b(k,214)
         b(k,85) = b(k,85) - lu(k,1512) * b(k,214)
         b(k,83) = b(k,83) - lu(k,1511) * b(k,214)
         b(k,82) = b(k,82) - lu(k,1510) * b(k,214)
         b(k,75) = b(k,75) - lu(k,1509) * b(k,214)
         b(k,74) = b(k,74) - lu(k,1508) * b(k,214)
         b(k,69) = b(k,69) - lu(k,1507) * b(k,214)
         b(k,68) = b(k,68) - lu(k,1506) * b(k,214)
         b(k,67) = b(k,67) - lu(k,1505) * b(k,214)
         b(k,66) = b(k,66) - lu(k,1504) * b(k,214)
         b(k,60) = b(k,60) - lu(k,1503) * b(k,214)
         b(k,59) = b(k,59) - lu(k,1502) * b(k,214)
         b(k,58) = b(k,58) - lu(k,1501) * b(k,214)
         b(k,56) = b(k,56) - lu(k,1500) * b(k,214)
         b(k,54) = b(k,54) - lu(k,1499) * b(k,214)
                                                                        
      end do
                                                                        
      end subroutine lu_slv09
                                                                        
      subroutine lu_slv10( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
         b(k,213) = b(k,213) * lu(k,1485)
         b(k,212) = b(k,212) - lu(k,1484) * b(k,213)
         b(k,211) = b(k,211) - lu(k,1483) * b(k,213)
         b(k,210) = b(k,210) - lu(k,1482) * b(k,213)
         b(k,209) = b(k,209) - lu(k,1481) * b(k,213)
         b(k,208) = b(k,208) - lu(k,1480) * b(k,213)
         b(k,192) = b(k,192) - lu(k,1479) * b(k,213)
         b(k,181) = b(k,181) - lu(k,1478) * b(k,213)
         b(k,169) = b(k,169) - lu(k,1477) * b(k,213)
         b(k,124) = b(k,124) - lu(k,1476) * b(k,213)
                                                                        
         b(k,212) = b(k,212) * lu(k,1463)
         b(k,211) = b(k,211) - lu(k,1462) * b(k,212)
         b(k,210) = b(k,210) - lu(k,1461) * b(k,212)
         b(k,208) = b(k,208) - lu(k,1460) * b(k,212)
                                                                        
         b(k,211) = b(k,211) * lu(k,1447)
         b(k,210) = b(k,210) - lu(k,1446) * b(k,211)
         b(k,208) = b(k,208) - lu(k,1445) * b(k,211)
         b(k,175) = b(k,175) - lu(k,1444) * b(k,211)
         b(k,97) = b(k,97) - lu(k,1443) * b(k,211)
                                                                        
         b(k,210) = b(k,210) * lu(k,1430)
         b(k,208) = b(k,208) - lu(k,1429) * b(k,210)
         b(k,192) = b(k,192) - lu(k,1428) * b(k,210)
         b(k,175) = b(k,175) - lu(k,1427) * b(k,210)
         b(k,168) = b(k,168) - lu(k,1426) * b(k,210)
         b(k,97) = b(k,97) - lu(k,1425) * b(k,210)
                                                                        
         b(k,209) = b(k,209) * lu(k,1415)
         b(k,192) = b(k,192) - lu(k,1414) * b(k,209)
         b(k,169) = b(k,169) - lu(k,1413) * b(k,209)
                                                                        
         b(k,208) = b(k,208) * lu(k,1401)
         b(k,175) = b(k,175) - lu(k,1400) * b(k,208)
         b(k,97) = b(k,97) - lu(k,1399) * b(k,208)
                                                                        
         b(k,207) = b(k,207) * lu(k,1386)
         b(k,206) = b(k,206) - lu(k,1385) * b(k,207)
         b(k,205) = b(k,205) - lu(k,1384) * b(k,207)
         b(k,204) = b(k,204) - lu(k,1383) * b(k,207)
         b(k,203) = b(k,203) - lu(k,1382) * b(k,207)
         b(k,202) = b(k,202) - lu(k,1381) * b(k,207)
         b(k,200) = b(k,200) - lu(k,1380) * b(k,207)
         b(k,199) = b(k,199) - lu(k,1379) * b(k,207)
         b(k,198) = b(k,198) - lu(k,1378) * b(k,207)
         b(k,197) = b(k,197) - lu(k,1377) * b(k,207)
         b(k,195) = b(k,195) - lu(k,1376) * b(k,207)
         b(k,192) = b(k,192) - lu(k,1375) * b(k,207)
         b(k,191) = b(k,191) - lu(k,1374) * b(k,207)
         b(k,189) = b(k,189) - lu(k,1373) * b(k,207)
         b(k,184) = b(k,184) - lu(k,1372) * b(k,207)
         b(k,168) = b(k,168) - lu(k,1371) * b(k,207)
         b(k,156) = b(k,156) - lu(k,1370) * b(k,207)
         b(k,148) = b(k,148) - lu(k,1369) * b(k,207)
         b(k,137) = b(k,137) - lu(k,1368) * b(k,207)
         b(k,104) = b(k,104) - lu(k,1367) * b(k,207)
                                                                        
         b(k,206) = b(k,206) * lu(k,1354)
         b(k,200) = b(k,200) - lu(k,1353) * b(k,206)
         b(k,195) = b(k,195) - lu(k,1352) * b(k,206)
         b(k,192) = b(k,192) - lu(k,1351) * b(k,206)
         b(k,168) = b(k,168) - lu(k,1350) * b(k,206)
         b(k,156) = b(k,156) - lu(k,1349) * b(k,206)
         b(k,148) = b(k,148) - lu(k,1348) * b(k,206)
         b(k,143) = b(k,143) - lu(k,1347) * b(k,206)
                                                                        
         b(k,205) = b(k,205) * lu(k,1332)
         b(k,204) = b(k,204) - lu(k,1331) * b(k,205)
         b(k,200) = b(k,200) - lu(k,1330) * b(k,205)
         b(k,195) = b(k,195) - lu(k,1329) * b(k,205)
         b(k,192) = b(k,192) - lu(k,1328) * b(k,205)
         b(k,190) = b(k,190) - lu(k,1327) * b(k,205)
         b(k,186) = b(k,186) - lu(k,1326) * b(k,205)
         b(k,181) = b(k,181) - lu(k,1325) * b(k,205)
         b(k,168) = b(k,168) - lu(k,1324) * b(k,205)
                                                                        
         b(k,204) = b(k,204) * lu(k,1311)
         b(k,200) = b(k,200) - lu(k,1310) * b(k,204)
         b(k,196) = b(k,196) - lu(k,1309) * b(k,204)
         b(k,195) = b(k,195) - lu(k,1308) * b(k,204)
         b(k,192) = b(k,192) - lu(k,1307) * b(k,204)
         b(k,191) = b(k,191) - lu(k,1306) * b(k,204)
         b(k,166) = b(k,166) - lu(k,1305) * b(k,204)
         b(k,102) = b(k,102) - lu(k,1304) * b(k,204)
                                                                        
         b(k,203) = b(k,203) * lu(k,1288)
         b(k,200) = b(k,200) - lu(k,1287) * b(k,203)
         b(k,199) = b(k,199) - lu(k,1286) * b(k,203)
         b(k,197) = b(k,197) - lu(k,1285) * b(k,203)
         b(k,196) = b(k,196) - lu(k,1284) * b(k,203)
         b(k,195) = b(k,195) - lu(k,1283) * b(k,203)
         b(k,192) = b(k,192) - lu(k,1282) * b(k,203)
         b(k,191) = b(k,191) - lu(k,1281) * b(k,203)
         b(k,184) = b(k,184) - lu(k,1280) * b(k,203)
         b(k,176) = b(k,176) - lu(k,1279) * b(k,203)
         b(k,174) = b(k,174) - lu(k,1278) * b(k,203)
         b(k,166) = b(k,166) - lu(k,1277) * b(k,203)
         b(k,157) = b(k,157) - lu(k,1276) * b(k,203)
         b(k,145) = b(k,145) - lu(k,1275) * b(k,203)
         b(k,141) = b(k,141) - lu(k,1274) * b(k,203)
         b(k,104) = b(k,104) - lu(k,1273) * b(k,203)
         b(k,84) = b(k,84) - lu(k,1272) * b(k,203)
                                                                        
         b(k,202) = b(k,202) * lu(k,1256)
         b(k,200) = b(k,200) - lu(k,1255) * b(k,202)
         b(k,199) = b(k,199) - lu(k,1254) * b(k,202)
         b(k,197) = b(k,197) - lu(k,1253) * b(k,202)
         b(k,196) = b(k,196) - lu(k,1252) * b(k,202)
         b(k,195) = b(k,195) - lu(k,1251) * b(k,202)
         b(k,192) = b(k,192) - lu(k,1250) * b(k,202)
         b(k,191) = b(k,191) - lu(k,1249) * b(k,202)
         b(k,168) = b(k,168) - lu(k,1248) * b(k,202)
         b(k,166) = b(k,166) - lu(k,1247) * b(k,202)
         b(k,157) = b(k,157) - lu(k,1246) * b(k,202)
         b(k,144) = b(k,144) - lu(k,1245) * b(k,202)
                                                                        
         b(k,201) = b(k,201) * lu(k,1232)
         b(k,172) = b(k,172) - lu(k,1231) * b(k,201)
         b(k,135) = b(k,135) - lu(k,1230) * b(k,201)
         b(k,106) = b(k,106) - lu(k,1229) * b(k,201)
                                                                        
         b(k,200) = b(k,200) * lu(k,1221)
         b(k,192) = b(k,192) - lu(k,1220) * b(k,200)
                                                                        
         b(k,199) = b(k,199) * lu(k,1209)
         b(k,192) = b(k,192) - lu(k,1208) * b(k,199)
         b(k,181) = b(k,181) - lu(k,1207) * b(k,199)
                                                                        
         b(k,198) = b(k,198) * lu(k,1193)
         b(k,197) = b(k,197) - lu(k,1192) * b(k,198)
         b(k,192) = b(k,192) - lu(k,1191) * b(k,198)
         b(k,191) = b(k,191) - lu(k,1190) * b(k,198)
         b(k,189) = b(k,189) - lu(k,1189) * b(k,198)
         b(k,174) = b(k,174) - lu(k,1188) * b(k,198)
         b(k,166) = b(k,166) - lu(k,1187) * b(k,198)
         b(k,157) = b(k,157) - lu(k,1186) * b(k,198)
         b(k,117) = b(k,117) - lu(k,1185) * b(k,198)
         b(k,113) = b(k,113) - lu(k,1184) * b(k,198)
                                                                        
         b(k,197) = b(k,197) * lu(k,1173)
         b(k,195) = b(k,195) - lu(k,1172) * b(k,197)
         b(k,192) = b(k,192) - lu(k,1171) * b(k,197)
         b(k,191) = b(k,191) - lu(k,1170) * b(k,197)
         b(k,184) = b(k,184) - lu(k,1169) * b(k,197)
         b(k,168) = b(k,168) - lu(k,1168) * b(k,197)
         b(k,166) = b(k,166) - lu(k,1167) * b(k,197)
         b(k,79) = b(k,79) - lu(k,1166) * b(k,197)
                                                                        
         b(k,196) = b(k,196) * lu(k,1153)
         b(k,195) = b(k,195) - lu(k,1152) * b(k,196)
         b(k,194) = b(k,194) - lu(k,1151) * b(k,196)
         b(k,192) = b(k,192) - lu(k,1150) * b(k,196)
         b(k,191) = b(k,191) - lu(k,1149) * b(k,196)
         b(k,190) = b(k,190) - lu(k,1148) * b(k,196)
         b(k,180) = b(k,180) - lu(k,1147) * b(k,196)
         b(k,88) = b(k,88) - lu(k,1146) * b(k,196)
                                                                        
         b(k,195) = b(k,195) * lu(k,1140)
                                                                        
         b(k,194) = b(k,194) * lu(k,1129)
         b(k,166) = b(k,166) - lu(k,1128) * b(k,194)
         b(k,119) = b(k,119) - lu(k,1127) * b(k,194)
                                                                        
         b(k,193) = b(k,193) * lu(k,1115)
         b(k,192) = b(k,192) - lu(k,1114) * b(k,193)
         b(k,191) = b(k,191) - lu(k,1113) * b(k,193)
         b(k,188) = b(k,188) - lu(k,1112) * b(k,193)
         b(k,180) = b(k,180) - lu(k,1111) * b(k,193)
         b(k,168) = b(k,168) - lu(k,1110) * b(k,193)
         b(k,166) = b(k,166) - lu(k,1109) * b(k,193)
         b(k,152) = b(k,152) - lu(k,1108) * b(k,193)
         b(k,88) = b(k,88) - lu(k,1107) * b(k,193)
                                                                        
         b(k,192) = b(k,192) * lu(k,1103)
         b(k,168) = b(k,168) - lu(k,1102) * b(k,192)
                                                                        
         b(k,191) = b(k,191) * lu(k,1096)
         b(k,184) = b(k,184) - lu(k,1095) * b(k,191)
         b(k,168) = b(k,168) - lu(k,1094) * b(k,191)
                                                                        
         b(k,190) = b(k,190) * lu(k,1085)
                                                                        
         b(k,189) = b(k,189) * lu(k,1073)
         b(k,184) = b(k,184) - lu(k,1072) * b(k,189)
         b(k,176) = b(k,176) - lu(k,1071) * b(k,189)
         b(k,174) = b(k,174) - lu(k,1070) * b(k,189)
         b(k,145) = b(k,145) - lu(k,1069) * b(k,189)
                                                                        
         b(k,188) = b(k,188) * lu(k,1059)
         b(k,180) = b(k,180) - lu(k,1058) * b(k,188)
         b(k,168) = b(k,168) - lu(k,1057) * b(k,188)
                                                                        
         b(k,187) = b(k,187) * lu(k,1045)
         b(k,183) = b(k,183) - lu(k,1044) * b(k,187)
         b(k,166) = b(k,166) - lu(k,1043) * b(k,187)
         b(k,149) = b(k,149) - lu(k,1042) * b(k,187)
         b(k,116) = b(k,116) - lu(k,1041) * b(k,187)
                                                                        
         b(k,186) = b(k,186) * lu(k,1024)
         b(k,181) = b(k,181) - lu(k,1023) * b(k,186)
         b(k,174) = b(k,174) - lu(k,1022) * b(k,186)
         b(k,168) = b(k,168) - lu(k,1021) * b(k,186)
         b(k,164) = b(k,164) - lu(k,1020) * b(k,186)
         b(k,156) = b(k,156) - lu(k,1019) * b(k,186)
                                                                        
         b(k,185) = b(k,185) * lu(k,999)
         b(k,184) = b(k,184) - lu(k,998) * b(k,185)
         b(k,183) = b(k,183) - lu(k,997) * b(k,185)
         b(k,181) = b(k,181) - lu(k,996) * b(k,185)
         b(k,180) = b(k,180) - lu(k,995) * b(k,185)
         b(k,179) = b(k,179) - lu(k,994) * b(k,185)
         b(k,178) = b(k,178) - lu(k,993) * b(k,185)
         b(k,168) = b(k,168) - lu(k,992) * b(k,185)
         b(k,114) = b(k,114) - lu(k,991) * b(k,185)
         b(k,86) = b(k,86) - lu(k,990) * b(k,185)
         b(k,44) = b(k,44) - lu(k,989) * b(k,185)
         b(k,41) = b(k,41) - lu(k,988) * b(k,185)
         b(k,40) = b(k,40) - lu(k,987) * b(k,185)
         b(k,39) = b(k,39) - lu(k,986) * b(k,185)
         b(k,38) = b(k,38) - lu(k,985) * b(k,185)
         b(k,37) = b(k,37) - lu(k,984) * b(k,185)
                                                                        
         b(k,184) = b(k,184) * lu(k,979)
         b(k,168) = b(k,168) - lu(k,978) * b(k,184)
         b(k,37) = b(k,37) - lu(k,977) * b(k,184)
                                                                        
         b(k,183) = b(k,183) * lu(k,969)
                                                                        
         b(k,182) = b(k,182) * lu(k,948)
         b(k,181) = b(k,181) - lu(k,947) * b(k,182)
         b(k,180) = b(k,180) - lu(k,946) * b(k,182)
         b(k,179) = b(k,179) - lu(k,945) * b(k,182)
         b(k,178) = b(k,178) - lu(k,944) * b(k,182)
         b(k,168) = b(k,168) - lu(k,943) * b(k,182)
         b(k,114) = b(k,114) - lu(k,942) * b(k,182)
         b(k,86) = b(k,86) - lu(k,941) * b(k,182)
         b(k,49) = b(k,49) - lu(k,940) * b(k,182)
         b(k,41) = b(k,41) - lu(k,939) * b(k,182)
         b(k,40) = b(k,40) - lu(k,938) * b(k,182)
         b(k,39) = b(k,39) - lu(k,937) * b(k,182)
         b(k,38) = b(k,38) - lu(k,936) * b(k,182)
         b(k,37) = b(k,37) - lu(k,935) * b(k,182)
                                                                        
      end do
                                                                        
      end subroutine lu_slv10
                                                                        
      subroutine lu_slv11( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
         b(k,181) = b(k,181) * lu(k,929)
         b(k,168) = b(k,168) - lu(k,928) * b(k,181)
                                                                        
         b(k,180) = b(k,180) * lu(k,922)
                                                                        
         b(k,179) = b(k,179) * lu(k,912)
         b(k,166) = b(k,166) - lu(k,911) * b(k,179)
         b(k,149) = b(k,149) - lu(k,910) * b(k,179)
         b(k,130) = b(k,130) - lu(k,909) * b(k,179)
                                                                        
         b(k,178) = b(k,178) * lu(k,899)
         b(k,170) = b(k,170) - lu(k,898) * b(k,178)
         b(k,155) = b(k,155) - lu(k,897) * b(k,178)
         b(k,154) = b(k,154) - lu(k,896) * b(k,178)
         b(k,151) = b(k,151) - lu(k,895) * b(k,178)
         b(k,134) = b(k,134) - lu(k,894) * b(k,178)
                                                                        
         b(k,177) = b(k,177) * lu(k,884)
         b(k,171) = b(k,171) - lu(k,883) * b(k,177)
         b(k,70) = b(k,70) - lu(k,882) * b(k,177)
                                                                        
         b(k,176) = b(k,176) * lu(k,873)
                                                                        
         b(k,175) = b(k,175) * lu(k,864)
         b(k,97) = b(k,97) - lu(k,863) * b(k,175)
                                                                        
         b(k,174) = b(k,174) * lu(k,857)
                                                                        
         b(k,173) = b(k,173) * lu(k,841)
         b(k,46) = b(k,46) - lu(k,840) * b(k,173)
         b(k,41) = b(k,41) - lu(k,839) * b(k,173)
         b(k,40) = b(k,40) - lu(k,838) * b(k,173)
                                                                        
         b(k,172) = b(k,172) * lu(k,832)
         b(k,91) = b(k,91) - lu(k,831) * b(k,172)
                                                                        
         b(k,171) = b(k,171) * lu(k,824)
         b(k,70) = b(k,70) - lu(k,823) * b(k,171)
                                                                        
         b(k,170) = b(k,170) * lu(k,815)
         b(k,166) = b(k,166) - lu(k,814) * b(k,170)
         b(k,109) = b(k,109) - lu(k,813) * b(k,170)
         b(k,92) = b(k,92) - lu(k,812) * b(k,170)
                                                                        
         b(k,169) = b(k,169) * lu(k,804)
                                                                        
         b(k,168) = b(k,168) * lu(k,801)
                                                                        
         b(k,167) = b(k,167) * lu(k,791)
         b(k,120) = b(k,120) - lu(k,790) * b(k,167)
                                                                        
         b(k,166) = b(k,166) * lu(k,786)
                                                                        
         b(k,165) = b(k,165) * lu(k,778)
         b(k,87) = b(k,87) - lu(k,777) * b(k,165)
                                                                        
         b(k,164) = b(k,164) * lu(k,768)
         b(k,138) = b(k,138) - lu(k,767) * b(k,164)
                                                                        
         b(k,163) = b(k,163) * lu(k,759)
                                                                        
         b(k,162) = b(k,162) * lu(k,748)
         b(k,160) = b(k,160) - lu(k,747) * b(k,162)
         b(k,158) = b(k,158) - lu(k,746) * b(k,162)
         b(k,145) = b(k,145) - lu(k,745) * b(k,162)
         b(k,127) = b(k,127) - lu(k,744) * b(k,162)
         b(k,105) = b(k,105) - lu(k,743) * b(k,162)
         b(k,96) = b(k,96) - lu(k,742) * b(k,162)
                                                                        
         b(k,161) = b(k,161) * lu(k,732)
         b(k,160) = b(k,160) - lu(k,731) * b(k,161)
         b(k,153) = b(k,153) - lu(k,730) * b(k,161)
         b(k,145) = b(k,145) - lu(k,729) * b(k,161)
         b(k,127) = b(k,127) - lu(k,728) * b(k,161)
         b(k,96) = b(k,96) - lu(k,727) * b(k,161)
                                                                        
         b(k,160) = b(k,160) * lu(k,721)
                                                                        
         b(k,159) = b(k,159) * lu(k,714)
         b(k,94) = b(k,94) - lu(k,713) * b(k,159)
         b(k,62) = b(k,62) - lu(k,712) * b(k,159)
                                                                        
         b(k,158) = b(k,158) * lu(k,701)
         b(k,145) = b(k,145) - lu(k,700) * b(k,158)
         b(k,127) = b(k,127) - lu(k,699) * b(k,158)
         b(k,105) = b(k,105) - lu(k,698) * b(k,158)
         b(k,96) = b(k,96) - lu(k,697) * b(k,158)
                                                                        
         b(k,157) = b(k,157) * lu(k,690)
         b(k,64) = b(k,64) - lu(k,689) * b(k,157)
                                                                        
         b(k,156) = b(k,156) * lu(k,684)
                                                                        
         b(k,155) = b(k,155) * lu(k,677)
         b(k,103) = b(k,103) - lu(k,676) * b(k,155)
                                                                        
         b(k,154) = b(k,154) * lu(k,666)
         b(k,134) = b(k,134) - lu(k,665) * b(k,154)
                                                                        
         b(k,153) = b(k,153) * lu(k,655)
         b(k,145) = b(k,145) - lu(k,654) * b(k,153)
         b(k,127) = b(k,127) - lu(k,653) * b(k,153)
         b(k,96) = b(k,96) - lu(k,652) * b(k,153)
                                                                        
         b(k,152) = b(k,152) * lu(k,642)
                                                                        
         b(k,151) = b(k,151) * lu(k,632)
         b(k,134) = b(k,134) - lu(k,631) * b(k,151)
                                                                        
         b(k,150) = b(k,150) * lu(k,625)
         b(k,128) = b(k,128) - lu(k,624) * b(k,150)
         b(k,93) = b(k,93) - lu(k,623) * b(k,150)
                                                                        
         b(k,149) = b(k,149) * lu(k,617)
                                                                        
         b(k,148) = b(k,148) * lu(k,610)
                                                                        
         b(k,147) = b(k,147) * lu(k,603)
                                                                        
         b(k,146) = b(k,146) * lu(k,594)
                                                                        
         b(k,145) = b(k,145) * lu(k,590)
                                                                        
         b(k,144) = b(k,144) * lu(k,581)
                                                                        
         b(k,143) = b(k,143) * lu(k,572)
                                                                        
         b(k,142) = b(k,142) * lu(k,564)
                                                                        
         b(k,141) = b(k,141) * lu(k,556)
                                                                        
         b(k,140) = b(k,140) * lu(k,548)
                                                                        
         b(k,139) = b(k,139) * lu(k,540)
                                                                        
         b(k,138) = b(k,138) * lu(k,532)
                                                                        
         b(k,137) = b(k,137) * lu(k,524)
                                                                        
         b(k,136) = b(k,136) * lu(k,518)
         b(k,65) = b(k,65) - lu(k,517) * b(k,136)
                                                                        
         b(k,135) = b(k,135) * lu(k,511)
                                                                        
         b(k,134) = b(k,134) * lu(k,506)
                                                                        
         b(k,133) = b(k,133) * lu(k,499)
         b(k,122) = b(k,122) - lu(k,498) * b(k,133)
                                                                        
         b(k,132) = b(k,132) * lu(k,491)
         b(k,75) = b(k,75) - lu(k,490) * b(k,132)
                                                                        
         b(k,131) = b(k,131) * lu(k,483)
         b(k,127) = b(k,127) - lu(k,482) * b(k,131)
         b(k,121) = b(k,121) - lu(k,481) * b(k,131)
                                                                        
         b(k,130) = b(k,130) * lu(k,474)
                                                                        
         b(k,129) = b(k,129) * lu(k,467)
                                                                        
         b(k,128) = b(k,128) * lu(k,463)
                                                                        
         b(k,127) = b(k,127) * lu(k,460)
                                                                        
         b(k,126) = b(k,126) * lu(k,454)
         b(k,107) = b(k,107) - lu(k,453) * b(k,126)
                                                                        
         b(k,125) = b(k,125) * lu(k,447)
                                                                        
         b(k,124) = b(k,124) * lu(k,441)
                                                                        
         b(k,123) = b(k,123) * lu(k,435)
         b(k,108) = b(k,108) - lu(k,434) * b(k,123)
         b(k,90) = b(k,90) - lu(k,433) * b(k,123)
                                                                        
         b(k,122) = b(k,122) * lu(k,427)
                                                                        
         b(k,121) = b(k,121) * lu(k,421)
                                                                        
         b(k,120) = b(k,120) * lu(k,415)
                                                                        
         b(k,119) = b(k,119) * lu(k,409)
                                                                        
         b(k,118) = b(k,118) * lu(k,403)
                                                                        
         b(k,117) = b(k,117) * lu(k,397)
                                                                        
         b(k,116) = b(k,116) * lu(k,391)
                                                                        
         b(k,115) = b(k,115) * lu(k,385)
                                                                        
         b(k,114) = b(k,114) * lu(k,379)
                                                                        
         b(k,113) = b(k,113) * lu(k,373)
                                                                        
         b(k,112) = b(k,112) * lu(k,365)
                                                                        
         b(k,111) = b(k,111) * lu(k,357)
                                                                        
         b(k,110) = b(k,110) * lu(k,349)
                                                                        
         b(k,109) = b(k,109) * lu(k,344)
                                                                        
         b(k,108) = b(k,108) * lu(k,339)
         b(k,90) = b(k,90) - lu(k,338) * b(k,108)
                                                                        
         b(k,107) = b(k,107) * lu(k,333)
                                                                        
         b(k,106) = b(k,106) * lu(k,328)
                                                                        
         b(k,105) = b(k,105) * lu(k,323)
                                                                        
         b(k,104) = b(k,104) * lu(k,320)
                                                                        
         b(k,103) = b(k,103) * lu(k,315)
                                                                        
         b(k,102) = b(k,102) * lu(k,310)
                                                                        
         b(k,101) = b(k,101) * lu(k,304)
         b(k,85) = b(k,85) - lu(k,303) * b(k,101)
                                                                        
         b(k,100) = b(k,100) * lu(k,297)
                                                                        
         b(k,99) = b(k,99) * lu(k,291)
                                                                        
         b(k,98) = b(k,98) * lu(k,285)
                                                                        
         b(k,97) = b(k,97) * lu(k,282)
                                                                        
         b(k,96) = b(k,96) * lu(k,279)
                                                                        
         b(k,95) = b(k,95) * lu(k,273)
                                                                        
         b(k,94) = b(k,94) * lu(k,269)
                                                                        
         b(k,93) = b(k,93) * lu(k,265)
                                                                        
         b(k,92) = b(k,92) * lu(k,261)
                                                                        
         b(k,91) = b(k,91) * lu(k,257)
         b(k,63) = b(k,63) - lu(k,256) * b(k,91)
                                                                        
         b(k,90) = b(k,90) * lu(k,253)
                                                                        
         b(k,89) = b(k,89) * lu(k,248)
         b(k,85) = b(k,85) - lu(k,247) * b(k,89)
                                                                        
         b(k,88) = b(k,88) * lu(k,244)
                                                                        
         b(k,87) = b(k,87) * lu(k,241)
                                                                        
         b(k,86) = b(k,86) * lu(k,238)
                                                                        
         b(k,85) = b(k,85) * lu(k,235)
                                                                        
         b(k,84) = b(k,84) * lu(k,230)
                                                                        
         b(k,83) = b(k,83) * lu(k,226)
                                                                        
         b(k,82) = b(k,82) * lu(k,221)
                                                                        
         b(k,81) = b(k,81) * lu(k,216)
                                                                        
         b(k,80) = b(k,80) * lu(k,208)
         b(k,78) = b(k,78) - lu(k,207) * b(k,80)
         b(k,51) = b(k,51) - lu(k,206) * b(k,80)
                                                                        
         b(k,79) = b(k,79) * lu(k,203)
                                                                        
         b(k,78) = b(k,78) * lu(k,199)
                                                                        
         b(k,77) = b(k,77) * lu(k,194)
                                                                        
         b(k,76) = b(k,76) * lu(k,187)
         b(k,50) = b(k,50) - lu(k,186) * b(k,76)
                                                                        
         b(k,75) = b(k,75) * lu(k,183)
                                                                        
         b(k,74) = b(k,74) * lu(k,179)
                                                                        
         b(k,73) = b(k,73) * lu(k,174)
                                                                        
         b(k,72) = b(k,72) * lu(k,170)
                                                                        
         b(k,71) = b(k,71) * lu(k,164)
         b(k,45) = b(k,45) - lu(k,163) * b(k,71)
                                                                        
         b(k,70) = b(k,70) * lu(k,161)
                                                                        
         b(k,69) = b(k,69) * lu(k,156)
                                                                        
         b(k,68) = b(k,68) * lu(k,151)
                                                                        
         b(k,67) = b(k,67) * lu(k,146)
                                                                        
         b(k,66) = b(k,66) * lu(k,141)
                                                                        
         b(k,65) = b(k,65) * lu(k,138)
                                                                        
         b(k,64) = b(k,64) * lu(k,135)
                                                                        
         b(k,63) = b(k,63) * lu(k,132)
                                                                        
         b(k,62) = b(k,62) * lu(k,129)
                                                                        
         b(k,61) = b(k,61) * lu(k,125)
                                                                        
         b(k,60) = b(k,60) * lu(k,121)
                                                                        
         b(k,59) = b(k,59) * lu(k,117)
                                                                        
         b(k,58) = b(k,58) * lu(k,113)
                                                                        
         b(k,57) = b(k,57) * lu(k,109)
                                                                        
         b(k,56) = b(k,56) * lu(k,105)
                                                                        
         b(k,55) = b(k,55) * lu(k,102)
                                                                        
         b(k,54) = b(k,54) * lu(k,99)
                                                                        
         b(k,53) = b(k,53) * lu(k,96)
                                                                        
         b(k,52) = b(k,52) * lu(k,93)
                                                                        
         b(k,51) = b(k,51) * lu(k,92)
         b(k,41) = b(k,41) - lu(k,91) * b(k,51)
         b(k,40) = b(k,40) - lu(k,90) * b(k,51)
         b(k,39) = b(k,39) - lu(k,89) * b(k,51)
         b(k,38) = b(k,38) - lu(k,88) * b(k,51)
         b(k,37) = b(k,37) - lu(k,87) * b(k,51)
                                                                        
      end do
                                                                        
      end subroutine lu_slv11
                                                                        
      subroutine lu_slv12( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,avec_len
         b(k,50) = b(k,50) * lu(k,86)
         b(k,41) = b(k,41) - lu(k,85) * b(k,50)
         b(k,40) = b(k,40) - lu(k,84) * b(k,50)
         b(k,39) = b(k,39) - lu(k,83) * b(k,50)
         b(k,38) = b(k,38) - lu(k,82) * b(k,50)
         b(k,37) = b(k,37) - lu(k,81) * b(k,50)
                                                                        
         b(k,49) = b(k,49) * lu(k,80)
         b(k,41) = b(k,41) - lu(k,79) * b(k,49)
         b(k,40) = b(k,40) - lu(k,78) * b(k,49)
         b(k,39) = b(k,39) - lu(k,77) * b(k,49)
         b(k,38) = b(k,38) - lu(k,76) * b(k,49)
         b(k,37) = b(k,37) - lu(k,75) * b(k,49)
                                                                        
         b(k,48) = b(k,48) * lu(k,74)
         b(k,47) = b(k,47) - lu(k,73) * b(k,48)
                                                                        
         b(k,47) = b(k,47) * lu(k,72)
         b(k,41) = b(k,41) - lu(k,71) * b(k,47)
         b(k,40) = b(k,40) - lu(k,70) * b(k,47)
         b(k,39) = b(k,39) - lu(k,69) * b(k,47)
         b(k,38) = b(k,38) - lu(k,68) * b(k,47)
         b(k,37) = b(k,37) - lu(k,67) * b(k,47)
                                                                        
         b(k,46) = b(k,46) * lu(k,66)
         b(k,41) = b(k,41) - lu(k,65) * b(k,46)
         b(k,40) = b(k,40) - lu(k,64) * b(k,46)
         b(k,39) = b(k,39) - lu(k,63) * b(k,46)
         b(k,38) = b(k,38) - lu(k,62) * b(k,46)
         b(k,37) = b(k,37) - lu(k,61) * b(k,46)
                                                                        
         b(k,45) = b(k,45) * lu(k,60)
         b(k,41) = b(k,41) - lu(k,59) * b(k,45)
         b(k,40) = b(k,40) - lu(k,58) * b(k,45)
         b(k,39) = b(k,39) - lu(k,57) * b(k,45)
         b(k,38) = b(k,38) - lu(k,56) * b(k,45)
         b(k,37) = b(k,37) - lu(k,55) * b(k,45)
                                                                        
         b(k,44) = b(k,44) * lu(k,54)
         b(k,41) = b(k,41) - lu(k,53) * b(k,44)
         b(k,40) = b(k,40) - lu(k,52) * b(k,44)
         b(k,39) = b(k,39) - lu(k,51) * b(k,44)
         b(k,38) = b(k,38) - lu(k,50) * b(k,44)
         b(k,37) = b(k,37) - lu(k,49) * b(k,44)
                                                                        
         b(k,43) = b(k,43) * lu(k,48)
         b(k,41) = b(k,41) - lu(k,47) * b(k,43)
         b(k,40) = b(k,40) - lu(k,46) * b(k,43)
         b(k,39) = b(k,39) - lu(k,45) * b(k,43)
         b(k,38) = b(k,38) - lu(k,44) * b(k,43)
         b(k,37) = b(k,37) - lu(k,43) * b(k,43)
                                                                        
         b(k,42) = b(k,42) * lu(k,42)
                                                                        
         b(k,41) = b(k,41) * lu(k,41)
                                                                        
         b(k,40) = b(k,40) * lu(k,40)
                                                                        
         b(k,39) = b(k,39) * lu(k,39)
                                                                        
         b(k,38) = b(k,38) * lu(k,38)
                                                                        
         b(k,37) = b(k,37) * lu(k,37)
                                                                        
         b(k,36) = b(k,36) * lu(k,36)
                                                                        
         b(k,35) = b(k,35) * lu(k,35)
                                                                        
         b(k,34) = b(k,34) * lu(k,34)
                                                                        
         b(k,33) = b(k,33) * lu(k,33)
                                                                        
         b(k,32) = b(k,32) * lu(k,32)
                                                                        
         b(k,31) = b(k,31) * lu(k,31)
                                                                        
         b(k,30) = b(k,30) * lu(k,30)
                                                                        
         b(k,29) = b(k,29) * lu(k,29)
                                                                        
         b(k,28) = b(k,28) * lu(k,28)
                                                                        
         b(k,27) = b(k,27) * lu(k,27)
                                                                        
         b(k,26) = b(k,26) * lu(k,26)
                                                                        
         b(k,25) = b(k,25) * lu(k,25)
                                                                        
         b(k,24) = b(k,24) * lu(k,24)
                                                                        
         b(k,23) = b(k,23) * lu(k,23)
                                                                        
         b(k,22) = b(k,22) * lu(k,22)
                                                                        
         b(k,21) = b(k,21) * lu(k,21)
                                                                        
         b(k,20) = b(k,20) * lu(k,20)
                                                                        
         b(k,19) = b(k,19) * lu(k,19)
                                                                        
         b(k,18) = b(k,18) * lu(k,18)
                                                                        
         b(k,17) = b(k,17) * lu(k,17)
                                                                        
         b(k,16) = b(k,16) * lu(k,16)
                                                                        
         b(k,15) = b(k,15) * lu(k,15)
                                                                        
         b(k,14) = b(k,14) * lu(k,14)
                                                                        
         b(k,13) = b(k,13) * lu(k,13)
                                                                        
         b(k,12) = b(k,12) * lu(k,12)
                                                                        
         b(k,11) = b(k,11) * lu(k,11)
                                                                        
         b(k,10) = b(k,10) * lu(k,10)
                                                                        
         b(k,9) = b(k,9) * lu(k,9)
                                                                        
         b(k,8) = b(k,8) * lu(k,8)
                                                                        
         b(k,7) = b(k,7) * lu(k,7)
                                                                        
         b(k,6) = b(k,6) * lu(k,6)
                                                                        
         b(k,5) = b(k,5) * lu(k,5)
                                                                        
         b(k,4) = b(k,4) * lu(k,4)
                                                                        
         b(k,3) = b(k,3) * lu(k,3)
                                                                        
         b(k,2) = b(k,2) * lu(k,2)
                                                                        
         b(k,1) = b(k,1) * lu(k,1)
                                                                        
      end do
                                                                        
      end subroutine lu_slv12
                                                                        
      subroutine lu_slv( avec_len, lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : clscnt4, nzcnt
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   avec_len
      real(r8), intent(in)    ::   lu(veclen,max(1,nzcnt))
      real(r8), intent(inout) ::   b(veclen,clscnt4)
                                                                        
      call lu_slv01( avec_len, lu, b )
      call lu_slv02( avec_len, lu, b )
      call lu_slv03( avec_len, lu, b )
      call lu_slv04( avec_len, lu, b )
      call lu_slv05( avec_len, lu, b )
      call lu_slv06( avec_len, lu, b )
      call lu_slv07( avec_len, lu, b )
      call lu_slv08( avec_len, lu, b )
      call lu_slv09( avec_len, lu, b )
      call lu_slv10( avec_len, lu, b )
      call lu_slv11( avec_len, lu, b )
      call lu_slv12( avec_len, lu, b )
                                                                        
      end subroutine lu_slv
                                                                        
      end module mo_lu_solve


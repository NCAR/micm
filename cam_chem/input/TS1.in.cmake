BEGSIM
output_unit_number = 7
output_file        = musica-perf.out
output_path        = @CAMCHEM_OUTPUT@
procout_path       = @CAMCHEM_MECHANISM_OUTPUT@
procfiles_path     = @CAMCHEM_PROCFILES_PATH@
src_path           = @CAMCHEM_BKEND_PATH@
temp_path          = @CAMCHEM_TMP_OUTPUT@
sim_dat_path       = @CAMCHEM_SIMDAT_PATH@

* Comments
* User-given Tag Description: TS1.2
* Tag database identifier : MZ327_TS1.2_20230307
* Tag created by : lke
* Tag created from branch : TS1.2
* Tag created on : 2023-03-07 13:48:22.37464-07
* Comments for this tag follow:
*     lke : 2023-03-07 : Updated with additional O1D+O3 reaction.

      SPECIES

      Solution
 ALKNIT -> C5H11ONO2,
 ALKOOH -> C5H12O2,
 AOA_NH -> CO,
 bc_a1 -> C,
 bc_a4 -> C,
 BCARY -> C15H24,
 BENZENE -> C6H6,
 BENZOOH -> C6H8O5,
 BEPOMUC -> C6H6O3,
 BIGALD -> C5H6O2,
 BIGALD1 -> C4H4O2,
 BIGALD2 -> C5H6O2,
 BIGALD3 -> C5H6O2,
 BIGALD4 -> C6H8O2,
 BIGALK -> C5H12,
 BIGENE -> C4H8,
 BR -> Br,
 BRCL -> BrCl,
 BRO -> BrO,
 BRONO2 -> BrONO2,
 BRY,
 BZALD -> C7H6O,
 BZOOH -> C7H8O2,
 C2H2,
 C2H4,
 C2H5OH,
 C2H5OOH,
 C2H6,
 C3H6,
 C3H7OOH,
 C3H8,
 C6H5OOH -> C6H5OOH,
 CCL4 -> CCl4,
 CF2CLBR -> CF2ClBr,
 CF3BR -> CF3Br,
 CFC11 -> CFCl3,
 CFC113 -> CCl2FCClF2,
 CFC114 -> CClF2CClF2,
 CFC115 -> CClF2CF3,
 CFC12 -> CF2Cl2,
 CH2BR2 -> CH2Br2,
 CH2O,
 CH3BR -> CH3Br,
 CH3CCL3 -> CH3CCl3,
 CH3CHO,
 CH3CL -> CH3Cl,
 CH3CN,
 CH3COCH3,
 CH3COCHO,
 CH3COOH,
 CH3COOOH,
 CH3OH,
 CH3OOH,
 CH4,
 CHBR3 -> CHBr3,
 CL -> Cl,
 CL2 -> Cl2,
 CL2O2 -> Cl2O2,
 CLO -> ClO,
 CLONO2 -> ClONO2,
 CLY,
 CO,
 CO2,
 COF2,
 COFCL -> COFCl,
 CRESOL -> C7H8O,
 DMS -> CH3SCH3,
 dst_a1 -> AlSiO5,
 dst_a2 -> AlSiO5,
 dst_a3 -> AlSiO5,
 E90 -> CO,
 EOOH -> HOCH2CH2OOH,
 F,
 GLYALD -> HOCH2CHO,
 GLYOXAL -> C2H2O2,
 H,
 H2,
 H2402 -> CBrF2CBrF2,
 H2O2,
 H2SO4 -> H2SO4,
 HBR -> HBr,
 HCFC141B -> CH3CCl2F,
 HCFC142B -> CH3CClF2,
 HCFC22 -> CHF2Cl,
 HCL -> HCl,
 HCN,
 HCOOH,
 HF,
 HNO3,
 HO2NO2,
 HOBR -> HOBr,
 HOCL -> HOCl,
 HONITR -> C4H9NO4,
 HPALD -> HOOCH2CCH3CHCHO,
 HYAC -> CH3COCH2OH,
 HYDRALD -> HOCH2CCH3CHCHO,
 IEPOX -> C5H10O3,
 ISOP -> C5H8,
 ISOPNITA -> C5H9NO4,
 ISOPNITB -> C5H9NO4,
 ISOPNO3 -> CH2CHCCH3OOCH2ONO2,
 ISOPNOOH -> C5H9NO5,
 ISOPOOH -> HOCH2COOHCH3CHCH2,
 IVOC -> C13H28,
 MACR -> CH2CCH3CHO,
 MACROOH -> CH3COCHOOHCH2OH,
 MEK -> C4H8O,
 MEKOOH -> C4H8O3,
 MPAN -> CH2CCH3CO3NO2,
 MTERP -> C10H16,
 MVK -> CH2CHCOCH3,
 N,
 N2O,
 N2O5,
 NC4CH2OH -> C5H9NO4,
 NC4CHO -> C5H7NO4,
 ncl_a1 -> NaCl,
 ncl_a2 -> NaCl,
 ncl_a3 -> NaCl,
 NH3,
 NH4,
 NH_5 -> CO,
 NH_50 -> CO,
 NO,
 NO2,
 NO3,
 NOA -> CH3COCH2ONO2,
 NTERPOOH -> C10H17NO5,
 num_a1 -> H,
 num_a2 -> H,
 num_a3 -> H,
 num_a4 -> H,
 O,
 O3,
 O3S -> O3,
 OCLO -> OClO,
 OCS -> OCS,
 ONITR -> C4H7NO4,
 PAN -> CH3CO3NO2,
 PBZNIT -> C7H5O3NO2,
 PHENO -> C6H5O,
 PHENOL -> C6H5OH,
 PHENOOH -> C6H8O6,
 pom_a1 -> C,
 pom_a4 -> C,
 POOH -> C3H6OHOOH,
 ROOH -> CH3COCH2OOH,
 S -> S,
 SF6,
 SO -> SO,
 SO2,
 SO3 -> SO3,
 so4_a1 -> NH4HSO4,
 so4_a2 -> NH4HSO4,
 so4_a3 -> NH4HSO4,
 soa1_a1 -> C15H38O2,
 soa1_a2 -> C15H38O2,
 soa2_a1 -> C15H38O2,
 soa2_a2 -> C15H38O2,
 soa3_a1 -> C15H38O2,
 soa3_a2 -> C15H38O2,
 soa4_a1 -> C15H38O2,
 soa4_a2 -> C15H38O2,
 soa5_a1 -> C15H38O2,
 soa5_a2 -> C15H38O2,
 SOAG0 -> C15H38O2,
 SOAG1 -> C15H38O2,
 SOAG2 -> C15H38O2,
 SOAG3 -> C15H38O2,
 SOAG4 -> C15H38O2,
 ST80_25 -> CO,
 SVOC -> C22H46,
 TEPOMUC -> C7H8O3,
 TERP2OOH -> C10H16O4,
 TERPNIT -> C10H17NO4,
 TERPOOH -> C10H18O3,
 TERPROD1 -> C10H16O2,
 TERPROD2 -> C9H14O2,
 TOLOOH -> C7H10O5,
 TOLUENE -> C7H8,
 XOOH -> HOCH2COOHCH3CHOHCHO,
 XYLENES -> C8H10,
 XYLENOOH -> C8H12O5,
 XYLOL -> C8H10O,
 XYLOLOOH -> C8H12O6,
 NHDEP -> N,
 NDEP -> N,
 ACBZO2 -> C7H5O3,
 ALKO2 -> C5H11O2,
 BCARYO2VBS -> C15H25O3,
 BENZO2 -> C6H7O5,
 BENZO2VBS -> C6H7O5,
 BZOO -> C7H7O2,
 C2H5O2,
 C3H7O2,
 C6H5O2,
 CH3CO3,
 CH3O2,
 DICARBO2 -> C5H5O4,
 ENEO2 -> C4H9O3,
 EO -> HOCH2CH2O,
 EO2 -> HOCH2CH2O2,
 HO2,
 HOCH2OO,
 ISOPAO2 -> HOC5H8O2,
 ISOPBO2 -> HOC5H8O2,
 ISOPO2VBS -> C5H9O3,
 IVOCO2VBS -> C13H29O3,
 MACRO2 -> CH3COCHO2CH2OH,
 MALO2 -> C4H3O4,
 MCO3 -> CH2CCH3CO3,
 MDIALO2 -> C4H5O4,
 MEKO2 -> C4H7O3,
 MTERPO2VBS -> C10H17O3,
 NTERPO2 -> C10H16NO5,
 O1D -> O,
 OH,
 PHENO2 -> C6H7O6,
 PO2 -> C3H6OHO2,
 RO2 -> CH3COCH2O2,
 TERP2O2 -> C10H15O4,
 TERPO2 -> C10H17O3,
 TOLO2 -> C7H9O5,
 TOLUO2VBS -> C7H9O5,
 XO2 -> HOCH2COOCH3CHOHCHO,
 XYLENO2 -> C8H11O5,
 XYLEO2VBS -> C8H11O5,
 XYLOLO2 -> C8H11O6,
 H2O

      End Solution


      Fixed
 M, O2, N2
      End Fixed

      Col-int
 O3 = 0.
 O2 = 0.
      End Col-int

      Not-Transported
 ACBZO2,
 ALKO2,
 BCARYO2VBS,
 BENZO2,
 BENZO2VBS,
 BZOO,
 C2H5O2,
 C3H7O2,
 C6H5O2,
 CH3CO3,
 CH3O2,
 DICARBO2,
 ENEO2,
 EO,
 EO2,
 HO2,
 HOCH2OO,
 ISOPAO2,
 ISOPBO2,
 ISOPO2VBS,
 IVOCO2VBS,
 MACRO2,
 MALO2,
 MCO3,
 MDIALO2,
 MEKO2,
 MTERPO2VBS,
 NTERPO2,
 O1D,
 OH,
 PHENO2,
 PO2,
 RO2,
 TERP2O2,
 TERPO2,
 TOLO2,
 TOLUO2VBS,
 XO2,
 XYLENO2,
 XYLEO2VBS,
 XYLOLO2
      End Not-Transported

   END Species


   Solution classes
      Explicit
 NHDEP
 NDEP
      End Explicit

      Implicit
 ALKNIT
 ALKOOH
 AOA_NH
 bc_a1
 bc_a4
 BCARY
 BENZENE
 BENZOOH
 BEPOMUC
 BIGALD
 BIGALD1
 BIGALD2
 BIGALD3
 BIGALD4
 BIGALK
 BIGENE
 BR
 BRCL
 BRO
 BRONO2
 BRY
 BZALD
 BZOOH
 C2H2
 C2H4
 C2H5OH
 C2H5OOH
 C2H6
 C3H6
 C3H7OOH
 C3H8
 C6H5OOH
 CCL4
 CF2CLBR
 CF3BR
 CFC11
 CFC113
 CFC114
 CFC115
 CFC12
 CH2BR2
 CH2O
 CH3BR
 CH3CCL3
 CH3CHO
 CH3CL
 CH3CN
 CH3COCH3
 CH3COCHO
 CH3COOH
 CH3COOOH
 CH3OH
 CH3OOH
 CH4
 CHBR3
 CL
 CL2
 CL2O2
 CLO
 CLONO2
 CLY
 CO
 CO2
 COF2
 COFCL
 CRESOL
 DMS
 dst_a1
 dst_a2
 dst_a3
 E90
 EOOH
 F
 GLYALD
 GLYOXAL
 H
 H2
 H2402
 H2O2
 H2SO4
 HBR
 HCFC141B
 HCFC142B
 HCFC22
 HCL
 HCN
 HCOOH
 HF
 HNO3
 HO2NO2
 HOBR
 HOCL
 HONITR
 HPALD
 HYAC
 HYDRALD
 IEPOX
 ISOP
 ISOPNITA
 ISOPNITB
 ISOPNO3
 ISOPNOOH
 ISOPOOH
 IVOC
 MACR
 MACROOH
 MEK
 MEKOOH
 MPAN
 MTERP
 MVK
 N
 N2O
 N2O5
 NC4CH2OH
 NC4CHO
 ncl_a1
 ncl_a2
 ncl_a3
 NH3
 NH4
 NH_5
 NH_50
 NO
 NO2
 NO3
 NOA
 NTERPOOH
 num_a1
 num_a2
 num_a3
 num_a4
 O
 O3
 O3S
 OCLO
 OCS
 ONITR
 PAN
 PBZNIT
 PHENO
 PHENOL
 PHENOOH
 pom_a1
 pom_a4
 POOH
 ROOH
 S
 SF6
 SO
 SO2
 SO3
 so4_a1
 so4_a2
 so4_a3
 soa1_a1
 soa1_a2
 soa2_a1
 soa2_a2
 soa3_a1
 soa3_a2
 soa4_a1
 soa4_a2
 soa5_a1
 soa5_a2
 SOAG0
 SOAG1
 SOAG2
 SOAG3
 SOAG4
 ST80_25
 SVOC
 TEPOMUC
 TERP2OOH
 TERPNIT
 TERPOOH
 TERPROD1
 TERPROD2
 TOLOOH
 TOLUENE
 XOOH
 XYLENES
 XYLENOOH
 XYLOL
 XYLOLOOH
 ACBZO2
 ALKO2
 BCARYO2VBS
 BENZO2
 BENZO2VBS
 BZOO
 C2H5O2
 C3H7O2
 C6H5O2
 CH3CO3
 CH3O2
 DICARBO2
 ENEO2
 EO
 EO2
 HO2
 HOCH2OO
 ISOPAO2
 ISOPBO2
 ISOPO2VBS
 IVOCO2VBS
 MACRO2
 MALO2
 MCO3
 MDIALO2
 MEKO2
 MTERPO2VBS
 NTERPO2
 O1D
 OH
 PHENO2
 PO2
 RO2
 TERP2O2
 TERPO2
 TOLO2
 TOLUO2VBS
 XO2
 XYLENO2
 XYLEO2VBS
 XYLOLO2
 H2O
      End Implicit

   End Solution classes


 CHEMISTRY
      Photolysis
*********************************
*** odd-oxygen
*********************************
[jh2o_b]                      H2O + hv -> H2 + O1D 
[jh2o_a]                      H2O + hv -> OH + H 
[jh2o_c]                      H2O + hv -> 2*H + O 
[jh2o2]                       H2O2 + hv -> 2*OH 
[jo2_a=userdefined,]          O2 + hv -> O + O1D 
[jo2_b=userdefined,]          O2 + hv -> 2*O 
[jo3_a]                       O3 + hv -> O1D + O2 
[jo3_b]                       O3 + hv -> O + O2 
*********************************
*** odd-nitrogen
*********************************
[jhno3]                       HNO3 + hv -> NO2 + OH 
[jho2no2_a]                   HO2NO2 + hv -> OH + NO3 
[jho2no2_b]                   HO2NO2 + hv -> NO2 + HO2 
[jn2o]                        N2O + hv -> O1D + N2 
[jn2o5_a]                     N2O5 + hv -> NO2 + NO3 
[jn2o5_b]                     N2O5 + hv -> NO + O + NO3 
[jno=userdefined,]            NO + hv -> N + O 
[jno2]                        NO2 + hv -> NO + O 
[jno3_b]                      NO3 + hv -> NO + O2 
[jno3_a]                      NO3 + hv -> NO2 + O 
*********************************
*** organics
*********************************
[jalknit->,jch3ooh]           ALKNIT + hv -> NO2 + 0.4*CH3CHO + 0.1*CH2O + 0.25*CH3COCH3 + HO2 + 0.8*MEK 
[jalkooh->,jch3ooh]           ALKOOH + hv -> 0.4*CH3CHO + 0.1*CH2O + 0.25*CH3COCH3 + 0.9*HO2 + 0.8*MEK + OH 
[jbenzooh->,jch3ooh]          BENZOOH + hv -> OH + GLYOXAL + 0.5*BIGALD1 + HO2 
[jbepomuc->,.10*jno2]         BEPOMUC + hv -> BIGALD1 + 1.5*HO2 + 1.5*CO 
[jbigald->,0.2*jno2]          BIGALD + hv -> 0.45*CO + 0.13*GLYOXAL + 0.56*HO2 + 0.13*CH3CO3 + 0.18*CH3COCHO 
[jbigald1->,.14*jno2]         BIGALD1 + hv -> 0.6*MALO2 + HO2 
[jbigald2->,.20*jno2]         BIGALD2 + hv -> 0.6*HO2 + 0.6*DICARBO2 
[jbigald3->,.20*jno2]         BIGALD3 + hv -> 0.6*HO2 + 0.6*CO + 0.6*MDIALO2 
[jbigald4->,.006*jno2]        BIGALD4 + hv -> HO2 + CO + CH3COCHO + CH3CO3 
[jbzooh->,jch3ooh]            BZOOH + hv -> BZALD + OH + HO2 
[jc2h5ooh->,jch3ooh]          C2H5OOH + hv -> CH3CHO + HO2 + OH 
[jc3h7ooh->,jch3ooh]          C3H7OOH + hv -> 0.82*CH3COCH3 + OH + HO2 
[jc6h5ooh->,jch3ooh]          C6H5OOH + hv -> PHENO + OH 
[jch2o_b]                     CH2O + hv -> CO + H2 
[jch2o_a]                     CH2O + hv -> CO + 2*H 
[jch3cho]                     CH3CHO + hv -> CH3O2 + CO + HO2 
[jacet]                       CH3COCH3 + hv -> CH3CO3 + CH3O2 
[jmgly]                       CH3COCHO + hv -> CH3CO3 + CO + HO2 
[jch3co3h->,0.28*jh2o2]       CH3COOOH + hv -> CH3O2 + OH + CO2 
[jch3ooh]                     CH3OOH + hv -> CH2O + H + OH 
[jch4_b]                      CH4 + hv -> 1.44*H2 + 0.18*CH2O + 0.18*O + 0.33*OH + 0.33*H + 0.44*CO2 + 0.38*CO + 0.05*H2O 
[jch4_a]                      CH4 + hv -> H + CH3O2 
[jco2]                        CO2 + hv -> CO + O 
[jeooh->,jch3ooh]             EOOH + hv -> EO + OH 
[jglyald]                     GLYALD + hv -> 2*HO2 + CO + CH2O 
[jglyoxal->,jmgly]            GLYOXAL + hv -> 2*CO + 2*HO2 
[jhonitr->,jch2o_a]           HONITR + hv -> NO2 + 0.67*HO2 + 0.33*CH3CHO + 0.33*CH2O + 0.33*CO + 0.33*GLYALD + 0.33*CH3CO3 + 0.17*HYAC + 0.17*CH3COCH3 
[jhpald->,.006*jno2]          HPALD + hv -> BIGALD3 + OH + HO2 
[jhyac]                       HYAC + hv -> CH3CO3 + HO2 + CH2O 
[jisopnooh->,jch3ooh]         ISOPNOOH + hv -> NO2 + HO2 + ISOPOOH 
[jisopooh->,jch3ooh]          ISOPOOH + hv -> 0.7*MVK + 0.3*MACR + OH + CH2O + HO2 
[jmacr_a]                     MACR + hv -> 1.34*HO2 + 0.66*MCO3 + 1.34*CH2O + 1.34*CH3CO3 
[jmacr_b]                     MACR + hv -> 0.66*HO2 + 1.34*CO 
[jmek->,jacet]                MEK + hv -> CH3CO3 + C2H5O2 
[jmekooh->,jch3ooh]           MEKOOH + hv -> OH + CH3CO3 + CH3CHO 
[jmpan->,jpan]                MPAN + hv -> MCO3 + NO2 
[jmvk]                        MVK + hv -> 0.7*C3H6 + 0.7*CO + 0.3*CH3O2 + 0.3*CH3CO3 
[jnc4cho->,jch2o_a]           NC4CHO + hv -> BIGALD3 + NO2 + HO2 
[jnoa->,jch2o_a]              NOA + hv -> NO2 + CH2O + CH3CO3 
[jnterpooh->,jch3ooh]         NTERPOOH + hv -> TERPROD1 + NO2 + OH 
[jonitr->,jch3cho]            ONITR + hv -> NO2 
[jpan]                        PAN + hv -> 0.6*CH3CO3 + 0.6*NO2 + 0.4*CH3O2 + 0.4*NO3 + 0.4*CO2 
[jphenooh->,jch3ooh]          PHENOOH + hv -> OH + HO2 + 0.7*GLYOXAL 
[jpooh->,jch3ooh]             POOH + hv -> CH3CHO + CH2O + HO2 + OH 
[jrooh->,jch3ooh]             ROOH + hv -> CH3CO3 + CH2O + OH 
[jtepomuc->,.10*jno2]         TEPOMUC + hv -> 0.5*CH3CO3 + HO2 + 1.5*CO 
[jterp2ooh->,jch3ooh]         TERP2OOH + hv -> OH + 0.375*CH2O + 0.3*CH3COCH3 + 0.25*CO + CO2 + TERPROD2 + HO2 + 0.25*GLYALD 
[jterpnit->,jch3ooh]          TERPNIT + hv -> TERPROD1 + NO2 + HO2 
[jterpooh->,jch3ooh]          TERPOOH + hv -> 0.4*CH2O + 0.05*CH3COCH3 + TERPROD1 + HO2 + OH 
[jterprd1->,jch3cho]          TERPROD1 + hv -> HO2 + CO + TERPROD2 
[jterprd2->,jch3cho]          TERPROD2 + hv -> 0.15*RO2 + 0.68*CH2O + 0.8*CO2 + 0.5*CH3COCH3 + 0.65*CH3CO3 + 1.2*HO2 + 1.7*CO 
[jtolooh->,jch3ooh]           TOLOOH + hv -> OH + 0.6*GLYOXAL + 0.4*CH3COCHO + HO2 + 0.2*BIGALD1 + 0.2*BIGALD2 + 0.2*BIGALD3 
[jxooh->,jch3ooh]             XOOH + hv -> OH 
[jxylenooh->,jch3ooh]         XYLENOOH + hv -> OH + HO2 + 0.34*GLYOXAL + 0.54*CH3COCHO + 0.06*BIGALD1 + 0.2*BIGALD2 + 0.15*BIGALD3 + 0.21*BIGALD4 
[jxylolooh->,jch3ooh]         XYLOLOOH + hv -> OH + 0.17*GLYOXAL + 0.51*CH3COCHO + HO2 
*********************************
*** halogens
*********************************
[jbrcl]                       BRCL + hv -> BR + CL 
[jbro]                        BRO + hv -> BR + O 
[jbrono2_b]                   BRONO2 + hv -> BRO + NO2 
[jbrono2_a]                   BRONO2 + hv -> BR + NO3 
[jccl4]                       CCL4 + hv -> 4*CL 
[jcf2clbr]                    CF2CLBR + hv -> BR + CL + COF2 
[jcf3br]                      CF3BR + hv -> BR + F + COF2 
[jcfcl3]                      CFC11 + hv -> 2*CL + COFCL 
[jcfc113]                     CFC113 + hv -> 2*CL + COFCL + COF2 
[jcfc114]                     CFC114 + hv -> 2*CL + 2*COF2 
[jcfc115]                     CFC115 + hv -> CL + F + 2*COF2 
[jcf2cl2]                     CFC12 + hv -> 2*CL + COF2 
[jch2br2]                     CH2BR2 + hv -> 2*BR 
[jch3br]                      CH3BR + hv -> BR + CH3O2 
[jch3ccl3]                    CH3CCL3 + hv -> 3*CL 
[jch3cl]                      CH3CL + hv -> CL + CH3O2 
[jchbr3]                      CHBR3 + hv -> 3*BR 
[jcl2]                        CL2 + hv -> 2*CL 
[jcl2o2]                      CL2O2 + hv -> 2*CL 
[jclo]                        CLO + hv -> CL + O 
[jclono2_a]                   CLONO2 + hv -> CL + NO3 
[jclono2_b]                   CLONO2 + hv -> CLO + NO2 
[jcof2]                       COF2 + hv -> 2*F 
[jcofcl]                      COFCL + hv -> F + CL 
[jh2402]                      H2402 + hv -> 2*BR + 2*COF2 
[jhbr]                        HBR + hv -> BR + H 
[jhcfc141b]                   HCFC141B + hv -> CL + COFCL 
[jhcfc142b]                   HCFC142B + hv -> CL + COF2 
[jhcfc22]                     HCFC22 + hv -> CL + COF2 
[jhcl]                        HCL + hv -> H + CL 
[jhf]                         HF + hv -> H + F 
[jhobr]                       HOBR + hv -> BR + OH 
[jhocl]                       HOCL + hv -> OH + CL 
[joclo]                       OCLO + hv -> O + CLO 
[jsf6]                        SF6 + hv -> sink 
*********************************
*** sulfur
*********************************
[jh2so4]                      H2SO4 + hv -> SO3 + H2O 
[jocs]                        OCS + hv -> S + CO 
[jso]                         SO + hv -> S + O 
[jso2]                        SO2 + hv -> SO + O 
[jso3]                        SO3 + hv -> SO2 + O 
*********************************
*** soa
*********************************
[jsoa1_a1->,.0004*jno2]       soa1_a1 + hv ->  
[jsoa1_a2->,.0004*jno2]       soa1_a2 + hv ->  
[jsoa2_a1->,.0004*jno2]       soa2_a1 + hv ->  
[jsoa2_a2->,.0004*jno2]       soa2_a2 + hv ->  
[jsoa3_a1->,.0004*jno2]       soa3_a1 + hv ->  
[jsoa3_a2->,.0004*jno2]       soa3_a2 + hv ->  
[jsoa4_a1->,.0004*jno2]       soa4_a1 + hv ->  
[jsoa4_a2->,.0004*jno2]       soa4_a2 + hv ->  
[jsoa5_a1->,.0004*jno2]       soa5_a1 + hv ->  
[jsoa5_a2->,.0004*jno2]       soa5_a2 + hv ->  
      End Photolysis

      Reactions
*********************************
*** Not Assigned to a Section
*********************************
[E90_tau]              E90  ->                                                  ; 1.29e-07 
*********************************
*** odd-oxygen
*********************************
[O1D_H2]               O1D + H2  -> H + OH                                      ; 1.2e-10 
[O1D_H2O]              O1D + H2O  -> 2*OH                                       ; 1.63e-10, 60 
[O1D_N2,cph=189.81]    O1D + N2  -> O + N2                                      ; 2.15e-11, 110 
[O1D_O2ab]             O1D + O2  -> O + O2                                      ; 3.3e-11, 55 
[O1D_O3]               O1D + O3  -> O2 + O2                                     ; 1.2e-10 
[O1D_O3a]              O1D + O3  -> O2 + 2*O                                    ; 1.2e-10 
[O_O3,cph=392.19]      O + O3  -> 2*O2                                          ; 8e-12, -2060 
[usr_O_O,cph=493.58]   O + O + M  -> O2 + M                                      
[usr_O_O2,cph=101.39]  O + O2 + M  -> O3 + M                                     
*********************************
*** odd-hydrogen
*********************************
[H2_O]                 H2 + O  -> OH + H                                        ; 1.6e-11, -4570 
[H2O2_O]               H2O2 + O  -> OH + HO2                                    ; 1.4e-12, -2000 
[H_HO2,cph=232.59]     H + HO2  -> H2 + O2                                      ; 6.9e-12 
[H_HO2a]               H + HO2  -> 2*OH                                         ; 7.2e-11 
[H_HO2b]               H + HO2  -> H2O + O                                      ; 1.6e-12 
[H_O2,cph=203.4]       H + O2 + M  -> HO2 + M                                   ; 5.3e-32, 1.8, 9.5e-11, -0.4, 0.6 
[HO2_O,cph=226.58]     HO2 + O  -> OH + O2                                      ; 3e-11, 200 
[HO2_O3,cph=120.1]     HO2 + O3  -> OH + 2*O2                                   ; 1e-14, -490 
[H_O3,cph=194.71]      H + O3  -> OH + O2                                       ; 1.4e-10, -470 
[OH_H2]                OH + H2  -> H2O + H                                      ; 2.8e-12, -1800 
[OH_H2O2]              OH + H2O2  -> H2O + HO2                                  ; 1.8e-12 
[OH_HO2,cph=293.62]    OH + HO2  -> H2O + O2                                    ; 4.8e-11, 250 
[OH_O,cph=67.67]       OH + O  -> H + O2                                        ; 1.8e-11, 180 
[OH_O3,cph=165.3]      OH + O3  -> HO2 + O2                                     ; 1.7e-12, -940 
[OH_OH]                OH + OH  -> H2O + O                                      ; 1.8e-12 
[OH_OH_M]              OH + OH + M  -> H2O2 + M                                 ; 6.9e-31, 1, 2.6e-11, 0, 0.6 
[usr_HO2_HO2,cph=165.51] HO2 + HO2  -> H2O2 + O2                                 
*********************************
*** odd-nitrogen
*********************************
[HO2NO2_OH]            HO2NO2 + OH  -> H2O + NO2 + O2                           ; 4.5e-13, 610 
[N_NO,cph=313.75]      N + NO  -> N2 + O                                        ; 2.1e-11, 100 
[N_NO2a]               N + NO2  -> N2O + O                                      ; 2.9e-12, 220 
[N_NO2b]               N + NO2  -> 2*NO                                         ; 1.45e-12, 220 
[N_NO2c]               N + NO2  -> N2 + O2                                      ; 1.45e-12, 220 
[N_O2,cph=133.75]      N + O2  -> NO + O                                        ; 3.3e-12, -3150 
[NO2_O,cph=193.02]     NO2 + O  -> NO + O2                                      ; 5.1e-12, 210 
[NO2_O3]               NO2 + O3  -> NO3 + O2                                    ; 1.2e-13, -2450 
[NO2_O_M]              NO2 + O + M  -> NO3 + M                                  ; 2.5e-31, 1.8, 2.2e-11, 0.7, 0.6 
[NO3_HO2]              NO3 + HO2  -> OH + NO2 + O2                              ; 3.5e-12 
[NO3_NO]               NO3 + NO  -> 2*NO2                                       ; 1.7e-11, 125 
[NO3_O]                NO3 + O  -> NO2 + O2                                     ; 1.3e-11 
[NO3_OH]               NO3 + OH  -> HO2 + NO2                                   ; 2.2e-11 
[N_OH]                 N + OH  -> NO + H                                        ; 5e-11 
[NO_HO2,cph=34.47]     NO + HO2  -> NO2 + OH                                    ; 3.44e-12, 260 
[NO_O3,cph=199.17]     NO + O3  -> NO2 + O2                                     ; 3e-12, -1500 
[NO_O_M]               NO + O + M  -> NO2 + M                                   ; 9e-32, 1.5, 3e-11, 0, 0.6 
[O1D_N2Oa]             O1D + N2O  -> 2*NO                                       ; 7.26e-11, 20 
[O1D_N2Ob]             O1D + N2O  -> N2 + O2                                    ; 4.64e-11, 20 
[tag_NO2_HO2]          NO2 + HO2 + M  -> HO2NO2 + M                             ; 1.9e-31, 3.4, 4e-12, 0.3, 0.6 
[tag_NO2_NO3]          NO2 + NO3 + M  -> N2O5 + M                               ; 2.4e-30, 3, 1.6e-12, -0.1, 0.6 
[tag_NO2_OH]           NO2 + OH + M  -> HNO3 + M                                ; 1.8e-30, 3, 2.8e-11, 0, 0.6 
[usr_HNO3_OH]          HNO3 + OH  -> NO3 + H2O                                   
[usr_HO2NO2_M]         HO2NO2 + M  -> HO2 + NO2 + M                              
[usr_N2O5_M]           N2O5 + M  -> NO2 + NO3 + M                                
*********************************
*** odd-chlorine
*********************************
[CL_CH2O]              CL + CH2O  -> HCL + HO2 + CO                             ; 8.1e-11, -30 
[CL_CH4]               CL + CH4  -> CH3O2 + HCL                                 ; 7.1e-12, -1270 
[CL_H2]                CL + H2  -> HCL + H                                      ; 3.05e-11, -2270 
[CL_H2O2]              CL + H2O2  -> HCL + HO2                                  ; 1.1e-11, -980 
[CL_HO2a]              CL + HO2  -> HCL + O2                                    ; 1.4e-11, 270 
[CL_HO2b]              CL + HO2  -> OH + CLO                                    ; 3.6e-11, -375 
[CL_O3]                CL + O3  -> CLO + O2                                     ; 2.3e-11, -200 
[CLO_CH3O2]            CLO + CH3O2  -> CL + HO2 + CH2O                          ; 3.3e-12, -115 
[CLO_CLOa]             CLO + CLO  -> 2*CL + O2                                  ; 3e-11, -2450 
[CLO_CLOb]             CLO + CLO  -> CL2 + O2                                   ; 1e-12, -1590 
[CLO_CLOc]             CLO + CLO  -> CL + OCLO                                  ; 3.5e-13, -1370 
[CLO_HO2]              CLO + HO2  -> O2 + HOCL                                  ; 2.6e-12, 290 
[CLO_NO]               CLO + NO  -> NO2 + CL                                    ; 6.4e-12, 290 
[CLONO2_CL]            CLONO2 + CL  -> CL2 + NO3                                ; 6.5e-12, 135 
[CLO_NO2_M]            CLO + NO2 + M  -> CLONO2 + M                             ; 1.8e-31, 3.4, 1.5e-11, 1.9, 0.6 
[CLONO2_O]             CLONO2 + O  -> CLO + NO3                                 ; 3.6e-12, -840 
[CLONO2_OH]            CLONO2 + OH  -> HOCL + NO3                               ; 1.2e-12, -330 
[CLO_O]                CLO + O  -> CL + O2                                      ; 2.8e-11, 85 
[CLO_OHa]              CLO + OH  -> CL + HO2                                    ; 7.4e-12, 270 
[CLO_OHb]              CLO + OH  -> HCL + O2                                    ; 6e-13, 230 
[HCL_O]                HCL + O  -> CL + OH                                      ; 1e-11, -3300 
[HCL_OH]               HCL + OH  -> H2O + CL                                    ; 1.8e-12, -250 
[HOCL_CL]              HOCL + CL  -> HCL + CLO                                  ; 3.4e-12, -130 
[HOCL_O]               HOCL + O  -> CLO + OH                                    ; 1.7e-13 
[HOCL_OH]              HOCL + OH  -> H2O + CLO                                  ; 3e-12, -500 
[O1D_CCL4]             O1D + CCL4  -> 4*CL                                      ; 2.607e-10 
[O1D_CF2CLBR]          O1D + CF2CLBR  -> CL + BR + COF2                         ; 9.75e-11 
[O1D_CFC11]            O1D + CFC11  -> 2*CL + COFCL                             ; 2.07e-10 
[O1D_CFC113]           O1D + CFC113  -> 2*CL + COFCL + COF2                     ; 2.088e-10 
[O1D_CFC114]           O1D + CFC114  -> 2*CL + 2*COF2                           ; 1.17e-10 
[O1D_CFC115]           O1D + CFC115  -> CL + F + 2*COF2                         ; 4.644e-11 
[O1D_CFC12]            O1D + CFC12  -> 2*CL + COF2                              ; 1.204e-10 
[O1D_HCLa]             O1D + HCL  -> CL + OH                                    ; 9.9e-11 
[O1D_HCLb]             O1D + HCL  -> CLO + H                                    ; 3.3e-12 
[tag_CLO_CLO_M]        CLO + CLO + M  -> CL2O2 + M                              ; 1.9e-32, 3.6, 3.7e-12, 1.6, 0.6 
[usr_CL2O2_M]          CL2O2 + M  -> CLO + CLO + M                               
*********************************
*** odd-bromine
*********************************
[BR_CH2O]              BR + CH2O  -> HBR + HO2 + CO                             ; 1.7e-11, -800 
[BR_HO2]               BR + HO2  -> HBR + O2                                    ; 4.8e-12, -310 
[BR_O3]                BR + O3  -> BRO + O2                                     ; 1.6e-11, -780 
[BRO_BRO]              BRO + BRO  -> 2*BR + O2                                  ; 1.5e-12, 230 
[BRO_CLOa]             BRO + CLO  -> BR + OCLO                                  ; 9.5e-13, 550 
[BRO_CLOb]             BRO + CLO  -> BR + CL + O2                               ; 2.3e-12, 260 
[BRO_CLOc]             BRO + CLO  -> BRCL + O2                                  ; 4.1e-13, 290 
[BRO_HO2]              BRO + HO2  -> HOBR + O2                                  ; 4.5e-12, 460 
[BRO_NO]               BRO + NO  -> BR + NO2                                    ; 8.8e-12, 260 
[BRO_NO2_M]            BRO + NO2 + M  -> BRONO2 + M                             ; 5.2e-31, 3.2, 6.9e-12, 2.9, 0.6 
[BRONO2_O]             BRONO2 + O  -> BRO + NO3                                 ; 1.9e-11, 215 
[BRO_O]                BRO + O  -> BR + O2                                      ; 1.9e-11, 230 
[BRO_OH]               BRO + OH  -> BR + HO2                                    ; 1.7e-11, 250 
[HBR_O]                HBR + O  -> BR + OH                                      ; 5.8e-12, -1500 
[HBR_OH]               HBR + OH  -> BR + H2O                                    ; 5.5e-12, 200 
[HOBR_O]               HOBR + O  -> BRO + OH                                    ; 1.2e-10, -430 
[O1D_CF3BR]            O1D + CF3BR  -> BR + F + COF2                            ; 4.5e-11 
[O1D_CHBR3]            O1D + CHBR3  -> 3*BR                                     ; 4.62e-10 
[O1D_H2402]            O1D + H2402  -> 2*BR + 2*COF2                            ; 1.2e-10 
[O1D_HBRa]             O1D + HBR  -> BR + OH                                    ; 9e-11 
[O1D_HBRb]             O1D + HBR  -> BRO + H                                    ; 3e-11 
*********************************
*** odd-fluorine
*********************************
[F_CH4]                F + CH4  -> HF + CH3O2                                   ; 1.6e-10, -260 
[F_H2]                 F + H2  -> HF + H                                        ; 1.4e-10, -500 
[F_H2O]                F + H2O  -> HF + OH                                      ; 1.4e-11, 0 
[F_HNO3]               F + HNO3  -> HF + NO3                                    ; 6e-12, 400 
[O1D_COF2]             O1D + COF2  -> 2*F                                       ; 2.14e-11 
[O1D_COFCL]            O1D + COFCL  -> F + CL                                   ; 1.9e-10 
*********************************
*** organic-halogens
*********************************
[CH2BR2_CL]            CH2BR2 + CL  -> 2*BR + HCL                               ; 6.3e-12, -800 
[CH2BR2_OH]            CH2BR2 + OH  -> 2*BR + H2O                               ; 2e-12, -840 
[CH3BR_CL]             CH3BR + CL  -> HCL + HO2 + BR                            ; 1.46e-11, -1040 
[CH3BR_OH]             CH3BR + OH  -> BR + H2O + HO2                            ; 1.42e-12, -1150 
[CH3CCL3_OH]           CH3CCL3 + OH  -> H2O + 3*CL                              ; 1.64e-12, -1520 
[CH3CL_CL]             CH3CL + CL  -> HO2 + CO + 2*HCL                          ; 2.03e-11, -1100 
[CH3CL_OH]             CH3CL + OH  -> CL + H2O + HO2                            ; 1.96e-12, -1200 
[CHBR3_CL]             CHBR3 + CL  -> 3*BR + HCL                                ; 4.85e-12, -850 
[CHBR3_OH]             CHBR3 + OH  -> 3*BR                                      ; 9e-13, -360 
[HCFC141B_OH]          HCFC141B + OH  -> CL + COFCL                             ; 1.25e-12, -1600 
[HCFC142B_OH]          HCFC142B + OH  -> CL + COF2                              ; 1.3e-12, -1770 
[HCFC22_OH]            HCFC22 + OH  -> H2O + CL + COF2                          ; 9.2e-13, -1560 
[O1D_CH2BR2]           O1D + CH2BR2  -> 2*BR                                    ; 2.57e-10 
[O1D_CH3BR]            O1D + CH3BR  -> BR                                       ; 1.8e-10 
[O1D_HCFC141B]         O1D + HCFC141B  -> CL + COFCL                            ; 1.794e-10 
[O1D_HCFC142B]         O1D + HCFC142B  -> CL + COF2                             ; 1.3e-10 
[O1D_HCFC22]           O1D + HCFC22  -> CL + COF2                               ; 7.65e-11 
*********************************
*** C1
*********************************
[CH2O_HO2]             CH2O + HO2  -> HOCH2OO                                   ; 9.7e-15, 625 
[CH2O_NO3]             CH2O + NO3  -> CO + HO2 + HNO3                           ; 6e-13, -2058 
[CH2O_O]               CH2O + O  -> HO2 + OH + CO                               ; 3.4e-11, -1600 
[CH2O_OH]              CH2O + OH  -> CO + H2O + H                               ; 5.5e-12, 125 
[CH3O2_CH3O2a]         CH3O2 + CH3O2  -> 2*CH2O + 2*HO2                         ; 5e-13, -424 
[CH3O2_CH3O2b]         CH3O2 + CH3O2  -> CH2O + CH3OH                           ; 1.9e-14, 706 
[CH3O2_HO2]            CH3O2 + HO2  -> CH3OOH + O2                              ; 4.1e-13, 750 
[CH3O2_NO]             CH3O2 + NO  -> CH2O + NO2 + HO2                          ; 2.8e-12, 300 
[CH3OH_OH]             CH3OH + OH  -> HO2 + CH2O                                ; 2.9e-12, -345 
[CH3OOH_OH]            CH3OOH + OH  -> 0.7*CH3O2 + 0.3*OH + 0.3*CH2O + H2O      ; 3.8e-12, 200 
[CH4_OH]               CH4 + OH  -> CH3O2 + H2O                                 ; 2.45e-12, -1775 
[HCN_OH]               HCN + OH + M  -> HO2 + M                                 ; 6.1e-33, 1.5, 9.8e-15, -4.6, 0.8 
[HCOOH_OH]             HCOOH + OH  -> HO2 + CO2 + H2O                           ; 4e-13 
[HOCH2OO_HO2]          HOCH2OO + HO2  -> HCOOH                                  ; 7.5e-13, 700 
[HOCH2OO_M]            HOCH2OO  -> CH2O + HO2                                   ; 2.4e+12, -7000 
[HOCH2OO_NO]           HOCH2OO + NO  -> HCOOH + NO2 + HO2                       ; 2.6e-12, 265 
[O1D_CH4a]             O1D + CH4  -> CH3O2 + OH                                 ; 1.31e-10 
[O1D_CH4b]             O1D + CH4  -> CH2O + H + HO2                             ; 3.5e-11 
[O1D_CH4c]             O1D + CH4  -> CH2O + H2                                  ; 9e-12 
[O1D_HCN]              O1D + HCN  -> OH                                         ; 1.08e-10, 105 
[usr_CO_OH]            CO + OH  -> CO2 + HO2                                     
*********************************
*** C2
*********************************
[C2H2_CL_M]            C2H2 + CL + M  -> CL + M                                 ; 5.2e-30, 2.4, 2.2e-10, 0.7, 0.6 
[C2H2_OH_M]            C2H2 + OH + M  -> 0.65*GLYOXAL + 0.65*OH + 0.35*HCOOH + 0.35*HO2 + 0.35*CO + M ; 5.5e-30, 0, 8.3e-13, -2, 0.6 
[C2H4_CL_M]            C2H4 + CL + M  -> CL + M                                 ; 1.6e-29, 3.3, 3.1e-10, 1, 0.6 
[C2H4_O3]              C2H4 + O3  -> 0.63*CO + 0.13*OH + 0.13*HO2 + 0.37*HCOOH + CH2O ; 1.2e-14, -2630 
[C2H5O2_C2H5O2]        C2H5O2 + C2H5O2  -> 1.6*CH3CHO + 1.2*HO2 + 0.4*C2H5OH    ; 6.8e-14 
[C2H5O2_CH3O2]         C2H5O2 + CH3O2  -> 0.7*CH2O + 0.8*CH3CHO + HO2 + 0.3*CH3OH + 0.2*C2H5OH ; 2e-13 
[C2H5O2_HO2]           C2H5O2 + HO2  -> C2H5OOH + O2                            ; 7.5e-13, 700 
[C2H5O2_NO]            C2H5O2 + NO  -> CH3CHO + HO2 + NO2                       ; 2.6e-12, 365 
[C2H5OH_OH]            C2H5OH + OH  -> HO2 + CH3CHO                             ; 6.9e-12, -230 
[C2H5OOH_OH]           C2H5OOH + OH  -> 0.5*C2H5O2 + 0.5*CH3CHO + 0.5*OH        ; 3.8e-12, 200 
[C2H6_CL]              C2H6 + CL  -> HCL + C2H5O2                               ; 7.2e-11, -70 
[C2H6_OH]              C2H6 + OH  -> C2H5O2 + H2O                               ; 7.66e-12, -1020 
[CH3CHO_NO3]           CH3CHO + NO3  -> CH3CO3 + HNO3                           ; 1.4e-12, -1900 
[CH3CHO_OH]            CH3CHO + OH  -> CH3CO3 + H2O                             ; 4.63e-12, 350 
[CH3CN_OH]             CH3CN + OH  -> HO2                                       ; 7.8e-13, -1050 
[CH3CO3_CH3CO3]        CH3CO3 + CH3CO3  -> 2*CH3O2 + 2*CO2                      ; 2.9e-12, 500 
[CH3CO3_CH3O2]         CH3CO3 + CH3O2  -> 0.9*CH3O2 + CH2O + 0.9*HO2 + 0.9*CO2 + 0.1*CH3COOH ; 2e-12, 500 
[CH3CO3_HO2]           CH3CO3 + HO2  -> 0.4*CH3COOOH + 0.15*CH3COOH + 0.15*O3 + 0.45*OH + 0.45*CH3O2 ; 4.3e-13, 1040 
[CH3CO3_NO]            CH3CO3 + NO  -> CH3O2 + CO2 + NO2                        ; 8.1e-12, 270 
[CH3COOH_OH]           CH3COOH + OH  -> CH3O2 + CO2 + H2O                       ; 3.15e-14, 920 
[CH3COOOH_OH]          CH3COOOH + OH  -> 0.5*CH3CO3 + 0.5*CH2O + 0.5*CO2 + H2O  ; 1e-12 
[EO2_HO2]              EO2 + HO2  -> EOOH                                       ; 7.5e-13, 700 
[EO2_NO]               EO2 + NO  -> 0.5*CH2O + 0.25*HO2 + 0.75*EO + NO2         ; 4.2e-12, 180 
[EO_M]                 EO  -> 2*CH2O + HO2                                      ; 1.6e+11, -4150 
[EO_O2]                EO + O2  -> GLYALD + HO2                                 ; 1e-14 
[GLYALD_OH]            GLYALD + OH  -> HO2 + 0.2*GLYOXAL + 0.8*CH2O + 0.8*CO2   ; 1e-11 
[GLYOXAL_OH]           GLYOXAL + OH  -> HO2 + CO + CO2                          ; 1.15e-11 
[PAN_OH]               PAN + OH  -> CH2O + NO3                                  ; 4e-14 
[tag_C2H4_OH]          C2H4 + OH + M  -> EO2 + M                                ; 8.6e-29, 3.1, 9e-12, 0.85, 0.48 
[tag_CH3CO3_NO2]       CH3CO3 + NO2 + M  -> PAN + M                             ; 7.3e-29, 4.1, 9.5e-12, 1.6, 0.6 
[usr_PAN_M]            PAN + M  -> CH3CO3 + NO2 + M                              
*********************************
*** C3
*********************************
[C3H6_NO3]             C3H6 + NO3  -> NOA                                       ; 4.6e-13, -1156 
[C3H6_O3]              C3H6 + O3  -> 0.5*CH2O + 0.12*HCOOH + 0.12*CH3COOH + 0.5*CH3CHO + 0.56*CO + 0.28*CH3O2 + 0.1*CH4 + 0.2*CO2 + 0.28*HO2 + 0.36*OH ; 6.5e-15, -1900 
[C3H7O2_CH3O2]         C3H7O2 + CH3O2  -> CH2O + HO2 + 0.82*CH3COCH3            ; 3.75e-13, -40 
[C3H7O2_HO2]           C3H7O2 + HO2  -> C3H7OOH + O2                            ; 7.5e-13, 700 
[C3H7O2_NO]            C3H7O2 + NO  -> 0.82*CH3COCH3 + NO2 + HO2 + 0.27*CH3CHO  ; 4.2e-12, 180 
[C3H7OOH_OH]           C3H7OOH + OH  -> H2O + C3H7O2                            ; 3.8e-12, 200 
[C3H8_OH]              C3H8 + OH  -> C3H7O2 + H2O                               ; 9.19e-12, -630 
[CH3COCHO_NO3]         CH3COCHO + NO3  -> HNO3 + CO + CH3CO3                    ; 1.4e-12, -1860 
[CH3COCHO_OH]          CH3COCHO + OH  -> CH3CO3 + CO + H2O                      ; 8.4e-13, 830 
[HYAC_OH]              HYAC + OH  -> CH3COCHO + HO2                             ; 3e-12 
[NOA_OH]               NOA + OH  -> NO2 + CH3COCHO                              ; 6.7e-13 
[PO2_HO2]              PO2 + HO2  -> POOH + O2                                  ; 7.5e-13, 700 
[PO2_NO]               PO2 + NO  -> CH3CHO + CH2O + HO2 + NO2                   ; 4.2e-12, 180 
[POOH_OH]              POOH + OH  -> 0.5*PO2 + 0.5*OH + 0.5*HYAC + H2O          ; 3.8e-12, 200 
[RO2_CH3O2]            RO2 + CH3O2  -> 0.3*CH3CO3 + 0.8*CH2O + 0.3*HO2 + 0.2*HYAC + 0.5*CH3COCHO + 0.5*CH3OH ; 7.1e-13, 500 
[RO2_HO2]              RO2 + HO2  -> 0.85*ROOH + 0.15*OH + 0.15*CH2O + 0.15*CH3CO3 ; 8.6e-13, 700 
[RO2_NO]               RO2 + NO  -> CH3CO3 + CH2O + NO2                         ; 2.9e-12, 300 
[ROOH_OH]              ROOH + OH  -> RO2 + H2O                                  ; 3.8e-12, 200 
[tag_C3H6_OH]          C3H6 + OH + M  -> PO2 + M                                ; 8e-27, 3.5, 3e-11, 0, 0.5 
[usr_CH3COCH3_OH]      CH3COCH3 + OH  -> RO2 + H2O                               
*********************************
*** C4
*********************************
[BIGENE_NO3]           BIGENE + NO3  -> NO2 + CH3CHO + 0.5*CH2O + 0.5*CH3COCH3  ; 3.5e-13 
[BIGENE_OH]            BIGENE + OH  -> ENEO2                                    ; 5.4e-11 
[ENEO2_NO]             ENEO2 + NO  -> CH3CHO + 0.5*CH2O + 0.5*CH3COCH3 + HO2 + NO2 ; 4.8e-12, 120 
[ENEO2_NOb]            ENEO2 + NO  -> HONITR                                    ; 5.1e-14, 693 
[HONITR_OH]            HONITR + OH  -> ONITR + HO2                              ; 2e-12 
[MACRO2_CH3CO3]        MACRO2 + CH3CO3  -> 0.25*CH3COCHO + CH3O2 + 0.22*CO + 0.47*HO2 + 0.53*GLYALD + 0.22*HYAC + 0.25*CH2O + 0.53*CH3CO3 ; 1.4e-11 
[MACRO2_CH3O2]         MACRO2 + CH3O2  -> 0.73*HO2 + 0.88*CH2O + 0.11*CO + 0.24*CH3COCHO + 0.26*GLYALD + 0.26*CH3CO3 + 0.25*CH3OH + 0.23*HYAC ; 5e-13, 400 
[MACRO2_HO2]           MACRO2 + HO2  -> MACROOH                                 ; 8e-13, 700 
[MACRO2_NO3]           MACRO2 + NO3  -> NO2 + 0.47*HO2 + 0.25*CH2O + 0.25*CH3COCHO + 0.22*CO + 0.53*GLYALD + 0.22*HYAC + 0.53*CH3CO3 ; 2.4e-12 
[MACRO2_NOa]           MACRO2 + NO  -> NO2 + 0.47*HO2 + 0.25*CH2O + 0.53*GLYALD + 0.25*CH3COCHO + 0.53*CH3CO3 + 0.22*HYAC + 0.22*CO ; 2.7e-12, 360 
[MACRO2_NOb]           MACRO2 + NO  -> HONITR                                   ; 1.3e-13, 360 
[MACR_O3]              MACR + O3  -> 0.12*CH2O + 0.24*OH + 0.65*CO + 0.1*CH3CO3 + 0.88*CH3COCHO + 0.33*HCOOH + 0.14*HO2 ; 1.5e-15, -2100 
[MACR_OH]              MACR + OH  -> 0.5*MACRO2 + 0.5*H2O + 0.5*MCO3            ; 9.6e-12, 360 
[MACROOH_OH]           MACROOH + OH  -> 0.5*MCO3 + 0.2*MACRO2 + 0.1*OH + 0.2*HO2 ; 2.3e-11, 200 
[MCO3_CH3CO3]          MCO3 + CH3CO3  -> 2*CO2 + CH3O2 + CH2O + CH3CO3          ; 4.6e-12, 530 
[MCO3_CH3O2]           MCO3 + CH3O2  -> 2*CH2O + HO2 + CO2 + CH3CO3             ; 2e-12, 500 
[MCO3_HO2]             MCO3 + HO2  -> 0.15*O3 + 0.15*CH3COOH + 0.4*CH3COOOH + 0.45*OH + 0.45*CO2 + 0.45*CH2O + 0.45*CH3CO3 ; 4.3e-13, 1040 
[MCO3_MCO3]            MCO3 + MCO3  -> 2*CO2 + 2*CH2O + 2*CH3CO3                ; 2.3e-12, 530 
[MCO3_NO]              MCO3 + NO  -> NO2 + CH2O + CH3CO3                        ; 5.3e-12, 360 
[MCO3_NO3]             MCO3 + NO3  -> NO2 + CH2O + CH3CO3                       ; 5e-12 
[MEKO2_HO2]            MEKO2 + HO2  -> 0.8*MEKOOH + 0.2*OH + 0.2*CH3CHO + 0.2*CH3CO3 ; 7.5e-13, 700 
[MEKO2_NO]             MEKO2 + NO  -> CH3CO3 + CH3CHO + NO2                     ; 4.2e-12, 180 
[MEK_OH]               MEK + OH  -> MEKO2                                       ; 2.3e-12, -170 
[MEKOOH_OH]            MEKOOH + OH  -> MEKO2                                    ; 3.8e-12, 200 
[MPAN_OH_M]            MPAN + OH + M  -> 0.5*HYAC + 0.5*NO3 + 0.5*CH2O + 0.5*HO2 + 0.5*CO2 + M + 0.5*NDEP ; 8e-27, 3.5, 3e-11, 0, 0.5 
[MVK_O3]               MVK + O3  -> 0.6*CH2O + 0.56*CO + 0.1*CH3CHO + 0.1*CO2 + 0.28*CH3CO3 + 0.5*CH3COCHO + 0.28*HO2 + 0.36*OH + 0.12*HCOOH ; 8.5e-16, -1520 
[MVK_OH]               MVK + OH  -> MACRO2                                      ; 4.13e-12, 452 
[tag_MCO3_NO2]         MCO3 + NO2 + M  -> MPAN + M                              ; 9.7e-29, 5.6, 9.3e-12, 1.5, 0.6 
[usr_MPAN_M]           MPAN + M  -> MCO3 + NO2 + M                               
*********************************
*** C5
*********************************
[ALKNIT_OH]            ALKNIT + OH  -> 0.4*CH2O + 0.8*CH3CHO + 0.8*CH3COCH3 + NO2 ; 1.6e-12 
[ALKO2_HO2]            ALKO2 + HO2  -> ALKOOH                                   ; 7.5e-13, 700 
[ALKO2_NO]             ALKO2 + NO  -> 0.4*CH3CHO + 0.1*CH2O + 0.25*CH3COCH3 + HO2 + 0.8*MEK + NO2 ; 6.7e-12 
[ALKO2_NOb]            ALKO2 + NO  -> ALKNIT                                    ; 5.4e-14, 870 
[ALKOOH_OH]            ALKOOH + OH  -> ALKO2                                    ; 3.8e-12, 200 
[BIGALK_OH]            BIGALK + OH  -> ALKO2                                    ; 3.5e-12 
[HPALD_OH]             HPALD + OH  -> XO2                                       ; 1.86e-11, 175 
[HYDRALD_OH]           HYDRALD + OH  -> XO2                                     ; 1.86e-11, 175 
[IEPOX_OH]             IEPOX + OH  -> XO2                                       ; 1.3e-11 
[ISOPAO2_CH3CO3]       ISOPAO2 + CH3CO3  -> CH3O2 + HO2 + CH2O + 0.39*MACR + 0.61*MVK + CO2 ; 1.4e-11 
[ISOPAO2_CH3O2]        ISOPAO2 + CH3O2  -> 0.25*CH3OH + HO2 + 1.5*CH2O + 0.31*MACR + 0.44*MVK ; 5e-13, 400 
[ISOPAO2_HO2]          ISOPAO2 + HO2  -> ISOPOOH                                ; 8e-13, 700 
[ISOPAO2_NO]           ISOPAO2 + NO  -> 0.08*ISOPNITA + 0.92*NO2 + 0.36*MACR + 0.56*MVK + 0.92*CH2O + 0.92*HO2 ; 4.4e-12, 180 
[ISOPAO2_NO3]          ISOPAO2 + NO3  -> NO2 + 0.4*MACR + 0.6*MVK + CH2O + HO2  ; 2.4e-12 
[ISOPBO2_CH3CO3]       ISOPBO2 + CH3CO3  -> HYDRALD + CH3O2 + HO2               ; 1.4e-11 
[ISOPBO2_CH3O2]        ISOPBO2 + CH3O2  -> 0.25*CH3OH + HO2 + 0.75*CH2O + 0.75*HYDRALD ; 5e-13, 400 
[ISOPBO2_HO2]          ISOPBO2 + HO2  -> ISOPOOH                                ; 8e-13, 700 
[ISOPBO2_M]            ISOPBO2  -> HPALD + HO2                                  ; 1.6e+09, -8300 
[ISOPBO2_NO]           ISOPBO2 + NO  -> 0.87*HYDRALD + 0.08*ISOPNITB + 0.92*NO2 + 0.92*HO2 + 0.05*GLYOXAL + 0.05*GLYALD + 0.05*CH3COCHO + 0.05*HYAC ; 4.4e-12, 180 
[ISOPBO2_NO3]          ISOPBO2 + NO3  -> NO2 + 0.95*HYDRALD + HO2 + 0.05*GLYOXAL + 0.05*GLYALD + 0.05*CH3COCHO + 0.05*HYAC ; 2.4e-12 
[ISOPNITA_OH]          ISOPNITA + OH  -> 0.7*HYAC + 0.7*GLYALD + 0.7*NO2 + 0.3*CH2O + 0.3*HONITR + 0.3*HO2 ; 4e-11 
[ISOPNITB_OH]          ISOPNITB + OH  -> 0.5*HYAC + 0.5*GLYALD + 0.5*NOA + HO2 + 0.5*HONITR ; 4e-11 
[ISOP_NO3]             ISOP + NO3  -> ISOPNO3                                   ; 3.03e-12, -446 
[ISOPNO3_CH3CO3]       ISOPNO3 + CH3CO3  -> NC4CHO + CH3O2 + HO2                ; 1.4e-11 
[ISOPNO3_CH3O2]        ISOPNO3 + CH3O2  -> 0.8*NC4CHO + 1.2*HO2 + 0.8*CH2O + 0.2*CH3OH + 0.2*NC4CH2OH ; 5e-13, 400 
[ISOPNO3_HO2]          ISOPNO3 + HO2  -> ISOPNOOH                               ; 8e-13, 700 
[ISOPNO3_NO]           ISOPNO3 + NO  -> NC4CHO + NO2 + HO2                      ; 2.7e-12, 360 
[ISOPNO3_NO3]          ISOPNO3 + NO3  -> NC4CHO + NO2 + HO2                     ; 2.4e-12 
[ISOPNOOH_OH]          ISOPNOOH + OH  -> NOA + HO2                              ; 4e-11 
[ISOP_O3]              ISOP + O3  -> 0.3*MACR + 0.2*MVK + 0.11*HCOOH + 0.62*CO + 0.32*OH + 0.37*HO2 + 0.91*CH2O + 0.08*CH3CO3 + 0.13*C3H6 + 0.05*CH3O2 ; 1.05e-14, -2000 
[ISOP_OH]              ISOP + OH  -> 0.6*ISOPAO2 + 0.4*ISOPBO2                  ; 2.54e-11, 410 
[ISOPOOH_OH]           ISOPOOH + OH  -> 0.4*XO2 + 0.6*IEPOX + 0.6*OH            ; 1.52e-11, 200 
[NC4CH2OH_OH]          NC4CH2OH + OH  -> GLYALD + NOA + HO2                     ; 7e-11 
[NC4CHO_OH]            NC4CHO + OH  -> GLYOXAL + NOA + HO2                      ; 1e-10 
[XO2_CH3CO3]           XO2 + CH3CO3  -> 0.25*CO + 0.25*CH2O + 0.25*GLYOXAL + CH3O2 + HO2 + 0.25*CH3COCHO + 0.25*HYAC + 0.25*GLYALD + CO2 ; 1.3e-12, 640 
[XO2_CH3O2]            XO2 + CH3O2  -> 0.3*CH3OH + 0.8*HO2 + 0.8*CH2O + 0.2*CO + 0.1*GLYOXAL + 0.1*CH3COCHO + 0.1*HYAC + 0.1*GLYALD ; 5e-13, 400 
[XO2_HO2]              XO2 + HO2  -> XOOH                                       ; 8e-13, 700 
[XO2_NO]               XO2 + NO  -> NO2 + HO2 + 0.25*CO + 0.25*CH2O + 0.25*GLYOXAL + 0.25*CH3COCHO + 0.25*HYAC + 0.25*GLYALD ; 2.7e-12, 360 
[XO2_NO3]              XO2 + NO3  -> NO2 + HO2 + 0.5*CO + 0.25*HYAC + 0.25*GLYOXAL + 0.25*CH3COCHO + 0.25*GLYALD ; 2.4e-12 
[XOOH_OH]              XOOH + OH  -> 0.5*XO2 + 0.5*OH                           ; 1.52e-12, 200 
*********************************
*** C7
*********************************
[ACBZO2_HO2]           ACBZO2 + HO2  -> 0.4*C6H5O2 + 0.4*OH                     ; 4.3e-13, 1040 
[ACBZO2_NO]            ACBZO2 + NO  -> C6H5O2 + NO2                             ; 7.5e-12, 290 
[BENZENE_OH]           BENZENE + OH  -> 0.53*PHENOL + 0.12*BEPOMUC + 0.65*HO2 + 0.35*BENZO2 ; 2.3e-12, -193 
[BENZO2_HO2]           BENZO2 + HO2  -> BENZOOH                                 ; 7.5e-13, 700 
[BENZO2_NO]            BENZO2 + NO  -> NO2 + GLYOXAL + 0.5*BIGALD1 + HO2        ; 2.6e-12, 365 
[BENZOOH_OH]           BENZOOH + OH  -> BENZO2                                  ; 3.8e-12, 200 
[BZALD_OH]             BZALD + OH  -> ACBZO2                                    ; 5.9e-12, 225 
[BZOO_HO2]             BZOO + HO2  -> BZOOH                                     ; 7.5e-13, 700 
[BZOOH_OH]             BZOOH + OH  -> BZOO                                      ; 3.8e-12, 200 
[BZOO_NO]              BZOO + NO  -> BZALD + NO2 + HO2                          ; 2.6e-12, 365 
[C6H5O2_HO2]           C6H5O2 + HO2  -> C6H5OOH                                 ; 7.5e-13, 700 
[C6H5O2_NO]            C6H5O2 + NO  -> PHENO + NO2                              ; 2.6e-12, 365 
[C6H5OOH_OH]           C6H5OOH + OH  -> C6H5O2                                  ; 3.8e-12, 200 
[CRESOL_OH]            CRESOL + OH  -> 0.2*PHENO2 + 0.73*HO2 + 0.07*PHENO       ; 4.7e-11 
[DICARBO2_HO2]         DICARBO2 + HO2  -> 0.4*OH + 0.07*HO2 + 0.07*CH3COCHO + 0.07*CO + 0.33*CH3O2 ; 4.3e-13, 1040 
[DICARBO2_NO]          DICARBO2 + NO  -> NO2 + 0.17*HO2 + 0.17*CH3COCHO + 0.17*CO + 0.83*CH3O2 ; 7.5e-12, 290 
[DICARBO2_NO2]         DICARBO2 + NO2 + M  -> M + 1*NDEP                        ; 9.7e-29, 5.6, 9.3e-12, 1.5, 0.6 
[MALO2_HO2]            MALO2 + HO2  -> 0.16*GLYOXAL + 0.16*HO2 + 0.16*CO        ; 4.3e-13, 1040 
[MALO2_NO]             MALO2 + NO  -> 0.4*GLYOXAL + 0.4*HO2 + 0.4*CO + NO2      ; 7.5e-12, 290 
[MALO2_NO2]            MALO2 + NO2 + M  -> M + 1*NDEP                           ; 9.7e-29, 5.6, 9.3e-12, 1.5, 0.6 
[MDIALO2_HO2]          MDIALO2 + HO2  -> 0.4*OH + 0.33*HO2 + 0.07*CH3COCHO + 0.14*CO + 0.07*CH3O2 + 0.07*GLYOXAL ; 4.3e-13, 1040 
[MDIALO2_NO]           MDIALO2 + NO  -> NO2 + 0.83*HO2 + 0.17*CH3COCHO + 0.35*CO + 0.17*CH3O2 + 0.17*GLYOXAL ; 7.5e-12, 290 
[MDIALO2_NO2]          MDIALO2 + NO2 + M  -> M + 1*NDEP                         ; 9.7e-29, 5.6, 9.3e-12, 1.5, 0.6 
[PHENO2_HO2]           PHENO2 + HO2  -> PHENOOH                                 ; 7.5e-13, 700 
[PHENO2_NO]            PHENO2 + NO  -> HO2 + 0.7*GLYOXAL + NO2                  ; 2.6e-12, 365 
[PHENOL_OH]            PHENOL + OH  -> 0.14*PHENO2 + 0.8*HO2 + 0.06*PHENO       ; 4.7e-13, 1220 
[PHENO_NO2]            PHENO + NO2  -> 1*NDEP                                   ; 2.1e-12 
[PHENO_O3]             PHENO + O3  -> C6H5O2                                    ; 2.8e-13 
[PHENOOH_OH]           PHENOOH + OH  -> PHENO2                                  ; 3.8e-12, 200 
[tag_ACBZO2_NO2]       ACBZO2 + NO2 + M  -> PBZNIT + M                          ; 9.7e-29, 5.6, 9.3e-12, 1.5, 0.6 
[TOLO2_HO2]            TOLO2 + HO2  -> TOLOOH                                   ; 7.5e-13, 700 
[TOLO2_NO]             TOLO2 + NO  -> NO2 + 0.6*GLYOXAL + 0.4*CH3COCHO + HO2 + 0.2*BIGALD1 + 0.2*BIGALD2 + 0.2*BIGALD3 ; 2.6e-12, 365 
[TOLOOH_OH]            TOLOOH + OH  -> TOLO2                                    ; 3.8e-12, 200 
[TOLUENE_OH]           TOLUENE + OH  -> 0.18*CRESOL + 0.1*TEPOMUC + 0.07*BZOO + 0.65*TOLO2 + 0.28*HO2 ; 1.7e-12, 352 
[usr_PBZNIT_M]         PBZNIT + M  -> ACBZO2 + NO2 + M                           
[XYLENES_OH]           XYLENES + OH  -> 0.15*XYLOL + 0.23*TEPOMUC + 0.06*BZOO + 0.56*XYLENO2 + 0.38*HO2 ; 1.7e-11 
[XYLENO2_HO2]          XYLENO2 + HO2  -> XYLENOOH                               ; 7.5e-13, 700 
[XYLENO2_NO]           XYLENO2 + NO  -> NO2 + HO2 + 0.34*GLYOXAL + 0.54*CH3COCHO + 0.06*BIGALD1 + 0.2*BIGALD2 + 0.15*BIGALD3 + 0.21*BIGALD4 ; 2.6e-12, 365 
[XYLENOOH_OH]          XYLENOOH + OH  -> XYLENO2                                ; 3.8e-12, 200 
[XYLOLO2_HO2]          XYLOLO2 + HO2  -> XYLOLOOH                               ; 7.5e-13, 700 
[XYLOLO2_NO]           XYLOLO2 + NO  -> HO2 + NO2 + 0.17*GLYOXAL + 0.51*CH3COCHO ; 2.6e-12, 365 
[XYLOL_OH]             XYLOL + OH  -> 0.3*XYLOLO2 + 0.63*HO2 + 0.07*PHENO       ; 8.4e-11 
[XYLOLOOH_OH]          XYLOLOOH + OH  -> XYLOLO2                                ; 3.8e-12, 200 
*********************************
*** C10
*********************************
[BCARY_NO3]            BCARY + NO3  -> NTERPO2                                  ; 1.9e-11 
[BCARY_O3]             BCARY + O3  -> 0.33*TERPROD1 + 0.3*TERPROD2 + 0.63*OH + 0.57*HO2 + 0.23*CO + 0.27*CO2 + 0.52*CH3COCH3 + 0.34*CH2O + 0.1*BIGALD + 0.05*HCOOH + 0.05*BIGALK + 0.06*CH3CO3 + 0.06*RO2 ; 1.2e-14 
[BCARY_OH]             BCARY + OH  -> TERPO2                                    ; 2e-10 
[MTERP_NO3]            MTERP + NO3  -> NTERPO2                                  ; 1.2e-12, 490 
[MTERP_O3]             MTERP + O3  -> 0.33*TERPROD1 + 0.3*TERPROD2 + 0.63*OH + 0.57*HO2 + 0.23*CO + 0.27*CO2 + 0.52*CH3COCH3 + 0.34*CH2O + 0.1*BIGALD + 0.05*HCOOH + 0.05*BIGALK + 0.06*CH3CO3 + 0.06*RO2 ; 6.3e-16, -580 
[MTERP_OH]             MTERP + OH  -> TERPO2                                    ; 1.2e-11, 440 
[NTERPO2_CH3O2]        NTERPO2 + CH3O2  -> 0.5*TERPNIT + 0.75*CH2O + 0.25*CH3OH + 0.5*HO2 + 0.5*TERPROD1 + 0.5*NO2 ; 2e-12, 500 
[NTERPO2_HO2]          NTERPO2 + HO2  -> NTERPOOH                               ; 7.5e-13, 700 
[NTERPO2_NO]           NTERPO2 + NO  -> 0.2*TERPNIT + 1.6*NO2 + 0.8*TERPROD1 + 0.2*NDEP ; 4.2e-12, 180 
[NTERPO2_NO3]          NTERPO2 + NO3  -> 2*NO2 + TERPROD1                       ; 2.4e-12 
[NTERPOOH_OH]          NTERPOOH + OH  -> NTERPO2                                ; 2e-11 
[TERP2O2_CH3O2]        TERP2O2 + CH3O2  -> TERPROD2 + 0.93*CH2O + 0.25*CH3OH + HO2 + 0.5*CO2 + 0.125*CO + 0.125*GLYALD + 0.15*CH3COCH3 ; 2e-12, 500 
[TERP2O2_HO2]          TERP2O2 + HO2  -> TERP2OOH                               ; 7.5e-13, 700 
[TERP2O2_NO]           TERP2O2 + NO  -> 0.1*ONITR + 0.9*NO2 + 0.34*CH2O + 0.27*CH3COCH3 + 0.225*CO + 0.9*CO2 + 0.9*TERPROD2 + 0.9*HO2 + 0.225*GLYALD ; 4.2e-12, 180 
[TERP2OOH_OH]          TERP2OOH + OH  -> TERP2O2                                ; 2.3e-11 
[TERPNIT_OH]           TERPNIT + OH  -> NO2 + TERPROD1                          ; 2e-11 
[TERPO2_CH3O2]         TERPO2 + CH3O2  -> TERPROD1 + 0.95*CH2O + 0.25*CH3OH + HO2 + 0.025*CH3COCH3 ; 2e-12, 500 
[TERPO2_HO2]           TERPO2 + HO2  -> TERPOOH                                 ; 7.5e-13, 700 
[TERPO2_NO]            TERPO2 + NO  -> 0.2*TERPNIT + 0.8*NO2 + 0.32*CH2O + 0.04*CH3COCH3 + 0.8*TERPROD1 + 0.8*HO2 ; 4.2e-12, 180 
[TERPOOH_OH]           TERPOOH + OH  -> TERPO2                                  ; 3.3e-11 
[TERPROD1_NO3]         TERPROD1 + NO3  -> 0.5*TERP2O2 + 0.5*NTERPO2 + 0.5*NDEP  ; 1e-12 
[TERPROD1_OH]          TERPROD1 + OH  -> TERP2O2                                ; 5.7e-11 
[TERPROD2_OH]          TERPROD2 + OH  -> 0.15*RO2 + 0.68*CH2O + 1.8*CO2 + 0.5*CH3COCH3 + 0.65*CH3CO3 + 0.2*HO2 + 0.7*CO ; 3.4e-11 
*********************************
*** Sulfur
*********************************
[DMS_NO3]              DMS + NO3  -> SO2 + HNO3                                 ; 1.9e-13, 520 
[DMS_OHa]              DMS + OH  -> SO2                                         ; 1.1e-11, -280 
[OCS_O]                OCS + O  -> SO + CO                                      ; 2.1e-11, -2200 
[OCS_OH]               OCS + OH  -> SO2 + CO + H                                ; 7.2e-14, -1070 
[S_O2]                 S + O2  -> SO + O                                        ; 2.3e-12 
[SO2_OH_M]             SO2 + OH + M  -> SO3 + HO2                               ; 2.9e-31, 4.1, 1.7e-12, -0.2, 0.6 
[S_O3]                 S + O3  -> SO + O2                                       ; 1.2e-11 
[SO_BRO]               SO + BRO  -> SO2 + BR                                    ; 5.7e-11 
[SO_CLO]               SO + CLO  -> SO2 + CL                                    ; 2.8e-11 
[S_OH]                 S + OH  -> SO + H                                        ; 6.6e-11 
[SO_NO2]               SO + NO2  -> SO2 + NO                                    ; 1.4e-11 
[SO_O2]                SO + O2  -> SO2 + O                                      ; 1.6e-13, -2280 
[SO_O3]                SO + O3  -> SO2 + O2                                     ; 3.4e-12, -1100 
[SO_OCLO]              SO + OCLO  -> SO2 + CLO                                  ; 1.9e-12 
[SO_OH]                SO + OH  -> SO2 + H                                      ; 2.6e-11, 330 
[usr_DMS_OH]           DMS + OH  -> 0.5*SO2 + 0.5*HO2                            
[usr_SO3_H2O]          SO3 + H2O  -> H2SO4                                       
*********************************
*** Tropospheric Aerosol
*********************************
[NH3_OH]               NH3 + OH  -> H2O + 1*NHDEP                               ; 1.7e-12, -710 
[usr_HO2_aer]          HO2  -> H2O                                               
[usr_HONITR_aer]       HONITR  -> HNO3                                           
[usr_ISOPNITA_aer]     ISOPNITA  -> HNO3                                         
[usr_ISOPNITB_aer]     ISOPNITB  -> HNO3                                         
[usr_N2O5_aer]         N2O5  -> 2*HNO3                                           
[usr_NC4CH2OH_aer]     NC4CH2OH  -> HNO3                                         
[usr_NC4CHO_aer]       NC4CHO  -> HNO3                                           
[usr_NH4_strat_tau]    NH4  -> 1*NHDEP                                          ; 6.34e-08 
[usr_NO2_aer]          NO2  -> 0.5*OH + 0.5*NO + 0.5*HNO3                        
[usr_NO3_aer]          NO3  -> HNO3                                              
[usr_NTERPOOH_aer]     NTERPOOH  -> HNO3                                         
[usr_ONITR_aer]        ONITR  -> HNO3                                            
[usr_TERPNIT_aer]      TERPNIT  -> HNO3                                          
*********************************
*** SOA
*********************************
[BCARY_NO3_vbs]        BCARY + NO3  -> BCARY + NO3 + 0.17493*SOAG3 + 0.59019*SOAG4 ; 1.9e-11 
[BCARYO2_HO2_vbs]      BCARYO2VBS + HO2  -> HO2 + 0.2202*SOAG0 + 0.2067*SOAG1 + 0.0653*SOAG2 + 0.1284*SOAG3 + 0.114*SOAG4 ; 2.75e-13, 1300 
[BCARYO2_NO_vbs]       BCARYO2VBS + NO  -> NO + 0.1279*SOAG0 + 0.1792*SOAG1 + 0.0676*SOAG2 + 0.079*SOAG3 + 0.1254*SOAG4 ; 2.7e-12, 360 
[BCARY_O3_vbs]         BCARY + O3  -> BCARY + O3 + 0.2202*SOAG0 + 0.2067*SOAG1 + 0.0653*SOAG2 + 0.1284*SOAG3 + 0.114*SOAG4 ; 1.2e-14 
[BCARY_OH_vbs]         BCARY + OH  -> BCARY + OH + BCARYO2VBS                   ; 2e-10 
[BENZENE_OH_vbs]       BENZENE + OH  -> BENZENE + OH + BENZO2VBS                ; 2.3e-12, -193 
[BENZO2_HO2_vbs]       BENZO2VBS + HO2  -> HO2 + 0.0023*SOAG0 + 0.0008*SOAG1 + 0.0843*SOAG2 + 0.0443*SOAG3 + 0.1621*SOAG4 ; 7.5e-13, 700 
[BENZO2_NO_vbs]        BENZO2VBS + NO  -> NO + 0.0097*SOAG0 + 0.0034*SOAG1 + 0.1579*SOAG2 + 0.0059*SOAG3 + 0.0536*SOAG4 ; 2.6e-12, 365 
[ISOP_NO3_vbs]         ISOP + NO3  -> ISOP + NO3 + 0.059024*SOAG3 + 0.025024*SOAG4 ; 3.03e-12, -446 
[ISOPO2_HO2_vbs]       ISOPO2VBS + HO2  -> HO2 + 0.0031*SOAG0 + 0.0035*SOAG1 + 0.0003*SOAG2 + 0.0271*SOAG3 + 0.0474*SOAG4 ; 2.12e-13, 1300 
[ISOPO2_NO_vbs]        ISOPO2VBS + NO  -> NO + 0.0003*SOAG0 + 0.0003*SOAG1 + 0.0073*SOAG2 + 0.0057*SOAG3 + 0.0623*SOAG4 ; 2.7e-12, 350 
[ISOP_O3_vbs]          ISOP + O3  -> ISOP + O3 + 0.0033*SOAG3                   ; 1.05e-14, -2000 
[ISOP_OH_vbs]          ISOP + OH  -> ISOP + OH + ISOPO2VBS                      ; 2.54e-11, 410 
[IVOCO2_HO2_vbs]       IVOCO2VBS + HO2  -> HO2 + 0.2381*SOAG0 + 0.1308*SOAG1 + 0.0348*SOAG2 + 0.0076*SOAG3 + 0.0113*SOAG4 ; 7.5e-13, 700 
[IVOCO2_NO_vbs]        IVOCO2VBS + NO  -> NO + 0.1056*SOAG0 + 0.1026*SOAG1 + 0.0521*SOAG2 + 0.0143*SOAG3 + 0.0166*SOAG4 ; 2.6e-12, 365 
[IVOC_OH_vbs]          IVOC + OH  -> OH + IVOCO2VBS                             ; 1.34e-11 
[MTERP_NO3_vbs]        MTERP + NO3  -> MTERP + NO3 + 0.17493*SOAG3 + 0.59019*SOAG4 ; 1.2e-12, 490 
[MTERPO2_HO2_vbs]      MTERPO2VBS + HO2  -> HO2 + 0.0508*SOAG0 + 0.1149*SOAG1 + 0.0348*SOAG2 + 0.0554*SOAG3 + 0.1278*SOAG4 ; 2.6e-13, 1300 
[MTERPO2_NO_vbs]       MTERPO2VBS + NO  -> NO + 0.0245*SOAG0 + 0.0082*SOAG1 + 0.0772*SOAG2 + 0.0332*SOAG3 + 0.13*SOAG4 ; 2.7e-12, 360 
[MTERP_O3_vbs]         MTERP + O3  -> MTERP + O3 + 0.0508*SOAG0 + 0.1149*SOAG1 + 0.0348*SOAG2 + 0.0554*SOAG3 + 0.1278*SOAG4 ; 6.3e-16, -580 
[MTERP_OH_vbs]         MTERP + OH  -> MTERP + OH + MTERPO2VBS                   ; 1.2e-11, 440 
[SVOC_OH]              SVOC + OH  -> OH + 0.5931*SOAG0 + 0.1534*SOAG1 + 0.0459*SOAG2 + 0.0085*SOAG3 + 0.0128*SOAG4 ; 1.34e-11 
[TOLUENE_OH_vbs]       TOLUENE + OH  -> TOLUENE + OH + TOLUO2VBS                ; 1.7e-12, 352 
[TOLUO2_HO2_vbs]       TOLUO2VBS + HO2  -> HO2 + 0.1364*SOAG0 + 0.0101*SOAG1 + 0.0763*SOAG2 + 0.2157*SOAG3 + 0.0738*SOAG4 ; 7.5e-13, 700 
[TOLUO2_NO_vbs]        TOLUO2VBS + NO  -> NO + 0.0154*SOAG0 + 0.0452*SOAG1 + 0.0966*SOAG2 + 0.0073*SOAG3 + 0.238*SOAG4 ; 2.6e-12, 365 
[usr_GLYOXAL_aer]      GLYOXAL  -> SOAG0                                         
[XYLENES_OH_vbs]       XYLENES + OH  -> XYLENES + OH + XYLEO2VBS                ; 1.7e-11 
[XYLEO2_HO2_vbs]       XYLEO2VBS + HO2  -> HO2 + 0.1677*SOAG0 + 0.0174*SOAG1 + 0.086*SOAG2 + 0.0512*SOAG3 + 0.1598*SOAG4 ; 7.5e-13, 700 
[XYLEO2_NO_vbs]        XYLEO2VBS + NO  -> NO + 0.0063*SOAG0 + 0.0237*SOAG1 + 0.0025*SOAG2 + 0.011*SOAG3 + 0.1185*SOAG4 ; 2.6e-12, 365 
*********************************
*** Stratospheric Aerosol
*********************************
[het1]                 N2O5  -> 2*HNO3                                           
[het10]                HOCL + HCL  -> CL2 + H2O                                  
[het11]                BRONO2  -> HOBR + HNO3                                    
[het12]                N2O5  -> 2*HNO3                                           
[het13]                CLONO2  -> HOCL + HNO3                                    
[het14]                BRONO2  -> HOBR + HNO3                                    
[het15]                CLONO2 + HCL  -> CL2 + HNO3                               
[het16]                HOCL + HCL  -> CL2 + H2O                                  
[het17]                HOBR + HCL  -> BRCL + H2O                                 
[het2]                 CLONO2  -> HOCL + HNO3                                    
[het3]                 BRONO2  -> HOBR + HNO3                                    
[het4]                 CLONO2 + HCL  -> CL2 + HNO3                               
[het5]                 HOCL + HCL  -> CL2 + H2O                                  
[het6]                 HOBR + HCL  -> BRCL + H2O                                 
[het7]                 N2O5  -> 2*HNO3                                           
[het8]                 CLONO2  -> HOCL + HNO3                                    
[het9]                 CLONO2 + HCL  -> CL2 + HNO3                               
*********************************
*** Tracers
*********************************
[NH_50_tau]            NH_50  ->                                                ; 2.31e-07 
[NH_5_tau]             NH_5  ->                                                 ; 2.31e-06 
[ST80_25_tau]          ST80_25  ->                                              ; 4.63e-07 
      End Reactions

      Ext Forcing
 num_a4 <- dataset 
 pom_a4 <- dataset 
 bc_a4 <- dataset 
 SVOC <- dataset 
 so4_a1 <- dataset 
 so4_a2 <- dataset 
 CO <- dataset 
 SO2 <- dataset 
 NO2 <- dataset 
 num_a1 <- dataset 
 num_a2 <- dataset 
 AOA_NH 
 NO 
 N 
      End Ext Forcing

      End Chemistry

      SIMULATION PARAMETERS

      Version Options
        machine = nec
        model   = cam
        model_architecture = VECTOR
        vector_length = 32
        architecture = hybrid
        namemod = on
      End Version Options


      End Simulation Parameters

ENDSIM
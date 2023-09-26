MODULE advec_interface
   IMPLICIT NONE
   CONTAINS
   SUBROUTINE advec_dwarf_interface(legacyoptmode,t_adv,   &
                           lupdatemulti, lvertsplit, ipoles)
   USE advec_driver, ONLY: advec_dwarf
   USE scratch_datafields, ONLY: xtracer,bcx,bcy
   USE scratch_datafields, ONLY: rhr,rhoadv
   USE scratch_datafields, ONLY: uadv,vadv,wadv
   USE mod_parameters, ONLY: ih,ibcx,ibcy,ibcz,ibcx,ibcy,ibcz
   USE mod_parameters, ONLY: np,mp,lp
   INTEGER,INTENT(IN) :: t_adv,legacyoptmode,ipoles
   LOGICAL,INTENT(IN) :: lupdatemulti,lvertsplit
   CALL advec_dwarf(legacyoptmode,t_adv,   &
                    lupdatemulti, lvertsplit,         &
                    np,mp,lp,ih,ipoles,ibcx,ibcy,ibcz,           &
                    xtracer,bcx,bcy,                             &
                    uadv,vadv,wadv,rhr,rhoadv) 
   END SUBROUTINE advec_dwarf_interface
END MODULE advec_interface
 

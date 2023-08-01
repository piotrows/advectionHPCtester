MODULE advec_interface_dp
   IMPLICIT NONE
   CONTAINS
   SUBROUTINE allocate_interface_dp(linitmpi, nprocx, nprocy, nprocz)
   USE advec_initialize, ONLY : allocate_and_initialize 
   LOGICAL, INTENT(IN) :: linitmpi
   INTEGER, INTENT(IN) :: nprocx, nprocy, nprocz
   CALL allocate_and_initialize(linitmpi, nprocx, nprocy, nprocz) 
   END SUBROUTINE allocate_interface_dp
   
   SUBROUTINE deallocate_interface_dp(itime_counter)
   USE advec_initialize, ONLY : deallocate_and_finalize 
   INTEGER, INTENT(IN) :: itime_counter 
   CALL deallocate_and_finalize(itime_counter) 
   END SUBROUTINE deallocate_interface_dp

   SUBROUTINE advec_dwarf_interface_dp(legacyoptmode,t_adv,   &
                           lupdatemulti, lvertdplit, ipoles,itime_counter)
   USE advec_driver, ONLY: advec_dwarf
   USE scratch_datafields, ONLY: xtracer,bcx,bcy
   USE scratch_datafields, ONLY: rhr,rhoadv
   USE scratch_datafields, ONLY: uadv,vadv,wadv
   USE parameters, ONLY: ih,ibcx,ibcy,ibcz,ibcx,ibcy,ibcz
   USE parameters, ONLY: np,mp,lp
   USE precisions
   INTEGER,INTENT(IN) :: t_adv,legacyoptmode,ipoles,itime_counter
   LOGICAL,INTENT(IN) :: lupdatemulti,lvertdplit
   IF(itime_counter.eq.1) print *,'Advec precision is',euwp,'bytes'
   CALL advec_dwarf(legacyoptmode,t_adv,   &
                    lupdatemulti, lvertdplit,         &
                    np,mp,lp,ih,ipoles,ibcx,ibcy,ibcz,           &
                    xtracer,bcx,bcy,                             &
                    uadv,vadv,wadv,rhr,rhoadv) 
   END SUBROUTINE advec_dwarf_interface_dp
END MODULE advec_interface_dp
 

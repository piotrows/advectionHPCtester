MODULE advec_interface_dpgpu
   USE, INTRINSIC :: iso_c_binding
   IMPLICIT NONE
   CONTAINS
   SUBROUTINE allocate_interface_dpgpu(linitmpi, nprocx, nprocy, nprocz)  bind(c, name='allocate_interface_dpgpu')
   USE advec_initialize, ONLY : allocate_and_initialize 
   LOGICAL, INTENT(IN) :: linitmpi
   INTEGER, INTENT(IN) :: nprocx, nprocy, nprocz
   CALL allocate_and_initialize(linitmpi, nprocx, nprocy, nprocz) 
   END SUBROUTINE allocate_interface_dpgpu
   
   SUBROUTINE deallocate_interface_dpgpu(itime_counter) bind(c, name='deallocate_interface_dpgpu')
   USE advec_initialize, ONLY : deallocate_and_finalize 
   INTEGER(C_INT), INTENT(IN) :: itime_counter 
   print *,'inside deallocate'
   CALL deallocate_and_finalize(INT(itime_counter))
   END SUBROUTINE deallocate_interface_dpgpu

   SUBROUTINE advec_dwarf_interface_dpgpu(legacyoptmode,t_adv,   &
                           lupdatemulti, lvertdplit, ipoles,itime_counter) bind(c, name='advec_dwarf_interface_dpgpu')
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
   END SUBROUTINE advec_dwarf_interface_dpgpu
END MODULE advec_interface_dpgpu
 

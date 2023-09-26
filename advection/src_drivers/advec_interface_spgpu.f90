#include  "../src_algorithms/renames.inc"
MODULE advec_interface_spgpu
   USE, INTRINSIC :: iso_c_binding
   IMPLICIT NONE
   CONTAINS
   SUBROUTINE allocate_interface_spgpu(linitmpi, nprocx, nprocy, nprocz)  bind(c, name='allocate_interface_spgpu')
   USE advec_initialize, ONLY : allocate_and_initialize 
   LOGICAL, INTENT(IN) :: linitmpi
   INTEGER, INTENT(IN) :: nprocx, nprocy, nprocz
   CALL allocate_and_initialize(linitmpi, nprocx, nprocy, nprocz) 
   END SUBROUTINE allocate_interface_spgpu
   
   SUBROUTINE deallocate_interface_spgpu(itime_counter) bind(c, name='deallocate_interface_spgpu')
   USE advec_initialize, ONLY : deallocate_and_finalize 
   INTEGER, INTENT(IN) :: itime_counter 
   CALL deallocate_and_finalize(itime_counter) 
   END SUBROUTINE deallocate_interface_spgpu

   SUBROUTINE advec_dwarf_interface_spgpu(legacyoptmode,t_adv,   &
                           lupdatemulti, lvertsplit, ipoles,itime_counter) bind(c, name='advec_dwarf_interface_spgpu')
   USE advec_driver, ONLY: advec_dwarf
   USE scratch_datafields, ONLY: xtracer,bcx,bcy
   USE scratch_datafields, ONLY: rhr,rhoadv
   USE scratch_datafields, ONLY: uadv,vadv,wadv
   USE mod_parameters, ONLY: ih,ibcx,ibcy,ibcz,ibcx,ibcy,ibcz
   USE mod_parameters, ONLY: np,mp,lp
   USE precisions
   INTEGER,INTENT(IN) :: t_adv,legacyoptmode,ipoles,itime_counter
   LOGICAL,INTENT(IN) :: lupdatemulti,lvertsplit
   IF(itime_counter.eq.1) print *,'Advec precision is',euwp,'bytes'
   CALL advec_dwarf(legacyoptmode,t_adv,   &
                    lupdatemulti, lvertsplit,         &
                    np,mp,lp,ih,ipoles,ibcx,ibcy,ibcz,           &
                    xtracer,bcx,bcy,                             &
                    uadv,vadv,wadv,rhr,rhoadv) 
   END SUBROUTINE advec_dwarf_interface_spgpu
END MODULE advec_interface_spgpu
 

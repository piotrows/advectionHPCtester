MODULE advec_init_interface_spgpu
   USE advec_initialize
   IMPLICIT NONE
CONTAINS
   SUBROUTINE allocate_interface_spgpu(linitmpi, nprocx, nprocy, nprocz)
   USE precisions
   LOGICAL, INTENT(IN) :: linitmpi
   INTEGER, INTENT(IN) :: nprocx, nprocy, nprocz 
   CALL allocate_and_initialize(linitmpi, nprocx, nprocy, nprocz)
   END SUBROUTINE allocate_interface_spgpu
   

   SUBROUTINE deallocate_interface_spgpu(itime_counter)
   INTEGER, INTENT(IN) :: itime_counter
   CALL deallocate_and_finalize(itime_counter)
   END SUBROUTINE deallocate_interface_spgpu
END MODULE advec_init_interface_spgpu

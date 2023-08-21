MODULE advec_init_interface_dpgpu
   USE advec_initialize
   IMPLICIT NONE
CONTAINS
   SUBROUTINE allocate_interface_dpgpu(linitmpi, nprocx, nprocy, nprocz)
   USE precisions
   LOGICAL, INTENT(IN) :: linitmpi
   INTEGER, INTENT(IN) :: nprocx, nprocy, nprocz 
   CALL allocate_and_initialize(linitmpi, nprocx, nprocy, nprocz)
   END SUBROUTINE allocate_interface_dpgpu
   

   SUBROUTINE deallocate_interface_dpgpu(itime_counter)
   INTEGER, INTENT(IN) :: itime_counter
   CALL deallocate_and_finalize(itime_counter)
   END SUBROUTINE deallocate_interface_dpgpu
END MODULE advec_init_interface_dpgpu

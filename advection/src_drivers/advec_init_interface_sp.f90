MODULE advec_init_interface_sp
   IMPLICIT NONE

CONTAINS
   SUBROUTINE allocate_interface_sp(linitmpi, nprocx, nprocy, nprocz)
   USE precisions
   LOGICAL, INTENT(IN) :: linitmpi
   INTEGER, INTENT(IN) :: nprocx, nprocy, nprocz 
   CALL allocate_and_initialize(linitmpi, nprocx, nprocy, nprocz)
   END SUBROUTINE allocate_interface_sp
   

   SUBROUTINE deallocate_interface_sp(itime_counter)
   INTEGER, INTENT(IN) :: itime_counter
   CALL deallocate_and_finalize(itime_counter)
   END SUBROUTINE deallocate_interface_sp
END MODULE advec_init_interface_sp

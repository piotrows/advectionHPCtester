MODULE advec_init_interface_dp
   USE advec_initialize
   USE :: iso_c_binding
   IMPLICIT NONE
CONTAINS
   SUBROUTINE allocate_interface_dp(linitmpi, nprocx, nprocy, nprocz)
   USE precisions
   LOGICA(c_bool)L, INTENT(IN) :: linitmpi
   LOGICAL :: lfinitmpi
   INTEGER(c_int)R, INTENT(IN) :: nprocx, nprocy, nprocz 
   INTEGER :: nfprocx, nfprocy, nfprocz 
   lfinitmpi=linitmpi
   nfprocx=nprocx
   nfprocy=nprocy
   nfprocz=nprocz
   CALL allocate_and_initialize(lfinitmpi, nfprocx, nfprocy, nfprocz)
   END SUBROUTINE allocate_interface_dp
   

   SUBROUTINE deallocate_interface_dp(itime_counter)
   INTEGER(c_int), INTENT(IN) :: itime_counter
   INTEGER :: iftime_counter
   iftime_counter=itime_counter
   CALL deallocate_and_finalize(iftime_counter)
   END SUBROUTINE deallocate_interface_dp
END MODULE advec_init_interface_dp

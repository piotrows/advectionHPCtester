#include  "../src_algorithms/renames.inc"
MODULE advec_init_interface_dpgpu
   USE advec_initialize
   USE :: iso_c_binding
   IMPLICIT NONE
CONTAINS
   SUBROUTINE allocate_interface_dpgpu(linitmpi, nprocx, nprocy, nprocz)
   USE precisions
   LOGICAL(c_bool), INTENT(IN) :: linitmpi
   LOGICAL :: lfinitmpi
   INTEGER(c_int), INTENT(IN) :: nprocx, nprocy, nprocz 
   INTEGER :: nfprocx, nfprocy, nfprocz 
   lfinitmpi=linitmpi
   nfprocx=nprocx
   nfprocy=nprocy
   nfprocz=nprocz
   CALL allocate_and_initialize(lfinitmpi, nfprocx, nfprocy, nfprocz)
   END SUBROUTINE allocate_interface_dpgpu
   

   SUBROUTINE deallocate_interface_dpgpu(itime_counter)
   INTEGER(c_int), INTENT(IN) :: itime_counter
   INTEGER :: iftime_counter
   iftime_counter=itime_counter
   CALL deallocate_and_finalize(iftime_counter)
   END SUBROUTINE deallocate_interface_dpgpu
END MODULE advec_init_interface_dpgpu

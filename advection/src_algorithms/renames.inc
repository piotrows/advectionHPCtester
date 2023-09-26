!#define PASTER(x,y) x ## _ ## y
!#define EVALUATOR(x,y)  PASTER(x,y)
#define PASTE2(x,y) x/**/y
#define EVALUATOR(x,y) PASTE2(PASTE2(x,_),y)
#define SUFFIXIZE(fun) EVALUATOR(fun,SYMBOLSUFFIX)
#define donor SUFFIXIZE(donor)
#define upwind3d_gpubc SUFFIXIZE(upwind3d_gpubc)
#define antidiff3d_gauge_gpubc SUFFIXIZE(antidiff3d_gauge_gpubc)
#define antidiff3d_standard_gpubc SUFFIXIZE(antidiff3d_standard_gpubc)
#define vcorra SUFFIXIZE(vcorra)
#define vdyfa SUFFIXIZE(vdyfa)
#define vcorr SUFFIXIZE(vcorr)
#define vdyf SUFFIXIZE(vdyf)
#define rat4 SUFFIXIZE(rat4)
#define rat2 SUFFIXIZE(rat2)
#define pn SUFFIXIZE(pn)
#define pp SUFFIXIZE(pp)
#define putbcinxh SUFFIXIZE(putbcinxh)
#define putbcinx SUFFIXIZE(putbcinx)
#define zerogradbcinx SUFFIXIZE(zerogradbcinx)
#define zerogradbcinxh SUFFIXIZE(zerogradbcinxh)
#define flip_halo_poles SUFFIXIZE(flip_halo_poles) 
#define enforce_cyclic SUFFIXIZE(enforce_cyclic) 
#define enforce_cyclic_z SUFFIXIZE(enforce_cyclic_z) 
#define enforce_cyclic_xy SUFFIXIZE(enforce_cyclic_xy) 
#define cp_last_to_halo_z SUFFIXIZE(cp_last_to_halo_z) 
#define cp_scnd_last_to_halo_z SUFFIXIZE(cp_scnd_last_to_halo_z) 
#define cp_scnd_last_to_halo_z_full SUFFIXIZE(cp_scnd_last_to_halo_z_full) 
#define cp_scnd_last_to_halo_xy SUFFIXIZE(cp_scnd_last_to_halo_xy) 
#define cp_scnd_last_to_halo_xyz SUFFIXIZE(cp_scnd_last_to_halo_xyz) 
#define cp_scnd_last_to_halo_xyz_full SUFFIXIZE(cp_scnd_last_to_halo_xyz_full) 
#define cp_last_to_halo_xy SUFFIXIZE(cp_last_to_halo_xy)
#define cp_last_to_halo_xyz_full SUFFIXIZE(cp_last_to_halo_xyz_full)
#define cp_last_to_halo_xyd SUFFIXIZE(cp_last_to_halo_xyd)
#define flux_x_bc SUFFIXIZE(flux_x_bc)
#define flux_y_bc SUFFIXIZE(flux_y_bc)
#define flux_z_bc SUFFIXIZE(flux_z_bc)
#define flux_x_bcd SUFFIXIZE(flux_x_bcd)
#define flux_y_bcd SUFFIXIZE(flux_y_bcd)
#define flux_z_bc_nohalo SUFFIXIZE(flux_z_bc_nohalo)
#define flux_z_bc_nohalod SUFFIXIZE(flux_z_bc_nohalod)
#define zero_halo_xy SUFFIXIZE(zero_halo_xy)
#define zero_halo_xyd SUFFIXIZE(zero_halo_xyd)
#define cp_scnd_last_to_halo_xy_fulld  SUFFIXIZE(cp_scnd_last_to_halo_xy_fulld)
#define cp_scnd_last_to_halo_xy_full  SUFFIXIZE(cp_scnd_last_to_halo_xy_full)
#define velprd_driver SUFFIXIZE(velprd_driver)
#define velprdA_to_C SUFFIXIZE(velprdA_to_C)
#define velprd_traj0 SUFFIXIZE(velprd_traj0)
#define velprd_traj0_adapt_dt SUFFIXIZE(velprd_traj0_adapt_dt)
#define advec_dwarf SUFFIXIZE(advec_dwarf)
#define upwind_firstep_gpubc_qg_multiker SUFFIXIZE(upwind_firstep_gpubc_qg_multiker)
#define upwind_firstep_gpubc_dry SUFFIXIZE(upwind_firstep_gpubc_dry)
#define upwind_firstep_gpubc_qv SUFFIXIZE(upwind_firstep_gpubc_qv)
#define upwind_firstep_gpubc_qc SUFFIXIZE(upwind_firstep_gpubc_qc)
#define upwind_firstep_gpubc_qr SUFFIXIZE(upwind_firstep_gpubc_qr)
#define upwind_firstep_gpubc_qi SUFFIXIZE(upwind_firstep_gpubc_qi)
#define upwind_firstep_gpubc_qs SUFFIXIZE(upwind_firstep_gpubc_qs)
#define upwind_firstep_gpubc_qg SUFFIXIZE(upwind_firstep_gpubc_qg)
#define initialize_for_gpu SUFFIXIZE(initialize_for_gpu)
#define fill_full_variable_set SUFFIXIZE(fill_full_variable_set)
#define fill_test_variable_set SUFFIXIZE(fill_test_variable_set)
#define update3dsplit  SUFFIXIZE(update3dsplit)
#define update3dsplit_gpud  SUFFIXIZE(update3dsplit_gpud)
#define compute_courlipsch_Agrid_full  SUFFIXIZE(compute_courlipsch_Agrid_full)
#define compute_courants  SUFFIXIZE(compute_courants)
#define advec_driver SUFFIXIZE(advec_driver)
#define module_upwind3d_gpubc SUFFIXIZE(module_upwind3d_gpubc)
#define module_antidiff3d_gauge_gpubc SUFFIXIZE(module_antidiff3d_gauge_gpubc)
#define module_antidiff3d_standard_gpubc SUFFIXIZE(module_antidiff3d_standard_gpubc)
#define mpdataoperators SUFFIXIZE(mpdataoperators)
#define mpi_parallel SUFFIXIZE(mpi_parallel)
#define nvtx SUFFIXIZE(nvtx)
#define advec_interface_dp SUFFIXIZE(advec_interface_dp)
#define advec_interface_sp SUFFIXIZE(advec_interface_sp)
#define advec_interface_dpgpu SUFFIXIZE(advec_interface_dpgpu)
#define advec_interface_spgpu SUFFIXIZE(advec_interface_spgpu)
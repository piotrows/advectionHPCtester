MODULE advec_initialize
   USE iso_c_binding
#ifdef CUDACODE
   USE cudafor
#endif
   USE mpi_parallel, ONLY: mype,geomset,end_code,end_mpi
   USE mpi_parallel, ONLY: ttini,ttprt,ttbeg,ttend
   USE mpi_parallel, ONLY: init_subdomains
   USE advec_driver, ONLY: advec_dwarf
   USE parameters, ONLY: ih,ibcx,ibcy,ibcz
   USE parameters, ONLY: n,m,l,np,mp,lp,nt
   USE parameters, ONLY: set_eulagcommon_mpi_parameters
   USE parameters, ONLY: set_eulagcommon_lib_domainsize
   USE geometry, ONLY: topolog,metryc
   USE scratch_datafields, ONLY: xtracer,bcx,bcy
   USE scratch_datafields, ONLY: rhr,rhoadv,rhr2,rhoadv2,rhr3,rhoadv3
   USE scratch_datafields, ONLY: uadv,vadv,wadv,wadv2
   USE testing, ONLY: set_initial_tracer,init_data, init_data_sphere
   USE testing, ONLY: initialize_for_gpu
#if STATICMEM == 0
   USE eulag_datafields,   ONLY:   allocate_eulag_datafields
   USE eulag_datafields,   ONLY: deallocate_eulag_datafields
   USE geometry_datafields,ONLY:   allocate_geometry_basic_datafields
   USE geometry_datafields,ONLY: deallocate_geometry_basic_datafields
   USE scratch_datafields, ONLY:   allocate_scratch_mpdata_datafields
   USE scratch_datafields, ONLY:   allocate_scratch_advdrv_datafields
   USE scratch_datafields, ONLY:   allocate_scratch_advtst_datafields
   USE scratch_datafields, ONLY:   allocate_scratch_geometry_datafields
   USE scratch_datafields, ONLY:   allocate_scratch_mpi_datafields
   USE scratch_datafields, ONLY:   allocate_scratch_aux_datafields
   USE scratch_datafields, ONLY: deallocate_scratch_mpdata_datafields
   USE scratch_datafields, ONLY: deallocate_scratch_advdrv_datafields
   USE scratch_datafields, ONLY: deallocate_scratch_advtst_datafields
   USE scratch_datafields, ONLY: deallocate_scratch_geometry_datafields
   USE scratch_datafields, ONLY: deallocate_scratch_mpi_datafields
   USE scratch_datafields, ONLY: deallocate_scratch_aux_datafields
#endif /*STATICMEM*/
   IMPLICIT NONE

CONTAINS
   SUBROUTINE allocate_and_initialize(linitmpi, nprocx, nprocy, nprocz)
   USE precisions
   LOGICAL, INTENT(IN) :: linitmpi
   INTEGER, INTENT(IN) :: nprocx, nprocy, nprocz 
   REAL(KIND=euwp) :: dt_loc, dt_orig
   INTEGER(KIND=iintegers) :: iscale, itmax, itmin, itstr, iprint
   INTEGER(KIND=iintegers) :: nx, ny, nz, nstarttest 
   INTEGER(KIND=iintegers), PARAMETER :: ipoles=1
   LOGICAL :: lupdatemulti, lvertsplit

#if (STATICMEM == 0)
!call get_args_from_commandline
#endif

#if (STATICMEM == 0)
!   CALL set_default_dynamicmem_parameters
#endif
!  CALL read_commandline_parameters
   print *,'Precision is',euwp,'bytes'
   itmin=1
   itmax=16
   itstr=1
!  IF(mype.eq.0) print *,itestset,itmin,itmax,itstr
   iscale=1
   nx=59*iscale
   dt_orig=(0.018_euwp*2._euwp*acos(-1._euwp))
   dt_loc=dt_orig/REAL(iscale)
   nt=((556+0))*iscale
   nt=((556+0))!*iscale
!  nt=500
!  nprocx=1
!  nprocy=1
!  nprocz=1
!  CALL set_advection_lib_parameters(59,59,59,59,59,59,1,1,1,0,0,0,2,3,1,dt_loc,561)
   CALL set_eulagcommon_mpi_parameters(n_in=nx,nprocx_in=nprocx,nprocy_in=nprocy,nprocz_in=nprocz,   &
                                                 ibcx_in=0,  ibcy_in=0,  ibcz_in=0,   &
                                              isphere_in=0,ipoles_in=ipoles)
   CALL set_eulagcommon_lib_domainsize(   n_in=nx,    m_in=nx,   l_in=nx, ih_in=1,    &
                                       dx00_in=100._euwp, dy00_in=100._euwp, dz00_in=100._euwp, &
                                         nt_in=nt, dt_in=dt_orig )
   CALL geomset (linitmpi)
   CALL init_subdomains
#if (STATICMEM == 0)
!   CALL init_dynamicmem_parameters
!   CALL allocate_eulag_datafields
   CALL allocate_geometry_basic_datafields
!  CALL allocate_scratch_datafields
   CALL allocate_scratch_mpdata_datafields
   CALL allocate_scratch_advdrv_datafields
   CALL allocate_scratch_advtst_datafields
!   CALL allocate_scratch_geometry_datafields
#ifdef PUREMPI
  CALL allocate_scratch_mpi_datafields
#endif
   CALL allocate_scratch_aux_datafields

#endif /*STATICMEM*/
   CALL ttini
   CALL init_data(rhr,rhoadv,bcx,bcy,np,mp,lp,ih)
   CALL set_initial_tracer(uadv,vadv,wadv,xtracer,rhoadv,np,mp,lp,ih)
   CALL set_eulagcommon_lib_domainsize(   n_in=nx,    m_in=nx,   l_in=nx, ih_in=1,    &
                                       dx00_in=100._euwp, dy00_in=100._euwp, dz00_in=100._euwp, &
                                         nt_in=nt, dt_in=dt_loc )
   CALL initialize_for_gpu(uadv,vadv,wadv,rhoadv,rhr,np,mp,lp,ih)
   lupdatemulti=.FALSE.
     lvertsplit=.FALSE.
        iprint =1
     nstarttest=nt-1

IF (mype == 0.AND.iprint==1) THEN
        PRINT *,'Begin advection cartesian dwarf testing'
        WRITE (*,120)  nt,n,m,l,nprocx,nprocy,nprocz
120   FORMAT('eulagdwarf will run for ',i0,' iterations on field with size ',i0,' x ',i0,' x ',i0, &
                                                        ' on processor grid',i2,' x ',i2,' x ',i2)
ENDIF
#ifdef TESTING
IF(mype.eq.0.and.n.eq.59.and.m.eq.59.and.l.eq.59) &
print *,'Reference published result for rotating sphere in 59^3 cubical mesh ("Error" L2) norm is 0.0028'
#endif /*TESTING*/

!Prepare advecting velocities and the advected field

   IF (mype == 0) PRINT *,'Begin MPDATA dwarf testing'


   CALL set_initial_tracer(uadv,vadv,wadv,xtracer,rhoadv,np,mp,lp,ih)
   CALL initialize_for_gpu(uadv,vadv,wadv,rhoadv,rhr,np,mp,lp,ih)
   END SUBROUTINE allocate_and_initialize
   

SUBROUTINE deallocate_and_finalize(itime_counter)
USE precisions
USE scratch_datafields, ONLY: xtest
USE testing, ONLY: test_solution_energy,test_solution_error
INTEGER, INTENT(IN) :: itime_counter

!#ifdef TESTING
!   IF(itimecnt.ge.nstarttest) THEN
!@cuf istat=cudaDeviceSynchronize()
      CALL test_solution_energy(xtracer,xtest,itime_counter,'a',np,mp,lp,ih)
      CALL test_solution_error(xtracer,xtest,itime_counter,'c',np,mp,lp,ih)
!   ENDIF
!#endif /*TESTING*/

#ifdef CUDACODE
   CALL cudaProfilerStop();
#endif
!-----------------------------------------------------------------------------------
IF (mype == 0) PRINT *,'End MPDATA dwarf executions'
IF (mype == 0) PRINT *,'Dwarf statistics summary'
CALL ttprt
#if (STATICMEM==0)
   CALL deallocate_geometry_basic_datafields
!  CALL allocate_scratch_datafields
   CALL deallocate_scratch_mpdata_datafields
   CALL deallocate_scratch_advdrv_datafields
   CALL deallocate_scratch_advtst_datafields
!   CALL allocate_scratch_geometry_datafields
#ifdef PUREMPI
  CALL deallocate_scratch_mpi_datafields
#endif
   CALL deallocate_scratch_aux_datafields
#endif /*STATICMEM*/
   END SUBROUTINE deallocate_and_finalize
END MODULE advec_initialize

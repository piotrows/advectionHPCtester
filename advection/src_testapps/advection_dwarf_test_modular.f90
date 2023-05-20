!#define DWARF0
!#define DWARF1
!#define DWARF2
#define DWARF3
!#define DWARF4
!#define DWARF5
!#define DWARF6
!#define DWARF7
!#define DWARF8
!#define DWARF9
!#define DWARF10
!#define DWARF11
!#define DWARF12
PROGRAM advection_dwarf_cartesian_test
   USE precisions
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
#ifdef TESTING
   USE parameters, ONLY: n,m,l
   USE scratch_datafields, ONLY: xtest
   USE testing, ONLY: test_solution_energy,test_solution_error
   USE testing, ONLY: test_solution_ERR2_91,test_solution_ERR1_91
   USE testing, ONLY: test_solution_ERR0_91,test_solution_ERRMIN_91
   USE testing, ONLY: test_solution_ERRMAX_91
#endif /*TESTING*/
#ifdef PNETCDF
   USE parameters, ONLY: npnetstep
   USE mpi_parallel, ONLY: pnet_out_chunk
#endif /*PNETCDF*/
   IMPLICIT NONE
!  TYPE(C_PTR)  :: xtracer_c_ptr,xant_c_ptr
#ifdef PNETCDF
   INTEGER(KIND=iintegers) npnettime 
#endif /*PNETCDF*/
   INTEGER(KIND=iintegers) itimecnt,iprint,itestcnt,maxtestcnt,nstarttest
   INTEGER(KIND=iintegers) nprocx,nprocy,nprocz
   INTEGER(KIND=iintegers),DIMENSION(100) :: opttype_list,algtype_list,tformat_list 
   INTEGER(KIND=iintegers),PARAMETER :: ipoles=0
   INTEGER(KIND=iintegers),PARAMETER :: iflip=0
   CHARACTER(LEN=22) :: pnetvar_list(100)
   LOGICAL lupdatemulti,lvertsplit
   REAL(KIND=euwp) dt_loc,dt_orig
   LOGICAL :: PrintHelp = .FALSE.
   INTEGER istat,iscale,nx,ier
   INTEGER itmin,itmax,itstr,itestset
#if (STATICMEM == 0)
 call get_args_from_commandline
#endif

#if (STATICMEM == 0)
!   CALL set_default_dynamicmem_parameters
#endif 
!  CALL read_commandline_parameters
   itmin=1
   itmax=16
   itstr=1
   testloop: DO itestset=itmin,itmax,itstr
!  IF(mype.eq.0) print *,itestset,itmin,itmax,itstr
   iscale=1
   iscale=itestset
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
   
   IF(itestset.eq.itmin) THEN 
     CALL geomset (.TRUE.)
   ELSE
     CALL geomset (.FALSE.)
   ENDIF 
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
#ifdef PNETCDF
     npnetstep=1
#endif
 
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
   CALL define_list_of_tests
   
   IF (mype == 0) PRINT *,'Begin MPDATA dwarf testing'

 
   CALL set_initial_tracer(uadv,vadv,wadv,xtracer,rhoadv,np,mp,lp,ih)
   CALL initialize_for_gpu(uadv,vadv,wadv,rhoadv,rhr,np,mp,lp,ih)
!@cuf istat=cudaDeviceSynchronize()
!  print *,'RHR',rhr
!  print *,'UADV',uadv
!  print *,'VADV',vadv
!  print *,'WADV',wadv
!  print *,'RHOADV',rhoadv
#ifdef CUDACODE
!   ier = cudaMemAdvise(uadv,sizeof(uadv/8),cudaMemAdviseSetReadMostly,0)
!   ier = cudaMemAdvise(vadv,sizeof(vadv/8),cudaMemAdviseSetReadMostly,0)
!   ier = cudaMemAdvise(wadv,sizeof(wadv/8),cudaMemAdviseSetReadMostly,0)
!   ier = cudaMemAdvise(rhoadv,sizeof(rhoadv/8),cudaMemAdviseSetReadMostly,0)
!   ier = cudaMemAdvise(rhr,sizeof(rhoadv/8),cudaMemAdviseSetReadMostly,0)
   CALL cudaProfilerStart();
#endif
   DO itestcnt=1,maxtestcnt
!@cuf istat=cudaDeviceSynchronize()
     DO itimecnt=1,nt
     IF(mype.eq.0.AND.iprint.eq.1.and.0.eq.1) THEN
       IF(tformat_list(itestcnt).eq.101) WRITE (*, 101,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.102) WRITE (*, 102,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.103) WRITE (*, 103,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.104) WRITE (*, 104,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.105) WRITE (*, 105,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.106) WRITE (*, 106,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.107) WRITE (*, 107,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.108) WRITE (*, 108,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.109) WRITE (*, 109,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.110) WRITE (*, 110,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.111) WRITE (*, 111,advance='no') achar(13),itimecnt
       IF(tformat_list(itestcnt).eq.112) WRITE (*, 112,advance='no') achar(13),itimecnt
     ENDIF
     CALL ttbeg(300)
     CALL advec_dwarf(opttype_list(itestcnt),algtype_list(itestcnt),              &
                    lupdatemulti,lvertsplit,                     &
                    np,mp,lp,ih,ipoles,ibcx,ibcy,ibcz,           &
                    xtracer,bcx,bcy,                             &
                    uadv,vadv,wadv,rhr,rhoadv) !,wadv2,rhr2,rhoadv2,rhr3,rhoadv3)
     CALL ttend(300)
    
!IF (mype == 0.AND.iprint==1) print *,'Completed:',itimecnt,'of',nt
!Optional parallel netcdf I/O
#ifdef PNETCDF
   npnettime=(itimecnt-1)/npnetstep
!IF (mype == 0.AND.iprint==1) print *,'About to write parallel netcdf'
!  IF(MODULO(itimecnt-1,npnetstep)==0)  CALL pnet_out_chunk(pnetvar_list(itestcnt),'mpdat.nc', &
!                                         1,1,1,1,npnettime,0,xtracer,np,mp,lp,ih)
#endif /*PNETCDF*/
#ifdef TESTING
   IF(itimecnt.ge.nstarttest) THEN
!@cuf istat=cudaDeviceSynchronize()
      CALL test_solution_energy(xtracer,xtest,itimecnt,'a',np,mp,lp,ih)
      CALL test_solution_error(xtracer,xtest,itimecnt,'c',np,mp,lp,ih)
   ENDIF
#endif /*TESTING*/
     ENDDO
   ENDDO

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

CALL end_code
   ENDDO testloop
IF (mype == 0) PRINT *,'Goodbye'
CALL end_mpi
100   FORMAT (a,'Upwind gpubc C execution ',i5,' ')
101   FORMAT (a,'MPDATA gauge legacy execution ',i5,' ')
102   FORMAT (a,'MPDATA gauge gpubc execution ',i5,' ')
103   FORMAT (a,'MPDATA twostep gauge gpubc execution ',i5,' ')
104   FORMAT (a,'MPDATA gauge halobc execution ',i5,' ')
105   FORMAT (a,'MPDATA twostep gauge halobc execution ',i5,' ')
106   FORMAT (a,'MPDATA standard legacy execution ',i5,' ')
107   FORMAT (a,'MPDATA standard gpubc execution ',i5,' ')
108   FORMAT (a,'MPDATA twostep standard gpubc execution ',i5,' ')
109   FORMAT (a,'MPDATA standard halobc execution ',i5,' ')
110   FORMAT (a,'MPDATA twostep standard halobc execution ',i5,' ')
111   FORMAT (a,'UPWIND halobc execution ',i5,' ')
112   FORMAT (a,'UPWIND gpubc execution ',i5,' ')

CONTAINS
 SUBROUTINE advec_dwarf_c
 USE, INTRINSIC :: ISO_C_BINDING
 END SUBROUTINE advec_dwarf_c 
 SUBROUTINE define_list_of_tests()
   maxtestcnt=0
#ifdef DWARF0  /*Upwind C*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=0 
   algtype_list(maxtestcnt)=1 
   tformat_list(maxtestcnt)=100 
   pnetvar_list(maxtestcnt)='x_upwind_gpubc'
#endif /*DWARF0*/
#ifdef DWARF1  /*MPDATA gauge legacy*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=0 
   algtype_list(maxtestcnt)=1 
   tformat_list(maxtestcnt)=101 
   pnetvar_list(maxtestcnt)='x_onestep_gauge_legacy'
#endif /*DWARF1*/
#ifdef DWARF2  /*MPDATA gauge HALOBC*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=1 
   algtype_list(maxtestcnt)=1 
   tformat_list(maxtestcnt)=102 
   pnetvar_list(maxtestcnt)='x_onestep_gauge_halobc'
#endif /*DWARF2*/
#ifdef DWARF3  /*TWOSTEP gauge GPUBC*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=4 
   algtype_list(maxtestcnt)=1 
   tformat_list(maxtestcnt)=103 
   pnetvar_list(maxtestcnt)='x_twostep_gauge__gpubc'
#endif /*DWARF3*/
#ifdef DWARF4  /*MPDATA gauge GPUBC*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=2 
   algtype_list(maxtestcnt)=1 
   tformat_list(maxtestcnt)=104 
   pnetvar_list(maxtestcnt)='x_onestep_gauge__gpubc'
#endif /*DWARF4*/
#ifdef DWARF5  /*TWOSTEP gauge HALOBC*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=3 
   algtype_list(maxtestcnt)=1 
   tformat_list(maxtestcnt)=105 
   pnetvar_list(maxtestcnt)='x_twostep_gauge_halobc'
#endif /*DWARF5*/
#ifdef DWARF6  /*MPDATA standard legacy*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=0 
   algtype_list(maxtestcnt)=2 
   tformat_list(maxtestcnt)=106 
   pnetvar_list(maxtestcnt)='x_onestep_stnd_legacy'
#endif /*DWARF6*/
#ifdef DWARF7  /*MPDATA standard GPUBC*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=1 
   algtype_list(maxtestcnt)=2 
   tformat_list(maxtestcnt)=107 
   pnetvar_list(maxtestcnt)='x_onestep_stnd__gpubc'
#endif /*DWARF7*/
#ifdef DWARF8  /*TWOSTEP standard GPUBC*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=4 
   algtype_list(maxtestcnt)=2 
   tformat_list(maxtestcnt)=108 
   pnetvar_list(maxtestcnt)='x_twostep_stnd__gpubc'
#endif /*DWARF8*/
#ifdef DWARF9  /*MPDATA standard HALOBC*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=2 
   algtype_list(maxtestcnt)=2 
   tformat_list(maxtestcnt)=109 
   pnetvar_list(maxtestcnt)='x_twostep_stnd_halobc'
#endif /*DWARF9*/
#ifdef DWARF10  /*TWOSTEP standard HALOBC*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=3 
   algtype_list(maxtestcnt)=2 
   tformat_list(maxtestcnt)=110 
   pnetvar_list(maxtestcnt)='x_twostep_stnd_halobc'
#endif /*DWARF10*/
#ifdef DWARF11  /*UPWIND HALOBC*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=5 
   algtype_list(maxtestcnt)=0 
   tformat_list(maxtestcnt)=111 
   pnetvar_list(maxtestcnt)='x_upwind_halobc'
#endif /*DWARF11*/
#ifdef DWARF12  /*UPWIND GPUBC*/
   maxtestcnt=maxtestcnt+1
   opttype_list(maxtestcnt)=6 
   algtype_list(maxtestcnt)=0 
   tformat_list(maxtestcnt)=112 
   pnetvar_list(maxtestcnt)='x_upwind_gpubc'
#endif /*DWARF12*/

 END SUBROUTINE define_list_of_tests

 SUBROUTINE get_args_from_commandline
 USE argsparser, ONLY: CheckForHelp, ParseArgumentInt, ParseArgumentReal, ParseArgumentLogical, ParseArgumentString
 INTEGER :: i
 INTEGER :: StatusCorrect = 0
 INTEGER :: StatusIncorrect = 0

   iprint=0
   IF (mype == 0.AND.iprint==1) THEN

      write(*,*) n,m,l,ih,nt
      write(*,*) nprocx, nprocy, nprocz
   ENDIF
      PrintHelp = CheckForHelp()
      DO i = 1,2
!     n = ParseArgumentInt("--sizex", .TRUE.,64, StatusCorrect, StatusIncorrect,PrintHelp,&
!             &"--sizex #: specify domain size in direction X")
!     m = ParseArgumentInt("--sizey", .TRUE.,64, StatusCorrect, StatusIncorrect,PrintHelp,&
!             &"--sizey #: specify domain size in direction Y")
!     l = ParseArgumentInt("--sizez", .TRUE.,64, StatusCorrect, StatusIncorrect,PrintHelp,&
!             &"--sizez #: specify domain size in direction Z")
!     nt = ParseArgumentInt("--maxiteration", .TRUE.,100, StatusCorrect, StatusIncorrect,PrintHelp,&
!             &"--maxiteration #: number of timesteps")

      nprocx = ParseArgumentInt("--nprocx", .TRUE.,1, StatusCorrect, StatusIncorrect,PrintHelp,&
              &"--nprocx #: specify number of processes in direction X")

      nprocy = ParseArgumentInt("--nprocy", .TRUE.,1, StatusCorrect, StatusIncorrect,PrintHelp,&
              &"--nprocy #: specify number of processes in direction Y")

      nprocz = ParseArgumentInt("--nprocz", .TRUE.,1, StatusCorrect, StatusIncorrect,PrintHelp,&
              &"--nprocz #: specify number of processes in direction Z")


      IF (PrintHelp) THEN
          STOP 0
      END IF
      IF (StatusIncorrect > 0) THEN
          PrintHelp = .TRUE.
      ELSE
          EXIT
      END IF
      END DO
  END SUBROUTINE get_args_from_commandline

END PROGRAM advection_dwarf_cartesian_test



#ifdef DWARF0
   itestcnt=0
!     IF(tformat_list(itestcnt).eq.100) WRITE (*, 100,advance='no') achar(13),itimecnt
   CALL set_initial_tracer(uadv,vadv,wadv,xtracer,rhoadv,np,mp,lp,ih)
      ALLOCATE (xtracer_c(np+2*ih  ,mp+2*ih  ,lp+2*ih))
      ALLOCATE (    rhr_c(np+2*ih  ,mp+2*ih  ,lp+2*ih))
      ALLOCATE (      h_c(np+2*ih  ,mp+2*ih  ,lp+2*ih))
      ALLOCATE (   uadv_c(np+2*ih+1,mp+2*ih  ,lp+2*ih))
      ALLOCATE (   vadv_c(np+2*ih  ,mp+2*ih+1,lp+2*ih))
      ALLOCATE (   wadv_c(np+2*ih  ,mp+2*ih  ,lp+2*ih+1))
     xtracer_c(:,:,:)=xtracer(:,:,:)
         rhr_c(:,:,:)=    rhr(:,:,:)
           h_c(:,:,:)= rhoadv(:,:,:)
         uadv_c(:,:,:)=  uadv(:,:,:)
         vadv_c(:,:,:)=  vadv(:,:,:)
         wadv_c(:,:,:)=  wadv(:,:,:)
      itimecnt=1
      CALL test_solution_error(xtracer,xtest,itimecnt,'c',np,mp,lp,ih)
     CALL UPWIND_C_TIMELOOP(C_LOC(xtracer_c),C_LOC(uadv_c),C_LOC(vadv_c),C_LOC(wadv_c),  &
                            C_LOC(rhr_c),C_LOC(h_c),   &
                            np,mp,lp,ih,nt)
      CALL test_solution_error(xtracer_c,xtest,itimecnt,'c',np,mp,lp,ih)
      DEALLOCATE (xtracer_c)
      DEALLOCATE (rhr_c)
      DEALLOCATE (h_c)
      DEALLOCATE (uadv_c)
      DEALLOCATE (vadv_c)
      DEALLOCATE (wadv_c)
  INTERFACE 
    SUBROUTINE UPWIND_C_TIMELOOP(xtracer,uadv,vadv,wadv,rhr,h,np,mp,lp,ih,nt) bind(c,name="upwind_c_timeloop") 
!   SUBROUTINE UPWIND_C(xtracer_ptr,xant_ptr) bind(c,name="upwind_c") 
    import :: c_ptr
    import :: c_int
    import :: c_double
    type(c_ptr), VALUE :: xtracer,uadv,vadv,wadv,rhr,h
    integer(c_int), VALUE :: np,mp,lp,ih,nt
!   REAL(C_DOUBLE), DIMENSION(:,:,:) :: xtracer,xant
    END SUBROUTINE UPWIND_C_TIMELOOP
  END INTERFACE
  REAL(C_DOUBLE), DIMENSION (:,:,:),  ALLOCATABLE, TARGET :: xtracer_c,uadv_c,vadv_c,wadv_c, rhr_c,h_c  
#endif

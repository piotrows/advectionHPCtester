!#define DWARF0
!#define DWARF1
!#define DWARF2
!#define DWARF3
!#define DWARF4
!#define DWARF5
!#define DWARF6
!#define DWARF7
!#define DWARF8
!#define DWARF9
!#define DWARF10
!#define DWARF11
#define DWARF12
PROGRAM advection_dwarf_cartesian_test 
   USE advec_interface_spgpu, ONLY: advec_dwarf_interface_spgpu
   USE advec_interface_spgpu, ONLY:   allocate_interface_spgpu  
   USE advec_interface_spgpu, ONLY: deallocate_interface_spgpu  
   INTEGER,DIMENSION(100) :: opttype_list,algtype_list,tformat_list 
   INTEGER,PARAMETER :: ipoles=0
   CHARACTER(LEN=22) :: pnetvar_list(100)
   LOGICAL lupdatemulti,lvertsplit
   LOGICAL :: linitmpi=.TRUE.
   INTEGER itimecnt,nt,itestcnt
   INTEGER  :: nprocx=1 
   INTEGER  :: nprocy=1 
   INTEGER  :: nprocz=1 

   CALL define_list_of_tests
   nt=((556+0))
   itestcnt=1
   lupdatemulti=.FALSE.
   lvertsplit=.FALSE.


   CALL allocate_interface_spgpu(.TRUE.,1,1,1)
   DO itimecnt=1,nt
     CALL advec_dwarf_interface_spgpu(opttype_list(itestcnt), &
                                algtype_list(itestcnt), &
                                lupdatemulti,lvertsplit, ipoles, &                    
                                itimecnt )                    
   ENDDO
   CALL deallocate_interface_spgpu(itimecnt)    
CONTAINS
 SUBROUTINE define_list_of_tests()
   INTEGER maxtestcnt
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


END PROGRAM advection_dwarf_cartesian_test




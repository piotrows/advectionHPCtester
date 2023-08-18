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
!  USE advec_interface_dp, ONLY: advec_dwarf_interface_dp
!  USE advec_interface_dp, ONLY:   allocate_interface_dp  
!  USE advec_interface_dp, ONLY: deallocate_interface_dp  
!  USE advec_interface_sp, ONLY: advec_dwarf_interface_sp
!  USE advec_interface_sp, ONLY:   allocate_interface_sp  
!  USE advec_interface_sp, ONLY: deallocate_interface_sp  
   USE :: iso_c_binding
   integer(c_int), parameter :: rtld_local=0 ! value extracte from the C header file
   integer(c_int), parameter :: rtld_lazy=1 ! value extracte from the C header file
    integer(c_int), parameter :: rtld_now=2 ! value extracte from the C header file
    !
    ! interface to linux API
    interface
        function dlopen(filename,mode) bind(c,name="dlopen")
            ! void *dlopen(const char *filename, int mode);
            use iso_c_binding
            implicit none
            type(c_ptr) :: dlopen
            character(c_char), intent(in) :: filename(*)
            integer(c_int), value :: mode
        end function

        function dlsym(handle,name) bind(c,name="dlsym")
            ! void *dlsym(void *handle, const char *name);
            use iso_c_binding
            implicit none
            type(c_funptr) :: dlsym
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
        end function

        function dlclose(handle) bind(c,name="dlclose")
            ! int dlclose(void *handle);
            use iso_c_binding
            implicit none
            integer(c_int) :: dlclose
            type(c_ptr), value :: handle
        end function
    end interface

    abstract interface
        subroutine compute_proc (iopt, ialg, lupdate, lvsplit, ipoles, itimecnt) bind(c)
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(in) ::iopt, ialg, ipoles, itimecnt 
            logical(c_bool), intent(in) ::lupdate, lvsplit 
        end subroutine compute_proc
    end interface
    abstract interface
        subroutine allocate_proc (lmpiini, nprocx, nprocy, nprocz) bind(c)
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(in) :: nprocx, nprocy, nprocz 
            logical(c_bool), intent(in) ::lmpiini 
        end subroutine allocate_proc
    end interface
    abstract interface
        subroutine deallocate_proc (itimecnt) bind(c)
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(in) :: itimecnt 
        end subroutine deallocate_proc
    end interface
    type(c_funptr) :: allocate_addr
    type(c_funptr) :: deallocate_addr
    type(c_funptr) :: compute_addr
    type(c_ptr)    :: libhandle
    character(256) :: pName, lName
    procedure(compute_proc), bind(c), pointer :: advec_dwarf_interface 
    procedure(allocate_proc), bind(c), pointer :: allocate_interface 
    procedure(deallocate_proc), bind(c), pointer :: deallocate_interface 

   INTEGER,DIMENSION(100) :: opttype_list,algtype_list,tformat_list 
   INTEGER,PARAMETER :: ipoles=0
   CHARACTER(LEN=22) :: pnetvar_list(100)
   LOGICAL(c_bool) :: linitmpi=.TRUE.
   LOGICAL(c_bool) :: lupdatemulti=.FALSE.
   LOGICAL(c_bool) :: lvertsplit=.FALSE.

   INTEGER(c_int) itimecnt,nt,itestcnt,ierr
   CALL define_list_of_tests
   nt=((556+0))
   itestcnt=1
   libhandle=dlopen("./lib/libadvection_interface_dp.so"//c_null_char, RTLD_LOCAL)
    if (.not. c_associated(libhandle))then
        print*, 'Unable to load DLL ./lib/libadvection_interface_dp.so'
        stop
    end if 
    compute_addr=dlsym(libhandle, "advec_dwarf_interface_dp"//c_null_char)
    if (.not. c_associated(compute_addr))then
        write(*,*) 'Unable to load the procedure advec_dwarf_interface_dp'
        stop
    end if
    call c_f_procpointer( compute_addr, advec_dwarf_interface )
    allocate_addr=dlsym(libhandle, "allocate_interface_dp"//c_null_char)
    if (.not. c_associated(allocate_addr))then
        write(*,*) 'Unable to load the procedure allocate_interface_dp'
        stop
    end if
    call c_f_procpointer( allocate_addr, allocate_interface )
    deallocate_addr=dlsym(libhandle, "deallocate_interface_dp"//c_null_char)
    if (.not. c_associated(deallocate_addr))then
        write(*,*) 'Unable to load the procedure deallocate_interface_dp'
        stop
    end if
    call c_f_procpointer( deallocate_addr, deallocate_interface )

   CALL allocate_interface(linitmpi,1,1,1)
   DO itimecnt=1,nt
     CALL advec_dwarf_interface(opttype_list(itestcnt), &
                                algtype_list(itestcnt), &
                                lupdatemulti,lvertsplit, ipoles, &
                                itimecnt )                    
   ENDDO
   CALL  deallocate_interface(itimecnt)    
   ierr=dlclose(libhandle)
   itestcnt=1
   libhandle=dlopen("./lib/libadvection_interface_sp.so"//c_null_char, RTLD_LOCAL)
    if (.not. c_associated(libhandle))then
        print*, 'Unable to load DLL ./lib/libadvection_interface_sp.so'
        stop
    end if 
    compute_addr=dlsym(libhandle, "advec_dwarf_interface_sp"//c_null_char)
    if (.not. c_associated(compute_addr))then
        write(*,*) 'Unable to load the procedure advec_dwarf_interface_sp'
        stop
    end if
    call c_f_procpointer( compute_addr, advec_dwarf_interface )
    allocate_addr=dlsym(libhandle, "allocate_interface_sp"//c_null_char)
    if (.not. c_associated(allocate_addr))then
        write(*,*) 'Unable to load the procedure allocate_interface_sp'
        stop
    end if
    call c_f_procpointer( allocate_addr, allocate_interface )
    deallocate_addr=dlsym(libhandle, "deallocate_interface_sp"//c_null_char)
    if (.not. c_associated(deallocate_addr))then
        write(*,*) 'Unable to load the procedure deallocate_interface_sp'
        stop
    end if
    call c_f_procpointer( deallocate_addr, deallocate_interface )
   CALL allocate_interface(linitmpi,1,1,1)
   DO itimecnt=1,nt
     CALL advec_dwarf_interface(opttype_list(itestcnt), &
                                algtype_list(itestcnt), &
                                lupdatemulti,lvertsplit, ipoles, &
                                itimecnt )                    
   ENDDO
   CALL deallocate_interface(itimecnt)    
!  itestcnt=1
!  CALL allocate_interface_sp(.FALSE.,1,1,1)
!  DO itimecnt=1,nt
!    CALL advec_dwarf_interface_sp(opttype_list(itestcnt), &
!                               algtype_list(itestcnt), &
!                               lupdatemulti,lvertsplit, ipoles, &                    
!                               itimecnt )                    
!  ENDDO
!  CALL deallocate_interface_sp(itimecnt)    
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




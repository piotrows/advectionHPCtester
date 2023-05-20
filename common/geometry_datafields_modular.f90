#if (STATICMEM == 1)
MODULE geometry_datafields
   USE precisions
   USE parameters, ONLY: n,m,l,np,mp,lp,ih
   IMPLICIT NONE
!Geometry: topological and metric information
   REAL (KIND=euwp) x(n),y(m),z(l)
   REAL (KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih) :: &
     xcr,ycr,sina,cosa,tnga,sinx,cosx
   REAL (KIND=euwp),DIMENSION(np,mp) :: &
     fcr2,fcr3
   REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
     zcr, &
     g11,g12,g13, &
     g21,g22,g23, &
     g33, &
     gac,gmm, &
     gf11, gf12, gf13, &
     gf21, gf32, gf23, &
     gf31, gf22, gf33
     gf11V, gf12V, gf13V, &
     gf21V, gf32V, gf23V, &
     gf31V, gf22V, gf33
   REAL(KIND=euwp),DIMENSION(np,mp,lp) :: &
     gmmi,  &
     gaci,  &
     g33i,  &
     g11g22Mg12g21i
   REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
     sg11,sg12,sg13, &
     sg21,sg22,sg23, &
     sg33


    REAL(KIND=euwp), DIMENSION(1:mp,1:lp,2) ::  tnxx, tnxy, tnxz, sfflx
    REAL(KIND=euwp), DIMENSION(1:np,1:lp,2) ::  tnyx, tnyy, tnyz, sffly
    REAL(KIND=euwp), DIMENSION(1:np,1:mp,2) ::  tnzx, tnzy, tnzz, sfflz

    REAL(KIND=euwp), DIMENSION(1:mp,1:lp,2) ::  tnxxV, tnxyV, tnxzV
    REAL(KIND=euwp), DIMENSION(1:np,1:lp,2) ::  tnyxV, tnyyV, tnyzV
    REAL(KIND=euwp), DIMENSION(1:np,1:mp,2) ::  tnzxV, tnzyV, tnzzV


END MODULE geometry_datafields
#else /*STATICMEM*/ 
MODULE geometry_datafields
   USE precisions
   IMPLICIT NONE 
#include "../advection/src_algorithms/defines.inc"

!Geometry: topological and metric information
   REAL_euwp,DIMENSION(:), ALLOCATABLE :: x,y,z 
   REAL_euwp,DIMENSION(:,:), ALLOCATABLE :: & !1-ih:np+ih,1-ih:mp+ih 
     xcr,ycr,sina,cosa,tnga,sinx,cosx,fcr2,fcr3
   REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
     zcr,g11,g12,g13,g21,g22,g23,g33,gac,gmm
   REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
     g33i,gaci,gmmi,g11g22Mg12g21i
   REAL(KIND=euwp), DIMENSION(:,:,:), ALLOCATABLE :: gf11, gf12, gf13
   REAL(KIND=euwp), DIMENSION(:,:,:), ALLOCATABLE :: gf21, gf22, gf23
   REAL(KIND=euwp), DIMENSION(:,:,:), ALLOCATABLE :: gf31, gf32, gf33
   REAL(KIND=euwp), DIMENSION(:,:,:), ALLOCATABLE :: gf11V, gf12V, gf13V
   REAL(KIND=euwp), DIMENSION(:,:,:), ALLOCATABLE :: gf21V, gf22V, gf23V
   REAL(KIND=euwp), DIMENSION(:,:,:), ALLOCATABLE :: gf31V, gf32V, gf33V
   REAL(KIND=euwp), DIMENSION(:,:,:), ALLOCATABLE :: sg11,sg12,sg13, &
                                                     sg21,sg22,sg23, &
                                                               sg33

   REAL(KIND=euwp), DIMENSION(:,:,:), ALLOCATABLE   :: tnxx, tnxy, tnxz, sfflx, &
                                                       tnyx, tnyy, tnyz, sffly, &
                                                       tnzx, tnzy, tnzz, sfflz, &
                                                       tnxxV, tnxyV, tnxzV, &
                                                       tnyxV, tnyyV, tnyzV, &
                                                       tnzxV, tnzyV, tnzzV


   CONTAINS

   SUBROUTINE allocate_geometry_basic_datafields
   USE parameters, ONLY: n,m,l,np,mp,lp,ih
   USE, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_signaling_NAN
   IMPLICIT NONE
   REAL(KIND=euwp) real_initial_value
     INTEGER ierr,ierrtot
!Geometry: topological and metric information
     ierrtot=0
!    print *,'EULAG allocate extent',np,mp,lp,ih
     ALLOCATE ( zcr(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g11(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g12(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g13(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g21(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g22(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g23(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g33(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g33i(np,mp,lp)                         ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gac(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gmm(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gaci(np,mp,lp)                         ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gmmi(np,mp,lp)                         ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g11g22Mg12g21i(np,mp,lp)               ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( xcr(1-ih:np+ih, 1-ih:mp+ih)            ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( ycr(1-ih:np+ih, 1-ih:mp+ih)            ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( sina(1-ih:np+ih, 1-ih:mp+ih)           ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( cosa(1-ih:np+ih, 1-ih:mp+ih)           ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( tnga(1-ih:np+ih, 1-ih:mp+ih)           ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( sinx(1-ih:np+ih, 1-ih:mp+ih)           ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( cosx(1-ih:np+ih, 1-ih:mp+ih)           ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( fcr2(np,mp)                            ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( fcr3(np,mp)                            ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( x(n)                                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( y(m)                                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( z(l)                                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr

     IF(ierrtot/=0) print *,'Allocation of geometry datafields status: total number of errors is: ',ierrtot
     IF(ierrtot/=0) STOP 'Lack of memory for geometry datafields' 

#ifdef INITIALIZEVARIABLES
  real_initial_value = IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
! real_initial_value = 0. !IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
     zcr(:,:,:)=real_initial_value 
     g11(:,:,:)=real_initial_value 
     g12(:,:,:)=real_initial_value 
     g13(:,:,:)=real_initial_value 
     g21(:,:,:)=real_initial_value 
     g22(:,:,:)=real_initial_value 
     g23(:,:,:)=real_initial_value 
     g33(:,:,:)=real_initial_value 
     gac(:,:,:)=real_initial_value 
    gaci(:,:,:)=real_initial_value
     gmm(:,:,:)=real_initial_value 
    gmmi(:,:,:)=real_initial_value 
    g11g22Mg12g21i(:,:,:)=real_initial_value
     xcr(:,:)  =real_initial_value
     ycr(:,:)  =real_initial_value
    sina(:,:)  =real_initial_value
    cosa(:,:)  =real_initial_value
    tnga(:,:)  =real_initial_value
    sinx(:,:)  =real_initial_value
    cosx(:,:)  =real_initial_value
    fcr2(:,:)  =real_initial_value
    fcr3(:,:)  =real_initial_value
         x(:)  =real_initial_value
         y(:)  =real_initial_value
         z(:)  =real_initial_value
#endif
   END SUBROUTINE allocate_geometry_basic_datafields
   SUBROUTINE allocate_geometry_sgs_datafields
   USE parameters, ONLY: n,m,l,np,mp,lp,ih
   USE, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_signaling_NAN
   IMPLICIT NONE
   REAL(KIND=euwp) real_initial_value
     INTEGER ierr,ierrtot
!Geometry: topological and metric information
     ierrtot=0
!    print *,'EULAG allocate extent',np,mp,lp,ih
     ALLOCATE ( sg11(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( sg12(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( sg13(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( sg21(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( sg22(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( sg23(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( sg33(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf11(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf12(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf13(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf21(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf22(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf23(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf31(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf32(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf33(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf11V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf12V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf13V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf21V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf22V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf23V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf31V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf32V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gf33V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr

     ALLOCATE (      tnxx(1:mp,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr; 
     ALLOCATE (      tnxy(1:mp,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (      tnxz(1:mp,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     sfflx(1:mp,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (      tnyx(1:np,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (      tnyy(1:np,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (      tnyz(1:np,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     sffly(1:np,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (      tnzx(1:np,1:mp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (      tnzy(1:np,1:mp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (      tnzz(1:np,1:mp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     sfflz(1:np,1:mp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
 
     ALLOCATE (     tnxxV(1:mp,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     tnxyV(1:mp,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     tnxzV(1:mp,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     tnyxV(1:np,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     tnyyV(1:np,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     tnyzV(1:np,1:lp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     tnzxV(1:np,1:mp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     tnzyV(1:np,1:mp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     ALLOCATE (     tnzzV(1:np,1:mp,2)                   ,STAT=ierr ) ;ierrtot=ierrtot+ierr;

     IF(ierrtot/=0) print *,'Allocation of geometry sgs datafields status: total number of errors is: ',ierrtot
     IF(ierrtot/=0) STOP 'Lack of memory for geometry sgs datafields' 

#ifdef INITIALIZEVARIABLES
  real_initial_value = IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
! real_initial_value = 0. !IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
    sg11=real_initial_value 
    sg12=real_initial_value 
    sg13=real_initial_value 
    sg21=real_initial_value 
    sg22=real_initial_value 
    sg23=real_initial_value 
    sg33=real_initial_value 
    gf11=real_initial_value 
    gf12=real_initial_value 
    gf13=real_initial_value 
    gf21=real_initial_value 
    gf22=real_initial_value 
    gf23=real_initial_value 
    gf31=real_initial_value 
    gf32=real_initial_value 
    gf33=real_initial_value  
    gf11V=real_initial_value 
    gf12V=real_initial_value 
    gf13V=real_initial_value 
    gf21V=real_initial_value 
    gf22V=real_initial_value 
    gf23V=real_initial_value 
    gf31V=real_initial_value 
    gf32V=real_initial_value 
    gf33V=real_initial_value  

    tnxx=real_initial_value
    tnxy=real_initial_value
    tnxz=real_initial_value
    sfflx=real_initial_value
    tnyx=real_initial_value
    tnyy=real_initial_value
    tnyz=real_initial_value
    sffly=real_initial_value
    tnzx=real_initial_value
    tnzy=real_initial_value
    tnzz=real_initial_value
    sfflz=real_initial_value
    tnxxV=real_initial_value
    tnxyV=real_initial_value
    tnxzV=real_initial_value
    tnyxV=real_initial_value
    tnyyV=real_initial_value
    tnyzV=real_initial_value
    tnzxV=real_initial_value
    tnzyV=real_initial_value
    tnzzV=real_initial_value
#endif
   END SUBROUTINE allocate_geometry_sgs_datafields

   SUBROUTINE deallocate_geometry_sgs_datafields
     INTEGER ierr,ierrtot
     ierrtot=0
     IF(ALLOCATED(sg11 ))  DEALLOCATE ( sg11, STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(sg12 ))  DEALLOCATE ( sg12, STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(sg13 ))  DEALLOCATE ( sg13, STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(sg21 ))  DEALLOCATE ( sg21, STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(sg22 ))  DEALLOCATE ( sg22, STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(sg23 ))  DEALLOCATE ( sg23, STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(sg33 ))  DEALLOCATE ( sg33, STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf11 ))  DEALLOCATE (  gf11,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf12 ))  DEALLOCATE (  gf12,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf13 ))  DEALLOCATE (  gf13,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf21 ))  DEALLOCATE (  gf21,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf22 ))  DEALLOCATE (  gf22,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf23 ))  DEALLOCATE (  gf23,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf31 ))  DEALLOCATE (  gf31,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf32 ))  DEALLOCATE (  gf32,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf33 ))  DEALLOCATE (  gf33,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf11V))  DEALLOCATE ( gf11V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf12V))  DEALLOCATE ( gf12V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf13V))  DEALLOCATE ( gf13V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf21V))  DEALLOCATE ( gf21V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf22V))  DEALLOCATE ( gf22V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf23V))  DEALLOCATE ( gf23V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf31V))  DEALLOCATE ( gf31V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf32V))  DEALLOCATE ( gf32V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gf33V))  DEALLOCATE ( gf33V,STAT=ierr ) ;ierrtot=ierrtot+ierr

     IF(ALLOCATED(tnxx )) DEALLOCATE (   tnxx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnxy )) DEALLOCATE (   tnxy,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnxz )) DEALLOCATE (   tnxz,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(sfflx)) DEALLOCATE (  sfflx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnyx )) DEALLOCATE (   tnyx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnyy )) DEALLOCATE (   tnyy,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnyz )) DEALLOCATE (   tnyz,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(sffly)) DEALLOCATE (  sffly,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnzx )) DEALLOCATE (   tnzx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnzy )) DEALLOCATE (   tnzy,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnzz )) DEALLOCATE (   tnzz,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(sfflz)) DEALLOCATE (  sfflz,STAT=ierr ) ;ierrtot=ierrtot+ierr;

     IF(ALLOCATED(tnxxV)) DEALLOCATE (  tnxxV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnxyV)) DEALLOCATE (  tnxyV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnxzV)) DEALLOCATE (  tnxzV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnyxV)) DEALLOCATE (  tnyxV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnyyV)) DEALLOCATE (  tnyyV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnyzV)) DEALLOCATE (  tnyzV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnzxV)) DEALLOCATE (  tnzxV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnzyV)) DEALLOCATE (  tnzyV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
     IF(ALLOCATED(tnzzV)) DEALLOCATE (  tnzzV,STAT=ierr ) ;ierrtot=ierrtot+ierr;

     IF(ierrtot/=0) print *,'Deallocation of geometry sgs datafields status: total number of errors is: ',ierrtot

   END SUBROUTINE deallocate_geometry_sgs_datafields

   SUBROUTINE deallocate_geometry_basic_datafields
     INTEGER ierr,ierrtot
     ierrtot=0
     IF(ALLOCATED(zcr  ))  DEALLOCATE ( zcr , STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(g11  ))  DEALLOCATE ( g11 , STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(g12  ))  DEALLOCATE ( g12 , STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(g13  ))  DEALLOCATE ( g13 , STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(g21  ))  DEALLOCATE ( g21 , STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(g22  ))  DEALLOCATE ( g22 , STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(g23  ))  DEALLOCATE ( g23 , STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(g33  ))  DEALLOCATE ( g33 , STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(g33i ))  DEALLOCATE ( g33i, STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gac  ))  DEALLOCATE ( gac,  STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gmm  ))  DEALLOCATE ( gmm,  STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gaci ))  DEALLOCATE ( gaci, STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(gmmi ))  DEALLOCATE ( gmmi, STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(g11g22Mg12g21i)) &
                           DEALLOCATE ( g11g22Mg12g21i,&
                                              STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(xcr  ))  DEALLOCATE (   xcr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(ycr  ))  DEALLOCATE (   ycr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(sina ))  DEALLOCATE (  sina,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(cosa ))  DEALLOCATE (  cosa,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(tnga ))  DEALLOCATE (  tnga,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(sinx ))  DEALLOCATE (  sinx,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(cosx ))  DEALLOCATE (  cosx,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(fcr2 ))  DEALLOCATE (  fcr2,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(fcr3 ))  DEALLOCATE (  fcr3,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(x    ))  DEALLOCATE (     x,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(y    ))  DEALLOCATE (     y,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(z    ))  DEALLOCATE (     z,STAT=ierr ) ;ierrtot=ierrtot+ierr

     IF(ierrtot/=0) print *,'Deallocation of geometry basic datafields status: total number of errors is: ',ierrtot
   END SUBROUTINE deallocate_geometry_basic_datafields

END MODULE geometry_datafields
#endif /*STATICMEM*/
#if(1==0)
   SUBROUTINE allocate_geometry_datafields
   USE parameters, ONLY: n,m,l,np,mp,lp,ih
   USE, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_signaling_NAN
   IMPLICIT NONE
   REAL(KIND=euwp) real_initial_value
     INTEGER ierr,ierrtot
!Geometry: topological and metric information
     ierrtot=0
!    print *,'EULAG allocate extent',np,mp,lp,ih
     ALLOCATE ( zcr(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g11(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g12(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g13(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g21(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g22(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g23(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g33(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g33i(np,mp,lp),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (sg11(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (sg12(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (sg13(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (sg21(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (sg22(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (sg23(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (sg33(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gac(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gmm(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gaci(np,mp,lp),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( gmmi(np,mp,lp),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( g11g22Mg12g21i(np,mp,lp),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( xcr(1-ih:np+ih, 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( ycr(1-ih:np+ih, 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( sina(1-ih:np+ih, 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( cosa(1-ih:np+ih, 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( tnga(1-ih:np+ih, 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( sinx(1-ih:np+ih, 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( cosx(1-ih:np+ih, 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( fcr2(np,mp),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( fcr3(np,mp),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( x(n),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( y(m),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( z(l),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf11(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf12(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf13(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf21(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf22(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf23(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf31(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf32(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf33(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf11V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf12V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf13V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf21V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf22V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf23V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf31V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf32V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( gf33V(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr

    ALLOCATE (      tnxx(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr; 
    ALLOCATE (      tnxy(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnxz(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnyx(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnyy(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnyz(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnzx(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnzy(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnzz(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;

    ALLOCATE (     tnxxV(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnxyV(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnxzV(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnyxV(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnyyV(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnyzV(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnzxV(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnzyV(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnzzV(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   IF(ierrtot/=0) print *,'initialize status: total number of errors is: ',ierrtot

#ifdef INITIALIZEVARIABLES
  real_initial_value = IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
! real_initial_value = 0. !IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
     x=real_initial_value
     y=real_initial_value
     z=real_initial_value
     zcr=real_initial_value 
     g11=real_initial_value 
     g12=real_initial_value 
     g13=real_initial_value 
     g21=real_initial_value 
     g22=real_initial_value 
     g23=real_initial_value 
     g33=real_initial_value 
    sg11=real_initial_value 
    sg12=real_initial_value 
    sg13=real_initial_value 
    sg21=real_initial_value 
    sg22=real_initial_value 
    sg23=real_initial_value 
    sg33=real_initial_value 
     gac=real_initial_value 
     gmm=real_initial_value 
     xcr=real_initial_value
     ycr=real_initial_value
    sina=real_initial_value
    cosa=real_initial_value
    tnga=real_initial_value
    sinx=real_initial_value
    cosx=real_initial_value
    fcr2=real_initial_value
    fcr3=real_initial_value
    g11g22Mg12g21i=real_initial_value
    gaci=real_initial_value
    gf11=real_initial_value 
    gf12=real_initial_value 
    gf13=real_initial_value 
    gf21=real_initial_value 
    gf22=real_initial_value 
    gf23=real_initial_value 
    gf31=real_initial_value 
    gf32=real_initial_value 
    gf33=real_initial_value  
    gf11V=real_initial_value 
    gf12V=real_initial_value 
    gf13V=real_initial_value 
    gf21V=real_initial_value 
    gf22V=real_initial_value 
    gf23V=real_initial_value 
    gf31V=real_initial_value 
    gf32V=real_initial_value 
    gf33V=real_initial_value  

    tnxx=real_initial_value
    tnxy=real_initial_value
    tnxz=real_initial_value
    tnyx=real_initial_value
    tnyy=real_initial_value
    tnyz=real_initial_value
    tnzx=real_initial_value
    tnzy=real_initial_value
    tnzz=real_initial_value
    tnxxV=real_initial_value
    tnxyV=real_initial_value
    tnxzV=real_initial_value
    tnyxV=real_initial_value
    tnyyV=real_initial_value
    tnyzV=real_initial_value
    tnzxV=real_initial_value
    tnzyV=real_initial_value
    tnzzV=real_initial_value
#endif
   END SUBROUTINE allocate_geometry_datafields

   SUBROUTINE deallocate_geometry_datafields
     INTEGER ierr,ierrtot
     ierrtot=0
     DEALLOCATE ( zcr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( g11,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( g12,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( g13,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( g21,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( g22,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( g23,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( g33,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( g33i ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( sg11,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( sg12,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( sg13,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( sg21,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( sg22,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( sg23,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( sg33,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gac,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gmm,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gaci,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gmmi ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( g11g22Mg12g21i ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( xcr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( ycr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( sina,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(cosa)) DEALLOCATE ( cosa,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(tnga)) DEALLOCATE ( tnga,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( sinx,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( cosx,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( fcr2,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( fcr3,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( x,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( y,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( z,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf11,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf12,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf13,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf21,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf22,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf23,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf31,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf32,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf33,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf11V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf12V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf13V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf21V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf22V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf23V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf31V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf32V,STAT=ierr ) ;ierrtot=ierrtot+ierr
     DEALLOCATE ( gf33V,STAT=ierr ) ;ierrtot=ierrtot+ierr

    DEALLOCATE (      tnxx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (      tnxy,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (      tnxz,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (      tnyx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (      tnyy,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (      tnyz,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (      tnzx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (      tnzy,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (      tnzz,STAT=ierr ) ;ierrtot=ierrtot+ierr;

    DEALLOCATE (     tnxxV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (     tnxyV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (     tnxzV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (     tnyxV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (     tnyyV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (     tnyzV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (     tnzxV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (     tnzyV,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE (     tnzzV,STAT=ierr ) ;ierrtot=ierrtot+ierr;

   IF(ierrtot/=0) print *,'geometry deallocate status: total number of errors is: ',ierrtot

   END SUBROUTINE deallocate_geometry_datafields
#endif

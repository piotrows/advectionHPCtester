#if (STATICMEM == 1)
MODULE scratch_datafields
   USE precisions
   USE parameters, ONLY: np,mp,lp,ih,msz,mucnt,lrd
!MPI scratch storage
!Buffers for I/O  here are np+1 instead of np, because at times
!in EULAG we use variables of size 1-ih:np+1+ih

      REAL(KIND=euwp) :: &
      tmpysnd1((np+1)*ih*(lp+1),0:msz),                                 &
      tmpysnd2((np+1)*ih*(lp+1),0:msz),                                 &
      tmpxsnd1((mp+1+2*ih)*ih*(lp+1),0:msz),                            &
      tmpxsnd2((mp+1+2*ih)*ih*(lp+1),0:msz),                            &
      tmpzsnd1((np+1+2*ih)*ih*(mp+1+2*ih),0:msz),                       &
      tmpzsnd2((np+1+2*ih)*ih*(mp+1+2*ih),0:msz),                       &
      tmpyrcv1((np+1)*ih*(lp+1),0:msz),                                 &
      tmpyrcv2((np+1)*ih*(lp+1),0:msz),                                 &
      tmpxrcv1((mp+1+2*ih)*ih*(lp+1),0:msz),                            &
      tmpxrcv2((mp+1+2*ih)*ih*(lp+1),0:msz),                            &
      tmpzrcv1((np+1+2*ih)*ih*(mp+1+2*ih),0:msz),                       &
      tmpzrcv2((np+1+2*ih)*ih*(mp+1+2*ih),0:msz)
      REAL(KIND=euwp) :: &
      tmpysnd1m(mucnt*(np+1)*ih*(lp+1)),                                 &
      tmpysnd2m(mucnt*(np+1)*ih*(lp+1)),                                 &
      tmpxsnd1m(mucnt*(mp+1+2*ih)*ih*(lp+1)),                            &
      tmpxsnd2m(mucnt*(mp+1+2*ih)*ih*(lp+1)),                            &
      tmpzsnd1m(mucnt*(np+1+2*ih)*ih*(mp+1+2*ih)),                       &
      tmpzsnd2m(mucnt*(np+1+2*ih)*ih*(mp+1+2*ih)),                       &
      tmpyrcv1m(mucnt*(np+1)*ih*(lp+1)),                                 &
      tmpyrcv2m(mucnt*(np+1)*ih*(lp+1)),                                 &
      tmpxrcv1m(mucnt*(mp+1+2*ih)*ih*(lp+1)),                            &
      tmpxrcv2m(mucnt*(mp+1+2*ih)*ih*(lp+1)),                            &
      tmpzrcv1m(mucnt*(np+1+2*ih)*ih*(mp+1+2*ih)),                       &
      tmpzrcv2m(mucnt*(np+1+2*ih)*ih*(mp+1+2*ih))
!MPDATA scratch storage
   REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
      cp,cn,xstore,xant,rhr,rhoadv,rhr2,rhoadv2,rhr3,rhoadv3
   REAL(KIND=euwp),DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+ih  ) :: uadv
   REAL(KIND=euwp),DIMENSION(1-ih:np+ih  , 1-ih:mp+1+ih, 1-ih:lp+ih  ) :: vadv
   REAL(KIND=euwp),DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+1+ih) :: wadv,wadv2
   REAL(KIND=euwp) :: &
      v1(1-ih:np+ih+1,1-ih:mp+ih,  1-ih:lp+ih  ), &
      v2(1-ih:np+ih,  1-ih:mp+ih+1,1-ih:lp+ih  ), &
      v3(1-ih:np+ih,  1-ih:mp+ih,  1-ih:lp+ih+1), &
      f1(1-ih:np+ih+1,1-ih:mp+ih,  1-ih:lp+ih  ), &
      f2(1-ih:np+ih,  1-ih:mp+ih+1,1-ih:lp+ih  ), &
      f3(1-ih:np+ih,  1-ih:mp+ih,  1-ih:lp+ih+1), &
      f1ns(1-ih:np+ih+1,1-ih:mp+ih,  1-ih:lp+ih  ), &
      f2ns(1-ih:np+ih,  1-ih:mp+ih+1,1-ih:lp+ih  ), &
      f3ns(1-ih:np+ih,  1-ih:mp+ih,  1-ih:lp+ih+1)
   REAL(KIND=euwp) :: &
      u1_out(1-ih:np+ih+1,1-ih:mp+ih,  1-ih:lp+ih  ), &
      u2_out(1-ih:np+ih,  1-ih:mp+ih+1,1-ih:lp+ih  ), &
   REAL(KIND=euwp),DIMENSION(np,mp,lp) :: &
      mx,mn
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih) :: &
  u3gnd_store,u3sky_store 
   REAL(KIND=euwp), DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)  :: &  !twostep scratch
     thant,uant,vant,want,thtant,pstrant,qvant,qcant,qrant,qiant,qsant,qgant 
!Geometry scratch storage
   REAL(KIND=euwp), DIMENSION(1-ih:np+ih,1-ih:mp+ih) :: &
    esxb, &
    esyb, &
    dsxb, &
    dsyb, &
    strxx, &
    strxy, &
    stryx, &
    stryy

  REAL(KIND=euwp), DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: &
    csxb, &
    csyb, &
    cszb, &
    strzx, &
    strzy, &
    strzz, &
    scr1
!solver scratch storage
  REAL(KIND=euwp), DIMENSION(np,mp,lp) :: &
   pexc,hlmc
  REAL(KIND=euwp), DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: &
    cf_etainv,fimy,su,sv,sw,px,py,pz,pfx,pfy,pfz,px0,py0,px1,py1,px2,py2
! GCR fields
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih, lrd) :: &
   x,ax
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
  r,qr,ar 
! Precon fields
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
  pr
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
  rax_ltmp,axax_ltmp,rlsqr_ltmp,axar_ltmp
!General scratch 
  REAL(KIND=euwp), DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),PARAMETER :: &
  scratch_ones=1._euwp
  REAL(KIND=euwp), DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),PARAMETER :: &
  scratch_zeros=0._euwp

END MODULE scratch_datafields
#else /*STATICMEM*/ 
MODULE scratch_datafields
#include "../advection/src_algorithms/defines.inc"
   USE precisions
!MPI scratch storage
!Buffers for I/O  here are np+1 instead of np, because at times
!in EULAG we use variables of size 1-ih:np+1+ih

      DEV_REAL_euwp, DIMENSION (:,:), ALLOCATABLE :: &
      tmpysnd1,tmpysnd2,tmpyrcv1,tmpyrcv2, & !((np+1)*ih*(lp+1),0:msz)
      tmpxsnd1,tmpxsnd2,tmpxrcv1,tmpxrcv2, & !((mp+1+2*ih)*ih*(lp+1),0:msz)
      tmpzsnd1,tmpzsnd2,tmpzrcv1,tmpzrcv2    !((np+1+2*ih)*ih*(mp+1+2*ih),0:msz)

      REAL(KIND=euwp), DIMENSION (:), ALLOCATABLE :: &
      tmpysnd1m,tmpysnd2m,tmpyrcv1m,tmpyrcv2m, & 
      tmpxsnd1m,tmpxsnd2m,tmpxrcv1m,tmpxrcv2m, &
      tmpzsnd1m,tmpzsnd2m,tmpzrcv1m,tmpzrcv2m   

!MPDATA scratch storage
   REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
      rhr,rhr2,rhr3
   DEV_REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
      cp,cn,xant !(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
   REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
      xstore !(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
   REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
      rhoadv,rhoadv2,rhoadv3
   REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
      uadv,vadv,wadv,wadv2 !(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
   DEV_REAL_euwp, DIMENSION(:,:,:),  ALLOCATABLE :: &
      v1,f1,f1ns, & !(1-ih:np+ih+1,1-ih:mp+ih,  1-ih:lp+ih  )
      v2,f2,f2ns, & !(1-ih:np+ih,  1-ih:mp+ih+1,1-ih:lp+ih  )
      v3,f3,f3ns    !(1-ih:np+ih,  1-ih:mp+ih,  1-ih:lp+ih+1)
   REAL_euwp, DIMENSION(:,:,:),  ALLOCATABLE :: &
      u1_out, & !(1-ih:np+ih+1,1-ih:mp+ih,  1-ih:lp+ih  )
      u2_out, & !(1-ih:np+ih,  1-ih:mp+ih+1,1-ih:lp+ih  )
      u3_out    !(1-ih:np+ih,  1-ih:mp+ih,  1-ih:lp+ih+1)
   REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
      mx,mn !(np,mp,lp)
   REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &  !twostep scratch
     thant,uant,vant,want,thtant,pstrant,qvant,qcant,qrant,qiant,qsant,qgant 

   REAL_euwp, DIMENSION(:,:), ALLOCATABLE :: &
    u3gnd_store,u3sky_store

  REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: hi,h

   REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: & 
                    xtracer,xtracer_out,xtest,&  !1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih
                          bcx                         ,&  !(mp,lp, 2)
                          bcy                             !(np,lp, 2)

   

!Geometry scratch storage
   REAL_euwp, DIMENSION(:,:), ALLOCATABLE :: &
    esxb,esyb,dsxb,dsyb, &
    strxx,strxy,stryx,stryy

  REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
    csxb,csyb, cszb,    &
    strzx, strzy,strzz, &
    scr1
!solver scratch storage
  REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
   pexc,hlmc
  REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
    cf_etainv,fimy,su,sv,sw,px,py,pz,pfx,pfy,pfz
  REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
    px0,py0,px1,py1,px2,py2,pr,f
  REAL_euwp,DIMENSION(:, :, :), ALLOCATABLE :: &
  r,qr,ar 
  REAL_euwp,DIMENSION(:, :, :), ALLOCATABLE :: &
  rax_ltmp,axax_ltmp,rlsqr_ltmp,axar_ltmp
  REAL_euwp,DIMENSION(:,:,:,:), ALLOCATABLE :: &
   x,ax
!General scratch storage
  REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
  scratch_zeros,scratch_ones

   CONTAINS
   SUBROUTINE allocate_scratch_datafields
   USE parameters, ONLY: np,mp,lp,ih,msz,mucnt,lrd
   USE, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_signaling_NAN
   IMPLICIT NONE
   REAL(KIND=euwp) real_initial_value
   INTEGER ierr,ierrtot
     ierrtot=0
   ALLOCATE ( tmpysnd1((np+1)*ih*(lp+1),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpysnd2((np+1)*ih*(lp+1),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpyrcv1((np+1)*ih*(lp+1),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpyrcv2((np+1)*ih*(lp+1),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   ALLOCATE ( tmpxsnd1((mp+1+2*ih)*ih*(lp+1),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpxsnd2((mp+1+2*ih)*ih*(lp+1),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpxrcv1((mp+1+2*ih)*ih*(lp+1),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpxrcv2((mp+1+2*ih)*ih*(lp+1),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   ALLOCATE ( tmpzsnd1((np+1+2*ih)*ih*(mp+1+2*ih),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpzsnd2((np+1+2*ih)*ih*(mp+1+2*ih),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpzrcv1((np+1+2*ih)*ih*(mp+1+2*ih),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpzrcv2((np+1+2*ih)*ih*(mp+1+2*ih),0:msz),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   ALLOCATE ( tmpysnd1m(mucnt*(np+1)*ih*(lp+1)),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpysnd2m(mucnt*(np+1)*ih*(lp+1)),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpyrcv1m(mucnt*(np+1)*ih*(lp+1)),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpyrcv2m(mucnt*(np+1)*ih*(lp+1)),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   ALLOCATE ( tmpxsnd1m(mucnt*(mp+1+2*ih)*ih*(lp+1)),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpxsnd2m(mucnt*(mp+1+2*ih)*ih*(lp+1)),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpxrcv1m(mucnt*(mp+1+2*ih)*ih*(lp+1)),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpxrcv2m(mucnt*(mp+1+2*ih)*ih*(lp+1)),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   ALLOCATE ( tmpzsnd1m(mucnt*(np+1+2*ih)*ih*(mp+1+2*ih)),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpzsnd2m(mucnt*(np+1+2*ih)*ih*(mp+1+2*ih)),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpzrcv1m(mucnt*(np+1+2*ih)*ih*(mp+1+2*ih)),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( tmpzrcv2m(mucnt*(np+1+2*ih)*ih*(mp+1+2*ih)),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   ALLOCATE (     cp(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (     cn(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( xstore(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   xant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( rhoadv(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (rhoadv2(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (rhoadv3(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (    rhr(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   rhr2(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   rhr3(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   uadv(1-ih:np+1+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   vadv(1-ih:np+ih, 1-ih:mp+1+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   wadv(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (  wadv2(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   hi(np,mp,lp),STAT=ierr); ierrtot=ierrtot+ierr
   ALLOCATE (    xtracer(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
   ALLOCATE (xtracer_out(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
   ALLOCATE (xtest(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
   ALLOCATE (    h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
   ALLOCATE (  bcx(mp,lp, 2),STAT=ierr) ;ierrtot=ierrtot+ierr
   ALLOCATE (  bcy(np,lp, 2),STAT=ierr) ;ierrtot=ierrtot+ierr


   ALLOCATE (  thant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   uant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   vant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   want(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( thtant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (pstrant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (  qvant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (  qcant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (  qrant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (  qsant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (  qiant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (  qgant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   ALLOCATE (   v1(1-ih:np+ih+1, 1-ih:mp+ih  , 1-ih:lp+ih  ),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   f1(1-ih:np+ih+1, 1-ih:mp+ih  , 1-ih:lp+ih  ),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( f1ns(1-ih:np+ih+1, 1-ih:mp+ih  , 1-ih:lp+ih  ),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (u1_out(1-ih:np+ih+1, 1-ih:mp+ih  , 1-ih:lp+ih  ),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   v2(1-ih:np+ih  , 1-ih:mp+ih+1, 1-ih:lp+ih  ),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   f2(1-ih:np+ih  , 1-ih:mp+ih+1, 1-ih:lp+ih  ),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( f2ns(1-ih:np+ih  , 1-ih:mp+ih+1, 1-ih:lp+ih  ),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (u2_out(1-ih:np+ih  , 1-ih:mp+ih+1, 1-ih:lp+ih  ),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   v3(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih+1),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (   f3(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih+1),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( f3ns(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih+1),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (u3_out(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih+1),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   ALLOCATE (  mx(np,mp,lp),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (  mn(np,mp,lp),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   ALLOCATE (u3gnd_store(1-ih:np+ih  , 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr; 
   ALLOCATE (u3sky_store(1-ih:np+ih  , 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr; 

   ALLOCATE (esxb(1-ih:np+ih  , 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr; 
   ALLOCATE (esyb(1-ih:np+ih  , 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr; 
   ALLOCATE (dsxb(1-ih:np+ih  , 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr; 
   ALLOCATE (dsyb(1-ih:np+ih  , 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr; 
   ALLOCATE (strxx(1-ih:np+ih  , 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr; 
   ALLOCATE (strxy(1-ih:np+ih  , 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr; 
   ALLOCATE (stryx(1-ih:np+ih  , 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr; 
   ALLOCATE (stryy(1-ih:np+ih  , 1-ih:mp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr; 

   ALLOCATE (csxb(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (csyb(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (cszb(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (strzx(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (strzy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE (strzz(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
   ALLOCATE ( scr1(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;

  ALLOCATE ( pexc(np,mp,lp),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE ( hlmc(np,mp,lp),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE (cf_etainv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE (fimy(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(su(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(sv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(sw(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(px(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(py(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(pz(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(pfx(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(pfy(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(pfz(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(px0(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(py0(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(px1(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(py1(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(px2(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(py2(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(pr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(f(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(r(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(qr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(ar(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE( rax_ltmp(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(axax_ltmp(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(axar_ltmp(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(rlsqr_ltmp(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE( x(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,lrd),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(ax(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,lrd),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(scratch_ones(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
  ALLOCATE(scratch_zeros(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
IF(ierrtot/=0) print *,'scratch allocate status: total number of errors is: ',ierrtot
!Initialize
#ifdef INITIALIZEVARIABLES
  real_initial_value = IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
     tmpysnd1  = real_initial_value
     tmpysnd2  = real_initial_value
     tmpyrcv1  = real_initial_value
     tmpyrcv2  = real_initial_value
     tmpxsnd1  = real_initial_value
     tmpxsnd2  = real_initial_value
     tmpxrcv1  = real_initial_value
     tmpxrcv2  = real_initial_value
     tmpzsnd1  = real_initial_value
     tmpzsnd2  = real_initial_value
     tmpzrcv1  = real_initial_value
     tmpzrcv2  = real_initial_value
     tmpysnd1m = real_initial_value
     tmpysnd2m = real_initial_value
     tmpyrcv1m = real_initial_value
     tmpyrcv2m = real_initial_value
     tmpxsnd1m = real_initial_value
     tmpxsnd2m = real_initial_value
     tmpxrcv1m = real_initial_value
     tmpxrcv2m = real_initial_value
     tmpzsnd1m = real_initial_value
     tmpzsnd2m = real_initial_value
     tmpzrcv1m = real_initial_value
     tmpzrcv2m = real_initial_value
     hi          = real_initial_value
     uadv        = real_initial_value
     vadv        = real_initial_value
     wadv        = real_initial_value
     xtracer     = real_initial_value
     xtracer_out = real_initial_value
     xtest       = real_initial_value
     h           = real_initial_value
     bcx         = real_initial_value
     bcy         = real_initial_value
      cp       = real_initial_value
      cn       = real_initial_value
     uadv      = real_initial_value
     vadv      = real_initial_value
     wadv      = real_initial_value
     wadv2     = real_initial_value
      xstore   = real_initial_value
      xant     = real_initial_value
     thant     = real_initial_value
      uant     = real_initial_value
      vant     = real_initial_value
      want     = real_initial_value
    thtant     = real_initial_value
     qvant     = real_initial_value
     qcant     = real_initial_value
     qrant     = real_initial_value
      v1       = real_initial_value
      f1       = real_initial_value
      f1ns     = real_initial_value
      u1_out   = real_initial_value
      v2       = real_initial_value
      f2       = real_initial_value
      f2ns     = real_initial_value
      u2_out   = real_initial_value
      v3       = real_initial_value
      f3       = real_initial_value
      f3ns     = real_initial_value
      u3_out   = real_initial_value
     rhr       = real_initial_value
     rhr2      = real_initial_value
     rhr3      = real_initial_value
  rhoadv       = real_initial_value
  rhoadv2      = real_initial_value
  rhoadv3      = real_initial_value
      mx       = real_initial_value
      mn       = real_initial_value
   u3gnd_store = real_initial_value
   u3sky_store = real_initial_value
    csxb       = real_initial_value
    csyb       = real_initial_value
    cszb       = real_initial_value 
    strzx      = real_initial_value
    strzy      = real_initial_value
    strzz      = real_initial_value
    scr1       = real_initial_value
    x          = real_initial_value
   ax          = real_initial_value
    r          = real_initial_value
   qr          = real_initial_value
   ar          = real_initial_value
  rax_ltmp     = real_initial_value
 axax_ltmp     = real_initial_value
rlsqr_ltmp     = real_initial_value
 axar_ltmp     = real_initial_value
 cf_etainv     = real_initial_value
 fimy          = real_initial_value
   su          = real_initial_value
   sv          = real_initial_value
   sw          = real_initial_value
  pfx          = real_initial_value
  pfy          = real_initial_value
  pfz          = real_initial_value
  px0          = real_initial_value
  py0          = real_initial_value
  px1          = real_initial_value
  py1          = real_initial_value
  px2          = real_initial_value
  py2          = real_initial_value
  pr           = real_initial_value
  f            = real_initial_value
   px          = real_initial_value
   py          = real_initial_value
   pz          = real_initial_value
     pexc      = real_initial_value
     hlmc      = real_initial_value
#endif
   scratch_ones=1._euwp
   scratch_zeros=0._euwp
   END SUBROUTINE allocate_scratch_datafields
   SUBROUTINE deallocate_scratch_datafields
     INTEGER ierr,ierrtot
     ierrtot=0
   DEALLOCATE ( tmpysnd1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpysnd2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpyrcv1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpyrcv2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpxsnd1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpxsnd2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpxrcv1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpxrcv2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpzsnd1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpzsnd2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpzrcv1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpzrcv2,STAT=ierr ) ;ierrtot=ierrtot+ierr;

   DEALLOCATE ( tmpysnd1m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpysnd2m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpyrcv1m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpyrcv2m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpxsnd1m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpxsnd2m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpxrcv1m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpxrcv2m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpzsnd1m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpzsnd2m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpzrcv1m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( tmpzrcv2m,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   hi,STAT=ierr ) ;ierrtot=ierrtot+ierr
   DEALLOCATE ( uadv,STAT=ierr ) ;ierrtot=ierrtot+ierr
   DEALLOCATE ( vadv,STAT=ierr ) ;ierrtot=ierrtot+ierr
   DEALLOCATE ( wadv,STAT=ierr ) ;ierrtot=ierrtot+ierr
   DEALLOCATE (xtracer,STAT=ierr ) ;ierrtot=ierrtot+ierr
   DEALLOCATE (xtracer_out,STAT=ierr ) ;ierrtot=ierrtot+ierr
   DEALLOCATE (xtest,STAT=ierr ) ;ierrtot=ierrtot+ierr
   DEALLOCATE (    h,STAT=ierr ) ;ierrtot=ierrtot+ierr
   DEALLOCATE (  bcx,STAT=ierr ) ;ierrtot=ierrtot+ierr
   DEALLOCATE (  bcy,STAT=ierr ) ;ierrtot=ierrtot+ierr


   DEALLOCATE (     cp,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (     cn,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( xstore,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   xant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (rhoadv ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (rhoadv2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (rhoadv3,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   rhr ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   rhr2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   rhr3,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   uadv,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   vadv,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   wadv,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (  wadv2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (  thant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   uant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   vant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   want,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( thtant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (pstrant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (  qvant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (  qcant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (  qrant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (  qsant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (  qiant,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (  qgant,STAT=ierr ) ;ierrtot=ierrtot+ierr;

   DEALLOCATE (   v1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   f1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( f1ns,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (u1_out,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   v2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   f2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( f2ns,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (u2_out,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   v3,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   f3,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( f3ns,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (u3_out,STAT=ierr ) ;ierrtot=ierrtot+ierr;

   DEALLOCATE (  mx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (  mn,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   u3gnd_store,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (   u3sky_store,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (esxb,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (esyb,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (dsxb,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (dsyb,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (strxx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (strxy,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (stryx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (stryy,STAT=ierr ) ;ierrtot=ierrtot+ierr;

   DEALLOCATE (csxb,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (csyb,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (cszb,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (strzx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (strzy,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE (strzz,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( scr1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( pexc ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
   DEALLOCATE ( hlmc ,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(cf_etainv,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(fimy,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(su,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(sv,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(sw,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(px,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(py,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(pz,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(pfx,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(pfy,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(pfz,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(px0,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(py0,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(px1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(py1,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(px2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(py2,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(pr,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(f,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(r,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(qr,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(ar,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(  rax_ltmp,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE( axax_ltmp,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(rlsqr_ltmp,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE( axar_ltmp,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE( x,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(ax,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(scratch_ones,STAT=ierr ) ;ierrtot=ierrtot+ierr;
  DEALLOCATE(scratch_zeros,STAT=ierr ) ;ierrtot=ierrtot+ierr;
IF(ierrtot/=0) print *,'scratch deallocate status: total number of errors is: ',ierrtot
   END SUBROUTINE deallocate_scratch_datafields
END MODULE scratch_datafields
#endif /*STATICMEM*/

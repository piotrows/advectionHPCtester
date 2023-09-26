#include  "../advection/src_algorithms/renames.inc"
#if (STATICMEM == 1)
MODULE eulag_datafields
   USE precisions
   USE mod_parameters, ONLY: n,m,l,np,mp,lp,ih
   IMPLICIT NONE
!Advection dwarf testing
   REAL(KIND=euwp),DIMENSION(np,mp,lp) :: hi
   REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
      h,xtest,xtracer, xtracer_out
   REAL(KIND=euwp),DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+ih  ) :: uadv
   REAL(KIND=euwp),DIMENSION(1-ih:np+ih  , 1-ih:mp+1+ih, 1-ih:lp+ih  ) :: vadv
   REAL(KIND=euwp),DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+1+ih) :: wadv
   REAL(KIND=euwp)  bcx(mp,lp, 2),bcy(np,lp, 2)
!Solver dwarf variables
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
   p,fd,pstr,pext
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
   th, tht
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
   ft,fx,fy,fz
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
   ftht_expl,ft_expl, fx_expl, fy_expl, fz_expl, fpext_expl

  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:1) :: &
  rho
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2) :: &
  u,v,w
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2) :: &
   ox,oy,oz
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
   tte,ppe,the,pexe,rhe
  REAL(KIND=euwp), DIMENSION(np,mp,lp)   ::  rhoi

  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
  uref,vref,thref,pextref

   REAL(KIND=euwp)  rh_a_lr(mp,lp, 2), rh_a_bt(np,lp, 2)
   REAL(KIND=euwp)  th_a_lr(mp,lp, 2), th_a_bt(np,lp, 2)
   REAL(KIND=euwp)   u_a_lr(mp,lp, 2),  u_a_bt(np,lp, 2)
   REAL(KIND=euwp)   v_a_lr(mp,lp, 2),  v_a_bt(np,lp, 2)
   REAL(KIND=euwp)   w_a_lr(mp,lp, 2),  w_a_bt(np,lp, 2)
   REAL(KIND=euwp) tht_a_lr(mp,lp, 2),tht_a_bt(np,lp, 2)
   REAL(KIND=euwp)  ex_a_lr(mp,lp, 2), ex_a_bt(np,lp, 2)
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


END MODULE eulag_datafields
#else /*STATICMEM*/ 
MODULE eulag_datafields
   USE precisions
   IMPLICIT NONE 
#include "../advection/src_algorithms/defines.inc"
  REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
   th, tht, &
   ft,fx,fy,fz, &
   ft_expl, ftht_expl, fx_expl, fy_expl, fz_expl, fpext_expl

  REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
    u_a_lr,v_a_lr,w_a_lr,th_a_lr,tht_a_lr,rh_a_lr,ex_a_lr, &
    u_a_bt,v_a_bt,w_a_bt,th_a_bt,tht_a_bt,rh_a_bt,ex_a_bt

  REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
   p,fd,pstr
  REAL_euwp,DIMENSION(:,:,:,:), ALLOCATABLE :: &
  rho,u,v,w

!Geometry: topological and metric information
  REAL_euwp,DIMENSION(:,:,:,:), ALLOCATABLE :: &
   ox,oy,oz
  REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
   tte,the,pexe,ppe,pext,rhe
  REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE   ::  rhoi
  REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE   ::  &
    uref,vref,thref,pextref 
  


   CONTAINS

   SUBROUTINE allocate_eulag_datafields
   USE mod_parameters, ONLY: np,mp,lp,ih
   USE, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_signaling_NAN
   REAL(KIND=euwp) real_initial_value
   INTEGER ierr,ierrtot
     ierr=0
     ierrtot=0
     ALLOCATE (  rho(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:1)  ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   th(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  tht(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    u(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)  ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    v(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)  ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    w(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)  ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   ox(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)  ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   oy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)  ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   oz(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)  ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    p(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( pstr(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( pext(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  tte(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  the(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( pexe(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  ppe(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  rhe(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
!    ALLOCATE (   fd(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( rhoi(np,mp,lp)                                ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fx(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fz(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   ft(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    fx_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    fy_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    fz_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    ft_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  ftht_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( fpext_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr

     IF(ierrtot/=0) print *,'Allocation eulag datafields status: total number of errors is: ',ierrtot
     IF(ierrtot/=0) STOP 'Lack of memory for eulag datafields' 
!Initialize
#ifdef INITIALIZEVARIABLES
  real_initial_value = IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
! real_initial_value = 0._euwp !IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
           rho(:,:,:,:) = real_initial_value
            th(:,:,:)   = real_initial_value
           tht(:,:,:)   = real_initial_value
             u(:,:,:,:) = real_initial_value
             v(:,:,:,:) = real_initial_value
             w(:,:,:,:) = real_initial_value
            ox(:,:,:,:) = real_initial_value
            oy(:,:,:,:) = real_initial_value
            oz(:,:,:,:) = real_initial_value
             p(:,:,:)   = real_initial_value
          pstr(:,:,:)   = real_initial_value
          pext(:,:,:)   = real_initial_value
           tte(:,:,:)   = real_initial_value
           the(:,:,:)   = real_initial_value
          pexe(:,:,:)   = real_initial_value
           ppe(:,:,:)   = real_initial_value
           rhe(:,:,:)   = real_initial_value
!           fd(:,:,:)   = real_initial_value
          rhoi(:,:,:)   = real_initial_value
            fx(:,:,:)   = real_initial_value
            fy(:,:,:)   = real_initial_value
            fz(:,:,:)   = real_initial_value
            ft(:,:,:)   = real_initial_value
       fx_expl(:,:,:)   = real_initial_value
       fy_expl(:,:,:)   = real_initial_value
       fz_expl(:,:,:)   = real_initial_value
       ft_expl(:,:,:)   = real_initial_value
       fpext_expl(:,:,:)= real_initial_value
#endif
   END SUBROUTINE allocate_eulag_datafields

   SUBROUTINE allocate_reference_datafields
   USE mod_parameters, ONLY: np,mp,lp,ih
   USE, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_signaling_NAN
   REAL(KIND=euwp) real_initial_value
   INTEGER ierr,ierrtot
     ierr=0
     ierrtot=0
!Solver dwarf variables
     ALLOCATE (  thref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  uref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  vref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  pextref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr

     IF(ierrtot/=0) print *,'Allocation of reference datafields status: total number of errors is: ',ierrtot
     IF(ierrtot/=0) STOP 'Lack of memory for reference datafields' 
!Initialize
#ifdef INITIALIZEVARIABLES
  real_initial_value = IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
! real_initial_value = 0._euwp !IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
       uref      = real_initial_value
       vref      = real_initial_value
       thref     = real_initial_value
       pextref   = real_initial_value
#endif
   END SUBROUTINE allocate_reference_datafields

   SUBROUTINE allocate_advbc_datafields
   USE mod_parameters, ONLY: np,mp,lp,ih
   USE, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_signaling_NAN
   REAL(KIND=euwp) real_initial_value
   INTEGER ierr,ierrtot
     ierr=0
     ierrtot=0

     ALLOCATE (   u_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   v_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   w_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  th_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( tht_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  rh_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  ex_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   u_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   v_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   w_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  th_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( tht_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  rh_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  ex_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr

     IF(ierrtot/=0) print *,'Allocation advbc status: total number of errors is: ',ierrtot
     IF(ierrtot/=0) STOP 'Lack of memory for advbc variables' 
!Initialize
#ifdef INITIALIZEVARIABLES
  real_initial_value = IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
! real_initial_value = 0._euwp !IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
        th_a_lr(:,:,:)   = real_initial_value
        th_a_bt(:,:,:)   = real_initial_value
       tht_a_lr(:,:,:)   = real_initial_value
       tht_a_bt(:,:,:)   = real_initial_value
         u_a_lr(:,:,:)   = real_initial_value
         u_a_bt(:,:,:)   = real_initial_value
         v_a_lr(:,:,:)   = real_initial_value
         v_a_bt(:,:,:)   = real_initial_value
         w_a_lr(:,:,:)   = real_initial_value
         w_a_bt(:,:,:)   = real_initial_value
        rh_a_lr(:,:,:)   = real_initial_value
        rh_a_bt(:,:,:)   = real_initial_value
        ex_a_lr(:,:,:)   = real_initial_value
        ex_a_bt(:,:,:)   = real_initial_value
#endif
   END SUBROUTINE allocate_advbc_datafields

   SUBROUTINE deallocate_eulag_datafields
     INTEGER ierrtot,ierr
     ierrtot=0
     IF(ALLOCATED(p       ))  DEALLOCATE (    p      ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(pstr    ))  DEALLOCATE (  pstr     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(pext    ))  DEALLOCATE (  pext     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
!     IF(ALLOCATED(fd)) DEALLOCATE (   fd ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(rho     ))  DEALLOCATE (   rho     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(u       ))  DEALLOCATE (     u     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(v       ))  DEALLOCATE (     v     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(w       ))  DEALLOCATE (     w     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(ox      ))  DEALLOCATE (    ox     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(oy      ))  DEALLOCATE (    oy     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(oz      ))  DEALLOCATE (    oz     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(tte     ))  DEALLOCATE (   tte     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(the     ))  DEALLOCATE (   the     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(pexe    ))  DEALLOCATE (   pexe    ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(ppe     ))  DEALLOCATE (   ppe     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(rhe     ))  DEALLOCATE (   rhe     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(rhoi    ))  DEALLOCATE (  rhoi     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(thref   ))  DEALLOCATE ( thref     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(uref    ))  DEALLOCATE (  uref     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(vref    ))  DEALLOCATE (  vref     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(pextref ))  DEALLOCATE ( pextref   ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(th      ))  DEALLOCATE (    th     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(tht     ))  DEALLOCATE (    tht    ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(fx      ))  DEALLOCATE (    fx     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(fy      ))  DEALLOCATE (    fy     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(fz      ))  DEALLOCATE (    fz     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(ft      ))  DEALLOCATE (    ft     ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(fx_expl ))  DEALLOCATE (    fx_expl,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(fy_expl ))  DEALLOCATE (    fy_expl,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(fz_expl ))  DEALLOCATE (    fz_expl,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(ft_expl ))  DEALLOCATE (    ft_expl,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED( ftht_expl)) DEALLOCATE ( ftht_expl,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(fpext_expl)) DEALLOCATE (fpext_expl,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(u_a_lr  ))   DEALLOCATE (    u_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(v_a_lr  ))   DEALLOCATE (    v_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(w_a_lr  ))   DEALLOCATE (    w_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(th_a_lr ))   DEALLOCATE (   th_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(tht_a_lr))   DEALLOCATE (  tht_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(rh_a_lr ))   DEALLOCATE (   rh_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(ex_a_lr ))   DEALLOCATE (   ex_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(u_a_bt  ))   DEALLOCATE (    u_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(v_a_bt  ))   DEALLOCATE (    v_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(w_a_bt  ))   DEALLOCATE (    w_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(th_a_bt ))   DEALLOCATE (   th_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(tht_a_bt))   DEALLOCATE (  tht_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(rh_a_bt ))   DEALLOCATE (   rh_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr
     IF(ALLOCATED(ex_a_bt ))   DEALLOCATE (   ex_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr

     IF(ierrtot/=0) print *,'EULAG datafields deallocate status: total number of errors is: ',ierrtot
   END SUBROUTINE deallocate_eulag_datafields

END MODULE eulag_datafields
#endif /*STATICMEM*/
#if(1==0)
   SUBROUTINE allocate_eulag_all_datafields
   USE mod_parameters, ONLY: np,mp,lp,ih
   USE, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_signaling_NAN
   REAL(KIND=euwp) real_initial_value
   INTEGER ierr,ierrtot
     ierr=0
     ierrtot=0
     ALLOCATE (   hi(np,mp,lp),STAT=ierr); ierrtot=ierrtot+ierr
     ALLOCATE ( uadv(1-ih:np+1+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
     ALLOCATE ( vadv(1-ih:np+ih, 1-ih:mp+1+ih, 1-ih:lp+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
     ALLOCATE ( wadv(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
     ALLOCATE (    xtracer(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
     ALLOCATE (xtracer_out(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
     ALLOCATE (xtest(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
     ALLOCATE (    h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr) ;ierrtot=ierrtot+ierr
     ALLOCATE (  bcx(mp,lp, 2),STAT=ierr) ;ierrtot=ierrtot+ierr
     ALLOCATE (  bcy(np,lp, 2),STAT=ierr) ;ierrtot=ierrtot+ierr
!Solver dwarf variables
     ALLOCATE (    p(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( pstr(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( pext(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fd(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  rho(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:1),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    u(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    v(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (    w(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   ox(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   oy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   oz(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   tte(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   the(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   pexe(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   ppe(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   rhe(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( rhoi(np,mp,lp)                              ,STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  thref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  uref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  vref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  pextref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
!Timestep dwaf variables
     ALLOCATE (   th(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (  tht(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fx(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fz(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   ft(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fx_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fy_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fz_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   ft_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE ( ftht_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
     ALLOCATE (   fpext_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (   u_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (   v_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (   w_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (  th_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( tht_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (  rh_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (  ex_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (   u_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (   v_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (   w_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (  th_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE ( tht_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (  rh_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr
    ALLOCATE (  ex_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr

     IF(ierrtot/=0) print *,'initialize status: total number of errors is: ',ierrtot
!Initialize
#ifdef INITIALIZEVARIABLES
  real_initial_value = IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
! real_initial_value = 0._euwp !IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
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
       p         = real_initial_value
     pext        = real_initial_value
       fd        = real_initial_value
      rho        = real_initial_value
        u        = real_initial_value
        v        = real_initial_value
        w        = real_initial_value
       ox        = real_initial_value
       oy        = real_initial_value
       oz        = real_initial_value
       tte       = real_initial_value
       the       = real_initial_value
       pexe      = real_initial_value
       ppe       = real_initial_value
       rhe       = real_initial_value
       pstr      = real_initial_value
     rhoi        = real_initial_value
       uref      = real_initial_value
       vref      = real_initial_value
       thref     = real_initial_value
       pextref   = real_initial_value
       th        = real_initial_value
      tht        = real_initial_value
       fx        = real_initial_value
       fy        = real_initial_value
       fz        = real_initial_value
       ft        = real_initial_value
       fx_expl   = real_initial_value
       fy_expl   = real_initial_value
       fz_expl   = real_initial_value
       ft_expl   = real_initial_value
       fpext_expl= real_initial_value
       th_a_lr   = real_initial_value
       th_a_bt   = real_initial_value
       tht_a_lr  = real_initial_value
       tht_a_bt  = real_initial_value
       u_a_lr    = real_initial_value
       u_a_bt    = real_initial_value
       v_a_lr    = real_initial_value
       v_a_bt    = real_initial_value
       w_a_lr    = real_initial_value
       w_a_bt    = real_initial_value
       rh_a_lr   = real_initial_value
       rh_a_bt   = real_initial_value
       ex_a_lr   = real_initial_value
       ex_a_bt   = real_initial_value
#endif
   END SUBROUTINE allocate_eulag_all_datafields
   SUBROUTINE allocate_geometry_all_datafields
   USE mod_parameters, ONLY: n,m,l,np,mp,lp,ih
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
    ALLOCATE (     sfflx(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnyx(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnyy(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnyz(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     sffly(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnzx(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnzy(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (      tnzz(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     sfflz(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;

    ALLOCATE (     tnxxV(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnxyV(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnxzV(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnyxV(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnyyV(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnyzV(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnzxV(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnzyV(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (     tnzzV(1:np,1:mp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;

     IF(ierrtot/=0) print *,'Allocation of geometry datafields status: total number of errors is: ',ierrtot
     IF(ierrtot/=0) STOP 'Lack of memory for geometry datafields' 

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
   END SUBROUTINE allocate_geometry_all_datafields
#endif /*1==0*/

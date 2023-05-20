MODULE advection_CEtimestep_twostep_GPUBC
   USE precisions
   USE mpi_parallel, ONLY: update,update_multi,update3,update3_multi
   USE mpi_parallel, ONLY: ttbeg,ttend,iup
   USE mpi_parallel, ONLY: leftedge,rightedge,botedge,topedge,gndedge,skyedge
   USE    module_antidiff3d_gauge_gpubc, ONLY: antidiff3d_gauge_gpubc
   USE module_antidiff3d_standard_gpubc, ONLY: antidiff3d_standard_gpubc
   USE            module_upwind3d_gpubc, ONLY:   upwind3d_gpubc 
#ifdef CUDACODE
   USE cudafor
   USE mpi_parallel, ONLY: istream1,istream2,istream3 
#endif
#ifdef PNETCDF
   USE mpi_parallel, ONLY: pnet_out_chunk 
#endif
#ifdef DIAGNOS 
!  USE eulag_diagnostics, ONLY: diagnos_actual
#endif

   USE scratch_datafields, ONLY: rhr,rhoadv
   IMPLICIT NONE
   INTEGER(KIND=iintegers) icnt,iprint
   CONTAINS
#include "../src_algorithms/defines.inc"
#include "../src_algorithms/mpdataoperators.inc"

   SUBROUTINE CE_advection(lupdatemulti,mpdata_alg_choice,liner, &
                            ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,  & 
                           do_serialization_in,do_serialization_out, &
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt,    &
                           qr,  qr_lr,  qr_bt,    &
                           qi,  qi_lr,  qi_bt,    &
                           qs,  qs_lr,  qs_bt,    &
                           qg,  qg_lr,  qg_bt) 

   USE scratch_datafields, ONLY: thant,uant,vant,want,thtant,pstrant
   USE scratch_datafields, ONLY: qvant,qcant,qrant,qiant,qsant,qgant
   USE mpi_parallel, ONLY: gndedge,skyedge

   INTEGER(KIND=iintegers),INTENT(IN) :: ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih
   INTEGER(KIND=iintegers),INTENT(IN) :: liner,mpdata_alg_choice
   INTEGER(KIND=iintegers),PARAMETER :: iflip1=1
   INTEGER(KIND=iintegers),PARAMETER :: iflip0=0
   INTEGER(KIND=iintegers) :: istat

   LOGICAL,INTENT(IN) :: lupdatemulti
   LOGICAL,OPTIONAL :: do_serialization_in
   LOGICAL,OPTIONAL :: do_serialization_out




   REAL_euwp, INTENT(IN) ::                               & 
                      rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),    & 
                      rhoadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(INOUT) ::                            & 
                    uadv(1-ih:np+1+ih,1-ih:mp+ih,1-ih:lp+ih),   & 
                    vadv(1-ih:np+ih,1-ih:mp+1+ih,1-ih:lp+ih),   &
                    wadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+1+ih),   &
                      th(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       u(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                       v(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       w(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                     tht(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    pstr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(INOUT), OPTIONAL  ::               &
                      qv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qc(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &   
                      qr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                      qs(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qi(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qg(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)   

   REAL_euwp, INTENT(IN) ::    &                
                       th_lr(mp,lp,2), &
                        u_lr(mp,lp,2), &
                        v_lr(mp,lp,2), &
                        w_lr(mp,lp,2), &
                      tht_lr(mp,lp,2), &
                     pstr_lr(mp,lp,2)

   REAL_euwp, INTENT(IN), OPTIONAL ::    & 
                       qv_lr(mp,lp,2), &
                       qc_lr(mp,lp,2), &
                       qr_lr(mp,lp,2), &
                       qs_lr(mp,lp,2), &
                       qi_lr(mp,lp,2), &
                       qg_lr(mp,lp,2)

   REAL_euwp, INTENT(IN) ::    &                
                       th_bt(np,lp,2), &
                        u_bt(np,lp,2), &
                        v_bt(np,lp,2), &
                        w_bt(np,lp,2), &
                      tht_bt(np,lp,2), &
                     pstr_bt(np,lp,2)   
   
   REAL_euwp, INTENT(INOUT), OPTIONAL ::    &            
                       qv_bt(np,lp,2), &
                       qc_bt(np,lp,2), &
                       qr_bt(np,lp,2), &
                       qs_bt(np,lp,2), &
                       qi_bt(np,lp,2), &
                       qg_bt(np,lp,2)
   REAL(KIND=euwp)    :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp 
   REAL(KIND=euwp)    :: hinv 
   INTEGER(KIND=iintegers) :: i,j,k
   PROCEDURE(antidiff3d_gauge_gpubc),POINTER :: mpdatm3d_antidiff 
   PROCEDURE(antidiff3d_gauge_gpubc),POINTER :: mpdata3d_antidiff 

 
#ifdef SERIALIZE
        IF(do_serialization_in) THEN !0th iteration
                 print *,'Serializing in fields...\n'

            !$ser savepoint mysnapshot.run-in
            !$ser mode write
            !$ser data uadv,vadv,wadv,rhoadv
            !$ser data u=u u_lr=u_lr u_bt=u_bt
            !$ser data v=v v_lr=v_lr v_bt=v_bt
            !$ser data w=w w_lr=w_lr w_bt=w_bt 
            !$ser data th=th th_lr=th_lr th_bt=th_bt
            !$ser data tht=tht tht_lr=tht_lr tht_bt=tht_bt 
            !$ser data pstr=pstr pstr_lr=pstr_lr pstr_bt=pstr_bt  
            !$ser data qv=qv qv_lr=qv_lr qv_bt=qv_bt 
            !$ser data qc=qc qc_lr=qc_lr qc_bt=qc_bt
            !$ser data qr=qr qr_lr=qr_lr qr_bt=qr_bt
            !$ser data qi=qi qi_lr=qi_lr qi_bt=qi_bt 
            !$ser data qs=qs qs_lr=qs_lr qs_bt=qs_bt
            !$ser data qg=qg qg_lr=qg_lr qg_bt=qg_bt
        ENDIF
#endif /*SERIALIZE*/
  IF(PRESENT(do_serialization_in).OR.PRESENT(do_serialization_out)) THEN
! Do nothing, just to remove compiler remark
  ENDIF

   CALL ttbeg(82) 
   CALL ttbeg(18) 
   CALL ttbeg(1018) 
   IF(mpdata_alg_choice.eq.0.OR.mpdata_alg_choice.gt.2) THEN
      mpdatm3d_antidiff =>antidiff3d_gauge_gpubc
      mpdata3d_antidiff =>antidiff3d_standard_gpubc
 ELSE IF(mpdata_alg_choice.eq.1) THEN
      mpdatm3d_antidiff =>antidiff3d_gauge_gpubc
      mpdata3d_antidiff =>antidiff3d_gauge_gpubc
 ELSE IF(mpdata_alg_choice.eq.2) THEN
      mpdatm3d_antidiff =>antidiff3d_standard_gpubc
      mpdata3d_antidiff =>antidiff3d_standard_gpubc
 ENDIF

   
   CALL putbcinxh(th  ,  th_lr,  th_bt,ibcx,ibcy,np,mp,lp,ih)
!   CALL zerogradbcinxh(th , np,mp,lp,ih)
   CALL putbcinxh(u   ,   u_lr,   u_bt,ibcx,ibcy,np,mp,lp,ih)
   CALL putbcinxh(v   ,   v_lr,   v_bt,ibcx,ibcy,np,mp,lp,ih)
   CALL putbcinxh(w   ,   w_lr,   w_bt,ibcx,ibcy,np,mp,lp,ih)
   CALL putbcinxh(tht , tht_lr, tht_bt,ibcx,ibcy,np,mp,lp,ih)
!   CALL zerogradbcinxh(tht , np,mp,lp,ih)
   CALL putbcinxh(pstr,pstr_lr,pstr_bt,ibcx,ibcy,np,mp,lp,ih)
!   CALL zerogradbcinxh(pstr , np,mp,lp,ih)
   IF (PRESENT(qv)) CALL putbcinxh(qv  ,  qv_lr,  qv_bt,ibcx,ibcy,np,mp,lp,ih)
   IF (PRESENT(qc)) CALL putbcinxh(qc  ,  qc_lr,  qc_bt,ibcx,ibcy,np,mp,lp,ih)
   IF (PRESENT(qr)) CALL putbcinxh(qr  ,  qr_lr,  qr_bt,ibcx,ibcy,np,mp,lp,ih)
   IF (PRESENT(qs)) CALL putbcinxh(qs  ,  qs_lr,  qs_bt,ibcx,ibcy,np,mp,lp,ih)
   IF (PRESENT(qi)) CALL putbcinxh(qi  ,  qi_lr,  qi_bt,ibcx,ibcy,np,mp,lp,ih)
   IF (PRESENT(qg)) CALL putbcinxh(qg  ,  qg_lr,  qg_bt,ibcx,ibcy,np,mp,lp,ih)
!   IF (PRESENT(qr)) CALL zerogradbcinxh(qr ,ibcx,ibcy, np,mp,lp,ih)
!   IF (PRESENT(qs)) CALL zerogradbcinxh(qs ,ibcx,ibcy, np,mp,lp,ih)
!   IF (PRESENT(qi)) CALL zerogradbcinxh(qi ,ibcx,ibcy, np,mp,lp,ih)
!   IF (PRESENT(qg)) CALL zerogradbcinxh(qg ,ibcx,ibcy, np,mp,lp,ih)
     CALL  update3(  th,np,mp,lp,np,mp,lp,iup,ih)
     CALL  update3(   u,np,mp,lp,np,mp,lp,iup,ih)
     CALL  update3(   v,np,mp,lp,np,mp,lp,iup,ih)
     CALL  update3(   w,np,mp,lp,np,mp,lp,iup,ih)
     CALL  update3( tht,np,mp,lp,np,mp,lp,iup,ih)
     CALL  update3(pstr,np,mp,lp,np,mp,lp,iup,ih)
     IF(PRESENT(qv)) CALL  update3(qv,np,mp,lp,np,mp,lp,iup,ih)
     IF(PRESENT(qc)) CALL  update3(qc,np,mp,lp,np,mp,lp,iup,ih)
     IF(PRESENT(qr)) CALL  update3(qr,np,mp,lp,np,mp,lp,iup,ih)
     IF(PRESENT(qi)) CALL  update3(qi,np,mp,lp,np,mp,lp,iup,ih)
     IF(PRESENT(qs)) CALL  update3(qs,np,mp,lp,np,mp,lp,iup,ih)
     IF(PRESENT(qg)) CALL  update3(qg,np,mp,lp,np,mp,lp,iup,ih)
   IF(PRESENT(qg)) THEN
     CALL upwind_firstep_GPUBC_qg( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt,    &
                           qr,  qr_lr,  qr_bt,    &
                           qi,  qi_lr,  qi_bt,    &
                           qs,  qs_lr,  qs_bt,    &
                           qg,  qg_lr,  qg_bt) 
   ELSE IF (PRESENT(qs)) THEN 
     CALL upwind_firstep_GPUBC_qs( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt,    &
                           qr,  qr_lr,  qr_bt,    &
                           qi,  qi_lr,  qi_bt,    &
                           qs,  qs_lr,  qs_bt)    
   ELSE IF (PRESENT(qi)) THEN 
     CALL upwind_firstep_GPUBC_qi( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt,    &
                           qr,  qr_lr,  qr_bt,    &
                           qi,  qi_lr,  qi_bt)     
   ELSE IF (PRESENT(qr)) THEN 
     CALL upwind_firstep_GPUBC_qr( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt,    &
                           qr,  qr_lr,  qr_bt)     
   ELSE IF (PRESENT(qc)) THEN 
     CALL upwind_firstep_GPUBC_qc( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt)     
   ELSE IF (PRESENT(qv)) THEN 
     CALL upwind_firstep_GPUBC_qv( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt)     
   ELSE !no presence of moist variables
     CALL upwind_firstep_GPUBC_dry( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt)     
   ENDIF
   CALL ttend(1018) 
!! @cuf istat=cudaDeviceSynchronize()
 !@cuf istat=cudaStreamSynchronize(istream3)
 !@cuf istat=cudaStreamSynchronize(istream2)
 !@cuf istat=cudaStreamSynchronize(istream1)

   IF(lupdatemulti) THEN 
     IF(PRESENT(qv)) THEN
      CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,thant,uant,vant,want,thtant,pstrant,qvant,qcant)
!      CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,thant,uant)
!      CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,vant,want)
!      CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,thtant,pstrant)
!      CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,qvant,qcant)
     ELSE
      CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,thant,uant,vant,want,thtant,pstrant)
     ENDIF
   ELSE
     CALL  update(  thant,np,mp,lp,np,mp,lp,iup,ih)
     CALL  update(   uant,np,mp,lp,np,mp,lp,iup,ih)
     CALL  update(   vant,np,mp,lp,np,mp,lp,iup,ih)
     CALL  update(   want,np,mp,lp,np,mp,lp,iup,ih)
     CALL  update( thtant,np,mp,lp,np,mp,lp,iup,ih)
     CALL  update(pstrant,np,mp,lp,np,mp,lp,iup,ih)
     IF(PRESENT(qv)) CALL  update(  qvant,np,mp,lp,np,mp,lp,iup,ih)
     IF(PRESENT(qc)) CALL  update(  qcant,np,mp,lp,np,mp,lp,iup,ih)
   ENDIF !lupdatemulti
  CALL ttend(18) 
#ifdef DIAGNOS 
!  CALL diagnos_actual(-1,u,v,w,thant,np,mp,lp,ih,thtant,pstr,rhoadv)
#endif
IF(liner == 0) THEN
   CALL mpdata3d_antidiff    (lupdatemulti,uadv,vadv,wadv,  thant,   th,rhr,rhoadv,iflip0,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
   CALL mpdatm3d_antidiff    (lupdatemulti,uadv,vadv,wadv,   uant,    u,rhr,rhoadv,iflip1,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
   CALL mpdatm3d_antidiff    (lupdatemulti,uadv,vadv,wadv,   vant,    v,rhr,rhoadv,iflip1,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
   CALL mpdatm3d_antidiff    (lupdatemulti,uadv,vadv,wadv,   want,    w,rhr,rhoadv,iflip0,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
   CALL mpdatm3d_antidiff    (lupdatemulti,uadv,vadv,wadv, thtant,  tht,rhr,rhoadv,iflip0,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
   CALL mpdatm3d_antidiff    (lupdatemulti,uadv,vadv,wadv,pstrant, pstr,rhr,rhoadv,iflip0,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
   IF(PRESENT(qv)) &
   CALL mpdata3d_antidiff    (lupdatemulti,uadv,vadv,wadv,  qvant,   qv,rhr,rhoadv,iflip0,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
   IF(PRESENT(qc)) &
   CALL mpdata3d_antidiff    (lupdatemulti,uadv,vadv,wadv,  qcant,   qc,rhr,rhoadv,iflip0,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)

 IF(PRESENT(qg)) THEN
   CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,qrant,qsant,qiant,qgant,qv,qc)
 !   CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,qrant,qsant)
 !   CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,qiant,qgant)
 !   CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,qv,qc)
 ELSE IF(PRESENT(qs)) THEN
    CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,qrant,qsant,qiant,qv,qc)
 ELSE IF(PRESENT(qr)) THEN
    CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,qrant,qv,qc)
 ELSE IF(PRESENT(qc)) THEN
    CALL  update_multi(np,mp,lp,np,mp,lp,iup,ih,qv,qc)
 ENDIF
 IF(PRESENT(qr)) &
 CALL mpdata3d_antidiff    (lupdatemulti,uadv,vadv,wadv,  qrant,   qr,rhr,rhoadv,iflip0,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
 IF(PRESENT(qs)) &
 CALL mpdata3d_antidiff    (lupdatemulti,uadv,vadv,wadv,  qsant,   qs,rhr,rhoadv,iflip0,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
 IF(PRESENT(qi)) &
 CALL mpdata3d_antidiff    (lupdatemulti,uadv,vadv,wadv,  qiant,   qi,rhr,rhoadv,iflip0,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
 IF(PRESENT(qg)) &
 CALL mpdata3d_antidiff    (lupdatemulti,uadv,vadv,wadv,  qgant,   qg,rhr,rhoadv,iflip0,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.) 
ELSE
   th(:,:,:)=  thant(:,:,:)
    u(:,:,:)=   uant(:,:,:)
    v(:,:,:)=   vant(:,:,:)
    w(:,:,:)=   want(:,:,:)
  tht(:,:,:)= thtant(:,:,:)
 pstr(:,:,:)=pstrant(:,:,:)
IF(PRESENT(qv)) qv(:,:,:)=qvant(:,:,:)
IF(PRESENT(qc)) qc(:,:,:)=qcant(:,:,:)
IF(PRESENT(qr)) qr(:,:,:)=qrant(:,:,:)
IF(PRESENT(qs)) qs(:,:,:)=qsant(:,:,:)
IF(PRESENT(qi)) qi(:,:,:)=qiant(:,:,:)
IF(PRESENT(qg)) qg(:,:,:)=qgant(:,:,:)
ENDIF
   CALL ttend(82)
#ifdef SERIALIZE
        IF(do_serialization_out) THEN !last iteration
                 print *,'Serializing out fields...\n'

            !$ser savepoint mysnapshot.run-out
            !$ser mode write
            !$ser data u1=uadv u2=vadv u3=wadv 
            !$ser data u=u u_lr=u_lr u_bt=u_bt
            !$ser data v=v v_lr=v_lr v_bt=v_bt 
            !$ser data w=w w_lr=w_lr w_bt=w_bt
            !$ser data th=th th_lr=th_lr th_bt=th_bt
            !$ser data tht=tht tht_lr=tht_lr tht_bt=tht_bt
            !$ser data pstr=pstr pstr_lr=pstr_lr pstr_bt=pstr_bt 
            !$ser data qv=qv qv_lr=qv_lr qv_bt=qv_bt 
            !$ser data qc=qc qc_lr=qc_lr qc_bt=qc_bt 
            !$ser data qr=qr qr_lr=qr_lr qr_bt=qr_bt 
            !$ser data qi=qi qi_lr=qi_lr qi_bt=qi_bt
            !$ser data qs=qs qs_lr=qs_lr qs_bt=qs_bt
            !$ser data qg=qg qg_lr=qg_lr qg_bt=qg_bt
        ENDIF
#endif /*SERIALIZE*/
   END SUBROUTINE CE_advection

   SUBROUTINE upwind_firstep_GPUBC_qg( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt,    &
                           qr,  qr_lr,  qr_bt,    &
                           qi,  qi_lr,  qi_bt,    &
                           qs,  qs_lr,  qs_bt,    &
                           qg,  qg_lr,  qg_bt) 

   USE scratch_datafields, ONLY: thant,uant,vant,want,thtant,pstrant
   USE scratch_datafields, ONLY: qvant,qcant,qrant,qiant,qsant,qgant
   USE mpi_parallel, ONLY: gndedge,skyedge

   INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
   INTEGER(KIND=iintegers) :: istat


   REAL_euwp, INTENT(IN) ::                               & 
                      rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),    & 
                      rhoadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::                            & 
                    uadv(1-ih:np+1+ih,1-ih:mp+ih,1-ih:lp+ih),   & 
                    vadv(1-ih:np+ih,1-ih:mp+1+ih,1-ih:lp+ih),   &
                    wadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+1+ih),   &
                      th(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       u(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                       v(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       w(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                     tht(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    pstr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN), OPTIONAL  ::               &
                      qv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qc(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &   
                      qr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                      qs(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qi(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qg(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)   

   REAL_euwp, INTENT(IN) ::    &                
                       th_lr(mp,lp,2), &
                        u_lr(mp,lp,2), &
                        v_lr(mp,lp,2), &
                        w_lr(mp,lp,2), &
                      tht_lr(mp,lp,2), &
                     pstr_lr(mp,lp,2)

   REAL_euwp, INTENT(IN), OPTIONAL ::    & 
                       qv_lr(mp,lp,2), &
                       qc_lr(mp,lp,2), &
                       qr_lr(mp,lp,2), &
                       qs_lr(mp,lp,2), &
                       qi_lr(mp,lp,2), &
                       qg_lr(mp,lp,2)

   REAL_euwp, INTENT(IN) ::    &                
                       th_bt(np,lp,2), &
                        u_bt(np,lp,2), &
                        v_bt(np,lp,2), &
                        w_bt(np,lp,2), &
                      tht_bt(np,lp,2), &
                     pstr_bt(np,lp,2)   
   
   REAL_euwp, INTENT(IN), OPTIONAL ::    &            
                       qv_bt(np,lp,2), &
                       qc_bt(np,lp,2), &
                       qr_bt(np,lp,2), &
                       qs_bt(np,lp,2), &
                       qi_bt(np,lp,2), &
                       qg_bt(np,lp,2)
   REAL(KIND=euwp)    :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp 
   REAL(KIND=euwp)    :: hinv 
   INTEGER(KIND=iintegers) :: i,j,k
        IF(gndedge == 1) THEN
          k=1
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qr(i,j,k  ), qr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qi(i,j,k  ), qi(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qs(i  ,j,k), qs(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qs(i-1,j,k), qs(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qs(i,j  ,k), qs(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qs(i,j-1,k), qs(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qs(i,j,k  ), qs(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qsant(i,j,k)=rhr(i,j,k)*                         &
                                    (qs(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qg(i  ,j,k), qg(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qg(i-1,j,k), qg(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qg(i,j  ,k), qg(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qg(i,j-1,k), qg(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qg(i,j,k  ), qg(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qgant(i,j,k)=rhr(i,j,k)*                         &
                                    (qg(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !gndedge

InnerZFullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
InnerZFullXYDomainLoopDC(
                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qr(i,j,k  ), qr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qr(i,j,k-1), qr(i,j,k  ), wadv(i  ,j  ,k  ));

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);


                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qi(i,j,k  ), qi(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qi(i,j,k-1), qi(i,j,k  ), wadv(i  ,j  ,k  ));

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qs(i  ,j,k), qs(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qs(i-1,j,k), qs(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qs(i,j  ,k), qs(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qs(i,j-1,k), qs(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qs(i,j,k  ), qs(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qs(i,j,k-1), qs(i,j,k  ), wadv(i  ,j  ,k  ));

                       qsant(i,j,k)=rhr(i,j,k)*                         &
                                    (qs(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qg(i  ,j,k), qg(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qg(i-1,j,k), qg(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qg(i,j  ,k), qg(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qg(i,j-1,k), qg(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qg(i,j,k  ), qg(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qg(i,j,k-1), qg(i,j,k  ), wadv(i  ,j  ,k  ));

                       qgant(i,j,k)=rhr(i,j,k)*                         &
                                    (qg(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)

      IF(skyedge==1) THEN
         k=lp
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qr(i,j,k-1), qr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qi(i,j,k-1), qi(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qs(i  ,j,k), qs(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qs(i-1,j,k), qs(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qs(i,j  ,k), qs(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qs(i,j-1,k), qs(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qs(i,j,k-1), qs(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qsant(i,j,k)=rhr(i,j,k)*                         &
                                    (qs(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qg(i  ,j,k), qg(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qg(i-1,j,k), qg(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qg(i,j  ,k), qg(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qg(i,j-1,k), qg(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qg(i,j,k-1), qg(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qgant(i,j,k)=rhr(i,j,k)*                         &
                                    (qg(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !skyedge

   END SUBROUTINE upwind_firstep_GPUBC_qg

   SUBROUTINE upwind_firstep_GPUBC_qs( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt,    &
                           qr,  qr_lr,  qr_bt,    &
                           qi,  qi_lr,  qi_bt,    &
                           qs,  qs_lr,  qs_bt)     

   USE scratch_datafields, ONLY: thant,uant,vant,want,thtant,pstrant
   USE scratch_datafields, ONLY: qvant,qcant,qrant,qiant,qsant
   USE mpi_parallel, ONLY: gndedge,skyedge

   INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
   INTEGER(KIND=iintegers) :: istat


   REAL_euwp, INTENT(IN) ::                               & 
                      rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),    & 
                      rhoadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::                            & 
                    uadv(1-ih:np+1+ih,1-ih:mp+ih,1-ih:lp+ih),   & 
                    vadv(1-ih:np+ih,1-ih:mp+1+ih,1-ih:lp+ih),   &
                    wadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+1+ih),   &
                      th(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       u(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                       v(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       w(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                     tht(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    pstr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN), OPTIONAL  ::               &
                      qv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qc(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &   
                      qr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                      qi(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qs(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::    &                
                       th_lr(mp,lp,2), &
                        u_lr(mp,lp,2), &
                        v_lr(mp,lp,2), &
                        w_lr(mp,lp,2), &
                      tht_lr(mp,lp,2), &
                     pstr_lr(mp,lp,2)

   REAL_euwp, INTENT(IN), OPTIONAL ::    & 
                       qv_lr(mp,lp,2), &
                       qc_lr(mp,lp,2), &
                       qr_lr(mp,lp,2), &
                       qi_lr(mp,lp,2), &
                       qs_lr(mp,lp,2)

   REAL_euwp, INTENT(IN) ::    &                
                       th_bt(np,lp,2), &
                        u_bt(np,lp,2), &
                        v_bt(np,lp,2), &
                        w_bt(np,lp,2), &
                      tht_bt(np,lp,2), &
                     pstr_bt(np,lp,2)   
   
   REAL_euwp, INTENT(IN), OPTIONAL ::    &            
                       qv_bt(np,lp,2), &
                       qc_bt(np,lp,2), &
                       qr_bt(np,lp,2), &
                       qi_bt(np,lp,2), &
                       qs_bt(np,lp,2)

   REAL(KIND=euwp)    :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp 
   REAL(KIND=euwp)    :: hinv 
   INTEGER(KIND=iintegers) :: i,j,k
        IF(gndedge == 1) THEN
          k=1
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qr(i,j,k  ), qr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qi(i,j,k  ), qi(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qs(i  ,j,k), qs(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qs(i-1,j,k), qs(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qs(i,j  ,k), qs(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qs(i,j-1,k), qs(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qs(i,j,k  ), qs(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qsant(i,j,k)=rhr(i,j,k)*                         &
                                    (qs(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !gndedge

InnerZFullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qr(i,j,k  ), qr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qr(i,j,k-1), qr(i,j,k  ), wadv(i  ,j  ,k  ));

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);


                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qi(i,j,k  ), qi(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qi(i,j,k-1), qi(i,j,k  ), wadv(i  ,j  ,k  ));

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qs(i  ,j,k), qs(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qs(i-1,j,k), qs(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qs(i,j  ,k), qs(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qs(i,j-1,k), qs(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qs(i,j,k  ), qs(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qs(i,j,k-1), qs(i,j,k  ), wadv(i  ,j  ,k  ));

                       qsant(i,j,k)=rhr(i,j,k)*                         &
                                    (qs(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)

      IF(skyedge==1) THEN
         k=lp
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qr(i,j,k-1), qr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qi(i,j,k-1), qi(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qs(i  ,j,k), qs(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qs(i-1,j,k), qs(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qs(i,j  ,k), qs(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qs(i,j-1,k), qs(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qs(i,j,k-1), qs(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qsant(i,j,k)=rhr(i,j,k)*                         &
                                    (qs(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
          ENDIF !skyedge

   END SUBROUTINE upwind_firstep_GPUBC_qs

   SUBROUTINE upwind_firstep_GPUBC_qi( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt,    &
                           qr,  qr_lr,  qr_bt,    &
                           qi,  qi_lr,  qi_bt)    

   USE scratch_datafields, ONLY: thant,uant,vant,want,thtant,pstrant
   USE scratch_datafields, ONLY: qvant,qcant,qrant,qiant
   USE mpi_parallel, ONLY: gndedge,skyedge

   INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
   INTEGER(KIND=iintegers) :: istat


   REAL_euwp, INTENT(IN) ::                               & 
                      rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),    & 
                      rhoadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::                            & 
                    uadv(1-ih:np+1+ih,1-ih:mp+ih,1-ih:lp+ih),   & 
                    vadv(1-ih:np+ih,1-ih:mp+1+ih,1-ih:lp+ih),   &
                    wadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+1+ih),   &
                      th(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       u(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                       v(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       w(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                     tht(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    pstr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN), OPTIONAL  ::               &
                      qv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qc(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &   
                      qr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                      qi(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) 

   REAL_euwp, INTENT(IN) ::    &                
                       th_lr(mp,lp,2), &
                        u_lr(mp,lp,2), &
                        v_lr(mp,lp,2), &
                        w_lr(mp,lp,2), &
                      tht_lr(mp,lp,2), &
                     pstr_lr(mp,lp,2)

   REAL_euwp, INTENT(IN), OPTIONAL ::    & 
                       qv_lr(mp,lp,2), &
                       qc_lr(mp,lp,2), &
                       qr_lr(mp,lp,2), &
                       qi_lr(mp,lp,2)

   REAL_euwp, INTENT(IN) ::    &                
                       th_bt(np,lp,2), &
                        u_bt(np,lp,2), &
                        v_bt(np,lp,2), &
                        w_bt(np,lp,2), &
                      tht_bt(np,lp,2), &
                     pstr_bt(np,lp,2)   
   
   REAL_euwp, INTENT(IN), OPTIONAL ::    &            
                       qv_bt(np,lp,2), &
                       qc_bt(np,lp,2), &
                       qr_bt(np,lp,2), &
                       qi_bt(np,lp,2)

   REAL(KIND=euwp)    :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp 
   REAL(KIND=euwp)    :: hinv 
   INTEGER(KIND=iintegers) :: i,j,k
        IF(gndedge == 1) THEN
          k=1
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qr(i,j,k  ), qr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qi(i,j,k  ), qi(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !gndedge

InnerZFullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qr(i,j,k  ), qr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qr(i,j,k-1), qr(i,j,k  ), wadv(i  ,j  ,k  ));

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);


                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qi(i,j,k  ), qi(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qi(i,j,k-1), qi(i,j,k  ), wadv(i  ,j  ,k  ));

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)

      IF(skyedge==1) THEN
         k=lp
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qr(i,j,k-1), qr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qi(i,j,k-1), qi(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
          ENDIF !skyedge

   END SUBROUTINE upwind_firstep_GPUBC_qi

   SUBROUTINE upwind_firstep_GPUBC_qr( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt,    &
                           qr,  qr_lr,  qr_bt)    

   USE scratch_datafields, ONLY: thant,uant,vant,want,thtant,pstrant
   USE scratch_datafields, ONLY: qvant,qcant,qrant
   USE mpi_parallel, ONLY: gndedge,skyedge

   INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
   INTEGER(KIND=iintegers) :: istat


   REAL_euwp, INTENT(IN) ::                               & 
                      rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),    & 
                      rhoadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::                            & 
                    uadv(1-ih:np+1+ih,1-ih:mp+ih,1-ih:lp+ih),   & 
                    vadv(1-ih:np+ih,1-ih:mp+1+ih,1-ih:lp+ih),   &
                    wadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+1+ih),   &
                      th(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       u(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                       v(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       w(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                     tht(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    pstr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN), OPTIONAL  ::               &
                      qv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qc(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &   
                      qr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::    &                
                       th_lr(mp,lp,2), &
                        u_lr(mp,lp,2), &
                        v_lr(mp,lp,2), &
                        w_lr(mp,lp,2), &
                      tht_lr(mp,lp,2), &
                     pstr_lr(mp,lp,2)

   REAL_euwp, INTENT(IN), OPTIONAL ::    & 
                       qv_lr(mp,lp,2), &
                       qc_lr(mp,lp,2), &
                       qr_lr(mp,lp,2)

   REAL_euwp, INTENT(IN) ::    &                
                       th_bt(np,lp,2), &
                        u_bt(np,lp,2), &
                        v_bt(np,lp,2), &
                        w_bt(np,lp,2), &
                      tht_bt(np,lp,2), &
                     pstr_bt(np,lp,2)   
   
   REAL_euwp, INTENT(IN), OPTIONAL ::    &            
                       qv_bt(np,lp,2), &
                       qc_bt(np,lp,2), &
                       qr_bt(np,lp,2)

   REAL(KIND=euwp)    :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp 
   REAL(KIND=euwp)    :: hinv 
   INTEGER(KIND=iintegers) :: i,j,k
        IF(gndedge == 1) THEN
          k=1
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qr(i,j,k  ), qr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !gndedge

InnerZFullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qr(i,j,k  ), qr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qr(i,j,k-1), qr(i,j,k  ), wadv(i  ,j  ,k  ));

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)

      IF(skyedge==1) THEN
         k=lp
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qr(i,j,k-1), qr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !skyedge

   END SUBROUTINE upwind_firstep_GPUBC_qr

   SUBROUTINE upwind_firstep_GPUBC_qc( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt)    

   USE scratch_datafields, ONLY: thant,uant,vant,want,thtant,pstrant
   USE scratch_datafields, ONLY: qvant,qcant
   USE mpi_parallel, ONLY: gndedge,skyedge

   INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
   INTEGER(KIND=iintegers) :: istat


   REAL_euwp, INTENT(IN) ::                               & 
                      rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),    & 
                      rhoadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::                            & 
                    uadv(1-ih:np+1+ih,1-ih:mp+ih,1-ih:lp+ih),   & 
                    vadv(1-ih:np+ih,1-ih:mp+1+ih,1-ih:lp+ih),   &
                    wadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+1+ih),   &
                      th(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       u(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                       v(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       w(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                     tht(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    pstr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN), OPTIONAL  ::               &
                      qv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qc(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::    &                
                       th_lr(mp,lp,2), &
                        u_lr(mp,lp,2), &
                        v_lr(mp,lp,2), &
                        w_lr(mp,lp,2), &
                      tht_lr(mp,lp,2), &
                     pstr_lr(mp,lp,2)

   REAL_euwp, INTENT(IN), OPTIONAL ::    & 
                       qv_lr(mp,lp,2), &
                       qc_lr(mp,lp,2)

   REAL_euwp, INTENT(IN) ::    &                
                       th_bt(np,lp,2), &
                        u_bt(np,lp,2), &
                        v_bt(np,lp,2), &
                        w_bt(np,lp,2), &
                      tht_bt(np,lp,2), &
                     pstr_bt(np,lp,2)   
   
   REAL_euwp, INTENT(IN), OPTIONAL ::    &            
                       qv_bt(np,lp,2), &
                       qc_bt(np,lp,2)

   REAL(KIND=euwp)    :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp 
   REAL(KIND=euwp)    :: hinv 
   INTEGER(KIND=iintegers) :: i,j,k
        IF(gndedge == 1) THEN
          k=1
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !gndedge

InnerZFullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)

      IF(skyedge==1) THEN
         k=lp
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !skyedge

   END SUBROUTINE upwind_firstep_GPUBC_qc

   SUBROUTINE upwind_firstep_GPUBC_qv( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt)    

   USE scratch_datafields, ONLY: thant,uant,vant,want,thtant,pstrant
   USE scratch_datafields, ONLY: qvant
   USE mpi_parallel, ONLY: gndedge,skyedge

   INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
   INTEGER(KIND=iintegers) :: istat


   REAL_euwp, INTENT(IN) ::                               & 
                      rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),    & 
                      rhoadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::                            & 
                    uadv(1-ih:np+1+ih,1-ih:mp+ih,1-ih:lp+ih),   & 
                    vadv(1-ih:np+ih,1-ih:mp+1+ih,1-ih:lp+ih),   &
                    wadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+1+ih),   &
                      th(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       u(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                       v(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       w(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                     tht(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    pstr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN), OPTIONAL  ::               &
                      qv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::    &                
                       th_lr(mp,lp,2), &
                        u_lr(mp,lp,2), &
                        v_lr(mp,lp,2), &
                        w_lr(mp,lp,2), &
                      tht_lr(mp,lp,2), &
                     pstr_lr(mp,lp,2)

   REAL_euwp, INTENT(IN), OPTIONAL ::    & 
                       qv_lr(mp,lp,2)

   REAL_euwp, INTENT(IN) ::    &                
                       th_bt(np,lp,2), &
                        u_bt(np,lp,2), &
                        v_bt(np,lp,2), &
                        w_bt(np,lp,2), &
                      tht_bt(np,lp,2), &
                     pstr_bt(np,lp,2)   
   
   REAL_euwp, INTENT(IN), OPTIONAL ::    &            
                       qv_bt(np,lp,2)

   REAL(KIND=euwp)    :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp 
   REAL(KIND=euwp)    :: hinv 
   INTEGER(KIND=iintegers) :: i,j,k
        IF(gndedge == 1) THEN
          k=1
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !gndedge

InnerZFullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)

      IF(skyedge==1) THEN
         k=lp
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
          ENDIF !skyedge

   END SUBROUTINE upwind_firstep_GPUBC_qv

   SUBROUTINE upwind_firstep_GPUBC_dry( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt)     

   USE scratch_datafields, ONLY: thant,uant,vant,want,thtant,pstrant
   USE scratch_datafields, ONLY: qvant
   USE mpi_parallel, ONLY: gndedge,skyedge

   INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
   INTEGER(KIND=iintegers) :: istat


   REAL_euwp, INTENT(IN) ::                               & 
                      rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),    & 
                      rhoadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::                            & 
                    uadv(1-ih:np+1+ih,1-ih:mp+ih,1-ih:lp+ih),   & 
                    vadv(1-ih:np+ih,1-ih:mp+1+ih,1-ih:lp+ih),   &
                    wadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+1+ih),   &
                      th(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       u(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                       v(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       w(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                     tht(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    pstr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)


   REAL_euwp, INTENT(IN) ::    &                
                       th_lr(mp,lp,2), &
                        u_lr(mp,lp,2), &
                        v_lr(mp,lp,2), &
                        w_lr(mp,lp,2), &
                      tht_lr(mp,lp,2), &
                     pstr_lr(mp,lp,2)

   REAL_euwp, INTENT(IN) ::    &                
                       th_bt(np,lp,2), &
                        u_bt(np,lp,2), &
                        v_bt(np,lp,2), &
                        w_bt(np,lp,2), &
                      tht_bt(np,lp,2), &
                     pstr_bt(np,lp,2)   

   REAL(KIND=euwp)    :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp 
   REAL(KIND=euwp)    :: hinv 
   INTEGER(KIND=iintegers) :: i,j,k
        IF(gndedge == 1) THEN
          k=1
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
          ENDIF !gndedge

InnerZFullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)

      IF(skyedge==1) THEN
         k=lp
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !skyedge

   END SUBROUTINE upwind_firstep_GPUBC_dry

   SUBROUTINE upwind_firstep_GPUBC_qg_multiker( np,mp,lp,ih,  & 
                           rhr,rhoadv,uadv,vadv,wadv, &     
                           th,  th_lr,  th_bt,    &
                            u,   u_lr,   u_bt,    &
                            v,   v_lr,   v_bt,    &
                            w,   w_lr,   w_bt,    &
                          tht, tht_lr, tht_bt,    &
                         pstr,pstr_lr,pstr_bt,    &
                           qv,  qv_lr,  qv_bt,    &
                           qc,  qc_lr,  qc_bt,    &
                           qr,  qr_lr,  qr_bt,    &
                           qi,  qi_lr,  qi_bt,    &
                           qs,  qs_lr,  qs_bt,    &
                           qg,  qg_lr,  qg_bt) 

   USE scratch_datafields, ONLY: thant,uant,vant,want,thtant,pstrant
   USE scratch_datafields, ONLY: qvant,qcant,qrant,qiant,qsant,qgant
   USE mpi_parallel, ONLY: gndedge,skyedge

   INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
   INTEGER(KIND=iintegers) :: istat


   REAL_euwp, INTENT(IN) ::                               & 
                      rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),    & 
                      rhoadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN) ::                            & 
                    uadv(1-ih:np+1+ih,1-ih:mp+ih,1-ih:lp+ih),   & 
                    vadv(1-ih:np+ih,1-ih:mp+1+ih,1-ih:lp+ih),   &
                    wadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+1+ih),   &
                      th(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       u(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                       v(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                       w(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                     tht(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    pstr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

   REAL_euwp, INTENT(IN), OPTIONAL  ::               &
                      qv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qc(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &   
                      qr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                      qs(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qi(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     & 
                      qg(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)   

   REAL_euwp, INTENT(IN) ::    &                
                       th_lr(mp,lp,2), &
                        u_lr(mp,lp,2), &
                        v_lr(mp,lp,2), &
                        w_lr(mp,lp,2), &
                      tht_lr(mp,lp,2), &
                     pstr_lr(mp,lp,2)

   REAL_euwp, INTENT(IN), OPTIONAL ::    & 
                       qv_lr(mp,lp,2), &
                       qc_lr(mp,lp,2), &
                       qr_lr(mp,lp,2), &
                       qs_lr(mp,lp,2), &
                       qi_lr(mp,lp,2), &
                       qg_lr(mp,lp,2)

   REAL_euwp, INTENT(IN) ::    &                
                       th_bt(np,lp,2), &
                        u_bt(np,lp,2), &
                        v_bt(np,lp,2), &
                        w_bt(np,lp,2), &
                      tht_bt(np,lp,2), &
                     pstr_bt(np,lp,2)   
   
   REAL_euwp, INTENT(IN), OPTIONAL ::    &            
                       qv_bt(np,lp,2), &
                       qc_bt(np,lp,2), &
                       qr_bt(np,lp,2), &
                       qs_bt(np,lp,2), &
                       qi_bt(np,lp,2), &
                       qg_bt(np,lp,2)
   REAL(KIND=euwp)    :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp 
   REAL(KIND=euwp)    :: hinv 
   INTEGER(KIND=iintegers) :: i,j,k
        IF(gndedge == 1) THEN
          k=1
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
FullXYDomainLoopDC(
                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
FullXYDomainLoopDC(
                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
FullXYDomainLoopDC(
                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
FullXYDomainLoopDC(
                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qr(i,j,k  ), qr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qi(i,j,k  ), qi(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(qs(i  ,j,k), qs(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qs(i-1,j,k), qs(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qs(i,j  ,k), qs(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qs(i,j-1,k), qs(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qs(i,j,k  ), qs(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qsant(i,j,k)=rhr(i,j,k)*                         &
                                    (qs(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(qg(i  ,j,k), qg(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qg(i-1,j,k), qg(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qg(i,j  ,k), qg(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qg(i,j-1,k), qg(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qg(i,j,k  ), qg(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =-f3ijkp;

                       qgant(i,j,k)=rhr(i,j,k)*                         &
                                    (qg(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !gndedge

InnerZFullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(th(i,j,k  ), th(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
InnerZFullXYDomainLoopDC(
                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( u(i,j,k  ),  u(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));

                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
InnerZFullXYDomainLoopDC(
                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( v(i,j,k  ),  v(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));

                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
InnerZFullXYDomainLoopDC(

                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor( w(i,j,k  ),  w(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));

                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
InnerZFullXYDomainLoopDC(

                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(tht(i,j,k  ), tht(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));

                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
InnerZFullXYDomainLoopDC(

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(pstr(i,j,k  ), pstr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));

                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
InnerZFullXYDomainLoopDC(
                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qv(i,j,k  ), qv(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));

                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
InnerZFullXYDomainLoopDC(

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qc(i,j,k  ), qc(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
InnerZFullXYDomainLoopDC(
                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qr(i,j,k  ), qr(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qr(i,j,k-1), qr(i,j,k  ), wadv(i  ,j  ,k  ));

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
InnerZFullXYDomainLoopDC(

                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qi(i,j,k  ), qi(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qi(i,j,k-1), qi(i,j,k  ), wadv(i  ,j  ,k  ));

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
InnerZFullXYDomainLoopDC(
                    f1ijkp=donor(qs(i  ,j,k), qs(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qs(i-1,j,k), qs(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qs(i,j  ,k), qs(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qs(i,j-1,k), qs(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qs(i,j,k  ), qs(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qs(i,j,k-1), qs(i,j,k  ), wadv(i  ,j  ,k  ));

                       qsant(i,j,k)=rhr(i,j,k)*                         &
                                    (qs(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
InnerZFullXYDomainLoopDC(

                    f1ijkp=donor(qg(i  ,j,k), qg(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qg(i-1,j,k), qg(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qg(i,j  ,k), qg(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qg(i,j-1,k), qg(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijkp=donor(qg(i,j,k  ), qg(i,j,k+1), wadv(i  ,j  ,k+1));
                    f3ijk =donor(qg(i,j,k-1), qg(i,j,k  ), wadv(i  ,j  ,k  ));

                       qgant(i,j,k)=rhr(i,j,k)*                         &
                                    (qg(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)

      IF(skyedge==1) THEN
         k=lp
FullXYDomainLoopDC(
                    hinv=1._euwp/rhoadv(i,j,k);

                    f1ijkp=donor(th(i  ,j,k), th(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(th(i-1,j,k), th(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(th(i,j  ,k), th(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(th(i,j-1,k), th(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(th(i,j,k-1), th(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       thant(i,j,k)=rhr(i,j,k)*                         &
                                    (th(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(
                    f1ijkp=donor( u(i  ,j,k),  u(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( u(i-1,j,k),  u(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( u(i,j  ,k),  u(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( u(i,j-1,k),  u(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( u(i,j,k-1),  u(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        uant(i,j,k)=rhr(i,j,k)*                         &
                                    ( u(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
FullXYDomainLoopDC(
                    f1ijkp=donor( v(i  ,j,k),  v(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( v(i-1,j,k),  v(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( v(i,j  ,k),  v(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( v(i,j-1,k),  v(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( v(i,j,k-1),  v(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        vant(i,j,k)=rhr(i,j,k)*                         &
                                    ( v(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
FullXYDomainLoopDC(
                    f1ijkp=donor( w(i  ,j,k),  w(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor( w(i-1,j,k),  w(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor( w(i,j  ,k),  w(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor( w(i,j-1,k),  w(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor( w(i,j,k-1),  w(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                        want(i,j,k)=rhr(i,j,k)*                         &
                                    ( w(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);

)
FullXYDomainLoopDC(
                    f1ijkp=donor(tht(i  ,j,k), tht(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(tht(i-1,j,k), tht(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(tht(i,j  ,k), tht(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(tht(i,j-1,k), tht(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(tht(i,j,k-1), tht(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      thtant(i,j,k)=rhr(i,j,k)*                         &
                                   (tht(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(pstr(i  ,j,k), pstr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(pstr(i-1,j,k), pstr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(pstr(i,j  ,k), pstr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(pstr(i,j-1,k), pstr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(pstr(i,j,k-1), pstr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                      pstrant(i,j,k)=rhr(i,j,k)*                        &
                                   (pstr(i,j,k)-                        &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(qv(i  ,j,k), qv(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qv(i-1,j,k), qv(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qv(i,j  ,k), qv(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qv(i,j-1,k), qv(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qv(i,j,k-1), qv(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;


                       qvant(i,j,k)=rhr(i,j,k)*                         &
                                    (qv(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(qc(i  ,j,k), qc(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qc(i-1,j,k), qc(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qc(i,j  ,k), qc(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qc(i,j-1,k), qc(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qc(i,j,k-1), qc(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qcant(i,j,k)=rhr(i,j,k)*                         &
                                    (qc(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(qr(i  ,j,k), qr(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qr(i-1,j,k), qr(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qr(i,j  ,k), qr(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qr(i,j-1,k), qr(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qr(i,j,k-1), qr(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qrant(i,j,k)=rhr(i,j,k)*                         &
                                    (qr(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(qi(i  ,j,k), qi(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qi(i-1,j,k), qi(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qi(i,j  ,k), qi(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qi(i,j-1,k), qi(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qi(i,j,k-1), qi(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qiant(i,j,k)=rhr(i,j,k)*                         &
                                    (qi(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(qs(i  ,j,k), qs(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qs(i-1,j,k), qs(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qs(i,j  ,k), qs(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qs(i,j-1,k), qs(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qs(i,j,k-1), qs(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qsant(i,j,k)=rhr(i,j,k)*                         &
                                    (qs(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
FullXYDomainLoopDC(

                    f1ijkp=donor(qg(i  ,j,k), qg(i+1,j,k), uadv(i+1,j,  k  ));
                    f1ijk =donor(qg(i-1,j,k), qg(i  ,j,k), uadv(i  ,j  ,k  ));
                    f2ijkp=donor(qg(i,j  ,k), qg(i,j+1,k), vadv(i  ,j+1,k  ));
                    f2ijk =donor(qg(i,j-1,k), qg(i,j  ,k), vadv(i  ,j  ,k  ));
                    f3ijk =donor(qg(i,j,k-1), qg(i,j,k  ), wadv(i  ,j  ,k  ));
                    f3ijkp=-f3ijk;

                       qgant(i,j,k)=rhr(i,j,k)*                         &
                                    (qg(i,j,k)-                         &
                                             ( f1ijkp-f1ijk             &
                                              +f2ijkp-f2ijk             &
                                              +f3ijkp-f3ijk )*hinv);
)
          ENDIF !skyedge

   END SUBROUTINE upwind_firstep_GPUBC_qg_multiker
END MODULE advection_CEtimestep_twostep_GPUBC


MODULE velprd_driver_module 
    USE precisions, ONLY  : iintegers,euwp
    USE mpi_parallel, ONLY: mype
#ifdef PNETCDF
    USE mpi_parallel, ONLY: pnet_out_chunk
#endif
 USE eulag_diagnostics, ONLY: checknans
    IMPLICIT NONE

CONTAINS
#include "defines.inc"
!#include "mpdataoperators.inc"

    SUBROUTINE velprd_driver(nsubsteps,nsubsteps_o,lsubstepping,       &
                             dtn,dtnold,ltimeadapt,                   &
                             lvertsplit,lcrdiag_nphlf,crmaxlim,dtmax,  &
                             ox,oy,oz,u,v,w,uadv,vadv,wadv,h,          &
                             ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih,    & 
                             dxi,dyi,dzi,dt)
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY: update,updatelr,updatebt,updategs,update3
        USE mpi_parallel, ONLY: update_multi,update3_multi
        USE mpi_parallel, ONLY: ttbeg,ttend
        USE eulag_diagnostics, ONLY: diagnos_advection, diagnos_advection_vz
!       USE eulag_diagnostics, ONLY: exam_var,courants_crash
        USE eulag_diagutils, ONLY: compute_courants
        USE eulag_diagutils, ONLY: print_courants
        USE module_velprd, ONLY: velprd_traj0,velprd_traj0_adapt_dt
        USE module_velprd, ONLY: velprdA_to_C 
        USE scratch_datafields, ONLY: pfx,pfy
        LOGICAL, INTENT(IN) :: lsubstepping,ltimeadapt,lvertsplit,lcrdiag_nphlf
        INTEGER(KIND=iintegers), INTENT(IN) :: np,mp,lp,ih
        INTEGER(KIND=iintegers), INTENT(IN) :: ipoles0,ibcx0,ibcy0,ibcz0 
        INTEGER(KIND=iintegers), INTENT(INOUT) :: nsubsteps,nsubsteps_o
        REAL(KIND=euwp), INTENT(INOUT) :: dtn,dtnold

        REAL_euwp, INTENT(INOUT) ::  &
                          ox(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),    &
                          oy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),    &
                          oz(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)

        REAL_euwp, INTENT(INOUT) ::  &
                          u(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),    &
                          v(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),    &
                          w(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)

        REAL_euwp, INTENT(OUT) ::  &
                          uadv(1-ih:np+1+ih, 1-ih:mp+ih, 1-ih:lp+ih),    &
                          vadv(1-ih:np+ih, 1-ih:mp+1+ih, 1-ih:lp+ih),    &
                          wadv(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih)


        REAL_euwp, INTENT(IN) ::  &
                             h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

        REAL(KIND=euwp),INTENT(IN) :: dxi,dyi,dzi,dt
        REAL(KIND=euwp),INTENT(IN) :: crmaxlim,dtmax 

        REAL(KIND=euwp), PARAMETER :: cr_fatal=10. 
        LOGICAL, PARAMETER :: lprint=.TRUE.
        LOGICAL ladapt_cr_ok,lsub_cr_ok
        REAL(KIND=euwp) :: cr1,cr2
        REAL(KIND=euwp) :: cr3d,crxy,crxz,cryz,crx,cry,crz 
        REAL(KIND=euwp) :: cr_for_dt
        REAL(KIND=euwp) :: subfac,dt_ratio
        REAL(KIND=euwp) :: dtntemp,dtnoldtemp,crmaxlimtemp 
        INTEGER(KIND=iintegers) :: i,j,k 

        CALL ttbeg(9)

        IF(lsubstepping.AND.ltimeadapt) STOP "substepping and adaptivity not supported together"

!---------------------------------------------------------------------
! Courant initial diagnostics
!---------------------------------------------------------------------

!1. Compute ratio of full timestep and previously used timestep (dt/number of substeps)
        dt_ratio=1._euwp
        dtntemp=dtn
!If first (potential) substep, then account for substepping in previous timestep
        IF(lsubstepping.AND.nsubsteps.eq.1) dt_ratio=(dt/(0.*nsubsteps+1))/(dt/nsubsteps_o)
!       IF(ltimeadapt)                      dt_ratio=dtn/dtnold

!2. Compute first order prediction of advecting velocities
!   Using  basic dt (substepping) or latest dt (adaptive)
        IF(lsubstepping) dtntemp=dt
        IF(ltimeadapt  ) dtntemp=dtn

        IF(lcrdiag_nphlf) THEN
! Diagnosing Courant at n+1/2                
          CALL velprd_traj0_adapt_dt(ox,oy,oz,h,dxi,dyi,dzi,dtntemp,dt_ratio,np,mp,lp,ih)
          CALL compute_courants(cr3d,crxy,crxz,cryz,crx,cry,crz,ox,oy,oz,h, &
                                2,dxi,dyi,dzi,dtntemp,np,mp,lp,ih)
          IF(mype.eq.0.AND.lprint) CALL  print_courants(3,cr3d,crxy,crxz,cryz,crx,cry,crz)

        ELSE
! Diagnosing Courant at n (cheaper)                
          CALL compute_courants(cr3d,crxy,crxz,cryz,crx,cry,crz,ox,oy,oz,h, &
                                1,dxi,dyi,dzi,dtntemp,np,mp,lp,ih)
          IF(mype.eq.0.AND.lprint) CALL  print_courants(4,cr3d,crxy,crxz,cryz,crx,cry,crz)
        ENDIF
! Store cr in aux. variable for further manipulations

        cr_for_dt=cr3d
!3. Account for the use of split schemes
        IF(lvertsplit) cr_for_dt=max(crxy,.5_euwp*crz)+TINY(cr_for_dt)
!        
!---------------------------------------------------------------------
! Algorithm to determine substepping 
!---------------------------------------------------------------------
    IF(lsubstepping) THEN
                                lsub_cr_ok=.TRUE.
       IF ((cr_for_dt > crmaxlim) .AND. nsubsteps == 1) THEN
!5a.  Initial subcycling guess to be later verified  
         nsubsteps=CEILING(cr_for_dt)
         IF (cr_for_dt/nsubsteps > crmaxlim) nsubsteps=nsubsteps+1 
         IF(mype.eq.0) PRINT *, "EULAG Courant numbers(before accounting for substepping):"
         IF(mype.eq.0) PRINT "(1x,'cour:',1(e23.16,2X))",cr_for_dt
                                lsub_cr_ok=.FALSE.

       ELSE IF ((cr_for_dt/nsubsteps > cr_fatal) .AND. nsubsteps > 1) THEN
!5b. If Courant number is too large and we are already in substepping mode, then crash
!       CALL courants_crash(ox,oy,oz,h,dxi,dyi,dzi,dt,np,mp,lp,ih)
         STOP 'Courant > 1 for nsubsteps >1'

       ELSE IF ((cr_for_dt/nsubsteps .LE. crmaxlim) .AND. nsubsteps > 1) THEN
         IF(mype.eq.0) PRINT *, "Looks like we are inside substepping and it works fine."
       ELSE 
         IF(mype.eq.0) PRINT *, "No substepping needed as cr stays below ",crmaxlim
         IF(mype.eq.0) PRINT "(1x,'cour:',1(f8.3,2X))",cr_for_dt
       ENDIF ! React to courant in substepping

       dtn=dt/nsubsteps

!6. Check if substepping adequate, or more steps is need. 
       IF(lvertsplit) THEN
         dt_ratio=dtn/dtnold
         CALL velprd_traj0(ox,oy,oz,uadv,vadv,wadv,h,dxi,dyi,dzi,dtn,dt_ratio,  &
                        ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
         CALL compute_courants(cr3d,crxy,crxz,cryz,crx,cry,crz,ox,oy,oz,h, &
                                 2,dxi,dyi,dzi,dtn,np,mp,lp,ih)
         cr_for_dt=cr3d
       ENDIF 
       subcrloop: DO WHILE(.NOT.lsub_cr_ok)
         dt_ratio=dtn/dtnold
         CALL velprd_traj0(ox,oy,oz,uadv,vadv,wadv,h,dxi,dyi,dzi,dtn,dt_ratio,  &
                        ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
         CALL compute_courants(cr3d,crxy,crxz,cryz,crx,cry,crz,ox,oy,oz,h, &
                                 2,dxi,dyi,dzi,dtn,np,mp,lp,ih)

         cr_for_dt=cr3d
         IF(lvertsplit) cr_for_dt=max(crxy,.5_euwp*crz)+TINY(cr_for_dt)

         IF(cr_for_dt.le.crmaxlim) THEN
            lsub_cr_ok=.TRUE.
         ELSE
            nsubsteps=nsubsteps+1
            dtn=dt/nsubsteps
         ENDIF
           
         IF(mype.eq.0.AND.lprint) THEN
            CALL  print_courants(2,cr3d,crxy,crxz,cryz,crx,cry,crz)
            PRINT "(1x,'         crmaxlim,dt,dt_ratio,nsubsteps:',3(f13.6,1X),I2)",          &
                                 crmaxlim,dtn,dt_ratio,nsubsteps
          ENDIF

                 ENDDO subcrloop
       IF(nsubsteps.ne.nsubsteps_o.AND.mype == 0)  &
          print *,'dtratio,nsub,nsub_o', dt_ratio,nsubsteps,nsubsteps_o

      nsubsteps_o=nsubsteps
      dtnold=dtn
      dtn=dt/nsubsteps
 
      IF (nsubsteps > 1 .OR. (.NOT.lcrdiag_nphlf)) THEN
!7a. We need to recompute prediction of advecting velocities, 
!   now we know the true dt_ratio between the previous and current timestep
       CALL velprd_traj0(ox,oy,oz,uadv,vadv,wadv,h,dxi,dyi,dzi,dtn,dt_ratio,  &
                      ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
      ELSE
!7b. or it is enough to compute C-grid advectors, since A-grid advectors are unchanged
        CALL velprdA_to_C(ox(:,:,:,1),oy(:,:,:,1),oz(:,:,:,1),uadv,vadv,wadv, &
               ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih )
      ENDIF !nsubsteps

!---------------------------------------------------------------------
! Algorithm to determine adaptive timestep 
!---------------------------------------------------------------------
!--------------------------------------------------------------------------------------
    ELSE IF(ltimeadapt) THEN
!--------------------------------------------------------------------------------------

!     dtn = min((crmaxlim/(abs(cr_for_dt)+TINY(cr_for_dt)))*dt,dtmax) ! proposed next timestep
       crmaxlimtemp=crmaxlim
       dtnold=dtn

                    ladapt_cr_ok=.FALSE.
       adaptcrloop: DO WHILE(.NOT.ladapt_cr_ok)
!Evaluate proposed next timestep
           dtntemp = min((crmaxlimtemp/(abs(cr_for_dt)+TINY(cr_for_dt)))*dtn,dtmax) 
!Optional avoid adaptation if the increase of dt is too small
        IF(dtntemp.gt.dtnold.AND.((dtntemp-dtnold)/dtnold).gt.0.1) dtn=dtntemp 
!Always lower the timestep if needed
        IF(dtntemp.le.dtnold) dtn=dtntemp 
        dt_ratio=dtn/dtnold

        CALL velprd_traj0(ox,oy,oz,uadv,vadv,wadv,h,dxi,dyi,dzi,dtn,dt_ratio,  &
                          ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
        CALL compute_courants(cr3d,crxy,crxz,cryz,crx,cry,crz,ox,oy,oz,h, &
                              2,dxi,dyi,dzi,dtn,np,mp,lp,ih)

        cr_for_dt=cr3d
        IF(lvertsplit) cr_for_dt=max(crxy,.5_euwp*crz)+TINY(cr_for_dt)

         IF(cr_for_dt.le.crmaxlim) THEN
              ladapt_cr_ok=.TRUE.
         ELSE
!Lower the courant limit temporarily, so the cr computed at n+1/2 hopefully stays below limit  
            crmaxlimtemp=0.99*crmaxlimtemp
         ENDIF
          IF(mype.eq.0.AND.lprint) THEN
            CALL  print_courants(1,cr3d,crxy,crxz,cryz,crx,cry,crz)
            PRINT "(1x,'         crmaxlim,crmaxlimtemp,dt,dt_ratio:',4(f11.3,1X))",crmaxlim,crmaxlimtemp,dtn,dt_ratio
          ENDIF
                 ENDDO adaptcrloop
      IF (mype == 0) THEN
        PRINT *, "EULAG timestep with adaptivity:"
        PRINT "(1x,'dt,dt_ratio: ',2(e11.3,1X))",dtn,dt_ratio
      ENDIF

    ELSE !not time apaptivity nor substepping

      IF ((.NOT.lcrdiag_nphlf)) THEN
!8a. We need to recompute prediction of advecting velocities, 
!   now we know the true dt_ratio between the previous and current timestep
       CALL velprd_traj0(ox,oy,oz,uadv,vadv,wadv,h,dxi,dyi,dzi,dtn,dt_ratio,  &
                      ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
      ELSE
!or it is enough to compute C-grid advectors, since A-grid advectors are unchanged
        CALL velprdA_to_C(ox(:,:,:,1),oy(:,:,:,1),oz(:,:,:,1),uadv,vadv,wadv, &
               ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih )
        IF (mype == 0) THEN
        PRINT *, "EULAG Courant number (no substepping nor adaptivity):"
!       PRINT "(1x,'cour,lipsh:',2(e23.16,2X))",cr_for_dt,cr2
        PRINT "(1x,'cour:',1(e23.16,2X))",cr_for_dt
        ENDIF
     ENDIF ! diagnostics at n+1/2
  ENDIF ! substepping or time adaptivitiy

  IF(lvertsplit) THEN
          FullXYZDomainWithHaloLoopDC( wadv(i,j,k)=0.5*wadv(i,j,k);) 
  ENDIF

FullXYZDomainLoopDC(
  u(i,j,k,1)=u(i,j,k,0)                          &
           +dt_ratio*(u(i,j,k,0)-u(i,j,k,2));
  v(i,j,k,1)=v(i,j,k,0)                          &
           +dt_ratio*(v(i,j,k,0)-v(i,j,k,2));
  w(i,j,k,1)=w(i,j,k,0)                          &
           +dt_ratio*(w(i,j,k,0)-w(i,j,k,2));
 
  u(i,j,k,2) =  u(i,j,k,0);
  v(i,j,k,2) =  v(i,j,k,0);
  w(i,j,k,2) =  w(i,j,k,0);

!-------- Store ox to be later used later as n-1; 
  ox(i,j,k,2)=ox(i,j,k,0)*h(i,j,k);
  oy(i,j,k,2)=oy(i,j,k,0)*h(i,j,k);
  oz(i,j,k,2)=oz(i,j,k,0)*h(i,j,k);
)
!---------------------------------------------
! IF (nsubsteps >1 ) THEN
!    subfac=1._euwp/nsubsteps
!Account for modification of the non-dimensional Courant number gc1,gc2,gc3 when smaller
!timestep   
!   uadv(:,:,:)  =subfac*uadv(:,:,:)
!   vadv(:,:,:)  =subfac*vadv(:,:,:)
!   wadv(:,:,:)  =subfac*wadv(:,:,:)
!   cr1=cr1*subfac
!   cr2=cr2*subfac
! ENDIF
  IF (mype == 0.AND.nsubsteps >1) THEN
    PRINT *, "----------------------------------------------------"
    PRINT *, "          EULAG Courant numbers(inside substepping):"
    PRINT "(1x,'          cour:',1(e23.16,2X))",cr_for_dt
!   PRINT *, "EULAG Courant and Lipschitz numbers(inside substepping):"
!   PRINT "(1x,'cour,lipsh:',2(e23.16,2X))",cr1,cr2
!   PRINT *, "----------------------------------------------------"
    PRINT "(1x,'           Dynamics is substepped times: ',1(i3,2X), "//&
            "'with timestep',1(f8.3,2X))",nsubsteps,dtn
!   PRINT *, "          Dynamic is substepped times: ",nsubsteps," with timestep", dtn
    PRINT *, "----------------------------------------------------"
  ENDIF

    END SUBROUTINE velprd_driver


END MODULE velprd_driver_module 

#include  "../advection/src_algorithms/renames.inc"
#include "../advection/src_algorithms/defines.inc"
#if (STATICMEM == 1)
MODULE eulag_moist_datafields
   USE precisions
   USE mod_parameters, ONLY: n,m,l,np,mp,lp,ih
   IMPLICIT NONE
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
  qv,qc,qr,qi,qs,qg
  REAL(KIND=euwp),DIMENSION(np,mp,lp) :: &
  qcrs
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
  qvref,qcref,qiref
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
  qv_eulag_old
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
  fqv_expl,fqc_expl,fqr_expl,fqi_expl
   REAL(KIND=euwp)  qv_a_lr(mp,lp, 2), qv_a_bt(np,lp, 2)
   REAL(KIND=euwp)  qc_a_lr(mp,lp, 2), qc_a_bt(np,lp, 2)
   REAL(KIND=euwp)  qr_a_lr(mp,lp, 2), qr_a_bt(np,lp, 2)
   REAL(KIND=euwp)  qi_a_lr(mp,lp, 2), qi_a_bt(np,lp, 2)
   REAL(KIND=euwp)  qs_a_lr(mp,lp, 2), qs_a_bt(np,lp, 2)
   REAL(KIND=euwp)  qg_a_lr(mp,lp, 2), qg_a_bt(np,lp, 2)
END MODULE eulag_moist_datafields
#else /*STATICMEM*/ 
MODULE eulag_moist_datafields
   USE precisions
  REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
  qv,qc,qr,qi,qs,qg,qcrs
   REAL_euwp, DIMENSION(:,:,:), ALLOCATABLE :: &
   qv_a_lr, qv_a_bt, qc_a_lr, qc_a_bt, qr_a_lr, qr_a_bt, &
   qi_a_lr, qi_a_bt, qs_a_lr, qs_a_bt, qg_a_lr, qg_a_bt
  REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
  qvref,qcref,qiref
  REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
  qv_eulag_old
  REAL_euwp,DIMENSION(:,:,:), ALLOCATABLE :: &
  fqv_expl,fqc_expl,fqr_expl,fqi_expl
CONTAINS
   SUBROUTINE allocate_eulag_moist_datafields
   USE mod_parameters, ONLY: np,mp,lp,ih
   USE, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_signaling_NAN
   IMPLICIT NONE
   REAL(KIND=euwp) real_initial_value
   INTEGER ierr,ierrtot
     ierrtot=0
    ALLOCATE ( qv(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qc(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qr(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qi(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qs(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qg(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qvref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qcref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qiref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qv_eulag_old(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( fqv_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( fqc_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( fqr_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( fqi_expl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),STAT=ierr ) ;ierrtot=ierrtot+ierr;

    ALLOCATE ( qv_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qc_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qr_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qi_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qs_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qg_a_lr(1:mp,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;

    ALLOCATE ( qv_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qc_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qr_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qi_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qs_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE ( qg_a_bt(1:np,1:lp,2),STAT=ierr ) ;ierrtot=ierrtot+ierr;
    ALLOCATE (qcrs(1:np,1:mp,1:lp),STAT=ierr ) ;ierrtot=ierrtot+ierr;

   IF(ierrtot/=0) print *,'moist allocate status: total number of errors is: ',ierrtot
#ifdef INITIALIZEVARIABLES
  real_initial_value = IEEE_VALUE(real_initial_value, IEEE_signaling_NAN)
    qv=real_initial_value
    qc=real_initial_value
    qr=real_initial_value
    qi=real_initial_value
    qs=real_initial_value
    qg=real_initial_value
    qvref=real_initial_value
    qcref=real_initial_value
    qiref=real_initial_value
    qcrs=real_initial_value
    qv_eulag_old=real_initial_value
    fqv_expl=real_initial_value
    fqc_expl=real_initial_value
    fqr_expl=real_initial_value
    fqi_expl=real_initial_value
    qv_a_lr=real_initial_value
    qc_a_lr=real_initial_value
    qr_a_lr=real_initial_value
    qi_a_lr=real_initial_value
    qg_a_lr=real_initial_value
    qs_a_lr=real_initial_value
    qv_a_bt=real_initial_value
    qc_a_bt=real_initial_value
    qr_a_bt=real_initial_value
    qi_a_bt=real_initial_value
    qg_a_bt=real_initial_value
    qs_a_bt=real_initial_value
#endif
   END SUBROUTINE allocate_eulag_moist_datafields
   SUBROUTINE deallocate_eulag_moist_datafields
     INTEGER ierr,ierrtot
     ierrtot=0
    DEALLOCATE ( qv,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qc,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qr,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qi,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qs,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qg,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qvref,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qcref,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qiref,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qv_eulag_old,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( fqv_expl,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( fqc_expl,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( fqr_expl,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    IF(ALLOCATED(fqi_expl)) DEALLOCATE ( fqi_expl,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qcrs,STAT=ierr ) ;ierrtot=ierrtot+ierr;

    DEALLOCATE ( qv_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qc_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qr_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qi_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qs_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qg_a_lr,STAT=ierr ) ;ierrtot=ierrtot+ierr;

    DEALLOCATE ( qv_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qc_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qr_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qi_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qs_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr;
    DEALLOCATE ( qg_a_bt,STAT=ierr ) ;ierrtot=ierrtot+ierr;
IF(ierrtot/=0) print *,'moist deallocate status: total number of errors is: ',ierrtot
   END SUBROUTINE deallocate_eulag_moist_datafields
END MODULE eulag_moist_datafields
#endif /*STATICMEM*/

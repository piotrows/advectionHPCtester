MODULE filters
USE precisions, ONLY: euwp,iintegers    
USE mpi_parallel, ONLY: leftedge,rightedge,botedge,topedge,gndedge,skyedge 
USE mpi_parallel, ONLY: updatelr,updatebt,update
IMPLICIT NONE
CONTAINS
#include "../advection/src_algorithms/defines.inc"
SUBROUTINE filtprf(a,ifl1,ifl2,np,mp,lp,ih)
        INTEGER(KIND=iintegers) :: np,mp,lp,ih,ifl1,ifl2,i,j,k
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih+1) :: a
        REAL(KIND=euwp),DIMENSION(1:lp) :: sz
  IF (ifl1 == 1) THEN
   DO j=1,mp
     DO i=1,np
       DO k=1+gndedge,lp-skyedge
         sz(k)=0.25_euwp*(a(i,j,k+1)+2._euwp*a(i,j,k)+a(i,j,k-1))
       END DO
       IF(gndedge == 1)  sz( 1)=a(i,j,1)
       IF(skyedge == 1)  sz(lp)=a(i,j,lp)
       DO k=1,lp
         a(i,j,k)=sz(k)
       END DO
     END DO
   END DO
  ENDIF
  IF (ifl2 == 1) THEN
    DO j=1,mp
      DO i=1,np
        DO k=1+2*gndedge,lp-2*skyedge
          sz(k)=0.25_euwp*(a(i,j,k+2)+2._euwp*a(i,j,k)+a(i,j,k-2))
        END DO
        IF(gndedge == 1)   sz(1 )=a(i,j,1)
        IF(gndedge == 1)   sz(2 )=a(i,j,2)
        IF(skyedge == 1)   sz(lp)=a(i,j,lp)
        IF(skyedge == 1) sz(lp-1)=a(i,j,lp-1)
        DO k=1,lp
          a(i,j,k)=sz(k)
       END DO
     END DO
   END DO
  ENDIF
  END SUBROUTINE filtprf
  !---------------------------------------------------------------------!
SUBROUTINE filstrP(a,nullfl,np,mp,lp,ih)
!---------------------------------------------------------------------!
  USE parameters, ONLY: &
    ibcx,ibcy,ibcz

  ! Subroutine arguments
  ! Local variables
  INTEGER(KIND=iintegers) :: np,mp,lp,ih
  REAL_euwp, INTENT(INOUT) :: a(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)


  REAL(KIND=euwp),DIMENSION(np) :: sx
  REAL(KIND=euwp),DIMENSION(mp) :: sy
  REAL(KIND=euwp),DIMENSION(lp) :: sz

  INTEGER(KIND=iintegers) :: i, j, k, nullfl

  INTEGER(KIND=iintegers) :: ierr
  INTEGER(KIND=iintegers) :: izerr

!-------------------------------------------------------------------------------
! Begin subroutine filstrP

  ierr = 0
  izerr = 0

  IF (nullfl == 0) THEN
    ierr = 0_iintegers
    RETURN
  ENDIF

  CALL updatelr(a,np,mp,lp,np,mp,lp,ih,ih) !izerr) ; ierr = ierr+izerr

  DO k=1,lp
    DO j=1,mp
      DO i=1+leftedge,np-rightedge
        sx(i)=0.25_euwp*(a(i+1,j,k)+2._euwp*a(i,j,k)+a(i-1,j,k))
      ENDDO
      IF (ibcx == 0) THEN
        IF (leftedge == 1) THEN
          sx(1)=0.75_euwp*a(1,j,k)+0.5_euwp*a(2,j,k)-0.25_euwp*a(3,j,k)
        END IF
        IF (rightedge == 1) THEN
          sx(np)=0.75_euwp*a(np,j,k)+0.5_euwp*a(np-1,j,k)-0.25_euwp*a(np-2,j,k)
        END IF
      ELSE IF(ibcx == 1) THEN
        IF (leftedge == 1) THEN
          sx(1)=0.25_euwp*(a(2,j,k)+2._euwp*a(1,j,k)+a(-1,j,k))
        END IF
        IF (rightedge == 1) THEN
          sx(np)=0.25_euwp*(a(np+2,j,k)+2._euwp*a(np+1,j,k)+a(np-1,j,k))
        END IF
      ENDIF !ibcx
      DO i=1,np
        a(i,j,k)=sx(i)
      ENDDO
    ENDDO
  ENDDO


    CALL updatebt(a,np,mp,lp,np,mp,lp,ih,ih) !izerr); ierr = ierr+izerr


    DO k=1,lp
      DO i=1,np
        DO j=1+botedge, mp-topedge
          sy(j)=0.25_euwp*(a(i,j+1,k)+2._euwp*a(i,j,k)+a(i,j-1,k))
        ENDDO
        IF (ibcy == 0) THEN
          IF (botedge == 1) THEN
            sy(1)=(0.75_euwp*a(i,1,k)+0.5_euwp*a(i,2,k)-0.25_euwp*a(i,3,k))
          END IF
          IF (topedge == 1) THEN
           sy(mp)=0.75_euwp*a(i,mp,k)+0.5_euwp*a(i,mp-1,k)-0.25_euwp*a(i,mp-2,k)
          END IF
        ELSE IF (ibcy == 1) THEN
          IF (botedge == 1) THEN
            sy(1)=0.25_euwp*(a(i,2,k)+2._euwp*a(i,1,k)+a(i,-1,k))
          END IF
          IF (topedge == 1) THEN
            sy(mp)=0.25_euwp*(a(i,mp+2,k)+2._euwp*a(i,mp+1,k)+a(i,mp-1,k))
          END IF
        ENDIF !ibcy
        DO j=1,mp
          a(i,j,k)=sy(j)
        ENDDO
      ENDDO
    ENDDO

!    CALL updategs(a,np,mp,lp,np,mp,lp,ih,ih) !izerr); ierr = ierr+izerr

  DO j=1,mp
    DO i=1,np
      DO k=1+gndedge,lp-skyedge
        sz(k)=0.25_euwp*(a(i,j,k+1)+2._euwp*a(i,j,k)+a(i,j,k-1))
      ENDDO
      IF (ibcz == 1) THEN
        IF (gndedge == 1) THEN
          sz(1)=0.25_euwp*(a(i,j,2)+2._euwp*a(i,j,1)+a(i,j,-1))
        ENDIF
        IF (skyedge == 1) THEN
          sz(lp)=0.25_euwp*(a(i,j,lp+2)+2._euwp*a(i,j,lp+1)+a(i,j,lp-1))
        ENDIF
      ELSE !Now ibcz=0
        IF (gndedge == 1) THEN
          sz(1)= (0.75_euwp*a(i,j,1)+0.5_euwp*a(i,j,2)-0.25_euwp*a(i,j,3))
        ENDIF
        IF (skyedge == 1) THEN
          sz(lp)=0.75_euwp*a(i,j,lp)+0.5_euwp*a(i,j,lp-1)-0.25_euwp*a(i,j,lp-2)
        ENDIF
      ENDIF
      DO k=1,lp
        a(i,j,k)=sz(k)
      ENDDO
    ENDDO
  ENDDO

  CALL update(a,np,mp,lp,np,mp,lp,ih,ih) !izerr) ;  ierr = ierr+izerr

END SUBROUTINE filstrP


END MODULE filters 

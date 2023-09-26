#include  "../advection/src_algorithms/renames.inc"
MODULE bconditions
   USE precisions
   USE parameters, ONLY: ibcx,ibcy,ibcz,ipoldiffmode 
   USE mpi_parallel, ONLY: leftedge,rightedge,botedge,topedge,gndedge,skyedge
   IMPLICIT NONE
CONTAINS
#include "../advection/src_algorithms/renames.inc"
SUBROUTINE halo_zero_z(x,np,mp,lp,ih)
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers) np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x

      !----------------
      ! C O N T E N T :
      !----------------
      IF (ibcz == 0) THEN
         IF (gndedge == 1) THEN
            x(1:np,1:mp,0   )=0._euwp
         ENDIF !gndedge
         IF (skyedge == 1) THEN
            x(1:np,1:mp,lp+1)=0._euwp
         ENDIF
      ELSE !ibcz
         IF (gndedge == 1) THEN
            x(1:np,1:mp,0   )= x(1:np,1:mp,  -1)
         ENDIF !gndedge
         IF (skyedge == 1) THEN
            x(1:np,1:mp,lp+1)= x(1:np,1:mp,lp+2)
         ENDIF
      ENDIF!ibcz
   END SUBROUTINE halo_zero_z


!---------------------------------------------------------------------!
SUBROUTINE halo_zero_y(x,np,mp,lp,ih)
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers) np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x

      !----------------
      ! C O N T E N T :
      !----------------
      IF (ibcy == 0) THEN
         IF (botedge == 1) THEN
            x(1:np,0   ,1:lp)=0._euwp
         ENDIF
         IF (topedge == 1) THEN
            x(1:np,mp+1,1:lp)=0._euwp
         ENDIF
      ELSE !ibcy
         IF (botedge == 1) THEN
            x(1:np,0   ,1:lp)= x(1:np,  -1,1:lp)
         ENDIF
         IF (topedge == 1) THEN
            x(1:np,mp+1,1:lp)= x(1:np,mp+2,1:lp)
         ENDIF
      ENDIF!ibcy
   END SUBROUTINE halo_zero_y


!---------------------------------------------------------------------!
! Extrapolate 3D field 'x' at X-edges
!---------------------------------------------------------------------!
SUBROUTINE halo_zero_x(x,np,mp,lp,ih)
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN):: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x

      !----------------
      ! C O N T E N T :
      !----------------
      IF (ibcx == 0) THEN
          IF (leftedge == 1) THEN
             x(   0,1:mp,1:lp)=0._euwp
          ENDIF
          IF (rightedge == 1) THEN
             x(np+1,1:mp,1:lp)=0._euwp
          ENDIF
      ELSE !ibcx
         IF (leftedge == 1) THEN
            x(   0,1:mp,1:lp)= x(  -1,1:mp,1:lp) 
         ENDIF
         IF (rightedge == 1) THEN
            x(np+1,1:mp,1:lp)= x(np+2,1:mp,1:lp)
         ENDIF
      ENDIF!ibcx
   END SUBROUTINE halo_zero_x

!---------------------------------------------------------------------!
!
!---------------------------------------------------------------------!
SUBROUTINE halo_one_sided_diff_z(x,np,mp,lp,ih)
      IMPLICIT NONE
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers) np,mp,lp,ih
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x

      !----------------
      ! C O N T E N T :
      !----------------
      IF (ibcz == 0) THEN
         IF (gndedge == 1) THEN
            x(1:np,1:mp,0   )=-x(1:np,1:mp,   2)+2._euwp*x(1:np,1:mp, 1)
         ENDIF !gndedge
         IF (skyedge == 1) THEN
            x(1:np,1:mp,lp+1)=-x(1:np,1:mp,lp-1)+2._euwp*x(1:np,1:mp,lp)
         ENDIF
      ELSE !ibcz
         IF (gndedge == 1) THEN
            x(1:np,1:mp,0   )= x(1:np,1:mp,  -1)
         ENDIF !gndedge
         IF (skyedge == 1) THEN
            x(1:np,1:mp,lp+1)= x(1:np,1:mp,lp+2)
         ENDIF
      ENDIF!ibcz
   END SUBROUTINE halo_one_sided_diff_z


!---------------------------------------------------------------------!
! Extrapolate 3D field 'x' at Y-edges
!---------------------------------------------------------------------!
SUBROUTINE halo_one_sided_diff_y(x,iflip,ipoles,np,mp,lp,ih)
      IMPLICIT NONE
!---------------------------------------------------------------------!
      INTEGER(KIND=iintegers),INTENT(IN) ::  ipoles,iflip,np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x

      !----------------
      ! C O N T E N T :
      !----------------
    IF(ipoles == 0) THEN
      IF (ibcy == 0) THEN
         IF (botedge == 1) THEN
            x(1:np,0   ,1:lp)=-x(1:np,2   ,1:lp)+2._euwp*x(1:np,1 ,1:lp)
         ENDIF
         IF (topedge == 1) THEN
            x(1:np,mp+1,1:lp)=-x(1:np,mp-1,1:lp)+2._euwp*x(1:np,mp,1:lp)
         ENDIF
      ELSE !ibcy
         IF (botedge == 1) THEN
            x(1:np,0   ,1:lp)= x(1:np,  -1,1:lp)
         ENDIF
         IF (topedge == 1) THEN
            x(1:np,mp+1,1:lp)= x(1:np,mp+2,1:lp)
         ENDIF
      ENDIF!ibcy
    ELSE
      if(iflip.eq.1) then
      if(ipoldiffmode.eq.1) then
       if(botedge.eq.1) then
          x(1:np,0   ,1:lp)=-x(1:np,1 ,1:lp)
       endif
       if(topedge.eq.1) then
          x(1:np,mp+1,1:lp)=-x(1:np,mp,1:lp)
       endif
      else !ipoldiffmode
       if(botedge.eq.1) then
          x(1:np,0   ,1:lp)=-x(1:np,0 ,1:lp)
       endif
       if(topedge.eq.1) then
          x(1:np,mp+1,1:lp)=-x(1:np,mp+1,1:lp)
       endif
      endif !ipoldiffmode
       endif !iflip

    ENDIF

   END SUBROUTINE halo_one_sided_diff_y


!---------------------------------------------------------------------!
! Extrapolate 3D field 'x' at X-edges
!---------------------------------------------------------------------!
SUBROUTINE halo_one_sided_diff_x(x,np,mp,lp,ih)
      IMPLICIT NONE
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN):: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x

      !----------------
      ! C O N T E N T :
      !----------------
      IF (ibcx == 0) THEN
          IF (leftedge == 1) THEN
             x(   0,1:mp,1:lp)=-x(2   ,1:mp,1:lp)+2._euwp*x(1 ,1:mp,1:lp)
          ENDIF
          IF (rightedge == 1) THEN
             x(np+1,1:mp,1:lp)=-x(np-1,1:mp,1:lp)+2._euwp*x(np,1:mp,1:lp)
          ENDIF
      ELSE !ibcx
         IF (leftedge == 1) THEN
            x(   0,1:mp,1:lp)= x(  -1,1:mp,1:lp) 
         ENDIF
         IF (rightedge == 1) THEN
            x(np+1,1:mp,1:lp)= x(np+2,1:mp,1:lp)
         ENDIF
      ENDIF!ibcx
   END SUBROUTINE halo_one_sided_diff_x


!---------------------------------------------------------------------!
! Extrapolate 3D field 'x' at Y-edges
!---------------------------------------------------------------------!
SUBROUTINE halo_one_sided_diff_prec_y(x,imode,ipoles,np,mp,lp,ih)
      IMPLICIT NONE
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN) ::  imode,ipoles,np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x
    IF(ipoles == 0) THEN
      IF (ibcy == 0) THEN
        IF(imode == 0) THEN
          IF (botedge == 1) THEN
             x(1:np,0   ,1:lp)=x(1:np,2   ,1:lp)
          ENDIF
          IF (topedge == 1) THEN
             x(1:np,mp+1,1:lp)=x(1:np,mp-1,1:lp)
          ENDIF
        ELSE
          IF (botedge == 1) THEN
             x(1:np,0   ,1:lp)=-x(1:np,2   ,1:lp)
          ENDIF
          IF (topedge == 1) THEN
             x(1:np,mp+1,1:lp)=-x(1:np,mp-1,1:lp)
          ENDIF
        ENDIF
      ELSE !ibcy
         IF (botedge == 1) THEN
            x(1:np,0   ,1:lp)= x(1:np,  -1,1:lp)
         ENDIF
         IF (topedge == 1) THEN
            x(1:np,mp+1,1:lp)= x(1:np,mp+2,1:lp)
         ENDIF
      ENDIF!ibcy
    ELSE
      if(imode.eq.1) then
        if(ipoldiffmode.eq.1) then
          if(botedge.eq.1) then
            x(1:np,0   ,1:lp)=-x(1:np,1 ,1:lp)
          endif
          if(topedge.eq.1) then
            x(1:np,mp+1,1:lp)=-x(1:np,mp,1:lp)
          endif
        else
          if(botedge.eq.1) then
            x(1:np,0   ,1:lp)=-x(1:np,0 ,1:lp)
          endif
          if(topedge.eq.1) then
            x(1:np,mp+1,1:lp)=-x(1:np,mp+1,1:lp)
          endif
        endif
      endif
    ENDIF
   END SUBROUTINE halo_one_sided_diff_prec_y


!---------------------------------------------------------------------!
! Extrapolate 3D field 'x' at X-edges
!---------------------------------------------------------------------!
SUBROUTINE halo_one_sided_diff_prec_x(x,imode,np,mp,lp,ih)
      IMPLICIT NONE
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN):: imode,np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x

      !----------------
      ! C O N T E N T :
      !----------------
      IF (ibcx == 0) THEN
        IF(imode == 0) THEN
          IF (leftedge == 1) THEN
             x(   0,1:mp,1:lp)=x(2   ,1:mp,1:lp)
          ENDIF
          IF (rightedge == 1) THEN
             x(np+1,1:mp,1:lp)=x(np-1,1:mp,1:lp)
          ENDIF
        ELSE
          IF (leftedge == 1) THEN
             x(   0,1:mp,1:lp)=-x(2   ,1:mp,1:lp)
          ENDIF
          IF (rightedge == 1) THEN
             x(np+1,1:mp,1:lp)=-x(np-1,1:mp,1:lp)
          ENDIF
        ENDIF
      ELSE !ibcx
         IF (leftedge == 1) THEN
            x(   0,1:mp,1:lp)= x(  -1,1:mp,1:lp)
         ENDIF
         IF (rightedge == 1) THEN
            x(np+1,1:mp,1:lp)= x(np+2,1:mp,1:lp)
         ENDIF
      ENDIF!ibcx
   END SUBROUTINE halo_one_sided_diff_prec_x
!---------------------------------------------------------------------!
! Extrapolate 2D field 'x' at Y-edges
!---------------------------------------------------------------------!
SUBROUTINE halo_one_sided_diff_flt_y(x,ipoles,np,mp,ih)
      IMPLICIT NONE
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers) ipoles,np,mp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih) :: x

      !----------------
      ! C O N T E N T :
      !----------------
    IF(ipoles == 0) THEN
      IF (ibcy == 0) THEN
         IF (botedge == 1) THEN
            x(1:np,0   )=-x(1:np,2   )+2*x(1:np,1 )
         ENDIF
         IF (topedge == 1) THEN
            x(1:np,mp+1)=-x(1:np,mp-1)+2*x(1:np,mp)
         ENDIF
      ELSE !ibcy
         IF (botedge == 1) THEN
            x(1:np,0   )= x(1:np,  -1) +x(1:np,1)-x(1:np,0)
         ENDIF
         IF (topedge == 1) THEN
            x(1:np,mp+1)= x(1:np,mp+2)-x(1:np,mp+1)+x(1:np,mp)
         ENDIF
      ENDIF!ibcy
    ELSE  !ipoles
    ENDIF !ipoles
   END SUBROUTINE halo_one_sided_diff_flt_y


!---------------------------------------------------------------------!
! Extrapolate 2D field 'x' at X-edges
!---------------------------------------------------------------------!
SUBROUTINE halo_one_sided_diff_flt_x(x,ipoles,np,mp,ih)
      IMPLICIT NONE
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers) ipoles,np,mp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih) :: x
    IF(ipoles == 0) THEN
      IF (ibcx == 0) THEN
         IF (leftedge == 1) THEN
      x(   0,1:mp)=-x(2   ,1:mp)+2._euwp*x( 1,1:mp)
         ENDIF
         IF (rightedge == 1) THEN
      x(np+1,1:mp)=-x(np-1,1:mp)+2._euwp*x(np,1:mp)
         ENDIF
      ELSE !ibcx
         IF (leftedge == 1) THEN
            x(   0,1:mp)= x(  -1,1:mp)+x(1   ,1:mp)-x(0   ,1:mp)
         ENDIF
         IF (rightedge == 1) THEN
            x(np+1,1:mp)= x(np+2,1:mp)-x(np+1,1:mp)+x(np  ,1:mp)
         ENDIF
      ENDIF!ibcx
    ENDIF
   END SUBROUTINE halo_one_sided_diff_flt_x

SUBROUTINE halo_one_sided_diff_flux_z(x,y,np,mp,lp,ih)
      IMPLICIT NONE
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers) np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x,y

      !----------------
      ! C O N T E N T :
      !----------------
      IF (ibcz == 0) THEN
         IF (gndedge == 1) THEN
            x(1:np,1:mp,0   )=2._euwp*x(1:np,1:mp,1)*y(1:np,1:mp,1)   &
                             -        x(1:np,1:mp,2)*y(1:np,1:mp,2)
                                   
            y(1:np,1:mp,0   )=1._euwp
         ENDIF !gndedge
         IF (skyedge == 1) THEN
            x(1:np,1:mp,lp+1)=2._euwp*x(1:np,1:mp,lp  )*y(1:np,1:mp,lp  ) &
                             -        x(1:np,1:mp,lp-1)*y(1:np,1:mp,lp-1)
            y(1:np,1:mp,lp+1)=1._euwp
         ENDIF
      ELSE !ibcz
         IF (gndedge == 1) THEN
            x(1:np,1:mp,0   )= x(1:np,1:mp,  -1)
            y(1:np,1:mp,0   )= y(1:np,1:mp,  -1)
         ENDIF !gndedge
         IF (skyedge == 1) THEN
            x(1:np,1:mp,lp+1)= x(1:np,1:mp,lp+2)
            y(1:np,1:mp,lp+1)= y(1:np,1:mp,lp+2)
         ENDIF
      ENDIF!ibcz
   END SUBROUTINE halo_one_sided_diff_flux_z


!---------------------------------------------------------------------!
! Extrapolate 3D field 'x' at Y-edges
!---------------------------------------------------------------------!
SUBROUTINE halo_one_sided_diff_flux_y(x,y,iflip,ipoles,np,mp,lp,ih)
      IMPLICIT NONE
!---------------------------------------------------------------------!
INTEGER(KIND=iintegers) iflip,ipoles,np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x,y

      !----------------
      ! C O N T E N T :
      !----------------
    IF(ipoles == 0) THEN
      IF (ibcy == 0) THEN
         IF (botedge == 1) THEN
            x(1:np,0   ,1:lp)=2._euwp*x(1:np,1   ,1:lp)*y(1:np,1   ,1:lp) &
                                   -  x(1:np,2   ,1:lp)*y(1:np,2   ,1:lp)
            y(1:np,0   ,1:lp)=1._euwp
         ENDIF
         IF (topedge == 1) THEN
            x(1:np,mp+1,1:lp)=2._euwp*x(1:np,mp  ,1:lp)*y(1:np,mp  ,1:lp) &
                                     -x(1:np,mp-1,1:lp)*y(1:np,mp-1,1:lp)

            y(1:np,mp+1,1:lp)=1._euwp
         ENDIF
      ELSE !ibcy
         IF (botedge == 1) THEN
            x(1:np,0   ,1:lp)= x(1:np,  -1,1:lp)
            y(1:np,0   ,1:lp)= y(1:np,  -1,1:lp)
         ENDIF
         IF (topedge == 1) THEN
            x(1:np,mp+1,1:lp)= x(1:np,mp+2,1:lp)
            y(1:np,mp+1,1:lp)= y(1:np,mp+2,1:lp)
         ENDIF
      ENDIF!ibcy
    ELSE
!Flip
      if(iflip.eq.1) then
       if(ipoldiffmode.eq.1) then
         if(botedge.eq.1) then
            x(1:np,0   ,1:lp)=-x(1:np,1 ,1:lp)
            y(1:np,0   ,1:lp)= y(1:np,1 ,1:lp)
         endif
         if(topedge.eq.1) then
            x(1:np,mp+1,1:lp)=-x(1:np,mp,1:lp)
            y(1:np,mp+1,1:lp)= y(1:np,mp,1:lp)
         endif
       else
         if(botedge.eq.1) then
            x(1:np,0   ,1:lp)=-x(1:np,0   ,1:lp)
         endif
         if(topedge.eq.1) then
            x(1:np,mp+1,1:lp)=-x(1:np,mp+1,1:lp)
         endif
       endif !ipoldiffmode
      endif
    ENDIF
   END SUBROUTINE halo_one_sided_diff_flux_y


!---------------------------------------------------------------------!
! Extrapolate 3D field 'x' at X-edges
!---------------------------------------------------------------------!
SUBROUTINE halo_one_sided_diff_flux_x(x,y,np,mp,lp,ih)
      IMPLICIT NONE
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN):: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: x,y

      !----------------
      ! C O N T E N T :
      !----------------
      IF (ibcx == 0) THEN
         IF (leftedge == 1) THEN
            x(   0,1:mp,1:lp)=2._euwp*x(1  ,1:mp,1:lp)*y( 1,1:mp,1:lp) &
                                    - x(2  ,1:mp,1:lp)*y( 2,1:mp,1:lp)
            y(   0,1:mp,1:lp)= 1._euwp
         ENDIF
         IF (rightedge == 1) THEN
            x(np+1,1:mp,1:lp)=2._euwp*x(np  ,1:mp,1:lp)*y(np  ,1:mp,1:lp) &
                                    - x(np-1,1:mp,1:lp)*y(np-1,1:mp,1:lp)
            y(np+1,1:mp,1:lp)= 1._euwp
         ENDIF
      ELSE !ibcx
         IF (leftedge == 1) THEN
            x(   0,1:mp,1:lp)= x(  -1,1:mp,1:lp)
            y(   0,1:mp,1:lp)= y(  -1,1:mp,1:lp)
         ENDIF
         IF (rightedge == 1) THEN
            x(np+1,1:mp,1:lp)= x(np+2,1:mp,1:lp)
            y(np+1,1:mp,1:lp)= y(np+2,1:mp,1:lp)
         ENDIF
      ENDIF!ibcx
   END SUBROUTINE halo_one_sided_diff_flux_x


!---------------------------------------------------------------------!
! Copy edge of domain to halo
!---------------------------------------------------------------------!
   SUBROUTINE halo_copy_edge_to_halo(x,np,mp,lp,ih)
     INTEGER(KIND=iintegers) :: np,mp,lp,ih
     INTEGER(KIND=iintegers) :: j,k
     REAL (KIND=euwp) :: x(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
     IF (ih == 3) THEN
       IF(leftedge.eq.1) THEN
         x( 0,1:mp,1:lp) = x(1,1:mp,1:lp)
         x(-1,1:mp,1:lp) = x(1,1:mp,1:lp)
         x(-2,1:mp,1:lp) = x(1,1:mp,1:lp)
       ENDIF
       IF(rightedge.eq.1) THEN
         x(np+1,1:mp,1:lp) = x(np,1:mp,1:lp)
         x(np+2,1:mp,1:lp) = x(np,1:mp,1:lp)
         x(np+3,1:mp,1:lp) = x(np,1:mp,1:lp)
       ENDIF

       IF(botedge.eq.1) THEN
         x(1:np, 0,1:lp) = x(1:np,1,1:lp)
         x(1:np,-1,1:lp) = x(1:np,1,1:lp)
         x(1:np,-2,1:lp) = x(1:np,1,1:lp)
       ENDIF
       IF(topedge.eq.1) THEN
         x(1:np,mp+1,1:lp) = x(1:np,mp,1:lp)
         x(1:np,mp+2,1:lp) = x(1:np,mp,1:lp)
         x(1:np,mp+3,1:lp) = x(1:np,mp,1:lp)
       ENDIF

       IF(gndedge.eq.1) THEN
         x(1:np,1:mp, 0) = x(1:np,1:mp,1)
         x(1:np,1:mp,-1) = x(1:np,1:mp,1)
         x(1:np,1:mp,-2) = x(1:np,1:mp,1)
       ENDIF
       IF(skyedge.eq.1) THEN
         x(1:np,1:mp,lp+1) = x(1:np,1:mp,lp)
         x(1:np,1:mp,lp+2) = x(1:np,1:mp,lp)
         x(1:np,1:mp,lp+3) = x(1:np,1:mp,lp)
       ENDIF
     ELSE IF (ih == 2) THEN
       IF(leftedge.eq.1) THEN
         x( 0,1:mp,1:lp) = x(1,1:mp,1:lp)
         x(-1,1:mp,1:lp) = x(1,1:mp,1:lp)
       ENDIF
       IF(rightedge.eq.1) THEN
         x(np+1,1:mp,1:lp) = x(np,1:mp,1:lp)
         x(np+2,1:mp,1:lp) = x(np,1:mp,1:lp)
       ENDIF

       IF(botedge.eq.1) THEN
         x(1:np, 0,1:lp) = x(1:np,1,1:lp)
         x(1:np,-1,1:lp) = x(1:np,1,1:lp)
       ENDIF
       IF(topedge.eq.1) THEN
         x(1:np,mp+1,1:lp) = x(1:np,mp,1:lp)
         x(1:np,mp+2,1:lp) = x(1:np,mp,1:lp)
       ENDIF

       IF(gndedge.eq.1) THEN
         x(1:np,1:mp, 0) = x(1:np,1:mp,1)
         x(1:np,1:mp,-1) = x(1:np,1:mp,1)
       ENDIF
       IF(skyedge.eq.1) THEN
         x(1:np,1:mp,lp+1) = x(1:np,1:mp,lp)
         x(1:np,1:mp,lp+2) = x(1:np,1:mp,lp)
       ENDIF
     ELSE IF (ih == 1) THEN
       IF(leftedge.eq.1) THEN
         DO k=1,lp
           DO j=1,mp 
         x(0   ,j,k) = x(1 ,j,k)
           ENDDO
         ENDDO
       ENDIF
       IF(rightedge.eq.1) THEN
         x(np+1,1:mp,1:lp) = x(np,1:mp,1:lp)
       ENDIF

       IF(botedge.eq.1) THEN
         x(1:np,0   ,1:lp) = x(1:np,1 ,1:lp)
       ENDIF
       IF(topedge.eq.1) THEN
         x(1:np,mp+1,1:lp) = x(1:np,mp,1:lp)
       ENDIF

       IF(gndedge.eq.1) THEN
         x(1:np,1:mp,0   ) = x(1:np,1:mp,1 )
       ENDIF
       IF(skyedge.eq.1) THEN
         x(1:np,1:mp,lp+1) = x(1:np,1:mp,lp)
       ENDIF
      ENDIF

    END SUBROUTINE halo_copy_edge_to_halo

!---------------------------------------------------------------------!
! Copy edge of domain to halo
!---------------------------------------------------------------------!

    SUBROUTINE cp_bcs_from3D_xy(x,x_a_lr,x_a_bt,np,mp,lp,ih)
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL (KIND=euwp), INTENT(IN)  :: x(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
      REAL (KIND=euwp), INTENT(OUT) :: x_a_lr(mp,lp, 2), x_a_bt(np,lp, 2) 
      IF( leftedge == 1) x_a_lr(1:mp,1:lp,1)= x(1 ,1:mp,1:lp)
      IF(rightedge == 1) x_a_lr(1:mp,1:lp,2)= x(np,1:mp,1:lp)
      IF(  botedge == 1) x_a_bt(1:np,1:lp,1)= x(1:np,1,1:lp)
      IF(  topedge == 1) x_a_bt(1:np,1:lp,2)= x(1:np,mp,1:lp)
    END SUBROUTINE cp_bcs_from3D_xy 

    SUBROUTINE cp_bcs_from1D_xy(x1d,x_a_lr,x_a_bt,np,mp,lp)
      INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp
      REAL (KIND=euwp), INTENT(IN)  :: x1d(1:lp)
      REAL (KIND=euwp), INTENT(OUT) :: x_a_lr(mp,lp, 2), x_a_bt(np,lp, 2)     
      INTEGER k
      DO k=1,lp
        IF( leftedge == 1) x_a_lr(1:mp,k,1)= x1d(k)
        IF(rightedge == 1) x_a_lr(1:mp,k,2)= x1d(k)
        IF(  botedge == 1) x_a_bt(1:np,k,1)= x1d(k)
        IF(  topedge == 1) x_a_bt(1:np,k,2)= x1d(k)
      ENDDO
    END SUBROUTINE cp_bcs_from1D_xy      

    SUBROUTINE cp_bcs_from0D_xy(x0d,x_a_lr,x_a_bt,np,mp,lp)
      INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp
      REAL (KIND=euwp), INTENT(IN)  :: x0d
      REAL (KIND=euwp), INTENT(OUT) :: x_a_lr(mp,lp, 2), x_a_bt(np,lp, 2)
      INTEGER k
      DO k=1,lp
        IF( leftedge == 1) x_a_lr(1:mp,k,1)= x0d
        IF(rightedge == 1) x_a_lr(1:mp,k,2)= x0d
        IF(  botedge == 1) x_a_bt(1:np,k,1)= x0d
        IF(  topedge == 1) x_a_bt(1:np,k,2)= x0d
      ENDDO
    END SUBROUTINE cp_bcs_from0D_xy
    SUBROUTINE remove_cyclic_offset_full(x,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
      USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge,gndedge,skyedge
      INTEGER(KIND=iintegers),INTENT(IN) :: ibcx0,ibcy0,ibcz0,ipoles0,np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: x
      IF(ipoles0 == 0) THEN
        IF(ibcx0.eq.1) THEN
          IF ( leftedge == 1) x(0   ,0:mp+1,0:lp+1)=x(-1  ,0:mp+1,0:lp+1)
          IF (rightedge == 1) x(np+1,0:mp+1,0:lp+1)=x(np+2,0:mp+1,0:lp+1)
        ENDIF
        IF(ibcy0.eq.1) THEN
          IF ( botedge == 1) x(0:np+1,0     ,0:lp+1)=x(0:np+1,-1  ,0:lp+1)
          IF ( topedge == 1) x(0:np+1,mp+1  ,0:lp+1)=x(0:np+1,mp+2,0:lp+1)
        ENDIF
      ENDIF
      IF(ibcz0.eq.1) THEN
        IF ( gndedge == 1) x(0:np+1,0:mp+1,0     )=x(0:np+1,0:mp+1,-1  )
        IF ( skyedge == 1) x(0:np+1,0:mp+1,  lp+1)=x(0:np+1,0:mp+1,lp+2)
      ENDIF
    END SUBROUTINE remove_cyclic_offset_full
    SUBROUTINE cp_scnd_last_to_halo_xyz(x,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge,gndedge,skyedge
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers),INTENT(IN) :: ibcx0,ibcy0,ibcz0,np,mp,lp,ih,ipoles0
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x

        IF(ibcz0.eq.0) THEN
          IF (gndedge == 1) x(1:np,1:mp,0   )=x(1:np,1:mp,2   )
          IF (skyedge == 1) x(1:np,1:mp,lp+1)=x(1:np,1:mp,lp-1)
        ENDIF
        IF(ipoles0.eq.0) THEN
          IF(ibcy0.eq.0) THEN
            IF (botedge == 1) x(1:np,0   ,1:lp)=x(1:np,2   ,1:lp)
            IF (topedge == 1) x(1:np,mp+1,1:lp)=x(1:np,mp-1,1:lp)
          ENDIF
          IF(ibcx0.eq.0) THEN
            IF (leftedge  == 1)  x(0   ,1:mp,1:lp)=x(2   ,1:mp,1:lp)
            IF (rightedge == 1)  x(np+1,1:mp,1:lp)=x(np-1,1:mp,1:lp)
          ENDIF
        ENDIF !ipoles
    END SUBROUTINE cp_scnd_last_to_halo_xyz
    SUBROUTINE cp_scnd_last_to_halo_xyz_full(x,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge,gndedge,skyedge
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers),INTENT(IN) :: ibcx0,ibcy0,ibcz0,ipoles0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x
        IF(ibcz0.eq.0) THEN
          IF (gndedge == 1)   x(0:np+1,0:mp+1,0     )=x(0:np+1,0:mp+1,2     )
          IF (skyedge == 1)   x(0:np+1,0:mp+1,  lp+1)=x(0:np+1,0:mp+1,  lp-1)
        ENDIF
        IF(ipoles0.eq.0) THEN
          IF(ibcx0.eq.0) THEN
            IF (leftedge == 1)  x(0     ,0:mp+1,0:lp+1)=x(2     ,0:mp+1,0:lp+1)
            IF (rightedge == 1) x(np+1  ,0:mp+1,0:lp+1)=x(np-1  ,0:mp+1,0:lp+1)
          ENDIF
          IF(ibcy0.eq.0) THEN
            IF (botedge == 1)   x(0:np+1,0     ,0:lp+1)=x(0:np+1,2     ,0:lp+1)
            IF (topedge == 1)   x(0:np+1,mp+1  ,0:lp+1)=x(0:np+1,  mp-1,0:lp+1)
          ENDIF
        ENDIF
    END SUBROUTINE cp_scnd_last_to_halo_xyz_full
END MODULE bconditions

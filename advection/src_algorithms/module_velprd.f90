!
! Module contains velocity predictor scheme in the form of linear extrapolation
! in time velocity predictor scheme in the form of linear extrapolation in time.
! First order accuracy is sufficient to achieve fully second-order integration
! of the equations of motion. More elaborate predictor based on Runge-Kutta
! method will be added in the future (itraj=1). 
!

MODULE module_velprd 
    USE precisions, ONLY  : iintegers,euwp
    USE mpi_parallel, ONLY: leftedge,rightedge,botedge,topedge,gndedge,skyedge
    IMPLICIT NONE
CONTAINS
#include "defines.inc"
    SUBROUTINE velprd_traj0(ox,oy,oz,u1,u2,u3,h,dxi,dyi,dzi,dt,dt_ratio, &
               ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih,do_serialization_in,do_serialization_out )
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY: ttbeg,ttend
        USE mpi_parallel, ONLY: iupx,iupy,iupz 
        USE mpi_parallel, ONLY: update,updatelr,updatebt,updategs
        LOGICAL,OPTIONAL :: do_serialization_in
        LOGICAL,OPTIONAL :: do_serialization_out
        INTEGER(KIND=iintegers), INTENT(IN) :: np,mp,lp,ih
        REAL_euwp, INTENT(INOUT) ::                            &
                      ox(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),     &
                      oy(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),     &
                      oz(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2)

        REAL_euwp, INTENT(OUT) ::  &
                          u1(1-ih:np+1+ih, 1-ih:mp+ih, 1-ih:lp+ih),    &
                          u2(1-ih:np+ih, 1-ih:mp+1+ih, 1-ih:lp+ih),    &
                          u3(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih)

        REAL_euwp, INTENT(IN) ::  &
                           h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
 
        REAL(KIND=euwp),INTENT(IN) :: dxi,dyi,dzi,dt,dt_ratio

        INTEGER(KIND=iintegers) :: i,j,k
        INTEGER(KIND=iintegers) :: ipoles0,ibcx0,ibcy0,ibcz0

        REAL(KIND=euwp) :: ox0 , ox1 , ox2,  oy0 , oy1 , oy2 , oz0 , oz1, oz2
        REAL(KIND=euwp) :: ox_i, oy_j, oz_k
        REAL(KIND=euwp) :: gc1,gc2,gc3 
  CALL ttbeg(40)
  CALL ttbeg(1040)
        IF(PRESENT(do_serialization_in).OR.PRESENT(do_serialization_out)) THEN
!Do nothing but remove compiler warning about unused dummy variables
        ENDIF

!------- Compute advecting velocities 
  gc1=dt*dxi
  gc2=dt*dyi
  gc3=dt*dzi
!  print *,'gc1,gc2,gc3',gc1,gc2,gc3
 !----------------------------------------------------------------------!
  ! Do linear extrapolation to time level n+1
  !----------------------------------------------------------------------!
  !-- x direction -----------
FullXYZDomainLoopDC(
        !First extrapolate the contravariant velocities to n+1;
        !Contravariant velocities(momentum) at level n:;
        ox0 =ox(i,j,k,0)*h(i,j,k);
        !Contravariant velocities(momentum) at level n-1:;
        ox2 =ox(i,j,k,2);
        !Construct n+1 contravariant velocities(momentum);
        ox1 =ox0+dt_ratio*(ox0-ox2);
        !Construct and store n+1/2 contravariant velocities(momentum);
        ox_i  =0.5_euwp*(ox1+ox0 )*gc1;
        ox(i,j,k,1)=ox_i;

        oy0   =oy(i,j  ,k,0)*h(i,j  ,k);
        oy2   =oy(i,j  ,k,2);
        oy1   =oy0+dt_ratio*(oy0    -oy2);
        oy_j  =0.5_euwp*(oy1 +oy0 )*gc2;
        oy(i,j,k,1)=oy_j;

        oz0 =oz(i,j,k  ,0)*h(i,j,k  );
        oz2 =oz(i,j,k  ,2);
        oz1 =oz0+dt_ratio*(oz0-oz2);
        oz_k  =0.5_euwp*(oz1 +oz0 )*gc3;
        oz(i,j,k,1)=oz_k;
)
    CALL ttend(1040)
    CALL velprdA_to_C(ox(:,:,:,1),oy(:,:,:,1),oz(:,:,:,1),u1,u2,u3, &
               ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih,do_serialization_in,do_serialization_out )

   CALL ttend(40)
    END SUBROUTINE velprd_traj0 

    SUBROUTINE velprd_traj0_adapt_dt(ox,oy,oz,h,dxi,dyi,dzi,dt,dt_ratio, &
               np,mp,lp,ih,do_serialization_in,do_serialization_out )
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY: ttbeg,ttend
        USE mpi_parallel, ONLY: iupx,iupy,iupz 
        USE mpi_parallel, ONLY: update,updatelr,updatebt,updategs
        LOGICAL,OPTIONAL :: do_serialization_in
        LOGICAL,OPTIONAL :: do_serialization_out
        INTEGER(KIND=iintegers), INTENT(IN) :: np,mp,lp,ih
        REAL_euwp, INTENT(INOUT) ::                            &
                      ox(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),     &
                      oy(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),     &
                      oz(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2)
        REAL_euwp, INTENT(IN) ::  &
                           h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
 
        REAL(KIND=euwp),INTENT(IN) :: dxi,dyi,dzi,dt,dt_ratio

        INTEGER(KIND=iintegers) :: i,j,k

        REAL(KIND=euwp) :: ox0 , ox1 , ox2,  oy0 , oy1 , oy2 , oz0 , oz1, oz2
        REAL(KIND=euwp) :: ox_i, oy_j, oz_k
        REAL(KIND=euwp) :: gc1,gc2,gc3 
        REAL(KIND=euwp) :: cr3d,crxy,crxz,cryz,crx,cry,crz

  CALL ttbeg(40)
  CALL ttbeg(1040)
        IF(PRESENT(do_serialization_in).OR.PRESENT(do_serialization_out)) THEN
!Do nothing but remove compiler warning about unused dummy variables
        ENDIF
! dtn=dt;
! CALL compute_courants(cr3d,crxy,crxz,cryz,crx,cry,crz,ox,oy,oz,h, &
!                       1,dxi,dyi,dzi,dt,np,mp,lp,ih);
! dtnplus = min(cour_max_allowed/(abs(cr3d)+TINY(cr3d))*dtn,dtmax) ! proposed next timestep;
! spex=spex*dtnplus_dtn;
! spexi=1./spex;

 
!------- Compute advecting velocities 
  gc1=dt*dxi
  gc2=dt*dyi
  gc3=dt*dzi
!  print *,'gc1,gc2,gc3',gc1,gc2,gc3;
 !----------------------------------------------------------------------!;
  ! Do linear extrapolation to time level n+1;
  !----------------------------------------------------------------------!;
  !-- x direction -----------;
FullXYZDomainLoopDC(
        !First extrapolate the contravariant velocities to n+1;
        !Contravariant velocities(momentum) at level n:;
        ox0 =ox(i,j,k,0)*h(i,j,k);
        !Contravariant velocities(momentum) at level n-1:;
        ox2 =ox(i,j,k,2);
        !Construct n+1 contravariant velocities(momentum);
        ox1 =ox0+dt_ratio*(ox0-ox2);
        !Construct and store n+1/2 contravariant velocities(momentum);
        ox_i  =0.5_euwp*(ox1+ox0 )*gc1;
        ox(i,j,k,1)=ox_i;

        oy0   =oy(i,j  ,k,0)*h(i,j  ,k);
        oy2   =oy(i,j  ,k,2);
        oy1   =oy0+dt_ratio*(oy0    -oy2);
        oy_j  =0.5_euwp*(oy1 +oy0 )*gc2;
        oy(i,j,k,1)=oy_j;

        oz0 =oz(i,j,k  ,0)*h(i,j,k  );
        oz2 =oz(i,j,k  ,2);
        oz1 =oz0+dt_ratio*(oz0-oz2);
        oz_k  =0.5_euwp*(oz1 +oz0 )*gc3;
        oz(i,j,k,1)=oz_k;
)
!velprdA_to_C intentionally missing, needs to be called separately
   CALL ttend(40)
    END SUBROUTINE velprd_traj0_adapt_dt 

    SUBROUTINE velprdA_to_C(ox,oy,oz,u1,u2,u3, &
               ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih,do_serialization_in,do_serialization_out )
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY: ttbeg,ttend
        USE mpi_parallel, ONLY: iupx,iupy,iupz 
        USE mpi_parallel, ONLY: updated,updatelrd,updatebtd,updategsd
        LOGICAL,OPTIONAL :: do_serialization_in
        LOGICAL,OPTIONAL :: do_serialization_out
        INTEGER(KIND=iintegers), INTENT(IN) :: np,mp,lp,ih
        REAL_euwp, INTENT(INOUT) ::                            &
                      ox(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                      oy(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                      oz(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

        REAL_euwp, INTENT(OUT) ::  &
                          u1(1-ih:np+1+ih, 1-ih:mp+ih, 1-ih:lp+ih),    &
                          u2(1-ih:np+ih, 1-ih:mp+1+ih, 1-ih:lp+ih),    &
                          u3(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih)


        INTEGER(KIND=iintegers),INTENT(IN) :: ipoles0,ibcx0,ibcy0,ibcz0
        INTEGER(KIND=iintegers) :: i,j,k

  CALL ttbeg(40)
  CALL ttbeg(1040)
        IF(PRESENT(do_serialization_in).OR.PRESENT(do_serialization_out)) THEN
!Do nothing but remove compiler warning about unused dummy variables
        ENDIF

  CALL ttend(1040)

  CALL updatelrd(ox(:,:,:),np,mp,lp,np,mp,lp,iupx,ih)
 !Prepare special values of ox,oy in halo for one-sided extrapolation in x0;
  !For example, x0(1,:,:)=1.5*ox(1,:,:)-0.5*ox(2,:,:);
  CALL ttbeg(1040)
  IF(ibcx0.eq.0) THEN
    IF (leftedge== 1) THEN
X2DWallFullYZDomainLoopDC(ox(0   ,j,k)=2._euwp*ox(1 ,j,k)-ox(2   ,j,k);)
    ENDIF
    IF (rightedge == 1) THEN 
X2DWallFullYZDomainLoopDC(ox(np+1,j,k)=2._euwp*ox(np,j,k)-ox(np-1,j,k);)
    ENDIF
  ENDIF

  CALL ttend(1040,.TRUE.)
  CALL updatebtd(oy(:,:,:),np,mp,lp,np,mp,lp,iupy,ih)
  !Prepare special values of ox,oy in halo for one-sided extrapolation in y0;
  CALL ttbeg(1040)
  IF(ibcy0.eq.0) THEN
    IF (botedge == 1) THEN 
 Y2DWallFullXZDomainLoopDC(oy(i,0    ,k)=2._euwp*oy(i, 1,k)-oy(i,2   ,k);)
    ENDIF
    IF (topedge == 1) THEN 
  Y2DWallFullXZDomainLoopDC(oy(i,mp+1,k)=2._euwp*oy(i,mp,k)-oy(i,mp-1,k);)
    ENDIF
  ENDIF

  CALL ttend(1040,.TRUE.)
  CALL updategsd(oz(:,:,:),np,mp,lp,np,mp,lp,iupz,ih)
  !Prepare special values of ox,oy in halo for one-sided extrapolation in z0
  CALL ttbeg(1040)
  IF(ibcz0.eq.0) THEN
    IF(gndedge == 1) THEN 
 Z2DWallFullXYDomainLoopDC(oz(i,j,0   )=2._euwp*oz(i,j,1 )-oz(i,j,2   );)
    ENDIF
    IF(skyedge == 1) THEN 
 Z2DWallFullXYDomainLoopDC(oz(i,j,lp+1)=2._euwp*oz(i,j,lp)-oz(i,j,lp-1);)
    ENDIF
  ENDIF
  CALL ttend(1040,.TRUE.)


  CALL ttbeg(1040)
  ! Create staggered C-grid advective non-dimensional velocities (momentum)

  ! Update for further use in MPDATA 
  IF(ipoles0.eq.0) THEN
    IF(ibcx0 == 1) THEN
CGRIDInnerXFullYZDomainLoopDC(
!     u1(1+leftedge:np,1:mp,1:lp)=0.5_euwp*(ox(1+leftedge:np,1:mp,1:lp,1)+ox(0+leftedge:np-1,1:mp,1:lp,1));
      u1(i,j,k)=0.5_euwp*(ox(i,j,k)+ox(i-1,j,k));
)
      IF( leftedge.eq.1) THEN
X2DWallFullYZDomainLoopDC(
      u1(   1,j,k)=0.5_euwp*(ox(   0,j,k)+ox(  -1,j,k));
)
      ENDIF
      IF(rightedge.eq.1) THEN 
X2DWallFullYZDomainLoopDC(
      u1(np+1,j,k)=0.5_euwp*(ox(np+1,j,k)+ox(np+2,j,k));
)
      ENDIF
      CALL ttend(1040,.TRUE.)
      CALL updated(u1,np+rightedge,mp,lp,np+1,mp,lp,1,ih)
      CALL ttbeg(1040)
    ELSE
CGRIDXFullYZDomainLoopDC(
      u1(i,j,k)=0.5_euwp*(ox(i,j,k)+ox(i-1,j,k));
)
      ! For non-periodic bc, it doesn`t matter what range is given at the
      ! rightedge since it is not communicated to leftedge anyway
      CALL ttend(1040,.TRUE.)
      CALL updated(u1,np+rightedge,mp,lp,np+1,mp,lp,1,ih)
      CALL ttbeg(1040)
    ENDIF 
  ELSE
CGRIDXFullYZDomainLoopDC(
      u1(i,j,k)=0.5_euwp*(ox(i,j,k)+ox(i-1,j,k));
)
      ! For non-periodic bc, it doesn`t matter what range is given at the
      ! rightedge since it is not communicated to leftedge anyway
      CALL ttend(1040,.TRUE.)
      CALL updated(u1,np,mp,lp,np+1,mp,lp,1,ih)
      CALL ttbeg(1040)
  ENDIF


  ! Update for further use in MPDATA 
  IF(ipoles0.eq.0) THEN
    IF(ibcy0 == 1) THEN
CGRIDInnerYFullXZDomainLoopDC(
       u2(i,j,k)=0.5_euwp*(oy(i,j,k)+oy(i,j-1,k));
)
      IF(botedge.eq.1) THEN
Y2DWallFullXZDomainLoopDC(
              u2(i,   1,k)=0.5_euwp*(oy(i,0   ,k)+oy(i  ,-1,k));
)
      ENDIF
      IF(topedge.eq.1) THEN
Y2DWallFullXZDomainLoopDC(
              u2(i,mp+1,k)=0.5_euwp*(oy(i,mp+1,k)+oy(i,mp+2,k));
)
      ENDIF
!     u2(1:np,1:mp,1:lp)=0.5_euwp*(oy(1:np,1:mp,1:lp,1)+oy(1:np,0:mp-1,1:lp,1))
      CALL ttend(1040,.TRUE.)
      CALL updated(u2,np,mp+topedge,lp,np,mp+1,lp,iupy,ih)
      CALL ttbeg(1040)
!     IF(botedge.eq.1) u2(1:np,   1,1:lp)=u2(1:np,   0,1:lp)
!     IF(topedge.eq.1) u2(1:np,mp+1,1:lp)=u2(1:np,mp+2,1:lp)
    ELSE
CGRIDYFullXZDomainLoopDC(
      u2(i,j,k)=0.5_euwp*(oy(i,j,k)+oy(i,j-1,k));
)
      CALL ttend(1040,.TRUE.)
      CALL updated(u2,np,mp+topedge,lp,np,mp+1,lp,1,ih)
      CALL ttbeg(1040)
    ENDIF
  ELSE
CGRIDInnerYFullXZDomainLoopDC(
    u2(i,j,k)=0.5_euwp*(oy(i,j,k)+oy(i,j-1,k));
)
    CALL ttend(1040,.TRUE.)
    CALL updated(u2,np,mp+topedge,lp,np,mp+1,lp,1,ih)
    CALL ttbeg(1040)
    IF(botedge.eq.1) THEN
Y2DWallFullXZDomainLoopDC(
            u2(i,   1,k)=0._euwp;
)
    ENDIF
    IF(topedge.eq.1) THEN
Y2DWallFullXZDomainLoopDC(
            u2(i,mp+1,k)=0._euwp;
)
    ENDIF
  ENDIF



  ! Update for further use in MPDATA 
  !CALL update(u3,np,mp,lp,np,mp,lp,1,ih)

  IF(ibcz0 == 1) THEN
CGRIDZInnerXYDomainLoopDC(
    u3(i,j,k)=0.5_euwp*(oz(i,j,k)+oz(i,j,k-1));
)
    IF(gndedge.eq.1) THEN
Z2DWallFullXYDomainLoopDC(
            u3(i,j,1   )=0.5_euwp*(oz(i,j,0   )+oz(i,j,-1  ));
)
    ENDIF
    IF(skyedge.eq.1) THEN 
Z2DWallFullXYDomainLoopDC(
            u3(i,j,lp+1)=0.5_euwp*(oz(i,j,lp+1)+oz(i,j,lp+2));
)
    ENDIF
    CALL ttend(1040,.TRUE.)
    CALL updated(u3,np,mp,lp+skyedge,np,mp,lp,iupz,ih)
    CALL ttbeg(1040)
  ELSE
CGRIDZFullXYDomainLoopDC(
    u3(i,j,k)=0.5_euwp*(oz(i,j,k)+oz(i,j,k-1));
)
    CALL ttend(1040,.TRUE.)
    CALL updated(u3,np,mp,lp+skyedge,np,mp,lp+1,1,ih)
    CALL ttbeg(1040)
  ENDIF 
!------- End of computing advecting velocities 
 CALL ttend(40)
    END SUBROUTINE velprdA_to_C 

    SUBROUTINE velprdA_to_C_ver2(uadv,vadv,wadv,u1,u2,u3, &
               ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih )
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY: ttbeg,ttend
        USE mpi_parallel, ONLY: iupx,iupy,iupz 
        USE mpi_parallel, ONLY: update,updatelr,updatebt,updategs
        INTEGER(KIND=iintegers), INTENT(IN) :: np,mp,lp,ih
        INTEGER(KIND=iintegers), INTENT(IN) :: ipoles0,ibcx0,ibcy0,ibcz0
        REAL_euwp, INTENT(INOUT) ::                            &
                      uadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                      vadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                      wadv(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

        REAL_euwp, INTENT(OUT) ::  &
                          u1(1-ih:np+1+ih, 1-ih:mp+ih, 1-ih:lp+ih),    &
                          u2(1-ih:np+ih, 1-ih:mp+1+ih, 1-ih:lp+ih),    &
                          u3(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih)


  CALL ttbeg(40)
  CALL updatelr(uadv,np,mp,lp,np,mp,lp,iupx,ih)
 !Prepare special values of ox,oy in halo for one-sided extrapolation in x0
  !For example, x0(1,:,:)=1.5*ox(1,:,:)-0.5*ox(2,:,:)
  CALL ttbeg(1040)
  IF(ibcx0.eq.0) THEN
    IF (leftedge== 1) &
      uadv(0   ,1:mp,1:lp)=2._euwp*uadv(1 ,1:mp,1:lp)-uadv(2   ,1:mp,1:lp)
    IF (rightedge == 1) & 
      uadv(np+1,1:mp,1:lp)=2._euwp*uadv(np,1:mp,1:lp)-uadv(np-1,1:mp,1:lp)
  ENDIF

  CALL ttend(1040)
  CALL updatebt(vadv,np,mp,lp,np,mp,lp,iupy,ih)
  !Prepare special values of ox,oy in halo for one-sided extrapolation in y0
  CALL ttbeg(1040)
  IF(ibcy0.eq.0) THEN
    IF (botedge == 1) &
      vadv(1:np,0   ,1:lp)=2._euwp*vadv(1:np, 1,1:lp)-vadv(1:np,2   ,1:lp)
    IF (topedge == 1) & 
      vadv(1:np,mp+1,1:lp)=2._euwp*vadv(1:np,mp,1:lp)-vadv(1:np,mp-1,1:lp)
  ENDIF

  CALL ttend(1040,.TRUE.)
  CALL updategs(wadv,np,mp,lp,np,mp,lp,iupz,ih)
  !Prepare special values of ox,oy in halo for one-sided extrapolation in z0
  CALL ttbeg(1040)
  IF(ibcz0.eq.0) THEN
    IF(gndedge == 1) &
      wadv(1:np,1:mp,0   )=2._euwp*wadv(1:np,1:mp,1 )-wadv(1:np,1:mp,2   )
    IF(skyedge == 1) &
      wadv(1:np,1:mp,lp+1)=2._euwp*wadv(1:np,1:mp,lp)-wadv(1:np,1:mp,lp-1)
  ENDIF
  CALL ttend(1040,.TRUE.)

  CALL ttbeg(1040)
  ! Create staggered C-grid advective non-dimensional velocities (momentum)

  ! Update for further use in MPDATA 
  IF(ipoles0.eq.0) THEN
    IF(ibcx0 == 1) THEN
      u1(1+leftedge:np,1:mp,1:lp)=0.5_euwp*( uadv(1+leftedge:np  ,1:mp,1:lp)   &
                                            +uadv(0+leftedge:np-1,1:mp,1:lp))
      IF( leftedge.eq.1)   &
      u1(            1,1:mp,1:lp)=0.5_euwp*( uadv(              0,1:mp,1:lp)   &
                                            +uadv(             -1,1:mp,1:lp))
      IF(rightedge.eq.1)   &
      u1(np+1         ,1:mp,1:lp)=0.5_euwp*( uadv(np+1           ,1:mp,1:lp)   &
                                            +uadv(np+2           ,1:mp,1:lp))
!      u1(1:np,1:mp,1:lp)=0.5_euwp*(ox(1:np,1:mp,1:lp,1)+ox(0:np-1,1:mp,1:lp,1))
!     u1(1:np+rightedge,1:mp,1:lp)=0.5_euwp*(ox(1:np+rightedge,1:mp,1:lp,1)+ox(0:np-1+rightedge,1:mp,1:lp,1))
      CALL ttend(1040,.TRUE.)
      CALL update(u1,np+rightedge,mp,lp,np+1,mp,lp,1,ih)
      CALL ttbeg(1040)
!     IF( leftedge.eq.1) u1(   1,1:mp,1:lp)=u1(   0,1:mp,1:lp)
!     IF(rightedge.eq.1) u1(np+1,1:mp,1:lp)=u1(np+2,1:mp,1:lp)
    ELSE
      u1(1:np+rightedge,1:mp,1:lp)=0.5_euwp*( uadv(1:np  +rightedge,1:mp,1:lp)    &
                                             +uadv(0:np-1+rightedge,1:mp,1:lp))
      ! For non-periodic bc, it doesn`t matter what range is given at the
      ! rightedge since it is not communicated to leftedge anyway
      CALL ttend(1040,.TRUE.)
      CALL update(u1,np+rightedge,mp,lp,np+1,mp,lp,1,ih)
      CALL ttbeg(1040)
    ENDIF 
  ELSE
      u1(1:np+rightedge,1:mp,1:lp)=0.5_euwp*( uadv(1:np  +rightedge,1:mp,1:lp)   &
                                             +uadv(0:np-1+rightedge,1:mp,1:lp))
      ! For non-periodic bc, it doesn`t matter what range is given at the
      ! rightedge since it is not communicated to leftedge anyway
      CALL ttend(1040,.TRUE.)
      CALL update(u1,np,mp,lp,np+1,mp,lp,1,ih)
      CALL ttbeg(1040)
  ENDIF

  ! Update for further use in MPDATA 
  IF(ipoles0.eq.0) THEN
    IF(ibcy0 == 1) THEN
      u2(1:np,1+botedge:mp,1:lp)=0.5_euwp*( vadv(1:np,1+botedge:mp  ,1:lp)     &
                                           +vadv(1:np,0+botedge:mp-1,1:lp))
      IF(botedge.eq.1) &
      u2(1:np,           1,1:lp)=0.5_euwp*( vadv(1:np,0              ,1:lp)    &
                                           +vadv(1:np,-1             ,1:lp))
      IF(topedge.eq.1) &
      u2(1:np,mp+1        ,1:lp)=0.5_euwp*( vadv(1:np,mp+1           ,1:lp)    &
                                           +vadv(1:np,mp+2           ,1:lp))
!     u2(1:np,1:mp,1:lp)=0.5_euwp*(vadv(1:np,1:mp,1:lp,1)+vadv(1:np,0:mp-1,1:lp,1))
      CALL ttend(1040,.TRUE.)
      CALL update(u2,np,mp+topedge,lp,np,mp+1,lp,iupy,ih)
      CALL ttbeg(1040)
!     IF(botedge.eq.1) u2(1:np,   1,1:lp)=u2(1:np,   0,1:lp)
!     IF(topedge.eq.1) u2(1:np,mp+1,1:lp)=u2(1:np,mp+2,1:lp)
    ELSE
      u2(1:np,1:mp+topedge,1:lp)=0.5_euwp*( vadv(1:np,1:mp  +topedge  ,1:lp)   &
                                           +vadv(1:np,0:mp-1+topedge  ,1:lp))
      CALL ttend(1040,.TRUE.)
      CALL update(u2,np,mp+topedge,lp,np,mp+1,lp,1,ih)
      CALL ttbeg(1040)
    ENDIF
  ELSE
      u2(1:np,1+botedge:mp,1:lp)=0.5_euwp*( vadv(1:np,1+botedge:mp    ,1:lp)  &
                                           +vadv(1:np,0+botedge:mp-1  ,1:lp))
    CALL ttend(1040,.TRUE.)
    CALL update(u2,np,mp+topedge,lp,np,mp+1,lp,1,ih)
    CALL ttbeg(1040)
    IF(botedge.eq.1)  u2(1:np,   1,1:lp)=0._euwp
    IF(topedge.eq.1)  u2(1:np,mp+1,1:lp)=0._euwp
  ENDIF



  ! Update for further use in MPDATA 
  !CALL update(u3,np,mp,lp,np,mp,lp,1,ih)

  IF(ibcz0 == 1) THEN
    u3(1:np,1:mp,1+gndedge:lp)=0.5_euwp*( wadv(1:np,1:mp,1+gndedge:lp  )     &
                                         +wadv(1:np,1:mp,0+gndedge:lp-1))
    IF(gndedge.eq.1)  &
    u3(1:np,1:mp,1           )=0.5_euwp*( wadv(1:np,1:mp,0             )     &
                                         +wadv(1:np,1:mp,-1            ))
    IF(skyedge.eq.1)  &
    u3(1:np,1:mp,lp+1        )=0.5_euwp*( wadv(1:np,1:mp,lp+1          )     &
                                         +wadv(1:np,1:mp,lp+2          ))
    CALL ttend(1040,.TRUE.)
    CALL update(u3,np,mp,lp+skyedge,np,mp,lp,iupz,ih)
    CALL ttbeg(1040)
  ELSE
    u3(1:np,1:mp,1:lp+skyedge)=0.5_euwp*( wadv(1:np,1:mp,1:lp  +skyedge)     &
                                         +wadv(1:np,1:mp,0:lp-1+skyedge))
    CALL ttend(1040,.TRUE.)
    CALL update(u3,np,mp,lp+skyedge,np,mp,lp+1,1,ih)
    CALL ttbeg(1040)
  ENDIF 
!------- End of computing advecting velocities 
  CALL ttend(40)
    END SUBROUTINE velprdA_to_C_ver2 
END MODULE module_velprd 


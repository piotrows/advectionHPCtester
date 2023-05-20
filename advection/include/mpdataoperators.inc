    PURE FUNCTION pp(y)
        REAL(KIND=euwp),INTENT(IN) :: y
        REAL(KIND=euwp) :: pp
        pp=max(0._euwp,y)
    END FUNCTION pp
    PURE FUNCTION pn(y)
        REAL(KIND=euwp),INTENT(IN) :: y
        REAL(KIND=euwp) :: pn
        pn=-min(0._euwp,y)
    END FUNCTION pn
    PURE FUNCTION donor(y1,y2,a)
        REAL(KIND=euwp),INTENT(IN) :: y1,y2,a
        REAL(KIND=euwp) :: donor
        donor = max(0._euwp,a)*y1 - (-min(0._euwp,a)*y2)
    END FUNCTION donor
    PURE FUNCTION rat2(z1,z2)
        REAL(KIND=euwp),INTENT(IN) :: z1,z2
        REAL(KIND=euwp) :: rat2
        rat2=(z2-z1)*.5_euwp
    END FUNCTION rat2
    PURE FUNCTION rat4(z0,z1,z2,z3)
        REAL(KIND=euwp),INTENT(IN) :: z0,z1,z2,z3
        REAL(KIND=euwp) :: rat4
        rat4=((z3+z2)-(z1+z0))*.25_euwp
    END FUNCTION rat4
    PURE FUNCTION vdyf(x1,x2,a,rinv)
        REAL(KIND=euwp),INTENT(IN) :: x1,x2,a,rinv
        REAL(KIND=euwp) :: vdyf
        vdyf=(abs(a)-a*a*rinv)*rat2(x1,x2)
    END FUNCTION vdyf
    PURE FUNCTION vcorr(a,b,y0,y1,y2,y3,rinv)
        REAL(KIND=euwp),INTENT(IN) :: y0,y1,y2,y3,a,b,rinv
        REAL(KIND=euwp) :: vcorr
        vcorr=-0.125_euwp*a*rinv*(b*rat4(y0,y1,y2,y3))
    END FUNCTION vcorr
    PURE FUNCTION vdyfa(x1,x2,a,rinv)
          USE epsilons, ONLY: ep => ep_nonos
        REAL(KIND=euwp),INTENT(IN) :: x1,x2,a,rinv
        REAL(KIND=euwp) :: vdyfa
        vdyfa=(abs(a)-a*a*rinv)*(abs(x2)-abs(x1)) &
            /(abs(x2)+abs(x1)+ep)
    END FUNCTION vdyfa
    PURE FUNCTION vcorra(a,b,y1,y2,rinv)
        REAL(KIND=euwp),INTENT(IN) :: y1,y2,a,b,rinv
        REAL(KIND=euwp) :: vcorra
        vcorra=-0.125_euwp*a*b*y1*rinv/(y2)
    END FUNCTION vcorra

#define PURE

    PURE FUNCTION eval_v1(xant,u1,u2,u3,hmx,i,j,k,np,mp,lp,ih)
        INTEGER(KIND=iintegers),INTENT(IN) ::i,j,k,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &
            u1, u2, u3, xant
        REAL(KIND=euwp),INTENT(IN) :: hmx
        REAL(KIND=euwp) :: eval_v1

        eval_v1 =vdyf(xant(i-1,j,k),xant(i,j,k),u1(i,j,k),hmx)       &
                +vcorr(u1(i  ,j,k),                                  &
                u2(i-1,j,k)+u2(i-1,j+1,k)+                           &
                u2(i,j+1,k)+u2(i  ,j  ,k),                           &
                xant(i-1,j-1,k),xant(i,j-1,k),                       &
                xant(i-1,j+1,k),xant(i,j+1,k),                       &
                hmx)                                                 &
                +vcorr(u1(i,j,k),                                    &
                u3(i-1,j,k  )+u3(i-1,j,k+1)+                         &
                u3(i  ,j,k+1)+u3(i  ,j,k  ),                         &
                xant(i-1,j,k-1),xant(i,j,k-1),                       &
                xant(i-1,j,k+1),xant(i,j,k+1),                       &
                hmx)
    END FUNCTION eval_v1
    PURE FUNCTION eval_v2(xant,u1,u2,u3,hmy,i,j,k,np,mp,lp,ih)
        INTEGER(KIND=iintegers),INTENT(IN) ::i,j,k,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &
            u1, u2, u3, xant
        REAL(KIND=euwp),INTENT(IN) :: hmy
        REAL(KIND=euwp) :: eval_v2
        eval_v2=vdyf(xant(i,j-1,k),xant(i,j,k),u2(i,j,k),hmy)       &
               +vcorr(u2(i,j,k),                                    &
               u1(i  ,j-1,k)+u1(i  ,j  ,k)+                         &
               u1(i+1,j  ,k)+u1(i+1,j-1,k),                         &
               xant(i-1,j-1,k),xant(i-1,j,k),                       &
               xant(i+1,j-1,k),xant(i+1,j,k),                       &
               hmy)                                                 &
               +vcorr(u2(i,j,k),                                    &
               u3(i,j-1,k  )+u3(i,j  ,k  )+                         &
               u3(i,j  ,k+1)+u3(i,j-1,k+1),                         &
               xant(i,j-1,k-1), xant(i,j  ,k-1),                    &
               xant(i,j-1,k+1), xant(i,j  ,k+1),                    &
                hmy)

    END FUNCTION eval_v2

    PURE FUNCTION eval_v3(xant,u1,u2,u3,hmz,i,j,k,np,mp,lp,ih)
        INTEGER(KIND=iintegers),INTENT(IN) ::i,j,k,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &
            u1, u2, u3, xant
        REAL(KIND=euwp),INTENT(IN) :: hmz
        REAL(KIND=euwp) :: eval_v3
        eval_v3=vdyf(xant(i,j,k-1),xant(i,j,k),u3(i,j,k),           &
               hmz)                                                 &
               +vcorr(u3(i,j,k),                                    &
               u1(i  ,j,k-1)+u1(i  ,j,k  )+                         &
               u1(i+1,j,k  )+u1(i+1,j,k-1),                         &
               xant(i-1,j,k-1), xant(i-1,j,k  ),                    &
               xant(i+1,j,k-1), xant(i+1,j,k  ),                    &
               hmz)                                                 &
               +vcorr(u3(i,j,k),                                    &
               u2(i,j  ,k-1)+u2(i,j+1,k-1)+                         &
               u2(i,j+1,k  )+u2(i,j  ,k  ),                         &
               xant(i,j-1,k-1), xant(i,j-1,k  ),                    &
               xant(i,j+1,k-1), xant(i,j+1,k  ),                    &
               hmz)
    END FUNCTION eval_v3
    PURE FUNCTION eval_v1_a(xant,u1,u2,u3,hmx,i,j,k,np,mp,lp,ih)
        INTEGER(KIND=iintegers),INTENT(IN) ::i,j,k,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &
            u1, u2, u3, xant
        REAL(KIND=euwp),INTENT(IN) :: hmx
        REAL(KIND=euwp) :: eval_v1_a
        eval_v1_a =vdyf(xant(i-1,j,k),xant(i,j,k),u1(i,j,k),hmx)        &
                   -0.125*u1(i  ,j,k)*hmx*                              &
                   ((u2(i-1,j,k)+u2(i-1,j+1,k)+                         &
                     u2(i,j+1,k)+u2(i  ,j  ,k))*                        &
                     rat4(xant(i-1,j-1,k),xant(i,j-1,k),                &
                          xant(i-1,j+1,k),xant(i,j+1,k))                &
                   +                                                    &
                    (u3(i-1,j,k  )+u3(i-1,j,k+1)+                       &
                     u3(i  ,j,k+1)+u3(i  ,j,k  ))*                      &
                     rat4(xant(i-1,j,k-1),xant(i,j,k-1),                &
                          xant(i-1,j,k+1),xant(i,j,k+1)))
    END FUNCTION eval_v1_a
    PURE FUNCTION eval_v2_a(xant,u1,u2,u3,hmy,i,j,k,np,mp,lp,ih)
        INTEGER(KIND=iintegers),INTENT(IN) ::i,j,k,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &
            u1, u2, u3, xant
        REAL(KIND=euwp),INTENT(IN) :: hmy
        REAL(KIND=euwp) :: eval_v2_a
        eval_v2_a=vdyf(xant(i,j-1,k),xant(i,j,k),u2(i,j,k),hmy)        &
                -0.125_euwp*u2(i,j,k)*hmy*                             &
                ((u1(i  ,j-1,k)+u1(i  ,j  ,k)+                         &
                  u1(i+1,j  ,k)+u1(i+1,j-1,k))*                        &
                  rat4(xant(i-1,j-1,k),xant(i-1,j,k),                  &
                       xant(i+1,j-1,k),xant(i+1,j,k))                  &
                  +                                                    &
                   (u3(i,j-1,k  )+u3(i,j  ,k  )+                       &
                    u3(i,j  ,k+1)+u3(i,j-1,k+1))*                      &
                    rat4(xant(i,j-1,k-1), xant(i,j  ,k-1),             &
                         xant(i,j-1,k+1), xant(i,j  ,k+1)))
    END FUNCTION eval_v2_a
    PURE FUNCTION eval_v3_a(xant,u1,u2,u3,hmz,i,j,k,np,mp,lp,ih)
        INTEGER(KIND=iintegers),INTENT(IN) ::i,j,k,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &
            u1, u2, u3, xant
        REAL(KIND=euwp),INTENT(IN) :: hmz
        REAL(KIND=euwp) :: eval_v3_a
        eval_v3_a=vdyf(xant(i,j,k-1),xant(i,j,k),u3(i,j,k),hmz)      &
                -0.125_euwp*u3(i,j,k)*hmz*                           &
                ((u1(i  ,j,k-1)+u1(i  ,j,k  )+                       &
                  u1(i+1,j,k  )+u1(i+1,j,k-1))*                      &
                 rat4(xant(i-1,j,k-1), xant(i-1,j,k  ),              &
                      xant(i+1,j,k-1), xant(i+1,j,k  ))              &
                 +                                                   &
                 (u2(i,j  ,k-1)+u2(i,j+1,k-1)+                       &
                  u2(i,j+1,k  )+u2(i,j  ,k  ))*                      &
                 rat4(xant(i,j-1,k-1), xant(i,j-1,k  ),              &
                      xant(i,j+1,k-1), xant(i,j+1,k  )))
    END FUNCTION eval_v3_a
#undef PURE 
    !---------------------------------------------------------------------!
    !
    !---------------------------------------------------------------------!
    SUBROUTINE putbcinx(x,bcx,bcy,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge,gndedge,skyedge 
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers) :: ibcx0,ibcy0,ibcz0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x
        REAL(KIND=euwp) :: &
            bcx(mp,lp, 2), &               ! local array
            bcy(np,lp, 2)                  ! local array

        IF(ibcx0.eq.0) THEN
          IF ( leftedge == 1) x(0   ,1:mp,1:lp)=bcx(1:mp,1:lp,1)
          IF (rightedge == 1) x(np+1,1:mp,1:lp)=bcx(1:mp,1:lp,2)
        ENDIF
        IF(ibcy0.eq.0) THEN
          IF (  botedge == 1) x(1:np,0   ,1:lp)=bcy(1:np,1:lp,1)
          IF (  topedge == 1) x(1:np,mp+1,1:lp)=bcy(1:np,1:lp,2)
        ENDIF
        IF(ibcz0.eq.0) THEN
          IF (  gndedge == 1) x(1:np,1:mp,0   )=  x(1:np,1:mp,2)
          IF (  skyedge == 1) x(1:np,1:mp,lp+1)=  x(1:np,1:mp,lp-1)
        ENDIF
    END SUBROUTINE putbcinx


    !---------------------------------------------------------------------!
    !
    !---------------------------------------------------------------------!
    SUBROUTINE putbcinxh(x,bcx,bcy,ibcx0,ibcy0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers) :: ibcx0,ibcy0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(OUT) :: &
            x
        REAL(KIND=euwp), INTENT(IN) :: &
            bcx(mp,lp, 2), &               ! local array
            bcy(np,lp, 2)                  ! local array
      IF(ibcx0.eq.0) THEN
        IF ( leftedge == 1) x(0   ,1:mp,1:lp)=bcx(1:mp,1:lp,1)
        IF (rightedge == 1) x(np+1,1:mp,1:lp)=bcx(1:mp,1:lp,2)
      ENDIF
      IF(ibcy0.eq.0) THEN
        IF (  botedge == 1) x(1:np,0   ,1:lp)=bcy(1:np,1:lp,1)
        IF (  topedge == 1) x(1:np,mp+1,1:lp)=bcy(1:np,1:lp,2)
      ENDIF
    END SUBROUTINE putbcinxh
    SUBROUTINE zerogradbcinxh(x,ibcx0,ibcy0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers) :: ibcx0,ibcy0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(INOUT) :: &
            x
        IF(ibcx0 == 0) THEN
          IF ( leftedge == 1) x(0   ,1:mp,1:lp)=x(1   ,1:mp,1:lp)
          IF (rightedge == 1) x(np+1,1:mp,1:lp)=x(  np,1:mp,1:lp)
        ENDIF
        IF(ibcy0 == 0) THEN
          IF (  botedge == 1) x(1:np,0   ,1:lp)=x(1:np,1   ,1:lp)
          IF (  topedge == 1) x(1:np,mp+1,1:lp)=x(1:np,  mp,1:lp)
        ENDIF
    END SUBROUTINE zerogradbcinxh
    SUBROUTINE zerogradbcinx(x,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge,gndedge,skyedge 
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers) :: ibcx0,ibcy0,ibcz0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(INOUT) :: &
            x
        IF(ibcx0 == 0) THEN
          IF ( leftedge == 1) x(0   ,1:mp,1:lp)=x(1   ,1:mp,1:lp)
          IF (rightedge == 1) x(np+1,1:mp,1:lp)=x(  np,1:mp,1:lp)
        ENDIF
        IF(ibcy0 == 0) THEN
          IF (  botedge == 1) x(1:np,0   ,1:lp)=x(1:np,1   ,1:lp)
          IF (  topedge == 1) x(1:np,mp+1,1:lp)=x(1:np,  mp,1:lp)
        ENDIF
        IF(ibcz0 == 0) THEN
          IF (  gndedge == 1) x(1:np,1:mp,0   )=  x(1:np,1:mp,2)
          IF (  skyedge == 1) x(1:np,1:mp,lp+1)=  x(1:np,1:mp,lp-1)
        ENDIF
    END SUBROUTINE zerogradbcinx

    !---------------------------------------------------------------------!
    !

    SUBROUTINE cp_scnd_last_to_halo_z(x,ibcz0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  gndedge,skyedge 
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers) :: ibcz0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x
      IF(ibcz0.eq.0) THEN
        IF (gndedge == 1) x(1:np,1:mp,0)=x(1:np,1:mp,2)
        IF (skyedge == 1) x(1:np,1:mp,lp+1)=x(1:np,1:mp,lp-1)
      ENDIF
    END SUBROUTINE cp_scnd_last_to_halo_z
    !---------------------------------------------------------------------!
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
    !---------------------------------------------------------------------!
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
    !---------------------------------------------------------------------!
    SUBROUTINE cp_scnd_last_to_halo_z_full(x,ibcz0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:gndedge,skyedge 
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers),INTENT(IN) :: ibcz0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x
          IF(ibcz0.eq.0) THEN
        IF (gndedge == 1)   x(0:np+1,0:mp+1,0     )=x(0:np+1,0:mp+1,2     )
        IF (skyedge == 1)   x(0:np+1,0:mp+1,  lp+1)=x(0:np+1,0:mp+1,  lp-1)
          ENDIF
    END SUBROUTINE cp_scnd_last_to_halo_z_full
    !---------------------------------------------------------------------!
    SUBROUTINE cp_last_to_halo_xyz_full(x,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge,gndedge,skyedge 
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers),INTENT(IN) :: ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x
        IF(ibcz0.eq.0) THEN
          IF (gndedge == 1)   x(0:np+1,0:mp+1,0     )=x(0:np+1,0:mp+1,1     )
          IF (skyedge == 1)   x(0:np+1,0:mp+1,  lp+1)=x(0:np+1,0:mp+1,  lp  )
        ENDIF
        IF(ipoles0.eq.0) THEN
          IF(ibcy0.eq.0) THEN
            IF (botedge == 1)   x(0:np+1,0     ,0:lp+1)=x(0:np+1,1     ,0:lp+1)
            IF (topedge == 1)   x(0:np+1,mp+1  ,0:lp+1)=x(0:np+1,  mp  ,0:lp+1)
          ENDIF
          IF(ibcx0.eq.0) THEN
            IF (leftedge == 1)  x(0     ,0:mp+1,0:lp+1)=x(1     ,0:mp+1,0:lp+1)
            IF (rightedge == 1) x(np+1  ,0:mp+1,0:lp+1)=x(np    ,0:mp+1,0:lp+1)
          ENDIF
        ENDIF
    END SUBROUTINE cp_last_to_halo_xyz_full
    !---------------------------------------------------------------------!
    !
    !---------------------------------------------------------------------!
    SUBROUTINE cp_scnd_last_to_halo_xy(x,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers) :: ipoles0,ibcx0,ibcy0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x

      IF( ipoles0  == 0) THEN
        IF(ibcy0.eq.0) THEN
          IF (botedge == 1)   x(1:np,0   ,1:lp)=x(1:np,2   ,1:lp)
          IF (topedge == 1)   x(1:np,mp+1,1:lp)=x(1:np,mp-1,1:lp)
        ENDIF
        IF(ibcx0.eq.0) THEN
          IF (leftedge == 1)  x(0   ,1:mp,1:lp)=x(2   ,1:mp,1:lp)
          IF (rightedge == 1) x(np+1,1:mp,1:lp)=x(np-1,1:mp,1:lp)
        ENDIF
      ENDIF
    END SUBROUTINE cp_scnd_last_to_halo_xy

    SUBROUTINE cp_scnd_last_to_halo_xy_full(x,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge
        !---------------------------------------------------------------------!
        INTEGER(KIND=iintegers) :: ipoles0,ibcx0,ibcy0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x

      IF( ipoles0  == 0) THEN
        IF(ibcy0.eq.0) THEN
          IF (botedge == 1)   x(0:np+1,0   ,1:lp)=x(0:np+1,2   ,1:lp)
          IF (topedge == 1)   x(0:np+1,mp+1,1:lp)=x(0:np+1,mp-1,1:lp)
        ENDIF
        IF(ibcx0.eq.0) THEN
          IF (leftedge == 1)  x(0   ,0:mp+1,1:lp)=x(2   ,0:mp+1,1:lp)
          IF (rightedge == 1) x(np+1,0:mp+1,1:lp)=x(np-1,0:mp+1,1:lp)
        ENDIF
      ENDIF
    END SUBROUTINE cp_scnd_last_to_halo_xy_full


    SUBROUTINE cp_last_to_halo_xyz(x,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge,gndedge,skyedge 
        INTEGER(KIND=iintegers), INTENT(IN) :: ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x
      IF(ibcz0.eq.0) THEN
        IF (gndedge == 1) x(1:np,1:mp,0   )=x(1:np,1:mp,1 )
        IF (skyedge == 1) x(1:np,1:mp,lp+1)=x(1:np,1:mp,lp)
      ENDIF
      IF( ipoles0  == 0) THEN
      IF(ibcy0.eq.0) THEN
        IF (botedge == 1)   x(1:np,0   ,1:lp)=x(1:np,1 ,1:lp)
        IF (topedge == 1)   x(1:np,mp+1,1:lp)=x(1:np,mp,1:lp)
      ENDIF
      IF(ibcx0.eq.0) THEN
        IF (leftedge == 1)  x(0   ,1:mp,1:lp)=x(1 ,1:mp,1:lp)
        IF (rightedge == 1) x(np+1,1:mp,1:lp)=x(np,1:mp,1:lp)
      ENDIF
      ENDIF
    END SUBROUTINE cp_last_to_halo_xyz
    SUBROUTINE cp_last_to_halo_z(x,ibcz0,np,mp,lp,ih)
        USE mpi_parallel, ONLY:  gndedge,skyedge 
        INTEGER(KIND=iintegers), INTENT(IN) :: ibcz0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x
      IF(ibcz0.eq.0) THEN
        IF (gndedge == 1) x(1:np,1:mp,0   )=x(1:np,1:mp,1 )
        IF (skyedge == 1) x(1:np,1:mp,lp+1)=x(1:np,1:mp,lp)
      ENDIF
    END SUBROUTINE cp_last_to_halo_z


    !---------------------------------------------------------------------!
    !
    !---------------------------------------------------------------------!
    SUBROUTINE cp_last_to_halo_xy(x,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge 
        INTEGER(KIND=iintegers) :: ipoles0,ibcx0,ibcy0,np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            x
        IF(ipoles0.eq.0) THEN
          IF(ibcy0.eq.0) THEN
            IF (botedge == 1)   x(1:np,0   ,1:lp)=x(1:np,1   ,1:lp)
            IF (topedge == 1)   x(1:np,mp+1,1:lp)=x(1:np,mp  ,1:lp)
          ENDIF   
          IF(ibcx0.eq.0) THEN
            IF (leftedge == 1)  x(0   ,1:mp,1:lp)=x(1   ,1:mp,1:lp)
            IF (rightedge == 1) x(np+1,1:mp,1:lp)=x(np  ,1:mp,1:lp)
          ENDIF   
        ENDIF
    END SUBROUTINE cp_last_to_halo_xy


    !---------------------------------------------------------------------!
    ! Calls in subroutine FLUX_X_BC:
    ! => updatelr (on line <3733>)
    SUBROUTINE flux_x_bc(v1,np,mp,lp,ih,ibctype)
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY:  updatelr
        USE mpi_parallel, ONLY:  leftedge,rightedge 


        !This SUBROUTINE imposes x boundary condition on C-grid variables in MPDATA,
        !i.e.
        !fluxes and antidiffusive velocities
        !Available bc types
        !-1 - DO nothing for non-periodic bcs
        !0 - value zero at bc
        !1 - "zero flux divergence" (e.g. f(1)=f(2))
        !2 - "zero flux           " (e.g. f(1)=-f(2))



        INTEGER(KIND=iintegers) :: np,mp,lp,ih,ibctype
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih+1, 1-ih:mp+ih, 1-ih:lp+ih) :: v1

        SELECT CASE (ibctype)
            CASE(-1)
              !DO nothing
            CASE(0)
                IF (leftedge == 1) THEN
                    v1(1   ,0:mp+1,0:lp+1)= 0._euwp
                END IF
                IF (rightedge == 1) THEN
                    v1(np+1,0:mp+1,0:lp+1)= 0._euwp
                END IF
            CASE(1)
                IF (leftedge == 1) THEN
                    v1(1   ,1:mp,1:lp)=v1(2 ,1:mp,1:lp)
                END IF
                IF (rightedge == 1) THEN
                    v1(np+1,1:mp,1:lp)=v1(np,1:mp,1:lp)
                ENDIF
            CASE(2)
                IF (leftedge == 1) THEN
                    v1(1   ,1:mp,1:lp)=-v1(2 ,1:mp,1:lp)
                END IF
                IF (rightedge == 1) THEN
                    v1(np+1,1:mp,1:lp)=-v1(np,1:mp,1:lp)
                ENDIF
        END SELECT
        CALL updatelr(v1,np+rightedge,mp,lp,np+1,mp,lp,1,ih)
    END SUBROUTINE flux_x_bc

    !---------------------------------------------------------------------!
    !
    !---------------------------------------------------------------------!
    ! Calls in subroutine FLUX_Y_BC:
    ! => updatebt (on line <3780>)
    SUBROUTINE flux_y_bc(v2,np,mp,lp,ih,ibctype)
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY:  updatebt
        USE mpi_parallel, ONLY:  botedge,topedge 
        !This SUBROUTINE imposes y boundary condition on C-grid variables in MPDATA, i.e.
        !fluxes and antidiffusive velocities
        !Available bc types
        !-1 - DO nothing for non-periodic bcs
        !0 - value zero at bc
        !1 - "zero flux divergence" (e.g. f(1)=f(2))
        !2 - "zero flux           " (e.g. f(1)=-f(2))
        INTEGER(KIND=iintegers) :: np,mp,lp,ih,ibctype
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih+1, 1-ih:lp+ih) :: v2


        SELECT CASE (ibctype)
            CASE(-1)
              !DO nothing
            CASE(0)
                IF (botedge == 1) THEN
                    v2(0:np+1,1   ,0:lp+1)= 0._euwp
                END IF
                IF (topedge == 1) THEN
                    v2(0:np+1,mp+1,0:lp+1)= 0._euwp
                END IF
            CASE(1)
                IF (botedge == 1) THEN
                    v2(1:np,1   ,1:lp)=v2(1:np,2 ,1:lp)
                END IF
                IF (topedge == 1) THEN
                    v2(1:np,mp+1,1:lp)=v2(1:np,mp,1:lp)
                END IF
            CASE(2)
                IF (botedge == 1) THEN
                    v2(1:np,1   ,1:lp)=-v2(1:np,2 ,1:lp)
                END IF
                IF (topedge == 1) THEN
                    v2(1:np,mp+1,1:lp)=-v2(1:np,mp,1:lp)
                END IF
        END SELECT
        CALL updatebt(v2,np,mp+topedge,lp,np,mp+1,lp,1,ih)
    END SUBROUTINE flux_y_bc



    !---------------------------------------------------------------------!
    !
    !---------------------------------------------------------------------!
    ! Calls in subroutine FLUX_Z_BC:
    ! => updategs (on line <3829>)
    SUBROUTINE flux_z_bc(v3,np,mp,lp,ih,ibctype)
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY:  updategs
        USE mpi_parallel, ONLY:  skyedge,gndedge 
        INTEGER(KIND=iintegers) :: np,mp,lp,ih,ibctype
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih+1) :: v3


        !This SUBROUTINE imposes z boundary condition on C-grid variables in MPDATA, i.e.
        !fluxes and antidiffusive velocities
        !Available bc types
        !-1 - DO nothing for non-periodic bcs
        !0 - value zero at bc
        !1 - "zero flux divergence" (e.g. f(1)=f(2))
        !2 - "zero flux           " (e.g. f(1)=-f(2))

        IF (ibctype == -1) THEN
           !DO nothing
        ELSE IF (ibctype == 0) THEN
            IF (gndedge == 1) THEN
                v3(0:np+1,0:mp+1, 1)= 0._euwp
            ENDIF
            IF (skyedge == 1) THEN
                v3(0:np+1,0:mp+1,lp+1)= 0._euwp
            ENDIF
        ELSE IF (ibctype == 1) THEN
            IF (gndedge == 1) THEN
                v3(0:np+1,0:mp+1, 1)=v3(0:np+1,0:mp+1, 2)
            ENDIF
            IF (skyedge == 1) THEN
                v3(0:np+1,0:mp+1,lp+1)=v3(0:np+1,0:mp+1,lp)
            ENDIF
        ELSE IF (ibctype == 2) THEN
            IF (gndedge == 1) THEN
                v3(0:np+1,0:mp+1, 1)=-v3(0:np+1,0:mp+1, 2)
            ENDIF
            IF (skyedge == 1) THEN
                v3(0:np+1,0:mp+1,lp+1)=-v3(0:np+1,0:mp+1,lp)
            ENDIF
        END IF !ibctype
        CALL updategs(v3,np,mp,lp+skyedge,np,mp,lp+1,1,ih)
    END SUBROUTINE flux_z_bc

    !---------------------------------------------------------------------!
    !
    !---------------------------------------------------------------------!
    ! Calls in subroutine FLUX_Z_BC_NOHALO:
    ! => updategs (on line <3876>)
    SUBROUTINE flux_z_bc_nohalo(v3,np,mp,lp,ih,ibctype)
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY:  updategs
        USE mpi_parallel, ONLY:  gndedge,skyedge 
        INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih,ibctype
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih+1),INTENT(INOUT) :: v3


        !This SUBROUTINE imposes z boundary condition on C-grid variables in MPDATA, i.e.
        !fluxes and antidiffusive velocities
        !Available bc types
        !-1 - DO nothing for non-periodic bcs
        !0 - value zero at bc
        !1 - "zero flux divergence" (e.g. f(1)=f(2))
        !2 - "zero flux           " (e.g. f(1)=-f(2))

        IF (ibctype == -1) THEN
           !DO nothing
        ELSE IF (ibctype == 0) THEN
            IF (gndedge == 1) THEN
                v3(1:np,1:mp, 1)= 0._euwp
            ENDIF
            IF (skyedge == 1) THEN
                v3(1:np,1:mp,lp+1)= 0._euwp
            ENDIF
        ELSE IF (ibctype == 1) THEN
            IF (gndedge == 1) THEN
                v3(1:np,1:mp, 1)=v3(1:np,1:mp, 2)
            ENDIF
            IF (skyedge == 1) THEN
                v3(1:np,1:mp,lp+1)=v3(1:np,1:mp,lp)
            ENDIF
        ELSE IF (ibctype == 2) THEN
            IF (gndedge == 1) THEN
                v3(1:np,1:mp, 1)=-v3(1:np,1:mp, 2)
            ENDIF
            IF (skyedge == 1) THEN
                v3(1:np,1:mp,lp+1)=-v3(1:np,1:mp,lp)
            ENDIF
        END IF !ibctype
        CALL updategs(v3,np,mp,lp+skyedge,np,mp,lp+1,1,ih)
    END SUBROUTINE flux_z_bc_nohalo

    SUBROUTINE flip_halo_poles(x,np,mp,lp,ih)
      USE mpi_parallel, ONLY:  botedge,topedge 
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: x 
      IF (botedge == 1)  x(0:np+1,0   ,0:lp+1)=-x(0:np+1,0   ,0:lp+1)
      IF (topedge == 1)  x(0:np+1,mp+1,0:lp+1)=-x(0:np+1,mp+1,0:lp+1)
    END SUBROUTINE flip_halo_poles 

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

    SUBROUTINE remove_cyclic_offset_xy_full(x,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)
      USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge
      INTEGER(KIND=iintegers),INTENT(IN) :: ibcx0,ibcy0,ipoles0,np,mp,lp,ih
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
    END SUBROUTINE remove_cyclic_offset_xy_full

    SUBROUTINE remove_cyclic_offset_z_full(x,ibcz0,np,mp,lp,ih)
      USE mpi_parallel, ONLY:  gndedge,skyedge 
      INTEGER(KIND=iintegers),INTENT(IN) :: ibcz0,np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: x 
      IF(ibcz0.eq.1) THEN
        IF ( gndedge == 1) x(0:np+1,0:mp+1,0     )=x(0:np+1,0:mp+1,-1  )
        IF ( skyedge == 1) x(0:np+1,0:mp+1,  lp+1)=x(0:np+1,0:mp+1,lp+2)
      ENDIF
    END SUBROUTINE remove_cyclic_offset_z_full

    SUBROUTINE zero_halo_xy(x,ibcx0,ibcy0,np,mp,lp,ih)
      USE mpi_parallel, ONLY:  leftedge,rightedge,botedge,topedge
      INTEGER(KIND=iintegers),INTENT(IN) :: ibcx0,ibcy0,np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: x
        IF(ibcx0.eq.0) THEN
          IF ( leftedge == 1) x(0   ,0:mp+1,0:lp+1)=0.
          IF (rightedge == 1) x(np+1,0:mp+1,0:lp+1)=0.
        ENDIF
        IF(ibcy0.eq.0) THEN
          IF ( botedge == 1) x(0:np+1,0     ,0:lp+1)=0.
          IF ( topedge == 1) x(0:np+1,mp+1  ,0:lp+1)=0.
        ENDIF
    END SUBROUTINE zero_halo_xy 


    SUBROUTINE enforce_cyclic(x,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
      USE mpi_parallel, ONLY:  rightedge,topedge,skyedge 
       USE mpi_parallel, ONLY:  updatelr,updatebt,updategs
      INTEGER(KIND=iintegers),INTENT(IN) :: ibcx0,ibcy0,ibcz0,ipoles0,np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: x 

       RETURN 

      IF(ipoles0 == 0) THEN
        IF(ibcx0.eq.1) THEN
          CALL updatelr(x,np,mp,lp,np,mp,lp,1,ih)
          IF (rightedge == 1) x(np,1:mp,1:lp)=x(np+1,1:mp,1:lp)
        ENDIF
        IF(ibcy0.eq.1) THEN
          CALL updatebt(x,np,mp,lp,np,mp,lp,1,ih)
          IF ( topedge == 1) x(1:np,mp,1:lp)=x(1:np,mp+1,1:lp)
        ENDIF
      ENDIF
      IF(ibcz0.eq.1) THEN
          CALL updategs(x,np,mp,lp,np,mp,lp,1,ih)
        IF ( skyedge == 1) x(1:np,1:mp,lp)=x(1:np,1:mp,lp+1)
      ENDIF
    END SUBROUTINE enforce_cyclic

    SUBROUTINE enforce_cyclic_z(x,ibcz0,np,mp,lp,ih)
      USE mpi_parallel, ONLY:  updategs,skyedge 
      INTEGER(KIND=iintegers),INTENT(IN) :: ibcz0,np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: x
      IF(ibcz0.eq.1) THEN
          CALL updategs(x,np,mp,lp,np,mp,lp,1,ih)
        IF ( skyedge == 1) x(1:np,1:mp,lp)=x(1:np,1:mp,lp+1)
      ENDIF
    END SUBROUTINE enforce_cyclic_z

    SUBROUTINE enforce_cyclic_xy(x,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)
      USE mpi_parallel, ONLY:  updatelr,updatebt,rightedge,topedge
      INTEGER(KIND=iintegers),INTENT(IN) :: ibcx0,ibcy0,ipoles0,np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: x

      IF(ipoles0 == 0) THEN
        IF(ibcx0.eq.1) THEN
          CALL updatelr(x,np,mp,lp,np,mp,lp,1,ih)
          IF (rightedge == 1) x(np,1:mp,1:lp)=x(np+1,1:mp,1:lp)
        ENDIF
        IF(ibcy0.eq.1) THEN
          CALL updatebt(x,np,mp,lp,np,mp,lp,1,ih)
          IF ( topedge == 1) x(1:np,mp,1:lp)=x(1:np,mp+1,1:lp)
        ENDIF
      ENDIF
    END SUBROUTINE enforce_cyclic_xy

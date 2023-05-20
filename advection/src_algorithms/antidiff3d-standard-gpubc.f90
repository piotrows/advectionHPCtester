MODULE module_antidiff3d_standard_gpubc
   USE precisions, ONLY  : iintegers,euwp
   USE mpi_parallel, ONLY: leftedge,rightedge,botedge,topedge,gndedge,skyedge
#ifdef CUDACODE 
   USE mpi_parallel, ONLY: istream1,istream2
#endif
   USE epsilons, ONLY: ep => ep_nonos

   IMPLICIT NONE
CONTAINS
#include "defines.inc"
#include "mpdataoperators.inc"

        SUBROUTINE antidiff3d_standard_gpubc(lupdatemulti,u1,u2,u3,xant,x_in,rhr,h, &
                                             iflip,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih, &
                                             do_serialization_out )
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY:  update,updatelr,updatebt,updategs,update3
        USE mpi_parallel, ONLY:  updated,updatelrd,updatebtd,updategsd,update3d,updated
        USE mpi_parallel, ONLY:  updatelrd
        USE mpi_parallel, ONLY:  update3_multi
        USE mpi_parallel, ONLY: ttbeg,ttend
        USE mpi_parallel, ONLY: iup,iupx,iupy,iupz 
        USE mpi_parallel, ONLY: enforce_cyclicg =>enforce_cyclic
#ifdef PNETCDF
        USE mpi_parallel, ONLY: pnet_out_chunk
#endif /*PNETCDF*/
        USE scratch_datafields, ONLY: v1,v2,v3,f1,f2,f3,f1ns,f2ns,f3ns,cp,cn
        LOGICAL, INTENT(IN) :: lupdatemulti
        LOGICAL,OPTIONAL :: do_serialization_out
        INTEGER(KIND=iintegers), INTENT(IN) :: iflip,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih
        REAL_euwp, INTENT(IN) ::  & 
                          u1(1-ih:np+1+ih, 1-ih:mp+ih, 1-ih:lp+ih),    &
                          u2(1-ih:np+ih, 1-ih:mp+1+ih, 1-ih:lp+ih),    &
                          u3(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih),    &
                           h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),    &
                         rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

        REAL_euwp, INTENT(INOUT) ::  &
                        x_in(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
        DEV_REAL_euwp, INTENT(INOUT) ::  &
                        xant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

        INTEGER(KIND=iintegers) :: &
            i,j,k


        REAL(KIND=euwp) :: hmx,hmy,hmz,hinv

        REAL(KIND=euwp) :: a1,a2

        REAL(KIND=euwp) ::  mxijk,mnijk,mxijk_o,mnijk_o,v1ijk,v2ijk,v3ijk

        CALL ttbeg(53)
        IF(PRESENT(do_serialization_out)) THEN
!Do nothing but remove compiler warning about unused dummy variables
        ENDIF
            !Full update needed as x is referenced in all edges
             CALL updated(xant,np,mp,lp,np,mp,lp,iup,ih)
            CALL ttbeg(1053)
!           CALL remove_cyclic_offset_full(xant,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
            !At the domain edges, copy the second last  meaningful domain point to x halo
            !This is needed to zero the corrective fluxes vcorr
            CALL cp_scnd_last_to_halo_xy_fulld(xant,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)
!           if(iflip.eq.1.and.ipoles0.eq.1)                    &
!           CALL flip_halo_poles(xant,np,mp,lp,ih)


            !For the first antidiffusive iteration, f1,f2,f3 in traditional MPDATA
            !can be effectively replaced with C-grid advective velocities u1,u2,u3
            !This would not be valid for further antidiffusive iterations.
            IF(gndedge == 1) THEN
              k=1
FullXYDomainLoopDC(
                        hmx= 1._euwp/(0.5_euwp*(h(i-1,j,k)+h(i,j,k)));
                        v1(i,j,k)=vdyfa(xant(i-1,j,k),xant(i,j,k),u1(i,j,k),hmx)      &
                                 +vcorra(u1(i  ,j,k),                                 &
                                  u2(i-1,j  ,k)+u2(i-1,j+1,k)+                        &
                                  u2(i  ,j+1,k)+u2(i  ,j  ,k),                        &
                                  (abs(xant(i-1,j+1,k))+abs(xant(i,j+1,k)))           &
                                 -(abs(xant(i-1,j-1,k))+abs(xant(i,j-1,k))),          &
                                   abs(xant(i-1,j+1,k))+abs(xant(i,j+1,k))            &
                                  +abs(xant(i-1,j-1,k))+abs(xant(i,j-1,k))+ep,        &
                                  hmx);   
                        f1ns(i,j,k)=donor(xant(i-1,j,k),xant(i,j,k),v1(i,j,k));
                        hmy= 1._euwp/(0.5_euwp*(h(i,j-1,k)+h(i,j,k)));
                        v2(i,j,k)=vdyfa(xant(i,j-1,k),xant(i,j,k),u2(i,j,k),hmy)      &
                                 +vcorra(u2(i,j,k),                                   &
                                  u1(i,j-1,k)+u1(i  ,j  ,k)+                          &
                                  u1(i+1,j,k)+u1(i+1,j-1,k),                          &
                                  (abs(xant(i+1,j-1,k))+abs(xant(i+1,j,k)))           &
                                 -(abs(xant(i-1,j-1,k))+abs(xant(i-1,j,k))),          &
                                   abs(xant(i+1,j-1,k))+abs(xant(i+1,j,k))            &
                                  +abs(xant(i-1,j-1,k))+abs(xant(i-1,j,k))+ep,        &
                                   hmy);   
                        f2ns(i,j,k)=donor(xant(i,j-1,k),xant(i,j,k),v2(i,j,k));
)
                 ENDIF !gndedge
InnerZFullXYDomainLoopDC(
                        hmx= 1._euwp/(0.5_euwp*(h(i-1,j,k)+h(i,j,k)));
                        v1(i,j,k)=vdyfa(xant(i-1,j,k),xant(i,j,k),u1(i,j,k),hmx)      &
                                 +vcorra(u1(i  ,j,k),                                 &
                                  u2(i-1,j  ,k)+u2(i-1,j+1,k)+                        &
                                  u2(i  ,j+1,k)+u2(i  ,j  ,k),                        &
                                  (abs(xant(i-1,j+1,k))+abs(xant(i,j+1,k)))           &
                                 -(abs(xant(i-1,j-1,k))+abs(xant(i,j-1,k))),          &
                                   abs(xant(i-1,j+1,k))+abs(xant(i,j+1,k))            &
                                  +abs(xant(i-1,j-1,k))+abs(xant(i,j-1,k))+ep,        &
                                  hmx)                                                &
                                 +vcorra(u1(i  ,j,k  ),                               &
                                  u3(i-1,j,k  )+u3(i-1,j,k+1)+                        &
                                  u3(i  ,j,k+1)+u3(i  ,j,k  ),                        &
                                  (abs(xant(i-1,j,k+1))+abs(xant(i,j,k+1)))           &
                                 -(abs(xant(i-1,j,k-1))+abs(xant(i,j,k-1))),          &
                                   abs(xant(i-1,j,k+1))+abs(xant(i,j,k+1))            &
                                  +abs(xant(i-1,j,k-1))+abs(xant(i,j,k-1))+ep,        &
                            hmx);
                        f1ns(i,j,k)=donor(xant(i-1,j,k),xant(i,j,k),v1(i,j,k));
                        hmy= 1._euwp/(0.5_euwp*(h(i,j-1,k)+h(i,j,k)));
                        v2(i,j,k)=vdyfa(xant(i,j-1,k),xant(i,j,k),u2(i,j,k),hmy)      &
                                 +vcorra(u2(i,j,k),                                   &
                                  u1(i,j-1,k)+u1(i  ,j  ,k)+                          &
                                  u1(i+1,j,k)+u1(i+1,j-1,k),                          &
                                  (abs(xant(i+1,j-1,k))+abs(xant(i+1,j,k)))           &
                                 -(abs(xant(i-1,j-1,k))+abs(xant(i-1,j,k))),          &
                                   abs(xant(i+1,j-1,k))+abs(xant(i+1,j,k))            &
                                  +abs(xant(i-1,j-1,k))+abs(xant(i-1,j,k))+ep,        &
                                   hmy)                                               &
                                 +vcorra(u2(i,j,k),                                   &
                                  u3(i,j-1,k  )+u3(i,j  ,k  )+                        &
                                  u3(i,j  ,k+1)+u3(i,j-1,k+1),                        &
                                  (abs(xant(i,j-1,k+1))+abs(xant(i,j,k+1)))           &
                                 -(abs(xant(i,j-1,k-1))+abs(xant(i,j,k-1))),          &
                                   abs(xant(i,j-1,k+1))+abs(xant(i,j,k+1))            &
                                  +abs(xant(i,j-1,k-1))+abs(xant(i,j,k-1))+ep,        &
                                   hmy);
                        f2ns(i,j,k)=donor(xant(i,j-1,k),xant(i,j,k),v2(i,j,k));
                        hmz= 1._euwp/(0.5_euwp*(h(i,j,k-1)+h(i,j,k)));
                        v3(i,j,k)=vdyfa(xant(i,j,k-1),xant(i,j,k),u3(i,j,k),hmz)      &
                                 +vcorra(u3(i,j,k),                                   &
                                  u1(i  ,j,k-1)+u1(i  ,j,k  )+                        &
                                  u1(i+1,j,k  )+u1(i+1,j,k-1),                        &
                                 (abs(xant(i+1,j,k-1))+abs(xant(i+1,j,k)))            &
                                -(abs(xant(i-1,j,k-1))+abs(xant(i-1,j,k))),           &
                                  abs(xant(i+1,j,k-1))+abs(xant(i+1,j,k))             &
                                 +abs(xant(i-1,j,k-1))+abs(xant(i-1,j,k))+ep,         &
                                  hmz)                                                &
                                 +vcorra(u3(i,j,k),                                   &
                                  u2(i,j,k-1)+u2(i,j+1,k-1)+                          &
                                  u2(i,j+1,k)+u2(i,j,k),                              &
                                 (abs(xant(i,j+1,k-1))+abs(xant(i,j+1,k)))            &
                                -(abs(xant(i,j-1,k-1))+abs(xant(i,j-1,k))),           &
                                  abs(xant(i,j+1,k-1))+abs(xant(i,j+1,k))             &
                                 +abs(xant(i,j-1,k-1))+abs(xant(i,j-1,k))+ep,         &
                                  hmz);
                        f3ns(i,j,k)=donor(xant(i,j,k-1),xant(i,j,k),v3(i,j,k));
)
            IF(skyedge == 1) THEN
              k=lp
FullXYDomainLoopDC(
                        hmx= 1._euwp/(0.5_euwp*(h(i-1,j,k)+h(i,j,k)));
                        v1(i,j,k)=vdyfa(xant(i-1,j,k),xant(i,j,k),u1(i,j,k),hmx)      &
                                 +vcorra(u1(i  ,j,k),                                 &
                                  u2(i-1,j  ,k)+u2(i-1,j+1,k)+                        &
                                  u2(i  ,j+1,k)+u2(i  ,j  ,k),                        &
                                  (abs(xant(i-1,j+1,k))+abs(xant(i,j+1,k)))           &
                                 -(abs(xant(i-1,j-1,k))+abs(xant(i,j-1,k))),          &
                                   abs(xant(i-1,j+1,k))+abs(xant(i,j+1,k))            &
                                  +abs(xant(i-1,j-1,k))+abs(xant(i,j-1,k))+ep,        &
                                  hmx);    
                        f1ns(i,j,k)=donor(xant(i-1,j,k),xant(i,j,k),v1(i,j,k));
                        hmy= 1._euwp/(0.5_euwp*(h(i,j-1,k)+h(i,j,k)));
                        v2(i,j,k)=vdyfa(xant(i,j-1,k),xant(i,j,k),u2(i,j,k),hmy)      &
                                 +vcorra(u2(i,j,k),                                   &
                                  u1(i,j-1,k)+u1(i  ,j  ,k)+                          &
                                  u1(i+1,j,k)+u1(i+1,j-1,k),                          &
                                  (abs(xant(i+1,j-1,k))+abs(xant(i+1,j,k)))           &
                                 -(abs(xant(i-1,j-1,k))+abs(xant(i-1,j,k))),          &
                                   abs(xant(i+1,j-1,k))+abs(xant(i+1,j,k))            &
                                  +abs(xant(i-1,j-1,k))+abs(xant(i-1,j,k))+ep,        &
                                   hmy);           
                        f2ns(i,j,k)=donor(xant(i,j-1,k),xant(i,j,k),v2(i,j,k));
                        hmz= 1._euwp/(0.5_euwp*(h(i,j,k-1)+h(i,j,k)));
                        v3(i,j,k)=vdyfa(xant(i,j,k-1),xant(i,j,k),u3(i,j,k),hmz)      &
                                 +vcorra(u3(i,j,k),                                   &
                                  u1(i  ,j,k-1)+u1(i  ,j,k  )+                        &
                                  u1(i+1,j,k  )+u1(i+1,j,k-1),                        &
                                  (abs(xant(i+1,j,k-1))+abs(xant(i+1,j,k)))           &
                                 -(abs(xant(i-1,j,k-1))+abs(xant(i-1,j,k))),          &
                                  abs(xant(i+1,j,k-1))+abs(xant(i+1,j,k))             &
                                 +abs(xant(i-1,j,k-1))+abs(xant(i-1,j,k))+ep,         &
                                  hmz)                                                &
                                 +vcorra(u3(i,j,k),                                   &
                                  u2(i,j,k-1)+u2(i,j+1,k-1)+                          &
                                  u2(i,j+1,k)+u2(i,j,k),                              &
                                 (abs(xant(i,j+1,k-1))+abs(xant(i,j+1,k)))            &
                                -(abs(xant(i,j-1,k-1))+abs(xant(i,j-1,k))),           &
                                  abs(xant(i,j+1,k-1))+abs(xant(i,j+1,k))             &
                                 +abs(xant(i,j-1,k-1))+abs(xant(i,j-1,k))+ep,         &
                                  hmz);
                        f3ns(i,j,k)=donor(xant(i,j,k-1),xant(i,j,k),v3(i,j,k));
)
            ENDIF

          IF(ipoles0.eq.0) THEN
            IF(ibcx0.eq.0) CALL flux_x_bcd(v1  ,np,mp,lp,ih, 0)
            IF(ibcx0.eq.0) CALL flux_x_bcd(f1ns,np,mp,lp,ih, 0)
            IF(ibcy0.eq.0) CALL flux_y_bcd(v2  ,np,mp,lp,ih, 0)
            IF(ibcy0.eq.0) CALL flux_y_bcd(f2ns,np,mp,lp,ih, 0)
          ENDIF
            IF(ibcz0.eq.0)  CALL flux_z_bc_nohalod(v3,np,mp,lp,ih, 2)
            IF(ibcz0.eq.0)  CALL flux_z_bc_nohalod(f3ns,np,mp,lp,ih, 2)
        CALL ttend(1053)
      CALL updatelrd(f1ns,np+rightedge*(1-ipoles0),mp,lp,np+1,mp  ,lp,iupx*(1-ipoles0)+ipoles0,ih)
      CALL updatebtd(f2ns,np,mp+topedge              ,lp,np  ,mp+1,lp,iupy,ih)
!      IF(ibcz0.eq.1) CALL updategs(v3  ,np,mp,lp+skyedge              ,np,mp,lp+1,iupz,ih)
      CALL updategsd(f3ns,np,mp,lp+skyedge              ,np,mp,lp+1,iupz,ih)
!     IF(ibcx0.eq.1) CALL updatelr(v1  ,np+rightedge*(1-ipoles0),mp,lp,np+1,mp  ,lp,iupx*(1-ipoles0)+ipoles0,ih)
!     IF(ibcy0.eq.1) CALL updatebt(v2  ,np,mp+topedge              ,lp,np  ,mp+1,lp,iupy,ih)

        CALL ttbeg(1053)
      IF(ipoles0 == 0) THEN
!      if(ibcx0.eq.1.and. leftedge.eq.1) then 
!              v1(1   ,1:mp,1:lp)=  v1(-1  ,1:mp,1:lp)
!            f1ns(1   ,1:mp,1:lp)=f1ns(-1  ,1:mp,1:lp)
!      endif
!      if(ibcx0.eq.1.and.rightedge.eq.1) then
!              v1(np+1,1:mp,1:lp)=  v1(np+3,1:mp,1:lp)
!            f1ns(np+1,1:mp,1:lp)=f1ns(np+3,1:mp,1:lp)
!      endif
!      if(ibcy0.eq.1.and. botedge.eq.1) then
!              v2(1:np,1   ,1:lp)=  v2(1:np,-1  ,1:lp)
!            f2ns(1:np,1   ,1:lp)=f2ns(1:np,-1  ,1:lp)
!      endif
!      if(ibcy0.eq.1.and. topedge.eq.1)  then
!              v2(1:np,mp+1,1:lp)=  v2(1:np,mp+3,1:lp)
!            f2ns(1:np,mp+1,1:lp)=f2ns(1:np,mp+3,1:lp)
!      endif
      ELSE
!              v2(1:np,1   ,1:lp)=0._euwp
!            f2ns(1:np,1   ,1:lp)=0._euwp
!              v2(1:np,mp+1,1:lp)=0._euwp 
!            f2ns(1:np,mp+1,1:lp)=0._euwp

      ENDIF
       if(ibcz0.eq.1.and.gndedge.eq.1) then
!              f3ns(1:np,1:mp,   1)=f3ns(1:np,1:mp,  -1)
!                v3(1:np,1:mp,   1)=  v3(1:np,1:mp,  -1)
       endif
       if(ibcz0.eq.1.and.skyedge.eq.1) then
!              f3ns(1:np,1:mp,lp+1)=f3ns(1:np,1:mp,lp+3)
!                v3(1:np,1:mp,lp+1)=  v3(1:np,1:mp,lp+3)
       endif

            !                 non-osscilatory option
           CALL cp_last_to_halo_xy(x_in,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)
        if(iflip.eq.1.and.ipoles0.eq.1)                    &
          CALL flip_halo_poles(x_in,np,mp,lp,ih)

            IF(gndedge == 1) THEN 
                k=1
FullXYDomainLoopDC(
                        mxijk_o=max(x_in(i-1,j  ,k  ),                                &
                                    x_in(i  ,j  ,k  ), x_in(i+1,j  ,k  ),             &
                                    x_in(i  ,j-1,k  ), x_in(i  ,j+1,k  ),             &
                                                       x_in(i  ,j  ,k+1));
                        mnijk_o=min(x_in(i-1,j  ,k  ),                                &
                                    x_in(i  ,j  ,k  ), x_in(i+1,j  ,k  ),             &
                                    x_in(i  ,j-1,k  ), x_in(i  ,j+1,k  ),             &
                                                       x_in(i  ,j  ,k+1));
                        mxijk  =max(xant(i-1,j  ,k  ), xant(i  ,j  ,k  ),             &
                                    xant(i+1,j  ,k  ),        mxijk_o,                &
                                    xant(i  ,j-1,k  ), xant(i  ,j+1,k  ),             &
                                                       xant(i  ,j  ,k+1));
                        mnijk  =min(xant(i  ,j  ,k  ),        mnijk_o,                &
                                    xant(i-1,j  ,k  ), xant(i+1,j  ,k  ),             &
                                    xant(i  ,j-1,k  ), xant(i  ,j+1,k  ),             &
                                                       xant(i  ,j  ,k+1));
                        cp(i,j,k)=(mxijk-xant(i,j,k))*h(i,j,k)/                       &
                                  (pn(f1ns(i+1,j,k))+pp( f1ns(i,j,k))                 &
                                  +pn(f2ns(i,j+1,k))+pp( f2ns(i,j,k))                 &
                                  +pn(f3ns(i,j,k+1))+pp(-f3ns(i,j,k+1))+ep);

                        cn(i,j,k)=(xant(i,j,k)-mnijk)*h(i,j,k)/                       &
                                  (pp(f1ns(i+1,j,k))+pn( f1ns(i,j,k))                 &
                                  +pp(f2ns(i,j+1,k))+pn( f2ns(i,j,k))                 &
                                  +pp(f3ns(i,j,k+1))+pn(-f3ns(i,j,k+1))+ep);
)
            ENDIF

InnerZFullXYDomainLoopDC(
                        mxijk_o=max(x_in(i-1,j  ,k  ),                                &
                                    x_in(i  ,j  ,k  ), x_in(i+1,j  ,k  ),             &
                                    x_in(i  ,j-1,k  ), x_in(i  ,j+1,k  ),             &
                                    x_in(i  ,j  ,k-1), x_in(i  ,j  ,k+1));
                        mnijk_o=min(x_in(i-1,j  ,k  ),                                &
                                    x_in(i  ,j  ,k  ), x_in(i+1,j  ,k  ),             &
                                    x_in(i  ,j-1,k  ), x_in(i  ,j+1,k  ),             &
                                    x_in(i  ,j  ,k-1), x_in(i  ,j  ,k+1));
                        mxijk  =max(xant(i-1,j  ,k  ), xant(i  ,j  ,k  ),             &
                                    xant(i+1,j  ,k  ),        mxijk_o,                &
                                    xant(i  ,j-1,k  ), xant(i  ,j+1,k  ),             &
                                    xant(i  ,j  ,k-1), xant(i  ,j  ,k+1));
                        mnijk  =min(xant(i  ,j  ,k  ),        mnijk_o,                &
                                    xant(i-1,j  ,k  ), xant(i+1,j  ,k  ),             &
                                    xant(i  ,j-1,k  ), xant(i  ,j+1,k  ),             &
                                    xant(i  ,j  ,k-1), xant(i  ,j  ,k+1));

                        cp(i,j,k)=(mxijk-xant(i,j,k))*h(i,j,k)/                       &
                                  (pn(f1ns(i+1,j,k))+pp(f1ns(i,j,k))                  &
                                  +pn(f2ns(i,j+1,k))+pp(f2ns(i,j,k))                  &
                                  +pn(f3ns(i,j,k+1))+pp(f3ns(i,j,k))+ep);

                        cn(i,j,k)=(xant(i,j,k)-mnijk)*h(i,j,k)/                       &
                                  (pp(f1ns(i+1,j,k))+pn(f1ns(i,j,k))                  &
                                  +pp(f2ns(i,j+1,k))+pn(f2ns(i,j,k))                  &
                                  +pp(f3ns(i,j,k+1))+pn(f3ns(i,j,k))+ep);
)
            IF(skyedge == 1) THEN 
                k=lp
FullXYDomainLoopDC(
                        mxijk_o=max(x_in(i-1,j  ,k  ),                                &
                                    x_in(i  ,j  ,k  ), x_in(i+1,j  ,k  ),             &
                                    x_in(i  ,j-1,k  ), x_in(i  ,j+1,k  ),             &
                                    x_in(i  ,j  ,k-1)                   );
                        mnijk_o=min(x_in(i-1,j  ,k  ),                                &
                                    x_in(i  ,j  ,k  ), x_in(i+1,j  ,k  ),             &
                                    x_in(i  ,j-1,k  ), x_in(i  ,j+1,k  ),             &
                                    x_in(i  ,j  ,k-1)                   );
                        mxijk  =max(xant(i-1,j  ,k  ), xant(i  ,j  ,k  ),             &
                                    xant(i+1,j  ,k  ),        mxijk_o,                &
                                    xant(i  ,j-1,k  ), xant(i  ,j+1,k  ),             &
                                    xant(i  ,j  ,k-1)                   );
                        mnijk  =min(xant(i  ,j  ,k  ),        mnijk_o,                &
                                    xant(i-1,j  ,k  ), xant(i+1,j  ,k  ),             &
                                    xant(i  ,j-1,k  ), xant(i  ,j+1,k  ),             &
                                    xant(i  ,j  ,k-1)                   );
                        cp(i,j,k)=(mxijk-xant(i,j,k))*h(i,j,k)/                       &
                                  (pn(f1ns(i+1,j,k))+pp(f1ns(i,j,k))                  &
                                  +pn(f2ns(i,j+1,k))+pp(f2ns(i,j,k))                  &
                                  +pn( -f3ns(i,j,k))+pp(f3ns(i,j,k))+ep);

                        cn(i,j,k)=(xant(i,j,k)-mnijk)*h(i,j,k)/                       &
                                  (pp(f1ns(i+1,j,k))+pn(f1ns(i,j,k))                  &
                                  +pp(f2ns(i,j+1,k))+pn(f2ns(i,j,k))                  &
                                  +pp( -f3ns(i,j,k))+pn(f3ns(i,j,k))+ep);
)
            ENDIF
#ifdef PNETCDF
            !    call pnet_out_chunk('f1_mpdta3','tsta1.nc',2,1,0,0,0,1,f1ns,np,mp,lp,ih)
            !    call pnet_out_chunk('f2_mpdta3','tsta2.nc',2,0,1,0,0,1,f2ns,np,mp,lp,ih)
            !    call pnet_out_chunk('f3_mpdta3','tsta3.nc',2,0,0,1,0,1,f3ns,np,mp,lp,ih)
            !    CALL pnet_out_chunk('cp       ','mpdat.nc',1,1,1,1,0,0,cp,np,mp,lp,ih)
            !    CALL pnet_out_chunk('cn       ','mpdat.nc',1,1,1,1,0,0,cn,np,mp,lp,ih)
            !    CALL pnet_out_chunk('xant     ','mpdat.nc',1,1,1,1,0,0,xant,np,mp,lp,ih)
#endif /*PNETCDF*/

            !
        CALL ttend(1053,.TRUE.)
!       IF(lupdatemulti) THEN
!         CALL update3_multi(np,mp,lp,np,mp,lp,iup,ih,cp,cn)
!       ELSE
          CALL update3d(cp,np,mp,lp,np,mp,lp,iup,ih)
          CALL update3d(cn,np,mp,lp,np,mp,lp,iup,ih)
!       ENDIF
        CALL ttbeg(1053)
!       CALL remove_cyclic_offset_full(cp,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
!       CALL remove_cyclic_offset_full(cn,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)

        CALL zero_halo_xyd(cp,ibcx0,ibcy0,np,mp,lp,ih) !Only works for open bc`s
        CALL zero_halo_xyd(cn,ibcx0,ibcy0,np,mp,lp,ih)
        InnerX_CGRID_POLES_FullYZDomainLoopDC(
                        a1=min(1.0_euwp,cp(i,j,k),cn(i-1,j,k));
                        a2=min(1.0_euwp,cp(i-1,j,k),cn(i,j,k));
                        v1ijk=                                                        &
                             pp(v1(i,j,k))*(a1*pp(sign(1.0_euwp, xant(i-1,j,k)))      &
                                          + a2*pp(sign(1.0_euwp,-xant(i-1,j,k))))     &
                            -pn(v1(i,j,k))*(a2*pp(sign(1.0_euwp, xant(i  ,j,k )))     &
                                          + a1*pp(sign(1.0_euwp,-xant(i  ,j,k ))));
                        f1(i,j,k)=donor(xant(i-1,j,k)/rhr(i-1,j,k),                   &
                                        xant(i  ,j,k)/rhr(i  ,j,k),                   &
                                        v1ijk);
)
         InnerY_CGRID_FullXZDomainLoopDC(
                        a1=min(1.0_euwp,cp(i,j,k),cn(i,j-1,k));
                        a2=min(1.0_euwp,cp(i,j-1,k),cn(i,j,k));
                        v2ijk=                                                        &
                             pp(v2(i,j,k))*(a1*pp(sign(1.0_euwp, xant(i,j-1,k)))      &
                                          + a2*pp(sign(1.0_euwp,-xant(i,j-1,k))) )    &
                            -pn(v2(i,j,k))*(a2*pp(sign(1.0_euwp, xant(i,j  ,k)))      &
                                          + a1*pp(sign(1.0_euwp,-xant(i,j  ,k))) );
                        f2(i,j,k)=donor(xant(i,j-1,k)/rhr(i,j-1,k),                   &
                                        xant(i,j  ,k)/rhr(i,j  ,k),                   &
                                        v2ijk);
)
          LInnerZFullXYDomainLoopDC(
                        a1=min(1.0_euwp,cp(i,j,k),cn(i,j,k-1));
                        a2=min(1.0_euwp,cp(i,j,k-1),cn(i,j,k));
                        v3ijk=                                                        &
                             pp(v3(i,j,k))*(a1*pp(sign(1.0_euwp, xant(i,j,k-1)))      &
                                          + a2*pp(sign(1.0_euwp,-xant(i,j,k-1))))     &
                            -pn(v3(i,j,k))*(a2*pp(sign(1.0_euwp, xant(i,j,k  )))      &
                                          + a1*pp(sign(1.0_euwp,-xant(i,j,k  ))));
                        f3(i,j,k)=donor(xant(i,j,k-1)/rhr(i,j,k-1),                   &
                                        xant(i,j,k  )/rhr(i,j,k  ) ,                  &
                                        v3ijk);
)
      IF(ipoles0 == 0) THEN
            IF(ibcx0.eq.0) CALL flux_x_bcd(f1,np,mp,lp,ih, 0)
            IF(ibcy0.eq.0) CALL flux_y_bcd(f2,np,mp,lp,ih, 0)
      ENDIF
            IF(ibcz0.eq.0) CALL flux_z_bc_nohalod(f3,np,mp,lp,ih, 2)

        CALL ttend(1053,.TRUE.)
            CALL updatelrd(f1,np+rightedge*(1-ipoles0),mp,lp,np+1,mp  ,lp  ,iupx*(1-ipoles0)+ipoles0,ih)
            CALL updatebtd(f2,np,mp+topedge,lp  ,np  ,mp+1,lp  ,iupy,ih)
            CALL updategsd(f3,np,mp,lp+skyedge  ,np  ,mp  ,lp+1,iupz,ih)
        CALL ttbeg(1053)
!     IF(ipoles0 == 0) THEN
!      if(ibcx0.eq.1.and. leftedge.eq.1) f1(1   ,1:mp,1:lp)=f1(-1  ,1:mp,1:lp)
!      if(ibcx0.eq.1.and.rightedge.eq.1) f1(np+1,1:mp,1:lp)=f1(np+3,1:mp,1:lp)
!      if(ibcy0.eq.1.and. botedge.eq.1)  f2(1:np,1   ,1:lp)=f2(1:np,-1  ,1:lp)
!      if(ibcy0.eq.1.and. topedge.eq.1)  f2(1:np,mp+1,1:lp)=f2(1:np,mp+3,1:lp)
!     ELSE
!      if(ibcy0.eq.1.and. botedge.eq.1)  f2(1:np,1   ,1:lp)=0.
!      if(ibcy0.eq.1.and. topedge.eq.1)  f2(1:np,mp+1,1:lp)=0.
!     ENDIF
!      if(ibcz0.eq.1.and.gndedge.eq.1)   f3(1:np,1:mp,   1)=  f3(1:np,1:mp,  -1)
!      if(ibcz0.eq.1.and.skyedge.eq.1)   f3(1:np,1:mp,lp+1)=  f3(1:np,1:mp,lp+3)

            IF(gndedge == 1) THEN
                k=1
FullXYDomainLoopDC(
                        hinv=1._euwp/h(i,j,k);
                        x_in(i,j,k)=rhr(i,j,k)*(xant(i,j,k)/rhr(i,j,k)-             &
                                   ( f1(i+1,j,k)-f1(i,j,k)                          &
                                    +f2(i,j+1,k)-f2(i,j,k)                          &
                                    +f3(i,j,k+1)+f3(i,j,k+1) )*hinv);
)
            ENDIF !gndedge

InnerZFullXYDomainLoopDC(
                        hinv=1._euwp/h(i,j,k);
                        x_in(i,j,k)=rhr(i,j,k)*(xant(i,j,k)/rhr(i,j,k)-             &
                                   ( f1(i+1,j,k)-f1(i,j,k)                          &
                                    +f2(i,j+1,k)-f2(i,j,k)                          &
                                    +f3(i,j,k+1)-f3(i,j,k) )*hinv);
)

            IF(skyedge == 1) THEN
                k=lp
FullXYDomainLoopDC(
                        hinv=1._euwp/h(i,j,k);
                        x_in(i,j,k)=rhr(i,j,k)*(xant(i,j,k)/rhr(i,j,k)-              &
                                    ( f1(i+1,j,k)-f1(i,j,k)                          &
                                     +f2(i,j+1,k)-f2(i,j,k)                          &
                                     -f3(i,j,k)  -f3(i,j,k) )*hinv);
)
            ENDIF
        CALL ttend(1053,.TRUE.)
        CALL enforce_cyclicg(x_in,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
        CALL update(x_in,np,mp,lp,np,mp,lp,iup,ih)
        CALL ttend(53)
    END SUBROUTINE antidiff3d_standard_gpubc
END MODULE module_antidiff3d_standard_gpubc

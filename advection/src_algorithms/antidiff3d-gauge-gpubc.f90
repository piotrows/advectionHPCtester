#include  "renames.inc"
MODULE module_antidiff3d_gauge_gpubc
   USE precisions, ONLY  : iintegers,euwp
   USE mpi_parallel, ONLY: leftedge,rightedge,botedge,topedge,gndedge,skyedge
#ifdef CUDACODE
   USE mpi_parallel, ONLY: istream1, istream2
   USE cudafor
#endif
   USE epsilons, ONLY: ep => ep_nonos
   USE mpdataoperators

   IMPLICIT NONE
CONTAINS
#include "renames.inc"
#include "defines.inc"
SUBROUTINE antidiff3d_gauge_gpubc(lupdatemulti,u1,u2,u3,xant,x_in,rhr,h,& 
                                   iflip,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih, &
                                           do_serialization_out )
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY:  update,updatelr,updatebt,updategs,update3
        USE mpi_parallel, ONLY:  updatelrd,updatebtd,updategsd,update3d,updated
        USE mpi_parallel, ONLY:  update3_multi
        USE mpi_parallel, ONLY: ttbeg,ttend
        USE mpi_parallel, ONLY: iup,iupx,iupy,iupz
        USE mpi_parallel, ONLY: enforce_cyclicg =>enforce_cyclic
        USE scratch_datafields, ONLY: v1,v2,v3,cp,cn
        LOGICAL, INTENT(IN) :: lupdatemulti
        LOGICAL,OPTIONAL :: do_serialization_out
        INTEGER(KIND=iintegers), INTENT(IN) :: iflip,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih
#ifdef CUDACODE
        REAL(KIND=euwp), INTENT(IN),MANAGED ::  &
                          u1(1-ih:np+1+ih, 1-ih:mp+ih, 1-ih:lp+ih),    &
                          u2(1-ih:np+ih, 1-ih:mp+1+ih, 1-ih:lp+ih),    &
                          u3(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih),    &
                           h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),    &
                         rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

        REAL(KIND=euwp), INTENT(INOUT),MANAGED ::                            &
                           x_in(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
        REAL(KIND=euwp), INTENT(INOUT),DEVICE ::                            &
                           xant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
#else
        REAL(KIND=euwp), INTENT(IN)::  &
                          u1(1-ih:np+1+ih, 1-ih:mp+ih, 1-ih:lp+ih),    &
                          u2(1-ih:np+ih, 1-ih:mp+1+ih, 1-ih:lp+ih),    &
                          u3(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih),    &
                           h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),    &
                         rhr(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

        REAL(KIND=euwp), INTENT(INOUT)::                            &
                           x_in(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),  &
                           xant(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
#endif

        INTEGER(KIND=iintegers) :: &
            i,j,k
        INTEGER :: istat 

        REAL(KIND=euwp) :: hmx,hmy,hmz

        REAL(KIND=euwp) :: mxijk,mnijk,mxijk_o,mnijk_o

        REAL(KIND=euwp) :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp
        CALL ttbeg(51)
        CALL ttbeg(1051)
        IF(PRESENT(do_serialization_out)) THEN
!Do nothing but remove compiler warning about unused dummy variables
        ENDIF

        !Full update needed as x is referenced in all edges
        !@cuf istat=cudaStreamSynchronize(istream1)

        CALL updated(xant,np,mp,lp,np,mp,lp,iup,ih)
!!        CALL remove_cyclic_offset_full(xant,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
        !At the domain edges, copy the second last  meaningful domain point to x halo
        !This is needed to zero the corrective fluxes vcorr
        !@cuf istat=cudaStreamSynchronize(istream2)
        CALL cp_scnd_last_to_halo_xy_fulld(xant,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)
!       if(iflip.eq.1.and.ipoles0.eq.1)                    &
!        CALL flip_halo_poles(xant,np,mp,lp,ih)


        !For the first antidiffusive iteration, f1,f2,f3 in traditional MPDATA
        !can be effectively replaced with C-grid advective velocities u1,u2,u3
        !This would not be valid for further antidiffusive iterations.
        !  rat2(z1,z2)=(z2-z1)*.5_euwp
        !  vdyf(x1,x2,a,rinv)=(abs(a)-a**2*rinv)*rat2(x1,x2)
        !  rat4(z0,z1,z2,z3)=(z3+z2-z1-z0)*.25_euwp
        !  vcorr(a,b,y0,y1,y2,y3,rinv)=-0.125_euwp*a*b*rinv*rat4(y0,y1,y2,y3)

        IF(gndedge == 1 ) THEN
            k=1
         FullXYDomainLoopDCs(stream=istream2,
                    !----------------------------------------------;
                    !compute antidiffusive velocities in x direction;
                    !----------------------------------------------;
                    hmx= 1._euwp/(0.5_euwp*(h(i-1,j,k)+h(i,j,k)));
                     v1(i,j,k)=vdyf(xant(i-1,j,k),xant(i,j,k),u1(i,j,k),hmx)    &
                          -0.125_euwp*u1(i  ,j,k)*hmx*                          &
                                    ((u2(i-1,j,k)+u2(i-1,j+1,k)                 &
                                     +u2(i,j+1,k)+u2(i  ,j  ,k) )               &
                              *rat4(xant(i-1,j-1,k),xant(i,j-1,k),              &
                                    xant(i-1,j+1,k),xant(i,j+1,k)));
                    !----------------------------------------------
                    !compute antidiffusive velocities in y direction
                    !----------------------------------------------;
                    hmy= 1._euwp/(0.5_euwp*(h(i,j-1,k)+h(i,j,k)));
!                    v2(i,j,k)=eval_v2_a(xant,u1,u2,u3,hmy,i,j,k,np,mp,lp,ih) ;
                     v2(i,j,k)=vdyf(xant(i,j-1,k),xant(i,j,k),u2(i,j,k),hmy)    &
                          -0.125_euwp*u2(i,j,k)*hmy*                            &
                                    ((u1(i  ,j-1,k)+u1(i  ,j  ,k)               &
                                     +u1(i+1,j  ,k)+u1(i+1,j-1,k))              &
                              *rat4(xant(i-1,j-1,k),xant(i-1,j,k),              &
                                    xant(i+1,j-1,k),xant(i+1,j,k)));
                     v3(i,j,k)=0._euwp ;
)
        ENDIF
        IF(skyedge == 1 ) THEN
            k=lp
         FullXYDomainLoopDCs(stream=istream2,

                    !----------------------------------------------
                    !compute antidiffusive velocities in x direction
                    !----------------------------------------------;
                    hmx= 1._euwp/(0.5_euwp*(h(i-1,j,k)+h(i,j,k)));
!                    v1(i,j,k)=eval_v1_a(xant,u1,u2,u3,hmx,i,j,k,np,mp,lp,ih) ;
                     v1(i,j,k)=vdyf(xant(i-1,j,k),xant(i,j,k),u1(i,j,k),hmx)    &
                          -0.125_euwp*u1(i  ,j,k)*hmx*                          &
                                    ((u2(i-1,j,k)+u2(i-1,j+1,k)                 &
                                     +u2(i,j+1,k)+u2(i  ,j  ,k) )               &
                              *rat4(xant(i-1,j-1,k),xant(i,j-1,k),              &
                                    xant(i-1,j+1,k),xant(i,j+1,k)));
                    !----------------------------------------------
                    !compute antidiffusive velocities in y direction
                    !----------------------------------------------;
                    hmy= 1._euwp/(0.5_euwp*(h(i,j-1,k)+h(i,j,k)));
!                    v2(i,j,k)=eval_v2_a(xant,u1,u2,u3,hmy,i,j,k,np,mp,lp,ih)  ;
                     v2(i,j,k)=vdyf(xant(i,j-1,k),xant(i,j,k),u2(i,j,k),hmy)    &
                          -0.125_euwp*u2(i,j,k)*hmy*                            &
                                    ((u1(i  ,j-1,k)+u1(i  ,j  ,k)               &
                                     +u1(i+1,j  ,k)+u1(i+1,j-1,k))              &
                              *rat4(xant(i-1,j-1,k),xant(i-1,j,k),              &
                                    xant(i+1,j-1,k),xant(i+1,j,k)));
                    hmz= 1._euwp/(0.5_euwp*(h(i,j,k-1)+h(i,j,k)));
                     v3(i,j,k)=vdyf(xant(i,j,k-1),xant(i,j,k),u3(i,j,k),hmz)    &
                          -0.125_euwp*u3(i,j,k)*hmz*                            &
                                    ((u1(i  ,j,k-1)+u1(i  ,j,k  )               &
                                     +u1(i+1,j,k  )+u1(i+1,j,k-1))              &
                              *rat4(xant(i-1,j,k-1), xant(i-1,j,k  ),           &
                                    xant(i+1,j,k-1), xant(i+1,j,k  ))           &
                 +                                                              &
                                     (u2(i,j  ,k-1)+u2(i,j+1,k-1)               &
                                     +u2(i,j+1,k  )+u2(i,j  ,k  ))              &
                              *rat4(xant(i,j-1,k-1), xant(i,j-1,k  ),           &
                                    xant(i,j+1,k-1), xant(i,j+1,k  )));
)
        ENDIF

     InnerZFullXYDomainLoopDCs(stream=istream1,
                    !----------------------------------------------
                    !compute antidiffusive velocities in x direction
                    !----------------------------------------------;
                    hmx= 1._euwp/(0.5_euwp*(h(i-1,j,k)+h(i,j,k)));
                     v1(i,j,k)=vdyf(xant(i-1,j,k),xant(i,j,k),u1(i,j,k),hmx)    &
                          -0.125_euwp*u1(i  ,j,k)*hmx*                          &
                                    ((u2(i-1,j,k)+u2(i-1,j+1,k)                 &
                                     +u2(i,j+1,k)+u2(i  ,j  ,k) )               &
                              *rat4(xant(i-1,j-1,k),xant(i,j-1,k),              &
                                    xant(i-1,j+1,k),xant(i,j+1,k))              &
                   +                                                            &
                                     (u3(i-1,j,k  )+u3(i-1,j,k+1)               &
                                     +u3(i  ,j,k+1)+u3(i  ,j,k  ))              &
                              *rat4(xant(i-1,j,k-1),xant(i,j,k-1),              &
                                    xant(i-1,j,k+1),xant(i,j,k+1)));
                    !----------------------------------------------
                    !compute antidiffusive velocities in y direction
                    !----------------------------------------------;
                    hmy= 1._euwp/(0.5_euwp*(h(i,j-1,k)+h(i,j,k)));
                     v2(i,j,k)=vdyf(xant(i,j-1,k),xant(i,j,k),u2(i,j,k),hmy)    &
                          -0.125_euwp*u2(i,j,k)*hmy*                            &
                                    ((u1(i  ,j-1,k)+u1(i  ,j  ,k)               &
                                     +u1(i+1,j  ,k)+u1(i+1,j-1,k))              &
                              *rat4(xant(i-1,j-1,k),xant(i-1,j,k),              &
                                    xant(i+1,j-1,k),xant(i+1,j,k))              &
                  +                                                             &
                                     (u3(i,j-1,k  )+u3(i,j  ,k  )               &
                                     +u3(i,j  ,k+1)+u3(i,j-1,k+1))              &
                              *rat4(xant(i,j-1,k-1), xant(i,j  ,k-1),           &
                                    xant(i,j-1,k+1), xant(i,j  ,k+1)));
                    !----------------------------------------------
                    !compute antidiffusive velocities in z direction
                    !----------------------------------------------;
                    hmz= 1._euwp/(0.5_euwp*(h(i,j,k-1)+h(i,j,k)));
                     v3(i,j,k)=vdyf(xant(i,j,k-1),xant(i,j,k),u3(i,j,k),hmz)    &
                          -0.125_euwp*u3(i,j,k)*hmz*                            &
                                    ((u1(i  ,j,k-1)+u1(i  ,j,k  )               &
                                     +u1(i+1,j,k  )+u1(i+1,j,k-1))              &
                              *rat4(xant(i-1,j,k-1), xant(i-1,j,k  ),           &
                                    xant(i+1,j,k-1), xant(i+1,j,k  ))           &
                 +                                                              &
                                     (u2(i,j  ,k-1)+u2(i,j+1,k-1)               &
                                     +u2(i,j+1,k  )+u2(i,j  ,k  ))              &
                              *rat4(xant(i,j-1,k-1), xant(i,j-1,k  ),           &
                                    xant(i,j+1,k-1), xant(i,j+1,k  )));
                            )

        !                 non-osscilatory option
 !@cuf istat=cudaStreamSynchronize(istream2)
 !@cuf istat=cudaStreamSynchronize(istream1)
        !
      IF(ipoles0 == 0) THEN 
        IF(ibcx0.eq.0) CALL flux_x_bcd(v1  ,np,mp,lp,ih, 0)
        IF(ibcy0.eq.0) CALL flux_y_bcd(v2  ,np,mp,lp,ih, 0) 
      ENDIF
        IF(ibcz0.eq.0) CALL flux_z_bc_nohalod(v3,np,mp,lp,ih, 2)

        CALL cp_last_to_halo_xy(x_in,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)
!       if(iflip.eq.1.and.ipoles0.eq.1)                    &
!        CALL flip_halo_poles(x_in,np,mp,lp,ih) ;

        CALL ttend(1051)
        
        CALL updatelrd (v1,np+rightedge*(1-ipoles0),mp,lp,np+1,mp  ,lp,iupx*(1-ipoles0)+ipoles0,ih)
        CALL updatebtd (v2,np,mp+topedge  ,lp,np  ,mp+1,lp,iupy,ih)
        CALL updategsd (v3,np,mp,lp+skyedge,np,mp,lp+1,iupz,ih)
        CALL ttbeg(1051)
      IF(ipoles0 == 0) THEN 
!       if(ibcx0.eq.1.and. leftedge.eq.1) v1(1   ,1:mp,1:lp)=v1(-1  ,1:mp,1:lp)
!       if(ibcx0.eq.1.and.rightedge.eq.1) v1(np+1,1:mp,1:lp)=v1(np+3,1:mp,1:lp)
!       if(ibcy0.eq.1.and. botedge.eq.1)  v2(1:np,1   ,1:lp)=v2(1:np,-1  ,1:lp)
!       if(ibcy0.eq.1.and. topedge.eq.1)  v2(1:np,mp+1,1:lp)=v2(1:np,mp+3,1:lp)
      ELSE
!      if(botedge.eq.1) v2(1:np,1   ,1:lp)=0.
!      if(topedge.eq.1) v2(1:np,mp+1,1:lp)=0.
      ENDIF
!      if(ibcz0.eq.1.and.gndedge.eq.1)   v3(1:np,1:mp,   1)=  v3(1:np,1:mp,  -1)
!      if(ibcz0.eq.1.and.skyedge.eq.1)   v3(1:np,1:mp,lp+1)=  v3(1:np,1:mp,lp+3)


        IF(gndedge == 1) THEN
          k=1;
          FullXYDomainLoopDCs(stream=istream2,
              mxijk_o=max(x_in(i-1,j  ,k  ),                                &
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

              cp(i,j,k)=(mxijk - xant(i,j,k))*h(i,j,k)/                     &
                        (pn(v1(i+1,j,k))+pp( v1(i,j,k))                     &
                        +pn(v2(i,j+1,k))+pp( v2(i,j,k))                     &
                        +pn(v3(i,j,k+1))+pp(-v3(i,j,k+1))+ep);

              cn(i,j,k)=(xant(i,j,k) - mnijk)*h(i,j,k)/                      &
                        (pp(v1(i+1,j,k))+pn( v1(i,j,k))                      &
                        +pp(v2(i,j+1,k))+pn( v2(i,j,k))                      &
                        +pp(v3(i,j,k+1))+pn(-v3(i,j,k+1))+ep);
                )
        ENDIF !gndedge

        IF(skyedge == 1) THEN
          k=lp
          FullXYDomainLoopDCs(stream=istream2,
              mxijk_o=max(x_in(i-1,j  ,k  ),                                &
                          x_in(i  ,j  ,k  ), x_in(i+1,j  ,k  ),             &
                          x_in(i  ,j-1,k  ), x_in(i  ,j+1,k  ),             &
                          x_in(i  ,j  ,k-1));
              mnijk_o=min(x_in(i-1,j  ,k  ),                                &
                          x_in(i  ,j  ,k  ), x_in(i+1,j  ,k  ),             &
                          x_in(i  ,j-1,k  ), x_in(i  ,j+1,k  ),             &
                          x_in(i  ,j  ,k-1));
              mxijk  =max(xant(i-1,j  ,k  ), xant(i  ,j  ,k  ),             &
                          xant(i+1,j  ,k  ),        mxijk_o,                &
                          xant(i  ,j-1,k  ), xant(i  ,j+1,k  ),             &
                          xant(i  ,j  ,k-1));
              mnijk  =min(xant(i  ,j  ,k  ),        mnijk_o,                &
                          xant(i-1,j  ,k  ), xant(i+1,j  ,k  ),             &
                          xant(i  ,j-1,k  ), xant(i  ,j+1,k  ),             &
                          xant(i  ,j  ,k-1));

              cp(i,j,k)=(mxijk -xant(i,j,k))*h(i,j,k)/                      &
                        (pn(v1(i+1,j,k))+pp(v1(i,j,k))                      &
                        +pn(v2(i,j+1,k))+pp(v2(i,j,k))                      &
                        +pn( -v3(i,j,k))+pp(v3(i,j,k))+ep);

              cn(i,j,k)=(xant(i,j,k) -mnijk)*h(i,j,k)/                      &
                        (pp(v1(i+1,j,k))+pn(v1(i,j,k))                      &
                        +pp(v2(i,j+1,k))+pn(v2(i,j,k))                      &
                        +pp( -v3(i,j,k))+pn(v3(i,j,k))+ep);
                )
        ENDIF

        InnerZFullXYDomainLoopDCs(stream=istream1,
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
              cp(i,j,k)=(mxijk -xant(i,j,k))*h(i,j,k)/                      &
                        (pn(v1(i+1,j,k))+pp(v1(i,j,k))                      &
                        +pn(v2(i,j+1,k))+pp(v2(i,j,k))                      &
                        +pn(v3(i,j,k+1))+pp(v3(i,j,k))+ep);

              cn(i,j,k)=(xant(i,j,k) -mnijk)*h(i,j,k)/                      &
                        (pp(v1(i+1,j,k))+pn(v1(i,j,k))                      &
                        +pp(v2(i,j+1,k))+pn(v2(i,j,k))                      &
                        +pp(v3(i,j,k+1))+pn(v3(i,j,k))+ep);
                )
        CALL ttend(1051,.TRUE.)
        IF(lupdatemulti) THEN
!         CALL update3_multi(np,mp,lp,np,mp,lp,iup,ih,cp,cn)
        ELSE
         !@cuf istat=cudaStreamSynchronize(istream1)
          CALL update3d(cp,np,mp,lp,np,mp,lp,iup ,ih)
          CALL update3d(cn,np,mp,lp,np,mp,lp,iup ,ih)
        ENDIF
        CALL ttbeg(1051)
!       CALL remove_cyclic_offset_full(cp,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
!       CALL remove_cyclic_offset_full(cn,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)

       CALL zero_halo_xyd(cp,ibcx0,ibcy0,np,mp,lp,ih) !Only works for open bc`s
       CALL zero_halo_xyd(cn,ibcx0,ibcy0,np,mp,lp,ih)  

       IF(gndedge == 1) THEN
            k=1;
          FullXYDomainLoopDCs(stream=istream2,
                    f1ijk=                                                        &
                         pp(v1(i,j,k  ))*min(1.0_euwp,cp(i  ,j,k),cn(i-1,j,k))    &
                        -pn(v1(i,j,k  ))*min(1.0_euwp,cp(i-1,j,k),cn(i  ,j,k));
                    f1ijkp=                                                       &
                         pp(v1(i+1,j,k))*min(1.0_euwp,cp(i+1,j,k),cn(i  ,j,k))    &
                        -pn(v1(i+1,j,k))*min(1.0_euwp,cp(i  ,j,k),cn(i+1,j,k));
                    f2ijk=                                                        &
                         pp(v2(i,j,k  ))*min(1.0_euwp,cp(i,j  ,k),cn(i,j-1,k))    &
                        -pn(v2(i,j,k  ))*min(1.0_euwp,cp(i,j-1,k),cn(i,j  ,k));
                    f2ijkp=                                                       &
                         pp(v2(i,j+1,k))*min(1.0_euwp,cp(i,j+1,k),cn(i,j  ,k))    &
                        -pn(v2(i,j+1,k))*min(1.0_euwp,cp(i,j  ,k),cn(i,j+1,k));
!                    f3ijk=                                                        &
!                         pp(v3(i,j,k  ))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k-1))    &
!                        -pn(v3(i,j,k  ))*min(1.0_euwp,cp(i,j,k-1),cn(i,j,k  ));
                    f3ijk=                                                        &
                         pp(-v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k+1))    &
                        -pn(-v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k+1),cn(i,j,k  ));
                    f3ijkp=                                                       &
                         pp( v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k+1),cn(i,j,k  ))    &
                        -pn( v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k+1));
                    x_in(i,j,k)=rhr(i,j,k)*                          &
                              (xant(i,j,k)/rhr(i,j,k)-               &
                                        ( f1ijkp-f1ijk               &
                                         +f2ijkp-f2ijk               &
                                         +f3ijkp-f3ijk )/h(i,j,k));
                                 )
        ENDIF

        IF(skyedge == 1) THEN
            k=lp;
         FullXYDomainLoopDCs(stream=istream2,
                    f1ijk=                                                        &
                         pp(v1(i,j,k  ))*min(1.0_euwp,cp(i  ,j,k),cn(i-1,j,k))    &
                        -pn(v1(i,j,k  ))*min(1.0_euwp,cp(i-1,j,k),cn(i  ,j,k));
                    f1ijkp=                                                       &
                         pp(v1(i+1,j,k))*min(1.0_euwp,cp(i+1,j,k),cn(i  ,j,k))    &
                        -pn(v1(i+1,j,k))*min(1.0_euwp,cp(i  ,j,k),cn(i+1,j,k));
                    f2ijk=                                                        &
                         pp(v2(i,j,k  ))*min(1.0_euwp,cp(i,j  ,k),cn(i,j-1,k))    &
                        -pn(v2(i,j,k  ))*min(1.0_euwp,cp(i,j-1,k),cn(i,j  ,k));
                    f2ijkp=                                                       &
                         pp(v2(i,j+1,k))*min(1.0_euwp,cp(i,j+1,k),cn(i,j  ,k))    &
                        -pn(v2(i,j+1,k))*min(1.0_euwp,cp(i,j  ,k),cn(i,j+1,k));
                    f3ijk=                                                        &
                         pp(v3(i,j,k  ))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k-1))    &
                        -pn(v3(i,j,k  ))*min(1.0_euwp,cp(i,j,k-1),cn(i,j,k  ));
!                    f3ijkp=                                                       &
!                         pp(v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k+1),cn(i,j,k  ))    &
!                        -pn(v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k+1));
                    f3ijkp=                                                       &
                         pp(-v3(i,j,k ))*min(1.0_euwp,cp(i,j,k-1),cn(i,j,k  ))    &
                        -pn(-v3(i,j,k ))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k-1));
                    x_in(i,j,k)=rhr(i,j,k)*                          &
                              (xant(i,j,k)/rhr(i,j,k)-               &
                                        ( f1ijkp-f1ijk               &
                                         +f2ijkp-f2ijk               &
                                         +f3ijkp-f3ijk )/h(i,j,k));
                                 )
        ENDIF

        InnerZFullXYDomainLoopDCs(stream=istream1,
                    f1ijk=                                                        &
                         pp(v1(i,j,k  ))*min(1.0_euwp,cp(i  ,j,k),cn(i-1,j,k))    &
                        -pn(v1(i,j,k  ))*min(1.0_euwp,cp(i-1,j,k),cn(i  ,j,k));
                    f1ijkp=                                                       &
                         pp(v1(i+1,j,k))*min(1.0_euwp,cp(i+1,j,k),cn(i  ,j,k))    &
                        -pn(v1(i+1,j,k))*min(1.0_euwp,cp(i  ,j,k),cn(i+1,j,k));
                    f2ijk=                                                        &
                         pp(v2(i,j,k  ))*min(1.0_euwp,cp(i,j  ,k),cn(i,j-1,k))    &
                        -pn(v2(i,j,k  ))*min(1.0_euwp,cp(i,j-1,k),cn(i,j  ,k));
                    f2ijkp=                                                       &
                         pp(v2(i,j+1,k))*min(1.0_euwp,cp(i,j+1,k),cn(i,j  ,k))    &
                        -pn(v2(i,j+1,k))*min(1.0_euwp,cp(i,j  ,k),cn(i,j+1,k));
                    f3ijk=                                                        &
                         pp(v3(i,j,k  ))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k-1))    &
                        -pn(v3(i,j,k  ))*min(1.0_euwp,cp(i,j,k-1),cn(i,j,k  ));
                    f3ijkp=                                                       &
                         pp(v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k+1),cn(i,j,k  ))    &
                        -pn(v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k+1));
                    x_in(i,j,k)=rhr(i,j,k)*                          &
                              (xant(i,j,k)/rhr(i,j,k)-               &
                                        ( f1ijkp-f1ijk               &
                                         +f2ijkp-f2ijk               &
                                         +f3ijkp-f3ijk )/h(i,j,k));
                                 )
        CALL ttend(1051,.TRUE.)
        CALL enforce_cyclicg(x_in,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
         !@cuf istat=cudaStreamSynchronize(istream2)
         !@cuf istat=cudaStreamSynchronize(istream1)
!       CALL update(x_in,np,mp,lp,np,mp,lp,iup,ih)
        CALL ttend(51)
    END SUBROUTINE antidiff3d_gauge_gpubc
END MODULE module_antidiff3d_gauge_gpubc

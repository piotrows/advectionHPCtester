#include  "../src_algorithms/renames.inc"
!#ifdef SERIALIZE
!   USE m_serialize, ONLY: fs_write_field, fs_create_savepoint, fs_read_field, &
!      fs_read_and_perturb_field, fs_add_savepoint_metainfo
!   USE utils_ppser, ONLY: ppser_get_mode, ppser_set_mode, ppser_savepoint, &
!      ppser_serializer, ppser_serializer_ref, ppser_intlength, ppser_reallength, &
!      ppser_realtype, ppser_zrperturb
!#endif
    USE precisions, ONLY  : iintegers,euwp
    USE mpi_parallel, ONLY: leftedge,rightedge,botedge,topedge,gndedge,skyedge
#ifdef CUDACODE 
    USE mpi_parallel, ONLY: istream1, istream2
#endif
#ifdef PNETCDF
    USE mpi_parallel, ONLY: pnet_out_chunk
#endif
    USE epsilons, ONLY: ep => ep_nonos
    USE noise, ONLY: perturb_2D_signal
    IMPLICIT NONE
CONTAINS
#include "defines.inc"
#include "mpdataoperators.inc"
SUBROUTINE mpdata3d_rho_gauge_gpubc(                                &
           lupdatemulti,nsubsteps,nsubsteps_o,                      &
           lsubstepping,                                            &
           dtn,dtnold,ltimeadapt,                                   &
           lvertsplit,lcrdiag_nphlf,crmaxlim,dtmax,                 &
        ox,oy,oz,u,v,w,u1,u2,u3,                                    &
        x_in,x_incomp,h,rhoadv,rhr,                                 &
        pexe,pexr1,pexr2,                                           &
        bcx,bcy,                                                    &
        ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih,dxi,dyi,dzi,dt,liner, &
        do_serialization_in,do_serialization_out )
        !---------------------------------------------------------------------!
        USE mpi_parallel, ONLY:  update,updatelr,updatebt,updategs,update3
        USE mpi_parallel, ONLY:  updated,updatelrd,updatebtd,updategsd,update3d
        USE mpi_parallel, ONLY: update_multi,update3_multi
        USE mpi_parallel, ONLY: ttbeg,ttend
        USE mpi_parallel, ONLY: mype,iup,iupx,iupy,iupz
        USE mpi_parallel, ONLY: globmax 
        USE scratch_datafields, ONLY: v1,v2,v3,cp,cn,xant,f1ns,f2ns,f3ns
        USE scratch_datafields, ONLY: pfx,pfy
        USE eulag_diagnostics, ONLY: diagnos_advection, diagnos_advection_vz
        USE eulag_diagnostics, ONLY: exam_var
        USE eulag_diagutils, ONLY: compute_courlipsch_Agrid
        USE eulag_diagutils, ONLY: compute_courlipsch_Agrid_full
        USE module_velprd, ONLY: velprd_traj0
        USE velprd_driver_module, ONLY:  velprd_driver
        LOGICAL, INTENT(IN) :: lupdatemulti,lsubstepping,ltimeadapt,lvertsplit,lcrdiag_nphlf
       LOGICAL, OPTIONAL :: do_serialization_in
       LOGICAL, OPTIONAL :: do_serialization_out
        INTEGER(KIND=iintegers), INTENT(IN) :: ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih
        INTEGER(KIND=iintegers), INTENT(INOUT) :: nsubsteps,nsubsteps_o,liner 
        REAL_euwp, INTENT(INOUT) ::                                     &
                          u(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),    &
                          v(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),    &
                          w(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)
        REAL_euwp, INTENT(INOUT) ::                                     &
                      ox(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),         &
                      oy(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),         &
                      oz(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2)

        REAL_euwp, INTENT(INOUT) ::  &
                        x_in(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),       &
                    x_incomp(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)     

        REAL_euwp, INTENT(OUT) ::  &
                       pexr1(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),       &
                       pexr2(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),       &
                      rhoadv(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),       &
                         rhr(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),       &
                          u1(1-ih:np+1+ih, 1-ih:mp+ih, 1-ih:lp+ih),     &
                          u2(1-ih:np+ih, 1-ih:mp+1+ih, 1-ih:lp+ih),     &
                          u3(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+1+ih)

        REAL_euwp, INTENT(IN) ::                                        &
                             h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),     &
                          pexe(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) 

        REAL_euwp, INTENT(IN) ::                                        &
                             bcx(mp,lp, 2),                             &      
                             bcy(np,lp, 2)        
        REAL(KIND=euwp),INTENT(IN) :: dxi,dyi,dzi,dt
        REAL(KIND=euwp), INTENT(INOUT) :: dtn,dtnold
        REAL(KIND=euwp),INTENT(IN) :: crmaxlim,dtmax

        INTEGER(KIND=iintegers) :: &
            i,j,k

        REAL(KIND=euwp) :: hmx,hmy,hmz

        REAL(KIND=euwp) ::           &
            mxijk,mnijk,mxijk_o,mnijk_o

        REAL(KIND=euwp) :: f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp
        REAL(KIND=euwp) :: cr1,cr2 
        REAL(KIND=euwp) :: subfac,dt_ratio

        CALL ttbeg(9)

        IF(PRESENT(do_serialization_in).OR.PRESENT(do_serialization_out)) THEN
!Do nothing but remove compiler warning about unused dummy variables
        ENDIF


      CALL  velprd_driver(nsubsteps,nsubsteps_o,lsubstepping, &
                       dtn,dtnold,ltimeadapt,                      &
                       lvertsplit,lcrdiag_nphlf,crmaxlim,dtmax,     &
                       ox,oy,oz,u,v,w,u1,u2,u3,h,  &
                       ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih,           &
                       dxi,dyi,dzi,dt)
        CALL ttbeg(1009)
!-------- Store ox to be later used later as n-1 

FullXYZDomainLoopDC(
      ox(i,j,k,2)=ox(i,j,k,0)*h(i,j,k);
      oy(i,j,k,2)=oy(i,j,k,0)*h(i,j,k);
      oz(i,j,k,2)=oz(i,j,k,0)*h(i,j,k);
)
        ! 0. Store computational density at time level n
FullXYZDomainLoopDC(
         rhoadv(i,j,k)=x_incomp(i,j,k);
)
        CALL putbcinxh(x_in,bcx,bcy,ibcx0,ibcy0,np,mp,lp,ih)

        CALL ttend(1009)
        IF(lupdatemulti) THEN
          CALL update(rhoadv,np,mp,lp,np,mp,lp,iup,ih)
          CALL update3(x_in,np,mp,lp,np,mp,lp,iup,ih)
        ELSE        
          CALL update_multi(np,mp,lp,np,mp,lp,iup,ih,rhoadv,x_in)
        ENDIF
        CALL remove_cyclic_offset_full(x_in,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)

        CALL ttbeg(1009)

      IF(gndedge == 1) THEN
         k=1
         FullXYDomainLoopDC(
              f1ijkp=donor(x_in(i  ,j,k), x_in(i+1,j,k), u1(i+1,j,  k  ));
              f1ijk =donor(x_in(i-1,j,k), x_in(i  ,j,k), u1(i  ,j  ,k  ));
              f2ijkp=donor(x_in(i,j  ,k), x_in(i,j+1,k), u2(i  ,j+1,k  ));
              f2ijk =donor(x_in(i,j-1,k), x_in(i,j  ,k), u2(i  ,j  ,k  ));
              f3ijkp=donor(x_in(i,j,k  ), x_in(i,j,k+1), u3(i  ,j  ,k+1));
              f3ijk =-f3ijkp;
              xant(i,j,k)=                                      &
                        (x_in(i,j,k)-                           &
                                     ( f1ijkp-f1ijk             &
                                      +f2ijkp-f2ijk             &
                                      +f3ijkp-f3ijk )/h(i,j,k));
              f1ns(i,j,k)=f1ijk;
              f2ns(i,j,k)=f2ijk;
              f3ns(i,j,k)=f3ijk;
)
        ENDIF !gndedge
      InnerZFullXYDomainLoopDC(
              f1ijkp=donor(x_in(i  ,j,k), x_in(i+1,j,k), u1(i+1,j,  k  ));
              f1ijk =donor(x_in(i-1,j,k), x_in(i  ,j,k), u1(i  ,j  ,k  ));
              f2ijkp=donor(x_in(i,j  ,k), x_in(i,j+1,k), u2(i  ,j+1,k  ));
              f2ijk =donor(x_in(i,j-1,k), x_in(i,j  ,k), u2(i  ,j  ,k  ));
              f3ijkp=donor(x_in(i,j,k  ), x_in(i,j,k+1), u3(i  ,j  ,k+1));
              f3ijk =donor(x_in(i,j,k-1), x_in(i,j,k  ), u3(i  ,j  ,k  ));
              xant(i,j,k)=                                      &
                        (x_in(i,j,k)-                           &
                                     ( f1ijkp-f1ijk             &
                                      +f2ijkp-f2ijk             &
                                      +f3ijkp-f3ijk )/h(i,j,k));
              f1ns(i,j,k)=f1ijk;
              f2ns(i,j,k)=f2ijk;
              f3ns(i,j,k)=f3ijk;
)

        IF(skyedge == 1) THEN
          k=lp
         FullXYDomainLoopDC(
              f1ijkp=donor(x_in(i  ,j,k), x_in(i+1,j,k), u1(i+1,j,  k  ));
              f1ijk =donor(x_in(i-1,j,k), x_in(i  ,j,k), u1(i  ,j  ,k  ));
              f2ijkp=donor(x_in(i,j  ,k), x_in(i,j+1,k), u2(i  ,j+1,k  ));
              f2ijk =donor(x_in(i,j-1,k), x_in(i,j  ,k), u2(i  ,j  ,k  ));
              f3ijk =donor(x_in(i,j,k-1), x_in(i,j,k  ), u3(i  ,j  ,k  ));
              f3ijkp=-f3ijk;
              xant(i,j,k)=                                      &
                        (x_in(i,j,k)-                           &
                                     ( f1ijkp-f1ijk             &
                                      +f2ijkp-f2ijk             &
                                      +f3ijkp-f3ijk )/h(i,j,k));
              f1ns(i,j,k)=f1ijk;
              f2ns(i,j,k)=f2ijk;
              f3ns(i,j,k)=f3ijk;
)
        ENDIF !skyedge
        IF (rightedge == 1) THEN
          i=np+1
X2DWallFullYZDomainLoopDC(
       f1ns(i,j,k)= donor(x_in(i-1,j,k), x_in(i,j,k  ), u1(i,j,k  ));
)
        ENDIF
        IF (topedge == 1) THEN
          j=mp+1
Y2DWallFullXZDomainLoopDC(
       f2ns(i,j,k)= donor(x_in(i,j-1,k), x_in(i,j,k  ), u2(i,j,k  ));
)
        ENDIF !topedge

        IF (skyedge == 1) THEN
          k=lp+1
Z2DWallFullXYDomainLoopDC(
       f3ns(i,j,k)=-donor(x_in(i,j,k-2), x_in(i,j,k-1), u3(i,j,k-1));
)
        ENDIF
        CALL ttend(1009,.TRUE.)
       IF(liner.eq.0) THEN
        !Full update needed as x is referenced in all edges
        CALL updated(xant,np,mp,lp,np,mp,lp,iup,ih)
        CALL ttbeg(1009)
!       CALL remove_cyclic_offset_full(xant,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
        !At the domain edges, copy the second last  meaningful domain point to x halo
        !This is needed to zero the corrective fluxes vcorr
        CALL cp_scnd_last_to_halo_xy_fulld(xant,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)

 
        !For the first antidiffusive iteration, f1,f2,f3 in traditional MPDATA
        !can be effectively replaced with C-grid advective velocities u1,u2,u3
        !This would not be valid for further antidiffusive iterations.
        !  rat2(z1,z2)=(z2-z1)*.5_euwp
        !  vdyf(x1,x2,a,rinv)=(abs(a)-a**2*rinv)*rat2(x1,x2)
        !  rat4(z0,z1,z2,z3)=(z3+z2-z1-z0)*.25_euwp
        !  vcorr(a,b,y0,y1,y2,y3,rinv)=-0.125_euwp*a*b*rinv*rat4(y0,y1,y2,y3)

        IF(gndedge == 1 ) THEN
            k=1
         FullXYDomainLoopDC(

                    hmx= 1._euwp/(0.5_euwp*(h(i-1,j,k)+h(i,j,k)));

                     v1(i,j,k)=vdyf(xant(i-1,j,k),xant(i,j,k),u1(i,j,k),hmx)    &
                          -0.125_euwp*u1(i  ,j,k)*hmx*                          &
                                    ((u2(i-1,j,k)+u2(i-1,j+1,k)                 &
                                     +u2(i,j+1,k)+u2(i  ,j  ,k) )               &
                              *rat4(xant(i-1,j-1,k),xant(i,j-1,k),              &
                                    xant(i-1,j+1,k),xant(i,j+1,k)));

                    hmy= 1._euwp/(0.5_euwp*(h(i,j-1,k)+h(i,j,k)));

                     v2(i,j,k)=vdyf(xant(i,j-1,k),xant(i,j,k),u2(i,j,k),hmy)    &
                          -0.125_euwp*u2(i,j,k)*hmy*                            &
                                    ((u1(i  ,j-1,k)+u1(i  ,j  ,k)               &
                                     +u1(i+1,j  ,k)+u1(i+1,j-1,k))              &
                              *rat4(xant(i-1,j-1,k),xant(i-1,j,k),              &
                                    xant(i+1,j-1,k),xant(i+1,j,k)));
                     v3(i,j,k)=0._euwp ;
)
        ENDIF
      InnerZFullXYDomainLoopDC(
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
        IF(skyedge == 1 ) THEN
            k=lp
         FullXYDomainLoopDC(
                    hmx= 1._euwp/(0.5_euwp*(h(i-1,j,k)+h(i,j,k)));

                     v1(i,j,k)=vdyf(xant(i-1,j,k),xant(i,j,k),u1(i,j,k),hmx)    &
                          -0.125_euwp*u1(i  ,j,k)*hmx*                          &
                                    ((u2(i-1,j,k)+u2(i-1,j+1,k)                 &
                                     +u2(i,j+1,k)+u2(i  ,j  ,k) )               &
                              *rat4(xant(i-1,j-1,k),xant(i,j-1,k),              &
                                    xant(i-1,j+1,k),xant(i,j+1,k)));

                    hmy= 1._euwp/(0.5_euwp*(h(i,j-1,k)+h(i,j,k)));

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

      IF(ipoles0 == 0) THEN 
        IF(ibcx0.eq.0) CALL flux_x_bcd(v1  ,np,mp,lp,ih, 0)
        IF(ibcy0.eq.0) CALL flux_y_bcd(v2  ,np,mp,lp,ih, 0) 
      ENDIF
        IF(ibcz0.eq.0) CALL flux_z_bc_nohalod(v3,np,mp,lp,ih, 2)

        CALL ttend(1009,.TRUE.)
        
        CALL updatelrd(v1,np+rightedge*(1-ipoles0),mp,lp,np+1,mp  ,lp,iupx*(1-ipoles0)+ipoles0,ih)
        CALL updatebtd(v2,np,mp+topedge  ,lp,np  ,mp+1,lp,iupy,ih)
        CALL updategsd(v3,np,mp,lp+skyedge,np,mp,lp+1,iupz,ih)
        CALL ttbeg(1009)
!     IF(ipoles0 == 0) THEN 
!       if(ibcx0.eq.1.and. leftedge.eq.1) v1(1   ,1:mp,1:lp)=v1(-1  ,1:mp,1:lp)
!       if(ibcx0.eq.1.and.rightedge.eq.1) v1(np+1,1:mp,1:lp)=v1(np+3,1:mp,1:lp)
!       if(ibcy0.eq.1.and. botedge.eq.1)  v2(1:np,1   ,1:lp)=v2(1:np,-1  ,1:lp)
!       if(ibcy0.eq.1.and. topedge.eq.1)  v2(1:np,mp+1,1:lp)=v2(1:np,mp+3,1:lp)
!     ELSE
!      if(botedge.eq.1) v2(1:np,1   ,1:lp)=0.
!      if(topedge.eq.1) v2(1:np,mp+1,1:lp)=0.
!     ENDIF
!      if(ibcz0.eq.1.and.gndedge.eq.1)   v3(1:np,1:mp,   1)=  v3(1:np,1:mp,  -1)
!      if(ibcz0.eq.1.and.skyedge.eq.1)   v3(1:np,1:mp,lp+1)=  v3(1:np,1:mp,lp+3)

        !                 non-osscilatory option
        CALL cp_last_to_halo_xy(x_in,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)

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
                        (pn(v1(i+1,j,k))+pp(v1(i,j,k))                      &
                        +pn(v2(i,j+1,k))+pp(v2(i,j,k))                      &
                        +pn(v3(i,j,k+1))+pp(-v3(i,j,k+1))+ep);

              cn(i,j,k)=(xant(i,j,k)-mnijk)*h(i,j,k)/                       &
                        (pp(v1(i+1,j,k))+pn(v1(i,j,k))                      &
                        +pp(v2(i,j+1,k))+pn(v2(i,j,k))                      &
                        +pp(v3(i,j,k+1))+pn(-v3(i,j,k+1))+ep);
)
        ENDIF !gndedge
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
                        (pn(v1(i+1,j,k))+pp(v1(i,j,k))                      &
                        +pn(v2(i,j+1,k))+pp(v2(i,j,k))                      &
                        +pn(v3(i,j,k+1))+pp(v3(i,j,k))+ep);

              cn(i,j,k)=(xant(i,j,k)-mnijk)*h(i,j,k)/                       &
                        (pp(v1(i+1,j,k))+pn(v1(i,j,k))                      &
                        +pp(v2(i,j+1,k))+pn(v2(i,j,k))                      &
                        +pp(v3(i,j,k+1))+pn(v3(i,j,k))+ep);
)
        IF(skyedge == 1) THEN
          k=lp
          FullXYDomainLoopDC(
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

              cp(i,j,k)=(mxijk-xant(i,j,k))*h(i,j,k)/                       &
                        (pn(v1(i+1,j,k))+pp(v1(i,j,k))                      &
                        +pn(v2(i,j+1,k))+pp(v2(i,j,k))                      &
                        +pn( -v3(i,j,k))+pp(v3(i,j,k))+ep);

              cn(i,j,k)=(xant(i,j,k)-mnijk)*h(i,j,k)/                       &
                        (pp(v1(i+1,j,k))+pn(v1(i,j,k))                      &
                        +pp(v2(i,j+1,k))+pn(v2(i,j,k))                      &
                        +pp( -v3(i,j,k))+pn(v3(i,j,k))+ep);
)
        ENDIF
        CALL ttend(1009,.TRUE.)
!       IF(lupdatemulti) THEN
!         CALL update3_multi(np,mp,lp,np,mp,lp,iup,ih,cp,cn)
!       ELSE
          CALL update3d(cp,np,mp,lp,np,mp,lp,iup,ih)
          CALL update3d(cn,np,mp,lp,np,mp,lp,iup,ih)
!       ENDIF

        CALL ttbeg(1009)
!       CALL remove_cyclic_offset_full(cp,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
!       CALL remove_cyclic_offset_full(cn,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)

        CALL zero_halo_xyd(cp,ibcx0,ibcy0,np,mp,lp,ih) !Only works for open bc`s
        CALL zero_halo_xyd(cn,ibcx0,ibcy0,np,mp,lp,ih)
        IF(gndedge == 1) THEN
          k=1
          FullXYDomainLoopDC(
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
                   pp(-v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k+1))   &
                  -pn(-v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k+1),cn(i,j,k  ));
              f3ijkp=                                                       &
                   pp(v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k+1),cn(i,j,k  ))    &
                  -pn(v3(i,j,k+1))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k+1));
              x_in(i,j,k)=(xant(i,j,k)-             &
                                   ( f1ijkp-f1ijk               &
                                    +f2ijkp-f2ijk               &
                                    +f3ijkp-f3ijk )/h(i,j,k));
              x_incomp(i,j,k)=x_in(i,j,k)*h(i,j,k);
              pexr1(i,j,k)=x_incomp(i,j,k)*pexe(i,j,k);
              pexr2(i,j,k)=x_incomp(i,j,k);
                   rhr(i,j,k)=rhoadv(i,j,k)/x_incomp(i,j,k);

              u1(i,j,k)=f1ns(i,j,k)+f1ijk;
              u2(i,j,k)=f2ns(i,j,k)+f2ijk;
              u3(i,j,k)=f3ns(i,j,k)+f3ijk;
)
        ENDIF

      InnerZFullXYDomainLoopDC(
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
              x_in(i,j,k)=(xant(i,j,k)-             &
                                   ( f1ijkp-f1ijk               &
                                    +f2ijkp-f2ijk               &
                                    +f3ijkp-f3ijk )/h(i,j,k));
              x_incomp(i,j,k)=x_in(i,j,k)*h(i,j,k);
              pexr1(i,j,k)=x_incomp(i,j,k)*pexe(i,j,k);
              pexr2(i,j,k)=x_incomp(i,j,k);
                   rhr(i,j,k)=rhoadv(i,j,k)/x_incomp(i,j,k);

              u1(i,j,k)=f1ns(i,j,k)+f1ijk;
              u2(i,j,k)=f2ns(i,j,k)+f2ijk;
              u3(i,j,k)=f3ns(i,j,k)+f3ijk;
)
        IF(skyedge == 1) THEN
          k=lp
         FullXYDomainLoopDC(
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
                   pp(-v3(i,j,k ))*min(1.0_euwp,cp(i,j,k-1),cn(i,j,k  ))    &
                  -pn(-v3(i,j,k ))*min(1.0_euwp,cp(i,j,k  ),cn(i,j,k-1));
              x_in(i,j,k)=(xant(i,j,k)-             &
                                   ( f1ijkp-f1ijk               &
                                    +f2ijkp-f2ijk               &
                                    +f3ijkp-f3ijk )/h(i,j,k));
              x_incomp(i,j,k)=x_in(i,j,k)*h(i,j,k);
              pexr1(i,j,k)=x_incomp(i,j,k)*pexe(i,j,k);
              pexr2(i,j,k)=x_incomp(i,j,k);
                   rhr(i,j,k)=rhoadv(i,j,k)/x_incomp(i,j,k);
              u1(i,j,k)=f1ns(i,j,k)+f1ijk;
              u2(i,j,k)=f2ns(i,j,k)+f2ijk;
              u3(i,j,k)=f3ns(i,j,k)+f3ijk;
)
        ENDIF
        CALL ttend(1009,.TRUE.)
          CALL update(u1,np+rightedge,mp,lp,np+1,mp,lp,iup,ih)
          CALL update(u2,np,mp+topedge,lp,np,mp+1,lp,iup,ih)
          CALL update(u3,np,mp,lp+skyedge,np,mp,lp+1,iup,ih)

        CALL ttbeg(1009)

            IF (rightedge == 1) THEN
               i=np+1
X2DWallFullYZDomainLoopDC(
                    u1(i,j,k)=f1ns(i,j,k)+                                 &
                     pp(v1(i,j,k  ))*min(1.0_euwp,cp(i  ,j,k),cn(i-1,j,k)) &
                    -pn(v1(i,j,k  ))*min(1.0_euwp,cp(i-1,j,k),cn(i  ,j,k));
)
            ENDIF
            IF (topedge == 1) THEN
                j=mp+1
Y2DWallFullXZDomainLoopDC(
                   u2(i,j,k)=f2ns(i,j,k) + &
                      pp(v2(i,j,k  ))*min(1.0_euwp,cp(i,j  ,k),cn(i,j-1,k)) &
                     -pn(v2(i,j,k  ))*min(1.0_euwp,cp(i,j-1,k),cn(i,j  ,k));
)
            ENDIF !topedge
        IF (skyedge == 1) THEN
            k=lp+1
Z2DWallFullXYDomainLoopDC(
                  u3(i,j,k)=f3ns(i,j,k)+ &
                   pp(-v3(i,j,k-1  ))*min(1.0_euwp,cp(i,j,k-2),cn(i,j,k-1)) &
                  -pn(-v3(i,j,k-1  ))*min(1.0_euwp,cp(i,j,k-1),cn(i,j,k-2));
)
        ENDIF
       ELSE !liner
CGRIDXFullYZDomainLoopDC(u1(i,j,k)=f1ns(i,j,k);)
CGRIDYFullXZDomainLoopDC(u2(i,j,k)=f2ns(i,j,k);)
CGRIDZFullXYDomainLoopDC(u3(i,j,k)=f3ns(i,j,k);)
FullXYZDomainLoopDC(
            x_in(i,j,k)=    xant(i,j,k);
        x_incomp(i,j,k)=    x_in(i,j,k)*       h(i,j,k);
           pexr1(i,j,k)=x_incomp(i,j,k)*    pexe(i,j,k);
           pexr2(i,j,k)=x_incomp(i,j,k);
             rhr(i,j,k)=  rhoadv(i,j,k)/x_incomp(i,j,k);
)

       ENDIF !liner
          CALL update(u1,np+rightedge*(1-ipoles0),mp,lp,np+1,mp,lp,1,ih)
          CALL update(u2,np              ,mp+topedge,lp,np,mp+1,lp,1,ih)
          CALL update(u3,np              ,mp,lp+skyedge,np,mp,lp+1,1,ih)

        CALL ttend(1009,.TRUE.)

        IF(lupdatemulti) THEN
          CALL update_multi(np,mp,lp,np,mp,lp,iup,ih,pexr1,pexr2,x_in)
        ELSE
          CALL update(pexr1,np,mp,lp,np,mp,lp,iup,ih)
          CALL update(pexr2,np,mp,lp,np,mp,lp,iup,ih)
          CALL update(x_in,np,mp,lp,np,mp,lp,iup,ih)
        ENDIF

        CALL ttend(9)
    END SUBROUTINE mpdata3d_rho_gauge_gpubc
END MODULE module_mpdata3d_rho_gauge_gpubc

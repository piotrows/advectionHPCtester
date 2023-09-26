#include "renames.inc"
MODULE module_upwind3d_gpubc
   USE precisions, ONLY  : iintegers,euwp
   USE mpi_parallel, ONLY: leftedge,rightedge,botedge,topedge,gndedge,skyedge
   USE mpdataoperators
#ifdef CUDACODE
   USE mpi_parallel, ONLY: istream1, istream2
#endif
   USE epsilons, ONLY: ep => ep_nonos
    USE noise, ONLY: perturb_2D_signal
#ifdef PNETCDF
   USE mpi_parallel, ONLY: pnet_out_chunk
#endif
#ifdef CUDACODE
   USE cudafor
#endif
   IMPLICIT NONE
#include "defines.inc"
!#include "mpdataoperators.inc"

CONTAINS
SUBROUTINE upwind3d_gpubc(u1,u2,u3, &
                          x_in,xant,&
                          xforc_impl,xforc_expl,&
                          rhr,h,&
                          bcx,bcy, &
                          ipoles0,& 
                          ibcx0,ibcy0,ibcz0,&
                          np,mp,lp,ih,dt, &
                          do_serialization_in,do_serialization_out )
!---------------------------------------------------------------------!
      USE mpi_parallel, ONLY:  update,update3,update3d
      USE mpi_parallel, ONLY: ttbeg,ttend
      USE mpi_parallel, ONLY: iup
        LOGICAL,OPTIONAL :: do_serialization_in
        LOGICAL,OPTIONAL :: do_serialization_out
        INTEGER(KIND=iintegers), INTENT(IN) :: ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih
        REAL(KIND=euwp), INTENT(IN) :: dt
#ifdef CUDACODE
        REAL(KIND=euwp),MANAGED, INTENT(INOUT) ::  &
                     x_in(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih)
        REAL(KIND=euwp),MANAGED, INTENT(IN)   ::   &
                       u1(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+ih),    &
                       u2(1-ih:np+ih  , 1-ih:mp+1+ih, 1-ih:lp+ih),    &
                       u3(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+1+ih)
        REAL(KIND=euwp), DEVICE, INTENT(OUT) ::  &
                     xant(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih)

        REAL(KIND=euwp), MANAGED, INTENT(IN) ::  &
                        h(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih), &
                      rhr(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih), &
               xforc_expl(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih)

        REAL(KIND=euwp), MANAGED, INTENT(IN), OPTIONAL ::  &
               xforc_impl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

        REAL(KIND=euwp), MANAGED, INTENT(IN) :: &
                             bcx(mp,lp, 2), &
                             bcy(np,lp, 2)
#else
        REAL(KIND=euwp), INTENT(INOUT) ::  &
                     x_in(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih)
        REAL(KIND=euwp), INTENT(IN)   ::   &
                       u1(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+ih),    &
                       u2(1-ih:np+ih  , 1-ih:mp+1+ih, 1-ih:lp+ih),    &
                       u3(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+1+ih)
        REAL(KIND=euwp),  INTENT(OUT) ::  &
                     xant(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih)

        REAL(KIND=euwp),  INTENT(IN) ::  &
                        h(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih), &
                      rhr(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih), &
               xforc_expl(1-ih:np+ih  , 1-ih:mp+ih  , 1-ih:lp+ih)

        REAL(KIND=euwp), INTENT(IN), OPTIONAL ::  &
               xforc_impl(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

        REAL(KIND=euwp), INTENT(IN) :: &
                             bcx(mp,lp, 2), &
                             bcy(np,lp, 2)
#endif


      INTEGER(KIND=iintegers) ::  i,j,k
      INTEGER istat
      REAL(KIND=euwp) :: hinv,f1ijk,f1ijkp,f2ijk,f2ijkp,f3ijk,f3ijkp
      REAL(KIND=euwp) :: x_ijk, x_im, x_ip, x_jm, x_jp, x_km, x_kp, u1_ijk, u1_ip, u2_ijk, u2_jp, u3_ijk, u3_kp

      CALL ttbeg(16)

!       IF(PRESENT(do_serialization_in).OR.PRESENT(do_serialization_out)) THEN
!Do nothing but remove compiler warning about unused dummy variables
!       ENDIF

!       IF(PRESENT(xforc_impl)) THEN
!          x_in(1:np,1:mp,1:lp) =           x_in(1:np,1:mp,1:lp)              &
!                                   + xforc_impl(1:np,1:mp,1:lp)*.5_euwp*dt   &
!                                   + xforc_expl(1:np,1:mp,1:lp)*dt
!       ELSE
!          x_in(1:np,1:mp,1:lp) =           x_in(1:np,1:mp,1:lp)              &
!                                   + xforc_expl(1:np,1:mp,1:lp)*dt
!       ENDIF
        CALL putbcinxh(x_in,bcx,bcy,ibcx0,ibcy0,np,mp,lp,ih)

        CALL update3(x_in,np,mp,lp,np,mp,lp,iup,ih)

!       !Put boundary conditions in x and y direction into the x_in halo
!       CALL remove_cyclic_offset_full(x_in,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)

      IF(gndedge == 1) THEN
         k=1
         FullXYDomainLoopDCs(stream=istream2,
               hinv=1._euwp/h(i,j,k);
              x_ijk=x_in(i  ,j  ,k  );
!              x_im=x_in(i-1,j  ,k  );
!              x_ip=x_in(i+1,j  ,k  );
!              x_jm=x_in(i  ,j-1,k  );
!              x_jp=x_in(i  ,j+1,k  );
!              x_km=x_in(i  ,j  ,k-1);
!              x_kp=x_in(i  ,j  ,k+1);
!            u1_ijk=  u1(i  ,j  ,k  );
!             u1_ip=  u1(i+1,j  ,k  );
!            u2_ijk=  u2(i  ,j  ,k  );
!             u2_jp=  u2(i  ,j+1,k  );
!            u3_ijk=  u3(i  ,j  ,k  );
!             u3_kp=  u3(i  ,j  ,k+1);
!              f1ijkp=max(0._euwp,u1_ip )*x_ijk - (-min(0._euwp,u1_ip ))*x_ip   !donor(x_in(i  ,j,k), x_in(i+1,j,k), u1(i+1,j,k));
!              f1ijk =max(0._euwp,u1_ijk)*x_im  - (-min(0._euwp,u1_ijk))*x_ijk  !donor(x_in(i-1,j,k), x_in(i  ,j,k), u1(i  ,j,k));
!              f2ijkp=max(0._euwp,u2_jp )*x_ijk - (-min(0._euwp,u2_jp ))*x_jp   !donor(x_in(i,j  ,k), x_in(i,j+1,k), u2(i,j+1,k));
!              f2ijk =max(0._euwp,u2_ijk)*x_jm  - (-min(0._euwp,u2_ijk))*x_ijk  !donor(x_in(i,j-1,k), x_in(i,j  ,k), u2(i,j  ,k));
!              f3ijkp=max(0._euwp,u3_kp )*x_ijk - (-min(0._euwp,u3_kp ))*x_kp   !donor(x_in(i,j,k  ), x_in(i,j,k+1), u3(i,j,k+1));
               f1ijkp=donor(x_in(i  ,j,k), x_in(i+1,j,k), u1(i+1,j,k));
               f1ijk =donor(x_in(i-1,j,k), x_in(i  ,j,k), u1(i  ,j,k));
               f2ijkp=donor(x_in(i,j  ,k), x_in(i,j+1,k), u2(i,j+1,k));
               f2ijk =donor(x_in(i,j-1,k), x_in(i,j  ,k), u2(i,j  ,k));
               f3ijkp=donor(x_in(i,j,k  ), x_in(i,j,k+1), u3(i,j,k+1));
               f3ijk=-f3ijkp;
               xant(i,j,k)=rhr(i,j,k)*(x_ijk -                       &
                                           ( f1ijkp-f1ijk            &
                                            +f2ijkp-f2ijk            &
                                            +f3ijkp-f3ijk )*hinv);
         )
      ENDIF

      IF(skyedge == 1) THEN
         k=lp
         FullXYDomainLoopDCs(stream=istream2,
               hinv=1._euwp/h(i,j,k);
              x_ijk=x_in(i  ,j  ,k  );
               x_im=x_in(i-1,j  ,k  );
               x_ip=x_in(i+1,j  ,k  );
               x_jm=x_in(i  ,j-1,k  );
               x_jp=x_in(i  ,j+1,k  );
               x_km=x_in(i  ,j  ,k-1);
               x_kp=x_in(i  ,j  ,k+1);
             u1_ijk=  u1(i  ,j  ,k  );
              u1_ip=  u1(i+1,j  ,k  );
             u2_ijk=  u2(i  ,j  ,k  );
              u2_jp=  u2(i  ,j+1,k  );
             u3_ijk=  u3(i  ,j  ,k  );
              u3_kp=  u3(i  ,j  ,k+1);
               f1ijkp=max(0._euwp,u1_ip )*x_ijk - (-min(0._euwp,u1_ip ))*x_ip   !donor(x_in(i  ,j,k), x_in(i+1,j,k), u1(i+1,j,k));
               f1ijk =max(0._euwp,u1_ijk)*x_im  - (-min(0._euwp,u1_ijk))*x_ijk  !donor(x_in(i-1,j,k), x_in(i  ,j,k), u1(i  ,j,k));
               f2ijkp=max(0._euwp,u2_jp )*x_ijk - (-min(0._euwp,u2_jp ))*x_jp   !donor(x_in(i,j  ,k), x_in(i,j+1,k), u2(i,j+1,k));
               f2ijk =max(0._euwp,u2_ijk)*x_jm  - (-min(0._euwp,u2_ijk))*x_ijk  !donor(x_in(i,j-1,k), x_in(i,j  ,k), u2(i,j  ,k));
               f3ijk =max(0._euwp,u3_ijk)*x_km  - (-min(0._euwp,u3_ijk))*x_ijk  !donor(x_in(i,j,k-1), x_in(i,j,k  ), u3(i,j,k  )) ; 
               f3ijkp=-f3ijk;
               xant(i,j,k)=rhr(i,j,k)*(x_ijk-                        &
                                           ( f1ijkp-f1ijk            &
                                            +f2ijkp-f2ijk            &
                                            +f3ijkp-f3ijk )*hinv);
         )
      ENDIF

      InnerZFullXYDomainLoopDCs(stream=istream1,
               hinv=1._euwp/h(i,j,k) ;
              x_ijk=x_in(i  ,j  ,k  );
               x_im=x_in(i-1,j  ,k  );
               x_ip=x_in(i+1,j  ,k  );
               x_jm=x_in(i  ,j-1,k  );
               x_jp=x_in(i  ,j+1,k  );
               x_km=x_in(i  ,j  ,k-1);
               x_kp=x_in(i  ,j  ,k+1);
             u1_ijk=  u1(i  ,j  ,k  );
              u1_ip=  u1(i+1,j  ,k  );
             u2_ijk=  u2(i  ,j  ,k  );
              u2_jp=  u2(i  ,j+1,k  );
             u3_ijk=  u3(i  ,j  ,k  );
              u3_kp=  u3(i  ,j  ,k+1);
               f1ijkp=max(0._euwp,u1_ip )*x_ijk - (-min(0._euwp,u1_ip ))*x_ip   !donor(x_in(i  ,j,k), x_in(i+1,j,k), u1(i+1,j,k));
               f1ijk =max(0._euwp,u1_ijk)*x_im  - (-min(0._euwp,u1_ijk))*x_ijk  !donor(x_in(i-1,j,k), x_in(i  ,j,k), u1(i  ,j,k));
               f2ijkp=max(0._euwp,u2_jp )*x_ijk - (-min(0._euwp,u2_jp ))*x_jp   !donor(x_in(i,j  ,k), x_in(i,j+1,k), u2(i,j+1,k));
               f2ijk =max(0._euwp,u2_ijk)*x_jm  - (-min(0._euwp,u2_ijk))*x_ijk  !donor(x_in(i,j-1,k), x_in(i,j  ,k), u2(i,j  ,k));
               f3ijkp=max(0._euwp,u3_kp )*x_ijk - (-min(0._euwp,u3_kp ))*x_kp   !donor(x_in(i,j,k  ), x_in(i,j,k+1), u3(i,j,k+1));
               f3ijk =max(0._euwp,u3_ijk)*x_km  - (-min(0._euwp,u3_ijk))*x_ijk  !donor(x_in(i,j,k-1), x_in(i,j,k  ), u3(i,j,k  )) ; 
                xant(i,j,k)=rhr(i,j,k)*(x_ijk-                        &
                                            ( f1ijkp-f1ijk            &
                                             +f2ijkp-f2ijk            &
                                             +f3ijkp-f3ijk )*hinv);
                             )

      CALL ttend(16)
   END SUBROUTINE upwind3d_gpubc
END MODULE module_upwind3d_gpubc

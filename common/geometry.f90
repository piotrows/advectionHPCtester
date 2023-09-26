#include "../advection/src_algorithms/renames.inc"
MODULE geometry
   USE precisions
   USE mod_parameters, ONLY: n,m,l,np,mp,lp,ih
   USE mod_parameters, ONLY: rds,rdsi,fcr0,btpl,ideep
   USE mod_parameters, ONLY: ibcx,ibcy,ibcz
   USE mpi_parallel, ONLY:  update,update3,updateflt,updategs,mybarrier
   USE mpi_parallel, ONLY: iup,npos,mpos,lpos,nsubpos,msubpos,lsubpos
   USE mpi_parallel, ONLY: mype,leftedge,rightedge,botedge,topedge,gndedge,skyedge
   USE bconditions, ONLY: halo_one_sided_diff_x,   &
                          halo_one_sided_diff_y,   &
                          halo_one_sided_diff_z
   USE bconditions, ONLY: halo_one_sided_diff_flt_x,   &
                          halo_one_sided_diff_flt_y
   USE bconditions, ONLY: halo_copy_edge_to_halo
USE geometry_datafields, ONLY: gf11, gf12, gf13
USE geometry_datafields, ONLY: gf21, gf22, gf23
USE geometry_datafields, ONLY: gf31, gf32, gf33
USE geometry_datafields, ONLY: x,y,z, &
                      xcr,ycr,sina,cosa,tnga,sinx,cosx,fcr2,fcr3, &
                      zcr,g11,g12 ,g13 ,g21 ,g22 ,g23 ,g33, &
                      gac,gmm,gmmi,g33i,g11g22Mg12g21i,gaci
USE geometry_datafields, ONLY: sg11,sg12 ,sg13 ,sg21 ,sg22 ,sg23 ,sg33      
USE geometry_datafields, ONLY: tnxx, tnxy, tnxz, sfflx, &
    tnyx, tnyy, tnyz, sffly, tnzx, tnzy, tnzz, sfflz, tnxxV, tnxyV, tnxzV, tnyxV, &
    tnyyV, tnyzV, tnzxV, tnzyV, tnzzV      
                      
#ifdef PNETCDF
   USE mpi_parallel, ONLY: pnet_out_chunk
#endif /*PNETCDF*/

   IMPLICIT NONE
   REAL (KIND=euwp) :: dxa,dya
   REAL (KIND=euwp) :: cx,cy,cz,gnr,pnx,pny,pnz,zsx,zsy
   INTEGER(KIND=iintegers) :: ii, jj, kk
CONTAINS
#include "../advection/src_algorithms/defines.inc"
!------------------------------------------------------------------
   SUBROUTINE topolog
!------------------------------------------------------------------
      USE mod_parameters, ONLY:  isphere, ipoles
      USE mod_parameters, ONLY: dx,dy,dz,dxa,dya
      USE mod_parameters, ONLY: pi,pi2,pih,btpl,yj0,ang
      INTEGER (KIND=IINTEGERS) :: i,j,k,ia,ja,ka

! ----  specify computational grid
      DO  i=1,n
         IF(isphere==1) THEN
!           x(i)=-pi+(i-1)*dxa                               ! sphere (radians)
            x(i)=(i-1)*dxa                               ! sphere (radians)
         ELSE
            x(i)=(i-1)*dx -.5*(n-1)*dx*(1-ibcx)           ! Cartesian (m)
         END IF
      ENDDO
!Symmetry correction?:
              if(isphere.eq.1.and.1.eq.0) then
       do  i=1,floor(n/2.)
        x(floor(n/2.)+1-i)=-i*dxa                           ! sphere (radians)
       enddo
        x(floor(n/2.)+1)=0._euwp
       do  i=2,floor(n/2.)
        x(floor(n/2.)+i)=(i-1)*dxa                           ! sphere (radians)
       enddo
       do  i=floor(n/2.)+2,n
            x(i)=x(i-1)+(x(i-floor(n/2.))-x(i-floor(n/2.)-1))
       enddo
       
              endif



      do  j=1,m
              if(isphere == 1) then
        y(j)=-pih+(j-0.5)*dya                        ! sphere (radians)
!special
!        y(j)=+(43.5/180.0)*pi+(j-0.5)*dya

       y(j)=-(90.0/180.0)*pi+(j-0.5)*dya
         ELSE
            y(j)=(j-1)*dy-.5*(m-1)*dy *(1-ibcy)           ! Cartesian (m)
         END IF
      ENDDO
!correction for exact symmetry for isphere=1
              if(isphere.eq.1.and.1.eq.0) then
       do  j=1,m
         y(m-j+1)=0.5*(y(m-j+1)+pih-(j-0.5)*dya)               ! sphere (radians)
       enddo
              end if

      do  k=1,l
         z(k)=(k-1)*dz
      ENDDO

!------------------------------------------------------------------
!--- compute horizontal grid point "physical" locations
!--- and compute trigonometric and Coriolis functions
!------------------------------------------------------------------
      if(isphere.eq.0) then  
        do j=1,mp
          ja=msubpos + j
          do i=1,np
            ia=nsubpos + i
            xcr(i,j)=x(ia)
            ycr(i,j)=y(ja)
            cosa(i,j)=1._euwp
            sina(i,j)=0.
            tnga(i,j)=0.
            fcr2(i,j)=fcr0*(cos(yj0)-btpl*sin(yj0)*y(ja))
            fcr3(i,j)=fcr0*(sin(yj0)+btpl*cos(yj0)*y(ja))
          enddo
        enddo
      else 
      do j=1,mp
        ja=msubpos + j
        do i=1,np
         ia=nsubpos + i
          xcr(i,j)=x(ia)*rds
          ycr(i,j)=y(ja)*rds
         cosa(i,j)=0.5_euwp*( cos(y(m-ja+1))+cos(y(ja)))
         cosa(i,j)=cos(y(ja))
         sina(i,j)=sin(y(ja)) !0.5_euwp*(-sin(y(m-ja+1))+sin(y(ja)))
         tnga(i,j)=tan(y(ja))
         cosx(i,j)=0.5_euwp*( cos(x(n-ia+1))+cos(x(ia)))
         sinx(i,j)=sin(x(ia))
        enddo
      enddo
      do j=1,mp
        ja=msubpos + j
        do i=1,np
         ia=nsubpos + i
        if(ia.lt.(floor(n/2.)+1)) then
          sinx(i,j)=-sin(x(floor(n/2.)+ia))
        endif
        enddo
      enddo
     
      do j=1,mp
        do i=1,np
          fcr2(i,j)=fcr0*cosa(i,j)*ideep
          fcr3(i,j)=fcr0*sina(i,j)
        enddo
      enddo
      endif
      call updateflt( xcr,iup,0,np,mp,ih)
      call updateflt( ycr,iup,0,np,mp,ih)

!------------------------------------------------------------------
!--- compute vertical grid point "physical" locations
!------------------------------------------------------------------

        do k=1,lp
          ka=lsubpos + k
          zcr(1:np,1:mp,k)=z(ka)  
        enddo
      call update3(zcr,np,mp,lp,np,mp,lp,iup,ih)

      call updateflt(sinx,iup,0,np,mp,ih)

      call updateflt(cosa,iup,0,np,mp,ih)
      call updateflt(sina,iup,0,np,mp,ih)
      call updateflt(tnga,iup,0,np,mp,ih)


!-----------------------------------------------------------c
!create values for metric terms at the boundary for POLES == 1
!-----------------------------------------------------------c
      IF(isphere==1.AND.ipoles==1) THEN
         IF (leftedge==1) THEN
            DO j=1,mp
              xcr(0,j)=xcr(0,j)-pi2*rds
            ENDDO
         ENDIF
         IF (rightedge==1) THEN
            DO j=1,mp
               xcr(np+1,j)=xcr(np+1,j)+pi2*rds
            ENDDO
         ENDIF
         IF (botedge==1) THEN
            DO i=1,np
               xcr(i,0)= xcr(i,1)
               ycr(i,0)= ycr(i,1)
            ENDDO
         ENDIF
         IF (topedge==1) THEN
            DO i=1,np
               xcr(i,mp+1)= xcr(i,mp)
               ycr(i,mp+1)= ycr(i,mp)
            ENDDO
         ENDIF
      ENDIF

   END SUBROUTINE topolog

    function fi(rad,width) result(r)
      REAL(KIND=euwp), intent(in) :: rad,width ! input
      REAL(KIND=euwp)             :: r ! output
      REAL(KIND=euwp)             :: pi
      pi=acos(-1._dp)
      r = 1.*exp(-rad**2/width**2) ! .5*(1.+cos(pi*rad))
      r = .5*(1.+cos(pi*rad))
    end function fi

    function rd(xx, yy) result(r)
      REAL(KIND=euwp), intent(in) :: xx,yy ! input
      REAL(KIND=euwp)             :: r ! output
      r = sqrt((xx/x((n/2)))**2+(yy/y((m/2)))**2)
    end function rd
!------------------------------------------------------------------
   SUBROUTINE topol_vert(amp,width,radius,np,mp,ih)
!------------------------------------------------------------------
      USE mod_parameters, ONLY: dz,l
      INTEGER (KIND=IINTEGERS) :: i,j,k,ia,ja,ka,np,mp,ih
      REAL(KIND=euwp)pi,pi2,pih
      REAL(KIND=euwp) angle,zb,cang,sang,xx,yy,rr 
      REAL(KIND=euwp),INTENT(IN)::amp,radius,width
      REAL(KIND=euwp), DIMENSION(np,mp) :: zs,gi
      ! fi(rad)=-.5*(1.+cos(pi*rad))
      ! rd(xx,yy)=sqrt((xx/x((n/2)))**2+(yy/y((m/2)))**2)

      pi=acos(-1._dp)
      pi2=2._euwp*pi
      pih=pi/2._euwp
!------------------------------------------------------------------
!--- compute vertical grid point "physical" locations
!------------------------------------------------------------------
            zb=dz*(l-1)
!          amp=zb/4._euwp
        angle=0._euwp
      cang=cos(pi*angle/180._euwp)
      sang=sin(pi*angle/180._euwp)
          do  j=1,mp
            ja=msubpos + j
            do  i=1,np
              ia= nsubpos + i
         xx= cang*(x(ia)-x((n/2)))+sang*(y(ja)-y((m/2)))
         yy=-sang*(x(ia)-x((n/2)))+cang*(y(ja)-y((m/2)))
         rr=rd(xx,yy)
             zs(i,j)= 0.
!             zs(i,j)=amp*fi(rr)
             zs(i,j)=amp*merge(0._euwp,fi(rr,width),(radius/rds-rr)<0._euwp)
            enddo
           enddo
!          print *,'zb is',zb
!          print *,'max zs  is',maxval(zs(1:np,1:mp))
!          print *,'min dif is',maxval(zb-zs(1:np,1:mp))
          gi(1:np,1:mp)=zb/(zb-zs(1:np,1:mp))

        do k=1,lp
          ka=lsubpos + k
          zcr(1:np,1:mp,k)=z(ka)/gi(1:np,1:mp)+zs(1:np,1:mp)  
!         IF(npos.eq.1.and.mpos.eq.1) print *,'zcr(1,1,k),z(ka)',zcr(1,1,k),z(ka)
!         CALL flush(6)
        enddo
!        STOP 'topol_vert'
      call update3(zcr,np,mp,lp,np,mp,lp,iup,ih)

   END SUBROUTINE topol_vert
!------------------------------------------------------------------
   SUBROUTINE topolog_lib( xcr_in,ycr_in,zcr_in, &
                         fcr2_in,fcr3_in,       &
                         cosa_in,sina_in,tnga_in, &
                         np,mp,lp,n,m,l,ih)
!------------------------------------------------------------------
      USE mod_parameters, ONLY:  isphere, ipoles
      USE mod_parameters, ONLY: dx,dy,dz,dxa,dya
      USE mod_parameters, ONLY: pi,pi2,pih
      INTEGER (KIND=IINTEGERS) :: i,j,k
      INTEGER (KIND=IINTEGERS) :: np,mp,lp,n,m,l,ih 
      REAL(KIND=euwp), DIMENSION(1-ih:np+ih,1-ih:mp+ih)   ::  xcr_in,ycr_in
      REAL(KIND=euwp), DIMENSION(np,mp)   :: fcr2_in,fcr3_in
      REAL(KIND=euwp), DIMENSION(1-ih:np+ih,1-ih:mp+ih)   :: cosa_in,sina_in,tnga_in
      REAL(KIND=euwp), DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) ::  zcr_in


! ----  specify computational grid
      DO  i=1,n
         IF(isphere==1) THEN
            x(i)=-pi+(i-1)*dxa                               ! sphere (radians)
         ELSE
            x(i)=(i-1)*dx -.5*(n-1)*dx*(1-ibcx)           ! Cartesian (m)
         END IF
      ENDDO

!Symmetry correction?:
              if(isphere.eq.1) then
       do  i=1,floor(n/2.)
        x(floor(n/2.)+1-i)=-i*dxa                           ! sphere (radians)
       enddo
        x(floor(n/2.)+1)=0._euwp
       do  i=2,floor(n/2.)
        x(floor(n/2.)+i)=(i-1)*dxa                           ! sphere (radians)
       enddo
       do  i=floor(n/2.)+2,n
            x(i)=x(i-1)+(x(i-floor(n/2.))-x(i-floor(n/2.)-1))
       enddo
       
              endif


      do  j=1,m
              if(isphere == 1) then
        y(j)=-pih+(j-0.5)*dya                        ! sphere (radians)
!       y(j)=-(90.0/180.0)*pi+(j-0.5)*dya
         ELSE
            y(j)=(j-1)*dy-.5*(m-1)*dy *(1-ibcy)           ! Cartesian (m)
         END IF
      ENDDO
!correction for exact symmetry for isphere=1
              if(isphere.eq.1) then
       do  j=1,m
         y(m-j+1)=0.5*(y(m-j+1)+pih-(j-0.5)*dya)               ! sphere (radians)
       enddo
              end if

      do  k=1,l
         z(k)=(k-1)*dz
      ENDDO
!------------------------------------------------------------------
!--- compute trigonometric functions
!------------------------------------------------------------------
       xcr(1:np,1:mp     )= xcr_in(1:np,1:mp)
       ycr(1:np,1:mp     )= ycr_in(1:np,1:mp)
       zcr(1:np,1:mp,1:lp)= zcr_in(1:np,1:mp,1:lp)
      fcr2(1:np,1:mp     )=fcr2_in(1:np,1:mp)
      fcr3(1:np,1:mp     )=fcr3_in(1:np,1:mp)
      cosa(1:np,1:mp     )=cosa_in(1:np,1:mp)
      sina(1:np,1:mp     )=sina_in(1:np,1:mp)
      tnga(1:np,1:mp     )=tnga_in(1:np,1:mp)
      call updateflt( xcr,iup,0,np,mp,ih)
      call updateflt( ycr,iup,0,np,mp,ih)
      call update3(zcr,np,mp,lp,np,mp,lp,iup,ih)
      call updateflt(cosa,iup,0,np,mp,ih)
      call updateflt(sina,iup,0,np,mp,ih)

   END SUBROUTINE topolog_lib

!----------------------------------------------------------------------!
   SUBROUTINE metryc_lib(g11_in,g12_in,g13_in, &
                         g21_in,g22_in,g23_in, &
                         g33_in,gac_in,gmm_in,np,mp,lp,ih)  
!----------------------------------------------------------------------!

      USE mod_parameters,  ONLY:                     &
         dxih, dyih, dzih, rds, isphere, ideep, ipoles
      USE mpi_parallel, ONLY: update,  updateflt, update3
      !----------------------------------------------------------------------!

      INTEGER (KIND=IINTEGERS) :: np,mp,lp,ih 
      REAL_euwp, DIMENSION(np,mp,lp), INTENT(IN) :: & 
      g11_in,g12_in,g13_in,g21_in,g22_in,g23_in,g33_in

      REAL_euwp, DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih), INTENT(IN) :: & 
      gac_in,gmm_in

      INTEGER(KIND=iintegers) :: iflip, i, j, k

      REAL(KIND=euwp) :: g110i, g220i, g330i

      !----------------------------------------------------------------------!


      DO k=1,lp
         DO j=1,mp
            DO i=1,np
              gmm(i,j,k)=gmm_in(i,j,k)
              gmmi(i,j,k)=1._euwp/gmm(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !CALL update(gmm,np,mp,lp,np,mp,lp,iup,ih)

      DO k=1,lp
         DO j=1,mp
            DO i=1,np
               G110i=(gmm(i,j,k)*cosa(i,j))
               G220i=gmm(i,j,k)
               G330i=1._euwp
               g11(i,j,k)=g11_in(i,j,k)
               g12(i,j,k)=g12_in(i,j,k)
               g13(i,j,k)=g13_in(i,j,k)
               g21(i,j,k)=g21_in(i,j,k)
               g22(i,j,k)=g22_in(i,j,k)
               g23(i,j,k)=g23_in(i,j,k)
               g33(i,j,k)=g33_in(i,j,k)
              sg11(i,j,k)=g11(i,j,k)*g110i
              sg12(i,j,k)=g12(i,j,k)*g110i
              sg13(i,j,k)=g13(i,j,k)*g110i
              sg21(i,j,k)=g21(i,j,k)*g220i
              sg22(i,j,k)=g22(i,j,k)*g220i
              sg23(i,j,k)=g23(i,j,k)*g220i
              sg33(i,j,k)=g33(i,j,k)*g330i
               gac(i,j,k)=gac_in(i,j,k)
              gaci(i,j,k)=1._euwp/gac(i,j,k)
    g11g22Mg12g21i(i,j,k)= 1._euwp/ &
                         (g11(i,j,k)*g22(i,j,k)-g12(i,j,k)*g21(i,j,k))
        g33i(i,j,k)= 1._euwp/g33(i,j,k)
            ENDDO
         ENDDO
      ENDDO 
      g11(0   ,1:mp,1:lp) = g11(1 ,1:mp,1:lp)
      CALL halo_copy_edge_to_halo(g11,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g12,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g13,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g21,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g22,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g23,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g33,np,mp,lp,ih)
      

      CALL update(zcr,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g11,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g12,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g13,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g21,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g22,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g23,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g33,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(gac,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(gmm,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg11,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg12,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg13,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg21,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg22,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg23,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg33,np,mp,lp,np,mp,lp,iup,ih)     

      CALL halo_copy_edge_to_halo(g11,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g12,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g13,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g21,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g22,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g23,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g33,np,mp,lp,ih)

      
#if(1==1)
      DO k=1,lp  
         DO j=1,mp  
            DO i=1,np  
            gf11(i,j,k)=.5*( ( g11(i+1,j,k)*g11(i+1,j,k) + g21(i+1,j,k)*g21(i+1,j,k) )&
                            +( g11(i  ,j,k)*g11(i  ,j,k) + g21(i  ,j,k)*g21(i  ,j,k) ) )
            gf12(i,j,k)=.5*( ( g11(i+1,j,k)*g12(i+1,j,k) + g21(i+1,j,k)*g22(i+1,j,k) )&
                            +( g11(i  ,j,k)*g12(i  ,j,k) + g21(i  ,j,k)*g22(i  ,j,k) ) )
            gf13(i,j,k)=.5*( ( g11(i+1,j,k)*g13(i+1,j,k) + g21(i+1,j,k)*g23(i+1,j,k) )&
                            +( g11(i  ,j,k)*g13(i  ,j,k) + g21(i  ,j,k)*g23(i  ,j,k) ) )
            gf21(i,j,k)=.5*( ( g11(i,j+1,k)*g12(i,j+1,k) + g21(i,j+1,k)*g22(i,j+1,k) )&
                            +( g11(i,j  ,k)*g12(i,j  ,k) + g21(i,j  ,k)*g22(i,j  ,k) ) )
            gf22(i,j,k)=.5*( ( g12(i,j+1,k)*g12(i,j+1,k) + g22(i,j+1,k)*g22(i,j+1,k) )&
                            +( g12(i,j  ,k)*g12(i,j  ,k) + g22(i,j  ,k)*g22(i,j  ,k) ) )
            gf23(i,j,k)=.5*( ( g12(i,j+1,k)*g13(i,j+1,k) + g22(i,j+1,k)*g23(i,j+1,k) )&
                            +( g12(i,j  ,k)*g13(i,j  ,k) + g22(i,j  ,k)*g23(i,j  ,k) ) )
            gf31(i,j,k)=.5*( ( g11(i,j,k+1)*g13(i,j,k+1) + g21(i,j,k+1)*g23(i,j,k+1) )&
                            +( g11(i,j,k  )*g13(i,j,k  ) + g21(i,j,k  )*g23(i,j,k  ) ) )
            gf32(i,j,k)=.5*( ( g12(i,j,k+1)*g13(i,j,k+1) + g22(i,j,k+1)*g23(i,j,k+1) )&
                            +( g12(i,j,k  )*g13(i,j,k  ) + g22(i,j,k  )*g23(i,j,k  ) ) )
            gf33(i,j,k)=.5*( ( g13(i,j,k+1)*g13(i,j,k+1) + g23(i,j,k+1)*g23(i,j,k+1) &
                             + g33(i,j,k+1)*g33(i,j,k+1) )&
                             +(g13(i,j,k  )*g13(i,j,k  ) + g23(i,j,k  )*g23(i,j,k  ) &
                             + g33(i,j,k  )*g33(i,j,k  ) ) )
            ENDDO
         ENDDO
      ENDDO

!     CALL update(gf11,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf12,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf13,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf21,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf22,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf23,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf31,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf32,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf33,np,mp,lp,np,mp,lp,iup,ih)
#endif
       DO i=1,np,np-1
          IF((i==leftedge).or.(i==(rightedge*np))) then
             ii=1+i/np
             DO k=1,lp
                DO j=1,mp
                   gnr=sqrt(1.+sg21(i,j,k)*sg21(i,j,k))
                   pnx=1./gnr
                   pny=sg21(i,j,k)*pnx
                   tnxx(j,k,ii)=pnx*sg11(i,j,k)+pny*sg21(i,j,k)
                   tnxy(j,k,ii)=pnx*sg12(i,j,k)+pny*sg22(i,j,k)
                   tnxz(j,k,ii)=pnx*sg13(i,j,k)+pny*sg23(i,j,k)
                ENDDO
             ENDDO
          ENDIF
       ENDDO

    !-------> transformed components of the normal to j=1,m surfaces
    !         (NOTE: boundary is defined not a function of z)
       DO j=1,mp,mp-1
          IF((j==botedge).or.(j==(topedge*mp))) then
             jj=1+j/mp
             DO k=1,lp
                DO i=1,np
                   gnr=sqrt(sg12(i,j,k)*sg12(i,j,k)+1.)
                   pny=1./gnr
                   pnx=sg12(i,j,k)*pny
                   tnyx(i,k,jj)=pnx*sg11(i,j,k)+pny*sg21(i,j,k)
                   tnyy(i,k,jj)=pnx*sg12(i,j,k)+pny*sg22(i,j,k)
                   tnyz(i,k,jj)=pnx*sg13(i,j,k)+pny*sg23(i,j,k)
                ENDDO
             ENDDO
          ENDIF
       ENDDO

    !-------> transformed components of the normal to k=1,lp surfaces
       DO k=1,lp,lp-1
          IF((k==gndedge).or.(k==(skyedge*lp))) THEN
             kk=1+k/lp
             DO j=1,mp
                DO i=1,np
                   cx=sg13(i,j,k)
                   cy=sg23(i,j,k)
                   cz= g33(i,j,k)
                   zsx=-cx/cz
                   zsy=-cy/cz
                   gnr=sqrt(1.+zsx**2+zsy**2)
                   pnz=  1./gnr
                   pnx=-zsx*pnz
                   pny=-zsy*pnz
                   tnzx(i,j,kk)=pnx*sg11(i,j,k)+pny*sg21(i,j,k)
                   tnzy(i,j,kk)=pnx*sg12(i,j,k)+pny*sg22(i,j,k)
                   tnzz(i,j,kk)=pnx*sg13(i,j,k)+pny*sg23(i,j,k)+pnz*g33(i,j,k)
                ENDDO
             ENDDO
          ENDIF
       ENDDO

   END SUBROUTINE metryc_lib

!----------------------------------------------------------------------!
   SUBROUTINE metryc
!----------------------------------------------------------------------!

  USE mod_parameters,  ONLY:                     &
         dxih, dyih, dzih, rds, isphere, ideep, ih, ipoles
  USE mpi_parallel, ONLY: update,  updateflt, update3
  USE scratch_datafields, ONLY: esxb,esyb,dsxb,dsyb,    &
                               strxx,strxy,stryx,stryy, &
                                csxb,csyb,cszb,         &
                               strzx,strzy,strzz,scr1 
      !----------------------------------------------------------------------!


      REAL(KIND=euwp) ::  den, &
         esxbij, esybij, dsxbij, dsybij, cszbinv, cszbinvden

      INTEGER(KIND=iintegers) :: iflip, i, j, k

      REAL(KIND=euwp) :: g110, g220, g330

      !----------------------------------------------------------------------!


      DO k=1,lp
         DO j=1,mp
            DO i=1,np
               gmm(i,j,k)=1._euwp+isphere*zcr(i,j,k)/rds*ideep
              gmmi(i,j,k)=1._euwp/gmm(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !CALL update(gmm,np,mp,lp,np,mp,lp,iup,ih)

      !----------------------------------------------------------------------!
      ! Compute Jacobi and inverse Jacobi matrix elements
      !----------------------------------------------------------------------!
      scr1(1:np,1:mp,1:lp)=zcr(1:np,1:mp,1:lp)
      CALL update3(scr1,np,mp,lp,np,mp,lp,iup,ih)
  CALL halo_one_sided_diff_flt_x(xcr,ipoles,np,mp,ih)
  CALL halo_one_sided_diff_flt_y(xcr,ipoles,np,mp,ih)
  CALL halo_one_sided_diff_flt_x(ycr,ipoles,np,mp,ih)
  CALL halo_one_sided_diff_flt_y(ycr,ipoles,np,mp,ih)
  CALL halo_one_sided_diff_x(scr1,np,mp,lp,ih)
  IF(ipoles.eq.0) THEN
    iflip=0
    CALL halo_one_sided_diff_y(scr1,iflip,ipoles,np,mp,lp,ih)
  ELSE
        if (botedge.eq.1) then
           ycr(1:np,0   )=-ycr(1:np,0   )-acos(-1.)*rds
        endif
        if (topedge.eq.1) then
           ycr(1:np,mp+1)=-ycr(1:np,mp+1)+acos(-1.)*rds
        endif
!        print *,'Define ycr for poles=1'
  ENDIF
  CALL halo_one_sided_diff_z(scr1,np,mp,lp,ih)


      !----------------------------------------------------------------------!
      ! Compute dX/dXb
      !----------------------------------------------------------------------!
      DO j=1,mp
         DO i=1,np
            esxb(i,j)=(xcr(i+1,j  )-xcr(i-1,j  ))*dxih  ! dx/d(xb)
            esyb(i,j)=(xcr(i  ,j+1)-xcr(i  ,j-1))*dyih  ! dx/d(yb)
            dsxb(i,j)=(ycr(i+1,j  )-ycr(i-1,j  ))*dxih  ! dy/d(xb)
            dsyb(i,j)=(ycr(i  ,j+1)-ycr(i  ,j-1))*dyih  ! dy/d(yb)
            IF(esxb(i,j).lt.0) THEN
              PRINT *,'esxb(i,j),dxih',esxb(i,j),dxih
              PRINT *,'xcr(i+1,j)-xcr(i-1,j)',xcr(i+1,j)-xcr(i-1,j)
              PRINT *,'xcr(i+1,j)',i,j,xcr(i+1,j)
              PRINT *,'xcr(i  ,j)',i,j,xcr(i,j)
              PRINT *,'xcr(i-1,j)',i,j,xcr(i-1,j)
              PRINT *,'xcr(i-2,j)',i,j,xcr(i-2,j)
              PRINT *,'xcr(1 ,j)',i,j,xcr( 1,j)
              PRINT *,'xcr(0 ,j)',i,j,xcr( 0,j)
              PRINT *,'xcr(-1,j)',i,j,xcr(-1,j)
              PRINT *,'xcr(np  ,j)',i,j,xcr(np,j)
              PRINT *,'xcr(np-1,j)',i,j,xcr(np-1,j)
              STOP 'Negative esxb'
            ENDIF
         ENDDO
      ENDDO
      DO k=1,lp
         DO j=1,mp
            DO i=1,np
               csxb(i,j,k)=(scr1(i+1,j  ,k  )-scr1(i-1,j  ,k  ))*dxih  ! dz/d(xb)
               csyb(i,j,k)=(scr1( i ,j+1,k  )-scr1(i  ,j-1,k  ))*dyih  ! dz/d(yb)
               cszb(i,j,k)=(scr1( i ,j  ,k+1)-scr1(i  ,j  ,k-1))*dzih  ! dz/d(zb)
            ENDDO
         ENDDO
      ENDDO

      DO j=1,mp
         DO i=1,np
!      estb=xcrd(i,j)                                   ! dx/d(tb)
!      dstb=ycrd(i,j)                                   ! dy/d(tb)
            !     use KD identities to invert derivatives -> dXb/dX
            den=1._euwp/(dsyb(i,j)*esxb(i,j)-dsxb(i,j)*esyb(i,j))
            strxx(i,j)= dsyb(i,j)*den                        ! d(xb)/dx
            strxy(i,j)=-esyb(i,j)*den                        ! d(xb)/dy
!      strxd(i,j)=(dstb*esyb(i,j)-dsyb(i,j)*estb)*den   ! d(xb)/dt
            stryx(i,j)=-dsxb(i,j)*den                        ! d(yb)/dx
            stryy(i,j)= esxb(i,j)*den                        ! d(yb)/dy
!      stryd(i,j)=(estb*dsxb(i,j)-esxb(i,j)*dstb)*den   ! d(yb)/dt
      IF(strxx(i,j).lt.0) THEN
        PRINT *,'strxx(i,j)',strxx(i,j)
        PRINT *,' dsyb(i,j)',dsyb(i,j)
        PRINT *,' esxb(i,j)',esxb(i,j)
        PRINT *,' dsyb*esxb(i,j)',dsyb(i,j)*esxb(i,j)
        PRINT *,' dsxb(i,j)',dsxb(i,j)
        PRINT *,' esyb(i,j)',esyb(i,j)
        PRINT *,' dsxb*esyb(i,j)',dsxb(i,j)*esyb(i,j)
        STOP
      ENDIF
         ENDDO
      ENDDO
      DO j=1,mp
         DO i=1,np
!      estb=xcrd(i,j)                                   ! dx/d(tb)
!      dstb=ycrd(i,j)                                   ! dy/d(tb)
            esxbij=esxb(i,j)
            esybij=esyb(i,j)
            dsxbij=dsxb(i,j)
            dsybij=dsyb(i,j)
            den=1._euwp/(dsyb(i,j)*esxb(i,j)-dsxb(i,j)*esyb(i,j))
            DO k=1,lp
!        cstb=zcrd(i,j,k)
               cszbinv=1._euwp/cszb(i,j,k)
               cszbinvden=cszbinv*den
               !       use KD identities to invert derivatives -> dXb/dX
               strzx(i,j,k)=(-csxb(i,j,k)*dsybij+dsxbij*csyb(i,j,k))*cszbinvden
               strzy(i,j,k)=(-esxbij*csyb(i,j,k)+esybij*csxb(i,j,k))*cszbinvden
               strzz(i,j,k)=cszbinv
               ! d(zb)/dt
!        strzd(i,j,k)=((esxbij*dsybij-dsxbij*esybij)*(-cstb) &
!                    + (dsxbij*estb  -esxbij*dstb  )*(-csyb(i,j,k)) &
!                    + (esybij*dstb  -dsybij*estb  )*(-csxb(i,j,k)))*cszbinvden
            ENDDO
         ENDDO
      ENDDO

      DO k=1,lp
         DO j=1,mp
            DO i=1,np
               G110=1._euwp/(gmm(i,j,k)*cosa(i,j))
               G220=1._euwp/gmm(i,j,k)
               G330=1._euwp
              sg11(i,j,k)=strxx(i,j)
              sg12(i,j,k)=stryx(i,j)
              sg13(i,j,k)=strzx(i,j,k)
              sg21(i,j,k)=strxy(i,j)
              sg22(i,j,k)=stryy(i,j)
              sg23(i,j,k)=strzy(i,j,k)
              sg33(i,j,k)=strzz(i,j,k)
               g11(i,j,k)=strxx(i,j)*g110
               g12(i,j,k)=stryx(i,j)*g110
               g13(i,j,k)=strzx(i,j,k)*g110
               g21(i,j,k)=strxy(i,j)*g220
               g22(i,j,k)=stryy(i,j)*g220
               g23(i,j,k)=strzy(i,j,k)*g220
               g33(i,j,k)=strzz(i,j,k)*g330
               gac(i,j,k)=1._euwp/strzz(i,j,k)*(gmm(i,j,k)**2*(cosa(i,j)))    &
                  /(stryy(i,j)*strxx(i,j)-stryx(i,j)*strxy(i,j))
              gaci(i,j,k)=1._euwp/gac(i,j,k)
    g11g22Mg12g21i(i,j,k)= 1._euwp/ &
                         (g11(i,j,k)*g22(i,j,k)-g12(i,j,k)*g21(i,j,k))
        g33i(i,j,k)= 1._euwp/g33(i,j,k)
      IF(cosa(i,j).lt.0) THEN
        PRINT *,'cosa(i,j)<0',i,j,cosa(i,j)
!        STOP
      ENDIF

            ENDDO
         ENDDO
      ENDDO 

      CALL halo_copy_edge_to_halo(g11,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g12,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g13,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g21,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g22,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g23,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g33,np,mp,lp,ih)
      

      CALL update(zcr,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g11,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g12,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g13,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g21,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g22,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g23,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(g33,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(gac,np,mp,lp,np,mp,lp,iup,ih)
      CALL update(gmm,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg11,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg12,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg13,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg21,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg22,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg23,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(sg33,np,mp,lp,np,mp,lp,iup,ih)     

      CALL halo_copy_edge_to_halo(g11,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g12,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g13,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g21,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g22,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g23,np,mp,lp,ih)
      CALL halo_copy_edge_to_halo(g33,np,mp,lp,ih)

      
#if(1==1)
      DO k=1,lp  
         DO j=1,mp  
            DO i=1,np  
            gf11(i,j,k)=.5*( ( g11(i+1,j,k)*g11(i+1,j,k) + g21(i+1,j,k)*g21(i+1,j,k) )&
                            +( g11(i  ,j,k)*g11(i  ,j,k) + g21(i  ,j,k)*g21(i  ,j,k) ) )
            gf12(i,j,k)=.5*( ( g11(i+1,j,k)*g12(i+1,j,k) + g21(i+1,j,k)*g22(i+1,j,k) )&
                            +( g11(i  ,j,k)*g12(i  ,j,k) + g21(i  ,j,k)*g22(i  ,j,k) ) )
            gf13(i,j,k)=.5*( ( g11(i+1,j,k)*g13(i+1,j,k) + g21(i+1,j,k)*g23(i+1,j,k) )&
                            +( g11(i  ,j,k)*g13(i  ,j,k) + g21(i  ,j,k)*g23(i  ,j,k) ) )
            gf21(i,j,k)=.5*( ( g11(i,j+1,k)*g12(i,j+1,k) + g21(i,j+1,k)*g22(i,j+1,k) )&
                            +( g11(i,j  ,k)*g12(i,j  ,k) + g21(i,j  ,k)*g22(i,j  ,k) ) )
            gf22(i,j,k)=.5*( ( g12(i,j+1,k)*g12(i,j+1,k) + g22(i,j+1,k)*g22(i,j+1,k) )&
                            +( g12(i,j  ,k)*g12(i,j  ,k) + g22(i,j  ,k)*g22(i,j  ,k) ) )
            gf23(i,j,k)=.5*( ( g12(i,j+1,k)*g13(i,j+1,k) + g22(i,j+1,k)*g23(i,j+1,k) )&
                            +( g12(i,j  ,k)*g13(i,j  ,k) + g22(i,j  ,k)*g23(i,j  ,k) ) )
            gf31(i,j,k)=.5*( ( g11(i,j,k+1)*g13(i,j,k+1) + g21(i,j,k+1)*g23(i,j,k+1) )&
                            +( g11(i,j,k  )*g13(i,j,k  ) + g21(i,j,k  )*g23(i,j,k  ) ) )
            gf32(i,j,k)=.5*( ( g12(i,j,k+1)*g13(i,j,k+1) + g22(i,j,k+1)*g23(i,j,k+1) )&
                            +( g12(i,j,k  )*g13(i,j,k  ) + g22(i,j,k  )*g23(i,j,k  ) ) )
            gf33(i,j,k)=.5*( ( g13(i,j,k+1)*g13(i,j,k+1) + g23(i,j,k+1)*g23(i,j,k+1) &
                             + g33(i,j,k+1)*g33(i,j,k+1) )&
                             +(g13(i,j,k  )*g13(i,j,k  ) + g23(i,j,k  )*g23(i,j,k  ) &
                             + g33(i,j,k  )*g33(i,j,k  ) ) )
            ENDDO
         ENDDO
      ENDDO

!     CALL update(gf11,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf12,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf13,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf21,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf22,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf23,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf31,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf32,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf33,np,mp,lp,np,mp,lp,iup,ih)
#endif
       DO i=1,np,np-1
          IF((i==leftedge).or.(i==(rightedge*np))) then
             ii=1+i/np
             DO k=1,lp
                DO j=1,mp
                   gnr=sqrt(1.+sg21(i,j,k)*sg21(i,j,k))
                   pnx=1./gnr
                   pny=sg21(i,j,k)*pnx
                   tnxx(j,k,ii)=pnx*sg11(i,j,k)+pny*sg21(i,j,k)
                   tnxy(j,k,ii)=pnx*sg12(i,j,k)+pny*sg22(i,j,k)
                   tnxz(j,k,ii)=pnx*sg13(i,j,k)+pny*sg23(i,j,k)
                ENDDO
             ENDDO
          ENDIF
       ENDDO

    !-------> transformed components of the normal to j=1,m surfaces
    !         (NOTE: boundary is defined not a function of z)
       DO j=1,mp,mp-1
          IF((j==botedge).or.(j==(topedge*mp))) then
             jj=1+j/mp
             DO k=1,lp
                DO i=1,np
                   gnr=sqrt(sg12(i,j,k)*sg12(i,j,k)+1.)
                   pny=1./gnr
                   pnx=sg12(i,j,k)*pny
                   tnyx(i,k,jj)=pnx*sg11(i,j,k)+pny*sg21(i,j,k)
                   tnyy(i,k,jj)=pnx*sg12(i,j,k)+pny*sg22(i,j,k)
                   tnyz(i,k,jj)=pnx*sg13(i,j,k)+pny*sg23(i,j,k)
                ENDDO
             ENDDO
          ENDIF
       ENDDO

    !-------> transformed components of the normal to k=1,lp surfaces
       DO k=1,lp,lp-1
          IF((k==gndedge).or.(k==(skyedge*lp))) THEN
             kk=1+k/lp
             DO j=1,mp
                DO i=1,np
                   cx=sg13(i,j,k)
                   cy=sg23(i,j,k)
                   cz= g33(i,j,k)
                   zsx=-cx/cz
                   zsy=-cy/cz
                   gnr=sqrt(1.+zsx**2+zsy**2)
                   pnz=  1./gnr
                   pnx=-zsx*pnz
                   pny=-zsy*pnz
                   tnzx(i,j,kk)=pnx*sg11(i,j,k)+pny*sg21(i,j,k)
                   tnzy(i,j,kk)=pnx*sg12(i,j,k)+pny*sg22(i,j,k)
                   tnzz(i,j,kk)=pnx*sg13(i,j,k)+pny*sg23(i,j,k)+pnz*g33(i,j,k)
                ENDDO
             ENDDO
          ENDIF
       ENDDO

   END SUBROUTINE metryc

   SUBROUTINE diffmetrics(eta) 
   USE geometry_datafields, ONLY: gf11V, gf12V, gf13V
   USE geometry_datafields, ONLY: gf21V, gf22V, gf23V
   USE geometry_datafields, ONLY: gf31V, gf32V, gf33V
   USE mpi_parallel, ONLY: update
   REAL(KIND=euwp), DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,3), INTENT(IN) :: eta 
      INTEGER(KIND=iintegers) :: i, j, k




      DO k=1,lp  
         DO j=1,mp  
            DO i=1,np  
            gf11V(i,j,k)=.5*( ( g11(i+1,j,k)**2*eta(i+1,j,k,1) + g21(i+1,j,k)**2*eta(i+1,j,k,2) )&
                             +( g11(i  ,j,k)**2*eta(i  ,j,k,1) + g21(i  ,j,k)**2*eta(i  ,j,k,2) ) )
            gf12V(i,j,k)=.5*( ( g11(i+1,j,k)*g12(i+1,j,k)*eta(i+1,j,k,1)  &
                               +g21(i+1,j,k)*g22(i+1,j,k)*eta(i+1,j,k,2) )&
                             +( g11(i  ,j,k)*g12(i  ,j,k)*eta(i  ,j,k,1)   &
                               +g21(i  ,j,k)*g22(i  ,j,k)*eta(i  ,j,k,2) ) )
            gf13V(i,j,k)=.5*( ( g11(i+1,j,k)*g13(i+1,j,k)*eta(i+1,j,k,1)  &
                               +g21(i+1,j,k)*g23(i+1,j,k)*eta(i+1,j,k,2) )&
                             +( g11(i  ,j,k)*g13(i  ,j,k)*eta(i  ,j,k,1)  &
                               +g21(i  ,j,k)*g23(i  ,j,k)*eta(i  ,j,k,2) ) )
            gf21V(i,j,k)=.5*( ( g11(i,j+1,k)*g12(i,j+1,k)*eta(i,j+1,k,1)  &
                               +g21(i,j+1,k)*g22(i,j+1,k)*eta(i,j+1,k,2) )&
                             +( g11(i,j  ,k)*g12(i,j  ,k)*eta(i,j  ,k,1)   &
                               +g21(i,j  ,k)*g22(i,j  ,k)*eta(i,j  ,k,2) ) )
            gf22V(i,j,k)=.5*( ( g12(i,j+1,k)**2*eta(i,j+1,k,1) + g22(i,j+1,k)**2*eta(i,j+1,k,2) )&
                             +( g12(i,j  ,k)**2*eta(i,j  ,k,1) + g22(i,j  ,k)**2*eta(i,j,  k,2) ) )
            gf23V(i,j,k)=.5*( ( g12(i,j+1,k)*g13(i,j+1,k)*eta(i,j+1,k,1)  &
                               +g22(i,j+1,k)*g23(i,j+1,k)*eta(i,j+1,k,2) )&
                             +( g12(i,j  ,k)*g13(i,j  ,k)*eta(i,j  ,k,1)  &
                               +g22(i,j  ,k)*g23(i,j  ,k)*eta(i,j  ,k,2) ) )
            gf31V(i,j,k)=.5*( ( g11(i,j,k+1)*g13(i,j,k+1)*eta(i,j,k+1,1) &
                               +g21(i,j,k+1)*g23(i,j,k+1)*eta(i,j,k+1,2) )&
                             +( g11(i,j,k  )*g13(i,j,k  )*eta(i,j,k  ,1)  &
                               +g21(i,j,k  )*g23(i,j,k  )*eta(i,j,k  ,2) ) )
            gf32V(i,j,k)=.5*( ( g12(i,j,k+1)*g13(i,j,k+1)*eta(i,j,k+1,1) &
                               +g22(i,j,k+1)*g23(i,j,k+1)*eta(i,j,k+1,2) )&
                             +( g12(i,j,k  )*g13(i,j,k  )*eta(i,j,k  ,1)  &
                               +g22(i,j,k  )*g23(i,j,k  )*eta(i,j,k  ,2) ) )
            gf33V(i,j,k)=.5*( ( g13(i,j,k+1)**2*eta(i,j,k+1,1)+ g23(i,j,k+1)**2*eta(i,j,k+1,2) &
                              + g33(i,j,k+1)**2*eta(i,j,k+1,3) )&
                              +(g13(i,j,k  )**2*eta(i,j,k  ,1)+ g23(i,j,k  )**2*eta(i,j,k  ,2) &
                              + g33(i,j,k  )**2*eta(i,j,k  ,3) ) )
            ENDDO
         ENDDO
      ENDDO

!     CALL update(gf11V,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf12V,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf13V,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf21V,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf22V,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf23V,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf31V,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf32V,np,mp,lp,np,mp,lp,iup,ih)
!     CALL update(gf33V,np,mp,lp,np,mp,lp,iup,ih)
   END SUBROUTINE diffmetrics 

   SUBROUTINE physical2contravariant(u_phys,v_phys,w_phys,ox_contr,oy_contr,oz_contr)
      INTEGER (KIND=IINTEGERS) :: i,j,k
   REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &
   u_phys,v_phys,w_phys
   REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(OUT) :: &
   ox_contr,oy_contr,oz_contr
      DO k=1,lp
        DO j=1,mp
          DO i=1,np
            ox_contr(i,j,k)=g11(i,j,k)*u_phys(i,j,k)+g21(i,j,k)*v_phys(i,j,k)
            oy_contr(i,j,k)=g12(i,j,k)*u_phys(i,j,k)+g22(i,j,k)*v_phys(i,j,k)
            oz_contr(i,j,k)=g13(i,j,k)*u_phys(i,j,k)+g23(i,j,k)*v_phys(i,j,k)+g33(i,j,k)*w_phys(i,j,k)
          ENDDO
        ENDDO
      ENDDO

   END SUBROUTINE physical2contravariant

   SUBROUTINE contravariant2physical(ox_contr,oy_contr,oz_contr,u_phys,v_phys,w_phys)
      INTEGER (KIND=IINTEGERS) :: i,j,k
   REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &
   ox_contr,oy_contr,oz_contr
   REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(OUT) :: &
   u_phys,v_phys,w_phys
      DO k=1,lp
        DO j=1,mp
          DO i=1,np
      u_phys(i,j,k)= (g22(i,j,k)*ox_contr(i,j,k)-g21(i,j,k)*oy_contr(i,j,k)) &
                    /(g22(i,j,k)*g11(i,j,k)-g12(i,j,k)*g21(i,j,k))
!                      *g11g22Mg12g21i(i,j,k)

      v_phys(i,j,k)=-(g12(i,j,k)*ox_contr(i,j,k)-g11(i,j,k)*oy_contr(i,j,k)) &
!                         *g11g22Mg12g21i(i,j,k)
                    /(g22(i,j,k)*g11(i,j,k)-g12(i,j,k)*g21(i,j,k))
      w_phys(i,j,k)=( oz_contr(i,j,k)-g13(i,j,k)*u_phys(i,j,k)  &
                                     -g23(i,j,k)*v_phys(i,j,k)) &
                                                            /g33(i,j,k)
!                          *g33i(i,j,k)

          ENDDO
        ENDDO
      ENDDO
   END SUBROUTINE contravariant2physical


   SUBROUTINE calc_dn(dn, i,j,k, ii,jj,kk)
       INTEGER(KIND=iintegers) :: i, ii
       INTEGER(KIND=iintegers) :: j, jj
       INTEGER(KIND=iintegers) :: k, kk
       REAL (KIND=euwp) :: dn

       dn=  tnxx(j,k,ii)*tnyy(i,k,jj)*tnzz(i,j,kk) &
           +tnxy(j,k,ii)*tnyz(i,k,jj)*tnzx(i,j,kk) &
           +tnxz(j,k,ii)*tnyx(i,k,jj)*tnzy(i,j,kk) &
           -tnxz(j,k,ii)*tnyy(i,k,jj)*tnzx(i,j,kk) &
           -tnxx(j,k,ii)*tnyz(i,k,jj)*tnzy(i,j,kk) &
           -tnxy(j,k,ii)*tnyx(i,k,jj)*tnzz(i,j,kk)
   END SUBROUTINE calc_dn

   SUBROUTINE calc_d1(d1, i,j,k, ii,jj,kk, qx,qy,qz)
       INTEGER(KIND=iintegers) :: i, ii
       INTEGER(KIND=iintegers) :: j, jj
       INTEGER(KIND=iintegers) :: k, kk
       REAL (KIND=euwp) :: qx, qy, qz
       REAL (KIND=euwp) :: d1

       d1= -          qx*tnyy(i,k,jj)*tnzz(i,j,kk)&
           -tnxy(j,k,ii)*tnyz(i,k,jj)*qz&
           -tnxz(j,k,ii)*          qy*tnzy(i,j,kk)&
           +tnxz(j,k,ii)*tnyy(i,k,jj)*qz&
           +          qx*tnyz(i,k,jj)*tnzy(i,j,kk)&
           +tnxy(j,k,ii)*          qy*tnzz(i,j,kk)
   END SUBROUTINE calc_d1
   
   SUBROUTINE calc_d2(d2, i,j,k, ii,jj,kk, qx,qy,qz)
       INTEGER(KIND=iintegers) :: i, ii
       INTEGER(KIND=iintegers) :: j, jj
       INTEGER(KIND=iintegers) :: k, kk
       REAL (KIND=euwp) :: qx, qy, qz
       REAL (KIND=euwp) :: d2

       d2= -tnxx(j,k,ii)*          qy*tnzz(i,j,kk)&
           -          qx*tnyz(i,k,jj)*tnzx(i,j,kk)&
           -tnxz(j,k,ii)*tnyx(i,k,jj)*qz&
           +tnxz(j,k,ii)*          qy*tnzx(i,j,kk)&
           +tnxx(j,k,ii)*tnyz(i,k,jj)*qz&
           +          qx*tnyx(i,k,jj)*tnzz(i,j,kk)
   END SUBROUTINE calc_d2

   SUBROUTINE calc_d3(d3, i,j,k, ii,jj,kk, qx,qy,qz)
       INTEGER(KIND=iintegers) :: i, ii
       INTEGER(KIND=iintegers) :: j, jj
       INTEGER(KIND=iintegers) :: k, kk
       REAL (KIND=euwp) :: qx, qy, qz
       REAL (KIND=euwp),INTENT(OUT) :: d3

       d3= -tnxx(j,k,ii)*tnyy(i,k,jj)*qz&
           -tnxy(j,k,ii)*          qy*tnzx(i,j,kk)&
           -          qx*tnyx(i,k,jj)*tnzy(i,j,kk)&
           +          qx*tnyy(i,k,jj)*tnzx(i,j,kk)&
           +tnxx(j,k,ii)*          qy*tnzy(i,j,kk)&
           +tnxy(j,k,ii)*tnyx(i,k,jj)*qz
   END SUBROUTINE calc_d3

END MODULE geometry



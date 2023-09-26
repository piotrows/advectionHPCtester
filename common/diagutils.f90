#include "../advection/src_algorithms/renames.inc"
MODULE eulag_diagutils
USE precisions
USE mod_parameters, ONLY: n,m,nm,nml
USE mpi_parallel, ONLY: ttbeg,ttend

IMPLICIT NONE

CONTAINS
#include "../advection/src_algorithms/defines.inc"
#if(1==0)
!---------------------------------------------------------------------!
! Filter EULAG`s meteorological field 'a', filtering is applied on    !
! the model levels
!---------------------------------------------------------------------!
SUBROUTINE filstrP(a,nullfl)
!---------------------------------------------------------------------!
  
  ! COSMO modules
  ! EULAG modules
  USE data_euconstants, ONLY: &
    dt,dx,dy,dz,dti,dxi,dyi,dzi,zb,igrid,j3, &
    ibcx,ibcy,ibcz,irlx,irly
  
  USE kind_mod_parameters, ONLY: euwp
  
  USE src_parallel, ONLY: updatelr, updatebt, update
  
  ! include 'param.misc'
  
  ! Subroutine arguments
  ! Local variables
  REAL(KIND=euwp) :: a(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
  
  INTEGER(KIND=iintegers) :: nullfl
  
  
  REAL(KIND=euwp),DIMENSION(np) :: sx
  REAL(KIND=euwp),DIMENSION(mp) :: sy
  REAL(KIND=euwp),DIMENSION(l) :: sz
  
  INTEGER(KIND=iintegers) :: illim, iulim, ibox, iboy, iboz, i, j, k, jllim, julim
  
  REAL(KIND=euwp) :: temp, tempo
  
  !----------------
  ! C O N T E N T :
  !----------------
  IF (nullfl == 0) RETURN
  
  CALL updatelr(a,np,mp,lp,np,mp,lp,ih)
  
  illim = 1  + 1*leftedge
  iulim = np - 1*rightedge
  ibox=1-ibcx
  iboy=1-ibcy
  iboz=1-ibcz
  
  DO k=1,l
    DO j=1,mp
      DO i=illim,iulim
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
  
  
  IF (j3 == 1) THEN
    CALL updatebt(a,np,mp,lp,np,mp,lp,ih)
    jllim=1  + j3*botedge
    julim=mp - j3*topedge
    
    DO k=1,l
      DO i=1,np
        DO j=jllim,julim
          sy(j)=0.25_euwp*(a(i,j+j3,k)+2._euwp*a(i,j,k)+a(i,j-j3,k))
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
  ENDIF
  
  DO j=1,mp
    DO i=1,np
      DO k=2,lp-1
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
  
  CALL update(a,np,mp,lp,np,mp,lp,ih)
  
END SUBROUTINE filstrP


!---------------------------------------------------------------------!
!
!---------------------------------------------------------------------!
SUBROUTINE ave_fld(f,fav)
!---------------------------------------------------------------------!
  USE data_eufields,    ONLY: zcr,zstr,bufsum
  USE src_parallel,     ONLY: globsumv 
  USE data_runcontrol,   ONLY:  ntstep
  
!---------------------------------------------------------------------!
  
  REAL(KIND=euwp) :: f(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
  REAL(KIND=euwp) :: fav(lp),zpnt(lp),wgt(lp),rcnt(lp)

!---------------------------------------------------------------------!

  REAL(KIND=euwp) coe2 
  INTEGER kk,i,j,k
  
!---------------------------------------------------------------------!

  fav(:)=0._euwp
  rcnt(:)=0._euwp

  DO i=1,np
    DO j=1,mp
      !Here set the treshold to cutoff the high topography influence
      !if (zs(i,j).lt.10.) then
      zpnt(1:lp)=zcr(i,j,1:lp) 
      
      DO k=1,lp
        IF (zstr(k) >= zcr(i,j,1)) THEN
          kk=2
          DO WHILE (zpnt(kk) < zstr(k).AND.kk < l)
            kk=kk+1 ! Find the nearest base level kk higher than local level
          END DO
          coe2=(zstr(k)-zpnt(kk-1))/(zpnt(kk)-zpnt(kk-1))
          fav(k)=fav(k)+f(i,j,kk-1) + coe2*(f(i,j,kk)-f(i,j,kk-1))
          rcnt(k)=rcnt(k)+1
        ENDIF
      ENDDO
      !endif !treshold
    END DO
  END DO
!      DO k=1,lp
!       fav(k)=SUM(f(1:np,1:mp,k))
!       rcnt(k)=np*mp
!      ENDDO 

  bufsum(1:lp)=fav(1:lp)
  bufsum(lp+1:2*lp)=rcnt(1:lp)
  CALL globsumv(2*lp)
  fav(1:lp)=bufsum(1:lp)
  fav(2:lp)=fav(2:lp)/bufsum(lp+2:2*lp)
  fav(1)=fav(2)
  bufsum(:)=0._euwp
!  IF(mype == 0.AND.ntstep == 0) print *,'Fav:', fav

END SUBROUTINE ave_fld


#endif /*1==0*/
!---------------------------------------------------------------------!
!
!---------------------------------------------------------------------!
SUBROUTINE calcavgmaxminflt(x,xmax,xmin,xave,np,mp,ih)
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,ih 
  REAL(KIND=euwp), INTENT(IN) :: x(1-ih:np+ih,1-ih:mp+ih)
  REAL(KIND=euwp), INTENT(OUT) :: xave,xmin,xmax
  REAL(KIND=euwp)  :: xloc
  INTEGER(KIND=iintegers) :: i,j
  xave=0._euwp
  xmax=-HUGE(xmax)
  xmin= HUGE(xmin)
  
  DO j=1,mp
    DO i=1,np
      xloc=x(i,j)
      !Computing x average
      xave=xave+xloc
      !Computing x maximum 
      xmax=max(xmax,xloc)
      !Computing x minimum 
      xmin=min(xmin,xloc)
    ENDDO
  END DO
  xave=xave/REAL(n*m)
END SUBROUTINE calcavgmaxminflt


!---------------------------------------------------------------------!
!
!---------------------------------------------------------------------!
SUBROUTINE calcavgmaxmin(x,xmax,xmin,xave,np,mp,lp,ih)
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih 
  REAL(KIND=euwp), INTENT(IN) :: &
                 x(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
  REAL(KIND=euwp) :: xave,xmin,xmax
  REAL(KIND=euwp) :: xloc
  INTEGER(KIND=iintegers) :: i,j,k
  xave=0._euwp
  xmax=-HUGE(xmax)
  xmin= HUGE(xmin) 
  
  DO k=1,lp
    DO j=1,mp
      DO i=1,np
        xloc=x(i,j,k)
        !Computing x average
        xave=xave+xloc
        !Computing x maximum 
        xmax=max(xmax,xloc)
        !Computing x minimum 
        xmin=min(xmin,xloc)
      ENDDO
    ENDDO
  ENDDO
  xave=xave/REAL(nml)
END SUBROUTINE calcavgmaxmin
SUBROUTINE calcavgmaxmin_vz(x,xmaxv,xminv,xavev,np,mp,lp,ih)
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih 
  REAL(KIND=euwp), INTENT(IN) :: &
                 x(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
  REAL(KIND=euwp), INTENT(OUT) :: xavev(lp),xminv(lp),xmaxv(lp)
  REAL(KIND=euwp) :: xloc,rnm
  INTEGER(KIND=iintegers) :: i,j,k
  rnm=1._euwp/REAL(n*m)
  DO k=1,lp
  xavev(k)=0._euwp
  xmaxv(k)=-HUGE(xmaxv(k))
  xminv(k)= HUGE(xminv(k)) 
    DO j=1,mp
      DO i=1,np
        xloc=x(i,j,k)
        !Computing x average
        xavev(k)=xavev(k)+xloc
        !Computing x maximum 
        xmaxv(k)=max(xmaxv(k),xloc)
        !Computing x minimum 
        xminv(k)=min(xminv(k),xloc)
      ENDDO
    ENDDO
  xavev(k)=xavev(k)*rnm
  ENDDO
END SUBROUTINE calcavgmaxmin_vz

!---------------------------------------------------------------------!
!
!---------------------------------------------------------------------!
SUBROUTINE wrttobuf(xmx,xmn,xav,bufavg,bufmax,bufmin,indx)
!---------------------------------------------------------------------!
  REAL(KIND=euwp),INTENT(OUT) :: bufavg(:),bufmax(:),bufmin(:)
  REAL(KIND=euwp),INTENT(IN) :: xmx,xmn,xav
  INTEGER, INTENT(IN) :: indx
  bufavg(indx)=xav
  bufmax(indx)=xmx
  bufmin(indx)=xmn
END SUBROUTINE wrttobuf


!---------------------------------------------------------------------!
!
!---------------------------------------------------------------------!
SUBROUTINE wrtfmbuf(xmx,xmn,xav,bufavg,bufmax,bufmin,indx)
!---------------------------------------------------------------------!
  REAL(KIND=euwp),INTENT(IN) :: bufavg(:),bufmax(:),bufmin(:)
  REAL(KIND=euwp),INTENT(OUT) :: xmx,xmn,xav
  INTEGER, INTENT(IN)  :: indx
  xav=bufavg(indx)
  xmx=bufmax(indx)
  xmn=bufmin(indx)
END SUBROUTINE wrtfmbuf

!---------------------------------------------------------------------!
!
!---------------------------------------------------------------------!
SUBROUTINE scalcavgmaxmin(x,xmax,xmin,xave,np,mp,lp)
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp 
  REAL(KIND=euwp), INTENT(IN) :: x(np,mp,lp)
  REAL(KIND=euwp) xave,xmin,xmax,xloc
  INTEGER(KIND=iintegers) :: i,j,k

  xave=0._euwp
  xmax=-HUGE(xmax)
  xmin= HUGE(xmin)
  
  DO k=1,lp
    DO j=1,mp
      DO i=1,np
        xloc=x(i,j,k)
        !Computing x average
        xave=xave+xloc
        !Computing x maximum 
        xmax=max(xmax,xloc)
        !Computing x minimum 
        xmin=min(xmin,xloc)
      ENDDO
    ENDDO
  ENDDO
  xave=xave/nml
END SUBROUTINE scalcavgmaxmin


!---------------------------------------------------------------------!
!
!---------------------------------------------------------------------!
SUBROUTINE scalcfltavgmaxmin(x,xmax,xmin,xave,np,mp)
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN) :: np,mp 
  REAL(KIND=euwp),INTENT(IN) :: x(np,mp) 
  REAL(KIND=euwp) xave,xmin,xmax,xloc
  INTEGER(KIND=iintegers) :: i,j

  xave=0._euwp
  xmax=-HUGE(xmax)
  xmin= HUGE(xmin) 
  
  DO j=1,mp
    DO i=1,np
      xloc=x(i,j)
      !Computing x average
      xave=xave+xloc
      !Computing x maximum 
      xmax=max(xmax,xloc)
      !Computing x minimum 
      xmin=min(xmin,xloc)
    ENDDO
  ENDDO
  xave=xave/nm
END SUBROUTINE scalcfltavgmaxmin

!---------------------------------------------------------------------!
SUBROUTINE calcfltavgmaxmin(x,xmax,xmin,xave,np,mp,ih)
!---------------------------------------------------------------------!
  INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,ih 
  REAL(KIND=euwp),INTENT(IN) :: x(1-ih:np+ih,1-ih:mp+ih) 
  REAL(KIND=euwp) xave,xmin,xmax,xloc
  INTEGER(KIND=iintegers) :: i,j

  xave=0._euwp
  xmax=-HUGE(xmax)
  xmin= HUGE(xmin) 
  
  DO j=1,mp
    DO i=1,np
      xloc=x(i,j)
      !Computing x average
      xave=xave+xloc
      !Computing x maximum 
      xmax=max(xmax,xloc)
      !Computing x minimum 
      xmin=min(xmin,xloc)
    ENDDO
  ENDDO
  xave=xave/nm
END SUBROUTINE calcfltavgmaxmin

!---------------------------------------------------------------------!
!
!---------------------------------------------------------------------!
SUBROUTINE calcsumaxminloc(x,varname,buffer,bufnames,indx,np,mp,lp,ih)
!---------------------------------------------------------------------!
  USE mpi_parallel, ONLY: nsubpos,msubpos,lsubpos
  INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih 
  REAL(KIND=euwp) x(1-ih:np+ih, 1-ih:mp+ih,1-ih:lp+ih)
  REAL(KIND=euwp)  xmax,xmin,xloc,xave,xsd,xmxcnt,xmncnt,xavcnt,xlocd
  INTEGER(KIND=iintegers)  ixyz(3)
  INTEGER(KIND=iintegers) :: i,j,k,indx
  INTEGER(KIND=iintegers),PARAMETER :: idmaxvars=40
  REAL(KIND=euwp) buffer(15,idmaxvars)
  character(5) bufnames(idmaxvars)
  character(5) varname
  xave= 0._euwp
  xsd= 0._euwp
  xmax=-HUGE(xmax)
  xmin= HUGE(xmin)
  xavcnt= 0._euwp
  xmxcnt= 0._euwp
  xmncnt= 0._euwp
  DO k=1,lp
    DO j=1,mp
      DO i=1,np
        xloc=x(i,j,k)
        !Computing x average
        xave=xave+xloc
        if(xmax == xloc) xmxcnt=xmxcnt+1
        if(xmin == xloc) xmncnt=xmncnt+1
        !Computing x maximum
        xmax=max(xmax,xloc)
        !Computing x minimum 
        xmin=min(xmin,xloc)
      END DO
    END DO
  END DO
  DO k=1,lp
    DO j=1,mp
      DO i=1,np
        xlocd=x(i,j,k)-xave
        xsd=xsd+xlocd*xlocd
      END DO
    END DO
  END DO
  xave=xave/nml
  buffer(1,indx)=xave
  buffer(2,indx)=xmax
  buffer(3,indx)=xmin
  ixyz=maxloc(x(1:np,1:mp,1:lp))
  buffer(4,indx)=FLOAT(nsubpos+ixyz(1))
  buffer(5,indx)=FLOAT(msubpos+ixyz(2))
  buffer(6,indx)=FLOAT(lsubpos+ixyz(3))
  buffer(7,indx)=FLOAT(ixyz(1))
  buffer(8,indx)=FLOAT(ixyz(2))
  buffer(9,indx)=FLOAT(ixyz(3))
  ixyz=minloc(x(1:np,1:mp,1:lp))
  buffer(10,indx)=FLOAT(nsubpos+ixyz(1))
  buffer(11,indx)=FLOAT(msubpos+ixyz(2))
  buffer(12,indx)=FLOAT(lsubpos+ixyz(3))
  buffer(13,indx)=FLOAT(ixyz(1))
  buffer(14,indx)=FLOAT(ixyz(2))
  buffer(15,indx)=FLOAT(ixyz(3))
  bufnames(indx)=varname
  !      if(mype.eq.0) print *, 'Checking symmetry for: ',varname
!  CALL mybarrier()
  !      CALL ckcycd (x,varname)
!  CALL mybarrier()
  !      if(mype.eq.0) print *, 'End of checking symmetry for: ',varname
  IF (varname /= 'nnnnn') THEN
    !      CALL check_symmetry(x,varname,'calcsumaxminloc')
  END IF
END SUBROUTINE calcsumaxminloc

SUBROUTINE compute_courlipsch_Agrid(cr1,cr2,ox,oy,oz,h,dxi,dyi,dzi,dt,np,mp,lp,ih) 
  USE mpi_parallel, ONLY:globmax,update3,iup
  USE mod_parameters, ONLY: ipoles,ibcx,ibcy,ibcz
  USE bconditions, ONLY: remove_cyclic_offset_full,cp_scnd_last_to_halo_xyz_full
  INTEGER(KIND=iintegers), INTENT(IN) :: np,mp,lp,ih
  REAL_euwp, INTENT(IN) ::                            &
                   ox(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),     &
                   oy(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),     &
                   oz(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2)
  REAL_euwp, INTENT(IN) ::  &
                    h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
  REAL(KIND=euwp),INTENT(IN) :: dxi,dyi,dzi,dt
  REAL(KIND=euwp),INTENT(OUT) :: cr1,cr2 
  REAL(KIND=euwp) :: bufmax_in(2),bufmax_out(2)
  REAL(KIND=euwp) :: gc1,gc2,gc3,dinv 
  INTEGER(KIND=iintegers) :: i,j,k

  CALL ttbeg(401)
  CALL ttbeg(1401)
  gc1=dt*dxi
  gc2=dt*dyi
  gc3=dt*dzi

    cr1=-HUGE(cr1)
    cr2=-HUGE(cr2)
  CALL update3(ox(:,:,:,1),np,mp,lp,np,mp,lp,iup,ih)
  CALL update3(oy(:,:,:,1),np,mp,lp,np,mp,lp,iup,ih)
  CALL update3(oz(:,:,:,1),np,mp,lp,np,mp,lp,iup,ih)
  CALL     remove_cyclic_offset_full(ox(:,:,:,1),ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
  CALL     remove_cyclic_offset_full(oy(:,:,:,1),ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
  CALL     remove_cyclic_offset_full(oz(:,:,:,1),ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
  CALL cp_scnd_last_to_halo_xyz_full(ox(:,:,:,1),ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
  CALL cp_scnd_last_to_halo_xyz_full(oy(:,:,:,1),ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
  CALL cp_scnd_last_to_halo_xyz_full(oz(:,:,:,1),ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)

    DO k=1,lp
      DO j=1,mp
        DO i=1,np
          dinv=1._euwp/h(i,j,k)
          cr1=  max(cr1,(abs(ox(i,j,k,1)) +            &
                         abs(oy(i,j,k,1)) +            &
                         abs(oz(i,j,k,1)))*dinv)
          cr2=  max(cr2,                                                     &
             .5_euwp*gc1/gc1*abs(ox(i+1,j  ,k  ,1)-ox(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc1*abs(ox(i  ,j+1,k  ,1)-ox(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc1*abs(ox(i  ,j  ,k+1,1)-ox(i  ,j  ,k-1,1))*dinv,  &
             .5_euwp*gc1/gc2*abs(oy(i+1,j  ,k  ,1)-oy(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc2*abs(oy(i  ,j+1,k  ,1)-oy(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc2*abs(oy(i  ,j  ,k+1,1)-oy(i  ,j  ,k-1,1))*dinv,  &
             .5_euwp*gc1/gc3*abs(oz(i+1,j  ,k  ,1)-oz(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc3*abs(oz(i  ,j+1,k  ,1)-oz(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc3*abs(oz(i  ,j  ,k+1,1)-oz(i  ,j  ,k-1,1))*dinv )
        END DO
      END DO
    END DO
  CALL ttend(1401)
  bufmax_in(1)=cr1
  bufmax_in(2)=cr2
  CALL globmax(bufmax_in,bufmax_out,2)
  cr1=bufmax_out(1)
  cr2=bufmax_out(2)

  CALL ttend(401)

END SUBROUTINE compute_courlipsch_Agrid

SUBROUTINE compute_courlipsch_Agrid_full(pfx,pfy,ox,oy,oz,h,dxi,dyi,dzi,dt,np,mp,lp,ih) 
  USE mpi_parallel, ONLY:globmax
  INTEGER(KIND=iintegers), INTENT(IN) :: np,mp,lp,ih
  REAL_euwp, INTENT(IN) ::                            &
                   ox(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),     &
                   oy(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),     &
                   oz(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2)
  REAL_euwp, INTENT(IN) ::  &
                    h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
  REAL_euwp, INTENT(OUT) ::  &
                    pfx(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),     &
                    pfy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
  REAL(KIND=euwp),INTENT(IN) :: dxi,dyi,dzi,dt
  REAL(KIND=euwp) :: gc1,gc2,gc3,dinv 
  INTEGER(KIND=iintegers) :: i,j,k


  gc1=dt*dxi
  gc2=dt*dyi
  gc3=dt*dzi

FullXYZDomainLoopDC(
          dinv=1._euwp/h(i,j,k);
          pfx(i,j,k)=  ((abs(ox(i,j,k,1)) +            &
                         abs(oy(i,j,k,1)) +            &
                         abs(oz(i,j,k,1)))*dinv);
          pfy(i,j,k)=  max(                                                  &
             .5_euwp*gc1/gc1*abs(ox(i+1,j  ,k  ,1)-ox(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc1*abs(ox(i  ,j+1,k  ,1)-ox(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc1*abs(ox(i  ,j  ,k+1,1)-ox(i  ,j  ,k-1,1))*dinv,  &
             .5_euwp*gc1/gc2*abs(oy(i+1,j  ,k  ,1)-oy(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc2*abs(oy(i  ,j+1,k  ,1)-oy(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc2*abs(oy(i  ,j  ,k+1,1)-oy(i  ,j  ,k-1,1))*dinv,  &
             .5_euwp*gc1/gc3*abs(oz(i+1,j  ,k  ,1)-oz(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc3*abs(oz(i  ,j+1,k  ,1)-oz(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc3*abs(oz(i  ,j  ,k+1,1)-oz(i  ,j  ,k-1,1))*dinv );
)
END SUBROUTINE compute_courlipsch_Agrid_full

SUBROUTINE compute_courants(cr3d,crxy,crxz,cryz,crx,cry,crz,ox,oy,oz,h,   &
                            ntimelevel,dxi,dyi,dzi,dt,np,mp,lp,ih) 
  USE mpi_parallel, ONLY:globmax,mype
  USE mod_parameters, ONLY: ipoles,ibcx,ibcy,ibcz
  INTEGER(KIND=iintegers), INTENT(IN) :: np,mp,lp,ih,ntimelevel
  REAL_euwp, INTENT(IN) ::                                 &
                    h(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),         &
                   ox(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),     &
                   oy(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),     &
                   oz(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2)
  REAL(KIND=euwp),INTENT(IN) :: dxi,dyi,dzi,dt
  REAL(KIND=euwp),INTENT(OUT) :: cr3d,crxy,crxz,cryz,crx,cry,crz
  REAL(KIND=euwp) :: bufmax_in(7),bufmax_out(7)
  REAL(KIND=euwp) :: gc1,gc2,gc3
  REAL(KIND=euwp) :: oxa,oya,oza
  REAL(KIND=euwp) :: dinv 
  INTEGER(KIND=iintegers) :: i,j,k

  CALL ttbeg(402)
  CALL ttbeg(1402)

    cr3d=-HUGE(cr3d)
    crxy=-HUGE(cr3d)
    crxz=-HUGE(cr3d)
    cryz=-HUGE(cr3d)
    crx=-HUGE(cr3d)
    cry=-HUGE(cr3d)
    crz=-HUGE(cr3d)
    SELECT CASE (ntimelevel)
    CASE (1)
      gc1=dt*dxi
      gc2=dt*dyi
      gc3=dt*dzi
!   IF(mype.eq.0)  print *,'Lib gc1,gc2,gc3',gc1,gc2,gc3,dxi,dyi,dzi,dt
!   CALL flush(6)
FullXYZDomainLoopDC(
            oxa=abs(gc1*ox(i,j,k,0));
            oya=abs(gc2*oy(i,j,k,0));
            oza=abs(gc3*oz(i,j,k,0));
            cr3d=  max(cr3d,(oxa+oya+oza));
            crxy=  max(crxy,(oxa+oya));
            crxz=  max(crxz,(oxa+oza));
            cryz=  max(cryz,(oya+oza));
            crx=  max(crx,(oxa));
            cry=  max(cry,(oya));
            crz=  max(crz,(oza));
)
    CASE(2)
FullXYZDomainLoopDC(
            dinv=1._euwp/h(i,j,k);
            oxa=abs(ox(i,j,k,1)*dinv);
            oya=abs(oy(i,j,k,1)*dinv);
            oza=abs(oz(i,j,k,1)*dinv);
            cr3d=  max(cr3d,(oxa+oya+oza));
            crxy=  max(crxy,(oxa+oya));
            crxz=  max(crxz,(oxa+oza));
            cryz=  max(cryz,(oya+oza));
            crx=  max(crx,(oxa));
            cry=  max(cry,(oya));
            crz=  max(crz,(oza));
)
    CASE DEFAULT
            STOP 'Wrong ox level in computecourants'
    END SELECT

  CALL ttend(1402)
  bufmax_in(1)=cr3d
  bufmax_in(2)=crxy
  bufmax_in(3)=crxz
  bufmax_in(4)=cryz
  bufmax_in(5)=crx
  bufmax_in(6)=cry
  bufmax_in(7)=crz
  CALL globmax(bufmax_in,bufmax_out,7)
  cr3d=bufmax_out(1)
  crxy=bufmax_out(2)
  crxz=bufmax_out(3)
  cryz=bufmax_out(4)
  crx=bufmax_out(5)
  cry=bufmax_out(6)
  crz=bufmax_out(7)

  CALL ttend(402)

END SUBROUTINE compute_courants

SUBROUTINE print_courants(imode,cr3d,crxy,crxz,cryz,crx,cry,crz)
  REAL(KIND=euwp),INTENT(IN) :: cr3d,crxy,crxz,cryz,crx,cry,crz
  INTEGER(KIND=iintegers) :: imode 

    SELECT CASE (imode)
    CASE(1)
           PRINT *, "----------------------------------------------------"
           PRINT *, "         Adaptive EULAG Courant numbers at n+1/2:"
    CASE(2)
            PRINT *, "----------------------------------------------------"
            PRINT *, "         Substep EULAG Courant numbers at n+1/2:"
    CASE(3)
            PRINT *, "----------------------------------------------------"
            PRINT *, "EULAG Courant numbers at n+1/2:"
    CASE(4)
            PRINT *, "----------------------------------------------------"
            PRINT *, "EULAG Courant numbers at n:"
    CASE DEFAULT
            STOP 'Wrong printing mode in in print_courants'
    END SELECT
           PRINT "(1x,'         cr3d,crxy,crz:',3(e23.16,2X))",cr3d,crxy,crz
           PRINT "(1x,'         cr3d,cryz,crx:',3(e23.16,2X))",cr3d,cryz,crx
           PRINT "(1x,'         cr3d,crxz,cry:',3(e23.16,2X))",cr3d,crxz,cry
           PRINT *, "----------------------------------------------------"

END SUBROUTINE print_courants 


#if(1==0)

!---------------------------------------------------------------------!
! Using vert. coordinate data for sea level (vcoord_in; COSMO),
! data profile (prof_in; COSMO) and local vert. coordinate data
! (prof_out; EULAG) linearly interpolate data profile to the local
! EULAG coordinates.
!---------------------------------------------------------------------!
SUBROUTINE acq_prof_lin(prof_in, vcoord_in, prof_out, vcoord_out)
!---------------------------------------------------------------------!
  
  
  ! Input data profile
  REAL(KIND=wp), DIMENSION(ke), INTENT(IN) :: prof_in
  
  ! Heights related to prof_in data  :: descending order
  REAL(KIND=wp), DIMENSION(ke), INTENT(IN) :: vcoord_in
  
  ! Heights related to prof_out data :: ascending order
  REAL(KIND=wp), DIMENSION(l), INTENT(IN) :: vcoord_out
  
  ! Subroutine output
  REAL(KIND=wp), DIMENSION(l), INTENT(OUT) :: prof_out
  
  REAL(KIND=wp) :: coe2
  INTEGER :: k, k1
  
  k1 = 1
  
  ! From surface to the stratoshere
  DO k=2,l
    DO WHILE (vcoord_in(ke-k1+1) < vcoord_out(k))
      k1 = k1 + 1
      IF (k1 == ke+1) EXIT
    ENDDO
    ! vcoord_in(ke-k1+1) >=vcoord_out(k) >= vcoord_in(ke-k1+2) 
    IF (k1 == 1) THEN
      prof_out(k) = prof_in(ke-k1+1)
    ELSE IF (k1 == ke+1) THEN
      prof_out(k) = prof_in(ke-k1+2)
    ELSE
      coe2 = (vcoord_out(k)-vcoord_in(ke-k1+2)) &
           / (vcoord_in(ke-k1+1)-vcoord_in(ke-k1+2))
      prof_out(k)  =               coe2  * prof_in(ke-k1+1)  &
                   +  (1._wp - coe2) * prof_in(ke-k1+2)
    ENDIF
  ENDDO
  
  
  prof_out(1) = prof_out(2)
  
  
END SUBROUTINE acq_prof_lin


!---------------------------------------------------------------------!
!
!---------------------------------------------------------------------!
SUBROUTINE printsummary(buffer,bufnames,maxvar,subname)
!---------------------------------------------------------------------!
  USE data_msg
  USE src_parallel
  INTEGER(KIND=iintegers),PARAMETER :: idmaxvars=40
  REAL(KIND=euwp) buffer(12,idmaxvars)
  INTEGER(KIND=iintegers),PARAMETER :: iflo=6
  CHARACTER(5) bufnames(idmaxvars)
  CHARACTER(len=*) subname
  INTEGER(KIND=iintegers) :: k,maxvar
  IF (mype == 0) THEN
    WRITE(iflo,*) '....................................'
    WRITE(iflo,*) 'Debugging subroutine ',subname
    WRITE(iflo,*) '....................................'
    
    WRITE(iflo,*) 'Variable name  ave  max at  ia,ja,ka  min at ia,ja,ka'
    DO k=1,maxvar
      WRITE(iflo,'(A5," ave",E12.4," max",E12.4," at",3I4," min",E12.4," at",3I4)') &
      bufnames(k),buffer(1,k),buffer(2,k), &
      NINT(buffer(4,k)),NINT(buffer(5,k)),NINT(buffer(6,k)), &
      buffer(3,k),NINT(buffer(7,k)),NINT(buffer(8,k)), &
      NINT(buffer(9,k))
    ENDDO
    WRITE(iflo,*) '....................................'
    WRITE(iflo,*) 'End debugging subroutine ',subname
    WRITE(iflo,*) '....................................'
  ENDIF
END SUBROUTINE   PRINTsummary


SUBROUTINE interpolate_moist_to_eulag(qv_eu,qc_eu, qr_eu, qs_eu, qi_eu, qg_eu,&
         fqv_eu, fqc_eu, fqr_eu, fqi_eu, ntime)

  USE src_tracer,        ONLY: trcr_get, trcr_errorstr
  USE data_modelconfig,  ONLY: idt_qv, idt_qc, idt_qr, idt_qi, idt_qs, idt_qg
  USE environment,       ONLY: model_abort
  USE data_param,        ONLY: imoist
  USE data_runcontrol,   ONLY: itype_gscp


  REAL (KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(OUT) :: &
    qv_eu,  & 
    qc_eu,  &
    qr_eu,  &
    qs_eu,  &
    qi_eu,  & 
    qg_eu,  &
    fqv_eu, &
    fqc_eu, &
    fqr_eu, &
    fqi_eu

  INTEGER   (KIND=iintegers)  :: ntime


  REAL (KIND = wp), POINTER     :: &                     ! TODO : check if OK
    qv_c(:,:,:)=> NULL(),       &
    qc_c(:,:,:)=> NULL(),       &
    qr_c(:,:,:)=> NULL(),       &
    qs_c(:,:,:)=> NULL(),       &
    qi_c(:,:,:)=> NULL(),       &
    qg_c(:,:,:)=> NULL()

  REAL (KIND = wp), POINTER     :: &
    qv_tend(:,:,:)=> NULL()  ,&
    qc_tend(:,:,:)=> NULL()  ,&
    qr_tend(:,:,:)=> NULL()  ,&
    qi_tend(:,:,:)=> NULL()

  REAL (KIND = wp),DIMENSION(ie,je,ke) :: &
    mr_mult

  CHARACTER (LEN=80) :: yzerrmsg
  CHARACTER (LEN=25) :: yzroutine
  INTEGER   (KIND=iintegers)  :: izerror

  IF(imoist == 1) THEN
  CALL trcr_get(izerror, idt_qv, ptr_tlev = ntime, ptr = qv_c)
  CALL trcr_get(izerror, idt_qc, ptr_tlev = ntime, ptr = qc_c)
  CALL trcr_get(izerror, idt_qr, ptr_tlev = ntime, ptr = qr_c)
  CALL trcr_get(izerror, idt_qs, ptr_tlev = ntime, ptr = qs_c)
    IF(itype_gscp > 2) THEN
      CALL trcr_get(izerror, idt_qi, ptr_tlev = ntime, ptr = qi_c)
      mr_mult(:,:,:) = 1._euwp &
                    / (1._euwp - (qv_c + qc_c + qr_c + qs_c +  qi_c))
    ELSE
      mr_mult(:,:,:) = 1._euwp &
                    / (1._euwp - (qv_c + qc_c + qr_c + qs_c))
    ENDIF
    IF(itype_gscp > 3) THEN
      CALL trcr_get(izerror, idt_qg, ptr_tlev = ntime, ptr = qg_c)
      mr_mult(:,:,:) = 1._euwp &
                    / (1._euwp - (qv_c + qc_c + qr_c + qs_c +  qi_c + qg_c))
    ENDIF
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  CALL trcr_get(izerror, 'QV', ptr_tens = qv_tend)
  CALL trcr_get(izerror, 'QC', ptr_tens = qc_tend)
  CALL trcr_get(izerror, 'QR', ptr_tens = qr_tend)
    IF(itype_gscp > 2) THEN
      CALL trcr_get(izerror, 'QI', ptr_tens = qi_tend)
    ENDIF
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  qv_c = qv_c*mr_mult
  qc_c = qc_c*mr_mult
  qr_c = qr_c*mr_mult
    IF(itype_gscp > 2) THEN
      qs_c = qs_c*mr_mult
      qi_c = qi_c*mr_mult
      IF(itype_gscp > 3) THEN
        qg_c = qg_c*mr_mult
      ENDIF
    ENDIF

  qv_tend = qv_tend*mr_mult
  qc_tend = qc_tend*mr_mult
  qr_tend = qr_tend*mr_mult
    IF(itype_gscp > 2) THEN
      qi_tend = qi_tend*mr_mult
    ENDIF
  ENDIF !imoist 
  ! Get mixing rations to EULAG variables : 
  CALL interpolate_scalar_to_eulag(qv_eu(:,:,:),  qv_c(:,:,:))
  CALL interpolate_scalar_to_eulag(qc_eu(:,:,:),  qc_c(:,:,:))
  CALL interpolate_scalar_to_eulag(qr_eu(:,:,:),  qr_c(:,:,:))
  IF(itype_gscp > 2) THEN
    CALL interpolate_scalar_to_eulag(qs_eu(:,:,:),  qs_c(:,:,:))
    CALL interpolate_scalar_to_eulag(qi_eu(:,:,:),  qi_c(:,:,:))
    IF(itype_gscp > 3) THEN
      CALL interpolate_scalar_to_eulag(qg_eu(:,:,:),  qg_c(:,:,:))
    ENDIF
  ENDIF

  CALL interpolate_scalar_to_eulag(fqv_eu(:,:,:),  qv_tend(:,:,:))
  CALL interpolate_scalar_to_eulag(fqc_eu(:,:,:),  qc_tend(:,:,:))
  CALL interpolate_scalar_to_eulag(fqr_eu(:,:,:),  qr_tend(:,:,:))
  IF(itype_gscp > 2) THEN
    CALL interpolate_scalar_to_eulag(fqi_eu(:,:,:),  qi_tend(:,:,:))
  ENDIF

! qv_eu(:,:,:) = 0._wp 
! qc_eu(:,:,:) = 0._wp 
! qr_eu(:,:,:) = 0._wp 
! qs_eu(:,:,:) = 0._wp 
! qi_eu(:,:,:) = 0._wp 
! qg_eu(:,:,:) = 0._wp 
! fqv_eu(:,:,:) = 0._wp 
! fqc_eu(:,:,:) = 0._wp 
! fqr_eu(:,:,:) = 0._wp 
! fqi_eu(:,:,:) = 0._wp 
  
END SUBROUTINE interpolate_moist_to_eulag



SUBROUTINE interpolate_moist_bc_to_eulag(qv_eu, qc_eu, qr_eu, qs_eu, &
  qi_eu, qg_eu, ntime)

  USE src_tracer,        ONLY : trcr_get, trcr_errorstr
  USE data_modelconfig,  ONLY:  idt_qv, idt_qc, idt_qr, idt_qi, idt_qs, idt_qg
  USE environment,       ONLY : model_abort
  USE data_runcontrol,   ONLY: itype_gscp


  REAL (KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(OUT) :: &
    qv_eu,  & 
    qc_eu,  &
    qr_eu,  &
    qs_eu,  &
    qi_eu,  & 
    qg_eu

  INTEGER   (KIND=iintegers)  :: ntime

  REAL (KIND = wp), POINTER     :: &                     ! TODO : check if OK
    qv_c(:,:,:)=> NULL(),       &
    qc_c(:,:,:)=> NULL(),       &
    qr_c(:,:,:)=> NULL(),       &
    qs_c(:,:,:)=> NULL(),       &
    qi_c(:,:,:)=> NULL(),       &
    qg_c(:,:,:)=> NULL()

  REAL (KIND = wp), POINTER     :: &
    qv_tend(:,:,:)=> NULL()  ,&
    qc_tend(:,:,:)=> NULL()  ,&
    qr_tend(:,:,:)=> NULL()  ,&
    qi_tend(:,:,:)=> NULL()

  REAL (KIND = wp),DIMENSION(ie,je,ke) :: &
    mr_mult

  CHARACTER (LEN=80) :: yzerrmsg
  CHARACTER (LEN=25) :: yzroutine
  INTEGER   (KIND=iintegers)  :: izerror

  CALL trcr_get(izerror, idt_qv, ptr_tlev = ntime, ptr = qv_c)
  CALL trcr_get(izerror, idt_qc, ptr_tlev = ntime, ptr = qc_c)
  CALL trcr_get(izerror, idt_qr, ptr_tlev = ntime, ptr = qr_c)
  CALL trcr_get(izerror, idt_qs, ptr_tlev = ntime, ptr = qs_c)
  IF(itype_gscp > 2) THEN
    CALL trcr_get(izerror, idt_qi, ptr_tlev = ntime, ptr = qi_c)
  mr_mult(:,:,:) = 1._euwp &
                 / (1._euwp - (qv_c + qc_c + qr_c + qs_c +  qi_c ))
  ELSE
  mr_mult(:,:,:) = 1._euwp &
                 / (1._euwp - (qv_c + qc_c + qr_c + qs_c ))
  ENDIF

  IF(itype_gscp > 3) THEN
      CALL trcr_get(izerror, idt_qg, ptr_tlev = ntime, ptr = qg_c)
  mr_mult(:,:,:) = 1._euwp &
                 / (1._euwp - (qv_c + qc_c + qr_c + qs_c +  qi_c + qg_c))
  ENDIF

  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF


  qv_c = qv_c*mr_mult
  qc_c = qc_c*mr_mult
  qr_c = qr_c*mr_mult
  qs_c = qs_c*mr_mult
  IF(itype_gscp > 2) THEN
    qi_c = qi_c*mr_mult
    IF(itype_gscp > 3) THEN
    qg_c = qg_c*mr_mult
    ENDIF
  ENDIF

  ! Get mixing rations to EULAG variables : 
  CALL interpolate_scalar_to_eulag(qv_eu(:,:,:),  qv_c(:,:,:))
  CALL interpolate_scalar_to_eulag(qc_eu(:,:,:),  qc_c(:,:,:))
  CALL interpolate_scalar_to_eulag(qr_eu(:,:,:),  qr_c(:,:,:))
  CALL interpolate_scalar_to_eulag(qs_eu(:,:,:),  qs_c(:,:,:))
  IF(itype_gscp > 2) THEN
    CALL interpolate_scalar_to_eulag(qi_eu(:,:,:),  qi_c(:,:,:))
    IF(itype_gscp > 3) THEN
      CALL interpolate_scalar_to_eulag(qg_eu(:,:,:),  qg_c(:,:,:))
    ENDIF
  ENDIF

!  qv_eu(:,:,:) = 0._wp 
!  qc_eu(:,:,:) = 0._wp 
!  qr_eu(:,:,:) = 0._wp 
!  qs_eu(:,:,:) = 0._wp 
!  qi_eu(:,:,:) = 0._wp 
!  qg_eu(:,:,:) = 0._wp 
  
END SUBROUTINE interpolate_moist_bc_to_eulag


!--------------------------------------------------------------------!
! 1D history package
!--------------------------------------------------------------------!
SUBROUTINE store_histo1d(hov_time,hov_z,hov_x_size,  & 
                         hov_y_size,hov_z_size,hov_halo_size,  &
                         hov_data,hov_varcnt,     &
                         pcosmo_nnew,peulag_term,peulag_dyn,  &
                         excosmo_nnew,exeulag_term,exeulag_dyn,  &
                         rho_cosmo,rho_eulag,t_cosmo,t_eulag,    &
                         w_cosmo,w_eulag,qv_cosmo,qv_eulag,th_cosmo,th_eulag,tht_eulag,fp)
  USE data_runcontrol,   ONLY:  ntstep,nstart,nstop
  USE data_param, ONLY: ipresdiag_int
  USE data_modelconfig, ONLY: dt 
  INTEGER(KIND=iintegers):: hov_x_size,hov_y_size,hov_z_size,hov_halo_size, &
                            hov_varcnt,hov_datacnt,icnt,nsnap,interv
  CHARACTER :: hov_name(13,hov_varcnt)
  CHARACTER,DIMENSION(13)  :: varstring
  REAL(KIND=euwp), DIMENSION(:,:,:) :: hov_data
  REAL(KIND=euwp), DIMENSION(:) :: hov_time,hov_z
  REAL(KIND=euwp), DIMENSION(hov_z_size)  :: aveprof
  REAL(KIND=euwp), DIMENSION(:,:,:) :: &
                         pcosmo_nnew,peulag_term,peulag_dyn,  &
                         excosmo_nnew,exeulag_term,exeulag_dyn,  &
                         rho_cosmo,rho_eulag,t_cosmo,t_eulag,    &
                          w_cosmo,w_eulag,qv_cosmo,qv_eulag,     &
                         th_cosmo,th_eulag,tht_eulag,fp 

  REAL(KIND=euwp), DIMENSION(1-hov_halo_size:hov_x_size+hov_halo_size,  &
                             1-hov_halo_size:hov_y_size+hov_halo_size,  &
                             1-hov_halo_size:hov_z_size+hov_halo_size)  ::  &
                             diftemp
  
   interv =ipresdiag_int
   IF((nstop-nstart+1)/interv > nhist) THEN
     STOP 'nhist is too low, as (nstop-nstart+1)/interv > nhist.'// &
           'Halting  now to prevent segmentation fault later'
   ENDIF
  IF(ntstep == nstart) nsnap = 0
  IF(modulo(ntstep,interv)==0) THEN
    nsnap=ntstep/interv+1
    hov_time(nsnap)=ntstep*dt/3600.
    icnt=0_iintegers 
!Diagnose dynamical, thermodynamical and COSMO large scale pressure
    icnt=icnt+1_iintegers
    CALL ave_fld(pcosmo_nnew,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('cosmo_p__nnew',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(peulag_term,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('eulag_p__term',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(peulag_dyn,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('eulag_p__dynm',varstring)

    icnt=icnt+1_iintegers
    diftemp=peulag_dyn-pcosmo_nnew 
    CALL ave_fld(diftemp,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('dif_p__dyn_co',varstring)

    icnt=icnt+1_iintegers
    diftemp=peulag_dyn-peulag_term 
    CALL ave_fld(diftemp,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('dif_p_dyn_ter',varstring)

!Diagnose dynamical, thermodynamical and COSMO large scale Exner 
    icnt=icnt+1_iintegers
    CALL ave_fld(excosmo_nnew,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('cosmo_ex_nnew',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(exeulag_term,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('eulag_ex_term',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(exeulag_dyn,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('eulag_ex_dynm',varstring)

    icnt=icnt+1_iintegers
    diftemp=exeulag_dyn-excosmo_nnew
    CALL ave_fld(diftemp,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('dif_ex_dyn_co',varstring)

    icnt=icnt+1_iintegers
    diftemp=exeulag_dyn-exeulag_term 
    CALL ave_fld(diftemp,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('dif_ex_dynter',varstring)

!Diagnose EULAG and COSMO density 
    icnt=icnt+1_iintegers
    CALL ave_fld(rho_cosmo,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('cosmo_rhonnew',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(rho_eulag,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('eulag_rhoterm',varstring)

    icnt=icnt+1_iintegers
    diftemp=rho_eulag-rho_cosmo 
    CALL ave_fld(diftemp,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('dif_rho_eu_co',varstring)


!Diagnose EULAG and COSMO temperature  
    icnt=icnt+1_iintegers
    CALL ave_fld(t_cosmo,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('cosmo_t__nnew',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(t_eulag,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('eulag_t__term',varstring)

    icnt=icnt+1_iintegers
    diftemp=t_eulag-t_cosmo 
    CALL ave_fld(diftemp,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('dif_t___eu_co',varstring)

!Diagnose EULAG and COSMO vertical velocity  
    icnt=icnt+1_iintegers
    CALL ave_fld(w_cosmo,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('cosmo_w__nnew',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(w_eulag,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('eulag_w__term',varstring)

    icnt=icnt+1_iintegers
    diftemp=w_eulag-w_cosmo
    CALL ave_fld(diftemp,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('dif_w___eu_co',varstring)

!Diagnose EULAG and COSMO water vapor  
    icnt=icnt+1_iintegers
    CALL ave_fld(qv_cosmo,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('cosmo_qv_nnew',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(qv_eulag,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('eulag_qv_term',varstring)

    icnt=icnt+1_iintegers
    diftemp=qv_eulag-qv_cosmo 
    CALL ave_fld(diftemp,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('dif_qv__eu_co',varstring)
!Diagnose EULAG and COSMO potential temperature  
    icnt=icnt+1_iintegers
    CALL ave_fld(th_cosmo,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('cosmo_th_nnew',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(th_eulag,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('eulag_th_term',varstring)

    icnt=icnt+1_iintegers
    diftemp=th_eulag-th_cosmo 
    CALL ave_fld(diftemp,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('dif_th__eu_co',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(fp,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('implicitprcor',varstring)

    icnt=icnt+1_iintegers
    CALL ave_fld(tht_eulag,aveprof)
    hov_data(1:hov_z_size,nsnap,icnt)=aveprof(1:hov_z_size)
    hov_name(:,icnt)=transfer('tht_eulag____',varstring)

      IF (ntstep == nstop.AND.mype == 0) THEN
        print *,'Varcnt,icnt',hov_varcnt,icnt
        print *,'Hov_name:',hov_name
        print *,'Number of snapshots nsnap',nsnap
        hov_datacnt=nsnap
        CALL output_histo1d_netcdf(hov_name,hov_time,zstr(1:lp),lp, &
                                         hov_data,nhisvars,hov_datacnt)
      ENDIF

  ENDIF ! ipresdiag_int
END SUBROUTINE store_histo1d

SUBROUTINE output_histo1d_netcdf(hov_name,hov_time,hov_z,hov_z_size, &
                                       hov_data,hov_varcnt,hov_datacnt)
  USE netcdf
  INTEGER(KIND=iintegers):: hov_z_size,hov_varcnt,hov_datacnt,i,ier
  INTEGER(KIND=iintegers):: filehandle,timevarhandle,zvarhandle
  INTEGER(KIND=iintegers):: varhandle(hov_varcnt)
  CHARACTER :: hov_name(13,hov_varcnt)
  REAL(KIND=euwp), DIMENSION(:,:,:) :: hov_data
  REAL(KIND=euwp), DIMENSION(:) :: hov_time
  REAL(KIND=euwp), DIMENSION(:)  :: hov_z 
  CALL net_init_histo1d_file('pressdiag.nc','time',filehandle,hov_datacnt,hov_z_size)

  DO i=1,hov_varcnt
    CALL net_init_histo1d_var(filehandle,hov_name(:,i),varhandle(i))
  END DO
  ier=nf90_enddef(filehandle)
  CALL net_out_histo1dtime(filehandle,hov_time,hov_datacnt)
  CALL net_out_histo1dz(filehandle,hov_z   ,hov_z_size)
  DO i=1,hov_varcnt
    CALL net_out_histo1d(filehandle,hov_data(:,:,i),varhandle(i),hov_datacnt,hov_z_size)
  END DO
  ier =  nf90_close(filehandle)
END SUBROUTINE output_histo1d_netcdf




SUBROUTINE net_init_histo1d_file(file_name,dim_name,ihfilehandle,timevarsize,zvarsize)
!--------------------------------------------------------------------!
  USE netcdf
  CHARACTER*4 dim_name
  CHARACTER(LEN=*) file_name
  INTEGER(KIND=iintegers) idim,ihfilehandle,itimevarid,ier,timevarsize,zvarsize
  INTEGER(KIND=iintegers) adim(1)
  
  ier = nf90_create(file_name,NF90_WRITE,ihfilehandle)
  IF (ier /= 0) PRINT *,'init_hist1di nf_create ERROR # ',nf90_strerror(ier)
  ier = nf90_def_dim(ihfilehandle,dim_name,timevarsize,idim)
  IF (ier /= 0) PRINT *,'init_hist1di def dimen ERROR # ',nf90_strerror(ier)
  adim(1)=idim
  ier = nf90_def_var(ihfilehandle,dim_name,NF90_DOUBLE,adim,itimevarid)
  IF (ier /= 0) PRINT *,'init_hist1di def var time ERROR # ',nf90_strerror(ier)
  ier = nf90_def_dim(ihfilehandle,'z',lp,idim)
  IF (ier /= 0) PRINT *,'init_hist1di def dimen ERROR # ',nf90_strerror(ier)
  adim(1)=idim
  ier = nf90_def_var(ihfilehandle,'z',NF90_DOUBLE,adim,itimevarid)
  IF (ier /= 0) PRINT *,'init_hist1di def var z ERROR # ',nf90_strerror(ier)
  !      ier = nf_def_dim(ihfilehandle,'sample',1,idim)
  !      if (ier.ne.NF_NOERR) print *,'init_hist1di def dimen ERROR # ',ier
END SUBROUTINE net_init_histo1d_file



!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE net_init_histo1d_var(ihfilehandle,var_name,varhandle)
!--------------------------------------------------------------------!
  USE netcdf
  CHARACTER :: var_name(13)
  CHARACTER(LEN=13) :: varstring
  INTEGER(KIND=iintegers) ihfilehandle,idim,idimv(2),ier,varhandle
  ier = nf90_inq_dimid(ihfilehandle,'time',idim)
  IF (ier /= 0) PRINT *,'init_hist1d defvudiminqERROR # ',nf90_strerror(ier)
  idimv(1)=idim
  ier = nf90_inq_dimid(ihfilehandle,'z',idim)
  IF (ier /= 0) PRINT *,'init_hist1d defv diminqERROR # ',nf90_strerror(ier)
  idimv(2)=idim
  !      ier = nf_inq_dimid(ihfilehandle,'sample',idim)
  !      if (ier.ne.nf_noerr) print *,'init_hist1dv def diminqERROR # ',ier
  !      idimv(2)=idim
  print *,'Creating pressdiag.nc variable: ',var_name
  ier=nf90_def_var(ihfilehandle,transfer(var_name,varstring),NF90_DOUBLE,idimv,varhandle)
  IF (ier /= 0) PRINT *,'init_hist1dv def var   ERROR # ',nf90_strerror(ier)
END SUBROUTINE net_init_histo1d_var


!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE net_out_histo1d(filehandle,var,varhandle,nsnap,zvarsize)
!--------------------------------------------------------------------!
  USE netcdf
  USE data_runcontrol,   ONLY: nstart, ntstep,nstop
  INTEGER(KIND=iintegers) filehandle,varhandle,istartv(2),icountv(2)
  INTEGER(KIND=iintegers) ier,nsnap,i,k,zvarsize 
  REAL(KIND=euwp) var(zvarsize,nhist)
  REAL(KIND=euwp) buffer(nsnap,zvarsize)
  DO i=1,nsnap
    DO k=1,zvarsize
      buffer(i,k)=var(k,i)
    ENDDO
  ENDDO
!  print *,'Fav_2:',buffer(1,1:lp)
  
  istartv(1)=1_iintegers
  istartv(2)=1_iintegers
  icountv(1)=nsnap
  icountv(2)=zvarsize
  PRINT *,'Writing dimension of size',nsnap,'times',zvarsize
  ier = nf90_put_var(filehandle,varhandle,buffer,istartv,icountv) 
    IF (ier /= 0) &
             print *,'out_history put  ERROR # ',nf90_strerror(ier)
END SUBROUTINE net_out_histo1d

!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE net_out_histo1dtime(filehandle,var,nsnap)
!--------------------------------------------------------------------!
  USE netcdf
  USE data_runcontrol,   ONLY: nstart, ntstep,nstop
  REAL(KIND=euwp) var(nhist)
  REAL(KIND=euwp) buffer(nsnap)
  INTEGER(KIND=iintegers) filehandle,varhandle,istartv(1),icountv(1)
  INTEGER(KIND=iintegers) ier,nsnap,i 
  ier = nf90_inq_varid(filehandle,'time',varhandle)
    IF (ier /= 0) &
             print *,'out_history inq  ERROR # ',nf90_strerror(ier)
  DO i=1,nsnap
    buffer(i)=var(i)
  ENDDO
  
  istartv(1)=1_iintegers
  icountv(1)=nsnap
  PRINT *,'Writing dimension of size',nsnap
  ier = nf90_put_var(filehandle,varhandle,buffer,istartv,icountv) 
    IF (ier /= 0) &
             print *,'out_history time put  ERROR # ',nf90_strerror(ier)
END SUBROUTINE net_out_histo1dtime

!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE net_out_histo1dz(filehandle,var,zvarsize)
!--------------------------------------------------------------------!
  USE netcdf
  USE data_runcontrol,   ONLY: nstart, ntstep,nstop
  INTEGER(KIND=iintegers) filehandle,varhandle,istartv(1),icountv(1)
  INTEGER(KIND=iintegers) ier,k,zvarsize 
  REAL(KIND=euwp) var(zvarsize)
  REAL(KIND=euwp) buffer(zvarsize)
  ier = nf90_inq_varid(filehandle,'z',varhandle)
    IF (ier /= 0) &
             print *,'out_history inq  ERROR # ',nf90_strerror(ier)
    DO k=1,zvarsize
      buffer(k)=var(k)
    ENDDO
  
  istartv(1)=1_iintegers
  icountv(1)=zvarsize
  PRINT *,'Writing dimension of size',zvarsize
  ier = nf90_put_var(filehandle,varhandle,buffer,istartv,icountv) 
    IF (ier /= 0) &
             print *,'out_history put  ERROR # ',nf90_strerror(ier)
END SUBROUTINE net_out_histo1dz
!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnet_init_histo2d_file(file_name,dim_name,ihfilehandle,idimvarid,timevarsize)
!--------------------------------------------------------------------!
  USE mpi
  USE pnetcdf
  USE data_parallel,   ONLY: icomm_cart
  character*13 var_name
  character*4 dim_name
  character*13 file_name
  INTEGER(KIND=iintegers) idim,ihfilehandle,idimvarid,ier,timevarsize
  INTEGER(KIND=iintegers) adim(1)
  INTEGER(KIND=MPI_OFFSET_KIND) pn,pm,pt
  pn=n
  pm=m
  pt=NF_UNLIMITED
  ier = nfmpi_create(icomm_cart,file_name,NF_WRITE+NF_64BIT_OFFSET,MPI_INFO_NULL,ihfilehandle)
  IF (ier /= 0 .AND. mype == 0) PRINT *,'init_hist2di nf_create ERROR # ',nfmpi_strerror(ier)
  ier = nfmpi_def_dim(ihfilehandle,dim_name,pt,idim)
  adim(1)=idim
  IF (ier /= 0 .AND. mype == 0) PRINT *,'init_hist2di def dimen ERROR # ',nfmpi_strerror(ier)
  ier = nfmpi_def_var(ihfilehandle,dim_name,NF_REAL,1,adim,idimvarid)
  IF (ier /= 0 .AND. mype == 0) PRINT *,'init_hist2di def vardimERROR # ',nfmpi_strerror(ier)
  ier = nfmpi_def_dim(ihfilehandle,'y',pm,idim)
  IF (ier /= 0 .AND. mype == 0) PRINT *,'init_hist2di def dimen ERROR # ',nfmpi_strerror(ier)
  !      ier = nf_def_dim(ihfilehandle,'sample',1,idim)
  ier = nfmpi_def_dim(ihfilehandle,'x',pn,idim)
  IF (ier /= 0 .AND. mype == 0) PRINT *,'init_hist2di def dimen ERROR # ',nfmpi_strerror(ier)
  !      if (ier.ne.NF_NOERR) print *,'init_hist1di def dimen ERROR # ',ier
END SUBROUTINE pnet_init_histo2d_file


!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnet_init_histo2d_var(ihfilehandle,var_name,idimvarid,varhandle)
!--------------------------------------------------------------------!
  USE mpi
  USE pnetcdf
  character*13 var_name
  INTEGER(KIND=iintegers) ihfilehandle,idim,idimv(3),idimvarid,ier,varhandle
  ier = nfmpi_inq_dimid(ihfilehandle,'time',idim)
  IF (ier /= 0 .AND. mype == 0) PRINT *,'init_hist2d defvudiminqERROR # ',nfmpi_strerror(ier)
  idimv(3)=idim
  ier = nfmpi_inq_dimid(ihfilehandle,'x',idim)
  IF (ier /= 0 .AND. mype == 0) PRINT *,'init_hist2d defv diminqERROR # ',nfmpi_strerror(ier)
  idimv(1)=idim
  ier = nfmpi_inq_dimid(ihfilehandle,'y',idim)
  IF (ier /= 0 .AND. mype == 0) PRINT *,'init_hist2d defv diminqERROR # ',nfmpi_strerror(ier)
  idimv(2)=idim
  !      ier = nf_inq_dimid(ihfilehandle,'sample',idim)
  !      if (ier.ne.nf_noerr) print *,'init_hist1dv def diminqERROR # ',ier
  !      idimv(2)=idim
  ier=nfmpi_def_var(ihfilehandle,var_name,NF_REAL,3,idimv,varhandle)
  IF (ier /= 0) PRINT *,'init_hist2dv def var   ERROR # ',nfmpi_strerror(ier)
END SUBROUTINE pnet_init_histo2d_var


!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnet_out_histo2d(filehandle,var,varhandle,nsnap)
!--------------------------------------------------------------------!
  USE mpi
  USE pnetcdf
  USE data_runcontrol,   ONLY: nstart, ntstep,nstop
  REAL(KIND=euwp) var(nhist,np,mp)
  REAL*4 buffer(np,mp,nsnap)
  INTEGER(KIND=iintegers) filehandle,varhandle,ier,nsnap,iter,i,j 
  INTEGER(KIND=MPI_OFFSET_KIND) istartv(3),icountv(3) 
  DO iter=1,nsnap
    DO i=1,np
      DO j=1,mp
        buffer(i,j,iter)=var(iter,i,j)
      ENDDO
    ENDDO
  ENDDO
  
  istartv(3)=1
  istartv(1)=nsubpos+1
  istartv(2)=msubpos+1
  icountv(3)=nsnap-1
  icountv(1)=np
  icountv(2)=mp
  !      icountv(2)=1
  ier = nfmpi_put_vara_real_all(filehandle,varhandle,  & 
      istartv,icountv,buffer) 
    IF (ier /= 0 .AND. mype == 0) &
             print *,'out_history2d put  ERROR # ',nfmpi_strerror(ier)
END SUBROUTINE pnet_out_histo2d



!--------------------------------------------------------------------!
!  Parallel Netcdf package for EULAG
!  For documentation on pnetcdf see http://trac.mcs.anl.gov/projects/parallel-netcdf  
!--------------------------------------------------------------------!
!  Subroutines:
!  pnet_create_common - Creates a parallel netcdf files and defines the dimensions 
!                       x,y,z,t
!  pnet_out_lng   - Writes data  to long restartable tape. Consists of full 
!                   model variable set and optional write of 
!                   second order restart data (u,v,w,ox,oy,oz additional levels). 
!                   Not intended to be changed often. 
!  pnet_out_std   - Writes data to standard tape - standard set of variables. 
!                   Intended to replace short tape for quick analysis
!                   of standard fields
!  pnet_out_chnk  - Subroutine writes a specific chunk 
!                   (vector, 2D plane, 3D variable) to arbitrary file
!  pnet_close     - Closing (or syncing) the pnetcdf files is essential, 
!                   files not closed properly are ureadable. 
!                   
!  pnettrans,pnettrans8,pnettransr - rewrites the matrix to a specified precision
!--------------------------------------------------------------------!
! IMPORTANT !!! Please note that Parallel Netcdf is sensitive to coding quality, 
!               especially to proper type of function arguments. Failure to use 
!               proper type of variables results in weird errors, use nfmpi_strerror 
!               for diagnostics. Long tape is double precision by default, other tapes 
!               assumed single precision for analysis. Precision change is possible
!               and easy, although precision switch is not implemented due to the 
!               coding difficulty yet. 
!--------------------------------------------------------------------!

!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnet_create_common(ipnind,nfhd,DID,tVID)
!--------------------------------------------------------------------!
  USE mpi
  USE pnetcdf 
  USE data_msg
  USE data_parallel,   ONLY: icomm_cart
  INTEGER(KIND=MPI_OFFSET_KIND) pn,pm,pl,pt
  INTEGER(KIND=iintegers) ier,iPRINT,ipnind,cmode,nfhd,tVID,DID(4)
  INTEGER(KIND=iintegers) nfp 
  !Cast domain size into OFFSET_KIND integers
  pn=n
  pm=m
  pl=l
  pt=NF_UNLIMITED
  nfp = NF_REAL
  cmode=NF_WRITE+NF_64BIT_OFFSET
  
  iPRINT=1 ! PRINT out all error messages including NF_NOERR
  iPRINT=0 ! PRINT out all real error messages
  
  !For compatibility with older pre Netcdf 3.6 (pre 2007) set cmode=NF_WRITE
  !This would limit the file size to 2 GB as permitted by CDF-1 format
  
  IF (ipnind == 1) THEN
    ier = nfmpi_create( icomm_cart,'tape.custom.nc',                &
        cmode, MPI_INFO_NULL, nfhd)
  ELSE IF(ipnind == 2) Then
    ier = nfmpi_create( icomm_cart,'tapes.nc',                      &
        cmode, MPI_INFO_NULL, nfhd)
  ELSE IF(ipnind == 3) THEN
    ier = nfmpi_create( icomm_cart,'tapef.nc',                      &
        cmode, MPI_INFO_NULL, nfhd)
  ELSE IF(ipnind == 4) THEN
    ier = nfmpi_create( icomm_cart,'tape.slice.nc',                 &
        cmode, MPI_INFO_NULL, nfhd)
  ENDIF
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           &print *,'PNCDF tape',ipnind,'creat stat ',nfmpi_strerror(ier)
  
  !Now defining dimensions of variables (x,y,z,t) or (x,z,t)
  ier=nfmpi_def_dim(nfhd,"x",pn,DID(1))
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF dim. x creat. stat ',nfmpi_strerror(ier)
  
  ier=nfmpi_def_dim(nfhd,"y",pm,DID(2))
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF dim. y creat. stat ',nfmpi_strerror(ier)
  
  ier=nfmpi_def_dim(nfhd,"z",pl,DID(3))
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF dim. z creat. stat ',nfmpi_strerror(ier)
  
  ier=nfmpi_def_dim(nfhd,"t",pt,DID(4))
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF dim. t creat. stat ',nfmpi_strerror(ier)
  
  ier=nfmpi_def_var(nfhd,"time",nfp,1,DID(4),tVID)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0)) &
           print *,'PNCDF var time creat. stat ',nfmpi_strerror(ier)
  
  ier = nfmpi_enddef(nfhd)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0)) &
           print *,'PNCDF enddef stat ',nfmpi_strerror(ier)
  
END SUBROUTINE pnet_create_common


!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnet_open_common(ipnind,nfhd,DID,tVID,iframe)
!--------------------------------------------------------------------!
  USE mpi
  USE pnetcdf 
  USE data_msg
  USE data_parallel,   ONLY: icomm_cart
  INTEGER(KIND=MPI_OFFSET_KIND) pn,pm,pl,pt
  INTEGER(KIND=iintegers) ier,iframe,iPRINT,ipnind,cmode,nfhd,tVID
  INTEGER(KIND=iintegers) DID(4)
  iPRINT=0
  cmode=NF_WRITE  !+NF_64BIT_OFFSET
  !For compatibility with older pre Netcdf 3.6 (pre 2007) set cmode=NF_WRITE
  !This would limit the file size to 2 GB as permitted by CDF-1 format
  ier=nf_noerr
  IF (IFrame > 0 .OR. iframe == -1) THEN
    IF (ipnind == 1) THEN
      ier = nfmpi_open( icomm_cart,'tape.custom.nc',                  &
          cmode, MPI_INFO_NULL, nfhd)
    ELSE IF(ipnind == 2) THEN
      ier = nfmpi_open( icomm_cart,'tapes.nc',                        &
          cmode, MPI_INFO_NULL, nfhd)
    ELSE IF(ipnind == 3) THEN
      ier = nfmpi_open( icomm_cart,'tapef.nc',                        &
          cmode, MPI_INFO_NULL, nfhd)
    ELSE IF(ipnind == 4) THEN
      ier = nfmpi_open( icomm_cart,'tape.slice.nc',                   &
          cmode, MPI_INFO_NULL, nfhd)
    ENDIF
  ENDIF
  
  !If file doesn''t exist, create it
  IF (ier /= nf_noerr  .OR.  IFrame == 0) THEN
    CALL pnet_create_common(ipnind,nfhd,DID,tVID) 
  ELSE
    !Now inquiring for dimension and dimension variables (x,y,z,t) handles
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))    &
             print *,'PNCDF tape',ipnind,'open stat ',nfmpi_strerror(ier)
    ier=nfmpi_inq_dimid(nfhd,"x",DID(1))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0)) &
             print *,'PNCDF dim. x inq. stat ',nfmpi_strerror(ier)
    ier=nfmpi_inq_dimid(nfhd,"y",DID(2))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0)) &
             print *,'PNCDF dim. y inq. stat ',nfmpi_strerror(ier)
    ier=nfmpi_inq_dimid(nfhd,"z",DID(3))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0)) &
             print *,'PNCDF dim. z inq. stat ',nfmpi_strerror(ier)
    ier=nfmpi_inq_dimid(nfhd,"t",DID(4))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0)) &
             print *,'PNCDF dim. t inq. stat ',nfmpi_strerror(ier)
    
    ier=nfmpi_inq_varid(nfhd,"time",tVID)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0)) &
             print *,'PNCDF time inq. stat ',nfmpi_strerror(ier)
  ENDIF
END SUBROUTINE pnet_open_common


!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnet_out_lng(iframe0, u,v,w,ox,oy,oz,rh,th,p,           &
                        fx,fy,fz,ft,qv,qc,qr,qs,qi,qg)     
!--------------------------------------------------------------------!
  USE mpi
  USE pnetcdf 
  USE data_msg
  USE data_eufields,     ONLY: zcr,the
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) ::  &
    &        th,rh,p,fx,fy,fz,ft,qv,qc,qr,qs,qi,qg
    REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2):: &
    &        u,v,w,ox,oy,oz
      
  INTEGER(KIND=iintegers) :: nfhd,ihnd,p3dim,p4dim,p5dim,iframe0,iframe,id
  INTEGER(KIND=iintegers) :: DLVID(40),tVID,spcDID
  
  integer(KIND=MPI_OFFSET_KIND) pspc,pframe,pone,ifrmask
  REAL(KIND=euwp) :: dane(np,mp,lp)        ! real*8
  
  integer(KIND=MPI_OFFSET_KIND) START3(3),COUNT3(3)
  integer(KIND=MPI_OFFSET_KIND) START4(4),COUNT4(4)
  INTEGER(KIND=iintegers) :: DID(4),DID5(5)
  INTEGER(KIND=iintegers) :: iPRINT,irsdta,nfp,ier,i
  !Switch for exact restart option
  !Default off, because at the moment not truly exact
  irsdta=0
  iframe=iframe0
  
  nfp = nf_double
  pone=1
  iPRINT=0
  
  
  CALL  pnet_open_common(3,nfhd,DID,tVID,iframe)
  !Automatic iframe setting
  IF (IFrame == 0) THEN
    ier=nfmpi_inq_unlimdim(nfhd,ihnd)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &   
             print *,'PNCDF chunk. inq unlim dim id ',nfmpi_strerror(ier)
    
    ier=nfmpi_inq_dimlen(nfhd,ihnd,ifrmask)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF chunk. inq dimlen ',nfmpi_strerror(ier)
    
    iframe=ifrmask+1
  ENDIF
  
  pframe=iframe
  
  START4(1)=nsubpos+1
  START4(2)=msubpos+1
  START4(3)=lsubpos+1
  START4(4)=iframe
  COUNT4(1)=np 
  COUNT4(2)=mp
  COUNT4(3)=lp
  COUNT4(4)=1
  p5dim=5    
  p4dim=4
  p3dim=3
  DO i=1,3
    START3(i)=START4(i)
    COUNT3(i)=COUNT4(i)
    DID5(i)=DID(i)
  ENDDO
  !Check if variables are already defined in the file
  ier=nfmpi_inq_varid(nfhd,'u',DLVID(1))
  !If not, we define them now
  IF (IFrame == 1 .AND. ier /= nf_noerr) THEN
    ier = nfmpi_redef(nfhd)
    id=1
    ier=nfmpi_def_var(nfhd,'u',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long u  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'v',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long v  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'w',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long w  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'ox',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long ox stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'oy',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long oy  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'oz',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long oz  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'rhc',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long rhc stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'th',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long th  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'p',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long  p stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'fx',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long fx  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'fy',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long fy  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'fz',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long fz  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'ft',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long ft  stat ',nfmpi_strerror(ier)
    
    id=id+1
    ier=nfmpi_def_var(nfhd,'qv',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qv  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'qc',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qc  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'qr',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qr  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'qs',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qs  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'qi',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qi  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'qg',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qg  stat ',nfmpi_strerror(ier)
    
    IF (IFrame == 1) THEN
      id=id+1
      ier=nfmpi_def_var(nfhd,"ELEVATION",nfp,p3dim,DID,DLVID(id)) 
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long ELEV stat ',nfmpi_strerror(ier)
    ENDIF
    IF (irsdta == 1) THEN
      id=id+1
      ier=nfmpi_def_var(nfhd,'u2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long u2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'v2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long v2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'w2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long w2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'ox2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long ox2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oy2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long oy2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oz2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long oz2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'u3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long u3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'v3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long v3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'w3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long w3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'ox3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long ox3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oy3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long oy3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oz3',nfp,p4dim,DID,DLVID(id))        
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long oz3  stat ',nfmpi_strerror(ier)
    ENDIF   !irsdta
    ier=nfmpi_enddef(nfhd)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF open_common enddef stat ',nfmpi_strerror(ier)
    
    
  ELSE   !iframe > 1 
    
    id=1
    ier=nfmpi_inq_varid(nfhd,'u',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long u  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'v',DLVID(id))                      
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long v  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'w',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long w  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'ox',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long ox inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'oy',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long oy  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'oz',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long oz  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'rhc',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long rhc inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'th',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long th  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'p',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long  p inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fx',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long fx  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fy',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long fy  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fz',DLVID(id)) 
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long fz  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'ft',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long ft  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qv',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qv  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qc',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qc  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qr',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qr  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qs',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qs  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qi',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qi  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qg',DLVID(id))
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
             print *,'PNCDF long qg  stat ',nfmpi_strerror(ier)
    
    IF (IFrame == 1) THEN
      id=id+1
      ier=nfmpi_inq_varid(nfhd,"ELEVATION",DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
               print *,'PNCDF long ELEV inq ',nfmpi_strerror(ier)
    ENDIF 
    IF (irsdta == 1) THEN
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'u2',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long u2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'v2',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long v2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'w2',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long w2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'ox2',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long ox2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oy2',DLVID(id)) 
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long oy2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oz2',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long oz2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'u3',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long u3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'v3',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long v3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'w3',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long w3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'ox3',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long ox3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oy3',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long oy3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oz3',DLVID(id))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))  &
               print *,'PNCDF long oz3  inq ',nfmpi_strerror(ier)
    ENDIF   ! irsdta
  ENDIF     !iframe
  !      ier=nfmpi_begin_indep_data(nfhd)
  !      ier = nfmpi_put_vara_real(nfhd,tVID,pframe,pone,time)
  !      IF((mype.eq.0.and.iprint.eq.1).or.(iprint.eq.0.and.ier.ne.0))
  !     &   print *,'PNCDF std var. time put stat ',nfmpi_strerror(ier)
  !      ier=nfmpi_end_indep_data(nfhd)
  
  id=1
  CALL pnettrans8(u(1-ih,1-ih,1-ih,0),dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long u wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(v(1-ih,1-ih,1-ih,0),dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long v wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(w(1-ih,1-ih,1-ih,0),dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long w wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(ox(1-ih,1-ih,1-ih,0),dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long ox wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(oy(1-ih,1-ih,1-ih,0),dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long oy wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(oz(1-ih,1-ih,1-ih,0),dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long oz wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(rh,dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long rhc wrt ',nfmpi_strerror(ier)
  id=id+1
  th(1:np,1:mp,1:lp)=th(1:np,1:mp,1:lp)+the(1:np,1:mp,1:lp)
  CALL pnettrans8(th,dane)
  th(1:np,1:mp,1:lp)=th(1:np,1:mp,1:lp)-the(1:np,1:mp,1:lp)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long th wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(p,dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long p wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(fx,dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long fx wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(fy,dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long fy wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(fz,dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long fz wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(ft,dane)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF long ft wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(qv,dane)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
           print *,'PNCDF long qv wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(qc,dane)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
           print *,'PNCDF long qc wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(qr,dane)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
           print *,'PNCDF long qr wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(qs,dane)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
           print *,'PNCDF long qs wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(qi,dane)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
           print *,'PNCDF long qi wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(qg,dane)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
           print *,'PNCDF long qg wrt ',nfmpi_strerror(ier)
  IF (IFrame == 1) THEN
    id=id+1
    dane(1:np,1:mp,1:lp)=zcr(1:np,1:mp,1:lp)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START3,COUNT3,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long ELEV wrt ',nfmpi_strerror(ier)
  ENDIF
  
  IF (irsdta == 1) THEN 
    id=id+1
    CALL pnettrans8(u(1-ih,1-ih,1-ih,1),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long u2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(v(1-ih,1-ih,1-ih,1),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long v2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(w(1-ih,1-ih,1-ih,1),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long w2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(ox(1-ih,1-ih,1-ih,1),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long ox2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oy(1-ih,1-ih,1-ih,1),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long oy2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oz(1-ih,1-ih,1-ih,1),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long oz2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(u(1-ih,1-ih,1-ih,2),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long u3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(v(1-ih,1-ih,1-ih,2),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long v3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(w(1-ih,1-ih,1-ih,2),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long w3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(ox(1-ih,1-ih,1-ih,2),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long ox3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oy(1-ih,1-ih,1-ih,2),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long oy3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oz(1-ih,1-ih,1-ih,2),dane)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF long oz3 wrt ',nfmpi_strerror(ier)
  ENDIF   ! irsdta
  
  ier = nfmpi_close(nfhd)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'pnet out lng close',nfmpi_strerror(ier)
  
END SUBROUTINE pnet_out_lng


!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnet_out_chunk(var_name,file_name,imod,ia1,ia2,ia3,ifrmd,imem,var) 
!--------------------------------------------------------------------!
  USE mpi
  USE pnetcdf 
  USE data_parallel,   ONLY: icomm_cart
  !Subroutine writes a specific chunk (vector, 2D plane, 3D variable)
  !Type of chunk is specified by imod:
  ! 0 - full 3D variable at time iframe in SINGLE PRECISION
  ! 1 - full 3D variable at time iframe in DOUBLE PRECISION
  ! 2 - full 3D variable at time iframe in DP of size np+ia1, mp+ia2,lp+ia3 with one larger dimension
  ! 3 - full 3D variable at time iframe in DP of size np+ia1, mp+ia2,lp+ia3 with two larger dimensions
  ! 4 - full 3D variable at time iframe in DP of size np+ia1, mp+ia2,lp+ia3 with three larger dimensions
  ! 6 - halo of  3D variable at time iframe in DP in LR direction
  ! 7 - halo of  3D variable at time iframe in DP in BT direction
  ! 8 - halo of  3D variable at time iframe in DP in GS direction
  ! 9 - halo of  3D variable at time iframe - 8 corners
  ! 10 - 2D xy plane at level ia3
  ! 11 - 2D xz plane at index ia2
  ! 12 - 2D yz plane at index ia1
  ! 13 - 2D xy plane at level ia3 - no halo in z
  ! 20 - 1D vector along z at point ia1,ia2
  ! 21 - 1D vector along y at point ia1,ia3
  ! 22 - 1D vector along x at point ia2,ia3
  ! 23 - 1D vector along z at point ia2,ia3 - no halo in x,y
  ! 30 - 0D point  
  ! 34 - full 3D variable without time dimension
  ! 35 - 2D xy plane without time dimension
  ! Write of 2D variables is realized in mode 10
  ! KEEP imem=0 EXCEPT WHEN IMOD=2-4 imem=1
  
  REAL(KIND=euwp) var(1-ih:np+ia1*imem+ih,                          &
           1-ih:mp+ia2*imem+ih,                          &
           1-ih:lp+ia3*imem+ih)
  REAL(KIND=euwp) :: v1d(lp)
  INTEGER(KIND=iintegers) :: p4dim,p3dim,p2dim,p1dim,iframe,ifrmd
  INTEGER(KIND=iintegers) :: DCHVID
  integer(KIND=MPI_OFFSET_KIND) pn,pm,pl,pt,ifrmask
  character(9) var_name
  character(8) file_name
  real*4 dane3d(np,mp,lp)
  real*8 dane83d(np,mp,lp)
  real*8 danea83d((np+ia1-1)*imem+1,                                &
  &            (mp+ia2-1)*imem+1,(lp+ia3-1)*imem+1)
  real*8 danex83d((np+ia1-1)*imem+1,mp,lp)
  real*8 danez83d(np,mp,(lp+ia3-1)*imem+1)
  real*8 halo8gs(np,mp,ih)
  real*8 halo8lr(ih,mp,lp)
  real*8 halo8bt(np,ih,lp)
  real*8 halo8cr(ih,ih,ih)
  REAL(KIND=euwp) dane2dxy(np,mp)
  real*4 dane2dxz(np,lp)
  real*4 dane2dyz(mp,lp)
  real*4 dane1dx(np)
  real*4 dane1dy(mp)
  real*4 dane1dz(lp)
  real*4 dane0d
  
  integer(KIND=MPI_OFFSET_KIND) START4(4),COUNT4(4)
  integer(KIND=MPI_OFFSET_KIND) START3(3),COUNT3(3)
  integer(KIND=MPI_OFFSET_KIND) START2(2),COUNT2(2)
  INTEGER(KIND=iintegers) :: DID(4)
  INTEGER(KIND=iintegers) :: DID3(3)
  INTEGER(KIND=iintegers) :: DID2(2)
  INTEGER(KIND=iintegers) :: DID1(1),iPRINT 
  INTEGER(KIND=iintegers) :: irsdta,imod,nfp,nfpd,nfhd,ihnd,i,j,k,ier
  INTEGER(KIND=iintegers) :: iedge,icorn,i1p,i2p,i3p,i1,i2,i3,ii,jj,kk
  INTEGER(KIND=iintegers) :: ia1,ia2,ia3,imem
  
  !      IF(imem.ne.0.and.imod.ne.2) 
  !     & STOP 'PROVIDE imem=0 in pnet_out_chunk UNLESS imod=2!!!!'  
  !      IF(imod.eq.2.and.imem.ne.1)
  !     & STOP 'PROVIDE imem=1 in pnet_out_chunk for imod=2!!!!'  
  p4dim=4
  p3dim=3
  p2dim=2 
  p1dim=1
  pn=n
  pm=m
  pl=l
  IF (imod >= 2 .AND. imod <= 4) THEN
    pn=n+ia1
    pm=m+ia2
    pl=l+ia3
  ELSE IF(imod == 6) THEN
    pn=2*ih
    pm=m
    pl=l
  ELSE IF(imod == 7) THEN
    pn=n
    pm=2*ih
    pl=l
  ELSE IF(imod == 8) THEN
    pn=n
    pm=m
    pl=2*ih
  ELSE IF(imod == 9) THEN
    pn=2*ih
    pm=2*ih
    pl=2*ih
  ENDIF
  pt=NF_UNLIMITED
  nfp = nf_real
  nfpd= nf_double
  iPRINT=0
  
  !-------------------------
  IF (ifrmd == 0) THEN         !if #1 create names and dims
  !-------------------------
    
    !==================================
    !Now create pnetcdf file if not exist
    ! nfmpi_open will generate MPI error
    ! perhaps implementing CALL system(test -e file_name) is a solution 
    !==================================
    
    ier = nfmpi_open( icomm_cart,file_name,                      &
        NF_WRITE, MPI_INFO_NULL, nfhd)
    IF (ier /= nf_noerr) THEN  !IF #2
      ier = nfmpi_create( icomm_cart,file_name,                  &
      !     &                    NF_WRITE, MPI_INFO_NULL, nfhd)               
          NF_WRITE+NF_64BIT_OFFSET, MPI_INFO_NULL, nfhd)
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
               print *,'PNCDF chunk '//file_name//' creat. stat ',        &
               nfmpi_strerror(ier)
    ENDIF   !create file on nf_noerr ENDIF #2
    
    ier = nfmpi_redef(nfhd)
    
    !==================================
    !Now inquire or define dimensions in the pnetcdf file
    !==================================
    
    IF (imod /= 12 .AND. imod /= 20 .AND. imod /= 21 .AND.  imod /= 30 .AND. imod /= 23 .AND. imod /= 33) THEN !IF #3
      ier=nfmpi_inq_dimid(nfhd,"x",DID(1))
      IF (ier /= nf_noerr) THEN   !IF #4
        ier=nfmpi_def_dim(nfhd,"x",pn,DID(1))
        IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))&
                 print *,'PNCDF dim. x creat. stat ',nfmpi_strerror(ier)
      ENDIF ! x dim create !ENDIF #4
    ENDIF  ! x dim inquire or create !ENDIF #3 
    
    IF (imod /= 11 .AND. imod /= 20 .AND. imod /= 22 .AND.  imod /= 30 .AND. imod /= 23 .AND. imod /= 33) THEN !IF #5
      ier=nfmpi_inq_dimid(nfhd,"y",DID(2))
      IF (ier /= nf_noerr) THEN  ! IF #6 
        ier=nfmpi_def_dim(nfhd,"y",pm,DID(2))
        IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))&
                 print *,'PNCDF dim. y creat. stat ',nfmpi_strerror(ier)
      ENDIF  ! y dim create  ENDIF #6 
    ENDIF  ! y dim inquire or create  ENDIF #5 
    
    IF (imod /= 10 .AND. imod /= 21 .AND. imod /= 22 .AND. imod /= 30 .AND. imod /= 13 .AND. imod /= 33) THEN !IF #7
      ier=nfmpi_inq_dimid(nfhd,"z",DID(3))
      IF (ier /= nf_noerr) THEN    !IF #8
        ier=nfmpi_def_dim(nfhd,"z",pl,DID(3))
        IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))&
                 print *,'PNCDF dim. z creat. stat ',nfmpi_strerror(ier)
      ENDIF  ! z dim create ENDIF #8
    ENDIF  ! z dim inquire or create #7 
    
    ier=nfmpi_inq_dimid(nfhd,"t",DID(4))
    IF (ier /= nf_noerr) THEN   !IF#9
      ier=nfmpi_def_dim(nfhd,"t",pt,DID(4))
      IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))&
               print *,'PNCDF dim. t creat. stat ',nfmpi_strerror(ier)
    ENDIF  !t dim create !ENDIF #9
    
    !==================================
    ! defining  dimensions done
    !==================================
    
    !==================================
    !Now define variables if necessary
    !==================================
    
    ier=nfmpi_inq_varid(nfhd,var_name,DCHVID)
    IF (ier /= nf_noerr) THEN  !IF #10 define variable if not defined
      
      IF(mype == 0) PRINT *,'PNCDF chunk '//file_name//' create var: '//var_name//' '
      
      
      IF (imod == 0) THEN      !if #11
        ier=nfmpi_def_var(nfhd,var_name,nfp,p4dim,DID,DCHVID)
        !           ELSE IF (imod.ge.1.and.imod.le.4.or.
        !     &              imod.ge.6.and.and.imod.le.8) THEN   !#11
      ELSE IF (imod >= 1 .AND. imod <= 9)THEN
        ier=nfmpi_def_var(nfhd,var_name,nfpd,p4dim,DID,DCHVID)
      ELSE IF (imod == 34) THEN
        DID3(1)=DID(1)
        DID3(2)=DID(2)
        DID3(3)=DID(3)
        ier=nfmpi_def_var(nfhd,var_name,nfp,p3dim,DID3,DCHVID)
      ELSE IF (imod < 20) THEN   !#11
        DID3(3)=DID(4)
        
        IF (imod == 10 .OR. imod == 13) THEN !IF #12
          DID3(1)=DID(1)
          DID3(2)=DID(2)
        ELSE IF (imod == 11) THEN  !#12
          DID3(1)=DID(1)
          DID3(2)=DID(3)
        ELSE IF (imod == 12) THEN !#12 
          DID3(1)=DID(2)
          DID3(2)=DID(3)
        ENDIF   !imod  10-13   !ENDIF #12
        
        ier=nfmpi_def_var(nfhd,var_name,nfpd,p3dim,DID3,DCHVID)
      ELSE IF (imod < 30) THEN !#11
        DID2(2)=DID(4)
        
        IF (imod == 20 .OR. imod == 23) THEN !IF #13
          DID2(1)=DID(3)
        ELSE IF (imod == 21) THEN !#13
          DID2(1)=DID(2)
        ELSE IF (imod == 22) THEN  !#13
          DID2(1)=DID(1)
        ENDIF !imod 20-23 !ENDIF #13 
        
        ier=nfmpi_def_var(nfhd,var_name,nfp,p2dim,DID2,DCHVID)
      ELSE IF (imod >= 30 .AND. imod <= 34) THEN
        DID1(1)=DID(4)
        ier=nfmpi_def_var(nfhd,var_name,nfp,p1dim,DID1,DCHVID)
      ELSE IF (imod == 35) THEN
        DID2(1)=DID(1)
        DID2(2)=DID(2)
        ier=nfmpi_def_var(nfhd,var_name,nfpd,p2dim,DID2,DCHVID)
      ENDIF  !defining dimension set and respective variable dep. on imod !#11
      
    ENDIF ! no existing variable definition ENDIF #10 
    
    !==================================
    !Now check again if variable exist
    !==================================
    ier=nfmpi_inq_varid(nfhd,var_name,DCHVID)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))    &
             &print *,'PNCDF chunk '//file_name//'_'//var_name//' inq ERROR: ', &
             nfmpi_strerror(ier)
    
    ier = nfmpi_enddef(nfhd)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF enddef chunk ',nfmpi_strerror(ier)
    
    ier = nfmpi_close(nfhd)
    IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
             print *,'PNCDF chunk '//file_name//' close stat',               &
             nfmpi_strerror(ier)
    
  !-------------------------
  ENDIF   !ifrmd == 0
  !-------------------------
  
  ier = nfmpi_open( icomm_cart,file_name,                       &
      NF_WRITE, MPI_INFO_NULL, nfhd)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
           &print *,'PNCDF chunk '//file_name//'_'//var_name//' open stat ',  &
           nfmpi_strerror(ier)
  
  ier=nfmpi_inq_varid(nfhd,var_name,DCHVID)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))   &
           &print *,'PNCDF chunk '//file_name//'_'//var_name//' inq. stat ',  &
           nfmpi_strerror(ier)
  
  ier=nfmpi_inq_unlimdim(nfhd,ihnd)       
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0)) &
           print *,'PNCDF chunk. inq unlim dim ',nfmpi_strerror(ier)
  
  ier=nfmpi_inq_dimlen(nfhd,ihnd,ifrmask)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0)) &
           print *,'PNCDF chunk. inq dimlen ',nfmpi_strerror(ier)
  
  !         iframe=ifrmask+1
  iframe=ifrmd+1
  
  !-------------------------
  !  different mods
  !-------------------------
  IF (imod <= 1) THEN
    START4(1)=nsubpos+1
    START4(2)=msubpos+1
    START4(3)=lsubpos+1
    START4(4)=iframe
    COUNT4(1)=np
    COUNT4(2)=mp
    COUNT4(3)=lp
    COUNT4(4)=1
  ELSE IF (imod >= 2 .AND. imod <= 4) THEN
    START4(1)=nsubpos+1
    START4(2)=msubpos+1
    START4(3)=lsubpos+1
    START4(4)=iframe
    COUNT4(1)=np+ia1*rightedge
    COUNT4(2)=mp+ia2*topedge
    COUNT4(3)=lp+ia3*skyedge
    COUNT4(4)=1
  ELSE IF (imod == 34) THEN
    START3(1)=nsubpos+1
    START3(2)=msubpos+1
    START3(3)=lsubpos+1
    COUNT3(1)=np
    COUNT3(2)=mp
    COUNT3(3)=lp
  ELSE IF (imod == 6) THEN
    START4(1)=1
    START4(2)=msubpos+1
    START4(3)=lsubpos+1
    START4(4)=iframe
    iedge=max(leftedge,rightedge)
    COUNT4(1)=ih*iedge
    COUNT4(2)=mp*iedge
    COUNT4(3)=lp*iedge
    COUNT4(4)=1 *iedge
  ELSE IF (imod == 7) THEN
    START4(1)=nsubpos+1
    START4(2)=1
    START4(3)=lsubpos+1
    START4(4)=iframe
    iedge=max(botedge,topedge)
    COUNT4(1)=np*iedge
    COUNT4(2)=ih*iedge
    COUNT4(3)=lp*iedge
    COUNT4(4)=1 *iedge
  ELSE IF (imod == 8) THEN
    START4(1)=nsubpos+1
    START4(2)=msubpos+1
    START4(3)=1
    START4(4)=iframe
    iedge=max(gndedge,skyedge)
    COUNT4(1)=np*iedge
    COUNT4(2)=mp*iedge
    COUNT4(3)=ih*iedge
    COUNT4(4)=1 *iedge
  ELSE IF (imod == 9) THEN
    START4(1)=1+ ih*rightedge 
    START4(2)=1+ ih*topedge
    START4(3)=1+ ih*skyedge  
    START4(4)=iframe
    icorn=0
    IF(gndedge == 1 .AND. leftedge  == 1 .AND. botedge == 1 .OR.          &
               gndedge.eq.1.and.rightedge.eq.1.and.botedge.eq.1.or.         &
               gndedge.eq.1.and.leftedge .eq.1.and.topedge.eq.1.or.         &
               gndedge.eq.1.and.rightedge.eq.1.and.topedge.eq.1.or.         &
               skyedge.eq.1.and.leftedge .eq.1.and.botedge.eq.1.or.         &
               skyedge.eq.1.and.rightedge.eq.1.and.botedge.eq.1.or.         &
               skyedge.eq.1.and.leftedge .eq.1.and.topedge.eq.1.or.         &
               skyedge.eq.1.and.rightedge.eq.1.and.topedge.eq.1) icorn=1
    COUNT4(1)=ih*icorn
    COUNT4(2)=ih*icorn
    COUNT4(3)=ih*icorn
    COUNT4(4)=1 *icorn
  ELSE IF (imod == 10 .OR. imod == 13) THEN
    START3(1)=nsubpos+1
    START3(2)=msubpos+1
    START3(3)=iframe
    COUNT3(1)=np
    COUNT3(2)=mp
    COUNT3(3)=1
  ELSE IF (imod == 11) THEN
    START3(1)=nsubpos+1
    START3(2)=lsubpos+1
    START3(3)=iframe
    COUNT3(1)=np
    COUNT3(2)=lp
    COUNT3(3)=1
  ELSE IF (imod == 12) THEN
    START3(1)=msubpos+1
    START3(2)=lsubpos+1
    START3(3)=iframe
    COUNT3(1)=mp
    COUNT3(2)=lp
    COUNT3(3)=1
  ELSE IF (imod == 20 .OR. imod == 23) THEN
    START2(1)=lsubpos+1
    START2(2)=iframe
    COUNT2(1)=lp
    COUNT2(2)=1
  ELSE IF (imod == 21) THEN
    START2(1)=msubpos+1
    START2(2)=iframe
    COUNT2(1)=mp
    COUNT2(2)=1
  ELSE IF (imod == 22) THEN
    START2(1)=nsubpos+1
    START2(2)=iframe
    COUNT2(1)=np
    COUNT2(2)=1
  ELSE IF (imod == 35) THEN
    START2(1)=nsubpos+1
    START2(2)=msubpos+1
    COUNT2(1)=np
    COUNT2(2)=mp
  ENDIF
  i1p=0 !initialize to remove compiler warning 
  i2p=0
  i3p=0
  i1=0
  i2=0
  i3=0
  IF (imod > 3) THEN
    i1p=(ia1-1)/np+1
    i2p=(ia2-1)/mp+1
    i3p=(ia3-1)/lp+1
    i1=ia1-(i1p-1)*np
    i2=ia2-(i2p-1)*mp
    i3=ia3-(i3p-1)*lp
  ENDIF
  IF (imod == 0) THEN
    CALL pnettrans(var,dane3d)
    ier=nfmpi_put_vara_real_all(nfhd,DCHVID,START4,COUNT4,dane3d)
  ELSE IF(imod == 1) THEN
    CALL pnettrans8(var,dane83d)
    ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,dane83d)
  ELSE IF(imod == 2) THEN
    IF ((rightedge == 1 .AND. ia1 > 0) .OR.    (  topedge == 1 .AND. ia2 > 0) .OR.        (skyedge   == 1 .AND. ia3 > 0)) THEN
      DO k=1,lp+ia3*skyedge
        DO j=1,mp+ia2*topedge
          DO i=1,np+ia1*rightedge
            danea83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,danea83d)
    ELSE
      DO k=1,lp
        DO j=1,mp
          DO i=1,np
            dane83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,dane83d)
    ENDIF
  ELSE IF(imod == 3) THEN
    
    IF ((rightedge == 1 .AND. ia1 > 0) .AND. (skyedge == 0 .AND. ia3 > 0)) THEN 
      DO k=1,lp
        DO j=1,mp+ia2*topedge
          DO i=1,np+ia1*rightedge
            danex83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,danex83d)
    ELSE IF((rightedge == 0 .AND. ia1 > 0) .AND. (skyedge == 1 .AND. ia3 > 0)) THEN
      DO k=1,lp+ia3*skyedge
        DO j=1,mp+ia2*topedge
          DO i=1,np
            danez83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,danez83d)
    ELSE IF((rightedge == 1 .AND. ia1 > 0) .AND. (skyedge == 1 .AND. ia3 > 0)) THEN
      DO k=1,lp+ia3*skyedge
        DO j=1,mp+ia2*topedge
          DO i=1,np+ia1*rightedge
            danea83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,danea83d)
    ELSE IF((ia2 > 0) .AND. (skyedge == 1 .AND. ia3 > 0)) THEN
      DO k=1,lp+ia3*skyedge
        DO j=1,mp+ia2*topedge
          DO i=1,np
            danez83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,danez83d)
    ELSE IF((ia2 > 0) .AND. (rightedge == 1 .AND. ia1 > 0)) THEN
      DO k=1,lp
        DO j=1,mp+ia2*topedge
          DO i=1,np+ia1*rightedge
            danex83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,danex83d)
    ELSE
      DO k=1,lp
        DO j=1,mp
          DO  i=1,np
            dane83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,dane83d)
    ENDIF
  ELSE IF(imod == 4) THEN
    IF (rightedge == 1 .AND. skyedge == 0) THEN
      DO k=1,lp
        DO j=1,mp
          DO i=1,np+ia1
            danex83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,danex83d)
    ELSE IF(rightedge == 0 .AND. skyedge == 1) THEN
      DO k=1,lp+ia3
        DO j=1,mp
          DO i=1,np
            danez83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,danez83d)
    ELSE IF(rightedge == 1 .AND. skyedge == 1) THEN
      DO  k=1,lp+ia3
        DO j=1,mp
          DO i=1,np+ia1
            danea83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,danea83d)
    ELSE IF(rightedge == 0 .AND. skyedge == 0) THEN
      DO k=1,lp
        DO  j=1,mp
          DO i=1,np
            dane83d(i,j,k)=var(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,dane83d)
    ENDIF
  ELSE IF(imod == 6) THEN
    IF (leftedge == 1) THEN
      DO  k=1,lp
        DO j=1,mp
          DO  i=1-ih,0
            halo8lr(ih+i,j,k)=var(i,j,k)    
          ENDDO
        ENDDO
      ENDDO
      !      print 99,mype,START4(1),COUNT4(1),START4(2),COUNT4(2),
      !     &                 START4(3),COUNT4(3)
      !   99 format('LWR',i2,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3)
    ENDIF
    IF (rightedge == 1) THEN
      DO  k=1,lp
        DO  j=1,mp
          DO  i=1-ih,0
            halo8lr(ih+i,j,k)=var(np+ih+i,j,k)    
          ENDDO
        ENDDO
      ENDDO
      START4(1)=1+ih
      !      print 98,mype,START4(1),COUNT4(1),START4(2),COUNT4(2),
      !     &                 START4(3),COUNT4(3)
      !   98 format('RWR',i2,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3)
    ENDIF
    CALL mybarrier()
    ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,halo8lr)
  ELSE IF(imod == 7) THEN
    IF (botedge == 1) THEN
      DO  k=1,lp
        DO  j=1-ih,0
          DO  i=1,np
            halo8bt(i,ih+j,k)=var(i,j,k)    
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF (topedge == 1) THEN
      DO  k=1,lp
        DO  j=1-ih,0
          DO  i=1,np
            halo8bt(i,ih+j,k)=var(i,mp+ih+j,k)    
          ENDDO
        ENDDO
      ENDDO
      START4(2)=1+ih
    ENDIF
    CALL mybarrier()
    ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,halo8bt)
  ELSE IF(imod == 8) THEN
    IF (gndedge == 1) THEN
      DO  k=1-ih,0
        DO  j=1,mp
          DO  i=1,np
            halo8gs(i,j,ih+k)=var(i,j,k)    
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF (skyedge == 1) THEN
      DO  k=1-ih,0
        DO  j=1,mp
          DO i=1,np
            halo8gs(i,j,ih+k)=var(i,j,lp+ih+k)    
          ENDDO
        ENDDO
      ENDDO
      START4(3)=1+ih
    ENDIF
    CALL mybarrier()
    ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,halo8gs)
  ELSE IF(imod == 9) THEN
    DO k=1,ih
      DO j=1,ih
        DO i=1,ih
          ii=i-ih*leftedge+np*rightedge
          jj=j-ih*botedge +mp*topedge
          kk=k-ih*gndedge +lp*skyedge
          halo8cr(i,j,k)=var(ii,jj,kk)
        ENDDO
      ENDDO
    ENDDO
    ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,halo8cr)
  ELSE IF(imod == 10 .AND. i3p == lpos) THEN
    DO j=1,mp
      DO i=1,np
        dane2dxy(i,j)=var(i,j,i3)
      ENDDO
    ENDDO
    ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START3,COUNT3,dane2dxy)
  ELSE IF(imod == 13 .AND. lpos == 1) THEN
    DO j=1,mp
      DO i=1,np
        dane2dxy(i,j)=var(i,j,1)
      ENDDO
    ENDDO
    ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START3,COUNT3,dane2dxy)
  ELSE IF(imod == 11) THEN
    DO k=1,lp
      DO i=1,np
        dane2dxz(i,k)=var(i,i2,k)
      ENDDO
    ENDDO
    IF (i2p /= mpos) THEN
      COUNT3(1)=0
      COUNT3(2)=0
      COUNT3(3)=0
    ENDIF
    ier = nfmpi_put_vara_real_all(nfhd,DCHVID,START3,COUNT3,dane2dxz)
  ELSE IF(imod == 12) THEN
    DO k=1,lp
      DO j=1,mp
        dane2dyz(j,k)=var(i1,j,k)
      ENDDO
    ENDDO
    IF (i1p /= npos) THEN 
      COUNT3(1)=0
      COUNT3(2)=0
      COUNT3(3)=0
    ENDIF  
    ier = nfmpi_put_vara_real_all(nfhd,DCHVID,START3,COUNT3,dane2dyz)
  ELSE IF(imod == 20) THEN
    DO k=1,lp
      dane1dz(k)=var(i1,i2,k)
    ENDDO
    IF (i1p /= npos .OR. i2p /= mpos) THEN
      COUNT2(1)=0
      COUNT2(2)=0
    ENDIF
    ier = nfmpi_put_vara_real_all(nfhd,DCHVID,START2,COUNT2,dane1dz)
  ELSE IF(imod == 23) THEN
    DO k=1,lp
      dane1dz(k)=v1d(k)
    ENDDO 
    IF (1 /= npos .OR. 1 /= mpos) THEN
      COUNT2(1)=0
      COUNT2(2)=0
    ENDIF
    
    ier = nfmpi_put_vara_real_all(nfhd,DCHVID,START2,COUNT2,dane1dz)
    
  ELSE IF(imod == 21) THEN
    DO j=1,mp
      dane1dy(j)=var(i1,j,i3)
    ENDDO
    IF (i1p /= npos .OR. i3p /= lpos) THEN
      COUNT2(1)=0
      COUNT2(2)=0
    ENDIF
    
    ier=nfmpi_put_vara_real_all(nfhd,DCHVID,START2,COUNT2,dane1dy)
  ELSE IF(imod == 22) THEN
    DO i=1,np
      dane1dx(i)=var(i,i2,i3)
    ENDDO
    IF (i2p /= mpos .OR. i3p /= lpos) THEN
      COUNT2(1)=0
      COUNT2(2)=0
    ENDIF
    ier=nfmpi_put_vara_real_all(nfhd,DCHVID,START2,COUNT2,dane1dx)
  ELSE IF(imod == 34) THEN
    CALL pnettrans(var,dane3d)
    ier=nfmpi_put_vara_real_all(nfhd,DCHVID,START3,COUNT3,dane3d)
  ELSE IF(imod == 35) THEN
    DO j=1,mp
      DO i=1,np
        dane2dxy(i,j)=var(i,j,1-ih)
      ENDDO
    ENDDO
    ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START2,COUNT2,dane2dxy)
    
  ENDIF
  
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF chunk '//var_name//' put stat ',               &
           nfmpi_strerror(ier)
  
  ier = nfmpi_close(nfhd)
  IF((mype == 0 .AND. iPRINT == 1) .OR. (iPRINT == 0 .AND. ier /= 0))     &
           print *,'PNCDF chunk '//file_name//' close stat',               &
           nfmpi_strerror(ier)
  
END SUBROUTINE pnet_out_chunk


!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnettrans(sour,dest)
!--------------------------------------------------------------------!
  REAL(KIND=euwp) sour(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
  real*4  dest(np,mp,lp)
  dest(1:np,1:mp,1:lp)=sour(1:np,1:mp,1:lp)
END SUBROUTINE pnettrans


!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnettrans8(sour,dest)
!--------------------------------------------------------------------!
  REAL(KIND=euwp) sour(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
  real*8  dest(np,mp,lp)
  dest(1:np,1:mp,1:lp)=sour(1:np,1:mp,1:lp)
END SUBROUTINE pnettrans8

!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnettransr8(sour,dest)
!--------------------------------------------------------------------!
  REAL(KIND=euwp) :: dest(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
  real*8  sour(np,mp,lp)
  dest(1:np,1:mp,1:lp)=sour(1:np,1:mp,1:lp)
END SUBROUTINE pnettransr8 


#endif /*1==0*/
END MODULE eulag_diagutils

MODULE testing
   USE precisions
   USE mpi_parallel, ONLY: mype,globsum,globmax,globmin
   USE mpi_parallel, ONLY: leftedge,rightedge,botedge,topedge
   USE mpi_parallel, ONLY: nsubpos,msubpos,lsubpos
   USE module_velprd, ONLY: velprdA_to_C_ver2
!   USE geometry, ONLY: physical2contravariant
   USE   parameters, ONLY: dx,dy,dz,dt,dx00,nt,nmli,timetot
#ifdef PNETCDF
   USE mpi_parallel, ONLY: pnet_out_chunk
#endif /*PNETCDF*/
   IMPLICIT NONE
   INTEGER ncnt

CONTAINS
#include "../src_algorithms/defines.inc"
   SUBROUTINE set_initial_tracer(uadv_c,vadv_c,wadv_c,x_data,rho,np,mp,lp,ih)
     INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL_euwp, DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: &
         uadv_c,vadv_c,wadv_c,x_data,rho
#if(TEST == 1)
      CALL set_initial84(uadv_c,vadv_c,wadv_c,x_data,rho,np,mp,lp,ih)
#elif(TEST == 10)
      CALL set_initial91(uadv_c,vadv_c,wadv_c,x_data,rho,np,mp,lp,ih)
#else
      CALL set_initial15(uadv_c,vadv_c,wadv_c,x_data,rho,np,mp,lp,ih)
#endif /*TEST CHOICE*/
   END SUBROUTINE set_initial_tracer

   ! Calls in subroutine SET_INITIAL84: 
   ! => velprd (on line <84>)
   SUBROUTINE set_initial84(uadv_c,vadv_c,wadv_c,x_data,rho,np,mp,lp,ih)
      USE scratch_datafields, ONLY: xtest
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      INTEGER(KIND=iintegers) :: i,j,k,ia,ja,ka
      REAL_euwp, DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: &
         uadv_a,vadv_a,wadv_a,uadv_c,vadv_c,wadv_c,x_data,rho
      REAL(KIND=euwp) :: d,r_squared,h_sphere,tmp,sqrt6i
      REAL(KIND=euwp) :: x0,y0,z0,xc,yc,zc,xp,yp,zp
      REAL(KIND=euwp) :: omega,omegax,omegay,omegaz
      sqrt6i=1._euwp/sqrt(6._euwp)
      d=7._euwp*sqrt6i*dx
      x0=.5_euwp*dx00
      y0=.5_euwp*dx00
      z0=.5_euwp*dx00
      xc=x0-d
      yc=y0-d
      zc=z0+2._euwp*d
      r_squared=7._euwp*dx
      r_squared=r_squared*r_squared
!  print *,'xc,yc,zc',xc,yc,zc
!  print *,'x0,y0,z0',x0,y0,z0
!  print *,'r_squared,dx00',r_squared,dx00
      h_sphere=4._euwp
      omega = 0.1_euwp
      omegax = 0.5_euwp*omega
      omegay = 0.5_euwp*omega
      omegaz = omega/sqrt(2._euwp)

      x_data(:,:,:)=0._euwp
      DO k=1,lp
         DO j=1,mp
            DO i=1,np
               ia=nsubpos+i 
               ja=msubpos+j 
               ka=lsubpos+k 
               xp=ia*dx
               yp=ja*dy
               zp=ka*dz
               tmp=(xp-xc)*(xp-xc)+(yp-yc)*(yp-yc)+(zp-zc)*(zp-zc)
               IF(tmp < r_squared) x_data(i,j,k)=h_sphere*sqrt(tmp/r_squared)
               uadv_a(i,j,k)=( - omegaz*(ja*dy - y0) + omegay*(ka * dz - z0))* dt / dx
               vadv_a(i,j,k)=(   omegaz*(ia*dx - x0) - omegax*(ka * dz - z0))* dt / dy
               wadv_a(i,j,k)=( - omegay*(ia*dx - x0) + omegax*(ja * dy - y0))* dt / dz
            ENDDO
         ENDDO
      ENDDO
! Store initial condition for testing of the result
      xtest(1:np,1:mp,1:lp)=x_data(1:np,1:mp,1:lp)
#ifdef PNETCDF
!    call pnet_out_chunk('xini     ','tinit.nc',1,1,1,1,0,0,x_data,np,mp,lp,ih)
#endif /*PNETCDF*/
! Evaluate staggered C-grid velocity components
      CALL velprdA_to_C_ver2(uadv=uadv_a,vadv=vadv_a,wadv=wadv_a, &
                             u1  =uadv_c,u2  =vadv_c,u3  =wadv_c, &
                             ipoles0=0,ibcx0=0,ibcy0=0,ibcz0=0,   &
                             np=np,mp=mp,lp=lp,ih=ih)
   END SUBROUTINE set_initial84

   ! Calls in subroutine SET_INITIAL91: 
   SUBROUTINE set_initial91(uadv_c,vadv_c,wadv_c,x_data,rho,np,mp,lp,ih)
      USE scratch_datafields, ONLY: xtest
      USE parameters, ONLY: rds,n,dxa
      USE geometry, ONLY: xcr,ycr,cosa,sina,sinx,cosx !,coscx,cosCy
      USE geometry, ONLY: g11,g12,g13,g21,g22,g23,g33
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      INTEGER(KIND=iintegers) :: i,j,k,ia,ja,ka
      REAL_euwp, DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: &
         uadv_a,vadv_a,wadv_a,uadv_c,vadv_c,wadv_c,x_data,rho
      REAL(KIND=euwp) u_phys,v_phys,w_phys
      REAL(KIND=euwp) :: r_squared,tmp
      REAL(KIND=euwp) :: xc,yc,zc,xp,yp
      
      REAL(KIND=euwp) Deg2Rad,radius,Umax, wind_angle
      REAL(KIND=euwp) zbeta,zvel,pi
      
      
!Definitions
      pi=acos(-1.)
      radius=rds
      Deg2Rad=pi/180.
!Initially signal centered at:
      xc= -0.5_euwp*pi
      yc= 0._euwp
      zc= 0._euwp

      r_squared=7._euwp*dxa
      r_squared=r_squared*r_squared

      Umax = pi/n
      wind_angle = 270._euwp
      zvel   = Umax!/radius
      zbeta  = wind_angle*Deg2Rad
x_data(:,:,:)=0._euwp
      ncnt=0
      DO k=1,lp
         DO j=1,mp
            DO i=1,np
               ia=nsubpos+i 
               ja=msubpos+j 
               ka=lsubpos+k 
               xp=xcr(i,j)/rds
               yp=ycr(i,j)/rds
! GMD 15 way
!  tmp=(xp-xc)*(xp-xc)+(yp-yc)*(yp-yc)
               tmp=2._euwp*(cosa(i,j)*cosa(i,j)*sin(.5_euwp*(xp-xc))**2 &
                  +sin(.5_euwp*(yp-yc))**2 )
!  print *,ia,ja,tmp/(rds*rds),r_squared/(rds*rds)
               IF(tmp <= r_squared) THEN
                  x_data(i,j,k)=1._euwp - sqrt(tmp/r_squared)
                  ncnt=ncnt+1
               ENDIF
! Square tracer:
!  xcrmx=tracer_lonmax*Deg2Rad*rds
!  xcrmn=tracer_lonmin*Deg2Rad*rds
!  ycrmx=tracer_latmax*Deg2Rad*rds
!  ycrmn=tracer_latmin*Deg2Rad*rds

!  if (  ycr(i,j)<=ycrmx .and. &
!        ycr(i,j)>=ycrmn .and. &
!        xcr(i,j)<=xcrmx .and. &
!        xcr(i,j)>=xcrmn ) then
!    x_data(i,j,k) = 1._euwp
!   endif
               u_phys = zvel*(cos(zbeta)*cosa(i,j)+sina(i,j)*cosx(i,j)*sin(zbeta))*radius
               v_phys = -zvel*sinx(i,j)*sin(zbeta)*radius
               w_phys =0._euwp

!Evaluate contravariant velocities from physical velocities
               uadv_a(i,j,k) =rho(i,j,k)*( g11(i,j,k)*u_phys + g21(i,j,k)*v_phys ) * dt/dx
               vadv_a(i,j,k) =rho(i,j,k)*( g12(i,j,k)*u_phys + g22(i,j,k)*v_phys ) * dt/dy
               wadv_a(i,j,k) =rho(i,j,k)*( g13(i,j,k)*u_phys + g23(i,j,k)*v_phys &
                  + g33(i,j,k)*w_phys ) * dt/dz
            ENDDO
         ENDDO
      ENDDO
! Store initial condition for testing of the result
      xtest(1:np,1:mp,1:lp)=x_data(1:np,1:mp,1:lp)
! Evaluate staggered C-grid velocity components
      CALL velprd_sphere(uadv_a,vadv_a,wadv_a,  &
                         uadv_c,vadv_c,wadv_c,rho,np,mp,lp,ih)
   END SUBROUTINE set_initial91

   SUBROUTINE set_initial15(uadv_c,vadv_c,wadv_c,x_data,rho,np,mp,lp,ih)
      USE scratch_datafields, ONLY: xtest
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      INTEGER(KIND=iintegers) :: i,j,k,ia,ja,ka
      REAL_euwp, DIMENSION(:,:,:),ALLOCATABLE :: &
         uadv_a,vadv_a,wadv_a
      REAL_euwp, DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: &
         x_data,rho
   REAL_euwp,DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+ih  ) :: uadv_c
   REAL_euwp,DIMENSION(1-ih:np+ih  , 1-ih:mp+1+ih, 1-ih:lp+ih  ) :: vadv_c
   REAL_euwp,DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+1+ih) :: wadv_c

      REAL(KIND=euwp) :: r_squared,h_sphere,tmp,sqrt3i
      REAL(KIND=euwp) :: x0,y0,z0,xc,yc,zc,xp,yp,zp
      REAL(KIND=euwp) :: omega,omegax,omegay,omegaz
      INTEGER iallocerr,ierr 
!      124   FORMAT('Driver will run on field with size ',i0,' x ',i0,' x ',i0)
!   WRITE (*,124)  np,mp,lp
! print *,'Allocation with np,mp,lp,ih:',np,mp,lp,ih
      ierr=0
      ALLOCATE(uadv_a(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih), STAT=iallocerr);ierr=ierr+iallocerr
      ALLOCATE(vadv_a(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih), STAT=iallocerr);ierr=ierr+iallocerr
      ALLOCATE(wadv_a(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih), STAT=iallocerr);ierr=ierr+iallocerr
      IF(ierr.ne.0) print *,'Allocation errors:',ierr

      
      sqrt3i=1._euwp/sqrt(3._euwp)
      x0=.5_euwp*dx00
      y0=.5_euwp*dx00
      z0=.5_euwp*dx00
      xc= x0-0.25*dx00*sqrt3i
      yc= y0+0.25*dx00*sqrt3i
      zc= z0+0.25*dx00*sqrt3i
      r_squared=15._euwp
      r_squared=r_squared*r_squared
! print *,'xc,yc,zc',xc,yc,zc
! print *,'x0,y0,z0',x0,y0,z0
! print *,'r_squared,dx00',r_squared,dx00
      h_sphere=4._euwp
      omega  = 0.1_euwp/sqrt(3._euwp)
      omegax = omega
      omegay = omega
      omegaz = omega

      x_data(:,:,:)=0._euwp
      ncnt=0
      DO k=1,lp
         DO j=1,mp
            DO i=1,np
               ia=nsubpos +i 
               ja=msubpos +j 
               ka=lsubpos +k 
               xp=ia*dx
               yp=ja*dy
               zp=ka*dz
               tmp=(xp-xc)*(xp-xc)+(yp-yc)*(yp-yc)+(zp-zc)*(zp-zc)
               IF(tmp <= r_squared) THEN
                  x_data(i,j,k)=h_sphere
                  ncnt=ncnt+1
               ENDIF
               uadv_a(i,j,k)=( - omegaz*(ja*dy - y0) + omegay*(ka * dz - z0))* dt / dx
               vadv_a(i,j,k)=(   omegaz*(ia*dx - x0) - omegax*(ka * dz - z0))* dt / dy
               wadv_a(i,j,k)=( - omegay*(ia*dx - x0) + omegax*(ja * dy - y0))* dt / dz
            ENDDO
         ENDDO
      ENDDO
! Store initial condition for testing of the result
      xtest(1:np,1:mp,1:lp)=x_data(1:np,1:mp,1:lp)
#ifdef PNETCDF
!ALL pnet_out_chunk('uadv','mpset.nc',1,1,1,1,0,0,uadv_a,np,mp,lp,ih)
#endif /*PNETCDF*/
! Evaluate staggered C-grid velocity components
!     CALL velprd(uadv_a,vadv_a,wadv_a,    &
!                 uadv_c,vadv_c,wadv_c,rho,np,mp,lp,ih)
      CALL velprdA_to_C_ver2(uadv=uadv_a,vadv=vadv_a,wadv=wadv_a, &
                             u1  =uadv_c,u2  =vadv_c,u3  =wadv_c, &
                             ipoles0=0,ibcx0=0,ibcy0=0,ibcz0=0,   &
                             np=np,mp=mp,lp=lp,ih=ih)
      DEALLOCATE (uadv_a)
      DEALLOCATE (vadv_a)
      DEALLOCATE (wadv_a)
   END SUBROUTINE set_initial15

   SUBROUTINE set_initial15_perturbed(noise_amp,noise_phase,uadv_c,vadv_c,wadv_c,x_data,rho,np,mp,lp,ih)
      USE scratch_datafields, ONLY: xtest
      USE parameters, ONLY: n,m,l,pi2 
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      INTEGER(KIND=iintegers) :: i,j,k,ia,ja,ka
      REAL_euwp, DIMENSION(:,:,:),ALLOCATABLE :: &
         uadv_a,vadv_a,wadv_a
      REAL_euwp, DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: &
         uadv_c,vadv_c,wadv_c,x_data,rho,noise
      REAL(KIND=euwp) :: r_squared,h_sphere,tmp,sqrt3i
      REAL(KIND=euwp) :: x0,y0,z0,xc,yc,zc,xp,yp,zp
      REAL(KIND=euwp) :: omega,omegax,omegay,omegaz
      REAL(KIND=euwp) :: noise_amp,noise_phase
      REAL(KIND=euwp),parameter ::  vel_noise_amp=0.0001_euwp 
      INTEGER iallocerr,ierr 
!      124   FORMAT('Driver will run on field with size ',i0,' x ',i0,' x ',i0)
!   WRITE (*,124)  np,mp,lp
  print *,'Allocation with np,mp,lp,ih:',np,mp,lp,ih
      ierr=0
      ALLOCATE(uadv_a(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih), STAT=iallocerr);ierr=ierr+iallocerr
      ALLOCATE(vadv_a(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih), STAT=iallocerr);ierr=ierr+iallocerr
      ALLOCATE(wadv_a(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih), STAT=iallocerr);ierr=ierr+iallocerr
  print *,'Allocation status:',ierr

      
      sqrt3i=1._euwp/sqrt(3._euwp)
      x0=.5_euwp*dx00
      y0=.5_euwp*dx00
      z0=.5_euwp*dx00
      xc= x0-0.25*dx00*sqrt3i
      yc= y0+0.25*dx00*sqrt3i
      zc= z0+0.25*dx00*sqrt3i
      r_squared=15._euwp
      r_squared=r_squared*r_squared
!  print *,'xc,yc,zc',xc,yc,zc
!  print *,'x0,y0,z0',x0,y0,z0
!  print *,'r_squared,dx00',r_squared,dx00
      h_sphere=4._euwp
      omega  = 0.1_euwp/sqrt(3._euwp)
      omegax = omega
      omegay = omega
      omegaz = omega

      x_data(:,:,:)=0._euwp
      ncnt=0
      DO k=1,lp
         DO j=1,mp
            DO i=1,np
               ia=nsubpos+i 
               ja=msubpos+j 
               ka=lsubpos+k 
               xp=ia*dx
               yp=ja*dy
               zp=ka*dz
               tmp=(xp-xc)*(xp-xc)+(yp-yc)*(yp-yc)+(zp-zc)*(zp-zc)
               IF(tmp <= r_squared) THEN
                  x_data(i,j,k)=h_sphere
                  ncnt=ncnt+1
               ENDIF
              ia=nsubpos+i
              ja=msubpos+j
              ka=lsubpos+k
               noise(i,j,k)=noise_amp*sin(pi2*ia/n+noise_phase*pi2) &
                                     *cos(pi2*ja/m+2.*noise_phase*pi2) &
                                     *sin(pi2*ka/l+3.*noise_phase*pi2)
               x_data(i,j,k)=x_data(i,j,k)+noise(i,j,k)
               uadv_a(i,j,k)=( - omegaz*(ja*dy - y0) + omegay*(ka * dz - z0))* dt / dx
               uadv_a(i,j,k)=uadv_a(i,j,k)+vel_noise_amp*noise(i,j,k)
               vadv_a(i,j,k)=(   omegaz*(ia*dx - x0) - omegax*(ka * dz - z0))* dt / dy
               vadv_a(i,j,k)=vadv_a(i,j,k)+vel_noise_amp*noise(i,j,k)
               wadv_a(i,j,k)=( - omegay*(ia*dx - x0) + omegax*(ja * dy - y0))* dt / dz
               wadv_a(i,j,k)=wadv_a(i,j,k)+vel_noise_amp*noise(i,j,k)
            ENDDO
         ENDDO
      ENDDO
! Store initial condition for testing of the result
      xtest(1:np,1:mp,1:lp)=x_data(1:np,1:mp,1:lp)
#ifdef PNETCDF
#endif /*PNETCDF*/
! Evaluate staggered C-grid velocity components
      CALL velprdA_to_C_ver2(uadv=uadv_a,vadv=vadv_a,wadv=wadv_a, &
                             u1  =uadv_c,u2  =vadv_c,u3  =wadv_c, &
                             ipoles0=0,ibcx0=0,ibcy0=0,ibcz0=0,   &
                             np=np,mp=mp,lp=lp,ih=ih)
      DEALLOCATE (uadv_a)
      DEALLOCATE (vadv_a)
      DEALLOCATE (wadv_a)
   END SUBROUTINE set_initial15_perturbed

   SUBROUTINE set_velocities_for_velprd(ox,oy,oz,rho,np,mp,lp,ih)
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      INTEGER(KIND=iintegers) :: i,j,k,ia,ja,ka
      REAL_euwp, DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2),INTENT(OUT) :: &
         ox,oy,oz 
      REAL_euwp, DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: &
         uadv_a,vadv_a,wadv_a,rho
      REAL(KIND=euwp) :: x0,y0,z0
      REAL(KIND=euwp) :: omega,omegax,omegay,omegaz
!      124   FORMAT('Driver will run on field with size ',i0,' x ',i0,' x ',i0)
!   WRITE (*,124)  np,mp,lp

      
      x0=.5_euwp*dx00
      y0=.5_euwp*dx00
      z0=.5_euwp*dx00
      omega  = 0.1_euwp/sqrt(3._euwp)
      omegax = omega
      omegay = omega
      omegaz = omega

      DO k=1,lp
         DO j=1,mp
            DO i=1,np
               ia=nsubpos+i 
               ja=msubpos+j 
               ka=lsubpos+k 
               uadv_a(i,j,k)=( - omegaz*(ja*dy - y0) + omegay*(ka * dz - z0))* dt / dx
               vadv_a(i,j,k)=(   omegaz*(ia*dx - x0) - omegax*(ka * dz - z0))* dt / dy
               wadv_a(i,j,k)=( - omegay*(ia*dx - x0) + omegax*(ja * dy - y0))* dt / dz
            ENDDO
         ENDDO
      ENDDO
      CALL physical2contravariant(  u_phys=uadv_a     ,  v_phys=vadv_a     ,  w_phys=wadv_a, &
                                  ox_contr=ox(:,:,:,0),oy_contr=oy(:,:,:,0),oz_contr=oz(:,:,:,0), &
                                  np=np,mp=mp,lp=lp,ih=ih)
      omega  = 0.1001_euwp/sqrt(3._euwp)
      omegax = omega
      omegay = omega
      omegaz = omega
      DO k=1,lp
         DO j=1,mp
            DO i=1,np
               ia=nsubpos+i 
               ja=msubpos+j 
               ka=lsubpos+k 
               uadv_a(i,j,k)=rho(i,j,k)*( - omegaz*(ja*dy - y0) + omegay*(ka * dz - z0))* dt / dx
               vadv_a(i,j,k)=rho(i,j,k)*(   omegaz*(ia*dx - x0) - omegax*(ka * dz - z0))* dt / dy
               wadv_a(i,j,k)=rho(i,j,k)*( - omegay*(ia*dx - x0) + omegax*(ja * dy - y0))* dt / dz
            ENDDO
         ENDDO
      ENDDO
      CALL physical2contravariant(  u_phys=uadv_a     ,  v_phys=vadv_a     ,  w_phys=wadv_a, &
                                  ox_contr=ox(:,:,:,2),oy_contr=oy(:,:,:,2),oz_contr=oz(:,:,:,2), &
                                  np=np,mp=mp,lp=lp,ih=ih)
   END SUBROUTINE set_velocities_for_velprd

   SUBROUTINE velprd_sphere(uadv_a,vadv_a,wadv_a, &
                            uadv_c,vadv_c,wadv_c,rho,np,mp,lp,ih)
      USE mpi_parallel
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL_euwp, DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: &
         uadv_a,vadv_a,wadv_a,uadv_c,vadv_c,wadv_c,rho
      REAL(KIND=euwp) cr1(1),cr1_tot(1)
!Definition of the advecting velocity on staggered C grid
      CALL updatelr(uadv_a,np,mp,lp,np,mp,lp,1,ih)
      uadv_c(1:np,1:mp,1:lp)=0.5_euwp*(uadv_a(1:np  ,1:mp,1:lp)  &
                                      +uadv_a(0:np-1,1:mp,1:lp))
      CALL update(uadv_c,np,mp,lp,np,mp,lp,1,ih)

      CALL updatebt(vadv_a,np,mp,lp,np,mp,lp,1,ih)
      vadv_c(1:np,1:mp,1:lp)=0.5_euwp*(vadv_a(1:np,1:mp  ,1:lp)  &
                                      +vadv_a(1:np,0:mp-1,1:lp))
      IF(botedge == 1) vadv_c(1:np,1   ,1:lp)=0._euwp !Jacobian vanishes at the pole where vadv_c is specified
      CALL update(vadv_c,np,mp,lp,np,mp,lp,1,ih)
      IF(topedge == 1)vadv_c(1:np,mp+1,1:lp)=0._euwp !Jacobian vanishes at the pole where vadv_c is specified

      CALL updategs(wadv_a,np,mp,lp,np,mp,lp,1,ih)
      IF(gndedge == 1) &
         wadv_a(1:np,1:mp,0   )=2._euwp*wadv_a(1:np,1:mp,   1) &
                                       -wadv_a(1:np,1:mp,   2)
      IF(skyedge == 1) &
         wadv_a(1:np,1:mp,lp+1)=2._euwp*wadv_a(1:np,1:mp,lp  ) &
                                       -wadv_a(1:np,1:mp,lp-1)

      wadv_c(1:np,1:mp,1:lp)=0.5_euwp*(wadv_a(1:np,1:mp,1:lp  )    &
                                      +wadv_a(1:np,1:mp,0:lp-1))
      CALL update(wadv_c,np,mp,lp,np,mp,lp,1,ih)

      IF(skyedge == 1) THEN
         wadv_c(1:np,1:mp,lp+1)=0.5_euwp*(wadv_a(1:np,1:mp,lp+1)  &
                                         +wadv_a(1:np,1:mp,lp  ))
      ENDIF
!compute courant
      cr1(1)=maxval((abs(uadv_a(1:np,1:mp,1:lp))+                &
                     abs(vadv_a(1:np,1:mp,1:lp))+                &
                     abs(wadv_a(1:np,1:mp,1:lp)))/rho(1:np,1:mp,1:lp))
 
      CALL globmax(cr1,cr1_tot,1)
      IF(mype==0) PRINT *,'Max courant number is',cr1_tot
   END SUBROUTINE velprd_sphere


   ! Calls in subroutine TEST_SOLUTION84: 
   ! => globsum (on line <398>)
   SUBROUTINE test_solution84(x_actual,x_exact,icnt,experiment,np,mp,lp,ih)
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) ::         &
         x_actual,x_exact
      REAL(KIND=euwp),DIMENSION(np,mp,lp) ::l2norm_a,l2norm_e
      REAL(KIND=euwp),DIMENSION(2) ::l2normsum,l2normtot
      REAL(KIND=euwp) ::l2print
      INTEGER :: icnt
      CHARACTER :: experiment
      l2norm_a(1:np,1:mp,1:lp)=x_actual(1:np,1:mp,1:lp)**2
      l2norm_e(1:np,1:mp,1:lp)= x_exact(1:np,1:mp,1:lp)**2
      l2normsum(1)=sum(l2norm_a)
      l2normsum(2)=sum(l2norm_e)
      CALL globsum(l2normsum,l2normtot,2)
      l2print=(l2normtot(1)-l2normtot(2))/l2normtot(2)
! l2print=l2print*nmli*dti
      IF(mype == 0) PRINT *,'L2 norm:',experiment,icnt,l2print
   END SUBROUTINE test_solution84

   ! Calls in subroutine TEST_SOLUTION_ENERGY: 
   ! => globsum (on line <418>)
   SUBROUTINE test_solution_energy(x_actual,x_exact,icnt,experiment,np,mp,lp,ih)
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) ::         &
         x_actual,x_exact
      REAL(KIND=euwp),DIMENSION(:,:,:),ALLOCATABLE ::l2norm_a,l2norm_e
      REAL(KIND=euwp),DIMENSION(2) ::l2normsum,l2normtot
      REAL(KIND=euwp) ::l2print
      INTEGER :: icnt,ierr
      CHARACTER :: experiment
      ALLOCATE (l2norm_a(np,mp,lp), STAT=ierr); IF(ierr.ne.0) PRINT *,'Error allocating lnorm_a'
      ALLOCATE (l2norm_e(np,mp,lp), STAT=ierr); IF(ierr.ne.0) PRINT *,'Error allocating lnorm_e'
      l2norm_a(1:np,1:mp,1:lp)=x_actual(1:np,1:mp,1:lp)*x_actual(1:np,1:mp,1:lp)
      l2norm_e(1:np,1:mp,1:lp)= x_exact(1:np,1:mp,1:lp)* x_exact(1:np,1:mp,1:lp)
      l2normsum(1)=sum(l2norm_a)
      l2normsum(2)=sum(l2norm_e)
      CALL globsum(l2normsum,l2normtot,2)
      l2print=abs((l2normtot(1)-l2normtot(2))) !/l2normtot(2))
      l2print=sqrt(l2print*nmli)/timetot
      IF(mype == 0) PRINT 101,achar(10),experiment,icnt,l2print
101   FORMAT (a,'Energy L2 norm: ',a,' ',i5,' ',1(F10.8,2x))
   END SUBROUTINE test_solution_energy
   ! Calls in subroutine TEST_SOLUTION_ERR2_91: 
   ! => globsum (on line <441>)
   SUBROUTINE test_solution_ERR2_91(x_actual,x_exact,np,mp,lp,ih)
      USE geometry, ONLY : gac
      INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) ::         &
         x_actual,x_exact
      REAL(KIND=euwp),DIMENSION(np,mp,lp) ::l2norm_a,l2norm_e
      REAL(KIND=euwp),DIMENSION(2) ::l2normsum,l2normtot
      REAL(KIND=euwp) ::l2print
      l2norm_a(:,:,:)=0._euwp
      l2norm_e(:,:,:)=0._euwp
      l2norm_a(1:np,1:mp,1)=x_actual(1:np,1:mp,1)**2*gac(1:np,1:mp,1)
      l2norm_e(1:np,1:mp,1)= x_exact(1:np,1:mp,1)**2*gac(1:np,1:mp,1)
      l2normsum(1)=sum(l2norm_a)
      l2normsum(2)=sum(l2norm_e)
      CALL globsum(l2normsum,l2normtot,2)
      l2print=((l2normtot(1)-l2normtot(2)))/l2normtot(2)
       IF(mype == 0) print *,'ERR2:',l2print
   END SUBROUTINE test_solution_ERR2_91

   SUBROUTINE test_solution_ERR1_91(x_actual,x_exact,np,mp,lp,ih)
      USE geometry, ONLY : gac
        INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) ::         &
         x_actual,x_exact
      REAL(KIND=euwp),DIMENSION(np,mp,lp) ::l1norm_a,l1norm_e
      REAL(KIND=euwp),DIMENSION(2) ::l1normsum,l1normtot
      REAL(KIND=euwp) ::l1print
      l1norm_a(:,:,:)=0._euwp
      l1norm_e(:,:,:)=0._euwp
      l1norm_a(1:np,1:mp,1)=x_actual(1:np,1:mp,1)*gac(1:np,1:mp,1)
      l1norm_e(1:np,1:mp,1)= x_exact(1:np,1:mp,1)*gac(1:np,1:mp,1)
      l1normsum(1)=sum(l1norm_a)
      l1normsum(2)=sum(l1norm_e)
      CALL globsum(l1normsum,l1normtot,2)
      l1print=((l1normtot(1)-l1normtot(2)))/l1normtot(2)
       IF(mype == 0) print *,'ERR1: ',l1print
101   FORMAT (a,'ERR1 L1 norm: ',a,' ',i5,' ',F10.8)
   END SUBROUTINE test_solution_ERR1_91
    SUBROUTINE test_solution_ERR0_91(x_actual,x_exact,np,mp,lp,ih)
      USE geometry, ONLY : gac
        INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) ::         &
         x_actual,x_exact
      REAL(KIND=euwp),DIMENSION(np,mp,lp) ::l2norm_a,l2norm_e
      REAL(KIND=euwp),DIMENSION(1) ::l2normsum,l2normstot
      REAL(KIND=euwp),DIMENSION(1) ::l2normmax,l2normmtot
      REAL(KIND=euwp) ::l2print
      l2norm_a(:,:,:)=0._euwp
      l2norm_e(:,:,:)=0._euwp
      l2norm_a(1:np,1:mp,1)=sqrt((x_actual(1:np,1:mp,1) &
                                  -x_exact(1:np,1:mp,1))**2*gac(1:np,1:mp,1))
      l2norm_e(1:np,1:mp,1)= x_exact(1:np,1:mp,1) !*gac(1:np,1:mp,1)
      l2normsum(1)=sum(l2norm_a)
      CALL globsum(l2normsum,l2normstot,1)
      l2normmax(1)=maxval(l2norm_e)
      CALL globmax(l2normmax,l2normmtot,1)
      l2print=l2normstot(1)/l2normmtot(1)
! l2print=l2print*nmli/timetot
   END SUBROUTINE test_solution_ERR0_91
    SUBROUTINE test_solution_ERRMIN_91(x_actual,x_exact,np,mp,lp,ih)
        INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) ::         &
         x_actual,x_exact
      REAL(KIND=euwp),DIMENSION(np,mp,1) ::lnorm_a,lnorm_e
      REAL(KIND=euwp),DIMENSION(2) ::lnormmin,lnorm_min_tot
      REAL(KIND=euwp),DIMENSION(2) ::lnormmax,lnorm_max_tot
      REAL(KIND=euwp) ::lprint
      lnorm_a(:,:,:)=0._euwp
      lnorm_e(:,:,:)=0._euwp
      lnorm_a(1:np,1:mp,1) = x_actual(1:np,1:mp,1)  
      lnorm_e(1:np,1:mp,1) = x_exact(1:np,1:mp,1) !*gac(1:np,1:mp,1)
      lnormmin(1)=minval(lnorm_a)
      lnormmin(2)=minval(lnorm_e)
      CALL globmin(lnormmin,lnorm_min_tot,2)
      lnormmax(1)=maxval(lnorm_e)
      CALL globmax(lnormmax,lnorm_max_tot,1)
      lprint=(lnorm_min_tot(1)-lnorm_min_tot(2))/lnorm_max_tot(1)
      IF(mype == 0) print *,'ERRMIN: ',lprint
101   FORMAT (a,'ERRMIN norm: ',a,' ',i5,' ',F10.8)
    END SUBROUTINE test_solution_ERRMIN_91
    SUBROUTINE test_solution_ERRMAX_91(x_actual,x_exact,np,mp,lp,ih)
        INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) ::         &
         x_actual,x_exact
      REAL(KIND=euwp),DIMENSION(np,mp,1) ::lnorm_a,lnorm_e
      
      REAL(KIND=euwp),DIMENSION(2) ::lnormmax,lnorm_max_tot
      REAL(KIND=euwp) ::lprint
      lnorm_a(:,:,:)=0._euwp
      lnorm_e(:,:,:)=0._euwp
      lnorm_a(1:np,1:mp,1) = x_actual(1:np,1:mp,1)  
      lnorm_e(1:np,1:mp,1) = x_exact(1:np,1:mp,1) !*gac(1:np,1:mp,1)
      lnormmax(1)=maxval(lnorm_a)
      lnormmax(2)=maxval(lnorm_e)
      CALL globmax(lnormmax,lnorm_max_tot,2)
      lprint=(lnorm_max_tot(1)-lnorm_max_tot(2))/lnorm_max_tot(1)
      IF(mype == 0) print *,'ERRMAX: ',lprint
   END SUBROUTINE test_solution_ERRMAX_91
   ! Calls in subroutine TEST_SOLUTION_ERROR: 
   ! => globsum (on line <462>)
   SUBROUTINE test_solution_error(x_actual,x_exact,icnt,experiment,np,mp,lp,ih)
        INTEGER(KIND=iintegers) :: np,mp,lp,ih
      REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) ::         &
         x_actual,x_exact
!     REAL(KIND=euwp),DIMENSION(np,mp,lp) ::l2norm,l2norm_e
      REAL(KIND=euwp),DIMENSION(:,:,:),ALLOCATABLE ::l2norm,l2norm_e
      REAL(KIND=euwp),DIMENSION(2) ::l2normsum,l2normtot
      REAL(KIND=euwp) ::l2print
      INTEGER :: icnt
      CHARACTER :: experiment
      ALLOCATE(l2norm(np,mp,lp))
      ALLOCATE(l2norm_e(np,mp,lp))
      l2norm(1:np,1:mp,1:lp)=(x_exact(1:np,1:mp,1:lp)-x_actual(1:np,1:mp,1:lp))**2
      l2norm_e(1:np,1:mp,1:lp)= x_exact(1:np,1:mp,1:lp)**2
      l2normsum(1)=sum(l2norm)
      l2normsum(2)=sum(l2norm_e)
      CALL globsum(l2normsum,l2normtot,2)
      l2print=sqrt(l2normtot(1)*nmli)/timetot
      IF(mype == 0) PRINT 101,experiment,icnt,l2print
      DEALLOCATE(l2norm)
      DEALLOCATE(l2norm_e)
101   FORMAT ('Error  L2 norm: ',a,' ',i5,' ',1(F10.8,2x))
   END SUBROUTINE test_solution_error

       SUBROUTINE init_data(rhr,h,bcx,bcy,np,mp,lp,ih)
        INTEGER(KIND=iintegers) :: np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            h,rhr
!       REAL(KIND=euwp),DIMENSION(np,mp,lp) :: &
!           hi
        REAL(KIND=euwp) :: &
            bcx(mp,lp, 2), &               ! local array
            bcy(np,lp, 2)                  ! local array

        !Set fluid densities
        rhr(:,:,:)=1._euwp
        h(:,:,:)=1._euwp
!       hi(:,:,:)=1._euwp
        !Set boundary conditions for the advected field
        bcx(:,:,:)=0._euwp
        bcy(:,:,:)=0._euwp

    END SUBROUTINE init_data
       SUBROUTINE init_forces(xforc_impl,xforc_expl,np,mp,lp,ih)
        INTEGER(KIND=iintegers) :: np,mp,lp,ih
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            xforc_impl, xforc_expl 
        !Set forces
        xforc_impl(1:np,1:mp,1:lp)=0.001_euwp
        xforc_expl(1:np,1:mp,1:lp)=0.003_euwp
    END SUBROUTINE init_forces
       SUBROUTINE init_var_perturbed(noise_amp,noise_phase,var,bcx,bcy,np,mp,lp,ih)
        USE parameters, ONLY: pi2,n,m,l
        INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
        INTEGER(KIND=iintegers) :: i,j,k,ia,ja,ka 
        REAL(KIND=euwp) :: noise_amp,noise_phase
        REAL(KIND=euwp),PARAMETER :: bc_amp=0.1_euwp 
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
           var 
        REAL(KIND=euwp),DIMENSION(np,mp,lp) :: &
            noise
        REAL(KIND=euwp) :: &
            bcx(mp,lp, 2), &               ! local array
            bcy(np,lp, 2)                  ! local array
        REAL(KIND=euwp) :: &
            bcxnoise(mp,lp), &               ! local array
            bcynoise(np,lp)                  ! local array
!        print *,'noise',noise
!        print *,'nsubpos,msubpos,lsubpos',nsubpos,msubpos,lsubpos
        !Set fluid densities
        DO k=1,lp
          DO j=1,mp
            DO i=1,np
              ia=nsubpos+i
              ja=msubpos+j
              ka=lsubpos+k
               noise(i,j,k)=noise_amp*sin(pi2*ia/n+        noise_phase*pi2/3._euwp) &
                                     *cos(pi2*ja/m+2._euwp*noise_phase*pi2/3._euwp) &
                                     *cos(pi2*ka/l+3._euwp*noise_phase*pi2/3._euwp)
              var(i,j,k)=1._euwp+noise(i,j,k)
            ENDDO
          ENDDO
        ENDDO
        !Set boundary conditions for the advected field
!        print *,'n,m,l,pi2',n,m,l,pi2
        DO k=1,lp
          DO j=1,mp
              ja=msubpos+j
              ka=lsubpos+k
               IF(leftedge == 1) THEN
               bcxnoise(j,k)=noise_amp*cos(pi2*ja/m+2._euwp*noise_phase*pi2/7._euwp) &
                                      *sin(pi2*ka/l+4._euwp*noise_phase*pi2/7._euwp)
        bcx(j,k,1)=var(1,j,k)+bc_amp*bcxnoise(j,k)
               ENDIF
               IF(rightedge == 1) THEN
               bcxnoise(j,k)=noise_amp*sin(pi2*ja/m+3._euwp*noise_phase*pi2/7._euwp) &
                                      *cos(pi2*ka/l+5._euwp*noise_phase*pi2/7._euwp)
        bcx(j,k,2)=var(np,j,k)+bc_amp*bcxnoise(j,k)
               ENDIF
          ENDDO
        ENDDO
        DO k=1,lp
          DO i=1,np
              ia=nsubpos+i
              ka=lsubpos+k
               IF(botedge == 1) THEN
               bcynoise(i,k)=noise_amp*sin(pi2*ia/n+2._euwp*noise_phase*pi2/5._euwp) &
                                      *cos(pi2*ka/l+4._euwp*noise_phase*pi2/5._euwp)
               bcy(i,k,1)=var(i,1,k)+bc_amp*bcynoise(i,k)
               ENDIF
               IF(topedge == 1) THEN
               bcynoise(i,k)=noise_amp*cos(pi2*ia/n+3._euwp*noise_phase*pi2/5._euwp) &
                                      *sin(pi2*ka/l+5._euwp*noise_phase*pi2/5._euwp)
               bcy(i,k,2)=var(i,mp,k)+bc_amp*bcynoise(i,k)
               ENDIF
          ENDDO
        ENDDO
!        print *,'bcxnoise',bcxnoise
!        print *,'bcynoise',bcynoise

    END SUBROUTINE init_var_perturbed
       SUBROUTINE init_data_perturbed(noise_amp,noise_phase,rhr,h,hi,bcx,bcy,np,mp,lp,ih)
        USE parameters, ONLY: pi2,n,m,l
        USE mpi_parallel, ONLY: nsubpos,msubpos,lsubpos 
        INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
        INTEGER(KIND=iintegers) :: i,j,k,ia,ja,ka 
        REAL(KIND=euwp) :: noise_amp,noise_phase
        REAL(KIND=euwp),PARAMETER :: bc_amp=0.1_euwp 
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            h,rhr
        REAL(KIND=euwp),DIMENSION(np,mp,lp) :: &
            hi,noise
        REAL(KIND=euwp) :: &
            bcx(mp,lp, 2), &               ! local array
            bcy(np,lp, 2)                  ! local array
        REAL(KIND=euwp) :: &
            bcxnoise(mp,lp), &               ! local array
            bcynoise(np,lp)                  ! local array
        DO k=1,lp
          DO j=1,mp
            DO i=1,np
              ia=nsubpos+i
              ja=msubpos+j
              ka=lsubpos+k
               noise(i,j,k)=noise_amp*cos(pi2*ia/n+   noise_phase*pi2) &
                                     *cos(pi2*ja/m+2.*noise_phase*pi2) &
                                     *cos(pi2*ka/l+3.*noise_phase*pi2)
              rhr(i,j,k)=1._euwp+noise(i,j,k)
            ENDDO
          ENDDO
        ENDDO
!        print *,'noise',noise
!        print *,'nsubpos,msubpos,lsubpos',nsubpos,msubpos,lsubpos
        !Set fluid densities
        DO k=1,lp
          DO j=1,mp
            DO i=1,np
              ia=nsubpos+i
              ja=msubpos+j
              ka=lsubpos+k
               noise(i,j,k)=noise_amp*sin(pi2*ia/n+        noise_phase*pi2/3._euwp) &
                                     *cos(pi2*ja/m+2._euwp*noise_phase*pi2/3._euwp) &
                                     *cos(pi2*ka/l+3._euwp*noise_phase*pi2/3._euwp)
              h(i,j,k)=1._euwp+noise(i,j,k)
            ENDDO
          ENDDO
        ENDDO
        hi(1:np,1:mp,1:lp)=1._euwp/h(1:np,1:mp,1:lp)
        !Set boundary conditions for the advected field
!        print *,'n,m,l,pi2',n,m,l,pi2
        DO k=1,lp
          DO j=1,mp
              ja=msubpos+j
              ka=lsubpos+k
               bcxnoise(j,k)=noise_amp*cos(pi2*ja/m+2._euwp*noise_phase*pi2/7._euwp) &
                                      *sin(pi2*ka/l+4._euwp*noise_phase*pi2/7._euwp)
        bcx(j,k,1)=h(1,j,k)+bc_amp*bcxnoise(j,k)
               bcxnoise(j,k)=noise_amp*sin(pi2*ja/m+3._euwp*noise_phase*pi2/7._euwp) &
                                      *cos(pi2*ka/l+5._euwp*noise_phase*pi2/7._euwp)
        bcx(j,k,2)=h(np,j,k)+bc_amp*bcxnoise(j,k)
          ENDDO
        ENDDO
        DO k=1,lp
          DO i=1,np
              ia=nsubpos+i
              ka=lsubpos+k
               bcynoise(i,k)=noise_amp*sin(pi2*ia/n+2._euwp*noise_phase*pi2/5._euwp) &
                                      *cos(pi2*ka/l+4._euwp*noise_phase*pi2/5._euwp)
               bcy(i,k,1)=h(i,1,k)+bc_amp*bcynoise(i,k)
               bcynoise(i,k)=noise_amp*cos(pi2*ia/n+3._euwp*noise_phase*pi2/5._euwp) &
                                      *sin(pi2*ka/l+5._euwp*noise_phase*pi2/5._euwp)
               bcy(i,k,2)=h(i,mp,k)+bc_amp*bcynoise(i,k)
          ENDDO
        ENDDO
!        print *,'bcxnoise',bcxnoise
!        print *,'bcynoise',bcynoise

    END SUBROUTINE init_data_perturbed
    SUBROUTINE init_data_sphere(h,hi,np,mp,lp,ih)
        USE mpi_parallel, ONLY :update
        USE geometry, ONLY : gac
        INTEGER(KIND=iintegers) :: np,mp,lp,ih
        REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
            h, &
            rhr
        REAL_euwp,DIMENSION(np,mp,lp) :: &
            hi
        REAL(KIND=euwp) rh00

        !Set fluid densities
        rh00=1._euwp
        rhr(:,:,:)=1._euwp
        h(1:np,1:mp,1:lp)=rh00*gac(1:np,1:mp,1:lp)
        CALL update(h,np,mp,lp,np,mp,lp,1,ih)
        hi(1:np,1:mp,1:lp)=1._euwp/h(1:np,1:mp,1:lp)
    END SUBROUTINE init_data_sphere
    SUBROUTINE fill_test_variable_set(xtracer,u,v,w,th,tht,pstr,qv,qc,qr,qi,qs,qg,  &
                                       u_a_lr,v_a_lr,w_a_lr,th_a_lr,tht_a_lr,ex_a_lr, &
                                       u_a_bt,v_a_bt,w_a_bt,th_a_bt,tht_a_bt,ex_a_bt, &
                                       qv_a_lr, qv_a_bt, qc_a_lr, qc_a_bt, qr_a_lr, qr_a_bt, &
                                       qi_a_lr, qi_a_bt, qs_a_lr, qs_a_bt, qg_a_lr, qg_a_bt, &
                                      np,mp,lp,ih)
        INTEGER(KIND=iintegers) :: np,mp,lp,ih
        INTEGER(KIND=iintegers) :: i,j,k
        REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &
           xtracer
        REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(OUT) :: &
           u,v,w,th,tht,pstr,qv,qc,qr,qi,qs,qg 
        REAL_euwp,DIMENSION(mp,lp,2),INTENT(OUT) ::   &
            u_a_lr,v_a_lr,w_a_lr,th_a_lr,tht_a_lr,ex_a_lr, &
            u_a_bt,v_a_bt,w_a_bt,th_a_bt,tht_a_bt,ex_a_bt
        REAL_euwp,DIMENSION(np,lp,2),INTENT(OUT) ::  &
            qv_a_lr, qv_a_bt, qc_a_lr, qc_a_bt, qr_a_lr, qr_a_bt, &
            qi_a_lr, qi_a_bt, qs_a_lr, qs_a_bt, qg_a_lr, qg_a_bt


!$cuf kernel do(3) <<< (*,*), (32,4) >>>
   DO k=1,lp
     DO j=1,mp
       DO i=1,np
        th(i,j,k)=xtracer(i,j,k)
         u(i,j,k)=xtracer(i,j,k) 
         v(i,j,k)=xtracer(i,j,k)  
         w(i,j,k)=xtracer(i,j,k)  
       tht(i,j,k)=xtracer(i,j,k)
      pstr(i,j,k)=xtracer(i,j,k)
        qv(i,j,k)=xtracer(i,j,k)
        qc(i,j,k)=xtracer(i,j,k)
        qr(i,j,k)=xtracer(i,j,k)
        qi(i,j,k)=xtracer(i,j,k)
        qs(i,j,k)=xtracer(i,j,k)
        qg(i,j,k)=xtracer(i,j,k)
        ENDDO
      ENDDO
    ENDDO
!$cuf kernel do(2) <<< (*,*), (32,4) >>>
   DO k=1,lp
     DO j=1,mp
        th_a_lr(j,k,1)=0.
         u_a_lr(j,k,1)=0.
         v_a_lr(j,k,1)=0.
         w_a_lr(j,k,1)=0.
       tht_a_lr(j,k,1)=0.
        ex_a_lr(j,k,1)=0.
        qv_a_lr(j,k,1)=0.
        qc_a_lr(j,k,1)=0.
        qr_a_lr(j,k,1)=0.
        qi_a_lr(j,k,1)=0.
        qs_a_lr(j,k,1)=0.
        qg_a_lr(j,k,1)=0.
        th_a_lr(j,k,2)=0.
         u_a_lr(j,k,2)=0.
         v_a_lr(j,k,2)=0.
         w_a_lr(j,k,2)=0.
       tht_a_lr(j,k,2)=0.
        ex_a_lr(j,k,2)=0.
        qv_a_lr(j,k,2)=0.
        qc_a_lr(j,k,2)=0.
        qr_a_lr(j,k,2)=0.
        qi_a_lr(j,k,2)=0.
        qs_a_lr(j,k,2)=0.
        qg_a_lr(j,k,2)=0.
      ENDDO
    ENDDO
!$cuf kernel do(2) <<< (*,*), (32,4) >>>
   DO k=1,lp
       DO i=1,np
        th_a_bt(i,k,1)=0.
         u_a_bt(i,k,1)=0.
         v_a_bt(i,k,1)=0.
         w_a_bt(i,k,1)=0.
       tht_a_bt(i,k,1)=0.
        ex_a_bt(i,k,1)=0.
        qv_a_bt(i,k,1)=0.
        qc_a_bt(i,k,1)=0.
        qr_a_bt(i,k,1)=0.
        qi_a_bt(i,k,1)=0.
        qs_a_bt(i,k,1)=0.
        qg_a_bt(i,k,1)=0.
        th_a_bt(i,k,2)=0.
         u_a_bt(i,k,2)=0.
         v_a_bt(i,k,2)=0.
         w_a_bt(i,k,2)=0.
       tht_a_bt(i,k,2)=0.
        ex_a_bt(i,k,2)=0.
        qv_a_bt(i,k,2)=0.
        qc_a_bt(i,k,2)=0.
        qr_a_bt(i,k,2)=0.
        qi_a_bt(i,k,2)=0.
        qs_a_bt(i,k,2)=0.
        qg_a_bt(i,k,2)=0.
        ENDDO
    ENDDO


    END SUBROUTINE fill_test_variable_set 
    SUBROUTINE fill_full_variable_set(xtracer,u,v,w,th,tht,pstr,qv,qc,qr,qi,qs,qg,  &
                                       u_a_lr,v_a_lr,w_a_lr,th_a_lr,tht_a_lr,ex_a_lr, &
                                       u_a_bt,v_a_bt,w_a_bt,th_a_bt,tht_a_bt,ex_a_bt, &
                                       qv_a_lr, qv_a_bt, qc_a_lr, qc_a_bt, qr_a_lr, qr_a_bt, &
                                       qi_a_lr, qi_a_bt, qs_a_lr, qs_a_bt, qg_a_lr, qg_a_bt, &
                                      np,mp,lp,ih)
        INTEGER(KIND=iintegers) :: np,mp,lp,ih
        INTEGER(KIND=iintegers) :: i,j,k
        REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &
           xtracer
        REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(OUT) :: &
           u,v,w,th,tht,pstr,qv,qc,qr,qi,qs,qg 
        REAL_euwp,DIMENSION(mp,lp,2),INTENT(OUT) ::   &
            u_a_lr,v_a_lr,w_a_lr,th_a_lr,tht_a_lr,ex_a_lr, &
            u_a_bt,v_a_bt,w_a_bt,th_a_bt,tht_a_bt,ex_a_bt
        REAL_euwp,DIMENSION(np,lp,2),INTENT(OUT) ::  &
            qv_a_lr, qv_a_bt, qc_a_lr, qc_a_bt, qr_a_lr, qr_a_bt, &
            qi_a_lr, qi_a_bt, qs_a_lr, qs_a_bt, qg_a_lr, qg_a_bt


!$cuf kernel do(3) <<< (*,*), (32,4) >>>
   DO k=1,lp
     DO j=1,mp
       DO i=1,np
        th(i,j,k)=xtracer(i,j,k)
         u(i,j,k)=xtracer(i,j,k) 
         v(i,j,k)=xtracer(i,j,k)  
         w(i,j,k)=xtracer(i,j,k)  
       tht(i,j,k)=xtracer(i,j,k)
      pstr(i,j,k)=xtracer(i,j,k)
        qv(i,j,k)=xtracer(i,j,k)
        qc(i,j,k)=xtracer(i,j,k)
        qr(i,j,k)=xtracer(i,j,k)
        qi(i,j,k)=xtracer(i,j,k)
        qs(i,j,k)=xtracer(i,j,k)
        qg(i,j,k)=xtracer(i,j,k)
        ENDDO
      ENDDO
    ENDDO
!$cuf kernel do(2) <<< (*,*), (32,4) >>>
   DO k=1,lp
     DO j=1,mp
        th_a_lr(j,k,1)=0.
         u_a_lr(j,k,1)=0.
         v_a_lr(j,k,1)=0.
         w_a_lr(j,k,1)=0.
       tht_a_lr(j,k,1)=0.
        ex_a_lr(j,k,1)=0.
        qv_a_lr(j,k,1)=0.
        qc_a_lr(j,k,1)=0.
        qr_a_lr(j,k,1)=0.
        qi_a_lr(j,k,1)=0.
        qs_a_lr(j,k,1)=0.
        qg_a_lr(j,k,1)=0.
        th_a_lr(j,k,2)=0.
         u_a_lr(j,k,2)=0.
         v_a_lr(j,k,2)=0.
         w_a_lr(j,k,2)=0.
       tht_a_lr(j,k,2)=0.
        ex_a_lr(j,k,2)=0.
        qv_a_lr(j,k,2)=0.
        qc_a_lr(j,k,2)=0.
        qr_a_lr(j,k,2)=0.
        qi_a_lr(j,k,2)=0.
        qs_a_lr(j,k,2)=0.
        qg_a_lr(j,k,2)=0.
      ENDDO
    ENDDO
!$cuf kernel do(2) <<< (*,*), (32,4) >>>
   DO k=1,lp
       DO i=1,np
        th_a_bt(i,k,1)=0.
         u_a_bt(i,k,1)=0.
         v_a_bt(i,k,1)=0.
         w_a_bt(i,k,1)=0.
       tht_a_bt(i,k,1)=0.
        ex_a_bt(i,k,1)=0.
        qv_a_bt(i,k,1)=0.
        qc_a_bt(i,k,1)=0.
        qr_a_bt(i,k,1)=0.
        qi_a_bt(i,k,1)=0.
        qs_a_bt(i,k,1)=0.
        qg_a_bt(i,k,1)=0.
        th_a_bt(i,k,2)=0.
         u_a_bt(i,k,2)=0.
         v_a_bt(i,k,2)=0.
         w_a_bt(i,k,2)=0.
       tht_a_bt(i,k,2)=0.
        ex_a_bt(i,k,2)=0.
        qv_a_bt(i,k,2)=0.
        qc_a_bt(i,k,2)=0.
        qr_a_bt(i,k,2)=0.
        qi_a_bt(i,k,2)=0.
        qs_a_bt(i,k,2)=0.
        qg_a_bt(i,k,2)=0.
        ENDDO
    ENDDO


    END SUBROUTINE fill_full_variable_set 

    SUBROUTINE initialize_for_gpu(uadv,vadv,wadv,rhoadv,rhr,np,mp,lp,ih)
      USE mpi_parallel, ONLY: leftedge,rightedge, botedge, topedge,gndedge, skyedge
        INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
        INTEGER(KIND=iintegers) :: i,j,k
        REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(INOUT) :: &
          rhoadv,rhr 
        REAL_euwp,DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+ih  ),INTENT(INOUT)  :: uadv
        REAL_euwp,DIMENSION(1-ih:np+ih  , 1-ih:mp+1+ih, 1-ih:lp+ih  ),INTENT(INOUT)  :: vadv
        REAL_euwp,DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+1+ih),INTENT(INOUT)  :: wadv
        REAL(KIND=euwp) :: numberone
        numberone=1.
        CGRIDXFullYZDomainLoopDC(uadv(i,j,k)=uadv(i,j,k)*numberone;)
        CGRIDYFullXZDomainLoopDC(vadv(i,j,k)=vadv(i,j,k)*numberone;)
        CGRIDZFullXYDomainLoopDC(wadv(i,j,k)=wadv(i,j,k)*numberone;)
             FullXYZDomainLoopDC( rhr(i,j,k)= rhr(i,j,k)*numberone;)
             FullXYZDomainLoopDC(rhoadv(i,j,k)=rhoadv(i,j,k)*numberone;)
    END SUBROUTINE initialize_for_gpu 
   SUBROUTINE physical2contravariant(u_phys,v_phys,w_phys,ox_contr,oy_contr,oz_contr,np,mp,lp,ih)
      USE geometry, ONLY: g11,g12,g13,g21,g22,g23,g33
      INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
      INTEGER (KIND=IINTEGERS) :: i,j,k
   REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(IN) :: &

   u_phys,v_phys,w_phys
   REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),INTENT(OUT) :: &

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

END MODULE testing

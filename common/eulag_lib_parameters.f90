#include  "../advection/src_algorithms/renames.inc"
MODULE mod_parameters
#ifdef PUREMPI 
#ifdef MPI90  
   USE mpi
#endif
#endif
   USE precisions
   IMPLICIT NONE
#ifdef PUREMPI 
#ifdef  MPI77
  INCLUDE "mpif.h"
#endif
#endif


! Physical constants
   REAL(KIND=euwp), PROTECTED :: rds,rdsi  !rds - sphere radius (m)      
   REAL(KIND=euwp), PROTECTED :: fcr0 !fcr0 - coriolis parameter (s^-1)
   REAL(KIND=euwp), PARAMETER ::  ang = 0 !f=fcr0*sin(ang) for f- or beta-plane approximation 
   REAL(KIND=euwp), PARAMETER :: btpl = 0 ! beta-plane approximation on or off                  
   REAL(KIND=dp), PARAMETER :: pi=         acos(-1._dp)
   REAL(KIND=dp), PARAMETER :: pi2=2._dp*acos(-1._dp)
   REAL(KIND=dp), PARAMETER :: pih=.5_dp*acos(-1._dp)
   REAL(KIND=euwp), PARAMETER :: yj0=pi/180.*ang 
   REAL(KIND=euwp) spexi ! 1.0/spex=1/(th00*cp*dth)
   REAL(KIND=euwp), PROTECTED ::  capi ! cp_d/r_d
   REAL(KIND=euwp), PROTECTED ::  cap ! cp_d/r_d
   REAL(KIND=euwp), PROTECTED ::    rg ! Dry gas constant, COSMO: r_d
   REAL(KIND=euwp), PROTECTED :: cmpex ! r_d/pr00
   REAL(KIND=euwp), PROTECTED :: wexnr ! r_d/(cp_d-r_d) 
   REAL(KIND=euwp), PROTECTED :: epsi ! rv/r_d 
   REAL(KIND=euwp), PROTECTED :: epsb ! rv=461.51 rv/r_d - 1 
   REAL(KIND=euwp), PROTECTED ::    g
   REAL(KIND=euwp), PROTECTED ::   st 
   REAL(KIND=euwp), PROTECTED :: th00 
   REAL(KIND=euwp), PROTECTED :: rh00 
   REAL(KIND=euwp), PROTECTED :: tt00 
   REAL(KIND=euwp), PROTECTED :: pr00 
   REAL(KIND=euwp), PROTECTED :: cp ! cp_d 


! Numerical parameters that can be set constant
   INTEGER(KIND=iintegers), PARAMETER :: j3 = 1
   INTEGER(KIND=iintegers), PARAMETER :: msz=6 
   INTEGER(KIND=iintegers), PARAMETER :: mucnt=8 
   INTEGER(KIND=iintegers), PARAMETER :: ipoldiffmode=1 
   INTEGER(KIND=iintegers) :: ipoles=0 
   INTEGER(KIND=iintegers) :: icmpex=0 
   INTEGER(KIND=iintegers) :: itp0=600 
   INTEGER(KIND=iintegers) :: itp1=100 
   INTEGER(KIND=iintegers) :: itmn=1 
   REAL(KIND=euwp) :: epp0=1.e-6 
   REAL(KIND=euwp) :: epp1=1.e-5  

   INTEGER(KIND=iintegers) :: ibcx,ibcy,ibcz
   INTEGER(KIND=iintegers),PARAMETER :: glob_red_size=10
   INTEGER(KIND=iintegers) :: debugflag=0 
#if(STATICMEM ==0) 
   INTEGER(KIND=iintegers), PROTECTED  :: ih
   INTEGER(KIND=iintegers) :: n,m,l
   INTEGER(KIND=iintegers) :: np,mp,lp,nt
   INTEGER(KIND=iintegers)  :: nprocx, nprocy, nprocz
   INTEGER(KIND=iintegers) :: lrd
#else
   INTEGER(KIND=iintegers),PARAMETER :: nprocx=1, nprocy=1, nprocz=1
   INTEGER(KIND=iintegers),PARAMETER :: n=64,m=64,l=64
   INTEGER(KIND=iintegers),PARAMETER :: np=n/nprocx,mp=m/nprocy,lp=l/nprocz
   INTEGER(KIND=iintegers),PARAMETER :: ih=1
   INTEGER(KIND=iintegers),PARAMETER :: lrd=2
#endif
   INTEGER(KIND=iintegers), PROTECTED :: nm,nml
   INTEGER(KIND=iintegers) :: initprs
   INTEGER(KIND=iintegers), PROTECTED :: isphere,ideep
   INTEGER(KIND=iintegers) :: itr 
   INTEGER(KIND=iintegers) :: iab,iabth,iabexn,iabqw
   INTEGER(KIND=iintegers), PROTECTED :: irlx,irly
   INTEGER(KIND=iintegers), PROTECTED :: icmprss !p(ysical property:solver use 
   INTEGER(KIND=iintegers), PROTECTED :: noslip !physical property:solver use 
   INTEGER(KIND=iintegers), PROTECTED :: impl_metric !physical property:solver use 
   INTEGER(KIND=iintegers), PROTECTED :: implgw !physical property:solver use 
   INTEGER(KIND=iintegers), PROTECTED :: mtrord !numerical property:solver use 
   INTEGER(KIND=iintegers), PROTECTED :: mtrimx !numerical property:solver use 
   INTEGER(KIND=iintegers), PROTECTED :: genbuo !numerical property:solver use 
   INTEGER(KIND=iintegers), PROTECTED :: implsh !numerical property:solver use 


   REAL(KIND=euwp), PROTECTED :: nmi,nmli 
   REAL(KIND=euwp), PROTECTED :: dx00
   REAL(KIND=euwp), PROTECTED :: dy00
   REAL(KIND=euwp), PROTECTED :: dz00
   REAL(KIND=euwp), PROTECTED :: dx
   REAL(KIND=euwp), PROTECTED :: dy
   REAL(KIND=euwp), PROTECTED :: dz
   REAL(KIND=euwp), PROTECTED :: dt
   REAL(KIND=euwp), PROTECTED :: dxa
   REAL(KIND=euwp), PROTECTED :: dya
   REAL(KIND=euwp), PROTECTED :: dxh
   REAL(KIND=euwp), PROTECTED :: dyh
   REAL(KIND=euwp), PROTECTED :: dzh
   REAL(KIND=euwp), PROTECTED, PRIVATE :: dth
   REAL(KIND=euwp), PROTECTED, PRIVATE :: dti
   REAL(KIND=euwp), PROTECTED :: dxd
   REAL(KIND=euwp), PROTECTED :: dyd
   REAL(KIND=euwp), PROTECTED :: dzd
   REAL(KIND=euwp), PROTECTED :: dxi
   REAL(KIND=euwp), PROTECTED :: dyi
   REAL(KIND=euwp), PROTECTED :: dzi
   REAL(KIND=euwp), PROTECTED :: dxih
   REAL(KIND=euwp), PROTECTED :: dyih
   REAL(KIND=euwp), PROTECTED :: dzih
   REAL(KIND=euwp), PROTECTED :: dtih
   REAL(KIND=euwp), PROTECTED :: timetot
   REAL(KIND=euwp), PROTECTED :: zab  
   REAL(KIND=euwp), PROTECTED :: zb  
   REAL(KIND=euwp), PROTECTED:: towx
   REAL(KIND=euwp), PROTECTED:: towy
   REAL(KIND=euwp), PROTECTED:: towz
   REAL(KIND=euwp), PROTECTED:: towth_rat
   REAL(KIND=euwp), PROTECTED:: towexn_rat
   REAL(KIND=euwp), PROTECTED:: towqw_rat
   REAL(KIND=euwp), PROTECTED:: dxabL
   REAL(KIND=euwp), PROTECTED:: dyabL
   REAL(KIND=euwp), PROTECTED:: dxabR
   REAL(KIND=euwp), PROTECTED:: dyabR
#ifdef PNETCDF
   INTEGER(KIND=iintegers) :: npnetstep
#endif
   INTEGER(KIND=iintegers) :: nstarttest

CONTAINS

SUBROUTINE set_eulagcommon_mpi_parameters(     n_in, & 
                                          nprocx_in,nprocy_in,nprocz_in,  &
                                            ibcx_in,  ibcy_in,  ibcz_in,  & 
                                         isphere_in,ipoles_in)
   INTEGER(KIND=iintegers),INTENT(IN) ::  n_in
   INTEGER(KIND=iintegers),INTENT(IN) :: nprocx_in,nprocy_in,nprocz_in
   INTEGER(KIND=iintegers),INTENT(IN) ::  ibcx_in, ibcy_in, ibcz_in
   INTEGER(KIND=iintegers),INTENT(IN) ::  isphere_in, ipoles_in
#if(STATICMEM == 0)
     n=n_in
     nprocx=nprocx_in
     nprocy=nprocy_in
     nprocz=nprocz_in
     ibcx=ibcx_in
     ibcy=ibcy_in
     ibcz=ibcz_in
#endif
    isphere=isphere_in
    ipoles = ipoles_in
END SUBROUTINE set_eulagcommon_mpi_parameters

SUBROUTINE set_eulagcommon_lib_parameters(     n_in,     m_in,     l_in,  &
                                              np_in,    mp_in,    lp_in,  &
                                          nprocx_in,nprocy_in,nprocz_in,  &
                                            ibcx_in,  ibcy_in,  ibcz_in,  & 
                                            dx00_in,  dy00_in,  dz00_in,  & 
                                            towx_in,  towy_in,  towz_in,  &
                                           dxabL_in,dxabR_in,             &
                                           dyabL_in,dyabR_in,             &
                                       towth_rat_in,towexn_rat_in,towqw_rat_in, &
                                   iab_in, iabth_in,iabexn_in, iabqw_in,  &
                                           irelx_in, irely_in,   zab_in,  &
                                              ih_in,   lrd_in,noslip_in,  &
                                             itr_in,     impl_metric_in,  & 
                                         icmprss_in,isphere_in, ideep_in, &
                                                         dt_in,    nt_in)
   REAL(KIND=euwp),INTENT(IN) :: dt_in 
   REAL(KIND=euwp),INTENT(IN) ::  dx00_in, dy00_in, dz00_in
   REAL(KIND=euwp),INTENT(IN) ::  towx_in, towy_in, towz_in
   REAL(KIND=euwp),INTENT(IN) ::  towth_rat_in, towexn_rat_in, towqw_rat_in
   REAL(KIND=euwp),INTENT(IN) ::  zab_in
   REAL(KIND=euwp),INTENT(IN) ::  dxabL_in,dxabR_in,dyabL_in,dyabR_in 
   INTEGER(KIND=iintegers),INTENT(IN) ::  n_in, m_in, l_in
   INTEGER(KIND=iintegers),INTENT(IN) :: np_in,mp_in,lp_in,nt_in
   INTEGER(KIND=iintegers),INTENT(IN) :: nprocx_in,nprocy_in,nprocz_in
   INTEGER(KIND=iintegers),INTENT(IN) ::  ibcx_in, ibcy_in, ibcz_in
   INTEGER(KIND=iintegers),INTENT(IN) ::  iab_in, iabth_in, iabexn_in, iabqw_in 
   INTEGER(KIND=iintegers),INTENT(IN) ::  irelx_in, irely_in 
   INTEGER(KIND=iintegers),INTENT(IN) ::  ih_in, lrd_in, noslip_in
   INTEGER(KIND=iintegers),INTENT(IN) ::  itr_in, impl_metric_in
   INTEGER(KIND=iintegers),INTENT(IN) ::  icmprss_in ,isphere_in, ideep_in
   INTEGER iprint,mype,ierr,iprr
     iprint=0
#if(STATICMEM == 0)
     n=n_in
     m=m_in
     l=l_in
     np=np_in
     mp=mp_in
     lp=lp_in
     nprocx=nprocx_in
     nprocy=nprocy_in
     nprocz=nprocz_in
#if(1==0)
      IF (iprint==1) THEN
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifdef PUREMPI
       CALL MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)  ! total numbers of PE''s
#endif
      IF(mype == 0)     PRINT 98,nprocx,nprocy,nprocz
      CALL flush(6)
         DO iprr=0,nprocx*nprocy*nprocz-1
            IF (mype==iprr) THEN
!               PRINT 99,mype,np,mp,lp,npos,mpos,lpos
               CALL flush(6)
            ENDIF
#ifdef PUREMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

         ENDDO
98       FORMAT('nprocx = ',i3,'  nprocy = ',i3,'  nprocz = ',i3/,         &
         &' mype  np  mp lp npos mpos lpos')
99       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3)
      ENDIF !iprint
#endif

     ibcx=ibcx_in
     ibcy=ibcy_in
     ibcz=ibcz_in
         ih= ih_in
        lrd=lrd_in
#endif
        itr=itr_in
 impl_metric=impl_metric_in
    isphere=isphere_in
    icmprss=icmprss_in
    ideep=ideep_in
   noslip=noslip_in
      iab=iab_in
    iabth=iabth_in
   iabexn=iabexn_in
     towx=towx_in
     towy=towy_in
     towz=towz_in
     dxabL=dxabL_in
     dxabR=dxabR_in
     dyabL=dyabL_in
     dyabR=dyabR_in
     towth_rat=towth_rat_in
     towexn_rat=towexn_rat_in
     towqw_rat=towqw_rat_in
     towy=towy_in
     towz=towz_in
    iabqw=iabqw_in
     irlx=irelx_in
     irly=irely_in
      zab=zab_in
       zb=dz00_in
  mtrord = 2_iintegers
  mtrimx = 2_iintegers
!  IF(impl_metric == 1)  mtrimx=1_iintegers
  IF (icmprss == 0) THEN
    mtrord=mtrord*(isphere) +(1-isphere)
   IF(mtrord == 1 .OR. mtrord == 3) mtrimx=1_iintegers
  ENDIF
     implgw = 1_iintegers
     nml=n*m*l
     nm=n*m
     nmli=1._euwp/nml
     nmi=1._euwp/nm
     nt=nt_in
     dt=dt_in
     dx00=dx00_in
     dy00=dy00_in
     dz00=dz00_in
     dz=dz00/(l-1)
     if(isphere.eq.1) then
       if(ipoles.eq.1) then
         dxa=pi2/dble(n)
         dya= pi/dble(m)
       else 
         dxa=(360._euwp/180._euwp)*pi/dble(n-1)   ! full zonal extent :radians
         dya=(180._euwp/180._euwp)*pi/dble(m)     ! full meridional extent
       endif
       dx=rds*dxa                      ! meters
       dy=rds*dya
     else !isphere
        dx=dx00/(n-1)
        dy=dy00/(m-1)
       dxa=dx
       dya=dy
     endif
!    dxa= pi/dble(n/2.)
!    dya= pi/dble(m)
!special
!    dxa=0.01*(pi/180.)
!    dya=0.01*(pi/180.) 
!    dx=rds*dxa                      ! meters
!    dy=rds*dya
!end of special
     dxi=1._euwp/dx
     dyi=1._euwp/dy
     dzi=1._euwp/dz
     dti=1._euwp/dt
     dxih=.5_euwp*dxi
     dyih=.5_euwp*dyi
     dzih=.5_euwp*dzi
     dxd=2._euwp*dx
     dyd=2._euwp*dy
     dzd=2._euwp*dz
     dth=.5*dt
#ifdef PNETCDF
      npnetstep  = nt-1
#endif
      nstarttest = nt-1
      timetot=nt*dt
END SUBROUTINE set_eulagcommon_lib_parameters

SUBROUTINE set_eulagcommon_lib_fracsphere_parameters(  &
                                               n_in,     m_in,     l_in,  &
                                              np_in,    mp_in,    lp_in,  &
                                          nprocx_in,nprocy_in,nprocz_in,  &
                                            ibcx_in,  ibcy_in,  ibcz_in,  & 
                                            dx00_in,  dy00_in,  dz00_in,  & 
                                            towx_in,  towy_in,  towz_in,  &
                                           dxabL_in,dxabR_in,             &
                                           dyabL_in,dyabR_in,             &
                                       towth_rat_in,towexn_rat_in,towqw_rat_in, &
                                   iab_in, iabth_in,iabexn_in, iabqw_in,  &
                                           irelx_in, irely_in,   zab_in,  &
                                              ih_in,   lrd_in,noslip_in,  &
                                             itr_in,     impl_metric_in,  & 
                                         icmprss_in,isphere_in, ideep_in, &
                                                         dt_in,    nt_in)
   REAL(KIND=euwp),INTENT(IN) :: dt_in 
   REAL(KIND=euwp),INTENT(IN) ::  dx00_in, dy00_in, dz00_in
   REAL(KIND=euwp),INTENT(IN) ::  towx_in, towy_in, towz_in
   REAL(KIND=euwp),INTENT(IN) ::  towth_rat_in, towexn_rat_in, towqw_rat_in
   REAL(KIND=euwp),INTENT(IN) ::  zab_in
   REAL(KIND=euwp),INTENT(IN) ::  dxabL_in,dxabR_in,dyabL_in,dyabR_in 
   INTEGER(KIND=iintegers),INTENT(IN) ::  n_in, m_in, l_in
   INTEGER(KIND=iintegers),INTENT(IN) :: np_in,mp_in,lp_in,nt_in
   INTEGER(KIND=iintegers),INTENT(IN) :: nprocx_in,nprocy_in,nprocz_in
   INTEGER(KIND=iintegers),INTENT(IN) ::  ibcx_in, ibcy_in, ibcz_in
   INTEGER(KIND=iintegers),INTENT(IN) ::  iab_in, iabth_in, iabexn_in, iabqw_in 
   INTEGER(KIND=iintegers),INTENT(IN) ::  irelx_in, irely_in 
   INTEGER(KIND=iintegers),INTENT(IN) ::  ih_in, lrd_in, noslip_in
   INTEGER(KIND=iintegers),INTENT(IN) ::  itr_in, impl_metric_in
   INTEGER(KIND=iintegers),INTENT(IN) ::  icmprss_in ,isphere_in, ideep_in
   INTEGER iprint,mype,ierr,iprr
     iprint=0
#if(STATICMEM == 0)
     n=n_in
     m=m_in
     l=l_in
     np=np_in
     mp=mp_in
     lp=lp_in
     nprocx=nprocx_in
     nprocy=nprocy_in
     nprocz=nprocz_in
#if(1==0)
      IF (iprint==1) THEN
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifdef PUREMPI
       CALL MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)  ! total numbers of PE''s
#endif
      IF(mype == 0)     PRINT 98,nprocx,nprocy,nprocz
      CALL flush(6)
         DO iprr=0,nprocx*nprocy*nprocz-1
            IF (mype==iprr) THEN
!               PRINT 99,mype,np,mp,lp,npos,mpos,lpos
               CALL flush(6)
            ENDIF
#ifdef PUREMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

         ENDDO
98       FORMAT('nprocx = ',i3,'  nprocy = ',i3,'  nprocz = ',i3/,         &
         &' mype  np  mp lp npos mpos lpos')
99       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3)
      ENDIF !iprint
#endif

     ibcx=ibcx_in
     ibcy=ibcy_in
     ibcz=ibcz_in
         ih= ih_in
        lrd=lrd_in
#endif
        itr=itr_in
 impl_metric=impl_metric_in
    isphere=isphere_in
    icmprss=icmprss_in
    ideep=ideep_in
   noslip=noslip_in
      iab=iab_in
    iabth=iabth_in
   iabexn=iabexn_in
     towx=towx_in
     towy=towy_in
     towz=towz_in
     dxabL=dxabL_in
     dxabR=dxabR_in
     dyabL=dyabL_in
     dyabR=dyabR_in
     towth_rat=towth_rat_in
     towexn_rat=towexn_rat_in
     towqw_rat=towqw_rat_in
     towy=towy_in
     towz=towz_in
    iabqw=iabqw_in
     irlx=irelx_in
     irly=irely_in
      zab=zab_in
       zb=dz00_in
  mtrord = 2_iintegers
  mtrimx = 2_iintegers
!  IF(impl_metric == 1)  mtrimx=1_iintegers
  IF (icmprss == 0) THEN
    mtrord=mtrord*(isphere) +(1-isphere)
   IF(mtrord == 1 .OR. mtrord == 3) mtrimx=1_iintegers
  ENDIF
     implgw = 1_iintegers
     nml=n*m*l
     nm=n*m
     nmli=1._euwp/nml
     nmi=1._euwp/nm
     nt=nt_in
     dt=dt_in
     dx00=dx00_in
     dy00=dy00_in
     dz00=dz00_in
     dz=dz00/(l-1)
     dx=dx00/(n-1)
     dy=dy00/(m-1)
!    dxa= pi/dble(n/2.)
!    dya= pi/dble(m)
!special
!    dxa=0.01*(pi/180.)
!    dya=0.01*(pi/180.) 
!    dx=rds*dxa                      ! meters
!    dy=rds*dya
!end of special
     dxi=1._euwp/dx
     dyi=1._euwp/dy
     dzi=1._euwp/dz
     dti=1._euwp/dt
     dxih=.5_euwp*dxi
     dyih=.5_euwp*dyi
     dzih=.5_euwp*dzi
     dxd=2._euwp*dx
     dyd=2._euwp*dy
     dzd=2._euwp*dz
     dth=.5*dt
#ifdef PNETCDF
      npnetstep  = nt-1
#endif
      nstarttest = nt-1
      timetot=nt*dt
END SUBROUTINE set_eulagcommon_lib_fracsphere_parameters

SUBROUTINE set_eulagcommon_lib_constants(fcr0_in,rds_in,spexi_in,rg_in,cp_in,         &
                                         epsi_in,g_in,th00_in,tt00_in,rh00_in,st_in)
   REAL(KIND=euwp),INTENT(IN) :: fcr0_in,rds_in,spexi_in,rg_in,  &
                                epsi_in,g_in,cp_in,th00_in,tt00_in,rh00_in,st_in 
!-----------------------------------------------------------------c
!Coriolis force specification                                     c
!-----------------------------------------------------------------c
!---> fcr0 - coriolis parameter (s^-1)                            c
!---> ang  - f=fcr0*sin(ang) for f- or beta-plane approximation   c
!---> btpl -  beta-plane approximation on or off                  c
!-----------------------------------------------------------------c
fcr0=fcr0_in
rds = rds_in
rdsi=1._euwp/rds
spexi = spexi_in
   rg =    rg_in
   cp =    cp_in
  cap = rg/cp
  capi= 1._euwp/cap
 epsi =  epsi_in
    g =     g_in
 epsb =  epsi - 1._euwp
 th00 =  th00_in
 tt00 =  tt00_in
 rh00 =  rh00_in
  pr00=rg*rh00*tt00
 cmpex=rg/pr00
   st =  st_in
 wexnr=rg/(cp-rg)!Rd/Cv
 cmpex=rg/pr00
END SUBROUTINE set_eulagcommon_lib_constants

SUBROUTINE set_eulagcommon_lib_domainsize(   n_in,    m_in,    l_in, ih_in, &
                                          dx00_in, dy00_in, dz00_in, nt_in, dt_in)  
   INTEGER(KIND=iintegers),INTENT(IN) ::  n_in, m_in, l_in, ih_in, nt_in
   REAL(KIND=euwp),INTENT(IN) ::  dx00_in, dy00_in, dz00_in, dt_in
#if(STATICMEM == 0)
    ih=ih_in
     n=n_in
     m=m_in
     l=l_in
     nml=n*m*l
     nm=n*m
     nmli=1._euwp/nml
     nmi=1._euwp/nm
#endif
     dx00=dx00_in
     dy00=dy00_in
     dz00=dz00_in
     dz=dz00/(l-1)
     dx=dx00/(n-1)
     dy=dy00/(m-1)
     dt=dt_in
     nt=nt_in
     timetot=nt*dt
END SUBROUTINE set_eulagcommon_lib_domainsize
        
#if(1==0)
SUBROUTINE header()
      real eps,ri00
      character(80) title 

! ---- computational parameters
      write (6,999)
      write (6,700) EXPER,FASTMPDATA,MPDATAMODE  
      write (6,701) EUCOMPILERE  
      write (6,702) COMPOPTIE  
      write (6,999)
      write (6,8992)
      write (6,999)
      write (6,8992)
         if(isphere.eq.1) then
#if (POLES == 0)
      write(6,*) 'Poles = 0'
#else
      write(6,*) 'Poles = 1'
      write(6,*) 'Boundary condition ipoldiffmode =',ipoldiffmode
#endif
      write (6,9000) (-90.+(m-0.5)*180./float(m)),rds/1000.
         end if
      write (6,901) n,m,L,lagr,ior
      write (6,914) nprocx,nprocy,nprocz
      write (6,902) dx,dy,dz,dt
      write (6,903) nt,nplot,nstore
      write (6,904) ibcx,ibcy,ibcz,irlx,irly
 !     write (6,9041) mabsmooth,mabsamb
!      write (6,9042) msmoother,nsmiter
      write (6,905) iab,iabth,iabqw
      write (6,906) zab,towz
      if(isphere.eq.1) then
          write (6,908) towxL,towxR,towy,dxabL/(rds*3.1416)*180.
     &                                  ,dxabR/(rds*3.1416)*180.
     &                                  ,dyab /(rds*3.1416)*180.
        else
          write (6,9080) towxL,towxR,towy,dxabL/1000.
     &                                   ,dxabR/1000.
     &                                   ,dyab/1000.
        end if


! ---- basic state
      write (6,999)
      write (6,8994)
      eps=1.e-15
      if((abs(u0z).lt.eps).and.(abs(v0z).le.eps)) then
        write (6,929) u00,v00
      else
        write (6,909) u00,u0z
      endif
      bv=sqrt(st*g)
      write (6,910) bv,lipps
      if((abs(u0z).lt.eps).and.(abs(v0z).le.eps)) then
        ri00=0.
        write (6,899)
        write (title,899)
      else
        ri00 = g*st/(u0z**2 + v0z**2)
        write (6,900) ri00
      endif

! ---- basic state
      write (6,999)
      write (6,8994)
      write (6,8999) th00,tt00,rh00,pr00
      eps=1.e-15
      if((abs(u0z).lt.eps).and.(abs(v0z).le.eps)) then
        write (6,929) u00,v00
      else
        write (6,909) u00,u0z
      endif
      bv=sqrt(st*g)
      write (6,910) bv,lipps
      if((abs(u0z).lt.eps).and.(abs(v0z).le.eps)) then
        ri00=0.
        write (6,899)
!       write (title,899)
      else
        ri00 = g*st/(u0z**2 + v0z**2)
        write (6,900) ri00
      endif

! ---- physical parameters
      write (6,999)
      write (6,8992)
      write (6,898) icmprss, itraj
      write (6,897) icont0, itras
      write (6,896) ipsinc, implgw
      write (6,892) implicit_metric, implsh, genbuo
      write (6,893) predict_pressure,lhelm2, prsoff
      write (6,894) mtrord,mtrimx
      write(6,980) ipreflag,itr,line 
      write(6,981) inihydrobal 
      write(6,907) epp1
      write (6,999)
      write (6,8996)
      if(isphere.eq.1) write (6,9001) rds/1000.
      write (6,912) icorio,ivis,itke
      write (6,913) moist,ice

! ---- transformation parameters
      write (6,999)
      if(timeadapt.eq.1) write (6,8995)
      if(mshadapt.eq.1.and.tstart.le.tt.and.tt.le.tend) then
           write (6,8998)
           write (6,9031) stime,tt,tend
           write (6,915) xb0i,yb0i,1./sxinvi,1./syinvi
           write (6,916) xb0f,yb0f,1./sxinvf,1./syinvf

      else
           write (6,8997)
           write (6,915) xb0i,yb0i,1./sxinvi,1./syinvi
      endif
      if(istr.eq.1) write (6,917) zhz/1000.
      write (6,911) xml,yml,amp
!     write (6,9031) stime,tt,tend
  700 format(1x,'Code optimization(EXPER, FASTMPDATA, MPDATAMODE) = ',1x
     .      ,3i2)
  701 format(1x,'Compiler choice (EUCOMPILER)          = ',i1)
  702 format(1x,'Compiler optimization level (COMPOPTI)= ',i2)
  892 format(1x,'impl_metr, implsh, genbuo  =',3i5)
  893 format(1x,'predict_pressure, lhelm2, prsoff',i5,1x,l,f8.3) 
  894 format(1x,'mtrord,mtrimx ',2i5) 
  898 format(1x,'icmprss, itraj   =',2i5)
  897 format(1x,'icont0,  itras   =',2i5)
  896 format(1x,'ipsinc,  implgw  =',2i5)
  899 format(1x,'Ri0 = infinity')
 8992 format(11x,'COMPUTATIONAL PARAMETERS')
 8994 format(11x,'BASIC STATE')
 8996 format(11x,'PHYSICAL PARAMETERS')
 8991 format(11x,'EQUATIONS/SOLVER PARAMETERS')
 8995 format(11x,'TIME ADAPTIVITY: active')
 8997 format(11x,'TRANSFORMATION PARAMETERS: stationary')
 8998 format(11x,'TRANSFORMATION PARAMETERS: adapting')
 8999 format(1x,'th00,tt00        = ',2(1x,e11.4)/
     &       1x,'rh00,pr00        = ',1x,f8.3,e12.4)
  900 format(1x,'Ri0 =',e11.4)
 9000 format(1x,'SPHERE: meridional grid = +/-',f8.3,' degrees'/
     &       11x,'radius = ',f8.1,' kms')
 9001 format(1x,'SPHERE:   radius = ',f8.1,' kms')
  901 format(1x,'n,m,l            =',3i5,5x,'lagr,ior =',2i5)
  902 format(1x,'dx,dy,dz,dt      =',4e11.4)
  903 format(1x,'nt,nplot,nstore  =',3i10)
 9031 format(1x,'time,tt,tend =',3e11.4)
  904 format(1x,'ibcx,ibcy,ibcz   =',3i5/
     &       1x,'irlx,irly        =',2i5)
 9041 format(1x,'mabsmooth,mabsamb=',2i5)
 9042 format(1x,'msmoother,nsmiter=',2i5)
  905 format(1x,'iab,iabth,iabqw  =',3i5)
  906 format(1x,'zab,towz         =',2e11.4)
  907 format(1x,'epp1             =',e11.4)
  908 format(1x,'towxL,towxR,towy =',3e11.4,/
     &       1x,'dxabL,dxabr,dyab =',3e11.4,' deg.')
  980 format(1x,'ipreflag,itr,line=',3i5)
  981 format(1x,'Hydro solver initialization type =',i5)
 9080 format(1x,'towxL,towxR,towy =',3e11.4/
     &       1x,'dxabL,dxabR,dyab        =',3e11.4,' kms')
 9081 format(1x,'towyp,dyabp,wmult=',e11.4,2x,f5.2,
     &      ' deg,  ',f5.2,'  (Ex.Polar Abs.)')
  909 format(1x,'Const shear profile:'/
     &       1x,'  U00,U00Z       =',2e11.4)
  915 format(1x,'initial horizonal grid parameters:'/
     &       1x,'  xb0i,yb0i      =',2e11.4/
     &       1x,'   sxi,syi       =',2e11.4)
  916 format(1x,'final horizonal grid parameters:'/
     &       1x,'  xb0f,yb0f      =',2e11.4/
     &       1x,'   sxf,syf       =',2e11.4)
  917 format(1x,'vertical grid stretching: scale ht = ',f6.2,' kms')
  929 format(1x,'Const wind profile:'/
     &       1x,'  U00,V00        =',2e11.4)
  910 format(1x,'Const stability profile:'/
     &       1x,'  N              =',1x,e11.4,' lipps =',i3)
  911 format(1x,'mountain scales:'/
     &       1x,'  Lx,Ly,h0       =',3e11.4)
  912 format(1x,'icorio,ivis,itke =',3i5)
  913 format(1x,'moist,ice        =',2i5)
  914 format(1x,'npr_x,npr_y,npr_z=',3i5)
  999 format(1x,' ****************** ')
      write (6,999)

      return
    
      end subroutine header
#endif /*(1==0)*/
END MODULE mod_parameters

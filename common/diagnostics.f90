#undef PNETCDF
#include "../advection/src_algorithms/renames.inc"
MODULE eulag_diagnostics 
USE precisions
USE mod_parameters, ONLY: spexi
USE mpi_parallel, ONLY: ttbeg,ttend
#ifdef PNETCDF
  USE mpi_parallel, ONLY: pnet_out_chunk
#endif

IMPLICIT NONE

CONTAINS
#include "../advection/src_algorithms/defines.inc"
#if(1==0)
!----------------------------------------------------------------------!
SUBROUTINE diagnos(u,v,w,ox,oy,oz,th,p,rho,rhoi,rh)
!----------------------------------------------------------------------!

  USE src_parallel, ONLY: update ,globsum,mybarrier,globmin, &
    globmax,globmaxv,globsumv,globsumaxminv


  !----------------------------------------------------------------------!

  REAL(KIND=euwp),INTENT(IN) :: &
     DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
    u, &
    v, &
    w, &
    ox, &
    oy, &
    oz, &
    th, &
    p, &
    x0, &
    y0, &
    z0, &
    rho, &
    rh, &
    div, &
    rhm, &
    tempp
  
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih) :: &
    psurf
   REAL(KIND=euwp),   INTENT(IN) ::                             &
                      ue(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                      ve(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                     the(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    pexe(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    uref(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                    vref(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                   thref(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),     &
                 pextref(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

  
  REAL(KIND=euwp), INTENT(IN) ::  thsum0, qws0, drgx, drgy     
  
  REAL(KIND=euwp) :: ptmp(1:np,1:mp,1:lp)  
  REAL(KIND=euwp),DIMENSION(np,mp,lp) ::  rhoi
  
  INTEGER(KIND=iintegers),SAVE :: ifirst = 0
  
  REAL (KIND=euwp) :: tsm0, xnorm, xnormj, divmx, divmn, divav, divsd, &
    deldf, ckmmx, ckmmn, ckmav, ckmsd, ommx, ommn, ummx, ummn, vmmx, vmmn, &
    tsm
  
  REAL (KIND=euwp) :: tmass, flinf, flout, dmass
  
  INTEGER(KIND=iintegers) :: igri3, illim, jllim, i, j, k, nitav, mitav, &
    ia, ja, irimin, jrimin, krimin, ibound, imass, iproc
  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2011 = &
  "(1x,' umax,  umin,  uave:',1X,3(e14.7,1x),"// &
  & "   1x,' vmax,  vmin,  vave:',1X,3(e14.7,1x),"// &
  & "   1x,' wmax,  wmin,  wave:',1X,3(e14.7,1x),"// &
  & "   1x,'thpmax,thpmin,thpave:',1X,3(e14.7,1x),"// &
  & "   1x,'thmax, thmin, thave:',1X,3(e14.7,1x),"// &
  & "   1x,'oxmax, oxmin, oxave:',1X,3(e14.7,1x),"// &
  & "   1x,'oymax, oymin, oyave:',1X,3(e14.7,1x),"// &
  & "   1x,'ozmax, ozmin, ozave:',1X,3(e14.7,1x),"// &
  & "   1x,' pmax,  pmin,  pave:',1X,3(e14.7,1x),"// &
  & "   1x,'psmax, psmin, psave:',1X,3(e14.7,1x))"
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
  "(1x,' uemax,  uemin,  ueave:',1X,3(e14.7,1x),"// &
  & "   1x,' vemax,  vemin,  veave:',1X,3(e14.7,1x),"// &
  & "   1x,'fcr2mx, fcr2mn, fcr2av:',1X,3(e14.7,1x),"// &
  & "   1x,'fcr3mx, fcr3mn, fcr3av:',1X,3(e14.7,1x),"// &
  & "   1x,'themax, themin, theave:',1X,3(e14.7,1x),"// &
  & "   1x,'ppemax, ppemin, ppeave:',1X,3(e14.7,1x),"// &
  & "   1x,' zsmax,  zsmin,  zsave:',1X,3(e14.7,1x),"// &
  & "   1x,' zhmax,  zhmin,  zhave:',1X,3(e14.7,1x),"// &
  & "   1x,'xcrmax, xcrmin, xcrave:',1X,3(e14.7,1x),"// &
  & "   1x,'ycrmax, ycrmin, ycrave:',1X,3(e14.7,1x),"// &
  & "   1x,'zcrmax, zcrmin, zcrave:',1X,3(e14.7,1x),"// &
  & "   1x,'gmmmax, gmmmin, gmmave:',1X,3(e14.7,1x),"// &
  & "   1x,'g11max, g11min, g11ave:',1X,3(e14.7,1x),"// &
  & "   1x,'g12max, g12min, g12ave:',1X,3(e14.7,1x),"// &
  & "   1x,'g13max, g13min, g13ave:',1X,3(e14.7,1x),"// &
  & "   1x,'g21max, g21min, g21ave:',1X,3(e14.7,1x),"// &
  & "   1x,'g22max, g22min, g22ave:',1X,3(e14.7,1x),"// &
  & "   1x,'g23max, g23min, g23ave:',1X,3(e14.7,1x),"// &
  & "   1x,'g33max, g33min, g33ave:',1X,3(e14.7,1x))"
  CHARACTER(LEN=*),PARAMETER :: fm2013 = &
  "(1x,'pextmax, pextmin, pextave:',1X,3(e14.7,1x),"// &
  &"    1x,' rhomax,  rhomin,  rhoave:',1X,3(e14.7,1x))"
  CHARACTER(LEN=*),PARAMETER :: fm2014 = &
  "(1x,' thtmax,  thtmin,  thtave:',1X,3(e14.7,1x))"
  
  CHARACTER(LEN=*),PARAMETER :: fm201 = &
  "(1x,'umax, umin:',I6,1X,2(e23.16,1x),/ 1x, 'vmax, vmin:',I6,1X,2(e23.16,1x),"// &
  &"1x,'wmax, wmin:',I6,1X,2(e23.16,1x),/ 1x, 'pmax, pmin:',I6,1X,2(e23.16,1x),"// &
  &"1x,'thmx, thmn:',I6,1X,2(e23.16,1x),2x,"// &
  &"1x,'zsmx, zsmn:',2e23.16,"// &
  &"1x,'zhmx, zhmn:',2e23.16,1x,/ 1x,'oxmx, oxmn:',2e23.16,"// &
  &"1x,'oymx, oymn:',2e23.16,1x,/ 1x,'ozmx, ozmn:',2e23.16)"
  
  CHARACTER(LEN=*),PARAMETER :: fm203 = &
  "(1x,'jcmx,jcmn,jcav,jcsd:',4f11.8)"
  
  
  !  CHARACTER(LEN=*),PARAMETER :: fm2022 = &
  !    "(1x,'ipsim=1, dftmx, dftmn, dftav, dftsd:',4e11.4)"
  
  CHARACTER(LEN=*),PARAMETER :: fm205 = &
  "(1x,'dvmx,dvmn:',2X,2(e20.13,1x),"//&
  &"1x,'dvav,dvsd:',2X,2(e20.13,1x),"//&
  &"1x,' eer, eem:',2X,2(e20.13,1x),"//&
  &"1x,'niter,nitav,miter,mitav:',4i4)"
  
  CHARACTER(LEN=*),PARAMETER :: fm207 = &
  "(1x,'k,ommx,ommn:',i5,2e11.4)"
  
  CHARACTER(LEN=*),PARAMETER :: fm208 = &
  "(1x,'i,ummx,ummn:',i5,2e11.4)"
  
  CHARACTER(LEN=*),PARAMETER :: fm209 = &
  "(1x,'j,vmmx,vmmn:',i5,2e11.4)"
  
  CHARACTER(LEN=*),PARAMETER :: fm206 = &
  "(1x,'uinf, vinf, oinf:',3e11.4,"// &
  &"1x,'uout, vout, oout:',3e11.4,"// &
  &"1x,'tflx, finf, fout:',3e11.4,"// &
  &"1x,'dmas, epsm, epsa:',3e11.4)"
  
  CHARACTER(LEN=*),PARAMETER :: fm2062 = &
  "(1x,'drgx, drgy, drgnorm:',3e11.4)"
  
  CHARACTER(LEN=*),PARAMETER :: fm317 = &
  "(1x,'max, min surf. precip (mm/day):',2e11.4)"
  
  CHARACTER(LEN=*),PARAMETER :: fm400 = &
  "(e11.5,2x,e11.5,2x,e11.5)"
  
  REAL(KIND=euwp) :: umx,umn,ums,vmx,vmn,vms,wmx,wmn,wms,  &
    pmx,pmn,pms,psmx,psmn,psms,            &
    pexmx,pexmn,pexms,pextmx,pextmn,pextms, &
    thmx,thmn,thms,thtmx,thtmn,thtms,rhomx,rhomn,rhoms,  &
    thpmx,thpmn,thpms, &
    zsmx,zsmn,zsms,zhmx,zhmn,zhms,         &
    oxmx,oxmn,oxms,oymx,oymn,oyms,ozmx,ozmn,ozms,    &
    uemx,uemn,uems,vemx,vemn,vems,themx,themn,thems, &
    rhemx,rhemn,rhems,ppemx,ppemn,ppems, &
    gmmx,gmmn,gmav, &
    g11mx,g11mn,g11ms,g12mx,g12mn,g12ms,g13mx,g13mn,g13ms, &
    g21mx,g21mn,g21ms,g22mx,g22mn,g22ms,g23mx,g23mn,g23ms, &
    g33mx,g33mn,g33ms,xcmx,xcmn,xcms,ycmx,ycmn,ycms,zcmx,zcmn,zcms
  
  REAL(KIND=euwp) :: fcr2mx, fcr2mn, fcr2ms, fcr3mx, fcr3mn, fcr3ms
  
  !----------------------------------------------------------------------!
  ! Get full potential temperature
  !----------------------------------------------------------------------!
  scr4(1:np,1:mp,1:lp) = th(1:np,1:mp,1:lp) + the(1:np,1:mp,1:lp)


  !----------------------------------------------------------------------!
  ! Prepare pressure in [Pa]
  !----------------------------------------------------------------------!
  CALL transform_exner_to_pascals(p,tempp,ppe,np,mp,lp,ih)
  ! CALL filstrP(tempp)
  IF (gndedge == 1) THEN
    psurf(1:np,1:mp)=tempp(1:np,1:mp,1)
  ENDIF

  !----------------------------------------------------------------------!
  ! TODO
  !----------------------------------------------------------------------!
  scr5(1:np,1:mp,1:lp)=rho(1:np,1:mp,1:lp)*gaci(1:np,1:mp,1:lp)

  
  !----------------------------------------------------------------------!
  ! Calc averages, max and mins of 3D variables 
  !----------------------------------------------------------------------!
  CALL calcavgmaxmin(u,umx,umn,ums)
  CALL calcavgmaxmin(v,vmx,vmn,vms)
  CALL calcavgmaxmin(w,wmx,wmn,wms)
!  CALL calcavgmaxmin(p,pmx,pmn,pms)
  CALL calcavgmaxmin(tempp,pmx,pmn,pms)
  CALL calcavgmaxmin(scr4,thmx,thmn,thms)
  CALL calcavgmaxmin(th,thpmx,thpmn,thpms)
  CALL calcavgmaxmin(ox,oxmx,oxmn,oxms)
  CALL calcavgmaxmin(oy,oymx,oymn,oyms)
  CALL calcavgmaxmin(oz,ozmx,ozmn,ozms)
  CALL calcavgmaxminflt(psurf,psmx,psmn,psms)
  IF (icmprss == 1) THEN
    CALL calcavgmaxmin(rho,rhomx,rhomn,rhoms)
    CALL calcavgmaxmin(spexi*pext,pextmx,pextmn,pextms)
  ENDIF
  IF (ipsim == 1) THEN
    CALL calcavgmaxmin(tht,thtmx,thtmn,thtms)
  ENDIF
  IF ( ntstep  == 0) THEN
    CALL calcavgmaxmin(ue,uemx,uemn,uems)
    CALL calcavgmaxmin(ve,vemx,vemn,vems)

    CALL scalcfltavgmaxmin(fcr2,fcr2mx,fcr2mn,fcr2ms)
    CALL scalcfltavgmaxmin(fcr3,fcr3mx,fcr3mn,fcr3ms)

    CALL calcavgmaxmin(the,themx,themn,thems)
    CALL calcavgmaxmin(ppe,ppemx,ppemn,ppems)
    IF ( icmprss == 1) THEN
      CALL calcavgmaxmin(rhe,rhemx,rhemn,rhems)
      CALL calcavgmaxmin(pexe,pexmx,pexmn,pexms)
    ENDIF
    CALL calcavgmaxmin(gmm,gmmx,gmmn,gmms)
    CALL scalcavgmaxmin(g11,g11mx,g11mn,g11ms)
    CALL scalcavgmaxmin(g12,g12mx,g12mn,g12ms)
    CALL scalcavgmaxmin(g13,g13mx,g13mn,g13ms)
    CALL scalcavgmaxmin(g21,g21mx,g21mn,g21ms)
    CALL scalcavgmaxmin(g22,g22mx,g22mn,g22ms)
    CALL scalcavgmaxmin(g23,g23mx,g23mn,g23ms)
    CALL scalcavgmaxmin(g33,g33mx,g33mn,g33ms)
    CALL calcavgmaxminflt(zs,zsmx,zsmn,zsms)
    CALL calcavgmaxminflt(zh,zhmx,zhmn,zhms)
    CALL calcavgmaxminflt(xcr,xcmx,xcmn,xcms)
    CALL calcavgmaxminflt(ycr,ycmx,ycmn,ycms)
    CALL calcavgmaxmin(zcr,zcmx,zcmn,zcms)
  ENDIF
  CALL calcavgmaxmin(div,divmx,divmn,divav)


  !Store local averages,maxima and minima to the buffer
  CALL wrttobuf(umx,umn,ums,bufavg,bufmax,bufmin,1)
  CALL wrttobuf(vmx,vmn,vms,bufavg,bufmax,bufmin,2)
  CALL wrttobuf(wmx,wmn,wms,bufavg,bufmax,bufmin,3)
  CALL wrttobuf(pmx,pmn,pms,bufavg,bufmax,bufmin,4)
  CALL wrttobuf(thmx,thmn,thms,bufavg,bufmax,bufmin,5)
  CALL wrttobuf(oxmx,oxmn,oxms,bufavg,bufmax,bufmin,6)
  CALL wrttobuf(oymx,oymn,oyms,bufavg,bufmax,bufmin,7)
  CALL wrttobuf(ozmx,ozmn,ozms,bufavg,bufmax,bufmin,8)
  CALL wrttobuf(psmx,psmn,psms,bufavg,bufmax,bufmin,9)
  IF (ipsim == 1) THEN
    CALL wrttobuf(thtmx,thtmn,thtms,bufavg,bufmax,bufmin,10)
  ENDIF
  IF (icmprss == 1) THEN
    CALL wrttobuf(rhomx,rhomn,rhoms,bufavg,bufmax,bufmin,11)
    CALL wrttobuf(pextmx,pextmn,pextms,bufavg,bufmax,bufmin,12)
  ENDIF
  IF ( ntstep  == 0) THEN
    CALL wrttobuf(uemx,uemn,uems,bufavg,bufmax,bufmin,13)
    CALL wrttobuf(vemx,vemn,vems,bufavg,bufmax,bufmin,14)
    CALL wrttobuf(themx,themn,thems,bufavg,bufmax,bufmin,15)
    CALL wrttobuf(ppemx,ppemn,ppems,bufavg,bufmax,bufmin,16)
    IF (icmprss == 1) THEN
      CALL wrttobuf(rhomx,rhomn,rhoms,bufavg,bufmax,bufmin,17)
      CALL wrttobuf(pexmx,pexmn,pexms,bufavg,bufmax,bufmin,18)
    ENDIF
    CALL wrttobuf(gmmx,gmmn,gmms,bufavg,bufmax,bufmin,19)
    CALL wrttobuf(g11mx,g11mn,g11ms,bufavg,bufmax,bufmin,20)
    CALL wrttobuf(g12mx,g12mn,g12ms,bufavg,bufmax,bufmin,21)
    CALL wrttobuf(g13mx,g13mn,g13ms,bufavg,bufmax,bufmin,22)
    CALL wrttobuf(g21mx,g21mn,g21ms,bufavg,bufmax,bufmin,23)
    CALL wrttobuf(g22mx,g22mn,g22ms,bufavg,bufmax,bufmin,24)
    CALL wrttobuf(g23mx,g23mn,g23ms,bufavg,bufmax,bufmin,25)
    CALL wrttobuf(g33mx,g33mn,g33ms,bufavg,bufmax,bufmin,26)
    CALL wrttobuf(zsmx,zsmn,zsms,bufavg,bufmax,bufmin,27)
    CALL wrttobuf(zhmx,zhmn,zhms,bufavg,bufmax,bufmin,28)
    CALL wrttobuf(xcmx,xcmn,xcms,bufavg,bufmax,bufmin,29)
    CALL wrttobuf(ycmx,ycmn,ycms,bufavg,bufmax,bufmin,30)
    CALL wrttobuf(zcmx,zcmn,zcms,bufavg,bufmax,bufmin,31)
  ENDIF
  CALL wrttobuf(divmx,divmn,divav,bufavg,bufmax,bufmin,32)
  CALL wrttobuf(thpmx,thpmn,thpms,bufavg,bufmax,bufmin,33)
  CALL wrttobuf(fcr2mx,fcr2mn,fcr2ms,bufavg,bufmax,bufmin,34)
  CALL wrttobuf(fcr3mx,fcr3mn,fcr3ms,bufavg,bufmax,bufmin,35)

  !Compute sums, maxs and mins in one global operation
  CALL globsumaxminv(35)

  !Retrive GLOBAL averages,maxima and minima from the buffer
  CALL wrtfmbuf(umx,umn,ums,bufavg,bufmax,bufmin,1)
  CALL wrtfmbuf(vmx,vmn,vms,bufavg,bufmax,bufmin,2)
  CALL wrtfmbuf(wmx,wmn,wms,bufavg,bufmax,bufmin,3)
  CALL wrtfmbuf(pmx,pmn,pms,bufavg,bufmax,bufmin,4)
  CALL wrtfmbuf(thmx,thmn,thms,bufavg,bufmax,bufmin,5)
  CALL wrtfmbuf(oxmx,oxmn,oxms,bufavg,bufmax,bufmin,6)
  CALL wrtfmbuf(oymx,oymn,oyms,bufavg,bufmax,bufmin,7)
  CALL wrtfmbuf(ozmx,ozmn,ozms,bufavg,bufmax,bufmin,8)
  CALL wrtfmbuf(psmx,psmn,psms,bufavg,bufmax,bufmin,9)
  IF (ipsim == 1) THEN
    CALL wrtfmbuf(thtmx,thtmn,thtms,bufavg,bufmax,bufmin,10)
  ENDIF
  IF (icmprss == 1) THEN
    CALL wrtfmbuf(rhomx,rhomn,rhoms,bufavg,bufmax,bufmin,11)
    CALL wrtfmbuf(pextmx,pextmn,pextms,bufavg,bufmax,bufmin,12)
  ENDIF
  IF ( ntstep  == 0) THEN
    CALL wrtfmbuf(uemx,uemn,uems,bufavg,bufmax,bufmin,13)
    CALL wrtfmbuf(vemx,vemn,vems,bufavg,bufmax,bufmin,14)
    CALL wrtfmbuf(themx,themn,thems,bufavg,bufmax,bufmin,15)
    CALL wrtfmbuf(ppemx,ppemn,ppems,bufavg,bufmax,bufmin,16)
    IF (icmprss == 1) THEN
      CALL wrtfmbuf(rhomx,rhomn,rhoms,bufavg,bufmax,bufmin,17)
      CALL wrtfmbuf(pexmx,pexmn,pexms,bufavg,bufmax,bufmin,18)
    ENDIF
    CALL wrtfmbuf(gmmx,gmmn,gmms,bufavg,bufmax,bufmin,19)
    CALL wrtfmbuf(g11mx,g11mn,g11ms,bufavg,bufmax,bufmin,20)
    CALL wrtfmbuf(g12mx,g12mn,g12ms,bufavg,bufmax,bufmin,21)
    CALL wrtfmbuf(g13mx,g13mn,g13ms,bufavg,bufmax,bufmin,22)
    CALL wrtfmbuf(g21mx,g21mn,g21ms,bufavg,bufmax,bufmin,23)
    CALL wrtfmbuf(g22mx,g22mn,g22ms,bufavg,bufmax,bufmin,24)
    CALL wrtfmbuf(g23mx,g23mn,g23ms,bufavg,bufmax,bufmin,25)
    CALL wrtfmbuf(g33mx,g33mn,g33ms,bufavg,bufmax,bufmin,26)
    CALL wrtfmbuf(zsmx,zsmn,zsms,bufavg,bufmax,bufmin,27)
    CALL wrtfmbuf(zhmx,zhmn,zhms,bufavg,bufmax,bufmin,28)
    CALL wrtfmbuf(xcmx,xcmn,xcms,bufavg,bufmax,bufmin,29)
    CALL wrtfmbuf(ycmx,ycmn,ycms,bufavg,bufmax,bufmin,30)
    CALL wrtfmbuf(zcmx,zcmn,zcms,bufavg,bufmax,bufmin,31)
  ENDIF
  CALL wrtfmbuf(divmx,divmn,divav,bufavg,bufmax,bufmin,32)
  CALL wrtfmbuf(thpmx,thpmn,thpms,bufavg,bufmax,bufmin,33)
  CALL wrtfmbuf(fcr2mx,fcr2mn,fcr2ms,bufavg,bufmax,bufmin,34)
  CALL wrtfmbuf(fcr3mx,fcr3mn,fcr3ms,bufavg,bufmax,bufmin,35)
  
  
  IF (mype == 0) THEN
    WRITE(6,fm2010)
    PRINT *,'TIMESTEP:',ntstep,' time',REAL(ntstep*dt/3600,4), 'hours'
    WRITE(6,fm2010)
    PRINT fm2011,umx,umn,ums,vmx,vmn,vms,wmx,wmn,wms, &
    thpmx,thpmn,thpms,&
    thmx,thmn,thms,oxmx,oxmn,oxms,oymx,oymn,oyms,ozmx,ozmn,ozms, &
    pmx,pmn,pms,psmx,psmn,psms
    WRITE(6,fm2010)
    IF (ntstep == 0) THEN
      PRINT fm2012,uemx,uemn ,uems ,vemx ,vemn ,vems ,&
      fcr2mx,fcr2mn,fcr2ms,fcr3mx,fcr3mn,fcr3ms, &
      themx,themn,thems,ppemx,ppemn,ppems,&
      zsmx,zsmn,zsms,zhmx,zhmn,zhms,      &
      xcmx,xcmn,xcms,ycmx,ycmn,ycms,      &
      zcmx,zcmn,zcms,gmmx,gmmn,gmms,      &
      g11mx,g11mn,g11ms,  &
      g12mx,g12mn,g13ms,  &
      g13mx,g13mn,g13ms,  &
      g21mx,g21mn,g21ms,  &
      g22mx,g22mn,g22ms,  &
      g23mx,g23mn,g23ms,  &
      g33mx,g33mn,g33ms
      WRITE(6,fm2010)
    ENDIF
    IF (icmprss == 1) THEN
      PRINT fm2013,pextmx,pextmn,pextms,&
      rhomx,rhomn,rhoms 
    ENDIF
    IF (ipsim == 1) THEN
      PRINT fm2014,thtmx,thtmn,thtms
      WRITE(6,fm2010)
    ENDIF
  ENDIF
  
  !----------------------------------------------------------------------!
  ! Check eulerian divergence
  !----------------------------------------------------------------------!
  divmx=-1.e15_euwp
  divmn= 1.e15_euwp
  divav=0._euwp
  
  
  divmx=globmax(div,1-ih,np+ih,1-ih,mp+ih,1-ih,lp+ih,1,np,1,mp,1,lp) 
  divmn=globmin(div,1-ih,np+ih,1-ih,mp+ih,1-ih,lp+ih,1,np,1,mp,1,lp) 
  divav=globsum(div,1-ih,np+ih,1-ih,mp+ih,1-ih,lp+ih,1,np,1,mp,1,lp) 
  divav=divav*rnmli
  divsd=0._euwp
  DO k=1,lp
    DO j=1,mp
      DO i=1,np
        scr4(i,j,k)=(div(i,j,k)-divav)**2
      END DO
    END DO
  END DO
  divsd=globsum(scr4,1-ih,np+ih,1-ih,mp+ih,1-ih,lp+ih,1,np,1,mp,1,lp) 
  divsd=sqrt(divsd*rnmli)
  
  divmx=divmx*dt
  divmn=divmn*dt
  divav=divav*dt
  divsd=divsd*dt
  nitav=nitsm/max(icount,1)
  mitav=mitsm/max(jcount,1)
  IF (mype == 0) THEN
    PRINT fm205, divmx,divmn,divav,divsd,eer,eem,niter,nitav,miter,mitav
    !    IF(ipsim == 1) THEN
    !      PRINT fm2022,dftmx,dftmn,dftav,dftsd
    !    ENDIF
  ENDIF
  
  
END SUBROUTINE diagnos
#endif
!----------------------------------------------------------------------!
SUBROUTINE diagnos_metrics(iwrite,np,mp,lp,ih)
!----------------------------------------------------------------------!

  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin,calcavgmaxminflt
  USE eulag_diagutils, ONLY: scalcfltavgmaxmin,scalcavgmaxmin
  USE eulag_diagutils, ONLY: calcfltavgmaxmin
  USE geometry_datafields, ONLY: g11, g12, g13
  USE geometry_datafields, ONLY: g21, g22, g23
  USE geometry_datafields, ONLY: g33
  USE geometry_datafields, ONLY: sg11, sg12, sg13
  USE geometry_datafields, ONLY: sg21, sg22, sg23
  USE geometry_datafields, ONLY: sg33
  USE geometry_datafields, ONLY: gf11, gf12, gf13
  USE geometry_datafields, ONLY: gf21, gf22, gf23
  USE geometry_datafields, ONLY: gf31, gf32, gf33
  USE geometry_datafields, ONLY: xcr,ycr,sina,cosa,tnga,fcr2,fcr3, &
                      zcr,gac,gmm,sinx,cosx !,g33i,g11g22Mg12g21i,gaci


  INTEGER(KIND=iintegers),INTENT(IN) ::  iwrite,np,mp,lp,ih

  REAL(KIND=euwp) :: &
      bufavg(30), &
      bufmax(30), &
      bufmin(30)
#ifdef PNETCDF
  REAL(KIND=euwp) :: flattemp(1-ih:np+ih,1-ih:mp+ih)
#endif
 

  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
  "(    1x,'fcr2mx, fcr2mn, fcr2av:',1X,3(e14.7,1x),"// &
    "   1x,'fcr3mx, fcr3mn, fcr3av:',1X,3(e14.7,1x),"// &
    "   1x,'xcrmax, xcrmin, xcrave:',1X,3(e14.7,1x),"// &
    "   1x,'ycrmax, ycrmin, ycrave:',1X,3(e14.7,1x),"// &
    "   1x,'zcrmax, zcrmin, zcrave:',1X,3(e14.7,1x),"// &
    "   1x,'gmmmax, gmmmin, gmmave:',1X,3(e14.7,1x),"// &
    "   1x,'gacmax, gacmin, gacave:',1X,3(e14.7,1x),"// &
    "   1x,'g11max, g11min, g11ave:',1X,3(e14.7,1x),"// &
    "   1x,'g12max, g12min, g12ave:',1X,3(e14.7,1x),"// &
    "   1x,'g13max, g13min, g13ave:',1X,3(e14.7,1x),"// &
    "   1x,'g21max, g21min, g21ave:',1X,3(e14.7,1x),"// &
    "   1x,'g22max, g22min, g22ave:',1X,3(e14.7,1x),"// &
    "   1x,'g23max, g23min, g23ave:',1X,3(e14.7,1x),"// &
    "   1x,'g33max, g33min, g33ave:',1X,3(e14.7,1x),"// &
    "   1x,'sg11max, sg11min, sg11ave:',1X,3(e14.7,1x),"// &
    "   1x,'sg12max, sg12min, sg12ave:',1X,3(e14.7,1x),"// &
    "   1x,'sg13max, sg13min, sg13ave:',1X,3(e14.7,1x),"// &
    "   1x,'sg21max, sg21min, sg21ave:',1X,3(e14.7,1x),"// &
    "   1x,'sg22max, sg22min, sg22ave:',1X,3(e14.7,1x),"// &
    "   1x,'sg23max, sg23min, sg23ave:',1X,3(e14.7,1x),"// &
    "   1x,'sg33max, sg33min, sg33ave:',1X,3(e14.7,1x),"// &
    "   1x,'gf11max, gf11min, gf11ave:',1X,3(e14.7,1x),"// &
    "   1x,'gf12max, gf12min, gf12ave:',1X,3(e14.7,1x),"// &
    "   1x,'gf13max, gf13min, gf13ave:',1X,3(e14.7,1x),"// &
    "   1x,'gf21max, gf21min, gf21ave:',1X,3(e14.7,1x),"// &
    "   1x,'gf22max, gf22min, gf22ave:',1X,3(e14.7,1x),"// &
    "   1x,'gf23max, gf23min, gf23ave:',1X,3(e14.7,1x),"// &
    "   1x,'gf31max, gf31min, gf31ave:',1X,3(e14.7,1x),"// &
    "   1x,'gf32max, gf32min, gf32ave:',1X,3(e14.7,1x),"// &
    "   1x,'gf33max, gf33min, gf33ave:',1X,3(e14.7,1x))"

  REAL(KIND=euwp) :: gmmx,gmmn,gmav,gcmx,gcmn,gcav, &  
    g11mx,g11mn,g11av,g12mx,g12mn,g12av,g13mx,g13mn,g13av, &
    g21mx,g21mn,g21av,g22mx,g22mn,g22av,g23mx,g23mn,g23av, &
    g33mx,g33mn,g33av, &
    sg11mx,sg11mn,sg11av,sg12mx,sg12mn,sg12av,sg13mx,sg13mn,sg13av, &
    sg21mx,sg21mn,sg21av,sg22mx,sg22mn,sg22av,sg23mx,sg23mn,sg23av, &
    sg33mx,sg33mn,sg33av, &
    gf11mx,gf11mn,gf11av,gf12mx,gf12mn,gf12av,gf13mx,gf13mn,gf13av, &
    gf21mx,gf21mn,gf21av,gf22mx,gf22mn,gf22av,gf23mx,gf23mn,gf23av, &
    gf31mx,gf31mn,gf31av,gf32mx,gf32mn,gf32av,gf33mx,gf33mn,gf33av, &
    xcmx,xcmn,xcav,ycmx,ycmn,ycav,zcmx,zcmn,zcav, & 
    cosmx,cosmn,cosav,sinmx,sinmn,sinav,tngmx,tngmn,tngav, &
    fcr2mx, fcr2mn, fcr2av, fcr3mx, fcr3mn, fcr3av

    CALL ttbeg(39)
  ! Calc averages, max and mins of 3D variables 
  !----------------------------------------------------------------------!
    CALL scalcfltavgmaxmin(fcr2,fcr2mx,fcr2mn,fcr2av,np,mp)
    CALL scalcfltavgmaxmin(fcr3,fcr3mx,fcr3mn,fcr3av,np,mp)
    CALL calcfltavgmaxmin(tnga,tngmx,tngmn,tngav,np,mp,ih)
    CALL calcavgmaxmin(gmm,gmmx,gmmn,gmav,np,mp,lp,ih)
    CALL calcavgmaxmin(gac,gcmx,gcmn,gcav,np,mp,lp,ih)
    CALL calcavgmaxmin(g11,g11mx,g11mn,g11av,np,mp,lp,ih)
    CALL calcavgmaxmin(g12,g12mx,g12mn,g12av,np,mp,lp,ih)
    CALL calcavgmaxmin(g13,g13mx,g13mn,g13av,np,mp,lp,ih)
    CALL calcavgmaxmin(g21,g21mx,g21mn,g21av,np,mp,lp,ih)
    CALL calcavgmaxmin(g22,g22mx,g22mn,g22av,np,mp,lp,ih)
    CALL calcavgmaxmin(g23,g23mx,g23mn,g23av,np,mp,lp,ih)
    CALL calcavgmaxmin(g33,g33mx,g33mn,g33av,np,mp,lp,ih)
    CALL calcavgmaxmin(gf11,gf11mx,gf11mn,gf11av,np,mp,lp,ih)
    CALL calcavgmaxmin(gf12,gf12mx,gf12mn,gf12av,np,mp,lp,ih)
    CALL calcavgmaxmin(gf13,gf13mx,gf13mn,gf13av,np,mp,lp,ih)
    CALL calcavgmaxmin(gf21,gf21mx,gf21mn,gf21av,np,mp,lp,ih)
    CALL calcavgmaxmin(gf22,gf22mx,gf22mn,gf22av,np,mp,lp,ih)
    CALL calcavgmaxmin(gf23,gf23mx,gf23mn,gf23av,np,mp,lp,ih)
    CALL calcavgmaxmin(gf31,gf31mx,gf31mn,gf31av,np,mp,lp,ih)
    CALL calcavgmaxmin(gf32,gf32mx,gf32mn,gf32av,np,mp,lp,ih)
    CALL calcavgmaxmin(gf33,gf33mx,gf33mn,gf33av,np,mp,lp,ih)
    CALL calcavgmaxmin(sg11,sg11mx,sg11mn,sg11av,np,mp,lp,ih)
    CALL calcavgmaxmin(sg12,sg12mx,sg12mn,sg12av,np,mp,lp,ih)
    CALL calcavgmaxmin(sg13,sg13mx,sg13mn,sg13av,np,mp,lp,ih)
    CALL calcavgmaxmin(sg21,sg21mx,sg21mn,sg21av,np,mp,lp,ih)
    CALL calcavgmaxmin(sg22,sg22mx,sg22mn,sg22av,np,mp,lp,ih)
    CALL calcavgmaxmin(sg23,sg23mx,sg23mn,sg23av,np,mp,lp,ih)
    CALL calcavgmaxmin(sg33,sg33mx,sg33mn,sg33av,np,mp,lp,ih)
    CALL calcavgmaxminflt(xcr,xcmx,xcmn,xcav,np,mp,ih)
    CALL calcavgmaxminflt(ycr,ycmx,ycmn,ycav,np,mp,ih)
    CALL calcavgmaxmin(zcr,zcmx,zcmn,zcav,np,mp,lp,ih)
    CALL calcavgmaxminflt(cosa,cosmx,cosmn,cosav,np,mp,ih)
    CALL calcavgmaxminflt(sina,sinmx,sinmn,sinav,np,mp,ih)

    CALL wrttobuf(g11mx,g11mn,g11av,bufavg,bufmax,bufmin,1)
    CALL wrttobuf(g12mx,g12mn,g12av,bufavg,bufmax,bufmin,2)
    CALL wrttobuf(g13mx,g13mn,g13av,bufavg,bufmax,bufmin,3)
    CALL wrttobuf(g21mx,g21mn,g21av,bufavg,bufmax,bufmin,4)
    CALL wrttobuf(g22mx,g22mn,g22av,bufavg,bufmax,bufmin,5)
    CALL wrttobuf(g23mx,g23mn,g23av,bufavg,bufmax,bufmin,6)
    CALL wrttobuf(g33mx,g33mn,g33av,bufavg,bufmax,bufmin,7)
    CALL wrttobuf(xcmx,xcmn,xcav,bufavg,bufmax,bufmin,8)
    CALL wrttobuf(ycmx,ycmn,ycav,bufavg,bufmax,bufmin,9)
    CALL wrttobuf(zcmx,zcmn,zcav,bufavg,bufmax,bufmin,10)
    CALL wrttobuf(sg11mx,sg11mn,sg11av,bufavg,bufmax,bufmin,11)
    CALL wrttobuf(sg12mx,sg12mn,sg12av,bufavg,bufmax,bufmin,12)
    CALL wrttobuf(sg13mx,sg13mn,sg13av,bufavg,bufmax,bufmin,13)
    CALL wrttobuf(sg21mx,sg21mn,sg21av,bufavg,bufmax,bufmin,14)
    CALL wrttobuf(sg22mx,sg22mn,sg22av,bufavg,bufmax,bufmin,15)
    CALL wrttobuf(sg23mx,sg23mn,sg23av,bufavg,bufmax,bufmin,16)
    CALL wrttobuf(sg33mx,sg33mn,sg33av,bufavg,bufmax,bufmin,17)
    CALL wrttobuf(gmmx,gmmn,gmav,bufavg,bufmax,bufmin,18)
    CALL wrttobuf(fcr2mx,fcr2mn,fcr2av,bufavg,bufmax,bufmin,19)
    CALL wrttobuf(fcr3mx,fcr3mn,fcr3av,bufavg,bufmax,bufmin,20)
    CALL wrttobuf(gf11mx,gf11mn,gf11av,bufavg,bufmax,bufmin,21)
    CALL wrttobuf(gf12mx,gf12mn,gf12av,bufavg,bufmax,bufmin,22)
    CALL wrttobuf(gf13mx,gf13mn,gf13av,bufavg,bufmax,bufmin,23)
    CALL wrttobuf(gf21mx,gf21mn,gf21av,bufavg,bufmax,bufmin,24)
    CALL wrttobuf(gf22mx,gf22mn,gf22av,bufavg,bufmax,bufmin,25)
    CALL wrttobuf(gf23mx,gf23mn,gf23av,bufavg,bufmax,bufmin,26)
    CALL wrttobuf(gf31mx,gf31mn,gf31av,bufavg,bufmax,bufmin,27)
    CALL wrttobuf(gf32mx,gf32mn,gf32av,bufavg,bufmax,bufmin,28)
    CALL wrttobuf(gf33mx,gf33mn,gf33av,bufavg,bufmax,bufmin,29)
    CALL wrttobuf(gcmx,gcmn,gcav,bufavg,bufmax,bufmin,30)
  !Compute sums, maxs and mins in one global operation
    CALL globsumaxminv(bufavg,bufmax,bufmin,30_iintegers)

    CALL wrtfmbuf(g11mx,g11mn,g11av,bufavg,bufmax,bufmin,1)
    CALL wrtfmbuf(g12mx,g12mn,g12av,bufavg,bufmax,bufmin,2)
    CALL wrtfmbuf(g13mx,g13mn,g13av,bufavg,bufmax,bufmin,3)
    CALL wrtfmbuf(g21mx,g21mn,g21av,bufavg,bufmax,bufmin,4)
    CALL wrtfmbuf(g22mx,g22mn,g22av,bufavg,bufmax,bufmin,5)
    CALL wrtfmbuf(g23mx,g23mn,g23av,bufavg,bufmax,bufmin,6)
    CALL wrtfmbuf(g33mx,g33mn,g33av,bufavg,bufmax,bufmin,7)
    CALL wrtfmbuf(xcmx,xcmn,xcav,bufavg,bufmax,bufmin,8)
    CALL wrtfmbuf(ycmx,ycmn,ycav,bufavg,bufmax,bufmin,9)
    CALL wrtfmbuf(zcmx,zcmn,zcav,bufavg,bufmax,bufmin,10)
    CALL wrtfmbuf(sg11mx,sg11mn,sg11av,bufavg,bufmax,bufmin,11)
    CALL wrtfmbuf(sg12mx,sg12mn,sg12av,bufavg,bufmax,bufmin,12)
    CALL wrtfmbuf(sg13mx,sg13mn,sg13av,bufavg,bufmax,bufmin,13)
    CALL wrtfmbuf(sg21mx,sg21mn,sg21av,bufavg,bufmax,bufmin,14)
    CALL wrtfmbuf(sg22mx,sg22mn,sg22av,bufavg,bufmax,bufmin,15)
    CALL wrtfmbuf(sg23mx,sg23mn,sg23av,bufavg,bufmax,bufmin,16)
    CALL wrtfmbuf(sg33mx,sg33mn,sg33av,bufavg,bufmax,bufmin,17)
    CALL wrtfmbuf(gmmx,gmmn,gmav,bufavg,bufmax,bufmin,18)
    CALL wrtfmbuf(fcr2mx,fcr2mn,fcr2av,bufavg,bufmax,bufmin,19)
    CALL wrtfmbuf(fcr3mx,fcr3mn,fcr3av,bufavg,bufmax,bufmin,20)
    CALL wrtfmbuf(gf11mx,gf11mn,gf11av,bufavg,bufmax,bufmin,21)
    CALL wrtfmbuf(gf12mx,gf12mn,gf12av,bufavg,bufmax,bufmin,22)
    CALL wrtfmbuf(gf13mx,gf13mn,gf13av,bufavg,bufmax,bufmin,23)
    CALL wrtfmbuf(gf21mx,gf21mn,gf21av,bufavg,bufmax,bufmin,24)
    CALL wrtfmbuf(gf22mx,gf22mn,gf22av,bufavg,bufmax,bufmin,25)
    CALL wrtfmbuf(gf23mx,gf23mn,gf23av,bufavg,bufmax,bufmin,26)
    CALL wrtfmbuf(gf31mx,gf31mn,gf31av,bufavg,bufmax,bufmin,27)
    CALL wrtfmbuf(gf32mx,gf32mn,gf32av,bufavg,bufmax,bufmin,28)
    CALL wrtfmbuf(gf33mx,gf33mn,gf33av,bufavg,bufmax,bufmin,29)
    CALL wrtfmbuf(gcmx,gcmn,gcav,bufavg,bufmax,bufmin,30)

  IF (mype == 0) THEN
  
      PRINT fm2012, fcr2mx,fcr2mn,fcr2av,fcr3mx,fcr3mn,fcr3av, &
      xcmx,xcmn,xcav,ycmx,ycmn,ycav,zcmx,zcmn,zcav, &
      gmmx,gmmn,gmav,gcmx,gcmn,gcav,      &
      g11mx,g11mn,g11av,  &
      g12mx,g12mn,g13av,  &
      g13mx,g13mn,g13av,  &
      g21mx,g21mn,g21av,  &
      g22mx,g22mn,g22av,  &
      g23mx,g23mn,g23av,  &
      g33mx,g33mn,g33av,  &
      sg11mx,sg11mn,sg11av,  &
      sg12mx,sg12mn,sg13av,  &
      sg13mx,sg13mn,sg13av,  &
      sg21mx,sg21mn,sg21av,  &
      sg22mx,sg22mn,sg22av,  &
      sg23mx,sg23mn,sg23av,  &
      sg33mx,sg33mn,sg33av,  &
      gf11mx,gf11mn,gf11av,  &
      gf12mx,gf12mn,gf13av,  &
      gf13mx,gf13mn,gf13av,  &
      gf21mx,gf21mn,gf21av,  &
      gf22mx,gf22mn,gf22av,  &
      gf23mx,gf23mn,gf23av,  &
      gf31mx,gf31mn,gf31av,  &
      gf32mx,gf32mn,gf32av,  &
      gf33mx,gf33mn,gf33av
      WRITE(6,fm2010)
  ENDIF
    CALL ttend(39)
  IF(iwrite == 1) THEN
#ifdef PNETCDF
  flattemp(1:np,1:mp)=fcr2(1:np,1:mp)
  CALL pnet_out_chunk('fcr2     ','tstng.nc',35,1,1,1,0,0,flattemp,np,mp,lp,ih)
  flattemp(1:np,1:mp)=fcr3(1:np,1:mp)
  CALL pnet_out_chunk('fcr3     ','tstng.nc',35,1,1,1,0,0,flattemp,np,mp,lp,ih)
  CALL pnet_out_chunk('xcr      ','tstng.nc',35,1,1,1,0,0,xcr,np,mp,lp,ih)
  CALL pnet_out_chunk('ycr      ','tstng.nc',35,1,1,1,0,0,ycr,np,mp,lp,ih)
  CALL pnet_out_chunk('sina     ','tstng.nc',35,1,1,1,0,0,sina,np,mp,lp,ih)
  CALL pnet_out_chunk('cosa     ','tstng.nc',35,1,1,1,0,0,cosa,np,mp,lp,ih)
  CALL pnet_out_chunk('tnga     ','tstng.nc',35,1,1,1,0,0,tnga,np,mp,lp,ih)
  CALL pnet_out_chunk('cosx     ','tstng.nc',35,1,1,1,0,0,cosx,np,mp,lp,ih)
  CALL pnet_out_chunk('sinx     ','tstng.nc',35,1,1,1,0,0,sinx,np,mp,lp,ih)
  CALL pnet_out_chunk('g11      ','tstn0.nc',1,1,1,1,0,0,g11,np,mp,lp,ih)
  CALL pnet_out_chunk('g12      ','tstn0.nc',1,1,1,1,0,0,g12,np,mp,lp,ih)
  CALL pnet_out_chunk('g13      ','tstn0.nc',1,1,1,1,0,0,g13,np,mp,lp,ih)
  CALL pnet_out_chunk('g21      ','tstn0.nc',1,1,1,1,0,0,g21,np,mp,lp,ih)
  CALL pnet_out_chunk('g22      ','tstn0.nc',1,1,1,1,0,0,g22,np,mp,lp,ih)
  CALL pnet_out_chunk('g23      ','tstn0.nc',1,1,1,1,0,0,g23,np,mp,lp,ih)
  CALL pnet_out_chunk('g33      ','tstn0.nc',1,1,1,1,0,0,g33,np,mp,lp,ih)
  CALL pnet_out_chunk('zcr      ','tstn0.nc',1,1,1,1,0,0,zcr,np,mp,lp,ih)
  CALL pnet_out_chunk('gac      ','tstn0.nc',1,1,1,1,0,0,gac,np,mp,lp,ih)
  CALL pnet_out_chunk('gmm      ','tstn0.nc',1,1,1,1,0,0,gmm,np,mp,lp,ih)
  CALL pnet_out_chunk('sg11     ','tstn0.nc',1,1,1,1,0,0,sg11,np,mp,lp,ih)
  CALL pnet_out_chunk('sg12     ','tstn0.nc',1,1,1,1,0,0,sg12,np,mp,lp,ih)
  CALL pnet_out_chunk('sg13     ','tstn0.nc',1,1,1,1,0,0,sg13,np,mp,lp,ih)
  CALL pnet_out_chunk('sg21     ','tstn0.nc',1,1,1,1,0,0,sg21,np,mp,lp,ih)
  CALL pnet_out_chunk('sg22     ','tstn0.nc',1,1,1,1,0,0,sg22,np,mp,lp,ih)
  CALL pnet_out_chunk('sg23     ','tstn0.nc',1,1,1,1,0,0,sg23,np,mp,lp,ih)
  CALL pnet_out_chunk('sg33     ','tstn0.nc',1,1,1,1,0,0,sg33,np,mp,lp,ih)
  CALL pnet_out_chunk('gf11     ','tstn0.nc',1,1,1,1,0,0,gf11,np,mp,lp,ih)
  CALL pnet_out_chunk('gf12     ','tstn0.nc',1,1,1,1,0,0,gf12,np,mp,lp,ih)
  CALL pnet_out_chunk('gf13     ','tstn0.nc',1,1,1,1,0,0,gf13,np,mp,lp,ih)
  CALL pnet_out_chunk('gf21     ','tstn0.nc',1,1,1,1,0,0,gf21,np,mp,lp,ih)
  CALL pnet_out_chunk('gf22     ','tstn0.nc',1,1,1,1,0,0,gf22,np,mp,lp,ih)
  CALL pnet_out_chunk('gf23     ','tstn0.nc',1,1,1,1,0,0,gf23,np,mp,lp,ih)
  CALL pnet_out_chunk('gf31     ','tstn0.nc',1,1,1,1,0,0,gf31,np,mp,lp,ih)
  CALL pnet_out_chunk('gf32     ','tstn0.nc',1,1,1,1,0,0,gf32,np,mp,lp,ih)
  CALL pnet_out_chunk('gf33     ','tstn0.nc',1,1,1,1,0,0,gf33,np,mp,lp,ih)
#endif
  ENDIF
  
END SUBROUTINE diagnos_metrics

!----------------------------------------------------------------------!
SUBROUTINE diagnos_metrics_vz(np,mp,lp,ih)
!----------------------------------------------------------------------!

  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin_vz
  USE geometry_datafields, ONLY: gf11, gf12, gf13
  USE geometry_datafields, ONLY: gf21, gf22, gf23
  USE geometry_datafields, ONLY: gf31, gf32, gf33


  INTEGER(KIND=iintegers),INTENT(IN) ::  np,mp,lp,ih
  INTEGER :: k

  REAL(KIND=euwp) :: &
      bufavg(lp*9), &
      bufmax(lp*9), &
      bufmin(lp*9)
 

  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
  "(    1x,'gf11max, gf11min, gf11ave:',1X,3(e14.7,1x),"// &
   "    1x,'gf12max, gf12min, gf12ave:',1X,3(e14.7,1x),"// &
   "    1x,'gf13max, gf13min, gf13ave:',1X,3(e14.7,1x),"// &
   "    1x,'gf21max, gf21min, gf21ave:',1X,3(e14.7,1x),"// &
   "    1x,'gf22max, gf22min, gf22ave:',1X,3(e14.7,1x),"// &
   "    1x,'gf23max, gf23min, gf23ave:',1X,3(e14.7,1x),"// &
   "    1x,'gf31max, gf31min, gf31ave:',1X,3(e14.7,1x),"// &
   "    1x,'gf32max, gf32min, gf32ave:',1X,3(e14.7,1x),"// &
   "    1x,'gf33max, gf33min, gf33ave:',1X,3(e14.7,1x))"

  REAL(KIND=euwp), DIMENSION(lp) :: &
                     gf11mx,gf11mn,gf11av,gf12mx,gf12mn,gf12av,gf13mx,gf13mn,gf13av, &
                     gf21mx,gf21mn,gf21av,gf22mx,gf22mn,gf22av,gf23mx,gf23mn,gf23av, &
                     gf31mx,gf31mn,gf31av,gf32mx,gf32mn,gf32av,gf33mx,gf33mn,gf33av

    CALL ttbeg(39)
  ! Calc averages, max and mins of 3D variables 
  !----------------------------------------------------------------------!
      CALL calcavgmaxmin_vz(gf11,gf11mx,gf11mn,gf11av,np,mp,lp,ih)
      CALL calcavgmaxmin_vz(gf12,gf12mx,gf12mn,gf12av,np,mp,lp,ih)
      CALL calcavgmaxmin_vz(gf13,gf13mx,gf13mn,gf13av,np,mp,lp,ih)
      CALL calcavgmaxmin_vz(gf21,gf21mx,gf21mn,gf21av,np,mp,lp,ih)
      CALL calcavgmaxmin_vz(gf22,gf22mx,gf22mn,gf22av,np,mp,lp,ih)
      CALL calcavgmaxmin_vz(gf23,gf23mx,gf23mn,gf23av,np,mp,lp,ih)
      CALL calcavgmaxmin_vz(gf31,gf31mx,gf31mn,gf31av,np,mp,lp,ih)
      CALL calcavgmaxmin_vz(gf32,gf32mx,gf32mn,gf32av,np,mp,lp,ih)
      CALL calcavgmaxmin_vz(gf33,gf33mx,gf33mn,gf33av,np,mp,lp,ih)
    DO k=1,lp
      CALL wrttobuf(gf11mx(k),gf11mn(k),gf11av(k),bufavg,bufmax,bufmin,(k-1)*9+1)
      CALL wrttobuf(gf12mx(k),gf12mn(k),gf12av(k),bufavg,bufmax,bufmin,(k-1)*9+2)
      CALL wrttobuf(gf13mx(k),gf13mn(k),gf13av(k),bufavg,bufmax,bufmin,(k-1)*9+3)
      CALL wrttobuf(gf21mx(k),gf21mn(k),gf21av(k),bufavg,bufmax,bufmin,(k-1)*9+4)
      CALL wrttobuf(gf22mx(k),gf22mn(k),gf22av(k),bufavg,bufmax,bufmin,(k-1)*9+5)
      CALL wrttobuf(gf23mx(k),gf23mn(k),gf23av(k),bufavg,bufmax,bufmin,(k-1)*9+6)
      CALL wrttobuf(gf31mx(k),gf31mn(k),gf31av(k),bufavg,bufmax,bufmin,(k-1)*9+7)
      CALL wrttobuf(gf32mx(k),gf32mn(k),gf32av(k),bufavg,bufmax,bufmin,(k-1)*9+8)
      CALL wrttobuf(gf33mx(k),gf33mn(k),gf33av(k),bufavg,bufmax,bufmin,(k-1)*9+9)
    ENDDO
  !Compute sums, maxs and mins in one global operation
    CALL globsumaxminv(bufavg,bufmax,bufmin,lp*9_iintegers)

    DO k=1,lp
      CALL wrtfmbuf(gf11mx(k),gf11mn(k),gf11av(k),bufavg,bufmax,bufmin,(k-1)*9+1)
      CALL wrtfmbuf(gf12mx(k),gf12mn(k),gf12av(k),bufavg,bufmax,bufmin,(k-1)*9+2)
      CALL wrtfmbuf(gf13mx(k),gf13mn(k),gf13av(k),bufavg,bufmax,bufmin,(k-1)*9+3)
      CALL wrtfmbuf(gf21mx(k),gf21mn(k),gf21av(k),bufavg,bufmax,bufmin,(k-1)*9+4)
      CALL wrtfmbuf(gf22mx(k),gf22mn(k),gf22av(k),bufavg,bufmax,bufmin,(k-1)*9+5)
      CALL wrtfmbuf(gf23mx(k),gf23mn(k),gf23av(k),bufavg,bufmax,bufmin,(k-1)*9+6)
      CALL wrtfmbuf(gf31mx(k),gf31mn(k),gf31av(k),bufavg,bufmax,bufmin,(k-1)*9+7)
      CALL wrtfmbuf(gf32mx(k),gf32mn(k),gf32av(k),bufavg,bufmax,bufmin,(k-1)*9+8)
      CALL wrtfmbuf(gf33mx(k),gf33mn(k),gf33av(k),bufavg,bufmax,bufmin,(k-1)*9+9)
    ENDDO

  IF (mype == 0) THEN
  
    DO k=1,lp
      PRINT *,'Level k=',k
      PRINT fm2012, & 
      gf11mx(k),gf11mn(k),gf11av(k),  &
      gf12mx(k),gf12mn(k),gf13av(k),  &
      gf13mx(k),gf13mn(k),gf13av(k),  &
      gf21mx(k),gf21mn(k),gf21av(k),  &
      gf22mx(k),gf22mn(k),gf22av(k),  &
      gf23mx(k),gf23mn(k),gf23av(k),  &
      gf31mx(k),gf31mn(k),gf31av(k),  &
      gf32mx(k),gf32mn(k),gf32av(k),  &
      gf33mx(k),gf33mn(k),gf33av(k)
      WRITE(6,fm2010)
    ENDDO
  ENDIF
    CALL ttend(39)
  
END SUBROUTINE diagnos_metrics_vz


!----------------------------------------------------------------------!
SUBROUTINE diagnos_reference(uref,vref,thref,np,mp,lp,ih,pextref,rhref)
!----------------------------------------------------------------------!

  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin


  !----------------------------------------------------------------------!

  INTEGER(KIND=iintegers) ::  np,mp,lp,ih
  REAL(KIND=euwp),INTENT(IN) :: &
     uref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     vref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
    thref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp),INTENT(IN),OPTIONAL :: &
    rhref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
   pextref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp) :: &
      bufavg(5), &
      bufmax(5), &
      bufmin(5)
 
  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
      "(1x,' urefmax,  urefmin,  urefave:',1X,3(e14.7,1x),"// &
      " 1x,' vrefmax,  vrefmin,  vrefave:',1X,3(e14.7,1x),"// &
      " 1x,'threfmax, threfmin, threfave:',1X,3(e14.7,1x))"

  CHARACTER(LEN=*),PARAMETER :: fm2013 = &
        "(1x,' rhrefmax,  rhrefmin,  rhrefave:',1X,3(e14.7,1x),"// &
        " 1x,'pexrefmax, pexrefmin, pexrefave:',1X,3(e14.7,1x))"

  REAL(KIND=euwp) ::     &
       urefmx,  urefmn, urefav, &
       vrefmx  ,vrefmn, vrefav, &
      threfmx, threfmn,threfav, &
      rhrefmx, rhrefmn,rhrefav, &
     pexrefmx,pexrefmn,pexrefav

    CALL ttbeg(39)
  
  ! Calc averages, max and mins of 3D variables 
  !----------------------------------------------------------------------!
    CALL calcavgmaxmin(uref,urefmx,urefmn,urefav,np,mp,lp,ih)
    CALL calcavgmaxmin(vref,vrefmx,vrefmn,vrefav,np,mp,lp,ih)
    CALL calcavgmaxmin(thref,threfmx,threfmn,threfav,np,mp,lp,ih)
    IF (PRESENT(rhref).AND.PRESENT(pextref)) THEN
      CALL calcavgmaxmin(rhref,rhrefmx,rhrefmn,rhrefav,np,mp,lp,ih)
      CALL calcavgmaxmin(pextref,pexrefmx,pexrefmn,pexrefav,np,mp,lp,ih)
    ENDIF

    CALL wrttobuf(urefmx,urefmn,urefav,bufavg,bufmax,bufmin,1)
    CALL wrttobuf(vrefmx,vrefmn,vrefav,bufavg,bufmax,bufmin,2)
    CALL wrttobuf(threfmx,threfmn,threfav,bufavg,bufmax,bufmin,3)
    IF (PRESENT(rhref).AND.PRESENT(pextref)) THEN
      CALL wrttobuf( rhrefmx, rhrefmn, rhrefav,bufavg,bufmax,bufmin,4)
      CALL wrttobuf(pexrefmx,pexrefmn,pexrefav,bufavg,bufmax,bufmin,5)
    ENDIF

  !Compute sums, maxs and mins in one global operation
  CALL globsumaxminv(bufavg,bufmax,bufmin,5)

  !Retrive GLOBAL averages,maxima and minima from the buffer
    CALL wrtfmbuf(urefmx,urefmn,urefav,bufavg,bufmax,bufmin,1)
    CALL wrtfmbuf(vrefmx,vrefmn,vrefav,bufavg,bufmax,bufmin,2)
    CALL wrtfmbuf(threfmx,threfmn,threfav,bufavg,bufmax,bufmin,3)
    IF (PRESENT(rhref).AND.PRESENT(pextref)) THEN
      CALL wrtfmbuf( rhrefmx, rhrefmn, rhrefav,bufavg,bufmax,bufmin,4)
      CALL wrtfmbuf(pexrefmx,pexrefmn,pexrefav,bufavg,bufmax,bufmin,5)
    ENDIF
  
  IF (mype == 0) THEN
    WRITE(6,fm2010)
      PRINT fm2012,       &
      urefmx, urefmn, urefav,   & 
      vrefmx ,vrefmn ,vrefav,   &
      threfmx,threfmn,threfav
      IF(PRESENT(rhref).AND.PRESENT(pextref)) THEN
        PRINT fm2013,     &
        rhrefmx,rhrefmn,rhrefav,&
        pexrefmx,pexrefmn,pexrefav
      ENDIF
    WRITE(6,fm2010)
  ENDIF
  
  
    CALL ttend(39)
  
END SUBROUTINE diagnos_reference
!----------------------------------------------------------------------!
SUBROUTINE diagnos_eulag_forcings(iwrite,fu,fv,fw,ft,np,mp,lp,ih)
!----------------------------------------------------------------------!

  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin


  !----------------------------------------------------------------------!

  INTEGER(KIND=iintegers),INTENT(IN) :: iwrite,np,mp,lp,ih
  REAL(KIND=euwp),INTENT(IN) :: &
     fu(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     fv(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     fw(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     ft(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp) :: &
      bufavg(4), &
      bufmax(4), &
      bufmin(4)
 
  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
      "(1x,' fxmax,  fumin,  fuave:',1X,3(e14.7,1x),"// &
      " 1x,' fvmax,  fvmin,  fvave:',1X,3(e14.7,1x),"// &
      " 1x,' fwmax,  fwmin,  fwave:',1X,3(e14.7,1x),"// &
      " 1x,' ftmax,  ftmin,  ftave:',1X,3(e14.7,1x))"

  REAL(KIND=euwp) ::     &
       fumx,  fumn, fuav, &
       fvmx,  fvmn, fvav, &
       fwmx,  fwmn, fwav, &
       ftmx,  ftmn, ftav

    CALL ttbeg(39)
  
  ! Calc averages, max and mins of 3D variables 
  !----------------------------------------------------------------------!
    CALL calcavgmaxmin(fu,fumx,fumn,fuav,np,mp,lp,ih)
    CALL calcavgmaxmin(fv,fvmx,fvmn,fvav,np,mp,lp,ih)
    CALL calcavgmaxmin(fw,fwmx,fwmn,fwav,np,mp,lp,ih)
    CALL calcavgmaxmin(ft,ftmx,ftmn,ftav,np,mp,lp,ih)

    CALL wrttobuf(fumx,fumn,fuav,bufavg,bufmax,bufmin,1)
    CALL wrttobuf(fvmx,fvmn,fvav,bufavg,bufmax,bufmin,2)
    CALL wrttobuf(fwmx,fwmn,fwav,bufavg,bufmax,bufmin,3)
    CALL wrttobuf(ftmx,ftmn,ftav,bufavg,bufmax,bufmin,4)

  !Compute sums, maxs and mins in one global operation
  CALL globsumaxminv(bufavg,bufmax,bufmin,4)

  !Retrive GLOBAL averages,maxima and minima from the buffer
    CALL wrtfmbuf(fumx,fumn,fuav,bufavg,bufmax,bufmin,1)
    CALL wrtfmbuf(fvmx,fvmn,fvav,bufavg,bufmax,bufmin,2)
    CALL wrtfmbuf(fwmx,fwmn,fwav,bufavg,bufmax,bufmin,3)
    CALL wrtfmbuf(ftmx,ftmn,ftav,bufavg,bufmax,bufmin,4)
  
  IF (mype == 0) THEN
    WRITE(6,fm2010)
      PRINT fm2012,       &
      fumx, fumn, fuav,   & 
      fvmx, fvmn, fvav,   & 
      fwmx, fwmn, fwav,   & 
      ftmx, ftmn, ftav      
    WRITE(6,fm2010)
  ENDIF
  IF (iwrite == 1) THEN
#ifdef PNETCDF
  CALL pnet_out_chunk('fu_eulag ','tstn1.nc',1,1,1,1,0,0,fu,np,mp,lp,ih)
  CALL pnet_out_chunk('fv_eulag ','tstn1.nc',1,1,1,1,0,0,fv,np,mp,lp,ih)
  CALL pnet_out_chunk('fw_eulag ','tstn1.nc',1,1,1,1,0,0,fw,np,mp,lp,ih)
  CALL pnet_out_chunk('ft_eulag ','tstn1.nc',1,1,1,1,0,0,ft,np,mp,lp,ih)
#endif
  ENDIF
  
  
    CALL ttend(39)
  
END SUBROUTINE diagnos_eulag_forcings
!----------------------------------------------------------------------!
SUBROUTINE diagnos_cosmo_forcings(iwrite,fu,fv,fw,ft,fqv,fqc,fqr,np,mp,lp,ih,fqi,fpext)
!----------------------------------------------------------------------!

  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin


  !----------------------------------------------------------------------!

  INTEGER(KIND=iintegers) ::  iwrite,np,mp,lp,ih
  REAL(KIND=euwp),INTENT(IN) :: &
     fu(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     fv(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     fw(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     ft(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
    fqv(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
    fqc(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
    fqr(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp),INTENT(IN),OPTIONAL :: &
    fqi(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),&
  fpext(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp) :: &
      bufavg(9), &
      bufmax(9), &
      bufmin(9)
 
  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
      "(1x,' fumax,  fumin,  fuave:',1X,3(e14.7,1x),"// &
      " 1x,' fvmax,  fvmin,  fvave:',1X,3(e14.7,1x),"// &
      " 1x,' fwmax,  fwmin,  fwave:',1X,3(e14.7,1x),"// &
      " 1x,' ftmax,  ftmin,  ftave:',1X,3(e14.7,1x),"// &
      " 1x,' fqvmax,  fqvmin,  fqvave:',1X,3(e14.7,1x),"// &
      " 1x,' fqcmax,  fqcmin,  fqcave:',1X,3(e14.7,1x),"// &
      " 1x,' fqrmax,  fqrmin,  fqrave:',1X,3(e14.7,1x))"

  CHARACTER(LEN=*),PARAMETER :: fm2013 = &
        "(1x,' fqimax,  fqimin,  fqiave:',1X,3(e14.7,1x))"
  CHARACTER(LEN=*),PARAMETER :: fm2014 = &
        "(1x,' fpemax,  fpemin,  fpeave:',1X,3(e14.7,1x))"

  REAL(KIND=euwp) ::     &
       fumx,  fumn, fuav, &
       fvmx,  fvmn, fvav, &
       fwmx,  fwmn, fwav, &
       ftmx,  ftmn, ftav, &
      fqvmx, fqvmn, fqvav, &
      fqcmx, fqcmn, fqcav, &
      fqrmx, fqrmn, fqrav, &
      fqimx, fqimn, fqiav, &
      fpemx, fpemn, fpeav

    CALL ttbeg(39)
  
  ! Calc averages, max and mins of 3D variables 
  !----------------------------------------------------------------------!
    CALL calcavgmaxmin(fu,fumx,fumn,fuav,np,mp,lp,ih)
    CALL calcavgmaxmin(fv,fvmx,fvmn,fvav,np,mp,lp,ih)
    CALL calcavgmaxmin(fw,fwmx,fwmn,fwav,np,mp,lp,ih)
    CALL calcavgmaxmin(ft,ftmx,ftmn,ftav,np,mp,lp,ih)
    CALL calcavgmaxmin(fqv,fqvmx,fqvmn,fqvav,np,mp,lp,ih)
    CALL calcavgmaxmin(fqc,fqcmx,fqcmn,fqcav,np,mp,lp,ih)
    CALL calcavgmaxmin(fqr,fqrmx,fqrmn,fqrav,np,mp,lp,ih)
    fqimx=0.;fqimn=0.;fqiav=0.
    IF (PRESENT(fqi)) CALL calcavgmaxmin(fqi,fqimx,fqimn,fqiav,np,mp,lp,ih)
    fpemx=0.;fpemn=0.;fpeav=0.
    IF (PRESENT(fpext)) CALL calcavgmaxmin(fpext,fpemx,fpemn,fpeav,np,mp,lp,ih)

    CALL wrttobuf(fumx,fumn,fuav,bufavg,bufmax,bufmin,1)
    CALL wrttobuf(fvmx,fvmn,fvav,bufavg,bufmax,bufmin,2)
    CALL wrttobuf(fwmx,fwmn,fwav,bufavg,bufmax,bufmin,3)
    CALL wrttobuf(ftmx,ftmn,ftav,bufavg,bufmax,bufmin,4)
    CALL wrttobuf(fqvmx,fqvmn,fqvav,bufavg,bufmax,bufmin,5)
    CALL wrttobuf(fqcmx,fqcmn,fqcav,bufavg,bufmax,bufmin,6)
    CALL wrttobuf(fqrmx,fqrmn,fqrav,bufavg,bufmax,bufmin,7)
    CALL wrttobuf(fqimx,fqimn,fqiav,bufavg,bufmax,bufmin,8)
    CALL wrttobuf(fpemx,fpemn,fpeav,bufavg,bufmax,bufmin,9)

  !Compute sums, maxs and mins in one global operation
  CALL globsumaxminv(bufavg,bufmax,bufmin,9)

  !Retrive GLOBAL averages,maxima and minima from the buffer
    CALL wrtfmbuf(fumx,fumn,fuav,bufavg,bufmax,bufmin,1)
    CALL wrtfmbuf(fvmx,fvmn,fvav,bufavg,bufmax,bufmin,2)
    CALL wrtfmbuf(fwmx,fwmn,fwav,bufavg,bufmax,bufmin,3)
    CALL wrtfmbuf(ftmx,ftmn,ftav,bufavg,bufmax,bufmin,4)
    CALL wrtfmbuf(fqvmx,fqvmn,fqvav,bufavg,bufmax,bufmin,5)
    CALL wrtfmbuf(fqcmx,fqcmn,fqcav,bufavg,bufmax,bufmin,6)
    CALL wrtfmbuf(fqrmx,fqrmn,fqrav,bufavg,bufmax,bufmin,7)
    CALL wrtfmbuf(fqimx,fqimn,fqiav,bufavg,bufmax,bufmin,8)
    CALL wrtfmbuf(fpemx,fpemn,fpeav,bufavg,bufmax,bufmin,9)
  
  IF (mype == 0) THEN
    WRITE(6,fm2010)
      PRINT fm2012,       &
      fumx, fumn, fuav,   & 
      fvmx, fvmn, fvav,   & 
      fwmx, fwmn, fwav,   & 
      ftmx, ftmn, ftav,   & 
      fqvmx, fqvmn, fqvav,   & 
      fqcmx, fqcmn, fqcav,   & 
      fqrmx, fqrmn, fqrav
      IF(PRESENT(fqi)) THEN
        PRINT fm2013,     &
        fqimx,fqimn,fqiav
      ENDIF
      IF(PRESENT(fpext)) THEN
        PRINT fm2014,     &
        fpemx,fpemn,fpeav
      ENDIF
    WRITE(6,fm2010)
  ENDIF
  IF (iwrite == 1) THEN
#ifdef PNETCDF
  CALL pnet_out_chunk('fu       ','tstn1.nc',1,1,1,1,0,0,fu,np,mp,lp,ih)
  CALL pnet_out_chunk('fv       ','tstn1.nc',1,1,1,1,0,0,fv,np,mp,lp,ih)
  CALL pnet_out_chunk('fw       ','tstn1.nc',1,1,1,1,0,0,fw,np,mp,lp,ih)
  CALL pnet_out_chunk('ft       ','tstn1.nc',1,1,1,1,0,0,ft,np,mp,lp,ih)
  CALL pnet_out_chunk('fqv      ','tstn1.nc',1,1,1,1,0,0,fqv,np,mp,lp,ih)
  CALL pnet_out_chunk('fqc      ','tstn1.nc',1,1,1,1,0,0,fqc,np,mp,lp,ih)
  CALL pnet_out_chunk('fpext    ','tstn1.nc',1,1,1,1,0,0,fpext,np,mp,lp,ih)
#endif
  ENDIF
  
  
    CALL ttend(39)
  
END SUBROUTINE diagnos_cosmo_forcings
!----------------------------------------------------------------------!
SUBROUTINE diagnos_actual(iwrite,u,v,w,th,np,mp,lp,ih,tht,pext,rho)
!----------------------------------------------------------------------!

  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin
  USE mod_parameters, ONLY : ibcx,ibcy,ibcz,ipoles


  !----------------------------------------------------------------------!

  INTEGER(KIND=iintegers),INTENT(IN) ::  iwrite,np,mp,lp,ih
  INTEGER(KIND=iintegers) :: iframe
  REAL(KIND=euwp),INTENT(IN) :: &
     u(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     v(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     w(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
    th(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp),INTENT(IN),OPTIONAL :: &
    tht(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
    rho(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
   pext(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp) :: &
      bufavg(6), &
      bufmax(6), &
      bufmin(6)
 
  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2011 = &
      "(1x,'LIB  umax,  umin,  uave:',1X,3(e14.7,1x),"// &
      " 1x,'LIB  vmax,  vmin,  vave:',1X,3(e14.7,1x),"// &
      " 1x,'LIB  wmax,  wmin,  wave:',1X,3(e14.7,1x),"// &
      " 1x,'LIB thmax, thmin, thave:',1X,3(e14.7,1x))"  
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
      "(1x,'LIB thtmax, thtmin, thtave:',1X,3(e14.7,1x))"

  CHARACTER(LEN=*),PARAMETER :: fm2013 = &
        "(1x,'LIB  rhomax,  rhomin,  rhoave:',1X,3(e14.7,1x),"// &
        " 1x,'LIB pextmax, pextmin, pextave:',1X,3(e14.7,1x))"

  REAL(KIND=euwp) ::     &
       umx,  umn, uav, &
       vmx  ,vmn, vav, &
       wmx  ,wmn, wav, &
      thmx, thmn,thav, &
      thtmx, thtmn,thtav, &
      rhmx, rhmn,rhav, &
     pexmx,pexmn,pexav
   INTEGER(KIND=iintegers) cntsum
    CALL ttbeg(39)

  ! Calc averages, max and mins of 3D variables 
  !----------------------------------------------------------------------!
    CALL calcavgmaxmin(u,umx,umn,uav,np,mp,lp,ih)
    CALL calcavgmaxmin(v,vmx,vmn,vav,np,mp,lp,ih)
    CALL calcavgmaxmin(w,wmx,wmn,wav,np,mp,lp,ih)
    CALL calcavgmaxmin(th ,thmx ,thmn ,thav,np,mp,lp,ih)
    IF (PRESENT(tht)) THEN
      CALL calcavgmaxmin(tht,thtmx,thtmn,thtav,np,mp,lp,ih)
    ENDIF
    IF (PRESENT(rho).AND.PRESENT(pext)) THEN
      CALL calcavgmaxmin(rho,rhmx,rhmn,rhav,np,mp,lp,ih)
      CALL calcavgmaxmin(pext*spexi,pexmx,pexmn,pexav,np,mp,lp,ih)
    ENDIF

    CALL wrttobuf(umx,umn,uav,bufavg,bufmax,bufmin,1)
    CALL wrttobuf(vmx,vmn,vav,bufavg,bufmax,bufmin,2)
    CALL wrttobuf(thmx,thmn,thav,bufavg,bufmax,bufmin,3)
      cntsum=3
    IF (PRESENT(tht)) THEN
      CALL wrttobuf(thtmx,thtmn,thtav,bufavg,bufmax,bufmin,4)
      cntsum=4
    ENDIF
    IF (PRESENT(rho).AND.PRESENT(pext)) THEN
      CALL wrttobuf( rhmx, rhmn, rhav,bufavg,bufmax,bufmin,5)
      CALL wrttobuf(pexmx,pexmn,pexav,bufavg,bufmax,bufmin,6)
      cntsum=6
    ENDIF

  !Compute sums, maxs and mins in one global operation
  CALL globsumaxminv(bufavg,bufmax,bufmin,cntsum)

  !Retrive GLOBAL averages,maxima and minima from the buffer
    CALL wrtfmbuf(umx,umn,uav,bufavg,bufmax,bufmin,1)
    CALL wrtfmbuf(vmx,vmn,vav,bufavg,bufmax,bufmin,2)
    CALL wrtfmbuf(thmx,thmn,thav,bufavg,bufmax,bufmin,3)
    IF (PRESENT(tht)) THEN
      CALL wrtfmbuf(thtmx,thtmn,thtav,bufavg,bufmax,bufmin,4)
    ENDIF
    IF (PRESENT(rho).AND.PRESENT(pext)) THEN
      CALL wrtfmbuf( rhmx, rhmn, rhav,bufavg,bufmax,bufmin,5)
      CALL wrtfmbuf(pexmx,pexmn,pexav,bufavg,bufmax,bufmin,6)
    ENDIF
  
  IF (mype == 0) THEN
    WRITE(6,fm2010)
      PRINT fm2011,       &
      umx, umn, uav,      & 
      vmx ,vmn ,vav,      &
      wmx ,wmn ,wav,      &
      thmx,thmn,thav
      IF (PRESENT(tht)) THEN
        PRINT fm2012,       &
        thtmx,thtmn,thtav
      ENDIF
      IF(PRESENT(rho).AND.PRESENT(pext)) THEN
        PRINT fm2013,        &
        rhmx,rhmn,rhav,      &
        pexmx,pexmn,pexav
      ENDIF
    WRITE(6,fm2010)
  ENDIF
  
  iframe=iwrite
  IF(iframe.ge.0) THEN
#ifdef PNETCDF
! If(mype.eq.0)  print *,'Nans before pnetcdf' 
! CALL flush(6)
!       CALL checknans (u   ,'u    ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
!       CALL checknans (v   ,'v    ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
!       CALL checknans (w   ,'w    ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
!       CALL checknans (th  ,'th   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
  CALL pnet_out_chunk('u       ','tstng.nc',1,1,1,1,iframe,0,u ,np,mp,lp,ih)
  CALL pnet_out_chunk('v       ','tstng.nc',1,1,1,1,iframe,0,v ,np,mp,lp,ih)
  CALL pnet_out_chunk('w       ','tstng.nc',1,1,1,1,iframe,0,w ,np,mp,lp,ih)
  CALL pnet_out_chunk('th      ','tstng.nc',1,1,1,1,iframe,0,th,np,mp,lp,ih)
! If(mype.eq.0)  print *,'Nans after pnetcdf' 
! CALL flush(6)
       CALL checknans (u   ,'u    ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
       CALL checknans (v   ,'v    ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
       CALL checknans (w   ,'w    ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
       CALL checknans (th  ,'th   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
       CALL checknans (tht ,'tht  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
       CALL checknans (pext,'pstr ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
       CALL checknans (rho ,'rho  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
       CALL checknans (th  ,'th   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
  IF(PRESENT(tht))  CALL pnet_out_chunk('tht     ','tstng.nc',1,1,1,1,iframe,0,tht ,np,mp,lp,ih)
  IF(PRESENT(rho))  CALL pnet_out_chunk('rho     ','tstng.nc',1,1,1,1,iframe,0,rho ,np,mp,lp,ih)
  IF(PRESENT(pext)) CALL pnet_out_chunk('pext    ','tstng.nc',1,1,1,1,iframe,0,pext,np,mp,lp,ih)
#endif
  ENDIF !iframe
  
    CALL ttend(39)
END SUBROUTINE diagnos_actual
!----------------------------------------------------------------------!
SUBROUTINE diagnos_actual_moist(iwrite,qv,qc,qr,qcrs,np,mp,lp,ih,isnap,qs,qi,qg)
!----------------------------------------------------------------------!

  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin


  !----------------------------------------------------------------------!

  INTEGER(KIND=iintegers) ::  iwrite,np,mp,lp,ih,isnap
  REAL(KIND=euwp),INTENT(IN) :: &
     qv(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     qc(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     qr(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     qcrs(np,mp,lp)

  REAL(KIND=euwp),INTENT(IN),OPTIONAL :: &
     qs(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     qi(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     qg(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp) :: &
     temp(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)

  REAL(KIND=euwp) :: &
      bufavg(7), &
      bufmax(7), &
      bufmin(7)
 
  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
      "(1x,' qvmax,  qvmin,  qvave:',1X,3(e14.7,1x),"// &
      " 1x,' qcmax,  qcmin,  qcave:',1X,3(e14.7,1x),"// &
      " 1x,' qrmax,  qrmin,  qrave:',1X,3(e14.7,1x),"// &
      " 1x,' qtmax,  qtmin,  qtave:',1X,3(e14.7,1x))"

  CHARACTER(LEN=*),PARAMETER :: fm2013 = &
      "(1x,' qsmax,  qsmin,  qsave:',1X,3(e14.7,1x),"// &
      " 1x,' qimax,  qimin,  qiave:',1X,3(e14.7,1x),"// &
      " 1x,' qgmax,  qgmin,  qgave:',1X,3(e14.7,1x))"

  REAL(KIND=euwp) ::     &
       qvmx, qvmn, qvav, &
       qcmx, qcmn, qcav, &
       qrmx, qrmn, qrav, &
       qtmx, qtmn, qtav, &
       qsmx, qsmn, qsav, &
       qimx, qimn, qiav, &
       qgmx, qgmn, qgav

    CALL ttbeg(39)

  ! Calc averages, max and mins of 3D variables 
  !----------------------------------------------------------------------!
    CALL calcavgmaxmin(qv,qvmx,qvmn,qvav,np,mp,lp,ih)
    CALL calcavgmaxmin(qc,qcmx,qcmn,qcav,np,mp,lp,ih)
    CALL calcavgmaxmin(qr,qrmx,qrmn,qrav,np,mp,lp,ih)
    temp(:,:,:)=0._euwp
    temp(1:np,1:mp,1:lp)=qcrs(1:np,1:mp,1:lp)
    CALL calcavgmaxmin(temp,qtmx,qtmn,qtav,np,mp,lp,ih)
    qsmx=0.
    qsmn=0.
    qsav=0.
    qimx=0.
    qimn=0.
    qiav=0.
    qgmx=0.
    qgmn=0.
    qgav=0.
    IF (PRESENT(qs)) THEN
      CALL calcavgmaxmin(qs,qsmx,qsmn,qsav,np,mp,lp,ih)
        IF (PRESENT(qi)) THEN
          CALL calcavgmaxmin(qi,qimx,qimn,qiav,np,mp,lp,ih)
            IF (PRESENT(qg)) THEN
              CALL calcavgmaxmin(qg,qgmx,qgmn,qgav,np,mp,lp,ih)
            ENDIF
        ENDIF
    ENDIF

    CALL wrttobuf(qvmx,qvmn,qvav,bufavg,bufmax,bufmin,1)
    CALL wrttobuf(qcmx,qcmn,qcav,bufavg,bufmax,bufmin,2)
    CALL wrttobuf(qrmx,qrmn,qrav,bufavg,bufmax,bufmin,3)
    CALL wrttobuf(qtmx,qtmn,qtav,bufavg,bufmax,bufmin,4)
    CALL wrttobuf(qsmx,qsmn,qsav,bufavg,bufmax,bufmin,5)
    CALL wrttobuf(qimx,qimn,qiav,bufavg,bufmax,bufmin,6)
    CALL wrttobuf(qgmx,qgmn,qgav,bufavg,bufmax,bufmin,7)

  !Compute sums, maxs and mins in one global operation
  CALL globsumaxminv(bufavg,bufmax,bufmin,7_iintegers)

  !Retrive GLOBAL averages,maxima and minima from the buffer
    CALL wrtfmbuf(qvmx,qvmn,qvav,bufavg,bufmax,bufmin,1)
    CALL wrtfmbuf(qcmx,qcmn,qcav,bufavg,bufmax,bufmin,2)
    CALL wrtfmbuf(qrmx,qrmn,qrav,bufavg,bufmax,bufmin,3)
    CALL wrtfmbuf(qtmx,qtmn,qtav,bufavg,bufmax,bufmin,4)
    CALL wrtfmbuf(qsmx,qsmn,qsav,bufavg,bufmax,bufmin,5)
    CALL wrtfmbuf(qimx,qimn,qiav,bufavg,bufmax,bufmin,6)
    CALL wrtfmbuf(qgmx,qgmn,qgav,bufavg,bufmax,bufmin,7)
  
  IF (mype == 0) THEN
    WRITE(6,fm2010)
      PRINT fm2012,       &
      qvmx, qvmn, qvav,      & 
      qcmx ,qcmn ,qcav,      &
      qrmx ,qrmn ,qrav,      &
      qtmx ,qtmn ,qtav
        PRINT fm2013,        &
        qsmx,qsmn,qsav,      &
        qimx,qimn,qiav,      &
        qgmx,qgmn,qgav
    WRITE(6,fm2010)
  ENDIF
  
   IF(iwrite.eq.1.and.isnap.eq.1) THEN
#ifdef PNETCDF
   CALL pnet_out_chunk('qv       ','tstn1.nc',1,1,1,1,0,0,qv,np,mp,lp,ih)
   CALL pnet_out_chunk('qc       ','tstn1.nc',1,1,1,1,0,0,qc,np,mp,lp,ih)
   CALL pnet_out_chunk('qr       ','tstn1.nc',1,1,1,1,0,0,qr,np,mp,lp,ih)
   CALL pnet_out_chunk('qt       ','tstn1.nc',1,1,1,1,0,0,temp,np,mp,lp,ih)
   IF (PRESENT(qs)) THEN
     CALL pnet_out_chunk('qs       ','tstn1.nc',1,1,1,1,0,0,qs,np,mp,lp,ih)
     IF (PRESENT(qi)) THEN
       CALL pnet_out_chunk('qi       ','tstn1.nc',1,1,1,1,0,0,qi,np,mp,lp,ih)
       IF (PRESENT(qg)) THEN
         CALL pnet_out_chunk('qg       ','tstn1.nc',1,1,1,1,0,0,qg,np,mp,lp,ih)
       ENDIF
     ENDIF
   ENDIF
        ELSE IF (isnap.eq.2) THEN
   CALL pnet_out_chunk('qv       ','tstn2.nc',1,1,1,1,0,0,qv,np,mp,lp,ih)
   CALL pnet_out_chunk('qc       ','tstn2.nc',1,1,1,1,0,0,qc,np,mp,lp,ih)
   CALL pnet_out_chunk('qr       ','tstn2.nc',1,1,1,1,0,0,qr,np,mp,lp,ih)
   IF (PRESENT(qs)) THEN
     CALL pnet_out_chunk('qs       ','tstn2.nc',1,1,1,1,0,0,qs,np,mp,lp,ih)
     IF (PRESENT(qi)) THEN
       CALL pnet_out_chunk('qi       ','tstn2.nc',1,1,1,1,0,0,qi,np,mp,lp,ih)
       IF (PRESENT(qg)) THEN
         CALL pnet_out_chunk('qg       ','tstn2.nc',1,1,1,1,0,0,qg,np,mp,lp,ih)
       ENDIF
     ENDIF
   ENDIF
        ELSE IF (isnap.eq.3) THEN
   CALL pnet_out_chunk('qv       ','tstn3.nc',1,1,1,1,0,0,qv,np,mp,lp,ih)
   CALL pnet_out_chunk('qc       ','tstn3.nc',1,1,1,1,0,0,qc,np,mp,lp,ih)
   CALL pnet_out_chunk('qr       ','tstn3.nc',1,1,1,1,0,0,qr,np,mp,lp,ih)
   IF (PRESENT(qs)) THEN
     CALL pnet_out_chunk('qs       ','tstn3.nc',1,1,1,1,0,0,qs,np,mp,lp,ih)
     IF (PRESENT(qi)) THEN
       CALL pnet_out_chunk('qi       ','tstn3.nc',1,1,1,1,0,0,qi,np,mp,lp,ih)
       IF (PRESENT(qg)) THEN
         CALL pnet_out_chunk('qg       ','tstn3.nc',1,1,1,1,0,0,qg,np,mp,lp,ih)
       ENDIF
     ENDIF
   ENDIF
#endif
   ENDIF !iwrite
  
    CALL ttend(39)
END SUBROUTINE diagnos_actual_moist
!----------------------------------------------------------------------!
SUBROUTINE diagnos_difference(iwrite,uref,vref,thref,the, &
                                 u,   v,th   ,tht, &
                                np,mp,lp,ih,       &
                           pextref,rhref,          &
                              pext,rho,isnap)
!----------------------------------------------------------------------!

  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv
  USE mod_parameters, ONLY: capi,rg,cmpex,spexi 
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin
  USE scratch_datafields, ONLY: scr1



  !----------------------------------------------------------------------!

  INTEGER(KIND=iintegers),INTENT(IN) ::  iwrite,np,mp,lp,ih
  REAL(KIND=euwp),INTENT(IN) :: &
     uref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     vref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
    thref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
      the(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
  REAL(KIND=euwp),INTENT(IN) :: &
     u(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
     v(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
    th(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
   tht(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp),INTENT(IN),OPTIONAL :: &
    rhref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
  pextref(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
  REAL(KIND=euwp),INTENT(IN),OPTIONAL :: &
    rho(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
   pext(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp) :: &
      bufavg(6), &
      bufmax(6), &
      bufmin(6)
 

  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
      "(1x,' udifmax,  udifmin,  udifave:',1X,3(e14.7,1x),"// &
      " 1x,' vdifmax,  vdifmin,  vdifave:',1X,3(e14.7,1x),"// &
      " 1x,'thdifmax, thdifmin, thdifave:',1X,3(e14.7,1x),"// &
      " 1x,'thtdifmax, thtdifmin, thtdifave:',1X,3(e14.7,1x))"

  CHARACTER(LEN=*),PARAMETER :: fm2013 = &
        "(1x,' rhdifmax,  rhdifmin,  rhdifave:',1X,3(e14.7,1x),"// &
        " 1x,'pexdifmax, pexdifmin, pexdifave:',1X,3(e14.7,1x))"

  REAL(KIND=euwp) ::     &
       udifmx,  udifmn,  udifav, &
       vdifmx  ,vdifmn,  vdifav, &
      thdifmx, thdifmn, thdifav, &
     thtdifmx,thtdifmn,thtdifav, &
      rhdifmx, rhdifmn, rhdifav, &
     pexdifmx,pexdifmn,pexdifav
  INTEGER,OPTIONAL :: isnap

    CALL ttbeg(39)
  ! Calc averages, max and mins of 3D variables 
  !----------------------------------------------------------------------!
    scr1(:,:,:)=uref(:,:,:)-u(:,:,:)
    CALL calcavgmaxmin(scr1,udifmx,udifmn,udifav,np,mp,lp,ih)
    scr1(:,:,:)=vref(:,:,:)-v(:,:,:)
    CALL calcavgmaxmin(scr1,vdifmx,vdifmn,vdifav,np,mp,lp,ih)
    scr1(:,:,:)=thref(:,:,:)-th(:,:,:)-the(:,:,:)
    CALL calcavgmaxmin(scr1,thdifmx,thdifmn,thdifav,np,mp,lp,ih)
    scr1(:,:,:)=thref(:,:,:)-tht(:,:,:)
    CALL calcavgmaxmin(scr1,thtdifmx,thtdifmn,thtdifav,np,mp,lp,ih)
    IF (PRESENT(rhref).AND.PRESENT(pextref)) THEN
      scr1(:,:,:)=rhref(:,:,:)-rho(:,:,:)
      CALL calcavgmaxmin(scr1,rhdifmx,rhdifmn,rhdifav,np,mp,lp,ih)
      scr1(:,:,:)=pextref(:,:,:)-pext(:,:,:)
      CALL calcavgmaxmin(scr1,pexdifmx,pexdifmn,pexdifav,np,mp,lp,ih)
    ENDIF

    CALL wrttobuf(udifmx,udifmn,udifav,bufavg,bufmax,bufmin,1)
    CALL wrttobuf(vdifmx,vdifmn,vdifav,bufavg,bufmax,bufmin,2)
    CALL wrttobuf(thdifmx,thdifmn,thdifav,bufavg,bufmax,bufmin,3)
    CALL wrttobuf(thtdifmx,thtdifmn,thtdifav,bufavg,bufmax,bufmin,4)
    IF (PRESENT(rhref).AND.PRESENT(pextref)) THEN
      CALL wrttobuf( rhdifmx, rhdifmn, rhdifav,bufavg,bufmax,bufmin,5)
      CALL wrttobuf(pexdifmx,pexdifmn,pexdifav,bufavg,bufmax,bufmin,6)
    ENDIF

  !Compute sums, maxs and mins in one global operation
  CALL globsumaxminv(bufavg,bufmax,bufmin,6_iintegers)

  !Retrive GLOBAL averages,maxima and minima from the buffer
    CALL wrtfmbuf(udifmx,udifmn,udifav,bufavg,bufmax,bufmin,1)
    CALL wrtfmbuf(vdifmx,vdifmn,vdifav,bufavg,bufmax,bufmin,2)
    CALL wrtfmbuf(thdifmx,thdifmn,thdifav,bufavg,bufmax,bufmin,3)
    CALL wrtfmbuf(thtdifmx,thtdifmn,thtdifav,bufavg,bufmax,bufmin,4)
    IF (PRESENT(rhref).AND.PRESENT(pextref)) THEN
      CALL wrtfmbuf( rhdifmx, rhdifmn, rhdifav,bufavg,bufmax,bufmin,5)
      CALL wrtfmbuf(pexdifmx,pexdifmn,pexdifav,bufavg,bufmax,bufmin,6)
    ENDIF
  
  IF (mype == 0) THEN
    WRITE(6,fm2010)
      PRINT fm2012,       &
      udifmx, udifmn, udifav,   & 
      vdifmx ,vdifmn ,vdifav,   &
      thdifmx,thdifmn,thdifav,  &
      thtdifmx,thtdifmn,thtdifav
      IF(PRESENT(rho).AND.PRESENT(pext)) THEN
        PRINT fm2013,     &
         rhdifmx, rhdifmn, rhdifav,&
        pexdifmx,pexdifmn,pexdifav
      ENDIF
    WRITE(6,fm2010)
  ENDIF
  IF(PRESENT(isnap).AND.iwrite == 1) THEN
#ifdef PNETCDF
        IF(isnap.eq.1) THEN
          CALL pnet_out_chunk('pextref  ','tstn1.nc',1,1,1,1,0,0,(rg/cmpex)*(spexi*pextref)**capi,np,mp,lp,ih)
          CALL pnet_out_chunk('thref    ','tstn1.nc',1,1,1,1,0,0,thref,np,mp,lp,ih)
          CALL pnet_out_chunk('th       ','tstn1.nc',1,1,1,1,0,0,th,np,mp,lp,ih)
          CALL pnet_out_chunk('tht      ','tstn1.nc',1,1,1,1,0,0,tht,np,mp,lp,ih)
          CALL pnet_out_chunk('rh       ','tstn1.nc',1,1,1,1,0,0,rho,np,mp,lp,ih)
          CALL pnet_out_chunk('u        ','tstn1.nc',1,1,1,1,0,0,u,np,mp,lp,ih)
          CALL pnet_out_chunk('v        ','tstn1.nc',1,1,1,1,0,0,v,np,mp,lp,ih)
          CALL pnet_out_chunk('pext     ','tstn1.nc',1,1,1,1,0,0,(rg/cmpex)*(spexi*pextref)**capi,np,mp,lp,ih)
        ELSE IF (isnap.eq.2) THEN
          CALL pnet_out_chunk('pextref  ','tstn2.nc',1,1,1,1,0,0,(rg/cmpex)*(spexi*pextref)**capi,np,mp,lp,ih)
          CALL pnet_out_chunk('thref    ','tstn2.nc',1,1,1,1,0,0,thref,np,mp,lp,ih)
          CALL pnet_out_chunk('th       ','tstn2.nc',1,1,1,1,0,0,th,np,mp,lp,ih)
          CALL pnet_out_chunk('tht      ','tstn2.nc',1,1,1,1,0,0,tht,np,mp,lp,ih)
          CALL pnet_out_chunk('rh       ','tstn2.nc',1,1,1,1,0,0,rho,np,mp,lp,ih)
          CALL pnet_out_chunk('u        ','tstn2.nc',1,1,1,1,0,0,u,np,mp,lp,ih)
          CALL pnet_out_chunk('v        ','tstn2.nc',1,1,1,1,0,0,v,np,mp,lp,ih)
          CALL pnet_out_chunk('pext     ','tstn2.nc',1,1,1,1,0,0,(rg/cmpex)*(spexi*pextref)**capi,np,mp,lp,ih)
        ELSE IF (isnap.eq.3) THEN
          CALL pnet_out_chunk('pextref ','tstn3.nc',1,1,1,1,0,0,(rg/cmpex)*(spexi*pextref)**capi,np,mp,lp,ih)
          CALL pnet_out_chunk('thref    ','tstn3.nc',1,1,1,1,0,0,thref,np,mp,lp,ih)
          CALL pnet_out_chunk('th       ','tstn3.nc',1,1,1,1,0,0,th,np,mp,lp,ih)
          CALL pnet_out_chunk('tht      ','tstn3.nc',1,1,1,1,0,0,tht,np,mp,lp,ih)
          CALL pnet_out_chunk('rh       ','tstn3.nc',1,1,1,1,0,0,rho,np,mp,lp,ih)
          CALL pnet_out_chunk('u        ','tstn3.nc',1,1,1,1,0,0,u,np,mp,lp,ih)
          CALL pnet_out_chunk('v        ','tstn3.nc',1,1,1,1,0,0,v,np,mp,lp,ih)
          CALL pnet_out_chunk('pext     ','tstn3.nc',1,1,1,1,0,0,(rg/cmpex)*(spexi*pextref)**capi,np,mp,lp,ih)
        ENDIF
#endif 
  ENDIF
  
    CALL ttend(39)
END SUBROUTINE diagnos_difference

!----------------------------------------------------------------------!
SUBROUTINE diagnos_advection(ox,oy,oz,rho,dxi,dyi,dzi,dt,np,mp,lp,ih)
!----------------------------------------------------------------------!
  USE,intrinsic :: ieee_arithmetic
  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv,update3,iup
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin
  USE bconditions, ONLY: remove_cyclic_offset_full,cp_scnd_last_to_halo_xyz_full
  USE mod_parameters, ONLY: nml,ibcx,ibcy,ibcz,ipoles


  !----------------------------------------------------------------------!

  INTEGER(KIND=iintegers) ::  np,mp,lp,ih
  REAL_euwp,INTENT(INOUT) :: &
     ox(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2), &
     oy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2), &
     oz(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)
  REAL_euwp,INTENT(IN) :: &
     rho(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp),INTENT(IN) :: &
      dxi,dyi,dzi,dt
  REAL(KIND=euwp) :: &
      bufavg(6), &
      bufmax(6), &
      bufmin(6)
  REAL(KIND=euwp) :: gc1,gc2,gc3
 
  INTEGER(KIND=iintegers) ::  i, j, k

  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
      "(1x,' Courant_h_max, Courant_h_min, Courant_h_ave:',1X,3(e14.7,1x),"// &
      " 1x,' Courant_v_max, Courant_v_min, Courant_v_ave:',1X,3(e14.7,1x),"// &
      " 1x,' Courant_3D_max, Courant_3D_min, Courant_3D_ave:',1X,3(e14.7,1x))"

  CHARACTER(LEN=*),PARAMETER :: fm2013 = &
      "(1x,' Lipschitz_h_max, Lipschitz_h_min, Lipschitz_h_ave:',1X,3(e14.7,1x),"// &
      " 1x,' Lipschitz_v_max, Lipschitz_v_min, Lipschitz_v_ave:',1X,3(e14.7,1x),"// &
      " 1x,' Lipschitz_3D_max, Lipschitz_3D_min, Lipschitz_3D_ave:',1X,3(e14.7,1x))"

  REAL(KIND=euwp) ::     &
       crhmx,  crhmn, crhav, &
       crvmx,  crvmn, crvav, &
       cr3mx,  cr3mn, cr3av, &
       lphmx,  lphmn, lphav, &
       lpvmx,  lpvmn, lpvav, &
       lp3mx,  lp3mn, lp3av, &
       lpmxloc, lpmnloc,lpavloc,crloc,dinv


    CALL ttbeg(39)
    gc1=dt*dxi
    gc2=dt*dyi
    gc3=dt*dzi

    crhmx=-HUGE(crhmx)
    crvmx=-HUGE(crvmx)
    cr3mx=-HUGE(cr3mx)
    crhmn=HUGE(crhmn)
    crvmn=HUGE(crvmn)
    cr3mn=HUGE(cr3mn)
    crhav=0._euwp
    crvav=0._euwp
    cr3av=0._euwp
    lphmx=-HUGE(lphmx)
    lpvmx=-HUGE(lpvmx)
    lp3mx=-HUGE(lp3mx)
    lphmn=HUGE(lphmn)
    lpvmn=HUGE(lpvmn)
    lp3mn=HUGE(lp3mn)
    lphav=0._euwp
    lpvav=0._euwp
    lp3av=0._euwp
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
          dinv=1._euwp/rho(i,j,k)
         
          crloc=  (abs(ox(i,j,k,1)) +           &
                   abs(oy(i,j,k,1)) )*dinv
          crhmx=  max(crhmx,crloc)
          crhmn=  min(crhmn,crloc)
          crhav=  crhav+crloc

          crloc=  (abs(oz(i,j,k,1)))*dinv
          crvmx=  max(crvmx,crloc)
          crvmn=  min(crvmn,crloc)
          crvav=  crvav+crloc

          crloc=  (abs(ox(i,j,k,1)) +            &
                   abs(oy(i,j,k,1)) +            &
                   abs(oz(i,j,k,1)))*dinv
          cr3mx=  max(cr3mx,crloc)
          cr3mn=  min(cr3mn,crloc)
          cr3av=  cr3av+crloc

          lpmxloc= max(                                                     &
             .5_euwp*gc1/gc1*abs(ox(i+1,j  ,k  ,1)-ox(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc1*abs(ox(i  ,j+1,k  ,1)-ox(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc1*abs(ox(i  ,j  ,k+1,1)-ox(i  ,j  ,k-1,1))*dinv,  &
             .5_euwp*gc1/gc2*abs(oy(i+1,j  ,k  ,1)-oy(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc2*abs(oy(i  ,j+1,k  ,1)-oy(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc2*abs(oy(i  ,j  ,k+1,1)-oy(i  ,j  ,k-1,1))*dinv)  
          lpmnloc= min(                                                     &
             .5_euwp*gc1/gc1*abs(ox(i+1,j  ,k  ,1)-ox(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc1*abs(ox(i  ,j+1,k  ,1)-ox(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc1*abs(ox(i  ,j  ,k+1,1)-ox(i  ,j  ,k-1,1))*dinv,  &
             .5_euwp*gc1/gc2*abs(oy(i+1,j  ,k  ,1)-oy(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc2*abs(oy(i  ,j+1,k  ,1)-oy(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc2*abs(oy(i  ,j  ,k+1,1)-oy(i  ,j  ,k-1,1))*dinv)  
          lpavloc= (1._euwp/6._euwp)*(                                      &
             .5_euwp*gc1/gc1*abs(ox(i+1,j  ,k  ,1)-ox(i-1,j  ,k  ,1))*dinv+  &
             .5_euwp*gc2/gc1*abs(ox(i  ,j+1,k  ,1)-ox(i  ,j-1,k  ,1))*dinv+  &
             .5_euwp*gc3/gc1*abs(ox(i  ,j  ,k+1,1)-ox(i  ,j  ,k-1,1))*dinv+  &
             .5_euwp*gc1/gc2*abs(oy(i+1,j  ,k  ,1)-oy(i-1,j  ,k  ,1))*dinv+  &
             .5_euwp*gc2/gc2*abs(oy(i  ,j+1,k  ,1)-oy(i  ,j-1,k  ,1))*dinv+  &
             .5_euwp*gc3/gc2*abs(oy(i  ,j  ,k+1,1)-oy(i  ,j  ,k-1,1))*dinv)  
          lphmx=  max(lphmx,lpmxloc)
          lphmn=  min(lphmn,lpmnloc)
          lphav=  lphav+lpavloc
          if(ieee_is_nan(lpavloc)) then
            print '(A,3I5)','NaN detected at i,j,k:',i,j,k
            print '(A,2F7.2)','oxip,oxim',ox(i+1,j  ,k  ,1),ox(i-1,j  ,k  ,1)
            print '(A,2F7.2)','oyip,oyim',oy(i+1,j  ,k  ,1),oy(i-1,j  ,k  ,1)
            print '(A,2F7.2)','oxjp,oxjm',ox(i  ,j+1,k  ,1),ox(i  ,j-1,k  ,1)
            print '(A,2F7.2)','oyjp,oyjm',oy(i  ,j+1,k  ,1),ox(i  ,j-1,k  ,1)
            print '(A,2F7.2)','oxkp,oxkm',ox(i  ,j  ,k+1,1),ox(i  ,j  ,k-1,1)
            print '(A,2F7.2)','oykp,oykm',oy(i  ,j  ,k+1,1),oy(i  ,j  ,k-1,1)
            STOP
          endif

          lpmxloc= max(                                                     &
             .5_euwp*gc1/gc3*abs(oz(i+1,j  ,k  ,1)-oz(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc3*abs(oz(i  ,j+1,k  ,1)-oz(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc3*abs(oz(i  ,j  ,k+1,1)-oz(i  ,j  ,k-1,1))*dinv )
          lpmnloc= min(                                                     &
             .5_euwp*gc1/gc3*abs(oz(i+1,j  ,k  ,1)-oz(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc3*abs(oz(i  ,j+1,k  ,1)-oz(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc3*abs(oz(i  ,j  ,k+1,1)-oz(i  ,j  ,k-1,1))*dinv )
          lpavloc= (1._euwp/3._euwp)*(                                      &
             .5_euwp*gc1/gc3*abs(oz(i+1,j  ,k  ,1)-oz(i-1,j  ,k  ,1))*dinv+  &
             .5_euwp*gc2/gc3*abs(oz(i  ,j+1,k  ,1)-oz(i  ,j-1,k  ,1))*dinv+  &
             .5_euwp*gc3/gc3*abs(oz(i  ,j  ,k+1,1)-oz(i  ,j  ,k-1,1))*dinv )

          lpvmx=  max(lpvmx,lpmxloc)
          lpvmn=  min(lpvmn,lpmnloc)
          lpvav=  lpvav+lpavloc

          lpmxloc= max(                                                     &
             .5_euwp*gc1/gc1*abs(ox(i+1,j  ,k  ,1)-ox(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc1*abs(ox(i  ,j+1,k  ,1)-ox(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc1*abs(ox(i  ,j  ,k+1,1)-ox(i  ,j  ,k-1,1))*dinv,  &
             .5_euwp*gc1/gc2*abs(oy(i+1,j  ,k  ,1)-oy(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc2*abs(oy(i  ,j+1,k  ,1)-oy(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc2*abs(oy(i  ,j  ,k+1,1)-oy(i  ,j  ,k-1,1))*dinv,  &
             .5_euwp*gc1/gc3*abs(oz(i+1,j  ,k  ,1)-oz(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc3*abs(oz(i  ,j+1,k  ,1)-oz(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc3*abs(oz(i  ,j  ,k+1,1)-oz(i  ,j  ,k-1,1))*dinv )
          lpmnloc= min(                                                     &
             .5_euwp*gc1/gc1*abs(ox(i+1,j  ,k  ,1)-ox(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc1*abs(ox(i  ,j+1,k  ,1)-ox(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc1*abs(ox(i  ,j  ,k+1,1)-ox(i  ,j  ,k-1,1))*dinv,  &
             .5_euwp*gc1/gc2*abs(oy(i+1,j  ,k  ,1)-oy(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc2*abs(oy(i  ,j+1,k  ,1)-oy(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc2*abs(oy(i  ,j  ,k+1,1)-oy(i  ,j  ,k-1,1))*dinv,  &
             .5_euwp*gc1/gc3*abs(oz(i+1,j  ,k  ,1)-oz(i-1,j  ,k  ,1))*dinv,  &
             .5_euwp*gc2/gc3*abs(oz(i  ,j+1,k  ,1)-oz(i  ,j-1,k  ,1))*dinv,  &
             .5_euwp*gc3/gc3*abs(oz(i  ,j  ,k+1,1)-oz(i  ,j  ,k-1,1))*dinv )
          lpavloc= (1._euwp/9._euwp)*(                                      &
             .5_euwp*gc1/gc1*abs(ox(i+1,j  ,k  ,1)-ox(i-1,j  ,k  ,1))*dinv+  &
             .5_euwp*gc2/gc1*abs(ox(i  ,j+1,k  ,1)-ox(i  ,j-1,k  ,1))*dinv+  &
             .5_euwp*gc3/gc1*abs(ox(i  ,j  ,k+1,1)-ox(i  ,j  ,k-1,1))*dinv+  &
             .5_euwp*gc1/gc2*abs(oy(i+1,j  ,k  ,1)-oy(i-1,j  ,k  ,1))*dinv+  &
             .5_euwp*gc2/gc2*abs(oy(i  ,j+1,k  ,1)-oy(i  ,j-1,k  ,1))*dinv+  &
             .5_euwp*gc3/gc2*abs(oy(i  ,j  ,k+1,1)-oy(i  ,j  ,k-1,1))*dinv+  &
             .5_euwp*gc1/gc3*abs(oz(i+1,j  ,k  ,1)-oz(i-1,j  ,k  ,1))*dinv+  &
             .5_euwp*gc2/gc3*abs(oz(i  ,j+1,k  ,1)-oz(i  ,j-1,k  ,1))*dinv+  &
             .5_euwp*gc3/gc3*abs(oz(i  ,j  ,k+1,1)-oz(i  ,j  ,k-1,1))*dinv )
          lp3mx=  max(lp3mx,lpmxloc)
          lp3mn=  min(lp3mn,lpmnloc)
          lp3av=  lp3av+lpavloc
        END DO
      END DO
    END DO
    crhav=crhav/nml
    crvav=crvav/nml
    cr3av=cr3av/nml
    lphav=lphav/nml
    lpvav=lpvav/nml
    lp3av=lp3av/nml

    CALL wrttobuf(crhmx,crhmn,crhav,bufavg,bufmax,bufmin,1)
    CALL wrttobuf(crvmx,crvmn,crvav,bufavg,bufmax,bufmin,2)
    CALL wrttobuf(cr3mx,cr3mn,cr3av,bufavg,bufmax,bufmin,3)
    CALL wrttobuf(lphmx,lphmn,lphav,bufavg,bufmax,bufmin,4)
    CALL wrttobuf(lpvmx,lpvmn,lpvav,bufavg,bufmax,bufmin,5)
    CALL wrttobuf(lp3mx,lp3mn,lp3av,bufavg,bufmax,bufmin,6)

  !Compute sums, maxs and mins in one global operation
  CALL globsumaxminv(bufavg,bufmax,bufmin,6_iintegers)

    CALL wrtfmbuf(crhmx,crhmn,crhav,bufavg,bufmax,bufmin,1)
    CALL wrtfmbuf(crvmx,crvmn,crvav,bufavg,bufmax,bufmin,2)
    CALL wrtfmbuf(cr3mx,cr3mn,cr3av,bufavg,bufmax,bufmin,3)
    CALL wrtfmbuf(lphmx,lphmn,lphav,bufavg,bufmax,bufmin,4)
    CALL wrtfmbuf(lpvmx,lpvmn,lpvav,bufavg,bufmax,bufmin,5)
    CALL wrtfmbuf(lp3mx,lp3mn,lp3av,bufavg,bufmax,bufmin,6)

  
  IF (mype == 0) THEN
    WRITE(6,fm2010)
      PRINT fm2012,       &
      crhmx, crhmn, crhav,   & 
      crvmx, crvmn, crvav,   & 
      cr3mx, cr3mn, cr3av

    WRITE(6,fm2010)

      PRINT fm2013,       &
      lphmx, lphmn, lphav,   & 
      lpvmx, lpvmn, lpvav,   & 
      lp3mx, lp3mn, lp3av
    WRITE(6,fm2010)
  ENDIF
  
  
  CALL ttend(39)
END SUBROUTINE diagnos_advection

!----------------------------------------------------------------------!
SUBROUTINE diagnos_advection_vz(ox,oy,oz,rho,np,mp,lp,ih)
!----------------------------------------------------------------------!

  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin
  USE mod_parameters, ONLY: n,m

  !----------------------------------------------------------------------!

  INTEGER(KIND=iintegers) ::  np,mp,lp,ih
  REAL(KIND=euwp),INTENT(INOUT) :: &
     ox(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2), &
     oy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2), &
     oz(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)
  REAL(KIND=euwp),INTENT(IN) :: &
     rho(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp) :: &
      bufavg(3*lp), &
      bufmax(3*lp), &
      bufmin(3*lp)
 
  REAL(KIND=euwp) cr_frac
  INTEGER(KIND=iintegers) ::  i, j, k

  CHARACTER(LEN=*),PARAMETER :: fm2010 = &
  "(1x,'---------------------------------------------------------')"
  CHARACTER(LEN=*),PARAMETER :: fm2012 = &
      "(1x,' Courant_h_max, Courant_h_min, Courant_h_ave:',1X,3(e14.7,1x),"// &
      " 1x,' Courant_v_max, Courant_v_min, Courant_v_ave:',1X,3(e14.7,1x),"// &
      "  1x,' Courant_3D_max, Courant_3D_min, Courant_3D_ave:',1X,3(e14.7,1x))"

  REAL(KIND=euwp) ::     &
       crhmx(lp),  crhmn(lp), crhav(lp), &
       crvmx(lp),  crvmn(lp), crvav(lp), &
       cr3mx(lp),  cr3mn(lp), cr3av(lp), &
       crloc,dinv


    CALL ttbeg(39)
    crhmx(:)=-HUGE(crhmx)
    crvmx(:)=-HUGE(crvmx)
    cr3mx(:)=-HUGE(cr3mx)
    crhmn(:)=HUGE(crhmn)
    crvmn(:)=HUGE(crvmn)
    cr3mn(:)=HUGE(cr3mn)
    crhav(:)=0._euwp
    crvav(:)=0._euwp
    cr3av(:)=0._euwp

    DO k=1,lp
      DO j=1,mp
        DO i=1,np
          dinv=1._euwp/rho(i,j,k)
         
          crloc=  (abs(ox(i,j,k,1)) +           &
                   abs(oy(i,j,k,1)) )*dinv
          crhmx(k)=  max(crhmx(k),crloc)
          crhmn(k)=  min(crhmn(k),crloc)
          crhav(k)=  crhav(k)+crloc

          crloc=  (abs(oz(i,j,k,1)))*dinv
          crvmx(k)=  max(crvmx(k),crloc)
          crvmn(k)=  min(crvmn(k),crloc)
          crvav(k)=  crvav(k)+crloc

          crloc=  (abs(ox(i,j,k,1)) +            &
                   abs(oy(i,j,k,1)) +            &
                   abs(oz(i,j,k,1)))*dinv
          cr3mx(k)=  max(cr3mx(k),crloc)
          cr3mn(k)=  min(cr3mn(k),crloc)
          cr3av(k)=  cr3av(k)+crloc

        END DO
      END DO
    END DO
    crhav(:)=crhav(:)/REAL(n*m)
    crvav(:)=crvav(:)/REAL(n*m)
    cr3av(:)=cr3av(:)/REAL(n*m)
    DO k=1,lp
      CALL wrttobuf(crhmx(k),crhmn(k),crhav(k),bufavg,bufmax,bufmin,(k-1)*3+1)
      CALL wrttobuf(crvmx(k),crvmn(k),crvav(k),bufavg,bufmax,bufmin,(k-1)*3+2)
      CALL wrttobuf(cr3mx(k),cr3mn(k),cr3av(k),bufavg,bufmax,bufmin,(k-1)*3+3)
    ENDDO

  !Compute sums, maxs and mins in one global operation
  CALL globsumaxminv(bufavg,bufmax,bufmin,lp*3_iintegers)

    DO k=1,lp
      CALL wrtfmbuf(crhmx(k),crhmn(k),crhav(k),bufavg,bufmax,bufmin,(k-1)*3+1)
      CALL wrtfmbuf(crvmx(k),crvmn(k),crvav(k),bufavg,bufmax,bufmin,(k-1)*3+2)
      CALL wrtfmbuf(cr3mx(k),cr3mn(k),cr3av(k),bufavg,bufmax,bufmin,(k-1)*3+3)
    ENDDO

  
  IF (mype == 0.AND.(MAXVAL(crhmx) > 1. .OR.&
                     MAXVAL(crvmx) > 1. .OR.&
                     MAXVAL(cr3mx) > 1       )) THEN
    WRITE(6,fm2010)
    DO k=1,lp
      cr_frac=0.85
      IF(crhmx(k).gt.cr_frac.OR.crvmx(k).GT.cr_frac.OR.cr3mx(k).GT.cr_frac) THEN
      PRINT *,'Level k= ',k
      PRINT fm2012,       &
      crhmx(k), crhmn(k), crhav(k),   & 
      crvmx(k), crvmn(k), crvav(k),   & 
      cr3mx(k), cr3mn(k), cr3av(k)
      ENDIF
    ENDDO

    WRITE(6,fm2010)

  ENDIF
  
  
  CALL ttend(39)
END SUBROUTINE diagnos_advection_vz
!----------------------------------------------------------------------!
! Transform solver pressure/exner perturbation to Pa
!----------------------------------------------------------------------!
SUBROUTINE transform_exner_to_pascals(pvar_in,pPa_out,ppe,np,mp,lp,ih)
USE mod_parameters, ONLY: spexi,rg,cmpex,capi
!------------------------------------------------------------------

  !--------------------------------------------------------------------!

  INTEGER(KIND=iintegers), INTENT(IN) :: np,mp,lp,ih
  REAL(KIND=euwp), INTENT(IN) :: &
   pvar_in(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), & 
   ppe(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) 

  REAL(KIND=euwp), INTENT(OUT) :: &
   pPa_out(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) 
  
!  REAL(KIND=euwp) pexf,pexamb,pprt,pful,thetot
  REAL(KIND=euwp) pexf,pprt,pful

  INTEGER(KIND=iintegers) :: i, j, k

  !--------------------------------------------------------------------!

  DO k=1,lp
    DO j=1,mp
      DO i=1,np
!        IF (icmprss == 1 .OR. icmpexpl == 1) THEN !convert exner to pressure
!          IF (icmprss*(1-icmpexpl) == 1) THEN
            pexf=pvar_in(i,j,k)*spexi
            !thermodynamic pressure
            ! pexf=(rho(i,j,k,0)*tht(i,j,k)*cmpex)**wexnr
!          ELSE
!            thetot=the(i,j,k)
!            ! thetot=thetot*(1.+qve(i,j,k)*epsi)
!            pexamb=(rhe(i,j,k)*thetot*cmpex)**wexnr !Exner ambient
!            pexf=p(i,j,k)*spexi+pexamb              !full Exner
!          END IF
          pful=(pexf**capi)*(rg/cmpex)
          pprt=pful-ppe(i,j,k)!-pprim(i,j,k)
          pPa_out(i,j,k)=pprt
!        ELSE
!          tempp(i,j,k)=2*dti*p(i,j,k)*rh0(i,j,k)
!        END IF !iexner
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE transform_exner_to_pascals

SUBROUTINE exam_var(var,vname,next,np,mp,lp,ih)
  USE mpi_parallel, ONLY: globsumaxminlocv,globsum,nsubpos,msubpos,lsubpos,mype
  USE mod_parameters, ONLY: nml
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin
  USE eulag_diagutils, ONLY: calcsumaxminloc
!  USE eulag_datafields, ONLY: zcr,g13,g23,g33
 
      INTEGER(KIND=iintegers),INTENT(IN) ::  np,mp,lp,ih
      INTEGER,PARAMETER ::iflo=6 ! Which file to output exam_var to, 6 - screen
!                         20 - symmetry, 21 - reports 
      REAL(KIND=euwp)   var(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      REAL(KIND=euwp) tempmx(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      REAL(KIND=euwp) tempmn(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      INTEGER(KIND=iintegers),PARAMETER :: idmaxvars=40
      REAL(KIND=euwp) buffer(15,idmaxvars)
      REAL(KIND=euwp) buf1din(idmaxvars)
      REAL(KIND=euwp) buf1dout(idmaxvars)
      CHARACTER(5) bufnames(100) 
      CHARACTER(5) vname5 
      CHARACTER(LEN=*) vname
      REAL(KIND=euwp) xsd,xmxcnt,xmncnt
      INTEGER i,j,k,next,ia,ja,ka,ilm,ium,jlm,jum,klm,kum
      buffer(1,:)=0
      buffer(2,:)=-HUGE(buffer)
      buffer(3,:)=HUGE(buffer)
      buffer(4:9,:)=0.
      vname5=TRIM(ADJUSTL(vname(1:5)))
      IF(mype.eq.0) then
        WRITE(iflo,*) '------------------------------------'
        WRITE(iflo,*) 'Examining variable ',TRIM(ADJUSTL(vname(1:5))) 
        WRITE(iflo,*) '------------------------------------'
      ENDIF
!************************
!First, examine variable average, maximum, minimum and their absolute locations
!************************

      CALL  calcsumaxminloc(var,vname5,buffer,bufnames,1,np,mp,lp,ih)
      CALL  globsumaxminlocv(buffer,1)

!************************
! Now compute standard deviation and number of occurences of maximum and minimum field value 
!************************
      xsd=0.
      xmxcnt=0.
      xmncnt=0.
      do k=1,lp  
      do j=1,mp  
      do i=1,np  
        xsd=xsd+(var(i,j,k)-buffer(1,1))**2
        if(var(i,j,k).eq.buffer(2,1)) xmxcnt=xmxcnt+1.
        if(var(i,j,k).eq.buffer(3,1)) xmncnt=xmncnt+1.
      enddo
      enddo
      enddo
      buf1din(1)=xsd
      buf1din(2)=xmxcnt
      buf1din(3)=xmncnt
      call  globsum(buf1din,buf1dout,3)
      xsd=sqrt(buf1dout(1)/nml)
      xmxcnt=buf1dout(2)
      xmncnt=buf1dout(3)

      IF(mype.eq.0) THEN
         WRITE(iflo,*)                                                &
        'Variable name  ave std  max at  ia,ja,ka  min at ia,ja,ka'
         WRITE(iflo,                                                  &
      '(A5," ave",E12.4," max",E12.4," std",E12.4,'//                 &
      '    " at",3I4," min",E12.4," at",3I4)')                        &
         bufnames(1),buffer(1,1),buffer(2,1),xsd,                     & ! variable name, average and maximum value
         NINT(buffer(4,1)),NINT(buffer(5,1)),NINT(buffer(6,1)),       &! absolute location of maximum value
         buffer(3,1),                                                 &   ! minimum value
         NINT(buffer(10,1)),NINT(buffer(11,1)),NINT(buffer(12,1))  ! absolute location of minimum value
         WRITE(iflo,'("Max value count: ",I0," Min value count ",I0)')  &
               NINT(buf1dout(2)),NINT(buf1dout(3)) 
         WRITE(iflo,*)                                                &
        'Variable name  max at  ia,ja,ka,i,j,k,npos,mpos,lpos'! topo zcr g13 g23 g33
         WRITE(iflo,                                                  &
      '(A5," max",E12.4," at",3I4,"   ",3I4,"   ",3I4)')                        &
         bufnames(1),buffer(2,1),                     & ! variable name, average and maximum value
         NINT(buffer(4,1)),NINT(buffer(5,1)),NINT(buffer(6,1)) ,       &! absolute location of maximum value
         NINT(buffer(7,1)),NINT(buffer(8,1)),NINT(buffer(9,1)) ,       &! relative location of maximum value
         NINT((buffer(4,1)-buffer(7,1))/np),                           &! npos
         NINT((buffer(5,1)-buffer(8,1))/mp),                           &! mpos
         NINT((buffer(6,1)-buffer(9,1))/lp)                             ! lpos

      ENDIF

!      IF((buffer(4,1).gt.nsubpos).AND.(buffer(4,1).le.nsubpos+np).AND.  &
!         (buffer(5,1).gt.msubpos).AND.(buffer(5,1).le.msubpos+mp)     ) THEN
!         WRITE(iflo,*)                                                &
!        'Variable name topo zcr g13 g23 g33'
!         WRITE (iflo,'(A5,5E12.4)')                                   & 
!         bufnames(1),                                                 &
!         zcr(NINT(buffer(4,1))-nsubpos,NINT(buffer(5,1))-msubpos,1),  &   
!         zcr(NINT(buffer(4,1))-nsubpos,NINT(buffer(5,1))-msubpos,NINT(buffer(6,1))-lsubpos),&   
!    atan(g13(NINT(buffer(4,1))-nsubpos,NINT(buffer(5,1))-msubpos,NINT(buffer(6,1))-lsubpos))*90./acos(-1.),&   
!    atan(g23(NINT(buffer(4,1))-nsubpos,NINT(buffer(5,1))-msubpos,NINT(buffer(6,1))-lsubpos))*90./acos(-1.),&
!         g33(NINT(buffer(4,1))-nsubpos,NINT(buffer(5,1))-msubpos,NINT(buffer(6,1))-lsubpos)
!         WRITE (iflo,'(A5,5E12.4)')                                   & 
!         bufnames(1),                                                 &
!         zcr(NINT(buffer(7,1)),NINT(buffer(8,1)),1),  &   
!         zcr(NINT(buffer(7,1)),NINT(buffer(8,1)),NINT(buffer(9,1))),&   
!    atan(g13(NINT(buffer(7,1)),NINT(buffer(8,1)),NINT(buffer(9,1))))*90./acos(-1.),&   
!    atan(g23(NINT(buffer(7,1)),NINT(buffer(8,1)),NINT(buffer(9,1))))*90./acos(-1.),&
!         g33(NINT(buffer(7,1)),NINT(buffer(8,1)),NINT(buffer(9,1)))
!      ENDIF
   
      vname5='nnnnn' !define special varname to disable symmetry check in calc.

!------------------------------------------------------------------
!Prepare two auxiliary variables to search for next-in-row max and min
!------------------------------------------------------------------

      CALL  calcsumaxminloc(var,vname5,buffer,bufnames,1,np,mp,lp,ih)
      CALL  calcsumaxminloc(var,vname5,buffer,bufnames,2,np,mp,lp,ih)
      CALL  globsumaxminlocv(buffer,2)
      tempmx(1:np,1:mp,1:lp)=var(1:np,1:mp,1:lp)
      tempmn(1:np,1:mp,1:lp)=var(1:np,1:mp,1:lp)

!------------------------------------------------------------------
!Loop over desired number of max/min values
!------------------------------------------------------------------
      DO i=1,next-1
        ia=NINT(buffer(4,1))
        ja=NINT(buffer(5,1))
        ka=NINT(buffer(6,1))
        ilm=nsubpos
        ium=nsubpos+np
        jlm=msubpos
        jum=msubpos+mp
        klm=lsubpos
        kum=lsubpos+lp
        IF(ia.gt.ilm.and.ia.le.ium.and.    & !Check if we are on the right processor
           ja.gt.jlm.and.ja.le.jum.and.    &
           ka.gt.klm.and.ka.le.kum)  then
           tempmx(ia-ilm,ja-jlm,ka-klm)=-1.e35  !remove maximum value to search for the next
        ENDIF
        CALL  calcsumaxminloc(tempmx,vname5,buffer,bufnames,1,np,mp,lp,ih) !determine next max value

        ia=NINT(buffer(7,2))
        ja=NINT(buffer(8,2))
        ka=NINT(buffer(9,2))
        if(ia.gt.ilm.and.ia.le.ium.and.    &
           ja.gt.jlm.and.ja.le.jum.and.    &
           ka.gt.klm.and.ka.le.kum)  then
           tempmn(ia-ilm,ja-jlm,ka-klm)=1.e35 !remove minimum value to search for the next
        endif
        CALL  calcsumaxminloc(tempmn,vname5,buffer,bufnames,2,np,mp,lp,ih) !determine next min value
        CALL  globsumaxminlocv(buffer,2)
        IF(mype.eq.0) THEN

         WRITE(iflo,                                                 &
      '("      followed by:   " ," max",E12.4," at",3I4,'//          &
      '                          " min",E12.4," at",3I4)')           &
         buffer(2,1),                                                &
         NINT(buffer(4,1)),NINT(buffer(5,1)),NINT(buffer(6,1)),      &
         buffer(3,2),                                                &
         NINT(buffer(7,2)),NINT(buffer(8,2)),NINT(buffer(9,2))      
         ENDIF
      enddo  !next

      IF(mype.eq.0) THEN
        WRITE(iflo,*) '------------------------------------'
        WRITE(iflo,*) 'End exam variable ',TRIM(ADJUSTL(vname(1:5))) 
        WRITE(iflo,*) '------------------------------------'
      ENDIF
      END SUBROUTINE exam_var

  SUBROUTINE exam_divergence(var,vname,dt,np,mp,lp,ih)
  USE mpi_parallel, ONLY: globsumaxminlocv,globsum,mype
  USE mod_parameters, ONLY: nml
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin
  USE eulag_diagutils, ONLY: calcsumaxminloc
!  USE eulag_datafields, ONLY: zcr,g13,g23,g33
 
      INTEGER(KIND=iintegers),INTENT(IN) ::  np,mp,lp,ih
      INTEGER,PARAMETER ::iflo=6 ! Which file to output exam_var to, 6 - screen
!                         20 - symmetry, 21 - reports 
      REAL(KIND=euwp)   var(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      INTEGER(KIND=iintegers),PARAMETER :: idmaxvars=40
      REAL(KIND=euwp) buffer(15,idmaxvars)
      REAL(KIND=euwp) buf1din(idmaxvars)
      REAL(KIND=euwp) buf1dout(idmaxvars)
      CHARACTER(5) bufnames(100) 
      CHARACTER(5) vname5 
      CHARACTER(LEN=*) vname
      REAL(KIND=euwp) xsd,xmxcnt,xmncnt
      REAL(KIND=euwp) dt 
      INTEGER i,j,k
      buffer(1,:)=0
      buffer(2,:)=-HUGE(buffer)
      buffer(3,:)=HUGE(buffer)
      buffer(4:9,:)=0.
      vname5=TRIM(ADJUSTL(vname(1:5)))
      IF(mype.eq.0) then
        WRITE(iflo,*) '------------------------------------'
        WRITE(iflo,*) 'Eulerian divergence                 ' 
        WRITE(iflo,*) '------------------------------------'
      ENDIF
!************************
!First, examine variable average, maximum, minimum and their absolute locations
!************************

      CALL  calcsumaxminloc(var,vname5,buffer,bufnames,1,np,mp,lp,ih)
      CALL  globsumaxminlocv(buffer,1)

!************************
! Now compute standard deviation and number of occurences of maximum and minimum field value 
!************************
      xsd=0.
      xmxcnt=0.
      xmncnt=0.
      do k=1,lp  
      do j=1,mp  
      do i=1,np  
        xsd=xsd+(var(i,j,k)-buffer(1,1))**2
        if(var(i,j,k).eq.buffer(2,1)) xmxcnt=xmxcnt+1.
        if(var(i,j,k).eq.buffer(3,1)) xmncnt=xmncnt+1.
      enddo
      enddo
      enddo
      buf1din(1)=xsd
      buf1din(2)=xmxcnt
      buf1din(3)=xmncnt
      call  globsum(buf1din,buf1dout,3)
      xsd=sqrt(buf1dout(1)/nml)
      xmxcnt=buf1dout(2)
      xmncnt=buf1dout(3)

      IF(mype.eq.0) THEN
         WRITE(iflo,*)                                                &
        'Variable name  ave std  max at  ia,ja,ka  min at ia,ja,ka'
         WRITE(iflo,                                                  &
      '(A5," ave",E12.4," max",E12.4," std",E12.4,'//                 &
      '    " at",3I4," min",E12.4," at",3I4)')                        &
         bufnames(1),buffer(1,1)*dt,buffer(2,1)*dt,xsd*dt,            & ! variable name, average and maximum value
         NINT(buffer(4,1)),NINT(buffer(5,1)),NINT(buffer(6,1)),       &! absolute location of maximum value
         buffer(3,1)*dt,                                              &   ! minimum value
         NINT(buffer(10,1)),NINT(buffer(11,1)),NINT(buffer(12,1))  ! absolute location of minimum value
         WRITE(iflo,'("Max value count: ",I0," Min value count ",I0)')  &
               NINT(buf1dout(2)),NINT(buf1dout(3)) 
         WRITE(iflo,*)                                                &
        'Variable name  max at  ia,ja,ka,i,j,k,npos,mpos,lpos'! topo zcr g13 g23 g33
         WRITE(iflo,                                                  &
      '(A5," max",E12.4," at",3I4,"   ",3I4,"   ",3I4)')                        &
         bufnames(1),buffer(2,1)*dt,                     & ! variable name, average and maximum value
         NINT(buffer(4,1)),NINT(buffer(5,1)),NINT(buffer(6,1)) ,       &! absolute location of maximum value
         NINT(buffer(7,1)),NINT(buffer(8,1)),NINT(buffer(9,1)) ,       &! relative location of maximum value
         NINT((buffer(4,1)-buffer(7,1))/np),                           &! npos
         NINT((buffer(5,1)-buffer(8,1))/mp),                           &! mpos
         NINT((buffer(6,1)-buffer(9,1))/lp)                             ! lpos

      ENDIF

      IF(mype.eq.0) THEN
        WRITE(iflo,*) '------------------------------------'
        WRITE(iflo,*) 'End Eulerian divergence             '
        WRITE(iflo,*) '------------------------------------'
      ENDIF
      END SUBROUTINE exam_divergence
  SUBROUTINE exam_errors(vara,vnamea, &
                         varb,vnameb, &
                         varc,vnamec, &
                         vard,vnamed, &
                         vare,vnamee, &
                          dt,np,mp,lp,ih)
  USE mpi_parallel, ONLY: globsumaxminlocv,globsum,mype
  USE mod_parameters, ONLY: nml
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin
  USE eulag_diagutils, ONLY: calcsumaxminloc
!  USE eulag_datafields, ONLY: zcr,g13,g23,g33
 
      INTEGER(KIND=iintegers),INTENT(IN) ::  np,mp,lp,ih
      INTEGER,PARAMETER ::iflo=6 ! Which file to output exam_var to, 6 - screen
!                         20 - symmetry, 21 - reports 
      REAL(KIND=euwp)   vara(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      REAL(KIND=euwp)   varb(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      REAL(KIND=euwp)   varc(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      REAL(KIND=euwp)   vard(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      REAL(KIND=euwp)   vare(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      INTEGER(KIND=iintegers),PARAMETER :: idmaxvars=40
      REAL(KIND=euwp) buffer(15,idmaxvars)
      REAL(KIND=euwp) buf1din(idmaxvars)
      REAL(KIND=euwp) buf1dout(idmaxvars)
      CHARACTER(5) bufnames(100) 
      CHARACTER(5) vname5a,vname5b,vname5c,vname5d,vname5e 
      CHARACTER(LEN=*) vnamea,vnameb,vnamec,vnamed,vnamee
      REAL(KIND=euwp) xsd(5)
      REAL(KIND=euwp) dt 
      INTEGER i,j,k
      INTEGER icnt
      buffer(1,:)=0
      buffer(2,:)=-HUGE(buffer)
      buffer(3,:)=HUGE(buffer)
      buffer(4:9,:)=0.
      vname5a=TRIM(ADJUSTL(vnamea(1:5)))
      vname5b=TRIM(ADJUSTL(vnameb(1:5)))
      vname5c=TRIM(ADJUSTL(vnamec(1:5)))
      vname5d=TRIM(ADJUSTL(vnamed(1:5)))
      vname5e=TRIM(ADJUSTL(vnamee(1:5)))
      IF(mype.eq.0) then
        WRITE(iflo,*) '------------------------------------'
        WRITE(iflo,*) ' GCRK errors diagnostics                 ' 
        WRITE(iflo,*) '------------------------------------'
      ENDIF
!************************
!First, examine variable average, maximum, minimum and their absolute locations
!************************

      CALL  calcsumaxminloc(vara,vname5a,buffer,bufnames,1,np,mp,lp,ih)
      CALL  calcsumaxminloc(varb,vname5b,buffer,bufnames,2,np,mp,lp,ih)
      CALL  calcsumaxminloc(varc,vname5c,buffer,bufnames,3,np,mp,lp,ih)
      CALL  calcsumaxminloc(vard,vname5d,buffer,bufnames,4,np,mp,lp,ih)
      CALL  calcsumaxminloc(vare,vname5e,buffer,bufnames,5,np,mp,lp,ih)
      CALL  globsumaxminlocv(buffer,5)

!************************
! Now compute standard deviation and number of occurences of maximum and minimum field value 
!************************
      xsd(:)=0.
      do k=1,lp  
      do j=1,mp  
      do i=1,np  
        xsd(1)=xsd(1)+(vara(i,j,k)-buffer(1,1))**2
        xsd(2)=xsd(2)+(varb(i,j,k)-buffer(1,2))**2
        xsd(3)=xsd(3)+(varc(i,j,k)-buffer(1,3))**2
        xsd(4)=xsd(4)+(vard(i,j,k)-buffer(1,4))**2
        xsd(5)=xsd(5)+(vare(i,j,k)-buffer(1,5))**2
      enddo
      enddo
      enddo
      buf1din(1:5)=xsd(1:5)
      call  globsum(buf1din,buf1dout,5)
      xsd(1:5)=sqrt(buf1dout(1:5)/nml)

      IF(mype.eq.0) THEN
         WRITE(iflo,*)                                                &
        'Variable name  ave std  max at  ia,ja,ka  min at ia,ja,ka'
        DO icnt=1,5
         WRITE(iflo,                                                  &
      '(A5," ave",E12.4," max",E12.4," std",E12.4,'//                 &
          '" at",3I4," min",E12.4," at",3I4)')                        &
         bufnames(icnt),buffer(1,icnt)*dt,buffer(2,icnt)*dt,xsd(icnt)*dt,     & ! variable name, average and maximum value
         NINT(buffer(4,icnt)),NINT(buffer(5,icnt)),NINT(buffer(6,icnt)),       &! absolute location of maximum value
         buffer(3,icnt)*dt,                                              &   ! minimum value
         NINT(buffer(10,icnt)),NINT(buffer(11,icnt)),NINT(buffer(12,icnt))  ! absolute location of minimum value
        ENDDO
      ENDIF

      IF(mype.eq.0) THEN
        WRITE(iflo,*) '------------------------------------'
        WRITE(iflo,*) 'End GCRK errors diagnostics             '
        WRITE(iflo,*) '------------------------------------'
      ENDIF
      END SUBROUTINE exam_errors

      SUBROUTINE ckcyc(a,varname,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
      USE mpi_parallel, ONLY: mype,globmaxminv,rightedge,topedge,skyedge 
      USE mpi_parallel, ONLY: updatelr, updatebt, updategs
      INTEGER(KIND=iintegers),INTENT(IN) ::  ipoles,np,mp,lp,ih
      REAL_euwp :: a(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
      CHARACTER(5) :: varname
      INTEGER(KIND=iintegers),INTENT(IN):: ibcx,ibcy,ibcz
      INTEGER(KIND=iintegers) iupdate,jx,kx,iy,ky,iz,jz,i,j,k
      REAL(KIND=euwp) ck,cmx,cmn,cm1,cm2,cmxg,cmng
      REAL(KIND=euwp) :: bufmax(1),bufmin(1)
      REAL(KIND=euwp),PARAMETER :: cktre=1.e-8 !14
      iupdate=1
      IF(ipoles.eq.0) THEN
        IF(ibcx.eq.1) THEN
          if(iupdate.eq.1) call updatelr(a,np,mp,lp,np,mp,lp,1,ih)
          cmx =-TINY(cmx)
          cmn = TINY(cmn) 
          cm1 = 0._euwp
          cm2 = 0._euwp
          jx=0
          kx=0
          if (rightedge.eq.1) then
             do k=1,lp
                do j=1,mp
                   ck=a(np,j,k)-a(np+1,j,k)
                   cmx=max(cmx,ck)
                   cmn=min(cmn,ck)
                   if((cmx.eq.ck).or.((cmn.eq.ck))) then
                     jx=j
                     kx=k
                     cm1=a(np  ,j,k)
                     cm2=a(np+1,j,k)
                   endif
                end do
             end do
          end if
          bufmax(1)=cmx
          bufmin(1)=cmn
          call globmaxminv(bufmax,bufmin,1)
          cmxg=bufmax(1)
          cmng=bufmin(1)

          if(cmxg.gt. cktre .or. cmng.lt.-cktre) then
            if((cmxg.eq.cmx).or.(cmng.eq.cmn)) then
              print 100,mype,varname,jx,kx,cmxg,cmng,cm1,cm2
            endif
          endif
        ENDIF !!ibcx

        if(ibcy.eq.1) then
          if(iupdate.eq.1) call updatebt(a,np,mp,lp,np,mp,lp,1,ih)
          cmx=-TINY(cmx)
          cmn= TINY(cmn) 
          iy=0
          ky=0
          cm1 = 0.
          cm2 = 0.
          if (topedge.eq.1) then
             do k=1,lp
                do i=1,np
                   ck=a(i,mp+1,k)-a(i,mp,k)
                   cmx=max(cmx,ck)
                   cmn=min(cmn,ck)
                   if((cmx.eq.ck).or.((cmn.eq.ck))) then
                     iy=i
                     ky=k
                     cm1=a(i,mp,k)
                     cm2=a(i,mp+1,k)
                   endif
                end do
             end do
          end if
          bufmax(1)=cmx
          bufmin(1)=cmn
          call globmaxminv(bufmax,bufmin,1)
          cmxg=bufmax(1)
          cmng=bufmin(1)


          if(cmxg.gt. cktre .or. cmng.lt.-cktre) then
            if((cmxg.eq.cmx).or.(cmng.eq.cmn)) then
              print 200,mype,varname,iy,ky,cmxg,cmng,cm1,cm2
            endif
          endif
        ENDIF !ibcy
      ENDIF !ipoles

      if((ibcz.eq.1)) then
        if(iupdate.eq.1) call updategs(a,np,mp,lp,np,mp,lp,1,ih)
        cmx=-TINY(cmx)
        cmn= TINY(cmn) 
        iz=0
        jz=0
        cm1 = 0.
        cm2 = 0.
        if(skyedge.eq.1) then
          do j=1,mp
             do i=1,np
                ck=a(i,j,lp+1)-a(i,j,lp)
                cmx=max(cmx,ck)
                cmn=min(cmn,ck)
                if((cmx.eq.ck).or.((cmn.eq.ck))) then
                  iz=i
                  jz=j
                  cm1=a(i,j,lp)
                  cm2=a(i,j,lp+1)
                endif
             end do
          end do
        endif
        bufmax(1)=cmx
        bufmin(1)=cmn
        call globmaxminv(bufmax,bufmin,1)
        cmxg=bufmax(1)
        cmng=bufmin(1)

        if(cmxg.gt. cktre .or. cmng.lt.-cktre) then
          if((cmxg.eq.cmx).or.(cmng.eq.cmn)) then
            print 300,mype,varname,iz,jz,cmxg,cmng,cm1,cm2
          endif
        endif
      endif !ibcz
  100 format(1x,i4,1x,a7,' ibcx j,k:',2i4,' mx,mn,a1,a2=',2e11.3,2e11.3)
  200 format(1x,i4,1x,a7,' ibcy i,k:',2i4,' mx,mn,a1,a2=',2e11.3,2e11.3)
  300 format(1x,i4,1x,a7,' ibcz i,j:',2i4,' mx,mn,a1,a2=',2e11.3,2e11.3)

      END SUBROUTINE ckcyc

      SUBROUTINE ckcyc_solver(u0,v0,w0,th,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,tht,pext)
        USE mpi_parallel, ONLY: mype
        INTEGER(KIND=iintegers),INTENT(IN) ::  ipoles,np,mp,lp,ih
        CHARACTER(5) :: varname
        INTEGER(KIND=iintegers),INTENT(IN):: ibcx,ibcy,ibcz
        REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: u0,v0,w0,th
        REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),OPTIONAL :: tht,pext
!       IF(mype == 0) print *,'Periodicity check fter solver'
        CALL flush(6)
        CALL  ckcyc(u0  ,'u0   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        CALL  ckcyc(v0  ,'v0   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        CALL  ckcyc(w0  ,'w0   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        CALL  ckcyc(th  ,'th   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        IF(PRESENT(tht))  CALL ckcyc(tht ,'tht  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        IF(PRESENT(pext)) CALL ckcyc(pext,'pext ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
!       IF(mype == 0) print *,'End of periodicity check after solver'
        CALL flush(6)
 
      END SUBROUTINE ckcyc_solver
      
      SUBROUTINE check_nans_prognostic(rho0,u0,v0,w0,th,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,tht,pext)
        USE mpi_parallel, ONLY: mype
        INTEGER(KIND=iintegers),INTENT(IN) ::  ipoles,np,mp,lp,ih
        CHARACTER(5) :: varname
        INTEGER(KIND=iintegers),INTENT(IN):: ibcx,ibcy,ibcz
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: rho0,u0,v0,w0,th
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih),OPTIONAL :: tht,pext
!       IF(mype == 0) print *,'Nans check'
        CALL flush(6)
        CALL checknans (rho0,'rho0 ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        CALL checknans (u0  ,'u0   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        CALL checknans (v0  ,'v0   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        CALL checknans (w0  ,'w0   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        CALL checknans (th  ,'th   ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        IF(PRESENT(tht)) CALL checknans (tht ,'tht  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
        IF(PRESENT(pext)) CALL checknans (pext,'pext ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
!       IF(mype == 0) print *,'End of NaNs check '
        CALL flush(6)
 
      END SUBROUTINE check_nans_prognostic

      SUBROUTINE check_nans_coef(c11,c12,c13,c21,c22,c23,c31,c32,c33,etainv,np,mp,lp,ih)
        USE mpi_parallel, ONLY: mype
        USE mod_parameters, ONLY: ibcx,ibcy,ibcz,ipoles
        INTEGER(KIND=iintegers),INTENT(IN) ::  np,mp,lp,ih
        CHARACTER(5) :: varname
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: etainv 
        REAL(KIND=euwp),DIMENSION(np,mp,lp) :: c11,c12,c13,c21,c22,c23,c31,c32,c33
!       IF(mype == 0) print *,'Nans coef check'
        CALL flush(6)
        CALL checknans (c11    ,'c11  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,0)
        CALL checknans (c12    ,'c12  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,0)
        CALL checknans (c13    ,'c13  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,0)
        CALL checknans (c21    ,'c21  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,0)
        CALL checknans (c22    ,'c22  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,0)
        CALL checknans (c23    ,'c23  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,0)
        CALL checknans (c31    ,'c31  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,0)
        CALL checknans (c32    ,'c32  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,0)
        CALL checknans (c33    ,'c33  ',ipoles,ibcx,ibcy,ibcz,np,mp,lp,0)
        CALL checknans (etainv ,'etain',ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
!       IF(mype == 0) print *,'End of NaNs coef check '
        CALL flush(6)
 
      END SUBROUTINE check_nans_coef

      SUBROUTINE checknans(a,varname,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih)
      USE mpi_parallel, ONLY: mype
      USE,intrinsic :: ieee_arithmetic
      INTEGER(KIND=iintegers),INTENT(IN) ::  ipoles,np,mp,lp,ih
      REAL(KIND=euwp) :: a(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
      CHARACTER(5) :: varname
      INTEGER(KIND=iintegers),INTENT(IN):: ibcx,ibcy,ibcz
      INTEGER(KIND=iintegers) :: i,j,k
      
      DO k=1,lp
        DO j=1,mp
          DO i=1,np
            if(ieee_is_nan(a(i,j,k))) then
              print '(A,3I5,2A)','NaN detected at i,j,k:',i,j,k,' for var: ',varname
              STOP
            endif
          ENDDO
        ENDDO
      ENDDO
      END SUBROUTINE checknans

      SUBROUTINE courants_crash(ox,oy,oz,h,dxi,dyi,dzi,dt,np,mp,lp,ih)
      USE eulag_diagutils, ONLY: compute_courlipsch_Agrid_full
      USE scratch_datafields, ONLY: scr1,scr2 
        INTEGER(KIND=iintegers), INTENT(IN) :: np,mp,lp,ih
        REAL_euwp, INTENT(INOUT) ::  &
                          ox(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),    &
                          oy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2),    &
                          oz(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2)
        REAL_euwp, INTENT(IN) ::  &
                           h(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)
        REAL(KIND=euwp),INTENT(IN) :: dxi,dyi,dzi,dt

        CALL compute_courlipsch_Agrid_full(scr1,scr2,ox,oy,oz,h,dxi,dyi,dzi,dt,np,mp,lp,ih)
!       CALL exam_var(scr1,'Cournt',10,np,mp,lp,ih)
!       CALL exam_var(scr2,'Lipsch',10,np,mp,lp,ih)
#ifdef PNETCDF
        CALL pnet_out_chunk('Courant  ','diagn.nc',0,1,1,1,0,0,scr1,np,mp,lp,ih)
        CALL pnet_out_chunk('Lipsch   ','diagn.nc',0,1,1,1,0,0,scr2,np,mp,lp,ih)
        CALL pnet_out_chunk('ox       ','diagn.nc',0,1,1,1,0,0,ox(:,:,:,1),np,mp,lp,ih)
        CALL pnet_out_chunk('oy       ','diagn.nc',0,1,1,1,0,0,oy(:,:,:,1),np,mp,lp,ih)
        CALL pnet_out_chunk('oz       ','diagn.nc',0,1,1,1,0,0,oz(:,:,:,1),np,mp,lp,ih)
        CALL pnet_out_chunk('rho      ','diagn.nc',0,1,1,1,0,0,h,np,mp,lp,ih)
#endif
      END SUBROUTINE courants_crash

  SUBROUTINE print_time(ltimeadapt,lsubstepping,it,isub,nsub,dt,stime,timescale)
   USE mpi_parallel, ONLY: mype
   REAL(KIND=euwp), INTENT(IN) :: dt,stime,timescale
   INTEGER(KIND=iintegers), INTENT(IN) :: it,isub,nsub
   LOGICAL  :: ltimeadapt,lsubstepping
      if(mype.eq.0) then
        if(ltimeadapt) then
               if(timescale.eq.86400*30) then; print 1211,it,dt,stime
          else if(timescale.eq.86400)    then; print 1212,it,dt,stime
          else if(timescale.eq.3600)     then; print 1213,it,dt,stime
          else if(timescale.eq.60)       then; print 1214,it,dt,stime
          else if(timescale.eq.1)        then; print 1215,it,dt,stime
          else if(timescale.eq.1.e-3)    then; print 1216,it,dt,stime
          else if(timescale.eq.1.e-6)    then; print 1217,it,dt,stime
          endif
       else if (lsubstepping) then
        if(nsub.gt.1) then
               if(timescale.eq.86400*30) then;print 2211,it,dt,isub+1,nsub,stime
          else if(timescale.eq.86400)    then;print 2212,it,dt,isub+1,nsub,stime
          else if(timescale.eq.3600)     then;print 2213,it,dt,isub+1,nsub,stime
          else if(timescale.eq.60)       then;print 2214,it,dt,isub+1,nsub,stime
          else if(timescale.eq.1)        then;print 2215,it,dt,isub+1,nsub,stime
          else if(timescale.eq.1.e-3)    then;print 2216,it,dt,isub+1,nsub,stime
          else if(timescale.eq.1.e-6)    then;print 2217,it,dt,isub+1,nsub,stime
               endif
        else
               if(timescale.eq.86400*30) then;print 3211,it,dt,stime
          else if(timescale.eq.86400)    then;print 3212,it,dt,stime
          else if(timescale.eq.3600)     then;print 3213,it,dt,stime
          else if(timescale.eq.60)       then;print 3214,it,dt,stime
          else if(timescale.eq.1)        then;print 3215,it,dt,stime
          else if(timescale.eq.1.e-3)    then;print 3216,it,dt,stime
          else if(timescale.eq.1.e-6)    then;print 3217,it,dt,stime
               endif
        endif !nsub
       else !no substepping or adaptivity
            if(timescale.eq.86400*30)    then;print 211,it,dt,stime
          else if(timescale.eq.86400)    then;print 212,it,dt,stime
          else if(timescale.eq.3600)     then;print 213,it,dt,stime
          else if(timescale.eq.60)       then;print 214,it,dt,stime
          else if(timescale.eq.1)        then;print 215,it,dt,stime
          else if(timescale.eq.1.e-3)    then;print 216,it,dt,stime
          else if(timescale.eq.1.e-6)    then;print 217,it,dt,stime
            endif
        endif
      endif !mype
  211   format(/,8x,' it=',i10,' dt(sec)=',f7.2,' time(s.days)=',f7.2)
  212   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' time(days)=',f7.2)
  213   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' time(hours)=',f7.2)
  214   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' time(mins)=',f7.2)
  215   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' time(secs)=',f7.2)
  216   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' time(milisecs)=',f7.2)
  217   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' time(micrsecs)=',f7.2)
 1211   format(/,8x,' it=',i10,' adaptive dt(sec)=',f7.2,' time(s.days)=',f7.2)
 1212   format(/,8x,' it=',i5, ' adaptive dt(sec)=',f7.2,' time(days)=',f7.2)
 1213   format(/,8x,' it=',i5, ' adaptive dt(sec)=',f7.2,' time(hours)=',f7.2)
 1214   format(/,8x,' it=',i5, ' adaptive dt(sec)=',f7.2,' time(mins)=',f7.2)
 1215   format(/,8x,' it=',i5, ' adaptive dt(sec)=',f7.2,' time(secs)=',f7.2)
 1216   format(/,8x,' it=',i5, ' adaptive dt(sec)=',f7.2,' time(milisecs)=',f7.2)
 1217   format(/,8x,' it=',i5, ' adaptive dt(sec)=',f7.2,' time(micrsecs)=',f7.2)
 2211   format(/,8x,' it=',i10,' dt(sec)=',f7.2,' substep: ',i3,' of',i3,' time(s.days)=',f7.2)
 2212   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' substep: ',i3,' of',i3,' time(days)=',f7.2)
 2213   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' substep: ',i3,' of',i3,' time(hours)=',f7.2)
 2214   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' substep: ',i3,' of',i3,' time(mins)=',f7.2)
 2215   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' substep: ',i3,' of',i3,' time(secs)=',f7.2)
 2216   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' substep: ',i3,' of',i3,' time(milisecs)=',f7.2)
 2217   format(/,8x,' it=',i5, ' dt(sec)=',f7.2,' substep: ',i3,' of',i3,' time(micrsecs)=',f7.2)
 3211   format(/,8x,' it=',i10,' first substep dt(sec)=',f7.2,' time(s.days)=',f7.2)
 3212   format(/,8x,' it=',i5, ' first substep dt(sec)=',f7.2,' time(days)=',f7.2)
 3213   format(/,8x,' it=',i5, ' first substep dt(sec)=',f7.2,' time(hours)=',f7.2)
 3214   format(/,8x,' it=',i5, ' first substep dt(sec)=',f7.2,' time(mins)=',f7.2)
 3215   format(/,8x,' it=',i5, ' first substep dt(sec)=',f7.2,' time(secs)=',f7.2)
 3216   format(/,8x,' it=',i5, ' first substep dt(sec)=',f7.2,' time(milisecs)=',f7.2)
 3217   format(/,8x,' it=',i5, ' first substep dt(sec)=',f7.2,' time(micrsecs)=',f7.2)
  END SUBROUTINE print_time

  SUBROUTINE diagnos_logics(lprintoutput,lsavelong,lsaveshrt,  &
                            nprintoutput,nsavelong,nsaveshrt,  &
                            ltimeadapt,lsubstepping,it,isub,noutp,nstore,nplot,dt00,tt)

    USE mpi_parallel, ONLY: mype
    REAL(KIND=euwp), INTENT(IN) :: dt00,tt
    INTEGER(KIND=iintegers), INTENT(IN) :: it,isub ,noutp,nstore,nplot
    INTEGER(KIND=iintegers), INTENT(INOUT) :: nprintoutput,nsavelong,nsaveshrt
    LOGICAL, INTENT(OUT)  :: lprintoutput,lsavelong,lsaveshrt
    LOGICAL, INTENT(IN)  :: ltimeadapt,lsubstepping 
      if(it.eq.1) then
        nprintoutput=1
        nsavelong=1
        nsaveshrt=1
      endif

      lprintoutput=.FALSE.
      if(it.eq.1) lprintoutput=.TRUE.
      if((ltimeadapt.AND.(tt.ge.nprintoutput*noutp*dt00)) &
         .or.  &
         (.NOT.ltimeadapt.AND.(modulo(it,noutp).eq.0))) then
         lprintoutput=.TRUE.
         nprintoutput=nprintoutput+1
         IF(mype.eq.0) print *,'nprintoutput',nprintoutput,tt,noutp*dt00
      endif

      lsavelong=.FALSE.
      if(it.eq.1) lsavelong=.TRUE.
      if((ltimeadapt.AND.(tt.ge.nsavelong*nstore*dt00))  &
         .or. &
         (.NOT.ltimeadapt.AND.((modulo(it,nstore).eq.0)  &
                               .and.it.gt.1))) then
         lsavelong=.TRUE.
         nsavelong=nsavelong+1
         IF(mype.eq.0) print *,'nsavelong',nsavelong,tt,nstore*dt00
      endif

      lsaveshrt=.FALSE.
      if(it.eq.1) lsaveshrt=.TRUE.
      if((ltimeadapt.AND.(tt.ge.nsaveshrt*nplot*dt00))  &
         .or. &
         (.NOT.ltimeadapt.AND.((modulo(it,nplot).eq.0)  &
                               .and.it.gt.1))) then
         lsaveshrt=.TRUE.
         nsaveshrt=nsaveshrt+1
         IF(mype.eq.0) print *,'nsaveshrt',nsaveshrt,tt,nplot*dt00
      endif

  END SUBROUTINE diagnos_logics 

!----------------------------------------------------------------------!
SUBROUTINE diagnos_driver(ltimeadapt,lsubstepping, &
                          it,isub,noutp,nstore,nplot,dt00,tt, &
                          u,v,w,ox,oy,oz,rh,th,np,mp,lp,ih,tht,pstr,pext,  &
                          fx,fy,fz,ft,fxexpl,fyexpl,fzexpl,ftexpl,fpextexpl, &
                          qv,qv_old,qc,qr,qs,qi,qg)
!----------------------------------------------------------------------!

  USE mpi_parallel, ONLY: mype
  USE mpi_parallel, ONLY: globsumaxminv
#ifdef PNETCDF
  USE mpi_parallel, ONLY: pnet_out_lng_anel,pnet_out_lng_cmpr 
#endif
  USE eulag_diagutils, ONLY: wrtfmbuf,wrttobuf,calcavgmaxmin
  USE mod_parameters, ONLY: icmprss 


  INTEGER(KIND=iintegers),INTENT(IN) ::  np,mp,lp,ih
  INTEGER(KIND=iintegers) :: iframe
  INTEGER(KIND=iintegers), SAVE :: nprintoutput,nsavelong,nsaveshrt
  INTEGER(KIND=iintegers), INTENT(IN) :: it,isub,noutp,nstore,nplot
  LOGICAL  :: lprintoutput,lsavelong,lsaveshrt
  LOGICAL, INTENT(IN)  :: ltimeadapt,lsubstepping 
  REAL(KIND=euwp), INTENT(IN) :: dt00,tt
  REAL(KIND=euwp),INTENT(IN) :: &
     u(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2), &
     v(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2), &
     w(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2), &
    ox(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2), &
    oy(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2), &
    oz(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:2), &
    rh(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih,0:1)

  REAL(KIND=euwp),INTENT(INOUT) :: &
    th(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

  REAL(KIND=euwp),INTENT(IN),OPTIONAL :: &
    tht(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
   pstr(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih), &
   pext(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih)

   REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(IN),OPTIONAL ::  &
             fx,fy,fz,ft,fxexpl,fyexpl,fzexpl,ftexpl,fpextexpl

  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(IN),OPTIONAL ::  &
              qv,qv_old,qc
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(IN),OPTIONAL ::  &
            qr,qs,qi,qg
  INTEGER(KIND=iintegers), PARAMETER:: itapetypeshrt=2
  INTEGER(KIND=iintegers), PARAMETER:: itapetypelong=3

  CALL diagnos_logics(lprintoutput,lsavelong,lsaveshrt,  &
                      nprintoutput,nsavelong,nsaveshrt,  &
                      ltimeadapt,lsubstepping,it,isub,noutp,nstore,nplot, &
                      dt00,tt)
  IF(lprintoutput) CALL diagnos_actual(nprintoutput,u,v,w,th,np,mp,lp,ih)
#ifdef PNETCDF
  IF(icmprss.eq.0) THEN
    IF(lsavelong) CALL pnet_out_lng_anel(itapetypelong,np,mp,lp,ih,nsavelong,u,v,w,ox,oy,oz,   &
                                         rh,th,    &
                                         fx,fy,fz,ft,fxexpl,fyexpl,fzexpl,ftexpl, &
                                         qv,qv_old,qc,qr,qs,qi,qg)
    IF(lsaveshrt) CALL pnet_out_lng_anel(itapetypeshrt,np,mp,lp,ih,nsaveshrt,u,v,w,ox,oy,oz,rh,th)
  ELSE
    IF(lsavelong) CALL pnet_out_lng_cmpr(itapetypelong,np,mp,lp,ih,nsavelong,u,v,w,ox,oy,oz,   &
                                         rh,th,tht,pstr,pext,  &
                                         fx,fy,fz,ft,fxexpl,fyexpl,fzexpl,ftexpl,fpextexpl, &
                                         qv,qv_old,qc,qr,qs,qi,qg)
    IF(lsaveshrt) CALL pnet_out_lng_cmpr(itapetypeshrt,np,mp,lp,ih,nsaveshrt,u,v,w,ox,oy,oz,rh,th,tht,pstr,pext)
  ENDIF
#endif
END SUBROUTINE diagnos_driver

END MODULE eulag_diagnostics


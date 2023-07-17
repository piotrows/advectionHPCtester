MODULE mpi_parallel
#ifdef CUDACODE
  USE cudafor
#endif
#ifdef TIMERSOMP 
   USE omp_lib
#endif
#ifdef PUREMPI
#ifdef MPI90
   USE mpi
#endif /*PUREMPI*/
#endif /*PUREMPI*/
   USE precisions, ONLY  : iintegers,euwp,sp,dp
   USE parameters, ONLY  : nprocx,nprocy,nprocz,msz,mucnt
   USE parameters, ONLY  : ibcx,ibcy,ibcz,isphere,ipoles
   USE parameters, ONLY  : n
   IMPLICIT NONE
#ifdef PUREMPI 
#ifdef MPI77
  INCLUDE "mpif.h"
#endif
#endif
#include "../advection/src_algorithms/defines.inc"

   INTEGER (KIND=iintegers) botedge,topedge,leftedge,rightedge,      &
      gndedge,skyedge,middle
   INTEGER  mype, mysize
   INTEGER  myce
   INTEGER  :: npe=1
   INTEGER  myple, mypme, mypne
   INTEGER  my_prcnm, my_prcnl, my_prcml
   INTEGER  myrow, mycol
!   INTEGER (KIND=iintegers) npos,mpos,lpos,                          &
   INTEGER  npos,mpos,lpos,nsubpos,msubpos,lsubpos,                   &
             peW,peE,peN,peS,peNE,peSE,peNW,peSW,peG,                 &
             peGW,peGE,peGN,peGS,peGNE,peGSE,peGNW,peGSW,peZ,         &
             peZW,peZE,peZN,peZS,peZNE,peZSE,peZNW,peZSW
   INTEGER mpi_cust_comm,DC_TYPE,DC_TYPE_DOUBLE,DC_TYPE_SINGLE
   INTEGER(KIND=IINTEGERS), DIMENSION(:),ALLOCATABLE :: ip !(0:n)
   INTEGER,PROTECTED :: iup,iupx,iupy,iupz
   INTEGER,PARAMETER :: isend=3
   INTEGER,PARAMETER :: updmulgran=1 

   INTEGER, PARAMETER :: maxtimers=2000
   INTEGER icounts(maxtimers)
#ifdef CUDACODE
   integer(kind=cuda_stream_kind) :: istream1,istream2,istream3,istream4
   type(cudaEvent) :: cttstart(maxtimers),cttend(maxtimers) 
!  DIMENSION(maxtimers) :: cttstart,cttend 
#endif
   REAL(KIND=dp) rtimerb(maxtimers)
   REAL(KIND=dp) rtimer0(maxtimers)
   REAL(KIND=dp) rtimerbav(maxtimers)
   REAL(KIND=dp) rtimerbmx(maxtimers)
   REAL(KIND=dp) rtimerbmn(maxtimers)
   CHARACTER clabels(maxtimers)*25
   CHARACTER(LEN=2) clabclr(maxtimers)
  INTEGER :: MPI_MAXMIN,MPI_SUMAX,MPI_SUMAXMIN, MPI_SUMAXMINLOC
!

CONTAINS
   SUBROUTINE end_code
   INTEGER i,ierr
#ifdef CUDACODE
     DO i=1,maxtimers
        ierr = cudaEventDestroy(cttstart(i))
        ierr = cudaEventDestroy(cttend(i))
     ENDDO
!       ierr= cudaDeviceReset()
#endif

     IF(ALLOCATED(ip)) DEALLOCATE (ip)
   END SUBROUTINE end_code

   SUBROUTINE end_mpi
   INTEGER ierr
#ifdef PUREMPI
      CALL MPI_Comm_free(my_prcnm,ierr)
      CALL MPI_Comm_free(my_prcnl,ierr)
      CALL MPI_Comm_free(my_prcml,ierr)
      CALL MPI_Op_free(MPI_MAXMIN,ierr)
      CALL MPI_Op_free(MPI_SUMAXMIN,ierr)
      CALL MPI_Op_free(MPI_SUMAXMINLOC,ierr)
      CALL MPI_Finalize(ierr)
#endif
   END SUBROUTINE end_mpi

   SUBROUTINE geomset(lmpiinit)
   USE parameters, ONLY: isphere,ipoles
      LOGICAL,OPTIONAL ::  lmpiinit
!
!     setup geometry information for each processor
!

      INTEGER iprr,ipri,i
#ifdef CUDACODE
  integer :: istat, nDevices, minpriority,maxpriority
  type (cudaDeviceProp) :: prop
  character(len=10) cudastr
  logical cudaflag
#endif
#ifdef PUREMPI
      INTEGER old_comm, new_comm, ndims!, reorder
      INTEGER dim_size(3)
      LOGICAL periods(3),reorder
      INTEGER coords(3),coordl(3)
      INTEGER ierr,ierr1,ierr1a,ierr2,ierr3,ierr4,rank
      ierr=0
      IF(lmpiinit) THEN
        CALL MPI_Init(ierr)
!************ Prepare copy of default communicator ************************
!       call MPI_Comm_dup(MPI_COMM_WORLD, mpi_cust_comm,ierr)
        mpi_cust_comm=MPI_COMM_WORLD
      ENDIF

#ifdef MPI90
!        CALL MPI_Type_create_f90_real( 15 , 307 , DC_TYPE_DOUBLE , ierr )
!        CALL MPI_Type_create_f90_real( 6  ,  37 , DC_TYPE_SINGLE , ierr )
         DC_TYPE_DOUBLE = MPI_REAL8 
         DC_TYPE_SINGLE = MPI_REAL4 
#endif
#ifdef MPI77
         DC_TYPE_DOUBLE = MPI_DOUBLE_PRECISION
         DC_TYPE_SINGLE = MPI_REAL
#endif
      IF(euwp == dp) THEN
         DC_TYPE = DC_TYPE_DOUBLE  
      ELSE
         DC_TYPE = DC_TYPE_SINGLE 
      ENDIF


!****** Check if no. of declared cores = no. of allocated cores ***********
       CALL MPI_Comm_size(mpi_cust_comm, mysize, ierr)  ! total numbers of PE''s
       CALL MPI_Comm_size(MPI_COMM_WORLD, mysize, ierr)  ! total numbers of PE''s
       npe=nprocx*nprocy*nprocz
      IF(mysize/=npe) THEN
      CALL MPI_Comm_rank(mpi_cust_comm, rank, ierr)  ! number of current PE
!     print *,'mype,mysize',rank,mysize
!     CALL flush(6)
      mype = rank
         IF(mype==0) THEN
            PRINT *,'!!!!!!  Prescribed processor number is not what is available from MPI !!!!!!'
            PRINT *,'!!!!!!  MPI_Comm_size .ne. nprocx*nprocy*nprocz  !!!!!!'
            PRINT *,'MPI size= ',mysize,'nprocx*nprocy*nprocz=',npe
            PRINT *,'nprocx ',nprocx,'nprocy',nprocy,'nprocz',nprocz
            PRINT *,'!!!!!!         EXIT             !!!!!!'
         ENDIF
         CALL MPI_Finalize(ierr)
         STOP 'GEOMSET'
      ENDIF

!****** If using MPI Cartesian topology, redefine communicator **************
      old_comm=mpi_cust_comm
      ndims = 3
      dim_size(1) = nprocx ! rows
      dim_size(2) = nprocy ! columns
      dim_size(3) = nprocz ! vertical
      periods(1) = .FALSE.
      periods(2) = .FALSE.
      periods(3) = .FALSE.
!        if(ibcx.eq.1) periods(1)=.TRUE. ! row periodicity
!        if(ibcy.eq.1) periods(2)=.TRUE. ! column periodicity
!        if(ibcz.eq.1) periods(3)=.TRUE. ! vertical periodicity
      reorder = .TRUE.  ! Allow for reordering cores, default 1
      reorder = .FALSE.
      IF(lmpiinit) THEN 
        CALL MPI_Cart_create(old_comm, ndims, dim_size,                 &
                periods, reorder, new_comm, ierr)
        mpi_cust_comm=new_comm
      ENDIF


!******* Define planes  ***************************************
        call MPI_Comm_split(mpi_cust_comm,myple,lpos-1,my_prcnm,ierr)
        call MPI_Comm_split(mpi_cust_comm,mypme,mpos-1,my_prcnl,ierr)
        call MPI_Comm_split(mpi_cust_comm,mypne,npos-1,my_prcml,ierr)
       myrow=my_prcml
       mycol=my_prcnl
!******* Ask for this core rank number ***************************************
      CALL MPI_Comm_rank(mpi_cust_comm, rank, ierr)  ! number of current PE
      mype = rank
#ifdef CUDACODE
  istat = cudaGetDeviceCount(nDevices)
! if(mype.eq.0) then
    do i = 0, nDevices-1
       istat = cudaGetDeviceProperties(prop, i)
       write(*,"(' Device Number: ',i0)") i
       write(*,"('   Device name: ',a)") trim(prop%name)
       write(*,"('   Memory Clock Rate (KHz): ', i0)") &
         prop%memoryClockRate
       write(*,"('   Memory Bus Width (bits): ', i0)") &
         prop%memoryBusWidth
       write(*,"('   Peak Memory Bandwidth (GB/s): ', f12.4)") &
         2.0*prop%memoryClockRate*(prop%memoryBusWidth/8)/10.0**6
       write(*,*)        
    enddo
! endif
!     call flush()
        ierr1a=cudaDeviceGetStreamPriorityRange(minpriority,maxpriority) 
        IF(mype.eq.0) print *,'Priority range for streams',minpriority,maxpriority
      DO iprr=0,npe
        IF(mype.eq.iprr) THEN
        myce=mod(mype,4)
        myce=0 
        ierr1=cudaSetDevice(myce)
        istream1=-10
        istream2=-10
        istream3=-10
        ierr2=cudaStreamCreateWithPriority(istream1,cudaStreamNonBlocking,maxpriority+2)
        ierr3=cudaStreamCreateWithPriority(istream2,cudaStreamNonBlocking,maxpriority+2)
        ierr4=cudaStreamCreateWithPriority(istream3,cudaStreamNonBlocking,maxpriority)
!       ierr2=cudaStreamCreate(istream1)!,cudaStreamNonBlocking,maxpriority+2)
!       ierr3=cudaStreamCreate(istream2)!,cudaStreamNonBlocking,maxpriority+2)
!       ierr4=cudaStreamCreate(istream3)!,cudaStreamNonBlocking,maxpriority)
!       print *,'CUDA set device from mype ', mype,'no of cuda device',myce,"streams:",istream1,istream2,istream3, &
!                                                                            'stat: ',cudaGetErrorString(ierr1), &
!                                                                                    cudaGetErrorString(ierr2), & 
!                                                                                    cudaGetErrorString(ierr3), &
!                                                                                    cudaGetErrorString(ierr4)
        call flush(6)
        cudaflag=.FALSE.
        CALL MPI_Info_get(MPI_INFO_ENV, "cuda_aware",10,cudastr , cudaflag,ierr);
        if (cudaflag) then
                print *,"The CUDA awareness in MPI is activated."
        else
                print *,"Parastation MPI seems no CUDA aware or not a Parastation MPI."
        endif
      
        call flush(6)
        ENDIF
        CALL MPI_BARRIER(mpi_cust_comm,ierr)
      ENDDO
#endif


      middle = 0
      rightedge = 0
      leftedge  = 0
      botedge   = 0
      topedge   = 0
      skyedge   = 0
      gndedge   = 0
      CALL MPI_Cart_coords(mpi_cust_comm, mype,ndims,coords, ierr)
      npos=coords(1)+1
      mpos=coords(2)+1
      lpos=coords(3)+1
      coordl(1)=coords(1)
      coordl(2)=coords(2)
      coordl(3)=coords(3)

!Define peS,peGS,peZS
      coordl(1)=coords(1)
      coordl(2)=coords(2)-1
      coordl(3)=coords(3)
      IF(coordl(2)==-1) coordl(2)=coords(2)+nprocy-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peS,ierr)

      coordl(3)=coords(3)-1
      IF(coordl(3)==-1) coordl(3)=nprocz-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGS,ierr)

      coordl(3)=coords(3)+1
      IF(coordl(3)==nprocz) coordl(3)=0
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZS,ierr)

!Define peN,peGN,peZN
      coordl(1)=coords(1)
      coordl(2)=coords(2)+1
      coordl(3)=coords(3)
      IF(coordl(2)==nprocy) coordl(2)=coords(2)-nprocy+1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peN,ierr)

      coordl(3)=coords(3)-1
      IF(coordl(3)==-1) coordl(3)=nprocz-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGN,ierr)

      coordl(3)=coords(3)+1
      IF(coordl(3)==nprocz) coordl(3)=0
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZN,ierr)

!Define peW,peGW,peZW
      coordl(1)=coords(1)-1
      coordl(2)=coords(2)
      coordl(3)=coords(3)
      IF(coordl(1)==-1) coordl(1)=coords(1)+nprocx-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peW,ierr)

      coordl(3)=coords(3)-1
      IF(coordl(3)==-1) coordl(3)=nprocz-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGW,ierr)

      coordl(3)=coords(3)+1
      IF(coordl(3)==nprocz) coordl(3)=0
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZW,ierr)

!Define peE,peGE,peZE
      coordl(1)=coords(1)+1
      coordl(2)=coords(2)
      coordl(3)=coords(3)
      IF(coordl(1)==nprocx) coordl(1)=coords(1)-nprocx+1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peE,ierr)

      coordl(3)=coords(3)-1
      IF(coordl(3)==-1) coordl(3)=nprocz-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGE,ierr)

      coordl(3)=coords(3)+1
      IF(coordl(3)==nprocz) coordl(3)=0
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZE,ierr)

      coordl(1)=coords(1)
      coordl(2)=coords(2)
      coordl(3)=coords(3)-1
      IF(coordl(3)==-1) coordl(3)=coords(3)+nprocz-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peG,ierr)

      coordl(3)=coords(3)+1
      IF(coordl(3)==nprocz) coordl(3)=coords(3)-nprocz+1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZ,ierr)

!Define peNE,GNE,ZNE
      coordl(1)=coords(1)+1
      coordl(2)=coords(2)+1
      coordl(3)=coords(3)
      IF(coordl(1)==nprocx) coordl(1)=0
      IF(coordl(2)==nprocy) coordl(2)=0
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peNE,ierr)
      coordl(3)=coords(3)-1
      IF(coordl(3)==-1) coordl(3)=nprocz-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGNE,ierr)
      coordl(3)=coords(3)+1
      IF(coordl(3)==nprocz) coordl(3)=0
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZNE,ierr)
!Define peSE,GSE,ZSE
      coordl(1)=coords(1)+1
      coordl(2)=coords(2)-1
      coordl(3)=coords(3)
      IF(coordl(1)==nprocx) coordl(1)=0
      IF(coordl(2)==-1) coordl(2)=nprocy-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peSE,ierr)
      coordl(3)=coords(3)-1
      IF(coordl(3)==-1) coordl(3)=nprocz-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGSE,ierr)
      coordl(3)=coords(3)+1
      IF(coordl(3)==nprocz) coordl(3)=0
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZSE,ierr)

!Define peSW,GSW,ZSW
      coordl(1)=coords(1)-1
      coordl(2)=coords(2)-1
      coordl(3)=coords(3)
      IF(coordl(1)==-1) coordl(1)=nprocx-1
      IF(coordl(2)==-1) coordl(2)=nprocy-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peSW,ierr)
      coordl(3)=coords(3)-1
      IF(coordl(3)==-1) coordl(3)=nprocz-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGSW,ierr)
      coordl(3)=coords(3)+1
      IF(coordl(3)==nprocz) coordl(3)=0
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZSW,ierr)

!Define peNW,GNW,ZNW
      coordl(1)=coords(1)-1
      coordl(2)=coords(2)+1
      coordl(3)=coords(3)
      IF(coordl(1)==-1) coordl(1)=nprocx-1
      IF(coordl(2)==nprocy) coordl(2)=0
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peNW,ierr)
      coordl(3)=coords(3)-1
      IF(coordl(3)==-1) coordl(3)=nprocz-1
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGNW,ierr)
      coordl(3)=coords(3)+1
      IF(coordl(3)==nprocz) coordl(3)=0
      CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZNW,ierr)
!Define edge marks
      IF (coords(1)==nprocx-1) rightedge = 1
      IF (coords(1)==0)    leftedge  = 1
      IF (coords(2)==0)    botedge   = 1
      IF (coords(2)==nprocy-1) topedge   = 1
      IF (coords(3)==0)    gndedge   = 1
      IF (coords(3)==nprocz-1) skyedge   = 1
      IF (rightedge==0 .AND. leftedge==0 .AND. botedge==0.AND.   &
            topedge==0 .AND.  gndedge==0 .AND. skyedge==0) middle = 1
      IF(isphere==1.AND.ipoles==1) CALL geomset_sphere(coords,coordl)

      ipri=0
      IF (ipri==1) THEN
         IF (mype==0) THEN
            PRINT *, 'Initializing processor setup'
            PRINT 96,nprocx,nprocy,nprocz
            PRINT 97,mype,middle,rightedge,leftedge,botedge,topedge,npos,mpos,&
                      peW,peE,peN,peS,peNE,peSE,peNW,peSW
         ENDIF

         CALL MPI_BARRIER(mpi_cust_comm,ierr)

         IF (mype/=0) THEN
            PRINT 97,mype,middle,rightedge,leftedge,botedge,topedge,npos,mpos,&
                          peW,peE,peN,peS,peNE,peSE,peNW,peSW
         ENDIF

         CALL MPI_BARRIER(mpi_cust_comm,ierr)

96       FORMAT('nprocx = ',i3,'  nprocy = ',i3,' nprocz = ',i3/,                         &
          ' my  mid  R   L   B   T   N   M  pW  pE  pN  pS pNE pSE pNW pSW')
97       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',   &
                 i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3)
         Call flush(6)
      ENDIF
        
    CALL MPI_OP_CREATE(maxmin,.TRUE.,MPI_MAXMIN,ierr)
    CALL MPI_OP_CREATE(sumaxmin,.TRUE.,MPI_SUMAXMIN,ierr)
    CALL MPI_OP_CREATE(sumaxminloc,.TRUE.,MPI_SUMAXMINLOC,ierr)
#else  /*PUREMPI*/
      IF(nprocx*nprocy*nprocz.ne.1) &
      STOP 'Wrong number of processors nprocx,nprocy,nprocz for single processor run'
      mype = 0
      middle = 0
      rightedge = 1
      leftedge = 1
      botedge = 1
      topedge = 1
      skyedge = 1
      gndedge = 1
      npos = 1
      mpos = 1
      lpos = 1
      peW = 0
      peE = 0
      peN = 0
      peS = 0
      peGW = 0
      peGE = 0
      peGN = 0
      peGS = 0
      peZW = 0
      peZE = 0
      peZN = 0
      peZS = 0
      peNE = 0
      peSE = 0
      peNW = 0
      peSW = 0
      peGNE = 0
      peGSE = 0
      peGNW = 0
      peGSW = 0
      peGNE = 0
      peGSE = 0
      peGNW = 0
      peGSW = 0
      npos=1
      mpos=1
      lpos=1
#endif /*PUREMPI*/
!
    iup=1+max(ibcx,ibcy,ibcz)
    iupx=1+ibcx
    iupy=1+ibcy
    iupz=1+ibcz
!   IF(mype == 0) print *,'EULAG update limits are: ',iup,iupx,iupy,iupz
!Vector to understand who is our neighbour at the other side of north/south pole
      IF(isphere*ipoles.eq.1) THEN
      ALLOCATE (ip(0:n))
      DO i=0,n
         ip(i)=mod(i+    nint(n/2.)-1,n  )+1
 !        print *,'i,ip(i)',i,ip(i)
      ENDDO
      ENDIF

!Diagnostics:
      ipri=0
      IF (ipri==1) THEN
         IF (mype==0) THEN
            PRINT 98,nprocx,nprocy,nprocz
         ENDIF
         npe=nprocx*nprocy*nprocz
         DO iprr=0,npe
            IF (mype==iprr) THEN
               PRINT 99,mype,rightedge,leftedge,                                 &
                         botedge,topedge,                                         &
                         gndedge,skyedge,npos,mpos,lpos,                          &
                         peW,peE,peN,peS,peNE,peSE,peNW,peSW,peG,                 &
                         peGW,peGE,peGN,peGS,peGNE,peGSE,peGNW,peGSW,peZ,         &
                         peZW,peZE,peZN,peZS,peZNE,peZSE,peZNW,peZSW
         Call flush(6)
            ENDIF
         ENDDO
98       FORMAT('nprocx = ',i3,'  nprocy = ',i3,'  nprocz = ',i3/,         &
          ' mype R   L   B   T   G   S   N   M   L  pW  pE pN pS'//        &
          ' pNE pSE pNW pSW  pG pGW pGE pGN pGS pGNE pGSE pGNW pGSW pZ pZW pZ'//&
          'E pZN pZS pZNE pZSE pZNW pZSW')
99       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',                               &
          i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,'',i3,'',i3,'',      &
          i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,'  ',  &
          i3,'  ',i3,                                                       &
          '  ',i3,'  ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,'  ',   &
          i3,'  ',i3,'  ',i3)
      ENDIF
   END SUBROUTINE geomset
#ifdef PUREMPI
   SUBROUTINE geomset_sphere(coords,coordl)
      INTEGER coords(3),coordl(3)
      INTEGER ierr,itemp
!Redefine processors at poles
      IF(botedge==1) THEN
!REDefine peS,peGS,peZS
         coordl(1)=coords(1)
         coordl(2)=coords(2)
         coordl(3)=coords(3)
         IF(coordl(1)*2<nprocx) THEN
            coordl(1)=coords(1)+floor(.5*nprocx)
         ELSE
            coordl(1)=coords(1)-floor(.5*nprocx)
         ENDIF
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peS,ierr)
         itemp=coordl(1)
         coordl(1)=itemp-1
         IF (coordl(1)==-1) coordl(1)=nprocx-1
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peSE,ierr)
         coordl(3)=coords(3)-1
         IF(coordl(3)==-1) coordl(3)=nprocz-1
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGSE,ierr)
         coordl(3)=coords(3)+1
         IF(coordl(3)==nprocz) coordl(3)=0
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZSE,ierr)

         coordl(1)=itemp+1
         coordl(3)=coords(3)
         IF (coordl(1)==nprocx) coordl(1)=0
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peSW,ierr)
         coordl(3)=coords(3)-1
         IF(coordl(3)==-1) coordl(3)=nprocz-1
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGSW,ierr)
         coordl(3)=coords(3)+1
         IF(coordl(3)==nprocz) coordl(3)=0
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZSW,ierr)

         coordl(1)=itemp
         coordl(3)=coords(3)-1
         IF(coordl(3)==-1) coordl(3)=nprocz-1
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGS,ierr)

         coordl(3)=coords(3)+1
         IF(coordl(3)==nprocz) coordl(3)=0
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZS,ierr)


      ENDIF !botedge

      IF(topedge==1) THEN
!Define peN,peGN,peZN
         coordl(1)=coords(1)
         coordl(2)=coords(2)
         coordl(3)=coords(3)
         IF(coordl(1)*2<nprocx) THEN
            coordl(1)=coords(1)+floor(.5*nprocx)
         ELSE
            coordl(1)=coords(1)-floor(.5*nprocx)
         ENDIF
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peN,ierr)
         itemp=coordl(1)
         coordl(1)=itemp-1
         IF (coordl(1)==-1) coordl(1)=nprocx-1
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peNE,ierr)
         coordl(3)=coords(3)-1
         IF(coordl(3)==-1) coordl(3)=nprocz-1
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGNE,ierr)
         coordl(3)=coords(3)+1
         IF(coordl(3)==nprocz) coordl(3)=0
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZNE,ierr)

         coordl(1)=itemp+1
         coordl(3)=coords(3)
         IF (coordl(1)==nprocx) coordl(1)=0
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peNW,ierr)
         coordl(3)=coords(3)-1
         IF(coordl(3)==-1) coordl(3)=nprocz-1
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGNW,ierr)
         coordl(3)=coords(3)+1
         IF(coordl(3)==nprocz) coordl(3)=0
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZNW,ierr)

         coordl(1)=itemp
         coordl(3)=coords(3)-1
         IF(coordl(3)==-1) coordl(3)=nprocz-1
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peGN,ierr)

         coordl(3)=coords(3)+1
         IF(coordl(3)==nprocz) coordl(3)=0
         CALL MPI_CART_RANK(mpi_cust_comm,coordl,peZN,ierr)
      ENDIF !topedge
   END SUBROUTINE geomset_sphere

   SUBROUTINE mybarrier
      INTEGER ierr
      CALL MPI_BARRIER(mpi_cust_comm,ierr)
   END SUBROUTINE mybarrier
#else
   SUBROUTINE mybarrier
   END SUBROUTINE mybarrier
#endif /*PUREMPI*/

#if (STATICMEM == 0)
   SUBROUTINE init_subdomains
   USE parameters, ONLY: np,mp,lp,n,m,l,nprocx,nprocy,nprocz 
   INTEGER(KIND=iintegers) npardiff,mpardiff,lpardiff
   INTEGER ierr,ipri
!    np=CEILING(REAL(n)/REAL(nprocx))
!    nsubpos=(npos-1)*np
!    IF(rightedge.eq.1) np=min(np,n-np*(npos-1))
!    mp=CEILING(REAL(m)/REAL(nprocy))
!    msubpos=(mpos-1)*mp
!    IF(topedge.eq.1) mp=min(mp,m-mp*(mpos-1))
!    lp=CEILING(REAL(l)/REAL(nprocz))
!    lsubpos=(lpos-1)*lp
!    IF(skyedge.eq.1) lp=min(lp,l-lp*(lpos-1))
   np=CEILING(REAL(n)/REAL(nprocx))
   mp=CEILING(REAL(m)/REAL(nprocy))
   lp=CEILING(REAL(l)/REAL(nprocz))
   npardiff=np*nprocx-n
    nsubpos=np*(npos-1)
   IF((npardiff-1).gt.(nprocx-npos)) nsubpos=np*(npos-1)-((npardiff-1)-(nprocx-npos))
   IF((nprocx-npos).lt.npardiff) np=np-1

    mpardiff=mp*nprocy-m
     msubpos=mp*(mpos-1)
   IF((mpardiff-1).gt.(nprocy-mpos)) msubpos=mp*(mpos-1)-((mpardiff-1)-(nprocy-mpos))
   IF((nprocy-mpos).lt.mpardiff) mp=mp-1

   lpardiff=lp*nprocz-l
    lsubpos=lp*(lpos-1)
   IF((lpardiff-1).gt.(nprocz-lpos)) lsubpos=lp*(lpos-1)-((lpardiff-1)-(nprocz-lpos))
   IF((nprocz-lpos).lt.lpardiff) lp=lp-1

   IF (np.le.0) THEN
     print *,'Original np',CEILING(REAL(n)/REAL(nprocx))
     print *,'New np',np,'nprocx-1',nprocx-1,'n',n
     STOP 'negative np'
   ENDIF
   IF (mp.le.0) THEN
     print *,'Original mp',CEILING(REAL(m)/REAL(nprocy))
     print *,'New mp',mp,'nprocy-1',nprocy-1,'m',m
     STOP 'negative mp'
   ENDIF
   IF (lp.le.0) THEN
     print *,'Original lp',CEILING(REAL(l)/REAL(nprocz))
     print *,'New lp',lp,'nprocz-1',nprocz-1,'l',l
     STOP 'negative lp'
   ENDIF
      ipri=0
      IF (ipri==1) THEN
         IF (mype==0) THEN
            PRINT *, 'Initializing subdomains'
            PRINT 96,nprocx,nprocy,nprocz,np,mp,lp
            PRINT 97,mype,npos,mpos,lpos,nsubpos,msubpos,lsubpos,&
                      peW,peE,peN,peS,peZ,peG
            CALL flush(9)
         ENDIF

#ifdef PUREMPI
        CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif

         IF (mype/=0) THEN
            PRINT 96,nprocx,nprocy,nprocz,np,mp,lp
            PRINT 97,mype,npos,mpos,lpos,nsubpos,msubpos,lsubpos,&
                      peW,peE,peN,peS,peZ,peG
            CALL flush(9)
         ENDIF

#ifdef PUREMPI
        CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif

96       FORMAT('nprocx = ',i3,'  nprocy = ',i3,' nprocz = ',i3,' np=',i3,' mp=',i3,' lp=',i3/,   &
          ' my   N   M   L  NS  MS  LS   pW  pE  pN  pS pZ  pG')
97       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',i3,' ',   &
                 i3,' ',i3,' ',i3,' ',i3,' ',i3)
      ENDIF
   END SUBROUTINE init_subdomains
#else
   SUBROUTINE init_subdomains
   USE parameters, ONLY: np,mp,lp,n,m,l 
     nsubpos=(npos-1)*np
     msubpos=(mpos-1)*mp
     lsubpos=(lpos-1)*lp
   END SUBROUTINE init_subdomains

#endif /*STATICMEM*/
   SUBROUTINE geomset_lib(mpi_comm_in,npos_in,mpos_in,lpos_in,nsubpos_in,msubpos_in,lsubpos_in, &
                          nprocx_in,nprocy_in,nprocz_in,                                        &
                          leftedge_in,rightedge_in,botedge_in,topedge_in,gndedge_in,skyedge_in, &
                           peE_in, peW_in, peS_in, peN_in, peG_in, peZ_in)

   USE parameters, ONLY  : ibcx,ibcy,ibcz
   INTEGER (KIND=iintegers),INTENT(IN) :: mpi_comm_in,                &
                                          leftedge_in,rightedge_in,   &
                                           botedge_in,  topedge_in,   &
                                           gndedge_in,  skyedge_in,   &
                                              npos_in,     mpos_in, lpos_in, &
                                              nprocx_in, nprocy_in, nprocz_in, &
                                           nsubpos_in,  msubpos_in, lsubpos_in
 

   INTEGER (KIND=iintegers),INTENT(IN) :: peE_in, peW_in, peS_in, peN_in, peG_in, peZ_in
#ifdef PUREMPI
   INTEGER ierr
#endif
   INTEGER rank,iprint,iprr,mysize

      mpi_cust_comm=mpi_comm_in 
!******* Ask for this core rank number ***************************************
#ifdef PUREMPI
      CALL MPI_Comm_rank(mpi_cust_comm, rank, ierr)  ! number of current PE
#else
      rank=0
#endif
      mype = rank
      middle = 0
      rightedge = rightedge_in 
      leftedge  =  leftedge_in 
      botedge   =   botedge_in
      topedge   =   topedge_in
      skyedge   =   skyedge_in
      gndedge   =   gndedge_in
      peE       = peE_in
      peW       = peW_in
      peS       = peS_in
      peN       = peN_in
      peG       = peG_in
      peZ       = peZ_in
      npos      = npos_in
      mpos      = mpos_in
      lpos      = lpos_in

#ifdef PUREMPI
      IF(euwp == dp) THEN
        DC_TYPE = MPI_DOUBLE_PRECISION
      ELSE IF (euwp == sp) THEN
        DC_TYPE = MPI_REAL
      ELSE 
        STOP  'Wrong precision in geomset_lib'
      ENDIF
#endif /*PUREMPI*/
    iup=1+max(ibcx,ibcy,ibcz)
    iupx=1+ibcx
    iupy=1+ibcy
    iupz=1+ibcz
!   IF(mype == 0) print *,'EULAG update limits are: ',iup,iupx,iupy,iupz
!     IF(rightedge.eq.0) THEN
!        nsubpos=(npos-1)*np
!     ELSE
!        nsubpos=n-np
!     ENDIF
!     IF(topedge.eq.0) THEN
!        msubpos=(mpos-1)*mp
!    ELSE
!        msubpos=m-mp
!     ENDIF
!     IF(skyedge.eq.0) THEN
!        lsubpos=(lpos-1)*lp
!     ELSE
!        lsubpos=l-lp
!     ENDIF
     nsubpos=nsubpos_in
     msubpos=msubpos_in
     lsubpos=lsubpos_in
#if (STATICMEM == 0)
     nprocx=nprocx_in
     nprocy=nprocy_in
     nprocz=nprocz_in
#endif
!******* Define planes  ***************************************
      myple = (mpos-1)*nprocx+npos-1
!Define extra planes for mpos and npos, only for communicator splitting purposes
      mypme = (lpos-1)*nprocx+npos-1
      mypne = (lpos-1)*nprocy+mpos-1

#ifdef PUREMPI
        call MPI_Comm_split(mpi_cust_comm,myple,lpos-1,my_prcnm,ierr)
        call MPI_Comm_split(mpi_cust_comm,mypme,mpos-1,my_prcnl,ierr)
        call MPI_Comm_split(mpi_cust_comm,mypne,npos-1,my_prcml,ierr)
#endif /*PUREMPI*/
       myrow=my_prcml
       mycol=my_prcnl
!******* Ask for this core rank number ***************************************
#ifdef PUREMPI
    CALL MPI_OP_CREATE(sumaxmin,.TRUE.,MPI_MAXMIN,ierr)
    CALL MPI_OP_CREATE(sumaxmin,.TRUE.,MPI_SUMAXMIN,ierr)
    CALL MPI_OP_CREATE(sumaxminloc,.TRUE.,MPI_SUMAXMINLOC,ierr)
#endif /*PUREMPI*/

      npe=nprocx*nprocy*nprocz
!Diagnostics:
      iprint=0
!     IF(mype == 0)     PRINT 98,nprocx,nprocy,nprocz
      IF (iprint==1) THEN
      mysize=1
#ifdef PUREMPI
       CALL MPI_Comm_size(mpi_cust_comm, mysize, ierr)  ! total numbers of PE''s
#endif /*PUREMPI*/
      IF(mype == 0)     PRINT *,'Mysize:',mysize
         DO iprr=0,mysize
            IF (mype==iprr) THEN
               PRINT 99,mype,rightedge,leftedge,                                 &
                         botedge,topedge,                                         &
                         npos,mpos,lpos,                          &
                         peW,peE,peN,peS
               CALL flush(6)
            ENDIF
#ifdef PUREMPI
         CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif /*PUREMPI*/
         ENDDO
98       FORMAT('nprocx = ',i3,'  nprocy = ',i3,'  nprocz = ',i3/,         &
          ' mype R   L   B   T   N   M   L    pW    pE    pN    pS')
99       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',                               &
                i3,' ',i3,' ',i3,' ',i3,' ', &
                i5,' ',i5,' ',i5,' ',i5)
      ENDIF
   END SUBROUTINE geomset_lib
   ! Calls in subroutine GLOBSUM: 
   ! => MPI_ALLREDUCE (on line <477>)
   SUBROUTINE globsum(mat_in,mat_out,mat_size)
      INTEGER(KIND=iintegers),INTENT(IN) :: mat_size
      REAL(KIND=euwp),DIMENSION(mat_size),INTENT(IN) :: mat_in
      REAL(KIND=euwp),DIMENSION(mat_size),INTENT(OUT) :: mat_out
#ifdef PUREMPI
      INTEGER :: ir
      CALL ttbeg(31)
      CALL ttbeg(32)
      CALL MPI_ALLREDUCE(mat_in,mat_out,mat_size,   &
         DC_TYPE,MPI_SUM,mpi_cust_comm,ir)
      IF(ir.ne.0) THEN
              print *,'globsum MPI error'
              STOP
      ENDIF
      CALL ttend(32)
      CALL ttend(31)
#else /*PUREMPI*/
      mat_out=mat_in
#endif

   END SUBROUTINE globsum
   SUBROUTINE globmax(mat_in,mat_out,mat_size)
      INTEGER(KIND=iintegers),INTENT(IN) :: mat_size
      REAL(KIND=euwp),DIMENSION(mat_size),INTENT(IN) :: mat_in
      REAL(KIND=euwp),DIMENSION(mat_size),INTENT(OUT) :: mat_out
#ifdef PUREMPI
      INTEGER :: ir
      CALL ttbeg(31)
      CALL ttbeg(33)
      CALL MPI_ALLREDUCE(mat_in,mat_out,mat_size,   &
         DC_TYPE,MPI_MAX,mpi_cust_comm,ir)
      CALL ttend(33)
      CALL ttend(31)
#else /*PUREMPI*/
      mat_out=mat_in
#endif
   END SUBROUTINE globmax

   SUBROUTINE globmin(mat_in,mat_out,mat_size)
      INTEGER,INTENT(IN) :: mat_size
      REAL(KIND=euwp),DIMENSION(mat_size),INTENT(IN) :: mat_in
      REAL(KIND=euwp),DIMENSION(mat_size),INTENT(OUT) :: mat_out
#ifdef PUREMPI
      INTEGER :: ir
      CALL ttbeg(31)
      CALL ttbeg(34)
      CALL MPI_ALLREDUCE(mat_in,mat_out,mat_size,   &
         DC_TYPE,MPI_MIN,mpi_cust_comm,ir)
      CALL ttend(34)
      CALL ttend(31)
#else /*PUREMPI*/
      mat_out=mat_in
#endif

   END SUBROUTINE globmin

SUBROUTINE globsumaxminv(bufavg,bufmax,bufmin,isize)

  INTEGER(KIND=iintegers), INTENT(IN) :: isize
  REAL(KIND=euwp),INTENT(INOUT) ::  &
           bufavg(isize), &
           bufmax(isize), &
           bufmin(isize) 
  REAL(KIND=euwp) :: array(3*isize),psumaxminv(3*isize)
#ifdef PUREMPI
  INTEGER(KIND=iintegers) :: ir
#endif

    array(        1:  isize)=bufavg(1:isize)
    array(  isize+1:2*isize)=bufmax(1:isize)
    array(2*isize+1:3*isize)=bufmin(1:isize)
#ifdef PUREMPI
    CALL MPI_ALLReduce(array,psumaxminv,3*isize,DC_TYPE,MPI_SUMAXMIN, &
                       mpi_cust_comm,ir)
#else
    psumaxminv=array
#endif
    bufavg(1:isize)=psumaxminv(        1:isize)
    bufmax(1:isize)=psumaxminv(  isize+1:2*isize)
    bufmin(1:isize)=psumaxminv(2*isize+1:3*isize)
END SUBROUTINE globsumaxminv
SUBROUTINE globmaxminv(bufmax,bufmin,isize)

  INTEGER(KIND=iintegers), INTENT(IN) :: isize
  REAL(KIND=euwp),INTENT(INOUT) ::  &
           bufmax(isize), &
           bufmin(isize) 
  REAL(KIND=euwp) :: array(2*isize),pmaxminv(2*isize)
#ifdef PUREMPI
  INTEGER(KIND=iintegers) :: ir
#endif

    array(        1:  isize)=bufmax(1:isize)
    array(  isize+1:2*isize)=bufmin(1:isize)
#ifdef PUREMPI
    CALL MPI_ALLReduce(array,pmaxminv,2*isize,DC_TYPE,MPI_MAXMIN, &
                       mpi_cust_comm,ir)
#else
    pmaxminv=array
#endif
    bufmax(1:isize)=pmaxminv(        1:isize)
    bufmin(1:isize)=pmaxminv(  isize+1:2*isize)
END SUBROUTINE globmaxminv

SUBROUTINE globsumaxminlocv(buffer,isize)

  INTEGER(KIND=iintegers), INTENT(IN) :: isize
  INTEGER(KIND=iintegers),PARAMETER :: idmaxvars=40
  REAL(KIND=euwp),INTENT(INOUT) ::  &
           buffer(15,idmaxvars)
  REAL(KIND=euwp) :: array(15*isize),psumaxminv(15*isize)
  INTEGER(KIND=iintegers) :: i,k
#ifdef PUREMPI
  INTEGER(KIND=iintegers) :: ir
#endif

  psumaxminv(:)=0.

  DO k=0,14
    DO i=1,isize
      array(i+isize*k)=buffer(1+k,i)
    ENDDO
  ENDDO

#ifdef PUREMPI
    CALL MPI_ALLReduce(array,psumaxminv,15*isize,DC_TYPE,MPI_SUMAXMINLOC, &
                       mpi_cust_comm,ir)
#else
    psumaxminv=array
#endif
  DO k=0,14
    DO i=1,isize
      buffer(1+k,i)=psumaxminv(i+isize*k)
    ENDDO
  ENDDO

END SUBROUTINE globsumaxminlocv
!---------------------------------------------------------------------!
! The operator associated with MPI_MAXMIN (sum,max,min)
!---------------------------------------------------------------------!
SUBROUTINE maxmin(wrkin,wrkinout,icnt)
!---------------------------------------------------------------------!

!  REAL(KIND=euwp)  :: sumaxmin
  INTEGER :: icnt
  REAL(KIND=euwp), INTENT(INOUT) :: wrkinout(icnt)
  REAL(KIND=euwp), INTENT(IN)    :: wrkin(icnt)

  INTEGER(KIND=iintegers) :: i,ihalf

  ihalf=icnt/2
  DO i=1,ihalf
    wrkinout(        i)=max(wrkinout(      i),wrkin(      i))
    wrkinout(  ihalf+i)=min(wrkinout(ihalf+i),wrkin(ihalf+i))
  ENDDO

END SUBROUTINE maxmin
!---------------------------------------------------------------------!
! The operator associated with MPI_SUMAXMIN (sum,max,min)
!---------------------------------------------------------------------!
SUBROUTINE sumaxmin(wrkin,wrkinout,icnt)
!---------------------------------------------------------------------!

!  REAL(KIND=euwp)  :: sumaxmin
  INTEGER :: icnt
  REAL(KIND=euwp), INTENT(INOUT) :: wrkinout(icnt)
  REAL(KIND=euwp), INTENT(IN)    :: wrkin(icnt)

  INTEGER(KIND=iintegers) :: i,ithrd

  ithrd=icnt/3
  DO i=1,ithrd
    wrkinout(        i)=wrkin(i)+wrkinout(i)
    wrkinout(  ithrd+i)=max(wrkinout(  ithrd+i),wrkin(  ithrd+i))
    wrkinout(2*ithrd+i)=min(wrkinout(2*ithrd+i),wrkin(2*ithrd+i))
  ENDDO

!  sumaxmin=0._euwp
END SUBROUTINE sumaxmin
!---------------------------------------------------------------------!
! The operator associated with MPI_SUMAXMINLOC (sum,max,min)
!---------------------------------------------------------------------!
SUBROUTINE sumaxminloc(wrkin,wrkinout,icnt)
!---------------------------------------------------------------------!

!  REAL(KIND=euwp)  :: sumaxmin
  INTEGER :: icnt
  REAL(KIND=euwp), INTENT(INOUT) :: wrkinout(icnt)
  REAL(KIND=euwp), INTENT(IN)    :: wrkin(icnt)

  INTEGER(KIND=iintegers) :: i,ithrd

  ithrd=icnt/15
  DO i=1,ithrd
    wrkinout(        i)=wrkin(i)+wrkinout(i)
      IF(wrkinout(ithrd+i) < wrkin(ithrd+i)) THEN
!      wrkinout(ithrd+i)=max(wrkinout(ithrd+i),wrkin(ithrd+i))
        wrkinout(1*ithrd+i)=wrkin(1*ithrd+i)
        wrkinout(3*ithrd+i)=wrkin(3*ithrd+i)
        wrkinout(4*ithrd+i)=wrkin(4*ithrd+i)
        wrkinout(5*ithrd+i)=wrkin(5*ithrd+i)
        wrkinout(6*ithrd+i)=wrkin(6*ithrd+i)
        wrkinout(7*ithrd+i)=wrkin(7*ithrd+i)
        wrkinout(8*ithrd+i)=wrkin(8*ithrd+i)
      ENDIF !storing max val location
      IF(wrkinout(2*ithrd+i) > wrkin(2*ithrd+i)) THEN
!      wrkinout(2*ithrd+i)=min(wrkinout(2*ithrd+i),wrkin(2*ithrd+i))
        wrkinout( 2*ithrd+i)=wrkin( 2*ithrd+i)
        wrkinout( 9*ithrd+i)=wrkin( 9*ithrd+i)
        wrkinout(10*ithrd+i)=wrkin(10*ithrd+i)
        wrkinout(11*ithrd+i)=wrkin(11*ithrd+i)
        wrkinout(12*ithrd+i)=wrkin(12*ithrd+i)
        wrkinout(13*ithrd+i)=wrkin(13*ithrd+i)
        wrkinout(14*ithrd+i)=wrkin(14*ithrd+i)
      ENDIF !storing min val location

  ENDDO

!  sumaxmin=0._euwp
END SUBROUTINE sumaxminloc


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     U P D A T E S    H A L L O    B E T W E E N    P R O C E S S O R S
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------------------------------------------------------------
!
!  Here definition of wrappers that call various setups of update3dsplit
!  Ideally, these should be inlined as they only translate the call
!  Beg and end respectively initiate and finish communication for given variable
!  Note that the number of buffer needs to be provided to avoid
!  accidental memory overwritting(inovar)
!
!  FOR NONPERIODIC BC BORDERS ARE NOT UPDATED !!!!!!!!!
!  THIS IS A SIGNIFICANT CHANGE FROM PREVIOUS VERSIONS AND NEEDS SPECIAL velbc SETUP !!!!!!!
!  CURRENTLY THE SOLUTION IS TO CALL UPDATE FOR INTERIOR AND THEN SEPARATELY FOR THE BORDER
!
!------------------------------------------------------------------
   ! Calls in subroutine UPDATE: 
   ! => ttbeg (on line <507>)
   ! ==> ttend (on line <515>)
   SUBROUTINE update(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(21)
      IF(isphere*ipoles==0) THEN
         CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,2,2,2,1,1,1,0,ih)
      ELSE
         CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
                 1,ldim,ihg,2,2,2,1,1,1,0,ih)
      ENDIF
      CALL ttend(21)
   END SUBROUTINE update
   SUBROUTINE updated(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      DEV_REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(21)
!     IF(isphere*ipoles==0) THEN
         CALL update3dsplit_GPUd(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,2,2,2,1,1,1,0,ih)
!     ELSE
!        CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
!                1,ldim,ihg,2,2,2,1,1,1,0,ih)
!     ENDIF
      CALL ttend(21)
   END SUBROUTINE updated

   SUBROUTINE update_multi(n1,m1,l1,ndim,mdim,ldim,ihg,ih,a,b,c,d,e,f,g,h)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      REAL_euwp,OPTIONAL::b(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      REAL_euwp,OPTIONAL::c(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      REAL_euwp,OPTIONAL::d(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      REAL_euwp,OPTIONAL::e(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      REAL_euwp,OPTIONAL::f(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      REAL_euwp,OPTIONAL::g(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      REAL_euwp,OPTIONAL::h(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(27)
      IF(isphere*ipoles==0) THEN
        ! IF (PRESENT(b)) THEN 
       !  ELSE
        SELECT CASE (updmulgran)
        CASE(1) 
          CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,1,0,ih)
          IF (PRESENT(b)) &
          CALL update3dsplit(b,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,1,0,ih)
          IF (PRESENT(c)) &
          CALL update3dsplit(c,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,1,0,ih)
          IF (PRESENT(d)) &
          CALL update3dsplit(d,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,1,0,ih)
          IF (PRESENT(e)) &
          CALL update3dsplit(e,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,1,0,ih)
          IF (PRESENT(f)) &
          CALL update3dsplit(f,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,1,0,ih)
          IF (PRESENT(g)) &
          CALL update3dsplit(g,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,1,0,ih)
          IF (PRESENT(h)) &
          CALL update3dsplit(h,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,1,0,ih)
       CASE(2) 
         IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b)
         ELSE 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,a)
         ENDIF
         IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,c)
         ENDIF 
 
         IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,e)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,g)
         ENDIF
       CASE(3) 
         IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,a)
         ENDIF
         IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,d,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,d,e)
         ELSE IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,d)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,g)
         ENDIF
       CASE(4) 
         IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,a)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,e,f,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,e,f,g)
         ELSE IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,1,0,ih,e)
         ENDIF 
       CASE(8) 
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b,c,d,e,f,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b,c,d,e,f,g)
         ELSE IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b,c,d,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b,c,d,e)
         ELSE IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,1,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,1,0,ih,a)
         ENDIF
       CASE DEFAULT  
          STOP 'Wrong updmulgran'
       END SELECT
      ELSE
         STOP 'isphere*ipoles not yet available'
      ENDIF
      CALL ttend(27)
   END SUBROUTINE update_multi

!   Now lr,bt,gs and lrbtgs (3D without corners and edges) update wrappers
!   Bor versions of variables update borders only, nothing is done for
!   middle processors
!   begf and endf SUBROUTINEs are wider updates to construct full update
!   with corners and edges in three steps
!   1. updatelrbegf/endf/begendf 2.updatebtbegf/endf 3.updategsf/endf
!
!------------------------------------------------------------------

   ! Calls in subroutine UPDATE3: 
   ! => ttbeg (on line <533>)
   ! ==> ttend (on line <541>)
   SUBROUTINE update3(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(22)
      IF(isphere*ipoles==0) THEN
         CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,2,2,2,1,1,0,0,ih)
      ELSE
         CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
                 1,ldim,ihg,2,2,2,1,1,0,0,ih)
      ENDIF
      CALL ttend(22)
   END SUBROUTINE update3

   SUBROUTINE update3d(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      DEV_REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(22)
!     IF(isphere*ipoles==0) THEN
         CALL update3dsplit_GPUd(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,2,2,2,1,1,0,0,ih)
!     ELSE
!        CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
!                1,ldim,ihg,2,2,2,1,1,0,0,ih)
!     ENDIF
      CALL ttend(22)
   END SUBROUTINE update3d

  SUBROUTINE update3_multi(n1,m1,l1,ndim,mdim,ldim,ihg,ih,a,b,c,d,e,f,g,h)
     INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
     REAL_euwp ::a(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::b(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::c(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::d(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::e(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::f(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::g(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::h(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     CALL ttbeg(28)
     IF(isphere*ipoles==0) THEN
        SELECT CASE (updmulgran)
        CASE(1) 
          CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,0,0,ih)
          IF (PRESENT(b)) &
          CALL update3dsplit(b,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,0,0,ih)
          IF (PRESENT(c)) &
          CALL update3dsplit(c,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,0,0,ih)
          IF (PRESENT(d)) &
          CALL update3dsplit(d,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,0,0,ih)
          IF (PRESENT(e)) &
          CALL update3dsplit(e,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,0,0,ih)
          IF (PRESENT(f)) &
          CALL update3dsplit(f,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,0,0,ih)
          IF (PRESENT(g)) &
          CALL update3dsplit(g,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,0,0,ih)
          IF (PRESENT(h)) &
          CALL update3dsplit(h,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,2,2,1,1,0,0,ih)
       CASE(2) 
         IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b)
         ELSE 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,c)
         ENDIF 
 
         IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,e)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,g)
         ENDIF
       CASE(3) 
         IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,d,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,d,e)
         ELSE IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,d)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,g)
         ENDIF
       CASE(4) 
         IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,e,f,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,e,f,g)
         ELSE IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,2,2,1,1,0,0,ih,e)
         ENDIF 
       CASE(8) 
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b,c,d,e,f,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b,c,d,e,f,g)
         ELSE IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b,c,d,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b,c,d,e)
         ELSE IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,2,2,1,1,0,0,ih,a)
         ENDIF
       CASE DEFAULT  
          STOP 'Wrong updmulgran'
       END SELECT
     ELSE
        STOP 'Not yet available'
     ENDIF
     CALL ttend(28)
  END SUBROUTINE update3_multi

   SUBROUTINE updatelr(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih, 1-ih:ldim+ih)
      CALL ttbeg(23)
      IF(isphere*ipoles==0) THEN
         CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,2,0,0,1,1,0,0,ih)
      ELSE
         CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
                 1,ldim,ihg,2,0,0,1,1,0,0,ih)
      ENDIF
      CALL ttend(23)
   END SUBROUTINE updatelr

   SUBROUTINE updatelrd(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      DEV_REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih, 1-ih:ldim+ih)
      CALL ttbeg(23)
!     IF(isphere*ipoles==0) THEN
         CALL update3dsplit_GPUd(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,2,0,0,1,1,0,0,ih)
!     ELSE
!        CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
!                1,ldim,ihg,2,0,0,1,1,0,0,ih)
!     ENDIF
      CALL ttend(23)
   END SUBROUTINE updatelrd

   SUBROUTINE updatelr_multi(n1,m1,l1,ndim,mdim,ldim,ihg,ih,a,b,c,d,e,f,g,h)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih, 1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::b(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::c(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::d(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::e(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::f(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::g(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::h(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(29)
      IF(isphere*ipoles==0) THEN
        SELECT CASE (updmulgran)
        CASE(1) 
          CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,0,0,1,1,0,0,ih)
          IF (PRESENT(b)) &
          CALL update3dsplit(b,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,0,0,1,1,0,0,ih)
          IF (PRESENT(c)) &
          CALL update3dsplit(c,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,0,0,1,1,0,0,ih)
          IF (PRESENT(d)) &
          CALL update3dsplit(d,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,0,0,1,1,0,0,ih)
          IF (PRESENT(e)) &
          CALL update3dsplit(e,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,0,0,1,1,0,0,ih)
          IF (PRESENT(f)) &
          CALL update3dsplit(f,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,0,0,1,1,0,0,ih)
          IF (PRESENT(g)) &
          CALL update3dsplit(g,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,0,0,1,1,0,0,ih)
          IF (PRESENT(h)) &
          CALL update3dsplit(h,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,2,0,0,1,1,0,0,ih)
       CASE(2) 
         IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b)
         ELSE 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,c)
         ENDIF 
 
         IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,e)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,g)
         ENDIF
       CASE(3) 
         IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,d,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,d,e)
         ELSE IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,d)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,g)
         ENDIF
       CASE(4) 
         IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,e,f,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,e,f,g)
         ELSE IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,2,0,0,1,1,0,0,ih,e)
         ENDIF 
       CASE(8) 
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b,c,d,e,f,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b,c,d,e,f,g)
         ELSE IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b,c,d,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b,c,d,e)
         ELSE IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,0,0,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,2,0,0,1,1,0,0,ih,a)
         ENDIF
       CASE DEFAULT  
          STOP 'Wrong updmulgran'
       END SELECT
      ELSE
        STOP 'Not yet available'
      ENDIF
      CALL ttend(29)
   END SUBROUTINE updatelr_multi

   SUBROUTINE updatelrw(a,n1,n2,m1,m2,l1,l2,ndim1,ndim2,mdim1,mdim2,ldim1,ldim2,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim1,mdim1,ldim1,ihg,ih
      INTEGER(KIND=iintegers),INTENT(IN)::n2,m2,l2,ndim2,mdim2,ldim2
      REAL_euwp :: a(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      CALL ttbeg(75)
      IF(isphere*ipoles==0) THEN
         CALL update3dsplit(a,n1,n2,m1,m2,l1,l2,ndim1,ndim2,mdim1,mdim2,&
             ldim1,ldim2,ihg,2,0,0,1,1,0,0,ih)
      ELSE
         CALL update3dsplit_sphere(a,n1,n2,m1,m2,l1,l2,ndim1,ndim2,mdim1,mdim2,&
             ldim1,ldim2,ihg,2,0,0,1,1,0,0,ih)
      ENDIF
      CALL ttend(75)
   END SUBROUTINE updatelrw

   SUBROUTINE updatebt(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(24)
      IF(isphere*ipoles==0) THEN
         CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,0,2,0,1,1,0,0,ih)
      ELSE
         CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
                 1,ldim,ihg,0,2,0,1,1,0,0,ih)
      ENDIF
      CALL ttend(24)
   END  SUBROUTINE updatebt

   SUBROUTINE updatebtd(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      DEV_REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(24)
!     IF(isphere*ipoles==0) THEN
         CALL update3dsplit_GPUd(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,0,2,0,1,1,0,0,ih)
!     ELSE
!        CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
!                1,ldim,ihg,0,2,0,1,1,0,0,ih)
!     ENDIF
      CALL ttend(24)
   END  SUBROUTINE updatebtd
   SUBROUTINE updatebt_multi(n1,m1,l1,ndim,mdim,ldim,ihg,ih,a,b,c,d,e,f,g,h)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih, 1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::b(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::c(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::d(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::e(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::f(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::g(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::h(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(29)
      IF(isphere*ipoles==0) THEN
        SELECT CASE (updmulgran)
        CASE(1) 
          CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,2,0,1,1,0,0,ih)
          IF (PRESENT(b)) &
          CALL update3dsplit(b,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,2,0,1,1,0,0,ih)
          IF (PRESENT(c)) &
          CALL update3dsplit(c,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,2,0,1,1,0,0,ih)
          IF (PRESENT(d)) &
          CALL update3dsplit(d,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,2,0,1,1,0,0,ih)
          IF (PRESENT(e)) &
          CALL update3dsplit(e,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,2,0,1,1,0,0,ih)
          IF (PRESENT(f)) &
          CALL update3dsplit(f,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,2,0,1,1,0,0,ih)
          IF (PRESENT(g)) &
          CALL update3dsplit(g,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,2,0,1,1,0,0,ih)
          IF (PRESENT(h)) &
          CALL update3dsplit(h,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,2,0,1,1,0,0,ih)
       CASE(2) 
         IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b)
         ELSE 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,c)
         ENDIF 
 
         IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,e)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,g)
         ENDIF
       CASE(3) 
         IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,d,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,d,e)
         ELSE IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,d)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,g)
         ENDIF
       CASE(4) 
         IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,e,f,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,e,f,g)
         ELSE IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,2,0,1,1,0,0,ih,e)
         ENDIF 
       CASE(8) 
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b,c,d,e,f,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b,c,d,e,f,g)
         ELSE IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b,c,d,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b,c,d,e)
         ELSE IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,2,0,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,2,0,1,1,0,0,ih,a)
         ENDIF
       CASE DEFAULT  
          STOP 'Wrong updmulgran'
       END SELECT
      ELSE
        STOP 'Not yet available'
      ENDIF
      CALL ttend(29)
   END SUBROUTINE updatebt_multi

   SUBROUTINE updatebtw(a,n1,n2,m1,m2,l1,l2,ndim1,ndim2,mdim1,mdim2,ldim1,ldim2,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim1,mdim1,ldim1,ihg,ih
      INTEGER(KIND=iintegers),INTENT(IN)::n2,m2,l2,ndim2,mdim2,ldim2
      REAL_euwp :: a(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      CALL ttbeg(75)
      IF(isphere*ipoles==0) THEN
         CALL update3dsplit(a,n1,n2,m1,m2,l1,l2,ndim1,ndim2,mdim1,mdim2,&
             ldim1,ldim2,ihg,0,2,0,1,1,0,0,ih)
      ELSE
         CALL update3dsplit_sphere(a,n1,n2,m1,m2,l1,l2,ndim1,ndim2,mdim1,mdim2,&
             ldim1,ldim2,ihg,0,2,0,1,1,0,0,ih)
      ENDIF
      CALL ttend(75)
   END  SUBROUTINE updatebtw

   SUBROUTINE updatelrbt(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih, 1-ih:ldim+ih)
      CALL ttbeg(23)
      IF(isphere*ipoles==0) THEN
         CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,2,2,0,1,1,0,0,ih)
      ELSE
         CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
                 1,ldim,ihg,2,2,0,1,1,0,0,ih)
      ENDIF
      CALL ttend(23)
   END SUBROUTINE updatelrbt

   SUBROUTINE updategs(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(25)
      IF(isphere*ipoles==0) THEN
         CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,0,0,2,1,1,0,0,ih)
      ELSE
         CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
                 1,ldim,ihg,0,0,2,1,1,0,0,ih)
      ENDIF
      CALL ttend(25)
   END SUBROUTINE updategs

   SUBROUTINE updategsd(a,n1,m1,l1,ndim,mdim,ldim,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      DEV_REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(25)
!     IF(isphere*ipoles==0) THEN
         CALL update3dsplit_GPUd(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                 1,ldim,ihg,0,0,2,1,1,0,0,ih)
!     ELSE
!        CALL update3dsplit_sphere(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,         &
!                1,ldim,ihg,0,0,2,1,1,0,0,ih)
!     ENDIF
      CALL ttend(25)
   END SUBROUTINE updategsd

   SUBROUTINE updategs_multi(n1,m1,l1,ndim,mdim,ldim,ihg,ih,a,b,c,d,e,f,g,h)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim,mdim,ldim,ihg,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih, 1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::b(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::c(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::d(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::e(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::f(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::g(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
     REAL_euwp,OPTIONAL::h(1-ih:ndim+ih, 1-ih:mdim+ih,1-ih:ldim+ih)
      CALL ttbeg(29)
      IF(isphere*ipoles==0) THEN
        SELECT CASE (updmulgran)
        CASE(1) 
          CALL update3dsplit(a,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,0,2,1,1,0,0,ih)
          IF (PRESENT(b)) &
          CALL update3dsplit(b,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,0,2,1,1,0,0,ih)
          IF (PRESENT(c)) &
          CALL update3dsplit(c,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,0,2,1,1,0,0,ih)
          IF (PRESENT(d)) &
          CALL update3dsplit(d,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,0,2,1,1,0,0,ih)
          IF (PRESENT(e)) &
          CALL update3dsplit(e,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,0,2,1,1,0,0,ih)
          IF (PRESENT(f)) &
          CALL update3dsplit(f,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,0,2,1,1,0,0,ih)
          IF (PRESENT(g)) &
          CALL update3dsplit(g,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,0,2,1,1,0,0,ih)
          IF (PRESENT(h)) &
          CALL update3dsplit(h,1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                  1,ldim,ihg,0,0,2,1,1,0,0,ih)
       CASE(2) 
         IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b)
         ELSE 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,c)
         ENDIF 
 
         IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,e)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,g)
         ENDIF
       CASE(3) 
         IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,d,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,d,e)
         ELSE IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,d)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,g)
         ENDIF
       CASE(4) 
         IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,a)
         ENDIF
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,e,f,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,e,f,g)
         ELSE IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                    1,ldim,ihg,0,0,2,1,1,0,0,ih,e)
         ENDIF 
       CASE(8) 
         IF (PRESENT(h)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b,c,d,e,f,g,h)
         ELSE IF (PRESENT(g)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b,c,d,e,f,g)
         ELSE IF (PRESENT(f)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b,c,d,e,f)
         ELSE IF (PRESENT(e)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b,c,d,e)
         ELSE IF (PRESENT(d)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b,c,d)
         ELSE IF (PRESENT(c)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b,c)
         ELSE IF (PRESENT(b)) THEN 
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,0,2,1,1,0,0,ih,a,b)
         ELSE  
           CALL update3dsplit_multi(1,n1,1,m1,1,l1,1,ndim,1,mdim,                &
                   1,ldim,ihg,0,0,2,1,1,0,0,ih,a)
         ENDIF
         CASE DEFAULT  
          STOP 'Wrong updmulgran'
         END SELECT
      ELSE
        STOP 'Not yet available'
      ENDIF
      CALL ttend(29)
   END SUBROUTINE updategs_multi

   SUBROUTINE updategsw(a,n1,n2,m1,m2,l1,l2,ndim1,ndim2,mdim1,mdim2,ldim1,ldim2,ihg,ih)
      INTEGER(KIND=iintegers),INTENT(IN)::n1,m1,l1,ndim1,mdim1,ldim1,ihg,ih
      INTEGER(KIND=iintegers),INTENT(IN)::n2,m2,l2,ndim2,mdim2,ldim2
      REAL_euwp :: a(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      CALL ttbeg(75)
      IF(isphere*ipoles==0) THEN
         CALL update3dsplit(a,n1,n2,m1,m2,l1,l2,ndim1,ndim2,mdim1,mdim2,&
             ldim1,ldim2,ihg,0,0,2,1,1,0,0,ih)
      ELSE
         CALL update3dsplit_sphere(a,n1,n2,m1,m2,l1,l2,ndim1,ndim2,mdim1,mdim2,&
             ldim1,ldim2,ihg,0,0,2,1,1,0,0,ih)
      ENDIF
      CALL ttend(75)
   END SUBROUTINE updategsw

!---------------------------------------------------------------------!
!       Description of iflg options
!       0 - Full 2D update
!       1 - EW 2D update
!       2 - NS 2D update
!----------------------------------------------------------------------!
   SUBROUTINE updateflt(a,ihg,iflg,ndim,mdim,ih)
!----------------------------------------------------------------------!

      INTEGER(KIND=iintegers), INTENT(IN) :: ihg,iflg,ndim,mdim,ih
      REAL_euwp :: a(1-ih:ndim+ih, 1-ih:mdim+ih)
      
      CALL ttbeg(26)

      IF(isphere*ipoles==0) THEN
        IF (IFlg == 0) THEN
           CALL update3dsplit(a,1_iintegers,ndim,1_iintegers,mdim,1_iintegers,1_iintegers,  &
                                1_iintegers,ndim,1_iintegers,mdim,1_iintegers+ih,1_iintegers-ih,ihg, &
                             2_iintegers,2_iintegers,0_iintegers,1_iintegers,1_iintegers,0_iintegers,0_iintegers,ih)
        ELSEIF(iflg == 1) THEN
           CALL update3dsplit(a,1,ndim,1,mdim,1,1,1,ndim,1,mdim, &
              1+ih,1-ih,ihg,0,2,0,1,1,0,0,ih)
        ELSEIF(iflg == 2) THEN
           CALL update3dsplit(a,1,ndim,1,mdim,1,1,1,ndim,1,mdim, &
              1+ih,1-ih,ihg,2,0,0,1,1,0,0,ih)
        ENDIF
      ELSE !isphere*ipoles
        IF (IFlg == 0) THEN
           CALL update3dsplit_sphere(a,1,ndim,1,mdim,1,1,1,ndim,1,mdim, &
              1+ih,1-ih,ihg,3,3,0,1,1,0,0,ih)
        ELSEIF(iflg == 1) THEN
           CALL update3dsplit_sphere(a,1,ndim,1,mdim,1,1,1,ndim,1,mdim, &
              1+ih,1-ih,ihg,3,3,0,1,1,0,0,ih)
        ELSEIF(iflg == 2) THEN
           CALL update3dsplit_sphere(a,1,ndim,1,mdim,1,1,1,ndim,1,mdim, &
              1+ih,1-ih,ihg,3,3,0,1,1,0,0,ih)
        ENDIF
      ENDIF

      CALL ttend(26)
   END SUBROUTINE updateflt
! End of update wrappers
#if(1==0)
   SUBROUTINE update3dsplit_GPUd                               &
    (a,nn1,nn2,mm1,mm2,ll1,ll2,ndim1,ndim2,mdim1,mdim2,               &
     ldim1,ldim2,ihg,icomx0,icomy0,icomz0,istart,ifinsh,              &
     isful0,inovar,ih) !iopt,ioptz,imode0,inovar)
      USE scratch_datafields, ONLY: &
      tmpxsnd1,tmpxsnd2,tmpxrcv1,tmpxrcv2, & 
      tmpysnd1,tmpysnd2,tmpyrcv1,tmpyrcv2, & 
      tmpzsnd1,tmpzsnd2,tmpzrcv1,tmpzrcv2  
!
!     this SUBROUTINE updates the halo (or ghost) cells surrounding
!     each processors subgrid
!
      INTEGER(KIND=iintegers) :: nn1,nn2,mm1,mm2,ll1,ll2,              &
        ndim1,ndim2,mdim1,mdim2,ldim1,ldim2,ihg,                         &
        icomx0,icomy0,icomz0,icomx,icomy,icomz,istart,ifinsh,isful,inovar
      INTEGER(KIND=iintegers) :: n1,n2,m1,m2,l1,l2,                    &
         nllim,mllim,lllim, nulim,mulim,lulim
      INTEGER nhl,mhl,lhl
      INTEGER(KIND=iintegers) :: ih,i,j,k,ihr
      INTEGER(KIND=iintegers) :: icnt,isful0
      INTEGER(KIND=iintegers) :: istat
      DEV_REAL_euwp,INTENT(INOUT) :: &
         a(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
#ifdef PUREMPI
      INTEGER ierr,itagoffset
      INTEGER statsx1(0:msz),statsx2(0:msz)
      INTEGER statrx1(0:msz),statrx2(0:msz)
      INTEGER statsy1(0:msz),statsy2(0:msz)
      INTEGER statry1(0:msz),statry2(0:msz)
      INTEGER statsz1(0:msz),statsz2(0:msz)
      INTEGER statrz1(0:msz),statrz2(0:msz)
      SAVE statsx1,statsx2,statrx1,statrx2
      SAVE statsy1,statsy2,statry1,statry2
      SAVE statsz1,statsz2,statrz1,statrz2
#endif /*PUREMPI*/
      LOGICAL commx,commy,commz
      LOGICAL readsenddataxL,recvwrtedataxL
      LOGICAL readsenddataxR,recvwrtedataxR
      LOGICAL readsenddatayB,recvwrtedatayB
      LOGICAL readsenddatayT,recvwrtedatayT
      LOGICAL readsenddatazG,recvwrtedatazG
      LOGICAL readsenddatazS,recvwrtedatazS
      LOGICAL timecounterdone 
#ifdef DEBUGMPI
      INTEGER iprint,iprr
#endif
#ifdef PUREMPI
      INTEGER status(MPI_STATUS_SIZE)
      timecounterdone = .FALSE.
      itagoffset=1000*isful0+icomx0*100+icomy0*300 
#endif /*PUREMPI*/
      IF(inovar.ne.0) STOP 'Wrong inovar'
      IF(ih<3.AND.ihg>=3)                                         &
        STOP 'Trying to update with too large halo extent 3'
      IF(ih<2.AND.ihg>=2)                                         &
        STOP 'Trying to update with too large halo extent 2'
!     ierr=cudaDeviceSynchronize()
!     ierr=cudaStreamSynchronize()
      ihr=min(ih,ihg)
      ihr=1
      IF(ih==0)                                         &
        STOP 'Calling update with zero halo size. Probably a bug'
      n1=nn1
      n2=nn2
      m1=mm1
      m2=mm2
      l1=ll1
      l2=ll2
      nllim=n1-1-ihr
      mllim=m1-1-ihr
      lllim=l1-1-ihr
      nulim=n2  -ihr
      mulim=m2  -ihr
      lulim=l2  -ihr
      isful=isful0*0+1
      nhl=(n2-n1+1)*ihr*(l2-l1+1)
      mhl=ihr*(m2-m1+1+2*ihr*isful)*(l2-l1+1)
      lhl=(n2-n1+1+2*ihr*isful)*(m2-m1+1+2*ihr*isful)*ihr

!     IF(mype == 0)     PRINT *,'Entering updateGPUd' 
!     CALL flush
!     CALL mybarierr()
 
#ifdef DEBUGMPI
      iprint=0
!      IF(mype == 0)     PRINT 98,nprocx,nprocy,nprocz,nhl,mhl,lhl
      CALL flush(6)
         CALL MPI_BARRIER(mpi_cust_comm,ierr)
      IF (iprint==1) THEN
       CALL MPI_Comm_size(mpi_cust_comm, mysize, ierr)  ! total numbers of PE''s
         DO iprr=0,mysize
            IF (mype==iprr) THEN
               PRINT 99,mype,rightedge,leftedge, &
                         botedge,topedge, &
                         npos,mpos,lpos,                          &
                         peW,peE,peN,peS,n2-n1+1,m2-m1+1,l2-l1+1
               CALL flush(6)
            ENDIF
         CALL MPI_BARRIER(mpi_cust_comm,ierr)
         ENDDO
98       FORMAT('nprocx = ',i3,'  nprocy = ',i3,'  nprocz = ',i3,' nhl =',i3,' mhl = ',i3,' lhl'/, &
          ' mype R   L   B   T   N   M   L    pW    pE    pN    pS  np  mp lp')
99       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',                               &
                i3,' ',i3,' ',i3,' ',i3,' ', &
                i5,' ',i5,' ',i5,' ',i5,' ',i5,' ',i5,' ',i5,' ')
      ENDIF
#endif

!     there are four border segments that must be sent and received.
!     first, we copy from the large array to the tmp arrays for
!     sending. next, we send our four border segments to surrounding
!     processors. then, we receive the four segments from the
!     surrounding processors and we copy that data to the big
!     array.

!Description of the communication choice logic:
!------------------------------------------------------------------
! Generally, we consider three main types of communication mainly for nonperiodic b.c.
! (For periodic b.c. we usually need to communicate everything, so there is less choice)
! 1. Communicate only between border outermost processor and update only
!    the rightmost (topmost, skymost) processor - ideal for cyclicity enforcement -icomx,y,z=1
! 2. Communication inside the domain for nonperiodic b.c., but do not
!      communicate by default between border outermost processors - icomx,y,z=2
! 3. Communication everywhere (as traditionally done in previous versions of EULAG) icomx,y,z=3
!------------------------------------------------------------------
      icomx=icomx0
      icomy=icomy0
      icomz=icomz0
      IF(ibcx==1.AND.icomx==2) icomx=3
      IF(ibcy==1.AND.icomy==2) icomy=3
      IF(ibcz==1.AND.icomz==2) icomz=3
      IF(isphere==1.AND.ipoles == 2.AND.icomx>0) icomx=3
      IF(isphere==1.AND.ipoles == 2.AND.icomy>0) icomy=3
!      if(icomx.gt.0) icomx=3
!      if(icomy.gt.0) icomy=3
!      if(icomz.gt.0) icomz=3
!Set boolean flags for communication in x,y and z directions respectively.
      readsenddataxL=(icomx==1.AND.leftedge==1).OR.                 &
                     (icomx==2.AND.leftedge/=1).OR.                 &
                      icomx==3
      readsenddataxR=(icomx/=1).AND.                                &
                    ((icomx==2 .AND.rightedge/=1).OR.               &
                      icomx==3)
      readsenddataZG=(icomz==1.AND.gndedge==1).OR.                  &
                     (icomz==2.AND.gndedge/=1).OR.                  &
                      icomz==3
      readsenddataZS=(icomz/=1).AND.                                &
                    ((icomz==2 .AND.skyedge/=1).OR.                 &
                      icomz==3)
      recvwrtedataxL=(icomx/=1).AND.                                &
                    ((icomx==2 .AND.leftedge/=1).OR.                &
                      icomx==3)
      recvwrtedataxR=(icomx==1.AND.rightedge==1).OR.                &
                     (icomx==2.AND.rightedge/=1).OR.                &
                      icomx==3
      recvwrtedataZG=(icomz/=1).AND.                                &
                    ((icomz==2 .AND.gndedge/=1).OR.                 &
                      icomz==3)
      recvwrtedataZS=(icomz==1.AND.skyedge==1).OR.                  &
                     (icomz==2.AND.skyedge/=1).OR.                  &
                      icomz==3
      readsenddatayB=(icomy==1.AND.botedge==1).OR.                  &
                     (icomy==2.AND.botedge/=1).OR.                  &
                      icomy==3
      readsenddatayT=(icomy/=1).AND.                                &
                    ((icomy==2 .AND.(topedge/=1)).OR.               &
                      icomy==3)
      recvwrtedatayB=(icomy/=1).AND.                                &
                    ((icomy==2 .AND.(botedge/=1)).OR.               &
                      icomy==3)
      recvwrtedatayT=(icomy==0.AND.topedge==1).OR.                  &
                     (icomy==2.AND.topedge/=1).OR.                  &
                      icomy==3
!Set boolean flags for initiating (sending) or finishing (receiving) message
!First flags to send message depending on the position on the proc grid and icomxyz flag
!For icomxyz=1 only the left,bottom and ground procs send messages
      commx=icomx>0
      commy=icomy>0
      commz=icomz>0
      IF(commy) THEN !icomy > 0 so we are communication in y direction`
!If starting communication
         IF(istart==1) THEN
            IF(readsenddatayB) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                        tmpysnd1(icnt,inovar)=a(i,      j,k)
                     END DO
                  END DO
               END DO
            ENDIF !bottom buffer y direction
            IF(readsenddatayT) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare top data segments for processor on the bottom:tmpysnd2
                        icnt=icnt+1
                        tmpysnd2(icnt,inovar)=a(i,mulim+j,k)
                     END DO
                  END DO
               END DO
            ENDIF ! top buffer y direction
         ENDIF !istart
!
!     exchange data now
!
         IF(nprocy>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN

               IF(istart==1) THEN
#ifdef DEBUGMPI
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
                  IF(readsenddatayB) THEN                                                 
                   CALL MPI_ISEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statsy1(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to snd ',nhl,'bytes to  ', peS, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayT) THEN                                                 
                   CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statry1(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peN, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayT) THEN                                                
                   CALL MPI_ISEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statsy2(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to snd ',nhl,'bytes to  ', peN, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN                                                
                   CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statry2(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peS, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
               ENDIF !istart
#ifdef DEBUGMPI
!               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
               IF(ifinsh==1) THEN
                  IF(readsenddatayB) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'waits to send',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                        CALL MPI_WAIT(statsy1(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'has sent',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayT) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'waits to send',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                         CALL MPI_WAIT(statsy2(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'has sent',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayT) THEN 
#ifdef DEBUGMPI
                    PRINT *,mype, 'expects ',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                        CALL MPI_WAIT(statry1(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'received',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'expects ',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                         CALL MPI_WAIT(statry2(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'received',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
               ENDIF !ifinsh
#ifdef DEBUGMPI
                CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedatayT) THEN 
                   CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statry1(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peN, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN  
                   CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statry2(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'prepared to rcv ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayB)  THEN
#ifdef DEBUGMPI
!                   PRINT *,mype,' sending ',nhl,'bytes to ', peS 
!                   CALL flush(6)
#endif
                   CALL  MPI_SEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,ierr)
#ifdef DEBUGMPI
!                  PRINT *,mype,' has sent ',nhl,'bytes to ', peS,' with ', ierr 
!                   CALL flush(6)
#endif
                   !PRINT *,mype,' sending err',ierr
                   !CALL flush(6)
                  ENDIF
                  IF(readsenddatayT) THEN 
#ifdef DEBUGMPI
!                   PRINT *,mype,' sending ',nhl,'bytes to ', peN 
!                   CALL flush(6)
#endif
                   CALL  MPI_SEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' has sent ',nhl,'bytes to ', peN,' with ', ierr 
!                   CALL flush(6)
#endif
                   !PRINT *,mype,' sending err',ierr
                   !CALL flush(6)
                  ENDIF
               ENDIF !istart
               IF(ifinsh==1) THEN
!                  print *,mype,'Before y wait'
                  IF(recvwrtedatayT) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry1(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry2(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
!                  print *,mype,'After y wait'
               ENDIF !ifinsh
            ENDIF
            CALL ttend(30)
            timecounterdone=.TRUE.
#endif /*PUREMPI*/
         ELSE ! nprocy=1
            IF(recvwrtedatayB) THEN
               DO i=1,nhl
                  tmpyrcv2(i,inovar)=tmpysnd2(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedatayT) THEN
               DO i=1,nhl
                  tmpyrcv1(i,inovar)=tmpysnd1(i,inovar)
               ENDDO
            ENDIF
         ENDIF ! nprocy=1
         IF(ifinsh==1) THEN
            IF(recvwrtedatayB) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        a(i,mllim+j,k)=tmpyrcv2(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatayB
            IF(recvwrtedatayT) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        a(i,   m2+j,k)=tmpyrcv1(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatayT
         ENDIF !ifinsh
      ENDIF !commy >0 Communication in y direction

      IF(commx) THEN ! commx >0 Communication in x direction
         IF(istart==1) THEN
            IF(readsenddataxL) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1(icnt,inovar)=a(      i,j,k)
                     END DO
                  END DO
               END DO
            ENDIF !readsenddataxL
            IF(readsenddataxR) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2(icnt,inovar)=a(nulim+i,j,k)
                     END DO
                  END DO
               END DO
            ENDIF !readsenddataxR
         ENDIF !istart
         IF(nprocx>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
               IF(istart==1) THEN
                  IF(readsenddataxL)                                                &
                   CALL MPI_ISEND(tmpxsnd1(1,inovar),mhl,DC_TYPE,peW,itagoffset+3+100*inovar    &
                                  ,mpi_cust_comm,statsx1(inovar),ierr)
                  IF(recvwrtedataxR)                                                &
                   CALL MPI_IRECV(tmpxrcv1(1,inovar),mhl,DC_TYPE,peE,itagoffset+3+100*inovar    &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxR)                                                &
                   CALL MPI_ISEND(tmpxsnd2(1,inovar),mhl,DC_TYPE,peE,itagoffset+4+100*inovar    &
                                  ,mpi_cust_comm,statsx2(inovar),ierr)
                  IF(recvwrtedataxL)                                                &
                   CALL MPI_IRECV(tmpxrcv2(1,inovar),mhl,DC_TYPE,peW,itagoffset+4+100*inovar    &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
               ENDIF !istart
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
               IF(ifinsh==1) THEN
                  IF(readsenddataxL) CALL MPI_WAIT(statsx1(inovar),status,ierr)
                  IF(readsenddataxR) CALL MPI_WAIT(statsx2(inovar),status,ierr)
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
               ENDIF !ifinsh
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedataxR)                                               &
                   CALL MPI_IRECV(tmpxrcv1(1,inovar),mhl,DC_TYPE,peE,itagoffset+3+100*inovar   &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxL)                                               &
                   CALL  MPI_SEND(tmpxsnd1(1,inovar),mhl,DC_TYPE,peW,itagoffset+3+100*inovar   &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedataxL)                                               &
                   CALL MPI_IRECV(tmpxrcv2(1,inovar),mhl,DC_TYPE,peW,itagoffset+4+100*inovar   &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
                  IF(readsenddataxR)                                               &
                   CALL  MPI_SEND(tmpxsnd2(1,inovar),mhl,DC_TYPE,peE,itagoffset+4+100*inovar   &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
!                  print *,mype,'Before x wait'
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
!                  print *,mype,'After x wait'

               ENDIF !ifinsh
            ENDIF
            CALL ttend(30,timecounterdone)
            timecounterdone=.TRUE.
#endif /*PURMPI*/
         ELSE !now nprocx.eq.1
            IF(recvwrtedataxR) THEN
               DO i=1,mhl
                  tmpxrcv1(i,inovar)=tmpxsnd1(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedataxL) THEN
               DO i=1,mhl
                  tmpxrcv2(i,inovar)=tmpxsnd2(i,inovar)
               ENDDO
            ENDIF
         ENDIF !nprocx

         IF(ifinsh==1) THEN
            IF(recvwrtedataxL) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        a(nllim+i,j,k)=tmpxrcv2(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedataxL
            IF(recvwrtedataxR) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        a(   n2+i,j,k)=tmpxrcv1(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedataxR
         ENDIF !ifinsh
      ENDIF !commx

      IF(commz) THEN ! commz >0 Communication in z direction

!        IF(istart==1) THEN
!           IF(readsenddatazG) THEN
!              icnt=1
!              DO k=1,ihr
!                 DO j=m1-ihr*isful,m2+ihr*isful
!                    DO i=n1-ihr*isful,n2+ihr*isful
!                       tmpzsnd1(icnt,inovar)=a(i,j,      k)
!                       icnt=icnt+1
!                    END DO
!                 END DO
!              END DO
!           ENDIF !readsenddatazG
!           IF(readsenddatazS) THEN
!              icnt=1
!              DO k=1,ihr
!                 DO j=m1-ihr*isful,m2+ihr*isful
!                    DO i=n1-ihr*isful,n2+ihr*isful
!                       tmpzsnd2(icnt,inovar)=a(i,j,lulim+k)
!                       icnt=icnt+1
!                    END DO
!                 END DO
!              END DO
!           ENDIF !readsenddatazS
!        ENDIF !istart

         IF(nprocz>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN
               IF(istart==1) THEN
                  IF(readsenddatazG)                                                &
                   CALL MPI_ISEND(a(1-ih,1-ih,1),lhl,DC_TYPE,peG,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statsz1(inovar),ierr)
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(a(1-ih,1-ih,l2+1),lhl,DC_TYPE,peZ,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL MPI_ISEND(a(1-ih,1-ih,lulim+1),lhl,DC_TYPE,peZ,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statsz2(inovar),ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(a(1-ih,1-ih,lllim+1),lhl,DC_TYPE,peG,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(readsenddatazG) CALL MPI_WAIT(statsz1(inovar),status,ierr)
                  IF(readsenddatazS) CALL MPI_WAIT(statsz2(inovar),status,ierr)
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh

            ENDIF
            IF (isend==2) THEN
               IF(istart==1) THEN
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(a(1-ih,1-ih,l2+1),lhl,DC_TYPE,peZ,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazG)                                                &
                   CALL  MPI_SEND(a(1-ih,1-ih,1),lhl,DC_TYPE,peG,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(a(1-ih,1-ih,lllim+1),lhl,DC_TYPE,peG,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL  MPI_SEND(a(1-ih,1-ih,lulim+1),lhl,DC_TYPE,peZ,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh
            ENDIF
            CALL ttend(30,timecounterdone)
#endif /*PUREMPI*/
         ELSE !now nprocz.eq.1
!           IF(recvwrtedatazS) THEN
!              DO i=1,lhl
!                 tmpzrcv1(i,inovar)=tmpzsnd1(i,inovar)
!              ENDDO
!           ENDIF
!           IF(recvwrtedatazG) THEN
!              DO i=1,lhl
!                 tmpzrcv2(i,inovar)=tmpzsnd2(i,inovar)
!              ENDDO
!           ENDIF
         ENDIF
!
!     store data in main array now
!


!        IF(ifinsh==1) THEN
!           IF(recvwrtedatazG) THEN
!              icnt=1
!              DO k=1,ihr
!                 DO j=m1-ihr*isful,m2+ihr*isful
!                    DO i=n1-ihr*isful,n2+ihr*isful
!                       a(i,j,lllim+k)=tmpzrcv2(icnt,inovar)
!                       icnt=icnt+1
!                    END DO
!                 END DO
!              END DO
!           ENDIF !recvwrtedatazG
!           IF(recvwrtedatazS) THEN
!              icnt=1
!              DO k=1,ihr
!                 DO j=m1-ihr*isful,m2+ihr*isful
!                    DO i=n1-ihr*isful,n2+ihr*isful
!                       a(i,j,   l2+k)=tmpzrcv1(icnt,inovar)
!                       icnt=icnt+1
!                    END DO
!                 END DO
!              END DO
!           ENDIF !recvwrtedatazZ
!        ENDIF !ifinsh
      ENDIF !commz
!  ierr=cudaDeviceSynchronize()
   END SUBROUTINE update3dsplit_GPUd                                      

   SUBROUTINE update3dsplit                                      &
    (a,nn1,nn2,mm1,mm2,ll1,ll2,ndim1,ndim2,mdim1,mdim2,               &
     ldim1,ldim2,ihg,icomx0,icomy0,icomz0,istart,ifinsh,              &
     isful,inovar,ih) !iopt,ioptz,imode0,inovar)
      USE scratch_datafields, ONLY: &
      tmpxsnd1,tmpxsnd2,tmpxrcv1,tmpxrcv2, & 
      tmpysnd1,tmpysnd2,tmpyrcv1,tmpyrcv2, & 
      tmpzsnd1,tmpzsnd2,tmpzrcv1,tmpzrcv2  
!
!     this SUBROUTINE updates the halo (or ghost) cells surrounding
!     each processors subgrid
!
      INTEGER(KIND=iintegers) :: nn1,nn2,mm1,mm2,ll1,ll2,              &
        ndim1,ndim2,mdim1,mdim2,ldim1,ldim2,ihg,                         &
        icomx0,icomy0,icomz0,icomx,icomy,icomz,istart,ifinsh,isful,inovar
      INTEGER(KIND=iintegers) :: n1,n2,m1,m2,l1,l2,                    &
         nllim,mllim,lllim, nulim,mulim,lulim
      INTEGER nhl,mhl,lhl
      INTEGER(KIND=iintegers) :: ih,i,j,k,ihr
      INTEGER(KIND=iintegers) :: icnt
      REAL_euwp,INTENT(INOUT) :: &
         a(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
#ifdef PUREMPI
      INTEGER ierr,itagoffset
      INTEGER statsx1(0:msz),statsx2(0:msz)
      INTEGER statrx1(0:msz),statrx2(0:msz)
      INTEGER statsy1(0:msz),statsy2(0:msz)
      INTEGER statry1(0:msz),statry2(0:msz)
      INTEGER statsz1(0:msz),statsz2(0:msz)
      INTEGER statrz1(0:msz),statrz2(0:msz)
      SAVE statsx1,statsx2,statrx1,statrx2
      SAVE statsy1,statsy2,statry1,statry2
      SAVE statsz1,statsz2,statrz1,statrz2
#endif /*PUREMPI*/
      LOGICAL commx,commy,commz
      LOGICAL readsenddataxL,recvwrtedataxL
      LOGICAL readsenddataxR,recvwrtedataxR
      LOGICAL readsenddatayB,recvwrtedatayB
      LOGICAL readsenddatayT,recvwrtedatayT
      LOGICAL readsenddatazG,recvwrtedatazG
      LOGICAL readsenddatazS,recvwrtedatazS
      LOGICAL timecounterdone 
#ifdef DEBUGMPI
      INTEGER iprint,iprr
#endif
#ifdef PUREMPI
      INTEGER status(MPI_STATUS_SIZE)
      timecounterdone = .FALSE.
      itagoffset=1000*isful+icomx0*100+icomy0*300 
#endif /*PUREMPI*/
      IF(inovar.ne.0) STOP 'Wrong inovar'
      IF(ih<3.AND.ihg>=3)                                         &
        STOP 'Trying to update with too large halo extent 3'
      IF(ih<2.AND.ihg>=2)                                         &
        STOP 'Trying to update with too large halo extent 2'
      ihr=min(ih,ihg)
      IF(ih==0)                                         &
        STOP 'Calling update with zero halo size. Probably a bug'
      n1=nn1
      n2=nn2
      m1=mm1
      m2=mm2
      l1=ll1
      l2=ll2
      nllim=n1-1-ihr
      mllim=m1-1-ihr
      lllim=l1-1-ihr
      nulim=n2  -ihr
      mulim=m2  -ihr
      lulim=l2  -ihr

      nhl=(n2-n1+1)*ihr*(l2-l1+1)
      mhl=ihr*(m2-m1+1+2*ihr*isful)*(l2-l1+1)
      lhl=(n2-n1+1+2*ihr*isful)*(m2-m1+1+2*ihr*isful)*ihr
#ifdef DEBUGMPI
      iprint=0
!      IF(mype == 0)     PRINT 98,nprocx,nprocy,nprocz,nhl,mhl,lhl
      CALL flush(6)
         CALL MPI_BARRIER(mpi_cust_comm,ierr)
      IF (iprint==1) THEN
       CALL MPI_Comm_size(mpi_cust_comm, mysize, ierr)  ! total numbers of PE''s
         DO iprr=0,mysize
            IF (mype==iprr) THEN
               PRINT 99,mype,rightedge,leftedge, &
                         botedge,topedge, &
                         npos,mpos,lpos,                          &
                         peW,peE,peN,peS,n2-n1+1,m2-m1+1,l2-l1+1
               CALL flush(6)
            ENDIF
         CALL MPI_BARRIER(mpi_cust_comm,ierr)
         ENDDO
98       FORMAT('nprocx = ',i3,'  nprocy = ',i3,'  nprocz = ',i3,' nhl =',i3,' mhl = ',i3,' lhl'/, &
          ' mype R   L   B   T   N   M   L    pW    pE    pN    pS  np  mp lp')
99       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',                               &
                i3,' ',i3,' ',i3,' ',i3,' ', &
                i5,' ',i5,' ',i5,' ',i5,' ',i5,' ',i5,' ',i5,' ')
      ENDIF
#endif

!     there are four border segments that must be sent and received.
!     first, we copy from the large array to the tmp arrays for
!     sending. next, we send our four border segments to surrounding
!     processors. then, we receive the four segments from the
!     surrounding processors and we copy that data to the big
!     array.

!Description of the communication choice logic:
!------------------------------------------------------------------
! Generally, we consider three main types of communication mainly for nonperiodic b.c.
! (For periodic b.c. we usually need to communicate everything, so there is less choice)
! 1. Communicate only between border outermost processor and update only
!    the rightmost (topmost, skymost) processor - ideal for cyclicity enforcement -icomx,y,z=1
! 2. Communication inside the domain for nonperiodic b.c., but do not
!      communicate by default between border outermost processors - icomx,y,z=2
! 3. Communication everywhere (as traditionally done in previous versions of EULAG) icomx,y,z=3
!------------------------------------------------------------------
      icomx=icomx0
      icomy=icomy0
      icomz=icomz0
      IF(ibcx==1.AND.icomx==2) icomx=3
      IF(ibcy==1.AND.icomy==2) icomy=3
      IF(ibcz==1.AND.icomz==2) icomz=3
      IF(isphere==1.AND.ipoles == 2.AND.icomx>0) icomx=3
      IF(isphere==1.AND.ipoles == 2.AND.icomy>0) icomy=3
!      if(icomx.gt.0) icomx=3
!      if(icomy.gt.0) icomy=3
!      if(icomz.gt.0) icomz=3
!Set boolean flags for communication in x,y and z directions respectively.
      readsenddataxL=(icomx==1.AND.leftedge==1).OR.                 &
                     (icomx==2.AND.leftedge/=1).OR.                 &
                      icomx==3
      readsenddataxR=(icomx/=1).AND.                                &
                    ((icomx==2 .AND.rightedge/=1).OR.               &
                      icomx==3)
      readsenddataZG=(icomz==1.AND.gndedge==1).OR.                  &
                     (icomz==2.AND.gndedge/=1).OR.                  &
                      icomz==3
      readsenddataZS=(icomz/=1).AND.                                &
                    ((icomz==2 .AND.skyedge/=1).OR.                 &
                      icomz==3)
      recvwrtedataxL=(icomx/=1).AND.                                &
                    ((icomx==2 .AND.leftedge/=1).OR.                &
                      icomx==3)
      recvwrtedataxR=(icomx==1.AND.rightedge==1).OR.                &
                     (icomx==2.AND.rightedge/=1).OR.                &
                      icomx==3
      recvwrtedataZG=(icomz/=1).AND.                                &
                    ((icomz==2 .AND.gndedge/=1).OR.                 &
                      icomz==3)
      recvwrtedataZS=(icomz==1.AND.skyedge==1).OR.                  &
                     (icomz==2.AND.skyedge/=1).OR.                  &
                      icomz==3
      readsenddatayB=(icomy==1.AND.botedge==1).OR.                  &
                     (icomy==2.AND.botedge/=1).OR.                  &
                      icomy==3
      readsenddatayT=(icomy/=1).AND.                                &
                    ((icomy==2 .AND.(topedge/=1)).OR.               &
                      icomy==3)
      recvwrtedatayB=(icomy/=1).AND.                                &
                    ((icomy==2 .AND.(botedge/=1)).OR.               &
                      icomy==3)
      recvwrtedatayT=(icomy==1.AND.topedge==1).OR.                  &
                     (icomy==2.AND.topedge/=1).OR.                  &
                      icomy==3
!Set boolean flags for initiating (sending) or finishing (receiving) message
!First flags to send message depending on the position on the proc grid and icomxyz flag
!For icomxyz=1 only the left,bottom and ground procs send messages
      commx=icomx>0
      commy=icomy>0
      commz=icomz>0

      IF(commy) THEN !icomy > 0 so we are communication in y direction`
!If starting communication
         IF(istart==1) THEN
            IF(readsenddatayB) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                        tmpysnd1(icnt,inovar)=a(i,      j,k)
                     END DO
                  END DO
               END DO
            ENDIF !bottom buffer y direction
            IF(readsenddatayT) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare top data segments for processor on the bottom:tmpysnd2
                        icnt=icnt+1
                        tmpysnd2(icnt,inovar)=a(i,mulim+j,k)
                     END DO
                  END DO
               END DO
            ENDIF ! top buffer y direction
         ENDIF !istart
!
!     exchange data now
!
         IF(nprocy>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN

               IF(istart==1) THEN
#ifdef DEBUGMPI
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
                  IF(readsenddatayB) THEN                                                 
                   CALL MPI_ISEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statsy1(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to snd ',nhl,'bytes to  ', peS, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayT) THEN                                                 
                   CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statry1(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peN, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayT) THEN                                                
                   CALL MPI_ISEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statsy2(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to snd ',nhl,'bytes to  ', peN, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN                                                
                   CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statry2(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peS, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
               ENDIF !istart
#ifdef DEBUGMPI
!               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
               IF(ifinsh==1) THEN
                  IF(readsenddatayB) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'waits to send',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                        CALL MPI_WAIT(statsy1(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'has sent',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayT) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'waits to send',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                         CALL MPI_WAIT(statsy2(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'has sent',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayT) THEN 
#ifdef DEBUGMPI
                    PRINT *,mype, 'expects ',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                        CALL MPI_WAIT(statry1(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'received',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'expects ',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                         CALL MPI_WAIT(statry2(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'received',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
               ENDIF !ifinsh
#ifdef DEBUGMPI
                CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedatayT) THEN 
                   CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statry1(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peN, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN  
                   CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statry2(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'prepared to rcv ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayB)  THEN
#ifdef DEBUGMPI
!                   PRINT *,mype,' sending ',nhl,'bytes to ', peS 
!                   CALL flush(6)
#endif
                   CALL  MPI_SEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,ierr)
#ifdef DEBUGMPI
!                  PRINT *,mype,' has sent ',nhl,'bytes to ', peS,' with ', ierr 
!                   CALL flush(6)
#endif
                   !PRINT *,mype,' sending err',ierr
                   !CALL flush(6)
                  ENDIF
                  IF(readsenddatayT) THEN 
#ifdef DEBUGMPI
!                   PRINT *,mype,' sending ',nhl,'bytes to ', peN 
!                   CALL flush(6)
#endif
                   CALL  MPI_SEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' has sent ',nhl,'bytes to ', peN,' with ', ierr 
!                   CALL flush(6)
#endif
                   !PRINT *,mype,' sending err',ierr
                   !CALL flush(6)
                  ENDIF
               ENDIF !istart
               IF(ifinsh==1) THEN
!                  print *,mype,'Before y wait'
                  IF(recvwrtedatayT) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry1(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry2(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
!                  print *,mype,'After y wait'
               ENDIF !ifinsh
            ENDIF
            CALL ttend(30)
            timecounterdone=.TRUE.
#endif /*PUREMPI*/
         ELSE ! nprocy=1
            IF(recvwrtedatayB) THEN
               DO i=1,nhl
                  tmpyrcv2(i,inovar)=tmpysnd2(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedatayT) THEN
               DO i=1,nhl
                  tmpyrcv1(i,inovar)=tmpysnd1(i,inovar)
               ENDDO
            ENDIF
         ENDIF ! nprocy=1
         IF(ifinsh==1) THEN
            IF(recvwrtedatayB) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        a(i,mllim+j,k)=tmpyrcv2(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatayB
            IF(recvwrtedatayT) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        a(i,   m2+j,k)=tmpyrcv1(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatayT
         ENDIF !ifinsh
      ENDIF !commy >0 Communication in y direction

      IF(commx) THEN ! commx >0 Communication in x direction
         IF(istart==1) THEN
            IF(readsenddataxL) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1(icnt,inovar)=a(      i,j,k)
                     END DO
                  END DO
               END DO
            ENDIF !readsenddataxL
            IF(readsenddataxR) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2(icnt,inovar)=a(nulim+i,j,k)
                     END DO
                  END DO
               END DO
            ENDIF !readsenddataxR
         ENDIF !istart
         IF(nprocx>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
               IF(istart==1) THEN
                  IF(readsenddataxL)                                                &
                   CALL MPI_ISEND(tmpxsnd1(1,inovar),mhl,DC_TYPE,peW,itagoffset+3+100*inovar    &
                                  ,mpi_cust_comm,statsx1(inovar),ierr)
                  IF(recvwrtedataxR)                                                &
                   CALL MPI_IRECV(tmpxrcv1(1,inovar),mhl,DC_TYPE,peE,itagoffset+3+100*inovar    &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxR)                                                &
                   CALL MPI_ISEND(tmpxsnd2(1,inovar),mhl,DC_TYPE,peE,itagoffset+4+100*inovar    &
                                  ,mpi_cust_comm,statsx2(inovar),ierr)
                  IF(recvwrtedataxL)                                                &
                   CALL MPI_IRECV(tmpxrcv2(1,inovar),mhl,DC_TYPE,peW,itagoffset+4+100*inovar    &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
               ENDIF !istart
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
               IF(ifinsh==1) THEN
                  IF(readsenddataxL) CALL MPI_WAIT(statsx1(inovar),status,ierr)
                  IF(readsenddataxR) CALL MPI_WAIT(statsx2(inovar),status,ierr)
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
               ENDIF !ifinsh
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedataxR)                                               &
                   CALL MPI_IRECV(tmpxrcv1(1,inovar),mhl,DC_TYPE,peE,itagoffset+3+100*inovar   &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxL)                                               &
                   CALL  MPI_SEND(tmpxsnd1(1,inovar),mhl,DC_TYPE,peW,itagoffset+3+100*inovar   &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedataxL)                                               &
                   CALL MPI_IRECV(tmpxrcv2(1,inovar),mhl,DC_TYPE,peW,itagoffset+4+100*inovar   &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
                  IF(readsenddataxR)                                               &
                   CALL  MPI_SEND(tmpxsnd2(1,inovar),mhl,DC_TYPE,peE,itagoffset+4+100*inovar   &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
!                  print *,mype,'Before x wait'
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
!                  print *,mype,'After x wait'

               ENDIF !ifinsh
            ENDIF
            CALL ttend(30,timecounterdone)
            timecounterdone=.TRUE.
#endif /*PURMPI*/
         ELSE !now nprocx.eq.1
            IF(recvwrtedataxR) THEN
               DO i=1,mhl
                  tmpxrcv1(i,inovar)=tmpxsnd1(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedataxL) THEN
               DO i=1,mhl
                  tmpxrcv2(i,inovar)=tmpxsnd2(i,inovar)
               ENDDO
            ENDIF
         ENDIF !nprocx

         IF(ifinsh==1) THEN
            IF(recvwrtedataxL) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        a(nllim+i,j,k)=tmpxrcv2(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedataxL
            IF(recvwrtedataxR) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        a(   n2+i,j,k)=tmpxrcv1(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedataxR
         ENDIF !ifinsh
      ENDIF !commx

      IF(commz) THEN ! commz >0 Communication in z direction

         IF(istart==1) THEN
            IF(readsenddatazG) THEN
               icnt=1
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        tmpzsnd1(icnt,inovar)=a(i,j,      k)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !readsenddatazG
            IF(readsenddatazS) THEN
               icnt=1
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        tmpzsnd2(icnt,inovar)=a(i,j,lulim+k)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !readsenddatazS
         ENDIF !istart

         IF(nprocz>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN

               IF(istart==1) THEN
                  IF(readsenddatazG)                                                &
                   CALL MPI_ISEND(tmpzsnd1(1,inovar),lhl,DC_TYPE,peG,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statsz1(inovar),ierr)
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(tmpzrcv1(1,inovar),lhl,DC_TYPE,peZ,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL MPI_ISEND(tmpzsnd2(1,inovar),lhl,DC_TYPE,peZ,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statsz2(inovar),ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(tmpzrcv2(1,inovar),lhl,DC_TYPE,peG,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(readsenddatazG) CALL MPI_WAIT(statsz1(inovar),status,ierr)
                  IF(readsenddatazS) CALL MPI_WAIT(statsz2(inovar),status,ierr)
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh

            ENDIF
            IF (isend==2) THEN
               IF(istart==1) THEN
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(tmpzrcv1(1,inovar),lhl,DC_TYPE,peZ,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazG)                                                &
                   CALL  MPI_SEND(tmpzsnd1(1,inovar),lhl,DC_TYPE,peG,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(tmpzrcv2(1,inovar),lhl,DC_TYPE,peG,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL  MPI_SEND(tmpzsnd2(1,inovar),lhl,DC_TYPE,peZ,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh
            ENDIF
            CALL ttend(30,timecounterdone)
#endif /*PUREMPI*/
         ELSE !now nprocz.eq.1
            IF(recvwrtedatazS) THEN
               DO i=1,lhl
                  tmpzrcv1(i,inovar)=tmpzsnd1(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedatazG) THEN
               DO i=1,lhl
                  tmpzrcv2(i,inovar)=tmpzsnd2(i,inovar)
               ENDDO
            ENDIF
         ENDIF
!
!     store data in main array now
!


         IF(ifinsh==1) THEN
            IF(recvwrtedatazG) THEN
               icnt=1
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        a(i,j,lllim+k)=tmpzrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatazG
            IF(recvwrtedatazS) THEN
               icnt=1
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        a(i,j,   l2+k)=tmpzrcv1(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatazZ
         ENDIF !ifinsh
      ENDIF !commz

   END SUBROUTINE update3dsplit
#endif


   SUBROUTINE update3dsplit &
    (a,nn1,nn2,mm1,mm2,ll1,ll2,ndim1,ndim2,mdim1,mdim2,               &
     ldim1,ldim2,ihg,icomx0,icomy0,icomz0,istart,ifinsh,              &
     isful,inovar,ih) !iopt,ioptz,imode0,inovar)
      USE scratch_datafields, ONLY: &
      tmpxsnd1,tmpxsnd2,tmpxrcv1,tmpxrcv2, & 
      tmpysnd1,tmpysnd2,tmpyrcv1,tmpyrcv2, & 
      tmpzsnd1,tmpzsnd2,tmpzrcv1,tmpzrcv2  
!
!     this SUBROUTINE updates the halo (or ghost) cells surrounding
!     each processors subgrid
!
      INTEGER(KIND=iintegers) :: nn1,nn2,mm1,mm2,ll1,ll2,              &
        ndim1,ndim2,mdim1,mdim2,ldim1,ldim2,ihg,                         &
        icomx0,icomy0,icomz0,icomx,icomy,icomz,istart,ifinsh,isful,inovar
      INTEGER(KIND=iintegers) :: n1,n2,m1,m2,l1,l2,                    &
         nllim,mllim,lllim, nulim,mulim,lulim
      INTEGER nhl,mhl,lhl
      INTEGER(KIND=iintegers) :: ih,i,j,k,ihr
      INTEGER(KIND=iintegers) :: icnt
      INTEGER(KIND=iintegers) :: ilow,iupp,jlow,jupp,klow,kupp,isiz,jsiz,ksiz
      REAL_euwp,INTENT(INOUT) :: &
         a(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
#ifdef PUREMPI
      INTEGER ierr,itagoffset
      INTEGER statsx1(0:msz),statsx2(0:msz)
      INTEGER statrx1(0:msz),statrx2(0:msz)
      INTEGER statsy1(0:msz),statsy2(0:msz)
      INTEGER statry1(0:msz),statry2(0:msz)
      INTEGER statsz1(0:msz),statsz2(0:msz)
      INTEGER statrz1(0:msz),statrz2(0:msz)
      SAVE statsx1,statsx2,statrx1,statrx2
      SAVE statsy1,statsy2,statry1,statry2
      SAVE statsz1,statsz2,statrz1,statrz2
#endif /*PUREMPI*/
      LOGICAL commx,commy,commz
      LOGICAL readsenddataxL,recvwrtedataxL
      LOGICAL readsenddataxR,recvwrtedataxR
      LOGICAL readsenddatayB,recvwrtedatayB
      LOGICAL readsenddatayT,recvwrtedatayT
      LOGICAL readsenddatazG,recvwrtedatazG
      LOGICAL readsenddatazS,recvwrtedatazS
      LOGICAL timecounterdone 
#ifdef DEBUGMPI
      INTEGER iprint,iprr
#endif
#ifdef PUREMPI
      INTEGER status(MPI_STATUS_SIZE)
      timecounterdone = .FALSE.
      itagoffset=1000*isful+icomx0*100+icomy0*300 
#endif /*PUREMPI*/
      IF(inovar.ne.0) STOP 'Wrong inovar'
      IF(ih<3.AND.ihg>=3)                                         &
        STOP 'Trying to update with too large halo extent 3'
      IF(ih<2.AND.ihg>=2)                                         &
        STOP 'Trying to update with too large halo extent 2'
      ihr=min(ih,ihg)
      IF(ih==0)                                         &
        STOP 'Calling update with zero halo size. Probably a bug'
!  ierr=cudaDeviceSynchronize()
      n1=nn1
      n2=nn2
      m1=mm1
      m2=mm2
      l1=ll1
      l2=ll2
      nllim=n1-1-ihr
      mllim=m1-1-ihr
      lllim=l1-1-ihr
      nulim=n2  -ihr
      mulim=m2  -ihr
      lulim=l2  -ihr

      nhl=(n2-n1+1)*ihr*(l2-l1+1)
      mhl=ihr*(m2-m1+1+2*ihr*isful)*(l2-l1+1)
      lhl=(n2-n1+1+2*ihr*isful)*(m2-m1+1+2*ihr*isful)*ihr
#ifdef DEBUGMPI
      iprint=0
!      IF(mype == 0)     PRINT 98,nprocx,nprocy,nprocz,nhl,mhl,lhl
      CALL flush(6)
         CALL MPI_BARRIER(mpi_cust_comm,ierr)
      IF (iprint==1) THEN
       CALL MPI_Comm_size(mpi_cust_comm, mysize, ierr)  ! total numbers of PE''s
         DO iprr=0,mysize
            IF (mype==iprr) THEN
               PRINT 99,mype,rightedge,leftedge, &
                         botedge,topedge, &
                         npos,mpos,lpos,                          &
                         peW,peE,peN,peS,n2-n1+1,m2-m1+1,l2-l1+1
               CALL flush(6)
            ENDIF
         CALL MPI_BARRIER(mpi_cust_comm,ierr)
         ENDDO
98       FORMAT('nprocx = ',i3,'  nprocy = ',i3,'  nprocz = ',i3,' nhl =',i3,' mhl = ',i3,' lhl'/, &
          ' mype R   L   B   T   N   M   L    pW    pE    pN    pS  np  mp lp')
99       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',                               &
                i3,' ',i3,' ',i3,' ',i3,' ', &
                i5,' ',i5,' ',i5,' ',i5,' ',i5,' ',i5,' ',i5,' ')
      ENDIF
#endif

!     there are four border segments that must be sent and received.
!     first, we copy from the large array to the tmp arrays for
!     sending. next, we send our four border segments to surrounding
!     processors. then, we receive the four segments from the
!     surrounding processors and we copy that data to the big
!     array.

!Description of the communication choice logic:
!------------------------------------------------------------------
! Generally, we consider three main types of communication mainly for nonperiodic b.c.
! (For periodic b.c. we usually need to communicate everything, so there is less choice)
! 1. Communicate only between border outermost processor and update only
!    the rightmost (topmost, skymost) processor - ideal for cyclicity enforcement -icomx,y,z=1
! 2. Communication inside the domain for nonperiodic b.c., but do not
!      communicate by default between border outermost processors - icomx,y,z=2
! 3. Communication everywhere (as traditionally done in previous versions of EULAG) icomx,y,z=3
!------------------------------------------------------------------
      icomx=icomx0
      icomy=icomy0
      icomz=icomz0
      IF(ibcx==1.AND.icomx==2) icomx=3
      IF(ibcy==1.AND.icomy==2) icomy=3
      IF(ibcz==1.AND.icomz==2) icomz=3
      IF(isphere==1.AND.ipoles == 2.AND.icomx>0) icomx=3
      IF(isphere==1.AND.ipoles == 2.AND.icomy>0) icomy=3
!      if(icomx.gt.0) icomx=3
!      if(icomy.gt.0) icomy=3
!      if(icomz.gt.0) icomz=3
!Set boolean flags for communication in x,y and z directions respectively.
      readsenddataxL=(icomx==1.AND.leftedge==1).OR.                 &
                     (icomx==2.AND.leftedge/=1).OR.                 &
                      icomx==3
      readsenddataxR=(icomx/=1).AND.                                &
                    ((icomx==2 .AND.rightedge/=1).OR.               &
                      icomx==3)
      readsenddataZG=(icomz==1.AND.gndedge==1).OR.                  &
                     (icomz==2.AND.gndedge/=1).OR.                  &
                      icomz==3
      readsenddataZS=(icomz/=1).AND.                                &
                    ((icomz==2 .AND.skyedge/=1).OR.                 &
                      icomz==3)
      recvwrtedataxL=(icomx/=1).AND.                                &
                    ((icomx==2 .AND.leftedge/=1).OR.                &
                      icomx==3)
      recvwrtedataxR=(icomx==1.AND.rightedge==1).OR.                &
                     (icomx==2.AND.rightedge/=1).OR.                &
                      icomx==3
      recvwrtedataZG=(icomz/=1).AND.                                &
                    ((icomz==2 .AND.gndedge/=1).OR.                 &
                      icomz==3)
      recvwrtedataZS=(icomz==1.AND.skyedge==1).OR.                  &
                     (icomz==2.AND.skyedge/=1).OR.                  &
                      icomz==3
      readsenddatayB=(icomy==1.AND.botedge==1).OR.                  &
                     (icomy==2.AND.botedge/=1).OR.                  &
                      icomy==3
      readsenddatayT=(icomy/=1).AND.                                &
                    ((icomy==2 .AND.(topedge/=1)).OR.               &
                      icomy==3)
      recvwrtedatayB=(icomy/=1).AND.                                &
                    ((icomy==2 .AND.(botedge/=1)).OR.               &
                      icomy==3)
      recvwrtedatayT=(icomy==1.AND.topedge==1).OR.                  &
                     (icomy==2.AND.topedge/=1).OR.                  &
                      icomy==3
!Set boolean flags for initiating (sending) or finishing (receiving) message
!First flags to send message depending on the position on the proc grid and icomxyz flag
!For icomxyz=1 only the left,bottom and ground procs send messages
      commx=icomx>0
      commy=icomy>0
      commz=icomz>0
      IF(commy) THEN !icomy > 0 so we are communication in y direction`
!If starting communication
         IF(istart==1) THEN
            ilow=n1
            iupp=n2
            jlow=1
            jupp=ihr
            klow=l1
            kupp=l2
            isiz=iupp-ilow+1
            jsiz=jupp-jlow+1
            ksiz=kupp-klow+1
            IF(readsenddatayB) THEN
#ifdef CUDACODE
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !1,ihr
                     DO i=ilow,iupp !n1,n2
                        tmpysnd1((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)=a(i,j,      k)
                     END DO
                  END DO
               END DO
#else
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                        tmpysnd1(icnt,inovar)=a(i,      j,k)
                     END DO
                  END DO
               END DO
#endif
            ENDIF !bottom buffer y direction
            IF(readsenddatayT) THEN
#ifdef CUDACODE
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !1,ihr
                     DO i=ilow,iupp !n1,n2
!     prepare top data segments for processor on the bottom:tmpysnd2
                        tmpysnd2((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)=a(i,mulim+j,k)
                     END DO
                  END DO
               END DO
#else
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare top data segments for processor on the bottom:tmpysnd2
                        icnt=icnt+1
                        tmpysnd2(icnt,inovar)=a(i,mulim+j,k)
                     END DO
                  END DO
               END DO
#endif
            ENDIF ! top buffer y direction
         ENDIF !istart
!
!     exchange data now
!
         IF(nprocy>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN

               IF(istart==1) THEN
#ifdef DEBUGMPI
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
                  IF(readsenddatayB) THEN                                                 
                   CALL MPI_ISEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statsy1(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to snd ',nhl,'bytes to  ', peS, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayT) THEN                                                 
                   CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statry1(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peN, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayT) THEN                                                
                   CALL MPI_ISEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statsy2(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to snd ',nhl,'bytes to  ', peN, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN                                                
                   CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statry2(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peS, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
               ENDIF !istart
#ifdef DEBUGMPI
!               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
               IF(ifinsh==1) THEN
                  IF(readsenddatayB) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'waits to send',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                        CALL MPI_WAIT(statsy1(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'has sent',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayT) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'waits to send',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                         CALL MPI_WAIT(statsy2(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'has sent',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayT) THEN 
#ifdef DEBUGMPI
                    PRINT *,mype, 'expects ',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                        CALL MPI_WAIT(statry1(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'received',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'expects ',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                         CALL MPI_WAIT(statry2(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'received',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
               ENDIF !ifinsh
#ifdef DEBUGMPI
                CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedatayT) THEN 
                   CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statry1(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peN, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN  
                   CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statry2(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'prepared to rcv ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayB)  THEN
#ifdef DEBUGMPI
!                   PRINT *,mype,' sending ',nhl,'bytes to ', peS 
!                   CALL flush(6)
#endif
                   CALL  MPI_SEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,ierr)
#ifdef DEBUGMPI
!                  PRINT *,mype,' has sent ',nhl,'bytes to ', peS,' with ', ierr 
!                   CALL flush(6)
#endif
                   !PRINT *,mype,' sending err',ierr
                   !CALL flush(6)
                  ENDIF
                  IF(readsenddatayT) THEN 
#ifdef DEBUGMPI
!                   PRINT *,mype,' sending ',nhl,'bytes to ', peN 
!                   CALL flush(6)
#endif
                   CALL  MPI_SEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' has sent ',nhl,'bytes to ', peN,' with ', ierr 
!                   CALL flush(6)
#endif
                   !PRINT *,mype,' sending err',ierr
                   !CALL flush(6)
                  ENDIF
               ENDIF !istart
               IF(ifinsh==1) THEN
!                  print *,mype,'Before y wait'
                  IF(recvwrtedatayT) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry1(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry2(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
!                  print *,mype,'After y wait'
               ENDIF !ifinsh
            ENDIF
            CALL ttend(30)
            timecounterdone=.TRUE.
#endif /*PUREMPI*/
         ELSE ! nprocy=1
            IF(recvwrtedatayB) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,nhl
                  tmpyrcv2(i,inovar)=tmpysnd2(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedatayT) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,nhl
                  tmpyrcv1(i,inovar)=tmpysnd1(i,inovar)
               ENDDO
            ENDIF
         ENDIF ! nprocy=1
         IF(ifinsh==1) THEN
               ilow=n1
               iupp=n2
               jlow=1
               jupp=ihr
               klow=l1
               kupp=l2
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1
            IF(recvwrtedatayB) THEN
#ifdef CUDACODE
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !1,ihr
                     DO i=ilow,iupp !n1,n2
                        a(i,mllim+j,k)=tmpyrcv2((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)
                     END DO
                  END DO
               END DO
#else
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        a(i,mllim+j,k)=tmpyrcv2(icnt,inovar)
                     END DO
                  END DO
               END DO
#endif
            ENDIF !recvwrtedatayB
            IF(recvwrtedatayT) THEN
#ifdef CUDACODE
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !1,ihr
                     DO i=ilow,iupp !n1,n2
                        a(i,   m2+j,k)=tmpyrcv1((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)
                     END DO
                  END DO
               END DO
#else
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        a(i,   m2+j,k)=tmpyrcv1(icnt,inovar)
                     END DO
                  END DO
               END DO
#endif
            ENDIF !recvwrtedatayT
         ENDIF !ifinsh
      ENDIF !commy >0 Communication in y direction

      IF(commx) THEN ! commx >0 Communication in x direction
         IF(istart==1) THEN
               ilow=1
               iupp=ihr
               jlow=m1-ihr*isful
               jupp=m2+ihr*isful
               klow=l1
               kupp=l2
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1
            IF(readsenddataxL) THEN
#ifdef CUDACODE
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !1,ihr
!     prepare top data segments for processor on the bottom:tmpysnd2
                        tmpxsnd1((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)=a(i,j,k)
                     END DO
                  END DO
               END DO
#else

               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1(icnt,inovar)=a(      i,j,k)
                     END DO
                  END DO
               END DO
#endif
            ENDIF !readsenddataxL
            IF(readsenddataxR) THEN
#ifdef CUDACODE
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !1,ihr
!     prepare top data segments for processor on the bottom:tmpysnd2
                        tmpxsnd2((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)=a(nulim+i,j,k)
                     END DO
                  END DO
               END DO
#else
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2(icnt,inovar)=a(nulim+i,j,k)
                     END DO
                  END DO
               END DO
#endif
            ENDIF !readsenddataxR
         ENDIF !istart
         IF(nprocx>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
               IF(istart==1) THEN
                  IF(readsenddataxL)                                                &
                   CALL MPI_ISEND(tmpxsnd1(1,inovar),mhl,DC_TYPE,peW,itagoffset+3+100*inovar    &
                                  ,mpi_cust_comm,statsx1(inovar),ierr)
                  IF(recvwrtedataxR)                                                &
                   CALL MPI_IRECV(tmpxrcv1(1,inovar),mhl,DC_TYPE,peE,itagoffset+3+100*inovar    &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxR)                                                &
                   CALL MPI_ISEND(tmpxsnd2(1,inovar),mhl,DC_TYPE,peE,itagoffset+4+100*inovar    &
                                  ,mpi_cust_comm,statsx2(inovar),ierr)
                  IF(recvwrtedataxL)                                                &
                   CALL MPI_IRECV(tmpxrcv2(1,inovar),mhl,DC_TYPE,peW,itagoffset+4+100*inovar    &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
               ENDIF !istart
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
               IF(ifinsh==1) THEN
                  IF(readsenddataxL) CALL MPI_WAIT(statsx1(inovar),status,ierr)
                  IF(readsenddataxR) CALL MPI_WAIT(statsx2(inovar),status,ierr)
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
               ENDIF !ifinsh
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedataxR)                                               &
                   CALL MPI_IRECV(tmpxrcv1(1,inovar),mhl,DC_TYPE,peE,itagoffset+3+100*inovar   &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxL)                                               &
                   CALL  MPI_SEND(tmpxsnd1(1,inovar),mhl,DC_TYPE,peW,itagoffset+3+100*inovar   &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedataxL)                                               &
                   CALL MPI_IRECV(tmpxrcv2(1,inovar),mhl,DC_TYPE,peW,itagoffset+4+100*inovar   &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
                  IF(readsenddataxR)                                               &
                   CALL  MPI_SEND(tmpxsnd2(1,inovar),mhl,DC_TYPE,peE,itagoffset+4+100*inovar   &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
!                  print *,mype,'Before x wait'
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
!                  print *,mype,'After x wait'

               ENDIF !ifinsh
            ENDIF
            CALL ttend(30,timecounterdone)
            timecounterdone=.TRUE.
#endif /*PURMPI*/
         ELSE !now nprocx.eq.1
            IF(recvwrtedataxR) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,mhl
                  tmpxrcv1(i,inovar)=tmpxsnd1(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedataxL) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,mhl
                  tmpxrcv2(i,inovar)=tmpxsnd2(i,inovar)
               ENDDO
            ENDIF
         ENDIF !nprocx

         IF(ifinsh==1) THEN
               ilow=1
               iupp=ihr
               jlow=m1-ihr*isful
               jupp=m2+ihr*isful
               klow=l1
               kupp=l2
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1

            IF(recvwrtedataxL) THEN
#ifdef CUDACODE
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !1,ihr
!     prepare top data segments for processor on the bottom:tmpysnd2
                       a(nllim+i,j,k)= tmpxrcv2((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)
                     END DO
                  END DO
               END DO
#else
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        a(nllim+i,j,k)=tmpxrcv2(icnt,inovar)
                     END DO
                  END DO
               END DO
#endif
            ENDIF !recvwrtedataxL
            IF(recvwrtedataxR) THEN
#ifdef CUDACODE
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !1,ihr
!     prepare top data segments for processor on the bottom:tmpysnd2
                        a(   n2+i,j,k)=tmpxrcv1((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)
                     END DO
                  END DO
               END DO

#else
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        a(   n2+i,j,k)=tmpxrcv1(icnt,inovar)
                     END DO
                  END DO
               END DO
#endif
            ENDIF !recvwrtedataxR
         ENDIF !ifinsh
      ENDIF !commx

      IF(commz) THEN ! commz >0 Communication in z direction

         IF(istart==1) THEN
               ilow=n1-ihr*isful
               iupp=n2+ihr*isful
               jlow=m1-ihr*isful
               jupp=m2+ihr*isful
               klow=1
               kupp=ihr
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1
            IF(readsenddatazG) THEN
               icnt=1
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !1,ihr
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !n1-ihr*isful,n2+ihr*isful
                        icnt=(k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1)
                        tmpzsnd1(icnt,inovar)=a(i,j,      k)
                     END DO
                  END DO
               END DO
            ENDIF !readsenddatazG
            IF(readsenddatazS) THEN
               icnt=1
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !1,ihr
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !n1-ihr*isful,n2+ihr*isful
                        icnt=(k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1)
                        tmpzsnd2(icnt,inovar)=a(i,j,lulim+k)
                     END DO
                  END DO
               END DO
            ENDIF !readsenddatazS
         ENDIF !istart

         IF(nprocz>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN

               IF(istart==1) THEN
                  IF(readsenddatazG)                                                &
                   CALL MPI_ISEND(tmpzsnd1(1,inovar),lhl,DC_TYPE,peG,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statsz1(inovar),ierr)
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(tmpzrcv1(1,inovar),lhl,DC_TYPE,peZ,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL MPI_ISEND(tmpzsnd2(1,inovar),lhl,DC_TYPE,peZ,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statsz2(inovar),ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(tmpzrcv2(1,inovar),lhl,DC_TYPE,peG,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(readsenddatazG) CALL MPI_WAIT(statsz1(inovar),status,ierr)
                  IF(readsenddatazS) CALL MPI_WAIT(statsz2(inovar),status,ierr)
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh

            ENDIF
            IF (isend==2) THEN
               IF(istart==1) THEN
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(tmpzrcv1(1,inovar),lhl,DC_TYPE,peZ,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazG)                                                &
                   CALL  MPI_SEND(tmpzsnd1(1,inovar),lhl,DC_TYPE,peG,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(tmpzrcv2(1,inovar),lhl,DC_TYPE,peG,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL  MPI_SEND(tmpzsnd2(1,inovar),lhl,DC_TYPE,peZ,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh
            ENDIF
            CALL ttend(30,timecounterdone)
#endif /*PUREMPI*/
         ELSE !now nprocz.eq.1
            IF(recvwrtedatazS) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,lhl
                  tmpzrcv1(i,inovar)=tmpzsnd1(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedatazG) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,lhl
                  tmpzrcv2(i,inovar)=tmpzsnd2(i,inovar)
               ENDDO
            ENDIF
         ENDIF
!
!     store data in main array now
!


         IF(ifinsh==1) THEN
               ilow=n1-ihr*isful
               iupp=n2+ihr*isful
               jlow=m1-ihr*isful
               jupp=m2+ihr*isful
               klow=1
               kupp=ihr
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1
            IF(recvwrtedatazG) THEN
               icnt=1
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !1,ihr
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !n1-ihr*isful,n2+ihr*isful
!              DO k=1,ihr
!                 DO j=m1-ihr*isful,m2+ihr*isful
!                    DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=(k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1)
                        a(i,j,lllim+k)=tmpzrcv2(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatazG
            IF(recvwrtedatazS) THEN
               icnt=1
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !1,ihr
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !n1-ihr*isful,n2+ihr*isful
!              DO k=1,ihr
!                 DO j=m1-ihr*isful,m2+ihr*isful
!                    DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=(k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1)
                        a(i,j,   l2+k)=tmpzrcv1(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatazZ
         ENDIF !ifinsh
      ENDIF !commz

!  ierr=cudaDeviceSynchronize()
   END SUBROUTINE update3dsplit    
   SUBROUTINE update3dsplit_GPUd &
    (a,nn1,nn2,mm1,mm2,ll1,ll2,ndim1,ndim2,mdim1,mdim2,               &
     ldim1,ldim2,ihg,icomx0,icomy0,icomz0,istart,ifinsh,              &
     isful,inovar,ih) !iopt,ioptz,imode0,inovar)
      USE scratch_datafields, ONLY: &
      tmpxsnd1,tmpxsnd2,tmpxrcv1,tmpxrcv2, & 
      tmpysnd1,tmpysnd2,tmpyrcv1,tmpyrcv2, & 
      tmpzsnd1,tmpzsnd2,tmpzrcv1,tmpzrcv2  
!
!     this SUBROUTINE updates the halo (or ghost) cells surrounding
!     each processors subgrid
!
      INTEGER(KIND=iintegers) :: nn1,nn2,mm1,mm2,ll1,ll2,              &
        ndim1,ndim2,mdim1,mdim2,ldim1,ldim2,ihg,                         &
        icomx0,icomy0,icomz0,icomx,icomy,icomz,istart,ifinsh,isful,inovar
      INTEGER(KIND=iintegers) :: n1,n2,m1,m2,l1,l2,                    &
         nllim,mllim,lllim, nulim,mulim,lulim
      INTEGER nhl,mhl,lhl
      INTEGER(KIND=iintegers) :: ih,i,j,k,ihr
      INTEGER(KIND=iintegers) :: icnt
      INTEGER(KIND=iintegers) :: ilow,iupp,jlow,jupp,klow,kupp,isiz,jsiz,ksiz
      DEV_REAL_euwp,INTENT(INOUT) :: &
         a(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
#ifdef PUREMPI
      INTEGER ierr,itagoffset
      INTEGER statsx1(0:msz),statsx2(0:msz)
      INTEGER statrx1(0:msz),statrx2(0:msz)
      INTEGER statsy1(0:msz),statsy2(0:msz)
      INTEGER statry1(0:msz),statry2(0:msz)
      INTEGER statsz1(0:msz),statsz2(0:msz)
      INTEGER statrz1(0:msz),statrz2(0:msz)
      SAVE statsx1,statsx2,statrx1,statrx2
      SAVE statsy1,statsy2,statry1,statry2
      SAVE statsz1,statsz2,statrz1,statrz2
#endif /*PUREMPI*/
      LOGICAL commx,commy,commz
      LOGICAL readsenddataxL,recvwrtedataxL
      LOGICAL readsenddataxR,recvwrtedataxR
      LOGICAL readsenddatayB,recvwrtedatayB
      LOGICAL readsenddatayT,recvwrtedatayT
      LOGICAL readsenddatazG,recvwrtedatazG
      LOGICAL readsenddatazS,recvwrtedatazS
      LOGICAL timecounterdone 
#ifdef DEBUGMPI
      INTEGER iprint,iprr
#endif
#ifdef PUREMPI
      INTEGER status(MPI_STATUS_SIZE)
      timecounterdone = .FALSE.
      itagoffset=1000*isful+icomx0*100+icomy0*300 
#endif /*PUREMPI*/
      IF(inovar.ne.0) STOP 'Wrong inovar'
      IF(ih<3.AND.ihg>=3)                                         &
        STOP 'Trying to update with too large halo extent 3'
      IF(ih<2.AND.ihg>=2)                                         &
        STOP 'Trying to update with too large halo extent 2'
      ihr=min(ih,ihg)
      IF(ih==0)                                         &
        STOP 'Calling update with zero halo size. Probably a bug'
!  ierr=cudaDeviceSynchronize()
      n1=nn1
      n2=nn2
      m1=mm1
      m2=mm2
      l1=ll1
      l2=ll2
      nllim=n1-1-ihr
      mllim=m1-1-ihr
      lllim=l1-1-ihr
      nulim=n2  -ihr
      mulim=m2  -ihr
      lulim=l2  -ihr

      nhl=(n2-n1+1)*ihr*(l2-l1+1)
      mhl=ihr*(m2-m1+1+2*ihr*isful)*(l2-l1+1)
      lhl=(n2-n1+1+2*ihr*isful)*(m2-m1+1+2*ihr*isful)*ihr
#ifdef DEBUGMPI
      iprint=0
!      IF(mype == 0)     PRINT 98,nprocx,nprocy,nprocz,nhl,mhl,lhl
!     CALL flush(6)
!        CALL MPI_BARRIER(mpi_cust_comm,ierr)
      IF (iprint==1) THEN
       CALL MPI_Comm_size(mpi_cust_comm, mysize, ierr)  ! total numbers of PE''s
         DO iprr=0,mysize
            IF (mype==iprr) THEN
               PRINT 99,mype,rightedge,leftedge, &
                         botedge,topedge, &
                         npos,mpos,lpos,                          &
                         peW,peE,peN,peS,n2-n1+1,m2-m1+1,l2-l1+1
               CALL flush(6)
            ENDIF
         CALL MPI_BARRIER(mpi_cust_comm,ierr)
         ENDDO
98       FORMAT('nprocx = ',i3,'  nprocy = ',i3,'  nprocz = ',i3,' nhl =',i3,' mhl = ',i3,' lhl'/, &
          ' mype R   L   B   T   N   M   L    pW    pE    pN    pS  np  mp lp')
99       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',                               &
                i3,' ',i3,' ',i3,' ',i3,' ', &
                i5,' ',i5,' ',i5,' ',i5,' ',i5,' ',i5,' ',i5,' ')
      ENDIF
#endif

!     there are four border segments that must be sent and received.
!     first, we copy from the large array to the tmp arrays for
!     sending. next, we send our four border segments to surrounding
!     processors. then, we receive the four segments from the
!     surrounding processors and we copy that data to the big
!     array.

!Description of the communication choice logic:
!------------------------------------------------------------------
! Generally, we consider three main types of communication mainly for nonperiodic b.c.
! (For periodic b.c. we usually need to communicate everything, so there is less choice)
! 1. Communicate only between border outermost processor and update only
!    the rightmost (topmost, skymost) processor - ideal for cyclicity enforcement -icomx,y,z=1
! 2. Communication inside the domain for nonperiodic b.c., but do not
!      communicate by default between border outermost processors - icomx,y,z=2
! 3. Communication everywhere (as traditionally done in previous versions of EULAG) icomx,y,z=3
!------------------------------------------------------------------
      icomx=icomx0
      icomy=icomy0
      icomz=icomz0
      IF(ibcx==1.AND.icomx==2) icomx=3
      IF(ibcy==1.AND.icomy==2) icomy=3
      IF(ibcz==1.AND.icomz==2) icomz=3
      IF(isphere==1.AND.ipoles == 2.AND.icomx>0) icomx=3
      IF(isphere==1.AND.ipoles == 2.AND.icomy>0) icomy=3
!      if(icomx.gt.0) icomx=3
!      if(icomy.gt.0) icomy=3
!      if(icomz.gt.0) icomz=3
!Set boolean flags for communication in x,y and z directions respectively.
      readsenddataxL=(icomx==1.AND.leftedge==1).OR.                 &
                     (icomx==2.AND.leftedge/=1).OR.                 &
                      icomx==3
      readsenddataxR=(icomx/=1).AND.                                &
                    ((icomx==2 .AND.rightedge/=1).OR.               &
                      icomx==3)
      readsenddataZG=(icomz==1.AND.gndedge==1).OR.                  &
                     (icomz==2.AND.gndedge/=1).OR.                  &
                      icomz==3
      readsenddataZS=(icomz/=1).AND.                                &
                    ((icomz==2 .AND.skyedge/=1).OR.                 &
                      icomz==3)
      recvwrtedataxL=(icomx/=1).AND.                                &
                    ((icomx==2 .AND.leftedge/=1).OR.                &
                      icomx==3)
      recvwrtedataxR=(icomx==1.AND.rightedge==1).OR.                &
                     (icomx==2.AND.rightedge/=1).OR.                &
                      icomx==3
      recvwrtedataZG=(icomz/=1).AND.                                &
                    ((icomz==2 .AND.gndedge/=1).OR.                 &
                      icomz==3)
      recvwrtedataZS=(icomz==1.AND.skyedge==1).OR.                  &
                     (icomz==2.AND.skyedge/=1).OR.                  &
                      icomz==3
      readsenddatayB=(icomy==1.AND.botedge==1).OR.                  &
                     (icomy==2.AND.botedge/=1).OR.                  &
                      icomy==3
      readsenddatayT=(icomy/=1).AND.                                &
                    ((icomy==2 .AND.(topedge/=1)).OR.               &
                      icomy==3)
      recvwrtedatayB=(icomy/=1).AND.                                &
                    ((icomy==2 .AND.(botedge/=1)).OR.               &
                      icomy==3)
      recvwrtedatayT=(icomy==1.AND.topedge==1).OR.                  &
                     (icomy==2.AND.topedge/=1).OR.                  &
                      icomy==3
!Set boolean flags for initiating (sending) or finishing (receiving) message
!First flags to send message depending on the position on the proc grid and icomxyz flag
!For icomxyz=1 only the left,bottom and ground procs send messages
      commx=icomx>0
      commy=icomy>0
      commz=icomz>0
      IF(commy) THEN !icomy > 0 so we are communication in y direction`
!If starting communication
         IF(istart==1) THEN
               ilow=n1
               iupp=n2
               jlow=1
               jupp=ihr
               klow=l1
               kupp=l2
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1
            IF(readsenddatayB) THEN
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !1,ihr
                     DO i=ilow,iupp !n1,n2
                        tmpysnd1((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)=a(i,j,      k)
                     END DO
                  END DO
               END DO
            ENDIF !bottom buffer y direction
            IF(readsenddatayT) THEN
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !1,ihr
                     DO i=ilow,iupp !n1,n2
!     prepare top data segments for processor on the bottom:tmpysnd2
                        tmpysnd2((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)=a(i,mulim+j,k)
                     END DO
                  END DO
               END DO
            ENDIF ! top buffer y direction
         ENDIF !istart
!
!     exchange data now
!
         IF(nprocy>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN

               IF(istart==1) THEN
#ifdef DEBUGMPI
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
                  IF(readsenddatayB) THEN                                                 
                   CALL MPI_ISEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statsy1(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to snd ',nhl,'bytes to  ', peS, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayT) THEN                                                 
                   CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statry1(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peN, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayT) THEN                                                
                   CALL MPI_ISEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statsy2(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to snd ',nhl,'bytes to  ', peN, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN                                                
                   CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statry2(inovar),ierr)
#ifdef DEBUGMPI
                    PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peS, 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
               ENDIF !istart
#ifdef DEBUGMPI
!               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
               IF(ifinsh==1) THEN
                  IF(readsenddatayB) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'waits to send',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                        CALL MPI_WAIT(statsy1(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'has sent',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayT) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'waits to send',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                         CALL MPI_WAIT(statsy2(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'has sent',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayT) THEN 
#ifdef DEBUGMPI
                    PRINT *,mype, 'expects ',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                        CALL MPI_WAIT(statry1(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'received',nhl,'bytes from', peN , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN
#ifdef DEBUGMPI
                    PRINT *,mype, 'expects ',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                         CALL MPI_WAIT(statry2(inovar),status,ierr)
#ifdef DEBUGMPI
                    PRINT *,mype, 'received',nhl,'bytes from', peS , 'with err',ierr
                    CALL flush(6)
#endif
                  ENDIF
               ENDIF !ifinsh
#ifdef DEBUGMPI
                CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedatayT) THEN 
                   CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statry1(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peN, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN  
                   CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statry2(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'prepared to rcv ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayB)  THEN
#ifdef DEBUGMPI
!                   PRINT *,mype,' sending ',nhl,'bytes to ', peS 
!                   CALL flush(6)
#endif
                   CALL  MPI_SEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,ierr)
#ifdef DEBUGMPI
!                  PRINT *,mype,' has sent ',nhl,'bytes to ', peS,' with ', ierr 
!                   CALL flush(6)
#endif
                   !PRINT *,mype,' sending err',ierr
                   !CALL flush(6)
                  ENDIF
                  IF(readsenddatayT) THEN 
#ifdef DEBUGMPI
!                   PRINT *,mype,' sending ',nhl,'bytes to ', peN 
!                   CALL flush(6)
#endif
                   CALL  MPI_SEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' has sent ',nhl,'bytes to ', peN,' with ', ierr 
!                   CALL flush(6)
#endif
                   !PRINT *,mype,' sending err',ierr
                   !CALL flush(6)
                  ENDIF
               ENDIF !istart
               IF(ifinsh==1) THEN
!                  print *,mype,'Before y wait'
                  IF(recvwrtedatayT) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry1(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry2(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
!                  print *,mype,'After y wait'
               ENDIF !ifinsh
            ENDIF
            CALL ttend(30)
            timecounterdone=.TRUE.
#endif /*PUREMPI*/
         ELSE ! nprocy=1
            IF(recvwrtedatayB) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,nhl
                  tmpyrcv2(i,inovar)=tmpysnd2(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedatayT) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,nhl
                  tmpyrcv1(i,inovar)=tmpysnd1(i,inovar)
               ENDDO
            ENDIF
         ENDIF ! nprocy=1
         IF(ifinsh==1) THEN
               ilow=n1
               iupp=n2
               jlow=1
               jupp=ihr
               klow=l1
               kupp=l2
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1
            IF(recvwrtedatayB) THEN
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !1,ihr
                     DO i=ilow,iupp !n1,n2
                        a(i,mllim+j,k)=tmpyrcv2((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatayB
            IF(recvwrtedatayT) THEN
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !1,ihr
                     DO i=ilow,iupp !n1,n2
                        a(i,   m2+j,k)=tmpyrcv1((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatayT
         ENDIF !ifinsh
      ENDIF !commy >0 Communication in y direction

      IF(commx) THEN ! commx >0 Communication in x direction
         IF(istart==1) THEN
               ilow=1
               iupp=ihr
               jlow=m1-ihr*isful
               jupp=m2+ihr*isful
               klow=l1
               kupp=l2
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1
            IF(readsenddataxL) THEN
               icnt=0
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !1,ihr
                       tmpxsnd1((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)=a(      i,j,k)
                     END DO
                  END DO
               END DO
            ENDIF !readsenddataxL
            IF(readsenddataxR) THEN
               icnt=0
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !1,ihr
                       tmpxsnd2((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)=a(nulim+i,j,k)
                     END DO
                  END DO
               END DO
            ENDIF !readsenddataxR
         ENDIF !istart
         IF(nprocx>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
               IF(istart==1) THEN
                  IF(readsenddataxL)                                                &
                   CALL MPI_ISEND(tmpxsnd1(1,inovar),mhl,DC_TYPE,peW,itagoffset+3+100*inovar    &
                                  ,mpi_cust_comm,statsx1(inovar),ierr)
                  IF(recvwrtedataxR)                                                &
                   CALL MPI_IRECV(tmpxrcv1(1,inovar),mhl,DC_TYPE,peE,itagoffset+3+100*inovar    &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxR)                                                &
                   CALL MPI_ISEND(tmpxsnd2(1,inovar),mhl,DC_TYPE,peE,itagoffset+4+100*inovar    &
                                  ,mpi_cust_comm,statsx2(inovar),ierr)
                  IF(recvwrtedataxL)                                                &
                   CALL MPI_IRECV(tmpxrcv2(1,inovar),mhl,DC_TYPE,peW,itagoffset+4+100*inovar    &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
               ENDIF !istart
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
               IF(ifinsh==1) THEN
                  IF(readsenddataxL) CALL MPI_WAIT(statsx1(inovar),status,ierr)
                  IF(readsenddataxR) CALL MPI_WAIT(statsx2(inovar),status,ierr)
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
               ENDIF !ifinsh
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedataxR)                                               &
                   CALL MPI_IRECV(tmpxrcv1(1,inovar),mhl,DC_TYPE,peE,itagoffset+3+100*inovar   &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxL)                                               &
                   CALL  MPI_SEND(tmpxsnd1(1,inovar),mhl,DC_TYPE,peW,itagoffset+3+100*inovar   &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedataxL)                                               &
                   CALL MPI_IRECV(tmpxrcv2(1,inovar),mhl,DC_TYPE,peW,itagoffset+4+100*inovar   &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
                  IF(readsenddataxR)                                               &
                   CALL  MPI_SEND(tmpxsnd2(1,inovar),mhl,DC_TYPE,peE,itagoffset+4+100*inovar   &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
!                  print *,mype,'Before x wait'
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
!                  print *,mype,'After x wait'

               ENDIF !ifinsh
            ENDIF
            CALL ttend(30,timecounterdone)
            timecounterdone=.TRUE.
#endif /*PURMPI*/
         ELSE !now nprocx.eq.1
            IF(recvwrtedataxR) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,mhl
                  tmpxrcv1(i,inovar)=tmpxsnd1(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedataxL) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,mhl
                  tmpxrcv2(i,inovar)=tmpxsnd2(i,inovar)
               ENDDO
            ENDIF
         ENDIF !nprocx

         IF(ifinsh==1) THEN
               ilow=1
               iupp=ihr
               jlow=m1-ihr*isful
               jupp=m2+ihr*isful
               klow=l1
               kupp=l2
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1
            IF(recvwrtedataxL) THEN
               icnt=0
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !1,ihr
                       a(nllim+i,j,k)=tmpxrcv2((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedataxL
            IF(recvwrtedataxR) THEN
               icnt=0
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !l1,l2
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !1,ihr
                       a(   n2+i,j,k)=tmpxrcv1((k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1),inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedataxR
         ENDIF !ifinsh
      ENDIF !commx

      IF(commz) THEN ! commz >0 Communication in z direction

         IF(istart==1) THEN
               ilow=n1-ihr*isful
               iupp=n2+ihr*isful
               jlow=m1-ihr*isful
               jupp=m2+ihr*isful
               klow=1
               kupp=ihr
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1
            IF(readsenddatazG) THEN
               icnt=1
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !1,ihr
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !n1-ihr*isful,n2+ihr*isful
                        icnt=(k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1)
                        tmpzsnd1(icnt,inovar)=a(i,j,      k)
                     END DO
                  END DO
               END DO
            ENDIF !readsenddatazG
            IF(readsenddatazS) THEN
               icnt=1
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !1,ihr
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !n1-ihr*isful,n2+ihr*isful
                        icnt=(k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1)
                        tmpzsnd2(icnt,inovar)=a(i,j,lulim+k)
                     END DO
                  END DO
               END DO
            ENDIF !readsenddatazS
         ENDIF !istart

         IF(nprocz>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN

               IF(istart==1) THEN
                  IF(readsenddatazG)                                                &
                   CALL MPI_ISEND(tmpzsnd1(1,inovar),lhl,DC_TYPE,peG,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statsz1(inovar),ierr)
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(tmpzrcv1(1,inovar),lhl,DC_TYPE,peZ,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL MPI_ISEND(tmpzsnd2(1,inovar),lhl,DC_TYPE,peZ,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statsz2(inovar),ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(tmpzrcv2(1,inovar),lhl,DC_TYPE,peG,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(readsenddatazG) CALL MPI_WAIT(statsz1(inovar),status,ierr)
                  IF(readsenddatazS) CALL MPI_WAIT(statsz2(inovar),status,ierr)
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh

            ENDIF
            IF (isend==2) THEN
               IF(istart==1) THEN
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(tmpzrcv1(1,inovar),lhl,DC_TYPE,peZ,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazG)                                                &
                   CALL  MPI_SEND(tmpzsnd1(1,inovar),lhl,DC_TYPE,peG,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(tmpzrcv2(1,inovar),lhl,DC_TYPE,peG,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL  MPI_SEND(tmpzsnd2(1,inovar),lhl,DC_TYPE,peZ,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh
            ENDIF
            CALL ttend(30,timecounterdone)
#endif /*PUREMPI*/
         ELSE !now nprocz.eq.1
            IF(recvwrtedatazS) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,lhl
                  tmpzrcv1(i,inovar)=tmpzsnd1(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedatazG) THEN
!$cuf kernel do(1) <<< (*,*), (32,4) >>>
               DO i=1,lhl
                  tmpzrcv2(i,inovar)=tmpzsnd2(i,inovar)
               ENDDO
            ENDIF
         ENDIF
!
!     store data in main array now
!


         IF(ifinsh==1) THEN
               ilow=n1-ihr*isful
               iupp=n2+ihr*isful
               jlow=m1-ihr*isful
               jupp=m2+ihr*isful
               klow=1
               kupp=ihr
               isiz=iupp-ilow+1
               jsiz=jupp-jlow+1
               ksiz=kupp-klow+1
            IF(recvwrtedatazG) THEN
               icnt=1
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !1,ihr
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !n1-ihr*isful,n2+ihr*isful
!              DO k=1,ihr
!                 DO j=m1-ihr*isful,m2+ihr*isful
!                    DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=(k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1)
                        a(i,j,lllim+k)=tmpzrcv2(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatazG
            IF(recvwrtedatazS) THEN
               icnt=1
!$cuf kernel do(3) <<< (*,*), (32,4) >>>
               DO k=klow,kupp !1,ihr
                  DO j=jlow,jupp !m1-ihr*isful,m2+ihr*isful
                     DO i=ilow,iupp !n1-ihr*isful,n2+ihr*isful
!              DO k=1,ihr
!                 DO j=m1-ihr*isful,m2+ihr*isful
!                    DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=(k-klow)*isiz*jsiz+(j-jlow)*isiz+(i-ilow+1)
                        a(i,j,   l2+k)=tmpzrcv1(icnt,inovar)
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatazZ
         ENDIF !ifinsh
      ENDIF !commz

!  ierr=cudaDeviceSynchronize()
   END SUBROUTINE update3dsplit_GPUd    

   SUBROUTINE update3dsplit_multi                                     &
    (nn1,nn2,mm1,mm2,ll1,ll2,ndim1,ndim2,mdim1,mdim2,               &
     ldim1,ldim2,ihg,icomx0,icomy0,icomz0,istart,ifinsh,              &
     isful,inovar,ih, &
     a,b,c,d,e,f,g,h ) !iopt,ioptz,imode0,inovar)
      USE scratch_datafields, ONLY: &
      tmpxsnd1m,tmpxsnd2m,tmpxrcv1m,tmpxrcv2m, & 
      tmpysnd1m,tmpysnd2m,tmpyrcv1m,tmpyrcv2m, & 
      tmpzsnd1m,tmpzsnd2m,tmpzrcv1m,tmpzrcv2m  
!
!     this SUBROUTINE updates the halo (or ghost) cells surrounding
!     each processors subgrid
!
      INTEGER(KIND=iintegers) :: nn1,nn2,mm1,mm2,ll1,ll2,              &
        ndim1,ndim2,mdim1,mdim2,ldim1,ldim2,ihg,                         &
        icomx0,icomy0,icomz0,icomx,icomy,icomz,istart,ifinsh,isful,inovar
      INTEGER(KIND=iintegers) :: n1,n2,m1,m2,l1,l2,                    &
         nllim,mllim,lllim, nulim,mulim,lulim
      INTEGER nhl,mhl,lhl
      INTEGER(KIND=iintegers) :: ih,i,j,k,ihr
      INTEGER(KIND=iintegers) :: icnt,vcnt
      REAL_euwp,INTENT(INOUT) :: &
         a(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      REAL_euwp,INTENT(INOUT),OPTIONAL :: &
         b(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      REAL_euwp,INTENT(INOUT),OPTIONAL :: &
         c(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      REAL_euwp,INTENT(INOUT),OPTIONAL :: &
         d(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      REAL_euwp,INTENT(INOUT),OPTIONAL :: &
         e(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      REAL_euwp,INTENT(INOUT),OPTIONAL :: &
         f(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      REAL_euwp,INTENT(INOUT),OPTIONAL :: &
         g(ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      REAL_euwp,INTENT(INOUT),OPTIONAL :: &
         h(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
#ifdef PUREMPI
      INTEGER ierr,itagoffset
      INTEGER statsx1(0:msz),statsx2(0:msz)
      INTEGER statrx1(0:msz),statrx2(0:msz)
      INTEGER statsy1(0:msz),statsy2(0:msz)
      INTEGER statry1(0:msz),statry2(0:msz)
      INTEGER statsz1(0:msz),statsz2(0:msz)
      INTEGER statrz1(0:msz),statrz2(0:msz)
      SAVE statsx1,statsx2,statrx1,statrx2
      SAVE statsy1,statsy2,statry1,statry2
      SAVE statsz1,statsz2,statrz1,statrz2
#endif /*PUREMPI*/
      LOGICAL commx,commy,commz
      LOGICAL readsenddataxL,recvwrtedataxL
      LOGICAL readsenddataxR,recvwrtedataxR
      LOGICAL readsenddatayB,recvwrtedatayB
      LOGICAL readsenddatayT,recvwrtedatayT
      LOGICAL readsenddatazG,recvwrtedatazG
      LOGICAL readsenddatazS,recvwrtedatazS
      LOGICAL timecounterdone
#ifdef DEBUGMPI
      INTEGER iprint,iprr
#endif
#ifdef PUREMPI
      INTEGER status(MPI_STATUS_SIZE)
      timecounterdone=.FALSE.
      itagoffset=1000*isful+icomx0*100+icomy0*300 
#endif /*PUREMPI*/
      IF(inovar.ne.0) STOP 'Wrong inovar'
      IF(ih<3.AND.ihg>=3)                                         &
        STOP 'Trying to update with too large halo extent 3'
      IF(ih<2.AND.ihg>=2)                                         &
        STOP 'Trying to update with too large halo extent 2'
      ihr=min(ih,ihg)
      IF(ih==0)                                         &
        STOP 'Calling update with zero halo size. Probably a bug'
      n1=nn1
      n2=nn2
      m1=mm1
      m2=mm2
      l1=ll1
      l2=ll2
      nllim=n1-1-ihr
      mllim=m1-1-ihr
      lllim=l1-1-ihr
      nulim=n2  -ihr
      mulim=m2  -ihr
      lulim=l2  -ihr

      nhl=(n2-n1+1)*ihr*(l2-l1+1)
      mhl=ihr*(m2-m1+1+2*ihr*isful)*(l2-l1+1)
      lhl=(n2-n1+1+2*ihr*isful)*(m2-m1+1+2*ihr*isful)*ihr
#ifdef DEBUGMPI
      iprint=0
      IF(mype == 0)     PRINT 98,nprocx,nprocy,nprocz,nhl,mhl,lhl
      CALL flush(6)
         CALL MPI_BARRIER(mpi_cust_comm,ierr)
      IF (iprint==1) THEN
       CALL MPI_Comm_size(mpi_cust_comm, mysize, ierr)  ! total numbers of PE''s
         DO iprr=0,mysize
            IF (mype==iprr) THEN
               PRINT 99,mype,rightedge,leftedge, &
                         botedge,topedge, &
                         npos,mpos,lpos,                          &
                         peW,peE,peN,peS,n2-n1+1,m2-m1+1,l2-l1+1
               CALL flush(6)
            ENDIF
         CALL MPI_BARRIER(mpi_cust_comm,ierr)
         ENDDO
98       FORMAT('nprocx = ',i3,'  nprocy = ',i3,'  nprocz = ',i3,' nhl =',i3,' mhl = ',i3,' lhl'/, &
          ' mype R   L   B   T   N   M   L    pW    pE    pN    pS  np  mp lp')
99       FORMAT(i3,' ',i3,' ',i3,' ',i3,' ',                               &
                i3,' ',i3,' ',i3,' ',i3,' ', &
                i5,' ',i5,' ',i5,' ',i5,' ',i5,' ',i5,' ',i5,' ')
      ENDIF
#endif

!     there are four border segments that must be sent and received.
!     first, we copy from the large array to the tmp arrays for
!     sending. next, we send our four border segments to surrounding
!     processors. then, we receive the four segments from the
!     surrounding processors and we copy that data to the big
!     array.

!Description of the communication choice logic:
!------------------------------------------------------------------
! Generally, we consider three main types of communication mainly for nonperiodic b.c.
! (For periodic b.c. we usually need to communicate everything, so there is less choice)
! 1. Communicate only between border outermost processor and update only
!    the rightmost (topmost, skymost) processor - ideal for cyclicity enforcement -icomx,y,z=1
! 2. Communication inside the domain for nonperiodic b.c., but do not
!      communicate by default between border outermost processors - icomx,y,z=2
! 3. Communication everywhere (as traditionally done in previous versions of EULAG) icomx,y,z=3
!------------------------------------------------------------------
      icomx=icomx0
      icomy=icomy0
      icomz=icomz0
      IF(ibcx==1.AND.icomx==2) icomx=3
      IF(ibcy==1.AND.icomy==2) icomy=3
      IF(ibcz==1.AND.icomz==2) icomz=3
      IF(isphere==1.AND.ipoles == 2.AND.icomx>0) icomx=3
      IF(isphere==1.AND.ipoles == 2.AND.icomy>0) icomy=3
!      if(icomx.gt.0) icomx=3
!      if(icomy.gt.0) icomy=3
!      if(icomz.gt.0) icomz=3
!Set boolean flags for communication in x,y and z directions respectively.
      readsenddataxL=(icomx==1.AND.leftedge==1).OR.                 &
                     (icomx==2.AND.leftedge/=1).OR.                 &
                      icomx==3
      readsenddataxR=(icomx/=1).AND.                                &
                    ((icomx==2 .AND.rightedge/=1).OR.               &
                      icomx==3)
      readsenddataZG=(icomz==1.AND.gndedge==1).OR.                  &
                     (icomz==2.AND.gndedge/=1).OR.                  &
                      icomz==3
      readsenddataZS=(icomz/=1).AND.                                &
                    ((icomz==2 .AND.skyedge/=1).OR.                 &
                      icomz==3)
      recvwrtedataxL=(icomx/=1).AND.                                &
                    ((icomx==2 .AND.leftedge/=1).OR.                &
                      icomx==3)
      recvwrtedataxR=(icomx==1.AND.rightedge==1).OR.                &
                     (icomx==2.AND.rightedge/=1).OR.                &
                      icomx==3
      recvwrtedataZG=(icomz/=1).AND.                                &
                    ((icomz==2 .AND.gndedge/=1).OR.                 &
                      icomz==3)
      recvwrtedataZS=(icomz==1.AND.skyedge==1).OR.                  &
                     (icomz==2.AND.skyedge/=1).OR.                  &
                      icomz==3
      readsenddatayB=(icomy==1.AND.botedge==1).OR.                  &
                     (icomy==2.AND.botedge/=1).OR.                  &
                      icomy==3
      readsenddatayT=(icomy/=1).AND.                                &
                    ((icomy==2 .AND.(topedge/=1)).OR.               &
                      icomy==3)
      recvwrtedatayB=(icomy/=1).AND.                                &
                    ((icomy==2 .AND.(botedge/=1)).OR.               &
                      icomy==3)
      recvwrtedatayT=(icomy==1.AND.topedge==1).OR.                  &
                     (icomy==2.AND.topedge/=1).OR.                  &
                      icomy==3
!Set boolean flags for initiating (sending) or finishing (receiving) message
!First flags to send message depending on the position on the proc grid and icomxyz flag
!For icomxyz=1 only the left,bottom and ground procs send messages
      commx=icomx>0
      commy=icomy>0
      commz=icomz>0
                      vcnt=8
               icnt=0
      IF(commy) THEN !icomy > 0 so we are communication in y direction`
!If starting communication
         IF(istart==1) THEN
            IF(readsenddatayB) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                        tmpysnd1m(icnt)=a(i,      j,k)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                         tmpysnd1m(icnt)=b(i,      j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                        tmpysnd1m(icnt)=c(i,      j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                        tmpysnd1m(icnt)=d(i,      j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                        tmpysnd1m(icnt)=e(i,      j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                        tmpysnd1m(icnt)=f(i,      j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                        tmpysnd1m(icnt)=g(i,      j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!     prepare bottom data segments for processor on the top:tmpysnd1
                        icnt=icnt+1
                        tmpysnd1m(icnt)=h(i,      j,k)
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !bottom buffer y direction
            IF(readsenddatayT) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!
!     prepare top data segments for processor on the bottom:tmpysnd2
!
                       icnt=icnt+1
                       tmpysnd2m(icnt)=a(i,mulim+j,k)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!
!     prepare top data segments for processor on the bottom:tmpysnd2
!
                       icnt=icnt+1
                       tmpysnd2m(icnt)=b(i,mulim+j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!
!     prepare top data segments for processor on the bottom:tmpysnd2
!
                        icnt=icnt+1
                        tmpysnd2m(icnt)=c(i,mulim+j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!
!     prepare top data segments for processor on the bottom:tmpysnd2
!
                        icnt=icnt+1
                        tmpysnd2m(icnt)=d(i,mulim+j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!
!     prepare top data segments for processor on the bottom:tmpysnd2
!
                        icnt=icnt+1
                        tmpysnd2m(icnt)=e(i,mulim+j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!
!     prepare top data segments for processor on the bottom:tmpysnd2
!
                        icnt=icnt+1
                        tmpysnd2m(icnt)=f(i,mulim+j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!
!     prepare top data segments for processor on the bottom:tmpysnd2
!
                        icnt=icnt+1
                        tmpysnd2m(icnt)=g(i,mulim+j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
!
!     prepare top data segments for processor on the bottom:tmpysnd2
!
                        icnt=icnt+1
                        tmpysnd2m(icnt)=h(i,mulim+j,k)
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF ! top buffer y direction
         ENDIF !istart
         
!     exchange data now
!
         IF(nprocy>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN

               IF(istart==1) THEN
#ifdef DEBUGMPI
!               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
                  IF(readsenddatayB) THEN                                                 
                   CALL MPI_ISEND(tmpysnd1m,icnt,DC_TYPE,peS,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statsy1(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' prepared to snd ',nhl,'bytes to  ', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayT) THEN                                                 
                   CALL MPI_IRECV(tmpyrcv1m,icnt,DC_TYPE,peN,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statry1(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peN, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayT) THEN                                                
                   CALL MPI_ISEND(tmpysnd2m,icnt,DC_TYPE,peN,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statsy2(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' prepared to snd ',nhl,'bytes to  ', peN, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN                                                
                   CALL MPI_IRECV(tmpyrcv2m,icnt,DC_TYPE,peS,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statry2(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
               ENDIF !istart
#ifdef DEBUGMPI
!               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
               IF(ifinsh==1) THEN
                  IF(readsenddatayB) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'waits to send',nhl,'bytes from', peS , 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statsy1(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'has sent',nhl,'bytes from', peS , 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayT) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'waits to send',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                         CALL MPI_WAIT(statsy2(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'has sent',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayT) THEN 
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry1(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peS , 'with err',ierr
!                   CALL flush(6)
#endif
                         CALL MPI_WAIT(statry2(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peS , 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
               ENDIF !ifinsh
#ifdef DEBUGMPI
!               CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedatayT) THEN 
                   CALL MPI_IRECV(tmpyrcv1m,icnt,DC_TYPE,peN,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,statry1(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' prepared to rcv ',nhl,'bytes from', peN, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN  
                   CALL MPI_IRECV(tmpyrcv2m,icnt,DC_TYPE,peS,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,statry2(inovar),ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'prepared to rcv ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(readsenddatayB)  THEN
#ifdef DEBUGMPI
!                   PRINT *,mype,' sending ',nhl,'bytes to ', peS 
!                   CALL flush(6)
#endif
                   CALL  MPI_SEND(tmpysnd1m,icnt,DC_TYPE,peS,itagoffset+1+100*inovar    &
                                  ,mpi_cust_comm,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' has sent ',nhl,'bytes to ', peS,' with ', ierr 
!                   CALL flush(6)
#endif
                   !PRINT *,mype,' sending err',ierr
                   !CALL flush(6)
                  ENDIF
                  IF(readsenddatayT) THEN 
#ifdef DEBUGMPI
!                   PRINT *,mype,' sending ',nhl,'bytes to ', peN 
!                   CALL flush(6)
#endif
                   CALL  MPI_SEND(tmpysnd2m,icnt,DC_TYPE,peN,itagoffset+2+100*inovar    &
                                  ,mpi_cust_comm,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype,' has sent ',nhl,'bytes to ', peN,' with ', ierr 
!                   CALL flush(6)
#endif
                   !PRINT *,mype,' sending err',ierr
                   !CALL flush(6)
                  ENDIF
               ENDIF !istart
               IF(ifinsh==1) THEN
!                  print *,mype,'Before y wait'
                  IF(recvwrtedatayT) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry1(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peN , 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
                  IF(recvwrtedatayB) THEN
#ifdef DEBUGMPI
!                   PRINT *,mype, 'expects ',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                        CALL MPI_WAIT(statry2(inovar),status,ierr)
#ifdef DEBUGMPI
!                   PRINT *,mype, 'received',nhl,'bytes from', peS, 'with err',ierr
!                   CALL flush(6)
#endif
                  ENDIF
!                  print *,mype,'After y wait'
               ENDIF !ifinsh
            ENDIF
            CALL ttend(30)
            timecounterdone=.TRUE.
#endif /*PUREMPI*/
         ELSE ! nprocy=1
            IF(recvwrtedatayB) THEN
                  tmpyrcv2m(1:icnt)=tmpysnd2m(1:icnt)
            ENDIF
            IF(recvwrtedatayT) THEN
                  tmpyrcv1m(1:icnt)=tmpysnd1m(1:icnt)
            ENDIF
         ENDIF ! nprocy=1
         IF(ifinsh==1) THEN
            IF(recvwrtedatayB) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        a(i,mllim+j,k)=tmpyrcv2m(icnt)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        b(i,mllim+j,k)=tmpyrcv2m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        c(i,mllim+j,k)=tmpyrcv2m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        d(i,mllim+j,k)=tmpyrcv2m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        e(i,mllim+j,k)=tmpyrcv2m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        f(i,mllim+j,k)=tmpyrcv2m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        g(i,mllim+j,k)=tmpyrcv2m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        h(i,mllim+j,k)=tmpyrcv2m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !recvwrtedatayB
            IF(recvwrtedatayT) THEN
               icnt=0
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        a(i,   m2+j,k)=tmpyrcv1m(icnt)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        b(i,   m2+j,k)=tmpyrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        c(i,   m2+j,k)=tmpyrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        d(i,   m2+j,k)=tmpyrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        e(i,   m2+j,k)=tmpyrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        f(i,   m2+j,k)=tmpyrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        g(i,   m2+j,k)=tmpyrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        icnt=icnt+1
                        h(i,   m2+j,k)=tmpyrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !recvwrtedatayT
         ENDIF !ifinsh
      ENDIF !commy >0 Communication in y direction

      IF(commx) THEN ! commx >0 Communication in x direction
         IF(istart==1) THEN
            IF(readsenddataxL) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1m(icnt)=a(      i,j,k)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1m(icnt)=b(      i,j,k) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1m(icnt)=c(      i,j,k) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1m(icnt)=d(      i,j,k) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1m(icnt)=e(      i,j,k) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1m(icnt)=f(      i,j,k) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1m(icnt)=g(      i,j,k) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd1m(icnt)=h(      i,j,k) 
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !readsenddataxL
            IF(readsenddataxR) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2m(icnt)=a(nulim+i,j,k)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2m(icnt)=b(nulim+i,j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2m(icnt)=c(nulim+i,j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2m(icnt)=d(nulim+i,j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2m(icnt)=e(nulim+i,j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2m(icnt)=f(nulim+i,j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2m(icnt)=g(nulim+i,j,k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        tmpxsnd2m(icnt)=h(nulim+i,j,k)
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !readsenddataxR
         ENDIF !istart
         IF(nprocx>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
               IF(istart==1) THEN
                  IF(readsenddataxL)                                                &
                   CALL MPI_ISEND(tmpxsnd1m,icnt,DC_TYPE,peW,itagoffset+3+100*inovar    &
                                  ,mpi_cust_comm,statsx1(inovar),ierr)
                  IF(recvwrtedataxR)                                                &
                   CALL MPI_IRECV(tmpxrcv1m,icnt,DC_TYPE,peE,itagoffset+3+100*inovar    &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxR)                                                &
                   CALL MPI_ISEND(tmpxsnd2m,icnt,DC_TYPE,peE,itagoffset+4+100*inovar    &
                                  ,mpi_cust_comm,statsx2(inovar),ierr)
                  IF(recvwrtedataxL)                                                &
                   CALL MPI_IRECV(tmpxrcv2m,icnt,DC_TYPE,peW,itagoffset+4+100*inovar    &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
               ENDIF !istart
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
               IF(ifinsh==1) THEN
                  IF(readsenddataxL) CALL MPI_WAIT(statsx1(inovar),status,ierr)
                  IF(readsenddataxR) CALL MPI_WAIT(statsx2(inovar),status,ierr)
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
               ENDIF !ifinsh
               CALL MPI_BARRIER(mpi_cust_comm,ierr)
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedataxR)                                               &
                   CALL MPI_IRECV(tmpxrcv1m,icnt,DC_TYPE,peE,itagoffset+3+100*inovar   &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxL)                                               &
                   CALL  MPI_SEND(tmpxsnd1m,icnt,DC_TYPE,peW,itagoffset+3+100*inovar   &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedataxL)                                               &
                   CALL MPI_IRECV(tmpxrcv2m,icnt,DC_TYPE,peW,itagoffset+4+100*inovar   &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
                  IF(readsenddataxR)                                               &
                   CALL  MPI_SEND(tmpxsnd2m,icnt,DC_TYPE,peE,itagoffset+4+100*inovar   &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
!                  print *,mype,'Before x wait'
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
!                  print *,mype,'After x wait'

               ENDIF !ifinsh
            ENDIF
            CALL ttend(30,timecounterdone)
            timecounterdone=.TRUE.
#endif /*PURMPI*/
         ELSE !now nprocx.eq.1
            IF(recvwrtedataxR) THEN
                  tmpxrcv1m(1:icnt)=tmpxsnd1m(1:icnt)
            ENDIF
            IF(recvwrtedataxL) THEN
                  tmpxrcv2m(1:icnt)=tmpxsnd2m(1:icnt)
            ENDIF
         ENDIF !nprocx

         IF(ifinsh==1) THEN
            IF(recvwrtedataxL) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        a(nllim+i,j,k)=tmpxrcv2m(icnt)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        b(nllim+i,j,k)=tmpxrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        c(nllim+i,j,k)=tmpxrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        d(nllim+i,j,k)=tmpxrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        e(nllim+i,j,k)=tmpxrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        f(nllim+i,j,k)=tmpxrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        g(nllim+i,j,k)=tmpxrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        h(nllim+i,j,k)=tmpxrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !recvwrtedataxL
            IF(recvwrtedataxR) THEN
               icnt=0
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        a(   n2+i,j,k)=tmpxrcv1m(icnt)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        b(   n2+i,j,k)=tmpxrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        c(   n2+i,j,k)=tmpxrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        d(   n2+i,j,k)=tmpxrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        e(   n2+i,j,k)=tmpxrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        f(   n2+i,j,k)=tmpxrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        g(   n2+i,j,k)=tmpxrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        icnt=icnt+1
                        h(   n2+i,j,k)=tmpxrcv1m(icnt) 
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !recvwrtedataxR
         ENDIF !ifinsh
      ENDIF !commx

      IF(commz) THEN ! commz >0 Communication in z direction

         IF(istart==1) THEN
            IF(readsenddatazG) THEN
               icnt=0
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd1m(icnt)=a(i,j,      k)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd1m(icnt)=b(i,j,      k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd1m(icnt)=c(i,j,      k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd1m(icnt)=d(i,j,      k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd1m(icnt)=e(i,j,      k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd1m(icnt)=f(i,j,      k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd1m(icnt)=g(i,j,      k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd1m(icnt)=h(i,j,      k)
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !readsenddatazG
            IF(readsenddatazS) THEN
               icnt=0
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd2m(icnt)=a(i,j,lulim+k)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd2m(icnt)=a(i,j,lulim+k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd2m(icnt)=a(i,j,lulim+k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd2m(icnt)=a(i,j,lulim+k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd2m(icnt)=a(i,j,lulim+k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd2m(icnt)=a(i,j,lulim+k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd2m(icnt)=a(i,j,lulim+k)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        tmpzsnd2m(icnt)=a(i,j,lulim+k)
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !readsenddatazS
         ENDIF !istart

         IF(nprocz>1) THEN
#ifdef PUREMPI
            CALL ttbeg(30)
            IF (isend==3) THEN

               IF(istart==1) THEN
                  IF(readsenddatazG)                                                &
                   CALL MPI_ISEND(tmpzsnd1m,icnt,DC_TYPE,peG,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statsz1(inovar),ierr)
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(tmpzrcv1m,icnt,DC_TYPE,peZ,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL MPI_ISEND(tmpzsnd2m,icnt,DC_TYPE,peZ,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statsz2(inovar),ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(tmpzrcv2m,icnt,DC_TYPE,peG,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(readsenddatazG) CALL MPI_WAIT(statsz1(inovar),status,ierr)
                  IF(readsenddatazS) CALL MPI_WAIT(statsz2(inovar),status,ierr)
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh

            ENDIF
            IF (isend==2) THEN
               IF(istart==1) THEN
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(tmpzrcv1m,icnt,DC_TYPE,peZ,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazG)                                                &
                   CALL  MPI_SEND(tmpzsnd1m,icnt,DC_TYPE,peG,itagoffset+5+100*inovar    &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(tmpzrcv2m,icnt,DC_TYPE,peG,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL  MPI_SEND(tmpzsnd2m,icnt,DC_TYPE,peZ,itagoffset+6+100*inovar    &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh
            ENDIF
            CALL ttend(30,timecounterdone)
#endif /*PUREMPI*/
         ELSE !now nprocz.eq.1
            IF(recvwrtedatazS) THEN
                  tmpzrcv1m(1:icnt)=tmpzsnd1m(1:icnt)
            ENDIF
            IF(recvwrtedatazG) THEN
                  tmpzrcv2m(1:icnt)=tmpzsnd2m(1:icnt)
            ENDIF
         ENDIF
!
!     store data in main array now
!


         IF(ifinsh==1) THEN
            IF(recvwrtedatazG) THEN
               icnt=0
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        a(i,j,lllim+k)=tmpzrcv2m(icnt)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        b(i,j,lllim+k)=tmpzrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        c(i,j,lllim+k)=tmpzrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        d(i,j,lllim+k)=tmpzrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        e(i,j,lllim+k)=tmpzrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        f(i,j,lllim+k)=tmpzrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        g(i,j,lllim+k)=tmpzrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        h(i,j,lllim+k)=tmpzrcv2m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !recvwrtedatazG
            IF(recvwrtedatazS) THEN
               icnt=0
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        a(i,j,   l2+k)=tmpzrcv1m(icnt)
                     END DO
                  END DO
               END DO
               IF(PRESENT(b)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        b(i,j,   l2+k)=tmpzrcv1m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(c)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        c(i,j,   l2+k)=tmpzrcv1m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(d)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        d(i,j,   l2+k)=tmpzrcv1m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(e)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        e(i,j,   l2+k)=tmpzrcv1m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(f)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        f(i,j,   l2+k)=tmpzrcv1m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(g)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        g(i,j,   l2+k)=tmpzrcv1m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
               IF(PRESENT(h)) THEN
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        icnt=icnt+1
                        h(i,j,   l2+k)=tmpzrcv1m(icnt)
                     END DO
                  END DO
               END DO
               ENDIF
            ENDIF !recvwrtedatazZ
         ENDIF !ifinsh
      ENDIF !commz

   END SUBROUTINE update3dsplit_multi

   SUBROUTINE update3dsplit_sphere                                   &
    (a,nn1,nn2,mm1,mm2,ll1,ll2,ndim1,ndim2,mdim1,mdim2,               &
     ldim1,ldim2,ihg,icomx0,icomy0,icomz0,istart,ifinsh,              &
     isful,inovar,ih) !iopt,ioptz,imode0,inovar)
      USE scratch_datafields, ONLY: &
      tmpxsnd1,tmpxsnd2,tmpxrcv1,tmpxrcv2, & 
      tmpysnd1,tmpysnd2,tmpyrcv1,tmpyrcv2, & 
      tmpzsnd1,tmpzsnd2,tmpzrcv1,tmpzrcv2  
!
!     this SUBROUTINE updates the halo (or ghost) cells surrounding
!     each processors subgrid
!
      INTEGER(KIND=iintegers) :: nn1,nn2,mm1,mm2,ll1,ll2,              &
        ndim1,ndim2,mdim1,mdim2,ldim1,ldim2,ihg,                         &
        icomx0,icomy0,icomz0,icomx,icomy,icomz,istart,ifinsh,isful,inovar
      INTEGER(KIND=iintegers) :: n1,n2,m1,m2,l1,l2,                    &
         nllim,mllim,lllim, nulim,mulim,lulim,nhl,mhl,lhl
      INTEGER(KIND=iintegers) :: ih,i,j,k,ihr
      INTEGER(KIND=iintegers) :: icnt
      REAL_euwp,INTENT(INOUT) :: &
         a(ndim1-ih:ndim2+ih,mdim1-ih:mdim2+ih,ldim1-ih:ldim2+ih)
      INTEGER,PARAMETER ::msz=0
#ifdef PUREMPI
      INTEGER ierr
      INTEGER statsx1(0:msz),statsx2(0:msz)
      INTEGER statrx1(0:msz),statrx2(0:msz)
      INTEGER statsy1(0:msz),statsy2(0:msz)
      INTEGER statry1(0:msz),statry2(0:msz)
      INTEGER statsz1(0:msz),statsz2(0:msz)
      INTEGER statrz1(0:msz),statrz2(0:msz)
      SAVE statsx1,statsx2,statrx1,statrx2
      SAVE statsy1,statsy2,statry1,statry2
      SAVE statsz1,statsz2,statrz1,statrz2
#endif /*PUREMPI*/
      LOGICAL commx,commy,commz
      LOGICAL readsenddataxL,recvwrtedataxL
      LOGICAL readsenddataxR,recvwrtedataxR
      LOGICAL readsenddatayB,recvwrtedatayB
      LOGICAL readsenddatayT,recvwrtedatayT
      LOGICAL readsenddatazG,recvwrtedatazG
      LOGICAL readsenddatazS,recvwrtedatazS
#ifdef PUREMPI
      INTEGER status(MPI_STATUS_SIZE)
#endif /*PUREMPI*/
#if(1==0)

      IF(ih<3.AND.ihg==3)                                         &
        STOP 'Trying to update with too large halo extent 3'
      IF(ih<2.AND.ihg==2)                                         &
        STOP 'Trying to update with too large halo extent 2'
      ihr=min(ih,ihg)
!     there are four border segments that must be sent and received.
!     first, we copy from the large array to the tmp arrays for
!     sending. next, we send our four border segments to surrounding
!     processors. then, we receive the four segments from the
!     surrounding processors and we copy that data to the big
!     array.

      n1=nn1
      n2=nn2
      m1=mm1
      m2=mm2
      l1=ll1
      l2=ll2
      nllim=n1-1-ihr
      mllim=m1-1-ihr
      lllim=l1-1-ihr
      nulim=n2  -ihr
      mulim=m2  -ihr
      lulim=l2  -ihr

      nhl=(n2-n1+1)*ihr*(l2-l1+1)
      mhl=ihr*(m2-m1+1+2*ihr*isful)*(l2-l1+1)
      lhl=(n2-n1+1+2*ihr*isful)*(m2-m1+1+2*ihr*isful)*ihr
!Description of the communication choice logic:
!------------------------------------------------------------------
! Generally, we consider three main types of communication mainly for nonperiodic b.c.
! (For periodic b.c. we usually need to communicate everything, so there is less choice)
! 1. Communicate only between border outermost processor and update only
!    the rightmost (topmost, skymost) processor - ideal for cyclicity enforcement -icomx,y,z=1
! 2. Communication inside the domain for nonperiodic b.c., but do not
!      communicate by default between border outermost processors - icomx,y,z=2
! 3. Communication everywhere (as traditionally done in previous versions of EULAG) icomx,y,z=3
!------------------------------------------------------------------
      icomx=icomx0
      icomy=icomy0
      icomz=icomz0
      IF(ibcx==1.AND.icomx==2) icomx=3
      IF(ibcy==1.AND.icomy==2) icomy=3
      IF(ibcz==1.AND.icomz==2) icomz=3
      IF(isphere==1.AND.icomx>0) icomx=3
      IF(isphere==1.AND.icomy>0) icomy=3
!      if(icomx.gt.0) icomx=3
!      if(icomy.gt.0) icomy=3
!      if(icomz.gt.0) icomz=3
!Set boolean flags for communication in x,y and z directions respectively.
      readsenddataxL=(icomx==1.AND.leftedge==1).OR.                 &
                      (icomx==2.AND.leftedge/=1).OR.                 &
                       icomx==3
      readsenddataxR=(icomx/=1).AND.                                  &
                     ((icomx==2 .AND.rightedge/=1).OR.               &
                       icomx==3)
      readsenddataZG=(icomz==1.AND.gndedge==1).OR.                  &
                      (icomz==2.AND.gndedge/=1).OR.                  &
                       icomz==3
      readsenddataZS=(icomz/=1).AND.                                  &
                     ((icomz==2 .AND.skyedge/=1).OR.                 &
                       icomz==3)
      recvwrtedataxL=(icomx/=1).AND.                                  &
                     ((icomx==2 .AND.leftedge/=1).OR.                &
                       icomx==3)
      recvwrtedataxR=(icomx==1.AND.rightedge==1).OR.                &
                      (icomx==2.AND.rightedge/=1).OR.                &
                       icomx==3
      recvwrtedataZG=(icomz/=1).AND.                                  &
                     ((icomz==2 .AND.gndedge/=1).OR.                 &
                       icomz==3)
      recvwrtedataZS=(icomz==1.AND.skyedge==1).OR.                  &
                      (icomz==2.AND.skyedge/=1).OR.                  &
                       icomz==3
      readsenddatayB=(icomy==1.AND.botedge==1).OR.                  &
                      (icomy==2.AND.botedge/=1).OR.                  &
                       icomy==3
      readsenddatayT=(icomy/=1).AND.                                  &
                     ((icomy==2 .AND.(topedge/=1)).OR.               &
                       icomy==3)
      recvwrtedatayB=(icomy/=1).AND.                                  &
                     ((icomy==2 .AND.(botedge/=1)).OR.               &
                       icomy==3)
      recvwrtedatayT=(icomy==1.AND.topedge==1).OR.                  &
                      (icomy==2.AND.topedge/=1).OR.                  &
                       icomy==3
!Set boolean flags for initiating (sending) or finishing (receiving) message
!First flags to send message depending on the position on the proc grid and icomxyz flag
!For icomxyz=1 only the left,bottom and ground procs send messages
      commx=icomx>0
      commy=icomy>0
      commz=icomz>0

      IF(commy) THEN !icomy > 0 so we are communication in y direction`
!If starting communication
         IF(istart==1) THEN
            icnt=1
            DO k=l1,l2
               DO j=1,ihr
                  DO i=n1,n2
                     tmpysnd1(icnt,inovar)=a(i,      j,k)
                     tmpysnd2(icnt,inovar)=a(i,mulim+j,k)
                     icnt=icnt+1
                  END DO
               END DO
            END DO
         ENDIF !istart

!
!     exchange data now
!
         IF(nprocy>1 .OR.(nprocy == 1.AND.nprocx>1)) THEN
#ifdef PUREMPI
            IF (isend==3) THEN

               IF(istart==1) THEN
                  IF(botedge==1) THEN
                     CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,11+100*inovar   &
                                     ,mpi_cust_comm,statry2(inovar),ierr)
                     CALL MPI_ISEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,11+100*inovar   &
                                     ,mpi_cust_comm,statsy1(inovar),ierr)
                    IF(nprocy > 1) THEN
                     CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,1+100*inovar    &
                                     ,mpi_cust_comm,statry1(inovar),ierr)
                     CALL MPI_ISEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,2+100*inovar    &
                                     ,mpi_cust_comm,statsy2(inovar),ierr)
                    ENDIF
                  ENDIF
                  IF (topedge==1) THEN
                     CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,12+100*inovar   &
                                     ,mpi_cust_comm,statry1(inovar),ierr)
                     CALL MPI_ISEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,12+100*inovar   &
                                     ,mpi_cust_comm,statsy2(inovar),ierr)
                    IF(nprocy > 1) THEN
                     CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,2+100*inovar    &
                                     ,mpi_cust_comm,statry2(inovar),ierr)
                     CALL MPI_ISEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,1+100*inovar    &
                                     ,mpi_cust_comm,statsy1(inovar),ierr)
                    ENDIF
                  ENDIF
                  IF(botedge==0.AND.topedge==0) THEN
                     CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,1+100*inovar    &
                                     ,mpi_cust_comm,statry1(inovar),ierr)
                     CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,2+100*inovar    &
                                     ,mpi_cust_comm,statry2(inovar),ierr)
                     CALL MPI_ISEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,1+100*inovar    &
                                     ,mpi_cust_comm,statsy1(inovar),ierr)
                     CALL MPI_ISEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,2+100*inovar    &
                                     ,mpi_cust_comm,statsy2(inovar),ierr)
                  ENDIF
               ENDIF !istart
               IF(ifinsh==1) THEN
                  CALL MPI_WAIT(statsy1(inovar),status,ierr)
                  CALL MPI_WAIT(statsy2(inovar),status,ierr)
                  CALL MPI_WAIT(statry1(inovar),status,ierr)
                  CALL MPI_WAIT(statry2(inovar),status,ierr)
               ENDIF !ifinsh
            ENDIF
            IF (isend==2) THEN
                    IF(nprocy == 1) STOP 'Change for isend=3 for nprocy 1.Blocking MPI not realizable'
               IF(istart==1) THEN
                  IF(botedge==1) THEN
                     CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,1+100*inovar    &
                                     ,mpi_cust_comm,statry1(inovar),ierr)
                     CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,11+100*inovar   &
                                     ,mpi_cust_comm,statry2(inovar),ierr)
                     CALL MPI_SEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,11+100*inovar    &
                                     ,mpi_cust_comm,ierr)
                     CALL MPI_SEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,2+100*inovar     &
                                     ,mpi_cust_comm,ierr)
                  ELSE IF (topedge==1) THEN
                     CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,12+100*inovar   &
                                     ,mpi_cust_comm,statry1(inovar),ierr)
                     CALL MPI_SEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,12+100*inovar    &
                                     ,mpi_cust_comm,ierr)
                     CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,2+100*inovar    &
                                     ,mpi_cust_comm,statry2(inovar),ierr)
                     CALL MPI_SEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,1+100*inovar     &
                                     ,mpi_cust_comm,ierr)
                  ELSE  IF(botedge==0.AND.topedge==0) THEN
                     CALL MPI_IRECV(tmpyrcv1(1,inovar),nhl,DC_TYPE,peN,1+100*inovar    &
                                     ,mpi_cust_comm,statry1(inovar),ierr)
                     CALL MPI_IRECV(tmpyrcv2(1,inovar),nhl,DC_TYPE,peS,2+100*inovar    &
                                     ,mpi_cust_comm,statry2(inovar),ierr)
                     CALL MPI_SEND(tmpysnd1(1,inovar),nhl,DC_TYPE,peS,1+100*inovar     &
                                     ,mpi_cust_comm,ierr)
                     CALL MPI_SEND(tmpysnd2(1,inovar),nhl,DC_TYPE,peN,2+100*inovar     &
                                     ,mpi_cust_comm,ierr)
                  ENDIF
               ENDIF !istart
               IF(ifinsh==1) THEN
                  CALL MPI_WAIT(statry1(inovar),status,ierr)
                  CALL MPI_WAIT(statry2(inovar),status,ierr)
               ENDIF !ifinsh
            ENDIF
#endif /*PUREMPI*/
         ELSE ! nprocy=1
            DO i=1,nhl
               tmpyrcv2(i,inovar)=tmpysnd1(i,inovar)
            ENDDO
            DO i=1,nhl
               tmpyrcv1(i,inovar)=tmpysnd2(i,inovar)
            ENDDO
         ENDIF ! nprocy=1
         IF(ifinsh==1.AND.nprocy>1.AND.nprocx>1) THEN
            icnt=1
            IF(botedge==1) THEN  !-------------------BOTTOM--------
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        a(i,   m2+j,k)=tmpyrcv1(icnt,inovar)
                        a(i,1-j,k)=tmpyrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF
            icnt=1
            IF (topedge==1) THEN  !-------------TOP----------
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        a(i,m2+1+ihr-j,k)=tmpyrcv1(icnt,inovar)  !north pole
                        a(i,    -ihr+j,k)=tmpyrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF
            icnt=1
            IF(topedge==0.AND.botedge==0) THEN  !---------MIDDLE---------
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        a(i,   m2+j,k)=tmpyrcv1(icnt,inovar)
                        a(i,mllim+j,k)=tmpyrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF
         ENDIF !ifinsh
         IF(ifinsh==1.AND.nprocy == 1.AND.nprocx>1) THEN
            icnt=1
!            IF(botedge==1) THEN  !-------------------BOTTOM--------
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        a(i,1-j,k)=tmpyrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
!            ENDIF
            icnt=1
!            IF (topedge==1) THEN  !-------------TOP----------
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        a(i,m2+1+ihr-j,k)=tmpyrcv1(icnt,inovar)  !north pole
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
!            ENDIF
          ENDIF !ifinsh
         IF(ifinsh==1.AND.nprocx==1.AND.nprocy>1) THEN
            icnt=1
            IF(botedge==1) THEN  !-------------------BOTTOM--------
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        a(i,   m2+j,k)=tmpyrcv1(icnt,inovar)
                        a(ip(i),1-j,k)=tmpyrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ELSE IF (topedge==1) THEN  !-------------TOP----------
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        a(ip(i),m2+1+ihr-j,k)=tmpyrcv1(icnt,inovar)  !north pole
                        a(i,    -ihr+j,k)=tmpyrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ELSE IF(topedge==0.AND.botedge==0) THEN  !---------MIDDLE---------
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        a(i,   m2+j,k)=tmpyrcv1(icnt,inovar)
                        a(i,mllim+j,k)=tmpyrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF
         ENDIF !ifinsh
         IF(ifinsh==1.AND.nprocx==1.AND.nprocy==1) THEN
            icnt=1
!            IF(botedge==1) THEN  !-------------------BOTTOM--------
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        a(ip(i),1-j,k)=tmpyrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
!            ELSE IF (topedge==1) THEN  !-------------TOP----------
            icnt=1
               DO k=l1,l2
                  DO j=1,ihr
                     DO i=n1,n2
                        a(ip(i),m2+1+ihr-j,k)=tmpyrcv1(icnt,inovar)  !north pole
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
!            ENDIF
         ENDIF !ifinsh
      ENDIF !commy >0 Communication in y direction

      IF(commx) THEN ! commx >0 Communication in x direction
         IF(istart==1) THEN
            IF(readsenddataxL) THEN
               icnt=1
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        tmpxsnd1(icnt,inovar)=a(      i,j,k)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !readsenddataxL
            IF(readsenddataxR) THEN
               icnt=1
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        tmpxsnd2(icnt,inovar)=a(nulim+i,j,k)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !readsenddataxR
         ENDIF !istart
         IF(nprocx>1) THEN
#ifdef PUREMPI
            IF (isend==3) THEN
               IF(istart==1) THEN
                  IF(readsenddataxL)                                                &
                   CALL MPI_ISEND(tmpxsnd1(1,inovar),mhl,DC_TYPE,peW,3+100*inovar    &
                                  ,mpi_cust_comm,statsx1(inovar),ierr)
                  IF(recvwrtedataxR)                                                &
                   CALL MPI_IRECV(tmpxrcv1(1,inovar),mhl,DC_TYPE,peE,3+100*inovar    &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxR)                                                &
                   CALL MPI_ISEND(tmpxsnd2(1,inovar),mhl,DC_TYPE,peE,4+100*inovar    &
                                  ,mpi_cust_comm,statsx2(inovar),ierr)
                  IF(recvwrtedataxL)                                                &
                   CALL MPI_IRECV(tmpxrcv2(1,inovar),mhl,DC_TYPE,peW,4+100*inovar    &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(readsenddataxL) CALL MPI_WAIT(statsx1(inovar),status,ierr)
                  IF(readsenddataxR) CALL MPI_WAIT(statsx2(inovar),status,ierr)
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
               ENDIF !ifinsh
            ENDIF
            IF (isend==2) THEN

               IF(istart==1) THEN
                  IF(recvwrtedataxR)                                               &
                   CALL MPI_IRECV(tmpxrcv1(1,inovar),mhl,DC_TYPE,peE,3+100*inovar   &
                                  ,mpi_cust_comm,statrx1(inovar),ierr)
                  IF(readsenddataxL)                                               &
                   CALL  MPI_SEND(tmpxsnd1(1,inovar),mhl,DC_TYPE,peW,3+100*inovar   &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedataxL)                                               &
                   CALL MPI_IRECV(tmpxrcv2(1,inovar),mhl,DC_TYPE,peW,4+100*inovar   &
                                  ,mpi_cust_comm,statrx2(inovar),ierr)
                  IF(readsenddataxR)                                               &
                   CALL  MPI_SEND(tmpxsnd2(1,inovar),mhl,DC_TYPE,peE,4+100*inovar   &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(recvwrtedataxR) CALL MPI_WAIT(statrx1(inovar),status,ierr)
                  IF(recvwrtedataxL) CALL MPI_WAIT(statrx2(inovar),status,ierr)
               ENDIF !ifinsh
            ENDIF
#endif /*PURMPI*/
         ELSE !now nprocx.eq.1
            IF(recvwrtedataxR) THEN
               DO i=1,mhl
                  tmpxrcv1(i,inovar)=tmpxsnd1(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedataxL) THEN
               DO i=1,mhl
                  tmpxrcv2(i,inovar)=tmpxsnd2(i,inovar)
               ENDDO
            ENDIF
         ENDIF !nprocx

         IF(ifinsh==1) THEN
            IF(recvwrtedataxL) THEN
               icnt=1
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        a(nllim+i,j,k)=tmpxrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedataxL
            IF(recvwrtedataxR) THEN
               icnt=1
               DO k=l1,l2
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=1,ihr
                        a(   n2+i,j,k)=tmpxrcv1(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedataxR
         ENDIF !ifinsh
      ENDIF !commx

      IF(commz) THEN ! commz >0 Communication in z direction

         IF(istart==1) THEN
            IF(readsenddatazG) THEN
               icnt=1
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        tmpzsnd1(icnt,inovar)=a(i,j,      k)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !readsenddatazG
            IF(readsenddatazS) THEN
               icnt=1
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        tmpzsnd2(icnt,inovar)=a(i,j,lulim+k)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !readsenddatazS
         ENDIF !istart

         IF(nprocz>1) THEN
#ifdef PUREMPI
            IF (isend==3) THEN

               IF(istart==1) THEN
                  IF(readsenddatazG)                                                &
                   CALL MPI_ISEND(tmpzsnd1(1,inovar),lhl,DC_TYPE,peG,5+100*inovar    &
                                  ,mpi_cust_comm,statsz1(inovar),ierr)
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(tmpzrcv1(1,inovar),lhl,DC_TYPE,peZ,5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL MPI_ISEND(tmpzsnd2(1,inovar),lhl,DC_TYPE,peZ,6+100*inovar    &
                                  ,mpi_cust_comm,statsz2(inovar),ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(tmpzrcv2(1,inovar),lhl,DC_TYPE,peG,6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(readsenddatazG) CALL MPI_WAIT(statsz1(inovar),status,ierr)
                  IF(readsenddatazS) CALL MPI_WAIT(statsz2(inovar),status,ierr)
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh

            ENDIF
            IF (isend==2) THEN
               IF(istart==1) THEN
                  IF(recvwrtedatazS)                                                &
                   CALL MPI_IRECV(tmpzrcv1(1,inovar),lhl,DC_TYPE,peZ,5+100*inovar    &
                                  ,mpi_cust_comm,statrz1(inovar),ierr)
                  IF(readsenddatazG)                                                &
                   CALL  MPI_SEND(tmpzsnd1(1,inovar),lhl,DC_TYPE,peG,5+100*inovar    &
                                  ,mpi_cust_comm,ierr)
                  IF(recvwrtedatazG)                                                &
                   CALL MPI_IRECV(tmpzrcv2(1,inovar),lhl,DC_TYPE,peG,6+100*inovar    &
                                  ,mpi_cust_comm,statrz2(inovar),ierr)
                  IF(readsenddatazS)                                                &
                   CALL  MPI_SEND(tmpzsnd2(1,inovar),lhl,DC_TYPE,peZ,6+100*inovar    &
                                  ,mpi_cust_comm,ierr)
               ENDIF !istart
               IF(ifinsh==1) THEN
                  IF(recvwrtedatazS) CALL MPI_WAIT(statrz1(inovar),status,ierr)
                  IF(recvwrtedatazG) CALL MPI_WAIT(statrz2(inovar),status,ierr)
               ENDIF !ifinsh
            ENDIF
#endif /*PUREMPI*/
         ELSE !now nprocz.eq.1
            IF(recvwrtedatazS) THEN
               DO i=1,lhl
                  tmpzrcv1(i,inovar)=tmpzsnd1(i,inovar)
               ENDDO
            ENDIF
            IF(recvwrtedatazG) THEN
               DO i=1,lhl
                  tmpzrcv2(i,inovar)=tmpzsnd2(i,inovar)
               ENDDO
            ENDIF
         ENDIF
!
!     store data in main array now
!


         IF(ifinsh==1) THEN
            IF(recvwrtedatazG) THEN
               icnt=1
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        a(i,j,lllim+k)=tmpzrcv2(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatazG
            IF(recvwrtedatazS) THEN
               icnt=1
               DO k=1,ihr
                  DO j=m1-ihr*isful,m2+ihr*isful
                     DO i=n1-ihr*isful,n2+ihr*isful
                        a(i,j,   l2+k)=tmpzrcv1(icnt,inovar)
                        icnt=icnt+1
                     END DO
                  END DO
               END DO
            ENDIF !recvwrtedatazZ
         ENDIF !ifinsh
      ENDIF !commz
#endif
   END SUBROUTINE update3dsplit_sphere

   SUBROUTINE ttini()
#ifdef CUDACODE
   integer i,ierr
     DO i=1,maxtimers
       ierr = cudaEventCreate(cttstart(i))
       ierr = cudaEventCreate(cttend(i))
     ENDDO
#endif
      clabclr(:)='37'
      clabels(1)='mpdatm_leg'
      clabels(2)='mpdatm_gbc'
      clabels(3)='mpdatm_hbc'
      clabels(4)='mpdatm_spl'
      clabels(5)='mpdatm_wr2'
      clabels(6)='mpdatm_rhospl'
      clabels(7)='mpdatm_wr3'
      clabels(8)='mpdatm_LAM'
      clabels(9)='mpdatm_rho'
      clabels(11)='mpdata_leg'
      clabels(12)='mpdata_gbc'
      clabels(13)='mpdata_hbc'
      clabels(14)='mpdata_spl'
      clabels(15)='upwind_expanded '
      clabels(16)='upwind_gpubc'
      clabels(17)='upwind_halobc'
      clabels(18)='upwinds_two_gbc'
      clabels(19)='upwinds_two_hbc'
      clabels(20)='upwinds_two_leg'
      clabclr(1:20)='32'
      clabels(50)='adif_gau_halobc'
      clabels(51)='adif_gau_gpubc'
      clabels(52)='adif_std_halobc'
      clabels(53)='adif_std_gpubc'
      clabels(54)='twos_gau_halobc'
      clabels(55)='twos_gau_gpubc'
      clabels(56)='twos_std_halobc'
      clabels(57)='twos_std_gpubc'
      clabclr(50:57)='32'
      clabels(21)='update'
      clabels(22)='update3'
      clabels(23)='updatelr'
      clabels(24)='updatebt'
      clabels(25)='updategs'
      clabels(26)='updateflt'
      clabels(27)='update_multi'
      clabels(28)='update3_multi'
      clabels(29)='updatelrbt_multi'
!      clabels(28)='mpdatm_sphere'
!      clabels(29)='upwind_sphere'
      clabels(30)='actualMPIp2peertime'
      clabels(31)='actualMPIreductiontime'
      clabels(32)='globsum'
      clabels(33)='globmax'
      clabels(34)='globmin'
      clabclr(21:34)='31'
      clabels(39)='diagnostic routines'
      clabclr(39)='33'
      clabels(40)='velprd traj0'
      clabels(41)='velprd traj1'
      clabclr(40:41)='32'
      clabels(70)='beforeGCRK-preparations'
      clabels(71)='afterGCRK-update integr'
      clabclr(70:71)='35'
      clabels(75)='wide updates lrwbtwgsw'
      clabclr(75)='31'
      clabels(79)='ADVcmprsunidriver'
      clabels(80)='ADVCEonestepGPUBC'
      clabels(81)='ADVCEonestepHBC'
      clabels(82)='ADVCEtwostepGPUBC'
      clabels(83)='ADVCEtwostepHBC'
      clabels(84)='ADVCElegacy'
      clabels(85)='ADVCEdriverinit'
      clabels(86)='ADVanelonestepGPUBC'
      clabels(87)='ADVanelonestepHBC'
      clabels(88)='ADVaneltwostepGPUBC'
      clabels(89)='ADVaneltwostepHBC'
      clabels(90)='ADVanellegacy'
      clabclr(79:90)='32'
      clabels(91)='GCRKsolverandscndhalf'
      clabclr(91)='35'
      clabels(92)='advec standalone'
      clabclr(92)='32'
      clabels(100)='gcrk'
      clabels(101)='gcrk_init'
      clabels(102)='prforc_legacy'
      clabels(103)='laplc_legacy'
      clabels(104)='coef0'
      clabels(105)='divlpc_legacy'
      clabels(106)='precon'
      clabels(107)='precon_coefs'
      clabels(108)='precon_coefs_cmp_z'
      clabels(109)='precon_bcz_init'
      clabels(110)='rhsdiv_legacy'
      clabels(111)='theimpl'
      clabels(112)='prforc_gpubc'
      clabels(113)='prforc_halobc'
      clabels(114)='laplc_gpubc'
      clabels(115)='laplc_halobc'
      clabels(116)='divlpc_gpubc'
      clabels(117)='divlpc_halobc'
      clabels(118)='rhsdiv_gpubc'
      clabels(119)='rhsdiv_halobc'
      clabels(120)='gcrk initial norms'
      clabels(121)='gcrk final norms'
      clabels(131)='precon_tridg_z'     
      clabels(132)='precon_tridg_y'     
      clabels(133)='precon_tridg_x'     
      clabels(134)='tdmapar_z'     
      clabels(135)='tdmapar_y'     
      clabels(136)='tdmapar_x'     
      clabclr(100:136)='35'
      clabels(200)='diff_fckflxdv'
      clabels(201)='diff_fckflxdv_fx'
      clabels(202)='diff_fckflxdv_fy'
      clabels(203)='diff_fckflxdv_fz'
      clabclr(200:203)='34'
      clabels(300)='eulagdwarf'
      clabels(301)='timeloop'
      clabels(401)='computecourlipsch'
      clabels(402)='computecourants'
      clabclr(401:402)='33'
      clabels(1002)='mpdt_gau_gbc_excl'
      clabels(1003)='mpdt_gau_hbc_excl'
      clabels(1009)='mpdt_gau_rho_excl'
      clabels(1012)='mpdt_std_gbc_excl'
      clabels(1013)='mpdt_std_hbc_excl'
      clabels(1018)='upwinds_tstg_excl'
      clabels(1019)='upwinds_tsth_excl'
      clabels(1030)='mpdt_gau_s_v_excl'
      clabels(1031)='mpdt_gau_s_h_excl'
      clabels(1040)='velprd_excl'
      clabels(1050)='adif_gau_hbc_excl'
      clabels(1051)='adif_gau_gbc_excl'
      clabels(1052)='adif_std_hbc_excl'
      clabels(1053)='adif_std_gbc_excl'
      clabclr(1002:1053)='32'
      clabels(1070)='beforeGCRK_excl'
      clabels(1071)='afterGCRK_excl'
      clabels(1100)='gcrk_exclusive'
      clabels(1102)='prforc_exclusive'
      clabels(1103)='laplc_exclusive'
      clabels(1105)='divlpc_exclusive'
      clabels(1106)='precon_exclusive'
      clabels(1109)='prc_bcz_init_exc'
      clabels(1110)='rhsdiv_exclusive'
      clabels(1120)='gcrk_innorm_excl'
      clabels(1121)='gcrk_fnnorm_excl'
      clabels(1131)='prc_tridg_z_excl'     
      clabels(1132)='prc_tridg_y_excl'     
      clabels(1133)='prc_tridg_x_excl'     
      clabclr(1100:1121)='35'
      clabels(1200)='total_excl'
      clabels(1201)='advection_excl'
      clabels(1202)='solver_excl'
      clabels(1401)='courlipsch_excl'
      clabels(1402)='courants_excl'
      clabclr(1401:1402)='33'

         rtimerb(:)=0.
         rtimer0(:)=0.
         icounts(:)=0
         rtimerbav(:)=0.
         rtimerbmx(:)=0.
         rtimerbmn(:)=0.
   END SUBROUTINE ttini

   ! Calls in subroutine TTBEG: 
   ! => cpu_time (on line <1717>)
   SUBROUTINE ttbeg(icountnr)
      INTEGER       icountnr,ierr
      REAL(KIND=dp) timeloc
      timeloc=0._dp
      ierr=0
#ifdef CUDACODE
!     CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif /*CUDACODE*/
#ifdef SCALTEST 
      if(icountnr.ne.300) return
#endif 
#ifdef DEBUGMPI
      IF(mype.eq.0) print *,'ttbeg ',icountnr,clabels(icountnr)
!      IF(mype.eq.0) print 101,icountnr,clabels(icountnr)
      call flush(6)
      CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
#ifdef TIMERSCPU
      CALL cpu_time(timeloc)
#endif /*TIMERSCPU*/
#ifdef TIMERSMPI
      timeloc= mpi_wtime()
#endif /*TIMERSMPI*/
#ifdef TIMERSCUDA
      ierr = cudaEventRecord(cttstart(icountnr),0) 
      timeloc =0.
#endif/* TIMERSCUDA*/ 
#ifdef TIMERSOMP
      timeloc=omp_get_wtime()
#endif /*TIMERSOMP*/
 101  FORMAT ('ttbeg: ',i3,A)
      rtimer0(icountnr)=timeloc
   END SUBROUTINE ttbeg
   ! Calls in subroutine TTEND: 
   ! => cpu_time (on line <1734>)
   SUBROUTINE ttend(icountnr,ignore_counter)
#ifdef CUDACODE
      INTEGER istat
#endif /*CUDACODE*/
      INTEGER   icountnr,ierr
      LOGICAL,OPTIONAL:: ignore_counter 
      REAL(KIND=dp) timeloc
      REAL timeloc_real
      timeloc=0.
      ierr=0
#ifdef CUDACODE
!     CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif /*CUDACODE*/
#ifdef SCALTEST 
      if(icountnr.ne.300) return
#endif 
#ifdef DEBUGMPI
      CALL MPI_BARRIER(mpi_cust_comm,ierr)
      IF(mype.eq.0) print *,'ttend ',icountnr,clabels(icountnr)
      call flush(6)
      CALL MPI_BARRIER(mpi_cust_comm,ierr)
#endif
#ifdef TIMERSCPU
      CALL cpu_time(timeloc)
#endif /*TIMERSCPU*/
#ifdef TIMERSMPI
      timeloc= mpi_wtime()
#endif /*TIMERSMPI*/
#ifdef TIMERSOMP
      timeloc=omp_get_wtime()
#endif /*TIMERSOMP*/
#ifdef TIMERSCUDA
      ierr =  cudaEventRecord(cttend(icountnr),0) 
      ierr = cudaEventSynchronize(cttend(icountnr))
      ierr = cudaEventElapsedTime(timeloc_real,cttstart(icountnr),cttend(icountnr))
      timeloc=timeloc_real/1000.
#endif/* TIMERSCUDA*/ 
      rtimerb(icountnr)=rtimerb(icountnr)+(timeloc-rtimer0(icountnr))
      icounts(icountnr)=icounts(icountnr)+1
         rtimer0(icountnr)=timeloc
     IF(PRESENT(ignore_counter)) THEN
       IF(ignore_counter) THEN
         icounts(icountnr)=icounts(icountnr)-1
       ENDIF
     ENDIF
   END SUBROUTINE ttend

   ! Calls in subroutine TTPRT: 
   ! => timer_reduce (on line <1750>)
   SUBROUTINE ttprt
      INTEGER i,j
      REAL(KIND=dp) runit,fmax,temp
      CHARACTER(LEN=:), ALLOCATABLE :: colbegin
      CHARACTER(LEN=:), ALLOCATABLE :: colend
      CALL timer_reduce
      IF(mype==0) THEN
!create artificial timers to measure computation time exclusive of MPI communication
!     clabels(1100)='total_excl'
!     clabels(1101)='advection_excl'
!     clabels(1102)='solver_excl'
      IF(icounts(1051).gt.0) THEN
         rtimerbav(1201)=rtimerbav(1018)   &
                        +rtimerbav(1009)   &
                        +rtimerbav(1040)   &
                        +rtimerbav(1401)   &
                        +rtimerbav(1051)   &!*icounts(1018)/icounts(1051)   &
                        +rtimerbav(1053)    !*icounts(1018)/icounts(1053)   
         rtimerbmx(1201)=rtimerbmx(1018)   &
                        +rtimerbmx(1009)   &
                        +rtimerbmx(1040)   &
                        +rtimerbmx(1401)   &
                        +rtimerbmx(1051)   & !*icounts(1018)/icounts(1051)   &
                        +rtimerbmx(1053)     !*icounts(1018)/icounts(1053)   
         rtimerbmn(1201)=rtimerbmn(1018)   &
                        +rtimerbmn(1009)   &
                        +rtimerbmn(1040)   &
                        +rtimerbmn(1401)   &
                        +rtimerbmn(1051)   & !*icounts(1018)/icounts(1051)   &
                        +rtimerbmn(1053)     !*icounts(1018)/icounts(1053)   

         rtimerbav(1202)=rtimerbav(1100)*icounts(1100)/max(icounts(1100),1)   &
                        +rtimerbav(1102)*icounts(1102)/max(icounts(1102),1)   &
                        +rtimerbav(1103)*icounts(1103)/max(icounts(1103),1)   &
                        +rtimerbav(1105)*icounts(1105)/max(icounts(1105),1)   &
                        +rtimerbav(1110)*icounts(1110)/max(icounts(1110),1)   &
                        +rtimerbav(1070)   &
                        +rtimerbav(1071)   & 
                        +rtimerbav(1120)   &
                        +rtimerbav(1121)   
         rtimerbmx(1202)=rtimerbmx(1100)*icounts(1100)/max(icounts(1100),1)   &
                        +rtimerbmx(1102)*icounts(1102)/max(icounts(1102),1)   &
                        +rtimerbmx(1103)*icounts(1103)/max(icounts(1103),1)   &
                        +rtimerbmx(1105)*icounts(1105)/max(icounts(1105),1)   &
                        +rtimerbmx(1110)*icounts(1110)/max(icounts(1110),1)   &
                        +rtimerbmx(1070)   &
                        +rtimerbmx(1071)   & 
                        +rtimerbmx(1120)   &
                        +rtimerbmx(1121)   
         rtimerbmn(1202)=rtimerbmn(1100)*icounts(1100)/max(icounts(1100),1)   &
                        +rtimerbmn(1102)*icounts(1102)/max(icounts(1102),1)   &
                        +rtimerbmn(1103)*icounts(1103)/max(icounts(1103),1)   &
                        +rtimerbmn(1105)*icounts(1105)/max(icounts(1105),1)   &
                        +rtimerbmn(1110)*icounts(1110)/max(icounts(1110),1)   &
                        +rtimerbmn(1070)   &
                        +rtimerbmn(1071)   & 
                        +rtimerbmn(1120)   &
                        +rtimerbmx(1121)   
          ENDIF
          rtimerbav(1200)=rtimerbav(1201)+rtimerbav(1202)
          rtimerbmx(1200)=rtimerbmx(1201)+rtimerbmx(1202)
          rtimerbmn(1200)=rtimerbmn(1201)+rtimerbmn(1202)
          icounts(1200)=icounts(1051)
          icounts(1201)=icounts(1051)
          icounts(1202)=icounts(1051)
         PRINT *, &
          'Pe #, Item #, Subrout,  # calls, .av timer,  mx timer,  mn timer, runit'
                 colend=achar(27)//'['//'0'//'m'
                 colend='    '
         DO i=1,maxtimers
            IF(icounts(i)>0) THEN
               colbegin=achar(27)//'['//clabclr(i)//'m'
                 colbegin='    '
               runit=rtimerbav(i)/icounts(i) 
               PRINT 99,colbegin,mype,i,clabels(i),icounts(i),rtimerbav(i), &
                          rtimerbmx(i),rtimerbmn(i),colend!,runit
            ENDIF
         ENDDO
         PRINT *,'-------------------------------------------------'
         PRINT *,'-----Sorted list of computing cost:--------------'
         PRINT *,'-------------------------------------------------'
        fmax=1.
        j=1
 ttsort:       DO WHILE (fmax.ne.0.AND.j.le.maxtimers)
          fmax=0.
          j=j+1
          DO i=1,maxtimers
            temp=rtimerbav(i)
            fmax=max(fmax,temp)
          ENDDO
          DO i=1,maxtimers
            IF(icounts(i).gt.0) THEN
               IF(fmax.eq.rtimerbav(i)) THEN
                 colbegin=achar(27)//'['//clabclr(i)//'m'
                 colbegin='    '
                 PRINT 99,colbegin,mype,i,clabels(i),icounts(i),rtimerbav(i), &
                          rtimerbmx(i),rtimerbmn(i),colend !,runit
                 rtimerbav(i)=0.
                 EXIT
               ENDIF      
            ENDIF
          ENDDO
              ENDDO ttsort 
      ENDIF


99    FORMAT(1x,a5,1x,i4,1x,i4,4x,a16,2x,i8,1x,f12.4,1x,f12.4,1x,f12.4,1x,a4) !,1x,f12.4)
   END SUBROUTINE ttprt

   SUBROUTINE timer_reduce
      INTEGER i
      REAL nthrdsi
#ifdef PUREMPI
      INTEGER ir
#endif /*PUREMPI*/

      IF(npe>1) THEN
#ifdef PUREMPI
         CALL MPI_ALLReduce(rtimerb,rtimerbav,maxtimers,DC_TYPE_DOUBLE,MPI_SUM,  &
                               mpi_cust_comm,ir)
         DO i=1,maxtimers
            rtimerbav(i)=rtimerbav(i)/npe
         ENDDO
         CALL MPI_ALLReduce(rtimerb,rtimerbmx,maxtimers,DC_TYPE_DOUBLE,MPI_MAX,  &
                               mpi_cust_comm,ir)
         CALL MPI_ALLReduce(rtimerb,rtimerbmn,maxtimers,DC_TYPE_DOUBLE,MPI_MIN,  &
                               mpi_cust_comm,ir)
#endif /*PUREMPI*/
      ELSE
        nthrdsi=1.!./4.
         DO i=1,maxtimers
            rtimerbav(i)=nthrdsi*rtimerb(i)/npe
            rtimerbmx(i)=nthrdsi*rtimerb(i)
            rtimerbmn(i)=nthrdsi*rtimerb(i)
         ENDDO
      ENDIF
   END SUBROUTINE timer_reduce


#ifdef PNETCDF
!--------------------------------------------------------------------!
   SUBROUTINE pnet_out_chunk(var_name,file_name,imode,ia1,ia2,ia3,ifrmd,imem,var,np,mp,lp,ih)
!--------------------------------------------------------------------!
      USE pnetcdf
      USE parameters, ONLY: n, l, m
      !Subroutine writes a specific chunk (vector, 2D plane, 3D variable)
      !Type of chunk is specified by imode:
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
      ! 30 - 0D point
      ! 34 - full 3D variable without time dimension
      ! 35 - 2D xy plane without time dimension
      ! Write of 2D variables is realized in mode 10
      ! KEEP imem=0 EXCEPT WHEN IMOD=2-4 imem=1
      INTEGER(KIND=iintegers) np,mp,lp,ih
      INTEGER(KIND=iintegers) :: ia1,ia2,ia3,imem
      REAL_euwp :: var(1-ih:np+ia1*imem+ih,                          &
                    1-ih:mp+ia2*imem+ih,                          &
                    1-ih:lp+ia3*imem+ih)
      INTEGER(KIND=iintegers) :: p4dim,p3dim,p2dim,p1dim,iframe,ifrmd
      INTEGER(KIND=iintegers) :: DCHVID
      INTEGER(KIND=MPI_OFFSET_KIND) pn,pm,pl,pt,ifrmask
      CHARACTER(len=*) var_name
      CHARACTER(len=*) file_name
      REAL(KIND=sp) dane3d(np,mp,lp)
      REAL(KIND=dp) dane83d(np,mp,lp)
      REAL(KIND=dp) danea83d((np+ia1-1)*imem+1,                                &
                   (mp+ia2-1)*imem+1,(lp+ia3-1)*imem+1)
      REAL(KIND=dp) danex83d((np+ia1-1)*imem+1,mp,lp)
      REAL(KIND=dp) danez83d(np,mp,(lp+ia3-1)*imem+1)
      REAL(KIND=dp) halo8gs(np,mp,ih)
      REAL(KIND=dp) halo8lr(ih,mp,lp)
      REAL(KIND=dp) halo8bt(np,ih,lp)
      REAL(KIND=dp) halo8cr(ih,ih,ih)
      REAL(KIND=dp) dane2dxy(np,mp)
      REAL(KIND=sp) dane2dxz(np,lp)
      REAL(KIND=sp) dane2dyz(mp,lp)
      REAL(KIND=sp) dane1dx(np)
      REAL(KIND=sp) dane1dy(mp)
      REAL(KIND=sp) dane1dz(lp)

      INTEGER(KIND=MPI_OFFSET_KIND) START4(4),COUNT4(4)
      INTEGER(KIND=MPI_OFFSET_KIND) START3(3),COUNT3(3)
      INTEGER(KIND=MPI_OFFSET_KIND) START2(2),COUNT2(2)
      INTEGER                 :: DID(4)
      INTEGER                 :: DID3(3)
      INTEGER                 :: DID2(2)
      INTEGER                 :: DID1(1),iprint
      INTEGER                 :: imode,nfp,nfpd,nfhd,ihnd,i,j,k,ier
      INTEGER                 :: iedge,icorn,i1p,i2p,i3p,i1,i2,i3,ii,jj,kk
  CHARACTER(LEN=*),PARAMETER :: fm2013 = "(1x,'npos,mpos,lpos ',3(1x,I2),"// &
                                             " ' nsubpos,msubpos,lsubpos ',3(1x,I4),"// &
                                             "  '  START4',4(1x,I4), "//& 
                                             "  '  START4+COUNT4',4(2x,I4))" 

      !      IF(imem.ne.0.and.imode.ne.2) &
      !      STOP 'PROVIDE imem=0 in pnet_out_chunk UNLESS imode=2!!!!'
      !      IF(imode.eq.2.and.imem.ne.1) &
      !      STOP 'PROVIDE imem=1 in pnet_out_chunk for imode=2!!!!'
      p4dim=4
      p3dim=3
      p2dim=2
      p1dim=1
      pn=n
      pm=m
      pl=l
      IF (imode >= 2 .AND. imode <= 4) THEN
         pn=n+ia1
         pm=m+ia2
         pl=l+ia3
      ELSE IF(imode == 6) THEN
         pn=2*ih
         pm=m
         pl=l
      ELSE IF(imode == 7) THEN
         pn=n
         pm=2*ih
         pl=l
      ELSE IF(imode == 8) THEN
         pn=n
         pm=m
         pl=2*ih
      ELSE IF(imode == 9) THEN
         pn=2*ih
         pm=2*ih
         pl=2*ih
      ENDIF
      pt=NF_UNLIMITED
      nfp = nf_real
      nfpd= nf_double
      iprint=0
!     IF(mype == 0) print *,'CHUNK var,iframe: ',var_name,ifrmd
!     CALL flush()
!     IF(mype == 0) print *,'ifrmd,np,mp,lp,ih',ifrmd,np,mp,lp,ih
!     CALL flush()
      !-------------------------
      IF (ifrmd == 0) THEN         !if #1 create names and dims
         !-------------------------

         !==================================
         !Now create pnetcdf file if not exist
         ! nfmpi_open will generate MPI error
         ! perhaps implementing CALL system(test -e file_name) is a solution
         !==================================

          ier = nfmpi_open( mpi_cust_comm,file_name,                      &
            NF_WRITE, MPI_INFO_NULL, nfhd)
         IF (ier /= nf_noerr) THEN  !IF #2
            ier = nfmpi_create( mpi_cust_comm,file_name,                  &
               NF_WRITE+NF_64BIT_OFFSET, MPI_INFO_NULL, nfhd)
            IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
               PRINT *,'PNCDF chunk '//file_name//' creat. stat ',        &
               nfmpi_strerror(ier)
         ENDIF   !create file on nf_noerr ENDIF #2

         ier = nfmpi_redef(nfhd)

         !==================================
         !Now inquire or define dimensions in the pnetcdf file
         !==================================

         IF (imode /= 12 .AND. imode /= 20 .AND. imode /= 21 .AND.  imode /= 30 .AND. imode /= 33) THEN !IF #3
            ier=nfmpi_inq_dimid(nfhd,"x",DID(1))
            IF (ier /= nf_noerr) THEN   !IF #4
               ier=nfmpi_def_dim(nfhd,"x",pn,DID(1))
               IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))&
                  PRINT *,'PNCDF dim. x creat. stat ',nfmpi_strerror(ier)
            ENDIF ! x dim create !ENDIF #4
         ENDIF  ! x dim inquire or create !ENDIF #3

         IF (imode /= 11 .AND. imode /= 20 .AND. imode /= 22 .AND.  imode /= 30 .AND. imode /= 33) THEN !IF #5
            ier=nfmpi_inq_dimid(nfhd,"y",DID(2))
            IF (ier /= nf_noerr) THEN  ! IF #6
               ier=nfmpi_def_dim(nfhd,"y",pm,DID(2))
               IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))&
                  PRINT *,'PNCDF dim. y creat. stat ',nfmpi_strerror(ier)
            ENDIF  ! y dim create  ENDIF #6
         ENDIF  ! y dim inquire or create  ENDIF #5

         IF (imode /= 10 .AND. imode /= 21 .AND. imode /= 22 .AND. imode /= 30 .AND. imode /= 13 .AND. imode /= 33) THEN !IF #7
            ier=nfmpi_inq_dimid(nfhd,"z",DID(3))
            IF (ier /= nf_noerr) THEN    !IF #8
               ier=nfmpi_def_dim(nfhd,"z",pl,DID(3))
               IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))&
                  PRINT *,'PNCDF dim. z creat. stat ',nfmpi_strerror(ier)
            ENDIF  ! z dim create ENDIF #8
         ENDIF  ! z dim inquire or create #7

         ier=nfmpi_inq_dimid(nfhd,"t",DID(4))
         IF (ier /= nf_noerr) THEN   !IF#9
            ier=nfmpi_def_dim(nfhd,"t",pt,DID(4))
            IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))&
               PRINT *,'PNCDF dim. t creat. stat ',nfmpi_strerror(ier)
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


            IF (imode == 0) THEN      !if #11
               ier=nfmpi_def_var(nfhd,var_name,nfp,p4dim,DID,DCHVID)
            ELSE IF (imode >= 1 .AND. imode <= 9)THEN
               ier=nfmpi_def_var(nfhd,var_name,nfpd,p4dim,DID,DCHVID)
            ELSE IF (imode == 34) THEN
               DID3(1)=DID(1)
               DID3(2)=DID(2)
               DID3(3)=DID(3)
               ier=nfmpi_def_var(nfhd,var_name,nfp,p3dim,DID3,DCHVID)
            ELSE IF (imode < 20) THEN   !#11
               DID3(3)=DID(4)

               IF (imode == 10 .OR. imode == 13) THEN !IF #12
                  DID3(1)=DID(1)
                  DID3(2)=DID(2)
               ELSE IF (imode == 11) THEN  !#12
                  DID3(1)=DID(1)
                  DID3(2)=DID(3)
               ELSE IF (imode == 12) THEN !#12
                  DID3(1)=DID(2)
                  DID3(2)=DID(3)
               ENDIF   !imode  10-13   !ENDIF #12

               ier=nfmpi_def_var(nfhd,var_name,nfpd,p3dim,DID3,DCHVID)
            ELSE IF (imode < 30) THEN !#11
               DID2(2)=DID(4)

               IF (imode == 20 ) THEN !IF #13
                  DID2(1)=DID(3)
               ELSE IF (imode == 21) THEN !#13
                  DID2(1)=DID(2)
               ELSE IF (imode == 22) THEN  !#13
                  DID2(1)=DID(1)
               ENDIF !imode 20-22 !ENDIF #13

               ier=nfmpi_def_var(nfhd,var_name,nfp,p2dim,DID2,DCHVID)
            ELSE IF (imode >= 30 .AND. imode <= 34) THEN
               DID1(1)=DID(4)
               ier=nfmpi_def_var(nfhd,var_name,nfp,p1dim,DID1,DCHVID)
            ELSE IF (imode == 35) THEN
               DID2(1)=DID(1)
               DID2(2)=DID(2)
               ier=nfmpi_def_var(nfhd,var_name,nfpd,p2dim,DID2,DCHVID)
            ENDIF  !defining dimension set and respective variable dep. on imode !#11

         ENDIF ! no existing variable definition ENDIF #10

         !==================================
         !Now check again if variable exist
         !==================================
         ier=nfmpi_inq_varid(nfhd,var_name,DCHVID)
         IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))    &
          PRINT *,'PNCDF chunk '//file_name//'_'//var_name//' inq ERROR: ', &
            nfmpi_strerror(ier)

         ier = nfmpi_enddef(nfhd)
         IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
            PRINT *,'PNCDF enddef chunk ',nfmpi_strerror(ier)

         ier = nfmpi_close(nfhd)
         IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
            PRINT *,'PNCDF chunk '//file_name//' close stat',               &
            nfmpi_strerror(ier)

         !-------------------------
      ENDIF   !ifrmd == 0
      !-------------------------

      ier = nfmpi_open( mpi_cust_comm,file_name,                       &
         NF_WRITE, MPI_INFO_NULL, nfhd)
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
       PRINT *,'PNCDF chunk '//file_name//'_'//var_name//' open stat ',  &
         nfmpi_strerror(ier)

      ier=nfmpi_inq_varid(nfhd,var_name,DCHVID)
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
       PRINT *,'PNCDF chunk '//file_name//'_'//var_name//' inq. stat ',  &
         nfmpi_strerror(ier)

      ier=nfmpi_inq_unlimdim(nfhd,ihnd)
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0)) &
         PRINT *,'PNCDF chunk. inq unlim dim ',nfmpi_strerror(ier)

      ier=nfmpi_inq_dimlen(nfhd,ihnd,ifrmask)
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0)) &
         PRINT *,'PNCDF chunk. inq dimlen ',nfmpi_strerror(ier)

      !         iframe=ifrmask+1
      iframe=ifrmd+1

      !-------------------------
      !  different mods
      !-------------------------
      IF (imode <= 1) THEN
         START4(1)=nsubpos+1
         START4(2)=msubpos+1
         START4(3)=lsubpos+1
         START4(4)=iframe
         COUNT4(1)=np
         COUNT4(2)=mp
         COUNT4(3)=lp
         COUNT4(4)=1
      ELSE IF (imode >= 2 .AND. imode <= 4) THEN
         START4(1)=nsubpos+1
         START4(2)=msubpos+1
         START4(3)=lsubpos+1
         START4(4)=iframe
         COUNT4(1)=np+ia1*rightedge
         COUNT4(2)=mp+ia2*topedge
         COUNT4(3)=lp+ia3*skyedge
         COUNT4(4)=1
      ELSE IF (imode == 34) THEN
         START3(1)=npos+1
         START3(2)=mpos+1
         START3(3)=lpos+1
         COUNT3(1)=np
         COUNT3(2)=mp
         COUNT3(3)=lp
      ELSE IF (imode == 6) THEN
         START4(1)=1
         START4(2)=mpos+1
         START4(3)=lpos+1
         START4(4)=iframe
         iedge=max(leftedge,rightedge)
         COUNT4(1)=ih*iedge
         COUNT4(2)=mp*iedge
         COUNT4(3)=lp*iedge
         COUNT4(4)=1 *iedge
      ELSE IF (imode == 7) THEN
         START4(1)=npos+1
         START4(2)=1
         START4(3)=lpos+1
         START4(4)=iframe
         iedge=max(botedge,topedge)
         COUNT4(1)=np*iedge
         COUNT4(2)=ih*iedge
         COUNT4(3)=lp*iedge
         COUNT4(4)=1 *iedge
      ELSE IF (imode == 8) THEN
         START4(1)=npos+1
         START4(2)=mpos+1
         START4(3)=1
         START4(4)=iframe
         iedge=max(gndedge,skyedge)
         COUNT4(1)=np*iedge
         COUNT4(2)=mp*iedge
         COUNT4(3)=ih*iedge
         COUNT4(4)=1 *iedge
      ELSE IF (imode == 9) THEN
         START4(1)=1+ ih*rightedge
         START4(2)=1+ ih*topedge
         START4(3)=1+ ih*skyedge
         START4(4)=iframe
         icorn=0
         IF(gndedge == 1 .AND. leftedge  == 1 .AND. botedge == 1 .OR.          &
            gndedge==1.AND.rightedge==1.AND.botedge==1.OR.         &
            gndedge==1.AND.leftedge ==1.AND.topedge==1.OR.         &
            gndedge==1.AND.rightedge==1.AND.topedge==1.OR.         &
            skyedge==1.AND.leftedge ==1.AND.botedge==1.OR.         &
            skyedge==1.AND.rightedge==1.AND.botedge==1.OR.         &
            skyedge==1.AND.leftedge ==1.AND.topedge==1.OR.         &
            skyedge==1.AND.rightedge==1.AND.topedge==1) icorn=1
         COUNT4(1)=ih*icorn
         COUNT4(2)=ih*icorn
         COUNT4(3)=ih*icorn
         COUNT4(4)=1 *icorn
      ELSE IF (imode == 10 .OR. imode == 13) THEN
         START3(1)=npos+1
         START3(2)=mpos+1
         START3(3)=iframe
         COUNT3(1)=np
         COUNT3(2)=mp
         COUNT3(3)=1
      ELSE IF (imode == 11) THEN
         START3(1)=npos+1
         START3(2)=lpos+1
         START3(3)=iframe
         COUNT3(1)=np
         COUNT3(2)=lp
         COUNT3(3)=1
      ELSE IF (imode == 12) THEN
         START3(1)=mpos+1
         START3(2)=lpos+1
         START3(3)=iframe
         COUNT3(1)=mp
         COUNT3(2)=lp
         COUNT3(3)=1
      ELSE IF (imode == 20 ) THEN
         START2(1)=lpos+1
         START2(2)=iframe
         COUNT2(1)=lp
         COUNT2(2)=1
      ELSE IF (imode == 21) THEN
         START2(1)=mpos+1
         START2(2)=iframe
         COUNT2(1)=mp
         COUNT2(2)=1
      ELSE IF (imode == 22) THEN
         START2(1)=npos+1
         START2(2)=iframe
         COUNT2(1)=np
         COUNT2(2)=1
      ELSE IF (imode == 35) THEN
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
      IF (imode > 3) THEN
         i1p=(ia1-1)/np+1
         i2p=(ia2-1)/mp+1
         i3p=(ia3-1)/lp+1
         i1=ia1-(i1p-1)*np
         i2=ia2-(i2p-1)*mp
         i3=ia3-(i3p-1)*lp
      ENDIF
      IF (imode == 0) THEN
         CALL pnettrans(var,dane3d,np,mp,lp,ih)
         ier=nfmpi_put_vara_real_all(nfhd,DCHVID,START4,COUNT4,dane3d)
      ELSE IF(imode == 1) THEN
         CALL pnettrans8(var,dane83d,np,mp,lp,ih)
!    print fm2013,npos,mpos,lpos,nsubpos,msubpos,lsubpos,START4,START4+COUNT4
         CALL flush(6)
         ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,dane83d)
      ELSE IF(imode == 2) THEN
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
!      print *,START4,COUNT4
         ENDIF
      ELSE IF(imode == 3) THEN

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
      ELSE IF(imode == 4) THEN
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
      ELSE IF(imode == 6) THEN
         IF (leftedge == 1) THEN
            DO  k=1,lp
               DO j=1,mp
                  DO  i=1-ih,0
                     halo8lr(ih+i,j,k)=var(i,j,k)
                  ENDDO
               ENDDO
            ENDDO
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
         ENDIF
!    CALL mybarrier()
         ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,halo8lr)
      ELSE IF(imode == 7) THEN
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
!    CALL mybarrier()
         ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,halo8bt)
      ELSE IF(imode == 8) THEN
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
!    CALL mybarrier()
         ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START4,COUNT4,halo8gs)
      ELSE IF(imode == 9) THEN
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
      ELSE IF(imode == 10 .AND. i3p == lpos) THEN
         DO j=1,mp
            DO i=1,np
               dane2dxy(i,j)=var(i,j,i3)
            ENDDO
         ENDDO
         ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START3,COUNT3,dane2dxy)
      ELSE IF(imode == 13 .AND. lpos == 1) THEN
         DO j=1,mp
            DO i=1,np
               dane2dxy(i,j)=var(i,j,1)
            ENDDO
         ENDDO
         ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START3,COUNT3,dane2dxy)
      ELSE IF(imode == 11) THEN
         DO k=1,lp
            DO i=1,np
               dane2dxz(i,k)=REAL(var(i,i2,k),sp)
            ENDDO
         ENDDO
         IF (i2p /= mpos) THEN
            COUNT3(1)=0
            COUNT3(2)=0
            COUNT3(3)=0
         ENDIF
         ier = nfmpi_put_vara_real_all(nfhd,DCHVID,START3,COUNT3,dane2dxz)
      ELSE IF(imode == 12) THEN
         DO k=1,lp
            DO j=1,mp
               dane2dyz(j,k)=REAL(var(i1,j,k),sp)
            ENDDO
         ENDDO
         IF (i1p /= npos) THEN
            COUNT3(1)=0
            COUNT3(2)=0
            COUNT3(3)=0
         ENDIF
         ier = nfmpi_put_vara_real_all(nfhd,DCHVID,START3,COUNT3,dane2dyz)
      ELSE IF(imode == 20) THEN
         DO k=1,lp
            dane1dz(k)=REAL(var(i1,i2,k),sp)
         ENDDO
         IF (i1p /= npos .OR. i2p /= mpos) THEN
            COUNT2(1)=0
            COUNT2(2)=0
         ENDIF
         ier = nfmpi_put_vara_real_all(nfhd,DCHVID,START2,COUNT2,dane1dz)

      ELSE IF(imode == 21) THEN
         DO j=1,mp
            dane1dy(j)=REAL(var(i1,j,i3),sp)
         ENDDO
         IF (i1p /= npos .OR. i3p /= lpos) THEN
            COUNT2(1)=0
            COUNT2(2)=0
         ENDIF

         ier=nfmpi_put_vara_real_all(nfhd,DCHVID,START2,COUNT2,dane1dy)
      ELSE IF(imode == 22) THEN
         DO i=1,np
            dane1dx(i)=REAL(var(i,i2,i3),sp)
         ENDDO
         IF (i2p /= mpos .OR. i3p /= lpos) THEN
            COUNT2(1)=0
            COUNT2(2)=0
         ENDIF
         ier=nfmpi_put_vara_real_all(nfhd,DCHVID,START2,COUNT2,dane1dx)
      ELSE IF(imode == 34) THEN
         CALL pnettrans(var,dane3d,np,mp,lp,ih)
         ier=nfmpi_put_vara_real_all(nfhd,DCHVID,START3,COUNT3,dane3d)
      ELSE IF(imode == 35) THEN
         DO j=1,mp
            DO i=1,np
               dane2dxy(i,j)=var(i,j,1-ih)
            ENDDO
         ENDDO
         ier=nfmpi_put_vara_double_all(nfhd,DCHVID,START2,COUNT2,dane2dxy)

      ENDIF

      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
         PRINT *,'PNCDF chunk '//var_name//' put stat ',               &
         nfmpi_strerror(ier)

      ier = nfmpi_close(nfhd)
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
         PRINT *,'PNCDF chunk '//file_name//' close stat',               &
         nfmpi_strerror(ier)

   END SUBROUTINE pnet_out_chunk

!--------------------------------------------------------------------!
   SUBROUTINE pnettrans(sour,dest,np,mp,lp,ih)
!--------------------------------------------------------------------!
      INTEGER(KIND=iintegers) np,mp,lp,ih
      REAL(KIND=euwp) sour(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      REAL(KIND=sp)  dest(np,mp,lp)
      dest(1:np,1:mp,1:lp)=REAL(sour(1:np,1:mp,1:lp),sp)
   END SUBROUTINE pnettrans


!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
   SUBROUTINE pnettrans8(sour,dest,np,mp,lp,ih)
!--------------------------------------------------------------------!
      INTEGER(KIND=iintegers) np,mp,lp,ih
      REAL_euwp :: sour(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
      REAL(KIND=dp) dest(np,mp,lp)
      dest(1:np,1:mp,1:lp)=REAL(sour(1:np,1:mp,1:lp),dp)
   END SUBROUTINE pnettrans8

!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
!--------------------------------------------------------------------!
SUBROUTINE pnet_create_common(ipnind,nfhd,DID,tVID)
!--------------------------------------------------------------------!
  USE pnetcdf
  USE parameters, ONLY: n,m,l
  INTEGER(KIND=MPI_OFFSET_KIND) pn,pm,pl,pt
  INTEGER(KIND=iintegers) ier,iprint,ipnind,cmode,nfhd,tVID,DID(4)
  INTEGER(KIND=iintegers) nfp
  !Cast domain size into OFFSET_KIND integers
  pn=n
  pm=m
  pl=l
  pt=NF_UNLIMITED
  nfp = NF_REAL
  cmode=NF_WRITE+NF_64BIT_OFFSET

! iprint=1 ! PRINT out all error messages including NF_NOERR
   iprint=0 ! PRINT out all real error messages

  !For compatibility with older pre Netcdf 3.6 (pre 2007) set cmode=NF_WRITE
  !This would limit the file size to 2 GB as permitted by CDF-1 format
  IF (ipnind == 1) THEN
    ier = nfmpi_create( mpi_cust_comm,'tape.custom.nc',                &
        cmode, MPI_INFO_NULL, nfhd)
  ELSE IF(ipnind == 2) Then
    ier = nfmpi_create( mpi_cust_comm,'tapes.nc',                      &
        cmode, MPI_INFO_NULL, nfhd)
  ELSE IF(ipnind == 3) THEN
    ier = nfmpi_create( mpi_cust_comm,'tapef.nc',                      &
        cmode, MPI_INFO_NULL, nfhd)
  ELSE IF(ipnind == 4) THEN
    ier = nfmpi_create( mpi_cust_comm,'tape.slice.nc',                 &
        cmode, MPI_INFO_NULL, nfhd)
  ENDIF
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
            print *,'PNCDF tape',ipnind,'creat stat ',nfmpi_strerror(ier)

  !Now defining dimensions of variables (x,y,z,t) or (x,z,t)
  ier=nfmpi_def_dim(nfhd,"x",pn,DID(1))
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF dim. x creat. stat ',nfmpi_strerror(ier)

  ier=nfmpi_def_dim(nfhd,"y",pm,DID(2))
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF dim. y creat. stat ',nfmpi_strerror(ier)

  ier=nfmpi_def_dim(nfhd,"z",pl,DID(3))
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF dim. z creat. stat ',nfmpi_strerror(ier)

  ier=nfmpi_def_dim(nfhd,"t",pt,DID(4))
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF dim. t creat. stat ',nfmpi_strerror(ier)

  ier=nfmpi_def_var(nfhd,"time",nfp,1,DID(4),tVID)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0)) &
           print *,'PNCDF var time creat. stat ',nfmpi_strerror(ier)

  ier = nfmpi_enddef(nfhd)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0)) &
           print *,'PNCDF enddef stat ',nfmpi_strerror(ier)

END SUBROUTINE pnet_create_common



!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
SUBROUTINE pnet_open_common(ipnind,nfhd,DID,tVID,iframe)
!--------------------------------------------------------------------!
  USE pnetcdf 
  INTEGER(KIND=iintegers) ier,iprint,ipnind,cmode,nfhd,tVID
  INTEGER(KIND=iintegers),INTENT(IN) :: iframe
  INTEGER(KIND=iintegers) DID(4)
  iprint=1
  cmode=NF_WRITE  !+NF_64BIT_OFFSET
  !For compatibility with older pre Netcdf 3.6 (pre 2007) set cmode=NF_WRITE
  !This would limit the file size to 2 GB as permitted by CDF-1 format
  ier=nf_noerr
  IF (iframe > 0 .OR. iframe == -1) THEN
    IF (ipnind == 1) THEN
      ier = nfmpi_open( mpi_cust_comm,'tape.custom.nc',                  &
          cmode, MPI_INFO_NULL, nfhd)
    ELSE IF(ipnind == 2) THEN
      ier = nfmpi_open( mpi_cust_comm,'tapes.nc',                        &
          cmode, MPI_INFO_NULL, nfhd)
    ELSE IF(ipnind == 3) THEN
      ier = nfmpi_open( mpi_cust_comm,'tapef.nc',                        &
          cmode, MPI_INFO_NULL, nfhd)
    ELSE IF(ipnind == 4) THEN
      ier = nfmpi_open( mpi_cust_comm,'tape.slice.nc',                   &
          cmode, MPI_INFO_NULL, nfhd)
    ENDIF
  ENDIF
  
  !If file doesn''t exist, create it
  IF (ier /= nf_noerr  .OR.  iframe == 0) THEN
    CALL pnet_create_common(ipnind,nfhd,DID,tVID) 
  ELSE
    !Now inquiring for dimension and dimension variables (x,y,z,t) handles
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))    &
             print *,'PNCDF tape',ipnind,'open stat ',nfmpi_strerror(ier)
    ier=nfmpi_inq_dimid(nfhd,"x",DID(1))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0)) &
             print *,'PNCDF dim. x inq. stat ',nfmpi_strerror(ier)
    ier=nfmpi_inq_dimid(nfhd,"y",DID(2))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0)) &
             print *,'PNCDF dim. y inq. stat ',nfmpi_strerror(ier)
    ier=nfmpi_inq_dimid(nfhd,"z",DID(3))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0)) &
             print *,'PNCDF dim. z inq. stat ',nfmpi_strerror(ier)
    ier=nfmpi_inq_dimid(nfhd,"t",DID(4))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0)) &
             print *,'PNCDF dim. t inq. stat ',nfmpi_strerror(ier)
    
    ier=nfmpi_inq_varid(nfhd,"time",tVID)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0)) &
             print *,'PNCDF time inq. stat ',nfmpi_strerror(ier)
  ENDIF
END SUBROUTINE pnet_open_common

#if(1==0)

SUBROUTINE pnet_out_lng_cmpr(itapetype,np,mp,lp,ih,iframe0, u,v,w,ox,oy,oz,rh,th,tht,pstr,pext,  &
                                fx,fy,fz,ft,fxexpl,fyexpl,fzexpl,ftexpl,fpextexpl, &
                                qv,qv_old,qc,qr,qs,qi,qg)     
!--------------------------------------------------------------------!
  USE pnetcdf 
  USE eulag_datafields,     ONLY: zcr
  USE eulag_datafields,     ONLY: the
  INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
  INTEGER(KIND=iintegers),INTENT(IN) :: itapetype 
  
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(IN) ::  &
             tht,pstr,pext
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(IN),OPTIONAL ::  &
             fx,fy,fz,ft,fxexpl,fyexpl,fzexpl,ftexpl,fpextexpl
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(INOUT) ::  &
             th
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(IN),OPTIONAL ::  &
              qv,qv_old,qc
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(IN),OPTIONAL ::  &
            qr,qs,qi,qg
    REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:1):: &
            rh 
    REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2):: &
             u,v,w,ox,oy,oz
      
  INTEGER(KIND=iintegers) :: nfhd,ihnd,p3dim,p4dim,p5dim,iframe0,iframe,id
  INTEGER(KIND=iintegers) :: DLVID(50),tVID
  
  integer(KIND=MPI_OFFSET_KIND) pframe,pone,ifrmask
  REAL(KIND=euwp) :: dane(np,mp,lp)        ! real*8
  
  integer(KIND=MPI_OFFSET_KIND) START3(3),COUNT3(3)
  integer(KIND=MPI_OFFSET_KIND) START4(4),COUNT4(4)
  INTEGER(KIND=iintegers) :: DID(4),DID5(5)
  INTEGER(KIND=iintegers) :: iprint,irsdta,nfp,ier,i
  !Switch for exact restart option
  !Default off, because at the moment not truly exact
  irsdta=1
  iframe=iframe0
  
  nfp = nf_double
  pone=1
  iprint=0
  
  
    IF(mype == 0 .AND. iprint == 1)   &   
   PRINT *,'iframe before open  is',iframe
  CALL  pnet_open_common(itapetype,nfhd,DID,tVID,iframe)
    IF(mype == 0 .AND. iprint == 1)   &   
   PRINT *,'iframe after open  is',iframe
  !Automatic iframe setting
  IF (iframe == 0) THEN
    ier=nfmpi_inq_unlimdim(nfhd,ihnd)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &   
             print *,'PNCDF chunk. inq unlim dim id ',nfmpi_strerror(ier)
    
    ier=nfmpi_inq_dimlen(nfhd,ihnd,ifrmask)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
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
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF first inq u  stat ',nfmpi_strerror(ier)
  !If not, we define them now
  IF (iframe == 1 .AND. ier /= nf_noerr) THEN
    ier = nfmpi_redef(nfhd)
    id=1
    ier=nfmpi_def_var(nfhd,'u',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape u  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'v',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape v  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'w',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape w  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'ox',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ox stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'oy',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape oy  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'oz',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape oz  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'rh0',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape rh0 stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'rh1',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape rh1 stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'th',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape th  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'tht',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape tht  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'pstr',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape  pstr stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'pext',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape  pext stat ',nfmpi_strerror(ier)
    IF(PRESENT(fx)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fx',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fx  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fy)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fy',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fy  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fz)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fz',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fz  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(ft)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'ft',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ft  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fxexpl)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fxexpl',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fxexpl  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fyexpl)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fyexpl',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fyexpl  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fzexpl)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fzexpl',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fzexpl  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(ftexpl)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'ftexpl',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ftexpl  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fpextexpl)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fpextexpl',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fpextexpl  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qv)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qv',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qv  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qc)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qv_old',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qv_old  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qc)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qc',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qc  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qr)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qr',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qr  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qs)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qs',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qs  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qi)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qi',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qi  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qg)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qg',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qg  stat ',nfmpi_strerror(ier)
    ENDIF
    
    IF (iframe == 1) THEN
      id=id+1
      ier=nfmpi_def_var(nfhd,"ELEVATION",nfp,p3dim,DID,DLVID(id)) 
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape ELEV stat ',nfmpi_strerror(ier)
    ENDIF
    IF (irsdta == 1) THEN
      id=id+1
      ier=nfmpi_def_var(nfhd,'u2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape u2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'v2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape v2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'w2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape w2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'ox2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape ox2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oy2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oy2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oz2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oz2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'u3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape u3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'v3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape v3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'w3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape w3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'ox3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape ox3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oy3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oy3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oz3',nfp,p4dim,DID,DLVID(id))        
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oz3  stat ',nfmpi_strerror(ier)
    ENDIF   !irsdta
    ier=nfmpi_enddef(nfhd)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF open_common enddef stat ',nfmpi_strerror(ier)
    
    
  ELSE   !iframe > 1 
    
    id=1
    ier=nfmpi_inq_varid(nfhd,'u',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape u  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'v',DLVID(id))                      
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape v  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'w',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape w  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'ox',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ox inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'oy',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape oy  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'oz',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape oz  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'rh0',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape rh0 inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'rh1',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape rh1 inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'th',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape th  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'tht',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape tht  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'pstr',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape  pstr inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'pext',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape  pext inq ',nfmpi_strerror(ier)
    IF(PRESENT(fx)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fx',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fx  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fy)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fy',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fy  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fz)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fz',DLVID(id)) 
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fz  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(ft)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'ft',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ft  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fxexpl)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fxexpl',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fxexpl  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fyexpl)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fyexpl',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fyexpl  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fzexpl)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fzexpl',DLVID(id)) 
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fzexpl  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(ftexpl)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'ftexpl',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ftexpl  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(ftexpl)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fpextexpl',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fpextexpl  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qv)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qv',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qv  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qv_old)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qv_old',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qv_old  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qc)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qc',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qc  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qr)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qr',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qr  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qs)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qs',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qs  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qi)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qi',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qi  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qg)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qg',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qg  stat ',nfmpi_strerror(ier)
    ENDIF
    
    IF (iframe == 1) THEN
      id=id+1
      ier=nfmpi_inq_varid(nfhd,"ELEVATION",DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
               print *,'PNCDF tape ELEV inq ',nfmpi_strerror(ier)
    ENDIF 
    IF (irsdta == 1) THEN
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'u2',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape u2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'v2',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape v2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'w2',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape w2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'ox2',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape ox2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oy2',DLVID(id)) 
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oy2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oz2',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oz2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'u3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape u3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'v3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape v3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'w3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape w3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'ox3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape ox3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oy3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oy3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oz3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oz3  inq ',nfmpi_strerror(ier)
    ENDIF   ! irsdta
  ENDIF     !iframe
  !      ier=nfmpi_begin_indep_data(nfhd)
  !      ier = nfmpi_put_vara_real(nfhd,tVID,pframe,pone,time)
  !      IF((mype.eq.0.and.iprint.eq.1).or.(iprint.eq.0.and.ier.ne.0))
  !     &   print *,'PNCDF std var. time put stat ',nfmpi_strerror(ier)
  !      ier=nfmpi_end_indep_data(nfhd)
  
  id=1
  CALL pnettrans8(u(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape u wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(v(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape v wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(w(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape w wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(ox(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape ox wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(oy(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape oy wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(oz(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape oz wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(rh(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape rh0 wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(rh(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape rh1 wrt ',nfmpi_strerror(ier)
  id=id+1
  th(1:np,1:mp,1:lp)=th(1:np,1:mp,1:lp)+the(1:np,1:mp,1:lp)
  CALL pnettrans8(th,dane,np,mp,lp,ih)
  th(1:np,1:mp,1:lp)=th(1:np,1:mp,1:lp)-the(1:np,1:mp,1:lp)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape th wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(tht,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape tht wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(pstr,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape pstr wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(pext,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape pext wrt ',nfmpi_strerror(ier)
  IF(PRESENT(fx)) THEN
  id=id+1
  CALL pnettrans8(fx,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fx wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fy)) THEN
  id=id+1
  CALL pnettrans8(fy,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fy wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fz)) THEN
  id=id+1
  CALL pnettrans8(fz,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fz wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(ft)) THEN
  id=id+1
  CALL pnettrans8(ft,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape ft wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fxexpl)) THEN
  id=id+1
  CALL pnettrans8(fxexpl,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fxexpl wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fyexpl)) THEN
  id=id+1
  CALL pnettrans8(fyexpl,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fyexpl wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fzexpl)) THEN
  id=id+1
  CALL pnettrans8(fzexpl,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fzexpl wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(ftexpl)) THEN
  id=id+1
  CALL pnettrans8(ftexpl,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape ftexpl wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fpextexpl)) THEN
  id=id+1
  CALL pnettrans8(fpextexpl,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fpextexpl wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qv)) THEN
  id=id+1
  CALL pnettrans8(qv,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qv wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qv_old)) THEN
  id=id+1
  CALL pnettrans8(qv_old,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qv_old wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qc)) THEN
  id=id+1
  CALL pnettrans8(qc,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qc wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qr)) THEN
  id=id+1
  CALL pnettrans8(qr,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qr wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qs)) THEN
  id=id+1
  CALL pnettrans8(qs,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qs wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qi)) THEN
  id=id+1
  CALL pnettrans8(qi,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qi wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qg)) THEN
  id=id+1
  CALL pnettrans8(qg,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qg wrt ',nfmpi_strerror(ier)
  ENDIF
  IF (iframe == 1) THEN
    id=id+1
    dane(1:np,1:mp,1:lp)=zcr(1:np,1:mp,1:lp)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START3,COUNT3,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape ELEV wrt ',nfmpi_strerror(ier)
  ENDIF
  
  IF (irsdta == 1) THEN 
    id=id+1
    CALL pnettrans8(u(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape u2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(v(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape v2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(w(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape w2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(ox(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape ox2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oy(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape oy2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oz(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape oz2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(u(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape u3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(v(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape v3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(w(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape w3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(ox(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape ox3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oy(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape oy3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oz(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape oz3 wrt ',nfmpi_strerror(ier)
  ENDIF   ! irsdta
  
  ier = nfmpi_close(nfhd)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'pnet out lng close',nfmpi_strerror(ier)
  
END SUBROUTINE pnet_out_lng_cmpr

SUBROUTINE pnet_out_lng_anel(itapetype,np,mp,lp,ih,iframe0, u,v,w,ox,oy,oz,rh,th,  &
                                fx,fy,fz,ft,fxexpl,fyexpl,fzexpl,ftexpl, &
                                qv,qv_old,qc,qr,qs,qi,qg)     
!--------------------------------------------------------------------!
  USE pnetcdf 
  USE eulag_datafields,     ONLY: zcr
  USE eulag_datafields,     ONLY: the
  INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
  INTEGER(KIND=iintegers),INTENT(IN) :: itapetype 
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(IN), OPTIONAL ::  &
             fx,fy,fz,ft,fxexpl,fyexpl,fzexpl,ftexpl
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(INOUT) ::  &
             th
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(IN),OPTIONAL ::  &
              qv,qv_old,qc
  REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih),INTENT(IN),OPTIONAL ::  &
            qr,qs,qi,qg
    REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:1):: &
            rh 
    REAL(KIND=euwp),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih,0:2):: &
             u,v,w,ox,oy,oz
      
  INTEGER(KIND=iintegers) :: nfhd,ihnd,p3dim,p4dim,p5dim,iframe0,iframe,id
  INTEGER(KIND=iintegers) :: DLVID(50),tVID
  
  integer(KIND=MPI_OFFSET_KIND) pframe,pone,ifrmask
  REAL(KIND=euwp) :: dane(np,mp,lp)        ! real*8
  
  integer(KIND=MPI_OFFSET_KIND) START3(3),COUNT3(3)
  integer(KIND=MPI_OFFSET_KIND) START4(4),COUNT4(4)
  INTEGER(KIND=iintegers) :: DID(4),DID5(5)
  INTEGER(KIND=iintegers) :: iprint,irsdta,nfp,ier,i
  !Switch for exact restart option
  !Default off, because at the moment not truly exact
  irsdta=1
  iframe=iframe0
  
  nfp = nf_double
  pone=1
  iprint=0
  
  
    IF(mype == 0 .AND. iprint == 1)   &   
   PRINT *,'iframe before open  is',iframe
  CALL  pnet_open_common(itapetype,nfhd,DID,tVID,iframe)
    IF(mype == 0 .AND. iprint == 1)   &   
   PRINT *,'iframe after open  is',iframe
  !Automatic iframe setting
  IF (iframe == 0) THEN
    ier=nfmpi_inq_unlimdim(nfhd,ihnd)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &   
             print *,'PNCDF chunk. inq unlim dim id ',nfmpi_strerror(ier)
    
    ier=nfmpi_inq_dimlen(nfhd,ihnd,ifrmask)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
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
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF first inq u  stat ',nfmpi_strerror(ier)
  !If not, we define them now
  IF (iframe == 1 .AND. ier /= nf_noerr) THEN
    ier = nfmpi_redef(nfhd)
    id=1
    ier=nfmpi_def_var(nfhd,'u',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape u  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'v',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape v  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'w',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape w  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'ox',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ox stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'oy',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape oy  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'oz',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape oz  stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'rh0',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape rh0 stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'rh1',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape rh1 stat ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_def_var(nfhd,'th',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape th  stat ',nfmpi_strerror(ier)
    IF(PRESENT(fx)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fx',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fx  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fy)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fy',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fy  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fz)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fz',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fz  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(ft)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'ft',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ft  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fxexpl)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fxexpl',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fxexpl  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fyexpl)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fyexpl',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fyexpl  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fzexpl)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'fzexpl',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fzexpl  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(ftexpl)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'ftexpl',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ftexpl  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qv)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qv',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qv  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qc)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qv_old',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qv_old  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qc)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qc',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qc  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qr)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qr',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qr  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qs)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qs',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qs  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qi)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qi',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qi  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qg)) THEN
    id=id+1
    ier=nfmpi_def_var(nfhd,'qg',nfp,p4dim,DID,DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qg  stat ',nfmpi_strerror(ier)
    ENDIF
    
    IF (iframe == 1) THEN
      id=id+1
      ier=nfmpi_def_var(nfhd,"ELEVATION",nfp,p3dim,DID,DLVID(id)) 
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape ELEV stat ',nfmpi_strerror(ier)
    ENDIF
    IF (irsdta == 1) THEN
      id=id+1
      ier=nfmpi_def_var(nfhd,'u2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape u2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'v2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape v2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'w2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape w2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'ox2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape ox2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oy2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oy2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oz2',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oz2  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'u3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape u3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'v3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape v3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'w3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape w3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'ox3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape ox3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oy3',nfp,p4dim,DID,DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oy3  stat ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_def_var(nfhd,'oz3',nfp,p4dim,DID,DLVID(id))        
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oz3  stat ',nfmpi_strerror(ier)
    ENDIF   !irsdta
    ier=nfmpi_enddef(nfhd)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF open_common enddef stat ',nfmpi_strerror(ier)
    
    
  ELSE   !iframe > 1 
    
    id=1
    ier=nfmpi_inq_varid(nfhd,'u',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape u  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'v',DLVID(id))                      
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape v  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'w',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape w  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'ox',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ox inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'oy',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape oy  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'oz',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape oz  inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'rh0',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape rh0 inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'rh1',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape rh1 inq ',nfmpi_strerror(ier)
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'th',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape th  inq ',nfmpi_strerror(ier)
    IF(PRESENT(fx)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fx',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fx  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fy)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fy',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fy  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fz)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fz',DLVID(id)) 
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fz  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(ft)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'ft',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ft  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fxexpl)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fxexpl',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fxexpl  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fyexpl)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fyexpl',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fyexpl  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(fzexpl)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'fzexpl',DLVID(id)) 
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape fzexpl  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(ftexpl)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'ftexpl',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape ftexpl  inq ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qv)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qv',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qv  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qv_old)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qv_old',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qv_old  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qc)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qc',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qc  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qr)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qr',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qr  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qs)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qs',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qs  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qi)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qi',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qi  stat ',nfmpi_strerror(ier)
    ENDIF
    IF(PRESENT(qg)) THEN
    id=id+1
    ier=nfmpi_inq_varid(nfhd,'qg',DLVID(id))
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
             print *,'PNCDF tape qg  stat ',nfmpi_strerror(ier)
    ENDIF
    
    IF (iframe == 1) THEN
      id=id+1
      ier=nfmpi_inq_varid(nfhd,"ELEVATION",DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
               print *,'PNCDF tape ELEV inq ',nfmpi_strerror(ier)
    ENDIF 
    IF (irsdta == 1) THEN
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'u2',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape u2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'v2',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape v2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'w2',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape w2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'ox2',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape ox2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oy2',DLVID(id)) 
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oy2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oz2',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oz2  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'u3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape u3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'v3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape v3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'w3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape w3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'ox3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape ox3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oy3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oy3  inq ',nfmpi_strerror(ier)
      id=id+1
      ier=nfmpi_inq_varid(nfhd,'oz3',DLVID(id))
      IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))  &
               print *,'PNCDF tape oz3  inq ',nfmpi_strerror(ier)
    ENDIF   ! irsdta
  ENDIF     !iframe
  !      ier=nfmpi_begin_indep_data(nfhd)
  !      ier = nfmpi_put_vara_real(nfhd,tVID,pframe,pone,time)
  !      IF((mype.eq.0.and.iprint.eq.1).or.(iprint.eq.0.and.ier.ne.0))
  !     &   print *,'PNCDF std var. time put stat ',nfmpi_strerror(ier)
  !      ier=nfmpi_end_indep_data(nfhd)
  
  id=1
  CALL pnettrans8(u(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape u wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(v(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape v wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(w(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape w wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(ox(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape ox wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(oy(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape oy wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(oz(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape oz wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(rh(1-ih,1-ih,1-ih,0),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape rh0 wrt ',nfmpi_strerror(ier)
  id=id+1
  CALL pnettrans8(rh(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape rh1 wrt ',nfmpi_strerror(ier)
  id=id+1
  th(1:np,1:mp,1:lp)=th(1:np,1:mp,1:lp)+the(1:np,1:mp,1:lp)
  CALL pnettrans8(th,dane,np,mp,lp,ih)
  th(1:np,1:mp,1:lp)=th(1:np,1:mp,1:lp)-the(1:np,1:mp,1:lp)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape th wrt ',nfmpi_strerror(ier)
  IF(PRESENT(fx)) THEN
  id=id+1
  CALL pnettrans8(fx,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fx wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fy)) THEN
  id=id+1
  CALL pnettrans8(fy,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fy wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fz)) THEN
  id=id+1
  CALL pnettrans8(fz,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fz wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(ft)) THEN
  id=id+1
  CALL pnettrans8(ft,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape ft wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fxexpl)) THEN
  id=id+1
  CALL pnettrans8(fxexpl,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fxexpl wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fyexpl)) THEN
  id=id+1
  CALL pnettrans8(fyexpl,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fyexpl wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(fzexpl)) THEN
  id=id+1
  CALL pnettrans8(fzexpl,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape fzexpl wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(ftexpl)) THEN
  id=id+1
  CALL pnettrans8(ftexpl,dane,np,mp,lp,ih)
  ier = nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'PNCDF tape ftexpl wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qv)) THEN
  id=id+1
  CALL pnettrans8(qv,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qv wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qv_old)) THEN
  id=id+1
  CALL pnettrans8(qv_old,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qv_old wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qc)) THEN
  id=id+1
  CALL pnettrans8(qc,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qc wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qr)) THEN
  id=id+1
  CALL pnettrans8(qr,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qr wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qs)) THEN
  id=id+1
  CALL pnettrans8(qs,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qs wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qi)) THEN
  id=id+1
  CALL pnettrans8(qi,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qi wrt ',nfmpi_strerror(ier)
  ENDIF
  IF(PRESENT(qg)) THEN
  id=id+1
  CALL pnettrans8(qg,dane,np,mp,lp,ih)
  ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
           print *,'PNCDF tape qg wrt ',nfmpi_strerror(ier)
  ENDIF
  IF (iframe == 1) THEN
    id=id+1
    dane(1:np,1:mp,1:lp)=zcr(1:np,1:mp,1:lp)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START3,COUNT3,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape ELEV wrt ',nfmpi_strerror(ier)
  ENDIF
  
  IF (irsdta == 1) THEN 
    id=id+1
    CALL pnettrans8(u(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape u2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(v(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape v2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(w(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape w2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(ox(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape ox2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oy(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape oy2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oz(1-ih,1-ih,1-ih,1),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape oz2 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(u(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape u3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(v(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape v3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(w(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape w3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(ox(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape ox3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oy(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape oy3 wrt ',nfmpi_strerror(ier)
    id=id+1
    CALL pnettrans8(oz(1-ih,1-ih,1-ih,2),dane,np,mp,lp,ih)
    ier=nfmpi_put_vara_double_all(nfhd,DLVID(id),START4,COUNT4,dane)
    IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))   &
             print *,'PNCDF tape oz3 wrt ',nfmpi_strerror(ier)
  ENDIF   ! irsdta
  
  ier = nfmpi_close(nfhd)
  IF((mype == 0 .AND. iprint == 1) .OR. (iprint == 0 .AND. ier /= 0))     &
           print *,'pnet out lng close',nfmpi_strerror(ier)
  
   END SUBROUTINE pnet_out_lng_anel
#endif /*1==0*/
#endif /*PNETCDF*/
    SUBROUTINE enforce_cyclic(x,ipoles0,ibcx0,ibcy0,ibcz0,np,mp,lp,ih)
      INTEGER(KIND=iintegers),INTENT(IN) :: ibcx0,ibcy0,ibcz0,ipoles0,np,mp,lp,ih
      REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: x

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
      INTEGER(KIND=iintegers),INTENT(IN) :: ibcz0,np,mp,lp,ih
      REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: x
      IF(ibcz0.eq.1) THEN
          CALL updategs(x,np,mp,lp,np,mp,lp,1,ih)
        IF ( skyedge == 1) x(1:np,1:mp,lp)=x(1:np,1:mp,lp+1)
      ENDIF
    END SUBROUTINE enforce_cyclic_z
    SUBROUTINE enforce_cyclic_xy(x,ipoles0,ibcx0,ibcy0,np,mp,lp,ih)
      INTEGER(KIND=iintegers),INTENT(IN) :: ibcx0,ibcy0,ipoles0,np,mp,lp,ih
      REAL_euwp,DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: x

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

      SUBROUTINE rcvbufferz(dest,i,j,k,ihx,ihy,ihz,icnt,lstep,np,mp,lp,ih)
      INTEGER,INTENT(IN) :: ihx,ihy,ihz,i,j,k,icnt,lstep
      INTEGER,INTENT(IN) :: np,mp,lp,ih 
      INTEGER :: ipercv,itgrcv,itgoff
      REAL(KIND=euwp) dest(1-ihx:np+ihx,1-ihy:mp+ihy,1-ihz:lp+ihz)
      LOGICAL rcv
      INTEGER :: request
#ifdef PUREMPI
      INTEGER :: status(MPI_STATUS_SIZE),ierr
#endif
      itgoff=1000*j !offset for MPI message tag to avoid conflicts with updates 
      rcv=.FALSE. !initialize to remove compiler warning
      if(lstep.gt.0) then
        rcv=gndedge.eq.0 !recv if not at gndedge
        ipercv=peG! recv from precessor lpos-2
        itgrcv=itgoff+npos*100+mpos*10+lpos-2 !message tag
      else if(lstep.lt.0) then
        rcv=skyedge.eq.0  !recv if not at skyedge
        ipercv=peZ !recv from precessor lpos
        itgrcv=itgoff+npos*100+mpos*10+lpos-1 !message tag
      endif
      IF(rcv) THEN
#ifdef PUREMPI
!      print *,'Before receive from',ipercv,'on',mype 
!     CALL flush(6)
         call MPI_IRecv(dest(i,j,k),icnt,DC_TYPE,ipercv,itgrcv,  &
                                      my_prcnm, request, ierr)
         call MPI_WAIT(request,status,ierr)
!      print *,'After receive from',ipercv,'on',mype 
!     CALL flush(6)
#endif
      ENDIF  ! rcv
      END SUBROUTINE rcvbufferz

      SUBROUTINE sndbufferz(dest,i,j,k,ihx,ihy,ihz,icnt,lstep,np,mp,lp,ih)
      INTEGER,INTENT(IN) :: ihx,ihy,ihz,i,j,k,icnt,lstep
      INTEGER,INTENT(IN) :: np,mp,lp,ih 
      INTEGER :: ierr,ipesnd,itgsnd,itgoff,istatus,status 
      REAL(KIND=euwp) dest(1-ihx:np+ihx,1-ihy:mp+ihy,1-ihz:lp+ihz)
      LOGICAL snd
      itgoff=1000*j !offset for MPI message tag to avoid conflicts with updates 
      snd=.FALSE. !initialize to remove compiler warning 
      if(lstep.gt.0) then
        snd=skyedge.eq.0 !send if not at skyedge
        ipesnd=peZ  ! send to processor lpos
        itgsnd=itgoff+npos*100+mpos*10+lpos-1 !message tag
      else if(lstep.lt.0) then
        snd=gndedge.eq.0 ! send if not at gndedge
        ipesnd=peG  !send to precessor lpos-2
        itgsnd=itgoff+npos*100+mpos*10+lpos-2 !message tag
      endif
      if(snd) then
#ifdef PUREMPI
!      print *,'Before send to',ipesnd,'on',mype 
!     CALL flush(6)
        call MPI_Send(dest(i,j,k),icnt,DC_TYPE,ipesnd,itgsnd,  &
                                              my_prcnm, ierr)
!      print *,'After send to',ipesnd,'on',mype 
!     CALL flush(6)
#endif
      endif ! snd 
      END SUBROUTINE sndbufferz

END MODULE mpi_parallel

#ifdef CUDACODE
module nvtx

use iso_c_binding
implicit none

integer,private :: col(7) = [ Z'0000ff00', Z'000000ff', Z'00ffff00', Z'00ff00ff', Z'0000ffff', Z'00ff0000', Z'00ffffff']
character,private,target :: tempName(256)

type, bind(C):: nvtxEventAttributes
  integer(C_INT16_T):: version=1
  integer(C_INT16_T):: size=48 !
  integer(C_INT):: category=0
  integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
  integer(C_INT):: color
  integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
  integer(C_INT):: reserved0
  integer(C_INT64_T):: payload   ! union uint,int,double
  integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1 
  type(C_PTR):: message  ! ascii char
end type

interface nvtxRangePush
  ! push range with custom label and standard color
  subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
  use iso_c_binding
  character(kind=C_CHAR) :: name(256)
  end subroutine

  ! push range with custom label and custom color
  subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
  use iso_c_binding
  import:: nvtxEventAttributes
  type(nvtxEventAttributes):: event
  end subroutine
end interface

interface nvtxRangePop
  subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
  end subroutine
end interface

contains

subroutine nvtxStartRange(name,id)
  character(kind=c_char,len=*) :: name
  integer, optional:: id
  type(nvtxEventAttributes):: event
  character(kind=c_char,len=256) :: trimmed_name
  integer:: i

  trimmed_name=trim(name)//c_null_char

  ! move scalar trimmed_name into character array tempName
  do i=1,LEN(trim(name)) + 1
     tempName(i) = trimmed_name(i:i)
  enddo


  if ( .not. present(id)) then
    call nvtxRangePush(tempName)
  else
    event%color=col(mod(id,7)+1)
    event%message=c_loc(tempName)
    call nvtxRangePushEx(event)
  end if
end subroutine

subroutine nvtxEndRange
  call nvtxRangePop
end subroutine

end module nvtx
#endif

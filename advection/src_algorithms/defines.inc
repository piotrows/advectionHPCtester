#ifdef CUDACODE
#define REAL_euwp  REAL(KIND=euwp),MANAGED  
#define REAL_sp    REAL(KIND=sp  ),MANAGED  
#define REAL_dp    REAL(KIND=dp  ),MANAGED  
#define DEV_REAL_euwp  REAL(KIND=euwp  ),DEVICE
#else
#define REAL_euwp  REAL(KIND=euwp)  
#define REAL_sp    REAL(KIND=sp  )  
#define REAL_dp    REAL(KIND=dp  )  
#define DEV_REAL_euwp  REAL(KIND=euwp  )
#endif

#define InnerZFullXYDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
      do k=1+gndedge,lp-skyedge;\
        do j=1,mp;\
          do i=1,np;\
            loop_body\
          end do;\
        end do;\
      end do;                                                                             
#define InnerZFullXYDomainLoopDCs(stream,loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4), stream >>>\
      do k=1+gndedge,lp-skyedge;\
        do j=1,mp;\
          do i=1,np;\
            loop_body\
          end do;\
        end do;\
      end do;                                                                             
#define FullXYDomainLoopDC(loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4) >>>\
      do j=1,mp;\
        do i=1,np;\
          loop_body\
        end do;\
      end do;
#define FullXYDomainLoopDCs(stream,loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4), stream >>>\
      do j=1,mp;\
        do i=1,np;\
          loop_body\
        end do;\
      end do;
#define FullXYZDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
    do k=1,lp;\
      do j=1,mp;\
        do i=1,np;\
          loop_body\
        end do;\
      end do;\
    end do;
#define FullXYZDomainWithHaloLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
    do k=1-ih,lp+ih;\
      do j=1-ih,mp+ih;\
        do i=1-ih,np+ih;\
          loop_body\
        end do;\
      end do;\
    end do;
#define InnerX_CGRID_POLES_FullYZDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
      do k=1,lp;\
        do j=1,mp;\
          do i=1+leftedge*(1-ipoles0), np+(1-rightedge);\
            loop_body\
          end do;\
        end do;\
      end do;
#define CGRIDXFullYZDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
      do k=1,lp;\
        do j=1,mp;\
          do i=1, np+rightedge;\
            loop_body\
          end do;\
        end do;\
      end do;
#define CGRIDInnerXFullYZDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
      do k=1,lp;\
        do j=1,mp;\
          do i=1+leftedge,np;\
            loop_body\
          end do;\
        end do;\
      end do;
#define InnerY_CGRID_FullXZDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
      do k=1,lp;\
        do j=1+botedge,mp+(1-topedge);\
          do i=1,np;\
            loop_body\
          end do;\
        end do;\
      end do;
#define CGRIDYFullXZDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
      do k=1,lp;\
        do j=1,mp+topedge;\
          do i=1,np;\
            loop_body\
          end do;\
        end do;\
      end do;
#define CGRIDInnerYFullXZDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
      do k=1,lp;\
        do j=1+botedge,mp;\
          do i=1,np;\
            loop_body\
          end do;\
        end do;\
      end do;
#define LInnerZFullXYDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
      do k=1+gndedge,lp;\
        do j=1,mp;\
          do i=1,np;\
            loop_body\
          end do;\
        end do;\
      end do;
#define CGRIDZFullXYDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
      do k=1,lp+skyedge;\
        do j=1,mp;\
          do i=1,np;\
            loop_body\
          end do;\
        end do;\
      end do;
#define CGRIDZInnerXYDomainLoopDC(loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4) >>>\
      do k=1+gndedge,lp;\
        do j=1,mp;\
          do i=1,np;\
            loop_body\
          end do;\
        end do;\
      end do;

#define X2DWallFullYZDomainWithHaloLoopDC(loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4) >>>\
      do k=1-ih,lp+ih;\
        do j=1-ih,mp+ih;\
          loop_body\
        end do;\
      end do; 
#define X2DWallFullYZDomainWithHaloLoopDCs(stream,loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4), stream >>>\
      do k=1-ih,lp+ih;\
        do j=1-ih,mp+ih;\
          loop_body\
        end do;\
      end do; 
#define X2DLeftWallFullYZDomainWithHaloLoopDCs(stream,loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4), stream >>>\
      do k=1-ih,lp+ih;\
        do j=1-ih,mp+ih;\
          do i=1,1;\
            loop_body\
          end do;\
        end do;\
      end do; 
#define X2DRightCGRIDWallFullYZDomainWithHaloLoopDCs(stream,loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4), stream >>>\
      do k=1-ih,lp+ih;\
        do j=1-ih,mp+ih;\
          do i=np+1,np+1;\
            loop_body\
          end do;\
        end do;\
      end do; 
#define X2DWallFullYZDomainLoopDC(loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4) >>>\
      do k=1,lp;\
        do j=1,mp;\
          loop_body\
        end do;\
      end do; 
#define X2DWallFullYZDomainLoopDCs(stream,loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4), stream>>>\
      do k=1,lp;\
        do j=1,mp;\
          loop_body\
        end do;\
      end do; 

#define Y2DWallFullXZDomainWithHaloLoopDCs(stream,loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4), stream >>>\
      do k=1-ih,lp+ih;\
        do i=1-ih,np+ih;\
          loop_body\
        end do;\
      end do; 
#define Y2DLeftWallFullXZDomainWithHaloLoopDCs(stream,loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4), stream >>>\
      do k=1-ih,lp+ih;\
        do j=1,1;\
          do i=1-ih,np+ih;\
            loop_body\
          end do;\
        end do;\
      end do; 
#define Y2DRightCGRIDWallFullXZDomainWithHaloLoopDCs(stream,loop_body)\
!$cuf kernel do(3) <<< (*,*), (32,4), stream >>>\
      do k=1-ih,lp+ih;\
        do j=np+1,np+1;\
          do i=1-ih,np+ih;\
            loop_body\
          end do;\
        end do;\
      end do; 
#define Y2DWallFullXZDomainWithHaloLoopDC(loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4) >>>\
      do k=1-ih,lp+ih;\
        do i=1-ih,np+ih;\
          loop_body\
        end do;\
      end do; 
#define Y2DWallFullXZDomainLoopDC(loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4) >>>\
      do k=1,lp;\
        do i=1,np;\
          loop_body\
        end do;\
      end do; 
#define Y2DWallFullXZDomainLoopDCs(stream,loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4), stream >>>\
      do k=1,lp;\
        do i=1,np;\
          loop_body\
        end do;\
      end do; 

#define Z2DWallFullXYDomainWithHaloLoopDC(loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4) >>>\
      do j=1-ih,mp+ih;\
        do i=1-ih,np+ih;\
          loop_body\
        end do;\
      end do; 
#define Z2DWallFullXYDomainLoopDC(loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4) >>>\
      do j=1,mp;\
        do i=1,np;\
          loop_body\
        end do;\
      end do; 
#define Z2DWallFullXYDomainLoopDCs(stream,loop_body)\
!$cuf kernel do(2) <<< (*,*), (32,4),stream >>>\
      do j=1,mp;\
        do i=1,np;\
          loop_body\
        end do;\
      end do; 

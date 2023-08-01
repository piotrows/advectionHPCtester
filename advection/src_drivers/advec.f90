#include "../src_algorithms/defines.inc"
MODULE advec_driver
   USE precisions
#ifdef CUDACODE
   USE cudafor
#endif
   USE mpi_parallel, ONLY: ttbeg,ttend,iup
   USE mpi_parallel, ONLY: mype
   USE module_upwind3d_gpubc,  ONLY:  upwind3d_gpubc
   USE module_antidiff3d_gauge_gpubc,  ONLY:  antidiff3d_gauge_gpubc
!  USE module_antidiff3d_standard_halobc, ONLY:  antidiff3d_standard_halobc
   USE module_antidiff3d_standard_gpubc,  ONLY:  antidiff3d_standard_gpubc
!  USE module_mpdata3d_gauge_halobc_split_compr_driver, ONLY:  mpdata3d_gauge_halobc_split_xy_z


#ifdef PNETCDF
   USE mpi_parallel, ONLY: pnet_out_chunk 
#endif
   IMPLICIT NONE
   INTEGER(KIND=iintegers) icnt,iprint
   CONTAINS

   SUBROUTINE advec_dwarf(legacyoptmode,t_adv,   &
                           lupdatemulti, lvertsplit,         &
                           np,mp,lp,ih,  &
                           ipoles,ibcx,ibcy,ibcz, &
                           x, x_lr, x_bt,  &
                           uadv,vadv,wadv,rhr,rhoadv, &
                           wadv2,rhr2,rhoadv2,rhr3,rhoadv3 )                   
   USE scratch_datafields, ONLY: scratch_zeros,xant
   INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
   INTEGER(KIND=iintegers),INTENT(IN) :: ipoles,ibcx,ibcy,ibcz 
   INTEGER(KIND=iintegers),INTENT(IN) :: t_adv,legacyoptmode
   REAL_euwp,  INTENT(INOUT) ::                            &
                       x(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih)
   REAL_euwp,  INTENT(INOUT) ::    &
                       x_lr(mp,lp,2)
   REAL_euwp,  INTENT(INOUT) ::    &
                       x_bt(np,lp,2)
   REAL_euwp,  INTENT(IN),DIMENSION(1-ih:np+ih,1-ih:mp+ih,1-ih:lp+ih) :: rhr,rhoadv                           
   REAL_euwp,  DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+ih  ),INTENT(IN) :: uadv
   REAL_euwp,  DIMENSION(1-ih:np+ih  , 1-ih:mp+1+ih, 1-ih:lp+ih  ),INTENT(IN) :: vadv
   REAL_euwp,  DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+1+ih) :: wadv
   REAL_euwp,DIMENSION(1-ih:np+1+ih, 1-ih:mp+ih  , 1-ih:lp+1+ih),OPTIONAL :: wadv2,rhr2,rhoadv2,rhr3,rhoadv3
   LOGICAL,INTENT(IN) :: lupdatemulti,lvertsplit
   REAL(KIND=euwp), PARAMETER :: dtloc=0._euwp   
   INTEGER(KIND=iintegers), PARAMETER :: iflip=0
   INTEGER(KIND=iintegers) :: counter_number 
   INTEGER i,j,k,istat
   PROCEDURE(upwind3d_gpubc),POINTER :: upwind3d_ptr
   PROCEDURE(antidiff3d_gauge_gpubc),POINTER :: antidiff3d_ptr

   SELECT CASE (legacyoptmode)
        CASE(4)
           SELECT CASE (t_adv)
                  CASE(1)
                      upwind3d_ptr =>upwind3d_GPUBC
                    antidiff3d_ptr =>antidiff3d_gauge_GPUBC
                    counter_number=55
                  CASE(2)
                      upwind3d_ptr =>upwind3d_GPUBC
                    antidiff3d_ptr =>antidiff3d_standard_GPUBC
                    counter_number=57
                  CASE DEFAULT
                    STOP 'Wrong legacyoptmode in advec_lib'
           END SELECT
        CASE(6)
                      upwind3d_ptr =>upwind3d_GPUBC
        CASE DEFAULT
              STOP 'Wrong legacyoptmode in advection driver'
   END SELECT

     IF(legacyoptmode.eq.5.OR.legacyoptmode.eq.6) THEN       
       CALL upwind3d_ptr(uadv,vadv,wadv,          x,xant,   scratch_zeros,scratch_zeros,       &
                       rhr,rhoadv,  x_lr, x_bt,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,dtloc,  &
                                  .FALSE.,.FALSE.)
!$cuf kernel do(3) <<<*,*>>>
     DO k=1,lp
       DO j=1,mp
         DO i=1,np
       x(i,j,k)=xant(i,j,k)
         ENDDO
       ENDDO
     ENDDO
     ELSE
       CALL ttbeg(counter_number)
       CALL upwind3d_ptr(uadv,vadv,wadv,          x,xant,   scratch_zeros,scratch_zeros,       &
                       rhr,rhoadv,  x_lr, x_bt,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,dtloc,  &
                                  .FALSE.,.FALSE.)
       CALL antidiff3d_ptr(lupdatemulti,uadv,vadv,wadv,xant,x,rhr,rhoadv,  &
                           iflip,ipoles,ibcx,ibcy,ibcz,np,mp,lp,ih,.FALSE.)
       CALL ttend(counter_number)
     ENDIF

     END SUBROUTINE advec_dwarf
END MODULE advec_driver
 

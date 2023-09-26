#include  "../advection/src_algorithms/renames.inc"
MODULE noise
    USE precisions, ONLY  : iintegers,euwp
    IMPLICIT NONE
CONTAINS
    SUBROUTINE perturb_2D_signal(field, mp, np, k, ih)
        INTEGER(KIND=iintegers), intent(in) :: np,mp,ih, k
        REAL(KIND=euwp), DIMENSION(1-ih:np+ih, 1-ih:mp+ih) :: field
        REAL(KIND=euwp) piq
        INTEGER(KIND=iintegers) :: i, j
        DO j=1,mp
            DO i=1,np
!                piq=acos(-1.)/4
                piq=1
                field(i,j)=field(i,j)+ cos(i*piq)*cos(j*piq)*cos(k*piq)*1.e-5
            END DO
        END DO
     END SUBROUTINE perturb_2D_signal
     SUBROUTINE noise_3d_wave(noise_amp,noise_phase,f1,np,mp,lp,ih)
        USE mod_parameters, ONLY: pi2,n,m,l
        USE mpi_parallel, ONLY: nsubpos,msubpos,lsubpos
        INTEGER(KIND=iintegers),INTENT(IN) :: np,mp,lp,ih
        INTEGER(KIND=iintegers) :: i,j,k,ia,ja,ka
        REAL(KIND=euwp) :: noise_amp,noise_phase
        REAL(KIND=euwp),DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: &
         f1 
        !REAL(KIND=euwp),DIMENSION(np,mp,lp) :: noise
!        print *,'noise',noise
!        print *,'nsubpos,msubpos,lsubpos',nsubpos,msubpos,lsubpos
        !Set fluid densities
        DO k=1,lp
          DO j=1,mp
            DO i=1,np
              ia=nsubpos+i
              ja=msubpos+j
              ka=lsubpos+k
                   f1(i,j,k)=noise_amp*sin(pi2*ia/n+        noise_phase*pi2/3._euwp) &
                                      *cos(pi2*ja/m+2._euwp*noise_phase*pi2/3._euwp) &
                                      *cos(pi2*ka/l+3._euwp*noise_phase*pi2/3._euwp)
            ENDDO
          ENDDO
        ENDDO

    END SUBROUTINE noise_3d_wave
 
     SUBROUTINE noise_elg_par_dflt(f1,f2,f3,zcr,ipoles,ibcx,ibcy,ibcz, np, mp, lp, ih, iflg)
     USE mpi_parallel, ONLY: mype
     INTEGER(KIND=iintegers), INTENT(IN) :: iflg,np,mp,lp,ih
        REAL(KIND=euwp), DIMENSION(1-ih:np+ih, 1-ih:mp+ih, 1-ih:lp+ih) :: zcr,f1,f2,f3
        REAL(KIND=euwp) x,randx,randomf1,randomf2,randomf3
        INTEGER(KIND=iintegers) :: i, j, k
        INTEGER(KIND=iintegers) :: ipoles, ibcx,ibcy,ibcz
        INTEGER(KIND=iintegers) :: im, ia, ic, iran
       randomf1(x)=(x-0.5_euwp)
       randomf2(x)=(x-0.5_euwp)
       randomf3(x,i,j,k)=(x-0.5_euwp)*max(0._euwp,1._euwp-zcr(i,j,k)/500._euwp)

!------------------------------------------------------------------
! Numerical Recipes in Fortran - Quick and Dirty Generators p.274-275
!       im         ia         ic           overflow
!      86436       1093      18254           2^27
!     117128       1277      24749           2^28
!     145800       3661      30809           2^29
!     139968       3877      29573           2^30
!     134456       8121      28411           2^31
!     233280       9301      49297           2^32
!------------------------------------------------------------------
      im=86436
      ia=1093
      ic=18254
      iran=mype+1+iflg             ! mype is a processor geometry parameter

        DO k=1,lp-ibcz*ipoles
          DO j=1,mp-ibcy*ipoles
            DO i=1,np-ibcx*ipoles
         iran=mod(iran*ia+ic,im)
         randx=float(iran)/float(im)
         f1(i,j,k)=randomf1(randx)
         iran=mod(iran*ia+ic,im)
         randx=float(iran)/float(im)
         f2(i,j,k)=randomf2(randx)
         iran=mod(iran*ia+ic,im)
         randx=float(iran)/float(im)
         f3(i,j,k)=randomf3(randx,i,j,k)
            END DO
          END DO
        END DO
     END SUBROUTINE noise_elg_par_dflt 
END MODULE noise

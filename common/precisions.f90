#include  "../advection/src_algorithms/renames.inc"
MODULE precisions
   USE, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INTEGER, PARAMETER :: iintegers = KIND  (1)
!  INTEGER, PARAMETER :: sp        = SELECTED_REAL_KIND( 6, 37) !< single precision
!  INTEGER, PARAMETER :: dp        = SELECTED_REAL_KIND(15,307) !< double precision
   INTEGER, PARAMETER :: sp = REAL32
   INTEGER, PARAMETER :: dp = REAL64
#ifdef FLOAT_PRECISION
#if FLOAT_PRECISION==4
   INTEGER, PARAMETER :: euwp = sp
#else 
   INTEGER, PARAMETER :: euwp = dp
#endif
#else
   INTEGER, PARAMETER :: euwp = dp
#endif
END MODULE precisions


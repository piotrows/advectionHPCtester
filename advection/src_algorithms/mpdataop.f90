#include "renames.inc"
MODULE mpdataoperators
   USE precisions, ONLY  : iintegers,euwp
#ifdef CUDACODE
   USE cudafor
   USE mpi_parallel, ONLY: istream1, istream2
#endif
IMPLICIT NONE
CONTAINS
#include"defines.inc"
#include"mpdataoperators.inc"

END MODULE mpdataoperators


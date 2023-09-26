#include  "../advection/src_algorithms/renames.inc"
MODULE epsilons
USE precisions, ONLY : euwp
REAL(KIND=euwp), PARAMETER:: ep_nonos=1.e-10_euwp
REAL(KIND=euwp), PARAMETER:: ep_abs=1.e-10_euwp
CONTAINS
END MODULE epsilons

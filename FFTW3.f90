MODULE FFTW3
      IMPLICIT NONE
      USE, INTRINSIC :: iso_c_binding
     
CONTAINS

      SUBROUTINE ffourier(u,fu)
        USE global_variables
        IMPLICIT NONE
        INCLUDE "fftw3-mpi.f03"



      END SUBROUTINE 

END MODULE

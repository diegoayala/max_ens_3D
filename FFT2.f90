MODULE FFT2
  IMPLICIT NONE  

CONTAINS

  SUBROUTINE ffourier(u,fu)
    !=======================================
    ! Forward Fourier transform
    !=======================================
    USE global_variables
    IMPLICIT NONE 
    INCLUDE "mpif.h"
  
    COMPLEX(pr), DIMENSION (1:n(1),1:n(2),1:local_nlast), INTENT(IN) :: u
    COMPLEX(pr), DIMENSION (1:n(1),1:n(2),1:local_nlast), INTENT(OUT) :: fu

    fu = u
    CALL fftwnd_fortran_mpi(fwdplan, 1, fu, fwork, FFTW_NORMAL_ORDER)
    fu = fu/SQRT(PRODUCT(REAL(n,pr))) 
  
  END SUBROUTINE ffourier

  SUBROUTINE bfourier(fu,u)
    !=======================================
    ! Inverse Fourier transform
    !=======================================
    USE global_variables
    IMPLICIT NONE 
    
    INCLUDE "mpif.h"

    COMPLEX(pr), DIMENSION (1:n(1),1:n(2),1:local_nlast), INTENT(OUT) :: u
    COMPLEX(pr), DIMENSION (1:n(1),1:n(2),1:local_nlast), INTENT(IN) :: fu
    COMPLEX(pr), DIMENSION (1:n(1),1:n(2),1:local_nlast) :: ftmp
        
    ftmp = fu
    CALL fftwnd_fortran_mpi(bwdplan, 1, ftmp, fwork, FFTW_NORMAL_ORDER)
    u = ftmp/SQRT(PRODUCT(REAL(n,pr)))

  END SUBROUTINE bfourier

!  SUBROUTINE ffourier_Fixres(u,fu)
!    !=======================================
!    ! Forward Fourier transform
!    !=======================================
!    USE global_variables
!    IMPLICIT NONE 
!
!    INCLUDE "mpif.h"
!  
!    COMPLEX(pr), DIMENSION (:,:), INTENT(in) :: u
!    COMPLEX(pr), DIMENSION (:,:), INTENT(out) :: fu
!    INTEGER :: i1,i2
!
!    fu = u
!    CALL fftwnd_fortran_mpi(fwdplan_Fixres, 1, fu, fwork_Fixres, FFTW_NORMAL_ORDER)
!    fu = fu/2048 
!  
!  END SUBROUTINE ffourier_Fixres
!
!  SUBROUTINE bfourier_Fixres(fu,u)
!    !=======================================
!    ! Inverse Fourier transform
!    !=======================================
!    USE global_variables
!    IMPLICIT NONE 
!    
!    INCLUDE "mpif.h"
!
!    COMPLEX(pr), DIMENSION (:,:), INTENT(out) :: u
!    COMPLEX(pr), DIMENSION (:,:), INTENT(in) :: fu
!    COMPLEX(pr), DIMENSION (1:2048, 1:local_nlastFixres) :: ftmp
!    INTEGER :: i1,i2
!        
!    ftmp = fu
!
!    CALL fftwnd_fortran_mpi(bwdplan_Fixres, 1, ftmp, fwork_Fixres, FFTW_NORMAL_ORDER)
!    u = ftmp/2048
!
!  END SUBROUTINE bfourier_Fixres
!
!
!  SUBROUTINE ffourier2(u,fu)
!    !=======================================
!    ! 2D Forward Fourier transform
!    !=======================================
!    USE global_variables
!    IMPLICIT NONE 
!
!    INCLUDE "mpif.h"
!  
!    COMPLEX(pr), DIMENSION (:,:,:), INTENT(in) :: u
!    COMPLEX(pr), DIMENSION (:,:,:), INTENT(out) :: fu
!    COMPLEX(pr), DIMENSION (1:n(1), 1:local_nlast) :: ftmp
!    INTEGER :: i1,i2
!
!    ftmp = u(:,:,1)
!    CALL fftwnd_fortran_mpi(fwdplan, 1, ftmp, fwork, FFTW_NORMAL_ORDER)
!    fu(:,:,1) = ftmp/SQRT(PRODUCT(REAL(n,pr)))
!
!    ftmp = u(:,:,2) 
!    CALL fftwnd_fortran_mpi(fwdplan, 1, ftmp, fwork, FFTW_NORMAL_ORDER)
!    fu(:,:,2) = ftmp/SQRT(PRODUCT(REAL(n,pr)))
! 
!  END SUBROUTINE ffourier2
!
!  SUBROUTINE bfourier2(fu,u)
!    !=======================================
!    ! 2D Inverse Fourier transform
!    !=======================================
!    USE global_variables
!    IMPLICIT NONE 
!    
!    INCLUDE "mpif.h"
!
!    COMPLEX(pr), DIMENSION (:,:,:), INTENT(out) :: u
!    COMPLEX(pr), DIMENSION (:,:,:), INTENT(in) :: fu
!    COMPLEX(pr), DIMENSION (1:n(1), 1:local_nlast) :: ftmp
!    INTEGER :: i1,i2
!        
!    ftmp = fu(:,:,1)
!    CALL fftwnd_fortran_mpi(bwdplan, 1, ftmp, fwork, FFTW_NORMAL_ORDER)
!    u(:,:,1) = ftmp/SQRT(PRODUCT(REAL(n,pr)))
!
!    ftmp = fu(:,:,2)
!    CALL fftwnd_fortran_mpi(bwdplan, 1, ftmp, fwork, FFTW_NORMAL_ORDER)
!    u(:,:,2) = ftmp/SQRT(PRODUCT(REAL(n,pr)))
!
!  END SUBROUTINE bfourier2
!
END MODULE

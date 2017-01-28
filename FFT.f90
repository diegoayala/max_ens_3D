module FFT
  implicit none  

  contains
	SUBROUTINE ffourier(u,fu)
    	!=======================================
    	! Forward Fourier transform
    	!=======================================
    	USE global_variables
    	IMPLICIT NONE 

    	INCLUDE "mpif.h"
  
     	real(pr), dimension (0:n(1)-1,local_last_start:local_last_start+local_nlast-1), intent(in) :: u
    	complex(pr), dimension (0:n(1)/2,local_last_start:local_last_start+local_nlast-1), intent(out) :: fu
    	integer :: i1,i2

    	DO i2=local_last_start,local_last_start+local_nlast-1
		DO i1 = 0,n(1)-1,2
			fu(i1/2,i2) = CMPLX(u(i1,i2),u(i1+1,i2))/SQRT (REAL(PRODUCT(n))) 
        	END DO
       	END DO

        CALL rfftwnd_fortran_mpi(fplan, 1, fu(:,:), fwork, FFTW_NORMAL_ORDER)

	END SUBROUTINE ffourier

	SUBROUTINE bfourier(fu,u)
    	!=======================================
    	! Inverse Fourier transform
    	!=======================================
    	USE global_variables
    	IMPLICIT NONE 

    	INCLUDE "mpif.h"

     	real(pr), dimension (0:n(1)-1,local_last_start:local_last_start+local_nlast-1), intent(out) :: u
    	complex(pr), dimension (0:n(1)/2,local_last_start:local_last_start+local_nlast-1), intent(in) :: fu
    	complex(pr), dimension (:,:), allocatable :: ftmp
    	integer :: i1,i2

	allocate (ftmp(0:n(1)/2, local_last_start:local_last_start+local_nlast-1))    
        ftmp = fu(:,:)

	CALL rfftwnd_fortran_mpi(bplan, 1, ftmp, fwork, FFTW_NORMAL_ORDER)

	DO i2=local_last_start,local_last_start+local_nlast-1
       		DO i1 = 0,n(1)-1,2
       			u(i1,i2) = REAL(ftmp(i1/2,i2))/SQRT (REAL(PRODUCT(n))) 
       			u(i1+1,i2) = AIMAG(ftmp(i1/2,i2))/SQRT (REAL(PRODUCT(n))) 
      		END DO
    	END DO
  	END SUBROUTINE bfourier
END MODULE

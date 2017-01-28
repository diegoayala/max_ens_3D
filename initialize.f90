SUBROUTINE initialize

  USE global_variables
  IMPLICIT NONE
  INCLUDE "mpif.h"

  INTEGER :: i,j,k

  !CALL rfftwnd_fortran_mpi_create_plan(fwdplan, MPI_COMM_WORLD, 3, n, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE)
  !CALL rfftwnd_fortran_mpi_create_plan(bwdplan, MPI_COMM_WORLD, 3, n, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE)
  !CALL rfftwnd_fortran_mpi_local_sizes(fwdplan, local_nlast, local_last_start,&
  !     local_nlast_after_trans, local_last_start_after_trans,total_local_size)

  CALL fftw3d_fortran_mpi_create_plan(fwdplan, MPI_COMM_WORLD, n(1), n(2), n(3), FFTW_FORWARD, FFTW_ESTIMATE)
  CALL fftw3d_fortran_mpi_create_plan(bwdplan, MPI_COMM_WORLD, n(1), n(2), n(3), FFTW_BACKWARD, FFTW_ESTIMATE)
  CALL fftwnd_fortran_mpi_local_sizes(fwdplan, local_nlast, local_last_start,&
       local_nlast_after_trans, local_last_start_after_trans, total_local_size)

  ALLOCATE( Uvec(1:n(1),1:n(2),1:local_nlast,1:3) )
  ALLOCATE( Wvec(1:n(1),1:n(2),1:local_nlast,1:3) )

  ALLOCATE( K1(1:n(1)) )
  ALLOCATE( K2(1:n(2)) )
  ALLOCATE( K3(1:n(3)) )
  ALLOCATE( fwork(0:total_local_size) )

  n_dim = 3*n(1)*n(2)*local_nlast
  dV = 1.0_pr/PRODUCT(REAL(n,pr))
  Kcut = 2.0_pr*PI*REAL(n(1),pr)/3.0_pr
  Kmax = Kcut/2.0_pr

  !--Set up wavenumbers
  DO i = 0, n(1)-1
    IF (i<=n(1)/2) THEN
       K1(i+1) = 2.0_pr*PI*REAL(i,pr)
    ELSE
       K1(i+1) = 2.0_pr*PI*REAL(i-n(1),pr)
    END IF
  END DO

  DO i = 0,n(2)-1
    IF (i <= n(2)/2) THEN
      K2(i+1) = 2.0_pr*PI*REAL(i,pr)
    ELSE
      K2(i+1) = 2.0_pr*PI*REAL(i-n(2),pr)
    END IF
  END DO

  DO i = 0, n(3)-1
    IF (i<=n(3)/2) THEN
       K3(i+1) = 2.0_pr*PI*REAL(i,pr)
    ELSE
       K3(i+1) = 2.0_pr*PI*REAL(i-n(3),pr)
    END IF
  END DO

  kappaTest = .FALSE.
  toDealias = .TRUE.
  add_pert = .FALSE.
  save_diag_NS = .FALSE.
  save_data_NS = .FALSE.
  save_diag_lineMin = .FALSE.
  save_data_lineMin = .FALSE.
  save_diag_Constr = .TRUE.
  save_data_Constr = .FALSE.
  save_diag_Optim = .TRUE.
  save_data_Optim = .TRUE.
 
  IF (n(1)<256) THEN
     parallel_data = .FALSE.
  ELSE
     parallel_data = .TRUE.
  END IF
 
END SUBROUTINE



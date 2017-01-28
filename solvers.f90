!======================================
! MODULE CONTAINING INTERFACES FOR 
! NAVIER-STOKES DIRECT SOLVER AND
! ADJOINT SOLVER
!
! (*) FWD_NS2D
! (*) ADJ_NS2D
! (*) D1_NS2D
! (*) D2_NS2D
! (*) EXACT_NS2D 
!
!======================================
MODULE solvers

  IMPLICIT NONE

  CONTAINS
     !==================================================
     ! Numerical solution of 2D Navier Stokes System
     ! using Krylov subspace method
     !================================================== 
     SUBROUTINE fwd_NS2D(U0, mysystem)
       USE global_variables
       USE krylov_routines
       USE data_ops
       IMPLICIT NONE
       
       INCLUDE "mpif.h" 

       REAL(pr), DIMENSION(1:n(1),1:local_nlast,1:2), INTENT(IN) :: U0
       CHARACTER(len=*), INTENT(IN) :: mysystem

       REAL(pr), DIMENSION(1:n_dim) :: u_kry       
       INTEGER, DIMENSION(8) :: t_ini, t_fin

       IF (rank==0 .AND. timing) THEN
          CALL date_and_time(VALUES=t_ini)
       END IF

       IF (mysystem == "nse") save_diag_NS = .TRUE.
 
       CALL matvec(U0, u_kry)
       IF (rank==0) PRINT *, " Solving direct problem..."
       CALL krylov_fwd(u_kry, mysystem)
       IF (rank==0) PRINT *, " Direct problem solved."

       IF (rank==0 .AND. timing) THEN
          CALL date_and_time(VALUES=t_fin)
          IF (Tcount==1 .AND. Pcount==1) THEN
             CALL save_elapsed_time("nse", 0, t_ini, t_fin)
          ELSE
             CALL save_elapsed_time("nse", 1, t_ini, t_fin)
          END IF
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)

       save_diag_NS = .FALSE.

     END SUBROUTINE fwd_NS2D

     !==================================================
     ! Numerical solution of 2D Adjoint System
     ! using Krylov subspace method
     !==================================================               
     SUBROUTINE adj_NS2D
       USE global_variables
       USE krylov_routines
       USE function_ops
       USE data_ops
       IMPLICIT NONE
       INCLUDE "mpif.h"

       REAL(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: u0
       REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: aux, w0
       REAL(pr), DIMENSION(1:n_dim) :: u_kry
       INTEGER :: ii, myindex
       REAL(pr) :: init_time
       LOGICAL :: FromFile

       CHARACTER(len = 100) :: filename
       CHARACTER(2) :: rang
       CHARACTER(4) :: resol
       CHARACTER(4) :: indexchar

       FromFile = .FALSE.
       myindex = 0
       WRITE(rang, '(i2.2)') rank
       WRITE(resol, '(i4.4)') n(1)
       WRITE(indexchar, '(i4.4)') myindex
 
       IF (FromFile) THEN
         IF (rank==0) THEN
            CALL read_init_time(init_time, myindex+1, "adj")
         END IF
         CALL MPI_BCAST(init_time, 1, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)
         filename = "/scratch/ayalada/MAXPALINS/DataGenerated/N_"//resol//"/Adj_U/U_"//indexchar//"_"//rang//".bin" 
         CALL read_field(u0, filename)
       ELSE
         init_time = T
         CALL read_velocity(u0, T, "nse")
         CALL bilaplacian(u0)
       END IF
 
       CALL matvec(u0, u_kry)
       
       IF (rank==0) PRINT *, " Solving adjoint problem..."
       CALL krylov_bwd2(u_kry, T - init_time, myindex)
       IF (rank==0) PRINT *, " Adjoint problem solved."
    
     END SUBROUTINE adj_NS2D

     !==================================================
     ! Numerical solution of Linearized 
     ! 2D Navier Stokes System
     ! using Krylov subspace method
     !================================================== 
     SUBROUTINE d1_NS2D(U0)
       USE global_variables
       USE krylov_routines
       USE data_ops
       IMPLICIT NONE
       
       INCLUDE "mpif.h" 

       REAL(pr), DIMENSION(:,:,:), INTENT(IN) :: U0

       REAL(pr), DIMENSION(1:n_dim) :: u_kry       

 
       CALL matvec(U0, u_kry)
       CALL krylov_linear_NS(u_kry)


     END SUBROUTINE d1_NS2D



     !==================================================
     ! Solution of the 2D heat equation
     !==================================================
     SUBROUTINE exact_NS2D(u0, J)
       USE global_variables
       USE FFT2
       USE function_ops
       IMPLICIT NONE

       REAL(pr), DIMENSION(:,:,:), INTENT(IN) :: u0
       REAL(pr), INTENT(OUT) :: J

       REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: w0
       COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: aux, faux
       INTEGER :: i1, i2
       REAL(pr) :: ksq
      
       CALL cal_vort(u0, w0)
       aux = CMPLX(w0, 0.0_pr)
       CALL ffourier(aux, faux)
       
       DO i2=1,local_nlast
         DO i1=1,n(1)
           ksq = K1(i1)**2 + K2(i2+local_last_start)**2
           faux(i1,i2) = faux(i1,i2)*EXP(-visc*ksq*T)
         END DO
       END DO
  
       CALL bfourier(faux, aux)
       w0 = REAL(aux)

       J=1

     END SUBROUTINE exact_NS2D

END MODULE


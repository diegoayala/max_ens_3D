MODULE optimization

  IMPLICIT NONE

  CONTAINS
    SUBROUTINE maxdEdt
      USE global_variables
      USE FFT2
      USE data_ops
      USE function_ops
      IMPLICIT NONE
      INCLUDE "mpif.h"

      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: gradJ0, gradJ1, diff_gradJ
      REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: f_scalar

      REAL(pr) :: J0, J1, deltaJ, ell, tau, beta, local_scalar_L2norm, divU_L2, divGradJ_L2
      REAL(pr), DIMENSION(1:3) :: local_field_L2norm, local_gradJK0, local_gradJK1, K, E, gradJ_K0, gradJ_K1
      REAL(pr), DIMENSION(1:2) :: tau_brack, dEdt

      INTEGER :: iter, gradType, mnbrak_flag, FixConstr_flag
      LOGICAL :: iterMAX = .FALSE.

      ALLOCATE( gradJ0(1:n(1),1:n(2),1:local_nlast,1:3) )
      ALLOCATE( gradJ1(1:n(1),1:n(2),1:local_nlast,1:3) )
      ALLOCATE( diff_gradJ(1:n(1),1:n(2),1:local_nlast,1:3) )
      ALLOCATE( f_scalar(1:n(1),1:n(2),1:local_nlast) )

      ell = lambda0

      CALL Fix_Constr(Uvec, FixConstr_flag)
      IF (FixConstr_flag /= 0) THEN
         CALL optim_msg_handle(13)
         RETURN
      END IF
 
      !====================================
      ! CALCULATE DIAGNOSTICS OF CONTROL
      !====================================
      local_field_L2norm = Energy(Uvec)
      CALL MPI_ALLREDUCE(local_field_L2norm, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      CALL divergence(Uvec, f_scalar)
      local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
      CALL MPI_ALLREDUCE(local_scalar_L2norm, divU_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

      local_field_L2norm = Enstrophy(Uvec)
      CALL MPI_ALLREDUCE(local_field_L2norm, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      dEdt = calc_dEdt(Uvec)

      J0 = eval_J(Uvec, "maxdEdt", "K0E0")
      J1 = 1.5_pr*J0
      deltaJ = ABS( (J1-J0)/J0 )    
      iter = 0

      IF (save_diag_Optim) THEN
         IF (rank==0) THEN
            CALL save_diagnostics_optim("maxdEdt", iter, 0.0_pr, 0.0_pr, J0, K, E, divU_L2, dEdt(1), dEdt(2))   
         END IF
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END IF
  
      IF (save_data_Optim) THEN
         CALL vel2vort(Uvec, Wvec)
         CALL diagnosticScalars(Uvec, Wvec, iter)
         CALL save_Ctrl(Uvec, Wvec, iter, "maxdEdt")
      END IF
 
      DO WHILE ( (ABS(deltaJ) > OPTIM_TOL) .AND. (iter<MAX_ITER) )
         !======================================================
         ! CALCULATE BASIC GRADIENT IN THE H^s TOPOLOGY
         !======================================================
         CALL eval_grad_J(Uvec, gradJ1, iter, "maxdEdt")
         IF (kappaTest) CALL kappa_test(Uvec, gradJ1, J0, "maxdEdt")
         kappaTest = .FALSE.
         CALL SobolevGradient(gradJ1, 1, ell)

         !======================================================
         ! DECIDE DIRECTION OF INCREASE
         !====================================================== 
         IF (MOD(iter,20)==0) THEN
            gradType = 0
         ELSE
            gradType = 0
         END IF
         
         SELECT CASE (gradType)
            CASE (0) !REGULAR STEEPEST ASCENT METHOD
               local_field_L2norm = 2.0_pr*Energy(gradJ1)
               CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K1, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
               beta = 0.0_pr
               gradJ_K0 = gradJ_K1
               diff_gradJ = gradJ1

            CASE (1) !FLETCHER-REEVES METHOD
               local_field_L2norm = 2.0_pr*Energy(gradJ1)
               CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K1, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
               beta = SUM(gradJ_K1)/SUM(gradJ_K0)
               gradJ_K0 = gradJ_K1
               diff_gradJ = gradJ1
               gradJ1 = gradJ1 + beta*gradJ0
 
            CASE (2) !POLAK-RIBIERE METHOD
               diff_gradJ = gradJ1 - diff_gradJ
               local_field_L2norm = field_inner_product(gradJ1, diff_gradJ, "L2")
               CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K1, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
               beta = SUM(gradJ_K1)/SUM(gradJ_K0)
               local_field_L2norm = 2.0_pr*Energy(gradJ1)
               CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K0, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
               diff_gradJ = gradJ1
               gradJ1 = gradJ1 + beta*gradJ0

         END SELECT

         !===========================================
         ! GET DIAGNOSTICS OF ASCENT DIRECTION AND SAVE
         !===========================================
         local_field_L2norm = Energy(gradJ1)
         CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K1, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
         CALL divergence(gradJ1, f_scalar)
         local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
         CALL MPI_ALLREDUCE(local_scalar_L2norm, divGradJ_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

         CALL save_gradient(gradJ1, iter, "maxdEdt")
       
         !======================================
         ! FIND OPTIMAL tau BY ARC OPTIMIZATION
         !====================================== 
         IF (iter==0) THEN 
            tau = alpha0
         END IF
         tau_brack(1) = 0.0_pr
 
         CALL optim_msg_handle(20) 
         tau_brack = mnbrak("maxdEdt", Uvec, gradJ1, tau_brack(1), tau, mnbrak_flag)
         IF (mnbrak_flag /= 0) THEN
            CALL optim_error_handle(mnbrak_flag)            
            IF (save_diag_Optim) THEN
               IF (rank==0) THEN
                  CALL save_diagnostics_optim("maxdEdt", iter+1, tau, beta, J0, K, E, divU_L2, dEdt(1), dEdt(2))   
               END IF
               CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
            END IF
            RETURN
         ELSE
            CALL optim_msg_handle(21)
         END IF

         CALL optim_msg_handle(30)         
         tau = brent("maxdEdt", Uvec, gradJ1, tau_brack)
         CALL optim_msg_handle(31)
         tau = MIN(tau, TAU_MAX)

         !======================================
         ! UPDATE CONTROL VARIABLE
         !======================================
         IF (tau == TAU_MAX) THEN
            CALL optim_msg_handle(32)
         END IF

         Uvec = Uvec + tau*gradJ1 
         CALL Fix_Constr(Uvec, FixConstr_flag)
         IF (FixConstr_flag /= 0) THEN
            CALL optim_msg_handle(13)
            RETURN
         END IF
 
         !====================================
         ! CALCULATE DIAGNOSTICS OF CONTROL
         !====================================
         local_field_L2norm = Energy(Uvec)
         CALL MPI_ALLREDUCE(local_field_L2norm, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

         CALL divergence(Uvec, f_scalar)
         local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
         CALL MPI_ALLREDUCE(local_scalar_L2norm, divU_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

         local_field_L2norm = Enstrophy(Uvec)
         CALL MPI_ALLREDUCE(local_field_L2norm, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
 
         dEdt = calc_dEdt(Uvec)

         J1 = eval_J(Uvec, "maxdEdt", "K0E0")
         deltaJ = (J1-J0)/ABS(J0)
         IF (deltaJ < 0.0_pr) THEN
            IF (save_diag_Optim) THEN
               IF (rank==0) THEN
                  CALL save_diagnostics_optim("maxdEdt", iter+1, tau, beta, J0, K, E, divU_L2, dEdt(1), dEdt(2))   
               END IF
               CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
            END IF
            CALL optim_msg_handle(0)
            save_data_optim = .FALSE.
            EXIT
         ELSE
            save_data_optim = .TRUE.
         END IF
 
         !===============================
         ! UPDATE VARIABLES
         !===============================     
         J0 = J1
         gradJ0 = gradJ1
         iter = iter + 1

         IF (save_diag_Optim) THEN
            IF (rank==0) THEN
               CALL save_diagnostics_optim("maxdEdt", iter, tau, beta, J1, K, E, divU_L2, dEdt(1), dEdt(2))   
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         END IF
  
         IF (save_data_Optim) THEN
            CALL vel2vort(Uvec, Wvec)
            CALL diagnosticScalars(Uvec, Wvec, iter)
            CALL save_Ctrl(Uvec, Wvec, iter, "maxdEdt")
         END IF
 
      END DO

      CALL optim_msg_handle(1)

      DEALLOCATE(gradJ0)
      DEALLOCATE(gradJ1)
      DEALLOCATE(diff_gradJ)
      DEALLOCATE(f_scalar)
 
    END SUBROUTINE maxdEdt

    !================================================= 
    ! SUBROUTINE THAT PROJECTS A GIVEN FIELD ONTO THE 
    ! CONSTRAINT MANIFOLD (CO-DIMENSION 1 or 2)
    !=================================================
    SUBROUTINE Fix_Constr(myfield, myflag)
      USE global_variables
      USE function_ops
      USE data_ops
      IMPLICIT NONE
      INCLUDE "mpif.h"

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: myfield
      INTEGER, INTENT(OUT) :: myflag

      INTEGER :: optim_msg_flag

      SELECT CASE (ConsType)
        CASE (1)
           CALL div_free(myfield)
           CALL Fix_E0H0(myfield)
           myflag = 0 
       
        CASE (2)
           CALL optim_msg_handle(10)
           CALL Fix_K0E0(myfield, optim_msg_flag)
           CALL optim_msg_handle(optim_msg_flag)
           IF (optim_msg_flag == 11) THEN
              myflag = 0
           ELSE
              myflag = 1
           END IF

      END SELECT

    END SUBROUTINE Fix_Constr

    !=================================================
    ! SUBROUTINE THAT FIXES HELICITY AND ENSTROPHY
    ! FOR A GIVEN VELOCITY FIELD
    !=================================================
    SUBROUTINE Fix_E0H0(U, FixH0E0_flag)
      USE global_variables
      USE data_ops
      USE FFT2
      USE function_ops
      IMPLICIT NONE
      INCLUDE "mpif.h"

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: U
      INTEGER, INTENT(OUT) :: FixK0E0_flag

      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: W, gradJ
      REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: f_scalar

      REAL(pr) :: J, J1, ell, tau, beta, local_scalar_L2norm, divU_L2, divGradJ_L2
      REAL(pr), DIMENSION(1:3) :: local_field_L2norm, K, E, gradJ_K
      REAL(pr), DIMENSION(1:2) :: tau_brack

      INTEGER :: iter, gradType, mnbrak_flag
      LOGICAL :: iterMAX

      ALLOCATE( W(1:n(1),1:n(2),1:local_nlast,1:3) )
      ALLOCATE( gradJ(1:n(1),1:n(2),1:local_nlast,1:3) )
      ALLOCATE( f_scalar(1:n(1),1:n(2),1:local_nlast) )

      ell = 1.0e-1
      beta = 0.0_pr
      FixH0E0_flag = 12

      CALL div_free(U)
      CALL Fix_E0(U)

      !====================================
      ! CALCULATE DIAGNOSTICS OF CONTROL
      !====================================
      local_field_L2norm = Energy(U)
      CALL MPI_ALLREDUCE(local_field_L2norm, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      CALL divergence(U, f_scalar)
      local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
      CALL MPI_ALLREDUCE(local_scalar_L2norm, divU_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

      local_field_L2norm = Enstrophy(U)
      CALL MPI_ALLREDUCE(local_field_L2norm, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      J = eval_J(U, "FixH0E0", "Enstrophy")     
      iter = 0

      IF (save_diag_Constr) THEN
         IF (rank==0) THEN
            CALL save_diagnostics_optim("FixK0E0", iter, 0.0_pr, ell, J, K, E, divU_L2, SUM(gradJ_K), divGradJ_L2)   
         END IF
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END IF
      
      IF (save_data_Constr) THEN
         CALL vel2vort(U, W)
         CALL save_Ctrl(U, W, iter, "FixK0E0")
      END IF
       
      DO WHILE ( (J > CONSTR_TOL) .AND. (iter<MAX_ITER_CONSTR) )
         !======================================================
         ! CALCULATE BASIC GRADIENT IN THE H^s TOPOLOGY
         !======================================================
         CALL eval_grad_J(U, gradJ, iter, "Enstrophy")
         IF (kappaTest) CALL kappa_test(U, gradJ, J, "FixK0E0")
         kappaTest = .FALSE.
         CALL SobolevGradient(gradJ, 1, ell)

         !===============================
         ! GET DIAGNOSTICS OF GRADIENT
         !===============================
         local_field_L2norm = 0.5_pr*field_inner_product(gradJ, gradJ, "L2")
         CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
         CALL divergence(gradJ, f_scalar)
         local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
         CALL MPI_ALLREDUCE(local_scalar_L2norm, divGradJ_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

         IF (save_data_Constr) THEN
            CALL save_gradient(gradJ, iter, "FixK0E0")
         END IF
                 
         !======================================
         ! FIND OPTIMAL tau BY ARC OPTIMIZATION
         !====================================== 
         IF (iter==0) THEN 
           tau = REAL(1.0, pr) 
         END IF
         tau_brack(1) = 0.0_pr
 
         mnbrak_flag = 1
         DO WHILE ( mnbrak_flag /= 0 .AND. ell > 1.0_pr/(2.0_pr*REAL(n(1),pr)) )
            tau_brack = mnbrak_FixK0E0(U, gradJ, tau_brack(1), tau, mnbrak_flag)         
            IF (mnbrak_flag /= 0) THEN
               CALL optim_error_handle(mnbrak_flag)                        
               ell = 0.95_pr*ell
               CALL eval_grad_J(U, gradJ, iter, "Enstrophy")
               CALL SobolevGradient(gradJ, 1, ell)
            ELSE
               EXIT
            END IF
         END DO
         
         IF (mnbrak_flag /= 0) THEN
            CALL optim_error_handle(15)
            tau = TAU_MAX
         ELSE
            tau = brent_FixK0E0(U, gradJ, tau_brack)
         END IF
         tau = MIN(tau,TAU_MAX)

         !======================================
         ! UPDATE CONTROL VARIABLE
         !====================================== 
         U = U + tau*gradJ
         CALL div_free(U) 
         CALL Fix_K0(U)

         !====================================
         ! CALCULATE DIAGNOSTICS OF CONTROL
         !====================================
         local_field_L2norm = Energy(U)
         CALL MPI_ALLREDUCE(local_field_L2norm, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

         CALL divergence(U, f_scalar)
         local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
         CALL MPI_ALLREDUCE(local_scalar_L2norm, divU_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

         local_field_L2norm = Enstrophy(U)
         CALL MPI_ALLREDUCE(local_field_L2norm, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
 
         J1 = eval_J(U, "FixK0E0", "Enstrophy")     
         iter = iter + 1

         IF (save_diag_Constr) THEN
            IF (rank==0) THEN
               CALL save_diagnostics_optim("FixK0E0", iter, tau, ell, J1, K, E, divU_L2, SUM(gradJ_K), divGradJ_L2)   
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         END IF

         IF (save_data_Constr) THEN
            CALL vel2vort(U, W)
            CALL save_Ctrl(U, W, iter, "FixK0E0")
         END IF

         J = J1
 
      END DO

      IF ( iter .GE. MAX_ITER_CONSTR ) THEN
         FixK0E0_flag = 12
      ELSE
         FixK0E0_flag = 11
      END IF

      DEALLOCATE(W)
      DEALLOCATE(gradJ)
      DEALLOCATE(f_scalar)
      
    END SUBROUTINE Fix_K0E0


    !=================================================
    ! SUBROUTINE THAT FIXES ENERGY AND ENSTROPHY
    ! FOR A GIVEN VELOCITY FIELD
    !=================================================
    SUBROUTINE Fix_K0E0(U, FixK0E0_flag)
      USE global_variables
      USE data_ops
      USE FFT2
      USE function_ops
      IMPLICIT NONE
      INCLUDE "mpif.h"

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: U
      INTEGER, INTENT(OUT) :: FixK0E0_flag

      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: W, gradJ
      REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: f_scalar

      REAL(pr) :: J, J1, ell, tau, beta, local_scalar_L2norm, divU_L2, divGradJ_L2
      REAL(pr), DIMENSION(1:3) :: local_field_L2norm, K, E, gradJ_K
      REAL(pr), DIMENSION(1:2) :: tau_brack

      INTEGER :: iter, gradType, mnbrak_flag
      LOGICAL :: iterMAX

      ALLOCATE( W(1:n(1),1:n(2),1:local_nlast,1:3) )
      ALLOCATE( gradJ(1:n(1),1:n(2),1:local_nlast,1:3) )
      ALLOCATE( f_scalar(1:n(1),1:n(2),1:local_nlast) )

      ell = 1.0e-1
      beta = 0.0_pr
      FixK0E0_flag = 12

      CALL div_free(U)
      CALL Fix_K0(U)

      !====================================
      ! CALCULATE DIAGNOSTICS OF CONTROL
      !====================================
      local_field_L2norm = Energy(U)
      CALL MPI_ALLREDUCE(local_field_L2norm, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      CALL divergence(U, f_scalar)
      local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
      CALL MPI_ALLREDUCE(local_scalar_L2norm, divU_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

      local_field_L2norm = Enstrophy(U)
      CALL MPI_ALLREDUCE(local_field_L2norm, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      J = eval_J(U, "FixK0E0", "Enstrophy")     
      iter = 0

      IF (save_diag_Constr) THEN
         IF (rank==0) THEN
            CALL save_diagnostics_optim("FixK0E0", iter, 0.0_pr, ell, J, K, E, divU_L2, SUM(gradJ_K), divGradJ_L2)   
         END IF
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END IF
      
      IF (save_data_Constr) THEN
         CALL vel2vort(U, W)
         CALL save_Ctrl(U, W, iter, "FixK0E0")
      END IF
       
      DO WHILE ( (J > CONSTR_TOL) .AND. (iter<MAX_ITER_CONSTR) )
         !======================================================
         ! CALCULATE BASIC GRADIENT IN THE H^s TOPOLOGY
         !======================================================
         CALL eval_grad_J(U, gradJ, iter, "Enstrophy")
         IF (kappaTest) CALL kappa_test(U, gradJ, J, "FixK0E0")
         kappaTest = .FALSE.
         CALL SobolevGradient(gradJ, 1, ell)

         !===============================
         ! GET DIAGNOSTICS OF GRADIENT
         !===============================
         local_field_L2norm = 0.5_pr*field_inner_product(gradJ, gradJ, "L2")
         CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
         CALL divergence(gradJ, f_scalar)
         local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
         CALL MPI_ALLREDUCE(local_scalar_L2norm, divGradJ_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

         IF (save_data_Constr) THEN
            CALL save_gradient(gradJ, iter, "FixK0E0")
         END IF
                 
         !======================================
         ! FIND OPTIMAL tau BY ARC OPTIMIZATION
         !====================================== 
         IF (iter==0) THEN 
           tau = REAL(1.0, pr) 
         END IF
         tau_brack(1) = 0.0_pr
 
         mnbrak_flag = 1
         DO WHILE ( mnbrak_flag /= 0 .AND. ell > 1.0_pr/(2.0_pr*REAL(n(1),pr)) )
            tau_brack = mnbrak_FixK0E0(U, gradJ, tau_brack(1), tau, mnbrak_flag)         
            IF (mnbrak_flag /= 0) THEN
               CALL optim_error_handle(mnbrak_flag)                        
               ell = 0.95_pr*ell
               CALL eval_grad_J(U, gradJ, iter, "Enstrophy")
               CALL SobolevGradient(gradJ, 1, ell)
            ELSE
               EXIT
            END IF
         END DO
         
         IF (mnbrak_flag /= 0) THEN
            CALL optim_error_handle(15)
            tau = TAU_MAX
         ELSE
            tau = brent_FixK0E0(U, gradJ, tau_brack)
         END IF
         tau = MIN(tau,TAU_MAX)

         !======================================
         ! UPDATE CONTROL VARIABLE
         !====================================== 
         U = U + tau*gradJ
         CALL div_free(U) 
         CALL Fix_K0(U)

         !====================================
         ! CALCULATE DIAGNOSTICS OF CONTROL
         !====================================
         local_field_L2norm = Energy(U)
         CALL MPI_ALLREDUCE(local_field_L2norm, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

         CALL divergence(U, f_scalar)
         local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
         CALL MPI_ALLREDUCE(local_scalar_L2norm, divU_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

         local_field_L2norm = Enstrophy(U)
         CALL MPI_ALLREDUCE(local_field_L2norm, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
 
         J1 = eval_J(U, "FixK0E0", "Enstrophy")     
         iter = iter + 1

         IF (save_diag_Constr) THEN
            IF (rank==0) THEN
               CALL save_diagnostics_optim("FixK0E0", iter, tau, ell, J1, K, E, divU_L2, SUM(gradJ_K), divGradJ_L2)   
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         END IF

         IF (save_data_Constr) THEN
            CALL vel2vort(U, W)
            CALL save_Ctrl(U, W, iter, "FixK0E0")
         END IF

         J = J1
 
      END DO

      IF ( iter .GE. MAX_ITER_CONSTR ) THEN
         FixK0E0_flag = 12
      ELSE
         FixK0E0_flag = 11
      END IF

      DEALLOCATE(W)
      DEALLOCATE(gradJ)
      DEALLOCATE(f_scalar)
      
    END SUBROUTINE Fix_K0E0

    !=================================================
    ! FUNCTION THAT CALCULATES COST FUNCTIONAL
    !=================================================
    RECURSIVE FUNCTION eval_J(myfield, mysystem, myJ) RESULT (J)
      USE global_variables
      USE data_ops
      USE FFT2
      USE function_ops
      IMPLICIT NONE
      INCLUDE "mpif.h"

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: myfield
      CHARACTER(len=*), INTENT(IN) :: mysystem
      CHARACTER(len=*), INTENT(IN) :: myJ
      REAL(pr) :: J

      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi, aux1, aux2, aux3
      REAL(pr), DIMENSION(1:3) :: local_VecField_norm, global_VecField_norm
      
      INTEGER :: constr_flag

      J = 0.0_pr

      SELECT CASE (mysystem)
        CASE ("maxdEdt")
           ALLOCATE( aux1(1:n(1),1:n(2),1:local_nlast,1:3) )
           ALLOCATE( aux2(1:n(1),1:n(2),1:local_nlast,1:3) )
           ALLOCATE( aux3(1:n(1),1:n(2),1:local_nlast,1:3) )

           aux1 = myfield 
           aux2 = myfield
           CALL div_free(aux1)
           CALL div_free(aux2)
           CALL advection(aux1, aux2, aux3)

           CALL laplacian(aux1)
           local_VecField_norm = Energy(aux1) 
           CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 

           J = -2.0_pr*visc*SUM(global_VecField_norm)

           local_VecField_norm = field_inner_product(aux3,aux1,"L2")           
           CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 

           J = J + SUM(global_VecField_norm)

           DEALLOCATE(aux1)
           DEALLOCATE(aux2)
           DEALLOCATE(aux3)
             
        CASE ("FixK0E0")
           SELECT CASE (myJ)
             CASE ("Enstrophy")
                local_vecField_norm = Enstrophy(myfield)
                CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 

                J = 0.5_pr*( ( SUM(global_VecField_norm) - E0 )/E0 )**2

             CASE ("Energy")
                local_VecField_norm = Energy(myfield) 
                CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 

                J = 0.5_pr*( ( SUM(global_VecField_norm) - K0 )/K0 )**2  

           END SELECT 
 
        CASE ("FixE0H0")
           SELECT CASE (myJ)
             CASE ("Enstrophy")
                local_vecField_norm = Enstrophy(myfield)
                CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 

                J = 0.5_pr*( ( SUM(global_VecField_norm) - E0 )/E0 )**2

             CASE ("Helicity")
                local_VecField_norm = Helicity(myfield) 
                CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 

                J = 0.5_pr*( SUM(global_VecField_norm) - H0 )**2  

           END SELECT 
             
        CASE ("LineMin")

           ALLOCATE( aux1(1:n(1),1:n(2),1:local_nlast,1:3) )
           aux1 = myfield

           SELECT CASE (myJ)
              CASE ("FixK0E0")

                 CALL div_free(aux1)
                 CALL Fix_K0(aux1)
                 J = eval_J(aux1, "FixK0E0", "Enstrophy")
                 
              CASE ("maxdEdt")

                 CALL div_free(aux1)
                 CALL Fix_Constr(aux1, constr_flag)
                 IF (constr_flag == 0) THEN
                    J = -1.0_pr*eval_J(aux1, "maxdEdt", "K0E0")
                 ELSE
                    CALL optim_msg_handle(14)
                    J = -1.0_pr*eval_J(aux1, "maxdEdt", "K0E0")
                 END IF

           END SELECT

           DEALLOCATE(aux1)
  
      END SELECT

    END FUNCTION eval_J
 
    !================================================
    ! Obtain gradient of cost functional
    !================================================
    SUBROUTINE eval_grad_J(phi, gradJ, iter_index, myJ)
      USE global_variables
      USE data_ops
      USE FFT2
      USE function_ops
      IMPLICIT NONE
      INCLUDE "mpif.h"

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: phi
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(OUT) :: gradJ
      INTEGER, INTENT(IN) :: iter_index
      CHARACTER(len=*), INTENT(IN) :: myJ

      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: aux1, aux2, aux3
      REAL(pr), DIMENSION(3) :: local_VecField_norm, global_VecField_norm

      INTEGER :: ii

      SELECT CASE (myJ)
        CASE ("Energy")

        CASE ("Enstrophy")
           DO ii=1,3
              gradJ(:,:,:,ii) = phi(:,:,:,ii)
           END DO
           local_VecField_norm = Enstrophy(phi)
           CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 

           CALL div_free(gradJ)
           CALL laplacian(gradJ)

           gradJ = (SUM(global_VecField_norm) - E0)*gradJ/E0**2

        CASE ("maxdEdt")
           ALLOCATE( aux1(1:n(1),1:n(2),1:local_nlast,1:3) )
           ALLOCATE( aux2(1:n(1),1:n(2),1:local_nlast,1:3) )
           ALLOCATE( aux3(1:n(1),1:n(2),1:local_nlast,1:3) )

           gradJ = phi
           CALL div_free(gradJ)
           aux1 = gradJ
           aux2 = gradJ
 

           CALL bilaplacian(gradJ)
           gradJ = -2.0_pr*visc*gradJ

           CALL advection(aux1,aux2,aux3)
           CALL laplacian(aux3)
           gradJ = gradJ + aux3

           CALL laplacian(aux2)           
           CALL advection(aux1,aux2,aux3)
           gradJ = gradJ - aux3
           
           CALL stretching(aux1,aux2,aux3)
           gradJ = gradJ + aux3

           DEALLOCATE( aux1 )
           DEALLOCATE( aux2 )
           DEALLOCATE( aux3 )

      END SELECT      

    END SUBROUTINE eval_grad_J

    !==================================================
    ! BRACKET THE LOCATION OF OPTIMAL TAU
    !==================================================

    FUNCTION mnbrak(mysystem, phi, grad, tA0, tB0, myflag) RESULT (tau_brack)
      USE global_variables
      USE data_ops
      USE function_ops
      IMPLICIT NONE
  
      CHARACTER(len=*), INTENT(IN) :: mysystem
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: phi
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: grad
      REAL(pr), INTENT(IN) :: tA0, tB0
      INTEGER, INTENT(INOUT) :: myflag
      REAL(pr), DIMENSION(1:2) :: tau_brack

      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi_bar
      REAL(pr) :: aux, tP, FP, Pmax, R, Q
      REAL(pr) :: FA, FB, FC, tA, tB, tC
 
      REAL(pr), PARAMETER :: GOLD = (1.0_pr + SQRT(5.0_pr))/2.0_pr
      REAL(pr), PARAMETER :: CGOLD = 1.0_pr/GOLD  
      REAL(pr), PARAMETER :: GLIMIT = 10.0_pr
      REAL(pr), PARAMETER :: tMAX = 10.0_pr
      INTEGER, PARAMETER :: ITMAX = 100
      INTEGER :: FuncEval, iter 
      LOGICAL :: saveLineMin

      saveLineMin = .TRUE.

      ALLOCATE( phi_bar(1:n(1),1:n(2),1:local_nlast,1:3) )

      FuncEval = 0
      iter = 0
 
      tA = tA0
      tB = MAX(tB0, MACH_EPSILON)

      phi_bar = phi + tA*grad
      FA = eval_J(phi_bar, "LineMin", mysystem)
      FuncEval = FuncEval+1

      phi_bar = phi + tB*grad
      FB = eval_J(phi_bar, "LineMin", mysystem)
      FuncEval = FuncEval+1

      IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, mysystem, "replace")
 
      DO WHILE (FB > FA .AND. tB > MACH_EPSILON) 
         tB = CGOLD*tB
         phi_bar = phi + tB*grad
         FB = eval_J(phi_bar, "LineMin", mysystem)
         FuncEval = FuncEval+1
         IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, mysystem, "append")
      END DO

      IF (tB .LE. MACH_EPSILON) THEN
         myflag = 1
         RETURN
      END IF

      tC = GOLD*tB
      phi_bar = phi + tC*grad
      FC = eval_J(phi_bar, "LineMin", mysystem)
      FuncEval = FuncEval+1

      IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, mysystem, "append")

      DO WHILE (FB>=FC .AND. iter<ITMAX)
         iter = iter+1
         tC = GOLD*tC
         phi_bar = phi + tC*grad
         FC = eval_J(phi_bar, "LineMin", mysystem)
         FuncEval = FuncEval+1
                  
         R = (tB-tA)*(FB-FC)
         Q = (tB-tC)*(FB-FA)
         tP = tB - 0.5_pr*((tB-tC)*Q - (tB-tA)*R)/( SIGN( MAX(ABS(Q-R),MACH_EPSILON), Q-R) )
    
         Pmax = tB + GLIMIT*(tC-tB)
    
         IF ( (tB-tP)*(tP-tC)>0 ) THEN
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", mysystem)
        
            IF (FP<FC) THEN
               tA = tB
               FA = FB
               tB = tP
               FB = FP
               EXIT 
            ELSEIF (FP>FB) THEN
               tC = tP
               FC = FP
               EXIT
            END IF
        
            tP = tC + GOLD*(tC-tB)
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", mysystem)
        
         ELSEIF ( (tC-tP)*(tP-Pmax)>0 ) THEN
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", mysystem)
       
            IF (FP<FC) THEN
               tB = tC
               tC = tP
               FB = FC
               FC = FP
               tP = tC+GOLD*(tC-tB)
               phi_bar = phi + tP*grad
               FP = eval_J(phi_bar, "LineMin", mysystem)
            END IF
        
        ELSEIF ( (tP-Pmax)*(Pmax-tC)>=0 ) THEN
            tP = Pmax
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", mysystem)
        
         ELSE
            tP = tC + GOLD*(tC-tB)
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", mysystem)
         END IF
    
         tA = tB
         tB = tC
         tC = tP
         FA = FB
         FB = FC
         FC = FP
        
         IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, mysystem, "append")
 
      END DO

      tau_brack(1) = tA
      tau_brack(2) = tC

      IF (iter .GE. ITMAX) THEN
         myflag = 2
      ELSE
         myflag = 0
      END IF

      DEALLOCATE( phi_bar )
    

    END FUNCTION mnbrak

  
    !============================================
    ! BRENT ALGORITHM FOR LINE OPTIMIZATION
    !============================================

    FUNCTION brent(mysystem, phi, grad, tau_brack) RESULT (X)
      USE global_variables
      USE data_ops
      USE function_ops
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: mysystem
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: phi, grad
      REAL(pr), DIMENSION(1:2), INTENT(IN) :: tau_brack
      REAL(pr) :: X
   
      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi_bar
      REAL(pr) :: D, A, B, V, W, E, ETEMP, P, Q, R, U, XM
      REAL(pr) :: FV, FW, FU, FX
      REAL(pr) :: TOL1, TOL2
      INTEGER :: FLAG, j


      INTEGER, PARAMETER :: ITMAX = 200
      REAL(pr), PARAMETER :: TOL = 1E-6
      REAL(pr), PARAMETER :: ZEPS = 1E-8
      REAL(pr), PARAMETER :: CGOLD = .381966


      ALLOCATE( phi_bar(1:n(1),1:n(2),1:local_nlast,1:3) )
 
      D = 0.0_pr
      A = MIN(tau_brack(1),tau_brack(2))
      B = MAX(tau_brack(1),tau_brack(2))
      V = tau_brack(2)*CGOLD 
      W = V
      X = V 
      E = 0.0_pr

      phi_bar = phi + X*grad
      FX = eval_J(phi_bar, "LineMin", mysystem)
 
      FV = FX 
      FW = FX

      DO j=1,ITMAX
        XM = 0.5_pr*(A+B)
        TOL1 = TOL*ABS(X)+ZEPS
        TOL2 = 2.0_pr*TOL1
    
        IF ( ABS(X-XM) <= (TOL2-0.5*(B-A)) ) EXIT
    
        FLAG = 1
        IF ( ABS(E) > TOL1 ) THEN
           R = (X-W)*(FX-FV)
           Q = (X-V)*(FX-FW)
           P = (X-V)*Q - (X-W)*R
           Q = 2.0_pr*(Q-R)
           IF ( Q > 0.0_pr ) P = -P 
      
           Q = ABS(Q)
           ETEMP = E
           E = D
        
           IF ( (ABS(P) >= ABS(0.5_pr*Q*ETEMP)) .OR. (P <= Q*(A-X)) .OR. (P >= Q*(B-X)) ) THEN
              FLAG = 1
           ELSE
              FLAG = 2
           END IF

        END IF
    
        SELECT CASE (FLAG)
          CASE (1)
            IF (X >= XM) THEN
                E = A-X
            ELSE
                E=B-X
            END IF
            D = CGOLD*E
          CASE (2)
            D = P/Q
            U = X+D
            IF ( (U-A < TOL2) .OR. (B-U < TOL2) ) D = SIGN(TOL1, XM-X)
            
        END SELECT
    
        IF ( ABS(D) >= TOL1 ) THEN
           U = X+D
        ELSE
           U = X + SIGN(TOL1,D)
        END IF
    
        phi_bar = phi + U*grad
        FU = eval_J(phi_bar, "LineMin", mysystem)
 
        IF ( FU <= FX ) THEN
           IF ( U >= X ) THEN
              A = X
           ELSE
              B = X 
           END IF
           V = W
           FV = FW
           W = X
           FW = FX
           X = U
           FX = FU
        ELSE
           IF ( U < X ) THEN
              A = U
           ELSE
              B = U
           END IF
        
           IF ( (FU <= FW) .OR. (W == X) ) THEN
              V = W
              FV = FW
              W = U
              FW = FU
           ELSEIF ( (FU <= FV) .OR. (V==X) .OR. (V==W)) THEN 
              V = U
              FV = FU
           END IF
        END IF
      END DO

      DEALLOCATE( phi_bar )

    END FUNCTION brent

    !==================================================
    ! BRACKET THE LOCATION OF OPTIMAL TAU
    !==================================================
    FUNCTION mnbrak_FixK0E0(phi, grad, tA0, tB0, myflag) RESULT (tau_brack)
      USE global_variables
      USE data_ops
      USE function_ops
      IMPLICIT NONE
  
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: phi
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: grad
      REAL(pr), INTENT(IN) :: tA0, tB0
      INTEGER, INTENT(INOUT) :: myflag
      REAL(pr), DIMENSION(1:2) :: tau_brack

      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi_bar
      REAL(pr) :: aux, tP, FP, Pmax, R, Q
      REAL(pr) :: FA, FB, FC, tA, tB, tC
 
      REAL(pr), PARAMETER :: GOLD = (1.0_pr + SQRT(5.0_pr))/2.0_pr
      REAL(pr), PARAMETER :: CGOLD = 1.0_pr/GOLD  
      REAL(pr), PARAMETER :: GLIMIT = 10.0_pr
      REAL(pr), PARAMETER :: tMAX = 10.0_pr
      INTEGER, PARAMETER :: ITMAX = 100
      INTEGER :: FuncEval, iter 
      LOGICAL :: saveLineMin

      saveLineMin = .TRUE.

      ALLOCATE( phi_bar(1:n(1),1:n(2),1:local_nlast,1:3) )

      FuncEval = 0
      iter = 0
 
      tA = tA0
      tB = MAX(tB0, MACH_EPSILON)

      phi_bar = phi + tA*grad
      FA = eval_J(phi_bar, "LineMin", "FixK0E0")
      FuncEval = FuncEval+1

      phi_bar = phi + tB*grad
      FB = eval_J(phi_bar, "LineMin", "FixK0E0")
      FuncEval = FuncEval+1

      IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, "FixK0E0", "replace")
 
      DO WHILE (FB > FA .AND. tB > MACH_EPSILON)
         tB = CGOLD*tB
         phi_bar = phi + tB*grad
         FB = eval_J(phi_bar, "LineMin", "FixK0E0")
         FuncEval = FuncEval+1
         IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, "FixK0E0", "append")
      END DO

      IF (tB .LE. 2.0_pr*MACH_EPSILON) THEN
         myflag = 11
         RETURN
      END IF

      tC = GOLD*tB
      phi_bar = phi + tC*grad
      FC = eval_J(phi_bar, "LineMin", "FixK0E0")
      FuncEval = FuncEval+1

      IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, "FixK0E0", "append")

      DO WHILE (FB>=FC .AND. iter<ITMAX)
         iter = iter+1
         !IF (FuncEval < ITMAX .AND. tC < tMAX ) THEN
         tC = GOLD*tC
         phi_bar = phi + tC*grad
         FC = eval_J(phi_bar, "LineMin", "FixK0E0")
         FuncEval = FuncEval+1
         !END IF
                  
         R = (tB-tA)*(FB-FC)
         Q = (tB-tC)*(FB-FA)
         tP = tB - 0.5_pr*((tB-tC)*Q - (tB-tA)*R)/( SIGN( MAX(ABS(Q-R),MACH_EPSILON), Q-R) )
    
         Pmax = tB + GLIMIT*(tC-tB)
    
         IF ( (tB-tP)*(tP-tC)>0 ) THEN
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", "FixK0E0")
        
            IF (FP<FC) THEN
               tA = tB
               FA = FB
               tB = tP
               FB = FP
               EXIT 
            ELSEIF (FP>FB) THEN
               tC = tP
               FC = FP
               EXIT
            END IF
        
            tP = tC + GOLD*(tC-tB)
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", "FixK0E0")
        
         ELSEIF ( (tC-tP)*(tP-Pmax)>0 ) THEN
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", "FixK0E0")
       
            IF (FP<FC) THEN
               tB = tC
               tC = tP
               FB = FC
               FC = FP
               tP = tC+GOLD*(tC-tB)
               phi_bar = phi + tP*grad
               FP = eval_J(phi_bar, "LineMin", "FixK0E0")
            END IF
        
         ELSEIF ( (tP-Pmax)*(Pmax-tC)>=0 ) THEN
            tP = Pmax
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", "FixK0E0")
        
         ELSE
            tP = tC + GOLD*(tC-tB)
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", "FixK0E0")
         END IF
    
         tA = tB
         tB = tC
         tC = tP
         FA = FB
         FB = FC
         FC = FP
        
         IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, "FixK0E0", "append")
 
      END DO

      tau_brack(1) = tA
      tau_brack(2) = tC

      IF (iter .GE. ITMAX) THEN
         myflag = 12
      ELSE
         myflag = 0
      END IF

      DEALLOCATE( phi_bar )
    
    END FUNCTION mnbrak_FixK0E0

  
    !============================================
    ! BRENT ALGORITHM FOR LINE OPTIMIZATION
    !============================================
    FUNCTION brent_FixK0E0(phi, grad, tau_brack) RESULT (X)
      USE global_variables
      USE data_ops
      USE function_ops
      IMPLICIT NONE

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: phi, grad
      REAL(pr), DIMENSION(1:2), INTENT(IN) :: tau_brack
      REAL(pr) :: X
   
      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi_bar
      REAL(pr) :: D, A, B, V, W, E, ETEMP, P, Q, R, U, XM
      REAL(pr) :: FV, FW, FU, FX
      REAL(pr) :: TOL1, TOL2
      INTEGER :: FLAG, j


      INTEGER, PARAMETER :: ITMAX = 300
      REAL(pr), PARAMETER :: TOL = 1E-6
      REAL(pr), PARAMETER :: ZEPS = 1E-8
      REAL(pr), PARAMETER :: CGOLD = .381966


      ALLOCATE( phi_bar(1:n(1),1:n(2),1:local_nlast,1:3) )
 
      D = 0.0_pr
      A = MIN(tau_brack(1),tau_brack(2))
      B = MAX(tau_brack(1),tau_brack(2))
      V = tau_brack(2)*CGOLD 
      W = V
      X = V 
      E = 0.0_pr

      phi_bar = phi + X*grad
      FX = eval_J(phi_bar, "LineMin", "FixK0E0")
 
      FV = FX 
      FW = FX

      DO j=1,ITMAX
        XM = 0.5_pr*(A+B)
        TOL1 = TOL*ABS(X)+ZEPS
        TOL2 = 2.0_pr*TOL1
    
        IF ( ABS(X-XM) <= (TOL2-0.5*(B-A)) ) EXIT
    
        FLAG = 1
        IF ( ABS(E) > TOL1 ) THEN
           R = (X-W)*(FX-FV)
           Q = (X-V)*(FX-FW)
           P = (X-V)*Q - (X-W)*R
           Q = 2.0_pr*(Q-R)
           IF ( Q > 0.0_pr ) P = -P 
      
           Q = ABS(Q)
           ETEMP = E
           E = D
        
           IF ( (ABS(P) >= ABS(0.5_pr*Q*ETEMP)) .OR. (P <= Q*(A-X)) .OR. (P >= Q*(B-X)) ) THEN
              FLAG = 1
           ELSE
              FLAG = 2
           END IF

        END IF
    
        SELECT CASE (FLAG)
          CASE (1)
            IF (X >= XM) THEN
                E = A-X
            ELSE
                E=B-X
            END IF
            D = CGOLD*E
          CASE (2)
            D = P/Q
            U = X+D
            IF ( (U-A < TOL2) .OR. (B-U < TOL2) ) D = SIGN(TOL1, XM-X)
            
        END SELECT
    
        IF ( ABS(D) >= TOL1 ) THEN
           U = X+D
        ELSE
           U = X + SIGN(TOL1,D)
        END IF
    
        phi_bar = phi + U*grad
        FU = eval_J(phi_bar, "LineMin", "FixK0E0")
 
        IF ( FU <= FX ) THEN
           IF ( U >= X ) THEN
              A = X
           ELSE
              B = X 
           END IF
           V = W
           FV = FW
           W = X
           FW = FX
           X = U
           FX = FU
        ELSE
           IF ( U < X ) THEN
              A = U
           ELSE
              B = U
           END IF
        
           IF ( (FU <= FW) .OR. (W == X) ) THEN
              V = W
              FV = FW
              W = U
              FW = FU
           ELSEIF ( (FU <= FV) .OR. (V==X) .OR. (V==W)) THEN 
              V = U
              FV = FU
           END IF
        END IF
      END DO

      DEALLOCATE( phi_bar )

    END FUNCTION brent_FixK0E0

    !============================================
    ! Perform kappa_test
    !============================================
    SUBROUTINE kappa_test(phi, gradJ, J0, mysystem)
      USE global_variables
      USE function_ops
      USE data_ops
      IMPLICIT NONE
      INCLUDE "mpif.h"

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: phi
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: gradJ
      REAL(pr), INTENT(IN) :: J0
      CHARACTER(len=*), INTENT(IN) :: mysystem

      REAL(pr) :: myepsilon, kappa, deltaJ
      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi_pert, phi_bar
      REAL(pr) :: myexp, J1
      REAL(pr), DIMENSION(1:3) :: dx, local_inner_prod, global_inner_prod
      INTEGER :: ii      
      CHARACTER(2) :: K0txt, E0txt, IGtxt

      WRITE(K0txt,'(i2.2)') K0_index
      WRITE(E0txt,'(i2.2)') E0_index
      WRITE(IGtxt,'(i2.2)') iguess
 
      ALLOCATE( phi_pert(1:n(1),1:n(2),1:local_nlast,1:3) )
      ALLOCATE( phi_bar(1:n(1),1:n(2),1:local_nlast,1:3) )

      CALL kappa_test_pert(phi_pert, "sine", 1.0_pr, 1.0_pr, -2.0_pr)
         
      local_inner_prod = field_inner_product(gradJ, phi_pert, "L2")
      CALL MPI_ALLREDUCE(local_inner_prod, global_inner_prod, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)  
 
      DO ii=1,21
         myexp = -10.0_pr + 0.5_pr*REAL(ii-1,pr)
         myepsilon = 10.0_pr**myexp         
         
         phi_bar = phi + myepsilon*phi_pert
         IF (rank==0) THEN
            OPEN(10, FILE="../LOGFILES/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_info.log", STATUS='OLD', POSITION='APPEND')
            WRITE(10,*) " =========================== "
            WRITE(10,*) "   Kappa index = ", ii
            WRITE(10,*) "   Epsilon = ", myepsilon
            CLOSE(10)
         END IF
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)         

         J1 = eval_J(phi_bar, mysystem, "Enstrophy")
         
         kappa = (J1-J0)/(myepsilon*SUM(global_inner_prod))
         deltaJ = J1-J0  
         IF (rank==0) THEN
            CALL save_kappa_test(myepsilon, SUM(global_inner_prod), deltaJ, kappa, ii, mysystem)
         END IF 
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END DO

      DEALLOCATE(phi_pert)
      DEALLOCATE(phi_bar)

    END SUBROUTINE kappa_test
 
END MODULE

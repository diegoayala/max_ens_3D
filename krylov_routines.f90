MODULE krylov_routines
       USE global_variables
       IMPLICIT NONE

       INTEGER, PARAMETER :: MAX_TIME_STEPS = 500
       INTEGER, PARAMETER :: kry_l = 25
       INTEGER, PARAMETER :: kry_m = 25, n_grid = 2, j_it = 2
       REAL(pr), SAVE :: dtwrite 
       REAL(pr), PARAMETER :: tol = 1.0e-4
    
     CONTAINS
       !============================================
       !  KRYLOV_FWD
       !============================================
       SUBROUTINE krylov_fwd(u_b, mysystem)
         USE global_variables
         USE RHS_routines
         USE data_ops
         USE function_ops
         IMPLICIT NONE
         INCLUDE "mpif.h"

         REAL(pr), DIMENSION(1:n_dim), INTENT (INOUT)  :: u_b
         CHARACTER(len=*), INTENT(IN) :: mysystem

         REAL(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: ub
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: wb

         INTEGER :: i, j, j_p, m
         REAL(pr) :: dt, dt_max, time, ub_local_max, ub_global_max
         REAL(pr) :: a, localerr, err_max, err_tot, mod_cj, mod_rhs_U, mod_Ub

         REAL(pr), DIMENSION (n_dim):: drhs, rhs, rhs_u, u_aux1, u_aux2, zsum, zsum_err
         REAL(pr), DIMENSION (n_grid) :: tau_p
         REAL(pr), DIMENSION (n_dim, n_grid) :: R_taup, u_err, u_tau
         REAL(pr), DIMENSION (n_dim, kry_l) :: v
         REAL(pr), DIMENSION (n_dim, j_it) :: c_j
         REAL(pr), DIMENSION (kry_l, kry_l) :: hess 
         REAL(pr), DIMENSION (n_dim, kry_m, j_it) :: v_m
         REAL(pr), DIMENSION (kry_m, kry_m, j_it) :: h_m 
         LOGICAL :: ireturn
         INTEGER :: iwrite

         REAL(pr) :: cfl, ens, palins, w_max
         REAL(pr), DIMENSION(1:2) :: ener, u_max

         !Determine CFL condition.
         ub_local_max = MAXVAL(ABS(u_b))
         CALL MPI_ALLREDUCE(ub_local_max, ub_global_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
         dt_max = 1.0_pr/(ub_global_max*n(1))
         dt = 0.1_pr*T
         DO WHILE (dt>dt_max)
            dt = 0.95_pr*dt
         END DO
         dtwrite = T/REAL(MAX_TIME_STEPS, pr)
         
         iwrite = 0
         time = T0
         CALL init_vars(u_b, ub, wb)
         
         CALL save_NS_velocity(u_b, iwrite, ".bin", mysystem)
         
         IF (save_diag_NS) THEN
            CALL save_NS_vorticity(wb, iwrite, ".nc", mysystem) 
            CALL diagnostics(time, ub, wb, 0, iwrite, mysystem)
            !CALL save_image(iwrite, wb, mysystem,1)
         END IF

         IF (rank==0) THEN
           CALL save_timevec(time, 0, mysystem)
         END IF
         CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
    
         iwrite = iwrite + 1
 
         DO WHILE (time <= T)  

            !--Generate orthonormal basis for Krylov subspace v and approximate operator h
            ub_local_max = MAXVAL(ABS(u_b))
            CALL MPI_ALLREDUCE(ub_local_max, ub_global_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
            a = DOT_PRODUCT (u_b, u_b)
            CALL MPI_ALLREDUCE(a, mod_Ub, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            mod_Ub = SQRT (mod_Ub)

            CALL rhs_kry(u_b, rhs_u)
               
            a = DOT_PRODUCT (rhs_u, rhs_u)
            CALL MPI_ALLREDUCE(a, mod_rhs_u, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            mod_rhs_u = SQRT(mod_rhs_u)
            CALL arnoldi (kry_l, rhs_u, mod_rhs_u, v, hess, u_b, ub, wb)

            ireturn = .TRUE.
            m = 1
            DO WHILE (ireturn)
               !--Chebychev grid     
               DO j_p = 1, n_grid
                  tau_p(j_p) = 0.5_pr*dt*(1.0_pr - COS (j_p * pi / n_grid))

                  !--Linear solution at each grid point
                  u_tau(:,j_p) = MATMUL (v, MATMUL (exp_mat (hess, tau_p(j_p), 1), u_vec (kry_l, mod_rhs_u)))
                  u_err(:,j_p) = MATMUL (v, MATMUL (exp_mat (hess, tau_p(j_p), 2), u_vec (kry_l, mod_rhs_u)))

                  !--Remainder term	
                  R_taup(:,j_p) = rmd_t (u_tau(:,j_p), u_b, rhs_u, ub, wb)  
               END DO

               !--Obtain the coefficients C_j
               CALL least_square (n_dim, n_grid, j_it, tau_p, R_taup, c_j)

               !--Form krylov subspace for each polynomial
               DO j = 1, j_it
                 a = DOT_PRODUCT (c_j(:,j), c_j(:,j))
                 CALL MPI_ALLREDUCE(a,mod_cj,1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                 mod_cj = SQRT(mod_cj)
                 CALL arnoldi (kry_m, c_j(:,j), mod_cj, v_m(:,:,j), h_m(:,:,j), u_b, ub, wb)
               END DO

               !--Calcul de la somme non lineaire pour chaque point de la grille
               DO j_p = 1, n_grid
                 zsum = u_tau(:, j_p)
                 zsum_err = u_err(:, j_p)
                 DO j = 1, j_it
                   a = DOT_PRODUCT (c_j(:,j), c_j(:,j))
                   CALL MPI_ALLREDUCE(a, mod_cj, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                   mod_cj = SQRT (mod_cj)
                   u_aux1 = MATMUL (v_m(:,:,j), MATMUL (integ (h_m(:,:,j), j+1, tau_p(j_p), 1), u_vec(kry_m, mod_cj)))
                   u_aux2 = MATMUL (v_m(:,:,j), MATMUL (integ (h_m(:,:,j), j+1, tau_p(j_p), 2), u_vec(kry_m, mod_cj))) &
                            + c_j(:,j) * tau_p(j_p)**(j+1)
                   zsum     = zsum     + u_aux1
                   zsum_err = zsum_err + u_aux2
                 END DO
                 !--Final form for u at each grid-point
                 u_tau(:,j_p) = zsum
               END DO

               !--Estimate linf error
               CALL rhs_kry (u_tau(:,n_grid)+u_b, rhs)
               localerr = MAXVAL(ABS(zsum_err - rhs))
               CALL MPI_ALLREDUCE (localerr, err_tot, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
               err_max = tau_p(n_grid) * err_tot/ub_global_max

               IF (err_max<tol) THEN

                 u_b = u_b + u_tau(:,n_grid)
                 time = time + tau_p(n_grid)
                 IF (err_max < tol*0.15_pr) THEN
                   dt = tau_p(n_grid) * 1.5_pr
                 ELSE
                   dt = tau_p(n_grid)
                 END IF
                 ireturn = .FALSE. 

                 !--Output analysis
                 CALL vecmat (u_b, ub)
                 CALL cal_vort (ub, wb)

                 IF (iwrite*dtwrite < time) THEN
                    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

                    CALL save_NS_velocity(u_b, iwrite, ".bin", mysystem)
                    IF (save_diag_NS) THEN
                       CALL save_NS_vorticity(wb, iwrite, ".nc", mysystem)
                       CALL diagnostics(time, ub, wb, 1, iwrite, mysystem)
                       !CALL save_image(iwrite, wb, mysystem,1)
                    END IF
                    iwrite = iwrite+1
                    IF (rank==0) THEN
                       CALL save_timevec(time, 1, mysystem)
                    END IF
                    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
                 END IF

               ELSE
                 dt = dt * 0.75_pr
                 m  = m+1
                 ireturn = .TRUE.
               END IF
            END DO
         END DO
         time_steps = iwrite

       END SUBROUTINE krylov_fwd

       !===============================================
       !  KRYLOV_BWD (Takes into account that adjoint
       !  is a linear system.
       !===============================================
       SUBROUTINE krylov_bwd(u_b)
         USE global_variables
         USE RHS_routines
         USE data_ops
         USE function_ops
         IMPLICIT NONE
         INCLUDE "mpif.h"

         REAL(pr), DIMENSION(1:n_dim), INTENT (INOUT)  :: u_b

         REAL(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: ub
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: wb, stream

         INTEGER :: i, j, j_p, m
         REAL(pr) :: dt, dt_max, time, ub_local_max, ub_global_max
         REAL(pr) :: a, localerr, err_max, err_tot, mod_cj, mod_rhs_U, mod_Ub

         REAL(pr), DIMENSION (n_dim):: drhs, rhs, rhs_u, u_aux1, u_aux2, zsum, zsum_err
         REAL(pr), DIMENSION (n_grid) :: tau_p
         REAL(pr), DIMENSION (n_dim, n_grid) :: R_taup, u_err, u_tau
         REAL(pr), DIMENSION (n_dim, kry_l) :: v
         REAL(pr), DIMENSION (n_dim, j_it) :: c_j
         REAL(pr), DIMENSION (kry_l, kry_l) :: hess 
         REAL(pr), DIMENSION (n_dim, kry_m, j_it) :: v_m
         REAL(pr), DIMENSION (kry_m, kry_m, j_it) :: h_m 
         LOGICAL :: ireturn
         INTEGER :: iwrite

         REAL(pr) :: cfl, ens, palins, w_max
         REAL(pr), DIMENSION(1:2) :: ener, u_max

         !Determine CFL condition.
         ub_local_max = MAXVAL(ABS(u_b))
         CALL MPI_ALLREDUCE(ub_local_max, ub_global_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
         dt_max = 1.0_pr/(ub_global_max*n(1))
         dt = 0.1_pr*T
         DO WHILE (dt>dt_max)
            dt = 0.95_pr*dt
         END DO
         dtwrite = T/REAL(MAX_TIME_STEPS, pr)
         
         iwrite = 0
         time = T
         CALL init_vars(u_b, ub, wb)
         CALL cal_stream(wb, stream)
         
         CALL save_adj_velocity(u_b, iwrite, ".bin")
         IF (save_diag_Adj) THEN
            CALL save_adj_vorticity(stream, iwrite, ".nc") 
            CALL diagnostics(time, ub, wb, 0, iwrite, "adj")
         END IF  
         IF (rank==0) THEN
            CALL save_adj_timevec(time, 0)
         END IF
         CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
    
         iwrite = iwrite + 1
 
         DO WHILE (time > 0.0_pr)  

            !--Generate orthonormal basis for Krylov subspace v and approximate operator h
            ub_local_max = MAXVAL(ABS(u_b))
            CALL MPI_ALLREDUCE(ub_local_max, ub_global_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
            a = DOT_PRODUCT (u_b, u_b)
            CALL MPI_ALLREDUCE(a, mod_Ub, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            mod_Ub = SQRT (mod_Ub)

            CALL rhs_adj_kry(u_b, rhs_u, time)
               
            a = DOT_PRODUCT (rhs_u, rhs_u)
            CALL MPI_ALLREDUCE(a, mod_rhs_u, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            mod_rhs_u = SQRT(mod_rhs_u)
            CALL arnoldi_adj2(kry_l, rhs_u, mod_rhs_u, v, hess, u_b, ub, wb, time)

            ireturn = .TRUE.
            m = 1
            DO WHILE (ireturn)
               !--Chebychev grid     
               DO j_p = 1, n_grid
                  tau_p(j_p) = -0.5_pr*dt*(1.0_pr - COS (j_p * pi / n_grid))

                  !--Linear solution at each grid point
                  u_tau(:,j_p) = MATMUL (v, MATMUL (exp_mat (hess, tau_p(j_p), 1), u_vec (kry_l, mod_rhs_u)))
                  u_err(:,j_p) = MATMUL (v, MATMUL (exp_mat (hess, tau_p(j_p), 2), u_vec (kry_l, mod_rhs_u)))

                  !--Remainder term	
                  R_taup(:,j_p) = rmd_t_adj(u_tau(:,j_p), u_b, rhs_u, ub, wb, time)  
               END DO

               !--Obtain the coefficients C_j
               CALL least_square (n_dim, n_grid, j_it, tau_p, R_taup, c_j)

               !--Form krylov subspace for each polynomial
               DO j = 1, j_it
                 a = DOT_PRODUCT (c_j(:,j), c_j(:,j))
                 CALL MPI_ALLREDUCE(a,mod_cj,1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                 mod_cj = SQRT(mod_cj)
                 CALL arnoldi_adj2(kry_m, c_j(:,j), mod_cj, v_m(:,:,j), h_m(:,:,j), u_b, ub, wb, time)
               END DO

               !--Calcul de la somme non lineaire pour chaque point de la grille
               DO j_p = 1, n_grid
                 zsum = u_tau(:, j_p)
                 zsum_err = u_err(:, j_p)
                 DO j = 1, j_it
                   a = DOT_PRODUCT (c_j(:,j), c_j(:,j))
                   CALL MPI_ALLREDUCE(a, mod_cj, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                   mod_cj = SQRT (mod_cj)
                   u_aux1 = MATMUL (v_m(:,:,j), MATMUL (integ (h_m(:,:,j), j+1, tau_p(j_p), 1), u_vec(kry_m, mod_cj)))
                   u_aux2 = MATMUL (v_m(:,:,j), MATMUL (integ (h_m(:,:,j), j+1, tau_p(j_p), 2), u_vec(kry_m, mod_cj))) &
                            + c_j(:,j) * tau_p(j_p)**(j+1)
                   zsum     = zsum     + u_aux1
                   zsum_err = zsum_err + u_aux2
                 END DO
                 !--Final form for u at each grid-point
                 u_tau(:,j_p) = zsum
               END DO

               !--Estimate linf error
               CALL rhs_adj_kry(u_tau(:,n_grid)+u_b, rhs, time)
               localerr = MAXVAL(ABS(zsum_err - rhs))
               CALL MPI_ALLREDUCE (localerr, err_tot, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
               err_max = tau_p(n_grid) * err_tot/ub_global_max

               IF (err_max<tol) THEN

                 u_b = u_b + u_tau(:,n_grid)
                 time = time + tau_p(n_grid)
                 IF (err_max < tol*0.15_pr) THEN
                   dt = tau_p(n_grid) * 1.5_pr
                 ELSE
                   dt = tau_p(n_grid)
                 END IF
                 ireturn = .FALSE. 

                 !--Output analysis
                 CALL vecmat (u_b, ub)
                 CALL cal_vort (ub, wb)
                 CALL cal_stream(wb, stream)
         
                 IF (iwrite*dtwrite < time) THEN
                    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

                    CALL save_adj_velocity(u_b, iwrite, ".bin")
                    IF (save_diag_Adj) THEN
                       CALL save_adj_vorticity(stream, iwrite, ".nc")
                       CALL diagnostics(time, ub, wb, 1, iwrite, "adj")
                    END IF
                    iwrite = iwrite+1
                    IF (rank==0) THEN
                       CALL save_adj_timevec(time, 1)
                    END IF
                    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
                 END IF

               ELSE
                 dt = dt * 0.75_pr
                 m  = m+1
                 ireturn = .TRUE.
               END IF
            END DO
         END DO
         adj_time_steps = iwrite


       END SUBROUTINE krylov_bwd



       !===============================================
       ! Another way of solving the adjoint equation
       !===============================================
       SUBROUTINE krylov_bwd2(u_b, init_time, myindex)
         USE global_variables
         USE RHS_routines
         USE data_ops
         USE function_ops
         IMPLICIT NONE
         INCLUDE "mpif.h"

         REAL(pr), DIMENSION(1:n_dim), INTENT (INOUT)  :: u_b
         REAL(pr), INTENT(IN) :: init_time
         INTEGER, INTENT(IN) :: myindex

         REAL(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: ub
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: wb, stream

         INTEGER :: i, j, j_p, m
         REAL(pr) :: dt, dt_max, time, ub_local_max, ub_global_max
         REAL(pr) :: a, err, err_max, err_tot, mod_cj, mod_rhs_U, mod_Ub

         REAL(pr), DIMENSION (n_dim):: drhs, rhs, rhs_u, u_aux1, u_aux2, zsum, zsum_err
         REAL(pr), DIMENSION (n_grid) :: tau_p
         REAL(pr), DIMENSION (n_dim, n_grid) :: R_taup, u_err, u_tau
         REAL(pr), DIMENSION (n_dim, kry_l) :: v
         REAL(pr), DIMENSION (n_dim, j_it) :: c_j
         REAL(pr), DIMENSION (kry_l, kry_l) :: hess 
         REAL(pr), DIMENSION (n_dim, kry_m, j_it) :: v_m
         REAL(pr), DIMENSION (kry_m, kry_m, j_it) :: h_m 
         LOGICAL :: ireturn
         INTEGER :: iwrite

         REAL(pr) :: cfl, ens, palins, w_max
         REAL(pr), DIMENSION(1:2) :: ener, u_max

         !Determine CFL condition.
         IF (rank==0) THEN
            CALL get_umax(ub_global_max)
         END IF
         CALL MPI_BCAST(ub_global_max, 1, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)
         dt_max = 1.0_pr/(ub_global_max*n(1))

         dt = 0.1_pr*T
         DO WHILE (dt>dt_max)
            dt = 0.95_pr*dt
         END DO
         dtwrite = T/REAL(MAX_TIME_STEPS, pr)
         
         IF (myindex > 0) THEN 
            iwrite = myindex
            time = init_time
         ELSE
            iwrite = 0
            time = 0.0_pr
            CALL init_vars(u_b, ub, wb)
            CALL cal_stream(wb, stream)
            CALL save_adj_velocity(u_b, iwrite, ".bin")
            IF (save_diag_Adj) THEN
               CALL save_adj_vorticity(stream, iwrite, ".nc") 
               CALL diagnostics(T-time, ub, wb, 0, iwrite, "adj")
               CALL save_image(iwrite, stream, "adj", 1)
            END IF  
            IF (rank==0) THEN
              CALL save_adj_timevec(T-time, 0)
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)

         END IF
         
         iwrite = iwrite + 1
 
         DO WHILE (time <= T)  

            !--Generate orthonormal basis for Krylov subspace v and approximate operator h
            ub_local_max = MAXVAL(ABS(u_b))
            CALL MPI_ALLREDUCE(ub_local_max, ub_global_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
            a = DOT_PRODUCT (u_b, u_b)
            CALL MPI_ALLREDUCE(a, mod_Ub, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            mod_Ub = SQRT (mod_Ub)

            CALL rhs_adj_kry(u_b, rhs_u, time)

            a = DOT_PRODUCT (rhs_u, rhs_u)
            CALL MPI_ALLREDUCE(a, mod_rhs_u, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            mod_rhs_u = SQRT(mod_rhs_u)
            CALL arnoldi_adj2(kry_l, rhs_u, mod_rhs_u, v, hess, u_b, ub, wb, time)

            ireturn = .TRUE.
            m = 1
            DO WHILE (ireturn)
               !--Chebychev grid     
               DO j_p = 1, n_grid
                  tau_p(j_p) = 0.5_pr*dt*(1.0_pr - COS (j_p * pi / n_grid))

                  !--Linear solution at each grid point
                  u_tau(:,j_p) = MATMUL (v, MATMUL (exp_mat (hess, tau_p(j_p), 1), u_vec (kry_l, mod_rhs_u)))
                  u_err(:,j_p) = MATMUL (v, MATMUL (exp_mat (hess, tau_p(j_p), 2), u_vec (kry_l, mod_rhs_u)))

                  !--Remainder term	
                  R_taup(:,j_p) = rmd_t_adj(u_tau(:,j_p), U_b, rhs_u, ub, wb, time)  
               END DO

               !--Obtain the coefficients C_j
               CALL least_square (n_dim, n_grid, j_it, tau_p, R_taup, c_j)

               !--Form krylov subspace for each polynomial
               DO j = 1, j_it
                 a = DOT_PRODUCT (c_j(:,j), c_j(:,j))
                 CALL MPI_ALLREDUCE(a,mod_cj,1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                 mod_cj = SQRT(mod_cj)
                 CALL arnoldi_adj2(kry_m, c_j(:,j), mod_cj, v_m(:,:,j), h_m(:,:,j), u_b, ub, wb, time)
               END DO

               !--Calcul de la somme non lineaire pour chaque point de la grille
               DO j_p = 1, n_grid
                 zsum = u_tau(:, j_p)
                 zsum_err = u_err(:, j_p)
                 DO j = 1, j_it
                   a = DOT_PRODUCT (c_j(:,j), c_j(:,j))
                   CALL MPI_ALLREDUCE(a, mod_cj, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                   mod_cj = SQRT (mod_cj)
                   u_aux1 = MATMUL (v_m(:,:,j), MATMUL (integ (h_m(:,:,j), j+1, tau_p(j_p), 1), u_vec(kry_m, mod_cj)))
                   u_aux2 = MATMUL (v_m(:,:,j), MATMUL (integ (h_m(:,:,j), j+1, tau_p(j_p), 2), u_vec(kry_m, mod_cj))) &
                            + c_j(:,j) * tau_p(j_p)**(j+1)
                   zsum     = zsum     + u_aux1
                   zsum_err = zsum_err + u_aux2
                 END DO
                 !--Final form for u at each grid-point
                 u_tau(:,j_p) = zsum
               END DO

               !--Estimate linf error
               CALL rhs_adj_kry (u_tau(:,n_grid)+u_b, rhs, time)
               err = MAXVAL(ABS(zsum_err - rhs))
               CALL MPI_ALLREDUCE (err, err_tot, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
               err_max = tau_p(n_grid) * err_tot/ub_global_max

               IF (err_max<tol) THEN

                 u_b = u_b + u_tau(:,n_grid)
                 time = time + tau_p(n_grid)
                 IF (err_max < tol*0.15_pr) THEN
                   dt = tau_p(n_grid) * 1.5_pr
                 ELSE
                   dt = tau_p(n_grid)
                 END IF
                 ireturn = .FALSE. 

                 !--Output analysis
                 CALL vecmat (u_b, ub)
                 CALL cal_vort (ub, wb)
                 CALL cal_stream(wb,stream)

                 IF (iwrite*dtwrite < time) THEN
                    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
                    CALL save_adj_velocity(u_b, iwrite, ".bin")
                    IF (save_diag_Adj) THEN
                       CALL save_adj_vorticity(stream, iwrite, ".nc")
                       CALL diagnostics(T-time, ub, wb, 1, iwrite, "adj")
                       CALL save_image(iwrite, stream, "adj", 1)
                    END IF
                    iwrite = iwrite+1
                    IF (rank==0) THEN
                       CALL save_adj_timevec(T-time, 1)
                    END IF
                    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
                 END IF

               ELSE
                 dt = dt * 0.75_pr
                 m  = m+1
                 ireturn = .TRUE.
               END IF
            END DO
         END DO
         adj_time_steps = iwrite

       END SUBROUTINE krylov_bwd2


       !===============================================
       ! SOLUTION OF THE PERTURBATION EQUATION
       !===============================================
       SUBROUTINE krylov_linear_NS(u_b)
         USE global_variables
         USE RHS_routines
         USE data_ops
         USE function_ops
         IMPLICIT NONE
         INCLUDE "mpif.h"

         REAL(pr), DIMENSION(:), INTENT (INOUT)  :: u_b

         REAL(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: ub
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: wb, stream

         INTEGER :: i, j, j_p, m
         REAL(pr) :: dt, dt_max, time, ub_local_max, ub_global_max
         REAL(pr) :: a, err, err_max, err_tot, mod_cj, mod_rhs_U, mod_Ub

         REAL(pr), DIMENSION (n_dim):: drhs, rhs, rhs_u, u_aux1, u_aux2, zsum, zsum_err
         REAL(pr), DIMENSION (n_grid) :: tau_p
         REAL(pr), DIMENSION (n_dim, n_grid) :: R_taup, u_err, u_tau
         REAL(pr), DIMENSION (n_dim, kry_l) :: v
         REAL(pr), DIMENSION (n_dim, j_it) :: c_j
         REAL(pr), DIMENSION (kry_l, kry_l) :: hess 
         REAL(pr), DIMENSION (n_dim, kry_m, j_it) :: v_m
         REAL(pr), DIMENSION (kry_m, kry_m, j_it) :: h_m 
         LOGICAL :: ireturn
         INTEGER :: iwrite

         REAL(pr) :: cfl, ens, palins, w_max
         REAL(pr), DIMENSION(1:2) :: ener, u_max

         !Determine CFL condition.
         IF (rank==0) THEN
            CALL get_umax(ub_global_max)
         END IF
         CALL MPI_BCAST(ub_global_max, 1, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)
         dt_max = 1.0_pr/(ub_global_max*n(1))

         dt = 0.1_pr*T
         DO WHILE (dt>dt_max)
            dt = 0.95_pr*dt
         END DO
         dtwrite = T/REAL(MAX_TIME_STEPS, pr)
         
         iwrite = 0
         time = 0.0_pr
         CALL init_vars(u_b, ub, wb)
         CALL save_NS_velocity(u_b, iwrite, ".bin", "d1_nse")
         IF (save_diag_D1NSE) THEN
            CALL save_NS_vorticity(wb, iwrite, ".nc", "d1_nse") 
            CALL diagnostics(time, ub, wb, 0, iwrite, "d1_nse")
         END IF  
         IF (rank==0) THEN
            CALL save_timevec(time, 0, "d1_nse")
         END IF
         CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)

         
         iwrite = iwrite + 1
 
         DO WHILE (time <= T)  

            !--Generate orthonormal basis for Krylov subspace v and approximate operator h
            ub_local_max = MAXVAL(ABS(u_b))
            CALL MPI_ALLREDUCE(ub_local_max, ub_global_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
            a = DOT_PRODUCT (u_b, u_b)
            CALL MPI_ALLREDUCE(a, mod_Ub, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            mod_Ub = SQRT (mod_Ub)

            CALL rhs_d1nse_kry(u_b, rhs_u, time)

            a = DOT_PRODUCT (rhs_u, rhs_u)
            CALL MPI_ALLREDUCE(a, mod_rhs_u, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            mod_rhs_u = SQRT(mod_rhs_u)
            CALL arnoldi_d1nse(kry_l, rhs_u, mod_rhs_u, v, hess, u_b, ub, wb, time)

            ireturn = .TRUE.
            m = 1
            DO WHILE (ireturn)
               !--Chebychev grid     
               DO j_p = 1, n_grid
                  tau_p(j_p) = 0.5_pr*dt*(1.0_pr - COS (j_p * pi / n_grid))

                  !--Linear solution at each grid point
                  u_tau(:,j_p) = MATMUL (v, MATMUL (exp_mat (hess, tau_p(j_p), 1), u_vec (kry_l, mod_rhs_u)))
                  u_err(:,j_p) = MATMUL (v, MATMUL (exp_mat (hess, tau_p(j_p), 2), u_vec (kry_l, mod_rhs_u)))

                  !--Remainder term	
                  R_taup(:,j_p) = rmd_t_adj (u_tau(:,j_p), U_b, rhs_u, ub, wb, time)  
               END DO

               !--Obtain the coefficients C_j
               CALL least_square (n_dim, n_grid, j_it, tau_p, R_taup, c_j)

               !--Form krylov subspace for each polynomial
               DO j = 1, j_it
                 a = DOT_PRODUCT (c_j(:,j), c_j(:,j))
                 CALL MPI_ALLREDUCE(a,mod_cj,1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                 mod_cj = SQRT(mod_cj)
                 CALL arnoldi_d1nse(kry_m, c_j(:,j), mod_cj, v_m(:,:,j), h_m(:,:,j), u_b, ub, wb, time)
               END DO

               !--Calcul de la somme non lineaire pour chaque point de la grille
               DO j_p = 1, n_grid
                 zsum = u_tau(:, j_p)
                 zsum_err = u_err(:, j_p)
                 DO j = 1, j_it
                   a = DOT_PRODUCT (c_j(:,j), c_j(:,j))
                   CALL MPI_ALLREDUCE(a, mod_cj, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                   mod_cj = SQRT (mod_cj)
                   u_aux1 = MATMUL (v_m(:,:,j), MATMUL (integ (h_m(:,:,j), j+1, tau_p(j_p), 1), u_vec(kry_m, mod_cj)))
                   u_aux2 = MATMUL (v_m(:,:,j), MATMUL (integ (h_m(:,:,j), j+1, tau_p(j_p), 2), u_vec(kry_m, mod_cj))) &
                            + c_j(:,j) * tau_p(j_p)**(j+1)
                   zsum     = zsum     + u_aux1
                   zsum_err = zsum_err + u_aux2
                 END DO
                 !--Final form for u at each grid-point
                 u_tau(:,j_p) = zsum
               END DO

               !--Estimate linf error
               CALL rhs_d1nse_kry (u_tau(:,n_grid)+u_b, rhs, time)
               err = MAXVAL(ABS(zsum_err - rhs))
               CALL MPI_ALLREDUCE (err, err_tot, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
               err_max = tau_p(n_grid) * err_tot/ub_global_max

               IF (err_max<tol) THEN

                 u_b = u_b + u_tau(:,n_grid)
                 time = time + tau_p(n_grid)
                 IF (err_max < tol*0.15_pr) THEN
                   dt = tau_p(n_grid) * 1.5_pr
                 ELSE
                   dt = tau_p(n_grid)
                 END IF
                 ireturn = .FALSE. 

                 !--Output analysis
                 CALL vecmat (u_b, ub)
                 CALL cal_vort (ub, wb)

                 IF (iwrite*dtwrite < time) THEN
                    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
                    CALL save_NS_velocity(u_b, iwrite, ".bin", "d1_nse")
                    IF (save_diag_D1NSE) THEN
                       CALL save_NS_vorticity(stream, iwrite, ".nc", "d1_nse")
                       CALL diagnostics(time, ub, wb, 1, iwrite, "d1_nse")
                    END IF
                    iwrite = iwrite+1
                    IF (rank==0) THEN
                       CALL save_timevec(time, 1, "d1_nse")
                    END IF
                    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
                 END IF

               ELSE
                 dt = dt * 0.75_pr
                 m  = m+1
                 ireturn = .TRUE.
               END IF
            END DO
         END DO
         d1NSE_time_steps = iwrite

       END SUBROUTINE krylov_linear_NS


       !============================================
       ! Initialize Krylov variables U0 and W0
       !============================================
       SUBROUTINE init_vars(u_kry, ub, wb) 
         USE global_variables
         USE function_ops
         USE data_ops
         USE FFT2
         IMPLICIT NONE
         INCLUDE 'mpif.h'
         REAL(pr), DIMENSION(n_dim) :: u_kry
         REAL(pr), DIMENSION(1:n(1),1:local_nlast,1:2), INTENT(OUT) :: ub
         REAL(pr), DIMENSION(1:n(1),1:local_nlast), INTENT(OUT) :: wb

         INTEGER :: i
         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: fu
         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: fw, aux, faux
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: stream

         !--Store the velocity in a matrix and save
         CALL vecmat (u_kry, ub)

         !--Save initial vorticity
         DO i = 1, 2
            aux = CMPLX(ub(:,:,i), 0.0_pr)
            CALL ffourier(aux, faux)
            fu(:,:,i) = faux
         END DO

         CALL cal_w(fu, fw)

         CALL bfourier (fw, aux)
         wb = REAL(aux)

       END SUBROUTINE init_vars


       !==============================================
       !   EXP_MAT
       !==============================================
       FUNCTION exp_mat (hess, dt, flag)

         USE global_variables
         IMPLICIT NONE

         INTEGER, INTENT (IN) :: flag
         REAL (pr), INTENT (IN) :: dt
         REAL (pr), DIMENSION (kry_l, kry_l), INTENT (IN) :: hess
  
         INTEGER :: i, j
         REAL (pr), DIMENSION (kry_l, kry_l) :: exp_mat
         COMPLEX (pr), DIMENSION (kry_l) :: eval
         COMPLEX (pr), DIMENSION (kry_l, kry_l) :: evec, evec_inv

         !--Variables for linear algebra routines
         INTEGER :: info
         INTEGER, DIMENSION (kry_l) :: ipvt
         REAL (pr), DIMENSION (34*kry_l) :: work
         REAL (pr), DIMENSION (kry_l) :: wr, wi
         REAL (pr), DIMENSION (kry_l, kry_l) :: h_tmp, vl, vr
         COMPLEX (pr), DIMENSION (2) :: det
         COMPLEX (pr), DIMENSION (kry_l) :: work1

         h_tmp = hess
         CALL dgeev('N', 'V', kry_l, h_tmp, kry_l, wr, wi, vl, 1, vr, kry_l, work, 34*kry_l, info)
         eval = CMPLX(wr, wi)

         DO i = 1, kry_l
           IF (wi(i) == 0.0_pr) THEN
             evec(:,i) = CMPLX (vr(:,i), 0.0_pr)
           ELSE IF (wi(i) > 0.0_pr) THEN
             evec(:,i)   = CMPLX (vr(:,i),  vr(:,i+1))
             evec(:,i+1) = CMPLX (vr(:,i), -vr(:,i+1))
           END IF
         END DO

         evec_inv = evec
         CALL zgetrf (kry_l, kry_l, evec_inv, kry_l, ipvt, info)
         CALL zgetri (kry_l, evec_inv, kry_l, ipvt, work1, 34*kry_l, info)

         !--Calculate exponential matrix   
         IF (flag == 1) exp_mat = REAL (MATMUL (evec, MATMUL (diag_mat(kry_l, (EXP (dt * eval) - 1.0_pr) / eval), evec_inv)), pr)  
         IF (flag == 2) exp_mat = REAL (MATMUL (evec, MATMUL (diag_mat(kry_l,  EXP (dt * eval)), evec_inv)), pr)
       END FUNCTION exp_mat


       !=============================================
       !  ARNOLDI METHOD
       !=============================================
       SUBROUTINE arnoldi(kry, u, mod_u, v, h, u_b, ub, wb)
         USE global_variables
         USE RHS_routines
         IMPLICIT NONE
         INCLUDE "mpif.h"

         INTEGER, INTENT(IN) :: kry
         REAL(pr), INTENT(IN) :: mod_u
         REAL(pr), DIMENSION(n_dim), INTENT(IN) :: u, u_b
         REAL(pr), DIMENSION(n_dim, kry), INTENT(OUT) :: v 
         REAL(pr), DIMENSION(kry, kry), INTENT (OUT) :: h
         REAL(pr), DIMENSION(1:n(1),1:local_nlast,1:2), INTENT(IN) :: ub
         REAL(pr), DIMENSION(1:n(1),1:local_nlast), INTENT(IN) :: wb

         INTEGER :: i, j
         REAL(pr), DIMENSION (n_dim) :: w
         REAL(pr) :: a

         h = 0.0_pr
         w = 0.0_pr
         v = 0.0_pr
         v(:,1) = u / mod_u 
         DO j = 1, kry
            CALL drhs_kry(v(:,j), w, ub, wb)
            DO i = 1, j
               a = DOT_PRODUCT (v(:,i), w) 
               CALL MPI_ALLREDUCE(a,h(i,j),1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)     
               w = w - v(:,i)*h(i,j) 
            END DO
            IF (j < kry) THEN
               a = DOT_PRODUCT(w, w)
               CALL MPI_ALLREDUCE(a,h(j+1,j),1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
               h(j+1,j)=SQRT(h(j+1,j))
               v(:,j+1) = w / h(j+1,j)
            END IF
         END DO
       END SUBROUTINE arnoldi

       !=============================================
       !  ARNOLDI METHOD FOR ADJOINT PROBLEM
       !=============================================
       SUBROUTINE arnoldi_adjoint(u, time, v, h)
         USE global_variables
         USE RHS_routines
         IMPLICIT NONE
         INCLUDE "mpif.h"

         REAL(pr), INTENT(IN) :: time
         REAL(pr), DIMENSION(n_dim), INTENT(IN) :: u
         REAL(pr), DIMENSION(n_dim, kry_l), INTENT(OUT) :: v 
         REAL(pr), DIMENSION(kry_l, kry_l), INTENT (OUT) :: h
         REAL(pr), DIMENSION(n_dim) :: w
         REAL(pr) :: mod_u, mod_u_local, a
         INTEGER :: i, j

         h = 0.0_pr
         v = 0.0_pr

         mod_u_local = DOT_PRODUCT(u,u)
         CALL MPI_ALLREDUCE(mod_u_local, mod_u, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
         mod_u = SQRT(mod_u)

         v(:,1) = u / mod_u 
         DO j = 1, kry_l
            CALL rhs_adj_kry(v(:,j), w, time)
            DO i = 1, j
               a = DOT_PRODUCT(v(:,i), w) 
               CALL MPI_ALLREDUCE(a, h(i,j), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)     
               w = w - v(:,i) * h(i,j) 
            END DO
            IF (j < kry_l) THEN
               a = DOT_PRODUCT(w, w)
               CALL MPI_ALLREDUCE(a,h(j+1,j),1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
               h(j+1,j)=SQRT(h(j+1,j))
               v(:,j+1) = w / h(j+1,j)
            END IF
         END DO
       END SUBROUTINE arnoldi_adjoint

       !===============================================
       !  ARNOLDI METHOD FOR ADJOINT (second version)
       !===============================================
       SUBROUTINE arnoldi_adj2(kry, u, mod_u, v, h, u_b, ub, wb, time)
         USE global_variables
         USE RHS_routines
         IMPLICIT NONE
         INCLUDE "mpif.h"

         INTEGER, INTENT(IN) :: kry
         REAL(pr), INTENT(IN) :: mod_u
         REAL(pr), DIMENSION(n_dim), INTENT(IN) :: u, u_b
         REAL(pr), DIMENSION(n_dim, kry), INTENT(OUT) :: v 
         REAL(pr), DIMENSION(kry, kry), INTENT (OUT) :: h
         REAL(pr), DIMENSION(1:n(1),1:local_nlast,1:2), INTENT(IN) :: ub
         REAL(pr), DIMENSION(1:n(1),1:local_nlast), INTENT(IN) :: wb
         REAL(pr), INTENT(IN) :: time

         INTEGER :: i, j
         REAL(pr), DIMENSION (n_dim) :: w
         REAL(pr) :: a

         h = 0.0_pr
         w = 0.0_pr
         v = 0.0_pr
         v(:,1) = u / mod_u 
         DO j = 1, kry
            CALL rhs_adj_kry(v(:,j), w, time)
            DO i = 1, j
               a = DOT_PRODUCT (v(:,i), w) 
               CALL MPI_ALLREDUCE(a,h(i,j),1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)     
               w = w - v(:,i) * h(i,j) 
            END DO
            IF (j < kry) THEN
               a = DOT_PRODUCT(w, w)
               CALL MPI_ALLREDUCE(a,h(j+1,j),1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
               h(j+1,j)=SQRT(h(j+1,j))
               v(:,j+1) = w / h(j+1,j)
            END IF
         END DO
       END SUBROUTINE arnoldi_adj2

       !===============================================
       !  ARNOLDI METHOD FOR PERTURBATION SYSTEM
       !===============================================
       SUBROUTINE arnoldi_d1nse(kry, u, mod_u, v, h, u_b, ub, wb, time)
         USE global_variables
         USE RHS_routines
         IMPLICIT NONE
         INCLUDE "mpif.h"

         INTEGER, INTENT(IN) :: kry
         REAL(pr), INTENT(IN) :: mod_u
         REAL(pr), DIMENSION(n_dim), INTENT(IN) :: u, u_b
         REAL(pr), DIMENSION(n_dim, kry), INTENT(OUT) :: v 
         REAL(pr), DIMENSION(kry, kry), INTENT (OUT) :: h
         REAL(pr), DIMENSION(:,:,:), INTENT(IN) :: ub
         REAL(pr), DIMENSION(:,:), INTENT(IN) :: wb
         REAL(pr), INTENT(IN) :: time

         INTEGER :: i, j
         REAL(pr), DIMENSION (n_dim) :: w
         REAL(pr) :: a

         h = 0.0_pr
         w = 0.0_pr
         v = 0.0_pr
         v(:,1) = u / mod_u 
         DO j = 1, kry
            CALL rhs_d1nse_kry(v(:,j), w, time)
            DO i = 1, j
               a = DOT_PRODUCT (v(:,i), w) 
               CALL MPI_ALLREDUCE(a,h(i,j),1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)     
               w = w - v(:,i) * h(i,j) 
            END DO
            IF (j < kry) THEN
               a = DOT_PRODUCT(w, w)
               CALL MPI_ALLREDUCE(a,h(j+1,j),1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
               h(j+1,j)=SQRT(h(j+1,j))
               v(:,j+1) = w / h(j+1,j)
            END IF
         END DO
       END SUBROUTINE arnoldi_d1nse


       !=================================================
       ! INTEGRAL
       !=================================================
       SUBROUTINE integral (myInt, lambda, j, dt)
         USE global_variables
         IMPLICIT NONE
         INTEGER :: i
         INTEGER , INTENT (IN) :: j
         REAL (pr), INTENT (IN) :: dt
         COMPLEX (pr), INTENT (OUT) :: myInt
         COMPLEX (pr),INTENT (IN) :: lambda
         COMPLEX (pr), DIMENSION (0:j) :: R

         R(0) = 1.0_pr/lambda * (EXP (lambda * dt) - 1.0_pr)
         DO i = 1, j
            R(i)= 1.0_pr/lambda * (REAL(i, pr) * R(i-1)- dt**i)
         END DO
         myInt = R(j)
       END SUBROUTINE integral


       !====================================
       ! LEAST_SQUARE
       !====================================
       SUBROUTINE  least_square (Ndim, Nobs, NCOF, tau, REM, BB)
         USE global_variables
         IMPLICIT NONE
         INTEGER, INTENT (IN)  :: Ndim, Nobs, NCOF
         REAL (pr), DIMENSION (Nobs),       INTENT (IN)  :: tau
         REAL (pr), DIMENSION (Ndim,Nobs),  INTENT (IN)  :: REM
         REAL (pr), DIMENSION (Ndim,Ncof),  INTENT (OUT) :: BB

         INTEGER                     :: i, nn, Ndeg
         REAL (pr), DIMENSION (Nobs) :: Xdata
         REAL (pr), DIMENSION (Nobs) :: Ydata
         REAL (pr)                   :: det, det1, det2

         Ndeg = ncof+1

         DO i = 1, Ndim
            DO nn = 1, Nobs
               Xdata(nn) = tau(nn)
               Ydata(nn) = REM(i,nn)
            END DO
            det  = xdata(1)**2 * xdata(2)**3 - xdata(2)**2 * xdata(1)**3
            det1 = ydata(1)    * xdata(2)**3 - ydata(2)    * xdata(1)**3
            det2 = xdata(1)**2 * ydata(2)    - xdata(2)**2 * ydata(1)

            BB(i,1) = det1/det
            BB(i,2) = det2/det
         END DO
       END SUBROUTINE least_square


       !=============================================
       ! INTEG
       !=============================================
       FUNCTION integ (hess, j, dt, flag)
         USE global_variables
         IMPLICIT NONE
         INTEGER :: i
         INTEGER, INTENT (IN) :: flag, j
         REAL (pr), INTENT (IN) :: dt
         REAL (pr), DIMENSION (kry_m, kry_m), INTENT (IN) :: hess
         REAL (pr), DIMENSION (kry_m, kry_m) :: integ
         COMPLEX (pr), DIMENSION (kry_m, kry_m) :: evec, evec_inv
         COMPLEX (pr), DIMENSION (kry_m)   :: eval, stok

         !--Variables for linear algebra routines
         INTEGER :: info
         INTEGER, DIMENSION (kry_m) :: ipvt
         REAL (pr), DIMENSION (34*kry_m) :: work
         REAL (pr), DIMENSION (kry_m) :: wr, wi
         REAL (pr), DIMENSION (kry_m, kry_m) :: h_tmp, vl, vr
         COMPLEX (pr), DIMENSION (2) :: det
         COMPLEX (pr), DIMENSION (kry_m) :: work1

         h_tmp = hess
         CALL dgeev ('N', 'V', kry_m, h_tmp, kry_m, wr, wi, vl, 1, vr, kry_m, work, 34*kry_m, info)
         eval = CMPLX (wr, wi)

         DO i = 1, kry_m
           IF (wi(i) == 0.0_pr) THEN
             evec(:,i) = CMPLX (vr(:,i), 0.0_pr)
           ELSE IF (wi(i) > 0.0_pr) THEN
             evec(:,i)   = CMPLX (vr(:,i),  vr(:,i+1))
             evec(:,i+1) = CMPLX (vr(:,i), -vr(:,i+1))
           END IF
         END DO  

         DO i = 1, kry_m
            CALL integral (stok(i), eval(i), j, dt)
         END DO

         evec_inv = evec
         CALL zgetrf (kry_m, kry_m, evec_inv, kry_m, ipvt, info)
         CALL zgetri (kry_m, evec_inv, kry_m, ipvt, work1, 34*kry_m, info)

         !--Calculate exponential matrix 
         IF (flag == 1) integ = REAL (MATMUL (evec, MATMUL (diag_mat (kry_m, stok), evec_inv)), pr)  
         IF (flag == 2) integ = REAL (MATMUL (evec, MATMUL (diag_mat (kry_m, eval*stok), evec_inv)), pr)
       END FUNCTION integ


       !========================================
       ! RMD_T
       !========================================
       FUNCTION rmd_t (u, u_b, rhs_u, ub, wb)
         USE RHS_routines    
         IMPLICIT NONE  
  
         REAL (pr), DIMENSION (n_dim), INTENT (IN) :: u, u_b
         REAL (pr), DIMENSION (n_dim), INTENT (IN) :: rhs_u
         REAL (pr), DIMENSION (n_dim)              :: rmd_t  
         REAL (pr), DIMENSION (1:n(1),1:local_nlast,1:2), INTENT (IN) :: ub
         REAL (pr), DIMENSION (1:n(1),1:local_nlast), INTENT (IN) :: wb
      
         REAL (pr), DIMENSION (n_dim) :: rhs1, drhs
  
         CALL drhs_kry (u, drhs, ub, wb)
         CALL rhs_kry  (u + u_b, rhs1)

         rmd_t = rhs1 - (rhs_u + drhs)   
       END FUNCTION rmd_t

       !========================================
       ! RMD_T ADJOINT
       !========================================
       FUNCTION rmd_t_adj (u, u_b, rhs_u, ub, wb, time)
         USE RHS_routines    
         IMPLICIT NONE  
  
         REAL (pr), DIMENSION (n_dim), INTENT (IN) :: u, u_b
         REAL (pr), DIMENSION (n_dim), INTENT (IN) :: rhs_u
         REAL (pr), DIMENSION (n_dim)              :: rmd_t_adj  
         REAL (pr), DIMENSION (1:n(1),1:local_nlast,1:2), INTENT (IN) :: ub
         REAL (pr), DIMENSION (1:n(1),1:local_nlast), INTENT (IN) :: wb
         REAL (pr), INTENT(IN) :: time
         REAL (pr), DIMENSION (n_dim) :: rhs1, drhs
  
         CALL rhs_adj_kry (u, drhs, time)
         CALL rhs_adj_kry (u + u_b, rhs1, time)

         rmd_t_adj = rhs1 - (rhs_u + drhs)   
       END FUNCTION rmd_t_adj


       !============================================
       ! DIAG_MAT
       !============================================
       FUNCTION diag_mat (m, val)
         USE global_variables
         IMPLICIT NONE
         INTEGER :: i
         INTEGER, INTENT (IN) :: m
         COMPLEX (pr), DIMENSION (m), INTENT (IN)   :: val
         COMPLEX (pr), DIMENSION (m,m) :: diag_mat

         diag_mat = CMPLX (0.0_pr, 0.0_pr)
         DO i = 1, m
            diag_mat(i,i) = val(i)
         END DO
       END FUNCTION diag_mat


       !===========================================
       ! U_VEC
       !===========================================
       FUNCTION u_vec (m, val)
         USE global_variables
         IMPLICIT NONE
         INTEGER :: i
         INTEGER, INTENT (IN) :: m
         REAL (pr), INTENT (IN) :: val
         REAL (pr), DIMENSION (m) :: u_vec

         u_vec (1) = val
         u_vec (2:m) = 0.0_pr
       END FUNCTION u_vec

END MODULE


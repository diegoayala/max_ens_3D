MODULE RHS_routines

     CONTAINS
       !==========================================
       !     RHS_KRY
       !==========================================
       SUBROUTINE rhs_kry (u_kry, rhs_u) 
         USE global_variables
         USE FFT2
         USE data_ops
         USE function_ops
         IMPLICIT NONE
         INCLUDE 'mpif.h'
         REAL(pr), DIMENSION (n_dim), INTENT(IN) :: u_kry
         REAL(pr), DIMENSION (n_dim), INTENT(OUT) :: rhs_u

         INTEGER :: i
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: rhs, u
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: w
         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: fdiss, frhs, fu
         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: fw

         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: aux, faux
 
         !--Store the velocity in a matrix
         CALL vecmat (u_kry, u)

         !--Transform to  Fourier space
         DO i=1,2
            aux = CMPLX(u(:,:,i),0.0_pr)
            CALL ffourier(aux, faux)
            fu(:,:,i) = faux
         END DO

         !--Calculate the dissipation in Fourier space
         CALL cal_diss(fu, fdiss)

         !--Calculate vorticity in Fourier space
         CALL cal_w(fu, fw)
 
         !--Transform to physical space
         CALL bfourier(fw, aux)
         w = REAL(aux)

         !--Calculate rhs uxw, and penalization terms in physical space
         rhs(:,:,1) = u(:,:,2)*w
         rhs(:,:,2) = -u(:,:,1)*w

         !--Transform to Fourier space
         DO i=1,2
            aux = CMPLX(rhs(:,:,i),0.0_pr)
            CALL ffourier(aux, faux)
            CALL dealiasing(faux)
            frhs(:,:,i) = faux
         END DO
         
         !--Add dissipation to rhs
         frhs = frhs + fdiss
 
         !--Make rhs divergence free
         CALL div_free(frhs)
 
         !--Transform rhs to physical space
         DO i = 1, 2
            faux = frhs(:,:,i)
            CALL bfourier(faux, aux)
            rhs(:,:,i) = REAL(aux)
         END DO

         !--Store the rhs in a vector
         CALL matvec (rhs, rhs_u) 
       END SUBROUTINE rhs_kry 
 
       !============================================
       !    DRHS_KRY Derivative of right hand side
       !============================================     
       SUBROUTINE drhs_kry(u_kry, drhs_u, u0, w0)
         USE global_variables
         USE data_ops
         USE function_ops
         USE FFT2
         IMPLICIT NONE
         INCLUDE "mpif.h"
         
         REAL(pr), DIMENSION(n_dim), INTENT(IN) :: u_kry
         REAL(pr), DIMENSION(n_dim), INTENT(OUT) :: drhs_u
         REAL(pr), DIMENSION(1:n(1),1:local_nlast,1:2), INTENT(IN) :: u0
         REAL(pr), DIMENSION(1:n(1),1:local_nlast), INTENT(IN) :: w0
 
         INTEGER :: i
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: drhs, u
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: w
         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: fdiss, fdrhs, fu 
         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: fwb, fw

         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: aux, faux  
   
         !--Store the velocity in a matrix
         CALL vecmat (u_kry, u)

         !--Transform to  Fourier space
         DO i = 1, 2
            aux = CMPLX(u(:,:,i),0.0_pr)
            CALL ffourier(aux, faux)
            fu(:,:,i) = faux
         END DO

         !--Calculate the dissipation in Fourier space
         CALL cal_diss(fu, fdiss)

         !--Calculate vorticity in Fourier space
         CALL cal_w(fu, fw)

         !--Transform vorticity to physical space
         CALL bfourier(fw, faux)
         w = REAL(faux)
       
         !--Calculate rhs: uxw and penalization terms in physical space
         drhs(:,:,1) = u(:,:,2)*w0 + u0(:,:,2)*w
         drhs(:,:,2) = -u(:,:,1)*w0 - u0(:,:,1)*w

         !--Transform to Fourier space
         DO i = 1, 2
            aux = drhs(:,:,i)
            CALL ffourier(aux, faux)
            CALL dealiasing(faux)
            fdrhs(:,:,i) = faux
         END DO
  
         !--Add dissipation to rhs
         fdrhs = fdrhs + fdiss
  
         !--Make drhs divergence free
         CALL div_free(fdrhs)
  
         !--Set mean of drhs = 0
         IF (rank==0) fdrhs(1,1,:) = CMPLX(0.0_pr,0.0_pr)

         !--Transform to physical space
         DO i = 1, 2
            faux = fdrhs(:,:,i)
            CALL bfourier (faux, aux)
            drhs(:,:,i) = REAL(aux)
         END DO
  
         !--Store the rhs in a vector
         CALL matvec (drhs, drhs_u) 
       END SUBROUTINE drhs_kry

       !==========================================
       !     RHS_ADJ_KRY
       !==========================================
       SUBROUTINE rhs_adj_kry (u_kry, rhs_u, time) 
         USE global_variables
         USE FFT2
         USE data_ops
         USE function_ops
         IMPLICIT NONE
         INCLUDE 'mpif.h'
         REAL(pr), DIMENSION (n_dim), INTENT(IN) :: u_kry
         REAL(pr), DIMENSION (n_dim), INTENT(OUT) :: rhs_u
         REAL(pr), INTENT(IN) :: time

         INTEGER :: i, i1, i2
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: rhs, u, NS_u, diss
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: w, NS_w, real_aux
         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: fdiss, frhs, fu
         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: fw

         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: aux, faux
         REAL(pr) :: mode 

         !--Store the velocity in a matrix
         CALL vecmat (u_kry, u)

         !--Transform to Fourier Space
         DO i=1,2
            aux = CMPLX(u(:,:,i),0.0_pr)
            CALL ffourier(aux,faux)
            fu(:,:,i) = faux
         END DO 

         !--Calculate the dissipation in fourier space
         CALL cal_diss(fu, fdiss)

         !--Calculate rhs rot(u x NS_u) + NS_w x u, and penalization terms in physical space
         CALL read_velocity(NS_u, T-time, "nse")
         CALL cal_vort(NS_u, NS_w)
         
         !CALL read_vorticity(NS_w, T-time, "nse")
         !aux = CMPLX(NS_w, 0.0_pr)
         !CALL ffourier(aux, faux)
         !CALL filter(faux)
         !CALL bfourier(faux, aux)
         !NS_w = REAL(aux)
         !DO i=1,2
         !   aux = CMPLX(NS_u(:,:,i),0.0_pr)
         !   CALL ffourier(aux, faux)
         !   CALL filter(faux)
         !   CALL bfourier(faux, aux)
         !   NS_u(:,:,i) = REAL(aux)
         !END DO
 
         real_aux = u(:,:,1)*NS_u(:,:,2) - u(:,:,2)*NS_u(:,:,1)
         CALL grad_perp(real_aux, rhs)
         rhs(:,:,1) = rhs(:,:,1) - u(:,:,2)*NS_w
         rhs(:,:,2) = rhs(:,:,2) + u(:,:,1)*NS_w

         !--Transform to Fourier space
         DO i=1,2
            aux = CMPLX(rhs(:,:,i),0.0_pr)
            CALL ffourier(aux, faux)
            CALL dealiasing(faux)
            frhs(:,:,i) = faux
         END DO
         
         !--Add dissipation to rhs
         frhs = frhs + fdiss
 
         !--Make rhs divergence free
         CALL div_free(frhs)
 
         !--Transform rhs to physical space
         DO i = 1, 2
            faux = frhs(:,:,i)
            CALL bfourier(faux, aux)
            rhs(:,:,i) = REAL(aux)
         END DO

         !--Store the rhs in a vector
         CALL matvec (rhs, rhs_u) 
 
      END SUBROUTINE rhs_adj_kry 

      !============================================
      ! Right hand side for perturbation system
      !============================================     
      SUBROUTINE rhs_d1nse_kry(u_kry, rhs_u, time)
         USE global_variables
         USE FFT2
         USE data_ops
         USE function_ops
         IMPLICIT NONE
         INCLUDE 'mpif.h'
         REAL(pr), DIMENSION (n_dim), INTENT(IN) :: u_kry
         REAL(pr), DIMENSION (n_dim), INTENT(OUT) :: rhs_u
         REAL(pr), INTENT(IN) :: time

         INTEGER :: i, i1, i2
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: rhs, u, NS_u, diss
         REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: w, NS_w, real_aux
         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: fdiss, frhs, fu
         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: fw

         COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: aux, faux
         REAL(pr) :: mode 

         !--Store the velocity in a matrix
         CALL vecmat (u_kry, u)

         !--Calculate vorticity
         CALL cal_vort(u,w)
         
         !--Transform to Fourier Space
         DO i=1,2
            aux = CMPLX(u(:,:,i),0.0_pr)
            CALL ffourier(aux,faux)
            fu(:,:,i) = faux
         END DO 

         !--Calculate the dissipation in fourier space
         CALL cal_diss(fu, fdiss)

         !--Calculate rhs rot(u x NS_u) + NS_w x u, and penalization terms in physical space
         CALL read_velocity(NS_u, time, "nse")
         CALL cal_vort(NS_u, NS_w)
         
         real_aux = u(:,:,1)*NS_u(:,:,1) + u(:,:,2)*NS_u(:,:,2)
         CALL gradient(real_aux, rhs)
         rhs(:,:,1) = u(:,:,2)*NS_w + NS_u(:,:,2)*w - rhs(:,:,1)
         rhs(:,:,2) = -u(:,:,1)*NS_w - NS_u(:,:,1)*w - rhs(:,:,2)

         !--Transform to Fourier space
         DO i=1,2
            aux = CMPLX(rhs(:,:,i),0.0_pr)
            CALL ffourier(aux, faux)
            CALL dealiasing(faux)
            frhs(:,:,i) = faux
         END DO
         
         !--Add dissipation to rhs
         frhs = frhs + fdiss
 
         !--Make rhs divergence free
         CALL div_free(frhs)
 
         !--Transform rhs to physical space
         DO i = 1, 2
            faux = frhs(:,:,i)
            CALL bfourier(faux, aux)
            rhs(:,:,i) = REAL(aux)
         END DO

         !--Store the rhs in a vector
         CALL matvec (rhs, rhs_u) 

      END SUBROUTINE rhs_d1nse_kry

END MODULE

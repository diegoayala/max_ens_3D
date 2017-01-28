!======================================================
! MODULE CONTAINS ROUTINES REPRESENTING OPERATIONS 
! APPLIED TO FUNCTIONS.
!
! (*) initial_guess
! (*) vort2vel
! (*) palinstrophy
! (*) cal_diss
! (*) cal_w
! (*) and more...
!=======================================================

MODULE function_ops
  
  IMPLICIT NONE

  CONTAINS
        !===================================
        !  Initial guess
        !===================================
        SUBROUTINE initial_guess
          USE global_variables
          USE data_ops
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:3) :: dx, local_K, K
          INTEGER :: ii, jj, kk, nn,  size_seed
          INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          REAL(pr) :: X, Y, Z

          LOGICAL :: read_from_file
          CHARACTER(100) :: filename
          CHARACTER(2) :: K0txt, E0txt, IGtxt, Fx_txt, Fy_txt, Fz_txt
          CHARACTER(4) :: Ntxt

          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: Ux, Uy, Uz
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: Wx, Wy, Wz

          dx = 1.0_pr/REAL(n, pr)
          ALLOCATE( Ux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( Uy(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( Uz(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( Wx(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( Wy(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( Wz(1:n(1),1:n(2),1:local_nlast) )
          
          SELECT CASE (iguess)

            CASE (0) ! Corresponds to solutions in the Taylor-Green Vortex family with |k|^2 = 3
              WRITE(K0txt,'(i2.2)') K0_index
              WRITE(E0txt,'(i2.2)') E0_index
              WRITE(IGtxt,'(i2.2)') iguess+1

              filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
              
              Fx_txt = "Ux"
              Fy_txt = "Uy"
              Fz_txt = "Uz" 

              read_from_file = .TRUE.

            CASE (1) ! Corresponds to solutions in the Taylor-Green Vortex family, from different value of Enstrophy
              WRITE(K0txt,'(i2.2)') K0_index
              WRITE(E0txt,'(i2.2)') E0_index+1
              WRITE(IGtxt,'(i2.2)') iguess

              filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
              
              Fx_txt = "Ux"
              Fy_txt = "Uy"
              Fz_txt = "Uz" 

              read_from_file = .TRUE.

            CASE (2) ! Corresponds to solutions in the Taylor-Green Vortex family with |k|^2 = 1 
              WRITE(K0txt,'(i2.2)') K0_index
              WRITE(E0txt,'(i2.2)') E0_index
              WRITE(IGtxt,'(i2.2)') iguess

              filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
              
              Fx_txt = "Ux"
              Fy_txt = "Uy"
              Fz_txt = "Uz" 

              read_from_file = .TRUE.

            CASE (3) ! Corresponds to silutions in the Ayala-Protas Vortex family, from different value of Enstrophy
              WRITE(K0txt,'(i2.2)') K0_index
              WRITE(E0txt,'(i2.2)') E0_index
              WRITE(IGtxt,'(i2.2)') iguess

              filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
              
              Fx_txt = "Ux"
              Fy_txt = "Uy"
              Fz_txt = "Uz" 

              read_from_file = .TRUE.

            CASE (4) ! Corresponds to solutions in the Taylor-Green Vortex family with |k|^2 = 2
              WRITE(K0txt,'(i2.2)') K0_index
              WRITE(E0txt,'(i2.2)') E0_index-1
              WRITE(IGtxt,'(i2.2)') iguess

              filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
              
              Fx_txt = "Ux"
              Fy_txt = "Uy"
              Fz_txt = "Uz" 

              read_from_file = .TRUE.

            CASE (5) ! Corresponds to solutions in the Taylor-Green Vortex family with |k|^2 = 2, from different value of enstrophy
              WRITE(K0txt,'(i2.2)') K0_index
              WRITE(E0txt,'(i2.2)') E0_index-1
              WRITE(IGtxt,'(i2.2)') iguess

              filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
              
              Fx_txt = "Ux"
              Fy_txt = "Uy"
              Fz_txt = "Uz" 

              read_from_file = .TRUE.


            CASE (10) !Initial guess is Taylor-Green Vortex
              DO kk=1, local_nlast 
                 DO jj=1, n(2)
                    DO ii=1, n(1)
                       X = REAL(ii-1,pr)*dx(1)
                       Y = REAL(jj-1,pr)*dx(2)
                       Z = REAL(kk-1,pr)*dx(3)
                       Ux(ii,jj,kk) = SIN(2.0_pr*PI*X)*COS(2.0_pr*PI*Y)*COS(2.0_pr*PI*Z)
                       Uy(ii,jj,kk) = -COS(2.0_pr*PI*X)*SIN(2.0_pr*PI*Y)*COS(2.0_pr*PI*Z)
                       Uz(ii,jj,kk) = 0.0_pr
                    END DO
                 END DO
              END DO

              Uvec(:,:,:,1) = Ux
              Uvec(:,:,:,2) = Uy
              Uvec(:,:,:,3) = Uz

              CALL vel2vort(Uvec,Wvec)
              Wx = Wvec(:,:,:,1)
              Wy = Wvec(:,:,:,2)
              Wz = Wvec(:,:,:,3)

              read_from_file = .FALSE.

              IF (save_data_Optim) THEN
                 WRITE(K0txt, '(i2.2)') K0_index
                 WRITE(E0txt, '(i2.2)') E0_index
                 WRITE(IGtxt, '(i2.2)') iguess
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
                 CALL save_field_R3toR3_ncdf(Ux, Uy, Uz, "Ux", "Uy", "Uz", filename, "netCDF")
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_w0.nc"
                 CALL save_field_R3toR3_ncdf(Wx, Wy, Wz, "Wx", "Wy", "Wz", filename, "netCDF")
              END IF 

            CASE (11) !Initial guess is a perturbation from Taylor-Green vortex
              DO kk=1, local_nlast
                 DO jj=1, n(2)
                    DO ii=1, n(1)
                       X = REAL(ii-1,pr)*dx(1)
                       Y = REAL(jj-1,pr)*dx(2)
                       Z = REAL(kk-1,pr)*dx(3)
                       Ux(ii,jj,kk) = SIN(2.0_pr*PI*X)*COS(2.0_pr*PI*Y)*COS(2.0_pr*PI*Z)
                       Uy(ii,jj,kk) = COS(2.0_pr*PI*X)*SIN(2.0_pr*PI*Y)*COS(2.0_pr*PI*Z)
                       Uz(ii,jj,kk) = -2.0_pr*COS(2.0_pr*PI*X)*COS(2.0_pr*PI*Y)*SIN(2.0_pr*PI*Z)
                    END DO
                 END DO
              END DO
             
              ! ADD PERTURBATION TO TAYLOR-GREEN VORTEX
              DO nn = 2,10
                 DO kk=1, local_nlast 
                    DO jj=1, n(2)
                       DO ii=1, n(1)
                          X = REAL(ii-1,pr)*dx(1)
                          Y = REAL(jj-1,pr)*dx(2)
                          Z = REAL(kk-1,pr)*dx(3)
                          Ux(ii,jj,kk) = Ux(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*SIN(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                          Uy(ii,jj,kk) = Uy(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*COS(2.0_pr*PI*REAL(nn,pr)*Z)
                          Uz(ii,jj,kk) = Uz(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*SIN(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                       END DO
                    END DO
                 END DO
              END DO
               
              Uvec(:,:,:,1) = Ux
              Uvec(:,:,:,2) = Uy
              Uvec(:,:,:,3) = Uz

              CALL div_free(Uvec)   
              
              CALL vel2vort(Uvec,Wvec)
              Wx = Wvec(:,:,:,1)
              Wy = Wvec(:,:,:,2)
              Wz = Wvec(:,:,:,3)

              read_from_file = .FALSE.

              IF (save_data_Optim) THEN
                 WRITE(K0txt, '(i2.2)') K0_index
                 WRITE(E0txt, '(i2.2)') E0_index
                 WRITE(IGtxt, '(i2.2)') iguess
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
                 CALL save_field_R3toR3_ncdf(Ux, Uy, Uz, "Ux", "Uy", "Uz", filename, "netCDF")
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_w0.nc"
                 CALL save_field_R3toR3_ncdf(Wx, Wy, Wz, "Wx", "Wy", "Wz", filename, "netCDF")
              END IF 

           CASE (20) 
              ! Initial guess is Ayala-Protas Vortex ( k_vec = (1,0,0) and
              ! permutations, i.e. |k|^2 = 1 )
              DO kk=1, local_nlast 
                 DO jj=1, n(2)
                    DO ii=1, n(1)
                       X = REAL(ii-1,pr)*dx(1)
                       Y = REAL(jj-1,pr)*dx(2)
                       Z = REAL(kk-1,pr)*dx(3)
                       Ux(ii,jj,kk) = 4.0_pr*COS(PI*(Y + Z))*COS(PI*(Y - Z))
                       Uy(ii,jj,kk) = 4.0_pr*COS(PI*(X + Z))*COS(PI*(X - Z)) 
                       Uz(ii,jj,kk) = 4.0_pr*COS(PI*(X + Y))*COS(PI*(X - Y))
                    END DO
                 END DO
              END DO

              Uvec(:,:,:,1) = Ux
              Uvec(:,:,:,2) = Uy
              Uvec(:,:,:,3) = Uz

              CALL vel2vort(Uvec,Wvec)
              Wx = Wvec(:,:,:,1)
              Wy = Wvec(:,:,:,2)
              Wz = Wvec(:,:,:,3)

              read_from_file = .FALSE.

              IF (save_data_Optim) THEN
                 WRITE(K0txt, '(i2.2)') K0_index
                 WRITE(E0txt, '(i2.2)') E0_index
                 WRITE(IGtxt, '(i2.2)') iguess
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
                 CALL save_field_R3toR3_ncdf(Ux, Uy, Uz, "Ux", "Uy", "Uz", filename, "netCDF")
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_w0.nc"
                 CALL save_field_R3toR3_ncdf(Wx, Wy, Wz, "Wx", "Wy", "Wz", filename, "netCDF")
              END IF 

            CASE (21) 
              ! Initial guess is a perturbation of Ayala-Protas Vortex 
              ! ( k_vec = (1,0,0) and permutations, i.e. |k|^2 = 1 )
              DO kk=1, local_nlast 
                 DO jj=1, n(2)
                    DO ii=1, n(1)
                       X = REAL(ii-1,pr)*dx(1)
                       Y = REAL(jj-1,pr)*dx(2)
                       Z = REAL(kk-1,pr)*dx(3)
                       Ux(ii,jj,kk) = 4.0_pr*COS(PI*(Y + Z))*COS(PI*(Y - Z))
                       Uy(ii,jj,kk) = 4.0_pr*COS(PI*(X + Z))*COS(PI*(X - Z)) 
                       Uz(ii,jj,kk) = 4.0_pr*COS(PI*(X + Y))*COS(PI*(X - Y))
                    END DO
                 END DO
              END DO
             
              ! ADD PERTURBATION TO AYALA-PROTAS VORTEX
              DO nn = 2,10
                 DO kk=1, local_nlast 
                    DO jj=1, n(2)
                       DO ii=1, n(1)
                          X = REAL(ii-1,pr)*dx(1)
                          Y = REAL(jj-1,pr)*dx(2)
                          Z = REAL(kk-1,pr)*dx(3)
                          Ux(ii,jj,kk) = Ux(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*SIN(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                          Uy(ii,jj,kk) = Uy(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*COS(2.0_pr*PI*REAL(nn,pr)*Z)
                          Uz(ii,jj,kk) = Uz(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*SIN(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                       END DO
                    END DO
                 END DO
              END DO
               
              Uvec(:,:,:,1) = Ux
              Uvec(:,:,:,2) = Uy
              Uvec(:,:,:,3) = Uz

              CALL div_free(Uvec)   
              
              CALL vel2vort(Uvec,Wvec)
              Wx = Wvec(:,:,:,1)
              Wy = Wvec(:,:,:,2)
              Wz = Wvec(:,:,:,3)

              read_from_file = .FALSE.

              IF (save_data_Optim) THEN
                 WRITE(K0txt, '(i2.2)') K0_index
                 WRITE(E0txt, '(i2.2)') E0_index
                 WRITE(IGtxt, '(i2.2)') iguess
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
                 CALL save_field_R3toR3_ncdf(Ux, Uy, Uz, "Ux", "Uy", "Uz", filename, "netCDF")
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_w0.nc"
                 CALL save_field_R3toR3_ncdf(Wx, Wy, Wz, "Wx", "Wy", "Wz", filename, "netCDF")
              END IF 

           CASE (30) 
              ! Initial guess is Ayala-Protas Vortex ( k_vec = (1,1,0) and
              ! permutations, i.e. |k|^2 = 2 )
              DO kk=1,local_nlast
                 DO jj=1, n(2)
                    DO ii=1, n(1)
                       X = REAL(ii-1,pr)*dx(1)
                       Y = REAL(jj-1,pr)*dx(2)
                       Z = REAL(kk-1,pr)*dx(3)
                       Ux(ii,jj,kk) = 4.0_pr*COS(2.0_pr*PI*Y)*COS(2.0_pr*PI*Z)
                       Uy(ii,jj,kk) = 4.0_pr*COS(2.0_pr*PI*X)*COS(2.0_pr*PI*Z) 
                       Uz(ii,jj,kk) = 4.0_pr*COS(2.0_pr*PI*X)*COS(2.0_pr*PI*Y)
                    END DO
                 END DO
              END DO

              Uvec(:,:,:,1) = Ux
              Uvec(:,:,:,2) = Uy
              Uvec(:,:,:,3) = Uz

              CALL vel2vort(Uvec,Wvec)
              Wx = Wvec(:,:,:,1)
              Wy = Wvec(:,:,:,2)
              Wz = Wvec(:,:,:,3)

              read_from_file = .FALSE.

              IF (save_data_Optim) THEN
                 WRITE(K0txt, '(i2.2)') K0_index
                 WRITE(E0txt, '(i2.2)') E0_index
                 WRITE(IGtxt, '(i2.2)') iguess
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
                 CALL save_field_R3toR3_ncdf(Ux, Uy, Uz, "Ux", "Uy", "Uz", filename, "netCDF")
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_w0.nc"
                 CALL save_field_R3toR3_ncdf(Wx, Wy, Wz, "Wx", "Wy", "Wz", filename, "netCDF")
              END IF 

            CASE (31) 
              ! Initial guess is a perturbation of Ayala-Protas Vortex 
              ! ( k_vec = (1,1,0) and permutations, i.e. |k|^2 = 2 )
              DO kk=1, local_nlast 
                 DO jj=1, n(2)
                    DO ii=1, n(1)
                       X = REAL(ii-1,pr)*dx(1)
                       Y = REAL(jj-1,pr)*dx(2)
                       Z = REAL(kk-1,pr)*dx(3)
                       Ux(ii,jj,kk) = 4.0_pr*COS(2.0_pr*PI*Y)*COS(2.0_pr*PI*Z)
                       Uy(ii,jj,kk) = 4.0_pr*COS(2.0_pr*PI*X)*COS(2.0_pr*PI*Z) 
                       Uz(ii,jj,kk) = 4.0_pr*COS(2.0_pr*PI*X)*COS(2.0_pr*PI*Y)
                    END DO
                 END DO
              END DO
             
              ! ADD PERTURBATION TO AYALA-PROTAS VORTEX
              DO nn = 2,10
                 DO kk=1, local_nlast 
                    DO jj=1, n(2)
                       DO ii=1, n(1)
                          X = REAL(ii-1,pr)*dx(1)
                          Y = REAL(jj-1,pr)*dx(2)
                          Z = REAL(kk-1,pr)*dx(3)
                          Ux(ii,jj,kk) = Ux(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*SIN(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                          Uy(ii,jj,kk) = Uy(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*COS(2.0_pr*PI*REAL(nn,pr)*Z)
                          Uz(ii,jj,kk) = Uz(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*SIN(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                       END DO
                    END DO
                 END DO
              END DO
               
              Uvec(:,:,:,1) = Ux
              Uvec(:,:,:,2) = Uy
              Uvec(:,:,:,3) = Uz

              CALL div_free(Uvec)   
              
              CALL vel2vort(Uvec,Wvec)
              Wx = Wvec(:,:,:,1)
              Wy = Wvec(:,:,:,2)
              Wz = Wvec(:,:,:,3)

              read_from_file = .FALSE.

              IF (save_data_Optim) THEN
                 WRITE(K0txt, '(i2.2)') K0_index
                 WRITE(E0txt, '(i2.2)') E0_index
                 WRITE(IGtxt, '(i2.2)') iguess
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
                 CALL save_field_R3toR3_ncdf(Ux, Uy, Uz, "Ux", "Uy", "Uz", filename, "netCDF")
                 filename = "/gwork/ayalada/MAX_ENS_3D/K"//K0txt//"/E"//E0txt//"/maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_w0.nc"
                 CALL save_field_R3toR3_ncdf(Wx, Wy, Wz, "Wx", "Wy", "Wz", filename, "netCDF")
              END IF 

            CASE DEFAULT
              DO kk=1, local_nlast 
                 DO jj=1, n(2)
                    DO ii=1, n(1)
                       Ux(ii,jj,kk) = 0.0_pr
                       Uy(ii,jj,kk) = 0.0_pr 
                       Uz(ii,jj,kk) = 0.0_pr
                    END DO
                 END DO
              END DO
              read_from_file = .FALSE.

          END SELECT

          IF (read_from_file) THEN

             CALL read_field_R3toR3_ncdf(Uvec, filename, Fx_txt, Fy_txt, Fz_txt) 
             !CALL vort2vel(Wvec, Uvec)            
             
             DO nn=1,3
                Ux = Uvec(:,:,:,nn)
                CALL dealiasing(Ux)
                Uvec(:,:,:,nn) = Ux
             END DO

             IF (add_pert) THEN
                DO nn = 10,20
                   DO kk=1, local_nlast 
                      DO jj=1, n(2)
                         DO ii=1, n(1)
                            X = REAL(ii-1,pr)*dx(1)
                            Y = REAL(jj-1,pr)*dx(2)
                            Z = REAL(kk-1,pr)*dx(3)
                            Ux(ii,jj,kk) = Ux(ii,jj,kk) + (1.0_pr/REAL(nn,pr))*SIN(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                            Uy(ii,jj,kk) = Uy(ii,jj,kk) + (1.0_pr/REAL(nn,pr))*COS(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*COS(2.0_pr*PI*REAL(nn,pr)*Z)
                            Uz(ii,jj,kk) = Uz(ii,jj,kk) + (1.0_pr/REAL(nn,pr))*COS(2.0_pr*PI*REAL(nn,pr)*X)*SIN(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                         END DO
                      END DO
                   END DO
                END DO
                Uvec(:,:,:,1) = Uvec(:,:,:,1) + Ux
                Uvec(:,:,:,2) = Uvec(:,:,:,2) + Uy
                Uvec(:,:,:,3) = Uvec(:,:,:,3) + Uz
             END IF
 
             CALL div_free(Uvec)   
             CALL vel2vort(Uvec, Wvec)

          END IF

          DEALLOCATE( Ux )
          DEALLOCATE( Uy )
          DEALLOCATE( Uz )
          DEALLOCATE( Wx )
          DEALLOCATE( Wy )
          DEALLOCATE( Wz )

        END SUBROUTINE initial_guess

        !===========================================================
        ! Initialize the random seed for RANDOM NUMBER GENERATOR
        !===========================================================
        SUBROUTINE init_random_seed()
           USE global_variables
           IMPLICIT NONE
           INCLUDE "mpif.h"

           INTEGER :: i, global_size_seed, size_seed, clock
           INTEGER, DIMENSION(:), ALLOCATABLE :: global_seed, local_seed
          
           CALL RANDOM_SEED(size = size_seed)
           global_size_seed = np*size_seed
           ALLOCATE(local_seed(size_seed))
           ALLOCATE(global_seed(global_size_seed))
           
           IF (rank==0) THEN
              OPEN(89, FILE='/dev/urandom', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
              READ(89) global_seed
              CLOSE(89)
           END IF
           CALL MPI_BCAST(global_seed, global_size_seed, MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo) 
           local_seed = global_seed(rank*size_seed+1:(rank+1)*size_seed)

           CALL RANDOM_SEED(put=local_seed)  
           
           DEALLOCATE(local_seed)
           DEALLOCATE(global_seed)

        END SUBROUTINE

        !===================================
        !  PERTURBATION FOR KAPPA TEST
        !===================================
        SUBROUTINE kappa_test_pert(phi_pert, mytype, m1, m2, m3)
          USE global_variables
          IMPLICIT NONE

          REAL(pr), DIMENSION(0:n(1)-1,0:n(2)-1,local_last_start:local_last_start+local_nlast-1,1:3), INTENT(OUT) :: phi_pert
          CHARACTER(len=*), INTENT(IN) :: mytype
          REAL(pr), INTENT(IN) :: m1, m2, m3
          
          REAL(pr), DIMENSION(1:3) :: dx
          INTEGER :: ii, jj, kk
          REAL(pr) :: X, Y, Z, ampl

          dx = 1.0/REAL(n, pr)
          SELECT CASE (mytype)

            CASE ("sine")
              DO kk=local_last_start, local_last_start+local_nlast-1 
                 DO jj=0, n(2)-1
                    DO ii=0, n(1)-1
                       X = REAL(ii,pr)*dx(1)
                       Y = REAL(jj,pr)*dx(2)
                       Z = REAL(kk,pr)*dx(3)
                       phi_pert(ii,jj,kk,1) = SIN(2.0_pr*PI*m1*X)*COS(2.0_pr*PI*m2*Y)*COS(2.0_pr*PI*m3*Z)
                       phi_pert(ii,jj,kk,2) = -COS(2.0_pr*PI*m1*X)*SIN(2.0_pr*PI*m2*Y)*COS(2.0_pr*PI*m3*Z)
                       phi_pert(ii,jj,kk,3) = 0.0_pr !COS(2.0_pr*PI*m1*X)*COS(2.0_pr*PI*m2*Y)*SIN(2.0_pr*PI*m3*Z)
                    END DO
                 END DO
              END DO

            !CASE ("gaussian")
            !DO j=local_last_start, local_last_start+local_nlast-1
            !   DO i=0,n(1)-1
            !      X = REAL(i,pr)*dx(1)
            !      Y = REAL(j,pr)*dx(2)
            !      w_pert(i,j) = EXP(-100.0_pr*(X-0.5_pr)**2 -100.0_pr*(Y-m1)**2 ) + &
            !                    EXP(-100.0_pr*(X-0.5_pr)**2 -100.0_pr*(Y-m2)**2 ) 
            !   END DO
            !END DO

           !CASE ("random")
            !CALL init_random_seed()
            !DO j=local_last_start, local_last_start+local_nlast-1
            !   DO i=0,n(1)-1
            !      !CALL init_random_seed()
            !      CALL RANDOM_NUMBER(X)
            !      Y = 2.0_pr*X - 1.0_pr
            !      vort0(i,j) = Y
            !   END DO
            !END DO
            !CALL filter(vort0, 2)

          
            !CASE DEFAULT
            !DO j=local_last_start, local_last_start+local_nlast-1
            !   DO i=0,n(1)-1
            !      w_pert(i,j) = 0.0
            !   END DO
            !END DO
          END SELECT

          !CALL vort2vel(w_pert, phi_pert)
          CALL div_free(phi_pert)
    
        END SUBROUTINE kappa_test_pert

        !=========================================================
        ! Calculate kinetic energy from function in physical space
        !=========================================================
        FUNCTION Energy(U) RESULT (kin_ener)
          USE global_variables
          IMPLICIT NONE
 
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: U
          INTEGER ::  i1, i2, i3
          REAL(pr), DIMENSION(1:3) :: kin_ener

          kin_ener = 0.0_pr

          DO i3=1,local_nlast
             DO i2=1,n(2) 
                DO i1=1,n(1)  
                   kin_ener(1) = kin_ener(1) + 0.5_pr*U(i1,i2,i3,1)**2*dV
                   kin_ener(2) = kin_ener(2) + 0.5_pr*U(i1,i2,i3,2)**2*dV
                   kin_ener(3) = kin_ener(3) + 0.5_pr*U(i1,i2,i3,3)**2*dV
                END DO
             END DO
          END DO
 
        END FUNCTION Energy

        !=========================================================
        ! Calculate enstrophy from function in physical space
        !=========================================================
        FUNCTION Enstrophy(U) RESULT (ens)
          USE global_variables
          IMPLICIT NONE
 
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: U

          REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: W

          REAL(pr), DIMENSION(1:3) :: ens
          INTEGER ::  i1, i2, i3
          REAL(pr) :: mode
          
          ALLOCATE( W(1:n(1),1:n(2),1:local_nlast,1:3) )

          CALL vel2vort(U,W)

          ens = 0.0_pr
          DO i3 = 1,local_nlast
             DO i2=1,n(2)
                DO i1 = 1,n(1)
                   ens(1) = ens(1) + 0.5_pr*W(i1,i2,i3,1)**2*dV
                   ens(2) = ens(2) + 0.5_pr*W(i1,i2,i3,2)**2*dV
                   ens(3) = ens(3) + 0.5_pr*W(i1,i2,i3,3)**2*dV
                END DO
             END DO
          END DO

          DEALLOCATE(W)

        END FUNCTION Enstrophy

        !=========================================================
        ! Calculate dEdt from velocity field in Physical space
        !=========================================================
        FUNCTION calc_dEdt(myfield) RESULT (R)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: myfield
          REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: aux1, aux2, aux3          

          REAL(pr), DIMENSION(1:3) :: local_VecField_norm, global_vecField_norm
          REAL(pr), DIMENSION(1:2) :: R

          R = 0.0_pr 
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

          R(1) = -2.0_pr*visc*SUM(global_VecField_norm)

          local_VecField_norm = field_inner_product(aux3,aux1,"L2")           
          CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 

          R(2) = SUM(global_VecField_norm)

          DEALLOCATE(aux1)
          DEALLOCATE(aux2)
          DEALLOCATE(aux3)
           
        END FUNCTION calc_dEdt


        !====================================================
        ! CALCULATE HELICITY FROM U AND W
        ! H = /int( U /cdot W)
        !====================================================
        FUNCTION Helicity(U,W) RESULT (H)
          USE global_variables
          IMPLICIT NONE
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: U,W    
          REAL(pr) :: H

          INTEGER :: ii, jj, kk

          H = 0.0_pr

          DO kk=1,local_nlast
             DO jj=1,n(2)
                DO ii=1,n(1)
                   H = H + ( U(ii,jj,kk,1)*W(ii,jj,kk,1) + U(ii,jj,kk,2)*W(ii,jj,kk,2) + U(ii,jj,kk,3)*W(ii,jj,kk,3) )*dV
                END DO
             END DO
          END DO

        END FUNCTION Helicity 

        !============================================================
        ! Diagnostics: |U|, |W|, stretching factor
        !============================================================
        SUBROUTINE diagnosticScalars(U, W, myIter)
          USE global_variables
          USE data_ops
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: U
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: W
          INTEGER, INTENT(IN) :: myIter
          
          REAL(pr), DIMENSION(:,:), ALLOCATABLE :: Spectrum
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux
          REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: auxVec, allDiagFields
          REAL(pr), DIMENSION(1:3) :: local_q3, K, E, Umax, Wmax, vorCoreData
          REAL(pr) :: local_q1, dEdt, H, magUmax, magWmax, maxHel, minHel, minStretch, maxStretch, divU_L2norm

          REAL(pr), DIMENSION(1:3) :: myPoint, myNormal

          INTEGER :: nn, i1, i2, i3 
          INTEGER, PARAMETER :: numDiagFields = 4 

          ALLOCATE( Spectrum(1:n(1)/2,1:2) )
          CALL calculate_spectral_data(U, Spectrum)

          !ALLOCATE( auxVec(1:n(1),1:n(2),1:local_nlast,1:3) )
          ALLOCATE( allDiagFields(1:n(1),1:n(2),1:local_nlast,1:numDiagFields) )
 
          !local_q3 = (/ 1.0_pr, -1.0_pr, 0.0_pr /) 
          !myPoint = FindMinPoint(magU,magW, local_q3 )
          !IF (myIter==0) THEN           
          !   myPoint = (/ REAL(120,pr)/REAL(n(1),pr), REAL(116,pr)/REAL(n(2),pr), REAL(120,pr)/REAL(n(3),pr) /)
          !   CALL shift_space(U, myPoint)
          !   CALL shift_space(W, myPoint)

             !myPoint = (/ PI/2.0_pr, PI/2.0_pr, PI/2.0_pr /)
             !CALL rotate_space(U, myPoint)
             !CALL vel2vort(U,W) 
          !END IF

          !CALL advection(W, U, auxVec)
          allDiagFields(:,:,:,1) = SQRT( U(:,:,:,1)**2 + U(:,:,:,2)**2 + U(:,:,:,3)**2 )
          allDiagFields(:,:,:,2) = SQRT( W(:,:,:,1)**2 + W(:,:,:,2)**2 + W(:,:,:,3)**2 )
          !allDiagFields(:,:,:,3) = ( W(:,:,:,1)*auxVec(:,:,:,1) + W(:,:,:,2)*auxVec(:,:,:,2) + W(:,:,:,3)*auxVec(:,:,:,3) )/( W(:,:,:,1)**2 + W(:,:,:,2)**2 + W(:,:,:,3)**2 ) 
          allDiagFields(:,:,:,3) = U(:,:,:,1)*W(:,:,:,1) + U(:,:,:,2)*W(:,:,:,2) + U(:,:,:,3)*W(:,:,:,3)
          !DEALLOCATE( auxVec )

          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) ) 
          CALL divergence(U, aux)          
          local_q1 = inner_product(aux, aux, "L2")
          CALL MPI_ALLREDUCE(local_q1, divU_L2norm, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

          local_q3 = Energy(U)
          CALL MPI_ALLREDUCE(local_q3, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
          local_q3 = Energy(W)
          CALL MPI_ALLREDUCE(local_q3, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

          local_q3(1) = MAXVAL(U(:,:,:,1))
          local_q3(2) = MAXVAL(U(:,:,:,2))
          local_q3(3) = MAXVAL(U(:,:,:,3))
          CALL MPI_ALLREDUCE(local_q3, Umax, 3, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
         
          local_q3(1) = MAXVAL(W(:,:,:,1))        
          local_q3(2) = MAXVAL(W(:,:,:,2))
          local_q3(3) = MAXVAL(W(:,:,:,3))
          CALL MPI_ALLREDUCE(local_q3, Wmax, 3, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
 
          local_q1 = Helicity(U, W)
          CALL MPI_ALLREDUCE(local_q1, H, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

          aux = allDiagFields(:,:,:,1)
          local_q1 = MAXVAL(aux)
          CALL MPI_ALLREDUCE(local_q1, magUmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
          aux = allDiagFields(:,:,:,2)
          local_q1 = MAXVAL(aux)
          CALL MPI_ALLREDUCE(local_q1, magWmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
          aux = allDiagFields(:,:,:,3)
          local_q1 = MAXVAL(aux)
          CALL MPI_ALLREDUCE(local_q1, maxHel, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
          local_q1 = MINVAL(aux)
          CALL MPI_ALLREDUCE(local_q1, minHel, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, Statinfo)

          CALL calculate_geometric_data(U, W, aux, "Q", vorCoreData)
          allDiagFields(:,:,:,4) = aux

          CALL calculate_ring_data(allDiagFields(:,:,:,2))

          CALL save_diagnosticScalars(allDiagFields, numDiagFields, "magU,magW,Helicity,VortexCore", "maxdEdt")
          IF (rank==0) THEN 
             CALL save_spectral_data(Spectrum)
             CALL save_diagnosticFields_global("maxdEdt", K, E, Umax, Wmax, magUmax, magWmax, H, maxHel, minHel, vorCoreData)
          END IF
          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

          DEALLOCATE( Spectrum )
          DEALLOCATE( aux )
          DEALLOCATE( allDiagFields )

        END SUBROUTINE diagnosticScalars 

!        !==========================================================
!        ! Calculate diagnostics for the current fields.
!        ! Energy, Enstrophy, Umax, Wmax, Helicity
!        !==========================================================
!        SUBROUTINE diagnosticFields_global(U, W, mysystem)
!          USE global_variables
!          USE data_ops
!          IMPLICIT NONE
!          INCLUDE "mpif.h" 
!
!
!          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: U, W
!          CHARACTER(len=*), INTENT(IN) :: mysystem
!
!          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: scalarField
!          REAL(pr), DIMENSION(1:3) :: local_q3, K, E, Umax, Wmax
!          REAL(pr) :: local_q1, H, magUmax, magWmax, maxHel, minHel, maxStretch, minStretch
!          INTEGER :: ii
!
!          local_q3 = Energy(U)
!          CALL MPI_ALLREDUCE(local_q3, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
!          local_q3 = Energy(W)
!          CALL MPI_ALLREDUCE(local_q3, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
!
!          local_q3(1) = MAXVAL(U(:,:,:,1))
!          local_q3(2) = MAXVAL(U(:,:,:,2))
!          local_q3(3) = MAXVAL(U(:,:,:,3))
!          CALL MPI_ALLREDUCE(local_q3, Umax, 3, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
!         
!          local_q3(1) = MAXVAL(W(:,:,:,1))        
!          local_q3(2) = MAXVAL(W(:,:,:,2))
!          local_q3(3) = MAXVAL(W(:,:,:,3))
!          CALL MPI_ALLREDUCE(local_q3, Wmax, 3, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
! 
!          local_q1 = Helicity(U, W)
!          CALL MPI_ALLREDUCE(local_q1, H, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 
!
!          ALLOCATE( scalarField(1:n(1),1:n(2),1:local_nlast) )
!
!          scalarField = 0.0_pr
!          DO ii=1,3
!             scalarField = scalarField + U(:,:,:,ii)**2
!          END DO
!          scalarField = SQRT(scalarField)
!          local_q1 = MAXVAL(scalarField)
!          CALL MPI_ALLREDUCE(local_q1, magUmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
!
!          scalarField = 0.0_pr
!          DO ii=1,3
!             scalarField = scalarField + W(:,:,:,ii)**2
!          END DO
!          scalarField = SQRT(scalarField)
!          local_q1 = MAXVAL(scalarField)
!          CALL MPI_ALLREDUCE(local_q1, magWmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
!
!          IF (rank==0) THEN 
!             CALL save_diagnosticFields_global(mysystem, K, E, Umax, Wmax, magUmax, magWmax, H, maxHel, minHel, maxStretch, minStretch)
!          END IF
!          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
!
!          DEALLOCATE( scalarField )
!
!
!        END SUBROUTINE diagnosticFields_global
!


!        !============================================================
!        ! Diagnostics of optimization
!        !============================================================
!        SUBROUTINE diagnostics_optim(iteration, tau, J, u, w, gradJ)
!          USE global_variables
!          USE data_ops
!          USE FFT2
!          IMPLICIT NONE
!          INCLUDE "mpif.h"
!
!          REAL(pr), DIMENSION(1:n(1),1:local_nlast,1:2), INTENT(IN) :: u, gradJ
!          REAL(pr), DIMENSION(1:n(1),1:local_nlast), INTENT(IN) :: w
!          REAL(pr), INTENT(IN) :: tau, J
!          INTEGER, INTENT(IN) :: iteration 
!
!          COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: fu, fgradJ
!          COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: fw, aux, faux
!          
!          REAL(pr) :: global_ens, global_wmax, L2div, L2div_gradJ
!          REAL(pr) :: local_ens, local_wmax, local_L2div
!          REAL(pr), DIMENSION(1:2) :: global_ener, global_umax, global_palins, avg_vel, global_ener_gradJ
!          REAL(pr), DIMENSION(1:2) :: local_ener, local_umax, local_palins
!          REAL(pr), DIMENSION(1:n(1)/2, 3) :: spectral_data
!          REAL(pr), DIMENSION(1:2) :: local_spectral_data, global_spectral_data
!          INTEGER :: i, i1, i2, mode_count
!          REAL(pr) :: kk_min, kk_max, norm_K
!
!          DO i=1,2
!             aux = CMPLX(u(:,:,i),0.0_pr)
!             CALL ffourier(aux, faux)
!             fu(:,:,i) = faux
!             
!             aux = CMPLX(gradJ(:,:,i),0.0_pr)
!             CALL ffourier(aux, faux)
!             fgradJ(:,:,i) = faux
!          END DO
!          aux = CMPLX(w,0.0_pr)
!          CALL ffourier(aux, fw)
! 
!          DO i=0,n(1)/2-1
!             kk_min = 2.0_pr*PI*REAL(i,pr)
!             kk_max = 2.0_pr*PI*REAL(i+1,pr)
!             spectral_data(i+1,1) = kk_max
!             local_spectral_data(1:2) = 0.0_pr
!             global_spectral_data(1:2) = 0.0_pr
!             mode_count = 0
!             DO i2=1,local_nlast
!                DO i1=1,n(1)
!                   norm_K = SQRT( K1(i1)**2 + K2(i2+local_last_start)**2 )
!                   IF (  kk_min < norm_K .AND. norm_K <= kk_max ) THEN
!                      mode_count = mode_count+1
!                      local_spectral_data(1) = local_spectral_data(1) + ABS(fu(i1,i2,1))**2*dV + ABS(fu(i1,i2,2))**2*dV
!                      local_spectral_data(2) = local_spectral_data(2) + ABS(fw(i1,i2))**2*dV
!                   END IF
!                END DO
!             END DO
!             IF ( mode_count /= 0) THEN
!                local_spectral_data(1) = local_spectral_data(1)/REAL(mode_count,pr)
!                local_spectral_data(2) = local_spectral_data(2)/REAL(mode_count,pr)
!             END IF
!
!             CALL MPI_ALLREDUCE(local_spectral_data, global_spectral_data, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
!             spectral_data(i+1,2) = global_spectral_data(1)
!             spectral_data(i+1,3) = global_spectral_data(2) 
!          END DO  
!          
!          IF (rank==0) THEN 
!             CALL save_spectral_data(spectral_data, iteration, "optim")
!          END IF
!          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)           
!          CALL get_max_mode("optim", iteration)
!
!          ! CALCULATE ENERGY OF PHI AND GRADJ
!          local_ener = kinetic_energy(fu)
!          CALL MPI_ALLREDUCE(local_ener, global_ener, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
!          local_ener = kinetic_energy(fgradJ)
!          CALL MPI_ALLREDUCE(local_ener, global_ener_gradJ, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
!          
!          ! CALCULATE DIVERGENCE OF PHI AND GRADJ
!          DO i2=1,local_nlast 
!             DO i1=1,n(1)
!                faux(i1,i2) = CMPLX(0.0_pr, K1(i1))*fu(i1,i2,1) + CMPLX(0.0_pr, K2(i2+local_last_start))*fu(i1,i2,2)
!             END DO
!          END DO
!          local_L2div = 2.0_pr*enstrophy(faux)
!          CALL MPI_ALLREDUCE(local_L2div, L2div, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
!          
!          DO i2=1,local_nlast 
!             DO i1=1,n(1)
!                faux(i1,i2) = CMPLX(0.0_pr, K1(i1))*fgradJ(i1,i2,1) + CMPLX(0.0_pr, K2(i2+local_last_start))*fgradJ(i1,i2,2)
!             END DO
!          END DO
!          local_L2div = 2.0_pr*enstrophy(faux)
!          CALL MPI_ALLREDUCE(local_L2div, L2div_gradJ, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
!          
!          !CALCULATE PALINSTROPHY
!          local_palins = palinstrophy(fu)
!          CALL MPI_ALLREDUCE(local_palins, global_palins, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
!        
!          !local_umax(1) = MAXVAL(ABS(u(:,:,1)))
!          !local_umax(2) = MAXVAL(ABS(u(:,:,2)))
!          !local_wmax = MAXVAL(ABS(w))
!          !CALL MPI_ALLREDUCE(local_umax, global_umax, 2, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
!          !CALL MPI_ALLREDUCE(local_wmax, global_wmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo) 
!         
!          IF (rank==0) THEN 
!             avg_vel(1) = fu(1,1,1)
!             avg_vel(2) = fu(1,1,2)
!             CALL save_diagnostics_optim(iteration, tau, J, global_ener, global_ener_gradJ, global_palins, L2div, L2div_gradJ, avg_vel)
!          END IF
!          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo) 
!        
!        END SUBROUTINE diagnostics_optim 

        !====================================================
        ! Calculate the spectrum of the velocity field
        ! using spherical shells
        !====================================================
        SUBROUTINE calculate_spectral_data(u, spectral_data)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: u
          REAL(pr), DIMENSION(1:n(1)/2, 1:2), INTENT(OUT) :: spectral_data

          COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: fu
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux

          INTEGER :: i, i1, i2, i3, mode_count
          REAL(pr) :: kk_min, kk_max, norm_K

          REAL(pr) :: local_spectral_data, global_spectral_data, spectral_Ener
          REAL(pr), DIMENSION(1:3) :: local_q3, Ener
 

          ALLOCATE( fu(1:n(1),1:n(2),1:local_nlast,1:3) )
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )

          local_q3 = Energy(u)
          CALL MPI_ALLREDUCE(local_q3, Ener, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
 
          DO i=1,3
             aux = CMPLX(u(:,:,:,i),0.0_pr)
             CALL ffourier(aux, faux)
             fu(:,:,:,i) = faux
          END DO
 
          DO i=0,n(1)/2-1
             kk_min = 2.0_pr*PI*REAL(i,pr)
             kk_max = 2.0_pr*PI*REAL(i+1,pr)
             spectral_data(i+1,1) = kk_max
             local_spectral_data = 0.0_pr
             global_spectral_data = 0.0_pr
             mode_count = 0
             DO i3=1,local_nlast
                DO i2=1,n(1)
                   DO i1=1,n(1)
                      norm_K = SQRT( K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2 )
                      IF (  kk_min < norm_K .AND. norm_K <= kk_max ) THEN
                         mode_count = mode_count+1
                         local_spectral_data = local_spectral_data + ABS(fu(i1,i2,i3,1))**2 + ABS(fu(i1,i2,i3,2))**2 + ABS(fu(i1,i2,i3,3))**2
                      END IF
                   END DO
                END DO
             END DO
             !IF ( mode_count /= 0) THEN
             !   local_spectral_data(1) = local_spectral_data(1)/REAL(mode_count,pr)
             !   local_spectral_data(2) = local_spectral_data(2)/REAL(mode_count,pr)
             !END IF

             CALL MPI_ALLREDUCE(local_spectral_data, global_spectral_data, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
             spectral_data(i+1,2) = global_spectral_data 
          END DO  

          spectral_Ener = 2.0_pr*PI*SUM(spectral_data(:,2))
          spectral_data(:,2) = SUM(Ener)/spectral_Ener*spectral_data(:,2)

          DEALLOCATE( fu )
          DEALLOCATE( aux )
          DEALLOCATE( faux )
 
        END SUBROUTINE calculate_spectral_data

        !====================================================================
        ! OBTAIN GEOMETRIC INFORMATION FROM VORTICITY FIELD
        !====================================================================
        SUBROUTINE calculate_geometric_data(U, W, vortexCore, vortexCriterion, diagScalars)
          USE global_variables
          USE data_ops
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: U, W
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast), INTENT(OUT) :: vortexCore
          CHARACTER(len=*), INTENT(IN) :: vortexCriterion
          REAL(pr), DIMENSION(1:3), INTENT(OUT) :: diagScalars          
 
          REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: grad_Ux, grad_Uy, grad_Uz
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux

          INTEGER, PARAMETER :: lwmax = 1000
 
          REAL(pr), DIMENSION(1:3,1:3) :: A, B
          REAL(pr), DIMENSION(1:3,1:3) :: dP, vR, vL
          REAL(pr), DIMENSION(1:3) :: wr, wi, lambda
          REAL(pr), DIMENSION(1:lwmax) :: work
          INTEGER :: info, lwork, ii, jj, kk, nn, mm

          REAL(pr) :: local_real, global_real
        
          ALLOCATE( grad_Ux(1:n(1),1:n(2),1:local_nlast,1:3) )  
          ALLOCATE( grad_Uy(1:n(1),1:n(2),1:local_nlast,1:3) )  
          ALLOCATE( grad_Uz(1:n(1),1:n(2),1:local_nlast,1:3) )  
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )

          aux = U(:,:,:,1)
          CALL gradient(aux, grad_Ux)
          aux = U(:,:,:,2)
          CALL gradient(aux, grad_Uy) 
          aux = U(:,:,:,3)
          CALL gradient(aux, grad_Uz)
          DEALLOCATE( aux ) 

          local_real = 0.0_pr
          global_real = 0.0_pr

          SELECT CASE (vortexCriterion)
             CASE ("EigVals")
                DO kk=1,local_nlast
                   DO jj=1,n(2)
                      DO ii=1,n(1)
                         A(1,1) = grad_Ux(ii,jj,kk,1)
                         A(1,2) = grad_Ux(ii,jj,kk,2)
                         A(1,3) = grad_Ux(ii,jj,kk,3)
                         A(2,1) = grad_Uy(ii,jj,kk,1)
                         A(2,2) = grad_Uy(ii,jj,kk,2)
                         A(2,3) = grad_Uy(ii,jj,kk,3)
                         A(3,1) = grad_Uz(ii,jj,kk,1)
                         A(3,2) = grad_Uz(ii,jj,kk,2)
                         A(3,3) = grad_Uz(ii,jj,kk,3)

                         B(1,1) = grad_Ux(ii,jj,kk,1)
                         B(1,2) = grad_Uy(ii,jj,kk,1)
                         B(1,3) = grad_Uz(ii,jj,kk,1)
                         B(2,1) = grad_Ux(ii,jj,kk,2)
                         B(2,2) = grad_Uy(ii,jj,kk,2)
                         B(2,3) = grad_Uz(ii,jj,kk,2)
                         B(3,1) = grad_Ux(ii,jj,kk,3)
                         B(3,2) = grad_Uy(ii,jj,kk,3)
                         B(3,3) = grad_Uz(ii,jj,kk,3)
                   
                         dP = 0.5_pr*( MATMUL(A,A) + MATMUL(B,B) )
                         CALL dgeev('N','V', 3, dP, 3, wr, wi, vL, 3, vR, 3, work, lwmax, info) 
              
                         DO nn=1,3
                            lambda(nn) = MAXVAL(wr)
                            mm = MAXLOC(wr,1)
                            wr(mm) = MINVAL(wr)-1.0_pr
                         END DO
                                             
                         vortexCore(ii,jj,kk) = SIGN(MIN(lambda(2), 0.0_pr),1.0_pr)

                         IF ( vortexCore(ii,jj,kk) > 0.0_pr ) THEN
                            local_real = local_real + dV
                         END IF 
 
                      END DO
                   END DO
                END DO 


             CASE ("Q")  
                DO kk=1,local_nlast
                   DO jj=1,n(2)
                      DO ii=1,n(1)
                         vortexCore(ii,jj,kk) = MAX(0.25_pr*( W(ii,jj,kk,1)**2 + W(ii,jj,kk,2)**2 + W(ii,jj,kk,3)**2 ) - &
                                                    0.25_pr*( grad_Ux(ii,jj,kk,2)+grad_Uy(ii,jj,kk,1) )**2 - &
                                                    0.25_pr*( grad_Ux(ii,jj,kk,3)+grad_Uz(ii,jj,kk,1) )**2 - &
                                                    0.25_pr*( grad_Uy(ii,jj,kk,3)+grad_Uz(ii,jj,kk,2) )**2 - &
                                                    0.50_pr*( grad_Ux(ii,jj,kk,1)**2 + grad_Uy(ii,jj,kk,2)**2 + grad_Uz(ii,jj,kk,3)**2 ), 0.0_pr)                         

                         IF ( vortexCore(ii,jj,kk) > 0.0_pr ) THEN
                            local_real = local_real + dV
                         END IF 

                      END DO
                   END DO
                END DO 

          END SELECT

          CALL MPI_ALLREDUCE(local_real, global_real, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 
          diagScalars(1) = global_real

          local_real = MAXVAL(vortexCore)
          CALL MPI_ALLREDUCE(local_real, global_real, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
          diagScalars(2) = global_real

          vortexCore = vortexCore/diagScalars(2)

          DEALLOCATE( grad_Ux )
          DEALLOCATE( grad_Uy )
          DEALLOCATE( grad_Uz )

        END SUBROUTINE calculate_geometric_data 

        !==========================================================
        ! Calculate location of the ring-like vorticity structure
        !==========================================================
        SUBROUTINE calculate_ring_data(magW)
          USE global_variables
          USE data_ops
          IMPLICIT NONE
          INCLUDE "mpif.h"
 
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast), INTENT(IN) :: magW

          INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ringMask
          REAL(pr), DIMENSION(:,:), ALLOCATABLE :: local_ringLoc, global_ringLoc
          REAL(pr), DIMENSION(1:3) :: dx 
          REAL(pr) :: local_maxWmag, global_maxWmag
          INTEGER :: mm, nn, ii, jj, kk, myflag
          INTEGER :: local_numPointsRing, global_numPointsRing

          local_maxWmag = MAXVAL(magW)
          CALL MPI_ALLREDUCE(local_maxWmag, global_maxWmag, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)

          dx = 1.0_pr/REAL(n,pr)

          ALLOCATE( ringMask(1:n(1),1:n(2),1:local_nlast) )

          DO kk=1,local_nlast
             DO jj=1,n(2)
                DO ii=1,n(1)
                   IF ( magW(ii,jj,kk)/global_maxWmag > 0.95_pr ) THEN
                      ringMask(ii,jj,kk) = 1
                   ELSE
                      ringMask(ii,jj,kk) = 0
                   END IF
                END DO
             END DO
          END DO
                      
          local_numPointsRing = SUM(ringMask)

          ALLOCATE( local_ringLoc(1:local_numPointsRing,3) )
          nn = 1

          DO kk=1,local_nlast
             DO jj=1,n(2)
                DO ii=1,n(1)
                   IF ( ringMask(ii,jj,kk) == 1 ) THEN
                      local_ringLoc(nn,1) = REAL(ii-1,pr)*dx(1)
                      local_ringLoc(nn,2) = REAL(jj-1,pr)*dx(2)
                      local_ringLoc(nn,3) = REAL(local_last_start+kk-1,pr)*dx(3)
               
                      nn = nn+1

                   END IF
                END DO
             END DO
          END DO
          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
          !CALL MPI_ALLREDUCE(local_numPointsRing, global_numPointsRing, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, Statinfo)

          myflag = 0
          DO mm=0,np-1 
             IF (rank==mm) THEN
                IF (local_numPointsRing > 10) THEN
                   CALL save_ring_data(local_ringLoc, local_numPointsRing, myflag)
                   myflag = 1
                END IF
             END IF
             CALL MPI_BCAST(myflag, 1, MPI_INTEGER, mm, MPI_COMM_WORLD, Statinfo) 
          END DO

          DEALLOCATE( ringMask )
          DEALLOCATE( local_ringLoc ) 
   
        END SUBROUTINE calculate_ring_data

        !=====================================================================
        ! OBTAINS THE DIV_FREE PORJECTION OF A GIVEN VECTOR FIELD
        !=====================================================================
        SUBROUTINE div_free(myfield)   
          USE global_variables
          USE FFT2
          IMPLICIT NONE
 
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: myfield

          REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: grad_phi
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: divU 
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux, divU_hat
          COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: grad_phi_hat, U_hat

          INTEGER  :: i1, i2, i3, ii
          REAL(pr) :: ksq
          REAL(pr), DIMENSION (1:3) :: k
          COMPLEX(pr), DIMENSION(1:3) :: tmp

          ALLOCATE( grad_phi(1:n(1),1:n(2),1:local_nlast,1:3) )
          ALLOCATE( divU(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( divU_hat(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( grad_phi_hat(1:n(1),1:n(2),1:local_nlast,1:3) )

          CALL divergence(myfield, divU)          
          aux = CMPLX(divU,0.0_pr)
          CALL ffourier(aux,faux)
          divU_hat = faux 

          DO i3 = 1, local_nlast
             DO i2 = 1, n(2)
                DO i1 = 1, n(1)        
                   ksq = K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2
                   IF (ksq > MACH_EPSILON) THEN
                      grad_phi_hat(i1,i2,i3,1) = -CMPLX(0.0_pr,K1(i1))*divU_hat(i1,i2,i3)/ksq
                      grad_phi_hat(i1,i2,i3,2) = -CMPLX(0.0_pr,K2(i2))*divU_hat(i1,i2,i3)/ksq
                      grad_phi_hat(i1,i2,i3,3) = -CMPLX(0.0_pr,K3(i3+local_last_start))*divU_hat(i1,i2,i3)/ksq 
                   ELSE
                      grad_phi_hat(i1,i2,i3,1) = 0.0_pr
                      grad_phi_hat(i1,i2,i3,2) = 0.0_pr
                      grad_phi_hat(i1,i2,i3,3) = 0.0_pr 
                   END IF
                END DO 
             END DO
          END DO

          DO ii=1,3
             faux = grad_phi_hat(:,:,:,ii)
             CALL bfourier(faux,aux)
             grad_phi(:,:,:,ii) = REAL(aux)
          END DO

          myfield = myfield - grad_phi

          DEALLOCATE( grad_phi )
          DEALLOCATE( divU )
          DEALLOCATE( aux )
          DEALLOCATE( faux )
          DEALLOCATE( divU_hat )
          DEALLOCATE( grad_phi_hat )


        END SUBROUTINE div_free

!        !===========================================================
!        ! Calculate bilinear form J(f,g) = det(f_x f_y; g_x g_y)
!        !===========================================================
!
!        SUBROUTINE J_bilinear_form(f,g,J)
!          USE global_variables
!          USE FFT2
!
!          REAL(pr), DIMENSION(1:n(1),1:local_nlast), INTENT(IN) :: f,g
!          REAL(pr), DIMENSION(1:n(1),1:local_nlast), INTENT(OUT) :: J   
!          
!          COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: aux, faux0, faux1
!          REAL(pr), DIMENSION(1:n(1), 1:local_nlast) :: fx, fy, gx, gy
!          INTEGER :: i1, i2
!
!          aux = CMPLX(f,0.0_pr)
!          CALL ffourier(aux, faux0)
!          
!          DO i2=1,local_nlast
!             DO i1=1,n(1)
!                faux1(i1,i2) = CMPLX(0.0_pr,K1(i1))*faux0(i1,i2)
!             END DO
!          END DO 
!          CALL bfourier(faux1, aux)
!          fx = REAL(aux)
!
!          DO i2=1,local_nlast
!             DO i1=1,n(1)
!                faux1(i1,i2) = CMPLX(0.0_pr,K2(i2+local_last_start))*faux0(i1,i2)
!             END DO
!          END DO 
!          CALL bfourier(faux1, aux)
!          fy = REAL(aux)
!
!          aux = CMPLX(g,0.0_pr)
!          CALL ffourier(aux, faux0)
!          
!          DO i2=1,local_nlast
!             DO i1=1,n(1)
!                faux1(i1,i2) = CMPLX(0.0_pr,K1(i1))*faux0(i1,i2)
!             END DO
!          END DO 
!          CALL bfourier(faux1, aux)
!          gx = REAL(aux)
!
!          DO i2=1,local_nlast
!             DO i1=1,n(1)
!                faux1(i1,i2) = CMPLX(0.0_pr,K2(i2+local_last_start))*faux0(i1,i2)
!             END DO
!          END DO 
!          CALL bfourier(faux1, aux)
!          gy = REAL(aux)
!
!          J = fx*gy - fy*gx
!
!        END SUBROUTINE J_bilinear_form
!
!
        !======================================================
        ! FIX THE ENERGY OF A GIVEN VELOCITY FIELD
        !======================================================
        SUBROUTINE Fix_K0(myfield)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1), 1:n(2), 1:local_nlast,1:3), INTENT(INOUT) :: myfield

          COMPLEX(pr), DIMENSION(1:n(1), 1:n(2), 1:local_nlast) :: aux, faux

          REAL(pr), DIMENSION(1:3) :: local_E, global_E
          REAL(pr), DIMENSION(1:3) :: local_K, global_K
          REAL(pr) :: alpha
          INTEGER :: ii
         
          local_K = Energy(myfield)
          CALL MPI_ALLREDUCE(local_K, global_K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
          alpha = SQRT(K0/SUM(global_K))

          myfield = alpha*myfield          

        END SUBROUTINE Fix_K0

        !======================================================
        ! FIX THE ENSTROPHY OF A GIVEN VELOCITY FIELD
        !======================================================
        SUBROUTINE Fix_E0(myfield)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1), 1:n(2), 1:local_nlast,1:3), INTENT(INOUT) :: myfield

          COMPLEX(pr), DIMENSION(1:n(1), 1:n(2), 1:local_nlast) :: aux, faux

          REAL(pr), DIMENSION(1:3) :: local_E, global_E
          REAL(pr), DIMENSION(1:3) :: local_K, global_K
          REAL(pr) :: alpha
          INTEGER :: ii
         
          local_E = Enstrophy(myfield)
          CALL MPI_ALLREDUCE(local_E, global_E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
          alpha = SQRT(E0/SUM(global_E))

          myfield = alpha*myfield          

        END SUBROUTINE Fix_E0

        !==============================================================
        ! Calculates the vorticity w in physical space given
        ! velocity u in physical space
        !==============================================================
        SUBROUTINE vel2vort(vel,vort)
          USE global_variables
          USE FFT2
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:n(1), 1:n(2), 1:local_nlast, 1:3), INTENT(IN) :: vel
          REAL(pr), DIMENSION(1:n(1), 1:n(2), 1:local_nlast, 1:3), INTENT(OUT) :: vort

          INTEGER  :: i1, i2, i3

          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: Ux_hat, Uy_hat, Uz_hat, Wx_hat, Wy_hat, Wz_hat
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux

          ALLOCATE( Ux_hat(1:n(1),1:n(2),1:local_nlast) ) 
          ALLOCATE( Uy_hat(1:n(1),1:n(2),1:local_nlast) ) 
          ALLOCATE( Uz_hat(1:n(1),1:n(2),1:local_nlast) ) 
          ALLOCATE( Wx_hat(1:n(1),1:n(2),1:local_nlast) ) 
          ALLOCATE( Wy_hat(1:n(1),1:n(2),1:local_nlast) ) 
          ALLOCATE( Wz_hat(1:n(1),1:n(2),1:local_nlast) ) 
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) ) 
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )

          aux = CMPLX(vel(:,:,:,1), 0.0_pr)
          CALL ffourier(aux, faux)
          Ux_hat = faux

          aux = CMPLX(vel(:,:,:,2), 0.0_pr)
          CALL ffourier(aux, faux)
          Uy_hat = faux

          aux = CMPLX(vel(:,:,:,3), 0.0_pr)
          CALL ffourier(aux, faux)
          Uz_hat = faux

         
          DO i3 = 1,local_nlast
             DO i2 = 1,n(2)
                DO i1 = 1,n(1)
                   Wx_hat(i1,i2,i3) = CMPLX(0.0_pr, K2(i2))*Uz_hat(i1,i2,i3) - CMPLX(0.0_pr, K3(i3+local_last_start))*Uy_hat(i1,i2,i3)
                   Wy_hat(i1,i2,i3) = CMPLX(0.0_pr, K3(i3+local_last_start))*Ux_hat(i1,i2,i3) - CMPLX(0.0_pr, K1(i1))*Uz_hat(i1,i2,i3)
                   Wz_hat(i1,i2,i3) = CMPLX(0.0_pr, K1(i1))*Uy_hat(i1,i2,i3) - CMPLX(0.0_pr, K2(i2))*Ux_hat(i1,i2,i3)
                END DO
             END DO
          END DO

          !--Transform back
          CALL bfourier (Wx_hat, aux)
          vort(:,:,:,1) = REAL(aux)
          CALL bfourier (Wy_hat, aux)
          vort(:,:,:,2) = REAL(aux)
          CALL bfourier (Wz_hat, aux)
          vort(:,:,:,3) = REAL(aux)
 
          DEALLOCATE( Ux_hat )
          DEALLOCATE( Uy_hat )
          DEALLOCATE( Uz_hat )
          DEALLOCATE( Wx_hat )
          DEALLOCATE( Wy_hat )
          DEALLOCATE( Wz_hat )
          DEALLOCATE( aux )
          DEALLOCATE( faux )

        END SUBROUTINE vel2vort

        !=====================================
        ! Transform vorticity to velocity
        !=====================================
        SUBROUTINE vort2vel(W, U)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: W
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(OUT) :: U
          
          COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: fU, fW
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
          INTEGER :: i1, i2, i3, ii
          REAL(pr) :: ksq

          ALLOCATE( fU(1:n(1),1:n(2),1:local_nlast,1:3) )
          ALLOCATE( fW(1:n(1),1:n(2),1:local_nlast,1:3) )
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )
          
          DO ii=1,3
             aux = CMPLX(W(:,:,:,ii),0.0_pr) 
             CALL ffourier(aux,faux)
             fW(:,:,:,ii) = faux
          END DO

          DO i3 = 1,local_nlast
             DO i2 = 1,n(2)
                DO i1 = 1,n(1)
                   ksq = K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2
                   IF (ksq > MACH_EPSILON) THEN
                      fU(i1,i2,i3,1) = ( CMPLX(0.0_pr,K2(i2))*fW(i1,i2,i3,3) - CMPLX(0.0_pr,K3(i3+local_last_start))*fW(i1,i2,i3,2) )/ksq
                      fU(i1,i2,i3,2) = ( CMPLX(0.0_pr,K3(i3+local_last_start))*fW(i1,i2,i3,1) - CMPLX(0.0_pr,K1(i1))*fW(i1,i2,i3,3) )/ksq
                      fU(i1,i2,i3,3) = ( CMPLX(0.0_pr,K1(i1))*fW(i1,i2,i3,2) - CMPLX(0.0_pr,K2(i2))*fW(i1,i2,i3,1) )/ksq
                   ELSE
                      fU(i1,i2,i3,1) = CMPLX(0.0_pr,0.0_pr)
                      fU(i1,i2,i3,2) = CMPLX(0.0_pr,0.0_pr)
                      fU(i1,i2,i3,3) = CMPLX(0.0_pr,0.0_pr)      
                   END IF
                END DO
             END DO
          END DO

          DO ii=1,3
             faux = fU(:,:,:,ii)
             CALL bfourier(faux,aux)
             U(:,:,:,ii) = REAL(aux)
          END DO    
     
          DEALLOCATE( fU )
          DEALLOCATE( fW )
          DEALLOCATE( aux )
          DEALLOCATE( faux )

        END SUBROUTINE vort2vel

        !========================================
        ! CALCULATE DERIVATIVE WRT ii-th VARIABLE
        !========================================
        SUBROUTINE derivative(u, ii)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast), INTENT(INOUT) :: u
          INTEGER, INTENT(IN) :: ii

          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: fu
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
          INTEGER :: i1, i2, i3
          REAL(pr) :: mode           
 
          ALLOCATE( fu(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )

          aux = CMPLX(u,0.0_pr)
          CALL ffourier(aux, fu)

          DO i3=1, local_nlast
             DO i2=1,n(2)
                DO i1=1, n(1)
                   IF (ii .EQ. 1) THEN
                      faux(i1,i2,i3) = CMPLX(0.0_pr, K1(i1))*fu(i1,i2,i3)
                   ELSEIF (ii .EQ. 2) THEN
                      faux(i1,i2,i3) = CMPLX(0.0_pr, K2(i2))*fu(i1,i2,i3)
                   ELSEIF (ii .EQ. 3) THEN
                      faux(i1,i2,i3) = CMPLX(0.0_pr, K3(i3+local_last_start))*fu(i1,i2,i3)
                   END IF
                END DO
             END DO
          END DO
          
          CALL bfourier(faux, aux)
          u = REAL(aux)

          DEALLOCATE(fu)
          DEALLOCATE(aux)
          DEALLOCATE(faux)
 
        END SUBROUTINE derivative

        !====================================
        ! CALCULATE GRADIENT OF FUNCTION
        !====================================
        SUBROUTINE gradient(u, grad_u)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast), INTENT(IN) :: u
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(OUT) :: grad_u

          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: fu
          COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: fgrad_u
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
          INTEGER :: i1, i2, i3, ii
          REAL(pr) :: mode           
 
          ALLOCATE( fu(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( fgrad_u(1:n(1),1:n(2),1:local_nlast,1:3) )
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )

          aux = CMPLX(u,0.0_pr)
          CALL ffourier(aux, fu)

          DO i3 = 1, local_nlast
             DO i2 = 1,n(2)
                DO i1 = 1, n(1)
                   fgrad_u(i1,i2,i3,1) = CMPLX(0.0_pr, K1(i1))*fu(i1,i2,i3)
                   fgrad_u(i1,i2,i3,2) = CMPLX(0.0_pr, K2(i2))*fu(i1,i2,i3) 
                   fgrad_u(i1,i2,i3,3) = CMPLX(0.0_pr, K3(i3+local_last_start))*fu(i1,i2,i3)
                END DO
             END DO
          END DO
          
          DO ii=1,3
             faux = fgrad_u(:,:,:,ii)
             CALL bfourier(faux, aux)
             grad_u(:,:,:,ii) = REAL(aux)
          END DO

          DEALLOCATE(fu)
          DEALLOCATE(fgrad_u)
          DEALLOCATE(aux)
          DEALLOCATE(faux)
          
        END SUBROUTINE gradient

!        !====================================
!        ! CALCULATE GRAD_PERP
!        !====================================
!        SUBROUTINE grad_perp(u, grad_u)
!          USE global_variables
!          USE FFT2
!          IMPLICIT NONE
!          INCLUDE "mpif.h"
!
!          REAL(pr), DIMENSION(1:n(1),1:local_nlast), INTENT(IN) :: u
!          REAL(pr), DIMENSION(1:n(1),1:local_nlast,1:2), INTENT(OUT) :: grad_u
!
!          COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: fu
!          COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast, 2) :: fgrad_u
!          COMPLEX(pr), DIMENSION(1:n(1), 1:local_nlast) :: aux, faux
!          REAL(pr), DIMENSION(1:n(1)/2, 1:3) :: spectral_data
!          INTEGER, DIMENSION(1:3) :: index_min
!          INTEGER :: i1, i2, i
!          REAL(pr) :: mode           
! 
!          aux = CMPLX(u,0.0_pr)
!          CALL ffourier(aux, fu)
!
!          DO i2 = 1, local_nlast
!             DO i1 = 1, n(1)
!                fgrad_u(i1,i2,1) = CMPLX(0.0_pr, K2(i2+local_last_start))*fu(i1,i2)
!                fgrad_u(i1,i2,2) = -CMPLX(0.0_pr, K1(i1))*fu(i1,i2) 
!             END DO
!          END DO
!          
!          DO i=1,2
!             faux = fgrad_u(:,:,i)
!             CALL bfourier(faux, aux)
!             grad_u(:,:,i) = REAL(aux)
!          END DO
!        
!        END SUBROUTINE grad_perp
!
        !===========================================
        ! CALCULATE DIVERGENCE OF A FIELD
        !===========================================
        SUBROUTINE divergence(u, divU)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: u
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast), INTENT(OUT) :: divU

          COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: u_hat
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: divU_hat
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
          INTEGER :: i1, i2, i3, ii
          REAL(pr) :: local_max_mode, max_mode, mode

          ALLOCATE( u_hat(1:n(1),1:n(2),1:local_nlast,1:3) )
          ALLOCATE( divU_hat(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )
          
          DO ii=1,3
             aux = CMPLX(u(:,:,:,ii), 0.0_pr)
             CALL ffourier(aux, faux)
             u_hat(:,:,:,ii) = faux
          END DO

          DO i3 = 1,local_nlast
             DO i2 = 1,n(2)
                DO i1 = 1,n(1)
                   divU_hat(i1,i2,i3) = CMPLX(0.0_pr, K1(i1))*u_hat(i1,i2,i3,1) + CMPLX(0.0_pr, K2(i2))*u_hat(i1,i2,i3,2) + CMPLX(0.0_pr, K3(i3+local_last_start))*u_hat(i1,i2,i3,3)
                END DO
             END DO
          END DO
          
          CALL bfourier(divU_hat, aux)
          divU = REAL(aux)
        
          DEALLOCATE( u_hat )
          DEALLOCATE( divU_hat )
          DEALLOCATE( aux )
          DEALLOCATE( faux )

        END SUBROUTINE divergence

        !=======================================
        ! CALCULATE LAPLACIAN OF A VECTOR FIELD
        !=======================================
        SUBROUTINE laplacian(f)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: f
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: f_hat, aux
          INTEGER :: i1, i2, i3, ii
          REAL(pr) :: ksq, mode, max_mode, local_max_mode
         
          ALLOCATE( f_hat(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )

          DO ii=1,3
             aux = CMPLX(f(:,:,:,ii), 0.0_pr)
             CALL ffourier(aux, f_hat)

             DO i3 = 1,local_nlast
                DO i2 = 1,n(2)
                   DO i1 = 1,n(1)
                      mode = SQRT(K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2)
                      ksq = K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2
                      f_hat(i1,i2,i3) = -ksq*f_hat(i1,i2,i3)*EXP(-36.0_pr*(mode/Kcut)**36)
                   END DO
                END DO
             END DO

             CALL bfourier(f_hat, aux)
             f(:,:,:,ii) = REAL(aux)
          END DO

          DEALLOCATE( f_hat )
          DEALLOCATE( aux )


        END SUBROUTINE laplacian  

        !==========================================
        ! Function that calculates bilaplacian
        !========================================== 
        SUBROUTINE bilaplacian(u)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: u
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: fu, aux, faux
          INTEGER :: i1, i2, i3, ii
          REAL(pr) :: k4, mode, max_mode, local_max_mode 

          ALLOCATE( fu(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )

          DO ii=1,3
             aux = CMPLX(u(:,:,:,ii), 0.0_pr)
             CALL ffourier(aux, faux)
            
             
             DO i3 = 1,local_nlast
                DO i2 = 1,n(2)
                   DO i1 = 1,n(1)
                      mode = SQRT(K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2) 
                      k4 = ( K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2 )**2
                      fu(i1,i2,i3) = k4*faux(i1,i2,i3)*EXP(-36.0_pr*(mode/Kcut)**36)
                   END DO
                END DO
             END DO

             CALL bfourier(fu, aux)
             u(:,:,:,ii) = REAL(aux)
          END DO 

          DEALLOCATE(fu)
          DEALLOCATE(aux)
          DEALLOCATE(faux) 

        END SUBROUTINE bilaplacian

        !========================================
        ! CALCULATE ADVECTION TERM
        !======================================== 
        SUBROUTINE advection(U, V, W)
          USE global_variables        
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: U, V
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(OUT) :: W

          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: F
          REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: gradF

          INTEGER :: ii

          ALLOCATE( F(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( gradF(1:n(1),1:n(2),1:local_nlast,1:3) )
 
          DO ii=1,3
             F = V(:,:,:,ii)
             CALL gradient(F, gradF) 
             F = U(:,:,:,1)*gradF(:,:,:,1) + U(:,:,:,2)*gradF(:,:,:,2) + U(:,:,:,3)*gradF(:,:,:,3)

             IF (toDealias) CALL dealiasing(F)
             W(:,:,:,ii) = F

          END DO

          DEALLOCATE(F)
          DEALLOCATE(gradF)

        END SUBROUTINE advection

        !===========================================
        ! CALCULATE (nabla U)^T * V = W
        !===========================================
        SUBROUTINE stretching(U,V,W)
          USE global_variables        
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: U, V
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(OUT) :: W

          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: F, G

          INTEGER :: ii, jj

          ALLOCATE( F(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( G(1:n(1),1:n(2),1:local_nlast) )
 
          W = 0.0_pr

          DO ii=1,3
             DO jj=1,3
                F = U(:,:,:,jj)
                CALL derivative(F, ii)
                G = F*V(:,:,:,jj)
                IF (toDealias) CALL dealiasing(G) 
                W(:,:,:,ii) = W(:,:,:,ii) + G
             END DO 
          END DO

          DEALLOCATE(F)
          DEALLOCATE(G)

        END SUBROUTINE stretching 


        !==========================================
        ! PERFORM DEALIASING
        !==========================================
        SUBROUTINE dealiasing(f)
          USE global_variables
          USE FFT2
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast), INTENT(INOUT) :: f

          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
          INTEGER :: i1, i2, i3
          REAL(pr), DIMENSION(1:3) :: k
          REAL(pr) :: mode 

          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )

          aux = CMPLX(f, 0.0_pr)
          CALL ffourier(aux, faux)

          DO i3 = 1,local_nlast
             DO i2 = 1,n(2)
                DO i1 = 1,n(1)
                   k(1) = K1(i1)
                   k(2) = K2(i2)
                   k(3) = K3(i3+local_last_start)
                   mode = SQRT(k(1)**2 + k(2)**2 + k(3)**2)
                   faux(i1,i2,i3) = faux(i1,i2,i3)*EXP(-36.0_pr*(mode/Kcut)**36)
                END DO
             END DO
          END DO

          CALL bfourier(faux,aux)
          f = REAL(aux)

          DEALLOCATE(aux)
          DEALLOCATE(faux)

        END SUBROUTINE dealiasing

        !==============================================
        ! SPATIAL SHIFT: mypoint \to (0.5, 0.5, 0.5)
        !==============================================
        SUBROUTINE shift_space(myfield, mypoint)
          USE global_variables
          USE FFT2
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: myfield
          REAL(pr), DIMENSION(1:3), INTENT(IN) :: mypoint

          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
          REAL(pr), DIMENSION(1:3) :: p_shift, k
          REAL(pr) :: theta
          INTEGER :: nn, i1, i2, i3

          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )

          p_shift(1) = 0.5_pr - mypoint(1)
          p_shift(2) = 0.5_pr - mypoint(2)
          p_shift(3) = 0.5_pr - mypoint(3)
          
          DO nn=1,3
             aux = CMPLX(myfield(:,:,:,nn), 0.0_pr)
             CALL ffourier(aux, faux)
             
             DO i3=1,local_nlast
                DO i2=1,n(2)
                   DO i1=1,n(1)
                      k(1) = K1(i1)
                      k(2) = K2(i2)
                      k(3) = K3(i3+local_last_start)
                      theta = k(1)*p_shift(1) + k(2)*p_shift(2) + k(3)*p_shift(3)
                      faux(i1,i2,i3) = faux(i1,i2,i3)*CMPLX(COS(theta),-SIN(theta))
                   END DO
                END DO
             END DO

             CALL bfourier(faux,aux)
             myfield(:,:,:,nn) = REAL(aux)
          END DO

          DEALLOCATE(aux)
          DEALLOCATE(faux)

        END SUBROUTINE shift_space
 
        !==============================================
        ! SPATIAL SHIFT: mypoint \to (0.5, 0.5, 0.5)
        !==============================================
        SUBROUTINE rotate_space(myfield, mypoint)
          USE global_variables
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: myfield
          REAL(pr), DIMENSION(1:3), INTENT(IN) :: mypoint

          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux_local, aux_global0, aux_global1
          COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: aux_field
          REAL(pr), DIMENSION(1:3) :: p_shift, k
          REAL(pr) :: theta
          INTEGER :: nn, i1, i2, i3

          IF (rank==0) THEN
             ALLOCATE( aux_global0(1:n(1),1:n(2),1:n(3)) )
             ALLOCATE( aux_global1(1:n(1),1:n(2),1:n(3)) )
          END IF
          ALLOCATE( aux_local(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( aux_field(1:n(1),1:n(2),1:local_nlast,1:3) )
          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

          !aux_field = myfield
          !myfield(:,:,:,1) = aux_field(:,:,:,3)
          !myfield(:,:,:,2) = -1.0_pr*aux_field(:,:,:,2)
          !myfield(:,:,:,3) = aux_field(:,:,:,1)

          DO nn=1,3
             aux_local = myfield(:,:,:,nn)
             CALL MPI_GATHER(aux_local, total_local_size, MPI_REAL8, aux_global0, total_local_size, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)   
             
             IF (rank==0) THEN
                DO i3=1,n(3)
                   DO i2=1,n(2)
                      DO i1=1,n(1)
                         aux_global1(i1,i2,i3) = aux_global0(i3, n(2)-i2+1, i1)
                      END DO
                   END DO
                END DO
             END IF
             CALL MPI_SCATTER(aux_global1, total_local_size, MPI_REAL8, aux_local, total_local_size, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)
             
             !IF (nn == 1) THEN
                myfield(:,:,:,nn) = aux_local
             !ELSE
             !   myfield(:,:,:,nn) = -1.0_pr*aux_local
             !END IF

          END DO

          DEALLOCATE(aux_field)
          DEALLOCATE(aux_local)
          IF (rank==0) THEN
             DEALLOCATE(aux_global0)
             DEALLOCATE(aux_global1)
          END IF
          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

        END SUBROUTINE rotate_space

        !========================================================
        ! Find the location in physical space of the first
        ! point where both magU and magW are zero
        ! Use the direction given in myNormal as a mask
        !========================================================
        FUNCTION FindMinPoint(magU,magW,myNormal) RESULT (myPoint)
          USE global_variables
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast), INTENT(IN) :: magU, magW
          REAL(pr), DIMENSION(1:3), INTENT(IN) :: myNormal

          REAL(pr), DIMENSION(1:3) :: myPoint, dx
          REAL(pr) :: res, X, Y, Z          

          INTEGER :: ii, jj, kk
          INTEGER, DIMENSION(1:3) :: min_index
          LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: myMask

          dx = 1.0_pr/REAL(n,pr)
          ALLOCATE( myMask(1:n(1),1:n(2),1:local_nlast) )

          DO kk=1,local_nlast
             DO jj=1,n(2)
                DO ii=1,n(1)
                   X = REAL(ii,pr)*dx(1)
                   Y = REAL(jj,pr)*dx(2)
                   Z = REAL(kk-1+local_last_start,pr)*dx(3)
                   res = myNormal(1)*(X-0.5_pr) + myNormal(2)*(Y-0.5_pr)+myNormal(3)*(Z-0.5_pr)
                   IF (ABS(res) < SQRT(MACH_EPSILON)) THEN
                      myMask(ii,jj,kk) = .TRUE.
                   ELSE
                      myMask(ii,jj,kk) = .FALSE.
                   END IF                     
                END DO
             END DO
          END DO
           
          min_index = MINLOC(magU, myMask)
                    
          DEALLOCATE( myMask )

        END FUNCTION FindMinPoint

        !=======================================================
        ! CALCULATE THE INNER PRODUCT BETWEEN TWO FUNCTIONS
        !=======================================================
        RECURSIVE FUNCTION inner_product(f, g, mytype) RESULT (inn_prod)
          USE global_variables
          USE FFT2
          IMPLICIT NONE
          !INCLUDE "mpif.h"
      
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast), INTENT(IN) :: f, g
          CHARACTER(len=*), INTENT(IN) :: mytype

          REAL(pr) :: real_loc_prod, inn_prod, ksq
          REAL(pr) :: local_inn_prod 

          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: fhat, ghat, aux
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: f_aux
          INTEGER :: i1, i2, i3

          ALLOCATE( fhat(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( ghat(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( f_aux(1:n(1),1:n(2),1:local_nlast) )

          SELECT CASE (mytype)
            CASE ("L2")
               local_inn_prod = 0.0_pr
               DO i3=1,local_nlast
                  DO i2=1,n(2)
                     DO i1=1,n(1)
                        local_inn_prod = local_inn_prod + f(i1,i2,i3)*g(i1,i2,i3)*dV
                     END DO
                  END DO
               END DO
               inn_prod = local_inn_prod
               !CALL MPI_ALLREDUCE(local_inn_prod, inn_prod, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

            CASE ("H1")
               aux = CMPLX(f,0.0_pr)
               CALL ffourier(aux, fhat)
               DO i3=1,local_nlast
                  DO i2=1,n(2)
                     DO i1=1,n(1)
                        ksq = K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2
                        fhat(i1,i2,i3) = (1.0_pr + lambda0**2*ksq)*fhat(i1,i2,i3)
                     END DO
                  END DO
               END DO
               CALL bfourier(fhat, aux)
               f_aux = REAL(aux)
               inn_prod = inner_product(f_aux,g,"L2")

            CASE ("semi_H1")
               aux = CMPLX(f,0.0_pr)
               CALL ffourier(aux, fhat)
               DO i3=1,local_nlast
                  DO i2=1,n(2)
                     DO i1=1,n(1)
                        ksq = K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2
                        fhat(i1,i2,i3) = lambda0**2*ksq*fhat(i1,i2,i3)
                     END DO
                  END DO
               END DO
               CALL bfourier(fhat, aux)
               f_aux = REAL(aux)
               inn_prod = inner_product(f_aux,g,"L2")

            CASE ("H2")
               aux = CMPLX(f,0.0_pr)
               CALL ffourier(aux, fhat)
               DO i3=1,local_nlast
                  DO i2=1,n(2)
                     DO i1=1,n(1)
                        ksq = K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2
                        fhat(i1,i2,i3) = (1.0_pr + lambda0**4*ksq**2)*fhat(i1,i2,i3)
                     END DO
                  END DO
               END DO
               CALL bfourier(fhat, aux)
               f_aux = REAL(aux)
               inn_prod = inner_product(f_aux,g,"L2")

            CASE ("semi_H2")
               aux = CMPLX(f,0.0_pr)
               CALL ffourier(aux, fhat)
               DO i3=1,local_nlast
                  DO i2=1,n(2)
                     DO i1=1,n(1)
                        ksq = K1(i1)**2 + K2(i2)**2 + K3(i3+local_last_start)**2
                        fhat(i1,i2,i3) = lambda0**4*ksq**2*fhat(i1,i2,i3)
                     END DO
                  END DO
               END DO
               CALL bfourier(fhat, aux)
               f_aux = REAL(aux)
               inn_prod = inner_product(f_aux,g,"L2")

          END SELECT

          DEALLOCATE( fhat )
          DEALLOCATE( ghat )
          DEALLOCATE( aux )
          DEALLOCATE( f_aux )

        END FUNCTION inner_product 

        !=======================================================
        ! CALCULATE THE INNER PRODUCT BETWEEN TWO FIELDS
        !=======================================================
        FUNCTION field_inner_product(f,g,mytype) RESULT (inn_prod)
          USE global_variables
          IMPLICIT NONE
      
          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(IN) :: f,g
          CHARACTER(len=*), INTENT(IN) :: mytype
         
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: faux, gaux
          REAL(pr), DIMENSION(1:3) :: inn_prod

          INTEGER :: ii
         
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( gaux(1:n(1),1:n(2),1:local_nlast) )
        
 
          inn_prod = 0.0_pr
          DO ii=1,3
             faux = f(:,:,:,ii)
             gaux = g(:,:,:,ii)
             inn_prod(ii) =  inner_product(faux, gaux, mytype)
          END DO
 
          DEALLOCATE(faux)
          DEALLOCATE(gaux)       
    
        END FUNCTION field_inner_product 

        !======================================================
        ! CALCULATE THE SOBOLEV GRADIENT OF ORDER P, GIVEN
        ! THE L2 GRADIENT
        !======================================================
        SUBROUTINE SobolevGradient(grad, order, ell)
          USE global_variables
          USE FFT2
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_nlast,1:3), INTENT(INOUT) :: grad
          INTEGER, INTENT(IN) :: order
          REAL(pr), INTENT(IN) :: ell

          COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: grad_hat
          COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
  
          REAL(pr) :: ksq
          INTEGER :: ii,jj,kk

          ALLOCATE( grad_hat(1:n(1),1:n(2),1:local_nlast,1:3) )
          ALLOCATE( aux(1:n(1),1:n(2),1:local_nlast) )
          ALLOCATE( faux(1:n(1),1:n(2),1:local_nlast) )

          DO ii=1,3
             aux = CMPLX(grad(:,:,:,ii),0.0_pr)
             CALL ffourier(aux,faux)
             grad_hat(:,:,:,ii) = faux
          END DO

          DO kk=1,local_nlast
             DO jj=1,n(2)
                DO ii=1,n(1)
                   ksq = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_last_start)**2 )
                   grad_hat(ii,jj,kk,:) = grad_hat(ii,jj,kk,:)/(1.0_pr + (ell)**(2*order)*ksq**(2*order) )
                END DO
             END DO
          END DO

          DO ii=1,3
             faux = grad_hat(:,:,:,ii)
             CALL bfourier(faux,aux)
             grad(:,:,:,ii) = REAL(aux)
          END DO 

          DEALLOCATE( grad_hat )
          DEALLOCATE( aux )
          DEALLOCATE( faux )

        END SUBROUTINE SobolevGradient


END MODULE

!-----------------------------------------------------!
! Program used to maximize dEdt                       !
! in the 3D Navier-Stokes system,                     !
! using 2 constraints.                                !
!                                                     !
! Parallel version, Complex-Complex FFT               !   
!                                                     !
! November, 2014.                                     !
!                                                     !
! Author: Diego Ayala                                 !
! Department of Mathematics and Statistics            !
! McMaster University                                 !
!-----------------------------------------------------!

PROGRAM maxdEdt_main

  USE global_variables
  USE data_ops
  USE function_ops
  USE optimization
  IMPLICIT NONE
  INCLUDE "mpif.h" 
  
  INTEGER :: RESOL, ii, jj, kk, E0index, E1index, Kindex, Nuindex, numPts_E0
  REAL(pr) :: aux
  REAL(pr), DIMENSION (:), ALLOCATABLE, SAVE :: K0_vec, E0_vec
  CHARACTER(2) :: K0txt, E0txt, E1txt, IGtxt, NUtxt

 !====================================================

  CALL MPI_INIT(Statinfo)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,Statinfo)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,Statinfo)

  ! Read in runtime parameters
  CALL getarg(1,K0txt)
  CALL getarg(2,E0txt)
  CALL getarg(3,E1txt)
  CALL getarg(4,NUtxt)

  READ(K0txt, '(I2.2)') Kindex 
  READ(E0txt, '(I2.2)') E0index
  READ(E1txt, '(I2.2)') E1index
  READ(NUtxt, '(I2.2)') Nuindex

  IF (Kindex==0) THEN
     ConsType = 1
  ELSE
     ConsType = 2
  END IF

  PI = 4.0_pr*ATAN2(1.0_pr,1.0_pr)
  !=============================================
  ! Create K0_vec and E0_vec
  !============================================= 
  SELECT CASE (ConsType)
     CASE (1)
        K0 = 0.0_pr
        H0 = 1.0_pr
        ALLOCATE(E0_vec(1:65))
        DO ii=1,65
           IF (ii .LE. 2) THEN
              aux = REAL(ii-4,pr)
              E0_vec(ii) = 1.0_pr*10**(aux)
           ELSEIF (ii .LE. 11) THEN
              aux = REAL(ii-2,pr)
              E0_vec(ii) = 0.1_pr*aux
           ELSEIF (ii .LE. 31) THEN
              aux = REAL(ii-11,pr)
              E0_vec(ii) = aux
           ELSEIF (ii .LE. 47) THEN
              aux = REAL(ii-27,pr)
              E0_vec(ii) = 5.0_pr*aux
           ELSEIF (ii .LE. 56) THEN
              aux = REAL(ii-46,pr)
              E0_vec(ii) = 100.0_pr*aux
           ELSE
              aux = REAL(ii-55,pr)
              E0_vec(ii) = 1000.0_pr*aux
           END IF
        END DO
     
     CASE (2)
        numPts_E0 = 12*(2**Kindex-1) + 1

        ALLOCATE(K0_vec(1:4))
        ALLOCATE(E0_vec(1:numPts_E0))

        DO ii=1,4 
           K0_vec(ii) = 10.0_pr**(ii-1)
        END DO
        K0 = K0_vec(Kindex)
 
        DO ii=1,numPts_E0
           aux = REAL(ii-1,pr)/REAL(numPts_E0-1, pr) 
           E0_vec(ii) = 1.1_pr*(2.0_pr*PI)**2.0_pr*K0*(10.0_pr**aux)
        END DO

  END SELECT

  !=========================================================
  ! Loop over different values of E0
  !=========================================================
  DO kk = E0index, E1index
     WRITE(E0txt,'(I2.2)') kk

     !--Read in parameters
     IF (rank==0) THEN
        OPEN (10, FILE="/gwork/ayalada/MAX_ENS_3D/Nu"//NUtxt//"/K"//K0txt//"/E"//E0txt//"/maxdEdt_params.dat", STATUS='OLD')
        READ (10, *) RESOL
        READ (10, *) K0_index
        READ (10, *) E0_index
        READ (10, *) iguess
        READ (10, *) NU_index
        READ (10, *) lambda0
        READ (10, *) alpha0
        CLOSE (10)
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
 
     CALL MPI_BCAST (RESOL, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo)
     CALL MPI_BCAST (K0_index, 1,  MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo)
     CALL MPI_BCAST (E0_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo)
     CALL MPI_BCAST (alpha0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)
     CALL MPI_BCAST (lambda0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo) 
     CALL MPI_BCAST (iguess, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo)
     CALL MPI_BCAST (NU_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo)

     n(1) = RESOL
     n(2) = RESOL
     n(3) = RESOL

     CALL initialize
     E0 = E0_vec(E0_index) 
 
     !=========================================================
     ! Define value of viscosity
     !=========================================================
     IF (NU_index==0) THEN
        visc = 0.0_pr
     ELSE
        visc = 10.0_pr**REAL(-Nuindex,pr)
     END IF

     WRITE(IGtxt,'(i2.2)') iguess
 
     IF (rank==0) THEN
        OPEN(10, FILE="~/RESEARCH/MAX_ENS_3D/LOGFILES/maxdEdt_Nu"//NUtxt//"_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_info.log", STATUS='REPLACE')
        WRITE(10,*) "======================================= "
        WRITE(10,*) "  Resolution N = ", n(1)
        WRITE(10,*) "  Energy K0 = ", K0
        WRITE(10,*) "  Enstrophy E0 = ", E0
        WRITE(10,*) "  Helicity H0 = ", H0
        WRITE(10,*) "  Viscosity = ", visc
        WRITE(10,*) "  Processors = ", np
        WRITE(10,*) "  Initial guess = ", iguess
        WRITE(10,*) "======================================= "
        CLOSE(10)
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
    
     CALL initial_guess 
     CALL maxdEdt
     
     CALL fftwnd_fortran_mpi_destroy_plan (fwdplan)
     CALL fftwnd_fortran_mpi_destroy_plan (bwdplan)
     DEALLOCATE(Uvec)
     DEALLOCATE(Wvec)
     DEALLOCATE(K1)
     DEALLOCATE(K2)
     DEALLOCATE(K3)
     DEALLOCATE(fwork)

  END DO
  CALL MPI_FINALIZE (Statinfo)

END PROGRAM maxdEdt_main


      PROGRAM testing_FFTW3

      USE, INTRINSIC :: iso_c_binding
      INCLUDE "fftw-mpi.f03"

      INTEGER(C_INTPTR_T), DIMENSION(1:3), PARAMETER :: n 

      TYPE(C_PTR) :: pwd_plan, bwd_plan, cdata
      COMPLEX(C_DOUBLE_COMPLEX), POINTER :: u(:,:,:)

      INTEGER(C_INTPTR_T) :: ii, jj, kk, alloc_local, local_nlast, local_last_start        
      REAL(C_DOUBLE) :: x, y, z      


      INTEGER :: rank, np, Statinfo

      CALL MPI_INIT(Statinfo) 
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,Statinfo)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,Statinfo)

      CALL fftw_mpi_init
      
      alloc_local = fftw_mpi_local_size_3d(n(3), n(2), n(1), MPI_COMM_WORLD, local_nlast, local_last_start)
      cdata = fftw_alloc_complex(alloc_local)

      CALL c_f_pointer(cdata, u, [n(1),n(2),local_nlast])

      fwd_plan = fftw_mpi_plan_dft_3d(n(3), n(2), n(1), u, u, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      
      DO kk=1,local_nlast
         DO jj=1,n(2)
            DO ii=1,n(1)
               x = REAL(ii-1)*dV(1)         
               u(ii,jj,kk) = 


      END PROGRAM

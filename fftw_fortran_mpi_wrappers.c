#include <fftw_mpi.h>
#include <fftw.h>
#include <rfftw_mpi.h>
#include <rfftw.h>
//
#if defined(HAVE_MPI_COMM_F2C)
#  define FFTW_MPI_COMM_F2C(comm) MPI_Comm_f2c(*((MPI_Comm *) comm))
#elif defined(FFTW_USE_F77_MPI_COMM)
#  define FFTW_MPI_COMM_F2C(comm) (* ((MPI_Comm *) comm))
#elif defined(FFTW_USE_F77_MPI_COMM_P)
#  define FFTW_MPI_COMM_F2C(comm) ((MPI_Comm) comm)
#else
#  define FFTW_MPI_COMM_F2C(comm) MPI_COMM_WORLD
#endif

/*
#define FFTW_FORTRANIZE_LOWERCASE
#define FFTW_FORTRANIZE_EXTRA_UNDERSCORE
*/

#ifdef FFTW_FORTRANIZE_LOWERCASE
#  ifdef FFTW_FORTRANIZE_EXTRA_UNDERSCORE
#    define FORTRAN_FUNC_(x,X) x ## _
#  else
#    define FORTRAN_FUNC_(x,X) x
#  endif
#endif

#ifdef FFTW_FORTRANIZE_LOWERCASE_UNDERSCORE
#  ifdef FFTW_FORTRANIZE_EXTRA_UNDERSCORE
#    define FORTRAN_FUNC_(x,X) x ## __
#  else
#    define FORTRAN_FUNC_(x,X) x ## _
#  endif
#endif

#ifdef FFTW_FORTRANIZE_UPPERCASE
#  ifdef FFTW_FORTRANIZE_EXTRA_UNDERSCORE
#    define FORTRAN_FUNC_(x,X) X ## _
#  else
#    define FORTRAN_FUNC_(x,X) X
#  endif
#endif

#ifdef FFTW_FORTRANIZE_UPPERCASE_UNDERSCORE
#  ifdef FFTW_FORTRANIZE_EXTRA_UNDERSCORE
#    define FORTRAN_FUNC_(x,X) X ## __
#  else
#    define FORTRAN_FUNC_(x,X) X ## _
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

/************************************************************************/

void FORTRAN_FUNC_(fftwnd_fortran_mpi_create_plan,FFTWND_FORTRAN_MPI_CREATE_PLAN)
(fftwnd_mpi_plan *p, void *comm, int *rank, int *n, int *idir, int *flags)
{
     fftw_direction dir = *idir;
     int i,tmp,nDim = *rank;

     for (i=nDim/2-1;i>=0;i--) {
	   tmp = n[i]; 
	   n[i] = n[nDim-1-i];
	   n[nDim-1-i] = tmp;
	 }
     *p = fftwnd_mpi_create_plan(FFTW_MPI_COMM_F2C(comm),
				  *rank, n, dir, *flags);
     for (i=nDim/2-1;i>=0;i--) {
	   tmp = n[i]; 
	   n[i] = n[nDim-1-i];
	   n[nDim-1-i] = tmp;
	 }
}

void FORTRAN_FUNC_(fftw2d_fortran_mpi_create_plan,FFTW2D_FORTRAN_MPI_CREATE_PLAN)
(fftwnd_mpi_plan *p, void *comm, int *nx, int *ny, int *idir, int *flags)
{
     fftw_direction dir = *idir;

     *p = fftw2d_mpi_create_plan(FFTW_MPI_COMM_F2C(comm), *ny,*nx,dir,*flags);
}

void FORTRAN_FUNC_(fftw3d_fortran_mpi_create_plan,FFTW3D_FORTRAN_MPI_CREATE_PLAN)
(fftwnd_mpi_plan *p, void *comm, 
 int *nx, int *ny, int *nz, int *idir, int *flags)
{
     fftw_direction dir = *idir;

     *p = fftw3d_mpi_create_plan(FFTW_MPI_COMM_F2C(comm), 
				  *nz,*ny,*nx,dir,*flags);
}

void FORTRAN_FUNC_(fftwnd_fortran_mpi_destroy_plan,FFTWND_FORTRAN_MPI_DESTROY_PLAN)
(fftwnd_mpi_plan *p)
{
     fftwnd_mpi_destroy_plan(*p);
}

void FORTRAN_FUNC_(fftwnd_fortran_mpi,FFTWND_FORTRAN_MPI)
(fftwnd_mpi_plan *p, int *n_fields, fftw_complex *local_data,
 fftw_real *work, int *ioutput_order)
{
     fftwnd_mpi_output_order output_order = *ioutput_order;

     fftwnd_mpi(*p, *n_fields, local_data, NULL, output_order);
}

void FORTRAN_FUNC_(fftwnd_fortran_mpi_local_sizes,FFTWND_FORTRAN_MPI_LOCAL_SIZES)
(fftwnd_mpi_plan *p,
 int *local_nx, int *local_x_start,
 int *local_ny_after_transform,
 int *local_y_start_after_transform,
 int *total_local_size)
{
     fftwnd_mpi_local_sizes(*p, local_nx, local_x_start,
			     local_ny_after_transform, 
			     local_y_start_after_transform,
			     total_local_size);
}

/****************************************************************************/

#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */


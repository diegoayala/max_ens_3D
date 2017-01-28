MODULE global_variables
  IMPLICIT NONE
  INCLUDE "mpif.h"

  INTEGER, PARAMETER :: pr = 8 !KIND(0.0d0)
  INTEGER, PARAMETER :: MAX_ITER = 1000
  INTEGER, PARAMETER :: KappaPoints = 16
  REAL, PARAMETER :: OPTIM_TOL = 1.0e-6
  REAL, PARAMETER :: MACH_EPSILON = 1.0e-16 
  REAL, PARAMETER :: visc = 1.0e-2
  REAL, PARAMETER :: TAU_MAX = 100

  LOGICAL :: kappaTest
  LOGICAL :: toDealias 
  LOGICAL :: save_diag_NS
  LOGICAL :: save_data_NS
  LOGICAL :: save_diag_Constr
  LOGICAL :: save_data_Constr
  LOGICAL :: save_diag_Optim
  LOGICAL :: save_data_Optim
  LOGICAL :: save_diag_lineMin
  LOGICAL :: save_data_lineMin
  LOGICAL :: parallel_data

 
  INTEGER, DIMENSION(3), SAVE :: n
  INTEGER, SAVE :: n_dim
  INTEGER, SAVE :: K0_index, E0_index, iguess, ConsType
  REAL(pr), SAVE :: E0, K0, PI, lambda0, alpha0, dV, Kcut, Kmax
  !--NOTE: Kcut = cut frequency used for dealiasing. 
  !-       Kmax = maximal frequency present in solution.

  REAL(pr), DIMENSION (:), ALLOCATABLE, SAVE :: K1, K2, K3
  REAL(pr), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Ux, Uy, Uz
  REAL(pr), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Wx0, Wy0, Wz0
  REAL(pr), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Wx1, Wy1, Wz1
  COMPLEX(pr), DIMENSION (:), ALLOCATABLE, SAVE :: fwork
  
  !========================================================================== 
  !                            MPI VARIABLES
  !==========================================================================
  INTEGER, SAVE :: rank, Statinfo, np
  INTEGER, SAVE :: local_nlast, local_last_start, local_nlast_after_trans 
  INTEGER, SAVE :: local_last_start_after_trans, total_local_size 
  INTEGER, SAVE :: local_start, local_n

  INTEGER, SAVE :: local_nlastFixres, local_last_startFixres, local_nlast_after_transFixres
  INTEGER, SAVE :: local_last_start_after_transFixres, total_local_sizeFixres, local_startFixres, local_nFixres

  INTEGER, SAVE :: local_nlastHighres, local_last_startHighres, local_nlast_after_transHighres
  INTEGER, SAVE :: local_last_start_after_transHighres, total_local_sizeHighres, local_startHighres, local_nHighres


  INTEGER, SAVE :: NF90_mytype
 
  !========================================================================== 
  !                            FFTW VARIABLES
  !==========================================================================
  INTEGER, PARAMETER :: FFTW_FORWARD=-1, FFTW_BACKWARD=1
  INTEGER, PARAMETER :: FFTW_REAL_TO_COMPLEX=-1, FFTW_COMPLEX_TO_REAL=1
  INTEGER, PARAMETER :: FFTW_ESTIMATE=0, FFTW_MEASURE=1
  INTEGER, PARAMETER :: FFTW_OUT_OF_PLACE=0, FFTW_IN_PLACE=8
  INTEGER, PARAMETER :: FFTW_USE_WISDOM=16, FFTW_THREADSAFE=128
  INTEGER, PARAMETER :: FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0
  INTEGER, PARAMETER :: FTW_SCRAMBLED_INPUT=8192, FFTW_SCRAMBLED_OUTPUT=16384
  INTEGER(8), SAVE  :: fwdplan, bwdplan, fwdplan_Fixres, bwdplan_Fixres, fwdplan_Highres, bwdplan_Highres

END MODULE

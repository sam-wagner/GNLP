PROGRAM DRIVER
USE CONSTANTS
USE MATH_MODULE
USE ORBITAL_FUNCTIONS
USE OPTIMIZATION_MODULE
use omp_lib
IMPLICIT NONE

INTEGER :: N_POP, N_GEN, N_INT, N_DOUBLE, N_RUNS, N1,N2, NCON, ITER_MAX_MBH, ITER_MAX_NLP, N_CON
INTEGER, ALLOCATABLE :: INTEGER_UPPER(:)
INTEGER, ALLOCATABLE :: INTEGER_LOWER(:)
DOUBLE PRECISION :: P_CROSS, P_REP,P_MUT 
DOUBLE PRECISION, ALLOCATABLE :: DOUBLE_UPPER(:), DOUBLE_LOWER(:), INPUT_ARRAY(:,:)
CHARACTER(LEN=30) :: CROSSOVER_TYPE, MUTATION_TYPE, SELECTION_TYPE, OPTIMIZATION_TYPE

INTEGER, ALLOCATABLE :: INTEGER_RUN(:,:), CHROM_INT(:)
DOUBLE PRECISION, ALLOCATABLE :: MIN_RUN(:), AVG_RUN(:), DOUBLE_RUN(:,:), LAST_HALF_CHANGE(:)

INTEGER :: SEED,  COUNT1, COUNT2, RATE, I, J, MIN_LOC, IPRINT, NGEN_CONVERGE
DOUBLE PRECISION :: TIME , FITNESS_INDV, MAX_TIME, TOL_CONVERGE
DOUBLE PRECISION, ALLOCATABLE :: FITNESS_AVG(:), FITNESS_MIN(:), DOUBLE_MIN(:,:), X0(:), X(:)
INTEGER, ALLOCATABLE :: INTEGER_MIN(:,:)


!SEED FOR FIRST RANDOM NUMBER GENERATION
SEED=-54525585



!IF AN INPUT ARRAY OF INFORMATION NEEDS TO BE PASSED TO THE COST FUNCTION USE
!INPUT_ARRAY WHICH HAS A SIZE OF N1, N2





!PROBLEM INPUTS
P_CROSS=0.9D0
P_REP=1.D0-P_CROSS
P_MUT=0.5D0
N_POP=10
N_GEN=1000000
N_INT=1
N_DOUBLE=6
N_CON=0
N_RUNS=10
IPRINT=1
NGEN_CONVERGE=50
MAX_TIME=1000.D0
TOL_CONVERGE=1.D-5


N1=n_runs
N2=n_int+1
ALLOCATE(INPUT_ARRAY(N1,N2))
INPUT_ARRAY=0.D0
!MAXIMUM NUMBER OF ITERATIONS EACH NLP SOLVER CALL IS ALLOWED
ITER_MAX_NLP=100
! MAXIMUM NUMBER OF BAD SOLUTIONS ALLOWED WHEN THE MONOTONIC BASIN HOPPING OPTIMIZATION METHOD IS USED
ITER_MAX_MBH=5000


!OPTIMIZATION TYPES INCLUDE:
!   GEN - PURE GENETIC ALGORITHM OPTIMIZATION
!   HYB - HYBRID ALGORITHM THAT ITERATES ON THE COST FUNCTION WITH A TRADITIONAL 
!         NLP SOLVER
!   ALT - ALTERNATES 10 GENERATIONS OF PURE GENETIC ALGORITHM AND 10 GENERATIONS 
!         OF THE HYBRID ALGORITHM
OPTIMIZATION_TYPE="HYB_COBYLA"
!OPTIMIZATION_TYPE="HYB_CONMIN"
!OPTIMIZATION_TYPE="HYB_UNCMIN"
!OPTIMIZATION_TYPE="GEN"

!OPTIMIZATION PARAMETERS TYPES
!CROSSOVER_TYPE="DOUBLE_POINT"
CROSSOVER_TYPE="ARITHMETIC"
!CROSSOVER_TYPE="UNIFORM"


MUTATION_TYPE="UNIFORM"
SELECTION_TYPE="ROULETTE"
!SELECTION_TYPE="TOURNAMENT"

!ALLOCATE VARIABLE LIMIT ARRAYS HERE
ALLOCATE(INTEGER_LOWER(N_INT))
ALLOCATE(INTEGER_UPPER(N_INT))

ALLOCATE(DOUBLE_LOWER(N_DOUBLE))
ALLOCATE(DOUBLE_UPPER(N_DOUBLE))


!INTEGER INPUT BOUNDS
INTEGER_UPPER(1)=3
INTEGER_LOWER(1)=3




!LAUNCH DATE
DOUBLE_UPPER=10.D0

DOUBLE_LOWER=0.D0





!********************************************************************************************!
! START TIMER
CALL SYSTEM_CLOCK(COUNT1,RATE)

ALLOCATE(FITNESS_AVG(N_GEN), FITNESS_MIN(N_GEN))
ALLOCATE(INTEGER_MIN(N_GEN, N_INT))
ALLOCATE(DOUBLE_MIN(N_GEN, N_DOUBLE))
ALLOCATE(INTEGER_RUN(N_RUNS, N_INT), MIN_RUN(N_RUNS), AVG_RUN(N_RUNS))
ALLOCATE(DOUBLE_RUN(N_RUNS,N_DOUBLE), LAST_HALF_CHANGE(N_RUNS))
ALLOCATE(X0(N_DOUBLE), X(N_DOUBLE), CHROM_INT(N_INT))

DO I=1,N_RUNS,1
    INPUT_ARRAY=0.D0
    
    ! START GENETIC ALGORITHM
CALL GNLP_DRIVER(IPRINT, N_POP, N_GEN, N_INT, N_DOUBLE, N1, N2, ITER_MAX_NLP, &
    N_CON, INTEGER_UPPER, INTEGER_LOWER, P_CROSS, P_REP,P_MUT, DOUBLE_UPPER, &
    DOUBLE_LOWER, INPUT_ARRAY, CROSSOVER_TYPE, MUTATION_TYPE, &
    SELECTION_TYPE, OPTIMIZATION_TYPE, SEED, FITNESS_MIN, FITNESS_AVG, &
    INTEGER_MIN, DOUBLE_MIN, MAX_TIME, NGEN_CONVERGE, TOL_CONVERGE)
    
    
    CALL SYSTEM_CLOCK(COUNT2)
    TIME=DBLE(COUNT2-COUNT1)/DBLE(RATE)
    MIN_RUN(I)=FITNESS_MIN(N_GEN)
    
    DOUBLE_RUN(I,1:N_DOUBLE)=DOUBLE_MIN(N_GEN,1:N_DOUBLE)
    INTEGER_RUN(I,1:N_INT)=INTEGER_MIN(N_GEN,1:N_INT)
    input_array(I,1:n_int)=dble(integer_run(I,1:n_int))
    input_array(I,n_int+1)=min_run(i)
    !WRITE(*,*) INPUT_ARRAY(I,1:N2)
    !WRITE(*,*) TIME, I, MIN_RUN(I), INTEGER_RUN(I,1:N_INT)
    OPEN(1,FILE='GENETIC_OUTPUT_SHEKEL.TXT',ACCESS='APPEND')
    WRITE(1,*) INT(INPUT_ARRAY(1,2)), INTEGER_RUN(i,1:N_INT), MIN_RUN(I), DOUBLE_RUN(i,1:N_DOUBLE)
    close(1)
END DO

!GET TIME

CALL SYSTEM_CLOCK(COUNT2)
TIME=DBLE(COUNT2-COUNT1)/DBLE(RATE)


!PERFORM CLEAN UP OPERATIONS HERE.  THIS ISN'T REQUIRED, BUT IT KEEPS THE INTEL ANALYZER HAPPY
DEALLOCATE(INTEGER_LOWER, INTEGER_UPPER, DOUBLE_LOWER, DOUBLE_UPPER)
DEALLOCATE(FITNESS_AVG, FITNESS_MIN, INTEGER_MIN, DOUBLE_MIN, INPUT_ARRAY)

END PROGRAM DRIVER
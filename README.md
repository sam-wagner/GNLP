Hybrid Genetic-Nonlinear Programming Optimization Algorithm
====
The GNLP optimization algorithm is meant solve optimization problems of the form:

<img src="https://github.com/sam-wagner/GNLP/blob/master/min.png" width="100px" height="25px" />

Subject to:

<img src="https://github.com/sam-wagner/GNLP/blob/master/constraints.png" width="100px" height="25px" />

However, the provided GNLP optimization routines can be used for both constrained and unconstrained optimization problems.  Multiple NLP solvers are included to flexibly allow this capability.

## Introduction

The GNLP optimization algorithm is meant to globally optimize mixed integer problems (it can be used for both integer and real valued variable problems).  This algorithm was originally design for space mission design problems, such as optimizing multiple gravity-assist and low-thrust interplanetary mission analysis.  In additional to being utilized for these types of problems the GNLP algorithms has been tested on a variety benchmark optimization problems.  The cost module and driver files for many of these problems can be found in the benchmark folder.  A variety of mission optimization problem cost and driver files are included in the astrodynamics folder.  In general the provided cost function files and driver files have corresponding names for each problem set.


Evolutionary algorithms, particularly genetic algorithms, are well suited for global optimization. However, when genetic algorithms are used to optimize impulse multiple gravity-assist and low thrust problems they often only find solutions near the global optimum. Alternatively, non-linear programming (NLP) solvers typically only converge to locally optimal solutions. By modifying the genetic algorithm to utilize an NLP solver to determines locally optimal solutions a robust global optimization algorithm can be developed. This hybrid algorithm, known as the GNLP optimization algorithm, is able determine near globally optimal solutions by combining the global convergence of the genetic algorithm with accuracy of the nonlinear programming algorithm. The proposed algorithm is able to efficiently solve complex problems by significantly reducing the population size and number of generations required to converge on near-globally optimal solutions.

The GNLP algorithm should only be used to optimize functions that are continuous and at least twice differentiable, at least in the neighborhood of the proposed solution. Because of this the NLP solver does not iterate on integer variables, which often introduce large discontinuities. The genetic algorithm is used exclusively to optimize integer variables. These properties make optimizing both high an low thrust mission design problems good candidates for the proposed
GNLP optimization algorithm.

The following freely available NLP solvers were used:
* UNCMIN
* CONMIN
* COBYLA

Each of the NLP solvers was modified to ensure that a consisted objective/constraint function calling method was used across each of the.  The modified files have been included, but it should be noted that the author did not write the original versions of any of these solvers.  If any of them are used, their individual references should be included in any publications.

More information the the GNLP solver, including details of the genetic algorithm, and the NLP solvers can be found in the listed references.


## Variables and Structure of the GNLP Driver

Values the user must set in the driver prior to calling the GNLP algorithm are listed below (inputs).

#### Integer input values:

IPRINT - 0 prints no output except the final solution to the screen, while anything greater will print the best solution found for each generation.

N_POP - Size of the population used by the algorithms.  This should be at least 4 (typically much larger than this) and must be an even number.  

N_GEN - Maximum number of generations the algorithm is allowed to process.

N_INT - Number of integer variables.  This should be set to at leats 1, even for strictly real valued problems.

N_DOUBLE - Number of real valued variables.  This should be set to atleast 1, even for strictly integer problems.

N1 - Width of  INPUT_ARRAY 

N2 -  Length of INPUT_ARRAY 

ITER_MAX_NLP - The maximum number of iterations the various NLP solvers are allowed to perform for any one single solution.  For problems with computationally expensive objective functions this value should be kept low, perhaps 10-20 or less iterations.  For many of the example problems provided it often works best to set the value to 150 or more so the NLP solver will converge on the best solution.

N_CON - Number of problem constraints. This value should be set to 0 if no constraints are used or if the UNCMIN solver is used.  One important note is that this variable must be the same as the number of constraint values returned with the objective function.

INTEGER_UPPER - Upper bounds for the integer variables. This array must have a length of N_INT.

INTEGER_LOWER - Lower bounds for the integer variables. This array must also have a length of N_INT. 

SEED - Seed valued for the random number generator, which should be initialed to a negative integer.

NGEN_CONVERGE - Number of generations the algorithm can stagnate on a solution before the optimization routine exits.  This number is problem specific, but as a rule of thumb 50 works well for many problems.

 

####Double precision input values:

P_CROSS - Probability that a crossover will occur, whould typically have a value close to 0.9.

P_REP - Probability that reproduction will occur (the two parents are simply inserted into the next generation).  This variable should have a value of 1.0-P_CROSS

P_MUT - Probability that a mutation will occur.  The recommended value in literature is often around 0.1, but higher values can sometimes lead to optimal solutions much more quickly (at the risk of making the search somewhat more "random").

DOUBLE_UPPER - Real valued upper bounds for the problem variables.  This array should have a length of N_DOUBLE

DOUBLE_LOWER - Real valued lower bounds for the problem variables.  This array should have a length of N_DOUBLE

INPUT_ARRAY   =   This array should be used for any additional inputs the cost function may require, for example for the astrodynamics based problems it is often used for the ephemeris data the cost function requires. For the MGA/MGA-DSM Galileo and Cassini type space missions this array was used as the "database" to ensure the same flyby trajectory isn't optimized more than once.  This array must have the size of N1xN2.

MAX_TIME - This is the maximum time any individual call of the optimization routine is allowed to run in seconds.  For unlimited time simply set it to a very large value.  If this time has been exceded at the end of a generation the current best solution will be return.

TOL_CONVERGE - The convergence tolerance the solution must find for NGEN_CONVERGENCE generations.  Once the solution hasn't changed for the set number of generations the algorithm exits and outputs the solution it obtained.  This tolerance will be problem specific, but a value of 1.d-5 typically works well.

#### Character Array Inputs:

All of the following inputs should have the character length set to 30.

CROSS_TYPE - Crossover type to be used.  The options are UNIFORM, SINGLE_POINT, DOUBLE_POINT, ARITHMETIC, and HEURISTIC.  As a default either the DOUBLE_POINT or ARITHMETIC crossover operators are recommended.

MUT_TYPE - Type of mutation to be used during the optimization process.  The options are: UNIFORM, SLIDING, and BOUNDARY.  As a default it is probably best to use the UNIFORM mutation operator.

SEL_TYP - Selection type to be used.  The options are ROULETTE and TOURNAMENT.  The TOURNAMENT selection method often works best, but it is a greedy overselection method that can result in early stagnation.  This this happens ROULETTE selection should be used.

OPT_TYPE - The type of local optimization to be used.  The options are: GEN, HYB_COBLYA, HYB_CONMIN, and HYB_UNCMIN.  For unconstrained minimization problems the HYB_UNCMIN optimization type should be used because it is typically faster than the conmin or cobyla NLP solvers.  To skip the local optimization step and rely strictly on the genetic algorithm for optimization the the GEN optimization type.  The NLP solvers utilized in this algorithm only operate on real values variables, so the GEN optimization option should be used for any combinatorial problems.

#### GNLP Outputs:

These output arrays must be allocated prior to call the GNLP_DRIVER subroutine.

FITNESS_MIN - Array of with a length of N_GEN that contains the best fitness values obtained for each generation.  This must  be a double precision array.

FITNESS_AVG - Array of with a length of N_GEN that contains the average fitness values obtained for each generation, which is useful to tell if stagnation has occurred.  This must be be a double precision array.

INTEGER_MIN - Integer array of size N_GENxN_INT that contains the integer part of the best solution found for each generation.

DOUBLE_MIN - Double precision array of size N_GENxN_DOUBLE that contains the real part of the best solution found for each generation.  The final best solution set will always be the last element of each of the output arrays.


After setting all of the problem input  variables and initializes the output variables the GNLP driver is called as follows:

```fortran
CALL GNLP_DRIVER(IPRINT, N_POP, N_GEN, N_INT, N_DOUBLE, N1, &
    N2, ITER_MAX_NLP, N_CON, INTEGER_UPPER, INTEGER_LOWER, P_CROSS, &
    P_REP,P_MUT, DOUBLE_UPPER, DOUBLE_LOWER, INPUT_ARRAY, CROSS_TYPE,&
    MUT_TYPE, SEL_TYPE, OPT_TYPE, SEED, FITNESS_MIN, FITNESS_AVG, &
    INTEGER_MIN, DOUBLE_MIN, MAX_TIME, NGEN_CONVERGE, TOL_CONVERGE)
```

It should be noted that the GNLP solver uses a stoichastic methods, so it should be called more than once.  In the example driver files there is a variable called N_RUNS, which is set at the same time as the problem definition variables.  The loop that runs the GNLP_DRIVER function is the best place to parallelize the solver (basically run the solver several times at once) with somthing like OpenMP.  There are also comment in the GNLP_DRIVER subroutine that can be used to parallelize individual runs of the routine with OpenMP.  This should be used with caution, because this can actually cause the algorithm to run slower if the objective function isn't properly set up for parallel operations. 

If there is any confusion as to how the GNLP driver is called, please examine one of the benchmark problem driver program files.

### Structure of the Cost Functions

The objective function must be in it's own module title "COST_MODULE".  The objective function and module should have the following structure: 

```fortran
MODULE COST_MODULE
! USE (INSERT ANY MODULES REQUIRED FOR THE OBJECTIVE FUNCTION HERE)
IMPLICIT NONE


CONTAINS

SUBROUTINE COST(N_DOUBLE, N_INT, N1, N2, CHROM_DOUBLE, CHROM_INT, FITNESS, ARRAY, G_CON, NCON)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N_INT, N_DOUBLE, N1, N2, CHROM_INT(N_INT), NCON
DOUBLE PRECISION, INTENT(IN) :: CHROM_DOUBLE(N_DOUBLE)
DOUBLE PRECISION, INTENT(INOUT) :: FITNESS, G_CON(NCON), ARRAY(N1,N2)

...
FITNESS= SOME FUNCTION OF THE INTEGER AND REAL VARIABLES

! IF THERE ARE CONSTRAINTS THE VALUE(S) MUST BE SET PRIOR TO EXITING THE COST FUNCTION
! IF THERE ARE NOT CONSTRAINTS A VALUE SHOULDN'T BE SET.

G_CON(1)= SOME FUNCTION OF THE INTEGER AND REAL VARIABLES
.
.
.
G_CON(NCON)= SOME FUNCTION OF THE INTEGER AND REAL VARIABLES
END SUBROUTINE COST

END MODULE COST_MODULE
```

The user should take extra care to ensure the COST subroutine never returns a NAN, as they will cause problems with the algorithm.  It's typically a good idea to check for an NAN at the end of the cost function.  If one is encountered simply set the fitness value to a large number and the GNLP algorithm will remove the solution from the population as the generations progress.


### Compiling and Running the Optimization Routine

The following files required for the GNLP Global Optimization Solver and should be compiled in the order they appear here:
* cost_module_XXX.f90 
* conmin_ifc.f or conmin_ifc.f90 (This fortran 90  version was recently converted to fortran 90 because I had trouble with the fortran 90 version on CONMIN provided on Alan Miller's fortran page, so it is still somewhat experimental but has worked properly for all the tests I've performed.)
* cobyla.f90 
* uncmin.f90
* optimization_module.f90
* driver_XXX.f90

Note that any other modules necessary for the cost module file should be compiled first.

The XXX indicates the user files to define the specific problem, which consist of the cost function (in it's own module) and the main function/driver where the GNLP input variables are defined a and the optimization process is initiated.






## Guarantee 
This code is provided as is with no guarantees.  However if you come across any major flaws please let me know by emailing details of the problem to gnlp_bugs@outlook.com.  If possible a minimum working example will help the debugging process.

## References
This code can be used and modified however you like, the author simply asks that you reference the following in any published works:

1. Wagner, Samuel, "Automated trajectory design for impulsive and low thrust interplanetary mission analysis" (2014), Graduate Theses and Dissertation, Paper XXXXX.

In addition, the following references for the modified NLP solvers should also be referenced in any published works.

UNCMIN:

2.  Schnabel, R., Koontz, J., and Weiss, B. A Modular System of Algorithms for Unconstrained Minimization. University of Colorado at Boulder: Department of Computer Science, 1985.

3. Kahaner, D., Moler, C., and Nash, S. Numerical Methods and Software. Prentice Hall Series in Computational Mathematics, Englewood Cliffs, New Jersey 07632, 2nd edition, 1989.

CONMIN:

4. Vanderplaats, G.N. Conmin User’s Manual . Technical Report X-62282, NASA Technical Memorandum , 1978.

5. Vanderplaats, G. N., and Moses, F. Structural Optimization by Methods of Feasible Directions. National Symposium on Computerized Structural Analysis and Design, March 1972.

COBYLA:

6. Powell, M.J.D. A Direct Search Optimization Method That Models the Objective and Constraint Functions by Linear Interpolation. In Advances in Optimization and Numerical Analysis, volume 275 of Mathematics and Its Applications, pages 51–67. Springer Netherlands, 1994.

7. Powell, M.J.D. A View of Algorithms for Optimization without Derivatives. Technical Report DAMTP 2007/NA03, Cambridge University, 2007.


Hybrid Genetic-Nonlinear Programming Optimization Algorithm
====



## Introduction

The GNLP optimization algorithm is meant to globally optimize mixed integer problems (it can be used for both integer and real valued variable problems).  This algorithm was originally design for space mission design problems, such as optimizing multiple gravity-assist and low-thrust interplanetary mission analysis.  In additional to being utilized for these types of problems the GNLP algorithms has been tested on a variety benchmark optimization problems.  The cost module and driver files for many of these problems can be found in the benchmark folder.  A variety of mission optimization problem cost and driver files are included in the astrodynamics folder.  In general the cost function files and driver files have corresponding names for each problem set.

[min]: https://github.com/sam-wagner/GNLP/min.png


Values the user must set in the driver prior to calling the GNLP algorithm:






Files required for the GNLP Global Optimization Solver:
* cost_module_XXX.f90 
* conmin_ifc.f (or conmin_ifc.f90)
* coblyla.f90 
* uncmin.f90
* optimization_module.f90
* driver_XXX.f90

Where the XXX indicates the user files to define the problem, which are the cost function (in it's own module) and the driver where the GNLP input variables are defined a and the optimization process is initiated.

```fortran
CALL GENETIC_DRIVER(IPRINT, N_POP, N_GEN, N_INT, N_DOUBLE, N1, &
    N2, ITER_MAX_NLP, N_CON, INTEGER_UPPER, INTEGER_LOWER, P_CROSS, &
    P_REP,P_MUT, DOUBLE_UPPER, DOUBLE_LOWER, INPUT_ARRAY, CROSS_TYPE,&
    MUT_TYPE, SEL_TYPE, OPT_TYPE, SEED, FITNESS_MIN, FITNESS_AVG, &
    INTEGER_MIN, DOUBLE_MIN, MAX_TIME, NGEN_CONVERGE, TOL_CONVERGE)
```





This code is provided as is with no guarantees.  However if you come across any major flaws please let me know by emailing details of the problem to gnlp_bugs@outlook.com.  If possible a minimum working example will help the debugging process.

This code can be used and modified however you like, the author simply asks that you reference the following in any published works:

1. Wagner, Samuel Wagner, "Automated trajectory design for impulsive and low thrust interplanetary mission analysis" (2014), Graduate Theses and Dissertation, Paper XXXXX.

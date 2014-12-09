GNLP
====

Values the user must set:






Files required for the GNLP Global Optimization Solver:
* cost_module_XXX.f90 
* conmin_ifc.f (or conmin_ifc.f90)
* coblyla.f90 
* uncmin.f90
* optimization_module.f90
* driver_XXX.f90

Where the XXX indicates the user files to define the problem, which are the cost function (in it's own module) and the driver where the GNLP input variables are defined a and the optimization process is initiated.

```fortran
SUBROUTINE GENETIC_DRIVER(IPRINT, N_POP, N_GEN, N_INT, N_DOUBLE, N1, &
    N2, ITER_MAX_NLP, N_CON, INTEGER_UPPER, INTEGER_LOWER, P_CROSS, &
    P_REP,P_MUT, DOUBLE_UPPER, DOUBLE_LOWER, INPUT_ARRAY, CROSS_TYPE,&
    MUT_TYPE, SEL_TYPE, OPT_TYPE, SEED, FITNESS_MIN, FITNESS_AVG, &
    INTEGER_MIN, DOUBLE_MIN, MAX_TIME, NGEN_CONVERGE, TOL_CONVERGE)
```





This code is provided as is with no guarantees.  However if you come across any major flaws please let me know by emailing me at gnlp_bugs@outlook.com.

This code can be used and modified however you like, the author simply asks that you reference the following in any published works:

Wagner, Samuel Wagner, "Automated trajectory design for impulsive and low thrust interplanetary mission analysis" (2014), Graduate Theses and Dissertation, Paper XXXXX.

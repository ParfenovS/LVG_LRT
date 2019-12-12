#pragma once

#define MIN_TAU -30.0                                   // minimum allowed optical depth
#define BETA_ACCURACY 1.e-8                             // relative accuracy that is used for integration of escape probability over angles in beta_LVG.h file
#define BETA_NUM_SPLINE_POINTS 100                      // a half of the number of points that will be used to construct the spline approximating the dependence of escape probability on optical depth (beta_LVG.h)
#define MIN_POP_FOR_DIFF_CALC 1.e-14                    // the levels with population below this limit will not be taken into account in calculations of the maximum difference between lev. pops. during local fixed point iterations
#define MIN_POP 1.e-30                                  // the populations can't be lower than this value
#define MAX_POP (1.0-1000*DBL_EPSILON)                  // the populations can't be higher than this value
#define Ng_order 2                                      // the order of Ng acceleration
#define Ng_start (2 * Ng_order)                         // the iteration number from which Ng acceleration will be switched on
#define Ng_step (Ng_order + 1)                          // Ng accelaration will be performed each Ng_step iteration
#define DoNg true                                       // whether to use Ng acceleration
#define UNDER_RELAX_FACTOR 0.3                          // underrelaxation factor used while solving statistical equilibrium equations
#define MAX_LOCAL_ACCURACY 1.e-9                        // maximum relative local error used at each step of integration over time
#define MAX_NUM_INNER_STEPS 10000                       // maximum number of local iterations to achive LOCAL_ACCURACY
#define INITIAL_TIME_STEP 1.e-2                         // the initial value of the time step used for integration over time, [s]
#define MIN_TIME_STEP 1.e-5                             // minimum time step used for itegration over time, [s]
#define MAX_TIME_STEP 10000.0                           // maximum time step used for itegration over time, [s]
#define TIME_STEP_DECREASE_FAC 0.5                      // the time step can be decreased by this factor
#define TIME_STEP_INCREASE_FAC 1.5                      // the time step can be increased by this factor
#define MON_FUN_LOWER_LIMIT 1.e-3                       // the lower limit for the monitor function used to control the time step for integration over time (see Jannelli & Fazio 2006, Journal of Computational and Applied Mathematics, 191, 246)
#define MON_FUN_UPPER_LIMIT 1.e-2                       // the upper limit for the monitor function used to control the time step for integration over time
#define TIME_INTEGRATION_OF_KINETIC_EQUATIONS false     // if false the code will solve statistical equilibrium equations; if true the code will integrate kinetic equations for populations over time
#define USE_PARTITION_FUNCTIONS_RATIO true              // if true then the sum of populations will not be equal to 1.0 but will be equal to the ratio of full and partial partition functions (see Sobolev & Deguchi 1994); the full partition function is defined below
//
//#define PARTITION_FUNCTION 1.0
//#define PARTITION_FUNCTION (0.64 * pow(modelPhysPars::Tks, 1.5))          // CH3OH full partition rotational function (this is an approximation based on the data from CDMS database)
//#define PARTITION_FUNCTION (0.167 * pow(modelPhysPars::Tks, 1.92))        // 13CH3OH full partition rotational function
//#define PARTITION_FUNCTION (0.38 * pow(modelPhysPars::Tks, 1.76))         // H2CO full partition rotational function
//#define PARTITION_FUNCTION (0.03087 * pow(modelPhysPars::Tks, 1.48437))   // HDO full partition rotational function (this is an approximation based on the data from JPL database)
#define PARTITION_FUNCTION (1.72916 * pow(modelPhysPars::Tks, 1.47515))     // HNCO full partition rotational function

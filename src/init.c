/* Copyright 2010 by Roger S. Bivand. */

#include "spdep.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[]  = {
    {"dfs", (DL_FUNC) &dfs, 5},
    {"compute_gabriel", (DL_FUNC) &compute_gabriel, 7},
    {"compute_relative", (DL_FUNC) &compute_relative, 7},
    {"prunemst", (DL_FUNC) &prunemst, 4},
    {"gcdist", (DL_FUNC) &gcdist, 5},
    {"knearneigh", (DL_FUNC) &knearneigh, 7},
    {NULL, NULL, 0}
};

static R_CallMethodDef CallEntries[] = {
    {"opt_error_free", (DL_FUNC) &opt_error_free, 1},
    {"hess_error_free", (DL_FUNC) &hess_error_free, 1},
    {"hess_lag_free", (DL_FUNC) &hess_lag_free, 1},
    {"opt_error_init", (DL_FUNC) &opt_error_init, 0},
    {"hess_error_init", (DL_FUNC) &hess_error_init, 0},
    {"hess_lag_init", (DL_FUNC) &hess_lag_init, 0},
    {"R_ml_sse_env", (DL_FUNC) &R_ml_sse_env, 2},
    {"R_ml1_sse_env", (DL_FUNC) &R_ml1_sse_env, 3},
    {"R_ml2_sse_env", (DL_FUNC) &R_ml2_sse_env, 3},
    {"card", (DL_FUNC) &card, 1},
    {"listw2dsT", (DL_FUNC) &listw2dsT, 4},
    {"listw2dgR", (DL_FUNC) &listw2dgR, 4},
    {"listw2sn", (DL_FUNC) &listw2sn, 4},
    {"dnearneigh", (DL_FUNC) &dnearneigh, 6},
    {"gearyw", (DL_FUNC) &gearyw, 6},
    {"gsymtest", (DL_FUNC) &gsymtest, 3},
    {"spInsiders", (DL_FUNC) &spInsiders, 2},
    {"jcintern", (DL_FUNC) &jcintern, 4},
    {"lagw", (DL_FUNC) &lagw, 6},
    {"nbdists", (DL_FUNC) &nbdists, 5},
    {"polypoly", (DL_FUNC) &polypoly, 5},
    {"spOverlap", (DL_FUNC) &spOverlap, 2},
/*    {"poly_loop", (DL_FUNC) &poly_loop, 8},*/
    {"poly_loop2", (DL_FUNC) &poly_loop2, 8},
    {"symtest", (DL_FUNC) &symtest, 3},
    {"g_components", (DL_FUNC) &g_components, 2},
    {"mom_calc_int2", (DL_FUNC) &mom_calc_int2, 5},
    {"lmin21", (DL_FUNC) &lmin21, 4},
    {"lmin22", (DL_FUNC) &lmin22, 5},
    {"lmin23", (DL_FUNC) &lmin23, 6},
    {"lmin3", (DL_FUNC) &lmin3, 6},
    {"lmin3S", (DL_FUNC) &lmin3S, 7},
    {NULL, NULL, 0}
};


void 
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_spdep(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

}




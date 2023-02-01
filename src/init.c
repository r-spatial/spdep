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
    {"card", (DL_FUNC) &card, 1},
    {"listw2sn", (DL_FUNC) &listw2sn, 4},
    {"dnearneigh", (DL_FUNC) &dnearneigh, 6},
    {"dnearneigh1", (DL_FUNC) &dnearneigh1, 5},
    {"gearyw", (DL_FUNC) &gearyw, 6},
    {"gsymtest", (DL_FUNC) &gsymtest, 3},
    {"spInsiders", (DL_FUNC) &spInsiders, 2},
    {"jcintern", (DL_FUNC) &jcintern, 4},
    {"lagw", (DL_FUNC) &lagw, 6},
    {"nbdists", (DL_FUNC) &nbdists, 5},
    {"polypoly", (DL_FUNC) &polypoly, 5},
    {"spOverlap", (DL_FUNC) &spOverlap, 2},
    {"poly_loop2", (DL_FUNC) &poly_loop2, 8},
    {"symtest", (DL_FUNC) &symtest, 3},
    {"g_components", (DL_FUNC) &g_components, 2},
    {"perm_no_replace", (DL_FUNC) &perm_no_replace, 3},
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




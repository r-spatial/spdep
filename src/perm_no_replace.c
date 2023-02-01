/* Copyright 2023 by Roger S. Bivand. */

#include "spdep.h"

/* nsim, n, k */

SEXP draw_no_replace(int n, int crdi);

SEXP perm_no_replace(SEXP nsim0, SEXP n0, SEXP crdi0) {
    SEXP y, yk;
    int nsim = INTEGER_POINTER(nsim0)[0];
    int n = INTEGER_POINTER(n0)[0];
    int crdi = INTEGER_POINTER(crdi0)[0];
    GetRNGstate();
    PROTECT(y = allocVector(INTSXP, crdi*nsim));
    for (int k = 0; k < nsim; k++) {
        yk = draw_no_replace(n, crdi);
        for (int i = 0; i < crdi; i++) {
            INTEGER_POINTER(y)[k + (i*nsim)] = INTEGER_POINTER(yk)[i];
        }
    }
    PutRNGstate();
    UNPROTECT(1);
    return y;
}

SEXP draw_no_replace(int n, int crdi) {
    SEXP y;
    PROTECT(y = allocVector(INTSXP, crdi));
    int *iy = INTEGER(y);
    int *x = (int *)R_alloc(n, sizeof(int));
    for (int i = 0; i < n; i++) x[i] = i;
    for (int i = 0; i < crdi; i++) {
        int j = (int)(R_unif_index(n));
        iy[i] = x[j] + ROFFSET;
        x[j] = x[--n];
    }
    UNPROTECT(1);
    return y;
} 

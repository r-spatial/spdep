/* Copyright 2010 by Roger S. Bivand. */

#include "spdep.h"

static int c__1 = 1;

SEXP mom_calc_int2(SEXP is, SEXP m, SEXP nb, SEXP weights, SEXP card) {
    SEXP Omega;
    int hm = INTEGER_POINTER(m)[0];
    int n = length(card);
    double *eta, *zeta, *omega, sum, res;
    int i, ii, j, k1, k2, k3;
    int iis = length(is);

    omega = (double *) R_alloc((size_t) hm, sizeof(double));
    eta = (double *) R_alloc((size_t) n, sizeof(double));
    zeta = (double *) R_alloc((size_t) n, sizeof(double));
    for (j=0; j<hm; j++) omega[j] = 0.0;

    for (ii=0; ii<iis; ii++) {
        R_CheckUserInterrupt();
        i = INTEGER_POINTER(is)[ii]-ROFFSET;
        for (j=0; j<n; j++) eta[j] = 0.0;
        eta[i] = 1.0;
        for (j=1; j<hm; j=j+2) {
            for (k1=0; k1<n; k1++) {
                k3 = INTEGER_POINTER(card)[k1];
                if (k3 == 0) {
                    zeta[k1] = 0.0;
                } else {
                    sum = 0.0;
                    for (k2=0; k2<k3; k2++) {
                        sum += eta[INTEGER_POINTER(VECTOR_ELT(nb, k1))[k2]
                            - ROFFSET] * NUMERIC_POINTER(VECTOR_ELT(weights,
                            k1))[k2];
                    }
                    zeta[k1] = sum;
                }
            }
            res = F77_CALL(ddot)(&n, zeta, &c__1, eta, &c__1);
            if (R_FINITE(res)) omega[(j-1)] += res;
            else error("non-finite dot product %d, %d", i, j);
            res = F77_CALL(ddot)(&n, zeta, &c__1, zeta, &c__1);
            if (R_FINITE(res)) omega[j] += res;
            else error("non-finite dot product %d, %d", i, j);
            for (k1=0; k1<n; k1++) eta[k1] = zeta[k1];
        }
    }

    PROTECT(Omega = NEW_NUMERIC(hm));
    for (j=0; j<hm; j++) NUMERIC_POINTER(Omega)[j] = omega[j];

    UNPROTECT(1);
    return(Omega);
}

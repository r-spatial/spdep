/* Copyright 2015 by Roger S. Bivand. */

#include "spdep.h"


SEXP lmin21(SEXP nb, SEXP y, SEXP cy, SEXP card) {
    int i, j, k, nswitch=0, n=length(card), pc=0;
    SEXP ans;
    double t1, t2, ytemp;
    double *Y, *CY;

    Y = (double *) R_alloc((size_t) n, sizeof(double));
    CY = (double *) R_alloc((size_t) n, sizeof(double));

    for (i=0; i<n; i++) {
        Y[i] = NUMERIC_POINTER(y)[i];
        CY[i] = NUMERIC_POINTER(cy)[i];
    }

    PROTECT(ans = NEW_LIST(2)); pc++;
    SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
    SET_VECTOR_ELT(ans, 1, NEW_INTEGER(1));

    for (i=0; i<n; i++) {
      if (INTEGER_POINTER(card)[i] > 0) {
        t1 = fabs(Y[i] - CY[i]);
        t2 = fabs(-2*CY[i]);
        for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
            k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j]-ROFFSET;
            t1 = t1 + fabs(Y[k] - CY[k]);
            t2 = t2 + fabs(Y[k] - (CY[k] - Y[i] - CY[i]));
        }
        if (t1 <= t2) {
            nswitch++;
            ytemp = Y[i];
            Y[i] = -CY[i];
            for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
                k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j]-ROFFSET;
                CY[k] = CY[k] - ytemp + Y[i];
            }
        }
      }
    }

    for (i=0; i<n; i++) {
        NUMERIC_POINTER(VECTOR_ELT(ans, 0))[i] = Y[i];
    }

    INTEGER_POINTER(VECTOR_ELT(ans, 1))[0] = nswitch;
    UNPROTECT(pc); /* ans */
    return(ans);
}

SEXP lmin22(SEXP nb, SEXP y, SEXP cy, SEXP card, SEXP beta) {
    int i, j, k, nswitch=0, n=length(card), pc=0;
    SEXP ans;
    double t1, t2, ytemp, yhat;
    double *Y, *CY, *B;

    Y = (double *) R_alloc((size_t) n, sizeof(double));
    CY = (double *) R_alloc((size_t) n, sizeof(double));
    B = (double *) R_alloc((size_t) length(beta), sizeof(double));

    for (i=0; i<n; i++) {
        Y[i] = NUMERIC_POINTER(y)[i];
        CY[i] = NUMERIC_POINTER(cy)[i];
    }
    for (i=0; i<length(beta); i++) {
        B[i] = NUMERIC_POINTER(beta)[i];
    }

    PROTECT(ans = NEW_LIST(2)); pc++;
    SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
    SET_VECTOR_ELT(ans, 1, NEW_INTEGER(1));

    for (i=0; i<n; i++) {
      if (INTEGER_POINTER(card)[i] > 0) {
        t1 = fabs(Y[i] - CY[i]);
        yhat = B[0] + B[1]*CY[i];
        t2 = fabs(yhat - CY[i]);
        for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
            k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j]-ROFFSET;
            t1 = t1 + fabs(Y[k] - CY[k]);
            t2 = t2 + fabs(Y[k] - (CY[k] - Y[i] + B[0] + B[1]*CY[i]));
        }
        if (t1 <= t2) {
            nswitch++;
            ytemp = Y[i];
            Y[i] = yhat;
            for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
                k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j]-ROFFSET;
                CY[k] = CY[k] - ytemp + Y[i];
            }
        }
      }
    }

    for (i=0; i<n; i++) {
        NUMERIC_POINTER(VECTOR_ELT(ans, 0))[i] = Y[i];
    }

    INTEGER_POINTER(VECTOR_ELT(ans, 1))[0] = nswitch;
    UNPROTECT(pc); /* ans */
    return(ans);
}

SEXP lmin23(SEXP nb, SEXP y, SEXP cy, SEXP card, SEXP beta, SEXP tol) {
    int i, j, k, nswitch=0, n=length(card), pc=0;
    SEXP ans;
    double tmp, var, yhat;
    double *Y, *CY, *B;

    Y = (double *) R_alloc((size_t) n, sizeof(double));
    CY = (double *) R_alloc((size_t) n, sizeof(double));
    B = (double *) R_alloc((size_t) length(beta), sizeof(double));

    for (i=0; i<n; i++) {
        Y[i] = NUMERIC_POINTER(y)[i];
        CY[i] = NUMERIC_POINTER(cy)[i];
    }
    for (i=0; i<length(beta); i++) {
        B[i] = NUMERIC_POINTER(beta)[i];
    }
    PROTECT(ans = NEW_LIST(2)); pc++;
    SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
    SET_VECTOR_ELT(ans, 1, NEW_INTEGER(1));

    for (i=0; i<n; i++) {
      if (INTEGER_POINTER(card)[i] > 0) {
        yhat = B[0] + B[1]*CY[i];
        var = fabs(Y[i] - yhat);
        if (var > NUMERIC_POINTER(tol)[0]) {
            nswitch++;
            tmp = Y[i];
            Y[i] = yhat;
            for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
                k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j]-ROFFSET;
                CY[k] = CY[k] - tmp + Y[i];
            }
        }
      }
    }

    for (i=0; i<n; i++) {
        NUMERIC_POINTER(VECTOR_ELT(ans, 0))[i] = Y[i];
    }

    INTEGER_POINTER(VECTOR_ELT(ans, 1))[0] = nswitch;
    UNPROTECT(pc); /* ans */
    return(ans);
}

SEXP lmin3(SEXP nb, SEXP ev1, SEXP ev1_lag, SEXP n_nei, SEXP beta, SEXP tol) {
    int i, j, k, nswitch=0, n=length(n_nei), pc=0;
    SEXP ans;
    double tmp, var, yhat, ntmp;
    double *Y, *CY, *B;

    Y = (double *) R_alloc((size_t) n, sizeof(double));
    CY = (double *) R_alloc((size_t) n, sizeof(double));
    B = (double *) R_alloc((size_t) length(beta), sizeof(double));

    for (i=0; i<n; i++) {
        Y[i] = NUMERIC_POINTER(ev1)[i];
        CY[i] = NUMERIC_POINTER(ev1_lag)[i];
    }
    for (i=0; i<length(beta); i++) {
        B[i] = NUMERIC_POINTER(beta)[i];
    }
    PROTECT(ans = NEW_LIST(2)); pc++;
    SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
    SET_VECTOR_ELT(ans, 1, NEW_INTEGER(1));

    for (i=0; i<n; i++) {
      if (INTEGER_POINTER(n_nei)[i] > 0) {
        yhat = B[0] + B[1]*CY[i];
        var = fabs(Y[i] - yhat);
        if (var > NUMERIC_POINTER(tol)[0]) {
            nswitch++;
            tmp = Y[i];
            Y[i] = yhat;
            for (j=0; j<INTEGER_POINTER(n_nei)[i]; j++) {
                k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j]-ROFFSET;
                ntmp = sqrt(INTEGER_POINTER(n_nei)[i] *
                    INTEGER_POINTER(n_nei)[k]);
                CY[k] = CY[k] - (tmp/ntmp) + (Y[i]/ntmp);
            }
        }
      }
    }

    for (i=0; i<n; i++) {
        NUMERIC_POINTER(VECTOR_ELT(ans, 0))[i] = Y[i];
    }

    INTEGER_POINTER(VECTOR_ELT(ans, 1))[0] = nswitch;
    UNPROTECT(pc); /* ans */
    return(ans);
}


SEXP lmin3S(SEXP nb, SEXP ev1, SEXP ev1_lag, SEXP n_nei, SEXP card, SEXP beta, SEXP tol) {
    int i, j, k, nswitch=0, n=length(card), pc=0;
    SEXP ans;
    double tmp, var, yhat, ntmp;
    double *Y, *CY, *B;

    Y = (double *) R_alloc((size_t) n, sizeof(double));
    CY = (double *) R_alloc((size_t) n, sizeof(double));
    B = (double *) R_alloc((size_t) length(beta), sizeof(double));

    for (i=0; i<n; i++) {
        Y[i] = NUMERIC_POINTER(ev1)[i];
        CY[i] = NUMERIC_POINTER(ev1_lag)[i];
    }
    for (i=0; i<length(beta); i++) {
        B[i] = NUMERIC_POINTER(beta)[i];
    }
    PROTECT(ans = NEW_LIST(2)); pc++;
    SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
    SET_VECTOR_ELT(ans, 1, NEW_INTEGER(1));

    for (i=0; i<n; i++) {
      if (INTEGER_POINTER(card)[i] > 0) {
        yhat = B[0] + B[1]*CY[i];
        var = fabs(Y[i] - yhat);
        if (var > NUMERIC_POINTER(tol)[0]) {
            nswitch++;
            tmp = Y[i];
            Y[i] = yhat;
            for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
                k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j]-ROFFSET;
                ntmp = sqrt(NUMERIC_POINTER(n_nei)[i] *
                    NUMERIC_POINTER(n_nei)[k]);
                CY[k] = CY[k] - (tmp/ntmp) + (Y[i]/ntmp);
            }
        }
      }
    }

    for (i=0; i<n; i++) {
        NUMERIC_POINTER(VECTOR_ELT(ans, 0))[i] = Y[i];
    }

    INTEGER_POINTER(VECTOR_ELT(ans, 1))[0] = nswitch;
    UNPROTECT(pc); /* ans */
    return(ans);
}

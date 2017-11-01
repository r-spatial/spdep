/* Copyright 2010 by Roger S. Bivand. */

#include "spdep.h"

/** */
static int c__1 = 1;

typedef struct opt_error_sse {
    double *y;
    double *x;
    double *yl;
    double *wy1;
    double *xlq;
    double *wx1;
    double *qy;
    double *xlqyl;
    double *work;
    double *qraux;
    int *jpvt;
    int set;
} OPT_ERROR_SSE;

typedef struct hess_error_sse {
    double *y;
    double *x;
    double *yl;
    double *wy1;
    double *xl;
    double *wx1;
    double *beta1;
    double *xlb;
    int set;
} HESS_ERROR_SSE;

typedef struct hess_lag_sse {
    double *y;
    double *x;
    double *yl;
    double *wy1;
    double *beta1;
    double *xb;
    int set;
} HESS_LAG_SSE;

SEXP opt_error_init() {

    OPT_ERROR_SSE *pt;
    SEXP ptr;

    pt = Calloc(1, OPT_ERROR_SSE);
    pt->set = FALSE;

    PROTECT(ptr = R_MakeExternalPtr(pt, R_NilValue, R_NilValue));

    UNPROTECT(1);
    return(ptr);

}


void opt_error_set(SEXP env) {

    OPT_ERROR_SSE *pt;
    SEXP y, x, wy, WX;
    int i, n, p, np, pc=0;

    n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
    p = INTEGER_POINTER(findVarInFrame(env, install("p")))[0];
    np = n*p;

    pt = (OPT_ERROR_SSE *) R_ExternalPtrAddr(findVarInFrame(env,
        install("ptr")));
    if (pt->set) error("opt_error_set: function called out of order");

    PROTECT(y = findVarInFrame(env, install("y"))); pc++;
    PROTECT(x = findVarInFrame(env, install("x"))); pc++;
    PROTECT(wy = findVarInFrame(env, install("wy"))); pc++;
    PROTECT(WX = findVarInFrame(env, install("WX"))); pc++;

    pt->y = Calloc(n, double);
    pt->x = Calloc(np, double);
    pt->yl = Calloc(n, double);
    pt->wy1 = Calloc(n, double);
    pt->xlq = Calloc(np, double);
    pt->wx1 = Calloc(np, double);
    pt->qy = Calloc(np, double);
    pt->xlqyl = Calloc(p, double);
    pt->jpvt = Calloc(p, int);
    pt->work = Calloc(p*2, double); 
/*    pt->work = Calloc(p, double); */
    pt->qraux = Calloc(p, double);

    for (i=0; i<n; i++) {
        pt->y[i] = NUMERIC_POINTER(y)[i];
        pt->wy1[i] = NUMERIC_POINTER(wy)[i];
    }
    for (i=0; i<np; i++) {
        pt->x[i] = NUMERIC_POINTER(x)[i];
        pt->wx1[i] = NUMERIC_POINTER(WX)[i];
    }
    pt->set = TRUE;
    UNPROTECT(pc);

    return;
}


SEXP opt_error_free(SEXP ptr) {

    OPT_ERROR_SSE *pt;

    pt = (OPT_ERROR_SSE *) R_ExternalPtrAddr(ptr);

    Free(pt->qraux);
    Free(pt->work);
    Free(pt->jpvt);
    Free(pt->xlqyl);
    Free(pt->qy);
    Free(pt->wx1);
    Free(pt->xlq);
    Free(pt->wy1);
    Free(pt->yl);
    Free(pt->x);
    Free(pt->y);

    Free(pt);
    R_ClearExternalPtr(ptr);
    return(R_NilValue);
}

SEXP hess_error_init() {

    HESS_ERROR_SSE *pt;
    SEXP ptr;

    pt = Calloc(1, HESS_ERROR_SSE);
    pt->set = FALSE;

    PROTECT(ptr = R_MakeExternalPtr(pt, R_NilValue, R_NilValue));

    UNPROTECT(1);
    return(ptr);

}


void hess_error_set(SEXP env) {

    HESS_ERROR_SSE *pt;
    SEXP y, x, wy, WX;
    int i, n, p, np, pc=0;

    n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
    p = INTEGER_POINTER(findVarInFrame(env, install("p")))[0];
    np = n*p;

    pt = (HESS_ERROR_SSE *) R_ExternalPtrAddr(findVarInFrame(env,
        install("ptr")));
    if (pt->set) error("hess_error_set: function called out of order");

    PROTECT(y = findVarInFrame(env, install("y"))); pc++;
    PROTECT(x = findVarInFrame(env, install("x"))); pc++;
    PROTECT(wy = findVarInFrame(env, install("wy"))); pc++;
    PROTECT(WX = findVarInFrame(env, install("WX"))); pc++;

    pt->y = Calloc(n, double);
    pt->x = Calloc(np, double);
    pt->yl = Calloc(n, double);
    pt->wy1 = Calloc(n, double);
    pt->xl = Calloc(np, double);
    pt->wx1 = Calloc(np, double);
    pt->beta1 = Calloc(p, double);
    pt->xlb = Calloc(n, double);

    for (i=0; i<n; i++) {
        pt->y[i] = NUMERIC_POINTER(y)[i];
        pt->wy1[i] = NUMERIC_POINTER(wy)[i];
    }
    for (i=0; i<np; i++) {
        pt->x[i] = NUMERIC_POINTER(x)[i];
        pt->wx1[i] = NUMERIC_POINTER(WX)[i];
    }
    pt->set = TRUE;
    UNPROTECT(pc);

    return;
}


SEXP hess_error_free(SEXP ptr) {

    HESS_ERROR_SSE *pt;

    pt = (HESS_ERROR_SSE *) R_ExternalPtrAddr(ptr);

    Free(pt->xlb);
    Free(pt->beta1);
    Free(pt->wx1);
    Free(pt->xl);
    Free(pt->wy1);
    Free(pt->yl);
    Free(pt->x);
    Free(pt->y);

    Free(pt);
    R_ClearExternalPtr(ptr);
    return(R_NilValue);
}

SEXP hess_lag_init() {

    HESS_LAG_SSE *pt;
    SEXP ptr;

    pt = Calloc(1, HESS_LAG_SSE);
    pt->set = FALSE;

    PROTECT(ptr = R_MakeExternalPtr(pt, R_NilValue, R_NilValue));

    UNPROTECT(1);
    return(ptr);

}


void hess_lag_set(SEXP env) {

    HESS_LAG_SSE *pt;
    SEXP y, x, wy;
    int i, n, p, np, pc=0;

    n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
    p = INTEGER_POINTER(findVarInFrame(env, install("m")))[0];
    np = n*p;

    pt = (HESS_LAG_SSE *) R_ExternalPtrAddr(findVarInFrame(env,
        install("ptr")));
    if (pt->set) error("hess_lag_set: function called out of order");

    PROTECT(y = findVarInFrame(env, install("y"))); pc++;
    PROTECT(x = findVarInFrame(env, install("x"))); pc++;
    PROTECT(wy = findVarInFrame(env, install("wy"))); pc++;

    pt->y = Calloc(n, double);
    pt->x = Calloc(np, double);
    pt->yl = Calloc(n, double);
    pt->wy1 = Calloc(n, double);
    pt->beta1 = Calloc(p, double);
    pt->xb = Calloc(n, double);

    for (i=0; i<n; i++) {
        pt->y[i] = NUMERIC_POINTER(y)[i];
        pt->wy1[i] = NUMERIC_POINTER(wy)[i];
    }
    for (i=0; i<np; i++) pt->x[i] = NUMERIC_POINTER(x)[i];
    pt->set = TRUE;
    UNPROTECT(pc);

    return;
}


SEXP hess_lag_free(SEXP ptr) {

    HESS_LAG_SSE *pt;

    pt = (HESS_LAG_SSE *) R_ExternalPtrAddr(ptr);

    Free(pt->xb);
    Free(pt->beta1);
    Free(pt->wy1);
    Free(pt->yl);
    Free(pt->x);
    Free(pt->y);

    Free(pt);
    R_ClearExternalPtr(ptr);
    return(R_NilValue);
}

/**
 * Calculate the sum of squared errors term for spatial regression
 * using an environment to hold data
 *
 * @param env pointer to an SEXP environment
 * @param coef current value of coefficient being optimzed
 * 
 * @return double, value of SSE for current coef
 *
 */
SEXP R_ml_sse_env(SEXP env, SEXP coef) {

  SEXP res;
//  SEXP y, x, wy, WX;
  int i, k, n, p, np;
  double tol=1e-7, cyl, cxlqyl, sse;
  char *trans = "T";
  double one = 1.0, zero = 0.0;
  double m_lambda = - NUMERIC_POINTER(coef)[0];
  int pc=0, first_time;
  OPT_ERROR_SSE *pt;

  first_time = LOGICAL_POINTER(findVarInFrame(env, install("first_time")))[0];
  if (first_time) {
    opt_error_set(env);
  }

  n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
  p = INTEGER_POINTER(findVarInFrame(env, install("p")))[0];
  np = n*p;
  pt = (OPT_ERROR_SSE *) R_ExternalPtrAddr(findVarInFrame(env,
        install("ptr")));

  for (i=0; i<n; i++) pt->yl[i] = pt->y[i];
  for (i=0; i<np; i++) pt->xlq[i] = pt->x[i];

  F77_CALL(daxpy)(&n, &m_lambda, pt->wy1, &c__1, pt->yl, &c__1);

  F77_CALL(daxpy)(&np, &m_lambda, pt->wx1, &c__1, pt->xlq, &c__1);

  F77_CALL(dqrdc2)(pt->xlq, &n, &n, &p, &tol, &k, pt->qraux, pt->jpvt,
    pt->work); 
  if (p != k) warning("Q looses full rank"); 
/*  k = 0;
  F77_CALL(dqrdc)(pt->xlq, &n, &n, &p, pt->qraux, pt->jpvt, pt->work, &k);*/

  for (i=0; i<n*k; i++) pt->qy[i] = 0.0;
  for (i=0; i<k; i++) pt->qy[(i +(n*i))] = 1.0;

  F77_CALL(dqrqy)(pt->xlq, &n, &k, pt->qraux, pt->qy, &k, pt->qy);

  F77_CALL(dgemv)(trans, &n, &k, &one, pt->qy, &n, pt->yl, &c__1, &zero,
    pt->xlqyl, &c__1);

  cyl = F77_CALL(ddot)(&n, pt->yl, &c__1, pt->yl, &c__1);

  cxlqyl = F77_CALL(ddot)(&k, pt->xlqyl, &c__1, pt->xlqyl, &c__1);

  sse = cyl - cxlqyl;

  PROTECT(res=NEW_NUMERIC(1)); pc++;
  NUMERIC_POINTER(res)[0] = sse;
  UNPROTECT(pc);

  return(res);

}

SEXP R_ml1_sse_env(SEXP env, SEXP lambda, SEXP beta) {

  SEXP res;
  int i, n, p, np;
  double sse;
  char *trans = "N";
  double one = 1.0, zero = 0.0, m_one = -1.0;
  double m_lambda = - NUMERIC_POINTER(lambda)[0];
  int pc=0, first_time;
  HESS_ERROR_SSE *pt;

  first_time = LOGICAL_POINTER(findVarInFrame(env, install("first_time")))[0];
  if (first_time) {
    hess_error_set(env);
  }

  n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
  p = INTEGER_POINTER(findVarInFrame(env, install("p")))[0];
  np = n*p;
  pt = (HESS_ERROR_SSE *) R_ExternalPtrAddr(findVarInFrame(env,
        install("ptr")));

  for (i=0; i<n; i++) pt->yl[i] = pt->y[i];
  for (i=0; i<np; i++) pt->xl[i] = pt->x[i];

  for (i=0; i<p; i++) pt->beta1[i] = NUMERIC_POINTER(beta)[i];

  F77_CALL(daxpy)(&n, &m_lambda, pt->wy1, &c__1, pt->yl, &c__1);

  F77_CALL(daxpy)(&np, &m_lambda, pt->wx1, &c__1, pt->xl, &c__1);

  F77_CALL(dgemv)(trans, &n, &p, &one, pt->xl, &n, pt->beta1, &c__1, &zero,
    pt->xlb, &c__1);

  F77_CALL(daxpy)(&n, &m_one, pt->xlb, &c__1, pt->yl, &c__1);

  sse = F77_CALL(ddot)(&n, pt->yl, &c__1, pt->yl, &c__1);

  PROTECT(res=NEW_NUMERIC(1)); pc++;
  NUMERIC_POINTER(res)[0] = sse;
  UNPROTECT(pc);

  return(res);

}

SEXP R_ml2_sse_env(SEXP env, SEXP rho, SEXP beta) {

  SEXP res;
  int i, n, p;
  double sse;
  char *trans = "N";
  double one = 1.0, zero = 0.0, m_one = -1.0;
  double m_rho1 = - NUMERIC_POINTER(rho)[0];
  int pc=0, first_time;
  HESS_LAG_SSE *pt;

  first_time = LOGICAL_POINTER(findVarInFrame(env, install("first_time")))[0];
  if (first_time) {
    hess_lag_set(env);
  }

  n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
  p = INTEGER_POINTER(findVarInFrame(env, install("m")))[0];
  pt = (HESS_LAG_SSE *) R_ExternalPtrAddr(findVarInFrame(env,
        install("ptr")));

  for (i=0; i<n; i++) pt->yl[i] = pt->y[i];
  for (i=0; i<p; i++) pt->beta1[i] = NUMERIC_POINTER(beta)[i];

  F77_CALL(daxpy)(&n, &m_rho1, pt->wy1, &c__1, pt->yl, &c__1);

  F77_CALL(dgemv)(trans, &n, &p, &one, pt->x, &n, pt->beta1, &c__1, &zero,
    pt->xb, &c__1);

  F77_CALL(daxpy)(&n, &m_one, pt->xb, &c__1, pt->yl, &c__1);

  sse = F77_CALL(ddot)(&n, pt->yl, &c__1, pt->yl, &c__1);

  PROTECT(res=NEW_NUMERIC(1)); pc++;
  NUMERIC_POINTER(res)[0] = sse;
  UNPROTECT(pc);

  return(res);

}


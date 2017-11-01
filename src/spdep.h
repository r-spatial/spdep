/* Copyright 2010 by Roger S. Bivand. */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Applic.h>
/* #include <R_ext/Linpack.h> */
#include <R_ext/Utils.h>
#define ROFFSET 1

SEXP opt_error_free(SEXP ptr);
SEXP hess_error_free(SEXP ptr);
SEXP hess_lag_free(SEXP ptr);
SEXP opt_error_init();
SEXP hess_error_init();
SEXP hess_lag_init();
SEXP R_ml_sse_env(SEXP env, SEXP coef);
SEXP R_ml1_sse_env(SEXP env, SEXP lambda, SEXP beta);
SEXP R_ml2_sse_env(SEXP env, SEXP rho, SEXP beta);
SEXP mom_calc_int2(SEXP is, SEXP m, SEXP nb, SEXP weights, SEXP card);

void opt_error_set(SEXP env);
void hess_error_set(SEXP env);
void hess_lag_set(SEXP env);

SEXP card(SEXP nb);
SEXP listw2dsT(SEXP nbs, SEXP wts, SEXP card, SEXP ncard2);
SEXP listw2dgR(SEXP nbs, SEXP wts, SEXP card, SEXP ncard);
SEXP listw2sn(SEXP nbs, SEXP wts, SEXP card, SEXP ncard);
SEXP dnearneigh(SEXP din1, SEXP din2, SEXP pnte, SEXP p, SEXP test, SEXP lonlat);
SEXP gearyw(SEXP nb, SEXP weights, SEXP x, SEXP card, SEXP zeropolicy, SEXP ftype);
SEXP gsymtest(SEXP nb, SEXP glist, SEXP card);
SEXP spInsiders(SEXP bbbi, SEXP bbbj);
SEXP jcintern(SEXP nb, SEXP weights, SEXP dum, SEXP card);
SEXP lagw(SEXP nb, SEXP weights, SEXP x, SEXP card, SEXP zeropolicy, SEXP naok);
SEXP nbdists(SEXP nb, SEXP x, SEXP np, SEXP dim, SEXP lonlat);
SEXP polypoly(SEXP p1, SEXP n01, SEXP p2, SEXP n02, SEXP snap);
SEXP spOverlap(SEXP bbbi, SEXP bbbj);
/*SEXP poly_loop(SEXP n, SEXP i_findInBox, SEXP bb, SEXP pl, SEXP nrs, SEXP dsnap, SEXP criterion, SEXP scale);*/
SEXP poly_loop2(SEXP n, SEXP i_findInBox, SEXP bb, SEXP pl, SEXP nrs, SEXP dsnap, SEXP criterion, SEXP scale);
SEXP symtest(SEXP nb, SEXP card, SEXP verbose);
SEXP g_components(SEXP nblst, SEXP cmpnm);
SEXP lmin21(SEXP nb, SEXP y, SEXP cy, SEXP card);
SEXP lmin22(SEXP nb, SEXP y, SEXP cy, SEXP card, SEXP beta);
SEXP lmin23(SEXP nb, SEXP y, SEXP cy, SEXP card, SEXP beta, SEXP tol);
SEXP lmin3(SEXP nb, SEXP ev1, SEXP ev1_lag, SEXP n_nei, SEXP beta, SEXP tol);
SEXP lmin3S(SEXP nb, SEXP ev1, SEXP ev1_lag, SEXP n_nei, SEXP card, SEXP beta, SEXP tol);

void dfs(SEXP nblst, SEXP cmpnm, SEXP visited, int curcmp, int nodeid);
void compute_gabriel(int *no_nodes, int *g1, int *g2, int *nogab, int *ngaballoc,  double *nodes_xd, double *nodes_yd);
void compute_relative(int *no_nodes, int *g1, int *g2, int *nogab, int *ngaballoc, double *nodes_xd, double *nodes_yd);
void prunemst(int *e1, int *e2, int *ne, int *gr);

void gcdist(double *lon1, double *lon2, double *lat1, double *lat2, double *dist);
void knearneigh(int *kin, int *pnte, int *p, double *test, int *res, double *dists, int *lonlat);








/* Copyright 2010-2021 by Roger S. Bivand. */

#define USE_FC_LEN_T
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif
#include <R_ext/Applic.h>
/* #include <R_ext/Linpack.h> */
#include <R_ext/Utils.h>
#define ROFFSET 1

SEXP card(SEXP nb);
SEXP listw2sn(SEXP nbs, SEXP wts, SEXP card, SEXP ncard);
SEXP dnearneigh(SEXP din1, SEXP din2, SEXP pnte, SEXP p, SEXP test, SEXP lonlat);
SEXP dnearneigh1(SEXP din1, SEXP din2, SEXP pnte, SEXP test, SEXP cands);

SEXP gearyw(SEXP nb, SEXP weights, SEXP x, SEXP card, SEXP zeropolicy, SEXP ftype);
SEXP gsymtest(SEXP nb, SEXP glist, SEXP card);
SEXP spInsiders(SEXP bbbi, SEXP bbbj);
SEXP jcintern(SEXP nb, SEXP weights, SEXP dum, SEXP card);
SEXP lagw(SEXP nb, SEXP weights, SEXP x, SEXP card, SEXP zeropolicy, SEXP naok);
SEXP nbdists(SEXP nb, SEXP x, SEXP np, SEXP dim, SEXP lonlat);
SEXP polypoly(SEXP p1, SEXP n01, SEXP p2, SEXP n02, SEXP snap);
SEXP spOverlap(SEXP bbbi, SEXP bbbj);
SEXP poly_loop2(SEXP n, SEXP i_findInBox, SEXP bb, SEXP pl, SEXP nrs, SEXP dsnap, SEXP criterion, SEXP scale);
SEXP symtest(SEXP nb, SEXP card, SEXP verbose);
SEXP g_components(SEXP nblst, SEXP cmpnm);
SEXP perm_no_replace(SEXP nsim0, SEXP n0, SEXP crsi0);

void dfs(SEXP nblst, SEXP cmpnm, SEXP visited, int curcmp, int nodeid);
void compute_gabriel(int *no_nodes, int *g1, int *g2, int *nogab, int *ngaballoc,  double *nodes_xd, double *nodes_yd);
void compute_relative(int *no_nodes, int *g1, int *g2, int *nogab, int *ngaballoc, double *nodes_xd, double *nodes_yd);
void prunemst(int *e1, int *e2, int *ne, int *gr);

void gcdist(double *lon1, double *lon2, double *lat1, double *lat2, double *dist);
void knearneigh(int *kin, int *pnte, int *p, double *test, int *res, double *dists, int *lonlat);








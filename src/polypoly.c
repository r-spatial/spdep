/* Copyright 2004 by Roger S. Bivand. */

#include "spdep.h"

SEXP polypoly(SEXP p1, SEXP n01, SEXP p2, SEXP n02, SEXP snap)
{
	int n1=INTEGER_POINTER(n01)[0], n2=INTEGER_POINTER(n02)[0], pc=0;
	int i, j, k=0;
	double sn=NUMERIC_POINTER(snap)[0], dist;
	double x1, x2, y1, y2, xd, yd;

	SEXP ans;
	PROTECT(ans = NEW_INTEGER(1)); pc++;

	for (i=0; (i < n1) && (k < 2); i++) {
		x1 = NUMERIC_POINTER(p1)[i];
		y1 = NUMERIC_POINTER(p1)[n1 + i];
		for (j=0; (j < n2) && (k < 2); j++) {
			x2 = NUMERIC_POINTER(p2)[j];
			y2 = NUMERIC_POINTER(p2)[n2 + j];
/*			dist = pythag((x1-x2), (y1-y2));
			if (dist < sn) k++;
			if (k > 1) break;*/
/* following lines Micah Altman 2010 */
			xd = x1-x2;
			if (fabs(xd)>sn) { continue; }
			yd = y1-y2;
			if (fabs(yd)>sn) { continue; }
			dist = hypot(xd, yd);
			if (dist <= sn) {
                            k++;
                        }
		}
	}
	
	INTEGER_POINTER(ans)[0] = k;

	UNPROTECT(pc); /* ans */
	return(ans);
}

/* function by Micah Altman */

SEXP spOverlap(SEXP bbbi, SEXP bbbj) {

	int pc=0,overlap=1;
	double bbi[4], bbj[4];
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(1)); pc++;
	bbi[0] = NUMERIC_POINTER(bbbi)[0];
	bbi[1] = NUMERIC_POINTER(bbbi)[1];
	bbi[2] = NUMERIC_POINTER(bbbi)[2];
	bbi[3] = NUMERIC_POINTER(bbbi)[3];
	bbj[0] = NUMERIC_POINTER(bbbj)[0];
	bbj[1] = NUMERIC_POINTER(bbbj)[1];
	bbj[2] = NUMERIC_POINTER(bbbj)[2];
	bbj[3] = NUMERIC_POINTER(bbbj)[3];

        if ((bbi[0]>bbj[2]) | (bbi[1]>bbj[3]) | 
		(bbi[2]<bbj[0]) | (bbi[3]<bbj[1]) ) {
		overlap=0;
	}

	INTEGER_POINTER(ans)[0] = overlap;		
	UNPROTECT(pc); /* ans */
	return(ans);
}

/* SEXP poly_loop(SEXP n, SEXP i_findInBox, SEXP bb, SEXP pl, SEXP nrs,
    SEXP dsnap, SEXP criterion, SEXP scale) {

    int nn = INTEGER_POINTER(n)[0];
    int crit = INTEGER_POINTER(criterion)[0];
    int Scale = INTEGER_POINTER(scale)[0];
    int uBound = nn*Scale;
    int i, j, jj, k, li, pc = 0;
    int ii = 0;
    int *card, *icard, *is, *jjs;

    SEXP bbi, bbj, jhit, khit, ans, pli, plj, nrsi, nrsj;

    int xx, yy, zz, ww;

    card = (int *) R_alloc((size_t) nn, sizeof(int));
    icard = (int *) R_alloc((size_t) nn, sizeof(int));
    is = (int *) R_alloc((size_t) uBound, sizeof(int));
    jjs = (int *) R_alloc((size_t) uBound, sizeof(int));

    for (i=0; i<nn; i++) {
        card[i] = 0;
        icard[i] = 0;
   
    }
    for (i=0; i<uBound; i++) {
        is[i] = 0;
        jjs[i] = 0;
    }

    PROTECT(bbi = NEW_NUMERIC(4)); pc++;
    PROTECT(bbj = NEW_NUMERIC(4)); pc++;
    PROTECT(jhit = NEW_INTEGER(1)); pc++;
    PROTECT(khit = NEW_INTEGER(1)); pc++;
    PROTECT(nrsi = NEW_INTEGER(1)); pc++;
    PROTECT(nrsj = NEW_INTEGER(1)); pc++;

    for (i=0; i<(nn-1); i++) {
        li = length(VECTOR_ELT(i_findInBox, i));
        INTEGER_POINTER(nrsi)[0] = INTEGER_POINTER(nrs)[i];
        for (k=0; k<4; k++) 
            NUMERIC_POINTER(bbi)[k] = NUMERIC_POINTER(bb)[i+(k*nn)];
        for (j=0; j<li; j++) {
            jj = INTEGER_POINTER(VECTOR_ELT(i_findInBox, i))[j] - ROFFSET;
            for (k=0; k<4; k++) 
                NUMERIC_POINTER(bbj)[k] = NUMERIC_POINTER(bb)[jj+(k*nn)];
            jhit = spOverlap(bbi, bbj);
            if (INTEGER_POINTER(jhit)[0] > 0) {
                INTEGER_POINTER(khit)[0] = 0;
                INTEGER_POINTER(nrsj)[0] = INTEGER_POINTER(nrs)[jj];
                if (INTEGER_POINTER(nrsi)[0]*INTEGER_POINTER(nrsj)[0] > 0){
                    khit = polypoly(VECTOR_ELT(pl, i), nrsi, VECTOR_ELT(pl, jj),
                        nrsj, dsnap);
                }
                if (INTEGER_POINTER(khit)[0] > crit) {
                    card[i]++;
                    card[jj]++;
                    is[ii] = i;
                    jjs[ii] = jj;
                    ii++;
                    if (ii == uBound) error("memory error, scale problem");
                }
            }
        }
    }

    PROTECT(ans = NEW_LIST(nn)); pc++;

    for (i=0; i<nn; i++) {
        if (card[i] == 0) {
            SET_VECTOR_ELT(ans, i, NEW_INTEGER(1));
            INTEGER_POINTER(VECTOR_ELT(ans, i))[0] = 0;
        } else {
            SET_VECTOR_ELT(ans, i, NEW_INTEGER(card[i]));
        }
    }

    for (i=0; i<ii; i++) {
        xx = is[i];
        yy = jjs[i];
        zz = icard[yy];
        ww = icard[xx];
        if (zz == card[yy]) error("memory error, overflow");
        if (ww == card[xx]) error("memory error, overflow");
        INTEGER_POINTER(VECTOR_ELT(ans, yy))[zz] = xx + ROFFSET;
        INTEGER_POINTER(VECTOR_ELT(ans, xx))[ww] = yy + ROFFSET;
        icard[yy]++;
        icard[xx]++;
    }

    for (i=0; i<nn; i++) {
        if ((li = length(VECTOR_ELT(ans, i))) > 1) {
            for (j=0; j<li; j++)
                icard[j] = INTEGER_POINTER(VECTOR_ELT(ans, i))[j];
            R_isort(icard, li);
            for (j=0; j<li; j++)
                INTEGER_POINTER(VECTOR_ELT(ans, i))[j] = icard[j];
        }
    }

    UNPROTECT(pc);
    return(ans);
} */


int spOverlapC(double bbi1, double bbi2, double bbi3, double bbi4, double bbj1, double bbj2, double bbj3, double bbj4) {

    int overlap=1;

    if ((bbi1>bbj3) || (bbi2>bbj4) || 
		(bbi3<bbj1) || (bbi4<bbj2) ) {
		overlap = 0;
    }

    return(overlap);
}

int polypolyC(double *px1, double *py1, int n1, double *px2, double *py2,
    int n2, double sn, int crit) {
	int i, j, k=0;
	double dist;
	double x1, x2, y1, y2, xd, yd;

	for (i=0; (i < n1) && (k < crit); i++) {
		x1 = px1[i];
		y1 = py1[i];
		for (j=0; (j < n2) && (k < crit); j++) {
			x2 = px2[j];
			y2 = py2[j];
			xd = x1-x2;
			if (fabs(xd)>sn) { continue; }
			yd = y1-y2;
			if (fabs(yd)>sn) { continue; }
			dist = hypot(xd, yd);
			if (dist <= sn) {
                            k++;
                        }
		}
	}
	
	return(k);
}


SEXP poly_loop2(SEXP n, SEXP i_findInBox, SEXP bb, SEXP pl, SEXP nrs,
    SEXP dsnap, SEXP criterion, SEXP nfIBB) {

    int nn = INTEGER_POINTER(n)[0];
    int crit = INTEGER_POINTER(criterion)[0];
/*    int Scale = INTEGER_POINTER(scale)[0];*/
    int uBound = (int) INTEGER_POINTER(nfIBB)[0]*2;
    int i, j, jj, li, pc = 0;
    int ii = 0;
    int *card, *icard, *is, *jjs, *NRS, *cNRS;
    double *bb1, *bb2, *bb3, *bb4, *plx, *ply;
    double Dsnap = NUMERIC_POINTER(dsnap)[0];

//    struct bbcontainer *bbs;

    SEXP ans;

    int jhit, khit, nrsi, nrsj;

    int xx, yy, zz, ww;

    card = (int *) R_alloc((size_t) nn, sizeof(int));
    icard = (int *) R_alloc((size_t) nn, sizeof(int));
    is = (int *) R_alloc((size_t) uBound, sizeof(int));
    jjs = (int *) R_alloc((size_t) uBound, sizeof(int));
    bb1 = (double *) R_alloc((size_t) nn, sizeof(double));
    bb2 = (double *) R_alloc((size_t) nn, sizeof(double));
    bb3 = (double *) R_alloc((size_t) nn, sizeof(double));
    bb4 = (double *) R_alloc((size_t) nn, sizeof(double));
    NRS = (int *) R_alloc((size_t) nn, sizeof(int));
    cNRS = (int *) R_alloc((size_t) nn, sizeof(int));

    for (i=0, li=0; i<nn; i++) {
        card[i] = 0;
        icard[i] = 0;
        bb1[i] = NUMERIC_POINTER(bb)[i];
        bb2[i] = NUMERIC_POINTER(bb)[i+(1*nn)];
        bb3[i] = NUMERIC_POINTER(bb)[i+(2*nn)];
        bb4[i] = NUMERIC_POINTER(bb)[i+(3*nn)];
        NRS[i] = INTEGER_POINTER(nrs)[i];
        li += NRS[i];
    }

    for (i=0; i<nn; i++) {
        if (i == 0) cNRS[i] = 0;
        else cNRS[i] = NRS[i-1] + cNRS[i-1];
    }

    for (i=0; i<uBound; i++) {
        is[i] = 0;
        jjs[i] = 0;
    }

    plx = (double *) R_alloc((size_t) li, sizeof(double));
    ply = (double *) R_alloc((size_t) li, sizeof(double));

    for (i=0, jj=0; i<nn; i++) {
        nrsi = NRS[i];
        for (j=0; j<nrsi; j++) {
            plx[jj] = NUMERIC_POINTER(VECTOR_ELT(pl, i))[j];
            ply[jj] = NUMERIC_POINTER(VECTOR_ELT(pl, i))[j+nrsi];
            jj++;
/*            if (i < (nn-1) && jj == li) error("polygon memory overflow");*/
        }
    }

    for (i=0; i<(nn-1); i++) {
        li = length(VECTOR_ELT(i_findInBox, i));
        nrsi = NRS[i];
        for (j=0; j<li; j++) {
            jj = INTEGER_POINTER(VECTOR_ELT(i_findInBox, i))[j] - ROFFSET;
            jhit = spOverlapC(bb1[i], bb2[i], bb3[i], bb4[i], bb1[jj],
                bb2[jj], bb3[jj], bb4[jj]);
            if (jhit > 0) {
                khit = 0;
                nrsj = NRS[jj];
                if (nrsi > 0 && nrsj > 0){
                    khit = polypolyC(&plx[cNRS[i]], &ply[cNRS[i]], nrsi,
                       &plx[cNRS[jj]], &ply[cNRS[jj]], nrsj, Dsnap, crit+1L);
                }
                if (khit > crit) {
                    card[i]++;
                    card[jj]++;
                    is[ii] = i;
                    jjs[ii] = jj;
                    ii++;
/*                    if (ii == uBound) error("memory error, scale problem");*/
                }
            }
        }
    }

    PROTECT(ans = NEW_LIST(nn)); pc++;

    for (i=0; i<nn; i++) {
        if (card[i] == 0) {
            SET_VECTOR_ELT(ans, i, NEW_INTEGER(1));
            INTEGER_POINTER(VECTOR_ELT(ans, i))[0] = 0;
        } else {
            SET_VECTOR_ELT(ans, i, NEW_INTEGER(card[i]));
        }
    }

    for (i=0; i<ii; i++) {
        xx = is[i];
        yy = jjs[i];
        zz = icard[yy];
        ww = icard[xx];
/*        if (zz == card[yy]) error("memory error, overflow");
        if (ww == card[xx]) error("memory error, overflow");*/
        INTEGER_POINTER(VECTOR_ELT(ans, yy))[zz] = xx + ROFFSET;
        INTEGER_POINTER(VECTOR_ELT(ans, xx))[ww] = yy + ROFFSET;
        icard[yy]++;
        icard[xx]++;
    }

    for (i=0; i<nn; i++) {
        if ((li = length(VECTOR_ELT(ans, i))) > 1) {
            for (j=0; j<li; j++)
                icard[j] = INTEGER_POINTER(VECTOR_ELT(ans, i))[j];
            R_isort(icard, li);
            for (j=0; j<li; j++)
                INTEGER_POINTER(VECTOR_ELT(ans, i))[j] = icard[j];
        }
    }

    UNPROTECT(pc);
    return(ans);
}



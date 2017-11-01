/* Copyright 2006 by Roger S. Bivand. */

#include "spdep.h"

SEXP listw2dsT(SEXP nbs, SEXP wts, SEXP card, SEXP ncard2)
{
	int i, ii, j, jj, n, pc=0;
	SEXP ans;
/*	double *card; */

	n = LENGTH(nbs);
	PROTECT(ans = NEW_LIST(3)); pc++;
	SET_VECTOR_ELT(ans, 0, NEW_INTEGER(INTEGER_POINTER(ncard2)[0]));
	SET_VECTOR_ELT(ans, 1, NEW_INTEGER(INTEGER_POINTER(ncard2)[0]));
	SET_VECTOR_ELT(ans, 2, NEW_NUMERIC(INTEGER_POINTER(ncard2)[0]));

	for (i=0, ii=0; i < n; i++) {
	    for (j=0; j < INTEGER_POINTER(card)[i]; j++) {
		jj = INTEGER_POINTER(VECTOR_ELT(nbs, i))[j];
		if (jj > i) {
		    INTEGER_POINTER(VECTOR_ELT(ans, 0))[ii] = i;
	            INTEGER_POINTER(VECTOR_ELT(ans, 1))[ii] = jj-1;
		    NUMERIC_POINTER(VECTOR_ELT(ans, 2))[ii] = 
		        NUMERIC_POINTER(VECTOR_ELT(wts, i))[j]; 
		    if (ii >= INTEGER_POINTER(ncard2)[0])
			error("ncard2 incorrectly given");
		    ii++;
                }
	    }
	}
	UNPROTECT(pc); 
	return(ans);
}

SEXP listw2dgR(SEXP nbs, SEXP wts, SEXP card, SEXP ncard)
{
	int i, ii, j, jj, n, pc=0;
	SEXP ans;
/*	double *card; */

	n = LENGTH(nbs);
	PROTECT(ans = NEW_LIST(2)); pc++;
	SET_VECTOR_ELT(ans, 0, NEW_INTEGER(INTEGER_POINTER(ncard)[0]));
	SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(INTEGER_POINTER(ncard)[0]));

	for (i=0, ii=0; i < n; i++) {
	    for (j=0; j < INTEGER_POINTER(card)[i]; j++) {
		jj = INTEGER_POINTER(VECTOR_ELT(nbs, i))[j];
		INTEGER_POINTER(VECTOR_ELT(ans, 0))[ii] = jj-1;
		NUMERIC_POINTER(VECTOR_ELT(ans, 1))[ii] = 
		    NUMERIC_POINTER(VECTOR_ELT(wts, i))[j]; 
		if (ii >= INTEGER_POINTER(ncard)[0])
		    error("ncard incorrectly given");
		ii++;
	    }
	}
	UNPROTECT(pc); 
	return(ans);
}


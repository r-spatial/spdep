/* Copyright 2000-4 by Roger S. Bivand. */

#include "spdep.h"

SEXP listw2sn(SEXP nbs, SEXP wts, SEXP card, SEXP ncard)
{
	int i, ii, j, n, pc=0;
	SEXP ans;
/*	double *card; */

	n = LENGTH(nbs);
	PROTECT(ans = NEW_LIST(3)); pc++;
	SET_VECTOR_ELT(ans, 0, NEW_INTEGER(INTEGER_POINTER(ncard)[0]));
	SET_VECTOR_ELT(ans, 1, NEW_INTEGER(INTEGER_POINTER(ncard)[0]));
	SET_VECTOR_ELT(ans, 2, NEW_NUMERIC(INTEGER_POINTER(ncard)[0]));

	for (i=0, ii=0; i < n; i++) {
	    for (j=0; j < INTEGER_POINTER(card)[i]; j++) {
		INTEGER_POINTER(VECTOR_ELT(ans, 0))[ii] = i+ROFFSET;
	        INTEGER_POINTER(VECTOR_ELT(ans, 1))[ii] = 
		    INTEGER_POINTER(VECTOR_ELT(nbs, i))[j]; 
		NUMERIC_POINTER(VECTOR_ELT(ans, 2))[ii] = 
		    NUMERIC_POINTER(VECTOR_ELT(wts, i))[j]; 
		ii++;
	    }
	}
	UNPROTECT(pc); 
	return(ans);
}



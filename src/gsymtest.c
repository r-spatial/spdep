/* Copyright 2001-8 by Roger S. Bivand. */

#include "spdep.h"

SEXP gsymtest(SEXP nb, SEXP glist, SEXP card)
{
	int i, icard, j, k, k1, n=length(nb), pc=0, l=TRUE;
	double g, g0, d=0.0, d1=0.0;
	SEXP ans;
	PROTECT(ans = NEW_LIST(2)); pc++;
	SET_VECTOR_ELT(ans, 0, NEW_LOGICAL(1));
	SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(1));

	for (i=0; i < n; i++) {
	    icard = INTEGER_POINTER(card)[i];
	    for (j=0; j<icard; j++) {
		k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j];
		g = NUMERIC_POINTER(VECTOR_ELT(glist, i))[j];
		if (k > 0 && k <= n) {
		    for (k1=0; k1<INTEGER_POINTER(card)[k-ROFFSET]; k1++) {
			if (i+ROFFSET == INTEGER_POINTER(VECTOR_ELT(nb,
			    k-ROFFSET))[k1]) {
			    g0 = NUMERIC_POINTER(VECTOR_ELT(glist,
			        k-ROFFSET))[k1];
			    d = fabs(g - g0);
/* Rprintf("%d %d %f %f %f\n", i, j, g, g0, d); */
			    if (d > 0.0) {
				l = FALSE;
				if (d > d1) d1 = d;
			    }
			}
		    }
		}
	    }
	}

	LOGICAL_POINTER(VECTOR_ELT(ans, 0))[0] = l;
	NUMERIC_POINTER(VECTOR_ELT(ans, 1))[0] = d1;
	UNPROTECT(pc); /* ans */
	return(ans);
}


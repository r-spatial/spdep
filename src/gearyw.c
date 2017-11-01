/* Copyright 2002 by Roger S. Bivand. */

#include "spdep.h"

SEXP gearyw(SEXP nb, SEXP weights, SEXP x, SEXP card, SEXP zeropolicy,
	SEXP ftype) {
	int i, j, k, n=length(card), pc=0;
	double sum, wt, diff, xi, res;
	SEXP ans;
	PROTECT(ans = NEW_NUMERIC(n)); pc++;

	for (i=0; i < n; i++) {
	    if (INTEGER_POINTER(card)[i] == 0) {
		if (LOGICAL_POINTER(zeropolicy)[0] == TRUE)
		    NUMERIC_POINTER(ans)[i] = 0;
		else
		    NUMERIC_POINTER(ans)[i] = NA_REAL;
	    }
	    else {
		sum = 0.0;
		xi = NUMERIC_POINTER(x)[i];
		for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
		    k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j];
		    wt = NUMERIC_POINTER(VECTOR_ELT(weights, i))[j];
		    diff = (xi - NUMERIC_POINTER(x)[k-ROFFSET]);
		    if (LOGICAL_POINTER(ftype)[0] == TRUE) 
			res = diff*diff;
		    else {
			res = diff;
			if (res < 0.0) res = (-1) * res;
		    }
		    sum += wt * res;
		}
		NUMERIC_POINTER(ans)[i] = sum;
	    }
        }

	UNPROTECT(pc); /* ans */
	return(ans);
}


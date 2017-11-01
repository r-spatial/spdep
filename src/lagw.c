/* Copyright 2001 by Roger S. Bivand. */

#include "spdep.h"

SEXP lagw(SEXP nb, SEXP weights, SEXP x, SEXP card, SEXP zeropolicy,
   SEXP naok) {
	int i, j, k, n=length(card), pc=0, naOK=LOGICAL_POINTER(naok)[0],
            nas;
	double sum, wt, tmp;
	SEXP ans;
	PROTECT(ans = NEW_NUMERIC(n)); pc++;

        if (naOK == FALSE) {
            for (i=0; i < n; i++) 
                if (!R_FINITE(NUMERIC_POINTER(x)[i]))
                    error("Variable contains non-finite values");
        }

	for (i=0; i < n; i++) {
	    if (INTEGER_POINTER(card)[i] == 0) {
		if (LOGICAL_POINTER(zeropolicy)[0] == TRUE)
		    NUMERIC_POINTER(ans)[i] = 0;
		else
		    NUMERIC_POINTER(ans)[i] = NA_REAL;
	    }
	    else {
		sum = 0.0;
                nas = 0;
		for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
		    k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j];
		    wt = NUMERIC_POINTER(VECTOR_ELT(weights, i))[j];
		    tmp = NUMERIC_POINTER(x)[k-ROFFSET];
		    if (R_FINITE(tmp)) sum += tmp * wt;
                    else nas++;
		}
		if (nas == 0) NUMERIC_POINTER(ans)[i] = sum;
                else {
                    NUMERIC_POINTER(ans)[i] = NA_REAL;
/*                    warning("NA in lagged variable");*/
                }
	    }
        }

	UNPROTECT(pc); /* ans */
	return(ans);
}



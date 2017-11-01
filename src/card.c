/* Copyright 2000-2 by Roger S. Bivand. */

#include "spdep.h"

SEXP card(SEXP nb)
{
	int i, n=length(nb), pc=0, first_value, li;
	SEXP ans;
	PROTECT(ans = NEW_INTEGER(n)); pc++;

	for (i=0; i < n; i++) {

            li = length(VECTOR_ELT(nb, i));
            if (li > 0) 
                first_value = INTEGER_POINTER(VECTOR_ELT(nb, i))[0];
            else
                error("zero length neighbour vector");
            
	    if (first_value == 0) 
		INTEGER_POINTER(ans)[i] = 0;
	    else
		INTEGER_POINTER(ans)[i] = li;
	}

	UNPROTECT(pc); /* ans */
	return(ans);
}


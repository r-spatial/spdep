/* Copyright 2000-3 by Roger S. Bivand. */

#include "spdep.h"

SEXP nbdists(SEXP nb, SEXP x, SEXP np, SEXP dim, SEXP lonlat)
{
	int i, j, j1, k, /* m,*/ n, d, pc=0, ll, first_value;
	SEXP ans;
        SEXP class;
	double lon1[1], lon2[1], lat1[1], lat2[1], gc[1];
	double tmp /*, tmp1*/;
	
	PROTECT(ans = NEW_LIST(1)); pc++;
	n = INTEGER_POINTER(np)[0];
	ll = INTEGER_POINTER(lonlat)[0];

	SET_VECTOR_ELT(ans, 0, NEW_LIST(n));
	d = INTEGER_POINTER(dim)[0];
	if (d > 2) error("Only 2D coordinates allowed");
	PROTECT(class = NEW_CHARACTER(1)); pc++;
	SET_STRING_ELT(class, 0, COPY_TO_USER_STRING("nbdist"));
	setAttrib(VECTOR_ELT(ans, 0), R_ClassSymbol, class);

	for (i=0; i < n; i++) {
                R_CheckUserInterrupt();
		k = length(VECTOR_ELT(nb, i));
/*		if (k == 1 && INTEGER_POINTER(VECTOR_ELT(nb, i))[0] == 0) {
			SET_VECTOR_ELT(VECTOR_ELT(ans, 0), i,
				NEW_NUMERIC(1));
			NUMERIC_POINTER(VECTOR_ELT(VECTOR_ELT(ans, 0), i))[0]
				= NA_REAL;
		} else { */

                if (k > 0) 
                    first_value = INTEGER_POINTER(VECTOR_ELT(nb, i))[0];
                else
                    error("zero length neighbour vector");

		if (first_value > 0) {
			SET_VECTOR_ELT(VECTOR_ELT(ans, 0), i,
				NEW_NUMERIC(k));
			for (j=0; j < k; j++) {
				j1 = INTEGER_POINTER(VECTOR_ELT(nb, i))[j]
					- ROFFSET;
			/*	tmp1 = 0;
				for (m=0; m < d; m++) {
					tmp = NUMERIC_POINTER(x)[i + m * n]
					- NUMERIC_POINTER(x)[j1 + m * n];
					tmp1 += tmp * tmp;
				} */
	    			lon1[0] = NUMERIC_POINTER(x)[i];
	    			lat1[0] = NUMERIC_POINTER(x)[i + n];
	    			lon2[0] = NUMERIC_POINTER(x)[j1];
	    			lat2[0] = NUMERIC_POINTER(x)[j1 + n];
	    			if (ll == 0) 
					tmp = hypot((lon1[0]-lon2[0]), 
							(lat1[0]-lat2[0]));
	    			else {
		    			gcdist(lon1, lon2, lat1, lat2, gc);
		    			tmp = gc[0];
	    			}

				NUMERIC_POINTER(VECTOR_ELT(VECTOR_ELT(ans, 0),
					i))[j] = tmp;
			}
		}
	}
	UNPROTECT(pc);
	return(ans);
}


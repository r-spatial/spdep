/*
 *  based on code taken from:
 *  class/class.c by W. N. Venables and B. D. Ripley  Copyright (C) 1994-9
 *  and written by Roger Bivand (C) 2001-3
 */

#include "spdep.h"

#define DOUBLE_XMAX DBL_MAX

#define EPS 1e-4		/* relative test of equality of distances */

#define MAX_TIES 1000
/* Not worth doing this dynamically -- limits k + # ties + fence, in fact */


void
knearneigh(int *kin, int *pnte, int *p, double *test, int *res, double *dists,
		int *lonlat)
{
    int   j, k, k1, kinit = *kin, kn, npat, nte = *pnte, ll = *lonlat;
    int   pos[MAX_TIES];
    double dist, /* tmp,*/ nndist[MAX_TIES];
    double lon1[1], lon2[1], lat1[1], lat2[1], gc[1];

/*
    Use a `fence' in the (k+1)st position to avoid special cases.
    Simple insertion sort will suffice since k will be small.
 */

    for (npat = 0; npat < nte; npat++) {
        R_CheckUserInterrupt();
	kn = kinit;
	for (k = 0; k < kn; k++)
	    nndist[k] = 0.99 * DOUBLE_XMAX;
	for (j = 0; j < nte; j++) {
	    if (j == npat) continue;
	    lon1[0] = test[npat];
	    lat1[0] = test[npat + nte];
	    lon2[0] = test[j];
	    lat2[0] = test[j + nte];
	    if (ll == 0) dist = hypot((lon1[0]-lon2[0]), (lat1[0]-lat2[0]));
	    else {
		    gcdist(lon1, lon2, lat1, lat2, gc);
		    dist = gc[0];
	    }
/*	    dist = 0.0;
	    for (k = 0; k < *p; k++) {
		tmp = test[npat + k * nte] - test[j + k * nte];
		dist += tmp * tmp;
	    }*/
	    
/* Use `fuzz' since distance computed could depend on order of coordinates */
	    if (dist <= nndist[kinit - 1] * (1 + EPS))
		for (k = 0; k <= kn; k++)
		    if (dist < nndist[k]) {
			for (k1 = kn; k1 > k; k1--) {
			    nndist[k1] = nndist[k1 - 1];
			    pos[k1] = pos[k1 - 1];
			}
			nndist[k] = dist;
			pos[k] = j;
/* Keep an extra distance if the largest current one ties with current kth */
			if (nndist[kn] <= nndist[kinit - 1])
			    if (++kn == MAX_TIES - 1)
				error("too many ties in knearneigh");
			break;
		    }
	    nndist[kn] = 0.99 * DOUBLE_XMAX;
	}
	for (k = 0; k < kinit; k++) {
	    res[k + (npat*kinit)] = pos[k]+1;
	    dists[k + (npat*kinit)] = nndist[k];
	}
    }
}



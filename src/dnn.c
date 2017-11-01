/*
 *  based on code taken from:
 *  class/class.c by W. N. Venables and B. D. Ripley  Copyright (C) 1994-9
 *  and written by Roger Bivand (C) 2001-2014
 */

#include "spdep.h"

SEXP
dnearneigh(SEXP din1, SEXP din2, SEXP pnte, SEXP p, SEXP test, SEXP lonlat)
{
    int   j, k, kn, npat, nte, pc=0, ll;
    int   *pos;
    double dist, /*tmp,*/ dn, dn0;
    int dn0_equal, dn_equal;
    double lon1[1], lon2[1], lat1[1], lat2[1], gc[1];
    SEXP ans;
    SEXP class;
    SEXP nbtype;
    SEXP dists;
    
    dn0 = NUMERIC_POINTER(din1)[0];
    dn = NUMERIC_POINTER(din2)[0];
    dn0_equal = LOGICAL_POINTER(getAttrib(din1, install("equal")))[0];
    dn_equal = LOGICAL_POINTER(getAttrib(din2, install("equal")))[0];
    nte = INTEGER_POINTER(pnte)[0];
/*    pdim = INTEGER_POINTER(p)[0];*/
    ll = INTEGER_POINTER(lonlat)[0];
    PROTECT(ans = NEW_LIST(1)); pc++;
    PROTECT(dists = NEW_NUMERIC(2)); pc++;
    NUMERIC_POINTER(dists)[0] = dn0;
    NUMERIC_POINTER(dists)[1] = dn;
    SET_VECTOR_ELT(ans, 0, NEW_LIST(nte));
    PROTECT(class = NEW_CHARACTER(1)); pc++;
    PROTECT(nbtype = NEW_CHARACTER(1)); pc++;
    SET_STRING_ELT(class, 0, COPY_TO_USER_STRING("nb"));
    SET_STRING_ELT(nbtype, 0, COPY_TO_USER_STRING("distance"));
    setAttrib(VECTOR_ELT(ans, 0), R_ClassSymbol, class);
    setAttrib(VECTOR_ELT(ans, 0), install("nbtype"), nbtype);
    setAttrib(VECTOR_ELT(ans, 0), install("distances"), dists);
/*    dn0 = dn0;
    dn = dn;*/
    pos = (int *) R_alloc((size_t) nte, sizeof(int));
    for (npat = 0; npat < nte; npat++) {
        R_CheckUserInterrupt();
	kn = 0;
	lon1[0] = NUMERIC_POINTER(test)[npat];
	lat1[0] = NUMERIC_POINTER(test)[npat + nte];
	for (j = 0; j < nte; j++) {
	    if (j == npat) continue;
/*	    dist = 0.0;
	    for (k = 0; k < pdim; k++) {
		tmp = NUMERIC_POINTER(test)[npat + k * nte]
			- NUMERIC_POINTER(test)[j + k * nte];
		dist += tmp * tmp;
	    } */
	    lon2[0] = NUMERIC_POINTER(test)[j];
	    lat2[0] = NUMERIC_POINTER(test)[j + nte];
//Rprintf("ll: %d\n", ll);
	    if (ll == 0) dist = hypot((lon1[0]-lon2[0]), (lat1[0]-lat2[0]));
	    else {
		    gcdist(lon1, lon2, lat1, lat2, gc);
		    dist = gc[0];
//Rprintf("npat: %d, j: %d, dist: %f\n", npat, j, dist);
	    }
	    if ((dn0_equal ? dist >= dn0 : dist > dn0)
                && (dn_equal ? dist <= dn: dist < dn)) {
		pos[kn] = j;
		if (++kn == nte - 1 && j == nte) {
			Rprintf("%d %d %d\n", kn, nte, j);
		    error("position array overrun");
		}
	    }
	}
	if (kn < 1) {
	    SET_VECTOR_ELT(VECTOR_ELT(ans, 0), npat, NEW_INTEGER(1));
	    INTEGER_POINTER(VECTOR_ELT(VECTOR_ELT(ans, 0), npat))[0] = 0;
	} else {
	    SET_VECTOR_ELT(VECTOR_ELT(ans, 0), npat, NEW_INTEGER(kn));
	    for (k = 0; k < kn; k++) {
	    	INTEGER_POINTER(VECTOR_ELT(VECTOR_ELT(ans, 0), npat))[k]
		    = pos[k]+ROFFSET;
	    }
	}
    }
    UNPROTECT(pc); 
    return(ans);
}

/* http://home.att.net/~srschmitt/greatcircle.html */

void gcdist(double *lon1, double *lon2, double *lat1, double *lat2, 
		double *dist) {
	
    double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C;
    double w, R, a, f, D, H1, H2;
    double lat1R, lat2R, lon1R, lon2R, DE2RA;

// Maeel Le Noc bug 2017-04-12

    if (fabs(lat1[0] - lat2[0]) < DOUBLE_EPS) {
        if (fabs(lon1[0] - lon2[0]) < DOUBLE_EPS) {
            dist[0] = 0.0;
            return;
/* Wouter Buytaert bug caught 100211 */
        } else if (fabs((fabs(lon1[0]) + fabs(lon2[0])) - 360.0) < DOUBLE_EPS) {
            dist[0] = 0.0;
            return;
        }
    }
/*    if (lat1[0] == lat2[0] && lon1[0] == lon2[0]) {
      dist[0] = 0.0;
    } else {*/
    
    DE2RA = M_PI/180;
    a = 6378.137;              /* WGS-84 equatorial radius in km */
    f = 1.0/298.257223563;     /* WGS-84 ellipsoid flattening factor */
    
    lat1R = lat1[0]*DE2RA;
    lat2R = lat2[0]*DE2RA;
    lon1R = lon1[0]*DE2RA;
    lon2R = lon2[0]*DE2RA;
    
    F = ( lat1R + lat2R )/2.0;
    G = ( lat1R - lat2R )/2.0;
    L = ( lon1R - lon2R )/2.0;

    sinG2 = R_pow_di( sin( G ), 2 );
    cosG2 = R_pow_di( cos( G ), 2 );
    sinF2 = R_pow_di( sin( F ), 2 );
    cosF2 = R_pow_di( cos( F ), 2 );
    sinL2 = R_pow_di( sin( L ), 2 );
    cosL2 = R_pow_di( cos( L ), 2 );

    S = sinG2*cosL2 + cosF2*sinL2;
    C = cosG2*cosL2 + sinF2*sinL2;

    w = atan( sqrt( S/C ) );
    R = sqrt( S*C )/w;

    D = 2*w*a;
    H1 = ( 3*R - 1 )/( 2*C );
// was H2 = ( 3*R + 2 )/( 2*S ); 161209
    H2 = ( 3*R + 1 )/( 2*S );

    dist[0] = D*( 1 + f*H1*sinF2*cosG2 - f*H2*cosF2*sinG2 ); 
//    }
}


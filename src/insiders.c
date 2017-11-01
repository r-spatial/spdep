#include "spdep.h"

int pipbb(double pt1, double pt2, double *bbs);

int between(double x, double low, double up); 

SEXP spInsiders(SEXP bbbi, SEXP bbbj) {

	int pc=0;
	int k1;
	double bbi[4], bbj[4];
	int jhit[8], hsum;
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(1)); pc++;
	bbi[0] = NUMERIC_POINTER(bbbi)[0];
	bbi[1] = NUMERIC_POINTER(bbbi)[1];
	bbi[2] = NUMERIC_POINTER(bbbi)[2];
	bbi[3] = NUMERIC_POINTER(bbbi)[3];

	hsum = 0;
	bbj[0] = NUMERIC_POINTER(bbbj)[0];
	bbj[1] = NUMERIC_POINTER(bbbj)[1];
	bbj[2] = NUMERIC_POINTER(bbbj)[2];
	bbj[3] = NUMERIC_POINTER(bbbj)[3];
	for (k1=0; k1 < 8; k1++) jhit[k1] = 0;
	jhit[0] = pipbb(bbi[2], bbi[3], bbj);
    	jhit[1] = pipbb(bbi[0], bbi[1], bbj);
	jhit[2] = pipbb(bbi[0], bbi[3], bbj);
	jhit[3] = pipbb(bbi[2], bbi[1], bbj);
	jhit[4] = pipbb(bbj[2], bbj[3], bbi);
    	jhit[5] = pipbb(bbj[0], bbj[1], bbi);
	jhit[6] = pipbb(bbj[0], bbj[3], bbi);
	jhit[7] = pipbb(bbj[2], bbj[1], bbi);

	for (k1=0; k1 < 8; k1++) hsum = hsum + jhit[k1];
	INTEGER_POINTER(ans)[0] = hsum;		
	UNPROTECT(pc); /* ans */
	return(ans);
}


int between(double x, double low, double up) {
	if (x >= low && x <= up) return(1);
	else return(0);
}

int pipbb(double pt1, double pt2, double *bbs) {
	if ((between(pt1, bbs[0], bbs[2]) == 1) && 
		(between(pt2, bbs[1], bbs[3]) == 1)) return(1);
	else return(0);
} 


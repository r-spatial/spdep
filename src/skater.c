#include "spdep.h"

void prunemst(int *e1, int *e2, int *ne, int *gr) {
  int i, j, no1[*ne], n1=1, li=0, ls=1;
  no1[0] = e1[0];
  for (i=0; i<*ne; i++)
    gr[i] = 0;
  do {
    for (i=li; i<ls; i++) {
      for (j=1; j<*ne; j++) {
	if (gr[j] == 0) {
	  if (no1[i]==e1[j]) {
	    gr[j] = 1;
	    no1[n1++] = e2[j];
	  }
	  if (no1[i]==e2[j]) {
	    gr[j] = 1;
	    no1[n1++] = e1[j];
	  }
	}
      }
    }
    li = ls;
    ls = n1;
  } while(li < ls);
}



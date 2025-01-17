/* Copyright 2001 by Nicholas Lewin-Koh.
 * interrupt added RSB 2024 */

#include "spdep.h"

#define BLACK 1
#define WHITE 0

void dfs(SEXP nblst, SEXP cmpnm, SEXP visited, int curcmp, int nodeid){
  int n,i,nvert;

  if (nodeid < 0 || nodeid >= Rf_length(nblst))
    Rf_error("n.comp.nb: recursive nodeid out of bounds\ninstall spatialreg and igraph to avoid legacy issue");
  

  INTEGER(cmpnm)[nodeid]=curcmp;
  INTEGER(visited)[nodeid]=BLACK;
  n=Rf_length(VECTOR_ELT(nblst,nodeid));

  for(i=0;i<n;i++){
    if(INTEGER(visited)[(INTEGER(VECTOR_ELT(nblst,nodeid))[i]-1)]==WHITE){ 
      dfs(nblst, cmpnm, visited, curcmp,
	  INTEGER(VECTOR_ELT(nblst,nodeid))[i]-1 );
    }
  }
}


SEXP g_components(SEXP nblst, SEXP cmpnm){
  int i, curcmp=1, nvert;
  SEXP visited;
  
  nvert=Rf_length(nblst);
  PROTECT(visited=Rf_allocVector(INTSXP,nvert));
  
  for(i=0; i < nvert; i++){
    INTEGER(visited)[i]=WHITE;
  }

  for(i=0; i < nvert; i++){
    R_CheckUserInterrupt();
    if(INTEGER(visited)[i]==WHITE){ 
      INTEGER(visited)[i]=BLACK;
      if(INTEGER(VECTOR_ELT(nblst,i))[0]==0){
	INTEGER(cmpnm)[i]=curcmp;
	curcmp++;
      }
      else{
	dfs(nblst, cmpnm, visited, curcmp, i);
	curcmp++;
      }
    }    
  }
  UNPROTECT(1);
  return(cmpnm);
}

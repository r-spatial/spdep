/* Copyright 2001 by Nicholas Lewin-Koh. */


#include "spdep.h"


static double distance(double x1, double y1, double x2, double y2){

  return(hypot(x1-x2,y1-y2));

}

void compute_relative(int *no_nodes, int *g1, int *g2, int *nogab,
      int *ngaballoc, double *nodes_xd, double *nodes_yd)
{
  int i,j,l, no_gab=0;
  double rad;

  for(i=0;i<*no_nodes;i++)
    {
      for(j=i+1;j<*no_nodes;j++)
	{
	  rad=distance(nodes_xd[i],nodes_yd[i],nodes_xd[j],nodes_yd[j]);
          /*Rprintf("hi \n");*/
	  for(l=0;l<*no_nodes;l++)
	    {
	      if((l!=i)&&(l!=j)&&
	      (distance(nodes_xd[i],nodes_yd[i],nodes_xd[l],nodes_yd[l])<rad)&&
              (distance(nodes_xd[j],nodes_yd[j],nodes_xd[l],nodes_yd[l])<rad))
		break;
	    }

/* bug Dan Putler 090121 */
	  if ((no_gab+1) > *ngaballoc) 
		error("number of neighbours overrun - increase nnmult");
	  if(l==*no_nodes)
	    {
	      g1[no_gab]=i+1; g2[no_gab++]=j+1;
	    }
	}
    }
  *nogab=no_gab;
  return;

}

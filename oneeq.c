/* ************************************************************************
 ** oneeq.c            Hilbert bases package version 0.1 (experimental)
 ** 
 ** Copyright (C) 1995 Dmitrii V. Pasechnik
 **                    RIACA, Amsterdam, The Netherlands
 **
 **                    front-end for hb.c
 ** (the solver of one (maybe inhomogeneous) linear Diophantine equation
 **
 ** Note1: injective companion trick is not used here; one day I'll 
 ** make it proper.
 **
 ** Note2: input format is unnecessarily restrictive here;
 **
* ************************************************************************ */

#include <time.h>
#include <stdio.h>
#include "hb.h"
main(int argc, char *argv[])
{
  int *x, *a, n, m, c, num, savesols, i;
  FILE *outfile;
  clock_t clock0;
  eq *e;

  if (argc==1) outfile=stdout;
  else outfile=fopen(argv[1],"w");
  clock0=clock();	
  scanf("%d%d%d%d%d",&prtlev,&savesols,&n,&m,&c);
  if (prtlev)
    {
      printf("\n inputting %d positive coefficients\n",n);
      printf(" and %d (assumed) negative coefficients\n",m);
      printf(" right-hand side %d\n",c);
      if (savesols) printf(" solutions will be saved\n");
      else printf(" just counting solutions\n");
    }
  if ((e=(eq *)calloc(1,sizeof(eq)))==NULL) oomem("oneeq(e)");
  e->a=getmat(stdin,1,n);
  e->b=getmat(stdin,1,m);
  e->c=c;
  e->n=n;
  e->m=m;

  if (c<0) {printf("\n right-hand side must be nonnegative. exit\n");
	  exit(1);
	  }
  for (i=0; i<n; i++) if (e->a[i]<=0)
    { printf("\n coefficient must be positive. exit \n"); exit(1); }
  for (i=0; i<m; i++) if (e->b[i]<=0)
    { printf("\n coefficient must be positive. exit \n"); exit(1); }

  if (prtlev) 
    { printf("\n positive coefficients:");
      outmat(stdout,1,n,e->a);
      printf("\n negative coefficients:");
      outmat(stdout,1,m,e->b);
    }

  printf("\n %d elements:\n",(num=hb(e,&x,savesols)));
  if (savesols && argc>1) printf(" writing the file %s\n",argv[1]);
/* the following doesnt work in SunOS */
/*  printf("\n used %d CPU sec.\n",(clock()-clock0)/CLOCKS_PER_SEC); */
  if (savesols && num) outmat(outfile,num,m+n,x);

}



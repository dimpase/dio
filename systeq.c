/* ************************************************************************
 ** systeq.c           Hilbert bases package version 0.1 (experimental)
 ** 
 ** Copyright (C) 1995 Dmitrii V. Pasechnik
 **                    RIACA, Amsterdam, The Netherlands
 **
 ** the front-end for the solver of systems of homogeneous linear
 ** Diopohantine equations
* ************************************************************************ */

#include <time.h>
#include <stdio.h>
#include "hb.h"
main(int argc, char *argv[])
{
  int *x, *a, n, m, num;
  FILE *outfile;
  clock_t clock0;

  if (argc==1) outfile=stdout;
  else outfile=fopen(argv[1],"w");
  clock0=clock();	
  scanf("%d%d%d",&prtlev,&n,&m);
  if (prtlev) printf("\n inputting %d x %d matrix\n",n,m);
  a=getmat(stdin,n,m);
  if (prtlev) outmat(stdout,n,m,a);
  printf("\n %d elements:\n",(num=hbs(n,m,a,&x)));
  if (argc>1) printf(" writing the file %s\n",argv[1]);
/* the following doesnt work in SunOS */
/*  printf("\n used %d CPU sec.\n",(clock()-clock0)/CLOCKS_PER_SEC); */
  if (num) outmat(outfile,num,m,x);
}


/* ************************************************************************
 ** hb_uts.c           Hilbert bases package version 0.1 (experimental)
 ** 
 ** Copyright (C) 1995 Dmitrii V. Pasechnik
 **                    RIACA, Amsterdam, The Netherlands
 **
 **                Utilites...
* ************************************************************************ */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hb.h"

/* *******************************************************************
   counting the number of k-subsets in n-element set
*/
int choose(int n, int k)
{
  int i,r;
  if (k<0 || n<k || n<0) return 0;
  if (!k) return 1;
  if (n-k<k) return choose(n,n-k);
  for (i=0,r=1; i<k; i++) { r*=n-i; r/=i+1; }
  return r;
}

/* *******************************************************************
   counting the number of nonneg. solutions the equation $\sum_1^n x_i = c $
*/
int nf(int n, int c)
{
  return choose(n+c-1,c);
}

/* *******************************************************************
   solving the equation $\sum_1^n x_i = c $
   given x, return the lex. next solution
*/
int nexpa(int n, int x[], int c) /* return 1 is success, 0 if no way to go */
{
  int i;
  
  for (i=n-1; i && !x[i]; i--); /* find 1st non-0 from the left */
  if (!i) return 0;
  x[i-1]++;
  x[n-1]=x[i]-1;
  if (i != n-1) x[i]=0;
  return 1;
}

/* *******************************************************************
   for the equation $\sum_1^n x_i = c $
   return the lex. 1st solution
*/
void firpa(int n, int x[], int c) /* return 1 is success, 0 if no way to go */
{
  int i;
  bzero(x,(n-1)*sizeof(int));
  x[n-1]=c;
}

/* ******************************************************************
   get rid of an equation
*/
void kieq(eq *e)
{
  free(e->a);
  free(e->b);
  free(e);
}

int comp_decr(const int *a, const int *b)
{ 
  register int r;
  r=*a-*b;
  if (r>0) return -1;
  if (r<0) return 1;
  return 0;
}


/* ******************************************************************
   the injective companion
   (assuming c>=0)
*/
eq *injcom(int n, int v[], int c)
{
  int i,j, *p, *d;
  eq *e;

  if ((d=(int *)calloc(n,sizeof(int)))==NULL) oomem("injcom");
  if ((e=(eq *)calloc(1,sizeof(eq)))==NULL) oomem("injcom");
  e->c=c; 

  for (i=0; i<n; i++) /* compute sizes of e.a and e.b */
    if (!d[i])        /* new element in v */
      {
	if (v[i]>0) e->n++;
	if (v[i]<0) e->m++;
	for (j=i; j<n; j++) if (v[j]==v[i]) d[j]=1;
      }

  for (i=0; i<n; i++) d[i]=0;

  if (e->n)
    if ((e->a=(int *)calloc(e->n,sizeof(int)))==NULL) oomem("injcom");
  if (e->m)
    if ((e->b=(int *)calloc(e->m,sizeof(int)))==NULL) oomem("injcom");

  for (i=0,e->n=0,e->m=0; i<n; i++) /* write e.a and e.b */
    if (!d[i]) /* new element in v */
      {
	if (v[i]>0) e->a[e->n++]=v[i];
	if (v[i]<0) e->b[e->m++]=-v[i];
	for (j=i; j<n; j++) if (v[j]==v[i]) d[j]=1;
      }
  /* sort e->a and e->b in decreasing order */
/*  qsort(e->a,e->n,sizeof(int),comp_decr); */
/*  qsort(e->b,e->m,sizeof(int),comp_decr); */
  free(d);
  return e;
}

/* *****************************************************************
   inputting an integer matrix from already opened file
*/
int *getmat(FILE *file, int n, int m)
{
  int i, j, *a, *p;
  if (n<0 || m<0)
    { fprintf(stderr,"\n wrong dimesions in getmat\n"); exit(1);}
  if (m*n==0) return NULL;
  if (NULL==(a=(int *)calloc(n*m,sizeof(int)))) oomem("getmat");

  for (i=0, p=a; i<n; i++)
    for (j=0; j<m; j++) 
      if (!fscanf(file,"%d",p++)) 
	{
	  fprintf(stderr,"\n input error in getmat\n");
	  exit(1);
	}
  return a;
}

/* *****************************************************************
   output of an integer matrix to already opened file
*/
void outmat(FILE *file, int n, int m, int *a)
{
  int i, j, *p;
  if (n<0 || m<0)
    { fprintf(stderr,"\n wrong dimesions in outmat\n"); exit(1);}
  if (n*m)
   {for (i=0, p=a; i<n; i++)
    {    
      fprintf(file,"\n");
      for (j=0; j<m; j++) fprintf(file,"%4d",*p++);
    }
   }
   fprintf(file,"\n");
} 


/* *******************************************************************
   solving the equation $\sum_1^n x_i = c $
   given x, return the lex. next solution
   x is not a usual array, but is given as an array of pointers
*/
int nex_x(int n, int *x[], int c) /* return 1 is success, 0 if no way to go */
{
  int i;
  
  for (i=n-1; i && !*(x[i]); i--); /* find 1st non-0 from the left */
  if (!i) return 0;
  (*(x[i-1]))++;
  (*(x[n-1]))=(*(x[i]))-1;
  if (i != n-1) (*(x[i]))=0;
  return 1;
}

/* *******************************************************************
   for the equation $\sum_1^n x_i = c $
   return the lex. 1st solution
   x is not a usual array, but is given as an array of pointers
*/
void fir_x(int n, int *x[], int c) /* return 1 is success, 0 if no way to go */
{
  int i;
  for (i=0; i<n; *(x[i++])=0);
  *(x[n-1])=c;
}

/* *******************************************************************
   reporting out of memory and quitting
*/
void oomem(char *t)
{
  fprintf(stderr,"out of dynamic memory at %s\n",t);
  exit(10);
}

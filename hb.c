/* ************************************************************************
 ** hb.c               Hilbert bases package version 0.1 (experimental)
 ** 
 ** Copyright (C) 1995 Dmitrii V. Pasechnik
 **                    RIACA, Amsterdam, The Netherlands
 **
 **   construct the Hilbert basis for the solutions of an inhomogeneous 
 **  diophantine equation
 **  (Clausen - Fortenbacher algorithm)
* *********************************************************************** */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hb.h"
struct step {int node; int edge;};

int hb(eq *e, int **sols, int savesols) /* return number elts in the basis; */
{                        /* or -1 if out of memory           */
  int nvert,amon,bmon,i,j,n,m,lsol,nnodes,
      *sol,*mark,amax,bmax,bcmax,nsol,lp;
  struct step *p, *p_t;
  int *mark_t, *sol_t, ntests;
  int w,*x,*pb,sign,nx,offset;
  int ismsol(eq *,int *,int *,struct step *,int *);
  void saves(int, int *, int, int, int *, int **), oomem(char *);
  int solmax;

  lsol=e->n + e->m;      /* length of the vectors */
  amax=getmax(e->n, e->a);
  bmax=getmax(e->m, e->b);
  if (e->c>bmax) bcmax=e->c;
  else bcmax=bmax;
  /* vertices of the graph are labelled by [-bcmax..amax] */
  nvert=amax+bcmax+1; /* number of vertices in the graph   */

  if (savesols) /* allocate inital chunk for solutions */
    {
      solmax=SOLCHUNK;
      if ((*sols=(int *)malloc(lsol*solmax*sizeof(int)))==NULL)
	oomem("hb(sols)");
    }

  if ((sol=(int *)calloc(lsol, sizeof(int)))==NULL) 
    oomem("hb"); /* to store (partial) solution */
  if ((p=(struct step *)calloc(nvert+1, sizeof(struct step)))==NULL)
    oomem("hb"); /* the path */
  if ((mark=(int *)calloc(nvert, sizeof(int)))==NULL) 
    oomem("hb"); /* mark visited vertices */
  if ((sol_t=(int *)calloc(lsol, sizeof(int)))==NULL) 
    oomem("hb"); /* to store (partial) solution */
  if ((p_t=(struct step *)calloc(nvert+1, sizeof(struct step)))==NULL)
    oomem("hb");  /* the path */
  if ((mark_t=(int *)calloc(nvert, sizeof(int)))==NULL)
    oomem("hb");  /* mark visited vertices */

  ntests=0;
  mark+=bcmax;             /* adjust to use the original node names */
  if (e->c) mark[-e->c]=1; /* mark the visited node (don't mark if 0) */
  lp=0;                    /* length of the path (# of edges) */
  nsol=0;                  /* number of elts in the basis */
  amon=0; bmon=0;          /* to mantain monotonousness */
  p[0].node=-e->c;         /* 1st vertex of the path is -c,trying to reach 0 */
  do
    {
      if (!p[lp].node && lp) /* found a solution */
	{
	    ntests++;
	  if (ismsol(e,sol,sol_t,p_t,mark_t)) 
	    saves(savesols, &solmax, lsol, ++nsol, sol, sols);
	  goto getnext;
	}
      else /* forward ! */
	{
	  if (p[lp].node<=0) 
	           { offset=0; nx=e->n; x=e->a; pb=&amon; sign=1;}
	  else  { offset=e->n; nx=e->m; x=e->b; pb=&bmon; sign=-1;}
	  for (w=*pb; w<nx && mark[p[lp].node+sign*x[w]]; w++);
	  if (w<nx) /* success */
	    {
	      p[lp].edge=w;
	      sol[offset+w]++;
	      *pb=w;
	      p[lp+1].node=p[lp].node+sign*x[w];
	      mark[p[++lp].node]=1;
	      ntests++;
	      if (!ismsol(e,sol,sol_t,p_t,mark_t)) goto getnext;
	    }
	  else 
	  getnext:
	    {
	      while (lp>=1)
		{
		  mark[p[lp].node]=0;
		  if (p[lp-1].node<=0) 
		    { offset=0; nx=e->n; x=e->a; pb=&amon; sign=1;}
		  else  { offset=e->n; nx=e->m; x=e->b; pb=&bmon; sign=-1;}
		  sol[offset+p[lp-1].edge]--;
		  for (w=p[lp-1].edge+1; 
		       w<nx && mark[p[lp-1].node+sign*x[w]]; w++);
		  if (w<nx) /* success */
		    {
		      p[lp-1].edge=w;
		      sol[offset+w]++;
		      *pb=w;
		      p[lp].node=p[lp-1].node+sign*x[w];
		      mark[p[lp].node]=1;
		      ntests++;
		      if (!ismsol(e,sol,sol_t,p_t,mark_t)) goto getnext;
		      goto cont;
		    }
		  else 
		    {
		      lp--;
		      /* reset *pb */
		      for (w=lp-1; w>=0 && (sig(p[w].node) != sign); w--);
		      if (w<0) *pb=0;
		      else *pb=p[w].edge;
		    }
		}
	    cont:
            ;
	    }
	}
    }
  while(lp);
  free(sol);
  free(p);
  free(mark-bcmax);
  free(sol_t);
  free(p_t);
  free(mark_t);
  if (prtlev>1) printf("\n number of nonmin. (partial) solutions found %d\n",
	   ntests - nsol);

  /* getting  rid of extra space */ 
  if (nsol)
    {
      if ((*sols=(int *)realloc(*sols,lsol*nsol*sizeof(int)))==NULL)
      oomem("hb(realloc)");
    }
  else free(*sols); /* no solutions found */
  return nsol;
}


/* ***************************************************************
   find maximum of an array of ints
*/
int getmax(int n, int a[])
{
  int i, maxa;
  for(maxa=0,i=0; i<n; i++)
    if (maxa<a[i]) maxa=a[i];
  return maxa;
}

/* ***************************************************************
   check minimality of a solution (1-minimal, 0-not minimal)
*/
int ismsol(eq *e, int bound[], int sol[], struct step p[], int mark[])
{                   
  int nvert,amon,bmon,i,j,n,m,lsol,nnodes,
      amax,bmax,bcmax,nsol,lp;
  int w,*x,*pb,sign,nx,offset;

  lsol=e->n + e->m;      /* length of the vectors */
  amax=getmax(e->n, e->a);
  bmax=getmax(e->m, e->b);
  if (e->c>bmax) bcmax=e->c;
  else bcmax=bmax;
  /* vertices of the graph are labelled by [-bcmax..amax] */
  nvert=amax+bcmax+1; /* number of vertces in the graph   */

  bzero(sol,lsol*sizeof(int));
  bzero(p,(nvert+1)*sizeof(struct step));
  bzero(mark,nvert*sizeof(int));
  
  mark+=bcmax;             /* adjust to use the original node names */
  if (e->c) mark[-e->c]=1; /* mark the visited node (don't mark if 0) */
  lp=0;                    /* length of the path (# of edges) */
  nsol=0;                  /* number of elts in the basis */
  amon=0; bmon=0;          /* to mantain monotonousness */
  p[0].node=-e->c;         /* 1st vertex of the path is -c,trying to reach 0 */
  do
    {
      if (!p[lp].node && lp) /* found a solution */
        {
/*	    if (nsol++) return 0; */
	    nsol++;
	    if (memcmp(bound,sol,lsol*sizeof(int))) return 0;
	    goto getnext;
	}
      else /* forward ! */
	{
	  if (p[lp].node<=0) 
	           { offset=0; nx=e->n; x=e->a; pb=&amon; sign=1;}
	  else  { offset=e->n; nx=e->m; x=e->b; pb=&bmon; sign=-1;}
	  for (w=*pb; w<nx && 
			  (mark[p[lp].node+sign*x[w]] ||
			  sol[offset+w]>=bound[offset+w]); w++);
	  if (w<nx) /* success */
	    {
	      p[lp].edge=w;
	      sol[offset+w]++;
	      *pb=w;
	      p[lp+1].node=p[lp].node+sign*x[w];
	      mark[p[++lp].node]=1;
	    }
	  else 
	  getnext:
	    {
	      while (lp>=1)
		{
		  mark[p[lp].node]=0;
		  if (p[lp-1].node<=0) 
		    { offset=0; nx=e->n; x=e->a; pb=&amon; sign=1;}
		  else  { offset=e->n; nx=e->m; x=e->b; pb=&bmon; sign=-1;}
		  sol[offset+p[lp-1].edge]--;
		  for (w=p[lp-1].edge+1; 
		       w<nx && (mark[p[lp-1].node+sign*x[w]] ||
				sol[offset+w]>=bound[offset+w]); w++);
		  if (w<nx) /* success */
		    {
		      p[lp-1].edge=w;
		      sol[offset+w]++;
		      *pb=w;
		      p[lp].node=p[lp-1].node+sign*x[w];
		      mark[p[lp].node]=1;
		      goto cont;
		    }
		  else 
		    {
		      lp--;
		      /* reset *pb */
		      for (w=lp-1; w>=0 && (sig(p[w].node) != sign); w--);
		      if (w<0) *pb=0;
		      else *pb=p[w].edge;
		    }
		}
	    cont:
            ;
	    }
	}
    }
  while(lp);
  return 1;
}

/* *************************************************************** 
*/
void saves(int savesol, int *solmax, int l, int num, int s[], int **sols)
{
  void oomem(char *);

  if (!(num%1000)) printf("\n %8d",num);
  if (savesol) /* storing the solution */
    {
      if (*solmax<num) /* get another chunk */
	{
	  *solmax+=SOLCHUNK;
	  if ((*sols=(int *)realloc(*sols, l*(*solmax)*sizeof(int)))==NULL)
	    oomem("savesol");
	}
      memcpy(*sols+l*(num-1),s,sizeof(int)*l);
    }
}

/* ***************************************************************
 */
inline int sig(int z)
{ 
    if (z<=0) return 1;
    return -1;
}



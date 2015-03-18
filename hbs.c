/* ************************************************************************
 ** hbs.c              Hilbert bases package version 0.1 (experimental)
 ** 
 ** Copyright (C) 1995 Dmitrii V. Pasechnik
 **                    RIACA, Amsterdam, The Netherlands
 **
 **  computing the Hilbert basis for a system of homogeneous linear
 **  Diophantine equations.
 **
 **  vsols performs one step of the algorithm; 
 **  hbs is the front-end.
* ************************************************************************ */

#include <stdio.h>
#include <string.h>
#include "hb.h"

typedef struct _sol {
  int *s; /* the vector itself */
  struct _sol *next; /* to link them */
  int deg;
  int mask;
} sol;
struct inj {int **p; int c;};
struct oss {int *vec; int sum;};

/* ********************************************************************
   return the number of elts in the Hilbert basis *s 
   of the n x m matrix a; the 4th parameter returns the pointer to the basis
   itself.
*/
int hbs(int n, int m, int *a, int **res)
{
  int i, j, *v, nn;
  sol *root;
  int *tomat(int, int,sol *);
  void delsls(sol **);
  int vsols(int,int *,int,sol **);

  if (NULL==(v=(int *)calloc(m,sizeof(int)))) oomem("hbs");
  for (i=0, root=NULL; i<m; i++) /* getting an identity matrix */
    {
      v[i]=1;
      if (i) v[i-1]=0;
      j=inssol(m,&root,v);
    }
  free(v);
  for (i=0, v=a, nn=m; i<n && (nn=vsols(m,v,nn,&root)); i++, v+=m);
  if (nn) *res=tomat(m,nn,root);
  delsls(&root);
  return nn;
}

int getmask(int n, int *v) /* bitwise encode nonzero positions */
{
  register int r,m;
  v+=n-1;
  r=0; m=1;
  while(n--)
    {
      if (*v--) r+=m;
      m=m<<1;
    }
  return r;
}
/* *******************************************************************
   lexicographic comparision of two vectors a, b
*/
inline int cmpvl(int n, int *a, int sa, int ma, int *b, int sb, int mb)
{
  if (sa != sb) return sa-sb;
  if (ma != mb) return ma-mb;
/*  while (n && *a==*b) {n--; a++; b++; }
  if (n) return *a-*b;
  return 0; */
  return memcmp(a,b,n);
}

/* *******************************************************************
   comparision of two vectors a, b 
   return >=1 if there is i such that a[i]>b[i]
   0 otherwise.
*/
inline int ismin(int n, int *a, int *b, int ma, int mb)
{
  if ((ma & mb) ^ ma) return 1; 
  while (n-- && *a++ <= *b++);
  return ++n;
}

int sumv(int n, int *v)
{ register int s,i;
  for (i=0,s=0; i<n; i++) s+=*(v++);
  return s;
}
/* *******************************************************************
   handling the list of lexicographically ordered solutions;
   given the list and a solution, check if it is minimal;
   if yes, insert it in the list and remove the non-min. 
   (against this new one) solutions.
   v is assumed to be created by malloc.
   return the number of solutions added (maybe <=0 in total )
*/
int inssol(int n, sol **root, int *v)
{
  int min, ctr, deg, nbytes, mask;
  sol *s, *s_prev, *p;
  s=*root; min=1;
  nbytes=n*sizeof(int);
  deg=sumv(n,v);
  mask=getmask(n,v);
  while (s != NULL && min && 
	 cmpvl(nbytes,s->s,s->deg,s->mask,v,deg,mask)<=0) 
    { /* while lexigraphically <= */
      min=ismin(n,s->s,v,s->mask,mask); 
      s_prev=s;
      s=s->next;
    }
  
  if (!min) return 0;
  /* insertion */
  if (NULL==(p=(sol *)malloc(sizeof(sol)))) oomem("inssol(1)");
  if (NULL==(p->s=(int *)malloc(n*sizeof(int)))) oomem("inssol(2)");
  memcpy(p->s,v,n*sizeof(int));
  p->deg=deg;
  p->mask=mask;
  p->next=s;
  if (*root==s) *root=p; /* insert in front */
  else s_prev->next=p;

  /* check the minimality of the rest, and remove if necessary */
  s_prev=p;
  ctr=1;
  while (s != NULL)
    {
      if (!ismin(n,v,s->s,mask,s->mask))
	{
	  s_prev->next=s->next;
	  free(s->s); free(s);
	  s=s_prev->next;
	  ctr--;
	}
      else 
	{
	  s_prev=s;
	  s=s->next;
	}
    }
  return ctr;
}


/* *******************************************************************
   checking the minimality of a partial solution;
   return 1 if minimal, 0 otherwise
   (basically cut from inssol).
   (note the different pointer level of root, since we're not gonna
   change it.)
*/
int ismps(int n, sol *root, int *v)
{
  int min,deg,nbytes,mask;
  sol *s;
  s=root; min=1;
  nbytes=n*sizeof(int);
  deg=sumv(n,v);
  mask=getmask(n,v);
  while (s != NULL && min && cmpvl(nbytes,s->s,s->deg,s->mask,v,deg,mask)<=0) 
    { /* while <= */
      min=ismin(n,s->s,v,s->mask,mask);
      s=s->next;
    }
  
  return min;
}

/* *******************************************************************
   deletion of a list of solutions
*/
void delsls(sol **root)
{
  sol *s, *p;
  s=*root; 
  while (s != NULL) 
    { 
      p=s;
      s=s->next; 
      free(p->s); free(p);
    }
  *root=NULL;
}
/* *******************************************************************
   making an array of pointers to solutions
*/
int **inds(int num, sol *s)
{
  int **w,i;
  if (NULL==(w=(int **)calloc(num, sizeof(int *)))) oomem("inds");
  for (i=0; i<num; i++)
    {
      w[i]=s->s;
      s=s->next;
    }
  return w;
}
/* *******************************************************************
   converting list into matrix
*/
int *tomat(int n, int num, sol *s)
{
  int nbytes, i, *w, *p;
  nbytes=n*sizeof(int);
  if (NULL==(w=(int *)malloc(nbytes*num))) oomem("tomat");
  for (i=0, p=w; i<num; i++)
    {
      memcpy(p,s->s,nbytes);
      s=s->next;
      p+=n;
    }
  return w;
}

/* *******************************************************************
   computing a*v
*/
int *mul1(int n, int len, sol *s, int *v)
{
  int i, j, *u, *p;

  if (NULL==(u=(int *)calloc(len, sizeof(int)))) oomem("mul1");
  for (i=0; i<len; i++, s=s->next) 
    for (j=0, p=s->s; j<n;) u[i]+=(*(p++))*v[j++]; 

  return u;
}

/* ********************************************************************
*/  
int vsols(int n, int v[],  /* the current equation                 */
	  int nsol0,       /* # of elts of the prev. H. basis      */
	  sol **sols      /* ptr to H. basis for the previous eqs  */
	  ) /* return size of the new H. basis, putting in in **sols */
{
  int *ord, **x, *w, j, **ind, n_inj, tot, i, *u, *p, *s, nsol1;
  eq *u_inj;
  extern int prtlev;
  struct inj *a;
  sol *root, *r;
  int nbytes;
  struct oss *os;
  int comp_oss(struct oss *, struct oss *);
  nbytes=n*sizeof(int);

  u=mul1(n,nsol0,*sols,v); /* computing (*sols)*v  */

  if (prtlev>1) { printf("\n the rewritten equation");
		  outmat(stdout,1,nsol0,u);}

  u_inj=injcom(nsol0, u, 0); /* computing inj. companion of u */

  nsol1=hb(u_inj,&s,1); /* finding the H. basis for u_inj */
  if (prtlev>1) 
    printf("\n found basis of size %d for the injective companion\n",nsol1);
  /* rewriting the H.basis for u_inj in the variables of u */

  root=NULL; /* set the beginning of the new solutions */

  /* the solutions corresponding to 0 entries in u are also OK for 
     the current equation; copy them. */
  for (i=0, r=*sols, tot=0; i<nsol0; i++, r=r->next) 
    if (!u[i]) tot+=inssol(n,&root,r->s);

  if (!nsol1) /* no new solutions except just copied */
    {
      kieq(u_inj);
      free(u);
      delsls(sols);
      *sols=root;
      return tot;
    }

  n_inj=u_inj->n+u_inj->m;

  /* sort the solutions of the inj.comp. in order of incr. degree */
  if (NULL==(os=(struct oss *)calloc(nsol1,sizeof(struct oss)))) 
    oomem("vsols(os)");
  for (i=0, p=s; i<nsol1; i++)
    {
      os[i].vec=p;
      for (j=0; j<n_inj; j++) os[i].sum+=*(p++);
    }
  qsort(os,nsol1,sizeof(struct oss),comp_oss);
  /* prepare pointers etc. */
  /* a[0].p contains the pointer to the beginning of the array of size nsol0 
     containing the positions of u_inj->a(or b)[i] in u;
     approriate places in this array are pointed by a[i].p;
     numbers of occurences of u_inj->a(b)[i] in u are in a[i].c; */

  if (NULL==(a=(struct inj *)calloc(n_inj,sizeof(struct inj)))) 
    oomem("vsols(a)");
  if (NULL==(a[0].p=(int **)calloc(nsol0,sizeof(int *))))
    oomem("vsols(a[0].p)");
    
  for (i=0, x=a[0].p; i<u_inj->n; i++)
    for (j=0, a[i].p=x, a[i].c=0; j<nsol0; j++) 
      if (u_inj->a[i]==u[j]) { *(x++)=&u[j]; a[i].c++; }
  
  for (i=0; i<u_inj->m; i++)
    for (j=0,a[i+u_inj->n].p=x, a[i+u_inj->n].c=0; j<nsol0; j++) 
      if (u_inj->b[i]==-u[j]) { *(x++)=&u[j]; a[i+u_inj->n].c++; }
  

  bzero(u,nsol0*sizeof(int));/*use u to store the current (partial) solution */
  if (NULL==(w=(int *)malloc(nbytes))) oomem("vsols(w)");
  if (NULL==(ord=(int *)malloc(n_inj*sizeof(int)))) oomem("vsols(ord)");
  ind=inds(nsol0,*sols);
  /* rewriting the H.basis for u back into the original n vars on v */

  for (i=0 /* , p=s */ /* beginning of the sols of inj.comp.*/; 
       i<nsol1; i++ /* , p+=n_inj */)
    {
      tot=addsls(n,u,n_inj,&root,a,ind,os[i].vec,w,tot,ord);
      if (prtlev && i && !(i%100)) printf("\n %6d %6d %6d",i,tot,os[i].sum);
    }

  free(a[0].p); free(a); 
  free(ord); free(os);
  kieq(u_inj);
  free(u);free(w);free(ind);
  delsls(sols);
  *sols=root;
  return tot;
}

/* ***********************************************************************
   handle the solutions corresponding to a solution p of inj. comp.
*/
int addsls(int n, int *u, int n_inj, sol **root, struct inj a[], 
	   int *ind[], int p[], int w[], int tot, int *ord)
{
  void mullev(int, int **, int *, int, struct inj *, int *, int *, int *);
  int lev;
  /* start with the original ordering of p */
  for (lev=0; lev<n_inj; lev++) ord[lev]=lev; 

  for (lev=0; lev<n_inj-1; lev++) /* sort in the desc. order with 0s first */
    { int j, max, maxj;
      for (max=p[ord[lev]], j=lev, maxj=lev; max && j<n_inj; j++)
	if (!p[ord[j]] || p[ord[j]]>max) 
	      {
		maxj=j;
		max=p[ord[j]];
	      }
      j=ord[lev];
      ord[lev]=ord[maxj];
      ord[maxj]=j;
    }

  lev=0;
  while(1)
    {
      int nn,ism;
      do
	{
	  while (lev<n_inj-1 && !p[ord[lev]]) lev++; /* skip 0s */
	  fir_x(a[ord[lev]].c,a[ord[lev]].p,p[ord[lev]]);
	  mullev(n,ind,w,lev,a,p,u,ord);
	  nn=1;
	  while (nn && !ismps(n,*root,w))
	    {
	      nn=nex_x(a[ord[lev]].c,a[ord[lev]].p,p[ord[lev]]);
	      if (nn) mullev(n,ind,w,lev,a,p,u,ord);
	    }
	}
      while(nn && ++lev<n_inj);
      if (lev==n_inj) tot+=inssol(n,root,w);
      do
	{
	  do /* skip 0s */
	    if (--lev<0) return tot;
	  while (!p[ord[lev]]); 
	  do
	    {
	      if ((nn=nex_x(a[ord[lev]].c,a[ord[lev]].p,p[ord[lev]])))
		{
		  mullev(n,ind,w,lev,a,p,u,ord);
		  if (lev==n_inj-1) tot+=inssol(n,root,w);
		  else ism=ismps(n,*root,w);
		}
	      else ism=0;
	    }
	  while (nn && (lev==n_inj-1 || !ism));
	}
      while(!ism);
      lev++;
    }
}

/* *******************************************************************
   computing u*sols for the first lev entries in the solution of the inj. comp.
*/
void mullev(int n, int *ind[], int *w, int lev, struct inj *a, int *p, 
	    int *u, int *ord)
{ 
  register int i;
  bzero(w,n*sizeof(int));
  for (i=0; i<=lev; i++ /* ,a++ */)
    /* if (*p++) */ /* nonzero entry in the solution of the inj. comp. */
    if (p[ord[i]])
      { 
	register int j, jjj, *r, x;
	for (j=0; /* j<a->c */ j<a[ord[i]].c; j++)
	 /* if ((x=*(a->p[j]))) */ /* nonzero entry in the solution itself */
	  if ((x=*(a[ord[i]].p[j])))
	    for (jjj=0, /* r=ind[a->p[j] - u] */ 
		 r=ind[a[ord[i]].p[j] - u]; jjj<n;) w[jjj++]+=*(r++)*x;
      }
}
		  
int comp_oss(struct oss *a, struct oss *b)
{
  if (a->sum > b->sum) return 1;
  if (a->sum < b->sum) return -1;
  return 0;
}

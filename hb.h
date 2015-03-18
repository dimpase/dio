/* ************************************************************************
 ** hb.h               Hilbert bases package version 0.1 (experimental)
 ** 
 ** Copyright (C) 1995 Dmitrii V. Pasechnik
 **                    RIACA, Amsterdam, The Netherlands
 **
 **                    headers file 
* ************************************************************************ */

typedef struct _eq { /* $\sum_1^n a_i + \sum_1^m b_i = c$ */
  int n;             /* $a_i, b_i >0$, $c>=0$             */
  int *a;
  int m;
  int *b;
  int c;
} eq;


/* each noneg. sol. of eq is given by $v+\sum_{i=1}^{ns_hom} t_i u_i$,
   where $t_i>=0$, $u_i\in s_hom$, $v\in s$ if $c\ne 0$, otherwise $v=0$ */

#define SOLCHUNK  1000 /* (re)allocate memory for solutions in these
			   bits (the number stands for # of solutions */

int hb(eq *, int **, int), nf(int, int), nexpa(int, int *, int);
eq *injcom(int, int *, int); 
void oomem(char *), firpa(int, int *, int), kieq(eq *);
int choose(int, int);
int hbs(int, int, int *, int **);
int *getmat(FILE *, int, int);
void outmat(FILE *, int, int, int *);
int prtlev;

/* Ising model in two dimensions, standard single-spin-flip algorithm */
/* using Metropolis flip rate min[1, exp(-Delta E/kT)]                */
/*                                Jian-Sheng Wang, November 1993      */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
                              /* macro definitions */
#define  L  16                /* lattice linear size */
#define  N  (L*L)             /* total number of spins */
#define  Z  4                 /* coordination number = 2*d */
#define  MCTOT 10000           /* total Monte Carlo steps */
#define  MCDIS 500            /* steps discarded in the beginning */
                              /* global variables */
int s[N];                     /* spin +1 or -1 */
double  T;            /* temperature */

			      /* funcition prototypes */
void neighbor(int i, int nn[ ]);
void monte_carlo();
void energy(double *e, double *e_sqr);

/*  The main program  */

void main()
{
   int i, mc;
   double e = 0, e_sqr = 0;
   double e_ave, e_sqr_ave;
   double e_var;
   double C;
   FILE *file;
   file = fopen("result2.txt","w");
   for(T=1; T<=4; T+=0.2)
   {
      for (i = 0; i < N; ++i)     /* initialize, all spin up */
         s[i] = 1;

      for(mc = 0; mc < MCTOT; ++ mc) {
         monte_carlo();
         if( mc >= MCDIS) 
            energy(&e, &e_sqr);
      }
      e_ave = e/(MCTOT-MCDIS)/N;
      e_sqr_ave = e_sqr/(MCTOT-MCDIS)/N;
      e_var = e_sqr_ave - e_ave*e_ave;
      C = e_var/(T*T);
      fprintf(file, "%lf, %lf, %lf\n", T, C, e_ave);
   }
   fclose(file);
}

/* This function monte_carlo performs one Monte Carlo step by trying to
flip spins L^2 times.  It picks a site at random and trys to flip it.
The flip is actually performed if energy change is negative, or if the 
random number is less than exp(-delta E/kT). N is the total number of spin,
a macro definition.  Spin s[], temperature T are passed globally. */

void monte_carlo()
{
   int i, j, k, e;                     /* i is the center site */
   int nn[Z];                          /* the name neighbors */

   for(k = 0; k < N; ++k) {
      i = drand48() * (double) N;      /* pick site at random */
      neighbor(i, nn);                 /* find neighbors of site i */
      for(e = 0, j = 0; j < Z; ++j)    /* go over the neighbors */
         e += s[nn[j]];                /* sum of the neighbor spins */
      e *= 2 * s[i];                   /* 2 times the center spin */
      if (e <= 0)                      /* when energy change is less */
         s[i] = - s[i];                /* than zero, spin is flipped */
      else if (drand48() < exp(-e/T))  /* othewise, it is flipped */
         s[i] = - s[i];                /* with probability less one */
   }
}

/* Neighbor returns in the array nn[ ] the neighbor sites of i.  The sites
are labelled sequentially, starting from 0.  It works for any hypercubic
lattice.  Z (=2*D) is the coordination number, passed as a macro
defintion.  L is linear size, also passed as a macro definition. */

void neighbor(int i, int nn[ ])
{
   int j, r, p, q;

   r = i;
   p = 1 - L;
   q = 1;

   for(j = 0; j < Z; j += 2) {
      nn[j] = (r + 1) % L == 0 ? i + p : i + q;
      nn[j+1]     = r % L == 0 ? i - p : i - q;
      r = r/L;
      p *= L;
      q *= L;
   }
}


/* This function energy add the energy of the current configuraion s[] 
to the argument e.  s[] and size information Z and N are passed globally */

void energy(double *e, double *e_sqr)
{
   int i, j, ie = 0; 
   int nn[Z];
   
   for(i = 0; i < N; ++i) {
      neighbor(i, nn);
      for(j = 0; j < Z; j += 2)    /* look at positive direction only */
         ie += s[i]*s[nn[j]];
   }
   assert(ie <= 2*N && ie >= -2*N);
   *e += ie;
   *e_sqr += ie*ie;
}

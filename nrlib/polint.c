//polint.c
#include <math.h>
#include "nrutil.h"

void polint(float xa[], float ya[], int n, float x, float *y, float *dy)

/*Given arrays xa[1..n] and ya[1..n], and given a value x, 
this routine returns a value y, and an error estimate dy.
If P(x) is the polynomial of degree N-1 such that P(xai) = yai,
then the returned value y=P(x).*/
{
	int i, m, ns=1;
	float den, dif, dift, ho, hp, w;
	float *c, *d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);

	for (i=1; i<=n; i++){       //Here we find the index ns of the closest table entry
		if( (dift=fabs(x-xa[i])) < dif){
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}

	*y=ya[ns--]; // y=ya[ns]; ns-=1; // This is the initial approximation to y.
	for (m=1; m<n; m++){        // For each column of the tableau,
		for (i=1; i<=n-m; i++){ // we looop over the current c's and d's and update them.
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			// This error can occur only if two input xa's (to within roundoff) are identical
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += ( *dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
		// if 2*ns < (n-m), that means upper path is better, vice versa.
	}
}


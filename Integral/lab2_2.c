#include <stdio.h>
#include <math.h>
#include "../nrlib/nrutil.c"
#include "../nrlib/nrutil.h"
#define EPS 1.0e-13
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5 //K for Roberg Integration

int main()
{
	long double test_func(long double x);
	long double qromb(long double (*func)(long double), long double a, long double b);
	long double s;
 	s = qromb(test_func, 0, 2);
	printf("I = %.16lf\n",s);
}

long double test_func(long double x)
{
	return pow(x,4)*log(x+sqrt(x*x+1));
}

void polint(long double xa[], long double ya[], int n, long double x, long double *y, long double *dy)

/*Given arrays xa[1..n] and ya[1..n], and given a value x, 
this routine returns a value y, and an error estimate dy.
If P(x) is the polynomial of degree N-1 such that P(xai) = yai,
then the returned value y=P(x).*/
{
	int i, m, ns=1;
	long double den, dif, dift, ho, hp, w;
	long double *c, *d;

	dif=fabsl(x-xa[1]);
	c=dvector(1,n);
	d=dvector(1,n);

	for (i=1; i<=n; i++){       //Here we find the index ns of the closest table entry
		if( (dift=fabsl(x-xa[i])) < dif){
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
	free_dvector(c, 1, 0);
	free_dvector(d, 1, 0);
}


// Basic trapzodial integration
#define FUNC(x) ((*func)(x))
long double trapzd(long double (*func)(long double), long double a, long double b, int n)
{
	long double x,tnm,sum,del;
	static long double s; //it will be initialize as 0 automatically
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

//Romberg Integration
long double qromb(long double (*func)(long double), long double a, long double b)
{
	void polint(long double xa[], long double ya[], int n, long double x, long double *y, long double *dy);	
	long double trapzd(long double (*func)(long double), long double a, long double b, int n);
	void nrerror(char error_text[]);
	long double ss, dss;
	long double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1; j<=JMAX; j++){
		s[j]=trapzd(func, a, b, j);
		if (j>=K) {
			polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
			if (fabsl(dss)<=EPS*fabsl(ss)) return ss;
		}
		h[j+1]=0.25*h[j]; // Here h is actually h^2 !!!!!
		/*This is a key step: The factor is 0.25 
		  even though the stepsize is decreased by only
		 0.5. This makes the extrapolation a polynomial 
		 in h^2 not just a polynomial in h*/ 
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;
}

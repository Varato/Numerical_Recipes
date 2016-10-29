#include <math.h>
#include "nrutil.h"
#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5 //K for Roberg Integration

// Basic trapzodial integration
#define FUNC(x) ((*func)(x))
double trapzd(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del;
	static double s; //it will be initialize as 0 automatically
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

// Trapzodial Integration
double qtrap(double (*func)(double), double a, double b)
{

	double trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	int j;
	double s, olds=0.0;

	for (j=1; j<=JMAX; j++){
		s = trapzd(func, a, b, j);
		if (j>5)
			if (abs(s-olds) < EPS*abs(olds) ||
				(s==0.0 && olds ==0.0)) return s;
		olds = s;
	}
	nrerror("Too many steps in routine qtrap");
	return 0.0;
}

// Simpson Integration
double qsimp(double (*func)(double), double a, double b)
{
	double trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	int j;
	double s, st, ost=0.0, os=0.0;

	for (j=1; j<=JMAX; j++){
		st = trapzd(func, a, b, j);
		s = (4.0*st-ost)/3.0; // Notice S2N and SN can be
		if (j>5)              // obtained by any two consecutive calls of trapzd()
			if (abs(s-os) < EPS*abs(os) || 
					(s == 0.0 && os == 0.0)) return s;
		os = s;
		ost = st;
	}
	nrerror("Too many steps in routine qsimp");
	return 0.0;
}


//Romberg Integration
double qromb(double (*func)(double), double a, double b)
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);	
	double trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	double ss, dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1; j<=JMAX; j++){
		s[j]=trapzd(func, a, b, j);
		if (j>=K) { // do the interpolation at least 5 rounds completed.
			polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
			if (abs(dss)<=EPS*abs(ss)) return ss;
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



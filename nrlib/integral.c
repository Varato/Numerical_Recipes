#include <math.h>
#include "nrutil.h"
#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5 //K for Roberg Integration

// Basic trapzodial integration
#define FUNC(x) ((*func)(x))
float trapzd(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del;
	static float s; //it will be initialize as 0 automatically
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
float qtrap(float (*func)(float), float a, float b)
{

	float trapzd(float (*func)(float), float a, float b, int n);
	void nrerror(char error_text[]);
	int j;
	float s, olds=0.0;

	for (j=1; j<=JMAX; j++){
		s = trapzd(func, a, b, j);
		if (j>5)
			if (fabs(s-olds) < EPS*fabs(olds) ||
				(s==0.0 && olds ==0.0)) return s;
		olds = s;
	}
	nrerror("Too many steps in routine qtrap");
	return 0.0;
}

// Simpson Integration
float qsimp(float (*func)(float), float a, float b)
{
	float trapzd(float (*func)(float), float a, float b, int n);
	void nrerror(char error_text[]);
	int j;
	float s, st, ost=0.0, os=0.0;

	for (j=1; j<=JMAX; j++){
		st = trapzd(func, a, b, j);
		s = (4.0*st-ost)/3.0; // Notice S2N and SN can be
		if (j>5)              // obtained by any two consecutive calls of trapzd()
			if (fabs(s-os) < EPS*fabs(os) || 
					(s == 0.0 && os == 0.0)) return s;
		os = s;
		ost = st;
	}
	nrerror("Too many steps in routine qsimp");
	return 0.0;
}


//Romberg Integration
float qromb(float (*func)(float), float a, float b)
{
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);	
	float trapzd(float (*func)(float), float a, float b, int n);
	void nrerror(char error_text[]);
	float ss, dss;
	float s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1; j<=JMAX; j++){
		s[j]=trapzd(func, a, b, j);
		if (j>=K) {
			polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
			if (fabs(dss)<=EPS*fabs(ss)) return ss;
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



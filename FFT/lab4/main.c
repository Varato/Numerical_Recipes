#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define x_step 0.01
#define p_step 
void show_cnum(double _Complex a);
int main()
{
	// fftw_complex I, q;
	// I[0]=0;
	// I[1]=1;
	show_cnum(I*I);
}

void show_cnum(double _Complex a)
{
	printf("%f + i%f\n", creal(a), cimag(a));
}
double wave_packet(x)
{
	return 1.;
}
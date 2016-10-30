#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define PI 3.1415926

int main()
{
	int N = 8, i, j;
	double xx[8], kk[8];	
	double complex x[8];
	double complex k[8];
	double h=PI/4;

	fftw_plan p, ip;

	p = fftw_plan_dft_1d(N, x, k, FFTW_FORWARD, FFTW_ESTIMATE);
	ip = fftw_plan_dft_1d(N, k, x, FFTW_BACKWARD, FFTW_ESTIMATE);
	for(i=0; i<8; i++){
		xx[i]=i*h;
		x[i]=cos(xx[i])+I*sin(xx[i]);
		if(i<N/2) kk[i] = i*2*PI/h/N;
		else kk[i] = (i-N)*2*PI/h/N; 	
	}
	for (i=0; i<N; i++){
		printf("%lf+%lfi ", creal(x[i]), cimag(x[i]));
	}
	fftw_execute(p);
	printf("\n");
	for (i=0; i<N; i++){
		printf("%lf+%lfi ", creal(k[i]), cimag(k[i]));
	}
	printf("\n");
	fftw_execute(ip);
	for (i=0; i<N; i++){
		printf("%lf+%lfi ", creal(x[i]), cimag(x[i]));
	}
	printf("\n");
	fftw_destroy_plan(p);
	fftw_destroy_plan(ip);

}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define PI 3.141592653589793
#define x_step 0.001
#define dt 0.1
#define L 40
#define sigma 1
#define N 40000
#define a 7.0710678118654755
#define x0 -6
#define p0 10
#define V0 1

//declaration of functions
double complex Gaussian_wp(double);
void normalize();
void initialize();
void write_file();
void evolve();

//global variables
double complex V_exp[N];
double complex T_exp[N];
double complex wave_packet[N];   // x representation of wave packet
double complex p_wave_packet[N]; // p representation of wave packet
fftw_plan fft, ifft;             // fft plan for fftw 

int main()
{
	FILE *f;
	int i;
	f=fopen("result", "w");
	initialize();
	fft = fftw_plan_dft_1d(N, wave_packet, p_wave_packet, FFTW_FORWARD, FFTW_ESTIMATE);
	ifft = fftw_plan_dft_1d(N, p_wave_packet, wave_packet, FFTW_BACKWARD, FFTW_ESTIMATE);
	write_file(f);
	for(i=0; i<=300; i++){
		printf("evolved step: %d\n", i);
		evolve();
		if(i%5==0)
			write_file(f);
	}
	fclose(f);
}

void write_file(FILE *f)
{
	int i;
	for(i=0; i<N; i++){
		fprintf(f, "%lf ", creal(wave_packet[i]));
		if(i==N-1) fprintf(f, "\n");
	}
	// for(i=0; i<N; i++){
	// 	fprintf(f, "%lf ", cimag(wave_packet[i]));
	// 	if(i==N-1) fprintf(f, "\n");
	// }
}
double complex Gaussian_wp(double x)
{
	return cexp(-(x-x0)*(x-x0)/(2*sigma*sigma)+I*p0*(x-x0)); 
}

void normalize()
// This function nomalizes the wave packet numerically
{
	double ww[N];
	double prod=0;
	int i;
	for(i=0; i<N; i++){
		ww[i] = wave_packet[i]*conj(wave_packet[i]);
	}
	//trapezoidal integral
	for(i=0; i<N-1; i++){
		prod += (ww[i]+ww[i+1])*x_step/2.0;
	}
	for(i=0; i<N; i++){
		wave_packet[i] = wave_packet[i]/prod;
	}

}

void initialize()
{
	int i;
	double x, p;
	for(i=0; i<N; i++){
		x = -L/2.0 + x_step*(i+1);
		if(0.<=x && x<=a) V_exp[i]=cexp(-I*V0*dt);
		else V_exp[i]=1;

		wave_packet[i] = Gaussian_wp(x);
	}

	for(i=0; i<N; i++){
		if(i<N/2) p = i*2.0*PI/x_step/N;
		else      p = (i-N)*2.0*PI/x_step/N;
		T_exp[i] = cexp(-I*p*p*dt/2.0);
	}

}

void evolve()
{
	int i;
	for(i=0; i<N; i++){
		wave_packet[i] *= V_exp[i];
	}
	fftw_execute(fft);
	for(i=0; i<N; i++){
		p_wave_packet[i] *= T_exp[i];
	}
	fftw_execute(ifft);
	normalize();
}












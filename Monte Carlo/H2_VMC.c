#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../nrlib/nrutil.h"
#include "../nrlib/nrutil.c"

// Global varables
int MCSteps = 5000;
int n_walkers = 150;
float step_size = 1.;
float E_sum, Esqrd_sum, nAccept;
float E_ave, E_var;
float **r10, **r20;

float R_norm = 1.398; // interval between protons
float R[3];

int main()
{
	float a;
	srand(2016);
	printf("%d\n%d", rand(),RAND_MAX);
}

float rand_uniform(float a, float b)
{
	if(b<=a) return 1;
	float rand_float = rand()/float(RAND_MAX);
	float rand_u;
	rand_u = rand_float*(b-a)+a;
	return rand_u;
}
float * rand_vector()
{
	float vector[3];
	for (int dd=0; dd<3; dd++)
	{
		vector[dd] = rand_uniform(-1, 1);
	}
	return vector;

}
float vector_norm(float *r)
{
	float sum=0;
	for(int dd=0; dd<3; dd++)
	{
		sum += r[i]*r[i]
	}
	return sqrtf(norm)
}
void initialize()
{
	float tmp1[3], tmp2[3];
	r10 = (float **) malloc(n_walkers*sizeof(float *));
	r20 = (float **) malloc(n_walkers*sizeof(float *));
	for (int i=0; i<3; i++)
	{
		r10[i] = (float *) malloc(3*sizeof(float));
		r20[i] = (float *) malloc(3*sizeof(float));
	}

	for(int i=0; i<n_walkers; i++)
	{
		for(int dd=0; dd<3; dd++)
		{
			tmp1[dd] = rand_uniform(-0.5, 0.5);
			tmp2[dd] = rand_uniform(-0.5, 0.5);
		}
		r10[i] = tmp1;
		r20[i] = tmp2;
	}
}

void zeroAccumulate()
{
	E_sum = 0;
	Esqrd_sum = 0;
}

float trial_wave_func(float *r1, float *r2)
{

}

float local_energy(float *r1, float *r2)
{

}

void Metropolis_step(int walker)
{

}

void oneMento_Carlo_step()
{

}














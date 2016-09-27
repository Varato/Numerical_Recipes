#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// Global varables
int MCSteps = 10000;
int n_walkers = 150;    //number of walkers
int nAccept;
double step_size = 1.;
double E_sum, Esqrd_sum;
double E_ave, E_var;
double **r10, **r20;

double R_norm; // interval between protons
double R[3];

int main()
{
	srand(1946); // set random seed
    
    // declaration of functions
	double rand_uniform(double a, double b);
	void initialize();
	void zeroAccumulate();
	void zeroAccumulate();
	double trial_wave_func(double *r1, double *r2);
	double local_energy(double *r1, double *r2);
	void Metropolis_step(int walker);
	void oneMonte_Carlo_step();
	void runMonte_Carlo();

	FILE *file;
	char *file_name;
	file_name = "result.txt";
	double Error;
	double rAB = 0.2;   // initial rAB
	int numPoints = 30; // going to compute 30 points
	double dr = 0.05;   // the interval between those points
	file = fopen(file_name,"w");
	for(int p=0; p<numPoints; p++)
	{
		printf("Current rAB = %f, %d/%d \n", rAB, p+1, numPoints);
		R_norm = rAB;
		initialize();
		runMonte_Carlo();
		Error = sqrt(E_var);
	        //write the result into a file
		fprintf(file, "%f, %f, %f\n", rAB, E_ave, Error);
		rAB += dr;
	printf("Acceptance ratio: %f\n", 1.0*nAccept/n_walkers/MCSteps);
	printf("Result:\nE_0 = %f (e/a0), E_var %f \n", E_ave, E_var);
	printf("---------------------------------------------------\n");
	}
	fclose(file);

}

double rand_uniform(double a, double b)
// this function gives a uniform random number between a and b
{
	if(b<=a) return 1;
	// A random float between 0 and 1
	double rand_float = rand()/(double)(RAND_MAX); 
	double rand_u;
	rand_u = rand_float*(b-a)+a;
	return rand_u;
}

void initialize()
{
	printf("initialize...\n");
	R[0] = 0; R[1]=0; R[2]=R_norm;
	r10 = (double **) malloc(n_walkers*sizeof(double *));
	r20 = (double **) malloc(n_walkers*sizeof(double *));
	for (int i=0; i<n_walkers; i++)
	{
		r10[i] = (double *) malloc(3*sizeof(double));
		r20[i] = (double *) malloc(3*sizeof(double));
		for (int d=0; d<3; d++)
		{
			r10[i][d] = rand_uniform(-0.5, 0.5);
			r20[i][d] = rand_uniform(-0.5, 0.5);
		}
	}
}

void zeroAccumulate()
{
	E_sum = 0;
	Esqrd_sum = 0;
}

double trial_wave_func(double *r1, double *r2)
{
	double rA1, rA2, rB1, rB2;
	double psi;
	rA1 = rA2 = rB1 = rB2 = 0;
	for(int d=0; d<3; d++) rA1 += pow(r1[d]-R[d],2);
	for(int d=0; d<3; d++) rA2 += pow(r2[d]-R[d],2);
	for(int d=0; d<3; d++) rB1 += pow(r1[d],2);
	for(int d=0; d<3; d++) rB2 += pow(r2[d],2);
	rA1 = sqrt(rA1);
	rA2 = sqrt(rA2);
	rB1 = sqrt(rB1);
	rB2 = sqrt(rB2);
	psi = exp(-rA1-rB2)+exp(-rA2-rB1);
	return psi;
}

double local_energy(double *r1, double *r2)
{

	double rA1, rA2, rB1, rB2, r12, rAB;
	rA1 = rA2 = rB1 = rB2 = r12 = 0;
	double T_up, T, V, E_loc;
	for(int d=0; d<3; d++) 
	{
		rA1 += pow(r1[d]-R[d],2);
		rA2 += pow(r2[d]-R[d],2);
		rB1 += pow(r1[d],2);
		rB2 += pow(r2[d],2);
		r12 += pow(r1[d]-r2[d],2);
	}
	rA1 = sqrt(rA1);
	rA2 = sqrt(rA2);
	rB1 = sqrt(rB1);
	rB2 = sqrt(rB2);
	r12 = sqrt(r12);
	rAB = R_norm;
	T_up = exp(-rA1-rB2)*(1./rA1+1./rB2-1) 
		   + exp(-rA2-rB1)*(1./rB1+1./rA2-1);
	T = T_up/trial_wave_func(r1, r2);
	V = -1./rA1-1./rB1-1./rA2-1./rB2+1./rAB+1./r12;
	E_loc = T+V;
	return E_loc;
}

void Metropolis_step(int walker)
{
	double trial_r1[3], trial_r2[3];
	double w0, w1;             // weight
	double p;                  // acceptence probability
	double E_loc;
	for(int d=0; d<3; d++)
	{
		trial_r1[d] = r10[walker][d] + step_size*rand_uniform(-1,1);
		trial_r2[d] = r20[walker][d] + step_size*rand_uniform(-1,1);
	}
	w0 = trial_wave_func(r10[walker], r20[walker]);
	w1 = trial_wave_func(trial_r1, trial_r2);
	p = (w1*w1)/(w0*w0);

	if(rand()/(double)(RAND_MAX) < p)
	{
		for(int d=0; d<3; d++)
		{
			r10[walker][d] = trial_r1[d];
			r20[walker][d] = trial_r2[d];
		}
		nAccept += 1;
	}
	E_loc = local_energy(r10[walker],r20[walker]);
	E_sum += E_loc;
	Esqrd_sum += E_loc*E_loc;
}

void oneMonte_Carlo_step()
{
	for(int walker=0; walker<n_walkers; walker++)
		Metropolis_step(walker);
}

void runMonte_Carlo()
{
	/* To ensure the acceptance ration is about 0.5,
       use 20% of MCSteps to do thermalization*/
	int thermSteps = (int)(0.2*MCSteps);
	int adjust_interval = (int)(0.1*thermSteps)+1;
	nAccept = 0;

	// Do thermalization to optimize the step_size
	printf("Thermalizing...\n");
	for(int i=0; i<thermSteps; i++)
	{
		oneMonte_Carlo_step();
		if( (i+1)%adjust_interval == 0 )
		{
			step_size*=nAccept/(0.5*n_walkers*adjust_interval);
			nAccept=0;
		}
	}
	printf("Adjusted step_size = %f\n", step_size);
	
	//After step_size-adjusting
	zeroAccumulate();
	nAccept=0;
	printf("Formal simulation:\n");
	for(int i=0; i<MCSteps; i++)
	{
		oneMonte_Carlo_step();
	}
	E_ave = E_sum/(double)(n_walkers)/(double)(MCSteps);
	E_var = Esqrd_sum/(double)(n_walkers)/(double)(MCSteps) 
	        - E_ave*E_ave;
}












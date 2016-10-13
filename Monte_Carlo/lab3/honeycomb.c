#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define Z  3
#define x  32
#define y  32 // y has to be even
#define N x*y
#define MCTOT  1000000
#define MCDIS  1000

int s[N+1]; // index used starts from 1
double T;   // temperature in unit of kT/J

void Monte_Carlo();
void initialize(double *e, double *e_sqr);
int neighbor(int i, int nn[Z]);
void energy(double *e, double *e_sqr);
void run_Monte_Carlo(double T1, double T2, double dT, FILE *file);

int main()
{

	FILE *file;
	char *file_name;
	file_name = "result_Lave1024.txt";
	file = fopen(file_name, "w");
	run_Monte_Carlo(0.5, 1.3, 0.2, file);
	run_Monte_Carlo(1.4, 1.6, 0.05, file);
	run_Monte_Carlo(1.8, 3, 0.2, file);
	fclose(file);
}
void run_Monte_Carlo(double T1, double T2, double dT, FILE *file)
{
	int mc;
	int n, n_ave = 10; //
	double c=0, C_ave; // capacity in unit of k
	double C_sum, C_sqr_sum, C_var;  
	double e;
	double e_sqr;
	double e_ave, e_sqr_ave;
	double e_var;	
	double t;

	for(t=T1; t<=T2; t += dT)
	{
		T=t;
		printf("Current T: %.3lf\n", T);
		C_sum = 0;
		C_sqr_sum = 0;
		for(n=0; n<n_ave; n++)
		{
			initialize(&e, &e_sqr); //set s[N] and zero e and e_sqr
			for(mc=0; mc<MCTOT; mc++)
			{
				Monte_Carlo();
				if(mc >= MCDIS) energy(&e, &e_sqr);
			}
			e_ave = e/((double)(MCTOT-MCDIS));
			e_sqr_ave = e_sqr/((double)(MCTOT-MCDIS));
			e_var = e_sqr_ave - e_ave*e_ave;
			c = e_var/(T*T);
			C_sum += c;
			C_sqr_sum += c*c;		
		}
		C_ave = C_sum/((double)n_ave);
		C_var = C_sqr_sum/((double)n_ave) - C_ave*C_ave;
		fprintf(file, "%lf, %lf, %lf, %lf, %lf\n", T, C_ave, sqrt(C_var), e_ave, sqrt(e_var));
	}
}

void initialize(double *e, double *e_sqr)
{
	for(int i=1; i<=N; i++)
		s[i]=1;
	*e = 0;
	*e_sqr = 0;

}
int neighbor(int i, int nn[Z])
{
	assert(1 <= i && i <= x*y );
	// upper nn
	if(i%x==1) nn[0] = i+x-1;
	else       nn[0] = i-1;
	// down nn
	if(i%x==0) nn[1] = i-x+1;
	else       nn[1] = i+1;
	// horizontal nn
	if( ( (i-1)%(2*x)<=(x-1) && (((i-1)%x)%2==0) ) ||
	    ( (i-1)%(2*x)>=x && (((i-1)%x)%2==1) ) )
	{
		if(x*y-x+1 <= i && i<= x*y) nn[2] = i-x*(y-1);
		else                        nn[2]=i+x;
	}	
	else
	{
		if(1 <= i && i<=x) nn[2] = i+x*(y-1);
		else               nn[2] = i-x;
	}
	return 0;
}
void Monte_Carlo()
{
	int k, i, j, de;
	int nn[Z];

	for(k=0; k<N; k++)
	{
		i = drand48()*((double) N)+1;
		neighbor(i,nn);
		for(de=0,j=0; j<Z; ++j)
			de += s[nn[j]];
		de *= 2*s[i];
		if(de <= 0 || drand48()<exp(-((double)de)/T))
			s[i] = - s[i];
    }
}
void energy(double *e, double *e_sqr)
{
	int i, j, ie = 0; 
	int nn[Z];

	for(i = 1; i <= N; ++i) {
		neighbor(i, nn);
		for(j = 0; j < Z; j++) 
		{
			if(i<nn[j])
				ie += -s[i]*s[nn[j]];
		}   		
	}
	assert((double)ie <= (3.*(double)(N)/2.) && ie >= -3.*(double)(N)/2.);
	*e += ie;
	*e_sqr += (ie)*(ie);
}

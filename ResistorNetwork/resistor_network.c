//author: Xin Chen
//update: 21 Aug 2016

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../nrlib/ludcmp.h"

//void show_2d_matrix(double **p, int n);

int main(int argc, char const *argv[])
{
	int L, n, *indx;
	int i, j;
	double **a, *b, d;
	double I, R; 


	if(argv[1]){
		L = atoi(argv[1]);
		n = (L+1)*(L+1);
		printf("L = %d, matrix_size = %d\n", L, n);
	} 
	else{
		printf("Usage:\n resistorSolver [L]\n");
		return 1;
	}
	indx = malloc((n+1)*sizeof(int));

	// allocate memory for the matrix and initialize it with 0
	a = (double **)malloc((n+1)*sizeof(double *));
	for (i=0; i<=n; i++){
		a[i] = (double *)malloc((n+1)*sizeof(double));
		for (j=0; j<=n; j++)
			a[i][j]=0;
	}
	b = (double *)malloc((n+1)*sizeof(double));

	// set the whole matrix
	for (i=1; i<=n; i++){
		if (i==1){         // left bottom vertice
			a[i][i]=1;
			b[i]=1;
		}
		else if(i==n){     // right top vertice
			a[i][i]=1;
			b[i]=0;
		}
		else if(i==L+1){   // right bottom vertice
			a[i][i]=2;
			a[i][i-1]=-1;
			a[i][i+L+1]=-1;
			b[i]=0;
		}
		else if(i==n-L){   // left top vertice
			a[i][i]=2;
			a[i][i+1]=-1;
			a[i][i-L-1]=-1;
			b[i]=0;
		}
		else if(1<i && i<L+1){ //botom side
			a[i][i]=3;
			a[i][i-1]=-1;
			a[i][i+1]=-1;
			a[i][i+L+1]=-1;
			b[i]=0;

		}
		else if(n-L<i && i<n){ // top side
			a[i][i]=3;
			a[i][i-1]=-1;
			a[i][i+1]=-1;
			a[i][i-L-1]=-1;
			b[i]=0;
		}
		else if(i%(L+1)==1){  //left side
			a[i][i]=3;
			a[i][i+1]=-1;
			a[i][i+L+1]=-1;
			a[i][i-L-1]=-1;
			b[i]=0;
		}
		else if(i%(L+1)==0){ //right side
			a[i][i]=3;
			a[i][i-1]=-1;
			a[i][i+L+1]=-1;
			a[i][i-L-1]=-1;
			b[i]=0;
		}
		else{                // inside
			a[i][i]=4;
			a[i][i+1]=-1;
			a[i][i-1]=-1;
			a[i][i+L+1]=-1;
			a[i][i-L-1]=-1;
			b[i]=0;	
		} //end if
	} //end for
	// test if the a matrix is correct
	// show_2d_matrix(a, n);
	ludcmp(a, n, indx, &d);
	lubksb(a, n, indx, b);
	I = 2*(b[1]-b[2]);
	R = 1./I;
	printf("---solution---\n");
	printf("I = %2.5lf\n", I);
	printf("R = %2.5lf\n", R);
	printf("--------------\n");
	free(a);
	free(b);
	free(indx);

	return 0;
}

void show_2d_matrix(double **p, int n)
{
	int i,j;
	for (i=0; i<=n; i++){
		for (j=0; j<=n; j++){
			printf("%2.0lf ", p[i][j]);
			if (j==n) printf("\n");
		}
	}
}

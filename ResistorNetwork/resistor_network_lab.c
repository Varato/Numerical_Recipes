#include <stdio.h>
#include <stdlib.h>
#include "../nrlib/linear_solve.h"

int main(int argc, char const *argv[])
{
	int L, n, *indx;
	int i, j;
	float **a, *b, d;
	float I, R; 

	if(argv[1]){
		L=atoi(argv[1]);
		n=(L+1)*(L+1);
		printf("L = %d, matrix_size = %d\n", L, n);
	} 
	else{
		printf("Usage:\n resistorSolver [L]\n");
		return 1;
	}

	// allocate memory for the matrix and initialize it with 0
	indx=malloc((n+1)*sizeof(int));
	b=(float *)malloc((n+1)*sizeof(float));
	a=(float **)malloc((n+1)*sizeof(float *));
	for (i=0; i<=n; i++){
		a[i] = (float *)malloc((n+1)*sizeof(float));
		for (j=0; j<=n; j++)
			a[i][j]=0;
	}

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
	}     //end for
	ludcmp(a, n, indx, &d);
	lubksb(a, n, indx, b);
	I=2*(b[1]-b[2]);
	R=1./I;
	printf("---solution---\n");
	printf("I = %2.5f\n", I);
	printf("R = %2.5f\n", R);
	printf("--------------\n");
	for (i=1;i<=n;i++) printf("%.2f ", b[i]);
	printf("\n");
	free(a);free(b);free(indx);
	return 0;
}

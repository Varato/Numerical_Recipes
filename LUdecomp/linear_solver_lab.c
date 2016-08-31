#include <stdio.h>
#include <stdlib.h>
#include "../nrlib/linear_solve.h"

int main()
{
	float **a, *b, d;
	int i, j, n=4, *indx;
	float bb[4]={0,2,3,-10};
	float aa[4][4]={{1, 3, 3, -5},
                        {2,-4, 7, -1},
                        {7,.5, 3, -6},
                        {9,-2, 3,  8}};
	// allocate memory for vectors and matrices
	b=malloc((n+1)*sizeof(float));
	indx=malloc((n+1)*sizeof(int));
	a=malloc((n+1)*sizeof(float *));
	for (i=0; i<=n; i++){
		a[i]=malloc((n+1)*sizeof(float));
		if (i>0){
			for(j=1; j<=n; j++) a[i][j]=aa[i-1][j-1];
			b[i]=bb[i-1];
		}
	}
	// Solve by using the routines
	ludcmp(a, n, indx, &d);
	lubksb(a, n, indx, b);
	printf("Solution:\n[");
	for (i=1; i<=n; i++)
		printf("%.3f ", b[i]);
	printf("]\n");
	free(a);free(b);free(indx);
	return 0;
}

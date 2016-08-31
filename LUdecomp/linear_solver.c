//author: Xin Chen
//update: 21 Aug 2016

#include <stdio.h>
#include <stdlib.h>
#include "../nrlib/linear_solve.h"

int main(int argc, char const *argv[])
{
	float **a, *b;
	float d;
	int n, *indx;
	int i,j;
	const char *file_name;
	FILE *file;

	// Load Matrix information from file
	if (argv[1]) file_name = argv[1];
	else {
		printf("Usage:\n LinearSolver [file_name.txt]\n");
		return 1;
	}

	file = fopen(file_name, "r");
	fscanf(file, "%d", &n);
	b = malloc((n+1)*sizeof(float));
	a = malloc((n+1)*sizeof(float *));
	indx = malloc((n+1)*sizeof(int));
	for (i=0; i<=n; i++){
		a[i] = malloc((n+1)*sizeof(float));
	}

	printf("============================\n");
	printf("Size: %d\n", n);
	printf("The matrix loaded:\n");
	for (i=1; i<=n; i++){
		for (j=1; j<=n; j++){
			if(file!=NULL){
				fscanf(file, "%f", &a[i][j]);
				if (i!=0 && j!=0){
					printf("%5.3f ", a[i][j]);
					if (j==n) printf("\n");
				}
			}
		}
	}
	printf("The b vector loaded:\n[");
	for (j=1; j<=n; j++){
		if(file!=NULL){
		fscanf(file, "%f", &b[j]);
		if(j!=0)
			printf("%5.3f ", b[j]);
		}
	}
	printf("]\n");
	fclose(file);
	printf("============================\n");

	printf("LU decomposition proccessing...\n");
	ludcmp(a, n, indx, &d);

	lubksb(a, n, indx, b);
	printf("Solution:\n");
	printf("[");
	for (i=1; i<=n; i++){
		printf("%.3f ", b[i]);
	}
	printf("]\n");

	return 0;
}

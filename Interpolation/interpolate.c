#include <stdio.h>
#include <math.h>
#include "../nrlib/polint.h"

int main(int argc, char const *argv[])
{
	int i,j,n=3;
	float xa[]={0,-1,0,1};
	float ya[]={0, 1,0,1};
	float x=0.5;
	float y;
	float dy;
	polint(xa, ya, n, x, &y, &dy);
	printf("y = %.3f @ %.3f\n", y, x);
	return 0;
}

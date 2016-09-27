#include <stdio.h>
#include <math.h>
#include "../nrlib/integral.h"
int main()
{
	float test_func(float x);
	float s;
 	s = qromb(test_func, 0, 2);
	printf("I = %f\n",s);
}

float test_func(float x)
{
	return pow(x,4)*logf(x+sqrt(x*x+1));
}



void rk4(float y[], float dydx[], int n, float x, float yout[], void (*derivs)(float, float [], float []))

/*Given values for variables y[1...n] and their derivatives dydx[1...n] known at x, use the 4th Runge-Kutta method to advance the solution over an interval h and return the incremented variables as yout[1...n], which need not to be a distinct array from y. The user supplies the routine derivs(x, y, dydx), which returns derivatives dydx at x.*/

{
	int i;
	float xh, hh, h6, *dym, *dyt, *yt;

	dym = vector(1,n);
	dyt = vector(1,n);
	yt = vector(1,n);
	hh = h*0.5;
	h6 = h/6.0;
	xh = x+hh;

	for (i=1; i<=n; i++) yt[i] = y[i] + hh*dydx[i];
	(*derivs)(xh, yt, dyt);

	for (i=1; i<=n; i++) yt[i] = y[i] + hh*dyt[i];
	(*derivs)(xh, yt, dym);

	for (i=1; i<=n; i++){
		yt[i] = y[i] + h*dym[i];
		dym[i] += dyt[i]; // k2+k3 
	}
	(*derivs)(x+h, yt, dyt); // dyt is reused

	for (i=1; i<=n; i++){
		yout[i] = y[i] + h6*(dydx[i] + dyt[i] + 2.0*dym[i]);
	}
	free_vector(yt,1,n);
	free_vector(dyt,1,n);
	free_vector(dym,1,n);



}

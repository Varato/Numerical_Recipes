#ifndef LUDMCP_H
#define LUDMCP_H
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
#endif //LUDMCP_H
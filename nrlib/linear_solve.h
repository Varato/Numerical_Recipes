#ifndef LUDMCP_H
#define LUDMCP_H
void ludcmp(float **a, int n, int *indx, float *d);
void lubksb(float **a, int n, int *indx, float b[]);
#endif //LUDMCP_H
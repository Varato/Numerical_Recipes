#ifndef INTEGRAL_H
#define INTEGRAL_H
float qtrap(float (*func)(float), float a, float b);
float qsimp(float (*func)(float), float a, float b);
float qromb(float (*func)(float), float a, float b);

#endif //INTEGRAL_H

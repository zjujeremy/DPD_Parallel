#ifndef RANDOMFUN_H
#define RANDOMFUN_H

#include "vector3D.h"

extern long int *iv;
extern long int idum, idum2, iy;
void ranils(unsigned long int iseed);
void randVec3(vector3D<double >& _e);
double ranuls();
double rangls();

#endif
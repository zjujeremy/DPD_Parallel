#include <iostream>
#include <cmath>
#include "vector3D.h"
#include "randomFun.h"
using namespace std;

long int *iv;
long int idum, idum2, iy;

void ranils(unsigned long int iseed){
	const int in = 2147483563;
	const int ik = 40014;
	const int iq = 53668;
	const int ir = 12211;
	const int ntab = 32;
	iv = new long int[ntab];
	if (iv == NULL){
		cout << "allocate space for iv error" << endl;
		throw 1;
	}

	idum = iseed + 123456789;
	idum2 = idum;

	for (int j = ntab + 7; j >= 0; j--){
		int k;
		k = idum / iq;
		idum = ik*(idum - k*iq) - k*ir;
		if (idum < 0)
			idum += in;
		if (j < ntab - 1)
			iv[j] = idum;
	}

	iy = iv[0];
}

void randVec3(vector3D<double> & _e){
	double x, y, r1, r2;
	double s = 2.0;
	while (s > 1){
		r1 = ranuls();
		r2 = ranuls();
		x = 2.*r1 - 1.0;
		y = 2.*r2 - 1.0;
		s = x*x + y*y;
	}

	_e[2] = 1 - 2 * s;
	s = 2 * sqrt(1 - s);
	_e[0] = s*x;
	_e[1] = s*y;
}

double ranuls(){
	const long int in1 = 2147483563;
	const long int ik1 = 40014;
	const long int iq1 = 53668;
	const long int ir1 = 12211;
	const long int in2 = 2147483399;
	const long int ik2 = 40692;
	const long int iq2 = 52774;
	const long int ir2 = 3791;
	const long int ntab = 32;
	const double an = 1. / in1;
	const int inm1 = in1 - 1;
	const double ndiv = 1 + inm1 / ntab;

	//Linear congruential generactor 1
	long double k = idum / iq1;
	idum = ik1*(idum - k*iq1) - k*ir1;
	if (idum < 0)
		idum = idum + in1;

	//Linear congruential generactor 2
	k = idum2 / iq2;
	idum2 = ik2*(idum2 - k*iq2) - k*ir2;
	if (idum2 < 0)
		idum2 = idum2 + in2;

	// Shuffling and subtracting
	int j = 1 + iy / ndiv;
	iy = iv[j-1] - idum2;
	iv[j-1] = idum;
	if (iy < 1)
		iy = iy + inm1;
	
	return (an*iy);
}

double rangls(){
	static int iflag = 0;
	static double gauss2 = 0;
	double x1, x2, xsq=0, aux;
	if (iflag == 0){
		do{
			//Pair of uniform random numbers in [-1,1]x[-1,1]
			x1 = 2*ranuls() - 1;
			x2 = 2*ranuls() - 1;
			//if not in the unit circle, try again
			xsq = x1*x1 + x2*x2;
		} while (xsq >= 1 || xsq <= 0);

		aux = sqrt(-2.*log(xsq) / xsq);
		gauss2 = x2*aux;
		iflag = 0;
		return (x1*aux);
	}
	else{
		iflag = 0;
		return gauss2;
	}
}



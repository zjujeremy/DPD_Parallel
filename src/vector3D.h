#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

template <class numtype>
class vector3D{
	public:
		vector3D(numtype _x = 0, numtype _y = 0, numtype _z = 0) : x(_x), y(_y), z(_z){ }
		vector3D(const vector3D<numtype> & _v){
			x = _v.x;
			y = _v.y;
			z = _v.z;
		}

		inline vector3D<numtype> operator + (const vector3D<numtype> & _v){
			return vector3D<numtype>(x + _v.x, y + _v.y, z + _v.z);
		}

		inline vector3D<numtype> operator += (const vector3D<numtype> & _v){
			return vector3D<numtype>(x + _v.x, y + _v.y, z + _v.z);
		}

		inline vector3D<numtype> operator - (const vector3D<numtype> & _v){
			return vector3D<numtype>(x - _v.x, y - _v.y, z - _v.z);
		}

		inline vector3D<numtype> operator -= (const vector3D<numtype> & _v){
			return vector3D<numtype>(x - _v.x, y - _v.y, z - _v.z);
		}

		inline vector3D<numtype> operator * (const double & _d){
			return vector3D<numtype>(x*_d, y*_d, z*_d);
		}

		inline vector3D<numtype> operator * (const int & _d){
			return vector3D<numtype>(x*_d, y*_d, z*_d);
		}

		inline numtype operator * (const vector3D<numtype> & _v){
			return (x*(_v.x) + y*(_v.y) + z*(_v.z));
		}

		inline vector3D<numtype> multi(const vector3D<numtype> & _v){
			return vector3D<numtype>(x*(_v.x), y*(_v.y), z*(_v.z));
		}

		inline vector3D<numtype> mixmulti(const vector3D<numtype> & _v){
			return vector3D<numtype>(x*(_v.y), x*(_v.z), y*(_v.z));
		}

		inline vector3D<numtype> operator *= (const double & _d){
			return vector3D<numtype>(x*_d, y*_d, z*_d);
		}

		inline vector3D<numtype> operator / (const double & _d){
			if (_d == 0) throw _d;
			return vector3D<numtype>(x / _d, y / _d, z / _d);
		}

		inline vector3D<numtype> operator /= (const double & _d){
			if (_d == 0) throw _d;
			return vector3D<numtype>(x / _d, y / _d, z / _d);
		}

		inline vector3D<numtype> operator /= (const int & _d){
			if (_d == 0) throw _d;
			return vector3D<numtype>(x / _d, y / _d, z / _d);
		}

		inline numtype& operator [] (const int &_i){
			if (_i == 0)
				return x;
			else{
				if (_i == 1)
					return y;
				if (_i == 2)
					return z;
				else
					throw _i;
			}	
		}

		inline numtype norm2(){
			return x*x + y*y + z*z;
		}

		inline numtype norm(){
			return sqrt(norm2());
		}

		template <class numtype>
		inline friend vector3D<numtype> operator * (const double &, const vector3D<numtype> &);
		template <class numtype>
		inline friend ostream& operator << (ostream&, const vector3D<numtype> &);
		template <class numtype>
		inline friend ifstream& operator >> (ifstream&, vector3D<numtype> &);
	private:
		numtype x, y, z;
};

template <class numtype>
ostream& operator << (ostream & oupt, const vector3D<numtype> & _v){
	oupt << "x = " << _v.x << ", " << "y = " << _v.y << ", " << "z = " << _v.z << endl;
	return oupt;
}

template <class numtype>
ifstream& operator >> (ifstream & input, vector3D<numtype> & _v){
	input >> _v.x >> _v.y >> _v.z;
	return input;
}

template <class numtype>
vector3D<numtype> operator * (const double & _d, const vector3D<numtype> & _v){
	return vector3D<numtype>(_d*_v.x, _d*_v.y, _d*_v.z);
}

#endif


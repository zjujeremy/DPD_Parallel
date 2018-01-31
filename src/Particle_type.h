#ifndef PARTICLE_TYPE_H
#define PARTICLE_TYPE_H

#include <iostream>
#include "DPD_particle.h"

using namespace std;
class Wall_particle :public DPD_particle{
public:
	Wall_particle(){
		for (int i = 0; i < 3; i++){
			position[i] = 0;
			velocity[i] = 0;
			acceleration[i] = 0;
			ncc[i] = 0;
			wn[i] = 0;
		}
	}
	// copy constructor
	Wall_particle(const Wall_particle & _wp){
		position = _wp.position;
		wn = _wp.wn;
	}
	inline void setwn(const vector3D<double> & _wn);
	inline vector3D<double> getwn() const;
	virtual inline double getsqrDist() const;
	virtual inline void setsqrDist(double _d);
	virtual inline int getatomID();
	virtual inline void setatomID(const int _i);
private:
	vector3D<double> wn;
};

inline void Wall_particle::setwn(const vector3D<double> & _wn){
	wn = _wn;
}

inline vector3D<double> Wall_particle::getwn() const{
	return wn;
}

inline double Wall_particle::getsqrDist() const{
	return 0; // wall particle not have sqrDist
}
inline void Wall_particle::setsqrDist(double _d){}

inline int Wall_particle::getatomID(){
	return -1; // wall particle not have atomID
}

inline void Wall_particle::setatomID(const int _i){}



class Drop_particle :public DPD_particle{
public:
	Drop_particle(){
		for (int i = 0; i < 3; i++){
			position[i] = 0;
			velocity[i] = 0;
			acceleration[i] = 0;
			ncc[i] = 0;
			centerDrop[i] = 0;
		}
		atomID = 0;
	}

	Drop_particle(const Drop_particle & _drop){
		position = _drop.position;
	}
	inline void setCenterDrop(const vector3D<double> &);
	inline vector3D<double> & getCenterDrop();
	virtual inline int getatomID();
	virtual inline void setatomID(const int _i);
	virtual inline double getsqrDist() const;
	virtual inline void setsqrDist(double _d);
private:
	vector3D<double> centerDrop;
	int atomID;
	double sqrDist;
};

inline void Drop_particle::setCenterDrop(const vector3D<double> &_c){
	centerDrop = _c;
}

inline vector3D<double> & Drop_particle::getCenterDrop(){
	return centerDrop;
}

inline int Drop_particle::getatomID(){
	return atomID;
}

inline void Drop_particle::setatomID(const int _i){
	atomID = _i;
}

inline double Drop_particle::getsqrDist() const{
	return sqrDist;
}
inline void Drop_particle::setsqrDist(double _d){
	sqrDist = _d;
}


class Chain_particle :public DPD_particle{
public:
	Chain_particle(){
		for (int i = 0; i < 3; i++){
			position[i] = 0;
			velocity[i] = 0;
			acceleration[i] = 0;
			ncc[i] = 0;
		}
		atomID = 0;
	}

	Chain_particle(const Chain_particle & _Chain){
		position = _Chain.position;
	}
	virtual inline int getatomID(){
		return atomID;
	}
	virtual inline void setatomID(const int _i){
		atomID = _i;
	}
	virtual inline double getsqrDist() const{
		return sqrDist;
	}
	virtual inline void setsqrDist(double _d){
		sqrDist = _d;
	}
private:
	int atomID;
	double sqrDist;
};

class Fluid_particle :public DPD_particle{
public:
	Fluid_particle() {
		for (int i = 0; i < 3; i++){
			position[i] = 0;
			velocity[i] = 0;
			acceleration[i] = 0;
			ncc[i] = 0;
		}
		atomID = 0;
	}

	Fluid_particle(const Fluid_particle & _fp){        // copy constructor
		position = _fp.position;
	}
	virtual inline int getatomID(){
		return atomID;
	}
	virtual inline void setatomID(const int _i){
		atomID = _i;
	}
	virtual inline double getsqrDist() const{
		return sqrDist;
	}
	virtual inline void setsqrDist(double _d){
		sqrDist = _d;
	}
private:
	int atomID;
	double sqrDist;
};

#endif
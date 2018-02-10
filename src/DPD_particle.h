#ifndef DPD_PARTICLE_H
#define DPD_PARTICLE_H

#include <iostream>
#include "vector3D.h"

using namespace std;

class DPD_particle{
	public:
		DPD_particle(){}
		DPD_particle(vector3D<double> _p, vector3D<double> _v, vector3D<double> _a) :position(_p), velocity(_v), acceleration(_a){ }
		//DPD_particle(const DPD_particle & _particle);
		/*inline void outSituation();*/

		inline void setPosition(const vector3D<double> &_r);
		inline void setVelocity(const vector3D<double> & _v);
		inline void setAcceleration(const vector3D<double> & _a);

		inline vector3D<double> & getPosition();
		inline vector3D<double> & getVelocity();
		inline vector3D<double> & getAcceleration();
		inline vector3D<int> & getncc();

		virtual double getsqrDist() const = 0;
		virtual void setsqrDist(double _d) = 0;
		virtual int getatomID() = 0;
		virtual void setatomID(const int _i) = 0;
	public:
		vector3D<double> accelCV, accelCR, accelDP, accelRD, accelSP, accelG, rvm;
		vector3D<double> rforce, rforce2;
		double mrho;
		//paramenters for repulsion potential in MDPD
		list<DPD_particle*> tag;
		list<double> weightcd, drij;
		list<vector3D<double>> ddr;
	protected:
		vector3D<double> position, velocity, acceleration;
		vector3D<int> ncc;
};

//inline void DPD_particle::outSituation(){
//
//}

inline void DPD_particle::setPosition(const vector3D<double> & _r){
	position = _r;
}

inline void DPD_particle::setVelocity(const vector3D<double> & _v){
	velocity = _v;
}

inline void DPD_particle::setAcceleration(const vector3D<double> & _a){
	acceleration = _a;
}

inline vector3D<double> & DPD_particle::getPosition() {
	return position;
}

inline vector3D<double> & DPD_particle::getVelocity(){
	return velocity;
}

inline vector3D<double> & DPD_particle::getAcceleration(){
	return acceleration;
}

inline vector3D<int> & DPD_particle::getncc(){
	return ncc;
}

#endif
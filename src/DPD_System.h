#ifndef DPD_SYSTEM_H
#define DPD_SYSTEM_H

#include <iostream>
#include <cmath>
#include <list>
#include <vector>
#include "Particle_type.h"
#include "randomFun.h"
using namespace std;

//runID:
//0 --no flow, periodic BC in 3 directions
//1 --simple shear flow or Poiseuille Flow, wall BC in Z - direction
//2 --periodic Poiseuille flow for viscosity test, periodic BC in 3 directions
//3 --elongation flow, periodic BC in 3 directions
//4 --pipe flow
//5 --simple shear flow with Lees - Edwards BC

class DPD_System{
	public:
		DPD_System(){
			moreCycle = 1;
			stepCount = 0;
			NDIM = 3;
			timeNow = 0;
		}
		//function for add  four (or less) type of the particle from the preprocess 
		template<typename T>
		inline void readEachParticle(ifstream & infile, unsigned int start, unsigned int end, vector<T> & _v);

		void initVels();    // initialize the velocity of particles
		void initAccels();  // initialize the acceleration of particles
		template<typename T>
		void setRandomVelocity(vector3D<double> & _e, vector3D<double> & _vTemp, vector<T> & _v);  // set random velocity for each type of particles 
		template<typename T>
		void normRandomVelocity(vector3D<double> & _vTemp, vector<T> & _v);  // eliminate the bias of random velocity to origin.
		template<typename T>
		void initEachTypeParticleAccels(T &); // using function template to initialize the each type of particles' acceleration

		void velocity_verlet(); // using velocity-verlet algorithm 
		template<typename T>
		void setNextTimeData(vector<T> & _v);
		void ApplyBoundaryCond();
		template<typename T>
		void applyBCtoEachParticle(vector<T> & _v);
		template<typename T>
		void setNearWallParicleVels(vector<T> & _v);

		void computeForces();
		void computeExternalForce();
		void updataAccelAndVelos();
		template<typename T>
		void updataAccelAndVelosub(vector<T> & _v);
		template<typename T>
		void comExForceforEachPart(vector<T> & _v);
		void computeParticleInEachCell();
		template<typename T>
		void handleEachTypeParticle(vector<T> & _v, const int & _ibase);
		template<typename T>
		void initForcePara(vector<T> & _v);
		void computeForceSub(DPD_particle* _p1, DPD_particle* _p2, const vector3D<double> & _shift);
		template<typename T>
		void computeRepForce(vector<T> & _v);
		

		void readIntiData();  // read the DPD particles' initial data
		void readConf();      // read the configuration in DPD system 
		void setParams();
		void setupJob();  

		void SingleStep(); // single step in DPD system
		void outputParticleSituation(); 

		inline int checkCycle();   // check the cycle 
		inline unsigned int getstepCount();
	private:
		//compute forces paraments
		int maxLists;
		list<int> *cellsList;
		vector3D<double> invWid;
		double alpha, cigama, gamma;
		double uSum, virSum;

		// DPD field parameters
		vector<Wall_particle> vectorWallPtc;
		vector<Drop_particle> vectorDropPtc;
		vector<Chain_particle> vectorChainPtc;
		vector<Fluid_particle> vectorFluidPtc;
		vector3D<int> initUcell, cells;
		vector3D<double> region, regionH, gap;
		unsigned int moreCycle, stepCount, stepLimit, NDIM, PNDIM, hSize;
		int nWallAtom, nDpEnd, nChainEnd, nAtom, nDp, DpFzn[10], nChain, ChainLen, nPDrop, nFreeAtom, nStartAtom;
		double RdsDp, wMingap, wLayer;
		double vMag, binVolm, sdtinv, rrfene;
		double timeNow;

		// DPD method paramenters 
		int runID;
		double shearRate, alphaf, alphafp, alphapp, alphaFD, alphaDD, alphaB;
		double alphaw, alphaw_top, alphawB;
		double rCut, rCut2, gammaF, gammaD, gammaFD, gammaW;
		double cigamaF, cigamaD, cigamaFD, cigamaW;
		double lambda, density, temperature, mass, WLadjust, Hfene, rmaxfene, reqfene;
		double deltaT, timeSteady;
		unsigned int stepAvg, stepEquil, startSample, stepSample, stepGrid, limitGrid;
		unsigned int stepChainProps, limitChainProps, nChainConf;
		vector3D<double> gravField;
		vector3D<int> sizeHistGrid;
};

inline int DPD_System::checkCycle(){
	if (stepCount == stepLimit)
		moreCycle = 0;
	return moreCycle;
}

inline unsigned int DPD_System::getstepCount(){
	return stepCount;
}

template<typename T>
inline void DPD_System::readEachParticle(ifstream & infile, unsigned int start, unsigned int end, vector<T> & _v){
	for (unsigned int n = start; n < end; n++) {
		T _p;
		vector3D<double> _pos;
		infile >> _pos;
		_p.setPosition(_pos);
		_v.push_back(_p);
	}
}

#endif
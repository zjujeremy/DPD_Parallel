#include <iostream>
#include <fstream>
#include <strstream>
#include <string>
#include <cmath>
#include <list>
#include "randomFun.h"
#include "DPD_System.h"
using namespace std;

void DPD_System::readIntiData(){
	ifstream infile("./data/initdpd.cnf", ios::in);
	if (!infile){
		cout << "open file('initdpd.cnf') error" << endl;
		throw 1;
	}
	infile >> nWallAtom >> nDpEnd >> nChainEnd >> nAtom >> nDp >> RdsDp;
	if (nDp > 0)
		for (int i = 0; i < nDp; i++)
			infile >> DpFzn[i];
	infile >> nChain >> ChainLen;
	infile >> initUcell >> region >> regionH >> gap;
	infile >> wMingap >> wLayer;
	
	readEachParticle(infile,         0, nWallAtom, vectorWallPtc);
	readEachParticle(infile, nWallAtom,    nDpEnd, vectorDropPtc);
	readEachParticle(infile,    nDpEnd, nChainEnd, vectorChainPtc);
	readEachParticle(infile, nChainEnd,     nAtom, vectorFluidPtc);

	for (auto & wp : vectorWallPtc){
		vector3D<double> _wn;
		infile >> _wn;
		wp.setwn(_wn);
	}

	infile.close();

}

void DPD_System::readConf(){
	ifstream infile("./config/config.dat", ios::in);
	string inform;
	if (!infile){
		cout << "open file('config.dat') error" << endl;
		exit(1);
	}
	int head = 1;
	while (getline(infile, inform)){
		switch (head){
			case 1: 
				infile >> runID;
				break;
			case 3: 
				infile >> shearRate;
				break;
			case 5: 
				infile >> alphaf >> alphafp >> alphapp >> alphaFD >> alphaDD >> alphaB;
				break;
			case 7: 
				infile >> alphaw >> alphaw_top >> alphawB;
				break;
			case 9:
				infile >> rCut >> rCut2;
				break;
			case 11:
				infile >> gammaF >> gammaD >> gammaFD;
				break;
			case 13:
				infile >> gammaW;
				break;
			case 15:
				infile >> lambda;
				break;
			case 17:
				infile >> density;
				break;
			case 19:
				infile >> temperature >> mass;
				break;
			case 21:
				infile >> gravField;
				break;
			case 23:
				infile >> WLadjust;
				break;
			case 25:
				infile >> Hfene >> rmaxfene >> reqfene;
				break;
			case 27:
				infile >> deltaT;
				break;
			case 29:
				infile >> stepAvg;
				break;
			case 31:
				infile >> stepEquil;
				break;
			case 33:
				infile >> startSample;
				break;
			case 35:
				infile >> stepSample;
				break;
			case 37:
				infile >> stepLimit;
				break;
			case 39:
				infile >> sizeHistGrid;
				break;
			case 41:
				infile >> stepGrid;
				break;
			case 43:
				infile >> limitGrid;
				break;
			case 45:
				infile >> stepChainProps;
				break;
			case 47:
				infile >> limitChainProps;
				break;
			case 49:
				infile >> nChainConf;
				break;
			case 51:
				infile >> timeSteady;
				break;
		}
		getline(infile, inform);
		head += 2;
	}

	infile.close();
}

void DPD_System::setParams(){
	if (runID == 1){
		nStartAtom = 1;
		PNDIM = NDIM - 1;
	}	
	else{
		nStartAtom = nWallAtom + 1;
		PNDIM = NDIM;
	}

	if (runID != 1){
		region[2] = initUcell[2] * gap[2];
		regionH[2] = region[2] / 2.0;
	}

	for (unsigned int i = 0; i < NDIM; i++){
		cells[i] = int(region[i] / rCut);
		if (cells[i] < 1)
			cells[i] = 1;
	}
	
	if (runID == 5)    // using simple shear flow with Lees-Edwards BC  
		cells[2] = cells[2] + 1;  

	if (nDp > 0)
		nPDrop = int((nDpEnd - nWallAtom) / nDp);
	else
		nPDrop = nDpEnd - nWallAtom;

	maxLists = cells[0] * cells[1] * cells[2];
	cellsList = new list<int>[maxLists];
	if (cellsList == NULL){
		cout << "Application heap memory failure for cellsList pointer" << endl;
		throw 1;
	}

	nFreeAtom = nAtom - nWallAtom;
	vMag = sqrt(NDIM*(1.0 - 1.0 / nFreeAtom)*temperature / mass);
	hSize = sizeHistGrid[0] * sizeHistGrid[1] * sizeHistGrid[2];    //divide the whole domain into many bins at three directions 
	                                                                //for compute density, average velocity, pressure and so on.
	binVolm = region[0] * region[1] * region[2] / hSize;
	cigamaF = sqrt(gammaF * 2 * temperature);
	cigamaD = sqrt(gammaD * 2 * temperature);
	cigamaFD = sqrt(gammaFD * 2 * temperature);
	cigamaW = sqrt(gammaW * 2 * temperature);

	sdtinv = 1.0 / sqrt(deltaT);
	wLayer = WLadjust*wLayer;

	rrfene = pow((rmaxfene - reqfene), 2);  // for Chain
}

void DPD_System::setupJob(){
	// generate a random number
	::ranils(290092405);  // using global scope function
	initVels();
	initAccels();
}

void DPD_System::SingleStep(){

	stepCount++;
	if (stepCount % 100 == 0)
		cout << "stepCount = " << stepCount << endl;
	if (stepCount > stepEquil)
		timeNow += deltaT;

	velocity_verlet(); 

}

void DPD_System::velocity_verlet(){
	setNextTimeData(vectorDropPtc);
	setNextTimeData(vectorChainPtc);
	setNextTimeData(vectorFluidPtc);

	ApplyBoundaryCond();

	setNearWallParicleVels(vectorDropPtc);
	setNearWallParicleVels(vectorChainPtc);
	setNearWallParicleVels(vectorFluidPtc);

	computeForces();
	computeExternalForce(); 
	
	updataAccelAndVelos();//compute sum of the accelration and update velocity for each particle
}


void DPD_System::initVels(){
	vector3D<double> e, vTemp;
	setRandomVelocity(e, vTemp, vectorDropPtc);
	setRandomVelocity(e, vTemp, vectorChainPtc);
	setRandomVelocity(e, vTemp, vectorFluidPtc);

	vTemp /= nFreeAtom;
	normRandomVelocity(vTemp, vectorDropPtc);
	normRandomVelocity(vTemp, vectorChainPtc);
	normRandomVelocity(vTemp, vectorFluidPtc);
}

void DPD_System::initAccels(){
	initEachTypeParticleAccels(vectorWallPtc);
	initEachTypeParticleAccels(vectorDropPtc);
	initEachTypeParticleAccels(vectorChainPtc);
	initEachTypeParticleAccels(vectorFluidPtc);
}

template<typename T>
void DPD_System::initEachTypeParticleAccels(T & _vp){
	for (auto & vp : _vp){
		for (unsigned int i = 0; i < NDIM; i++){
			vp.setAcceleration(vector3D<double>(0, 0, 0));  
			vp.accelCV[i] = 0;
			vp.accelCR[i] = 0;
			vp.accelDP[i] = 0;
			vp.accelRD[i] = 0;
			vp.accelSP[i] = 0;
			vp.accelG[i] = 0;
			vp.rforce[i] = 0;
			vp.rforce2[i] = 0;
		}
	}
}

template<typename T>
void DPD_System::setRandomVelocity(vector3D<double> & _e, vector3D<double> & _vTemp, vector<T> & _v){
	for (auto & _p : _v){
		randVec3(_e);
		_p.setVelocity(vMag*_e);
		_vTemp += _p.getVelocity();
	}
}

template<typename T>
void DPD_System::normRandomVelocity(vector3D<double> & _vTemp, vector<T> & _v){
	for (auto & _p : _v)
		_p.setVelocity(_p.getVelocity() - _vTemp);
}

template<typename T>
void DPD_System::setNextTimeData(vector<T> & _v){
	vector3D<double> temp;
	for (auto & _vp : _v){
		temp = deltaT*_vp.getAcceleration();
		_vp.rvm = _vp.getVelocity();
		_vp.setPosition(_vp.getPosition() + deltaT*(_vp.getVelocity() + 0.5*temp));
		_vp.setVelocity(_vp.getVelocity() + lambda*temp);
		_vp.rvm += 0.5*temp;
	}
}

void DPD_System::ApplyBoundaryCond(){
	applyBCtoEachParticle(vectorDropPtc);
	applyBCtoEachParticle(vectorChainPtc);
	applyBCtoEachParticle(vectorFluidPtc);
}

template<typename T>
void DPD_System::applyBCtoEachParticle(vector<T> & _v){
	for (auto & _vp : _v){
		for (unsigned int i = 0; i < PNDIM; i++){
			if (_vp.getPosition()[i] > region[i] || _vp.getPosition()[i] < -region[i]){
				cout << "molcule moves too far, reduce time step" << endl;
				throw 1;
			}

			if (_vp.getPosition()[i] >= regionH[i]){
				_vp.getPosition()[i] -= region[i];
				_vp.getncc()[i] += 1;
			}
			
			if (_vp.getPosition()[i] < -regionH[i]){
				_vp.getPosition()[i] += region[i];
				_vp.getncc()[i] -= 1;
			}
		}
	}
}

template<typename T>
void DPD_System::setNearWallParicleVels(vector<T> & _v){
	vector3D<double> e;
	for (auto & _vp : _v){
		if (_vp.getatomID() > 0){
			randVec3(e);
			if (stepCount <= stepEquil){
				_vp.setVelocity(0.5*vMag*e);
				_vp.rvm = _vp.getVelocity();
			}
			else{
				_vp.getVelocity()[0] = 1.2*vMag*e[0] + _vp.getPosition()[2] * shearRate;
				_vp.rvm[0] = _vp.getVelocity()[0];
				_vp.getVelocity()[1] = 1.2*vMag*e[1];
				_vp.rvm[1] = _vp.getVelocity()[1];
				_vp.getVelocity()[2] = 1.2*vMag*e[2];
				_vp.rvm[2] = _vp.getVelocity()[2];
			}
		}
	}
}

void DPD_System::computeForces(){
	
	static const int cellOffset[14][3] = {
		{ 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { -1, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, 
		{ 1, 1, 1 }, { 0, 1, 1 }, { -1, 1, 1 }, { -1, 0, 1 }, { -1, -1, 1 }, { 0, -1, 1 }, { 1, -1, 1 }
	};

	for (unsigned int k = 0; k < NDIM; k++)
		invWid[k] = double(cells[k] / region[k]);

	if (runID == 5)
		invWid[NDIM - 1] = double((cells[NDIM - 1] - 1) / region[NDIM - 1]);

	for (int i = 0; i < maxLists; i++)
		cellsList[i].clear();

	initAccels(); 

	computeParticleInEachCell();

	initForcePara(vectorWallPtc);
	initForcePara(vectorDropPtc);
	initForcePara(vectorChainPtc);
	initForcePara(vectorFluidPtc);

	int cellZ;
	if (runID == 5)
		cellZ = cells[2] - 1;
	else
		cellZ = cells[2];

	uSum = 0;
	virSum = 0;

	int m1, m2;
	vector3D<double> shift;
	list<int>::iterator iter1, iter2;
	for (int m1Z = 0; m1Z < cellZ; m1Z++){
		for (int m1Y = 0; m1Y < cells[1]; m1Y++){
			for (int m1X = 0; m1X < cells[0]; m1X++){

				m1 = (m1Z*cells[1] + m1Y)*cells[0] + m1X;

				for (int offset = 0; offset < 14; offset++){

					int m2X = m1X + cellOffset[offset][0];
					shift[0] = 0;
					if (m2X >= cells[0]){
						m2X = 0;
						shift[0] = region[0];
					}
					else{
						if (m2X < 0){
							m2X = cells[0] - 1;
							shift[0] = -region[0];
						}
					}
						
					int m2Y = m1Y + cellOffset[offset][1];
					shift[1] = 0;
					if (m2Y >= cells[1]){
						m2Y = 0;
						shift[1] = region[1];
					}
					else{
						if (m2Y < 0){
							m2Y = cells[1] - 1;
							shift[1] = -region[1];
						}
					}

					int m2Z = m1Z + cellOffset[offset][2];
					shift[2] = 0;
					if (runID == 1 || runID == 5){
						if (m2Z < 0 || m2Z >= cells[2])
							continue;
					}
					else{
						if (m2Z >= cells[2]){
							m2Z = 0;
							shift[2] = region[2];
						}
						else{
							if (m2Z < 0){
								m2Z = cells[2] - 1;
								shift[2] = -region[2];
							}
						}
					}

					m2 = ((m2Z*cells[1] + m2Y)*cells[0]) + m2X;

					for (iter1 = cellsList[m1].begin(); iter1 != cellsList[m1].end(); iter1++){  
						for (iter2 = cellsList[m2].begin(); iter2 != cellsList[m2].end(); iter2++){
							if (m1 != m2 || *iter1 > *iter2){
								const int j1 = *iter1;
								const int j2 = *iter2;

								//both j1 and j2 are wall particles
								if (j1 <= nWallAtom&&j2 <= nWallAtom){
									Wall_particle *wp1 = &vectorWallPtc[j1 - 1];
									Wall_particle *wp2 = &vectorWallPtc[j2 - 1];
									alpha = 0;
									cigama = 0;
									gamma = 0;
									computeForceSub(wp1, wp2, shift);
								}

								//one of j1 or j2 is a wall particle and another is a chain(or drop) bead or solvent particle
								if ((j1 <= nWallAtom && j2 > nWallAtom) || (j1 > nWallAtom && j2 <= nWallAtom)){ 
									DPD_particle *dp = NULL;
									Wall_particle *wp = NULL;
									bool ord;
									vector3D<double> dr;
									if (j1 < j2){  //j1 is a wall particle
										ord = true;
										wp = &vectorWallPtc[j1 - 1];
										if (j2 <= nDpEnd){
											dp = &vectorDropPtc[j2 - nWallAtom - 1];
										}
										if (j2 > nDpEnd&&j2 <= nChainEnd){
											dp = &vectorDropPtc[j2 - nDpEnd - 1];
										}
										if (j2 > nChainEnd){
											dp = &vectorFluidPtc[j2 - nChainEnd - 1];
										}
										dr = wp->getPosition() - dp->getVelocity() - shift;
									}
									else{ //j2 is a wall particle
										ord = false;
										wp = &vectorWallPtc[j2 - 1];
										if (j1 <= nDpEnd){
											dp = &vectorDropPtc[j1 - nWallAtom - 1];
										}
										if (j1 > nDpEnd&&j1 <= nChainEnd){
											dp = &vectorDropPtc[j1 - nDpEnd - 1];
										}
										if (j1 > nChainEnd){
											dp = &vectorFluidPtc[j1 - nChainEnd - 1];
										}	
										dr = dp->getPosition() - wp->getVelocity() - shift;
									}

									double rr = dr.norm2();
									double rij = dr.norm();
									if (rij < rCut){
										int tag1 = j1 > j2 ? j1 : j2; 
										int tag2 = j1 < j2 ? j1 : j2;  // wall particle
										double dn = dr*(wp->getwn());
										if (abs(dn) <= wLayer && rr < dp->getsqrDist()){
											dp->setsqrDist(rr);
											dp->setatomID(tag2);
										}
									}

									if (ord)
										computeForceSub(wp, dp, shift);
									else
										computeForceSub(dp, wp, shift);
								}


								//both j1 and j2 are droplet particle
								if (j1 > nWallAtom && j1 <= nDpEnd && j2 > nWallAtom && j2 <= nDpEnd){
									DPD_particle* dp1 = &vectorDropPtc[j1 - nWallAtom - 1];
									DPD_particle* dp2 = &vectorDropPtc[j2 - nWallAtom - 1];
									alpha = alphaDD;
									cigama = cigamaD;
									gamma = gammaD;
									computeForceSub(dp1, dp2, shift);
								}

								//one of j1 and j2 is a drop particle and another is a solvent particle or chain particle
								if ((j1 > nWallAtom && j1 <= nDpEnd && j2 > nDpEnd) || (j2 > nWallAtom && j2 <= nDpEnd && j1 > nDpEnd)){
									DPD_particle *dp = NULL;
									DPD_particle *fp = NULL;
									bool ord;
									if (j1 < j2){
										ord = true;
										dp = &vectorDropPtc[j1 - nWallAtom - 1];
										if (j2 <= nChainEnd)
											fp = &vectorChainPtc[j2 - nDpEnd - 1];
										else
											fp = &vectorFluidPtc[j2 - nChainEnd - 1];
									}
									else{
										ord = false;
										dp = &vectorDropPtc[j2 - nWallAtom - 1];
										if (j1 <= nChainEnd)
											fp = &vectorChainPtc[j1 - nDpEnd - 1];
										else
											fp = &vectorFluidPtc[j1 - nChainEnd - 1];
									}
									alpha = alphaFD;
									cigama = cigamaFD;
									gamma = gammaFD;
									if (ord)
										computeForceSub(dp, fp, shift);
									else
										computeForceSub(fp, dp, shift);	
								}

								//both j1 and j2 are the beads of chains
								if (j1 > nDpEnd && j1 <= nChainEnd && j2 > nDpEnd && j2 <= nChainEnd){
									DPD_particle* cp1 = &vectorChainPtc[j1 - nDpEnd - 1];
									DPD_particle* cp2 = &vectorChainPtc[j2 - nDpEnd - 1];
									alpha = alphapp;
									cigama = cigamaF;
									gamma = gammaF;
									computeForceSub(cp1, cp2, shift);
								}

								//one of j1 and j2 is a bead of chains and another is a solvent particle
								if ((j1 > nDpEnd && j1 <= nChainEnd && j2 > nChainEnd) || (j2 > nDpEnd && j2 <= nChainEnd && j1 > nChainEnd)){
									DPD_particle *cp = NULL;
									DPD_particle *fp = NULL;
									bool ord;
									if (j1 < j2){
										ord = true;
										cp = &vectorChainPtc[j1 - nDpEnd - 1];
										fp = &vectorFluidPtc[j2 - nChainEnd - 1];
									}
									else{
										ord = false;
										cp = &vectorChainPtc[j2 - nDpEnd - 1];
										fp = &vectorFluidPtc[j1 - nChainEnd - 1];
									}
									alpha = alphafp;
									cigama = cigamaF;
									gamma = gammaF;
									if (ord)
										computeForceSub(cp, fp, shift);
									else
										computeForceSub(fp, cp, shift);
									
								}

								//both j1 and j2 are solvent particles
								if (j1 > nChainEnd && j2 > nChainEnd){
									DPD_particle* fp1 = &vectorFluidPtc[j1 - nChainEnd - 1];
									DPD_particle* fp2 = &vectorFluidPtc[j2 - nChainEnd - 1];
									alpha = alphaf;
									cigama = cigamaF;
									gamma = gammaF;
									computeForceSub(fp1, fp2, shift);
								}
							}
						}
					}

				}
			}
		}
	}
	// compute repulsive potential in MDPD
	computeRepForce(vectorDropPtc);
	computeRepForce(vectorChainPtc);
	computeRepForce(vectorFluidPtc);
}

void DPD_System::computeParticleInEachCell(){
	handleEachTypeParticle(vectorWallPtc, 0);
	handleEachTypeParticle(vectorDropPtc, nWallAtom);
	handleEachTypeParticle(vectorChainPtc, nDpEnd);
	handleEachTypeParticle(vectorFluidPtc, nChainEnd);
}

template<typename T>
void DPD_System::handleEachTypeParticle(vector<T> & _v, const int & _ibase){
	int c, n = _ibase;
	for (auto & _vp : _v){
		n += 1;
		c = (int((_vp.getPosition()[2] + regionH[2])*invWid[2])*cells[1] + \
			int((_vp.getPosition()[1] + regionH[1])*invWid[1]))*cells[0] + \
			int((_vp.getPosition()[0] + regionH[0])*invWid[0]);

		if (c < 0){
			cout << "in cf, stepCount = " << stepCount << endl;
			throw 1;
		}

		if (c >= maxLists){
			cout << " stepCount = "<< stepCount << endl;
			cout << " c = " << c << endl;
			throw 1;
		}
		cellsList[c].push_front(n);
	}
}

template<typename T>
void DPD_System::initForcePara(vector<T> & _v){
	for (auto & _vp : _v){
		_vp.setsqrDist(100);
		_vp.setatomID(0);
		_vp.tag.clear();
		_vp.weightcd.clear();
		_vp.drij.clear();
		_vp.ddr.clear();
		_vp.mrho = 0;
	}
}

void DPD_System::computeForceSub(DPD_particle* _p1, DPD_particle* _p2, const vector3D<double> & _shift){
	const double pi = 3.14159265358979;
	vector3D<double> dr = _p1->getPosition() - _p2->getVelocity() - _shift;
	vector3D<double> dv = _p1->getVelocity() - _p2->getVelocity();

	double rr = dr.norm2();
	double rij = dr.norm();
	dr = dr / rij;
	double rdvij = dr*dv;
	if (rij<rCut){
		double weight = 1 - rij / rCut;
		double weightc = weight;
		double weightd = pow(weight, 2);
		double weightr = weight;

		double fcVal = alpha*weightc;
		double uVal = -alpha*(rij - 0.5*rr);

		double fdVal = -gamma*weightd*rdvij;
		double frVal = cigama*weightr*rangls()*sdtinv;
		double ftot = fcVal + fdVal + frVal;
		double fstl = ftot*rij;

		vector3D<double> fCV = dr*fcVal;
		vector3D<double> fDP = dr*fdVal;
		vector3D<double> fRD = dr*frVal;

		_p1->accelCV += fCV / mass;
		_p2->accelCV -= fCV / mass;

		_p1->accelDP += fDP / mass;
		_p2->accelDP -= fDP / mass;

		_p1->accelRD += fRD / mass;
		_p2->accelRD -= fRD / mass;

		vector3D<double> fr = (dr.multi(dr))*fstl;
		_p1->rforce += fr;
		_p2->rforce += fr;

		fr = (dr.mixmulti(dr))*fstl;
		_p1->rforce2 += fr;
		_p2->rforce2 += fr;

		uSum += uVal;
		virSum += fcVal*rij;

		// set paramenters for repulsion potential in MDPD
		if (rij < rCut2){
			double weightp = 15 / (2 * pi*pow(rCut2, 3))*pow((1 - rij / rCut2), 2);
			_p1->mrho += weightp;
			_p2->mrho += weightp;
			(_p1->tag).push_back(_p2);
			(_p2->tag).push_back(_p1);
			(_p1->weightcd).push_back(1 - rij / rCut2);
			(_p2->weightcd).push_back(1 - rij / rCut2);
			(_p1->drij).push_back(rij);
			(_p2->drij).push_back(rij);
			(_p1->ddr).push_back(dr);
			(_p2->ddr).push_back(-1 * dr);
		}
	}
}

template<typename T>
void DPD_System::computeRepForce(vector<T> & _v){
	for (auto & _vp : _v){
		list<DPD_particle*>::iterator iter_tag;
		list<vector3D<double>>::iterator iter_ddr;
		list<double>::iterator iter_weightcd, iter_drij;
		iter_weightcd = _vp.weightcd.begin();
		iter_drij = _vp.drij.begin();
		iter_tag = _vp.tag.begin();
		iter_ddr = _vp.ddr.begin();

		for (; iter_tag != _vp.tag.end(); iter_tag++, iter_weightcd++, iter_drij++, iter_ddr++){
			double fcVal = (*iter_weightcd)*alphaB*(_vp.mrho + (*iter_tag)->mrho);
			double fstl = fcVal*(*iter_drij);
			vector3D<double> fCR = (*iter_ddr)*fcVal;
			_vp.accelCR += fCR / mass;

			vector3D<double> fr = (iter_ddr->multi(*iter_ddr))*fstl;
			_vp.rforce += fr;
			fr = (iter_ddr->mixmulti(*iter_ddr))*fstl;
			_vp.rforce2 += fr;
		}
	}
}

void DPD_System::computeExternalForce(){
	comExForceforEachPart(vectorDropPtc);
	comExForceforEachPart(vectorChainPtc);
	comExForceforEachPart(vectorFluidPtc);
}

template<typename T>
void DPD_System::comExForceforEachPart(vector<T> & _v){
	if ((runID == 1 || runID == 4) && stepCount >= stepEquil&&gravField.norm2() != 0){
		for (auto & _vp : _v){
			_vp.accelG = gravField;
		}
	}

	if (runID == 2 && gravField.norm2() != 0 && stepCount >= stepEquil){
		for (auto & _vp : _v){
			_vp.accelG[0] = _vp.getPosition()[0] > 0 ? gravField[0] : -gravField[0];
		}
	}
}

void DPD_System::updataAccelAndVelos(){
	updataAccelAndVelosub(vectorDropPtc);
	updataAccelAndVelosub(vectorChainPtc);
	updataAccelAndVelosub(vectorFluidPtc);
}
template<typename T>
void DPD_System::updataAccelAndVelosub(vector<T> & _v){
	for (auto & _vp : _v){
		_vp.setAcceleration(_vp.accelCV + _vp.accelCR + _vp.accelDP + _vp.accelRD + _vp.accelSP + _vp.accelG);
		int n = _vp.getatomID();
		if (n <= 0){
			_vp.setVelocity(_vp.rvm + 0.5*deltaT*_vp.getAcceleration());
		}
		else{
			Wall_particle & wp = vectorWallPtc[n - 1];
			vector3D<double > e, vr;
			double vn = 0;
			randVec3(e);
			if (stepCount >= stepEquil)
				vr[0] = 1.2*vMag*e[0] + _vp.getPosition()[2] * shearRate;
			else
				vr[0] = 1.3*vMag*e[0];

			vn = vn + vr[0] * wp.getwn()[0];
			vr[1] = 1.2*vMag*e[1];
			vn = vn + vr[1] * wp.getwn()[1];
			vr[2] = 1.2*vMag*e[2];
			vn = vn + vr[2] * wp.getwn()[2];
			
			_vp.setVelocity(vr + (abs(vn) - vn)*wp.getwn());
		}
	}
}

void DPD_System::outputParticleSituation() {
	char filename[60];
	sprintf_s(filename, "./data/instantSituation/Situation_Step_%d.plt", stepCount);
	ofstream outFile(filename, ios::out);
	if (!outFile){
		cout << "open outfile error" << endl;
		throw 1;
	}
	outFile << "TITLE = \"All Particles Coordinates\"" << endl;
	outFile << "Variables= \"X\",\"Y\",\"Z\",\"Vx\",\"Vy\",\"Vz\",\"ax\",\"ay\",\"az\"" << endl;

	if (vectorWallPtc.size() > 0 && runID != 0 && runID != 6){
		outFile << "ZONE" << endl;
		outFile << "T=\"WallPariticles\"" << endl;
		outFile << "I = " << vectorWallPtc.size() << ", F = POINT" << endl;
		for (auto & _p : vectorWallPtc){
			outFile << _p.getPosition() << _p.getVelocity() << _p.getAcceleration() << endl;
		}
	}
	if (vectorDropPtc.size() > 0){
		outFile << "ZONE" << endl;
		outFile << "T=\"DropPariticles\"" << endl;
		outFile << "I = " << vectorDropPtc.size() << ", F = POINT" << endl;
		for (auto & _p : vectorDropPtc){
			outFile << _p.getPosition() << _p.getVelocity() << _p.getAcceleration() << endl;
		}
	}
	if (vectorChainPtc.size() > 0){
		outFile << "ZONE" << endl;
		outFile << "T=\"ChainPariticles\"" << endl;
		outFile << "I = " << vectorChainPtc.size() << ", F = POINT" << endl;
		for (auto & _p : vectorChainPtc){
			outFile << _p.getPosition() << _p.getVelocity() << _p.getAcceleration() << endl;
		}
	}
	if (vectorFluidPtc.size() > 0){
		outFile << "ZONE" << endl;
		outFile << "T=\"FluidPariticles\"" << endl;
		outFile << "I = " << vectorFluidPtc.size() << ", F = POINT" << endl;
		for (auto & _p : vectorFluidPtc){
			outFile << _p.getPosition() << _p.getVelocity() << _p.getAcceleration() << endl;
		}
	}
	outFile.close();
}
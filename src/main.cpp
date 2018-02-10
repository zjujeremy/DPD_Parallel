#include <iostream>
#include <omp.h>
#include <cstdlib>
#include "DPD_System.h"

using namespace std;

DPD_System dpdFlow;
int main(){
	try{
		dpdFlow.readIntiData();
		dpdFlow.readConf();
		dpdFlow.setParams();
		dpdFlow.setupJob();
		while (dpdFlow.checkCycle()){
			if (dpdFlow.getstepCount() % 500 == 0)
				dpdFlow.outputParticleSituation();
			dpdFlow.SingleStep();
		}
	}
	catch (...){
		cout << "error" << endl;
	}
	system("pause");
	return 0;
}


//#pragma omp parallel sections//定义以下的代码块用4个线程同时处理
//   {
//	#pragma omp section//第一块并行部分
//	{
//		int i = omp_get_thread_num();//获取每个线程的序号
//		printf_s("from thread %d\n", i);//结果打印四条序号不同的hello...
//	}
//	#pragma omp section//第二块并行部分
//	{
//		int i = omp_get_thread_num();//获取每个线程的序号
//		printf_s("from thread %d\n", i);//结果打印四条序号不同的hello...
//	}
//   }
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


//#pragma omp parallel sections//�������µĴ������4���߳�ͬʱ����
//   {
//	#pragma omp section//��һ�鲢�в���
//	{
//		int i = omp_get_thread_num();//��ȡÿ���̵߳����
//		printf_s("from thread %d\n", i);//�����ӡ������Ų�ͬ��hello...
//	}
//	#pragma omp section//�ڶ��鲢�в���
//	{
//		int i = omp_get_thread_num();//��ȡÿ���̵߳����
//		printf_s("from thread %d\n", i);//�����ӡ������Ų�ͬ��hello...
//	}
//   }
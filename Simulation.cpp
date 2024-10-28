// Service.cpp : �������̨Ӧ�ó������ڵ㡣
//
#include <stdio.h>
#include<iostream>
#include <iomanip> //��ȷС�����λ��

#include "Param.h"
#include "Simulation.h"
#include <stdlib.h>   //���������ͷ�ļ�
#include <map>
#include <vector>
#include <string>
#include <math.h>
using namespace std;
void Simulation::minMatrix()
{
	m_database.Suurballe_phyMatrix(1);
	m_database.linkStateInitial();

}

void Simulation::mapping(int ithTimeSFC,int SFC)
{
	int ithSFC = 0;
	sfcSet getASFC;
	vector<SFCMapPhy> aSFCMapSet;
	for (int i = 0; i < MaxNodeNum;i++)
	{
		m_database.phyVnfPair[i] = -1;
	}

	pq_getPhyNodeSer = pq_phyNodeSer;
	clock_t timeStart = clock();
	while(!pq_sfcSet.empty())
	{
		//clock_t timeStart = clock();//time �ĵ�λ�Ǻ��� //clock_t ���� ��Ļ���Ҫ����1000

		getASFC = pq_sfcSet.top();
		m_database.getaSFCMapOrder(getASFC); 

		map<int,vector<SFCMapPhy> >::iterator SFCMap = m_database.SFCMap.find(getASFC.SFC_id);
		aSFCMapSet = SFCMap->second;
		m_database.aSFCMapping(ithTimeSFC,getASFC,aSFCMapSet);
		pq_sfcSet.pop();
		if (pq_sfcSet.size() % 5 == 0 || pq_sfcSet.size() == 0)
		{
			blocka[figNodeNum - pq_sfcSet.size() / 5 - 1] += blockTime2;
			clock_t timeEnd = clock();
			double duration = double(timeEnd - timeStart) / CLOCKS_PER_SEC;
			run_Time[figNodeNum - pq_sfcSet.size() / 5 - 1] += duration;
		}
		m_database.SFCMap.clear();
		ithSFC++;
	}

	double bp=(double)blockTime1/totalSFC;//������
	blockTime1=0;
	blockTime2 = 0;
	ofstream outfile;
	outfile.open("BP-N-kun-a=0.8-ILP.xls", ios::app);//ios::out

	if(!outfile)
	{
		cerr << "��ȡ���ݵ��ļ��򲻿�" << endl;
	}
	else
	{
		outfile<<bp<<"\t"<<endl;
		outfile.close();
	}
	while(!pq_getPhyNodeSer.empty())
	{
		pq_getPhyNodeSer.top();
		pq_getPhyNodeSer.pop();
	}
	while(!pq_phyNodeSer.empty())
	{
		pq_phyNodeSer.top();
		pq_phyNodeSer.pop();
	}
}

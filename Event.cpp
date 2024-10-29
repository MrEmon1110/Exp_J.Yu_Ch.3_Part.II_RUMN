//Service.cpp : 定义控制台应用程序的入口点。
//
#include <fstream>
#include <iostream>
#include<stdio.h>
#include<time.h>
#include<algorithm>
using namespace std;
#include "Param.h"


#include <math.h>
#include <stdlib.h>   //随机数产生头文件
#include <map>
#include <vector>

#include "Database.h"
#include "Event.h"
vector<int> GenerateEvent::randperm(int Num)
{
	vector<int>temp;
	for (int i = 0; i < Num; i++)
	{
		temp.push_back(i);
	}
	random_shuffle(temp.begin(), temp.end());
	return temp;
}
int GenerateEvent::frand(int ia,int ib)
{
	double fr;	
	for(int i = 0; ; ++i)
	{	
		fr = (double)rand()*1.0/(double)RAND_MAX;
		if (fr != 0 && fr != 1)
		{
			break;
		}
		else
		{
			continue;
		}
	}
	//cout<<" frand = "<<int((ib+1-ia)*fr+ia)<<endl;
	return int((ib+1-ia)*fr+ia);	
}


double GenerateEvent::arrive_delta_time_gen (double lambda)
{
	double u=0.0, x =0.0, ln = 0.0;
	for(int i = 0; ; ++i)
	{	
		u = (double)rand()*1.0/(double)RAND_MAX;
		if (u != 0 && u != 1)
		{
			break;
		}
		else
		{
			continue;
		}
	}
	ln = log(u);
	x = ln/lambda;
	x= -1*x;
	return x;
}
double GenerateEvent::hold_delta_time_gen(double mu)
{
	double u=0.0, x =0.0, ln = 0.0;
	for(int i = 0; ; ++i)
	{	
		u = (double)rand()*1.0/(double)RAND_MAX;
		if (u != 0 && u != 1)
		{
			break;
		}
		else
		{
			continue;
		}
	}
	ln = log(u);
	x = ln/mu;
	x= -1*x;
	return x;
}

void GenerateEvent::randomSrcDst(int &src, int &des)
{
	src = frand(0, MaxNodeNum-1);
	do
	{
		des = frand(0, MaxNodeNum-1);
	}
	while (des == src);	
}

int GenerateEvent::GBPS()//每个业务的SLOTS需求数
{
	//return frand(Slot1,Slot2);

	int slot1 = GBLow;
	int slot2 = GBHigh;
	int gbps = 0;
	double fr = 0;
	for(int i = 0; ; ++i)
	{	
		fr = (double)rand()/RAND_MAX;
		if (fr != 0 && fr != 1)
		{
			gbps = (int)(slot2+1-slot1)*fr+slot1;
			break;
		}
	}
	//cout<< slotNum + GB <<endl;
	return gbps;	
}

EventMember GenerateEvent::generateIdTimeType(int idConReqNum)
{
	//double serviceArriveTime ;
	//double serviceFinishTime;

	if (idConReqNum == 0)
	{
		m_database.serviceArriveTime = 0;
		m_database.serviceFinishTime = 0;
		//srand(1);//随机种子
	}

	int sourceNode;
	int destNode;
	EventMember eachEvent;
	//srand(1);//随机种子
	double arriveDeltaTime = arrive_delta_time_gen(lambda); //每一个连接请求的时间间隔
	double holdDeltaTime = hold_delta_time_gen(mu); //每一个连接请求的时间间隔

	m_database.serviceArriveTime += arriveDeltaTime;
	m_database.serviceFinishTime = m_database.serviceArriveTime + holdDeltaTime;

	eachEvent.m_eventTime = m_database.serviceArriveTime;
	
	eachEvent.m_arriveTime = m_database.serviceArriveTime;
	eachEvent.m_finishTime = m_database.serviceFinishTime;

	eachEvent.m_eventType = Arrival;

	eachEvent.m_id = idConReqNum;

	randomSrcDst(sourceNode,destNode);

	eachEvent.m_sourceNode = sourceNode;
	eachEvent.m_destNode = destNode;

	eachEvent.m_GBps = GBPS();
	
/*
    int gbps = eachEvent.m_GBps;
	int slotNum = 0;
	int val = (int)gbps/lineRate;
	if (gbps/lineRate > (double)val )
	{
		slotNum = (int)gbps/lineRate + 1;
	}
	else
	{
		slotNum = (int)gbps/lineRate;
	}
*/

	eachEvent.failProReq = 0;

	double random = 0.0;
	double failPr = 0.0;
	while(1)
	{
		random = (double)rand()/RAND_MAX;
		//if (random >= 0.025 && random <= 0.035)
		if (random >= CRFP1 && random <= CRFP2)
		{
			failPr = random;
			break;
		} 
	}
	//cout<<"failPro "<<failPr<<endl;
	eachEvent.failProReq = -log10(1.0-failPr);
	
	//cout<<"logEachEven "<<eachEvent0.failProReq<<endl;

	Event = eachEvent;
	return Event;
}

void GenerateEvent::ERtopology(int SFC)
{
	int i = 0, j = 0;
	int n[MaxNodeNum];//度
	int k[MaxNodeNum];
	double p[MaxNodeNum];//度分布

	int totalEdges = 0;
	int addEdges = 0;
	int leftEdges = 0;
	int position = 0;
	int getEdges = 0;

	set<int> selectEdges;
	vector<int> fuc;
	if ((int)selectEdges.size() != 0)
	{
		selectEdges.clear();
	}


	for (i = 0; i < MaxNodeNum; i++)
	{
		k[i] = 0;
		n[i] = 0;
		p[i] = 0.0;
	}
	for (int i = 0; i < MaxNodeNum; i++)
	{
		for (int j = i; j < MaxNodeNum; j++)
		{
			if (i != j)
			{
				virWeight[i][j] = virWeight[j][i] = INFINITE;
			}
			else
			{
				virWeight[i][j] = virWeight[j][i] = 0;
			}
		}
	}

	//可重构
	if (m0 > 2)
	{
		totalEdges = (virMaxNode - m0) * m;
	}
	else if (m0 == 1)
	{
		totalEdges = (virMaxNode - m0) * m;
	}
	else if (m0 == 2)
	{
		totalEdges = 1 + (virMaxNode - m0) * m;
	}


	if (totalEdges < virMaxNode - 1)	//很疑惑，，表示生成的边数不足以形成一个连通图（至少需要 virMaxNode - 1 条边才能保证连通性）
	{
		cout << "The total edges is less than virMaxNode - 1. Tis VON failures." << endl;
	}
	//cout<<"totalEdges "<<totalEdges<<endl;
	map<int, vector<int> > ithNodePair;
	position = 0;
	for (int i = 0; i < virMaxNode; i++)
	{
		for (int j = i; j < virMaxNode; j++)
		{
			if (j == i + 1)		//表示这两个节点之间有一条边
			{
				virWeight[i][j] = 1;
				//weight[j][i] = 1;
				addEdges += 1;  //加了两条
			}
		}
	}
	//设置为1的距离为Gbps
	getEdges = 0;
	for (int i = 0; i < virMaxNode; i++)
	{
		for (int j = i; j < virMaxNode; j++)
		{
			if (i != j && virWeight[i][j] == 1)		//随机生成一个边权重
			{
				virWeight[j][i] = virWeight[i][j] = frand(GBLow, GBHigh);
				getEdges++;
			}
		}
	}
	if (getEdges != totalEdges)
	{
		cout << "The generated VON failures!!" << endl;
	}
	virLinkSeral NS;
	virNodeSeral virReqComRes, aVir;
	sfcSet ithSFC;

	int indexLink = 0;
	fuc = randperm(5);
	for (int i = 0; i < virMaxNode; i++)
	{
		for (int j = 0; j < virMaxNode; j++)
		{
			if (i != j && virWeight[i][j] < INFINITE)
			{
				NS.sfc_id = SFC;
				NS.virWeight = virWeight[i][j];
				NS.link_id = indexLink++;
				NS.source = i;
				NS.destination = j;
				pq_virLinkSer.push(NS);
			}
		}
		int Req = 0;
		virReqComRes.sfc_id = SFC;
		virReqComRes.virNetFuc_id = i;
		virReqComRes.reqComRes = frand(virNodeReqComResLow, virNodeReqComResHigh);
		virReqComRes.funC = fuc[i];
		switch (fuc[i])
		{
		case 0:
		{
			virReqComRes.instantRes = 20;
			break;
		}
		case 1:
		{
			virReqComRes.instantRes = 25;
			break;
		}
		case 2:
		{
			virReqComRes.instantRes = 30;
			break;
		}
		case 3:
		{
			virReqComRes.instantRes = 35;
			break;
		}
		case 4:
		{
			virReqComRes.instantRes = 40;
			break;
		}
		/*case 5:
		{
			virReqComRes.instantRes = 6;
			break;
		}
		case 6:
		{
			virReqComRes.instantRes = 7;
			break;
		}
		case 7:
		{
			virReqComRes.instantRes = 8;
			break;
		}
		case 8:
		{
			virReqComRes.instantRes = 9;
			break;
		}
		case 9:
		{
			virReqComRes.instantRes = 10;
			break;
		}*/
		}
		pq_virNodeSer.push(virReqComRes);
	}
	ithSFC.SFC_id = SFC;
	ithSFC.pq_sfcLinkSet = pq_virLinkSer;
	ithSFC.pq_sfcNodeSet = pq_virNodeSer;
	int NodeRev = 0;
	int LinkRev = 0;
	int SFCRev = 0;
	while (!pq_virLinkSer.empty())
	{
		LinkRev += 1 * pq_virLinkSer.top().virWeight;
		pq_virLinkSer.pop();
	}
	while (!pq_virNodeSer.empty())
	{
		NodeRev += 1 * pq_virNodeSer.top().reqComRes;
		pq_virNodeSer.pop();
	}
	SFCRev = NodeRev + LinkRev / 2;
	ithSFC.SFCRev = SFCRev;
	pq_sfcSet.push(ithSFC);
}


#ifndef _EVENT_H_
#define _EVENT_H_

#include<iostream>
using namespace std;
#include <stdlib.h>   //随机数产生头文件
#include <map>
#include <vector>
#include "Param.h"

class GenerateEvent
{
public:
	Database m_database;
	//Simulation m_simulation;
	EventMember Event;
	//pair<EventMember, EventMember> Event;
	int frand(int,int);
	int slots();//每个业务的SLOTS需求数
	int GBPS();//GBps
	void ERtopology(int);
	double arrive_delta_time_gen (double lambda); //每一个业务到达时间间隔，服从泊松流即到达时间为负指数分布
	double hold_delta_time_gen (double mu); //每一个业务保持时间间隔，即服务的时间间隔为负指数分布（离开时间间隔）
	void randomSrcDst(int &src, int &des);
	EventMember generateIdTimeType(int);
	vector<int> randperm(int Num);
	
	//pair<EventMember, EventMember> generateIdTimeType(int);
    //priority_queue<EventMember> Event_pq, Event_pq_tem, Event_pq_deal;
	//priority_queue<EventMember> Event_pq_tem;
};

#endif

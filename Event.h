#ifndef _EVENT_H_
#define _EVENT_H_

#include<iostream>
using namespace std;
#include <stdlib.h>   //���������ͷ�ļ�
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
	int slots();//ÿ��ҵ���SLOTS������
	int GBPS();//GBps
	void ERtopology(int);
	double arrive_delta_time_gen (double lambda); //ÿһ��ҵ�񵽴�ʱ���������Ӳ�����������ʱ��Ϊ��ָ���ֲ�
	double hold_delta_time_gen (double mu); //ÿһ��ҵ�񱣳�ʱ�������������ʱ����Ϊ��ָ���ֲ����뿪ʱ������
	void randomSrcDst(int &src, int &des);
	EventMember generateIdTimeType(int);
	vector<int> randperm(int Num);
	
	//pair<EventMember, EventMember> generateIdTimeType(int);
    //priority_queue<EventMember> Event_pq, Event_pq_tem, Event_pq_deal;
	//priority_queue<EventMember> Event_pq_tem;
};

#endif

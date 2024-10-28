#ifndef _Simulation_H_
#define _Simulation_H_


#include "Database.h"
#include "Event.h"
#include <vector>
#include <queue>
#include <map>
#include <time.h>


class Simulation
{
public:
	double totalGB;
	double energy;
	double priEnergy;
	double proEnergy;
	double tranEnergy;
	double regenEnergy;

	int spectrumNum;
	int priSpectrumNum;
	int proSpectrumNum;

	int tranNum;
	int regenNum;
	int priRegenNum;
	int proRegenNum;
	int hopNum;
	int priHopNum;
	int proHopNum;


	Suurballe m_suurballe;
	Graph m_graph;

	vector<int> mapVirNode, mapPhyNode;
	//map<nodePair, nodePairElement> allNodePair;

	//map<nodePair, pathElement> allPathElement;


	//ILPMode m_ILPMode;
	Database m_database;
	GenerateEvent m_generateEvent;
	//EventMember m_eventMember;
	//EventPair Event;
	void initialize();
	bool dealWithEvent1(EventMember);
	bool dealWithEvent(EventMember);
	bool getTraffic();
	void run(int);
	void minMatrix();
	void mapping(int,int);
	int m_sumOfFailedService;
private:
    int ithConRes;
	//double serviceArriveTime;
	//double serviceFinishTime;

	int m_sumOfBreak;
	int m_sumOfService;
	int m_sumOfBrokenService;
	int m_sumOfTerminatedService;
	int m_failedByRouting;
	int m_failedByRerouting;
	int m_sumOfWorkSlotNumber;
	int m_sumOfBackupSlotNumber;
	double m_utilizationRate;
	int m_statisticNumber;
};


#endif
#ifndef __PARAM__H__
#define __PARAM__H__

#include <queue>
#include <map>
#include <string>

using namespace std;
//event

enum Event_Type { Arrival, End, Fail, Fecover};

const int m0 = 1;
const int m = 1;

//main
const int MaxNodeNum = 14;
const int WaveNumber = 50;
const int figNodeNum = 15;
const int K_SP_Time = 1;
//const int runTime = 100;
const int runTime = 400;
const int GBNum = 1;
const int EACHCRNum = 50;
const int CRNUM = 10;
const int GAP = 10;
const int iniSFC =5;
const int SFCGap = 5;
//const int totalVON = 65;
const int totalSFC = 75;
const int CoreNumber = 7;
//const double a=1.0/3;
const double b=0.8;
const double c=0.2;
const int  computingResourceLow = 400;
const int  computingResourceHigh = 700;
//const int  computingResource = 100;
//graph
const int INFINITE = 100000;
const double INFValue = 10000000000000;
const double LFP = 0.001;
const double CRFP1 = 0.015;
const double CRFP2 = 0.025;
//database
extern int  printer; 
extern int eachCRNum;
extern int ConnectionRequestNum;
extern int rejectConnectionRequestNum;
extern int acceptNumID;
extern int acceptNum;
extern double lambda;
const double mu = 0.01;

const int Free = 0;
const int Busy = 1;
const int BusyProtection = 2;
const int BusyPrimary = 1;
const int InitializeVale = -1;

const int GBLow = 300;
const int GBHigh = 600;

const int virNodeReqComResLow = 5;
const int virNodeReqComResHigh = 10;

extern double minNodesWeight;
extern int virMaxNode;
extern int virTopoTime;
extern double minEnergy;
extern int SFCnum;

//const int GB = 1;
extern double maxReach_KM;
extern double tranUnitEnergy;
extern double regenUnitEnergy;
extern double tranUnitCost;
extern double regenUnitCost;

extern int specNum;//12.5GHZ
extern int lineRate;
//extern double channelWidth;
extern double distaceUnite;
extern int channelWidth;

extern int K_SP;
extern int conNumBreakEvent;
extern int ithTime;
extern int totalSpectrumResource;
//extern double totalSpectrum;
extern int temRegenPri;
extern int TRANumPRI;
extern int TRANumPRO;
extern int REGNumPRI;
extern int REGNumPRO;
extern int TRANum;
extern int REGNum;
extern double WASTEBAND;

extern int HOP;
extern double WEIGHT;

extern double TRANum400_1500;
extern double TRANum100_1700;
extern double TRANum40_1800;

extern int REGNum400_1500;
extern int REGNum100_1700;
extern int REGNum40_1800;

extern double GBbs;
extern double totalGBbs;

extern double receiveGBbs;
extern double RejectGBbs;

extern int priSpectrum;
extern int proSpectrum;
extern double priEnergy;
extern double proEnergy;
extern double priCost;
extern double proCost;

extern int totalSpectrum;
extern int totalSpectrum1;
extern double totalEnergy;
extern int utilizationTime;
extern double spectrUtilization;

extern int start_ksp;
extern int start_conNum;
extern int GB;

extern int blockTime1;
extern int blockTime2;
extern int blockTime3;
extern int blockTime4;
extern int blocka[figNodeNum];
extern double run_Time[figNodeNum];

extern int start_ksp;
extern int start_conNum;
extern int GB;


extern double virWeight[MaxNodeNum][MaxNodeNum];
extern double virUpdateWeight[MaxNodeNum][MaxNodeNum];
extern double virUpdateWeight[MaxNodeNum][MaxNodeNum];
extern double virNodes[MaxNodeNum][MaxNodeNum];
extern double weight[MaxNodeNum][MaxNodeNum]; 
extern double updateWeight[MaxNodeNum][MaxNodeNum];
extern double NodeInf[MaxNodeNum];

extern double RegeUpdateWeight[MaxNodeNum][MaxNodeNum];

class eachWaveLength
{
public:
	int workConReqId;
	int protectConReqId;
	int lineRate;
	int waveLengthState;
};
class  VLink
{
public:
	int slots;
	int avail_WaveLength;
	int selectCore;
	vector<int> primaryPath;
};

class nodePairElement
{
public:
	//pair<int,int> nodePair;
	
	double totalWeigth;
	int source;
	int destination;
	double priWeigth;
	double proWeigth;
	vector<int> primaryPath;
	vector<int> protectPath;
	int slots;
	map<int, vector<eachWaveLength>> CoreWaveLength;
	bool operator<(const nodePairElement& other) const
	{
		if(totalWeigth > other.totalWeigth)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};
class nodeSlots
{
public:
	int slots;
	bool operator<(const nodeSlots& other) const
	{
		if(slots > other.slots)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

};
class EventMember
{
public:
	double m_eventTime;
	double m_arriveTime;
	double m_finishTime;
	Event_Type m_eventType;
	int m_id;
	int m_sourceNode;
	int m_destNode;
	int m_GBps;
	int m_slots;
	double failProReq;

	bool operator<(const EventMember& other) const
	{
		if(m_id > other.m_id)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

class virLinkSeral
{
public:
	double virWeight;
	int sfc_id;
	int link_id;
	int source;
	int destination;

	bool operator<(const virLinkSeral& other) const
	{
		if(virWeight < other.virWeight)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

};

class VnfNodeIm
{
public:
	int totalComRes;
	int virNetFuc_id;
	int sfc_id;
	double vnfNodeIm;
	int funC;
	bool operator<(const VnfNodeIm& other) const
	{
		if(vnfNodeIm < other.vnfNodeIm )
		{
			return true;
		}
		else
		{
			return false;
		}
	}

};
class vnfNodeSeral
{
public:
	int totalComRes;
	int reqComRes;
	int instantRes;
	int virNetFuc_id;
	int sfc_id;
	int funC;
	bool operator<(const vnfNodeSeral& other)const
	{
		if (totalComRes < other.totalComRes)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};
class virNodeSeral
{
public:
	int reqComRes;
	int virNetFuc_id;
	int sfc_id;
	int instantRes;
	int funC;
	bool operator<(const virNodeSeral& other) const
	{
		if(reqComRes < other.reqComRes)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

};


class sfcSet
{
public:
	int SFC_id;
	int SFCRev;
	//virLinkSeral vonLinkSet;
	//virNodeSeral vonNodeSet;
	priority_queue<virLinkSeral> pq_sfcLinkSet;
	priority_queue<virNodeSeral> pq_sfcNodeSet;

	bool operator<(const sfcSet& other) const
	{
		if(SFC_id > other.SFC_id)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

};

class phyNodeSeral
{
public:
	int comRes;
	int phyNode_id;
	int totalComRes;
	int pfunC;
	map<int, virNodeSeral> vonSer;
	bool operator<(const phyNodeSeral& other) const
	{
		if(comRes < other.comRes)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

};

class phyNodeIm
{
public:
	int comRes;
	int phyNode_id;
	double phynodeIm;
	map<int, virNodeSeral> vonSer;
	bool operator<(const phyNodeIm& other) const
	{
		if(phynodeIm < other.phynodeIm)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

};


extern priority_queue<EventMember> Event_pq_para,Event_pq_tem, Event_pq_deal, ILP_pq_del;
extern map<int, const char *> nodeString;
extern priority_queue<virLinkSeral> pq_virLinkSer;
extern priority_queue<virNodeSeral> pq_virNodeSer;
extern priority_queue<sfcSet> pq_sfcSet, pq_eachVonSolution;
extern priority_queue<phyNodeSeral> pq_phyNodeSer, pq_phyNodeVirSer,pq_getPhyNodeSer,pq_getPhyNodeSer1;



#endif


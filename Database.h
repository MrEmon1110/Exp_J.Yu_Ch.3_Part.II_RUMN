#ifndef _DB_H_
#define _DB_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <queue>

#include "Graph.h"
#include "Suurballe.h"
#include "Param.h"

using namespace lch::graph;

using namespace std;

class paraWavelength
{
public:
	int waveLengthNum0;
	int waveLengthNum1;
	int waveLengthNum2;
	bool operator<(const paraWavelength& other) const
	{
		if(waveLengthNum0< other.waveLengthNum0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

class paraEnergy
{
public:
	int regNum;
	int transNum;
	int specNum;
	int lineRate;
	double channelWidth;
	int waveLengthNum;
	double energy;
	bool operator<(const paraEnergy& other) const
	{
		if(lineRate < other.lineRate)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

class chooseWaveLength
{
public:

	int waveLengthNum;
	int waveLengthNum400;
	int waveLengthNum100;
	int waveLengthNum40;

	bool operator<(const chooseWaveLength& other) const
	{
		if(waveLengthNum > other.waveLengthNum)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};




class NodeSer
{
public:
	int nodeNum;
	int SFCRev;

	int regenNum;
	int regenNum400;
	int regenNum100;
	int regenNum40;

	int transNum;
	int transNum400;
	int transNum100;
	int transNum40;

	double energy;
	double energy400;
	double energy100;
	double energy40;

	int specNum;
	int specNum400;
	int specNum100;
	int specNum40;

	int slots400;
	int slots100;
	int slots40;

	int lineRate400;
	int lineRate100;
	int lineRate40;

	double hop;

	int waveLengthNum;
	int waveLengthNum400;
	int waveLengthNum100;
	int waveLengthNum40;

	double channelWidth400;
	double channelWidth100;
	double channelWidth40;

	double totalWeigth;
	double updatetotalWeight;
	double UseRatio;
	double spectrUtilization;
	vector<int> aPath, primaryPath,protectPath;

	bool operator<(const NodeSer& other) const
	{
		if(energy > other.energy)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

//typedef pair<virNodeSeral,phyNodeSeral> virNodeMapPhyNode;

class SFCMapPhy
{
public:
	int ithSFC;
	//VirLinkIm virLink;
	vnfNodeSeral vnfNode;
	int avail_WaveLength;
	int selectCore;

	//vonSet VONSet;
	phyNodeSeral phyNode;
	//map<int, virNodeMapPhyNode> virNodeMapToPhyNode;
};


typedef pair<int,int> Link;



//typedef list<int> KS_Path;

//class eachWaveLength
//{
//public:
//	int workConReqId;
//	int protectConReqId;
//	int lineRate;
//	int waveLengthState;
//	int protectWaveLengthTimes;
//};


class LinkState
{
public:
	map< Link, vector<eachWaveLength> > linkState;
	map< Link, map<int, vector<eachWaveLength>> > linkState2;
	vector<int> workPath; //从源节点到宿节点的路径
	vector<int> protectPath; //从源节点到宿节点的路径
	map<int,vector<int> > conReqId_workPath;  //给每一条工作路径进行编号及其对应的路径
	map<int,vector<int> > conReqId_protectPath; //给每一保护条路径进行编号及其对应的路径

	//map<int（编号第几个连接请求）,vector<int>（存路径，用整数）>
	vector<int> findAvailWaveLength;

	vector<int> eachNeiCore;
	map<int,vector<int>> eachCoreInfor;
	map<int,Link> allLink;
	map<int,Link> failLink;
	~LinkState();
};


class pathElement
{
public:
	NodeSer priNodeSer;
	NodeSer proNodeSer;
};



class Database
{
public:

	Suurballe m_suurballe;

	LinkState m_linkstate;	
	Graph graph;

	map<nodePair, nodePairElement> allNodePair;
	map<nodePair, pathElement> allPathElement;

	map <int ,VLink> vLink,vLink1;

	map<int,vector<SFCMapPhy> > SFCMap;
	map<int,map<int, vector<SFCMapPhy> > > VONMap_AllSolution;
	int chooseIthSlot(int*, int);
	int Slots(int);

	//map<int, vector<NodeSer>> allResults;//map<<ith RUN, index of ith traffic results>, ith traffic results>;
	map<nodePair, NodeSer> allResults;//map<<ith RUN, index of ith traffic results>, ith traffic results>;
    map<nodePair, double> ithRunithVONTime;//map<<ith RUN, index of ith traffic results>, ith traffic results>;
	//void InitialWeigth();
	void topology();
	void linkStateInitial();
	vector<int> Dijkstra(int,int);
	int  regen(vector<int>);
	bool FPslots(vector<int>,int,sfcSet,int,int,bool);
	NodeSer energyFuction(vector<int>, int);

	bool Suurballe_disjoitPaths(int,int,int,int);
	bool Suurballe_phyMatrix(int);
	void aSFCMapping(int, sfcSet, vector<SFCMapPhy>);
	double aVONMapping_minSolution(sfcSet, vector<SFCMapPhy>);
	void aVONAllSolutios(sfcSet);
	void getaSFCMapOrder(sfcSet);

	bool setupPrimaryPath(int,int,int,int);
	bool reservePrimaryPathSource(int,double);

	void linkStateUpdatePrimaryPath(int);
	void linkStateUpdateProtectionPath(int);

	bool setupProtectionPath(int,int,int,int);
	bool reserveProtectionPathSource(int,double);

	void relievePrimaryPathSource(int);
	void relieveProtectionPathSource(int);

	void relievePrimaryProctionSource(int);
	
	void resourceUtilizationRatio();
	void statisticResult();

	//int vertexNum;  // nodes' number in a graph.
	//double weight[MaxNodeNum][MaxNodeNum]; // The adjoint matrix 
	//double updateWeight[MaxNodeNum][MaxNodeNum]; // The adjoint matrix
	
	int ithConRes;
	double serviceArriveTime;
	double serviceFinishTime;

	int connectionRequestNum;
	//int rejectConnectionRequestNum;
	int protectionResource;
	int primaryResource;
	int RUR;//resourceUtilizationRatio
	int StaTimeNum; //statisticTotalNumber
	int phyVnfPair[MaxNodeNum];
	int InstantCost;
	//int connetionPrimaryProtection;//代表工作路径分配资源和保护路径的计算及分配资源。0表示工作路径分配资源成功，可是计算保护路径；1表示工作路径分配资源失败，不分配保护路径的资源
	priority_queue< NodeSer > Path_pq,Path_pq_Pri,Path_pq_Pro;
};

#endif

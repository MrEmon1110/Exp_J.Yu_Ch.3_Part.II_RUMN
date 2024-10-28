#include "Param.h"
#include "Database.h"


int virMaxNode;
int virTopoTime = 100;
double minNodesWeight = INFValue;
double minEnergy = 10000000000000;
int SFCnum = 0;
//database
int zww = 0;
int  printer = 0; 
int eachCRNum = 1;
int ConnectionRequestNum = 1;
int rejectConnectionRequestNum = 0;
int acceptNumID = 0;
int acceptNum = 0;
int K_SP = 3;
double lambda = 0.5;
int conNumBreakEvent = 5;
//int conNumBreakGraph = 0;
int ithTime = 0;

int totalSpectrumResource = 0;

double maxReach_KM = 1800.0;
double tranUnitEnergy = 167.0;
double regenUnitEnergy = 334.0;
double tranUnitCost = 2.5;
double regenUnitCost = 5.0;

int lineRate = 40;
//int channelWidth = 2;//12.5GHZ
double distaceUnite = 0.0001;

int temRegenPri = 0;
int TRANumPRI = 0;
int TRANumPRO = 0;
int REGNumPRI = 0;
int REGNumPRO = 0;
int TRANum = 0;
int REGNum = 0;
double WASTEBAND = 0;

int HOP = 0;
double WEIGHT = 0;

double TRANum400_1500;
double TRANum100_1700;
double TRANum40_1800;

int REGNum400_1500;
int REGNum100_1700;
int REGNum40_1800;

double GBbs = 1.0;
double totalGBbs = 0;
int channelWidth=2;///////////////////////////

double receiveGBbs = 0;
double RejectGBbs = 0;

int priSpectrum = 0;
int proSpectrum = 0;
double priEnergy = 0;
double proEnergy = 0;
double priCost = 0;
double proCost = 0;
int totalSpectrum = 0;
int totalSpectrum1 = 0;
double totalEnergy = 0;
double totalCost = 0;
int utilizationTime = 0;
double spectrUtilization = 0;


int start_ksp = 0;
int start_conNum = 0;
int GB = 0;

int blockTime1 = 0;
int blockTime2 = 0;
int blockTime3 = 0;
int blockTime4 = 0;
int blocka[figNodeNum] = { 0 };
double run_Time[figNodeNum] = { 0 };


double virWeight[MaxNodeNum][MaxNodeNum];
double virUpdateWeight[MaxNodeNum][MaxNodeNum];
double virNodes[MaxNodeNum][MaxNodeNum];
priority_queue<virLinkSeral> pq_NodeSer;

double weight[MaxNodeNum][MaxNodeNum]; 
double updateWeight[MaxNodeNum][MaxNodeNum];
double NodeInf[MaxNodeNum];
double RegeUpdateWeight[MaxNodeNum][MaxNodeNum];

priority_queue<EventMember> Event_pq_para,Event_pq_tem, Event_pq_deal, ILP_pq_del;
map<int, const char *> nodeString;

priority_queue<virLinkSeral> pq_virLinkSer;
priority_queue<virNodeSeral> pq_virNodeSer;
priority_queue<sfcSet> pq_sfcSet, pq_eachVonSolution;

priority_queue<phyNodeSeral> pq_phyNodeSer,pq_phyNodeVirSer,pq_getPhyNodeSer,pq_getPhyNodeSer1;



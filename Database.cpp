// Service.cpp : 定义控制台应用程序的入口点。
//
//#include <stdio.h>
#include "Param.h"

#include <fstream>
#include <stdio.h> //C中的FILE *fp头文件
#include <iomanip>
//#include <map>
#include <vector>
#include <set>
#include <map>
#include<math.h>

//#include <iomanip>
#include "Database.h"

#include <iostream>
using namespace std;

typedef list<int> pathList;			// the top is the first node
typedef vector<pathList> pathVector;		// the paths are sorted by ASC
int o = 0;
LinkState::~LinkState()
{
	linkState.clear();
	workPath.clear();
	protectPath.clear();
	conReqId_workPath.clear();
	//conReqId_protectPath.clear;
	findAvailWaveLength.clear();
	allLink.clear();
	failLink.clear();
}
void Database::linkStateInitial()
{
	int node1, nodeConWithNode1;
	int core, neicore;
	vector<eachWaveLength> eachLinkWaveLength;
	eachWaveLength waveLengthInfor;
	map<int, vector<eachWaveLength>> ithCoreWaveLength;
	int eachWaveNum = 0;
	while (eachWaveNum < WaveNumber)//频谱隙状态初始化
	{
		waveLengthInfor.workConReqId = InitializeVale;
		waveLengthInfor.waveLengthState = Free;
		//waveLengthInfor.protectWaveLengthTimes = 0;
		eachLinkWaveLength.push_back(waveLengthInfor);
		eachWaveNum++;
	}
	m_linkstate.linkState2.clear();
	for (core = 0; core < CoreNumber; core++)
	{
		for (neicore = 0; neicore < CoreNumber; neicore++)
		{
			if (core != CoreNumber - 1 && (abs(core - neicore) == 1 || neicore == CoreNumber - 1 || abs(core - neicore) == CoreNumber - 1 - 1))
			{
				m_linkstate.eachNeiCore.push_back(neicore);
			}
			if (core == CoreNumber - 1 && core != neicore)
			{
				m_linkstate.eachNeiCore.push_back(neicore);

			}
		}
		while (m_linkstate.eachNeiCore.size() != 0)
		{
			m_linkstate.eachCoreInfor.insert(make_pair(core, m_linkstate.eachNeiCore));//每个核心相邻的核的索引
			m_linkstate.eachNeiCore.clear();
		}
	}
	for (core = 0; core < CoreNumber; core++)
	{
		ithCoreWaveLength.insert(make_pair(core, eachLinkWaveLength));
	}

	for (node1 = 0; node1 < MaxNodeNum; node1++)
	{
		for (nodeConWithNode1 = 0; nodeConWithNode1 < MaxNodeNum; nodeConWithNode1++)
		{
			if (weight[node1][nodeConWithNode1] != 0 && weight[node1][nodeConWithNode1] < INFINITE)//如果两个节点间直接相连
			{
				//m_linkstate.linkState.insert(make_pair(make_pair(node1, nodeConWithNode1),core1));////////////////
				m_linkstate.linkState2.insert(make_pair(make_pair(node1, nodeConWithNode1), ithCoreWaveLength));//42个链路，每个链路上都有ithCoreWaveLength（7个光芯，每个芯有50个slots,每个slot有5个属性）

			}
		}
	}
	eachLinkWaveLength.clear();
}

vector<int> Database::Dijkstra(int source, int destination)
{
	//cout<<source<<" "<<destination<<endl;
	int routeTable_hopNum = 0;
	int path[MaxNodeNum];
	int routeTable_nodeId[MaxNodeNum];
	double S[MaxNodeNum];
	double dist[MaxNodeNum];

	for (int i = 0; i < MaxNodeNum; i++)
	{
		path[i] = -1;
		routeTable_nodeId[i] = -1;
		S[i] = -1;
		dist[i] = -1;
	}

	//nodePairElement pairElement;

	for (int i = 0; i < MaxNodeNum; i++)//initialize dist[n], S[source], path[n]
	{
		dist[i] = RegeUpdateWeight[source][i];
		//cout<<updateWeight[source][i]<<endl;
		S[i] = 0;
		if (i != source && dist[i] < INFINITE)
		{
			path[i] = source;
		}
		else
		{
			path[i] = -1;
		}
	}


	S[source] = 1;      //顶点source加入已求出最短路径的顶点集合
	dist[source] = 0;      //标记顶点source为源点
	//cout << source <<"  ";
	//求出其其余n-1个顶点的最短路径
	for (int num = 0; num < MaxNodeNum - 1; num++)	 //当数组的顶点小于图的顶点数时循环
	{
		double min = INFINITE;
		int u = source;
		for (int i = 0; i < MaxNodeNum; i++) //选择当前不在集合source中具有最短路径的顶点u
		{
			if (!S[i] && dist[i] < min)
			{
				u = i;
				min = dist[i];
			}
		}
		//cout << u <<"  ";
		S[u] = 1;  //顶点u加入已求出最短路径的顶点集合
		for (int i = 0; i < MaxNodeNum; i++)//修改其余顶点的当前最短路径
		{
			//if ((!S[i]) && updateWeight[u][i]!=0 && updateWeight[u][i]< INFINITE && (dist[u]+updateWeight[u][i]< dist[i]))
			if ((!S[i]) && RegeUpdateWeight[u][i] < INFINITE && (dist[u] + RegeUpdateWeight[u][i] < dist[i]))
			{
				dist[i] = dist[u] + RegeUpdateWeight[u][i];
				path[i] = u;
			}
		}
	}
	//cout << endl;
	// 下面是得到逆向的算路径结果，例如4<-3<-0,它表示源节点为0,目的节点为4, 对应的路径为4->3->0,是034吧
	if (destination == source)
	{
		routeTable_hopNum = 0;
		cout << "path computation fails because destination==source. " << endl << endl;
		//return false;
	}
	else
	{
		routeTable_hopNum = 0;    // hop number in ith path
		int key = destination;
		routeTable_nodeId[0] = destination;
		//cout<< routeTable_nodeId[0] << " <- ";
		vector<int> pathList;
		//map<int,nodePairElement>::iterator insertPriPath = dis_paths.find(0);

		for (int k = 0; k < MaxNodeNum; k++)
		{
			if (path[key] == source)
			{
				routeTable_hopNum = k + 1;
				routeTable_nodeId[k + 1] = source;
				//cout<<routeTable_nodeId[k+1]<<endl;
				//cout << "hopNum = "<<routeTable_hopNum[i] <<endl;
				//从源节点到目的节点连接成一条路径
				for (int num = routeTable_hopNum; num >= 0; num--)
				{
					//cout << routeTable_nodeId[num]<<" ";
					int node = routeTable_nodeId[num];
					pathList.push_back(node);
				}
				//cout<<endl;
				//cout<<"Setup path2 success"<<endl;
				return pathList;
				break;
			}
			else
			{
				routeTable_nodeId[k + 1] = path[key];
				key = path[key];
				//cout<<routeTable_nodeId[k+1]<<endl;
				if (routeTable_nodeId[k + 1] == -1)
				{
					cout << " Setup path21 fail " << endl << endl;
					//return false;
					//break;
				}
			}
		}
	}
}

int Database::Slots(int gbps)
{
	int waveLengthNum = 0;
	int val = (int)gbps / lineRate;
	if ((double)gbps / lineRate > val)
	{
		waveLengthNum = (int)gbps / lineRate + 1;
	}
	else
	{
		waveLengthNum = (int)gbps / lineRate;
	}
	int slots = waveLengthNum * channelWidth;
	return slots;
}
int Database::chooseIthSlot(int* P, int totalSlot)//首次命中选频隙
{
	int iTime = 0;
	int totlSpeGap = 0;
	int eachWaveNum;//每一个波长编号
	int avail_WaveLength = -1;
	int totalFreeSlots = 0;
	int chooseWay = 0;
	//int priPath[MaxNodeNum];
	int eachLinkWavaLength[WaveNumber];
	for (int i = 0; i < WaveNumber; i++)
	{
		eachLinkWavaLength[i] = *P;
		P++;
	}

	eachWaveNum = 0;
	while (eachWaveNum < WaveNumber)
	{
		totalFreeSlots = 0;
		for (int i = 0; i < totalSlot; i++)
		{
			if (eachLinkWavaLength[eachWaveNum + i] == Free)
			{
				totalFreeSlots++;
			}
		}
		if (totalFreeSlots == totalSlot)
		{
			avail_WaveLength = eachWaveNum;
			break;
		}
		eachWaveNum++;
	}

	return avail_WaveLength;
}
bool Database::FPslots(vector<int> primaryPath, int slots, sfcSet getASFC, int LinkNum, int Size, bool pd)
{//分配频谱隙
	int avail_WaveLength = -1;
	/////chuanraozengjia///////////////////
	/*int nodeA,nodeB,N;
	int BF=0,FF=0;
	double μ,h,L;
	double XT1,XT;
	double η=0.8;
	double ε=0.5;
	double K=3.16e-5;
	double R=5.5*10e-5;
	int β=4e6;
	double W=4.5e-8;
	h=(2*K*K*R)/(β*W);
	ofstream outfile;
	double TH=-32.0,MXT=0.0;
	int logal=1;
	int avail_WaveLengthcr=-1;
	int flag=0,flag1=0;
	int selectCorecr=-1;
	int forxunhuandecanshu = 0;*/
	/////////end////////////////////////////////
	int selectCore = -1;
	int z = 0;
	VLink vlink, vlink1;
	bool lp = false;
	map<int, vector<eachWaveLength>> findIthCoreAndWaveLength;
	if (!primaryPath.empty())
	{
		int priPath[MaxNodeNum];
		for (int i = 0; i < MaxNodeNum; i++)
		{
			priPath[i] = -1;
		}
		int node = 0;
		for (int i = 0; i < primaryPath.size(); i++)
		{
			priPath[node] = primaryPath.at(i);
			node++;
		}
		int eachLinkWavaLength[WaveNumber];
		for (int core = 0; core < CoreNumber; core++)
		{
			for (int eachWave = 0; eachWave < WaveNumber; eachWave++)
			{
				eachLinkWavaLength[eachWave] = Free;
			}
			for (int i = 0; i < (int)primaryPath.size() - 1; i++)//更新当前纤芯的频谱隙占用状态
			{
				int node1 = priPath[i];
				int nodeConWithNode1 = priPath[i + 1];
				if (nodeConWithNode1 == -1)
				{
					break;
				}
				map< Link, map<int, vector<eachWaveLength>> >::iterator resPriPathNeighborNode = m_linkstate.linkState2.find(make_pair(node1, nodeConWithNode1));
				findIthCoreAndWaveLength = resPriPathNeighborNode->second;
				map<int, vector<eachWaveLength>>::iterator findWaveLengthInIthCore = findIthCoreAndWaveLength.find(core);
				int eachWaveNum = 0;

				while (eachWaveNum < WaveNumber)
				{
					if (findWaveLengthInIthCore->second[eachWaveNum].waveLengthState == BusyPrimary)
					{
						eachLinkWavaLength[eachWaveNum] = BusyPrimary;
					}
					eachWaveNum++;
				}
			}
			avail_WaveLength = chooseIthSlot(eachLinkWavaLength, slots + GB);//chooseIthSlot是用首次命中法选择频谱隙
			selectCore = core;//选择当前遍历的纤芯作为准备选纤芯，后续判断是否满足交叉串扰值约束
			//
/////////////////////////////////////串扰值计算begin//////////////////////////////////////////////
//			map<int,vector<int>> :: iterator findEachneiCore= m_linkstate.eachCoreInfor.find(core);
//			//取出当前纤芯的相邻纤芯，比如和0号纤芯相邻的1.5.6
//			N=findEachneiCore->second.size();//N代表共有几个相邻纤芯
//			while(avail_WaveLength != -1 && avail_WaveLength+slots-1<WaveNumber)
//			{//当备选频谱隙在有意义的范围内时（备选频谱隙没有超出一共提供的频谱隙范围）
//				for (int j = 0; j < (int)primaryPath.size()-1; j++)
//				{//对于工作链路上的每段链路计算交叉串扰，并准备进行下一步操作
//					int node1 = priPath[j];
//					int node2 = priPath[j+1];//取出正在处理的工作链路的源目节点，用于取出源目节点距离，计算串扰值
//					if (node2 == -1)
//					{
//						break;
//					}
//					if(j==0)
//					{
//						nodeA=priPath[j];
//					}
//					if(j==primaryPath.size()-2)
//					{
//						nodeB=priPath[j+1];//nodeA、nodeB代表该工作链路的源目节点
//					}
//					for(int i=0;i<N;i++)//统计相邻纤芯准备选频谱隙总占用状态，存入BF中，用于计算串扰值
//					{
//						int neiCore = findEachneiCore->second[i];
//						map<Link,map<int,vector<eachWaveLength>>> :: iterator preFindneiCoreInfor = m_linkstate.linkState2.find(make_pair(node1,node2));
//						map<int,vector<eachWaveLength>> :: iterator findNeiCoreInfor =  preFindneiCoreInfor->second.find(neiCore);
//						for(int k=avail_WaveLength;k<avail_WaveLength+slots;k++)
//						{
//							if(findNeiCoreInfor->second[k].waveLengthState == BusyPrimary)
//							{
//								BF++;
//							}
//						}
//					}
//					//计算交叉串扰值
//					L=weight[node1][node2];
//					FF=N*slots-BF;
//					XT1=(N-N*exp(-(N-1)*2*h*L))/(1+N*exp(-(N+1)*2*h*L));
//					if(BF==0)
//					{
//						XT=0;
//					}
//					else
//					{
//						μ=(N*BF+ε)/(slots*FF+η);
//						XT=4.4*log10(XT1*μ);//得到交叉串扰值
//					}
//					BF=0;
//					//输出串扰值
//					if(XT==0)
//					{
//						XT=-44;
//					}
//					//outfile.open("chaunraozhi-xtunsetth.xls", ios::app);//ios::out
//					//if(!outfile)
//					//{
//					//	cerr << "存取数据的文件打不开" << endl;
//					//}
//					//else
//					//{
//					//	outfile<<getAVON.VON_id<<"\t"<<XT<<endl;
//					//	outfile.close();
//					//}
//				}//结束for循环，退出到while
//				break;
//			}//结束while,回到core层循环
//////////////////////////////////////////串扰值计算end////////////////////////////////////////////////
			//
			if (avail_WaveLength != -1 && selectCore != -1)//如果已经
			{
				lp = true;
				int n = 0;
				while (n < (slots + GB))
				{
					int availWaveLegth = (avail_WaveLength + n);
					for (int pathLength = 0; pathLength < (int)primaryPath.size() - 1; pathLength++)//分配波长
					{//修改链路状态m_linkstate，即给链路上的各条路径将要占用的频谱隙置为BusyPrimary
						int node1 = primaryPath.at(pathLength);
						int nodeConWithNode1 = primaryPath.at(pathLength + 1);
						map< Link, map<int, vector<eachWaveLength>> >::iterator findCoreWave = m_linkstate.linkState2.find(make_pair(node1, nodeConWithNode1));
						map<int, vector<eachWaveLength>>::iterator findWaveLengthInIthCore = findCoreWave->second.find(selectCore);
						findWaveLengthInIthCore->second[availWaveLegth].waveLengthState = BusyPrimary;
						findWaveLengthInIthCore->second[availWaveLegth].workConReqId = getASFC.SFC_id;
					}
					n++;
				}
				//shanchudedifang //
				z = LinkNum;
				vlink.avail_WaveLength = avail_WaveLength;
				vlink.selectCore = selectCore;
				vlink.slots = slots;
				vlink.primaryPath = primaryPath;
				vLink.insert(make_pair(z, vlink));//当前映射节点与已经映射过的节点，存在链路相连的数组
				vlink1.avail_WaveLength = avail_WaveLength;
				vlink1.primaryPath = primaryPath;
				vlink1.selectCore = selectCore;
				vlink1.slots = slots;
				vLink1.insert(make_pair(o, vlink1));//总的链路数组
				o++;
				if (LinkNum == Size)
				{
					vLink.clear();
				}
				break;
			}
		}

		if (lp == false)
		{
			for (int i = 0; i < vLink.size(); i++)
			{
				for (int j = 0; j < vLink1.size(); j++)
				{
					map <int, VLink>::iterator linkinfor = vLink1.find(j);
					if (linkinfor->second.primaryPath == vLink.find(vLink.size() - i - 1)->second.primaryPath)
					{
						slots = vLink1.find(j)->second.slots;
						avail_WaveLength = vLink1.find(j)->second.avail_WaveLength;
						selectCore = vLink1.find(j)->second.selectCore;
						primaryPath = vLink1.find(j)->second.primaryPath;
						int n = 0;
						while (n < (slots + GB))
						{
							int availWaveLegth = (avail_WaveLength + n);
							for (int pathLength = 0; pathLength < (int)primaryPath.size() - 1; pathLength++)//分配波长
							{
								int node1 = primaryPath.at(pathLength);
								int nodeConWithNode1 = primaryPath.at(pathLength + 1);
								map< Link, map<int, vector<eachWaveLength>> >::iterator findCoreWave = m_linkstate.linkState2.find(make_pair(node1, nodeConWithNode1));
								map<int, vector<eachWaveLength>>::iterator findWaveLengthInIthCore = findCoreWave->second.find(selectCore);
								findWaveLengthInIthCore->second[availWaveLegth].waveLengthState = Free;
								findWaveLengthInIthCore->second[availWaveLegth].workConReqId = -1;
							}
							n++;
						}
						vLink1.erase(j);
						o--;
					}
				}
			}
			vLink.clear();
		}
		if (vLink1.size() != 2 && pd)
		{
			lp = false;
			for (int i = 0; i < vLink1.size(); i++)
			{
				slots = vLink1.find(i)->second.slots;
				avail_WaveLength = vLink1.find(i)->second.avail_WaveLength;
				selectCore = vLink1.find(i)->second.selectCore;
				primaryPath = vLink1.find(i)->second.primaryPath;
				int n = 0;
				while (n < (slots + GB))
				{
					int availWaveLegth = (avail_WaveLength + n);
					for (int pathLength = 0; pathLength < (int)primaryPath.size() - 1; pathLength++)//分配波长
					{
						int node1 = primaryPath.at(pathLength);
						int nodeConWithNode1 = primaryPath.at(pathLength + 1);
						map< Link, map<int, vector<eachWaveLength>> >::iterator findCoreWave = m_linkstate.linkState2.find(make_pair(node1, nodeConWithNode1));
						map<int, vector<eachWaveLength>>::iterator findWaveLengthInIthCore = findCoreWave->second.find(selectCore);
						findWaveLengthInIthCore->second[availWaveLegth].waveLengthState = Free;
						findWaveLengthInIthCore->second[availWaveLegth].workConReqId = -1;
					}
					n++;
				}
			}
			vLink1.clear();
		}
	}
	return lp;
}

int Database::regen(vector<int> pathVector)
{
	int node1 = -1, node2 = -1;
	int regenNum = -1;
	vector<int> rePath;
	int Regen = 0;
	double Energy = 0;

	int source = (int)pathVector.at(0);
	int pathSize = pathVector.size();
	int destination = (int)pathVector.at(pathSize - 1);

	for (int i = 0; i < MaxNodeNum; i++)
	{
		for (int j = 0; j < MaxNodeNum; j++)
		{
			if (i != j)
			{
				RegeUpdateWeight[j][i] = RegeUpdateWeight[i][j] = INFINITE;
			}
			else
			{
				RegeUpdateWeight[i][j] = 0;
			}
		}
	}

	for (int i = 0; i < (int)pathVector.size(); i++)
	{
		node1 = pathVector.at(i);
		for (int j = i + 1; j < (int)pathVector.size(); j++)
		{
			Energy = 0;
			node2 = pathVector.at(j);
			for (int k = i; k < j; k++)
			{
				int node3 = pathVector.at(k);
				int node4 = pathVector.at(k + 1);
				if (weight[node3][node4] != INFINITE)
				{
					Energy += weight[node3][node4];
				}
			}
			if (Energy <= maxReach_KM)
			{
				RegeUpdateWeight[node1][node2] = 1;
			}
		}
	}
	rePath = Dijkstra(source, destination);
	Regen = rePath.size() - 2;//去头去尾，即是放置再生器的数目

	return Regen;
}
NodeSer Database::energyFuction(vector<int> pathVector, int gbps)
{
	NodeSer ithPathEnergy;
	double chooseEnergy = INFINITE;
	double varEnergy = 0;
	int lineRateType = 2;
	int chooseType = 0;

	int aLinRate = 0;
	int channelWidth = 0;
	double energy = 0;
	/*int hop=0;*/
	int spectrum = 0;
	int regNum = 0;
	int totalRegNum = 0;
	int tranNum = 0;
	int waveLengthNum = 0;
	double pathWeight = 0;

	priority_queue< chooseWaveLength > waveLengthDecending;

	chooseWaveLength waveNum;

	vector<int> R400;
	vector<int> R100;
	vector<int> R40;

	int lineRateType0 = 0;
	int lineRateType1 = 0;
	int lineRateType2 = 0;

	int waveLengthNum0 = 0;
	int waveLengthNum1 = 0;
	int waveLengthNum2 = 0;

	int waveLengthNum400 = 0;
	int waveLengthNum100 = 0;
	int waveLengthNum40 = 0;

	int bandwidthGap = 5000;


	for (int node = 0; node < (int)pathVector.size() - 1; node++)
	{
		int i = pathVector.at(node);
		int j = pathVector.at(node + 1);
		pathWeight += weight[i][j];
	}

	for (lineRateType = 2; lineRateType < 3; lineRateType++)
	{
		switch (lineRateType)
		{
		case 0:
		{
			lineRate = 400;
			int val = (int)gbps / lineRate;
			if ((double)gbps / lineRate > val)
			{
				waveLengthNum0 = (int)gbps / lineRate + 1;
			}
			else
			{
				waveLengthNum0 = (int)gbps / lineRate;
			}
			for (int i = 0; i <= waveLengthNum0; i++)
			{
				R400.push_back(i);
			}
			break;
		}
		case 1:
		{
			lineRate = 100;
			int val = (int)gbps / lineRate;
			if ((double)gbps / lineRate > val)
			{
				waveLengthNum1 = (int)gbps / lineRate + 1;
			}
			else
			{
				waveLengthNum1 = (int)gbps / lineRate;
			}
			for (int i = 0; i <= waveLengthNum1; i++)
			{
				R100.push_back(i);
			}
			break;
		}
		case 2:
		{
			lineRate = 40;
			int val = (int)gbps / lineRate;
			if ((double)gbps / lineRate > val)
			{
				waveLengthNum2 = (int)gbps / lineRate + 1;
			}
			else
			{
				waveLengthNum2 = (int)gbps / lineRate;
			}
			for (int i = 0; i <= waveLengthNum2; i++)
			{
				R40.push_back(i);
			}
			break;
		}
		}
	}
	waveLengthNum40 = waveLengthNum2;
	for (lineRateType = 2; lineRateType < 3; lineRateType++)
	{
		//lineRateType = cases[i][j];
		switch (lineRateType)
		{
		case 0:
		{
			maxReach_KM = 1500;
			tranUnitEnergy = 450.0;
			regenUnitEnergy = 900.0;
			tranUnitCost = 5.5;
			regenUnitCost = 11.0;
			channelWidth = 10;//12.5GHZ
			lineRate = 400;
			waveLengthNum = waveLengthNum400;
			break;
		}
		case 1:
		{
			maxReach_KM = 1700;
			tranUnitEnergy = 300.0;
			regenUnitEnergy = 600.0;
			tranUnitCost = 3.75;
			regenUnitCost = 7.5;
			channelWidth = 3;//12.5GHZ
			lineRate = 100;
			waveLengthNum = waveLengthNum100;
			break;
		}
		case 2:
		{
			maxReach_KM = 1800;
			tranUnitEnergy = 167.0;
			regenUnitEnergy = 334.0;
			tranUnitCost = 2.5;
			regenUnitCost = 5.0;
			channelWidth = 2;//12.5GHZ
			lineRate = 40;
			waveLengthNum = waveLengthNum40;
			break;
		}
		}

		regNum = regen(pathVector);

		switch (lineRateType)
		{
		case 0:
		{
			ithPathEnergy.lineRate400 = 400;
			ithPathEnergy.channelWidth400 = 10;
			ithPathEnergy.waveLengthNum400 = waveLengthNum;
			ithPathEnergy.energy400 = (waveLengthNum * 2.0 * tranUnitCost + waveLengthNum * regNum * regenUnitCost);
			//ithPathEnergy.lineRateType = chooseType;
			ithPathEnergy.specNum400 = ((int)pathVector.size() - 1) * waveLengthNum * channelWidth;
			ithPathEnergy.transNum400 = waveLengthNum * 2;
			ithPathEnergy.regenNum400 = waveLengthNum * regNum;
			break;
		}
		case 1:
		{
			ithPathEnergy.lineRate100 = 100;
			ithPathEnergy.channelWidth100 = 3;
			ithPathEnergy.waveLengthNum100 = waveLengthNum;
			ithPathEnergy.energy100 = (waveLengthNum * 2.0 * tranUnitCost + waveLengthNum * regNum * regenUnitCost);
			//ithPathEnergy.lineRateType = chooseType;
			ithPathEnergy.specNum100 = ((int)pathVector.size() - 1) * waveLengthNum * channelWidth;
			ithPathEnergy.transNum100 = waveLengthNum * 2;
			ithPathEnergy.regenNum100 = waveLengthNum * regNum;
			break;
		}
		case 2:
		{
			ithPathEnergy.lineRate40 = 40;
			ithPathEnergy.channelWidth40 = 2;
			ithPathEnergy.waveLengthNum40 = waveLengthNum;//      gbps/linerate 
			ithPathEnergy.energy40 = (waveLengthNum * 2.0 * tranUnitCost + waveLengthNum * regNum * regenUnitCost);
			//ithPathEnergy.lineRateType = chooseType;
			ithPathEnergy.specNum40 = ((int)pathVector.size() - 1) * waveLengthNum * channelWidth;
			ithPathEnergy.transNum40 = waveLengthNum * 2;
			ithPathEnergy.regenNum40 = waveLengthNum * regNum;
			break;
		}
		}

		/*hop+=(int)pathVector.size()-1;*/
		energy += (waveLengthNum * 2.0 * tranUnitCost + waveLengthNum * regNum * regenUnitCost);
		//ithPathEnergy.lineRateType = chooseType;
		spectrum += ((int)pathVector.size() - 1) * waveLengthNum * channelWidth;
		tranNum += waveLengthNum * 2;
		totalRegNum += (waveLengthNum * regNum);
	}
	ithPathEnergy.energy = energy;
	ithPathEnergy.specNum = spectrum;
	ithPathEnergy.transNum = tranNum;
	ithPathEnergy.regenNum = totalRegNum;
	/*ithPathEnergy.hop=hop;*/

	ithPathEnergy.totalWeigth = pathWeight;
	ithPathEnergy.nodeNum = (int)pathVector.size();
	ithPathEnergy.aPath = pathVector;
	return ithPathEnergy;
}
bool Database::Suurballe_phyMatrix(int idConnectionRequestNum) //基于dijkstra算法计算工作路径
{
	NodeSer nodeSer;
	vector<int> primaryPath, protectPath;

	allNodePair = m_suurballe.allDisPath(idConnectionRequestNum);

	/*pathElement pathPairElement;
	for (int source = 0; source<MaxNodeNum; source++)
	{
		for (int destination = 0; destination<MaxNodeNum; destination++)
		{
			if (source != destination )
			{
				map<nodePair, nodePairElement>::iterator findDisPath = allNodePair.find(make_pair(source,destination));
				if (findDisPath != allNodePair.end())
				{
					primaryPath = findDisPath->second.primaryPath;
					protectPath = findDisPath->second.protectPath;

					NodeSer ithPriPathEnergy = energyFuction(primaryPath,);
					NodeSer ithProPathEnergy = energyFuction(protectPath);
					pathPairElement.priNodeSer = ithPriPathEnergy;
					pathPairElement.proNodeSer = ithProPathEnergy;
					allPathElement.insert(make_pair(make_pair(source,destination), pathPairElement));
					//cout<<"source "<<source<<endl;
				}
			}
		}
	}*/
	return true;
}
void Database::getaSFCMapOrder(sfcSet getASFC)
{
	int ithSFC = 0;
	bool status;
	sfcSet getSFC;
	virNodeSeral virNode;
	phyNodeSeral phyNode;
	virLinkSeral virLink;
	vnfNodeSeral vnfNode;
	VnfNodeIm vnfNode1;
	SFCMapPhy sfcMapToPhy;
	vector<SFCMapPhy> aSFCMapSet;
	vector<int> vnfNodeVector, eachNodeVect, vnfNodeVector1, vnfNodeVector2;
	priority_queue<virNodeSeral>  getVirNode1;
	priority_queue<virLinkSeral>  getVirLink1;
	priority_queue<VnfNodeIm>  getVnfNodeIm;
	priority_queue<vnfNodeSeral>getVnfNode, getVnfNode1, getVnfNode2, getVnfNode3, getVnfNode4;
	map<nodePair, virLinkSeral> allVirLink;

	map<int, double> VirNodeAB;
	map<int, virNodeSeral> ithVirNode;
	map<int, vnfNodeSeral>ithVnfNode, vnfNodeRes;
	map<int, virLinkSeral> ithvirLink;
	map<int, VnfNodeIm> virNodeMapOrder;
	map<int, int> vnfNodeOrder;
	map<int, vector<int>> funCphyVnf, phyfunCvnf;
	vector<int>unfunCphyVnf, unfunCphyVnf1;
	vector<int>funCphyNode;
	/////////////////////////////////////////////////////////////////////
	vLink.clear();
	vLink1.clear();
	o = 0;//o是vlink1总的链路数组前面的链路编号 
	InstantCost = 0;
	int avail_WaveLength = -1;
	int core = 0;
	int selectCore = -1;
	int ith = 0;
	int ithLink = 0;
	int ithMapNum = 0;

	bool lp = false;//  lp是干什么用的?
	int availWaveLegth;
	//////////////////////////////////////////////////////////////////////////


	int eachLinkWavaLength[WaveNumber];
	int priPath[MaxNodeNum];
	map<int, vector<eachWaveLength>> findIthCoreAndWaveLength;
	map<int, virLinkSeral> virLinkMapOrder1;
	priority_queue<vnfNodeSeral> vnfVONSet;

	nodePairElement phyNode1;
	phyNodeIm phyNodeim;
	nodeSlots phyNode2;
	priority_queue<nodePairElement> fixNodePair1;
	priority_queue<nodeSlots> fixNodePair2;
	priority_queue<phyNodeIm> eachNodePair, fixNodePair;
	map<int, nodePairElement>eachNodePair1;
	map<int, nodeSlots>eachNodePair2;
	map<int, phyNodeSeral> phyNodeRes, phyNodeRes1;
	map<int, SFCMapPhy> mapPhyNode;//每个SFC业务进来的时候，都会重新生成一遍，所以不用clear
	map<int, phyNodeSeral> unMapPhyNode, haveMapPhyNode;

	map<int, vector<SFCMapPhy> > SFCMap1, chooseSFCMap;


	//////////////////////////////////////////////////////////////////////////
	if ((int)aSFCMapSet.size() != 0)
	{
		aSFCMapSet.clear();
	}
	if ((int)vnfNodeVector.size() != 0)
	{
		vnfNodeVector.clear();
	}

	while (!getVnfNodeIm.empty())
	{
		getVnfNodeIm.top();
		getVnfNodeIm.pop();
	}
	if (!virNodeMapOrder.empty())
	{
		virNodeMapOrder.clear();
	}
	if (!vnfNodeOrder.empty())
	{
		vnfNodeOrder.clear();
	}
	pq_phyNodeVirSer = pq_phyNodeSer;
	sfcMapToPhy.ithSFC = getASFC.SFC_id;

	getSFC = getASFC;
	int m = 0;
	while (!getSFC.pq_sfcLinkSet.empty())
	{
		virLink = getSFC.pq_sfcLinkSet.top();
		ithvirLink.insert(make_pair(m, virLink));
		m++;
		getSFC.pq_sfcLinkSet.pop();
	}///////////////////////////////////////////////在虚拟链路前加入key值
	int n = 0;
	while (!getSFC.pq_sfcNodeSet.empty())
	{
		virNode = getSFC.pq_sfcNodeSet.top();
		ithVirNode.insert(make_pair(n, virNode));
		n++;
		getSFC.pq_sfcNodeSet.pop();
	}///////////////////////////////////////////////在虚拟节点前加入key值
	getSFC = getASFC;
	while (!getSFC.pq_sfcLinkSet.empty())
	{
		int virSou;
		int virDes;
		//virLinkSeral virLink;
		virLink = getSFC.pq_sfcLinkSet.top();
		virSou = virLink.source;
		virDes = virLink.destination;
		allVirLink.insert(make_pair(make_pair(virSou, virDes), virLink));//虚拟链路前添加了一个key值
		getSFC.pq_sfcLinkSet.pop();
	}
	getSFC = getASFC;
	while (!getSFC.pq_sfcNodeSet.empty())
	{
		virNode = getSFC.pq_sfcNodeSet.top();
		vnfNode.sfc_id = virNode.sfc_id;
		vnfNode.totalComRes = virNode.reqComRes + virNode.instantRes;
		vnfNode.reqComRes = virNode.reqComRes;
		vnfNode.instantRes = virNode.instantRes;
		vnfNode.virNetFuc_id = virNode.virNetFuc_id;
		vnfNode.funC = virNode.funC;
		getVnfNode.push(vnfNode);
		getSFC.pq_sfcNodeSet.pop();
	}
	int y = 0;
	getVnfNode1 = getVnfNode;
	getVnfNode2 = getVnfNode;

	while (!getVnfNode.empty())
	{
		vnfNode = getVnfNode.top();
		ithVnfNode.insert(make_pair(y, vnfNode));
		y++;
		getVnfNode.pop();
	}
	getSFC = getASFC;
	double ABmax = -1;
	double ABmin = 10000;
	for (int i = 0; i < (int)getSFC.pq_sfcNodeSet.size(); i++)//该虚节点周围虚链路的带宽均值
	{
		double virnodeB = 0;
		int linkn = 0;
		for (int j = 0; j < (int)getSFC.pq_sfcNodeSet.size(); j++)
		{
			map<nodePair, virLinkSeral>::iterator findVirLink = allVirLink.find(make_pair(i, j));
			if (findVirLink != allVirLink.end())
			{
				linkn++;
				virnodeB += findVirLink->second.virWeight;
			}
			map<nodePair, virLinkSeral>::iterator findVirLink2 = allVirLink.find(make_pair(j, i));
			if (findVirLink2 != allVirLink.end())
			{
				linkn++;
				virnodeB += findVirLink2->second.virWeight;
			}
		}
		double inf = (double)virnodeB / linkn;
		if (inf < ABmin)
			ABmin = inf;
		if (inf > ABmax)
			ABmax = inf;
		VirNodeAB.insert(make_pair(i, inf));
	}

	int Nodemax = -1;
	int Nodemin = -1;

	Nodemax = ithVnfNode.find(0)->second.totalComRes;
	Nodemin = ithVnfNode.find(ithVirNode.size() - 1)->second.totalComRes;
	int Nodecha = Nodemax - Nodemin;
	double ABcha = ABmax - ABmin;


	getVirNode1 = getSFC.pq_sfcNodeSet;

	while (!getVnfNode1.empty())
	{
		double Nodeuniform = 0;
		double ABuniform = 0;
		vnfNode = getVnfNode1.top();
		ABuniform = (double)(VirNodeAB.find(vnfNode.virNetFuc_id)->second - ABmin) / ABcha;
		Nodeuniform = (double)(vnfNode.totalComRes - Nodemin) / Nodecha;

		vnfNode1.vnfNodeIm = c * Nodeuniform + b * ABuniform;//陈琪第五章的0.8和0.2
		vnfNode1.totalComRes = vnfNode.totalComRes;
		vnfNode1.virNetFuc_id = vnfNode.virNetFuc_id;
		vnfNode1.sfc_id = vnfNode.sfc_id;
		vnfNode1.funC = vnfNode.funC;

		getVnfNodeIm.push(vnfNode1);
		getVnfNode1.pop();
	}

	int k;
	while (!getVnfNodeIm.empty())
	{
		vnfNode1 = getVnfNodeIm.top();
		k = vnfNode1.virNetFuc_id;
		vnfNodeVector.push_back(k);
		getVnfNodeIm.pop();
	}

	for (int i = 0; i < (int)vnfNodeVector.size(); i++)
	{
		vnfNodeOrder.insert(make_pair(i, vnfNodeVector.at(i)));
	}
	///////////////////////////////////////////////虚拟节点映射顺序

	//////////////////////////////////////////////////////////////////////////
	sfcSet aSFC = getASFC;
	for (int i = 0; i < (int)getVnfNode2.size(); i++)
	{
		getVnfNode3 = getVnfNode2;
		while (!getVnfNode3.empty())
		{
			vnfNode = getVnfNode3.top();
			if (i == vnfNode.virNetFuc_id)
			{
				vnfNodeRes.insert(make_pair(i, vnfNode));//虚拟节点按序号排序
				unfunCphyVnf1.push_back(i);
				break;
			}
			getVnfNode3.pop();
		}
	}
	for (int i = 0; i < vnfNodeRes.size(); i++)
	{
		int fun_ction;
		fun_ction = vnfNodeRes.find(i)->second.funC;
		vector<int> candiphyNode;
		if (candiphyNode.size() != 0)
		{
			candiphyNode.clear();
		}
		for (int j = 0; j < MaxNodeNum; j++)
		{
			if (phyVnfPair[j] == fun_ction)
			{
				candiphyNode.push_back(j);
			}
		}
		if (candiphyNode.size() != 0)
		{
			funCphyVnf[i] = candiphyNode;
		}
	}

	for (map<int, vector<int>>::iterator it = funCphyVnf.begin(); it != funCphyVnf.end(); it++)
	{
		int Wjr;
		Wjr = it->first;
		int j = 0;
		for (int i = 0; i < unfunCphyVnf1.size(); i++)
		{
			if (unfunCphyVnf1[i] == Wjr)
			{
				swap(unfunCphyVnf1[i], unfunCphyVnf1[unfunCphyVnf1.size() - 1]);
				unfunCphyVnf1.pop_back();
				i--;
			}
		}
	}
	unfunCphyVnf = unfunCphyVnf1;
	for (int i = 0; i < vnfNodeVector.size(); i++)
	{
		int fiona = vnfNodeVector[i];
		for (int j = 0; j < unfunCphyVnf.size(); j++)
		{
			if (fiona == unfunCphyVnf[j])
			{
				vnfNodeVector1.push_back(fiona);//业务中存在未实例化的vnf，这些vnf的映射顺序
			}
		}
	}
	for (int i = 0; i < vnfNodeVector.size(); i++)
	{
		int kenny = vnfNodeVector[i];
		for (map<int, vector<int>>::iterator it = funCphyVnf.begin(); it != funCphyVnf.end(); it++)
		{
			if (kenny == it->first)
			{
				vnfNodeVector2.push_back(kenny);
			}
		}
	}


	for (int i = 0; i < (int)pq_getPhyNodeSer.size(); i++)
	{
		pq_phyNodeVirSer = pq_getPhyNodeSer;
		while (!pq_phyNodeVirSer.empty())
		{
			phyNode = pq_phyNodeVirSer.top();
			if (i == phyNode.phyNode_id)
			{
				phyNodeRes.insert(make_pair(i, phyNode));//物理节点按序号排序
				break;
			}
			pq_phyNodeVirSer.pop();
		}
	}
	int ithNode = 0;
	pq_phyNodeVirSer = pq_getPhyNodeSer;
	while (!pq_phyNodeVirSer.empty())//物理节点计算资源从大到小，且带标号
	{
		phyNode = pq_phyNodeVirSer.top();
		phyNodeRes1.insert(make_pair(ithNode, phyNode));
		ithNode++;
		pq_phyNodeVirSer.pop();
	}

	for (int j = 0; j < MaxNodeNum; j++)
	{
		int slots = 0;
		double Weight = 0;
		for (int l = 0; l < MaxNodeNum; l++)
		{
			if (l != j)
			{
				map<nodePair, nodePairElement>::iterator ithNodePair = allNodePair.find(make_pair(j, l));
				if (ithNodePair != allNodePair.end())
				{
					for (int i = 0; i < MaxNodeNum; i++)
					{
						priPath[i] = -1;
					}
					int node = 0;
					phyNode1.source = j;
					phyNode1.primaryPath = ithNodePair->second.primaryPath;
					if (ithNodePair->second.primaryPath.size() == 2)
					{
						Weight += ithNodePair->second.totalWeigth;
					}
					for (int i = 0; i < (int)phyNode1.primaryPath.size(); i++)
					{
						priPath[node] = phyNode1.primaryPath.at(i);
						node++;
					}
					if (ithNodePair->second.primaryPath.size() == 2)
					{
						for (int core = 0; core < CoreNumber; core++)
						{
							for (int eachWave = 0; eachWave < WaveNumber; eachWave++)
							{
								eachLinkWavaLength[eachWave] = Free;
							}
							for (int i = 0; i < (int)phyNode1.primaryPath.size() - 1; i++)
							{
								int node1 = priPath[i];
								int nodeConWithNode1 = priPath[i + 1];
								if (nodeConWithNode1 == -1)
								{
									break;
								}
								map< Link, map<int, vector<eachWaveLength>> >::iterator resPriPathNeighborNode = m_linkstate.linkState2.find(make_pair(node1, nodeConWithNode1));
								findIthCoreAndWaveLength = resPriPathNeighborNode->second;
								map<int, vector<eachWaveLength>>::iterator findWaveLengthInIthCore = findIthCoreAndWaveLength.find(core);
								int eachWaveNum = 0;
								while (eachWaveNum < WaveNumber)
								{
									if (findWaveLengthInIthCore->second[eachWaveNum].waveLengthState == BusyPrimary)
									{
										eachLinkWavaLength[eachWaveNum] = BusyPrimary;
									}
									eachWaveNum++;
								}
							}
							for (int eachWaveNum = 0; eachWaveNum < WaveNumber; eachWaveNum++)
							{
								if (eachLinkWavaLength[eachWaveNum] == 0)
									slots++;
							}
						}
					}
				}
			}
		}
		phyNode1.totalWeigth = Weight;
		phyNode1.slots = slots;
		phyNode2.slots = slots;
		fixNodePair1.push(phyNode1);
		fixNodePair2.push(phyNode2);
	}

	int num = 0;
	while (fixNodePair1.size() != 0)
	{
		phyNode1 = fixNodePair1.top();
		phyNode2 = fixNodePair2.top();
		eachNodePair1.insert(make_pair(num, phyNode1));//带标号，物理节点按相邻物理链路距离，从小到达排序
		eachNodePair2.insert(make_pair(num, phyNode2));//带标号，物理节点按相邻物理链路上的可用频隙数量，从小到大排序
		num++;
		fixNodePair1.pop();
		fixNodePair2.pop();
	}

	int phyNodemax = -1;
	int phyNodemin = -1;
	double phyWeightmax = -1;
	double phyWeightmin = -1;
	int phySlotsmax = -1;
	int phySlotsmin = -1;

	phyNodemax = phyNodeRes1.find(0)->second.comRes;
	phyNodemin = phyNodeRes1.find(phyNodeRes1.size() - 1)->second.comRes;

	phyWeightmin = eachNodePair1.find(0)->second.totalWeigth;
	phyWeightmax = eachNodePair1.find(eachNodePair1.size() - 1)->second.totalWeigth;

	phySlotsmin = eachNodePair2.find(0)->second.slots;
	phySlotsmax = eachNodePair2.find(eachNodePair2.size() - 1)->second.slots;

	int phyNodecha = phyNodemax - phyNodemin;
	double phyWeightcha = phyWeightmax - phyWeightmin;
	double phySlotscha = phySlotsmax - phySlotsmin + 0.01;
	for (int i = 0; i < (int)eachNodePair1.size(); i++)
	{
		double phynodeIm = 0;
		double phyNodeuniform = 0;
		double phyWeightuniform = 0;
		double phySlotsuniform = 0;
		phyNode1 = eachNodePair1.find(i)->second;
		phyWeightuniform = (double)(phyWeightmax - phyNode1.totalWeigth) / phyWeightcha;
		phyNodeuniform = (double)(phyNodeRes.find(phyNode1.source)->second.comRes - phyNodemin) / phyNodecha;
		phySlotsuniform = (double)(phyNode1.slots - phySlotsmin) / phySlotscha;
		phynodeIm = 2 * phyNodeuniform + 2 * phySlotsuniform + 1 * phyWeightuniform;//计算资源越多，频谱隙越多，距离周围节点距离越短

		phyNodeim.phynodeIm = phynodeIm;
		phyNodeim.phyNode_id = phyNode1.source;
		phyNodeim.comRes = phyNodeRes.find(phyNode1.source)->second.comRes;
		fixNodePair.push(phyNodeim);//所有节点按重要度从大到小
	}///////////////////////////////////////////////////////////////////////////物理节点的顺序；
	eachNodePair = fixNodePair;
	unMapPhyNode = phyNodeRes;//物理节点按序号排列，带标号
	int NodePair[3] = { -1,-1,-1 };
	bool fond = false;


	//开始映射虚拟vnf
	if (funCphyVnf.size() == 0)
	{
		for (int i = 0; i < (int)vnfNodeVector.size(); i++)
		{
			bool pp = false;
			status = false;
			if (i == vnfNodeOrder.size() - 1)
			{
				pp = true;
			}
			int N;
			map<int, vnfNodeSeral>::iterator vir = vnfNodeRes.find(vnfNodeVector.at(i));
			while (mapPhyNode.size() == 0)
			{
				if (eachNodePair.size() != 0)
				{
					N = eachNodePair.top().phyNode_id;
					if (phyVnfPair[N] == -1)
					{
						map<int, phyNodeSeral>::iterator  phy = phyNodeRes.find(N);
						if (phy->second.comRes >= vir->second.totalComRes)
						{
							NodePair[vir->second.virNetFuc_id] = N;
							sfcMapToPhy.ithSFC = getASFC.SFC_id;
							sfcMapToPhy.phyNode = phy->second;
							sfcMapToPhy.vnfNode = vir->second;//这里有点问题，我暂时先把sfcMapToPhy的virNodeSeral属性换成vnfNodeSeral
							aSFCMapSet.push_back(sfcMapToPhy);
							mapPhyNode.insert(make_pair(ithMapNum, sfcMapToPhy));
							ithMapNum++;
							map<int, phyNodeSeral> deleMapPhyNode;
							int indexNode = 0;
							for (int i = 0; i < (int)unMapPhyNode.size(); i++)
							{
								map<int, phyNodeSeral>::iterator phy_nodeId = unMapPhyNode.find(i);
								if (phy_nodeId->second.phyNode_id != sfcMapToPhy.phyNode.phyNode_id)
								{
									deleMapPhyNode.insert(make_pair(indexNode, phy_nodeId->second));
									indexNode++;
								}
								if (i == (int)unMapPhyNode.size() - 1)
								{
									unMapPhyNode.clear();
									unMapPhyNode = deleMapPhyNode;
									deleMapPhyNode.clear();
								}
							}
							status = true;
							break;
						}
					}
					eachNodePair.pop();
				}
				else
				{
					status = false;
					break;
				}
			}
			if (i > 0)
			{
				eachNodePair = fixNodePair;
				while (!eachNodePair.empty())
				{
					fond = false;
					lp = false;
					bool dd = false;
					bool oo = false;
					if (eachNodePair.size() == 1)
					{
						dd = true;
						if (dd & pp)
							oo = true;

					}
					int v = 0;
					N = eachNodePair.top().phyNode_id;
					if (phyVnfPair[N] == -1)
					{
						for (int j = 0; j < (int)unMapPhyNode.size(); j++)
						{
							map<int, phyNodeSeral>::iterator Mapnode = unMapPhyNode.find(j);
							if (Mapnode->second.phyNode_id == N)
							{
								fond = true;
								break;
							}

						}
					}
					if (fond == true)
					{
						map<int, phyNodeSeral>::iterator  phy = phyNodeRes.find(N);
						if (phy->second.comRes >= vir->second.totalComRes)
						{
							if (virLinkMapOrder1.size() != 0)
							{
								virLinkMapOrder1.clear();
							}
							for (int i = 0; i < (int)aSFC.pq_sfcNodeSet.size(); i++)
							{
								for (int j = 0; j < (int)aSFC.pq_sfcNodeSet.size(); j++)
								{
									if (i != j)
									{
										map<nodePair, virLinkSeral>::iterator findVirLink = allVirLink.find(make_pair(i, j));
										if (findVirLink != allVirLink.end())
										{
											for (int l = 0; l < mapPhyNode.size(); l++)
											{
												if (i == mapPhyNode.find(l)->second.vnfNode.virNetFuc_id && j == vir->second.virNetFuc_id)
												{

													virLinkSeral virLink = findVirLink->second;
													virLinkMapOrder1.insert(make_pair(v, virLink));
													v++;
												}
											}
										}
									}
								}
							}
							if (virLinkMapOrder1.size() == 0)
								lp = true;
							for (int m = 0; m < virLinkMapOrder1.size(); m++)
							{
								int virs;
								int vird;
								int ps;
								int pd;
								int GBPs;
								int slots;
								vector<int> primaryPath;
								virLinkSeral virLink = virLinkMapOrder1.find(m)->second;
								GBPs = virLink.virWeight;
								slots = Slots(GBPs);
								virs = virLink.source;
								vird = virLink.destination;
								ps = NodePair[virs];
								pd = eachNodePair.top().phyNode_id;

								primaryPath = allNodePair.find(make_pair(ps, pd))->second.primaryPath;
								lp = FPslots(primaryPath, slots, aSFC, m, virLinkMapOrder1.size() - 1, oo);
								if (lp == false)
								{
									break;
								}
							}
						}

						if (lp == true)
						{
							NodePair[vir->second.virNetFuc_id] = N;
							sfcMapToPhy.ithSFC = getASFC.SFC_id;
							sfcMapToPhy.vnfNode = vir->second;
							sfcMapToPhy.phyNode = phy->second;
							aSFCMapSet.push_back(sfcMapToPhy);
							mapPhyNode.insert(make_pair(ithMapNum, sfcMapToPhy));
							ithMapNum++;

							map<int, phyNodeSeral> deleMapPhyNode;
							int indexNode = 0;
							for (int i = 0; i < (int)unMapPhyNode.size(); i++)
							{
								map<int, phyNodeSeral>::iterator phy_nodeId = unMapPhyNode.find(i);
								if (phy_nodeId->second.phyNode_id != N)
								{
									deleMapPhyNode.insert(make_pair(indexNode, phy_nodeId->second));
									indexNode++;
								}
								if (i == (int)unMapPhyNode.size() - 1)
								{
									unMapPhyNode.clear();
									unMapPhyNode = deleMapPhyNode;
									deleMapPhyNode.clear();
								}
							}
							status = true;
						}//该节点映射过程中，产生的所有相连链路均可以映射，则虚实节点对配对
						eachNodePair.pop();
						if (status == true)
						{
							break;
						}
					}
					if (fond == false)
					{
						eachNodePair.pop();
					}
				}

			}
			if (status == false)
			{
				break;
			}

		}
	}
	else
	{
		status = false;
		for (int i = 0; i < (int)funCphyVnf.size(); i++)
		{
			bool pp = false;
			status = false;
			if (i == vnfNodeVector.size() - 1)
			{
				pp = true;
			}
			vector<int> M1;
			int M2;
			vector<int> M3;
			int M4;
			int M5;
			int Bill = vnfNodeVector2[i];
			map<int, vnfNodeSeral>::iterator vir = vnfNodeRes.find(Bill);
			while (mapPhyNode.size() == 0)
			{
				M1 = funCphyVnf.find(Bill)->second;
				for (int i = 0; i < M1.size(); i++)
				{
					int M11;
					M11 = M1[i];
					map<int, phyNodeSeral>::iterator phy1 = phyNodeRes.find(M11);
					if (phy1->second.comRes >= vir->second.reqComRes)
					{
						NodePair[vir->second.virNetFuc_id] = M11;
						sfcMapToPhy.ithSFC = getASFC.SFC_id;
						sfcMapToPhy.phyNode = phy1->second;
						sfcMapToPhy.vnfNode = vir->second;
						aSFCMapSet.push_back(sfcMapToPhy);
						mapPhyNode.insert(make_pair(ithMapNum, sfcMapToPhy));
						ithMapNum++;
						map<int, phyNodeSeral>deleMapPhyNode;
						int indexNode = 0;
						for (int i = 0; i < (int)unMapPhyNode.size(); i++)
						{
							map<int, phyNodeSeral>::iterator phy_nodeId = unMapPhyNode.find(i);
							if (phy_nodeId->second.phyNode_id != sfcMapToPhy.phyNode.phyNode_id)
							{
								deleMapPhyNode.insert(make_pair(indexNode, phy_nodeId->second));
								indexNode++;
							}
							if (i == (int)unMapPhyNode.size() - 1)
							{
								unMapPhyNode.clear();
								unMapPhyNode = deleMapPhyNode;
								deleMapPhyNode.clear();
							}
						}
						status = true;
						break;
					}
				}

				if (status == false)//已经实例化过的vnf中的，第一个点，没有lp这一步，有备用的点就用，没有就直接status=false;
				{
					eachNodePair = fixNodePair;
					while (!eachNodePair.empty())
					{
						M2 = eachNodePair.top().phyNode_id;
						if (phyVnfPair[M2] == -1)
						{
							map<int, phyNodeSeral>::iterator phy2 = phyNodeRes.find(M2);
							if (phy2->second.comRes >= vir->second.reqComRes)
							{
								NodePair[vir->second.virNetFuc_id] = M2;
								sfcMapToPhy.ithSFC = getASFC.SFC_id;
								sfcMapToPhy.phyNode = phy2->second;
								sfcMapToPhy.vnfNode = vir->second;
								aSFCMapSet.push_back(sfcMapToPhy);
								mapPhyNode.insert(make_pair(ithMapNum, sfcMapToPhy));
								ithMapNum++;

								map<int, phyNodeSeral>deleMapPhyNode;
								int indexNode = 0;
								for (int i = 0; i < (int)unMapPhyNode.size(); i++)
								{
									map<int, phyNodeSeral>::iterator phy_nodeId = unMapPhyNode.find(i);
									if (phy_nodeId->second.phyNode_id != sfcMapToPhy.phyNode.phyNode_id)
									{
										deleMapPhyNode.insert(make_pair(indexNode, phy_nodeId->second));
										indexNode++;
									}
									if (i == (int)unMapPhyNode.size() - 1)
									{
										unMapPhyNode.clear();
										unMapPhyNode = deleMapPhyNode;
										deleMapPhyNode.clear();
									}
								}
								status = true;
								break;
							}
						}
						eachNodePair.pop();
					}

				}
				if (status == false)
				{
					break;
				}
			}
			if (i > 0)
			{
				M3 = funCphyVnf.find(Bill)->second;
				for (int i = 0; i < M3.size(); i++)
				{
					int M33;
					M33 = M3[i];
					int v = 0;
					bool oo = false;
					map<int, phyNodeSeral>::iterator phy3 = phyNodeRes.find(M33);
					if (phy3->second.comRes >= vir->second.reqComRes)
					{
						bool lp = false;
						if (virLinkMapOrder1.size() != 0)
						{
							virLinkMapOrder1.clear();
						}
						for (int i = 0; i < (int)aSFC.pq_sfcNodeSet.size(); i++)
						{
							for (int j = 0; j < (int)aSFC.pq_sfcNodeSet.size(); j++)
							{
								if (i != j)
								{
									map<nodePair, virLinkSeral>::iterator findVirLink1 = allVirLink.find(make_pair(i, j));
									if (findVirLink1 != allVirLink.end())
									{
										for (int l = 0; l < mapPhyNode.size(); l++)
										{
											if (i == mapPhyNode.find(l)->second.vnfNode.virNetFuc_id && j == vir->second.virNetFuc_id)
											{

												virLinkSeral virLink1 = findVirLink1->second;
												virLinkMapOrder1.insert(make_pair(v, virLink1));
												v++;
											}
										}
									}
								}
							}
						}
						if (virLinkMapOrder1.size() == 0)
							lp = true;
						for (int m = 0; m < virLinkMapOrder1.size(); m++)
						{
							int virs;
							int vird;
							int ps;
							int pd;
							int GBPs;
							int slots;
							vector<int> primaryPath;
							virLinkSeral virLink = virLinkMapOrder1.find(m)->second;
							GBPs = virLink.virWeight;
							slots = Slots(GBPs);
							virs = virLink.source;
							vird = virLink.destination;
							ps = NodePair[virs];
							pd = M33;

							primaryPath = allNodePair.find(make_pair(ps, pd))->second.primaryPath;
							lp = FPslots(primaryPath, slots, aSFC, m, virLinkMapOrder1.size() - 1, oo);
							if (lp == false)
							{
								break;
							}
						}
						if (lp == true)
						{
							NodePair[vir->second.virNetFuc_id] = M33;
							sfcMapToPhy.ithSFC = getASFC.SFC_id;
							sfcMapToPhy.vnfNode = vir->second;
							sfcMapToPhy.phyNode = phy3->second;
							aSFCMapSet.push_back(sfcMapToPhy);
							mapPhyNode.insert(make_pair(ithMapNum, sfcMapToPhy));
							ithMapNum++;

							map<int, phyNodeSeral> deleMapPhyNode;
							int indexNode = 0;
							for (int i = 0; i < (int)unMapPhyNode.size(); i++)
							{
								map<int, phyNodeSeral>::iterator phy_nodeId = unMapPhyNode.find(i);
								if (phy_nodeId->second.phyNode_id != M33)
								{
									deleMapPhyNode.insert(make_pair(indexNode, phy_nodeId->second));
									indexNode++;
								}
								if (i == (int)unMapPhyNode.size() - 1)
								{
									unMapPhyNode.clear();
									unMapPhyNode = deleMapPhyNode;
									deleMapPhyNode.clear();
								}
							}
							status = true;
							break;
						}
					}
				}
				if (status == false)
				{
					eachNodePair = fixNodePair;
					while (!eachNodePair.empty())
					{
						fond = false;
						lp = false;
						bool dd = false;
						bool oo = false;
						if (eachNodePair.size() == 1)
						{
							dd = true;
							if (dd & pp)
								oo = true;
						}
						int v = 0;
						M4 = eachNodePair.top().phyNode_id;
						if (phyVnfPair[M4] == -1)
						{
							for (int j = 0; j < (int)unMapPhyNode.size(); j++)
							{
								map<int, phyNodeSeral>::iterator Mapnode = unMapPhyNode.find(j);
								if (Mapnode->second.phyNode_id == M4)
								{
									fond = true;
									break;
								}
							}
						}
						if (fond == true)
						{
							map<int, phyNodeSeral>::iterator phy = phyNodeRes.find(M4);
							if (phy->second.comRes >= vir->second.reqComRes)
							{
								if (virLinkMapOrder1.size() != 0)
								{
									virLinkMapOrder1.clear();
								}
								for (int i = 0; i < (int)aSFC.pq_sfcNodeSet.size(); i++)
								{
									for (int j = 0; j < (int)aSFC.pq_sfcNodeSet.size(); j++)
									{
										if (i != j)
										{
											map<nodePair, virLinkSeral>::iterator findVirLink = allVirLink.find(make_pair(i, j));
											if (findVirLink != allVirLink.end())
											{
												for (int l = 0; l < mapPhyNode.size(); l++)
												{
													if (i == mapPhyNode.find(l)->second.vnfNode.virNetFuc_id && j == vir->second.virNetFuc_id)
													{

														virLinkSeral virLink = findVirLink->second;
														virLinkMapOrder1.insert(make_pair(v, virLink));
														v++;
													}
												}
											}
										}
									}
								}
								if (virLinkMapOrder1.size() == 0)
									lp = true;
								for (int m = 0; m < virLinkMapOrder1.size(); m++)
								{
									int virs;
									int vird;
									int ps;
									int pd;
									int GBPs;
									int slots;
									vector<int>primaryPath;
									virLinkSeral virLink = virLinkMapOrder1.find(m)->second;
									GBPs = virLink.virWeight;
									slots = Slots(GBPs);
									virs = virLink.source;
									vird = virLink.destination;
									ps = NodePair[virs];
									pd = M4;
									primaryPath = allNodePair.find(make_pair(ps, pd))->second.primaryPath;
									lp = FPslots(primaryPath, slots, aSFC, m, virLinkMapOrder1.size() - 1, oo);
									if (lp == false)
									{
										break;
									}
								}
							}
							if (lp == true)
							{
								NodePair[vir->second.virNetFuc_id] = M4;
								sfcMapToPhy.ithSFC = getASFC.SFC_id;
								sfcMapToPhy.vnfNode = vir->second;
								sfcMapToPhy.phyNode = phy->second;
								aSFCMapSet.push_back(sfcMapToPhy);
								mapPhyNode.insert(make_pair(ithMapNum, sfcMapToPhy));
								ithMapNum++;

								map<int, phyNodeSeral> deleMapPhyNode;
								int indexNode = 0;
								for (int i = 0; i < (int)unMapPhyNode.size(); i++)
								{
									map<int, phyNodeSeral>::iterator phy_nodeId = unMapPhyNode.find(i);
									if (phy_nodeId->second.phyNode_id != M4)
									{
										deleMapPhyNode.insert(make_pair(indexNode, phy_nodeId->second));
										indexNode++;
									}
									if (i == (int)unMapPhyNode.size() - 1)
									{
										unMapPhyNode.clear();
										unMapPhyNode = deleMapPhyNode;
										deleMapPhyNode.clear();
									}
								}
								status = true;
							}
							eachNodePair.pop();
							if (status == true)
							{
								break;
							}
						}
						if (fond == false)
						{
							eachNodePair.pop();
						}
					}
				}
			}
			if (status == false)
			{
				break;
			}
		}
		if (status == true)
		{
			for (int i = 0; i < (int)vnfNodeVector1.size(); i++)
			{
				bool pp = false;
				status = false;
				if (i == vnfNodeVector.size() - funCphyVnf.size() - 1)
				{
					pp = true;
				}
				int N1;
				map<int, vnfNodeSeral>::iterator vir = vnfNodeRes.find(vnfNodeVector1.at(i));
				eachNodePair = fixNodePair;
				while (!eachNodePair.empty())
				{
					fond = false;
					lp = false;
					bool dd = false;
					bool oo = false;
					if (eachNodePair.size() == 1)
					{
						dd = true;
						if (dd & pp)
							oo = true;
					}
					int v = 0;
					N1 = eachNodePair.top().phyNode_id;
					if (phyVnfPair[N1] == -1)
					{
						for (int j = 0; j < (int)unMapPhyNode.size(); j++)
						{
							map<int, phyNodeSeral>::iterator Mapnode = unMapPhyNode.find(j);
							if (Mapnode->second.phyNode_id == N1)
							{
								fond = true;
								break;
							}
						}
					}
					if (fond == true)
					{
						map<int, phyNodeSeral>::iterator phy = phyNodeRes.find(N1);
						if (phy->second.comRes >= vir->second.totalComRes)
						{
							if (virLinkMapOrder1.size() != 0)
							{
								virLinkMapOrder1.clear();
							}
							for (int i = 0; i < (int)aSFC.pq_sfcNodeSet.size(); i++)
							{
								for (int j = 0; j < (int)aSFC.pq_sfcNodeSet.size(); j++)
								{
									if (i != j)
									{
										map<nodePair, virLinkSeral>::iterator findVirLink = allVirLink.find(make_pair(i, j));
										if (findVirLink != allVirLink.end())
										{
											for (int l = 0; l < mapPhyNode.size(); l++)
											{
												if (i == mapPhyNode.find(l)->second.vnfNode.virNetFuc_id && j == vir->second.virNetFuc_id)
												{

													virLinkSeral virLink = findVirLink->second;
													virLinkMapOrder1.insert(make_pair(v, virLink));
													v++;
												}
											}
										}
									}
								}
							}
							if (virLinkMapOrder1.size() == 0)
								lp = true;
							for (int m = 0; m < virLinkMapOrder1.size(); m++)
							{
								int virs;
								int vird;
								int ps;
								int pd;
								int GBPs;
								int slots;
								vector<int> primaryPath;
								virLinkSeral virLink = virLinkMapOrder1.find(m)->second;
								GBPs = virLink.virWeight;
								slots = Slots(GBPs);
								virs = virLink.source;
								vird = virLink.destination;
								ps = NodePair[virs];
								pd = N1;

								primaryPath = allNodePair.find(make_pair(ps, pd))->second.primaryPath;
								lp = FPslots(primaryPath, slots, aSFC, m, virLinkMapOrder1.size() - 1, oo);
								if (lp == false)
								{
									break;
								}
							}
						}
						if (lp == true)
						{
							NodePair[vir->second.virNetFuc_id] = N1;
							sfcMapToPhy.ithSFC = getASFC.SFC_id;
							sfcMapToPhy.vnfNode = vir->second;
							sfcMapToPhy.phyNode = phy->second;
							aSFCMapSet.push_back(sfcMapToPhy);
							mapPhyNode.insert(make_pair(ithMapNum, sfcMapToPhy));
							ithMapNum++;

							map<int, phyNodeSeral> deleMapPhyNode;
							int indexNode = 0;
							for (int i = 0; i < (int)unMapPhyNode.size(); i++)
							{
								map<int, phyNodeSeral>::iterator phy_nodeId = unMapPhyNode.find(i);
								if (phy_nodeId->second.phyNode_id != N1)
								{
									deleMapPhyNode.insert(make_pair(indexNode, phy_nodeId->second));
									indexNode++;
								}
								if (i == (int)unMapPhyNode.size() - 1)
								{
									unMapPhyNode.clear();
									unMapPhyNode = deleMapPhyNode;
									deleMapPhyNode.clear();
								}
							}
							status = true;
						}
						eachNodePair.pop();
						if (status == true)
						{
							break;
						}
					}
					if (fond == false)
					{
						eachNodePair.pop();
					}
				}
				if (status == false)
				{
					break;
				}
			}
		}
	}

	//	//////////////////////////////////////////////////////////////////////////
	SFCMap.insert(make_pair(getASFC.SFC_id, aSFCMapSet));
	if (status == true)
	{
		vector<SFCMapPhy>aPhyNode1;
		map<int, vector<SFCMapPhy>>::iterator getPhyNode1 = SFCMap.find(getASFC.SFC_id);
		if (getPhyNode1 != SFCMap.end())
		{
			aPhyNode1 = getPhyNode1->second;
			for (int i = 0; i < (int)aPhyNode1.size(); i++)
			{
				phyVnfPair[aPhyNode1.at(i).phyNode.phyNode_id] = aPhyNode1.at(i).vnfNode.funC;
			}
		}
	}
	if (status == true)
	{
		/*VONMap.insert(make_pair(getAVON.VON_id,aVONMapSet));*/
		while (!pq_phyNodeVirSer.empty())
		{
			pq_phyNodeVirSer.top();
			pq_phyNodeVirSer.pop();
		}

		while (!pq_getPhyNodeSer.empty())
		{
			phyNodeSeral node = pq_getPhyNodeSer.top();
			vector<SFCMapPhy> aPhyNode;
			map<int, vector<SFCMapPhy> >::iterator getPhyNode = SFCMap.find(getASFC.SFC_id);
			if (getPhyNode != SFCMap.end())
			{
				aPhyNode = getPhyNode->second;
				if (funCphyVnf.size() == 0)
				{
					for (int i = 0; i < (int)aPhyNode.size(); i++)
					{
						if (node.phyNode_id == aPhyNode.at(i).phyNode.phyNode_id)
						{
							node.comRes = node.comRes - aPhyNode.at(i).vnfNode.totalComRes;
							InstantCost += aPhyNode.at(i).vnfNode.instantRes;
							break;
						}
					}
				}
				else
				{
					for (map<int, vector<int>>::iterator it = funCphyVnf.begin(); it != funCphyVnf.end(); it++)
					{
						vector<int> ocrPhy = it->second;
						for (int i = 0; i < (int)ocrPhy.size(); i++)
						{
							int p = ocrPhy[i];
							for (int j = 0; j < (int)aPhyNode.size(); j++)
							{
								if (p == aPhyNode.at(j).phyNode.phyNode_id)
								{
									if (node.phyNode_id == aPhyNode.at(j).phyNode.phyNode_id)
									{
										node.comRes = node.comRes - aPhyNode.at(j).vnfNode.reqComRes;
									}
								}
							}
						}
					}
					for (int i = 0; i < (int)unfunCphyVnf.size(); i++)
					{
						for (int j = 0; j < (int)aPhyNode.size(); j++)
						{
							if (unfunCphyVnf[i] == aPhyNode.at(j).vnfNode.virNetFuc_id)
							{
								if (node.phyNode_id == aPhyNode.at(j).phyNode.phyNode_id)
								{
									node.comRes = node.comRes - aPhyNode.at(j).vnfNode.totalComRes;
									InstantCost += aPhyNode.at(j).vnfNode.instantRes;
								}
							}
						}
					}
				}
			}
			pq_phyNodeVirSer.push(node);////////////////////更新后的物理节点的计算资源
			pq_getPhyNodeSer.pop();
		}
		pq_getPhyNodeSer = pq_phyNodeVirSer;
	}
	if (SFCMap.find(getASFC.SFC_id)->second.size() != virMaxNode)
	{
		cout << getASFC.SFC_id << " is failure" << endl;
		blockTime1++;
		blockTime2 = blockTime1;
	}
	else
	{
		pq_getPhyNodeSer = pq_phyNodeVirSer;
	}
	////////////////////////////////////////////////////////////////////////////
}

void Database::aSFCMapping(int ithTimeVON, sfcSet getaSFC, vector<SFCMapPhy> aSFCSet)
{
	virLinkSeral virLink1;
	virLinkSeral virLink;
	map<nodePair, virLinkSeral> allVirLink;
	sfcSet getVON;
	priority_queue<virLinkSeral>  getVirLink1;
	priority_queue<virLinkSeral>  getVirLink2;
	map<int, int> virPhyNodes;
	int tra = 0;
	bool status = true;
	NodeSer nodeSer, aResult;
	vector<int> primaryPath, protectPath;

	if ((int)virPhyNodes.size() != 0)
	{
		virPhyNodes.clear();
	}

	if ((int)getaSFC.pq_sfcNodeSet.size() != (int)aSFCSet.size())
	{
		status = false;
	}
	if (!allVirLink.empty())
	{
		allVirLink.clear();
	}
	while (!getVirLink1.empty())
	{
		getVirLink1.top();
		getVirLink1.pop();
	}

	for (int i = 0; i < (int)aSFCSet.size(); i++)
	{
		SFCMapPhy von = aSFCSet.at(i);
		int virNode = von.vnfNode.virNetFuc_id;
		int phyNode = von.phyNode.phyNode_id;
		virPhyNodes.insert(make_pair(virNode, phyNode));//虚实节点，配对信息
	}

	aResult.SFCRev = 0;
	aResult.nodeNum = 0;
	aResult.energy = 0;
	aResult.specNum = 0;
	aResult.transNum = 0;
	aResult.regenNum = 0;

	aResult.transNum400 = 0;
	aResult.transNum100 = 0;
	aResult.transNum40 = 0;

	aResult.regenNum400 = 0;
	aResult.regenNum100 = 0;
	aResult.regenNum40 = 0;

	aResult.energy400 = 0;
	aResult.energy100 = 0;
	aResult.energy40 = 0;

	aResult.specNum400 = 0;
	aResult.specNum100 = 0;
	aResult.specNum40 = 0;

	aResult.slots400 = 0;
	aResult.slots100 = 0;
	aResult.slots40 = 0;

	aResult.waveLengthNum = 0;
	aResult.waveLengthNum400 = 0;
	aResult.waveLengthNum100 = 0;
	aResult.waveLengthNum40 = 0;

	aResult.channelWidth400 = 0;
	aResult.channelWidth100 = 0;
	aResult.channelWidth40 = 0;

	aResult.totalWeigth = 0;
	aResult.hop = 0;
	while (!getaSFC.pq_sfcLinkSet.empty())//是双向的，因为映射当前节点与已映射节点的时候，链接的链路可能是回头的比如021顺序时，0-1和2-1，2-1就是回头的
	{
		int virSou;
		int virDes;
		//virLinkSeral virLink;
		virLink1 = getaSFC.pq_sfcLinkSet.top();
		virSou = virLink1.source;
		virDes = virLink1.destination;
		allVirLink.insert(make_pair(make_pair(virSou, virDes), virLink1));//双向的虚拟链路，前面带（源，宿）
		getaSFC.pq_sfcLinkSet.pop();
	}
	for (int i = 0; i < (int)getaSFC.pq_sfcNodeSet.size(); i++)
	{
		for (int j = i; j < (int)getaSFC.pq_sfcNodeSet.size(); j++)
		{
			if (i != j)
			{
				map<nodePair, virLinkSeral>::iterator findVirLink1 = allVirLink.find(make_pair(i, j));
				if (findVirLink1 != allVirLink.end())
				{
					getVirLink2.push(findVirLink1->second);//单向的虚拟链路
				}
			}
		}
	}
	while (!getVirLink2.empty())
	{
		if (status == false)
		{
			cout << getaSFC.SFC_id << ": This SFC mapping is failure" << endl;
			break;
		}
		virLinkSeral virLinkSet = getVirLink2.top();//取出单向链路里权重最大的
		int virSource = virLinkSet.source;
		int virDestination = virLinkSet.destination;
		int Gbps = virLinkSet.virWeight;

		map<int, int>::iterator sourNode = virPhyNodes.find(virSource);
		map<int, int>::iterator destinationNode = virPhyNodes.find(virDestination);
		if (sourNode != virPhyNodes.end() && destinationNode != virPhyNodes.end())
		{
			int phySource = sourNode->second;
			int phyDestination = destinationNode->second;

			map<nodePair, nodePairElement>::iterator findDisPath = allNodePair.find(make_pair(phySource, phyDestination));//单向链路源宿   对应的  物理节点源宿
			if (findDisPath != allNodePair.end())
			{
				primaryPath = findDisPath->second.primaryPath;

				NodeSer ithPriPathEnergy = energyFuction(primaryPath, Gbps);//工作路径，工作路径上的业务带宽

				aResult.hop += (double)ithPriPathEnergy.hop;
				aResult.energy += ithPriPathEnergy.energy;
				aResult.specNum += ithPriPathEnergy.specNum;
				aResult.transNum += ithPriPathEnergy.transNum;
				aResult.regenNum += ithPriPathEnergy.regenNum;

				aResult.transNum400 += ithPriPathEnergy.transNum400;
				aResult.transNum100 += ithPriPathEnergy.transNum100;
				aResult.transNum40 += ithPriPathEnergy.transNum40;

				aResult.regenNum400 += ithPriPathEnergy.regenNum400;
				aResult.regenNum100 += ithPriPathEnergy.regenNum100;
				aResult.regenNum40 += ithPriPathEnergy.regenNum40;

				aResult.energy400 += ithPriPathEnergy.energy400;
				aResult.energy100 += ithPriPathEnergy.energy100;
				aResult.energy40 += ithPriPathEnergy.energy40;

				aResult.specNum400 += ithPriPathEnergy.specNum400;
				aResult.specNum100 += ithPriPathEnergy.specNum100;
				aResult.specNum40 += ithPriPathEnergy.specNum40;

				aResult.slots400 += ithPriPathEnergy.slots400;
				aResult.slots100 += ithPriPathEnergy.slots100;
				aResult.slots40 += ithPriPathEnergy.slots40;

				aResult.waveLengthNum400 += ithPriPathEnergy.waveLengthNum400;
				aResult.waveLengthNum100 += ithPriPathEnergy.waveLengthNum100;
				aResult.waveLengthNum40 += ithPriPathEnergy.waveLengthNum40;
				aResult.waveLengthNum = aResult.waveLengthNum40;

				aResult.totalWeigth += ithPriPathEnergy.totalWeigth;

				totalSpectrum += ithPriPathEnergy.specNum;
				totalEnergy += ithPriPathEnergy.energy;
				TRANum += ithPriPathEnergy.transNum;
				REGNum += ithPriPathEnergy.regenNum;
				//tra = (ithPriPathEnergy.transNum + ithProPathEnergy.transNum);

				TRANum400_1500 += ithPriPathEnergy.transNum400;
				TRANum100_1700 += ithPriPathEnergy.transNum100;
				TRANum40_1800 += ithPriPathEnergy.transNum40;

				REGNum400_1500 += ithPriPathEnergy.regenNum400;
				REGNum100_1700 += ithPriPathEnergy.regenNum100;
				REGNum40_1800 += ithPriPathEnergy.regenNum40;

				//cout<<tra<<endl;
			}
		}
		getVirLink2.pop();
		if (getVirLink2.empty())
		{
			//ithVONResult.push_back(aResult);
			aResult.energy += InstantCost;
			aResult.SFCRev = getaSFC.SFCRev;
			int totalS = m_linkstate.linkState2.size() * CoreNumber * WaveNumber;
			aResult.spectrUtilization = (double)aResult.specNum / totalS;
			allResults.insert(make_pair(make_pair(ithTimeVON, getaSFC.SFC_id), aResult));
		}
	}
	if ((int)virPhyNodes.size() != 0)
	{
		virPhyNodes.clear();
	}
	//cout<<getaVON.VON_id<<" "<<TRANum<<endl;
}
void Database::resourceUtilizationRatio()//频谱碎片
{
	long int primaryLadaNum = 0, freeLadaNum = 0;
	int eachLinkWavaLength[WaveNumber];
	map<int, vector<eachWaveLength>> findIthCoreAndWaveLength;
	for (int node1 = 0; node1 < MaxNodeNum; node1++)
	{
		for (int nodeConWithNode1 = 0; nodeConWithNode1 < MaxNodeNum; nodeConWithNode1++)
		{
			for (int core = 0; core < CoreNumber; core++)
			{
				map< Link, map<int, vector<eachWaveLength>> >::iterator resPriPathNeighborNode = m_linkstate.linkState2.find(make_pair(node1, nodeConWithNode1));
				if (resPriPathNeighborNode != m_linkstate.linkState2.end())
				{
					findIthCoreAndWaveLength = resPriPathNeighborNode->second;
					map<int, vector<eachWaveLength>>::iterator resUtiRatio = findIthCoreAndWaveLength.find(core);
					int g = 0;
					if (resUtiRatio != findIthCoreAndWaveLength.end())//防止越界
					{
						if (resUtiRatio->second[0].waveLengthState == Free)
						{
							eachLinkWavaLength[0] = 0;
							g++;
						}
						for (int i = 0; i < (WaveNumber - 1); i++)
						{
							if (resUtiRatio->second[i].waveLengthState == BusyPrimary && resUtiRatio->second[i + 1].waveLengthState == Free)
							{
								eachLinkWavaLength[g] = i + 1;
								g++;
								/*primaryLadaNum++;*/
							}
						}
						for (int l = 0; l < g; l++)//eachLinkWaveLength[l]里放的是连续空闲频隙缝的首位，
						{
							freeLadaNum = 0;
							for (int j = eachLinkWavaLength[l]; j < WaveNumber; j++)
							{
								if (resUtiRatio->second[j].waveLengthState == Free)
									freeLadaNum++;
								else
									break;
							}//之后从第一位空闲频隙缝的统计，一共空闲了几个slot
							if (freeLadaNum < 16)//如果空闲频谱隙的个数小于16的话，我就认为它是频谱碎片
								primaryLadaNum += freeLadaNum;
						}
						for (int k = 0; k < WaveNumber; k++)
						{
							eachLinkWavaLength[k] = -1;

						}

					}

				}


			}

		}
	}
	int totalS = m_linkstate.linkState2.size() * CoreNumber * WaveNumber;//链路总频隙slot数7*80=560
	double SP;
	SP = (double)primaryLadaNum / totalS;//频谱利用率
	ofstream outfile;
	outfile.open("SP-NID-kun.xls", ios::app);//ios::out

	if (!outfile)
	{
		cerr << "存取数据的文件打不开" << endl;
	}
	else
	{
		outfile << SP << endl;
		outfile.close();
	}

}
void Database::statisticResult()
{
	double totalTime = 0;
	int sfcNum = 0;
	NodeSer aveResult;
	//double ave_hop;
	double toCost = 0;
	double cost400 = 0;
	double cost100 = 0;
	double cost40 = 0;

	double toSpec = 0;
	double spec400 = 0;
	double spec100 = 0;
	double spec40 = 0;

	double toTran = 0;
	double tran400 = 0;
	double tran100 = 0;
	double tran40 = 0;

	double toRege = 0;
	double rege400 = 0;
	double rege100 = 0;
	double rege40 = 0;
	double UseRatio = 0;
	double toRev = 0;
	double tocha = 0;
	double tospectrUtilization = 0;

	sfcNum = iniSFC;
	for (int figNode = 0; figNode < figNodeNum; figNode++)
	{
		totalTime = 0;
		aveResult.UseRatio = 0;
		aveResult.SFCRev = 0;

		aveResult.nodeNum = 0;
		aveResult.energy = 0;
		aveResult.specNum = 0;
		aveResult.transNum = 0;
		aveResult.regenNum = 0;

		aveResult.transNum400 = 0;
		aveResult.transNum100 = 0;
		aveResult.transNum40 = 0;

		aveResult.regenNum400 = 0;
		aveResult.regenNum100 = 0;
		aveResult.regenNum40 = 0;

		aveResult.energy400 = 0;
		aveResult.energy100 = 0;
		aveResult.energy40 = 0;

		aveResult.specNum400 = 0;
		aveResult.specNum100 = 0;
		aveResult.specNum40 = 0;

		aveResult.slots400 = 0;
		aveResult.slots100 = 0;
		aveResult.slots40 = 0;

		aveResult.waveLengthNum = 0;
		aveResult.waveLengthNum400 = 0;
		aveResult.waveLengthNum100 = 0;
		aveResult.waveLengthNum40 = 0;

		aveResult.channelWidth400 = 0;
		aveResult.channelWidth100 = 0;
		aveResult.channelWidth40 = 0;

		aveResult.totalWeigth = 0;
		aveResult.spectrUtilization = 0;
		//aveResult.hop = 0;

		for (int i = 0; i < runTime; i++)
		{
			for (int j = 0; j < sfcNum; j++)
			{

				//map<nodePair, double>::iterator timeResult = ithRunithVONTime.find(make_pair(i,j));
				//if (timeResult != ithRunithVONTime.end())
				//{
				//	totalTime += timeResult->second;
				//}

				map<nodePair, NodeSer>::iterator results = allResults.find(make_pair(i, j));
				if (results != allResults.end())
				{
					aveResult.spectrUtilization += results->second.spectrUtilization;
					aveResult.energy += results->second.energy;
					aveResult.energy400 += results->second.energy400;
					aveResult.energy100 += results->second.energy100;
					aveResult.energy40 += results->second.energy40;

					aveResult.specNum += results->second.specNum;
					aveResult.specNum400 += results->second.specNum400;
					aveResult.specNum100 += results->second.specNum100;
					aveResult.specNum40 += results->second.specNum40;

					aveResult.transNum += results->second.transNum;
					aveResult.transNum400 += results->second.transNum400;
					aveResult.transNum100 += results->second.transNum100;
					aveResult.transNum40 += results->second.transNum40;

					aveResult.regenNum += results->second.regenNum;
					aveResult.regenNum400 += results->second.regenNum400;
					aveResult.regenNum100 += results->second.regenNum100;
					aveResult.regenNum40 += results->second.regenNum40;
					aveResult.SFCRev += results->second.SFCRev;
					aveResult.hop += results->second.hop;
					totalTime++;
				}
			}
		}

		ofstream outfile;
		outfile.open("NID_B_ILP.xls", ios::app);
		//outfile.open("Result_6_3.txt", ios::out);
		if (!outfile)
		{
			cerr << "存取数据的文件打不开" << endl;
		}
		else
		{
			tospectrUtilization = (double)aveResult.spectrUtilization / runTime;
			toRev = aveResult.SFCRev / runTime;
			UseRatio = aveResult.UseRatio / runTime;
			toCost = aveResult.energy / runTime;
			tocha = toRev - toCost;
			cost400 = aveResult.energy400 / runTime;
			cost100 = aveResult.energy100 / runTime;
			cost40 = aveResult.energy40 / runTime;

			toSpec = (double)aveResult.specNum / runTime;
			spec400 = (double)aveResult.specNum400 / runTime;
			spec100 = (double)aveResult.specNum100 / runTime;
			spec40 = (double)aveResult.specNum40 / runTime;

			toTran = (double)aveResult.transNum / runTime;
			tran400 = (double)aveResult.transNum400 / runTime;
			tran100 = (double)aveResult.transNum100 / runTime;
			tran40 = (double)aveResult.transNum40 / runTime;

			toRege = (double)aveResult.regenNum / runTime;
			rege400 = (double)aveResult.regenNum400 / runTime;
			rege100 = (double)aveResult.regenNum100 / runTime;
			rege40 = (double)aveResult.regenNum40 / runTime;
			//ave_hop=(double)aveResult.hop/(runTime*totalTime);
			double bp = (double)blocka[figNode] / sfcNum;
			double bp_rt = bp / runTime;
			double time_rt = (double)run_Time[figNode] / runTime;

			outfile << sfcNum << "\t" << toRev << "\t" << toCost << "\t" << tocha << "\t" << toSpec << "\t" << toTran << "\t" << toRege << "\t\t"
				<< tospectrUtilization << "\t" << bp_rt << "\t" << time_rt << "\t" << endl;

		}
		sfcNum += SFCGap;
	}
}

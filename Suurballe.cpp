//本程序是从文件3420中读取数据,计算3420条序列的相似数以及它们的位置
//把数据存入pipe2-85中.

#include <math.h>
#include <stdlib.h>   //随机数产生头文件

#include "Suurballe.h"
#include <iostream>
#include <fstream>
using namespace std;

int Suurballe::frand(int ia,int ib)
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

void Suurballe::readTopology(int size, char* name)
{
	int computingResource=0;
	double weightValue;
	phyNodeSeral comRes;
	ifstream infile;
	infile.open(name, ios::in);
	if(!infile)
	{
		cerr << "open error!" << endl;
		exit(-2);
	}
	for(int i=0;i<size;i++)
	{
		int node=0;
		int link=0;
		for(int j=0;j<size;j++)
		{
			infile >> weightValue;
			if(weightValue!=INFINITE&&weightValue!=0)
			{
				link++;
				node+=weightValue;
			}
			weight[i][j] = weightValue;
			if (i != j)
			{
				updateWeight[i][j] = INFINITE;
			}
			else
			{
				updateWeight[i][j] = 0;
			}
		}
		NodeInf[i]=node/link;
		computingResource=frand(computingResourceLow,computingResourceHigh);
		comRes.comRes = computingResource;
		comRes.totalComRes = computingResource;
		comRes.phyNode_id = i;
		pq_phyNodeSer.push(comRes);
	}
	infile.close();
}

void Suurballe::Prim()
{
	int i,k,j,m,t;
	double w = 0, min = 0;
	for(i=0;i<MaxNodeNum-1;i++)
	{
		CT[i].fromvex=0;
		CT[i].endvex=i+1;
		CT[i].weight=weight[0][i+1];
	}
	for(k=1;k<MaxNodeNum;k++)
	{
		min = (INFINITE  + 10);
		m=k-1;
		for(j=k-1;j<MaxNodeNum-1;j++)
		{
			if(CT[j].weight< min)
			{   
				min=CT[j].weight;
				m=j;
			}
		}
		temp=CT[k-1];
		CT[k-1]=CT[m];
		CT[m]=temp;
		j=CT[k-1].endvex;
		for(i=k;i<MaxNodeNum-1;i++)
		{
			t=CT[i].endvex;
			w=weight[j][t];
			if(w<CT[i].weight)
			{
				CT[i].weight=w;
				CT[i].fromvex=j;
			}
		}
	}
	for(i=0;i<MaxNodeNum-1;i++)
	{
		int k = -1;
		int l = -1;
		k = CT[i].fromvex;
		l = CT[i].endvex;
		updateWeight[l][k] = updateWeight[k][l] = CT[i].weight;
		//updateWeight[k][l] = CT[i].weight;
	}
/*
	for(i=0;i<MaxNodeNum-1;i++)
	{
		cout<<CT[i].fromvex<<"\t";
	}
	cout<<endl;
	for(i=0;i<MaxNodeNum-1;i++)
	{
		cout<<CT[i].endvex<<"\t";
	}
	cout<<endl;
	for(i=0;i<MaxNodeNum-1;i++)
	{
		cout<<CT[i].weight<<"\t";
	}
	cout<<endl;*/
}


vector<int> Suurballe::Dijkstra(int source,int destination)
{
	int routeTable_hopNum = 0;
	int path[MaxNodeNum];
	int routeTable_nodeId[MaxNodeNum];
	double S[MaxNodeNum];
	double dist[MaxNodeNum];

	for (int i = 0; i<MaxNodeNum; i++)
	{
		path[i] = -1;
		routeTable_nodeId[i] = -1;
		S[i] = -1;
		dist[i] = -1;
	}

	nodePairElement pairElement;

	for (int i = 0; i < MaxNodeNum; i++)//initialize dist[n], S[source], path[n]
	{
		dist[i] = updateWeight[source][i];
		//cout<<updateWeight[source][i]<<endl;
		S[i] = 0;
		if(i != source && dist[i] < INFINITE)
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
	for (int num = 0; num<MaxNodeNum-1; num++)	 //当数组的顶点小于图的顶点数时循环
	{
		double min = INFINITE;
		int u = source;
		for (int i = 0; i<MaxNodeNum; i++) //选择当前不在集合source中具有最短路径的顶点u
		{
			if (!S[i] && dist[i]<min)
			{
				u = i;
				min = dist[i];
			}
		}
		//cout << u <<"  ";
		S[u] = 1;  //顶点u加入已求出最短路径的顶点集合
		for (int i = 0; i<MaxNodeNum; i++)//修改其余顶点的当前最短路径
		{
			//if ((!S[i]) && updateWeight[u][i]!=0 && updateWeight[u][i]< INFINITE && (dist[u]+updateWeight[u][i]< dist[i]))
			if ((!S[i]) && updateWeight[u][i]< INFINITE && (dist[u]+updateWeight[u][i]< dist[i]))
			{
				dist[i] = dist[u] + updateWeight[u][i];
				path[i] = u;
			}
		}
	}
	//cout << endl;
	// 下面是得到逆向的算路径结果，例如4<-3<-0,它表示源节点为0,目的节点为4, 对应的路径为4->3->0,是034吧
	if(destination==source)
	{
		routeTable_hopNum = 0;
		cout<< "path computation fails because destination==source. "<<endl<<endl;
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

		for(int k=0;k<MaxNodeNum;k++)
		{
			if(path[key]==source)
			{
				routeTable_hopNum = k+1;
				routeTable_nodeId[k+1] = source;
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
				routeTable_nodeId[k+1] = path[key];
				key = path[key];
				//cout<<routeTable_nodeId[k+1]<<endl;
				if (routeTable_nodeId[k+1] == -1)
				{
					cout<<" Setup path21 fail "<<endl<<endl;
					//return false;
					//break;
				}
			}
		}
	}
}

bool Suurballe::disjointPaths(int idConnectionRequestNum,int source,int destination) //基于dijkstra算法计算工作路径
{ 
	Prim();//用Prim算法计算拓扑的最小生成树,目的是在后面计算源节点	//到另外所有点的距离,和从源到宿的最路径

	//计算从源点到所有点的最短距离,并保存在costUpdate中,并且它包含了从源点到宿点的最短路径,
	// 存在	double priWeigth; vector<int> primaryPath; 并且double proWeigth;	 vector<int> protectPath没有值;
	for (int i = 0; i< MaxNodeNum; i++)
	{
		if (i != source)
		{
			nodePairElement pairElement;
			vector<int> pathList  = Dijkstra(source,i);

			if (pathList.size()>0)
			{
				pairElement.primaryPath = pathList;

				double pathWeigth = 0;
				for (int m = 0; m< (int)pairElement.primaryPath.size()-1;m++)
				{
					int k = pairElement.primaryPath.at(m);
					int l = pairElement.primaryPath.at(m+1);
					pathWeigth += weight[k][l];
				}
				pairElement.priWeigth = pathWeigth;
				pair<int, int> element = make_pair(source,i);
				costUpdate.insert(make_pair(element, pairElement));
				//cout<<"Setup path1 success"<<endl;
			}
		}
	}

	//以下部分是更新整个拓扑网络的链路的权值,方便另一条保护路径的计算
	double d[MaxNodeNum][MaxNodeNum];
	for (int i = 0; i < MaxNodeNum; i++)
	{
		for (int j = i; j < MaxNodeNum; j++)
		{
			d[j][i] = d[i][j] = 0;
		}
	}
	for (int i = 0; i<MaxNodeNum; i++)
	{
		if (source != i)
		{
			map<nodePair,nodePairElement>::iterator d_svu = costUpdate.find(make_pair(source,i));
			if (d_svu != costUpdate.end())
			{
				double dis = d_svu->second.priWeigth;
				d[source][i] = dis;
			}
		}
	}

	for (int i = 0; i < MaxNodeNum; i++)
	{
		for (int j = 0; j < MaxNodeNum; j++)
		{
			updateWeight[i][j] = weight[i][j];
		}
	}
	for (int u = 0; u<MaxNodeNum;u++)
	{
		for (int v = 0; v<MaxNodeNum; v++)
		{
			if (weight[u][v] != INFINITE )
			{
				updateWeight[u][v] = weight[u][v] - d[source][v] + d[source][u];
			}
		}
	}
	map<nodePair,nodePairElement>::iterator findPriPath = costUpdate.find(make_pair(source,destination));
	if (findPriPath != costUpdate.end())
	{
		nodePairElement elementPrim;
		elementPrim.primaryPath = findPriPath->second.primaryPath;
		elementPrim.priWeigth = findPriPath->second.priWeigth;
		dis_paths.insert(make_pair(idConnectionRequestNum,elementPrim));
		
		//cout<<elementPrim.priWeigth<<endl;

		for (int m = 0; m< (int)elementPrim.primaryPath.size()-1;m++)
		{
			int k = elementPrim.primaryPath.at(m);
			int l = elementPrim.primaryPath.at(m+1);
			updateWeight[k][l] = INFINITE;
			updateWeight[l][k] = 0;
		}
	}

	// 另一条路径的计算
	map<int,nodePairElement>::iterator insertPriPath = dis_paths.find(idConnectionRequestNum);
	vector<int> pathList = Dijkstra(source,destination);
	if (pathList.size()>0)
	{
		insertPriPath->second.protectPath = pathList;

		double proWeigth = 0;
		for (int i = 0; i< (int)pathList.size()-1;i++)
		{
			int k = pathList.at(i);
			int l = pathList.at(i+1);
			proWeigth += weight[k][l];
		}
		insertPriPath->second.proWeigth = proWeigth;
	}
	map<int,nodePairElement>::iterator findPath = dis_paths.find(idConnectionRequestNum);
	if (findPath != dis_paths.end())
	{
		return true;
	}
	else
	{
		cout<<"There are no disjoint paths"<<endl;
		return false;
	}
}


map<int, nodePairElement> Suurballe::combinePath(int idConnectionRequestNum,int source,int destination) //基于dijkstra算法计算工作路径
{ 
	vector<int> primaryPath;
	vector<int> protectPath;
	for(int i=0;i<MaxNodeNum;i++)
	{
		for(int j=i;j<MaxNodeNum;j++)
		{
			if (i != j)
			{
				updateWeight[j][i] = updateWeight[i][j] = INFINITE;
			}
		}	
	}
	map<int,nodePairElement>::iterator findPriProPath = dis_paths.find(idConnectionRequestNum);
	if (findPriProPath != dis_paths.end())
	{
		primaryPath = findPriProPath->second.primaryPath;
		protectPath = findPriProPath->second.protectPath;
		dis_paths.clear();
		for (int i = 0; i< (int)primaryPath.size()-1;i++)
		{
			int k = primaryPath.at(i);
			int l = primaryPath.at(i+1);
			updateWeight[k][l] = weight[k][l];
		}

		for (int i = 0; i< (int)protectPath.size()-1;i++)
		{
			int k = protectPath.at(i);
			int l = protectPath.at(i+1);
			updateWeight[k][l] = weight[k][l];
		}

		for(int i=0;i<MaxNodeNum;i++)
		{
			for(int j=i;j<MaxNodeNum;j++)
			{
				if (updateWeight[j][i] == updateWeight[i][j] && i != j)
				{
					updateWeight[j][i] = updateWeight[i][j] = INFINITE;
				}
			}	
		}
	}
	nodePairElement pairElement;
	pairElement.priWeigth = 0;
	pairElement.proWeigth = 0;
	pairElement.source = source;
	pairElement.destination = destination;

	int pathIndex = 0;
	while (pathIndex < 2)
	{
		//cout<<"pathIndex "<<pathIndex<<endl;

		vector<int> pathList = Dijkstra(source,destination);

		double  totalWeight = 0;
		for (int m = 0; m< (int)pathList.size()-1;m++)
		{
			int k = pathList.at(m);
			int l = pathList.at(m+1);
			totalWeight += weight[k][l];
			if (pathIndex == 0)
			{
				updateWeight[k][l] = INFINITE;
			}
		}
		if (pathIndex == 0)
		{
			pairElement.primaryPath = pathList;
			pairElement.priWeigth = totalWeight;
		}
		else
		{
			pairElement.protectPath = pathList;
			pairElement.proWeigth = totalWeight;
			pairElement.totalWeigth = pairElement.priWeigth;// + pairElement.proWeigth;
		}

		if (pathIndex == 1)
		{
			dis_paths.insert(make_pair(idConnectionRequestNum, pairElement));
		}

		pathIndex++;
	}

	map<int,nodePairElement>::iterator findPath = dis_paths.find(idConnectionRequestNum);

	if (findPath != dis_paths.end())
	{
/*
		cout<<"Primay path ";
		for(int i = 0; i<(int)findPath->second.primaryPath.size();i++)
		{
			cout<<findPath->second.primaryPath.at(i)<<" ";
		}
		cout<<endl;
		cout<<"Backup path ";
		for(int i = 0; i<(int)findPath->second.protectPath.size();i++)
		{
			cout<<findPath->second.protectPath.at(i)<<" ";
		}
		cout<<endl;*/
		return dis_paths;
	}
	else
	{
		cout<<"There are no disjoint paths"<<endl;
	}
}

map<nodePair, nodePairElement> Suurballe::allDisPath(int idConnectionRequestNum)
{
	//readTopology(MaxNodeNum, "6nodes-KM.txt");
	readTopology(MaxNodeNum, "NSF_KM.txt");
	//graph.readTopology(MaxNodeNum, "USNET_KM.txt");
	//graph.readTopology(MaxNodeNum, "testNet.txt");

	map<int, nodePairElement> oneNodePair;
	map<nodePair, nodePairElement> allNodePair;
	nodePairElement pairElement;

	if ((int)oneNodePair.size() != 0)
	{
		allNodePair.clear();
	}
	if ((int) allNodePair.size() != 0)
	{
		allNodePair.clear();
	}
	for (int source = 0; source<MaxNodeNum; source++)
	{
		for (int destination = 0; destination<MaxNodeNum; destination++)
		{
			if (source != destination )
			{
				disjointPaths(idConnectionRequestNum,source,destination);
				oneNodePair = combinePath(idConnectionRequestNum,source,destination);
				dis_paths.clear();
				map<int,nodePairElement>::iterator findDisPath = oneNodePair.find(idConnectionRequestNum);
				if (findDisPath != oneNodePair.end())
				{
					pairElement = findDisPath->second;
					allNodePair.insert(make_pair(make_pair(source,destination), pairElement));
					//cout<<"source "<<source<<endl;
				}
				oneNodePair.clear();
			}
		}
	}
	return allNodePair;
}
#ifndef _Suurballe_H_
#define _Suurballe_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <queue>

#include "Param.h"

using namespace std;



struct edge
{
	int fromvex;
	int endvex;
	double weight;
};
//struct edge CT[MaxNodeNum],temp;

typedef pair<int,int> nodePair;


class Suurballe
{
public:

	struct edge CT[MaxNodeNum],temp;
	map<nodePair,nodePairElement> costUpdate;
	map<int, nodePairElement> dis_paths;

	void readTopology(int, char*);

	int frand(int,int);
	void randomSrcDst(int &, int &);
	void Prim();
	vector<int> Dijkstra(int,int);
	bool dijkstra_1(int,int);
	bool disjointPaths(int,int,int); 
	map<int, nodePairElement> combinePath(int,int,int);
	map<nodePair, nodePairElement> allDisPath(int);
};

#endif
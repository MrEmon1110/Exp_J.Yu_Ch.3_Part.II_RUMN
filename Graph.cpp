#include "Param.h"

#include "Graph.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
#include <math.h>
#include <stdlib.h>   //随机数产生头文件
#include <set>
using namespace std;
using namespace lch::graph;



Graph::Graph() : m_initSize(0), m_size(0)
{
}

Graph::Graph(const Adjacency_list& adjMap)
{
	m_initSize = m_size = adjMap.size();
	m_adj_map = adjMap;
}


Graph::~Graph()
{
	delete_adj_list(m_adj_map);
	delete_adj_list(m_in_map);
	delete_adj_list(m_init_adj_map);
	delete_adj_list(m_init_in_map);
}

void Graph::readTopology(unsigned size, char* name)
{
	for(int i=0;i<(int)size;i++)
	{
		for(int j=0;j<(int)size;j++)
		{
			updateWeight[i][j] = weight[i][j] = 0;
		}	
	}

	ifstream infile;
	infile.open(name, ios::in);
	if(!infile)
	{
		cerr << "open error!" << endl;
		exit(-2);
	}
	for(int i=0;i<(int)size;i++)
	{
		for(int j=0;j<(int)size;j++)
		{
			WEIGHTVALUE  weightValue;
			infile >> weightValue;
			updateWeight[i][j] = weight[i][j] = weightValue;
		}	
	}
	infile.close();
}

/*
void Graph::failProbability()
{
	double failPr = 0;
	for(int i=0;i<MaxNodeNum;i++)
	{
		for(int j= 0;j<MaxNodeNum;j++)
		{
			if (updateWeight[i][j] != 0 && updateWeight[i][j] < INFINITE)
			{
				double random = 0.0;
				double failPr = 0.0;
				while(1)//每条链路故障的概率为[0，FP]
				{

					random = (double)rand()/RAND_MAX;
					//if (random >= LFP1 && random <= LFP2)
					if (random <= LFP)
					{
						failPr = random;
						break;
					} 
				}
				//updateWeight[j][i] = weight[j][i] = updateWeight[i][j] = weight[i][j] = -log(1.0-failPr);
				updateWeight[i][j] = weight[i][j]= -log(1.0-failPr);
			}
		}	
	}
	ofstream outfile;
	outfile.open("Arcs.txt", ios::app);
	//outfile.open("Arcs.txt", ios::out);

	if(!outfile)
	{
		cerr << "存取数据的文件打不开" << endl;
	}
	else
	{
		outfile <<"K_SP  "<<K_SP <<"  conNum "<<ithTime<<endl;
		for(int i=0;i<MaxNodeNum;i++)
		{
			for(int j= 0;j<MaxNodeNum;j++)
			{
				if (updateWeight[i][j] != 0 && updateWeight[i][j] < INFINITE)
				{
					//<"1","4", 2,5.5>
					outfile <<"<\""<<i<<"\", \""<<j<<"\",  "<<fixed<<setprecision(8)<<updateWeight[i][j]<<">"<<endl;
					//cout<<updateWeight[j][i]<<" ";
				}
			}
			//cout<<endl;
		}
		outfile <<endl;
	}
	outfile.close();
}*/
void Graph::createLinkTab()
{
	unsigned size;
	m_size = m_initSize = size = MaxNodeNum;
	for (int i=0;i<(int)size;i++)
	{
		Vertex ver(i, 0);
		m_adj_map.push_back(ver);
		m_in_map.push_back(ver);
	}
	for(int i=0;i<(int)size;i++)
	{
		for(int j=0;j<(int)size;j++)
		{
			WEIGHTVALUE weightValue;
			weightValue = updateWeight[i][j];
			if (weightValue > 0 && weightValue < INFINITE - 1 )
			{
				Vertex* vertex1 = new Vertex(j, weightValue);
				Vertex* vertex2 = new Vertex(i, weightValue);
				vertex1->m_next = m_adj_map[i].m_next;
				m_adj_map[i].m_next = vertex1;
				vertex2->m_next = m_in_map[j].m_next;
				m_in_map[j].m_next = vertex2;
			}				
		}	
	}
	m_init_in_map = m_in_map;
	m_init_adj_map = m_adj_map;
}


void Graph::delete_adj_list(Adjacency_list& map)
{
	for (int i=0;i<(int)map.size();i++)
	{
		Vertex* v = map.at(i).m_next;
		Vertex* p = NULL;
		while(v!=NULL)
		{
			p = v;
			v = v->m_next;
			delete p;
		}
	}
}

void Graph::dijkstra(int srcId, WEIGHTVALUE* dist, int* path)
{
	//cout<<"srcId "<< srcId <<"dist: "<<*dist<<endl;
	for(int i=0;i<(int)m_size;i++)
	{
		dist[i] = INFINITE;
		path[i] = -1;
	}
	dist[srcId] = 0;
	bool *done = new bool[m_size];
	for (int i=0;i<(int)m_size;i++)
	{
		done[i] = false;
	}
	//memset(done, 0, m_size);
	priority_queue < Edge > pq;
	Edge edge(srcId, 0);
	pq.push(edge);
	while(!pq.empty())
	{
		Edge out = pq.top();
		pq.pop();
		if(done[out.m_id]==true)
		{
			continue;
		}
		done[out.m_id] = true;
		dist[out.m_id] = out.m_metric;
		for(Vertex* vertex=m_adj_map[out.m_id].m_next;vertex!=NULL;vertex=vertex->m_next)
		{
			if(done[vertex->m_id]==false)
			{
				if(vertex->m_weightValue + out.m_metric < dist[vertex->m_id])
				{
					dist[vertex->m_id] = vertex->m_weightValue + out.m_metric;
					path[vertex->m_id] = out.m_id;
					Edge in(vertex->m_id, dist[vertex->m_id]);
					pq.push(in);
				}
			}
		}
	}
	delete[] done;
}

double Graph::dijkstra( int start , int end , Path_list& path_list )
{
	int* paths = new int[m_size];
	WEIGHTVALUE* dist = new WEIGHTVALUE[m_size];
	dijkstra(start, dist, paths);
	if(dist[end] > INFINITE -1) 
	{ 
		delete[] paths;
		delete[] dist;
		return -1;
	}
	// parse the shortest path
	int i = end;	
	while(i>=0)
	{		
		path_list.push_front(i);
		i=paths[i];			
	}
	WEIGHTVALUE min = dist[end];
	delete []paths;
	delete []dist;
	return min;
}
int Graph::ksp( int src, int dest, unsigned k , Path_vector& kpaths )
{	
	/* Initialize */	
	vector<int> Path, Prime, Base;
	vector<double> Dist;
	Path.resize(m_size); 
	Prime.resize(m_size);
	Base.resize(m_size);
	Dist.resize(m_size);
	for(int i=0;i<(int)m_size;i++)
	{ 
		Prime[i] = -1;
		Base[i] = i; 
	}

	/* Find the shortest path firstly */	
	dijkstra(src, &Dist[0], &Path[0]);
	if( Dist[dest] > INFINITE -1 ) return 0;
	Path_list path; 
	int j = Path[dest];	
	path.push_back(dest);
	while(j >= 0)
	{		
		path.push_front(j); j=Path[j];			
	}
	kpaths.push_back(path); // store the shortest path

	/* Find the 2th - kth shortest paths */
	int ki = 1,kj = 1;
	while(ki < (int)k && kj < 50)
	{
		/* Find the first node with more than a single incoming arc */
		unsigned int nh = -1;
		while( path.size() )
		{
			unsigned node = path.front(); 
			path.pop_front();
			int count = 0;
			Vertex* vertex = &m_in_map[node];
			while(vertex->m_next!=NULL)
			{
				count++;
				vertex = vertex->m_next;
				if( count > 1 ) break;
			}			
			if( count > 1 ) 
			{ 
				nh = node; 
				break; 
			}
		}

		if( nh == -1 ) break; // there is NOT an alternative path, exit!

		int ni = -1;
		/* Add the first prime node to graph */
		if( Prime[nh] < 0 )
		{
			unsigned nh1 = addNode(nh,Path[nh]);

			/* compute the minimal distance from node 0 to nh1 */
			double min_dist = INFINITE;
			int min_node = -1;
			for(Vertex* ver = m_in_map[nh1].m_next; ver != NULL ; ver = ver->m_next)
			{
				//cout << Dist[i] << " " << m_adj_map[i][nh1] << endl; // for debug
				int id = ver->m_id;
				WEIGHTVALUE wei = ver->m_weightValue;
				if( Dist[id] + wei < min_dist )
				{
					min_dist = Dist[id] + wei;
					min_node = id;
				}
			}			
			Dist.push_back(min_dist);
			Path.push_back(min_node);
			Prime.push_back(-1);
			Prime[nh] = nh1;

			/* record the base node */
			unsigned basei = nh;
			while(basei != Base[basei])
				basei = Base[basei];
			Base.push_back(basei);

			if(path.size())
			{ 
				ni = path.front(); 
				path.pop_front(); 
			}
		}
		/*  Get node ni, it must meet it's the first node following nh in path, but its prime node ni` is NOT in graph */
		else
		{
			while( path.size() )
			{
				ni = path.front(); path.pop_front();
				if(Prime[ni] < 0) break;				
			}
		}		

		/* Add the other prime nodes to graph */
		while(true && ni != -1)//修改
		{
			unsigned ni1 = addNode(ni,Path[ni]);
			int temp1 = Path[ni];
			int temp2 = Prime[temp1];
			if( temp2 >= 0 )	
			{
				WEIGHTVALUE wei = Dist[ni] - Dist[temp1];
				Vertex* ver1 = new Vertex(ni1, wei);
				Vertex* ver2 = new Vertex(temp2, wei);
				ver1->m_next = m_adj_map[temp2].m_next;
				m_adj_map[temp2].m_next = ver1;
				ver2->m_next = m_in_map[ni1].m_next;
				m_in_map[ni1].m_next = ver2;
			}
			/* compute the minimal distance from node 0 to ni1 */
			double min_dist = (unsigned)INFINITE;
			int min_node = -1;

			for(Vertex* ver = m_in_map[ni1].m_next; ver != NULL ; ver = ver->m_next)
			{
				//cout << Dist[i] << " " << m_adj_map[i][nh1] << endl; // for debug
				int id = ver->m_id;
				WEIGHTVALUE wei = ver->m_weightValue;
				if( Dist[id] + wei < min_dist )
				{
					min_dist = Dist[id] + wei;
					min_node = id;
				}
			}

			Dist.push_back(min_dist);
			Path.push_back(min_node);
			Prime.push_back(-1);
			Prime[ni] = ni1;

			/* record the base node */
			unsigned basei = ni;
			while(basei != Base[basei])
				basei = Base[basei];
			Base.push_back(basei);

			if( !path.size() ) break;
			ni = path.front(); 
			path.pop_front();			
		}

		/* get the kth shortest path */			
		if( ni==-1 ) ni = nh; // if nh is just the end node.
		Path_list temp;		
		int j = Prime[ni];
		while(j>=0)
		{		
			path.push_front(j); temp.push_front(Base[j]); j=Path[j];
		}
		if(temp.size()<2) break; 
		Path_list::iterator it1, it2;
		bool reg = true;
		for (it1=temp.begin();it1!=temp.end();it1++)
		{
			for (it2=it1,it2++;it2!=temp.end();it2++)
			{
				if (*it1 == *it2)
				{
					reg = false;
					break;
				}
			}
			if (!reg)
			{
				break;
			}
			//cout<<*it1<<" ";
		}
		//cout<<endl;
		if (reg)
		{
			kpaths.push_back(temp);  // store the kth shortest path
			ki++;
		}
		kj++;
		//this->Output(); // for debug
	}
	return ki;
}


/*
Look out: Dijkstra algorithm doesn't assure the
generated tree from graph is minimal generated tree.
*/
/// Look out: the returned paths doesn'n include the start node!


/*
Here, using vector is more efficient than 2-dimension dynamic array.
*/
unsigned Graph::addNode(unsigned ni,int preni)
{
	Vertex vertex(m_size, 0);
	m_adj_map.push_back(vertex);
	m_in_map.push_back(vertex);
	for(Vertex* ver = m_in_map[ni].m_next;ver != NULL;ver = ver->m_next)
	{
		int id = ver->m_id;
		WEIGHTVALUE wei = ver->m_weightValue;
		if(id != preni)
		{
			Vertex* vertex2 = new Vertex(m_size, wei);
			vertex2->m_next = m_adj_map[id].m_next;
			m_adj_map[id].m_next = vertex2;

			Vertex* vertex3 = new Vertex(id, wei);
			vertex3->m_next = m_in_map[m_size].m_next;
			m_in_map[m_size].m_next = vertex3;
		}
	}
	return m_size++;
}

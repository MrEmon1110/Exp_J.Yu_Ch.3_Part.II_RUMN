/*===============================================================================
<summary>
<filename> Graph.h </filename>
<copyright> Copyright (c) 2006 David D. All Rights Reserved. </copyright> 
<guide>
Graph Data Structure.			

Look out: If the input graph isn't directed acyclic
graph, the MSAforKSP maybe throw some exceptions.
What's more, the input graph must have only one START
node. If the origin graph has more than one START node,
you need add a 'virtual' START node to input graph 
which should be linked to the real START nodes. To the 
END node, do the same as START node.

注意：在使用MSAforKSP求解前k条最短路径时，必须确保
输入的有向图是简单有向图（没有环，路径权重非负），
并且必须有且只有一个开始节点，如果原始图有多个开始
节点的话，就需要添加一个虚的开始节点来连接原始多个
开始节点。同样，必须添加一个虚的结束节点连接所有原
始结束节点。

* Any suggestion, contact with SharperDavid@hotmail.com.
</guide>
<version> 1.0 </version>
<author> David D. </author>
<date> 2006.2 </date>
<summary>
=================================================================================*/

#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <list>

typedef double WEIGHTVALUE;
using namespace std;


//extern double weight[MaxNodeNum][MaxNodeNum]; 
//extern double updateWeight[MaxNodeNum][MaxNodeNum];


namespace lch
{
	namespace graph
	{
		class Vertex
		{
		public:
			Vertex& operator = (const Vertex& _graph)
			{
				m_id = _graph.m_id;
				m_weightValue = _graph.m_weightValue;
				if (_graph.m_next != NULL)
				{
					m_next = new Vertex;
					*m_next = *_graph.m_next;
				}
				else
				{
					m_next = NULL;
				}
				return *this;
			}
			Vertex(const Vertex& ver)
			{
				*this = ver;
			}
		private:
			friend class Edge;
			friend class Graph;
			Vertex()
			{
				m_next = NULL;
			}
			Vertex(int id, WEIGHTVALUE w)
			{
				m_next = NULL;
				m_id = id;
				m_weightValue = w;
			}

			int m_id;
			Vertex* m_next;
			WEIGHTVALUE m_weightValue;
		};
		class Edge
		{
		public:
			Edge(){}
			Edge(int id, WEIGHTVALUE w)
			{
				m_id = id;
				m_metric = w;
			}
			bool operator < (const Edge& edge) const
			{
				return this->m_metric > edge.m_metric;
			}
			int m_id;
			WEIGHTVALUE m_metric;
		};

		class Graph
		{
		public:
			typedef list<int> Path_list;			// the top is the first node
			typedef vector<Path_list>  Path_vector;		// the paths are sorted by ASC
			typedef vector<Vertex> Adjacency_list; // the two dimension array storing a DAG
			void clear()
			{
				delete_adj_list(m_adj_map);
				delete_adj_list(m_in_map);
				delete_adj_list(m_init_adj_map);
				delete_adj_list(m_init_in_map);
				m_adj_map.clear();
				m_in_map.clear();
				m_init_adj_map.clear();
				m_init_in_map.clear();
			}

			/// Input array: there must be only a start node!
			///	if the value of one element of array < 0, indicates there are no arc.
			Graph();
			Graph( const Adjacency_list& );
			Graph( const double **array , unsigned N );
			~Graph();
			/// Reset the directed acyclic graph (DAG).
			void Restart( const Adjacency_list& array );
			void Restart( const double **array , unsigned N );
			void readTopology(unsigned, char*);
			//void ERtopology();
			//void failProbability();
			void createLinkTab();
			/// Find the shortest path from node 0 to node N-1.
			/// Dijkstra Algorithm Implementation.			
			/// Parameter <path> return the shortest path.
			///	If fails, return -1, or return the shortest distance.			
			double Dijkstra( Path_list& path );
			/// Find the shortest path from the named start node to the named end node.
			/// Parameter <start> - start node.
			/// Parameter <end> - end node. 
			/// Make sure <end> != <start>.			
			double Dijkstra( unsigned start , unsigned end , Path_list& path );

			/// Find the k shortest paths (KSP) from node 0 to N-1.						
			/// Martins' Algorithm (deletion algorithm) Implementation.
			/// Parameter <paths> return all the shortest paths.			
			/// If fails, return 0, or return the real number of all the shortest paths.								
			int ksp(int src, int dest, unsigned k , Path_vector& kpaths );
			void dijkstra( int srcId, WEIGHTVALUE* dist, int* path);
			WEIGHTVALUE dijkstra( int start , int end , Path_list& path_list );

			/// Output the content of "m_adj_map" for debug.
			//void Output( ostream& out = cerr );

		private:
			/// Default is to compute the shortest distance from node 0 to node N-1.

			/// Add a node to graph.
			/// Return the number of new node.
			unsigned addNode( unsigned ni , int preni );		
			void delete_adj_list(Adjacency_list&);
		private:
			size_t m_initSize; // original size of "m_adj_map", it's fixed.
			size_t m_size; // size of "m_adj_map", because the "m_adj_map" maybe be reallocated. 
			Adjacency_list m_adj_map; // 
			Adjacency_list m_init_adj_map;
			Adjacency_list m_in_map;
			Adjacency_list m_init_in_map;				
		};
	}
}

#endif	// GRAPH_H

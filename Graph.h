#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <iostream>
#include <list>
#include <stack>
#include <queue>

class Graph {
public:
	Graph(int size);
	void addEdge(int, int);
	void print();
//	void transpose() ;

	void dfs(int);
	void dfsRec(int, std::vector<bool>&);
	void dfsIterative(int);
	void dfs_for_getComponent(int, std::vector<bool>&);

	int getComponentCount();

	void bfsRecursive(int);
	void bfsRec(std::vector<bool>&, std::queue<int>);
	void bfsIter(int);

	int shortestDistenceInUnweightedGraph(int, int);

	std::vector<std::vector<int>> allPaths(int u, int v);

    void dfs_for_allPath(int u, int v, std::vector<bool>& visited,
                         std::vector<std::vector<int>>& allPaths,
                         std::vector<int>& path);
	

	//undirect graph with parenting 
	bool hasCycle();
	bool dfsCycle(int, std::vector<bool>& visited,int);

	bool hasCycleInDirectedGraph();
	bool dfsCycleInDirectedGraph(int, std::vector<bool>& visited,std::vector<bool>&);

	std::vector<int> toplogicalSorting();
	void topSortRec(std::vector<int>&, std::vector<bool>&, int); 
	
	//Khan's Algorithm
	std::vector<int> KhansAlgorithm();

	//Kosoraju's Algorithm
	std::vector<std::vector<int>> Kosoraju();
	void dfs_kosoraju(std::vector<bool>& , int, std::stack<int>&);
	void dfs_transpose(int , std::vector<bool>&, std::vector<int>& component, std::vector<std::vector<int>>& transposed);
	std::vector<std::vector<int>> transpose();

	//Tarjan's Algorithm
	std::vector<std::vector<int>> Tarjan();
	void dfs_Tarjan(int u, std::vector<int>& disc, std::vector<int>& low, std::stack<int>& st, std::vector<bool>& onStack, std::vector<std::vector<int>>& result) ;

private:
	std::vector<std::vector<int>> adjMatrix;
};
#endif //GRAPH_H

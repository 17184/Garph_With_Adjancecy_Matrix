#include <queue>
#include <algorithm>
#include "Graph.h"

Graph::Graph(int size) {
	adjMatrix = std::vector<std::vector<int>>(size, std::vector<int>(size, 0));
}
void Graph::addEdge(int u, int v) {
	adjMatrix[u][v] = 1;
//	adjMatriv[u][u] = 1; for undirected graph
}
void Graph::print() {
		for(const auto & row : adjMatrix) {
			for(int val : row) {
				std::cout << val << " ";
			}
			std::cout << std::endl;
		}
}
/*
 //i comment it for else transpose in Kosoraju's algorithm
void Graph::transpose() {
	int currSize = adjMatrix.size();	
	//same as matrix transpose row will be column and column will be row
	std::vector<std::vector<int>> tmpMatrix(currSize, std::vector<int>(currSize, 0));
	for(int i = 0; i < currSize; ++i) {
		for(int j = 0; j < currSize; ++j) {
			tmpMatrix[i][j] = adjMatrix[j][i];
		}
	}
	adjMatrix = tmpMatrix;
}
*/
void Graph::dfs(int startVertex) {
	std::vector<bool> visited(adjMatrix.size(), false);
	dfsRec(startVertex, visited);
}
void Graph::dfsRec(int startVertex, std::vector<bool>& visited) {
	int size = adjMatrix.size();
	visited[startVertex] = true;
	std::cout << startVertex << " ";
	for(int i = 0; i < size; ++i) {
		if(visited[i] == false && adjMatrix[startVertex][i] == 1) {
			//it is adjucent and not visited
			dfsRec(i,visited);
		}
	}
}

void Graph::dfsIterative(int startVertex) {
	std::stack<int> s;
	std::vector<bool> visited(adjMatrix.size(), false);
	s.push(startVertex);
	visited[startVertex] = true;
	while(!s.empty()) {
		int curr = s.top();
		s.pop();
		std::cout << curr << " ";
		for(int i = 0; i < adjMatrix[curr].size(); ++i) {
			if(!visited[i] && adjMatrix[curr][i] == 1) {
				s.push(i);
				visited[i] = true;
			}
		}	
	}
}

int Graph::getComponentCount() {
    int currSize = adjMatrix.size();
    int componentCount = 0;
    std::vector<bool> visited(currSize, false);

    for (int i = 0; i < currSize; ++i) {
        if (!visited[i]) {
            dfs_for_getComponent(i, visited);
            componentCount++;
        }
    }

    return componentCount;
}
void Graph::dfs_for_getComponent(int startNode, std::vector<bool>& visited) {
    std::stack<int> s;
    s.push(startNode);

    while (!s.empty()) {
        int curr = s.top();
        s.pop();

        if (!visited[curr]) {
            visited[curr] = true;

            for (int i = 0; i < adjMatrix.size(); ++i) {
                if (!visited[i] && adjMatrix[curr][i] == 1) {
                    s.push(i);
                }
            }
        }
    }
}

void Graph::bfsRecursive(int u) {
	std::vector<bool> vis(adjMatrix.size(), false);
	std::queue<int> q;
	q.push(u);
	vis[u] = true;
		bfsRec(vis, q);
	
}

void Graph::bfsRec(std::vector<bool>& vis, std::queue<int> q) {
	if(q.empty()) { return; }
	int curr = q.front();
	q.pop();
	std::cout << curr << " ";

	for(int i = 0; i < adjMatrix[curr].size(); ++i) {
			if(!vis[i] && adjMatrix[curr][i]) {
				q.push(i);
				vis[i] = true;
		}
	}
	bfsRec(vis, q);
} 

void Graph::bfsIter(int startVertex) { 
	std::queue<int> q;
	std::vector<bool> visited(adjMatrix.size(), false);
	q.push(startVertex);

	while(!q.empty()) {
		int curr = q.front();
		std::cout << curr << " ";
		q.pop();
		visited[curr] = true;
		for(int i = 0; i < adjMatrix.size(); ++i) {
			if(!visited[i] && adjMatrix[curr][i]) {
				q.push(i);	
				visited[i] = true;
			}
		}
	}
}

int Graph::shortestDistenceInUnweightedGraph(int u, int v) {
	std::vector<bool> visited(adjMatrix.size(), false);
	std::vector<int> distance(adjMatrix.size(), -1);
	visited[u] = true;
	distance[u] = 0;
	
	std::queue<int> q;
	q.push(u);
	
	while(!q.empty()) {
		int curr = q.front();
		q.pop();
		for(int i = 0; i < adjMatrix.size(); ++i) {
			if(!visited[i] && adjMatrix[curr][i]) {
				q.push(i);
				visited[i] = true;
				distance[i] = distance[curr] + 1;

			if(i == v) {
					return distance[i];
				} 
			}
		}
	}
	return -1;
}
std::vector<std::vector<int>> Graph::allPaths(int u, int v) {
    std::vector<bool> visited(adjMatrix.size(), false);
    std::vector<int> path;
    std::vector<std::vector<int>> allPaths;
    dfs_for_allPath(u, v, visited, allPaths, path);
    return allPaths;
}

void Graph::dfs_for_allPath(int u, int v, std::vector<bool>& visited,
                             std::vector<std::vector<int>>& allPaths,
                             std::vector<int>& path) {
    visited[u] = true;
    path.push_back(u);

    if (u == v) {
        allPaths.push_back(path);
    } else {
        for (int i = 0; i < adjMatrix.size(); ++i) {
            if (!visited[i] && adjMatrix[u][i]) {
                dfs_for_allPath(i, v, visited, allPaths, path);
            }
        }
    }

    visited[u] = false;
    path.pop_back();
}

bool Graph::hasCycle() {
	std::vector<bool> visited(adjMatrix.size() , false);

	for(int i = 0; i < adjMatrix.size(); ++i) {
		if(!visited[i]) {
			if(dfsCycle(i, visited, -1)) {
			return true;
		}
	}
}
	return false;

}
//for undirected graph
bool Graph::dfsCycle(int vertex, std::vector<bool>& visited, int parent) {
	visited[vertex] = true;
	for(int neight = 0; neight < adjMatrix.size(); ++neight) {
		if(adjMatrix[vertex][neight]) {
			if(!visited[neight]) {
				if(dfsCycle(neight, visited, vertex)) {
					return true;
				}
			} else if(parent != neight) {
				return true;
		}
	}
	
}
	return false;
}
bool Graph::hasCycleInDirectedGraph() {
	std::vector<bool> visited(adjMatrix.size(), false);
	std::vector<bool> onStack(adjMatrix.size(), false);
	for(int i = 0; i < adjMatrix.size(); ++i) {
		if(!visited[i] && dfsCycleInDirectedGraph(i, visited, onStack)) {
			return true;
		}
	}
	return false;
}
bool Graph::dfsCycleInDirectedGraph(int vertex, std::vector<bool>& visited, std::vector<bool>& onStack) {
	visited[vertex] = true;
	onStack[vertex] = true;
	for(int i = 0; i < adjMatrix.size(); ++i) {
		if(adjMatrix[vertex][i]) {
			if(!visited[i]) {
				if(dfsCycleInDirectedGraph(i, visited, onStack)) {
					return true;
				}
			}
			 else if(onStack[i]) {
				return true;
		} 		}
	}
	onStack[vertex] = false;
	return false;
}

std::vector<int> Graph::toplogicalSorting() {
	std::vector<int> topOrder;


	std::vector<bool> vis(adjMatrix.size(), false);
	std::vector<int> stack;
	for(int i = 0; i < adjMatrix.size(); ++i) {
		if(!vis[i]) {
			topSortRec(stack, vis, i);
		}	
	}
//		std::vector<int> topOrder;
			while(!stack.empty()) {
				topOrder.push_back(stack.back());
				stack.pop_back();
		}
	return topOrder;

//	std::cout << "TopSort working only with DAG\n";
//	return topOrder;
	//}
}
void Graph::topSortRec(std::vector<int>& stack, std::vector<bool>& vis, int u) {
	vis[u] = true;
	
	for(int i = 0; i < adjMatrix[u].size(); ++i) {
		if(!vis[adjMatrix[u][i]]) {
			topSortRec(stack, vis, adjMatrix[u][i]);
		}
	}
	stack.push_back(u);
}

//topSOrting with Khan's Algorithm
std::vector<int> Graph::KhansAlgorithm() {
	std::vector<int> inDegree(adjMatrix.size(), 0);
	std::vector<int> topOrder;
	
	for(const auto & neight :  adjMatrix) {
		for(int i  : neight){
			inDegree[i]++;
	}
	}

	std::queue<int> q;
	
	for(int i = 0; i < adjMatrix.size(); ++i) {
		if(inDegree[i] == 0) {
			q.push(i);
		}
	}

	while(!q.empty()) {
		int curr = q.front();
		q.pop();
		topOrder.push_back(curr);

		for(int i : adjMatrix[curr]) {
			inDegree[i]--;
			if(inDegree[i] == 0) {
				q.push(i);
			}
		}
	}

	if(topOrder.size() != inDegree.size()) {
		std::cout << "Its not DAG top sorting can be only with DAG\n";
		return {};
	}
	return topOrder;
}



//Kosoraju's Algorithm
std::vector<std::vector<int>> Graph::Kosoraju() {
	std::vector<bool> vis(adjMatrix.size(), false);
	std::stack<int> st;
	
	for(int i = 0; i < adjMatrix.size(); ++i) {
		if(!vis[i]) {
			dfs_kosoraju(vis, i, st);
		}
	}
	std::fill(vis.begin(), vis.end(), false);
	std::vector<std::vector<int>> gt = transpose();
	std::vector<std::vector<int>> sccs;

//	std::fill(vis.begin(), vis.end(), false);

	while(!st.empty()) {
		int tp = st.top();
		st.pop();
		if(!vis[tp]) {
			std::vector<int> comp;
			dfs_transpose(tp, vis, comp, gt);
			sccs.push_back(comp);
		}
		}
	return sccs;
}
void Graph::dfs_kosoraju(std::vector<bool>& vis, int v, std::stack<int>& s) {
	vis[v] = 1;
	for(int i = 0; i < adjMatrix.size(); ++i) {
		if(!vis[i] && adjMatrix[v][i]) {
			dfs_kosoraju(vis, i, s);	
		}
	}	
	s.push(v);
}

void Graph::dfs_transpose(int v, std::vector<bool>& vis, std::vector<int>& component, std::vector<std::vector<int>>& transposed) {
	vis[v] = 1;
	component.push_back(v);	
	for(int i = 0; i < adjMatrix.size(); ++i) {
		if(!vis[i] && transposed[v][i]) {
			dfs_transpose(i, vis, component, transposed);
		}
	}
}
std::vector<std::vector<int>> Graph::transpose() {
	int v = adjMatrix.size();
	std::vector<std::vector<int>> t(v, std::vector<int>(v, 0));	
	for(int i = 0; i < adjMatrix.size(); ++i) {
		for(int j = 0; j < adjMatrix.size(); ++j) {
			t[i][j] = adjMatrix[j][i];
		}		
	}
	return t;
}

//Tarjan's Algorihtm
std::vector<std::vector<int>> Graph::Tarjan() {
	int numVertices = adjMatrix.size(); 
     std::vector<int> disc(numVertices, -1);  // Discovery time
    std::vector<int> low(numVertices, -1);   // Low link value
    std::stack<int> st;  
	                    // Stack for DFS traversal
    std::vector<bool> onStack(numVertices, false); // To check if a vertex is on the stack
    std::vector<std::vector<int>> result;    // Resulting strongly connected components

    for (int i = 0; i < numVertices; ++i) {
        if (disc[i] == -1) {
            dfs_Tarjan(i, disc, low, st, onStack, result);
        }
    }

    return result;
}

void Graph::dfs_Tarjan(int u, std::vector<int>& disc, std::vector<int>& low, std::stack<int>& st, std::vector<bool>& onStack, std::vector<std::vector<int>>& result) {
    static int time = 0;  
	int numVertices = adjMatrix.size();
    disc[u] = low[u] = ++time;
    st.push(u);
    onStack[u] = true;

    for (int v = 0; v < numVertices; ++v) {
        if (adjMatrix[u][v]) {
            if (disc[v] == -1) {
                dfs_Tarjan(v, disc, low, st, onStack, result);
                low[u] = std::min(low[u], low[v]);
            } else if (onStack[v]) {
                low[u] = std::min(low[u], disc[v]);
            }
        }
    }

    if (low[u] == disc[u]) {
        std::vector<int> component;
        for(int i = st.top(); ; i = st.top()) {
            int v = st.top();
            st.pop();
            onStack[v] = false;
            component.push_back(v);
            if (u == v) {
                break;
            }
        }
    //    std::reverse(component.begin(), component.end());  // Reverse to get correct order
        result.push_back(component);
    }
}

int main(){

	Graph graph(7);

    graph.addEdge(1,2);
    graph.addEdge(2,3);
    graph.addEdge(3,2);
    graph.addEdge(3,4);
	graph.addEdge(1,4);
	graph.addEdge(0,1);
	graph.addEdge(6,0);
 	graph.addEdge(1,6);
    graph.addEdge(6,2);
	graph.addEdge(4,5);
	graph.addEdge(5,4);
	graph.addEdge(3,5);

	std::cout << "Kosoraju's Algorihtm \n";
	std::vector<std::vector<int>> stronglyConnectedComponents = graph.Kosoraju();
	 for (const auto& component : stronglyConnectedComponents) {
        std::cout << "Component: ";
        for (int vertex : component) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }

	std::cout << std::endl;
	std::cout << "Tarjan's Algorithm \n";
	std::vector<std::vector<int>> TarjanSccs = graph.Tarjan();
	for(const auto & i : TarjanSccs) {
		std::cout << "Component ";
		for(int v : i) {
			std::cout << v << " ";
		}
	std::cout << std::endl;
	}

	
	
/*    std::cout << "BFS recursive starting from vertex 0: ";
    graph.bfsRecursive(0);
    std::cout << std::endl;

  	std::cout << "BFS ierative starting from vertex 0: ";
    graph.bfsIter(0);
    std::cout << std::endl;

 	std::cout << "DFS Recursive starting from vertex 0: ";
    graph.dfs(0);
    std::cout << std::endl;

	std::cout << "DFS Iterative starting from vertex 0: ";
    graph.dfsIterative(0);
    std::cout << std::endl;

	Graph g(6);
    g.addEdge(5, 2);
    g.addEdge(5, 0);
    g.addEdge(4, 0);
    g.addEdge(4, 1);
    g.addEdge(2, 3);
    g.addEdge(3, 1); 
*/
/*	std::vector<int> resTopSOrt = graph.toplogicalSorting();
	std::cout << "Topological Sorting result - ";
    std::vector<int> order = graph.toplogicalSorting();
	for(int i = 0; i < order.size(); ++i) {
		std::cout << order[i] << " ";
	}

	
	std::cout << "\nKhan's ALgorithm result - ";
  	std::vector<int> orderKhans = g.toplogicalSorting();
	 if (!orderKhans.empty()) {
        std::cout << "Topological Order: ";
        for (int node : orderKhans) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }

*/

/*	Graph g1(3);
	g1.addEdge(1, 0);
	g1.addEdge(2, 1);
	g1.addEdge(0,2);
*/
/*	g.addEdge(1,2);
	g.addEdge(2,4);
	g.addEdge(3,4);*/
//	g.addEdge(1,4);
/*	g.addEdge(1,3);
	g.addEdge(2,4);*/
/*    std::vector<std::vector<int>> allPaths = g.allPaths(2, 4);
    std::cout << "All paths from 2 to 4 -\n ";
    for (const auto& path : allPaths) {
        for (int vertex : path) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }*/
	
//	std::cout << (g1.hasCycle()   ? "True\n" : "False \n");
/*	g1.hasCycleInDirectedGraph() ? std::cout << "Graph contains cycle\n" : std::cout << "Graph doesn't contain cycle\n";
	Graph g2(3);
	g2.addEdge(0,1);
	g2.addEdge(0,2);
	g2.hasCycle() ? std::cout << "Graph contains cycle\n" : std::cout << "Graph doesn't contain cycle\n";
	*/
	/*	std::cout << "Transpose Matrix  \n";
	g.transpose();
	g.print();
*/
/*	std::cout << "Result of dfsRec startig from vertex of 1- ";
	g.dfs(1);
	std::cout << std::endl;

	std::cout << "Result of dfsItrative startig from vertex of 1- ";
	g.dfsIterative(1);
	std::cout << std::endl;

	std::cout << "Result of BfsRec starting from 1 vertex - ";
	g.bfs(1);
	std::cout << std::endl;
	
	std::cout << "Result of BfsIter starting vertex of 1 - ";
	g.bfsIter(1);
	std::cout << std::endl;
	
    Graph myGraph(6);

    myGraph.addEdge(0, 1);
    myGraph.addEdge(1, 2);
    myGraph.addEdge(2, 0);

    myGraph.addEdge(3, 4);
    myGraph.addEdge(4, 5);

	std::cout << "Get component count - " << myGraph.getComponentCount();
	std::cout << std::endl;

	int shortDist = myGraph.shortestDistenceInUnweightedGraph(0,2);
	if(shortDist == -1) {
		std::cout << "There is not distance betweeen two vertex\n";
	} else {
		std::cout << "SHort dist is - " << shortDist << std::endl;
	}
	
	*/
	
}

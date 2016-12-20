#pragma once
#include <iostream>
#include <vector>

using std::cin;
using std::cout;
using std::endl;
using std::vector;

typedef int DataType;
const int SUSCEPTIBLE = 0;
const int INFECTED = 1;

class Node {
public:
	int state;
	DataType data;
	vector<int> edges;
	Node(int s = SUSCEPTIBLE, DataType d = 0) :state(s), data(d) {}  // default, every node is susceptible, data is 0

	bool operator== (Node &y) {
		if (data == y.data && edges == y.edges)
			return true;
		else
			return false;
	}

	bool operator!=(Node &y) {
		!(*this == y);
	}
};

class Graph {
public:
	int vexNum;
	int edgeNum;
	vector<Node> Vertex;

	Graph() {
		vexNum = 0;
		edgeNum = 0;
	}

	Graph(vector<Node> &v);

	int GetVexNum();
	int GetEdgeNum();

	void InsertVertex(Node &v);

	bool AddEdge(int x, int y);

	~Graph() {};
};


Graph::Graph(vector<Node>& v)
{
	Vertex.assign(v.begin(), v.end());
	vexNum = Vertex.size();
	edgeNum = GetEdgeNum();
}

int Graph::GetVexNum()
{
	return Vertex.size();
}

int Graph::GetEdgeNum()
{
	edgeNum = 0;
	for (auto iter = Vertex.begin(); iter != Vertex.end(); ++iter) {
		edgeNum += (iter->edges).size();
	}
	return edgeNum;
}

void Graph::InsertVertex(Node & v)
{
	Vertex.push_back(v);
	++vexNum;
}

// add the edge x->y to the graph, no multiple-dege allowed. FOR UNDIRECTED GRAPH, use AddEdge(x,y) and AddEdge(y,x)
bool Graph::AddEdge(int x, int y)
{
	if ((x >= Vertex.size() || x < 0) || (y >= Vertex.size()) || y < 0) {
		cout << "AddEdge problem : the vertex id is out of range !" << endl;
		return false;
	}
	int size = Vertex[x].edges.size();
	for (auto iter = Vertex[x].edges.begin(); iter != Vertex[x].edges.end(); ++iter) {
		if (Vertex[*iter] == Vertex[y]) {
			//cout << "The edge " << x << " ---> " << y << "is already existed !" << endl;
			return false;
		}
	}
	Vertex[x].edges.push_back(y);
	++edgeNum;
	return true;
}

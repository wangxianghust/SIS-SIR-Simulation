#include "Graph.h"
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <map>
#include <utility>

using std::ifstream;
using std::ofstream;
using std::istream;
using std::ostream;
using std::string;
using std::map;
using std::pair;

// customed split function for string process
vector<string> split(const string &s, const string &delim) {
	vector<string> res;
	string::size_type front = 0;
	string::size_type last = s.find_first_of(delim, front);
	while (last != string::npos) {
		if (last > front) {
			string tmp = s.substr(front, last - front);
			res.push_back(tmp);
		}
		front = last + 1;
		last = s.find_first_of(delim, front);
	}
	if (last > front) {
		res.push_back(s.substr(front, last - front));
	}
	return res;
}

/*
Input file required:
from(int) \t	to(int)
from..	  \t    to..
*/
auto formatData(string filePath) {
	vector<int> edges;
	//edges for saved the edges information, 
	//edge0: 0,1 edge1:2,3 etc..the adjacent element forms a edge
	ifstream input(filePath);
	if (input) {
		string line;
		while (getline(input, line)) {
			auto from = split(line, "\t")[0];
			edges.push_back(stoi(from));
			auto to = split(line, "\t")[1];
			edges.push_back(stoi(to));
		}
	}
	return edges;
}

Graph constructGraph(vector<int> edges) {
	Graph G;
	auto max = max_element(edges.begin(), edges.end());
	auto vertexNum = *max + 1;
	vector<Node> v(vertexNum);
	// add graph data
	for (int i = 0; i < vertexNum; ++i) {
		v[i].data = i;
		G.InsertVertex(v[i]);
	}
	for (int i = 0, j = 1; j < edges.size(); i += 2, j += 2) {
		G.AddEdge(edges[i], edges[j]);
		G.AddEdge(edges[j], edges[i]);  // for undirected-graph
	}
	return G;
}

void testInfo(const string& info) {
	cout << "-----" << info << "-----" << endl;
}

template<typename T>
void print(const vector<T> &v) {
	for (auto i : v) {
		cout << i << " - ";
	}
	cout << endl;
}

/*/////////////////////////////
The functions for the SIS model. 
*//////////////////////////////


// input: the total vertex number, totalNum,,,the infected node, idex
// output: a vector which only the correspond node is infected.
vector<int> initialState(int totalNum, int index) {
	vector<int> initial(totalNum, 0);
	initial[index] = 1;
	return initial;
}


// calculate for the probability, is it be infected or not.
// the ratio is the probability, and the accuracy which means, eg, 1000 denots the accuracy is 0.001
bool guessTrue(double ratio) {
	int accuracy = 1000;
	int guess = rand() % accuracy;
	double guess_format = guess*1.0 / accuracy;
	//cout << guess_format << "  " << ratio << " ";
	if (guess_format <= ratio) {
		//cout << "true" << endl;
		return true;
	}
	else {
		//cout << "false" << endl;
		return false;
	}
}

// for recovery process, but we can not change the state of this step, instead, we just save the temp, and update it 
// at suitable time.
// input: is I or S, state,,, recovery_ratio
vector<int> recoveryProcess(double recovery_ratio, Graph& graph) {
	vector<int> ret;
	for (int index = 0; index < graph.GetVexNum(); ++index) {
		int node_state = graph.Vertex[index].state;
		if (node_state) {     // we only choose the infected node to choose recovery or not.
			if (guessTrue(recovery_ratio)) {
				node_state = 0;
			} 
		}
		ret.push_back(node_state);
	}
	return ret;
}

// we only consider the S state node, see if it will be infected by his neighbors; 
// the processing of the neighbor infecte is dependentely.
// input :  is I or S? state,,, infected_ratio,,, and graph
// output : the new state
vector<int> infectedProcess(vector<int> &state, double infected_ratio, const Graph& graph) {  
	vector<int> ret;
	for (int index = 0; index < state.size(); ++index) {
		int node_state = state[index];
		if (!node_state) {
			// get this node all neighbors and then guess the changing state
			auto neighbors = graph.Vertex[index].edges;
			for (auto neighbor : neighbors) {   // the index of the edges
				if (state[neighbor] == 1) {      // if the neighbor is infected
					if (guessTrue(infected_ratio)) node_state = 1;  // guessTrue, then infect it.
				}
			}
		}
		ret.push_back(node_state);
	}
	return ret;
}

// update the graph nodes state by new_state
// Input: new_s,,, infected_ratio,,, and graph
void update(vector<int> &new_state, Graph& graph) {
	for (int index = 0; index < new_state.size(); ++index) {
		graph.Vertex[index].state = new_state[index];
	}
}


// run SIS model
// Input
//
vector<double> runSIS(double infect_ratio, double recovery_ratio, int iterator_num, int step_num, Graph& graph) {
	int vertexNum = graph.GetVexNum();
	int edgesNum = graph.GetEdgeNum();

	vector<double> ret(vertexNum, 0); // for save average score of all iterations. 
	vector<vector<double>> scores; // scores are all average score, score is denoted by the infected numbers.

	for (int i = 0; i < iterator_num; ++i) {

		cout << endl;
		cout << "-----This the " << i << " iterator-----" << endl;
		cout << endl;

		vector<double> score; // save every node's score in the vector, every iteration initial these to 0.
							  // now we calculate the one iterator,and get the every node's infected
		for (int index = 0; index < vertexNum; ++index) {
			int infected_index = index; // from the first node to be infected
			auto initial_state = initialState(vertexNum, infected_index);
			testInfo("Initial state");
			print(initial_state);
			update(initial_state, graph);

			vector<int> recovery_state(vertexNum, 0);
			vector<int> infected_state(vertexNum, 0);
			int infectedNum = 0; // denote the index_th node have infected numbers.
			for (int j = 0; j < step_num; ++j) {
				cout << " STEP " << j << " " << endl;
				testInfo("recovery processing");
				recovery_state = recoveryProcess(recovery_ratio, graph);
				print(recovery_state);
				update(recovery_state, graph);

				testInfo("infected processing");
				infected_state = infectedProcess(recovery_state, infect_ratio, graph);
				print(infected_state);
				update(infected_state, graph);

				// static the infected numbers
				infectedNum = count(infected_state.begin(), infected_state.end(), 1);
				cout << "infected number is " << infectedNum << endl;
			}
			score.push_back(infectedNum);
		}
		for (auto i : score) cout << i << endl;
		scores.push_back(score);   // every iteration result is saved in the scores.
	}

	for (int iter = 0; iter < iterator_num; ++iter){
		for (int index = 0; index < vertexNum; ++index) {
			auto iter_socre = scores[iter][index];
			ret[index] += iter_socre;
			cout << iter_socre << " ";
		}
		cout << endl;
	}

	for (auto &s : ret) {   // get the average score
		s /= iterator_num;
	}
	 
	return ret;
}

////////////// for rank the scores with correspond index //////////////////
/*//////////////////////////////////////////
//A function to rank the vector number.
*///////////////////////////////////////////
// a compare function for sort, then we can use it the sort the map by value 
struct CmpByValue {
	bool operator()(const pair<int, double> & lhs, const pair<int, double> & rhs)
	{
		return lhs.second > rhs.second;
	}
};

void importanceRank(vector<double>& index_value, vector<pair<int, double>> &sorted_index) {
	map<int, double> tmp;
	int nodeID = 0;
	for (auto value : index_value) {
		tmp.insert({ nodeID++, value });
		//tmp.insert({ ++nodeID, value });
	}
	for (auto i : tmp) {
		sorted_index.push_back(i);
	}
	sort(sorted_index.begin(), sorted_index.end(), CmpByValue());
}

void formatOutput(vector<pair<int, double>>& sorted_index, const string& filePath) {
	ofstream output(filePath);
	for (auto ele : sorted_index) {
		cout << ele.first << "\t" << ele.second << endl;
		output << ele.first << "\t" << ele.second << endl;
	}
}

int main() {
	// construct the graph from the edge datas.
	string filePath("data/test.txt");
	auto data = formatData(filePath);
	Graph graph = constructGraph(data);

	// some configunation of the SIS, and we can modify it later.
	double infect_ratio = 0.6;
	double recovery_ratio = 0.1;
	int iterator_num = 50; 
	int step_num = 5;

	auto result = runSIS(infect_ratio, recovery_ratio, iterator_num, step_num, graph);
	print(result);

	// rank the result and save to result folder
	vector<pair<int, double>> rank;
	importanceRank(result, rank);
	formatOutput(rank, "./result/sis_score.txt");
	
	system("PAUSE");
	return 0;
}
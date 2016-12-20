#include "Graph.h"
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <map>
#include <utility>
#include <random>
#include <cmath>
#include <stdlib.h>

using std::ifstream;
using std::ofstream;
using std::istream;
using std::ostream;
using std::string;
using std::map;
using std::pair;
using std::default_random_engine;
using std::uniform_real_distribution;

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
vector<int> formatData(string filePath) {
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

void formatOutput(vector<pair<int, double>>& sorted_index, const string& filePath) {
	ofstream output(filePath);
	for (auto ele : sorted_index) {
		cout << ele.first << "\t" << ele.second << endl;
		output << ele.first << "\t" << ele.second << endl;
	}
}

template<typename T>
void formatOutput(vector<T>& data, const string& filePath){
    ofstream output(filePath);
    int size = data.size();
    for(int index = 0;  index < size; ++index){
        output << index << "\t" << data[index] << endl;
    }
}

template<typename T, typename N>
void formatOutput(vector<T>& data, vector<N>& other, const string& filePath) {
	ofstream output(filePath);
	int size = data.size();
	for (int index = 0; index < size; ++index) {
		output << index << "\t" << data[index] << "\t" << other[index] << endl;
	}
}

/*/////////////////////////////
The functions for the SIR model.
*//////////////////////////////


// input: the total vertex number, totalNum,,,the infected node, idex
// output: a vector which only the correspond node is infected.
vector<int> initialState(int totalNum, int index) {
	vector<int> initial(totalNum, 0);
	initial[index] = 1;
	return initial;
}

// input: graph
// output: the newest state from the graph
vector<int> getStateFromGraph(Graph &graph) {
	vector<int> ret;
	for (int index = 0; index < graph.Vertex.size(); ++index) {
		ret.push_back(graph.Vertex[index].state);
	}
	return ret;
}


// calculate for the probability, is it be infected or not.
// the ratio is the probability, and the accuracy which means, eg, 1000 denots the accuracy is 0.001
//bool guessTrue(double ratio) {
//	int accuracy = 1000;
//	int guess = rand() % accuracy;
//	double guess_format = guess*1.0 / accuracy;
//	//cout << guess_format << "  " << ratio << " ";
//	if (guess_format <= ratio) {
//		//cout << "true" << endl;
//		return true;
//	}
//	else {
//		//cout << "false" << endl;
//		return false;
//	}
//}

bool guessTrue(double ratio){
    static default_random_engine e;
    static uniform_real_distribution<double> u(0,1);
    double guess = u(e);
    if(guess <= ratio){
        return true;
    } else {
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
		if (node_state == 1) {     // we only choose the infected node to choose recovery or not.
			if (guessTrue(recovery_ratio)) {
				node_state = 2;
			}
		}
		ret.push_back(node_state);
	}
	return ret;
}

// we only consider the S state node, see if it will be infected by his neighbors; 
// the processing of the neighbor infecte is dependentely.
// input :  is I or S? state, we can get it from the graph; infected_ratio; graph
// output : the new state
vector<int> infectedProcess(double infected_ratio, const Graph& graph) {
	vector<int> ret;
	int size = graph.Vertex.size();
	for (int index = 0; index < size; ++index) {
		int node_state = graph.Vertex[index].state;
		if (node_state == 0) {   // we only get the S state node
			// get this node all neighbors and then guess the changing state
			auto neighbors = graph.Vertex[index].edges;
			//for (auto neighbor : neighbors) {   // the index of the edges
			//	if (graph.Vertex[neighbor].state == 1) {      // if the neighbor is infected
			//		if (guessTrue(infected_ratio)) node_state = 1;  // guessTrue, then infect it.
			//	}
			//}
            int infected_neighbor_num = count(neighbors.begin(), neighbors.end(), 1);
            double infected_probability = 1 - pow(1 - infected_ratio, infected_neighbor_num);
            if(guessTrue(infected_probability)) node_state = 1;
		}
		ret.push_back(node_state);
	}
	return ret;
}

// update the graph nodes state by new_state, just assign the new state to the old state.
// Input: new_s,,, infected_ratio,,, and graph
void update(vector<int> &new_state, Graph& graph) {
	for (int index = 0; index < new_state.size(); ++index) {
		graph.Vertex[index].state = new_state[index];
	}
}

// In this SIR, in each step, each infected node will infecte each of its neighbors 
// Then each infected start the recovery state.
// Input: initial_state; infected_state:the first action of each step of SIR; recovery_state:the second action of the step of SIR
// Then update the state, write it to the Graph.
void update(vector<int> &initial_state, vector<int> &infected_state, vector<int> &recovery_state, Graph &graph) {
	for (int index = 0; index < initial_state.size(); ++index) {
		if (initial_state[index] == 0) {	// if the initial is S, we see is it be infected
			graph.Vertex[index].state = infected_state[index];
		}
		else if (initial_state[index] == 1) {
			graph.Vertex[index].state = recovery_state[index];	// if the initial is I, we see is it recovery
		}
		else {
			graph.Vertex[index].state = 2;		// if the initial state is R, then we it will be R all the time.
		}
	}
}



// test is the end condition satisfied? if there is any node is infected, then go on.
// Input : node_state
// Output : true or false
bool isEnd(vector<int> &state) {
	for (auto s : state) {
		if (s == 1) {
			return false;
		}
	}
	return true;
}


// run SIR model
// Input
//
vector<double> runSIR(double infect_ratio, double recovery_ratio, int iterator_num, Graph& graph) {
	int vertexNum = graph.GetVexNum();
	int edgesNum = graph.GetEdgeNum();

	vector<double> ret(vertexNum, 0); // for save average score of all iterations. 
	vector<vector<double>> scores; // scores are all average score, score is denoted by the infected numbers.

	for (int iter = 0; iter < iterator_num; ++iter) {

		cout << endl;
		cout << "#iterator " << iter << endl;
		cout << endl;

		vector<double> score; // save every node's score in the vector, every iteration initial these to 0.
							  // now we calculate the one iterator,and get the every node's infected
		vector<int> infected_step(vertexNum, 0);   // calculate the step of each infected process
		for (int index = 0; index < vertexNum; ++index) {
			int infected_index = index; // from the first node to be infected
			auto initial_state = initialState(vertexNum, infected_index);
			//print(initial_state);
			update(initial_state, graph); // assign the initial state to the graph

			vector<int> recovery_state(vertexNum, 0);
			vector<int> infected_state(vertexNum, 0);
			int infectedNum = 0; // denote the index_th node have infected numbers.
			int recoveryNum = 0;
			int step = -1;
			do {
				++step;
				//cout << " STEP " << ++step << " " << endl;
				//testInfo("Initial state");
				initial_state = getStateFromGraph(graph);
				//print(initial_state);

				//testInfo("infected");
				infected_state = infectedProcess(infect_ratio, graph);
				//print(infected_state);

				//testInfo("recovery");
				recovery_state = recoveryProcess(recovery_ratio, graph);
				//print(recovery_state);

				update(initial_state, infected_state, recovery_state, graph);

				// static the infected numbers
				infectedNum = count(infected_state.begin(), infected_state.end(), 1); // 注意，这里有一个问题，每次会多迭代一步，但是对最后结果没有影响，下个版本修改 
				//cout << "infected number is " << infectedNum << endl;
			} while (infectedNum != 0);
			infected_step.push_back(step);
			recoveryNum = count(recovery_state.begin(), recovery_state.end(), 2);
			score.push_back(recoveryNum);
			int degree = graph.Vertex[index].edges.size();
			cout << "#R " << recoveryNum << "   #degree " << degree << "   #step " << step << endl;
		}
        formatOutput(score,infected_step, "./results/score_" + std::to_string(iter));
		for (auto i : score) cout << i << " ";
		scores.push_back(score);   // every iteration result is saved in the scores.
	}

	for (int iter = 0; iter < iterator_num; ++iter) {
		for (int index = 0; index < vertexNum; ++index) {
			auto iter_socre = scores[iter][index];
			ret[index] += iter_socre;
			//cout << iter_socre << " ";
		}
		//cout << endl;
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


int main() {
    cout << "Program start..." << endl;
	// construct the graph from the edges datas.
	//string filePath("data/test.txt");
	string filePath("data/Email_Enron.txt");
	auto data = formatData(filePath);
	Graph graph = constructGraph(data);

	// some configunation of the SIS, and we can modify it later.
	double infect_ratio = 0.0105;
	double recovery_ratio = 1;
	int iterator_num = 50;

	auto result = runSIR(infect_ratio, recovery_ratio, iterator_num, graph);
	print(result);

	// rank the result and save to result folder
	vector<pair<int, double>> rank;
	importanceRank(result, rank);
	formatOutput(rank, "./results/sir_score.txt");

	cout << "done !" << endl;

	system("PAUSE");
	return 0;
}

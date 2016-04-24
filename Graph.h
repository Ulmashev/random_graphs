#ifndef GRAPH_H
#define GRAPH_H

#include <stack>
#include <queue>
#include <iostream>
#include <time.h>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Set_of_degrees.h"
#include "walker.h"

using namespace std;

#define maxn 500
#define e 2.7182818284590452
#define eps 10e-7
#define d 10000
#define INF 10e7
#define TP 10

struct Graph { //structure for the graph
	double maxdeg, tau_in, tau_out, kappa; // parametrs of the distribution
	Set_of_degrees ind, outd, sumd; //oblects for keeping in-degrees, out-degrees and common degrees
	int number_of_edges;
	int delta;

	int indeg[maxn], outdeg[maxn], sumdeg[maxn]; //arrays for in-degrees, out-degrees and common degrees (temporary arrays, data from them will be written in corresponding structures)

	vector <int> g[maxn]; //adjacency list for the graph
	bool adj[maxn][maxn];

	Graph(int MD, double T_in, double T_out, double K) {
		maxdeg = MD;
		tau_in = T_in;
		tau_out = T_out;
		kappa = K;
		double sum = 0;
		number_of_edges = 0;

		srand(time(NULL));
		this->simulate(); //modeling graph until sums of in-degress and out-degress are equal
		
	//	return;

		ind = Set_of_degrees();
		outd = Set_of_degrees();
		sumd = Set_of_degrees();


		for (int i = 0; i < n; i++) { // transfer data from arrays to structures
			outd.values[i] = outdeg[i];
			ind.values[i] = indeg[i];
			sumd.values[i] = outdeg[i] + indeg[i];
			number_of_edges += indeg[i];
			for (int j = 0; j < n; j++)
				adj[i][j] = 0;

		}

		int edges_left = number_of_edges; // the variable pointing on the number of unused edges 

		while (edges_left != 0) {
			vector <pair <int, int> > available;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					if (i != j && !adj[i][j] && outdeg[i] > 0 && indeg[j] > 0) {
						available.push_back(make_pair(i, j));
					}
			if (available.size() == 0) {
				break;
			}
			srand(time(NULL));
			int edge = rand() % available.size();
			int out = available[edge].first;
			int in = available[edge].second;
			adj[out][in] = 1;
			outdeg[out]--;
			indeg[in]--;
			g[out].push_back(in);
			edges_left--;
		}
		cout << "Edges left: " << edges_left << endl;
		if (edges_left != 0) {
			for (int i = 0; i < n; i++) { // recounting if there are edges left
				ind.values[i] = 0;
				outd.values[i] = 0;
			}
			for (int i = 0; i < n; i++) {
				outd.values[i] = g[i].size();
				for (int j = 0; j < g[i].size(); j++)
					ind.values[g[i][j]]++;
			}
		}
		/*while (edges_left != 0) {  // matching nodes until there are any unused edges 
			vector <pair<int, int> > v_in; //vector for unused nodes and their in-degrees (for quick match)
			vector <pair<int, int> > v_out;
			for (int i = 0; i < n; i++) { //filling vectors
				if (indeg[i] != 0) v_in.push_back(make_pair(i, indeg[i]));
				if (outdeg[i] != 0) v_out.push_back(make_pair(i, outdeg[i]));
			}
			srand(time(NULL));
			int in = -1;
			int out = -1;
			if (v_in.size() == 1 && v_out.size() == 1 && v_in[0].first == v_out[0].first) { //the case when we can add only one edge but it is a self-loop
				in = v_in[0].first;
				out = v_out[0].first;
				indeg[in]--;
				outdeg[out]--;
				edges_left--;
			}
			else { // else we need to choose a couple of nodes from v_in and v_out
				while (in == out) {
					in = v_in[rand() % v_in.size()].first; // taking an in-node randomly
					out = v_out[rand() % v_out.size()].first; // taking an out-node randomly
				}
				g[out].push_back(in); 
				indeg[in]--;
				outdeg[out]--;
				edges_left--;

			}
		}
		for (int i = 0; i < n; i++) { //if there were some self-loops in the end, then we need to recount degrees
			ind.values[i] = 0;
			outd.values[i] = 0;
		}
		for (int i = 0; i < n; i++) {
			outd.values[i] = g[i].size();
			for (int j = 0; j < g[i].size(); j++)
				ind.values[g[i][j]]++;
		}*/


	}

	bool existing_edge(int out, int in) {
		for (int i = 0; i < g[out].size(); i++)
			if (g[out][i] == in) return true;
		return false;
	}

	void show_edges() {
		cout << "List of edges:" << endl;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < g[i].size(); j++)
				cout << i + 1 << " " << g[i][j] + 1 << endl;
	}

	void statistics() { // method for statistics outputing  
		cout << "|V| = " << n << endl;
		cout << "|E| = " << this->number_of_edges << endl;

		ind.findmax(); ind.findmin(); ind.findaver(); ind.findmed();
		outd.findmax(); outd.findmin(); outd.findaver(); outd.findmed();
		sumd.findmax(); sumd.findmin(); sumd.findaver(); sumd.findmed();

		cout << "Max: " << endl << "In: " << ind.mx << " Out: " << outd.mx << " Sum: " << sumd.mx << endl;
		cout << "Min: " << endl << "In: " << ind.mn << " Out: " << outd.mn << " Sum: " << sumd.mn << endl;
		cout << "Average: " << endl << "In: " << ind.average << " Out: " << outd.average << " Sum: " << sumd.average << endl;
		cout << "Median: " << endl << "In: " << ind.med << " Out: " << outd.med << " Sum: " << sumd.med << endl;

		cout << endl;

		cout << "In ";
		ind.top(TP);
		cout << "Out ";
		outd.top(TP);
		cout << "Sum ";
		sumd.top(TP);

		int node_start = 1;
		int node_end = 2;
		pair <int, int> pair_ind = ind.find_vert_from_top(TP, node_start, node_end);
		pair <int, int> pair_outd = outd.find_vert_from_top(TP, node_start, node_end);
		pair <int, int> pair_sumd = sumd.find_vert_from_top(TP, node_start, node_end);

		cout << "Shortest path between nodes with top in-degrees (number " << node_start << " and " << node_end << " in rate list) is: ";
		cout << find_path(pair_ind.first, pair_ind.second) << endl;
		cout << "Shortest path between nodes with top out-degrees (number " << node_start << " and " << node_end << " in rate list) is: ";
		cout << find_path(pair_outd.first, pair_ind.second) << endl;
		cout << "Shortest path between nodes with top sum-degrees (number " << node_start << " and " << node_end << " in rate list) is: ";
		cout << find_path(pair_sumd.first, pair_ind.second) << endl;


	}

	int find_path(int a, int b) { // Dejkstra
		int ans;

		int dist[maxn], prev[maxn];
		bool u[maxn];
		for (int i = 0; i < n; i++) {
			dist[i] = INF;
			prev[i] = -1;
			u[i] = false;
		}
		dist[a] = 0;
		for (int i = 0; i < n; i++) {
			int v = -1;
			for (int j = 0; j < n; j++)
				if (!u[j] && (v == -1 || dist[j] < dist[v]))
					v = j;
			if (dist[v] == INF)
				break;
			u[v] = true;

			for (int j = 0; j < g[v].size(); j++) {
				int end = g[v][j];
				if (dist[v] + 1 < dist[end]) {
					dist[end] = dist[v] + 1;
					prev[end] = v;
				}
			}
		}

		return dist[b];
	}
	void simulate() { // method for modeling pair of degrees 
		int check_in = 0;
		int check_out = 0;

		for (int i = 0; i < n; i++) {
			indeg[i] = 0;
			outdeg[i] = 0;
		}

		Random_variable in = Random_variable(maxdeg, tau_in, kappa);
		Random_variable out = Random_variable(maxdeg, tau_out, kappa);


		for (int i = 0; i < n; i++) {
			indeg[i] = in.walker();
			outdeg[i] = out.walker();
			check_in += indeg[i];
			check_out += outdeg[i];
		}

		if (check_in != check_out) {
			srand(time(NULL));
			int node = rand() % n;
			check_in -= indeg[node];
			check_out -= outdeg[node];
			if (check_in == check_out) {
				indeg[node] = 1;
				outdeg[node] = 1;
			}
			else if (check_in > check_out) {
				outdeg[node] = check_in - check_out + 1;
				indeg[node] = 1;
			}
			else {
				outdeg[node] = 1;
				indeg[node] = check_out - check_in + 1;
			}



		}
	}
	

};

#endif 
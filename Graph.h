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
#include "edge.h"

using namespace std;

#define maxn 500
#define e 2.7182818284590452
#define eps 10e-7
#define d 10000
#define INF 10e7
#define TP 10
#define gamma 0.04

struct Graph { //structure for the graph
	double maxdeg, tau_in, tau_out, kappa; // parametrs of the distribution
	Set_of_degrees ind, outd, sumd; //oblects for keeping in-degrees, out-degrees and common degrees
	int number_of_edges;
	int delta;

	int indeg[maxn], outdeg[maxn], sumdeg[maxn]; //arrays for in-degrees, out-degrees and common degrees (temporary arrays, data from them will be written in corresponding structures)
	double l[maxn], a_ib[maxn], a_m[maxn], dep[maxn];

	vector <Edge> g[maxn]; //adjacency list for the graph
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


		for (int i = 0; i < n; i++) { // transfer data from arrays to structures and initializating
			outd.values[i] = outdeg[i];
			ind.values[i] = indeg[i];
			sumd.values[i] = outdeg[i] + indeg[i];
			number_of_edges += indeg[i];
			for (int j = 0; j < n; j++)
				adj[i][j] = 0;
			a_m[i] = 0.8;
			a_ib[i] = 0;
			l[i] = 0;
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
			g[out].push_back(Edge(in, 0));
			edges_left--;
		}
		cout << "Edges left: " << edges_left << endl;
		if (edges_left != 0) { // recounting degrees if there are edges left
			for (int i = 0; i < n; i++) { 
				ind.values[i] = 0;
				outd.values[i] = 0;
			}
			for (int i = 0; i < n; i++) {
				outd.values[i] = g[i].size();
				for (int j = 0; j < g[i].size(); j++)
					ind.values[g[i][j].u]++;
			}
		}
		for (int i = 0; i < n; i++) { //counting weights
			for (int j = 0; j < g[i].size(); j++)
				g[i][j].weight = 0.2 / ind.values[g[i][j].u];
		}
		for (int i = 0; i < n; i++) { //counting financial parametrs
			for (int j = 0; j < g[i].size(); j++) {
				a_ib[g[i][j].u] += g[i][j].weight;
				l[i] += g[i][j].weight;
			}
		}
		for (int i = 0; i < n; i++) {
			dep[i] = 1 - l[i] - gamma;
		}

	}

	bool existing_edge(int out, int in) {
		for (int i = 0; i < g[out].size(); i++)
			if (g[out][i].u == in) return true;
		return false;
	}

	void show_edges() {
		cout << "List of edges:" << endl;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < g[i].size(); j++)
				cout << i + 1 << " " << g[i][j].u + 1 << endl;
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
				int end = g[v][j].u;
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
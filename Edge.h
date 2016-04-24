#ifndef EDGE_H
#define EDGE_H

using namespace std;

struct Edge { // struct for edge in the adjacency list
	int u;
	double weight;
	Edge(int U, double W) {
		u = U;
		weight = W;
	}
};

#endif 
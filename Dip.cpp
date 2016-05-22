#define _CRT_SECURE_NO_WARNINGS

#include <stack>
#include <queue>
#include <iostream>
#include <time.h>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Graph.h"

using namespace std;

#define maxn 500
#define e 2.7182818284590452
#define eps 10e-7
#define d 10000
#define INF 10e7
#define TP 3
#define Iter 10

vector <int> size;
vector <int> edges_left;
vector <int> indeg;
vector <int> outdeg;
vector <int> answers_edges;
vector <int> answers_faied;
vector <double> difference_in;
vector <double> difference_out;
vector <double> answers_dif_in;
vector <double> answers_dif_out;

int n;


double average(vector <int> v) { // counting average of <int> vector
	double sum = 0;
	for (int i = 0; i < v.size(); i++)
		sum += v[i];
	return sum / v.size();
}

double average_double(vector <double > v) { // counting average if <double> vector
	double sum = 0;
	for (int i = 0; i < v.size(); i++)
		sum += v[i];
	return sum / v.size();
}

int main() {
	freopen("dip.out", "w", stdout);

	n = 100;

	for (int i = 0; i < Iter; i++) {
		vector <int> a;
		vector <int> b;
		for (int j = 0; j < 30; j++) {
			Graph gr = Graph(n, n - 1, 1.92, 2.64, 2.5);
		//	cout << "n = " << n << " Expectation in = " << gr.z_in << " Expectation out = " << gr.z_out << " Average = " << gr.average_half_deg << endl;
			difference_in.push_back(abs(gr.z_in - gr.average_half_deg));
			difference_out.push_back(abs(gr.z_out - gr.average_half_deg));
		/*	srand(time(NULL)); // simulating defaults
			gr.fail(rand() % (n));
			a.push_back(gr.number_of_edges);
			b.push_back(gr.defaulted_banks);
			cout << gr.defaulted_banks << endl; */
		}
		n += 5;
		answers_dif_in.push_back(average_double(difference_in));
		answers_dif_out.push_back(average_double(difference_out));
		cout << answers_dif_in[i] << " " << answers_dif_out[i] << endl;
	/*	answers_edges.push_back(average(a));
		answers_faied.push_back(average(b));
		cout << "n = " << n << " " << answers_edges[i] << " " << answers_faied[i] << endl;
		cout << gr.number_of_edges << " " << gr.ind.average << endl;
		size.push_back(gr.number_of_edges);
		edges_left.push_back(gr.edges_left);
		indeg.push_back(gr.ind.average);
		outdeg.push_back(gr.outd.average); */
	}
	
	return 0;
}
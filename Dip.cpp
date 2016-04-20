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
//#include "Set_of_degrees.h"
#include "Graph.h"

using namespace std;

#define maxn 1000
#define e 2.7182818284590452
#define eps 10e-7
#define d 10000
#define INF 10e7
#define TP 10
#define n 30

int main() {

	freopen("dip.out", "w", stdout);

	Graph gr = Graph(7, 0.52, 2.5);
	gr.show_edges();
	gr.statistics();
	cout << endl;
	cout << "-----------------------------------------------------------";
	cout << endl; 

	return 0;
}
#ifndef WALKER_H
#define WALKER_H

#include <stack>
#include <queue>
#include <iostream>
#include <time.h>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

#define maxn 1000
#define e 2.7182818284590452
#define eps 10e-7
#define d 10000
#define INF 10e7
#define TP 10
#define n 30

struct Random_variable {
	double uniform, tau, kappa, maxdeg;
	double threshold[maxn]; //array for thresholds
	double p[maxn]; //array for current probabilities
	double sum, c;

	stack <pair <double, int> > donor;
	stack <pair <double, int> > recepient;

	vector <pair <double, int> > v; //final vector of thresholds


	Random_variable(int MD, double T, double K) {
		maxdeg = MD;
		tau = T;
		kappa = K;
		uniform = 1 / maxdeg;
		sum = 0;
		for (int i = 0; i < maxn; i++) { //initialization
			p[i] = 0;
			threshold[i] = 0;
		//	indeg[i] = 0;
		//	outdeg[i] = 0;
		//	sumdeg[i] = 0;
		}
		for (int i = 1; i <= maxdeg; i++) { //counting of non-normalized probabilities
			p[i] = 1 / (pow((double)i, tau) * pow(e, (double)i / kappa));
			sum += p[i];
		}
		c = 1 / sum;
		for (int i = 1; i <= maxdeg; i++) { //normalization
			p[i] = c / (pow((double)i, tau));// * pow(e, (double)i / kappa));
			if (p[i] > uniform && abs(p[i] - uniform) > eps) //including recepients in the corresponding stack
				recepient.push(make_pair(p[i], i));
			else if (p[i] < uniform && abs(p[i] - uniform) > eps) //including donors in the corresponding stack 
				donor.push(make_pair(p[i], i));
			threshold[i] = uniform * (i - 1); // setting initial current thresholds on values corresponding to the uniform distribution
			if (p[i] < uniform)
				threshold[i] += p[i];

			v.push_back(make_pair(uniform * (i - 1), i));
		}
		while (donor.empty() == false && recepient.empty() == false) { // transfer excess probabilities from donors to recepients 
			double a = recepient.top().first;
			double b = donor.top().first;
			if (a - uniform >= uniform - b) {
				recepient.top().first -= uniform - b;
				v.push_back(make_pair(threshold[donor.top().second], recepient.top().second));
				threshold[donor.top().second] += uniform - b;
				donor.pop();
				if (a - uniform == uniform - b)
					recepient.pop();
			}
			else if (a - uniform < uniform - b) {
				donor.top().first += a - uniform;
				v.push_back(make_pair(threshold[donor.top().second], recepient.top().second));
				threshold[donor.top().second] += a - uniform;
				recepient.pop();
			}
		}

		sort(v.begin(), v.end()); //sorting of final thresholds 
	}
	int walker() { // Walker
		int deg = 0;
		double x = rand() % (int)d;
		x = (double)x / d;
		for (int i = 0; i < v.size(); i++) {
			if (i == v.size() - 1 || v[i].first <= x && v[i + 1].first > x) {
				return v[i].second;
			}
		}
	}

};


#endif

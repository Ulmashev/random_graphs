#ifndef SET_OF_DEGREES_H
#define SET_OF_DEGREES_H

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

struct Set_of_degrees {    //struct for arrays of degrees and their statistical characteristics  
	double values[maxn];
	double mx, mn, average, med;
	Set_of_degrees() {
		for (int i = 0; i < n; i++)
			values[i] = 0;
		mx = -1;
		mn = maxn + 1;

	}
	void top(int k) { //output of k best elements
		cout << "Top " << k << endl;
		vector <int> v;
		for (int i = 0; i < n; i++)
			v.push_back(-values[i]);
		sort(v.begin(), v.end());
		for (int i = 0; i < k; i++)
			cout << -v[i] << " ";
		cout << endl;
	}
	pair <int, int> find_vert_from_top(int k, int a, int b) { //search of a-th and b-th element in the top rank
		vector <pair <int, int> > v;
		for (int i = 0; i < n; i++)
			v.push_back(make_pair(-values[i], i));
		sort(v.begin(), v.end());
		pair <int, int> ans = make_pair(v[a - 1].second, v[b - 1].second);
		return ans;
	}
	void findmax() {
		for (int i = 0; i < n; i++)
			if (values[i] > mx)
				mx = values[i];
	}
	void findmin() {
		for (int i = 0; i < n; i++)
			if (values[i] < mn)
				mn = values[i];
	}
	void findaver() { //search of the average element in the array
		double sum = 0;
		for (int i = 0; i < n; i++)
			sum += values[i];
		average = sum / n;
	}
	void findmed() { //search of the median in the array
		vector <int> v;
		for (int i = 0; i < n; i++)
			v.push_back(values[i]);
		sort(v.begin(), v.end());
		med = v[(n + 1) / 2];
	}
};

#endif 
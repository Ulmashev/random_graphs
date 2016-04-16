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

using namespace std;

#define maxn 1000
#define e 2.7182818284590452
#define eps 10e-7
#define d 10000
#define INF 10e7
#define TP 10

int n;

struct Set_of_degrees {    //��������� ��� �������� ������� �������� � �������� �� �������������� �������������
	double values[maxn];
	double mx, mn, average, med;
	Set_of_degrees() {
		for (int i = 0; i < n; i++)
			values[i] = 0;
		mx = -1;
		mn = maxn + 1;

	}
	void top(int k) { //����� k ��������� ���������
		cout << "Top " << k << endl;
		vector <int> v;
		for (int i = 0; i < n; i++)
			v.push_back(-values[i]);
		sort(v.begin(), v.end());
		for (int i = 0; i < k; i++)
			cout << -v[i] << " ";
		cout << endl;
	}
	pair <int, int> find_vert_from_top(int k, int a, int b) { //����� �-�� � b-�� �������� � ����
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
	void findaver() {
		double sum = 0;
		for (int i = 0; i < n; i++)
			sum += values[i];
		average = sum / n;
	}
	void findmed() {
		vector <int> v;
		for (int i = 0; i < n; i++)
			v.push_back(values[i]);
		sort(v.begin(), v.end());
		med = v[(n + 1) / 2];
	}
};


struct Graph { //��������� ��� �����
	double maxdeg, tau, kappa, uniform, c; // ��������� �������������
	Set_of_degrees ind, outd, sumd; //������� ��� �������� ��������, ��������� � ��������� ��������
	int number_of_edges; //���������� ����� � �����
	int delta;
	
	stack <pair <double, int> > donor; 
	stack <pair <double, int> > recepient;

	vector <pair <double, int> > v; //�������� ������ ������� ��� �������
	
	double p[maxn], threshold[maxn]; //������� ��� �������� ������������ � ������ ������� �������
	int indeg[maxn], outdeg[maxn], sumdeg[maxn]; //������� ��� �������� ��������, ��������� � ��������� �������� (�������� �� ������� ����, �������� ����� ����� ��������� �������� �� ��������� ������)

	vector <int> g[maxn]; //������ ��������� �����
		
	Graph(int MD, double T, double K) {
		tau = T;
		kappa = K;
		//n = N;
		maxdeg = MD;
		uniform = 1 / maxdeg;
		double sum = 0;
		number_of_edges = 0;
		for (int i = 0; i < maxn; i++) { //������������� 
			p[i] = 0;
			threshold[i] = 0;
			indeg[i] = 0;
			outdeg[i] = 0;
			sumdeg[i] = 0;
		}
		for (int i = 1; i <= maxdeg; i++) { //������� ��������������� ������������
			p[i] = 1 / (pow((double)i, tau) * pow(e, (double)i / kappa));
			sum += p[i];
		}
		c = 1 / sum;
		for (int i = 1; i <= maxdeg; i++) { //����������
			p[i] = c / (pow((double)i, tau) * pow(e, (double)i / kappa));
			if (p[i] > uniform && abs(p[i] - uniform) > eps) //���������� ����������� �����
				recepient.push(make_pair(p[i], i));
			else if (p[i] < uniform && abs(p[i] - uniform) > eps) //��������� ������� � ����
				donor.push(make_pair(p[i], i)); 
			threshold[i] = uniform * (i - 1); // ��������� �������������� ������� ������� �� ��������, ������������� ������������ ������������� 
			if (p[i] < uniform)
				threshold[i] += p[i];

			v.push_back(make_pair(uniform * (i - 1), i));
		}
		while (donor.empty() == false && recepient.empty() == false) { // ������������� ������ ������������ �� ������� � ������������
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

		sort(v.begin(), v.end()); //���������� ����������� ������� ������� �� ����������� ������������ 

		srand(time(NULL));

		while (!this->simulate()) {} //���������� ���� �� ��� ���, ���� ��������� ������� �������� � ��������� �� ����� ���������

		ind = Set_of_degrees();
		outd = Set_of_degrees();
		sumd = Set_of_degrees();


		for (int i = 0; i < n; i++) { // ��������� ������ � �������� � ���������
			outd.values[i] = outdeg[i];
			ind.values[i] = indeg[i];
			sumd.values[i] = outdeg[i] + indeg[i];
			number_of_edges += indeg[i];
		}

		int edges_left = number_of_edges; //���������� ��������� �� ���������� ��������������� �����

		while (edges_left != 0) {  // �������� ������� �� ��� ���, ���� �� �������� �����
			vector <pair<int, int> > v_in; //������ ��� ���������������� ������ � �������� �� �������� (����� ������� ��������)
			vector <pair<int, int> > v_out;
			for (int i = 0; i < n; i++) { //���������� ��������
				if (indeg[i] != 0) v_in.push_back(make_pair(i, indeg[i]));
				if (outdeg[i] != 0) v_out.push_back(make_pair(i, outdeg[i]));
			}
			srand(time(NULL));
			int in = -1;
			int out = -1;
			if (v_in.size() == 1 && v_out.size() == 1 && v_in[0].first == v_out[0].first) { //���� �������� ���� �����, ������� ����� ���� ������ ������
				in = v_in[0].first;
				out = v_out[0].first;
				indeg[in]--;
				outdeg[out]--;
				edges_left--;
			}
			else {
				while (in == out) {
					in = v_in[rand() % v_in.size()].first; //��������� ������� �������� �������� �������
					out = v_out[rand() % v_out.size()].first; //��������� ������� �������� ��������� �������
				}
				//if (!existing_edge(out, in)) 
				g[out].push_back(in); //���� ����� ����� ��� ����, �� ��� ������ ��� ��������� � ����
				indeg[in]--;
				outdeg[out]--;
				edges_left--;

			}
		}
		for (int i = 0; i < n; i++) { //��������� � ����� ����� �� �������� �����, �� ����� ����������� �������
			ind.values[i] = 0;
			outd.values[i] = 0;
		}
		for (int i = 0; i < n; i++) {
			outd.values[i] = g[i].size();
			for (int j = 0; j < g[i].size(); j++)
				ind.values[g[i][j]]++;
		}


	}

	bool existing_edge(int out, int in) { //����� ��� �������� ������� ����� � ����� 
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

	void statistics() { //����� ��� ������ ���������� �� �����
		cout << "|V| = " << n << endl;
		cout << "|E| = " << this->number_of_edges << endl;

		ind.findmax(); ind.findmin(); ind.findaver(); ind.findmed();
		outd.findmax(); outd.findmin(); outd.findaver(); outd.findmed();
		sumd.findmax(); sumd.findmin(); sumd.findaver(); sumd.findmed();

		cout << "Max: " << endl <<  "In: " << ind.mx << " Out: " << outd.mx << " Sum: " << sumd.mx << endl;
		cout << "Min: " << endl << "In: " << ind.mn << " Out: " << outd.mn << " Sum: " << sumd.mn << endl;
		cout << "Average: "<< endl << "In: " << ind.average << " Out: " << outd.average << " Sum: " << sumd.average << endl;
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

	int find_path(int a, int b) { // ��������
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
	bool simulate() { //����� ��� ����������� ��� �������� 
		int check_in = 0;
		int check_out = 0;

		for (int i = 0; i < n; i++) {
			indeg[i] = 0;
			outdeg[i] = 0;
		}

		for (int i = 0; i < n; i++) {
			indeg[i] = this->walker();
			outdeg[i] = this->walker();
			check_in += indeg[i];
			check_out += outdeg[i];
		}

		if (check_in == check_out)
			return true;
		else
			return false;
	}
	int walker() { //����� �������
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


int main() {

	freopen("dip.out", "w", stdout);

	n = 10;
	
	Graph gr = Graph(5, 0.52, 2.5);
	gr.show_edges();
	gr.statistics();
	cout << endl;
	cout << "-----------------------------------------------------------";
	cout << endl; 
	//Graph pp = Graph(21, 0.51, 2.5);
	//pp.statistics();

	return 0;
}
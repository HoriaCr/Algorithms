#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <queue>
#include <limits>

using namespace std;

void dijkstra(const vector< vector< pair<int, int> > >& graph, vector<int>& dist) {
	priority_queue< pair<int, int> > q;
	const int src = 0;
	dist[src] = 0;
	q.push({0, src});
	while (!q.empty()) {
		int v = q.top().second;
		int cost = -q.top().first;
		q.pop();
		if (cost != dist[v]) continue;
		
		for (const auto& w : graph[v]) {
			if (dist[w.first] > dist[v] + w.second) {
				dist[w.first] = dist[v] + w.second;
				q.push({-dist[w.first], w.first});
			}
		}
	}
}

int main() {
	ifstream cin("dijkstra.in");
	ofstream cout("dijkstra.out");
	int n, m;
	cin >> n >> m;
	vector< vector< pair<int, int> > > graph(n);
	for (int i = 0; i < m; i++) {
		int x, y, z;
		cin >> x >> y >> z;
		x--, y--;
		graph[x].push_back({y, z});
	}
	
	vector<int> dist(n, numeric_limits<int>::max());
	dijkstra(graph, dist);
	for (int i = 1; i < n; i++) {
		cout << (numeric_limits<int>::max() == dist[i] ? 0 : dist[i]) << " ";
	}
	return 0;
}
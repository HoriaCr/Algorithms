#include <fstream>
#include <vector>

using namespace std;

class DSU {
		vector<int> parent;
		int size;
		
		void compressPaths(const int& x,const int& root) {
			for (int v = x; v != parent[v]; ) {
				int aux = parent[v];
				parent[v] = root;
				v = aux;
			}
		}
		
		int find(const int& x) {
			return x != parent[x] ? find(parent[x]) : x;
		}
		
		int find(const int& x,const bool& compress) {
			int root = find(x);
			if (compress) {
				compressPaths(x, root);
			}
			return root;
		}
		
		
	public:
	
		DSU(int size_) {
			size = size_;
			parent.resize(size);
			
			for (int i = 0; i < size; i++) {
				parent[i] = i;
			}
		}
	
		void unite(const int& a,const int& b) {
			parent[find(b, true)] = find(a, true);
		}
		
		bool inSameSet(const int& a,const int& b) {
			return find(b, true) == find(a, true);
		}
	
};

int main() {
	ifstream cin("disjoint.in");
	ofstream cout("disjoint.out");
	int n, m;
	cin >> n >> m;
	DSU sets(n);
	int t, a, b;
	while (m--) {
		cin >> t >> a >> b;
		a--, b--;
		if (t == 1) {
			sets.unite(a, b);
		} else
		if (t == 2) {
			cout << (sets.inSameSet(a, b) ? "DA\n" : "NU\n"); 
		}
	}
	return 0;
}
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

template<class DataType>
class SegmentTree {
		int size;
		vector< DataType > data;
		
		void init() {
			int p = 1;
            while ( p <= size) {
                p <<= 1;
            }
            p <<= 1;
			data.resize(p);
		}
		
		void buildTree(int node, int l,int r, const vector< DataType >& a) {
			if (l == r) {
				data[node] = a[l - 1];
				return;
			}
			
			int mid = (l + r) / 2;
			buildTree(node * 2, l, mid, a);
			buildTree(node * 2 + 1, mid + 1, r, a);
			data[node] = max(data[node * 2], data[node * 2 + 1]);
		}
		
		int query(int node,int l,int r, int queryLeft,int queryRight) {
			if (queryLeft <= l && r <= queryRight) {
				return data[node];
			}
			int mid = (l + r) / 2;
			int ret = 0;
			if (queryLeft <= mid) {
				ret = max(ret, query(node * 2,l, mid, queryLeft, queryRight));
			}
			if (queryRight > mid) {
				ret = max(ret, query(node * 2 + 1,mid + 1, r, queryLeft, queryRight));
			}
		
			return ret;
		}
		
		void update(int node, int l, int r, const int& position, const int& value) {
			if (l == r) {
				data[node] = value;
				return;
			}
			int mid = (l + r) / 2;
			if (position <= mid) {
				update(node * 2, l, mid, position, value);
			} else {
				update(node * 2 + 1, mid + 1, r, position, value);
			}
			data[node] = max(data[node * 2], data[node * 2 + 1]);
		}
	
	public:
	
		SegmentTree(int size_) {
			size = size_;
			init();
		}
		
		SegmentTree(const vector< DataType >& a) {
			size = static_cast<int>(a.size());
			init();
			buildTree(1, 1, size, a);
		}
		
		void update(const int& position, const int& value) {
			update(1, 1, size, position, value);
		}
		
		int query(const int& l,const int& r) {
			return query(1, 1, size, l, r);
		}
};

int main() {
	ifstream cin("arbint.in");
	ofstream cout("arbint.out");
	int n, m;
	cin >> n >> m;
	vector<int> v(n);
	for (int i = 0; i < n; i++) {
		cin >> v[i];
	}	
	
	SegmentTree<int> tree(v);
	int t, a, b;
	while (m--) {
		cin >> t >> a >> b;
		if (t == 0) {
			cout << tree.query(a, b) << "\n";
		} else
		if (t == 1) {
			tree.update(a, b);
		}
	}
	return 0;
}
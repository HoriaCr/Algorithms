#include <fstream>
#include <vector>

using namespace std;

class BitTree {
	vector<int> data;
	int size;
	
	public:
	
		BitTree(int size_) {
			size = size_;
			data.resize(size_ + 1);
		}
		
		int getMinPosition(int sum) {
			int ret = 0;
			for (int step = 1 << 20; step > 0; step >>= 1) {
				if (ret + step <= size && data[ret + step] <= sum) {
					sum -= data[ret + step];
					ret += step;
					if (!sum) {
						return ret;
					}
				}
			}
			return -1;
		}
		
		void update(const int& position,const int& value) {
			for (int i = position; i <= size; i += (i & -i)) {
				data[i] += value;
			}
		}
		
		int query(const int& position) {
			int ret = 0;
			for (int i = position; i > 0; i -= (i & -i)) {
				ret += data[i];
			}
			return ret;
		}
		
		int query(const int& a,const int& b) {
			if (a > b) return 0;
			return query(b) - query(a - 1);
		}
		
};

int main() {
	ifstream cin("aib.in");
	ofstream cout("aib.out");
	int n, m;
	cin >> n >> m;
	BitTree tree(n);
	for (int i = 1; i <= n; i++) {
		int val;
		cin >> val;
		tree.update(i, val);
	}
	
	for (int i = 1; i <= m; i++) {
		int t, a, b;
		cin >> t;
		if (t == 0) {
			cin >> a >> b;
			tree.update(a, b);
		} else
		if (t == 1) {
			cin >> a >> b;
			cout << tree.query(a, b) << "\n";
		} else
		if (t == 2) {
			cin >> a;
			cout << tree.getMinPosition(a) << "\n";
		}
	}
	return 0;
}
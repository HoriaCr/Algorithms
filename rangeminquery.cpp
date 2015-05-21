#include <fstream>
#include <vector>

using namespace std;

template<class DataType>
class RangeMinimumQuery {
		vector< vector<int> > rmq;
		vector< DataType > data;
		vector< int > logCeil;
		int size;
		int rows;
		
		int getLogCeil(int x) {
			int p = 0;
			while ( (1 << p) < x) {
				p++;
			}
			return p;
		}
		
	public:
		
		
		RangeMinimumQuery(const vector< DataType >& data_) {
			data = data_;
			size = static_cast<int>(data.size());
			rows = getLogCeil(size) + 1;
			rmq.resize(rows + 1);
			logCeil.resize(size + 1);
			rmq[0].resize(size + 1);
			rmq[0][1] = 1;
			for (int i = 2; i <= size; i++) {
				logCeil[i] = logCeil[i >> 1] + 1;
				rmq[0][i] = i;
			}
			
		}
		
		void build() {
			for (int i = 1; i <= rows; i++) {
				rmq[i].resize(size);
				for (int j = 1; j <= size - (1 << i) + 1; j++)  {
					rmq[i][j] = data[rmq[i - 1][j] - 1] < data[rmq[i - 1][j + (1 << (i - 1))] - 1] ? rmq[i - 1][j] : rmq[i - 1][j + (1 << (i - 1))];
				}
			}
		}
		
		DataType query(const int& a,const int& b) {
			int l = logCeil[b - a + 1];
			return min(data[rmq[l][a] - 1], data[rmq[l][b - (1 << l) + 1] - 1]);
		}
};

int main() {
	ifstream cin("rmq.in");
	ofstream cout("rmq.out");
	int n, m;
	cin >> n >> m;
	vector<int> a(n);
	for (int i = 0; i < n; i++) {
		cin >> a[i];
	}
	
	RangeMinimumQuery<int> rmq(a);
	
	rmq.build();
	
	while (m--) {
		int x, y;
		cin >> x >> y;
		cout << rmq.query(x, y) << "\n";
	}
	return 0;
}
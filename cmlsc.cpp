#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

vector<int> lcs(const vector<int>& a,const vector<int>& b) {
	vector< vector<int> > best(a.size() + 1, vector<int>(b.size() + 1));
	int n = static_cast<int>(a.size());
	int m = static_cast<int>(b.size());
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			if (a[i - 1] == b[j - 1]) {
				best[i][j] = 1 + best[i - 1][j - 1];
			} else {
				best[i][j] = max(best[i - 1][j], best[i][j - 1]);
			}
		}
	}
	
	vector<int> ans;
	
	for (int j = n, k = m; best[j][k] > 0 ; ) {
		if (a[j - 1] == b[k - 1]) {
			ans.push_back(a[j - 1]);
			j--;
			k--;
		} else
		if (j > 0 && best[j - 1][k] == best[j][k]) {
			j--;
		} else {
			k--;
		}
	}
	reverse(begin(ans), end(ans));
	
	return ans;
}

int main() {
	ifstream cin("cmlsc.in");
	ofstream cout("cmlsc.out");
	int n, m;
	cin >> n >> m;
	vector<int> a(n), b(m);
	for (int i = 0; i < n; i++) {
		cin >> a[i];
	}
	
	for (int i = 0; i < m; i++) {
		cin >> b[i];
	}
	
	vector<int> ans = lcs(a, b);
	
	cout << ans.size() << "\n";
	for (const int& x : ans) {
		cout << x << " ";
	}
	
	return 0;
}
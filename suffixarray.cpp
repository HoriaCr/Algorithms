#include <iostream>
#include <algorithm>
#include <vector>
#include <tuple>
#include <string>

using namespace std;

int main() {
	const int maxn = 1 << 18;
	vector< vector<int> > sa(18, vector<int>(maxn));
	string str;
	cin >> str;
	int n = static_cast<int>(str.size());
	vector< tuple<int, int, int> > L(n);
	for (int i = 0; i < n; i++) {
		sa[0][i] = str[i] - 'a';
	}

	int step;
	for (step = 1; (1 << (step - 1)) < n; step++) {
		for (int i = 0; i < n; i++) {
			L[i] = make_tuple(sa[step - 1][i], (i + (1 << (step - 1)) < n ? sa[step - 1][i + (1 << (step - 1))] : -1), i);
		}
		sort(L.begin(), L.end());
		for (int i = 1; i < n; i++) {
			sa[step][get<2>(L[i])] = get<0>(L[i]) == get<0>(L[i - 1]) && get<1>(L[i]) == get<1>(L[i - 1]) ? sa[step][get<2>(L[i - 1])] : i;
		}
	}

	cout << "\n";
	for (int i = 0; i < n; i++) {
		cout << sa[step - 1][i] << " ";
	}
	return 0;
}
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <iostream>
#include <cassert>

using namespace std;

class SuffixAutomaton {

	struct Node {
		string str;
		map<char, int> next;
		int link;
		int substringsNumber;
		unsigned short length;
		bool isClone;
		Node() : link(-1), substringsNumber(0), length(0), isClone(false) {}

	};
	vector< Node > state;
	vector<int> cnt;
	int size;
	int last;

	void extend(const string& data) {
		int current;
		int q, p;
		int clone;
		for (const char& c : data) {
			current = size++;
			state[current].length = state[last].length + 1;
			for (p = last; p != -1 && !state[p].next.count(c); p = state[p].link) {
				state[p].next[c] = current;
			}

			if (p == -1) {
				state[current].link = 0;
			}
			else {
				q = state[p].next[c];
				if (state[p].length + 1 == state[q].length) {
					state[current].link = q;
				}
				else {
					clone = size++;
					state[clone].link = state[q].link;
					state[clone].next = state[q].next;
					state[clone].length = state[p].length + 1;
					state[clone].isClone = true;

					for (; p != -1 && state[p].next[c] == q; p = state[p].link) {
						state[p].next[c] = clone;
					}
					state[q].link = state[current].link = clone;
				}
			}
			last = current;
		}
	}

	int getState(const string& str) {
		int v = 0;
		for (size_t i = 0; i < str.size(); i++) {
			auto it = state[v].next.find(str[i]);
			if (it == state[v].next.end()) return -1;
			v = it->second;
		}
		return v;
	}

public:

	SuffixAutomaton(const string &data = "") {
		state.resize(data.size() << 1);
		state[0].link = -1;
		state[0].substringsNumber = 0;
		size = 1;
		last = 0;
		extend(data);
	}

	int df(int v) {
		if (state[v].substringsNumber) {
			return state[v].substringsNumber;
		}

		state[v].substringsNumber = 1;
		for (const auto& w : state[v].next) {
			//if (state[w.second].str.size() < 1) state[w.second].str = state[v].str + w.first;
			//else cout << state[w.second].str << " " << state[v].str + w.first << " kkk\n";
			state[v].substringsNumber += df(w.second);

			// translation
			//cout << state[v].str << " -> " << state[w.second].str << "\n";

		}
		return state[v].substringsNumber;
	}

	void go(int v) {
		//	cout << state[v].str << " : " << state[v].substringsNumber << "\n";
		//if (state[v].link != -1 && v != state[v].link) cout << state[v].str << " --> " << state[state[v].link].str << "\n";
		for (const auto& w : state[v].next) {
			go(w.second);
		}
	}

	bool isIn(const string &str) {
		int v = 0;
		for (size_t i = 0; i < str.size(); i++) {
			auto it = state[v].next.find(str[i]);
			if (it == state[v].next.end()) return false;
			v = it->second;
		}

		return true;
	}

	string getMin(int v, int n) {
		if (!n) return "";
		if (!state[v].next.empty()) {
			return state[v].next.begin()->first + getMin(state[v].next.begin()->second, n - 1);
		}
		return "";
	}

	string getKth(int v, int k) {
		for (auto w = rbegin(state[v].next); w != rend(state[v].next); w++) {
			if (state[w->second].substringsNumber <= k) {
				k -= state[w->second].substringsNumber;
			} else {
				return w->first + getKth(w->second, k);
			}
		}
		return "";
	}

	string getKth(int k) {
		//  substringTotal - substringIndex because of the way the tree is constructed
		// substringTotal is state[0].substringsNumber - 1
		const int substringTotal = state[0].substringsNumber - 1;
		return getKth(0, substringTotal - k - 1);
	}

	void computeCnt() {
		vector<int> p(size);
		vector< Node >& state_ = state;
		cnt.resize(size);
		for (int i = 0; i < size; i++){
			p[i] = i;
			if (state[i].isClone == false) {
				cnt[i] = 1;
			}
		}
		sort(p.begin(), p.end(), [&state_](const int& i, const int& j) -> bool {
			return state_[i].length > state_[j].length;
		});

		for (int i = 0; i < size; i++) {
			if (state[p[i]].link != -1) {
				cnt[state[p[i]].link] += cnt[p[i]];
			}
		}
	}

	int getStrCnt(const string& s) {
		int v = getState(s);
		if (v == -1) return 0;
		return cnt[v];
	}
};

vector<string> brute(string str, ofstream& cout) {
	set<string> s;
	for (int i = 0; i < (int)str.size(); i++) {
		string aux;
		for (int j = i; j < (int)str.size(); j++) {
			aux += str[j];
			s.insert(aux);
		}
	}

	return vector<string>(s.begin(), s.end());
}

string lexminShift(const string& str) {
	SuffixAutomaton sa(str + str);
	return sa.getMin(0, str.size());
}


class SuffixArray {

	vector< vector<int> > P;
	int n, m;
	int getLog(int n) {
		return static_cast<int>(ceil(log2(n))) + 1;
	}

public:

	SuffixArray(const string& str) {
		n = static_cast<int>(str.size());
		vector< pair< pair<int,int>, int> > L(n);
		P.resize(getLog(n) + 1, vector<int>(n));
		for (int i = 0; i < n; i++) {
			P[0][i] = str[i] - 'a';
		}
		int step, cnt;
		for (step = 1, cnt = 1; cnt >> 1 < n; ++step, cnt <<= 1) {
			for (int i = 0; i < n; i++) {
				L[i] = { {P[step - 1][i], i + cnt < n ? P[step - 1][i + cnt] : -1 }, i};
			}
			sort(L.begin(), L.end());
		
			for (int i = 1; i < n; i++) {
				P[step][L[i].second] = L[i].first == L[i - 1].first ? P[step][L[i - 1].second] : i;
			}
			for (int i = 0; i < n; i++) cout << P[step][i] << " "; cout << "\n";
		}

	}

	int lcp(int a, int b) {
		int ret = 0;
		for (int k = m - 1; k >= 0 && a < n && b < n; k--) {
			if (P[k][a] == P[k][b]) {
				ret += 1 << k;
				a += 1 << k;
				b += 1 << k;
			}
		}

		return ret;
	}
};

int getLongestPalindrome(const string& str) {
	string s = str + "$" + str;
	reverse(s.begin() + str.size() + 1, s.end());
	SuffixArray sa(s);
	int n = (int)str.size();
	int ans = 1;
	//baabx xbaab
	//01234 56789
	cout << s << "\n";
	for (int i = 1; i < n - 1; i++) {
		//cout << n - i + 1 << " ";
		int val = sa.lcp(i, 2 * n - i + 1 );
		cout << i << " " << 2 * n - i + 1 << " val :" << val << "\n";
		cout << string(s.begin() + i, s.begin() + str.size()) << 
		"\n" << string(s.begin() + 2 * n - i + 1, s.end()) << "\n";
		
	
	}
	return ans;
}


void testKth(ifstream& cin, ofstream& cout) {
	SuffixAutomaton sa("axadyax");
	vector<string> good = brute("axadyax", cout);
	cout << sa.df(0) - 1 << "\n";
	cout << good.size() << "\n";
	for (auto& x : good) {
		cout << x << " ";
	}

	cout << "\n";
	for (int i = 0; i < (int) good.size(); i++) {
		cout << sa.getKth(i) << " ";
	}
}

int main()
{
	ifstream cin("test.in"); 
	ofstream cout("test.out");
	//testKth(cin, cout);
	cout << getLongestPalindrome("xybaab");
	return 0;
}
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

map<string,int> num;
int index = 0;


class SuffixAutomaton {
	public :
	SuffixAutomaton(const string &data = "") {
		state.resize(data.size() << 1);
		state[0].length = 0;
		state[0].link = -1;
		state[0].stateCount = 0;
		size = 1;
		last = 0;
		for (const char &c : data) {
			extend(c);
		}
	}

	void compute() {

		vector<int> p(state.size());
		for (int i = 0;i < (int)p.size();i++) {
			p[i] = i;
		}

		sort (p.begin(),p.end(),[&](const int &i,const int &j) { 
			return state[i].length < state[j].length;
		});


		for (int i = 0;i < (int)p.size();i++) {
			if (state[p[i]].link != -1 ) {
				state[state[p[i]].link].stateCount += state[p[i]].stateCount;
			}
		}
	}

	void df(int v,string curr) {
		num[curr] += state[v].stateCount;
		//cout << curr <<  ": " << (++index) << "\n";
		cout << curr << " : " << state[v].stateCount << "\n";
		string aux = "";
		for (const auto w : state[v].next) {
			df (w.second,curr + w.first);
		}
	}

	void printResult() {
		cout << state[0].stateCount - 1 << "\n";
	}

	private:
	struct State {
		int length;
		int link;
		int stateCount;
		bool isClone;
		map<char,int> next;
		State() : length(0), link(-1), stateCount(0), isClone(false) {}
	};
	vector< State > state;
	int size, last;

	void extend (char c) {
		int current = size++;
		state[current].length = state[last].length + 1;
		state[current].stateCount = 1;
		int p;
		for (p = last; p != -1 && !state[p].next.count(c); p = state[p].link) {
			state[p].next[c] = current;
		}
	
		if (p == -1) {
			state[current].link = 0;
		} else {
			int q = state[p].next[c];
			if (state[p].length + 1 == state[q].length) {
				state[current].link = q;
			} else {
				int clone = size++;
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

};

int main()
{
	string str = "abcab";
	SuffixAutomaton sa(str);
	sa.compute();
	sa.df(0,"");
	auto f = [](string s,string aux) {
		int ret = 0;
		size_t pos, i = 0;
		while ( (pos = s.find(aux,i)) != string :: npos) {
			ret++;
			i = pos + 1;
		}
		return ret;
	};
	cout << num.size() << " ";
//	for (const auto &x : num) cout << x.first << " -> " << x.second << " ? " << f(str,x.first) << "\n";
	return 0;
}

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

string str;


class SuffixAutomaton {
	public :
	SuffixAutomaton(const string &data = "") {
		state.resize(data.size() << 1);
		state[0].length = 0;
		state[0].link = -1;
		state[0].substringsNumber = 0;
		size = 1;
		last = 0;
		for (const char &c : data) {
			extend(c);
		}
	}

	void df(int v,string curr) {
		state[v].substringsNumber = 1;
		state[v].substringsLength = -1;
		string aux = "";
		for (const auto w : state[v].next) {
			df (w.second,curr + w.first);
			state[v].substringsLength += state[w.second].substringsLength;
			state[v].substringsLength += state[w.second].substringsNumber;
			state[v].substringsNumber += state[w.second].substringsNumber;
		}

		cout << curr << " : " << state[v].substringsLength << "\n";
	}

	void printResult() {
		cout << state[0].substringsNumber - 1 << "\n";
		cout << state[0].substringsLength + state[0].substringsNumber << "\n";
	}

	private:
	struct State {
		int length;
		int link;
		int substringsNumber;
		int substringsLength;
		bool isClone;
		map<char,int> next;
		State() : length(0), link(-1), substringsNumber(0), substringsLength(0), isClone(false) {}
	};
	vector< State > state;
	int size, last;

	void extend (char c) {
		int current = size++;
		state[current].length = state[last].length + 1;
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
	str = "abcassasdcidicjicdjia";
	SuffixAutomaton sa(str);
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

	map<string,int> m;
	for (int i = 0;i < str.size();i++) {
		for (int j = 1;j <= str.size() - i;j++) {
			m[str.substr(i,j)] = 1;
		}
	}
	int sum = 0;
	for (auto x : m) sum += x.first.size();

	cout << m.size() << "\n";
	cout << sum << "\n";
	sa.printResult();
	return 0;
}

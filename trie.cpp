#include <fstream>
#include <map>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

class Trie {

	struct TrieNode {
		map<char, TrieNode*> next;
		int wordCount;
		int wordsInTree;
		TrieNode() {
			wordsInTree = 0;
			wordCount = 0;
		}
	};

	TrieNode* root;

	void updateCount(const string& str, const int& value) {
		TrieNode *node = root;
		for (size_t i = 0; i < str.size(); i++) {
			node->wordsInTree += value;
			if (node->next.find(str[i]) == node->next.end()) {
				node->next[str[i]] = new TrieNode();
			}
			node = node->next[str[i]];
		}

		node->wordCount += value;
		node->wordsInTree += value;
	}

	void print(TrieNode* node, string str) {
		if (node->wordCount > 0) {
			cout << str << "\n";
		}
		for (auto& w : node->next) {
			print(w.second, str + w.first);
		}
	}

public:

	Trie() {
		root = new TrieNode();
	}

	void insert(const string& str) {
		updateCount(str, 1);
	}

	void remove(const string& str) {
		updateCount(str, -1);
	}

	int queryLcp(const string& str) {
		TrieNode *node = root;
		int ans = 0;
		for (size_t i = 0; i < str.size(); i++) {
			if (node->next.find(str[i]) == node->next.end()) {
				break;
			}
			node = node->next[str[i]];
			if (node->wordsInTree > 0)
				ans = static_cast<int>(i + 1);			
		}

		return ans;
	}

	int queryWordCount(const string& str) {
		TrieNode *node = root;
		for (size_t i = 0; i < str.size(); i++) {
			auto it = node->next.find(str[i]);
			if (it == node->next.end()) {
				return 0;
			}
			node = node->next[str[i]];
		}
		return node->wordCount;
	}

	void printWords() {
		print(root, "");
		cout << "###\n";
	}

};

int main() {
	ifstream cin("trie.in");
	ofstream cout("trie.out");
	Trie trie;
	string str;
	while (getline(cin, str)) {
		string w = string(str.begin() + 2, str.end());
		if (str[0] == '0') {
			trie.insert(w);
		}
		else
		if (str[0] == '1') {
			trie.remove(w);
		}
		else
		if (str[0] == '2') {
			cout << trie.queryWordCount(w) << "\n";
		}
		else
		if (str[0] == '3') {
			cout << trie.queryLcp(w) << "\n";
		}
	}
	return 0;
}
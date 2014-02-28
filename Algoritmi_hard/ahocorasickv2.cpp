#include <fstream>
#include <vector>
#include <queue>
#include <cstring>
#include <algorithm>
#define ch(s) (s - 'a')

using namespace std; 

ifstream cin("test.in");
ofstream cout("test.out");

struct trie {
	int val;
	trie *son[26], *next;
	vector<trie*> M;
	trie() {
		memset(son,0,sizeof(son));
		next = 0;
		val = 0;
	}
};

const int AMAX = int(1e6) + 2, LDMAX = int(1e4) + 2, NMAX = 102, SIGMA = 26;
int N;
char A[AMAX], S[NMAX][LDMAX];
trie *T, *leaf[NMAX];
queue<trie*> Q;

trie *add(trie *t,char *str) {
	if(*str == 0) {
		return t;
	}
	if(t->son[ch(*str)] == 0) {
		t->son[ch(*str)] = new trie;
	}
	t = t->son[ch(*str)];
	str++;
	return add(t,str);
}

void readData() {
	cin.getline(A,AMAX);
	cin>>N;
	cin.get();
	for(int i = 0;i < N;i++) {
		cin.getline(S[i],LDMAX);
	}
}

void buildTrie() {
	T = new trie;
	for(int i = 0;i < N;i++) {
		leaf[i] = add(T,S[i]);
	}
}

void bfs() {
	Q.push(T);
	trie *v, *w;
	while(!Q.empty()) {
		v = Q.front();
		Q.pop();
		for(int i = 0;i < SIGMA;i++) {
			if(v->son[i] != 0) {
				for(w = v->next;w && w->son[i] == 0;) {
					w = w->next;
				}
				if(w) {
					v->son[i]->next = w->son[i];
					w->son[i]->M.push_back(v->son[i]);
				} else {
					v->son[i]->next = T;
					T->M.push_back(v->son[i]);
				}
				Q.push(v->son[i]);
			}
		}
	}
}

void df(trie *v) {
	for(vector<trie*>::const_iterator w = v->M.begin();w != v->M.end();w++) {
		df(*w);
		v->val += (*w)->val;
	}
}

void ahoCorasick() {
	buildTrie();
	bfs();

	trie *v = T;

	for(int i = 0;A[i] != 0;i++) {
		while(v && v->son[ch(A[i])] == 0) {
			v = v->next;
		}
		if(v == 0) {
			v = T;
		} else {
			v = v->son[ch(A[i])];
			v->val++;
		}
	}

	df(T);

	for(int i = 0;i < N;i++) {
		cout<<leaf[i]->val<<"\n";
	}
}

int main()
{
	readData();
	ahoCorasick();
	return 0;
}

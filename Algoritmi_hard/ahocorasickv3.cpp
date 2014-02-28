#include <fstream>
#include <vector>
#include <cstring>
#define ch(c) (c - 'a')

using namespace std;

ifstream cin("ahocorasick.in");
ofstream cout("ahocorasick.out");

struct trie {
	trie *son[26], *Next;
	int val;
	trie() {
		Next = 0;
		memset(son,0,sizeof(son));
		val = 0;
	}
};

const int SMAX = 1000002, LMAX = 10002, NMAX = 102;
char A[SMAX], S[NMAX][LMAX];
int N;
trie *T, *leaf[NMAX];
vector<trie*> Q;

trie *insert(trie *t,char *str) {
	if(*str == 0) return t;
	if(t->son[ch(*str)] == 0) t->son[ch(*str)] = new trie;
	t = t->son[ch(*(str++))];
	return insert(t,str);
}

void readData() {
	cin.getline(A,SMAX);
	cin>>N;
	cin.get();
	for(int i = 0;i < N;i++) {
		cin.getline(S[i],LMAX);
	}
}

void buildTrie() {
	T = new trie;
	for(int i = 0;i < N;i++) {
		leaf[i] = insert(T,S[i]);
	}
}

void ahoCorasick() {
	buildTrie();
	Q.push_back(T);
	trie *v, *w;
	int Left = 0;
	while(Left < (int)Q.size()) {
		v = Q[Left++];
		for(int i = 0;i < 26;i++) {
			if(v->son[i] != 0) {
				for(w = v->Next;w && w->son[i] == 0;) w = w->Next;
				v->son[i]->Next = w ? w->son[i] : T;
				Q.push_back(v->son[i]);
			}
		}
	}
	v = T;
	for(int i = 0;A[i] != 0;i++) {
		while(v && v->son[ch(A[i])] == 0) v = v->Next;
		if(!v) {
			v = T;
		} else {
			v = v->son[ch(A[i])];
			v->val++;
		}
	}
	for(int i = (int)Q.size() - 1;i >= 0;i--) {
		if(Q[i]->Next != 0)	 {
			Q[i]->Next->val += Q[i]->val;
		}
	}
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

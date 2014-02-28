#include <cstdio>
#include <vector>

using namespace std;
 
const int NMax = int(1e5);
int N, M;
vector<int> Val, G[NMax], Lev, weight, Chain, First, Bos;
vector< vector<int> > CSet, STree;

inline int position(const int &x){
    return Lev[x] - Lev[First[Chain[x]]];
}
 
int queryTree(const int &chain, const int &node, const int &left, const int &right, const int &a, const int &b) {
    if (a <= left && right <= b)
        return STree[chain][node];
 
    int mid = (left + right) >> 1;
    int lv = -1, rv = -1;
    if (a <= mid) lv = queryTree(chain, node << 1 | 1, left, mid, a, b);
    if (b > mid) rv = queryTree(chain, (node + 1) << 1, mid + 1, right, a, b);
    return max(lv, rv);
}

int query(const int &x, const int &y) {
    if (x == y) return Val[x];
    int cx = Chain[x], cy = Chain[y];

    if (cx == cy) {
        int px = position(x), py = position(y);
        return queryTree(cx, 0, 0, CSet[cx].size() - 1, min(px, py), max(px, py));
    }
    int qquery = Lev[Bos[cx]] > Lev[Bos[cy]] ? query(Bos[cx], y) : query(Bos[cy], x);
    int tquery = Lev[Bos[cx]] > Lev[Bos[cy]] ? 
		queryTree(cx, 0, 0, CSet[cx].size() - 1, 0, position(x)):
		queryTree(cy, 0, 0, CSet[cy].size() - 1, 0, position(y));
 
    return max(qquery, tquery);
}
 
void updateTree(const int &chain, const int &node, const int &left, const int &right, const int &a, const int &b, const int &v) {
    if (a <= left && right <= b) {
        STree[chain][node] = v;
        return;
    }
 
    int mid = (left + right) >> 1;
    if (a <= mid) updateTree(chain, (node << 1) | 1, left, mid, a, b, v);
    if (b > mid) updateTree(chain, (node + 1) << 1, mid + 1, right, a, b, v);
    STree[chain][node] = max(STree[chain][(node << 1) | 1], STree[chain][(node + 1) << 1]);
}
 
void update(const int &x, const int &y) {
    int c = Chain[x];
    Val[x] = y;
    updateTree(c, 0, 0, CSet[c].size() - 1, position(x), position(x), y);
}
 
void dfs(const int &v, const int &parent, const int &currentLevel) {
    int c, s = 0, m = 0, r = 0, d = 0;
    Lev[v] = currentLevel, weight[v] = 1;
 
	//leaf
    if (G[v].size() == 1 && G[v][0] == parent) { 
        CSet.push_back(vector<int>(1, v));
        Chain[v] = (int)(CSet.size() - 1);
        First.push_back(v);
        Bos.push_back(N);
        return;
    }

    for (auto w : G[v]) {
        if (w != parent) {
            dfs(w, v, currentLevel + 1);
            s += (d = weight[w]);
            if (d > m)
                m = d, r = Chain[w];
        }
    }
 
    CSet[r].push_back(v);
    weight[v] += s;
    First[r] = v;
    Chain[v] = r;
 
    for (auto w : G[v]) {
        if (w != parent && Chain[w] != r)
            Bos[Chain[w]] = v;
    }
}
 
void buildTree(const int &c, const int &v, const int &l, const int &r) {
    if (l == r) {
        STree[c][v] = Val[CSet[c][l]];
        return;
    }

    int m = (l + r) >> 1;
    buildTree(c, (v << 1) | 1, l, m);
    buildTree(c, (v + 1) << 1, m + 1, r);
    STree[c][v] = max(STree[c][(v << 1) | 1], STree[c][(v + 1) << 1]);
}
 
void buildSegmentTree(const int &c) {
    int s = CSet[c].size()<<1;
    s & (s - 1)?STree.push_back(vector<int>((s & (s - 1)) << 1, 0)):STree.push_back(vector<int>(s, 0));
    buildTree(c, 0, 0, CSet[c].size() - 1);
}
 
int main()
{
    int i, a, b;
    freopen("test.in", "r", stdin);
    freopen("test.out", "w", stdout);

    scanf("%d%d", &N, &M);
    Val.assign(N, -1);
    for (i = 0; i < N; ++i)
        scanf("%d", &Val[i]);
 
    for (i = 1; i < N; ++i) {
        scanf("%d%d", &a, &b);
        --a, --b;
        G[a].push_back(b);
        G[b].push_back(a);
    }
 
    Lev.assign(N, 0);
    weight.assign(N, 0);
    Chain.assign(N, -1);
 
    dfs(0, -1, 0);
    Lev.push_back(-1);
 

    for (i = 0; i < CSet.size(); ++i) {
        reverse(CSet[i].begin(), CSet[i].end());
        buildSegmentTree(i);
    }

    while (M--) {
        scanf("%d%d%d", &i, &a, &b);
        i? (void)printf("%d\v", query(a - 1, b - 1)) : update(a - 1, b);
    }  
 
    return 0;
}


#include <fstream>
#include <vector>
#include <algorithm>

#define REP(i,n) for(int i = 0;i < (int)n;i++)

using namespace std;

ifstream cin("heavypath.in");
ofstream cout("heavypath.out");

const int nmax = int(1e5) + 2;
int n, m;
int value[nmax];
int level[nmax];
int weight[nmax];
int nodeIndex[nmax];
int pathParent[nmax];
int nodePath[nmax];
int pathNumber;
vector< vector<int> > paths, segmentTrees;
vector<int> G[nmax];


inline void readData() {
	cin>>n>>m;
	for(int i = 1;i <= n;i++) {
		cin>>value[i];
	}
	int a, b;
	for(int i = 1;i < n;i++) {
		cin>>a>>b;
		G[a].push_back(b);
		G[b].push_back(a);
	}
}

void df(int v,int vp) {
	weight[v] = 1;
	int heavyestSon = -1;
	for(auto w : G[v]) {
		if(w == vp) continue;
		level[w] = level[v] + 1;
		df(w,v); 
		weight[v] += weight[w];
		if(heavyestSon == -1 || weight[heavyestSon] < weight[w]) {
			heavyestSon = w;
		}
	}

	if(heavyestSon == -1) {
		paths.push_back(vector<int>());
		paths.back().push_back(v);
		nodePath[v] = pathNumber++;
	} else {
		nodePath[v] = nodePath[heavyestSon];
		paths[nodePath[v]].push_back(v);
		for(auto w : G[v]) {
			if(w == vp || w == heavyestSon) continue;
			pathParent[nodePath[w]] = v;
		}
	}
}

void buildTree(const int &treeIndex,int node,int l,int r) {
	if(l > r) return;
	if(l == r) {
		segmentTrees[treeIndex][node] = value[paths[treeIndex][l - 1]];
		return;
	}
	int mid = (l + r)>>1;
	buildTree(treeIndex,node<<1,l,mid);
	buildTree(treeIndex,node<<1 | 1,mid + 1,r);
	segmentTrees[treeIndex][node] = max(segmentTrees[treeIndex][node<<1],segmentTrees[treeIndex][node<<1 | 1]);
}

inline void buildTrees() {
	REP(i,pathNumber) {
		reverse(paths[i].begin(),paths[i].end());
		REP(j,paths[i].size()) {
			nodeIndex[paths[i][j]] = j + 1;
		}
		segmentTrees.push_back(vector<int>(paths[i].size()<<2,0));
		buildTree(i,1,1,paths[i].size());
	}
}

void updateTree(const int &treeIndex,int node,int l,int r,int pos) {
	if(l == r) {
		segmentTrees[treeIndex][node] = value[paths[treeIndex][l - 1]];
		return;
	}
	int mid = (l + r)>>1;
	pos <= mid ? updateTree(treeIndex,node<<1,l,mid,pos) : updateTree(treeIndex,node<<1 | 1,mid + 1,r,pos);
	segmentTrees[treeIndex][node] = max(segmentTrees[treeIndex][node<<1],segmentTrees[treeIndex][node<<1 | 1]);
}

inline void update(const int &v,const int &newValue) {
	value[v] = newValue;
	updateTree(nodePath[v],1,1,paths[nodePath[v]].size(),nodeIndex[v]);
}

int queryTree(const int &treeIndex,int node,int l,int r,const int &ql,const int &qr) {
	if(r < ql || l > qr ) return 0;
	if(ql <= l && r <= qr) {
		return segmentTrees[treeIndex][node];
	}
	int mid = (l + r)>>1;
	return max(queryTree(treeIndex,node<<1,l,mid,ql,qr),queryTree(treeIndex,node<<1 | 1,mid + 1,r,ql,qr));
}

int query(int a,int b) {
	if(a == b) {
		return value[a];
	}
	
	if(nodePath[a] == nodePath[b]) {
		int posA = nodeIndex[a];
		int posB = nodeIndex[b];
		if(posA > posB) swap(posA,posB);
		return queryTree(nodePath[a],1,1,paths[nodePath[a]].size(),posA,posB);
	}

	if(level[pathParent[nodePath[a]]] > level[pathParent[nodePath[b]]]) {
		return max(query(pathParent[nodePath[a]],b),queryTree(nodePath[a],1,1,paths[nodePath[a]].size(),1,nodeIndex[a]));
	}

	return max(query(pathParent[nodePath[b]],a),queryTree(nodePath[b],1,1,paths[nodePath[b]].size(),1,nodeIndex[b]));
}

int main()
{
	readData();

	level[1] = 1;
	df(1,1);

	buildTrees();

	for(int t, a, b;m;m--) {
		cin>>t>>a>>b;
		if(!t) {
			update(a,b);
		} else {
			cout<<query(a,b)<<"\n";
		}
	}
    return 0;
}


#include <iostream>
#include <vector>
#include <algorithm>
#include <bitset>

using namespace std; 

const int NMAX = 1024 , INF = int(2e9);
vector<int> G[NMAX];
bitset<NMAX> visited;
vector<int> nodes;
int cost[NMAX][NMAX]

void dfs(const int &v,vector<int> G[NMAX]) {
   visited[v] = true;
   for(vector<int>::const_iterator w = G[v].begin();w != G[v].end();++w) {
	if(!visited[*w]) 
		dfs(*w,G);
   }  
   nodes.push_back(v);	
}


int bfs(const int &S,const int &D,vector<int> G[NMAX]) {
  int *Q , *dist , st ,  dr;
  Q = new int[NMAX + 2];
  dist = new int[NMAX + 2];
  memset(dist,INF,sizeof(dist));
  dist[S] = 0;
  st = dr = 0;
  Q[++dr] = S;
  while(st < dr) {
     int v = Q[++st];
     for(vector<int>::const_iterator w = G[v].begin();w != G[v].end();++w) {
	  if(dist[*w] > dist[v] + cost[v][*w]) {
		dist[*w] = dist[v] + cost[v][*w];
		Q[++dr] = *w;
		}
   	}  		
  }
  return dist[D];
}
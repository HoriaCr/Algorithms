#include <cstdio>
#include <algorithm>
 
using namespace std;
 
const int NMAX = 305, inf = int(1e9);
 
int n, w[NMAX][NMAX], ax[NMAX], ay[NMAX];
int edgeIndex[NMAX][NMAX];
int edges[NMAX], minCost;
bool vx[NMAX], vy[NMAX];
int mch[NMAX];
int N, M, E;
 
bool dfs(int x){
    vx[x] = 1;
    for(int y = 0; y < n; y++){
        int t = ax[x] + ay[y] - w[x][y];
        if(!vy[y] && t == 0){
            vy[y] = 1;
            if(mch[y] < 0 || dfs(mch[y])) {
                mch[y] = x;
                return 1;
            }
        }
    }
    return 0;
}
 
void khun(){
    fill(mch, mch + n, -1);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            ax[i] = max(ax[i], w[i][j]);
        }
    }
    for(int i = 0; i < n; i++){
        while(!dfs(i)){
            int d = inf;
            for(int j = 0; j < n; j++) {
                if(!vy[j]) {
                    for(int k = 0; k < n; k++) {
                        if(vx[k]) {
						 d = min(d, ax[k] + ay[j] - w[k][j]);
                        }
                    }
                }
            }
            for(int j = 0; j < n; j++) {
				if(vx[j]) ax[j] -= d, vx[j] = 0;
				if(vy[j]) ay[j] += d, vy[j] = 0;
			}
        }
    }
}
 
int main ()
{
    freopen("cmcm.in","r",stdin);
    freopen("cmcm.out","w",stdout);
    scanf("%d %d %d",&N,&M,&E);
    n = max(N,M);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            w[i][j] = -inf;
        }
    }
 
    for(int a, b, c, i = 0;i < E;i++) {
        scanf("%d %d %d",&a,&b,&c);
        a--;b--;
        w[a][b] = -c;
        edgeIndex[a][b] = i + 1;
    }
 
    khun();
 
    for(int i = 0; i < n; i++) {
        if(w[mch[i]][i] != -inf) {
            minCost += w[mch[i]][i];
            edges[++edges[0]] = edgeIndex[mch[i]][i];;
        }
    }
    printf("%d %d\n",edges[0],-minCost);
    for(int i = 1;i <= edges[0];i++) {
        printf("%d ",edges[i]);
    }
    return 0;
}

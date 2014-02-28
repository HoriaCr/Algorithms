/*
    * All pairs max flow
    * Algorithm : Gusfield's algorithm
    * Order : N^3 + cost of N max flow
*/
     
#include<stdio.h>
#include<string.h>
#include<vector>
#include<algorithm>
using namespace std;
     
#define MAX_V
#define MAX_E
#define INF 7777777
     
struct EDGE{
        long v,c;
};
     
long nV;
long SRC,TNK;
long eI;
EDGE Init[MAX_E+7]; // to keep the initial edge cap
EDGE Edge[MAX_E+7];     // edge list
long Next[MAX_E+7]; // next pointer of vertex v
long Last[MAX_V+7]; // last index of adj edge of vertex v
long Dist[MAX_V+7];     // level from src
long Start[MAX_V+7];// temporary used for last
     
long Par[MAX_V+7];
long Flow[MAX_V+7][MAX_V+7]; // to keep the max flow between i th node and j th node
     
inline void SetEdge( long u,long v,long c )
{
        Edge[eI].v = v;
        Edge[eI].c = c;
        Next[eI] = Last[u];
        Last[u] = eI++;
        Edge[eI].v = u;
        Edge[eI].c = 0;
        Next[eI] = Last[v];
        Last[v] = eI++;
}
     
long Q[MAX_V+7];
long Frnt,End;
     
bool Bfs( void )
{
        Frnt = End = 0;
        Q[End++] = SRC;
        long i;
        for( i=0;i<nV;i++){
                Dist[i] = INF;
        }
        Dist[SRC] = 0;
        long u,v;
        while( Frnt<End ){
                u = Q[Frnt++];
                for( i=Last[u];i!=-1;i=Next[i]){
                        v = Edge[i].v;
                        if( !Edge[i].c || Dist[v]<INF ) continue;
                        Dist[v] = Dist[u] + 1;
                        Q[End++] = v;
                }
        }
        return Dist[TNK] < INF;
}
     
long AugmentPath( long u,long f )
{
        if( u==TNK ) return f;
        long Tot = 0;
        for( long &i = Start[u];i!=-1;i=Next[i] ){
                long v = Edge[i].v;
                if( !Edge[i].c ) continue;
                if( Dist[v] != Dist[u]+1 ) continue;
                long Tmp = AugmentPath( v,min( f,Edge[i].c ));
                Edge[i].c -= Tmp;
                Edge[i ^ 1].c += Tmp;
                Tot += Tmp;
                f -= Tmp;
                if( !f ) break;
        }
        return Tot;
}
     
long MaxFlow( void )
{
        long Flw = 0;
        while( Bfs()){
                memcpy( Start,Last,(nV)*sizeof(long));
                Flw += AugmentPath( SRC,2*INF );
        }
        return Flw;
}
     
int main( void )
{
    long i,j;
    eI = 0;
    memset( Last,-1,nV*sizeof(long));
    memcpy( Init,Edge,eI*sizeof(EDGE));
    /* set parent of each vertex to 0 for o based indexin
        if 1th based index is used then parent will be 1 */
    memset( Par,0,sizeof(Par));
    /* set Flow to INF for all vertex pair */
    for( i=0;i<nV;i++ ){
        for( j=0;j<nV;j++ ){
            Flow[i][j] = INF;
        }
    }
     
    for( i=1;i<nV;i++ ){
    /* Compute the minimum cut between i and parent[i].
        Let the i-side of the min cut be S, and the value of the min-cut be F */
        SRC = i;
        TNK = Par[i];
        memcpy( Edge,Init,eI*sizeof(EDGE));
        long F = MaxFlow();
        //Bfs();
        for( j=i+1;j<nV;j++ ){
            /*if ((j is in S) && parent[j]==parent[i])
                parent[j]=i;
            */
            if( Dist[j]<INF && Par[j]==Par[i] ){
                Par[j] = i;
            }
        }
        Flow[i][Par[i]] = Flow[Par[i]][i] = F;
        for( j=0;j<i;j++ ){
            Flow[i][j] = Flow [j][i] = min( F,Flow[Par[i]][j] );
        }
    }
    /*
        Flow[i][j] with contain all pairs of flow
    */
        return 0;
}

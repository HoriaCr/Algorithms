/* 
 * Strongly Connected Component
 * Algorithm : Tarjan's Algorithm
 * Order : O( V+E )
 */
#include<stdio.h>
#include<string.h>
#include<vector>
#include<stack>
#include<algorithm>
using namespace std;
#define MAX 4007
#define pb push_back
long N,M;
vector<long> Edge[MAX+7];
bool Visit[MAX+7];
bool InStk[MAX+7];
long Low[MAX+7],I;
long Ind[MAX+7];
stack<long> Stk;
void SCC( long u )
{
    Visit[u] = true;
    InStk[u] = true;
    Ind[u] = ++I;
    Low[u] = I;
    Stk.push( u );
    long i;
    for( i=0;i<Edge[u].size();i++){
        long v = Edge[u][i];
        if( !Visit[v] ){
            SCC( v );
            Low[u] = min( Low[u],Low[v] );
        }
        else if( InStk[v] ){
            Low[u] = min( Low[u],Ind[v] );
        }
    }
    if( Low[u]!=Ind[u] ) return;
    // found new component
    while( Stk.top()!=u ){
        long v = Stk.top();
        Stk.pop();
        InStk[v] = false;
    }
    Stk.pop();
    InStk[u] = false;
}

int main( void )
{
    long i,j,u,v,Icase,k = 0;
    //freopen("text1.txt","r",stdin );
    scanf("%ld%ld",&N,&M );
    for( i=1;i<=M;i++){
        scanf("%ld%ld",&u,&v );
        Edge[u].pb( v );
    }
    for( i=1;i<=N;i++){
        if( Visit[i] ) continue;
        SCC( i );
    }
    return 0;
}

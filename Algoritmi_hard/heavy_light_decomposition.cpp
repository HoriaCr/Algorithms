/*
	Algorithm: heavy-light decompositon
	Order: O(N) / O(nlogn) if any update is needed
*/

#include<stdio.h>
#include<string.h>
#include<vector>
#include<algorithm>
using namespace std;

#define MAX 30007

long N; // number of node in tree
vector<long> Edge[MAX+7];
long SubT[MAX+7]; // subtree size
long Par[MAX+7]; // parent of a node
long Level[MAX+7]; // level of a node

long nC; // number of chain
long ChainLdr[MAX+7]; // chainleadr of a node
					  // for light edge chainldr is that node
long Chain[MAX+7];	  // node v in is which chain
long nP;			// number of position , obviously == N
long Pos[MAX+7];    // Pos of a node in chain/dfs order

/* find subtree size and level */
long Explore( long u,long p,long l )
{
	SubT[u] = 1;
	Par[u] = p;
	Level[u] = l;
	long i;
	for( i=0;i<Edge[u].size();i++){
		long v = Edge[u][i];
		if( p==v ) continue;
		SubT[u] += Explore( v,u,l+1 );
	}
	return SubT[u];
}
/* if IsL make this node a chainledr of new chain */
void HeavyLight( long u,long k,bool IsL )
{
	if( IsL ){
		k = ++nC;
		ChainLdr[k] = u;
	}
	Chain[u] = k;
	Pos[u] = ++nP;
	//Update( nP,W[u] ); if query is need can b updated here 	
	long i,mx = -1; // max subtree size child is mx
	for( i=0;i<Edge[u].size();i++){
		long v = Edge[u][i];
		if( Par[u]==v ) continue;
		if( mx==-1 ) mx = v;
		else if( SubT[v] > SubT[mx] ) mx = v;
	}
	if( mx==-1 ) return;
	HeavyLight( mx,k,false );
	for( i=0;i<Edge[u].size();i++){
		long v = Edge[u][i];
		if( Par[u]==v || mx == v ) continue;
		HeavyLight( v,0,true );
	}
}

long LCA( long u,long v )
{
    while( Chain[u]!=Chain[v] ){
        if( Level[ChainLdr[Chain[u]]] < Level[ChainLdr[Chain[v]]] ){
            v = Par[ChainLdr[Chain[v]]];
        }
        else{
            u = Par[ChainLdr[Chain[u]]];
        }
    }
    if( Level[u] < Level[v] ) return u;
    else return v;
}


int main( void )
{
	//freopen("text1.txt","r",stdin );
	// node is index from 0
	scanf("%ld",&N );
	Explore( 0,0,0 );
	HeavyLight( 0,0,true );
	return 0;
}

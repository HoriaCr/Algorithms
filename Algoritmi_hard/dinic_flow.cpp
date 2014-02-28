/*
	Algorithm: Dinic'c max flow
			   using Edge List 
*/

#include<stdio.h>
#include<string.h>

#define MAX_V 
#define MAX_E 
#define INF 7777777

struct EDGE{
	long v,c;
};
long nV;
long SRC,TNK;
long eI;
EDGE Edge[MAX_E+7];	// edge list
long Next[MAX_E+7]; // next pointer of vertex v
long Last[MAX_V+7]; // last index of adj edge of vertex v
long Dist[MAX_V+7];	// level from src 
long Start[MAX_V+7];// temporary used for last 

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
	for( i=1;i<=TNK;i++){
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
	return Dist[TNK] < INF;;
}

#define MIN( a,b ) a<b ? a:b

long AugmentPath( long u,long f )
{
	if( u==TNK ) return f;
	long Tot = 0;
	for( long &i = Start[u];i!=-1;i=Next[i] ){
		long v = Edge[i].v;
		if( !Edge[i].c ) continue;
		if( Dist[v] != Dist[u]+1 ) continue;
		long Tmp = AugmentPath( v,MIN( f,Edge[i].c ));
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

long main( void )
{
	//freopen("text1.txt","r",stdin );
	SRC = 0;
	TNK = N+E;
	nV = TNK+1;
	eI = 0;
	memset( Last,-1,nV*sizeof(long));
	MaxFlow()){
	return 0;
}

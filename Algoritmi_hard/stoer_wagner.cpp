////////////////////////////////////////////////////////////////////////////////
// All pairs min cut
////////////////////////////////////////////////////////////////////////////////
#include<stdio.h>
#include<algorithm>
using namespace std;

#define INF 777777777
#define MAX_V 1007

long Cap[ MAX_V+7][MAX_V+7];
long Dist[ MAX_V+7],Done[MAX_V+7];
long nVertex,nEdge;

long Best, Ind;

int main( void)
{
	long Icase,k=0;
	long i,j;
	long u,v,C;

	scanf("%ld",&Icase);
	while(Icase--){

		scanf("%ld%ld",&nVertex, &nEdge);
		memset(Cap, 0, sizeof(Cap));

		while(nEdge--){
			scanf("%ld%ld%ld", &u, &v, &C);
			Cap[u][v] += C;
			Cap[v][u] += C;
		}

		long Ans = INF;
		while(nVertex > 1){

			memset(Dist, 0, sizeof(Dist));
			memset(Done, 0, sizeof(Done));
			for( i=1;i<=nVertex; i++){

				Best = Ind = -1;
				for( j=1;j<=nVertex; j++){
					if(!Done[j] && Dist[j] > Best){
						Best = Dist[j];
						Ind = j;
					}
				}

				if( i+1 == nVertex){
					u = Ind;
				}

				if( i== nVertex){
					v = Ind;
					//Ans <?= Best;
					Ans =Best<Ans ? Best:Ans;
				}

				Done[Ind] = 1;
				for( j=1; j<=nVertex; j++){
					Dist[j] += Cap[Ind][j];
				}

			}

			if(u > v){
				swap( u,v);
			}

			for( i=1;i<=nVertex;i++){
				Cap[u][i] += Cap[v][i];
				Cap[i][u] += Cap[i][v];
			}
			for( i=1;i<=nVertex;i++){
				Cap[v][i] = Cap[nVertex][i];
				Cap[i][v] = Cap[i][nVertex];
			}
			nVertex--;

		}

		printf("Case #%ld: %ld\n", ++k, Ans);

	}
	return 0;
}

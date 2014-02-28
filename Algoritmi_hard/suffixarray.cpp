/*
 *	Algorithm: Suffix array
 *	Order: nlogn
*/

#include<stdio.h>
#include<string.h>
#include<vector>
#include<algorithm>
using namespace std;

#define MAX 200007

struct DATA{
	long Val[2];
	long Ind;
};
/* is needed if qsort is used */
bool operator<( const DATA &a,const DATA &b )
{
	if( a.Val[0] != b.Val[0] ) return a.Val[0] < b.Val[0];
	else return a.Val[1] < b.Val[1];
}
char Str[MAX+7];
long Len;
long Buc[22][MAX+7];
long SuffA[MAX+7]; // will contain index of all suffix in sorted order
long Stp;         

long Cnt[MAX+7],Cum[MAX+7];

void CountSort( vector<DATA> &Inf,long I )
{
	// max is used to count max val is used
	// +1 adding is needed to handle -1 
	long i,Max = 0;
	for( i=0;i<Inf.size();i++){
		if( Max < Inf[i].Val[I] + 1 ) Max = Inf[i].Val[I] + 1;
	}

	memset( Cnt,0,(Max+1)*sizeof(long));
	for( i=0;i<Inf.size();i++){
		Cnt[Inf[i].Val[I] + 1 ]++;
	}

	Cum[0] = Cnt[0];
	for( i=1;i<=Max;i++){
		Cum[i] = Cum[i-1] + Cnt[i];
	}

	vector<DATA> Tmp;
	Tmp.resize( Inf.size());
	for( i=Inf.size()-1;i>=0;i--){
	    Cum[ Inf[i].Val[I] + 1 ]--;
		Tmp[ Cum[ Inf[i].Val[I] +1 ] ] = Inf[i];
		
	}
	Inf = Tmp;
}

void SuffArrCreate( void )
{
	long i,j;

	for( i=0;i<Len;i++){
		Buc[0][i] = Str[i];
	}
	vector<DATA> Inf;
	Inf.resize( Len );
	long c;
	for( c=1,Stp=1;c < Len;c<<=1,Stp++ ){

		for( j=0;j<Len;j++){
			Inf[j].Val[0] = Buc[Stp-1][j];
			Inf[j].Val[1] = j+c < Len ? Buc[Stp-1][j+c]:-1;
			Inf[j].Ind = j;
		}
		// if O( nlogn^2 ) pos 
		//sort( Inf.begin(),Inf.end());
		CountSort( Inf , 1 );
		CountSort( Inf , 0 );

		Buc[Stp][Inf[0].Ind] = 0;
		for( j=1;j<Inf.size();j++){
			if( Inf[j].Val[0]==Inf[j-1].Val[0] && Inf[j].Val[1]==Inf[j-1].Val[1] ) Buc[Stp][Inf[j].Ind] = Buc[Stp][Inf[j-1].Ind];
			else{
				Buc[Stp][Inf[j].Ind] = j;
			}
		}
	}
	Stp--;
	for( i=0;i<Len;i++){
		SuffA[Buc[Stp][i]] = i;
	}

}
// fining longest common pref lenth of two suffix in O(longn)
long Lcp( long i,long j )
{
	long k,Tot = 0;
	for( k=Stp;k>=0 && i<Len && j<Len;k--){
		if( Buc[k][i]==Buc[k][j] ){
			i += 1<<k;
			j += 1<<k;
			Tot += 1<<k;
		}
	}
	return Tot;
}

int main( void )
{
	scanf("%s",Str );
	Len = strlen(Str);
	SuffArrCreate();
	return 0;
}

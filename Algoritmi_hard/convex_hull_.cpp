/*
	Algorithm: Convex Hull( Graham Scan )
	Compexity: O( nlogn )
*/

#include<stdio.h>
#include<string.h>
#include<vector>
#include<algorithm>
using namespace std;

struct POINT{
	long x,y;
	long I;
	POINT( long x = 0,long y = 0,long I = 0 )
	{
		this->x = x;
		this->y = y;
		this->I = I;
	}
};
bool operator<( const POINT &a,const POINT &b )
{
	if( a.x != b.x ) return a.x < b.x;
	else return a.y < b.y;
}

long N;
vector<POINT> Pt;
vector<bool> Del;

long Area2( POINT a,POINT b,POINT c )
{
	return a.x*b.y - a.y*b.x + a.y*c.x - a.x*c.y + b.x*c.y - c.x*b.y;
}

bool Cmp( const POINT &a,const POINT &b )
{
	long Ar = Area2( Pt[0],a,b );
	if( Ar ) return Ar>0;
	long Dx = labs( Pt[0].x-a.x ) - labs( Pt[0].x-b.x );
	long Dy = labs( Pt[0].y-a.y ) - labs( Pt[0].y-b.y );
	if( Dx<0 || Dy<0 ){
		Del[a.I] = true;
		return true;
	}
	else if( Dx>0 || Dy>0 ){
		Del[b.I] = true;
		return false;
	}
	return true;
}
/* to del all linear point */
void Squash( void )
{
	long i,j;
	for( i=j=0;i<N;i++){
		if( Del[Pt[i].I] ) continue;
		Pt[j] = Pt[i];
		j++;
	}
	Pt.resize( j );
	N = j;
}
void ConvexHull( vector<POINT> &Hull )
{
	sort( Pt.begin(),Pt.end()); // Pt[0] wiil contain leftmst-lowest point 
	sort( Pt.begin()+1,Pt.end(),Cmp );
	Squash();
	if( Pt.size()>=1 ) Hull.push_back( Pt[0] );
	if( Pt.size()>=2 ) Hull.push_back( Pt[1] );
	long i = 2;
	while( i<N ){
		long s = Hull.size();
		if( Area2( Hull[s-2],Hull[s-1],Pt[i] )>0 ){
			Hull.push_back( Pt[i] );
			i++;
		}
		else Hull.pop_back();
	}
}

int main( void )
{
	//freopen("text1.txt","r",stdin );
	Pt.resize(N); // repeated point is needed to omit by set
	Del.resize(N);
	vector<POINT> Hull;
	ConvexHull( Hull );
	return 0;
}

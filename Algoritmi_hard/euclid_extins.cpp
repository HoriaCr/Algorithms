#include<stdio.h>
#include<string.h>
#include<math.h>

typedef long long Long;

struct TRIPLE{
	Long x,y;
	Long d;
	TRIPLE( Long x=0,Long y=0,Long d=0 )
	{
		this->x = x;
		this->y = y;
		this->d = d;
	}
};


Long A,B,C;
Long X1,X2;
Long Y1,Y2;

TRIPLE EGCD( Long a,Long b )
{
	if( !b ){
		return TRIPLE( 1,0,a );
	}
	else{
		TRIPLE t = EGCD( b,a%b );
		return TRIPLE( t.y,t.x - (a/b)*t.y ,t.d );
	}
}

inline Long max( Long a,Long b )
{
	return a>b ? a:b;
}

inline Long min( Long a,Long b )
{
	return a<b ? a:b;
}


void FindAns( void )
{
	TRIPLE t = EGCD( A,B );
	if( C%t.d ){
		printf("%lld\n",0 );
		return;
	}
	A /= t.d;
	B /= t.d;
	C /= t.d;
	Long X = t.x*C;
	Long Y = t.y*C;

	X1 -= X;
	X2 -= X;
	Long LX,RX;
	if( B > 0 ){
		LX = (Long)ceil( (double)X1/B );
		RX = (Long)floor( (double)X2/B );
	}
	else{
		LX = (Long)ceil( (double)X2/B );
		RX = (Long)floor( (double)X1/B );
	}

	Y1 -= Y;
	Y2 -= Y;
	Long LY,RY;
	A *= -1;
	if( A > 0 ){
		LY = (Long)ceil( (double)Y1/A );
		RY = (Long)floor( (double)Y2/A );
	}
	else{
		LY = (Long)ceil( (double)Y2/A );
		RY = (Long)floor( (double)Y1/A );
	}

	Long Ans = min( RX,RY) - max( LX,LY ) + 1 ;
	Ans = max( Ans,0 );
	printf("%lld\n",Ans );
}



int main( void )
{
	Long i,Icase,k = 0;

	//freopen("text1.txt","r",stdin );

	scanf("%lld",&Icase );
	while( Icase-- ){
		scanf("%lld%lld%lld",&A,&B,&C );
		C *= -1;
		scanf("%lld%lld%lld%lld",&X1,&X2,&Y1,&Y2 );
		printf("Case %lld: ",++k );
		if( !A && !B ){
			if( C ) printf("0\n" );
			else printf("%lld\n",(X2-X1+1)*(Y2-Y1+1));
		}
		else if( !A && B ){
			if( C%B ) printf("0\n");
			else if( Y1<=(C/B) && (C/B)<=Y2 ) printf("%lld\n",(X2-X1+1) );
			else printf("0\n");
		}
		else if( A && !B ){
			if( C%A ) printf("0\n");
			else if( X1<=(C/A) && (C/A)<=X2 ) printf("%lld\n",(Y2-Y1+1) );
			else printf("0\n");
		}
		else FindAns();
	}

	return 0;
}

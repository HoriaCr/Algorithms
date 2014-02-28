#include <cstdio>



/*
*
*
    if n < 1,373,653, it is enough to test a = 2 and 3;
    if n < 9,080,191, it is enough to test a = 31 and 73;
    if n < 4,759,123,141, it is enough to test a = 2, 7, and 61;
    if n < 1,122,004,669,633, it is enough to test a = 2, 13, 23, and 1662803;
    if n < 2,152,302,898,747, it is enough to test a = 2, 3, 5, 7, and 11;
    if n < 3,474,749,660,383, it is enough to test a = 2, 3, 5, 7, 11, and 13;
    if n < 341,550,071,728,321, it is enough to test a = 2, 3, 5, 7, 11, 13, and 17.

*
*
*/
int expo(int a,int b,int mod)
{
	int res;
	for(res = 1;b;b>>=1,a = 1ll*a*a%mod) if(b & 1) res = 1ll*res*a%mod;
	return res;
}

bool miller_rabin(int a, int n) 
{
	int d = n - 1 , s = 0;
	while((d & 1) == 0) d>>=1 , s++;
	a = expo(a, d, n);
    if (a == 1) return 1;
    for (;s > 1;s--) {
        if (a == n-1) return 1;
        a = expo(a, 2, n);
    }
    if (a == n-1) return 1;
    return 0;
}

bool test(int n) 
{
    if(n == 2) return 1;
	if(n <= 1 || n%2 == 0) return 0;
	return (miller_rabin( 2, n) && (n <= 7  || miller_rabin( 7, n)) && (n <= 61 || miller_rabin(61, n)));
}

int main()
{
	freopen("test.in","r",stdin);
	freopen("test.out","w",stdout);
	int N , ans = 0, k;
	for(scanf("%d",&N);N;N--)
	{
		scanf("%d",&k);
		ans+=test(k);
	}
	printf("%d\n",ans);
	return 0;
}


template<class T> T MAX(T a,T b) { return a > b ? a : b;}

template<class T> int lcs(T a,T b) {
   int N = sizeof(a);
   int M = sizeof(b);
   int *dp[2]; 
   dp[0] = new int[M + 2];
   dp[1] = new int[M + 2];
   for(int i = 0;i <= M;++i) 
	dp[0][i] = dp[1][i] = 0;

   int L = 0;
   for(int i = 1;i <= N;++i,L = 1 - L)
	for(int j = 1;j <= M;++j)
	  dp[1 - L][j] = MAX(MAX(dp[1 - L][j - 1],dp[L][j]),dp[L][j - 1] + (a[i - 1] == b[j - 1] ? 1 : 0) );

  return dp[L][M];
}

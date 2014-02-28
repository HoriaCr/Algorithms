
const int NMAX = 1<<10;
int D[NMAX][NMAX];

void roy_floyd(const int &n)
{
	for(int k = 1;k <= n;++k)
		for(int i = 1;i <= n ;++i)
			for(int j = 1;j <= n;++j)
				if(i != j && D[i][k] && D[k][j] && (D[i][j] > D[i][k] + D[k][j] || !D[i][j]))
					D[i][j] = D[i][k] + D[k][j];
}

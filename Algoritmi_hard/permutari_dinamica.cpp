#include <fstream>

using namespace std;

ifstream fin("go.in");
ofstream fout("go.out");

int ans[64];

void getpermutation(int N,int rang)
{
	int i , fact , P[50] , O[50] , lg = 0;
	O[0] = 0 , fact = 1;
	for(i = 1; i <= N;++i)
		fact*= P[i] = i , O[i] = 0;
	
	for(int n = N , ic , k;rang != 1;n--)
	{
		fact/=n; 
		k = (rang - 1)/fact + 1;
		rang-=fact*(k - 1);
		for(ic = 0 , i = 1;ic < k;i++)
			if(O[i] == 0) ic++;

		O[i - 1] = 1;
		ans[++lg] = P[i - 1];
		//fout<<P[i - 1]<<" ";
	}
	for(int i = 1;i <= N;++i)
			if(O[i] == 0) ans[++lg] = i; 
				//fout<<i<<" ";
	//fout<<'\n';
}

void getpermutation_rang(int N,int p[])
{
	int fact = 1 , ans = 0;
	for(int i = 2;i <= N;++i) fact*=i;
	for(int i = 1;i < N;++i)
	{
		fact/=(N - i + 1);
		ans+=(p[i] - 1)*fact;
		for(int j = i + 1;j <= N;++j)
			if(p[j] > p[i]) p[j]--;
	}
	fout<<ans<<'\n';
}

int main()
{
	for(int i = 1; i <= 10;++i)
		getpermutation(5,i) , getpermutation_rang(5,ans);
	return 0;
}

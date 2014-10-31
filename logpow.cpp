#include <fstream>

using namespace std;

int logpow(int a,int b, const int& mod) {
	if (b < 0) {
		return logpow( logpow(a, -b), mod - 2, mod);
	}
	if (mod > 2 && b >= mod - 1) {
		b %= (mod - 1);
	}
	int ret = 1;
	while (b > 0) {
		if (b & 1) {
			ret = 1LL * ret * a % mod;
		}
		a = 1LL * a * a % mod;
		b >>= 1;
	}
	
	return ret;
}

int main() {
	ifstream cin("lgput.in");
	ofstream cout("lgput.out");
	const int mod = 1999999973;
	int N, P;
	cin >> N >> P;
	cout << logpow(N, P, mod);
	return 0;
}
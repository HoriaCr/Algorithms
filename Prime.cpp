#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

class Primes {

	public:
	static bool isPrime(const int& x) {
		if (x == 1) return 0;
		if (x == 2) return 1;
		if (!(x & 1)) return 0;
		for (int i = 3;i * i <= x; i += 2) {
			if (x % i == 0) {
				return false;
			}
		}
		
		return true;
	}

	static vector< pair<int,int> > getPrimeFactors(int x) {
		vector< pair<int,int> > ret;
		
		if (!(x & 1)) {
			ret.push_back( make_pair(2,0));
			do {
				x >>= 1;
				ret.back().second++;
			} while (!(x & 1));
		}

		int c;

		for (int d = 3;d * d <= x; d += 2) {
			if (x == (c = x / d) * d) {
				ret.push_back(make_pair(d,1));
				while ( (c = x / d) * d == x) {
					x = c;
					ret.back().second++;
				}
			}
		}

		if (x > 1) {
			ret.push_back(make_pair(x,1));
		}

		return ret;
	}

	static vector<int> generatePrimes(const int& n) {
		vector<int> ret;
		char *p = new char[(n >> 4) + 1]();
		for (int i = 1;((i * i) << 1) + (i << 1) <= n;i++) {
			if ( ((p[i >> 3] >> (i & 7)) & 1) == 0) {
				for (int j = ((i * i) << 1) + (i << 1); (j << 1) + 1 <= n;j += (i << 1) + 1) {
					p[j >> 3] |= 1 << (j & 7);
				}
			}
		}

		ret.push_back(2);
		for (int i = 1;(i << 1) + 1 <= n;i++) { 
			if ((p[i >> 3] & (1 << (i & 7))) == 0) {
				ret.push_back( (i << 1) + 1);
			}
		}

		delete[] p;

		return ret;
	}

};


class Combinatorics {
	
	vector<int> factorial;
	vector<int> inverseFactorial;
	int modulo;
	int n;
	

	public:

	int pow(int a,int b) {
		int ret = 1;
		for (;b;b >>= 1, a = 1LL * a * a % modulo) {
			if (b & 1) {
				ret = 1LL * a * a % modulo;
			}
		}

		return ret;
	}
	
	Combinatorics(const int& n_,const int& mod) {
		n = n_;
		modulo = mod;
		factorial.resize(n);
		factorial[0] = 1;
		inverseFactorial = factorial;
		for (int i = 1; i <= n; i++) {
			factorial[i] = 1LL * factorial[i - 1] * i % modulo;
		}

		inverseFactorial[n] = pow(factorial[n],modulo - 2);
		for (int i = n - 1;i >= 1;i--) {
			inverseFactorial[i] = 1LL * inverseFactorial[i + 1] * (i + 1) % modulo;
		}

	}


	int binomialCoefficient(int n,int k) {
		return 1LL * (1LL * factorial[n] * inverseFactorial[n - k] % modulo) * inverseFactorial[k] % modulo;
	}

	int P(int n,int k) {
		return 1LL * factorial[n] * inverseFactorial[n - k] % modulo;
	}
	


};

int main()
{

	return 0;
}


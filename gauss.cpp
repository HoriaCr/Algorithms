#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <queue>
#include <sstream>
#include <tuple>
#include <iostream>
#include <map>
#include <bitset>
#include <cstring>

using namespace std;

template<class DataType>
class Rational {
	DataType p, q;

		long long gcd(long long a, long long b) {
			return !b ? a : gcd(b, a % b);
		}

		long long lcm(long long a, long long b) {
			return a * b / gcd(a, b);
		}

		void simplify() {
			DataType g = gcd(p, q);
			p /= g;
			q /= g;
		}

	public:

		Rational(DataType P = DataType(), DataType Q = DataType(1)) :
			p(P), q(Q) {
		
			}

		Rational& operator = (const Rational& b) {
			p = b.p;
			q = b.q;
			return *this;
		}

		Rational& operator = (const int& b) {
			p = b;
			q = 1;
			return *this;
		}

		void operator += (const Rational& b) {
			if (b.q == DataType()) return;
			if (q == 0) {
				*this = b;
				return;
			}
			long long y = lcm(q, b.q);
			long long x = p * y / q + b.p * y / b.q;
			p = static_cast<DataType>(x);
			q = static_cast<DataType>(y);
			if (p == 0 || q == 0) {
				p = q = 0;
			}
		}

		Rational operator + (const Rational& b) const {
			Rational ret = *this;
			ret += b;
			return ret;
		}

		void operator -= (const Rational& b) {
			if (b.q == DataType()) return;
			if (q == 0) {
				*this = b;
				return;
			}
			long long y = lcm(q, b.q);
			long long x = p * y / q - b.p * y / b.q;
			p = static_cast<DataType>(x);
			q = static_cast<DataType>(y);
			if (p == 0 || q == 0) {
				p = q = 0;
			}
		}

		Rational operator - (const Rational& b) const {
			Rational ret = *this;
			ret -= b;
			return ret;
		}

		void operator *= (const Rational& b) {
			p *= b.p;
			q *= b.q;
			if (p == 0 || q == 0) {
				p = q = 0;
				return;
			}
		}

		Rational operator * (const Rational& b) const {
			Rational ret = *this;
			ret *= b;
			return ret;
		}

		void operator /= (const Rational& b) {
			p *= b.q;
			q *= b.p;
			if (p == 0 || q == 0) {
				p = q = 0;
				return;
			}
		}

		Rational operator / (const Rational& b) const {
			Rational ret = *this;
			ret /= b;
			return ret;
		}

		friend istream& operator >> (istream& in, Rational& f) {
			in >> f.p >> f.q;
			return in;
		}

		friend ostream& operator << (ostream& out,Rational& f) {
			if (!f.q) {
				out << 0;
			} else {
				out << 1.0 * f.p / f.q;
			}
			return out;
		}

		bool operator == (const Rational& b) {
			return (p == b.p && q == b.q);
		}

		bool operator < (const Rational& b) {
			return 1.0 * p / q < 1.0* b.p / b.q;
		}

		bool operator > (const Rational& b) {
			return 1.0 * p / q > 1.0* b.p / b.q;
		}
};

template<class DataType>
void printMatrix(vector< vector<DataType> >& A, ostream& out) {
	for (auto& x : A) {
		for (auto& y : x) {
			out << y << " ";
		}
		out << "\n";
	}
	out << "\n";
}

template<class DataType>
void gauss(vector< vector<DataType> >& A,ostream& out) {
	size_t n = A.size();
	if (!n) {
		return;
	}

	DataType det = 1;
	size_t m = A[0].size();
	printMatrix(A, out);
	for (size_t i = 0, j = 0; i < n  && j < m;) {
		size_t k = i;
		while (k < n && A[k][j] == DataType()) {
			k++;
		}
		
		if (k == n) {
			j++;
			continue;
		}

		if (i != k) {
			swap(A[i], A[k]);
			det *= (-1);
		}

		for (size_t w = j + 1; w < m; w++) {
			A[i][w] /= A[i][j];
		}

		det *= A[i][j];
		A[i][j] = 1;
		for (k = i + 1; k < n; k++) {
			for (size_t w = j + 1; w < m; w++) {
				A[k][w] -= A[i][w] * A[k][j];
			}
			A[k][j] = 0;
		}
		printMatrix(A, out);
		i++;
		j++;
	}

	out << "det : "<< det << "\n";
}

int main() {
	ifstream cin("test.in");
	ofstream cout("test.out");
	
	vector< vector< Rational<int> > > A = {
		{ { 1, 1 }, { 3, 1 }, { 1, 1 }, { 9, 1 } },
		{ { 1, 1 }, { 1, 1 }, { -1, 1 }, { 1, 1 } },
		{ { 3, 1 }, { 11, 1 }, { 5, 1 }, { 35, 1 } },
	};

	//gauss(A,cout);
	

	vector< vector< Rational<int> > > Z = {
		{ { 1, 1 }, { 2, 1 }, { 3, 1 } },
		{ { 4, 1 }, { 5, 1 }, { 6, 1 } },
		{ { 4, 1 }, { 2, 1 }, { 3, 1 } },
	};

	gauss(Z, cout);
	vector< vector< Rational<int> > > V = {
		{ { 0, 1 }, { 1, 1 }, { 3, 1 } },
		{ { 1, 1 }, { 2, 1 }, { 0, 1 } },
		{ { 0, 1 }, { 3, 1 }, { 4, 1 } },
	};
	
	gauss(V, cout);

	vector< vector< Rational<int> > > W = {
		{ { 4, 1 }, { 6, 1 } },
		{ { 3, 1 }, { 8, 1 } }	
	};
	gauss(W, cout);
	return 0;
}
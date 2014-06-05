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
#include <complex>

using namespace std;

template<class DataType>
class Polynomial {

	vector< DataType > data;

	void fft(vector< complex<double> >& A, const int& direction, const int NMAX = 1 << 16) {
		const long double PI = 3.14159265358979323846;
		for (int i = 1, j = 0; i < NMAX; i++) {
			int b = (NMAX >> 1);
			for (; j >= b; b >>= 1) j -= b;
			j += b;
			if (i < j) swap(A[i], A[j]);
		}

		for (int len = 2; len <= NMAX; len <<= 1) {
			long double theta = 2 * PI / len * direction;
			complex<double> z(cos(theta), sin(theta));
			for (int i = 0; i < NMAX; i += len) {
				complex<double> x(1, 0);
				for (int j = 0; j < (len >> 1); j++) {
					complex<double> u = A[i + j];
					complex<double> v = A[i + j + (len >> 1)] * x;
					A[i + j] = u + v;
					A[i + j + (len >> 1)] = u - v;
					x = x * z;
				}
			}
		}
	}

	Polynomial fastMuliply(const Polynomial& B) {
		const int NMAX = 1 << 15;
		Polynomial ret(data.size() + B.data().size());
		vector< complex<double> > a(NMAX);
		vector< complex<double> > b(NMAX);

		for (int i = 0; i < data.size(); i++) {
			a[i] = complex<double>(data[i], 0);
		}

		for (int i = 0; i < B.data.size(); i++) {
			b[i] = complex<double>(B.data[i], 0);
		}

		fft(a, 1);
		fft(b, 1);

		for (int i = 0; i < NMAX; i++) {
			a[i] = a[i] * b[i];
		}

		fft(a, -1);

		complex<double> y = complex<double>(1.0 / NMAX, 0);
		for (int i = 0; i < NMAX; i++) {
			a[i] = a[i] * y;
			ret[i] = static_cast<DataType>(round(a[i].real()));

		}

		return ret;
	}

public:

	Polynomial() {
	}

	Polynomial(const unsigned int& Size) {
		data.resize(Size);
	}

	Polynomial(const Polynomial& p) {
		data = p.data;
	}

	Polynomial(const vector<DataType>& p) {
		data = p;
	}

	Polynomial& operator = (const Polynomial& p) {
		data = p.data;
		return *this;
	}

	DataType& operator[] (const unsigned int& index) {
		try {
			if (index < 0 || index >= data.size()) throw out_of_range("Index  out of bounds");
		}
		catch (out_of_range& ex) {
			cout << ex.what() << "\n";
		}

		return data[index];
	}

	void operator += (const Polynomial& b) {
		if (b.data.size() > data.size()) {
			data.resize(b.data.size());
		}

		for (unsigned int i = 0; i < data.size() && i < b.data.size(); i++) {
			data[i] += b.data[i];
		}
	}

	Polynomial operator + (const Polynomial& b) const {
		Polynomial ret = *this;
		ret += b;
		return ret;
	}

	void operator -= (const Polynomial& b) {
		if (b.data.size() > data.size()) {
			data.resize(b.data.size());
		}

		for (unsigned int i = 0; i < data.size() && i < b.data.size(); i++) {
			data[i] -= b.data[i];
		}
	}

	Polynomial operator - (const Polynomial& b) const {
		Polynomial ret = *this;
		ret -= b;
		return ret;
	}

	void operator *= (const Polynomial& b) {

		Polynomial ret(data.size() + b.data.size());
		for (unsigned int i = 0; i < data.size(); i++) {
			for (unsigned int j = 0; j < b.data.size(); j++) {
				ret[i + j] += data[i] * b.data[j];
			}
		}

		data = ret.data;
	}

	Polynomial operator * (const Polynomial& b) const {
		Polynomial ret = *this;
		ret *= b;
		return ret;
	}

	void operator / (const Polynomial& b) {

	}

	friend ostream& operator << (ostream& out, Polynomial& p) {
		for (size_t i = 0; i < p.data.size(); i++) {
			if (p.data[i] != 0) {
				out << p.data[i];
				if (i) out << "x^" << i;
				if (i < p.data.size() - 1) out << " + ";
			}
		}
		return out;
	}
};

template<class DataType = int>
class SquareMatrix {
	unsigned int Size;
	vector< vector<DataType> > matrix;

	public:

	size_t size() {
		return Size;
	}

	SquareMatrix(unsigned int Size_ = 0) {
		Size = Size_;
		matrix.resize(Size, vector<DataType>(Size, DataType()));
	}

	SquareMatrix(const vector< vector<DataType> >& data) {
		if (data.size() == data[0].size()) {
			Size = data.size();
			matrix = data;
		}
	}

	vector<DataType>& operator [](const int& x) {
		return matrix[x];
	}

	SquareMatrix getIdentityMatrix() {
		SquareMatrix ret(matrix.Size);
		for (size_t i = 0; i < Size; i++) {
			for (size_t j = 0; j < Size; j++) {
				ret[i][j] = DataType();
			}
		}

		for (size_t i = 0; i < Size; i++) {
			ret[i][i] = static_cast<DataType>(1);
		}

		return ret;
	}

	SquareMatrix& operator = (SquareMatrix B) {
		matrix.resize(B.Size);
		for (size_t i = 0; i < Size; i++) {
			for (size_t j = 0; j < Size; j++) {
				matrix[i][j] = B.matrix[i][j];
			}
		}
		return *this;
	}

	void operator += (const SquareMatrix &B) {
		for (size_t i = 0; i < Size; i++) {
			for (size_t j = 0; j < Size; j++) {
				matrix[i][j] += B.matrix[i][j];
			}
		}
	}

	SquareMatrix operator + (const SquareMatrix &B) const {
		SquareMatrix ret(matrix);
		ret += B;
		return ret;
	}

	SquareMatrix operator * (const SquareMatrix &B) const {
		SquareMatrix ret(matrix.Size);
		for (size_t i = 0; i < Size; i++) {
			for (size_t j = 0; j < Size; j++) {
				for (size_t k = 0; k < Size; k++) {
					ret[i][j] += matrix[i][k] * B.matrix[k][j];
				}
			}
		}
		return ret;
	}

	SquareMatrix operator * (const DataType &scalar) const {
		SquareMatrix ret(matrix.Size);
		for (size_t i = 0; i < Size; i++) {
			for (size_t j = 0; j < Size; j++) {
				ret[i][j] = matrix[i][j] * scalar;
			}
		}
		return ret;
	}


	void operator *= (const DataType &scalar) {
		for (size_t i = 0; i < Size; i++) {
			for (size_t j = 0; j < Size; j++) {
				matrix[i][j] *= scalar;
			}
		}
	}


	void operator *= (const SquareMatrix &B) {
		matrix = matrix*B;
	}

	bool operator == (const SquareMatrix& B) const {
		if (Size != B.Size) return false;
		for (size_t i = 0; i < Size; i++) {
			for (size_t j = 0; j < Size; j++) {
				if (matrix[i][j] != B.matrix[i][j])
					return false;
			}
		}
		return true;
	}

	bool operator != (const SquareMatrix& B) const {
		return !(*this == B);
	}

	template<class powType>
	SquareMatrix pow(powType K) {
		SquareMatrix B = *this;
		SquareMatrix ret = B.getIdentityMatrix();

		for (; K > 0; K >>= 1) {
			if (K & 1)
				ret = ret * B;

			B = B * B;
		}

		return ret;
	}

	DataType getDeterminant() {
		size_t n = Size;
		size_t m = n;
		vector< vector< DataType > > A = matrix;
		DataType det = 1;
		for (size_t i = 0, j = 0; i < n && j < m;) {
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
				det *= static_cast<DataType>(-1);
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

			i++;
			j++;
		}

		return det;
	}

	friend istream& operator >> (istream& in, SquareMatrix& A) {
		in >> A.Size;
		A.matrix.resize(A.Size, vector<DataType>(A.Size));
		for (size_t i = 0; i < A.Size; i++) {
			for (size_t j = 0; j < A.Size; j++) {
				in >> A.matrix[i][j];
			}
		}
		return in;
	}

	friend ostream& operator << (ostream& out, SquareMatrix& A) {
		for (size_t i = 0; i < A.Size; i++) {
			for (size_t j = 0; j < A.Size; j++) {
				out << A.matrix[i][j] << " ";
			}
			out << "\n";
		}
		return out;
	}

};

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

	friend ostream& operator << (ostream& out, Rational& f) {
		out << f.p << "/" << f.q;
		return out;
	}

	bool operator == (const Rational& b) {
		return (p == b.p && q == b.q);
	}

	bool operator < (const Rational& b) {
		return 1.0 * p / q < 1.0* b.p / b.q;
	}

	bool operator >(const Rational& b) {
		return 1.0 * p / q > 1.0* b.p / b.q;
	}
};

int main() {
	ifstream cin("test.in");
	ofstream cout("test.out");

	SquareMatrix< Rational< Polynomial<int> > > A = vector< vector< Rational< Polynomial<int> > > > {
		
		
		{ Rational< Polynomial<int> >{ vector<int> { 1, 2, 3 }, vector<int> { 2, 4, 5 } }, Rational< Polynomial<int> >{ vector<int> { 1, 2, 3 }, vector<int> { 2, 4, 5 } } },
		{ Rational< Polynomial<int> >{ vector<int> { 1, 2, 3 }, vector<int> { 2, 4, 5 } }, Rational< Polynomial<int> >{ vector<int> { 1, 2, 3 }, vector<int> { 2, 4, 5 } } }

	};

	cout << A << "\n";
	
	return 0;
}
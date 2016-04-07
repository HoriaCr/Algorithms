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
#include <cassert>

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

	friend ostream& operator << (ostream& out, const Rational& f) {
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

	friend ostream& operator << (ostream& out, const Polynomial& p) {
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
	const static unsigned int minMatrixSize = 64;
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
        for (const auto& x : data) {
            assert(x.size() == data.size());
        }
        Size = data.size();
        matrix = data;
	}

	vector<DataType>& operator [](const int& x) {
		return matrix[x];
	}

	SquareMatrix getIdentityMatrix() const {
		SquareMatrix ret(Size);
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

	SquareMatrix& operator = (const SquareMatrix& B) {
		Size = B.Size;
		matrix = B.matrix;
        return *this;
	}

	void operator += (const SquareMatrix &B) {
		for (size_t i = 0; i < Size; i++) {
			for (size_t j = 0; j < Size; j++) {
				matrix[i][j] += B.matrix[i][j];
			}
		}
	}

    void operator -= (const SquareMatrix &B) {
		for (size_t i = 0; i < Size; i++) {
			for (size_t j = 0; j < Size; j++) {
				matrix[i][j] -= B.matrix[i][j];
			}
		}
    }

	SquareMatrix operator + (const SquareMatrix &B) const {
		SquareMatrix ret(matrix);
		ret += B;
		return ret;
	}

	SquareMatrix operator - (const SquareMatrix &B) const {
		SquareMatrix ret(matrix);
		ret -= B;
		return ret;
	}

	SquareMatrix operator * (const SquareMatrix &B) const {
		SquareMatrix ret(B.Size);
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

    SquareMatrix getTranspose() const { 
		SquareMatrix ret(Size);
		for (size_t i = 0; i < Size; i++) {
			for (size_t j = 0; j < Size; j++) {
                ret.matrix[i][j] = matrix[j][i];
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

    vector<DataType> multiplyByVector(const vector<DataType>& b) const {
       vector<DataType> ret(Size);
       for (size_t i = 0; i < Size; i++) {
            for (size_t j = 0; j < Size; j++) {
               ret[i] += b[j] * matrix[i][j]; 
            }
       }

       return ret;
    }


    SquareMatrix multiplyByConstant(const DataType& c) const {
        SquareMatrix ret(Size);
        for (size_t i = 0; i < Size; i++) {
            for (size_t j = 0; j < Size; j++) {
                ret.matrix[i][j] = matrix[i][j] * c;
            }
        }
        return ret;
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



    void strassenR(const SquareMatrix &A, 
                   const SquareMatrix &B, 
                         SquareMatrix &C) {
        if (A.Size <= minMatrixSize) {
            C = A * B;
            return;
        }   
        // other cases are treated here:
        else {
            int newSize = A.Size / 2;
            SquareMatrix
                a11(newSize), a12(newSize), a21(newSize), a22(newSize),
                b11(newSize), b12(newSize), b21(newSize), b22(newSize),
                c11(newSize), c12(newSize), c21(newSize), c22(newSize),
                p1(newSize), p2(newSize), p3(newSize), p4(newSize), 
                p5(newSize), p6(newSize), p7(newSize),
                aResult(newSize), bResult(newSize);
     
            int i, j;
     
            //dividing the matrices in 4 sub-matrices:
            for (i = 0; i < newSize; i++) {
                for (j = 0; j < newSize; j++) {
                    a11.matrix[i][j] = A.matrix[i][j];
                    a12.matrix[i][j] = A.matrix[i][j + newSize];
                    a21.matrix[i][j] = A.matrix[i + newSize][j];
                    a22.matrix[i][j] = A.matrix[i + newSize][j + newSize];
     
                    b11.matrix[i][j] = B.matrix[i][j];
                    b12.matrix[i][j] = B.matrix[i][j + newSize];
                    b21.matrix[i][j] = B.matrix[i + newSize][j];
                    b22.matrix[i][j] = B.matrix[i + newSize][j + newSize];
                }
            }
            // Calculating p1 to p7:
     
            aResult = a11 + a22;// a11 + a22
            bResult = b11 + b22; // b11 + b22
            strassenR(aResult, bResult, p1); // p1 = (a11+a22) * (b11+b22)
     
            aResult = a21 + a22; // a21 + a22
            strassenR(aResult, b11, p2); // p2 = (a21+a22) * (b11)
     
            bResult = b12 - b22; // b12 - b22
            strassenR(a11, bResult, p3); // p3 = (a11) * (b12 - b22)
     
            bResult = b21 - b11; // b21 - b11
            strassenR(a22, bResult, p4); // p4 = (a22) * (b21 - b11)
     
            aResult = a11 + a12; // a11 + a12
            strassenR(aResult, b22, p5); // p5 = (a11+a12) * (b22)   
     
            aResult = a21 - a11; // a21 - a11
            bResult = b11 + b12; // b11 + b12
            strassenR(aResult, bResult, p6); // p6 = (a21-a11) * (b11+b12)
     
            aResult = a12 - a22; // a12 - a22
            bResult  = b21 + b22; // b21 + b22
            strassenR(aResult, bResult, p7); // p7 = (a12-a22) * (b21+b22)
     
     
            c12 = p3 + p5; // c12 = p3 + p5
            c21 = p2 + p4; // c21 = p2 + p4
     
            aResult = p1 + p4; // p1 + p4
            bResult = aResult + p7; // p1 + p4 + p7
            c11 = bResult - p5; // c11 = p1 + p4 - p5 + p7
     
            aResult = p1 + p3; // p1 + p3
            bResult = aResult + p6; // p1 + p3 + p6
            c22 =  bResult - p2; // c22 = p1 + p3 - p2 + p6
    

            // Grouping the results obtained in a single matrix:
            for (i = 0; i < newSize ; i++) {
                for (j = 0 ; j < newSize ; j++) {
                    C.matrix[i][j] = c11.matrix[i][j];
                    C.matrix[i][j + newSize] = c12.matrix[i][j];
                    C.matrix[i + newSize][j] = c21.matrix[i][j];
                    C.matrix[i + newSize][j + newSize] = c22.matrix[i][j];
                }
            }
        }
    }


    unsigned int nextPowerOfTwo(int n) {
        return 1 << int(ceil(log2(n)));
    }


    void strassen(const SquareMatrix &A, 
                  const SquareMatrix &B, 
                  SquareMatrix &C) {
        unsigned int n = A.Size;
        unsigned int m = nextPowerOfTwo(n);
        SquareMatrix APrep(m), BPrep(m), CPrep(m);
        for(unsigned int i=0; i<n; i++) {
            for (unsigned int j=0; j<n; j++) {
                APrep.matrix[i][j] = A.matrix[i][j];
                BPrep.matrix[i][j] = B.matrix[i][j];
            }
        }

        strassenR(APrep, BPrep, CPrep);
        C = A;
        for(unsigned int i=0; i<n; i++) {
            for (unsigned int j=0; j<n; j++) {
                C.matrix[i][j] = CPrep.matrix[i][j];
            }
        }
    }

    DataType getNormInf() const {
        DataType ret = numeric_limits<DataType>::min();
        for (size_t i = 0; i < Size; i++) {
           ret = max(ret, accumulate(begin(matrix[i]), end(matrix[i]),
                       DataType(), [](const DataType& s, DataType b) {
               return s + abs(b);     
            }));
        }
        return ret;
    }

    DataType getNorm1() const {
        SquareMatrix m = getTranspose();
        return m.getNormInf();
    }

    SquareMatrix getInverse() const {
        
		size_t n = Size;
		size_t m = n;
        assert(n == m);
        SquareMatrix ret(Size);
        ret = ret.getIdentityMatrix();

		
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
                swap(ret.matrix[i], ret.matrix[k]);
				det *= static_cast<DataType>(-1);
			}

			for (size_t w = j + 1; w < m; w++) {
				A[i][w] /= A[i][j];
			}

            for (size_t w = 0; w < m; w++) {
                ret.matrix[i][w] /= A[i][j];
            }

			det *= A[i][j];
			A[i][j] = 1;
			for (k = i + 1; k < n; k++) { 
                for (size_t w = 0; w < m; w++) {
                    ret.matrix[k][w] -= ret.matrix[i][w]* A[k][j]; 
                }
                for (size_t w = j + 1; w < m; w++) {
					A[k][w] -= A[i][w] * A[k][j];
				}

				A[k][j] = 0;
			}

			i++;
			j++;
		}

        for (int i = (int)n - 2; i >= 0; i--) {
            for (int j = 0; j <= i; j++) {
                DataType val = A[j][i + 1];
                for (size_t k = 0; k < n; k++) {
                    ret.matrix[j][k] -= ret.matrix[i + 1][k] * A[j][i + 1]; 
                }

                for (size_t k = 0; k < n; k++) {
                    A[j][k] -= A[i + 1][k] * val;
                }
            }
        }
        return ret;
    }

	DataType getDeterminant() const {
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

    static vector< vector<DataType> > computeBinomial(uint32_t n, uint32_t k) {
        vector< vector<DataType> > C(n + 1, vector<DataType>(k + 1));
        C[0][0] = 1;
        for (size_t i = 1; i <= n; i++) {
            C[i][0] = 1;
            for (size_t j = 1; j <= k; j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];  
            }
        }
        return C;
    }

    pair<SquareMatrix, SquareMatrix> LU() const {
        SquareMatrix L(Size), U(Size);

        for (size_t i = 0; i < Size; i++) {
            L.matrix[i][0] = matrix[i][0];
        }
        U.matrix[0][0] = 1;

        for (size_t j = 1; j < Size; j++) {
            U.matrix[0][j] = matrix[0][j] / L.matrix[0][0];
        }

        for (size_t k = 1; k < Size; k++) {
            for (size_t i = k; i < Size; i++) {
                L.matrix[i][k] = matrix[i][k];
                for (size_t p = 0; p < k; p++) {
                    L.matrix[i][k] -= L.matrix[i][p] * U.matrix[p][k];
                }
            }
            U.matrix[k][k] = 1;
            for (size_t j = k + 1; j < Size; j++) {
                U.matrix[k][j] = matrix[k][j];
                for (size_t p = 0; p < k; p++) {
                    U.matrix[k][j] -= L.matrix[k][p] * U.matrix[p][j];
                }
                U.matrix[k][j] /= L.matrix[k][k];
            }
        }
        return make_pair(L, U); 
    }


    // Ax = b
    vector<DataType> solveWithLU(const vector<DataType>& b) const {
        vector<DataType> x(Size), y(Size);
        pair<SquareMatrix, SquareMatrix> lu = LU();
        for (size_t i = 0; i < Size; i++) {
           y[i] = b[i];
           for (size_t k = 0; k < i; k++) {
                y[i] -= lu.first.matrix[i][k] * y[k];
           }
           y[i] /= lu.first.matrix[i][i];
        }
        for (int i = Size - 1; i >= 0; i--) {
            x[i] = y[i];
            for (size_t k = (size_t)i + 1; k < Size; k++) {
                x[i] -= lu.second.matrix[i][k] * x[k];
            }
        }
        return x;
    }



    pair<SquareMatrix, SquareMatrix> cholesky() const {
        SquareMatrix L(Size), Lt(Size);

        for (size_t j = 0; j < Size; j++) {
            
            L.matrix[j][j] = matrix[j][j];
            for (size_t k = 0; k < j; k++) {
                L.matrix[j][j] -= L.matrix[j][k] * L.matrix[j][k];
            }
            L.matrix[j][j] = sqrt(L.matrix[j][j]);
            for (size_t i = j + 1; i < Size; i++) {
                L.matrix[i][j] = matrix[i][j];
                for (size_t k = 0; k < j; k++) {
                    L.matrix[i][j] -= L.matrix[i][k] * L.matrix[j][k];
                }
                L.matrix[i][j] /= L.matrix[j][j];
            }


        }
        Lt = L.getTranspose();
        return make_pair(L, Lt); 
    }

    vector<DataType> solveWithCholesky(const vector<DataType>& b) const {
        vector<DataType> x(Size), y(Size);
        pair<SquareMatrix, SquareMatrix> llt = cholesky();
        for (size_t i = 0; i < Size; i++) {
           y[i] = b[i];
           for (size_t k = 0; k < i; k++) {
                y[i] -= llt.first.matrix[i][k] * y[k];
           }
           y[i] /= llt.first.matrix[i][i];
        }
        for (int i = Size - 1; i >= 0; i--) {
            x[i] = y[i];
            for (size_t k = (size_t)i + 1; k < Size; k++) {
                x[i] -= llt.first.matrix[k][i] * x[k];
            }
            x[i] /= llt.first.matrix[i][i]; 
        }
        return x;
    }


    pair<SquareMatrix, SquareMatrix> QR() const {
        SquareMatrix Q(Size), R(Size);
        for (size_t i = 0; i < Size; i++) {
            R.matrix[0][0] += matrix[i][0] * matrix[i][0];

        }
        R.matrix[0][0] = sqrt(R.matrix[0][0]);

        for (size_t i = 0; i < Size; i++) {
            Q[i][0] = matrix[i][0] / R.matrix[0][0];
        }

        for (size_t k = 1; k < Size; k++) { 
            for (size_t j = 0; j < k; j++) { 
                for (size_t i = 0; i < Size; i++) {
                    R.matrix[j][k] += matrix[i][k] * Q.matrix[i][j];
                }
            }
            for (size_t j = 0; j < k; j++) { 
                R.matrix[k][k] -= R.matrix[j][k] * R.matrix[j][k];
            }
            for (size_t i = 0; i < Size; i++) {
                R.matrix[k][k] += matrix[i][k] * matrix[i][k];
            } 
            R.matrix[k][k] = sqrt(R.matrix[k][k]);
            for (size_t i = 0; i < Size; i++) {
                Q.matrix[i][k] = matrix[i][k];
                for (size_t j = 0; j < k; j++) {
                    Q.matrix[i][k] -= R.matrix[j][k] * Q.matrix[i][j];
                }
                Q.matrix[i][k] /= R.matrix[k][k];
            }
            
        }

        return make_pair(Q, R); 
    }

    vector<DataType> solveWithQR(const vector<DataType>& b) const {
        vector<DataType> x(Size), y(Size);
        pair<SquareMatrix, SquareMatrix> qr = QR();
        for (size_t i = 0; i < Size; i++) {
            for (size_t j = 0; j < Size; j++) {
                y[i] += qr.first.matrix[j][i] * b[j];
            }
        }
        x[Size - 1] = y[Size - 1] / qr.second.matrix[Size - 1][Size -1];
        for (int i = Size - 1; i >= 0; i--) {
            x[i] = y[i];
            for (size_t j = (size_t)i + 1; j < Size; j++) {
                x[i] -= qr.second.matrix[i][j] * x[j];
            }
            x[i] /= qr.second.matrix[i][i];
        }
        return x;
    }

    SquareMatrix getDiag() const {
        SquareMatrix ret(Size);
        for (size_t i = 0; i < Size; i++) {
            ret.matrix[i][i] = matrix[i][i]; 
        }
        return ret;
    }

    tuple<vector<DataType>, DataType, size_t> jacobiRelaxata(const 
            size_t& iterations, long double epsilon,  const vector<DataType>& a)
        const {
        SquareMatrix I = getIdentityMatrix();

        SquareMatrix B(Size);
        vector<DataType> b(Size);
        for (size_t i = 0; i < Size; i++) {
            for (size_t j = 0; j < Size; j++) {
                B.matrix[i][j] = matrix[i][j] / matrix[i][i];
            }
            b[i] = a[i] / matrix[i][i];
        }
        
        DataType t = 2.0 / B.getNormInf();
        DataType h = t / iterations;
        vector<DataType> ret;
        size_t minSteps = numeric_limits<size_t>::max();
        DataType sigmaMin = 0.0;
        for (size_t k = 1; k < iterations; k++) {
            DataType sigma = k * h;
            // I - sigma * D^(-1)*A
            SquareMatrix C = I - B.multiplyByConstant(sigma);
            vector<long double> c = b;
            for (size_t i = 0; i < Size; i++) {
                c[i] = sigma * b[i];
            }
            DataType err = 0.0;
            size_t steps = 0; 
            vector<DataType> x(Size), y;
            do {
                y = C.multiplyByVector(x);
                for (size_t i = 0; i < Size; i++) {
                    y[i] += c[i];
                }
                
                DataType norm = 0.0;
                for (size_t i = 0; i < Size; i++) {
                    norm += matrix[i][i] * abs(y[i] - x[i]) * abs(y[i] - x[i]);
                }

                norm = sqrt(norm);
                err = norm;
                swap(x, y);
                steps++;
            } while (err >= epsilon);
            
            if (steps < minSteps) {
                minSteps = steps;
                sigmaMin = sigma;
                ret = x;
            }
        }
        
        return make_tuple(ret, sigmaMin, minSteps);
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

	friend ostream& operator << (ostream& out, const SquareMatrix& A) {
		for (size_t i = 0; i < A.Size; i++) {
			for (size_t j = 0; j < A.Size; j++) {
				out << A.matrix[i][j] << " ";
			}
			out << "\n";
		}
		return out;
	}

    vector<DataType> metodaGradientuluiConjugat(long double epsilon,  const
            vector<DataType>& b) const {
        vector<DataType> x(Size);
        vector<DataType> r = b;
        vector<DataType> aux = (*this).multiplyByVector(x);
        auto getDot = [](const vector<DataType>& a, const vector<DataType>& b) {
            DataType ret = DataType();
            for (size_t i = 0; i < a.size(); i++) {
                ret += a[i] * b[i];
            }
            return ret;
        };

        for (size_t i = 0; i < Size; i++) {
            r[i] -= aux[i];
        }
        
        vector<DataType> p = r;
        DataType rsold = getDot(r, r);
        for (size_t i = 0; i < Size; i++) {
            vector<DataType> ap = (*this).multiplyByVector(p);
            DataType alpha = rsold / getDot(p, ap);
            for (size_t j = 0; j < Size; j++) {
                x[j] += alpha * p[j];
                r[j] -= alpha * ap[j];
            }
            DataType rsnew = getDot(r, r);
            if (sqrt(rsnew) < epsilon) {
                break;
            }
            DataType f = rsnew / rsold;
            for (size_t j = 0; j < Size; j++) {
                p[j] = r[j] + f * p[j]; 
            }
            rsold = rsnew;
        }
        return x;

    }

    tuple<vector<DataType>, DataType, size_t> GaussSeidelRelaxata(const 
            size_t& iterations, long double epsilon,  const vector<DataType>&
            b) const {
        
        DataType t = 2.0;
        DataType h = t / iterations;
        vector<DataType> ret;
        size_t minSteps = numeric_limits<size_t>::max();
        DataType sigmaMin = 0.0;
        for (size_t k = 1; k < iterations; k++) {
            DataType sigma = k * h;
            DataType err = 0.0;
            size_t steps = 0; 
            vector<DataType> x(Size), y(Size);
            do {
                for (size_t i = 0; i < Size; i++) {
                    DataType s = -b[i];
                    for (size_t j = 0; j < i; j++) {
                        s += matrix[i][j] * y[j]; 
                    }
                    for (size_t j = i + 1; j < Size; j++) {
                        s += matrix[i][j] * x[j]; 
                    }
                    y[i] = (1 - sigma) * x[i] - (s * sigma) / matrix[i][i];
                }
                
                DataType norm = 0.0;
                for (size_t i = 0; i < Size; i++) {
                    for (size_t j = 0; j < Size; j++) {
                        norm += matrix[i][j] * abs(y[i] - x[i]) * abs(y[j] -
                                x[j]);
                    }
                }

                norm = sqrt(norm);
                err = norm;
                swap(x, y);
                steps++;
            } while (err >= epsilon);
            
            if (steps < minSteps) {
                minSteps = steps;
                sigmaMin = sigma;
                ret = x;
            }
        }
        
        return make_tuple(ret, sigmaMin, minSteps);

    }
    static SquareMatrix getAp(uint32_t n, uint32_t p) {
        SquareMatrix ret(n);
        vector< vector<DataType> > C = computeBinomial(n + p, n + p);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                ret.matrix[i][j] = static_cast<DataType>(C[p + j][i]);
            }
        }
        return ret;
    }

    static pair<SquareMatrix, vector<DataType> > getAt3(uint32_t m) {
        SquareMatrix A(m);
        vector<DataType> b(m);
        for (size_t i = 0; i < m; i++) {
            b[i] = 1.0 / (m * m);
            A.matrix[i][i] = 2 + b[i];
            if (i != m - 1) {
                A.matrix[i][i + 1] = A.matrix[i + 1][i] = -1;
            }
        }
        return make_pair(A, b);
    }

};


template<class T> void printV(T v, const string sep = " ") {
    cout.precision(16);
    for (auto x : v) {
        cout << x << sep;
    }
}

void rezolvaGaussSeidel() {
    cout << "Metoda gauss seidel\n";
    size_t m = 10;
    long double eps = 1e-10;
    size_t intervale = 100;
    pair<SquareMatrix<long double>, vector<long double> > p = SquareMatrix<long
        double>::getAt3(m);

    vector<long double> x;
    long double sigma;
    size_t steps;
    tie(x, sigma, steps) = p.first.GaussSeidelRelaxata(intervale, eps, p.second);
    cout << "Solutia : \n";
    printV(x, "\n");
    cout << "\nSigma :" << sigma << "\n";
    cout << "Numar pasi: " << steps << "\n";
}

void rezolvaJacobi() {
    cout << "Metoda jacobi\n";
    size_t m = 10;
    long double eps = 1e-10;
    size_t intervale = 100;
    pair<SquareMatrix<long double>, vector<long double> > p = SquareMatrix<long
        double>::getAt3(m);

    /*
    cout << p.first << "\n";
    printV(p.second);
    cout << "\n";
    */
    vector<long double> x;
    long double sigma;
    size_t steps;
    tie(x, sigma, steps) = p.first.jacobiRelaxata(intervale, eps, p.second);
    cout << "Solutia : \n";
    printV(x, "\n");
    cout << "\nSigma :" << sigma << "\n";
    cout << "Numar pasi: " << steps << "\n\n\n";
}

void rezolvaMetGConj() {
    size_t m = 10;
    long double eps = 1e-10;
    pair<SquareMatrix<long double>, vector<long double> > p = SquareMatrix<long
        double>::getAt3(m);
    vector<long double> x = p.first.metodaGradientuluiConjugat(eps, p.second);

    cout << "\n\nMetoda Gradientului conjugat\n";
    printV(x, "\n");
}

int main() {
    rezolvaJacobi();
    rezolvaGaussSeidel();
    rezolvaMetGConj();
    return 0;
}

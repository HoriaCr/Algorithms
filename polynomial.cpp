#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <complex>
#include <stdexcept>

using namespace std;


template<class DataType>
class Poly {

        vector< DataType > data;

        void fft(vector< complex<double> >& A,const int& direction,const int NMAX = 1 << 16) {
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

        Poly fastMuliply(const Poly& B) {
            const int NMAX = 1 << 15;
            Poly ret(data.size() + B.data().size());
            vector< complex<double> > a(NMAX);
            vector< complex<double> > b(NMAX);

            for (int i = 0; i < data.size(); i++) {
                a[i] = complex<double>(data[i],0);
            }

            for (int i = 0; i < B.data.size(); i++) {
                b[i] = complex<double>(B.data[i],0);
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

        Poly() {
        }

        Poly(const unsigned int& Size) {
            data.resize(Size);
        }

        Poly(const Poly& p) {
            data = p.data;
        }

        Poly& operator = (const Poly& p) {
            data = p.data;
            return *this;
        }

        DataType& operator[] (const unsigned int& index) {
            try {
                if (index < 0 || index >= data.size()) throw out_of_range("Index  out of bounds");
            }
            catch(out_of_range& ex) {
                cout << ex.what() << "\n";
            }

            return data[index];
        }

        void operator += (const Poly& b) {
            if (b.data.size() > data.size()) {
                data.resize(b.data.size());
            }

            for (unsigned int i = 0; i < data.size() && i < b.data.size(); i++) {
                data[i] += b.data[i];
            }
        }

        Poly operator + (const Poly& b) const {
            Poly ret = *this;
            ret += b;
            return ret;
        }

        void operator -= (const Poly& b) {
            if (b.data.size() > data.size()) {
                data.resize(b.data.size());
            }

            for (unsigned int i = 0; i < data.size() && i < b.data.size(); i++) {
                data[i] -= b.data[i];
            }
        }

        Poly operator - (const Poly& b) const {
            Poly ret = *this;
            ret -= b;
            return ret;
        }

        void operator *= (const Poly& b) {

            Poly ret(data.size() + b.data.size());
            for (unsigned int i = 0; i < data.size() ; i++) {
                for (unsigned int j = 0; j < b.data.size(); j++) {
                    ret[i + j] += data[i] * b.data[j];
                }
            }

            data = ret.data;
        }

        Poly operator * (const Poly& b) const {
            Poly ret = *this;
            ret *= b;
            return ret;
        }

        void operator / (const Poly& b) {

        }

        friend ostream& operator << (ostream& out,Poly& p) {
            for (int i = 0; i < (int)p.data.size(); i++) {
                out << p.data[i] << " ";
            }
            return out;
        }
};


int main()
{


    return 0;
}

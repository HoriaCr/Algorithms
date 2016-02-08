#include <iostream>
#include <cstring>
#include <unordered_set>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cstring>
 
using namespace std;
 
 
const int NMAX = 1 << 16;
const int root = 576953005, mod = 983826433;
//mod = NMAX * K +1, primitive root root^NMAX % mod = 1
 
inline int addMod(const int& x, const int& y) {
    return x + y < mod ? x + y : x + y - mod;
}
 
inline int subtractMod(const int& x, const int& y) {
    return x - y < 0 ? x - y + mod : x - y;
}
 
inline int expo(int a, int b, const int& mod) {
    int ret = 1;
    for (; b; b >>= 1, a = 1LL * a * a % mod) {
        if (b & 1) 
            ret = 1LL * ret * a % mod;
    }
    return ret;
}
 
inline void ntt(int* A, int direction) {
 
    for (int i = 1, j = 0; i < NMAX; i++) {
        int b = (NMAX >> 1);
        for (; j >= b; b >>= 1) j -= b;
        j += b;
        if (i < j) swap(A[i], A[j]);
    }
 
    int k = root;
    if (direction == -1) k = expo(root, mod - 2, mod); // root^(-1)
    for (int len = 2; len <= NMAX; len <<= 1) {
        int z = expo(k, NMAX / len, mod);
        for (int i = 0; i < NMAX; i += len) {
            int x = 1;
            for (int j = 0; j < (len >> 1); j++) {
                int u = A[i + j];
                int v = 1LL * A[i + j + (len >> 1)] * x % mod;
                A[i + j] = addMod(u, v);
                A[i + j + (len >> 1)] = subtractMod(u, v);
                x = 1LL * x * z % mod;
            }
        }
    }
 
    if (direction == -1) {
        int invN = expo(NMAX, mod - 2, mod);
        for (int i = 0; i < NMAX; i++) {
            A[i] = 1LL * A[i] * invN % mod;
            if (A[i]) cout << i << " " << A[i] << "\n";
            A[i] = 0;
        }
 
        cout << "\n";
    }
}
 
 
void print(int n, int* x, int m, int* y) {
    cout << n << "\n";
    for (int i = 0; i < n; i++) cout << x[i] << " ";
    cout << m << "\n";
    for (int j = 0; j < m; j++) cout << y[j] << " ";
    cout << "\n";
}
 
vector< int > brute(int n, int* x, int m, int* y) {
    vector<int> a(NMAX);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            a[x[i] + y[j]]++;
        }
    }
 
    return a;
}
 
int T, N, M;
int a[NMAX], b[NMAX];
int X[NMAX], Y[NMAX];
 
void generateData() {
    int t = 1000;
    srand(static_cast<unsigned int>(time(0)));
    cout << t << "\n";
    while (t--) {
        int n = rand() % 100 + 1;
        int m = rand() % 100 + 1;
        for (int i = 0; i < n; i++) X[i] = rand() % 30000 + 1;
        for (int i = 0; i < m; i++) Y[i] = rand() % 30000 + 1;
        print(n, X, m, Y);
    }
}
 
int gcd(int a, int b) { return !b ? a : gcd(b, a % b); }
 
void findPrimitive() {
    for (int i = 15000; i < 6000; i++) {
        int j = (i << 16) + 1;
        if (j % 2 == 0) continue;
        bool ok = 1;
        for (int k = 3; k * k <= j && ok; k += 2) if (j %k == 0) ok = 0;
        if (ok) {
            cout << j << "\n";
            return;
        }
    }
 
    int pr = 11;
    int z = 0;
    int y = pr;
    for (int i = 1; i < mod && z < 10; i++)  {
        if (gcd(y, mod - 1) == 1) {
            if (expo(y, NMAX, mod) == 1)  {
                cout << y << "\n";
                z++;
            }
        }
 
        y = 1LL*y * pr % mod;
    }
}
 
int main() {
 
    cin >> N;
    int x;
    for (int i = 0; i < N; i++) {
        cin >> x;
        a[x]++;
    }

    cin >> M;
    for (int i = 0; i < M; i++) {
        cin >> x;
        b[x]++;
    }

    ntt(a, 1);
    ntt(b, 1);

    for (int i = 0; i < NMAX; i++) {
        a[i] = 1LL * a[i] * b[i] % mod;
        b[i] = 0;
    }

    ntt(a, -1);

    return 0;
}

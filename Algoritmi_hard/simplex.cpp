//nya's library of simplex
//A*X <= B
//minimize C*X

#include <utility>
#include <limits>
#include <iostream>
#include <algorithm>
#include <vector>
#include <numeric>

using namespace std;

#define REP(i,n) for(int i = 0; i < (int)(n); i++) 
#define FOR(i,c) for(__typeof((c).begin()) i = (c).begin(); i != (c).end(); ++i) 
#define ALLOF(c) ((c).begin()), ((c).end()) 

const double EPS = 1.0e-10; 
const double INF = numeric_limits<double>::infinity(); 
  
typedef vector<double> vector_t; 
typedef vector<vector_t> matrix_t; 
  
vector_t simplex(matrix_t A, vector_t b, vector_t c) { 
  
    const int n = c.size(), m = b.size(); 
  
    // modify b to non-negative 
    REP(i, m) if (b[i] < 0) { 
        REP(j, n) 
            A[i][j] *= -1; 
        b[i] *= -1; 
    } 
  
    // list of base/independent variable ids 
    vector<int> bx(m), nx(n); 
    REP(i, m) 
        bx[i] = n+i; 
    REP(i, n) 
        nx[i] = i; 
  
    // extend A, b 
    A.resize(m+2); 
    REP(i, m+2) 
        A[i].resize(n+m, 0); 
    REP(i, m) 
        A[i][n+i] = 1; 
    REP(i, m) REP(j, n) 
        A[m][j] += A[i][j]; 
    b.push_back(accumulate(ALLOF(b), (double)0.0)); 
    REP(j, n) 
        A[m+1][j] = -c[j]; 
    REP(i, m) 
        A[m+1][n+i] = -INF; 
    b.push_back(0); 
  
    // main optimization 
    REP(phase, 2) { 
  
        for(;;) { 
  
            // select an independent variable 
            int ni = -1; 
            REP(i, n) 
                if (A[m][nx[i]] > EPS && (ni < 0 || nx[i] < nx[ni])) 
                    ni = i; 
            if (ni < 0) 
                break; 
  
            int nv = nx[ni]; 
  
            // select a base variable 
            vector_t bound(m); 
            REP(i, m) 
                bound[i] = (A[i][nv] < EPS ? INF : b[i] / A[i][nv]); 
            if (!(*min_element(ALLOF(bound)) < INF)) 
                return vector_t(); // -infinity 
  
            int bi = 0; 
            REP(i, m) 
                if (bound[i] < bound[bi]-EPS || (bound[i] < bound[bi]+EPS && bx[i] < bx[bi])) 
                    bi = i; 
  
            // pivot 
            double pd = A[bi][nv]; 
            REP(j, n+m) 
                A[bi][j] /= pd; 
            b[bi] /= pd; 
  
            REP(i, m+2) if (i != bi) { 
                double pn = A[i][nv]; 
                REP(j, n+m) 
                    A[i][j] -= A[bi][j] * pn; 
                b[i] -= b[bi] * pn; 
            } 
  
            swap(nx[ni], bx[bi]); 
        } 
  
        if (phase == 0 && abs(b[m]) > EPS) 
            return vector_t(); // no solution 
  
        A[m].swap(A[m+1]); 
        swap(b[m], b[m+1]); 
    } 
  
    vector_t x(n+m, 0); 
    REP(i, m) 
        x[bx[i]] = b[i]; 
    x.resize(n); 
  
    return x; 
}
/*
// Simon Lo's
// Simplex algorithm on augmented matrix a of dimension (m+1)x(n+1)
// returns 1 if feasible, 0 if not feasible, -1 if unbounded
// returns solution in b[] in original var order, max(f) in ret
// form: maximize sum_j(a_mj*x_j)-a_mn s.t. sum_j(a_ij*x_j)<=a_in
// in standard form.
// To convert into standard form:
// 1. if exists equality constraint, then replace by both >= and <=
// 2. if variable x doesn't have nonnegativity constraint, then replace by
// difference of 2 variables like x1-x2, where x1>=0, x2>=0
// 3. for a>=b constraints, convert to -a<=-b
// note: watch out for -0.0 in the solution, algorithm may cycle
// eps = 1e-7 may give wrong answer, 1e-10 is better
void pivot(int m, int n, double a[maxm][maxn], int B[maxm], int N[maxn], int r, int c) {
int i, j;
swap(N[c], B[r]);
a[r][c]=1/a[r][c];
for (j=0; j<=n; j++)if (j!=c) a[r][j]*=a[r][c];
for (i=0; i<=m; i++)if (i!=r) {
for (j=0; j<=n; j++)if (j!=c)
a[i][j]-=a[i][c]*a[r][j];
a[i][c] = -a[i][c]*a[r][c];
}
}
int feasible(int m, int n, double a[maxm][maxn], int B[maxm], int N[maxn]) {
int r, c, i; double p, v;
while (1) {
for (p=inf, i=0; i<m; i++) if (a[i][n]<p) p=a[r=i][n];
if (p>-eps) return 1;
for (p=0, i=0; i<n; i++) if (a[r][i]<p) p=a[r][c=i];
if (p>-eps) return 0;
p = a[r][n]/a[r][c];
for (i=r+1; i<m; i++) if (a[i][c]>eps) {
v = a[i][n]/a[i][c];
if (v<p) r=i, p=v;
}
pivot(m, n, a, B, N, r, c);
}
}
int simplex(int m, int n, double a[maxm][maxn], double b[maxn], double& ret) {
int B[maxm], N[maxn], r, c, i; double p, v;
for (i=0; i<n; i++) N[i]=i;
for (i=0; i<m; i++) B[i]=n+i;
if (!feasible(m, n, a, B, N)) return 0;
while (1) {
for (p=0, i=0; i<n; i++) if (a[m][i]>p)
p=a[m][c=i];
if (p<eps) {
for (i=0; i<n; i++) if (N[i]<n)
b[N[i]]=0;
for (i=0; i<m; i++) if (B[i]<n)
b[B[i]]=a[i][n];
ret = -a[m][n];
return 1;
}
for (p=inf, i=0; i<m; i++) if (a[i][c]>eps) {
v = a[i][n]/a[i][c];
if (v<p) p=v, r=i;
}
if (p==inf) return -1;
pivot(m, n, a, B, N, r, c);
}
}
*/

int main()
{

	return 0;
}
#include <fstream>
#include <vector>
#include <vector>
#include <limits>
#include <cstring>
#include <algorithm>
 
using namespace std;
 
ifstream cin("test.in");
ofstream cout("test.out");


int getSgn(vector<int> &a) {
	int ret = 0;
	for (size_t i = 0;i < a.size();i++) {
		for (size_t j = i + 1;j < a.size();j++) {
			if (a[i] > a[j]) {
				ret ^= 1;
			}
		}
	}
	
	return ret ? -1 : 1;
}

template<class T> 
T computeDeterminant(vector< vector<T> > A) {
	T ret = 0;
	vector<int> a(A.size());
	for (size_t i = 0;i < a.size();i++) {
		a[i] = i;
	}

	do {
		T prod = getSgn(a);
		for (size_t i = 0;i < a.size();i++) {
			prod *= A[i][a[i]];
		}

		ret += prod;

	} while (next_permutation(a.begin(),a.end()));
	return ret;
}



int main()
{
	int n;
	cin >> n;
	vector< vector< double> > A(n,vector<double>(n));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cin >> A[i][j];
		}
	}

	cout << fixed << computeDeterminant(A) << "\n";
    return 0;
}
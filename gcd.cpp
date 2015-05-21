#include <fstream>

using namespace std;

int gcd(int a,int b) {
	return !b ? a : gcd(b, a % b);
}

int main() {
	ifstream cin("gcd.in");
	ofstream cout("gcd.out");
	int T;
	for (cin >> T; T; T--) {
		int a, b;
		cin >> a >> b;
		cout << gcd(a, b) << "\n";
	}
	return 0;
}
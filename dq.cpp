#include <fstream>
#include <deque>
#include <vector>

using namespace std;

int main() {
	ifstream cin("deque.in");
	ofstream cout("deque.out");
	int N, K;
	cin >> N >> K;
	vector<int> a(N);
	deque<int> deq;
	long long sum = 0;
	for (int i = 0; i < N; i++) {
		cin >> a[i];
		while (!deq.empty() && a[deq.back()] >= a[i]) {
			deq.pop_back();
		}
		
		deq.push_back(i);
		
		if (deq.front() == i - K) {
			deq.pop_front();
		}
		
		if (i >= K - 1) {
			sum += a[deq.front()];
		}
	}
	
	cout << sum;
}
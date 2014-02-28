#include <vector>
#include <fstream>

using namespace std;

ifstream fin("lis.in");
ofstream fout("lis.out");

template <typename T>void print (T const& coll){ for (typename T::const_iterator pos=coll.begin() ,end(coll.end()); pos!=end; ++pos) fout << *pos << ' ';   fout <<'\n';}
vector<int> getLis(const vector<int> &a) {
	vector<int> b;
	vector<int> p(a.size());
	if(a.empty()) return b;
	b.push_back(0);
	int u , v , mid;
	for(size_t i = 1;i < a.size();++i) {
		
		if(a[i] > a[b.back()]) {
			p[i] = b.back();
			b.push_back(i);
			continue;
		}

		for(u = 0 , v = b.size() - 1;u < v;) {
			mid = u + ((v - u)>>1);
			if(a[b[mid]] < a[i]) u = mid + 1; else v = mid;
		}

		if(a[i] < a[b[u]]) {
			b[u] = i;
			if(u > 0) p[i] = b[u - 1]; 
		}
	}
	for(u = b.size() , v = b.back();u--;v = p[v]) b[u] = a[v];
	return b;
}

int main()
{
	int N;
	fin>>N;
	vector<int> v(N) , ans;
	for(int i = 0;i < N;i++) {
		fin>>v[i];
	}
	ans = getLis(v);
	fout<<ans.size()<<'\n';
	print(ans);
	return 0;
}

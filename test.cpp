#include <iostream>
#include <functional>

using namespace std;

namespace Unpack {

	template<int ...> struct Seq {};
	template<int N, int ...S> struct Gens : Gens<N - 1, N - 1, S...> {};
	template<int ...S> struct Gens<0, S...>{ typedef Seq<S...> type; };

	template <typename returnType, typename ...Args>
	struct saveIt {
		tuple<Args...> params;
		returnType (*func)(Args...);

		returnType delayedDispatch() {
			return callFunc(typename Gens<sizeof...(Args)>::type());
		}

		template<int ...S>
		returnType callFunc(Seq<S...>) {
			return func(get<S>(params) ...);
		}
	};


	template<typename returnType,typename...Args>
	static returnType call(tuple<Args...>& myTuple,returnType (*func)(Args...)){
		saveIt<returnType,Args...> saved = { myTuple, func };
		return saved.delayedDispatch();
	}
};


int solve(int a, int b, int c) {
	return a * a + b - 2 + c * 5;
}
int main()
{
	tuple<int, int, int> t = make_tuple( 1, 2, 3);
	cout << Unpack::call(t, solve) << "\n";
	return 0;
}


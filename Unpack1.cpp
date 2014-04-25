#include <iostream>
#include <tuple>

using namespace std;

template<size_t...> struct index_tuple{};

template<size_t I, typename IndexTuple, typename... Types>
struct make_indices_impl;

template<size_t I, size_t... Indices, typename T, typename... Types>
struct make_indices_impl<I, index_tuple<Indices...>, T, Types...> {
	typedef typename make_indices_impl<I + 1, index_tuple<Indices...,
	I>, Types...>::type type;
};

template<size_t I, size_t... Indices>
struct make_indices_impl<I, index_tuple<Indices...> > {
	typedef index_tuple<Indices...> type;
};

template<typename... Types>
struct make_indices : make_indices_impl<0, index_tuple<>, Types...>{};


template <class returnType,class... Args>
class Unpacker {
		template <size_t... Indices>
		static returnType forward0(returnType(*func)(Args...), index_tuple<Indices...>, tuple<Args...>& args) {
			return go(std::forward<Args>(get<Indices>(args))...);
		}
		tuple<Args...> members;

	public:

		Unpacker(Args... args) : members(args...) {}

		returnType forward(returnType (*func)(Args...)) {
			typedef typename make_indices<Args...>::type Indices;
			return forward0(func, Indices(), members);
		}

};

int go(int x, int y, int z) {
	return x + (int)y + (int)z;
}

int main()
{
	cout << Unpacker<int, int, int, int>(42, 5, 3).forward(go);
}


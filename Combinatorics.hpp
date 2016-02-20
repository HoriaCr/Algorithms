#ifndef COMBINATORICS_HPP_
#define COMBINATORICS_HPP_

#include <cassert>
#include <vector>
#include <algorithm>

using namespace std;

template<class DataType>
class Combinatorics {

        vector< vector<DataType> > C;
        vector<DataType> F;
    
    public:

        void precomputeBinomial(uint32_t n, uint32_t k) {
            C = computeBinomial(n, k); 
        }

        void precomputeFactorial(uint32_t n) {
           F = computeFactorial(n); 
        }

        bool nextCombination(uint32_t n, uint32_t k, vector<uint32_t>& c) {
            uint32_t t = k - 1;
            while (t != 0 && c[t] == n - k + t) {
                t--;
            }
            c[t]++;
            for (size_t i = t + 1; i < k; i++) {
                c[i] = c[i - 1] + 1;
            }

            if (c[k - 1] == n)
                return false;
            return true;
        }

        bool nextPermutation(vector<uint32_t>& p) {
            return next_permutation(begin(p), end(p));
        }


        vector< vector<DataType> > computeBinomial(uint32_t n, uint32_t k) {
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


        vector<DataType> computeFactorial(uint32_t n) {
            vector<DataType> f(n + 1);
            f[0] = 1;
            for (uint32_t i = 1; i <= n; i++) {
                f[i] = f[i - 1] * i;
            }
            return f;
        }

        DataType comb(uint32_t n, uint32_t k) {
            // return 0 if undefined
            if (n < k) return 0;
            return C[n][k];
        }

        DataType combFact(uint32_t n, uint32_t k) {
            assert (n >= k);
            return F[n] / (F[n - k] * F[k]); 
        }

        DataType partialP(uint32_t n, uint32_t k) {
            assert (n >= k);
            return F[n] / F[n - k]; 
        }
        
        DataType factorial(uint32_t n) {
            return F[n];
        }

        int getIndex(const vector<int>& cards, const int& n, const int& k) {
            int ret = 0;
            for (int i = 0; i < k; i++) {
                ret += comb(n - cards[i], k - i);
            }
            
            return comb(n, k) - ret - 1;
        }

        vector<int> indexToVector(int index, const int& n, const int& k) {
            vector<int> cards(k);
            int cardIndex = 1;
            index = comb(n, k) - index - 1; 
            for (int i = 0; i < k; i++) {
                while (
                    static_cast<int>(comb(n - cardIndex, k - i)) > index) {
                    cardIndex += 1;
                }
                cards[i] = cardIndex;
                index -= comb(n - cardIndex, k - i);
                cardIndex += 1;
            
            }
            return cards;
        }

        bool nextPartialPermutation(uint32_t n, uint32_t k,
                vector<uint32_t>& a) {
            int edge = static_cast<int>(k) - 1;
            uint32_t j = k;
            while (j < n and a[edge] >= a[j])
                ++j;
            if (j < n) {
                swap(a[edge], a[j]);
            } else {
                reverse(a.begin() + k, a.end());
                // find rightmost ascent to left of edge
                int i = edge - 1;
                while (i >= 0 and a[i] >= a[i+1])
                      --i;

                if (i < 0) 
                    return false;

                j = n - 1;
                while (static_cast<int>(j) > i and a[i] >= a[j])
                      --j;

                swap (a[i], a[j]);
                reverse(a.begin() + i + 1, a.end());
            }
            return true;
        }

        vector< vector<uint32_t> > generatePartialPermutations(
                uint32_t n, uint32_t k) {
            vector<uint32_t> a(n);
            vector< vector<uint32_t> > ans;
            for (uint32_t i = 0; i < n; i++) {
                a[i] = i;
            }
            do {
                ans.push_back(vector<uint32_t>(a.begin(), a.begin() + k));
            } while (nextPartialPermutation(n, k, a));
            return ans;
        }

        vector< vector<uint32_t> > generatePermutations(uint32_t n) {
            return generatePartialPermutations(n, n);
        }

        vector< vector<uint32_t> > generateCombinations(
                uint32_t n, uint32_t k) {

            vector< vector<uint32_t> > ans;
            vector<uint32_t> a(k);
            for (uint32_t i = 0; i < k; i++) {
                a[i] = i;
            }

            do {
                ans.push_back(a);
            } while (nextCombination(n, k, a));

            return ans;
        }
        
        vector< vector<int> > getPermutationCycles(vector<int> a) {
            vector< vector<int> > ans;
            int N = (int)a.size(), i, k;
            for(int j = 0; j < N; j++) {
                if (a[j] == -1)
                    continue;
                k = j;
                vector<int> cycle;
                while( (i = a[k] + 1) != 0) { 
                    cycle.push_back(k);
                    a[k] = -1;
                    k = --i;
                }
                ans.push_back(cycle);
            }
            return ans;
        }

        vector<int> getPermutationInverse(vector<int> a) {
            vector<int> ans(a.size());
            for (size_t i = 0; i < a.size(); i++) {
                ans[a[i]] = (int)i;
            }
            return ans;
        }
};

#endif // COMBINATORICS_HPP_

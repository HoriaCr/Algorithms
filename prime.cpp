#include <iostream>

using namespace std;

int main()
{

    int n;
    cin >> n;
    for (int i = 2;i * i <= n; i++) {
        if (n % i == 0) {
            cout << n << " is not prime \n";
            return 0;
        }
    }

    cout << n << " is prime \n";
    
    return 0;
}

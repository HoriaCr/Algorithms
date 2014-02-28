/*
 * Algorithm: z-algorithm
 * Order : linear O(n)
 */
#include<stdio.h>
#include<string.h>
#define MAX 100007
char Str[MAX+7],N; // N Str len
long z[MAX+7];
void Z_Alg( void )
{
    long i,L = 0,R = 0;
    z[0] = 0;   // set on ur own
    for( i=1;i<N;i++ ){
        if (i > R) {
            L = R = i;
            while( R < N && Str[R-L] == Str[R] ){
                R++;
            }
            z[i] = R-L; R--;
        }
        else{
            long k = i-L;
            if( z[k] < R-i+1) z[i] = z[k];
            else{
                L = i;
                while( R < N && Str[R-L] == Str[R] ){
                    R++;
                }
                z[i] = R-L; R--;
            }
        }
    }
}

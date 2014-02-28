/*
 * Algorithm: Algorithm of Aho-Corasick
 * Order : O( n )
 * Note: Root is in 0th node
 */
#include<stdio.h>
#include<string.h>
#include<vector>
#include<algorithm>
using namespace std;
#define MAX_ALPH 26
struct NODE{
    char Ch;                // char with which this node reaches
    long Par;               // Parent node
    long Child[MAX_ALPH];   // Child node for coreesponding char
    long State[MAX_ALPH];   // state for coreesponding char
    long Prefix;            // node which is also the longest suffix of Node
    NODE( char Ch = 0,long Par = 0 )
    {
        this->Ch = Ch;
        this->Par = Par;
        memset( Child,-1,sizeof(Child));
        memset( State,-1,sizeof(State));
        Prefix = -1;
    }
};
vector<NODE> Node;

void AddString( long I,char *s )
{
    long i,n = 0;
    for( i=0;s[i];i++){
        char Ch = s[i]-'A';
        if( Node[n].Child[Ch]==-1 ){
            Node.push_back( NODE( Ch,n ) );
            Node[n].Child[Ch] = Node.size()-1;
        }
        n = Node[n].Child[Ch];
    }
}
long GoState( long I,char Ch );
/* find the node having path of longest prefix which is also suffiz of I*/
long FindPrefix( long I )
{
    if( Node[I].Prefix==-1 ){
        if( I==0 || Node[I].Par==0 ){
            Node[I].Prefix = 0;
        }
        else{
            Node[I].Prefix = GoState( FindPrefix( Node[I].Par ) , Node[I].Ch );
        }
    }
    return Node[I].Prefix;
}
/* will return state from I using char CH */
long GoState( long I,char Ch )
{
    if( Node[I].State[Ch]==-1){
        if( Node[I].Child[Ch]!=-1 ){
            Node[I].State[Ch] = Node[I].Child[Ch];
        }
        else{
            Node[I].State[Ch] = I==0 ? 0:GoState( FindPrefix( I ) , Ch );
        }
    }
    return Node[I].State[Ch];
}

int main( void )
{
    Node.resize( 1 );
    return 0;
}

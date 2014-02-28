#include <iostream>
#include <fstream>
#include <string>
#include <direct.h>
#include <sys/stat.h>
 
using namespace std;

bool createProject(string problemName) {

	string directoryName = "F:\\VC surse\\" + problemName;
    struct stat st;
    if (stat(directoryName.c_str(), &st) == 0)
    {
        cout << "The project exists.\n";
    }
    else
    {
        int mkdirResult = _mkdir(directoryName.c_str());
        if (mkdirResult == 0)
        {
			directoryName += "\\" + problemName;
			ofstream cppFile(directoryName + ".cpp",ios::out | ios::app);
			string content = 
			"#include <algorithm>\n"
			"#include <string>\n"
			"#include <set>\n"
			"#include <map>\n"
			"#include <vector>\n"
			"#include <queue>\n"
			"#include <iostream>\n"
			"#include <iterator>\n"
			"#include <math.h>\n"
			"#include <cstdio>\n"
			"#include <cstdlib>\n"
			"#include <sstream>\n\n"

			"using namespace std;\n\n"

			"typedef pair<int,int> pii;\n"  
			"typedef long long i64;\n"
			"typedef vector<int> vi;\n\n"

			"#define UN(v) SORT(v),v.erase(unique(v.begin(),v.end()),v.end())\n"
			"#define SORT(c) sort((c).begin(),(c).end())\n"
			"#define FOR(i,a,b) for (int i = (a); i < (b); i++)\n"
			"#define REP(i,n) FOR(i,0,n)\n"
			"#define CL(a,b) memset(a,b,sizeof(a))\n"
			"#define pb push_back\n"
			"\n"
			;
			content += "ifstream fin(\"" +  problemName + ".in\");\n";
			content += "ofstream fout(\"" + problemName + ".out\");\n\n";

			content += "int main()\n" 
					   "{\n\n\n\n"
					   "	return 0;\n}\n"
					;	
			cppFile<<content;
			ofstream inFile(directoryName + ".in",ios::out | ios::app);
			ofstream outFile(directoryName + ".out",ios::out | ios::app);
			cppFile.close();
			inFile.close();
			outFile.close();
            cout << "The project is created.\n";
			return true;
        }
        else
        {
            cout << "The project creation failed with error: " + mkdirResult << endl;
			return false;
        }
    }
	return false;
}
 
int main ()
{
	string problemName;
	cout<<"Type problem name and press Enter!\n";
	cin>>problemName;
	string input;
	cout<<"If you want to continue type yes (lower letters).\n";
	cin>>input;
	if(input == "yes") {
		createProject(problemName);
	} else {
		cout<<"Thanks!\n";
	}
	return 0;
}

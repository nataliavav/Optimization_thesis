//Schaffer function N. 1
//(-100,100)

#include <cmath>
#include <iostream>
#include <fstream>
#include <valarray>
#include <tgmath.h>

using namespace std;

int main()
{
	//
	// Read design variables
	std::ifstream in ("task.dat",std::ios::in);
	int n; in >> n;
	valarray<double> x; x.resize(n); x=0;
	for (int i=0;i<n;i++) in >> x[i];
	in.close();
	//

	
	double f1=x[0]*x[0];

	double f2=(x[0]-2)*(x[0]-2);
	
/*	std::ofstream out1 ("task.cns",std::ios::out);
	out1<<std::scientific<<f1<<endl;
	out1.close();*/

	// Compute & Return
	std::ofstream out("task.res",std::ios::out);
	out<<std::scientific<<f1<<std::endl;
	out<<std::scientific<<f2<<std::endl;
	out.close();
	//
	
	
	return 0;
}



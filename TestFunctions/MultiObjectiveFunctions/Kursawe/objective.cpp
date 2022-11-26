//Kursawe function

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

	
	double f1=0;
	for (int i=0;i<n-1;i++){
		f1+=-10.0*exp(-.2*sqrt(x[i]*x[i]+x[i+1]*x[i+1]));
	}

	double f2=0;
	for (int i=0;i<n;i++){
		f2+=pow(abs(x[i]),.8)+5.0*sin(x[i]*x[i]*x[i]);
	}
	
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



//Rastrigin. min=0//   
#include <cmath>
#include <iostream>
#include <fstream>
#include <valarray>

int main()
{
	//
	// Read design variables
	std::ifstream in ("task.dat",std::ios::in);
	int n; in >> n;
	std::valarray<double> x; x.resize(n); x=0;
	for (int i=0;i<n;i++) in >> x[i];
	in.close();
	//
	// Compute f
	const double pi   = 4.0*atan(1.0);	
	const double alfa = 10.0;	
	double f=alfa*n;
	for (int i=0; i<n; i++){
  		const double xi = x[i];
  		f+=xi*xi-alfa*cos(2.0*pi*xi);
 	}	
	//
	// Return
	std::ofstream out("task.res",std::ios::out);
	out.precision(14);
	out << std::scientific << f << std::endl;
	out.close();
	//
	return 0;
}

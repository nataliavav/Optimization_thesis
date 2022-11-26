//Sx^2
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
	double f=0;
	for (int i=0; i<n; i++){
  		const double xi = x[i];
  		f+=xi*xi;
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

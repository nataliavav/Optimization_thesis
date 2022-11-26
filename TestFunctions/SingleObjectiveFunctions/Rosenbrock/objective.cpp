//Rosenbrock's logarithmic banana function. Min=0, x=1//
//(-10,10) (DEN KSERW POY THN BRHKA)
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
	double f=log(1+(1-x[0])*(1-x[0])+100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0]));

	//
	// Return
	std::ofstream out("task.res",std::ios::out);
	out.precision(14);
	out << std::scientific << f << std::endl;
	out.close();
	//
	return 0;
}

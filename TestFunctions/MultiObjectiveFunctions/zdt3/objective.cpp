#include <cmath>
#include <iostream>
#include <fstream>
#include <valarray>


using namespace std;

const double pi=4.0*atan(1.0);

const double ffun(const valarray<double>& x)
{
	return x[0];
}

const double gfun(const valarray<double>& x)
{
	const int nn=x.size();
	double sum=0;
	if (nn==1) return 1;
	for (int i=1;i<nn;i++)sum+=x[i];
	return 1.0+sum*9.0/double(nn-1);
}

const double hfun(const double f1, const double g)
{
	return 1-sqrt(f1/g)-f1/g*sin(10.0*pi*f1);
}


const double f1(const valarray<double>& x)	
{
	return ffun(x);
}

const double f2(const valarray<double>& x)
{
	return gfun(x)*hfun(ffun(x),gfun(x));
}


int main()
{
	//
	// Read design variables
	std::ifstream in ("F:\\task.dat",std::ios::in);
	int n; in >> n;
	valarray<double> x; x.resize(n); x=0;
	for (int i=0;i<n;i++) in >> x[i];
	in.close();
	//
	// Compute & Return
	std::ofstream out("F:\\task.res",std::ios::out);
	out<<std::scientific<<f1(x)<<std::endl;
	out<<std::scientific<<f2(x)<<std::endl;
	out.close();
	//
	return 0;
}

//Poloni's two objective function:
//(-ð,ð)
#include <cmath>
#include <iostream>
#include <fstream>
#include <valarray>


using namespace std;



const double B1(const valarray<double>& x)
{
	double B=0.5*sin(x[0])-2.0*cos(x[0])+sin(x[1])-1.5*cos(x[1]);
	return B;
}

const double B2(const valarray<double>& x)
{
	double B=1.5*sin(x[0])-cos(x[0])+2.0*sin(x[1])-0.5*cos(x[1]);
	return B;
}

/*const double Á_1()
{
	double A=0.5*sin(1.0)-2.0*cos(1.0)+sin(2.0)-1.5*cos(2.0);
	return A;
}

const double Á_2()
{
	double A=1.5*sin(1.0)-cos(1.0)+2.0*sin(2.0)-0.5*cos(2.0);
	return A;
}*/

const double f1(const valarray<double>& x)	
{
	double F=1+(0.5*sin(1.0)-2.0*cos(1.0)+sin(2.0)-1.5*cos(2.0)-B1(x))*(0.5*sin(1.0)-2.0*cos(1.0)+sin(2.0)-1.5*cos(2.0)-B1(x))+(1.5*sin(1.0)-cos(1.0)+2.0*sin(2.0)-0.5*cos(2.0)-B2(x))*(1.5*sin(1.0)-cos(1.0)+2.0*sin(2.0)-0.5*cos(2.0)-B2(x));
	return F;
}

const double f2(const valarray<double>& x)
{
	double F=(x[0]+3)*(x[0]+3)+(x[1]+1)*(x[1]+1);
	return F;
}


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
	// Compute & Return
	std::ofstream out("task.res",std::ios::out);
	out<<std::scientific<<f1(x)<<std::endl;
	out<<std::scientific<<f2(x)<<std::endl;
	out.close();
	//
	return 0;
}

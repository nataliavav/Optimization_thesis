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
    double f;
    const double Q = x[0];
    const double d = x[1];
    double c;
    if (Q<500){
       c=0.3*Q;
    }
    else if (Q<1000){
        c=0.29+5/Q;
    }
    else{
        c=0.28+15/Q;
    }
    f=600.0*c+8.0*600.0/Q+0.2*c/2.0*(Q-d)*(Q-d)/Q+1/2.0*5.0*d*d/Q;
    //
    // Return
    std::ofstream out("task.res",std::ios::out);
    out.precision(14);
    out << std::scientific << f << std::endl;
    out.close();
    //
    return 0;
}
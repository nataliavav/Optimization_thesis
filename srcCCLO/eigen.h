

#ifndef _EIGEN_H
#define _EIGEN_H


#include "matrix.h"

bool RealSymEigen(const Matrix<double>& A, Vector<double>& l,
		  Matrix<double>& v, const int maxiter);
void PCACompress(const Vector<double>& v, const Matrix<double>& eiv,
		 Vector<double>& x);
void PCADecompress(const Vector<double>& v, const Matrix<double>& eiv,
		 Vector<double>& x);

namespace NEIGEN
{
	double pythag(const double a, const double b);
	void tred2(double **a, int n, double *d, double *e);
	bool tqli(double *d, double *e, int n, double **z, const int maxiter);

};

#endif



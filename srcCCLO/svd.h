#include "matrix.h"

#ifndef _SVD_H
#define _SVD_H


bool SVDSolve(const Matrix<double>& A, Matrix<double>& X,
	      const Matrix<double>& B, const int maxiter,
	      const double wtol);

namespace NSVD
{
	double pythag(const double a, const double b);
	bool svdcmp(double **a, int m, int n, double *w, 
		    double **v, int maxiter);
	void dsvbksb(double **u, double *w, double **v, int m, int n, 
	     	     double *b, double *x);
};


#endif

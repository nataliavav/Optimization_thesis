/////////////////////////////////////////////////////////////////////////////
//
//
//                 My basic linear algebra subroutines
//
//     10-01-2005    I.C. Kampolis
//
//
/////////////////////////////////////////////////////////////////////////////
//


#ifndef MYBLAS_H
#define MYBLAS_H

#include "matrix.h"
#include "mymath.h"
#include "svd.h"

template<class T>
void FullSearch(const Matrix<T>& A, long int& ir, long int& ic,
		const long int n, const long int row);


template <class T> inline void SwapCol(Matrix<T>& A, long int* X,
		long int c1, long int c2, long int N);

template <class T, class R> inline void SwapRow(Matrix<T>& A, Matrix<R>& B,
		long int r1, long int r2, long int N, long int M);

template <class T, class R> bool Gauss_Elim(Matrix<T>& A, Matrix<R>& X, 
		Matrix<R>& B);

template <class T> bool QPSolve(Matrix<T>& Q, Matrix<T>& c, 
				Matrix<T>& A, Matrix<T>& b, Matrix<T>& x);

template <class T> bool QPSolveI(const Matrix<T>& G, const Matrix<T>& d, 
				 const Matrix<T>& A, const Matrix<T>& b, 
				 Matrix<T>& x, const int maxiter=100,
				 const bool bdisp=false);

bool QPSolveI_D(const Matrix<double>& G, const Matrix<double>& d, 
		const Matrix<double>& A, const Matrix<double>& b, 
		Matrix<double>& x,const int maxiter,
		const bool bdisp);

#endif

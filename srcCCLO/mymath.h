// 
// Contains basic math routines.
// 
// dmax1, dmin1, imax1, imin1 (as fortran)
// Declaration of dmin,dmax,deps,pi
// min1, max1 Templates
// Determinant, Cofactor, Transpose, Adjoint, Inverse for matrices
// IsFinite() -- Checks for NaN  (Warning, do not forget to compile with -DWIN32 for windows)
// tri_diag -- Solves a tridiagonal system
// Gauss_Elim -- Gauss Elimination with Full pivoting
// 
// Ioannis C. Kampolis LTT/NTUA 2005




#ifndef MYMATH_H
#define MYMATH_H


#include <float.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>
	#ifdef WIN32
		#define IsFinite _finite
	#else //LINUX
		#define IsFinite __finite
	#endif


//Fortran like routines 
inline const double& dmax1(const double& a, const double& b) 
		{if (a>b) return a; return b;}
inline const double& dmin1(const double& a, const double& b) 
		{if (a<b) return a; return b;}
inline const double  dsgn1(const double& a) 		
		{if (a>=0) return +1; return -1;}
//
inline const long& imax1(const long& a, const long& b) 
		{if (a>b) return a; return b;}
inline const long& imin1(const long& a, const long& b) 
		{if (a<b) return a; return b;}
//
//General Templates min,max,absolute value
template<class T> inline const T& max1(const T& a, const T& b) 
		{if (a>b) return a; return b;}
template<class T> inline const T& min1(const T& a, const T& b) 
		{if (a<b) return a; return b;}
//
template<class T> inline const T abs1(const T& a) 
		{if (a<0) return -a; return a;}
template<class T> inline const T dabs(const T& a) 
		{if (a<0) return -a; return a;}
inline const long factorial(const long a) 
		{if (a==0) return 1;  return a*factorial(a-1);}
inline const double dsgn2(const double& a, const double& b)
		{if (b>=0) return dabs(a); return -dabs(a);}
//
//Some constants
//const double deps=std::numeric_limits<double>::epsilon(); 
const double deps=1.e-13;
const double dmax=std::numeric_limits<double>::max()/1000;
const double dmin=std::numeric_limits<double>::min();
const double pi=4.0*atan(1.0);

inline double Hc(const double& x, const double& c)
{
	if ((x-c)<=0)	return 0;
	return 1;
}

inline double Dc(const double& x, const double& c, const double ee=.1)
{
	const double d=x-c;
//	return ee*ee/(d*d+ee*ee);
	return ee/(d*d+ee*ee);
//	if (dabs(x-c)>deps) return 0;
//	return id;
}

inline double CumNormDist(double x)
{
	int neg = (x<0);
	if(neg) x *= -1;
	double k(1/(1+0.2316419*x));
	double y=((((1.330274429*k-1.821255978)*k+1.781477937)*k
				-0.356563782)*k+0.319381530)*k;
	y = 1.0 - 0.398942280401*exp(-0.5*x*x)*y;
	return (1-neg)*y + neg*(1-y);
}

namespace MATH
{
	//Solution of a tridiagonal system of linear equations
	//using the Thomas algorithm.
	//
	//It will ruin the diag vector
	// 
	//n   : number of equations
	//m   : number of columns  
	//bef : element before the diag.    bef[n]
	//diag: diagonal elements          diag[n]
	//aft : element after the diag.     aft[n]
	//rh  : right hand side              rh[n][m]   //In the end it contains the solution
	inline bool tri_diag(const long n,const long m,double* bef,double* diag,double* aft,double** rh)
	{
		std::cout <<"TRIDIAG NOT OK!\n";
		exit(0);
		long i,j; //the only indexes used
		//Elimination
		for (i=1;i<n;i++)
		{
			if (dabs(diag[i-1])<deps) {return false;} //divide by zero protection
			double r=bef[i]/diag[i-1];
			diag[i]-=r*aft[i-1];
			for (j=0;j<m;j++) rh[i][j]-=r*rh[i-1][j];
		}
		//Compute rhs[n-1][]
		for (j=0;j<m;j++) rh[n-1][j]/=diag[n-1];
		//Back substitution
		for (i=n-2;i>=0;i--)
		{
			for (j=0;j<m;j++) rh[i][j]=(rh[i][j]-aft[i]*rh[i+1][j])/diag[i];
		}
		return true;
	}


	
	inline double Determinant(double** a,long n)
	{
		long i,j,j1,j2;
		double det=0;
		double **m=0;

		if (n<1) {}
		else if (n==1) {det=a[0][0];}
		else if (n==2) {det=a[0][0]*a[1][1]-a[1][0]*a[0][1];}
		else
		{
			det=0;
			for (j1=0;j1<n;j1++)
			{
	  			m= new double*[n-1];
				for (i=0;i<n-1;i++) m[i]=new double[n-1];
				for (i=1;i<n;i++)
				{
					j2=0;
					for (j=0;j<n;j++)
					{
						if (j==j1) continue;
						m[i-1][j2]=a[i][j];
						j2++;
					}
				}
				det+=pow(-1.0,j1+2.0)*a[0][j1]*Determinant(m,n-1);
				for (i=0;i<n-1;i++) delete[] m[i];
				delete[] m;
			}
		}
		return det;
	}


	inline void Transpose(double **a ,long n)
	{
		long i,j;
		double tmp;
		for (i=1;i<n;i++)
		{
			for (j=0;j<i;j++)
			{
				tmp=a[i][j];
				a[i][j]=a[j][i];
				a[j][i]=tmp;
			}
		}
	}
	

	inline void CoFactor (double **a, long n, double** b)
	{
		long i,j,ii,jj,i1,j1;
		double det;
		double** c;
		c= new double*[n-1];
		for (i=0;i<n-1;i++) c[i]=new double[n-1];
		for (j=0;j<n;j++)
		{
			for (i=0;i<n;i++)
			{
				i1=0;
				for (ii=0;ii<n;ii++)
				{
					if (ii==i) continue;
					j1=0;
					for (jj=0;jj<n;jj++)
					{
						if (jj==j) continue;
						c[i1][j1]=a[ii][jj];
						j1++;
					}
					i1++;
				}
				det=Determinant(c,n-1);
				b[i][j]=pow(-1.0,i+j+2.0)*det;
			}
		}
		for (i=0;i<n-1;i++) delete[] c[i];
		delete[] c;
	}

	inline void Adjoint (double **a, long n, double** b)
	{
		CoFactor(a,n,b);
		Transpose(b,n);
	
	}

	inline bool Inverse (double **a, long n, double **ainv)
	{
		double Det=Determinant(a,n);
		if (dabs(Det)<deps) return false;
		Adjoint(a,n,ainv);
		for (long i=0;i<n;i++) for (long j=0;j<n;j++) ainv[i][j]/=Det;
		return true;
	}



	//GAUSS ELIMINATION ROUTINE +Helpers
	inline void FullSearch(double ** A,long& ir, long& ic,const long n, const long row)
	{
		double max=-1;
		ir=0;ic=0;
		for (long i=row;i<n;i++)
		{
			for (long j=row;j<n;j++)
			{
				if (dabs(A[i][j])>max) 
				{
					max=dabs(A[i][j]);
					ir=i; ic=j;
				}
			}
		}
	}

	template<class T> inline void SwapCol(double** A, T* X, long c1, long c2,long N)
	{
		double c;
		for (long i=0;i<N;i++)
		{
			c=A[i][c1]; A[i][c1]=A[i][c2]; A[i][c2]=c;
		}

		T x2;
		x2=X[c1]; X[c1]=X[c2]; X[c2]=x2;
	}

	template<class T> inline void SwapRow(double** A,T* B, long r1, long r2,long N)
	{
		double c;
		for (long i=0;i<N;i++)
		{
			c=A[r1][i]; A[r1][i]=A[r2][i]; A[r2][i]=c;
		}

		T c2;
		c2=B[r1]; B[r1]=B[r2]; B[r2]=c2;
	}

	template<class T> inline bool Gauss_Elim(double** A, T* X, T* B, long N)
	{
		long *x2=new long[N];
		long rowout,columnout;
		for (long i=0;i<N;i++) x2[i]=i;

		//Elimination
		for (long k=0;k<N-1;k++)
		{
			if (dabs(A[k][k])<deps)
			{
				FullSearch(A,rowout,columnout,N,k);
				SwapRow(A,B,rowout,k,N);
				SwapCol(A,x2,columnout,k,N);
			}
			if (dabs(A[k][k])<deps)
			{
				delete[] x2; x2=0; return false;
			}
			for (long i=k+1;i<N;i++)
			{
				double factor=A[i][k]/A[k][k];
				for (long j=k+1;j<N;j++)
				{
					A[i][j]-=factor*A[k][j];
				}
				B[i]-=factor*B[k];
			}
		}
		//Back Substitution
		X[x2[N-1]]=B[N-1]/A[N-1][N-1];
		for (long i=N-2;i>=0;i--)
		{
			T sum=0;
			for (long j=i+1;j<N;j++) sum+=A[i][j]*X[x2[j]];
			X[x2[i]]=(B[i]-sum)/A[i][i];
		}
		delete[] x2; x2=0;
		return true;
	}	

		
};


	

		

				








#endif

#include "myblas.h"
#include <cmath>
#include "svd.h"

template <class T> bool QPSolveI(const Matrix<T>& G, const Matrix<T>& d, 
				const Matrix<T>& A, const Matrix<T>& b, 
				Matrix<T>& x,const int maxiter,
				const bool bdisp)
{
	std::cerr<<"TEMPORARILY DISABLED!!!!\n";
	exit(0);
	//.5xGx +dx : Ax>=b
	////For how it is solved look: 
	//	J. Nocedal & S.J. Wright, Numerical Optimization, Ch.16.7,
	//	
	//				
	//Solved using an interior point method
	const int varN=G.dim1();
	const int conN=A.dim1();
	const double sigma=.5;	//Algorithm specific
	//
	bool converged=false;
	int iter=0;
	Matrix<T> l; l.Resize(conN,1,1.0);
	Matrix<T> y; y.Resize(conN,1,0.0);
	Matrix<T> rd; rd.Resize(varN,1,0.0);
	Matrix<T> rb; rb.Resize(conN,1,0.0);
	Matrix<T> yil; yil.Resize(conN,conN,0.0);
	T mi;
	Matrix<T> lhs;
	Matrix<T> rhs;
	Matrix<T> delx; delx.Resize(varN,1,0);
	Matrix<T> dell; dell.Resize(conN,1,0);
	Matrix<T> dely; dely.Resize(conN,1,0);
	Matrix<T> ep; ep.Resize(conN,1,0.0);
	Vector<T> xold; xold.resize(varN);
	Matrix<T> At; At.Resize(A.dim2(),A.dim1(),0);
	for (int k=0;k<A.dim1();k++)
		for (int j=0;j<A.dim2();j++)
			At[j][k]=A[k][j];
	while (iter<maxiter && !converged)
	{
		y=(A*x) -b;	
		rd=(G*x)-(At*l)+d;
		rb=(A*x)-y-b;
		yil=0;
		T mi=0;
		for (int k=0;k<conN;k++)
		{
			yil[k][k]=l[k][0]/(y[k][0]+deps);
			ep[k][0]=sigma/(l[k][0]+deps);
			mi+=y[k][0]*l[k][0];
		}
		for (int k=0;k<varN;k++) xold[k]=x[k][0];
		mi/=T(conN);
		for (int k=0;k<conN;k++) ep[k][0]=ep[k][0]*mi+deps;
		lhs=G+(At*yil*A);
		rhs=((At*yil)*(ep-y-rb))-rd;
		delx=0;
		const double wtol=deps*1000.0;
		const int svditer=100;
		if (!Gauss_Elim(lhs,delx,rhs)) 
		{
			if (bdisp) std::cout <<"GE: Error 1\n";
			return false;
		}
		//Now compute the lagrange multipliers
		dell=(A*delx) +rb; 
		dely=(A*delx);// -b;
		//I have found delx and dell
		T dx=0;
		//Compute the step and correct
		double alphamax=1;
		for (int k=0;k<conN;k++)
		{
			double al=-l[k][0]/dell[k][0];
			double ay=-y[k][0]/dely[k][0];
			if (al<0) al=1;
			if (ay<0) ay=1;
			if (dabs(l[k][0])<deps && dell[k][0]>0) al=1;
			if (dabs(y[k][0])<deps && dely[k][0]>0) ay=1;
			if (k==0) alphamax=al;
			alphamax=min1(alphamax,al);
			alphamax=min1(alphamax,ay);
		}
		alphamax=min1(alphamax,1.0);
		
			
		for (int k=0;k<varN;k++)
		{
			x[k][0]+=alphamax*delx[k][0];
			dx+= (x[k][0]-xold[k])*(x[k][0]-xold[k]);
		}
		for (int k=0;k<conN;k++)
			l[k][0]+=alphamax*dell[k][0];

			
		dx/=T(varN);
		if (bdisp) std::cout	<<iter
					<<"   \t   "<<dx
					<<"   \t   "<<log10(dx)
					<<"   \t   "<<alphamax<<std::endl;
		if (dx<deps) converged=true;
		iter++;
	}
	return converged;
}
	
	



template <class T> bool QPSolve(Matrix<T>& Q, Matrix<T>& c, 
				Matrix<T>& A, Matrix<T>& b, Matrix<T>& x)
{
	//xQx +cx : Ax=b
	//
	Matrix<T> I,Qinv,t1,t1_inv;
	Matrix<T> h=A*x-b; 
	//Compute the inverse of Q
	const long int n=Q.dim1();
	I.Resize(n,n,0);
	for (long int i=0;i<n;i++) I[i][i]=1.0;
	//
	if (!Gauss_Elim(Q,Qinv,I)) return false;
	t1=A*Qinv*A.Transpose();
	const long int nt=t1.dim1();
	/////
	I.Resize(nt,nt,0);
	for (long int i=0;i<nt;i++) I[i][i]=1.0;
	if (!Gauss_Elim(t1,t1_inv,I)) return false;
	x=t1_inv*(h-A*Qinv*c); 
	////
	x=Qinv*(A.Transpose()*x+c);  x*=-1;
	return true;
}
	
	
	


template<class T>
void FullSearch(const Matrix<T>& A, long int& ir, long int& ic,
		const long int n, const long int row)
{
	double max=-1;
	ir=0; ic=0;
	for (long int i=row;i<n;i++)
		for (long int j=row;j<n;j++)
		{
			if (dabs(A[i][j])>max)
			{
				max=dabs(A[i][j]);
				ir=i; ic=j;
			}
		}
}

template <class T> inline void SwapCol(Matrix<T>& A, long int* X,
		long int c1, long int c2, long int N)
{
	T c;
	for (long int i=0;i<N;i++)
	{
		c=A[i][c1]; A[i][c1]=A[i][c2]; A[i][c2]=c;
	}
	long int x2;
	x2=X[c1]; X[c1]=X[c2]; X[c2]=x2;
}

template <class T, class R> inline void SwapRow(Matrix<T>& A, Matrix<R>& B,
		long int r1, long int r2, long int N, long int M)
{
	T c;
	for (long int i=0;i<N;i++)
	{
		c=A[r1][i]; A[r1][i]=A[r2][i]; A[r2][i]=c;
	}
	R c2;
	for (long int j=0;j<M;j++)
	{
		c2=B[r1][j]; B[r1][j]=B[r2][j]; B[r2][j]=c2;
	}
}


template <class T, class R> bool Gauss_Elim(Matrix<T>& A, Matrix<R>& X, 
		Matrix<R>& B)
{
	const long int N=A.dim1();
	const long int M=B.dim2();
	X.Resize(N,M);
	if (N!=A.dim2()||N!=B.dim1()) return false;
	long int rowout,columnout;
	long int *x2=new long int[N];
	for (long int i=0;i<N;i++) x2[i]=i;
	//Elimination
	for (long int k=0;k<N-1;k++)
	{
		//if (dabs(A[k][k])<deps)
		{
			FullSearch(A,rowout,columnout,N,k);
			SwapRow(A,B,rowout,k,N,M);
			SwapCol(A,x2,columnout,k,N);
		}
		if (dabs(A[k][k])<deps)
		{
			delete[] x2; x2=0; return false;
		}
		for (long int i=k+1;i<N;i++)
		{
			T factor=A[i][k]/A[k][k];
			for (long int j=k+1;j<N;j++)
			{
				A[i][j]-=factor*A[k][j];
			}
			for (long int j=0;j<M;j++)
				B[i][j]-=factor*B[k][j];
		}
	}
	//Back Substitution
	for (long int m=0;m<M;m++)
	{
		X[x2[N-1]][m]=B[N-1][m]/A[N-1][N-1];
		for (long int i=N-2;i>=0;i--)
		{
			T sum=0;
			for (long int j=i+1;j<N;j++)
				sum+=A[i][j]*X[x2[j]][m];
			X[x2[i]][m]=(B[i][m]-sum)/A[i][i];
		}
	}
	delete[] x2; x2=0;
	return true;
}	
	

#ifndef INCLUDE_TEMPLATES_ONLY
bool QPSolveI_D(const Matrix<double>& G, const Matrix<double>& d, 
		const Matrix<double>& A, const Matrix<double>& b, 
		Matrix<double>& x,const int maxiter,
		const bool bdisp)
{
	//.5xGx +dx : Ax>=b
	////For how it is solved look: 
	//	J. Nocedal & S.J. Wright, Numerical Optimization, Ch.16.7,
	//	
	//				
	//Solved using an interior point method
	const int varN=G.dim1();
	const int conN=A.dim1();
	const double sigma=.5;	//Algorithm specific
	//
	bool converged=false;
	int iter=0;
	Matrix<double> l; l.Resize(conN,1,1.0);
	Matrix<double> y; y.Resize(conN,1,0.0);
	Matrix<double> rd; rd.Resize(varN,1,0.0);
	Matrix<double> rb; rb.Resize(conN,1,0.0);
	Matrix<double> yil; yil.Resize(conN,conN,0.0);
	double mi;
	Matrix<double> lhs;
	Matrix<double> rhs;
	Matrix<double> delx; delx.Resize(varN,1,0);
	Matrix<double> dell; dell.Resize(conN,1,0);
	Matrix<double> dely; dely.Resize(conN,1,0);
	Matrix<double> ep; ep.Resize(conN,1,0.0);
	Vector<double> xold; xold.resize(varN);
	Matrix<double> At; At.Resize(A.dim2(),A.dim1(),0);
	Matrix<double> sme; sme.Resize(conN,1,0.0);
	Matrix<double> Y2; Y2.Resize(conN,conN,0.0);
	for (int k=0;k<A.dim1();k++)
		for (int j=0;j<A.dim2();j++)
			At[j][k]=A[k][j];
	y=(A*x) -b;	
	while (iter<maxiter && !converged)
	{
		rb=(A*x)-y-b;
		rd=(G*x)-(At*l)+d;
		//
		yil=0;
		sme=0;
		Y2=0;
		double mi=0;
		for (int k=0;k<conN;k++)
		{
			yil[k][k]=l[k][0]/(y[k][0]+deps);
			ep[k][0]=sigma/(l[k][0]+deps);
			mi+=y[k][0]*l[k][0];
			Y2[k][k]=1.0/(y[k][0]+deps);
		}
		for (int k=0;k<varN;k++) xold[k]=x[k][0];
		mi/=double(conN);
		//
		for (int k=0;k<conN;k++) ep[k][0]=ep[k][0]*mi+deps;
		lhs=G+(At*yil*A);
		rhs=((At*yil)*(ep-y-rb))-rd;
		delx=0;
		const double wtol=deps*1000.0;
		const int svditer=100;
		if (!SVDSolve(lhs,delx,rhs,svditer,wtol)) 
		{
			if (bdisp) std::cout <<"GE: Error 1\n";
			return false;
		}
		//Now compute the lagrange multipliers
		dely=(A*delx)+rb;
		for (int k=0;k<conN;k++)
			sme[k][0]=-l[k][0]*y[k][0]+sigma*mi-l[k][0]*dely[k][0];
		dell=Y2*sme;

//		dell=(A*delx) +rb; 
//		dely=(A*delx);//-b;
		//I have found delx and dell
		double dx=0;
		//Compute the step and correct
		double alphamax=1;
		for (int k=0;k<conN;k++)
		{
			double al=-l[k][0]/(dell[k][0]+deps);
			double ay=-y[k][0]/(dely[k][0]+deps);
			if (al<0) al=1;
			if (ay<0) ay=1;
			if (dabs(l[k][0])<deps && dell[k][0]>0) al=1;
			if (dabs(y[k][0])<deps && dely[k][0]>0) ay=1;
			if (k==0) alphamax=al;
			alphamax=min1(alphamax,al);
			alphamax=min1(alphamax,ay);
		}
		alphamax=min1(alphamax,1.0);
		
			
		for (int k=0;k<varN;k++)
		{
			x[k][0]+=alphamax*delx[k][0];
			dx+= (x[k][0]-xold[k])*(x[k][0]-xold[k]);
		}
		for (int k=0;k<conN;k++)
		{
			l[k][0]+=alphamax*dell[k][0];
			y[k][0]+=alphamax*dely[k][0];
		}

			
		dx/=double(varN);
		if (bdisp) std::cout	<<iter
					<<"   \t   "<<dx
					<<"   \t   "<<log10(dx)
					<<"   \t   "<<alphamax<<std::endl;
		if (dx<deps) converged=true;
		iter++;
	}
	return converged;
}
#endif


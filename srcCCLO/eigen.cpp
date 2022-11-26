#include "mymath.h"
#include "eigen.h"
#include "util.h"


bool RealSymEigen(const Matrix<double>& A, Vector<double>& l,
		  Matrix<double>& v, const int maxiter)
{
	const int n=A.dim1();
	if (A.dim2()!=n) return false;
	l.resize(n);
	v.Resize(n,n,0.0);

	double** aa = new double*[n+1];
	double*   d = new double[n+1];
	double*   e = new double[n+1];
	for (int i=0;i<=n;i++)
	{
		aa[i]=new double[n+1];
		if (i>0)
			for (int j=1;j<=n;j++) aa[i][j]=A[i-1][j-1];
		d[i]=e[i]=0.0;
	}
	//Rotate to tridiag
	NEIGEN::tred2(aa,n,d,e);
	const bool b1=NEIGEN::tqli(d,e,n,aa,maxiter);
	if (b1)
	{
		for (int i=0;i<n;i++)
		{
			l[i]=d[i+1];
			//Eigenvectors in v are stored in ROWS
			//while in aa are stored in COLUMNS
			for (int j=0;j<n;j++)
				v[i][j]=aa[j+1][i+1];
		}
	}
	for (int i=0;i<=n;i++) delete[] aa[i];
	delete[] aa; aa=0;
	delete[] d;  d=0;
	delete[] e;  e=0;
	return b1;
}
	
	
void PCACompress(const Vector<double>& v, const Matrix<double>& eiv,
		 Vector<double>& x)
{
	const int n=v.size();
	if (n!=eiv.dim1())
	{
		x.assign(v);
		return;
	}
	//
	x.resize(n);
	//
	for (int i=0;i<n;i++)
	{
		x[i]=0;
		for (int j=0;j<n;j++)
			x[i]+=eiv[j][i]*v[j];
	}

}

void PCADecompress(const Vector<double>& v, const Matrix<double>& eiv,
		 Vector<double>& x)
{
	const int n=v.size();
	if (n!=eiv.dim1())
	{
		x.assign(v);
		return;
	}
	//
	x.resize(n);
	//
	for (int i=0;i<n;i++)
	{
		x[i]=0;
		for (int j=0;j<n;j++)
			x[i]+=eiv[i][j]*v[j];
	}
}


double NEIGEN::pythag(const double a, const double b)
{
	//return sqrt(a*a+b*b);
	//Basically it returns the above expression protecting it
	//from underflows and overflows
	double p,r,s,t,u;
	p=dmax1(dabs(a),dabs(b));
	if (p==0.0) return p;
	r=dmin1(dabs(a),dabs(b))/p;	r*=r;
	while(true)
	{
		t=4.0+r;
		if (t==4.0) return p;
		s=r/t;
		u=1.0+2.0*s;
		p*=u;
		r*=(s/u)*(s/u);
	}
}

void NEIGEN::tred2(double **a, int n, double *d, double *e)
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += dabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=1;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[1]=0.0;
	e[1]=0.0;
	// Contents of this loop can be omitted if eigenvectors not
	// wanted except for statement d[i]=a[i][i]; 
	for (i=1;i<=n;i++) {
		l=i-1;
		if (d[i]) {
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}


bool NEIGEN::tqli(double *d, double *e, int n, double **z, const int maxiter)
{
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i];
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) {
				dd=dabs(d[m])+dabs(d[m+1]);
				if ((dabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == maxiter) return false;
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+dsgn2(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=1;k<=n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
	return true;
}

#include "matrix.h"

#ifndef INCLUDE_TEMPLATES_ONLY
estring Matrix_exception::label="Matrix exception:";
#endif


#define MATBOUNDCHECK

/////////////////////////////---- MATRIX ROW
template <class T>
T& Matrix_row<T>::operator[](const size_t& i)
{
#ifdef MATBOUNDCHECK
	if (i>n2||i<0) throw Matrix_exception
		("Matrix row out of bounds\n");
#endif
	return m[i];
}

template <class T>
const T& Matrix_row<T>::operator[](const size_t& i) const
{
#ifdef MATBOUNDCHECK
	if (i>n2||i<0) throw Matrix_exception
		("Matrix row out of bounds\n");
#endif
	return m[i];
}

template <class T>
Matrix_row<T>::Matrix_row(const Matrix_row<T>& c)
{
	n2=c.n2;
	m=new T[n2];
	for (size_t i=0;i<n2;i++) m[i]=c.m[i];
}

template <class T>
Matrix_row<T>::Matrix_row(const size_t n)
{
	n2=n;
	m=new T[n2];
}
/////////////////////////////---- MATRIX REPRESENTATION
///// Matrix access!!
template <class T> 
T& Matrix_repg<T>::operator()(const size_t& i, const size_t& j) const
{
#ifdef MATBOUNDCHECK
	if (i>=n1||j>=n2) throw Matrix_exception
		("Matrix out of bounds!\n");
#endif
	return m[i][j];
}

template <class T> 
Matrix_row<T>& Matrix_repg<T>::operator[](const size_t& j) 
{
#ifdef MATBOUNDCHECK
	if (j>=n1) throw Matrix_exception
		("Matrix out of bounds!\n");
#endif
	return m[j];
}


template <class T> 
const Matrix_row<T>& Matrix_repg<T>::operator[](const size_t& j) const
{
#ifdef MATBOUNDCHECK
	if (j>=n1) throw Matrix_exception
		("Matrix out of bounds!\n");
#endif
	return m[j];
}


///// General utilites for matrix representation
template <class T> void Matrix_repg<T>::Destroy()
{
	if (m)
	{
//		for (size_t i=0;i<n1;i++)
//			delete[] m[i];
		delete[] m;
		m=0; n1=0; n2=0;
	}
}

template <class T> void Matrix_repg<T>::Resize(size_t _n1, size_t _n2, 
						const T& val)
{
	Destroy();
	Create(_n1,_n2,val);
}

template <class T> void Matrix_repg<T>::Create(size_t _n1, size_t _n2, 
						const T& val)
{
//	if (_n1<1||_n2<1) throw Matrix_exception
//		("Matrix Create/Resize called with invalid dimensions!\n");
	n1=_n1; n2=_n2;
	//
	m =new Matrix_row<T>[n1];
	for (size_t i=0;i<n1;i++)
	{
		m[i].m=new T[n2];
		m[i].n2=n2;
		for (size_t j=0;j<n2;j++) m[i][j]=val;
	}
}

template <class T> void Matrix_repg<T>::Set(const T& val)
{
	for (size_t i=0;i<n1;i++)
		for (size_t j=0;j<n2;j++)
			m[i][j]=val;
}
		

template <class T> Matrix_repg<T>::Matrix_repg(const Matrix_repg<T>& c)
{
	n1=0; n2=0; m=0; ncopies=0;
	*this=c;
}

template <class T> bool Matrix_repg<T>::operator=(const Matrix_repg<T>& c)
{
	if (this!=&c)
	{
		Destroy();
		Create(c.n1,c.n2);
		Assign(c);
	}
	return true;
}

template<class T> void Matrix_repg<T>::Assign(const Matrix_repg<T>& c)
{
	for (size_t i=0;i<n1;i++)
		for (size_t j=0;j<n2;j++)
			m[i][j]=c.m[i][j];
}

template<class T> long Matrix_repg<T>::GetCopies() const
{
	return ncopies;
}
	
template<class T> long Matrix_repg<T>::IncCopies() 
{
	return ++ncopies;
}

template<class T> long Matrix_repg<T>::DecCopies() 
{
	return --ncopies;
}


template <class T> Matrix_repg<T>::Matrix_repg(size_t _n1,size_t _n2)
{
	n1=0; n2=0; ncopies=0; m=0; Create(_n1,_n2);
}

/////////////////////////////---- MATRIX
template <class T> Matrix<T>::Matrix()
{
	m=new Matrix_repg<T>;
	m->IncCopies();
}

template <class T> Matrix<T>::Matrix(size_t _n1, size_t _n2)
{
	m=new Matrix_repg<T>(_n1,_n2);
	m->IncCopies();
}

template <class T> Matrix<T>::Matrix(const Matrix<T>& _m)
{
	m=_m.m;
	m->IncCopies();
}

template <class T> Matrix<T>::~Matrix()
{
	m->DecCopies();
	if (m->GetCopies()==0) delete m; 
	m=0;
}
	
	
template <class T> inline void Matrix<T>::Clone_representation()
{
	if (m->GetCopies()==1) return;
	//Else copy it
	Matrix_repg<T>* old=m;
	m=0;
	m=new Matrix_repg<T>(*old); m->IncCopies();
	old->DecCopies(); //Unregister
	old=0;
}

template <class T> bool Matrix<T>::Resize(size_t _n1, size_t _n2, const T& val)
{
	Clone_representation();
	m->Resize(_n1,_n2,val);
	return true;
}
	

template <class T> const Matrix_row<T>& 
		Matrix<T>::operator[] (const size_t& j) const
{
	return (*m)[j];
}

template <class T> Matrix_row<T>& 
		Matrix<T>::operator[] (const size_t& j) 
{
	Clone_representation(); //Passing out reference to a non-constant
				//object
	return (*m)[j];
}


template <class T> T& 
		Matrix<T>::operator() (const size_t& i,const size_t& j) 
{
	Clone_representation(); //Passing out reference to a non-constant
				//object
	return (*m)(i,j);
}


template<class T> bool Matrix<T>::operator==(const Matrix<T>& c) const
{
	if (m->n1!=c.m->n1||m->n2!=c.m->n2) return false;
	for (size_t i=0;i<m->n1;i++)
		for (size_t j=0;j<m->n2;j++)
			if ((*this)[i][j]!=c[i][j]) return false;
	return true;
}

template<class T> bool Matrix<T>::operator!=(const Matrix<T>& c) const
{
	return !(*this==c);
}

template <class T> void Matrix<T>::Unregister()
{
	m->DecCopies();
	if (m->GetCopies()==0) delete m; 
	m=0;
}

template <class T> Matrix<T>& Matrix<T>::operator=(const Matrix<T>& c)
{
	Unregister();
	m=c.m; m->IncCopies();
	return *this;
}

template <class T> void Matrix<T>::Dump(std::ostream& out) const
{
	out<<std::endl;
	for (size_t i=0;i<m->n1;i++)
	{
		for (size_t j=0;j<m->n2;j++)
			out <<m->m[i][j]<<"   ";
		out<<std::endl;
	}
}

template<class T> Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& c) 
{
	if (m->n1!=c.m->n1||m->n2!=c.m->n2) throw
		Matrix_exception
		("Operation += applied to a matrix with different dims!\n");
	Clone_representation();
	for (size_t i=0;i<m->n1;i++)
		for (size_t j=0;j<m->n2;j++)
			m->m[i][j]+=c[i][j]; //It is const so the Clone_rep
					     //method will not be invoked!
	return *this;
}

template<class T> Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& c) 
{
	if (m->n1!=c.m->n1||m->n2!=c.m->n2) throw
		Matrix_exception
		("Operation += applied to a matrix with different dims!\n");
	Clone_representation();
	for (size_t i=0;i<m->n1;i++)
		for (size_t j=0;j<m->n2;j++)
			m->m[i][j]-=c[i][j]; //It is const so the Clone_rep
					     //method will not be invoked!
	return *this;
}


template<class T> Matrix<T>& Matrix<T>::operator=(const T& val)
{
	Clone_representation();
	m->Set(val); return *this;
}

template<class T> Matrix<T>& Matrix<T>::operator+=(const T& c) 
{
	Clone_representation();
	for (size_t i=0;i<m->n1;i++)
		for (size_t j=0;j<m->n2;j++)
			m->m[i][j]+=c;
	return *this;
}

template<class T> Matrix<T>& Matrix<T>::operator-=(const T& c) 
{
	Clone_representation();
	for (size_t i=0;i<m->n1;i++)
		for (size_t j=0;j<m->n2;j++)
			m->m[i][j]-=c;
	return *this;
}

template<class T> Matrix<T>& Matrix<T>::operator*=(const T& c) 
{
	Clone_representation();
	for (size_t i=0;i<m->n1;i++)
		for (size_t j=0;j<m->n2;j++)
			m->m[i][j]*=c;
	return *this;
}

template<class T> Matrix<T>& Matrix<T>::operator/=(const T& c) 
{
	Clone_representation();
	for (size_t i=0;i<m->n1;i++)
		for (size_t j=0;j<m->n2;j++)
			m->m[i][j]/=c;
	return *this;
}

template <class T> Matrix<T> Matrix<T>::Transpose() const
{
	Matrix<T> temp;
	temp.Resize(m->n2,m->n1);
	for (size_t i=0;i<m->n1;i++)
		for (size_t j=0;j<m->n2;j++)
			temp.m->m[j][i]=m->m[i][j];
	return temp;
}
	


template <class T> Matrix<T> operator+ 
	(const Matrix<T>& A, const Matrix<T>& B)
{
	if ( A.dim1()!=B.dim1() ||
	     A.dim2()!=B.dim2() ) throw Matrix_exception
		("Attempted to add incomatible matrices!\n");
	const size_t n1=A.dim1();
	const size_t n2=A.dim2();
	Matrix<T> C;
	C.Resize(n1,n2);
	for (size_t i=0;i<n1;i++)
		for (size_t j=0;j<n2;j++)
			C[i][j]=A[i][j]+B[i][j];
	return C;
}
			

template <class T> Matrix<T> operator- 
	(const Matrix<T>& A, const Matrix<T>& B)
{
	if ( A.dim1()!=B.dim1() ||
	     A.dim2()!=B.dim2() ) throw Matrix_exception
		("Attempted to add incomatible matrices!\n");
	const size_t n1=A.dim1();
	const size_t n2=A.dim2();
	Matrix<T> C;
	C.Resize(n1,n2);
	for (size_t i=0;i<n1;i++)
		for (size_t j=0;j<n2;j++)
			C[i][j]=A[i][j]-B[i][j];
	return C;
}

template <class T> Matrix<T> operator*
	(const Matrix<T>& A, const Matrix<T>& B)
{
	if ( A.dim2()!=B.dim1()) throw Matrix_exception
		("Attempted to add incomatible matrices!\n");
	const size_t n1=A.dim1();
	const size_t n2=B.dim2();
	const size_t nc=A.dim2();
	Matrix<T> C;
	C.Resize(n1,n2);
	T sum;
	for (size_t i=0;i<n1;i++)
		for (size_t j=0;j<n2;j++)
		{
			sum=0;
			for (size_t k=0;k<nc;k++)
				sum+=A[i][k]*B[k][j];
			C[i][j]=sum;
		}
	return C;
}
		

#ifndef INCLUDE_TEMPLATES_ONLY
//Explicit instantiation
//double
template Matrix<double> operator+(const Matrix<double>&,const Matrix<double>&);
template Matrix<double> operator-(const Matrix<double>&,const Matrix<double>&);
template Matrix<double> operator*(const Matrix<double>&,const Matrix<double>&);
template class Matrix<double>;
//int
template Matrix<int> operator+(const Matrix<int>&,const Matrix<int>&);
template Matrix<int> operator-(const Matrix<int>&,const Matrix<int>&);
template Matrix<int> operator*(const Matrix<int>&,const Matrix<int>&);
template class Matrix<int>;
//long
template Matrix<long> operator+(const Matrix<long>&,
				    const Matrix<long>&);
template Matrix<long> operator-(const Matrix<long>&,
				    const Matrix<long>&);
template Matrix<long> operator*(const Matrix<long>&,
				    const Matrix<long>&);
template class Matrix<long>;
#endif
	


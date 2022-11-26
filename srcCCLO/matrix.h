/////////////////////////////////////////////////////////////////////////////
/////////////			MATRIX Class			/////////////
//
//
//   09 Jan 2006  Ioannis C. Kampolis
//
//
//	- New matrix class with representation handle that
//	  supports copy_on_write
/////////////////////////////////////////////////////////////////////////////

#ifndef MATRIX_H
#define MATRIX_H

#include "util.h"
#include <iostream>

class Matrix_exception: public Default_exception
{
	public:
		explicit Matrix_exception(const estring& _msg)
			: Default_exception(""), msg(_msg) {}
		virtual ~Matrix_exception() {}
		virtual const char* what() const { return join_str(label,msg);}
	private:
		static estring label;
		estring msg;
};


//////////////////////////////////////////////////////////////////////////////
//                            MATRIX

// All the declarations -- for friendships
template <class T> class Matrix;
template <class T> class Matrix_repg;
template <class T> class Matrix_row;
template <class T> Matrix<T> operator+ (const Matrix<T>&, const Matrix<T>&);
template <class T> Matrix<T> operator- (const Matrix<T>&, const Matrix<T>&);
template <class T> Matrix<T> operator* (const Matrix<T>&, const Matrix<T>&);

///////////////////////////   MATRIX REPRESENATION
template <class T> class Matrix_row
{
	//
	//-Friends
	//
	friend class Matrix<T>;
	friend class Matrix_repg<T>;
	friend Matrix<T> operator+<> (const Matrix<T>&, const Matrix<T>&);
	friend Matrix<T> operator-<> (const Matrix<T>&, const Matrix<T>&);
	friend Matrix<T> operator*<> (const Matrix<T>&, const Matrix<T>&);
	
	//==============================================
	protected: // Members
	//==============================================
	//
	T* m;
	size_t n2;
	//
	//==============================================
	public: //Constructor/Destructor only
	//==============================================
	//
	virtual ~Matrix_row() {delete m; m=0;}
	Matrix_row() {m=0;}
	Matrix_row(const Matrix_row&) ;
	Matrix_row(const size_t) ;
	//
	//==============================================
	public:	// METHODS
	//==============================================
	//
	virtual T& operator[](const size_t&);
	virtual const T& operator[](const size_t&) const;
	
};

template <class T> class Matrix_repg
{
	//
	//-Friends
	//
	friend class Matrix<T>;
	friend Matrix<T> operator+<> (const Matrix<T>&, const Matrix<T>&);
	friend Matrix<T> operator-<> (const Matrix<T>&, const Matrix<T>&);
	friend Matrix<T> operator*<> (const Matrix<T>&, const Matrix<T>&);
	
	//==============================================
	protected: // Members
	//==============================================
	Matrix_row<T>* m;
	size_t n1,n2;
	private:
	long ncopies;
	//==============================================
	public: //Constructor/Destructor only
	//==============================================
	////
	Matrix_repg(){n1=n2=0; m=0; ncopies=0;} 
	Matrix_repg(size_t,size_t);
	Matrix_repg(const Matrix_repg<T>&);
	virtual ~Matrix_repg() {Destroy();}
	/////
	//==============================================
	public:	// METHODS
	//==============================================
	//
	//- Access
	//
	virtual T& operator()(const size_t& i, const size_t& j) const;
	virtual Matrix_row<T>& operator[](const size_t& j);
	virtual const Matrix_row<T>& operator[](const size_t& j) const;
	//
	//- Manipulation
	//
	virtual void Resize(size_t _n1, size_t _n2, const T& val=T(0));
	virtual size_t dim1() {return n1;}
	virtual size_t dim2() {return n2;}
	//
	//- Operators
	//
	bool operator=(const Matrix_repg<T>&);
	//==============================================
	protected: // Protected Helper Methods
	//==============================================
	virtual void Destroy();
	virtual void Create(size_t _n1, size_t _n2, const T& val=T(0));
	virtual void Set(const T& val=0);
	virtual void Assign(const Matrix_repg<T>&);
	long GetCopies() const;
	long IncCopies();
	long DecCopies();
	//virtual T** GetPointer(){return m;}
		
		
};



///////////////////////////////  MATRIX
//
//
//
// The Matrix
template <class T> class Matrix
{
	//Friends
	friend Matrix<T> operator+<> (const Matrix<T>&, const Matrix<T>&);
	friend Matrix<T> operator-<> (const Matrix<T>&, const Matrix<T>&);
	friend Matrix<T> operator*<> (const Matrix<T>&, const Matrix<T>&);
	
	protected:
	//Protected members
	Matrix_repg<T>* m;
	//==============================================
	public: //Constructor/Destructor only
	//==============================================
	Matrix(); 
	Matrix(size_t _n1, size_t _n2);
	Matrix(const Matrix& _m);
	virtual ~Matrix();
	//==============================================
	public: //METHODS
	//==============================================	
	//
	//- Dimension handling
	//
	bool Resize(size_t _n1, size_t _n2, const T& c=T(0));
	size_t dim1() const { return m->n1; }
	size_t dim2() const { return m->n2; }
	//
	//- Access to elements
	//
	Matrix_row<T>& operator[](const size_t&);  
	T& operator()(const size_t&,const size_t&);  
	const Matrix_row<T>& operator[](const size_t&) const;
	void Dump(std::ostream&) const;
	long RepresentationCopies(){return m->GetCopies();}
	//T** Array2D(){return m->m;}
	//==============================================
	public:	// Operators
	//==============================================	
	//
	bool operator==(const Matrix<T>&) const;
	bool operator!=(const Matrix<T>&) const;
	//- operations with Matrices
	//
	Matrix<T>& operator=(const Matrix<T>&);	//set equal to another matrix
	Matrix<T>& operator+=(const Matrix&);
	Matrix<T>& operator-=(const Matrix&);
	//
	//- operations with values
	//
	Matrix<T>& operator=(const T&);		//set all elements = val
	Matrix<T>& operator+=(const T&);
	Matrix<T>& operator-=(const T&);
	Matrix<T>& operator*=(const T&);
	Matrix<T>& operator/=(const T&);
	//
	//- Various
	//
	Matrix<T> Transpose() const;
	//
	//==============================================
	private: // METHODS -- Private -- mostly helpers
	//==============================================
	void Clone_representation();
	void Unregister();
	
};





#endif

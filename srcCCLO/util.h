/////////////////////////////////////////////////////////////////////////////
//////////                      Utilities                         ///////////
//
// Generic templates
// Global exceptions, assertions
//
// 23 Apr 2002, agiotis
// 10 Oct 2005, gkampolis
/////////////////////////////////////////////////////////////////////////////


#ifndef UTIL_H
#define UTIL_H

#include <valarray>
#include <algorithm>
#include <string.h>
#include <vector>
#include <stdexcept>
#include "mymath.h"

using namespace std;

class Default_exception;

// Range checked Vector 
template <class T> class Vector : public valarray<T>
{
public:
	Vector():valarray<T>() {}
	Vector(int s) : valarray<T>(s) {}
	Vector(T t,int sz): valarray<T>(t,sz) {};
	
	Vector<T>& operator=(const T& _Val)
	{
		valarray<T>::operator=(_Val);
		return *this;
	}

	T& operator[](int i) 
	{	if(i<0 || (static_cast<unsigned int>(i))>=this->size())
			throw Default_exception("Caught a Vector range error");
		return valarray<T>::operator[](i);
	}
	T operator[](int i) const
	{	if(i<0 || (static_cast<unsigned int>(i))>=this->size())
			throw Default_exception("Caught a Vector range error");
		return valarray<T>::operator[](i);
	}
	Vector<T>& assign(const valarray<T>& _c)
	{
		valarray<T>::resize(_c.size());
		for (size_t i=0;i<_c.size();i++)
			(*this)[i]=_c[i];
		return *this;
	}

	Vector<T>& assign(const std::vector<T>& _c)
	{
		valarray<T>::resize(_c.size());
		for (size_t i=0;i<_c.size();i++)
			(*this)[i]=_c[i];
		return *this;
	}

};




/////////////////////
// Default Exception
typedef std::string estring;

class Default_exception //:public exception// General exception
{ 
public:
	Default_exception()	{ msg = "Unknown default exception";}
	Default_exception(const estring& _msg): msg(_msg) {}
	virtual ~Default_exception() {}
	virtual const char *what() const { return msg.c_str(); }
protected:
	const char* join_str(const estring& a, const estring& b) const
	{
		estring err=a+b;
		char *p=new char[err.size()+1];
		strcpy(p,err.c_str()); //deprecated
		return p;
	}
		
private:
	estring msg;
};

/////////////
// Assertion
#ifdef NDEBUG
const bool ARG_CHECK = false;	// we are not debugging: disable checks
#else
const bool ARG_CHECK = true;	// we are debugging
#endif

class Bad_arg{};	// thrown as exception by Assert
template<class X, class A> inline void Assert(A assertion)
{
	if(!assertion) throw X();
}




template<class T> T mul(const valarray<T>&v1, const valarray<T>& v2)
{
	T res = 0;
	for(int i=0; i<v1.size(); ++i)
		res += v1[i]*v2[i];
	return res;
}


class indexx
{
public:
	void operator() (const int n, const vector<double>& arrin, vector<int>& indx);
};


#endif



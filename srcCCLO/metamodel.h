/////////////////////////////////////////////////////////////////////////////
///////////////           METAMODELS				/////////////
//
//
//
//	Generic Abstract Class
//	I. Kampolis  04-01-2005
//
//
/////////////////////////////////////////////////////////////////////////////

#ifndef METAMODEL_H
#define METAMODEL_H

#include "util.h"
#include "matrix.h"
#include <string>

namespace METAMODELS
{
	enum  Generalization {Low=0, Normal, High};
	const Generalization Lowg=Low;
	const Generalization Normalg=Normal;
	const Generalization Highg=High;

	class Mmodel_exception: public Default_exception
	{
		public:
			explicit Mmodel_exception(const estring& _msg)
				: Default_exception(""), msg(_msg) {}
			virtual ~Mmodel_exception() {}
			virtual const char* what() const 
				{ return join_str(label,msg); }
		private:
			static estring label;
			estring msg;
	};

  	class CMetamodel
  	{
		public: //Virtual Destructor
			virtual ~CMetamodel(){};
			CMetamodel(const std::string p): mm_desc(p)
			{Reset();}
		protected:
			Matrix<double> x;		//Input Patterns
			Matrix<double> y;		//Output values
			const std::string mm_desc;	//Metamodel Description
			bool is_built;
			bool is_trained;
  		public:
			virtual void Reset() {is_built=is_trained=false;}
			virtual bool Build( const Matrix<double>&,
					    const Matrix<double>& )=0;
			virtual bool Train()=0;
			virtual bool Eval ( const Vector<double>&,
					          Vector<double>&)=0;
			virtual bool Evalgrad( const Vector<double>&,
						     Matrix<double>&)=0;
			virtual std::string Description()
				{return mm_desc;}
			virtual bool SetParams( const Vector<double>& )=0;
			
	};
};

#endif

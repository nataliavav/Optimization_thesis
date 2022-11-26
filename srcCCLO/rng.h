/////////////////////////////////////////////////////////////////////////////
//////////          Random Number Generator                       ///////////
//
// 20 Jun 2006, gkampolis
/////////////////////////////////////////////////////////////////////////////

#ifndef RNG_H_MY
#define RNG_H_MY

#include <iostream>

const char RNG_UNI_VER[]="RNG_UNIv1.0";
const char RNG_GAU_VER[]="RNG_GAUv1.0";

class RandUnif
{
	public:
		RandUnif(); 
		virtual ~RandUnif(); 
	private:
		void Initialize();
	public:
		void Reset();
		double operator()();
		int operator()(const int min, const int max);
		long Calls() const;
		void Save(std::ostream&) const;
		bool Restore(std::istream&);
	private:
		long callN;	// Mainly to debug
		//
		int inext,inextp;
		long *ma;
		//
		const long MBIG	,MSEED,MZ;
		const double FAC;
};

class RandGauss
{
	private:
		RandUnif& myrand;
		int iset;
		double gset;
	public:
		RandGauss(RandUnif&);
		double operator()();
		void Save(std::ostream&) const;
		bool Restore(std::istream&);
};


#endif

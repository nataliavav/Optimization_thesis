#include "rng.h"
#include "util.h"
#include "filesystem.h"


RandUnif::RandUnif () : MBIG(1000000000), MSEED(161803398),MZ(0),
			FAC(1.0/double(MBIG))
{
	ma=new long[56]; //Allocate memory
	Reset();
}

RandUnif::~RandUnif()
{
	delete[] ma; ma=0;
}

void RandUnif::Initialize()
{
	
	long mj,mk;
	int i,ii,k;
	mj=MSEED -1;	//-( *idum < 0 ? -*idum :*idum);
	// Initialize ma[55] using the seed idum and the large number MSEED
	mj %= MBIG;
	ma[55]=mj;
	mk=1;
	for (i=1;i<=54; i++ ) 
	{
		// Now initialize the rest of the
		// table,in a slightly random order,
		// with numbers not especially random
		ii=(21*i) % 55;
		ma[ii]=mk;  
		mk=mj-mk; 
		if (mk < MZ) mk += MBIG;
		mj=ma[ii];
	} 
	for (k=1;k<=4;k ++)	// We randomize them by
				//   "warming up the generator"
		for (i=1;i<=55 ;i++) 
		{
			ma[i] -= ma[1+(i+30) % 55];
			if (ma[i] < MZ) ma[i] += MBIG;
		}
	inext=0;	// Prepare indices for our first generated number
	inextp=31;	// The constant 31 is special; see Knuth
}

void RandUnif::Reset()
{
	callN=0; Initialize();
}

double RandUnif::operator()()
{
	// Returns a uniform random deviate between 0.0 and 1.0.
//	if(++callN >= LONG_MAX)
	if(++callN >= 10000000)
		throw Default_exception(
		"The internal counter of RandUnif overflowed");
	long mj;

	// Here is where we start, except on initialization.
	if (++inext == 56) inext=1;	// Increment inext and inextp,
					//   wrapping around 56 to 1.
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];	// Generate a new random number 
					//    subtractively
	if (mj < MZ) mj += MBIG;	// Besure that it is in range.
	ma[inext]=mj;			// Store it,
	return mj*FAC;			// and output the derived uniform 
					//    deviate. 
}


int RandUnif::operator()(const int min, const int max) 
{
	if( (max-min) <=1 ) return min;
	else return static_cast<int>(min+operator()()*(max-min));
}

long RandUnif::Calls() const
{
	return callN;
}
	
void RandUnif::Save(std::ostream& out) const
{
	const int nsize=sizeof(RNG_UNI_VER)/sizeof(char);
	out.write(RNG_UNI_VER,nsize);
	//
	bin_write(out,callN);
	bin_write(out,inext);
	bin_write(out,inextp);
	for (int i=0;i<56;i++) bin_write(out,ma[i]);
}

bool RandUnif::Restore(std::istream& in)
{
	const int _len=sizeof(RNG_UNI_VER)/sizeof(char);
	char _VER[_len];
	in.read(_VER,_len);
	if (strcmp(_VER,RNG_UNI_VER)) return false;
	//
	bin_read(in,callN);
	bin_read(in,inext);
	bin_read(in,inextp);
	for (int i=0;i<56;i++) bin_read(in,ma[i]);
	return istream_operational(in);
}

//////////
//

RandGauss::RandGauss(RandUnif& _myrand) :myrand(_myrand)
{
	iset=0;
}

double RandGauss::operator()()
{
	// Returns a gauss distributed N(0,1) random number

	// Computational time: more than twice that of uniform !

	// From NUMERICAL RECIPES
	double fac,rsq,v1,v2;
   
	if  (iset == 0) 
	{
		do {
			v1=2.0*myrand()-1.0;
			v2=2.0*myrand()-1.0;

			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} 
	else 
	{
		iset=0;
		return gset;
	}
}

void RandGauss::Save(std::ostream& out) const
{
	bin_write(out,RNG_GAU_VER);
	//
	bin_write(out,iset);
	bin_write(out,gset);
}

bool RandGauss::Restore(std::istream& in)
{
	const int _len=sizeof(RNG_GAU_VER)/sizeof(char);
	char _VER[_len];
	bin_read(in,_VER);
	if (strcmp(_VER,RNG_GAU_VER)) return false;
	//
	bin_read(in,iset);
	bin_read(in,gset);
	return istream_operational(in);
}

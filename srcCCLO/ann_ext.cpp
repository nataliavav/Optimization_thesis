#include "ann_ext.h"
#include <fstream>
#include <string>

using namespace METAMODELS;

void ANNExternal::Reset()
{
	remove(dbfile);
	remove(patfile);
	remove(resfile);
}

bool ANNExternal::Build(TrainingPattern *_tp, const int _npat,
			const int _varN, const int _objN)
{
	npat=_npat;
	varN=_varN;
	objN=_objN;
	
	std::ofstream out (dbfile,std::ios::out);
	out <<_npat<<"  "<<_varN<<"  "<<_objN<<"\n";
	out.precision(14);
	for (int i=0;i<_npat;i++)
	{
		for (int j=0;j<_varN;j++)out<<_tp[i].var[j]<<"  ";
		out <<"  ";
		for (int j=0;j<_objN;j++)out<<_tp[i].obj[j]<<"  ";
		out <<"\n";
	}
	out.close();
	return true;
}

bool ANNExternal::CheckTrainingSet() const
{
	return true;
}

bool ANNExternal::Predict(const Vector<double>& dat,
			  Vector<double>& res)
{
	if (dat.size()!=varN) return false;
	if (res.size()!=objN) return false;
	///////////////////////////////////////////////
	std::ofstream out (patfile,std::ios::out);
	out <<varN<<"\n";
	out.precision(14);
	for (int i=0;i<varN;i++) out<<dat[i]<<"\n";
	out.close();
	///////////////////////////////////////////////
	SystemCall(std::string(usefile));
	///////////////////////////////////////////////
	std::ifstream in(resfile,std::ios::in);
	for (int i=0;i<objN;i++) in>>res[i];
	if (!istream_operational(in))
	{
		in.close();
		res=dmax;
		return false;
	}
	in.close();
	return true;
}
	
bool ANNExternal::PredictGrad(const Vector<double>& dat,
			      Matrix<double>& res)
{
	return false;
}

bool ANNExternal::Train()
{
	SystemCall(std::string(trfile));
	return true;
}
			      

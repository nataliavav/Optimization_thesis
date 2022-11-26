#include "param.h"
#include "util.h"
#include "mymath.h"
#include "filesystem.h"

Parameterization::Parameterization()
{
	Reset();
}

Parameterization::~Parameterization(){}

void Parameterization::Reset()
{
	dvars.resize(0);
	nvar=ncon=0;
	con.resize(0);
}

int Parameterization::getVarN() const
{
	return nvar;
}

int Parameterization::getConN() const
{
	return ncon;
}

void Parameterization::addVariable(const double min, const double max,
				   const int bits)
{
	VariableInfo nv;
	nv.min=min; nv.max=max; nv.ext=max-min; nv.bits=bits;
	dvars.push_back(nv);
	nvar++;
}

void Parameterization::addConstraint(const double accv, const double fval,
					const double factor)
{
	if (fval<=accv) throw Default_exception
		("addConstraint: Nominal Thr.>=Relaxed Thr. error\n");
	if (factor<=0) throw Default_exception
		("addConstraint: Constraint ampl. factor <=0 error\n");

	ConstraintInfo newcon;
	newcon.acceptval=accv;
	newcon.failval=fval;
	newcon.factor=factor;
	con.push_back(newcon);
	ncon++;
}


bool Parameterization::ReadParameterization(std::istream& in)
{
	int nvars;
	in>>nvars; istream_nl(in);
	int bits; double b1,b2; 
	for (int i=0;i<nvars;i++)
	{
		in>>bits>>b1>>b2; istream_nl(in);
		addVariable(b1,b2,bits);
	}
	return true; 
}

bool Parameterization::ReadConstraints(std::istream& in)
{
	int nc;
	in>> nc; istream_nl(in);
	double v1,v2,v3;
	for (int i=0;i<nc;i++)
	{
		in>>v1>>v2>>v3; istream_nl(in);
		addConstraint(v1,v2,v3);
	}
	return true;
	//return istream_operational(in); //if other data follow
}

double Parameterization::ComputePenalty(const vector<double>& c,
		 			bool& pen, bool& death) const
{
	if (ncon==0) 
	{	
		death=pen=false;
		return 0;
	}

	if (int(c.size())!=ncon) throw Default_exception(
			"EvaluatePenalty: Invalid Constraints no\n");
	double p1=1.0; double s1=-.5;
	death=pen=false;
	for (int i=0;i<ncon;i++)
	{
		const double den=(con[i].failval-con[i].acceptval);
		const double cdr=(c[i]-con[i].acceptval)/den;
		s1+=Hc(cdr,0);
		p1*=exp(con[i].factor*cdr*Hc(cdr,0));
		if (c[i]>con[i].failval) death=true;
	}
	if (s1>0.1) pen=true;
	return Hc(s1,0)*p1;
}

void Parameterization::IsPenalized(const vector<double>& c,
				   bool& penalized, bool &death) const
{
	ComputePenalty(c,penalized,death);
}

void Parameterization::Penalize(const vector<double>& c, 
		const vector<double>& o, vector<double>& po,
		bool& penalized, bool& death)
{
	const double pen=ComputePenalty(c, penalized, death);
	po.resize(o.size());
	const double oblow=dmax/1000.0;
	for (int i=0;i<int(o.size());i++)
	{
		po[i]=o[i]+pen;
		if (penalized && o[i]>oblow && po[i]<o[i]) po[i]=dmax;
	}
}


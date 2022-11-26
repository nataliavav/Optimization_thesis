#ifndef PARAM_H
#define PARAM_H

#include<iostream>
#include<vector>

struct VariableInfo
{
	double min;	//lower bound
	double max;	//upper bound
	double ext;	//extension 
	int bits;	//bits 
};

struct ConstraintInfo
{
	double acceptval;	//Maximum accepted value
	double failval;		//Death Penalty Value
	double factor;		//Exponential factor
};

class Parameterization
{
	public: // --- Constructor - Destructor 
		Parameterization();
		virtual ~Parameterization();
	public:
	std::vector<VariableInfo>   dvars; 		
	std::vector<ConstraintInfo> con; // Constraints

	protected:
	int nvar;	// Number of Design Variables
	int ncon;	// Number of Constraints

	void Reset();
	//
	public:
	int getVarN() const;
	void addVariable(const double min, const double max,
			 const int bits);
	bool ReadParameterization(std::istream& in);
	//
	int getConN() const ;
	void addConstraint(const double, const double,
			   const double);
	bool ReadConstraints(std::istream& in);

	double ComputePenalty(const std::vector<double>&,
			      bool&,bool&) const;

	void IsPenalized(const std::vector<double>&,
			 bool&,bool&) const;

	void Penalize(const std::vector<double>& ,   // constraint 
		      const std::vector<double>& o,  // objective
		      std::vector<double>& po,	     // penalized obj
		      bool& penalized, bool& death); //penalized, Ded
	
};

#endif


#ifndef INDIV_H
#define INDIV_H

#include <vector>

class Individual
{
	public:
	std::vector <double> var;	//Variables
	std::vector <double> obj;	//Objectives
	std::vector <double> pobj;	//Penalized Objectives
	std::vector <double> con;	//Constraints
	double phi;		//SPEA fitness assignement
	double raw;		//SPEA Raw assignement
	int flag;		//flag
	int flagEval;		//flag
	void copy(Individual ind)
	{
		phi	 = ind.phi;
		raw	 = ind.raw;
		flag 	 = ind.flag;
		flagEval = ind.flagEval;
		var = ind.var;
		obj = ind.obj;
		pobj = ind.pobj;
		con = ind.con;
	};
};

#endif

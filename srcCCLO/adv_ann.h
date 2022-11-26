

#ifndef ADV_ANN_H
#define ADV_ANN_H

#include "util.h"
#include <set>
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <queue>
#include "matrix.h"


namespace METAMODELS
{

	enum ActivationType {GaussianActivationTypeI=0,
			     MultiQuadActivationTypeI,
			     InvMultiQuadActivationTypeI};
	const ActivationType GaussianActivationType=GaussianActivationTypeI;
	const ActivationType MultiQuadActivationType=MultiQuadActivationTypeI;
	const ActivationType InvMultiQuadActivationType=
					InvMultiQuadActivationTypeI;



	class ActivationFunction
	{
		public:
			virtual const double activ  (const double, 
				  	             const double) const=0;
			virtual const double dactiv (const double, 
				  	             const double) const=0;
			virtual const double ddactiv(const double,
						     const double) const=0;
	};

	class GaussianActivation : public ActivationFunction
	{	
		public:
			virtual const double activ  (const double,
						     const double) const;
			virtual const double dactiv (const double,
						     const double) const;
			virtual const double ddactiv(const double,
						     const double) const;
	};
		

	class MultiQuadActivation : public ActivationFunction
	{
		public:
			virtual const double activ  (const double,
						     const double) const;
			virtual const double dactiv (const double,
						     const double) const;
			virtual const double ddactiv(const double,
						     const double) const;
	};
		
	class InvMultiQuadActivation : public ActivationFunction
	{
		public:
			virtual const double activ  (const double,
						     const double) const;
			virtual const double dactiv (const double,
						     const double) const;
			virtual const double ddactiv(const double,
						     const double) const;
	};

	
	class Neuron
	{
		protected:
			Vector<double> c;	//The center of the neuron
		public:
			Neuron();
			Neuron(const Neuron&);
			Neuron(const Vector<double>&);
			virtual ~Neuron();
		public:
			void SetEqualTo(const Neuron&);
			const int GetDimension() const;
			void SetDimension(const int);
			Vector<double>& GetCenter();
			const Vector<double>& GetCenterConst() const;
	};

	class RBFNeuron : virtual public Neuron
	{
		protected:
			double r;
			ActivationFunction* a;
		public:
			RBFNeuron();
			RBFNeuron(const RBFNeuron&);
			RBFNeuron(const Vector<double>&);
			RBFNeuron(const Vector<double>&, const double);
			virtual ~RBFNeuron();
		public:
			const double Activate(const Vector<double>&) const;
			const double ActivateDer(const Vector<double>&,
						 Vector<double>& ) const;
			const double ActivateIF(const Vector<double>&,
					        const Matrix<double>&,
						const int) const;
			const double ActivateDerIF(const Vector<double>&,
					           const Matrix<double>&,
						   const int,
 						   Vector<double>& ) const;
			const double GetR() const;
			void SetR(const double);
	};

	class GNGNeuron: virtual public Neuron
	{
		protected:
			std::set<int> connect;
		public:
			GNGNeuron();
			GNGNeuron(const GNGNeuron&);
			GNGNeuron(const Vector<double>&);
			virtual ~GNGNeuron();
		public:
			std::set<int>& GetConnectivity();
			const std::set<int>& GetConnectivityConst() const;
	};

	struct TrainingPattern
	{
		Vector<double> var;
		Vector<double> obj;
		bool bfailed;
	};


	class AdvancedANN
	{
		protected:
			int varN;
			int objN;
			int npat;
			TrainingPattern* tp;
			bool is_built;
			bool is_trained;
		public:
			virtual ~AdvancedANN();
			virtual void Reset()=0;
		public:
			virtual bool Build(TrainingPattern*, const int,
					   const int, const int)=0;
			virtual bool CheckTrainingSet() const=0;
			virtual bool Predict(const Vector<double>&,
					     Vector<double>&)=0;
			virtual bool PredictGrad(const Vector<double>&,
					     Matrix<double>&)=0;
			virtual bool Train()=0;
	};


	//Edge of an MST
	struct edge_t 
	{
	
		edge_t(const int src_, const int dest_, const double& w_)
			: src(src_), dest(dest_), w(w_)
		{};
		int src;
		int dest;
		double w;
		bool operator < (const edge_t& e) const
		{
			// The smaller the w, the higher the priority:
			return w > e.w;
		}
	};


	
	class CMST
	{
		private:
			TrainingPattern* tp;
			const int npat;
			const int varN;
			const int objN;
			double d_ave;
			double d_dev;
			double** dist;
			std::vector< std::list<int> > tree;
		public:				//npat     //varN
			CMST(TrainingPattern*, const int, const int, 
						//objN
					       const int);
			virtual ~CMST();
		protected:
			const double Norm2Var(const int, const int) const;
			const double Norm2Obj(const int, const int) const;
			void CalcDistances();
			void ConstructMST();
			void WritePath(std::ostream&,const int, 
						     const int) const;
		public:
			const double GetAveDist() const;
			const double GetStDev() const;
					//   start node  % of pats to keep
			void ConditionedWalk(const int, const double,
					//    times stdev to keep
					     const double,
					//    number I kept
					     int&,
					//    the ids I kept
				             std::vector<int>&,
					//   traversed ave dist & std_dev
					     double&, double&) const;
					// start node, stream
			void WriteTree(const int, std::ostream&) const;
					// start node, fname
			void WriteTree(const int, std::string) const;
		private:
			CMST(const CMST&);	//do not allow copy-constructor
			
	};

	class ARBFN_Neuron :	virtual public RBFNeuron, 
				virtual public GNGNeuron,
				virtual public Neuron
	{
		public:
			ARBFN_Neuron ();
			ARBFN_Neuron (const ARBFN_Neuron&);
			ARBFN_Neuron (const Vector<double>&);
			ARBFN_Neuron (const Vector<double>&,const double);
			~ARBFN_Neuron();
		public:
			virtual void SetEqualTo(const ARBFN_Neuron&);
			virtual Vector<double>& GetCenter();
			virtual const Vector<double>& GetCenterConst() const;
			virtual void SetDimension(const int);
			virtual void Save(std::ostream&) const;
			virtual bool Restore(std::istream&);
			
			
	};


	
			

};
		

		



#endif

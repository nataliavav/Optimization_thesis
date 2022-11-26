//////////////////////////////////////////////////////////////////////////////
////////////  Advanced RBFN Network (ARBFN)		//////////////////////
//
//
//
//	 3 Oct 2006 Ioannis C. Kampolis  
//
//
//////////////////////////////////////////////////////////////////////////////
//
//


#ifndef A_RBFN_H
#define A_RBFN_H

#include "adv_ann.h"
#include "gng.h"
#include "matrix.h"
#include <sstream>

 
namespace METAMODELS
{

	const double ERR_BLOWUP = 1.e10;
	//
	class ARBFN: virtual public GNGCluster, virtual public AdvancedANN
	{
		protected:
			ARBFN_Neuron** neurons;		//The neurons	
			int n_neuron;			//Number of neurons
			//	
			bool b_interp;			//Inter/ Aprox
			//
			int* active;			//Active indices
			int nactive;			//Number of active
			//	
			double** w;			//Weights
			//	
			TrainingPattern* testpat;	//Test patterns
			int ntest;			//number of test pats
			//
			double test_to_tot_ratio;
			int min_cen,max_cen;
			double r_multfac;			//
			//
			double errtest_min;
			//
			double *lms_corr_trace;
			double **ls_h;
			int lms_Ntrain;
			double lms_lrfac;	 		//
			//
			ARBFN_Neuron* best_neuron;
			int idle_n_growth;                      //
			int nbest;
			GNGConfig* gngconfig;
			Matrix<double> IF;			// Imp.Factor
			bool bIF;
			std::ostringstream gngout;
		public:
			ARBFN();
			virtual ~ARBFN();
		public:
							//min cen, max cen
			void SetTrainingApproximation(const int ,const int,
							//test to total ratio
						      const double,
						      	//GNG config
						      const GNGConfig* const,
						      	//r_fac
						      const double,
						        //learning rate fac
						      const double,
						        //idle growth
						      const int);
			void SetTrainingInterpolation();
			void SetIFUsage(const bool);
			bool SetIFs(const Matrix<double>&);
			//From AdvancedANN
			virtual void Reset();
			virtual bool Build(TrainingPattern*, const int, const
					   int, const int);
			virtual bool CheckTrainingSet() const;
			virtual bool Train();
			virtual bool Predict(const Vector<double>&,
					     Vector<double>&);
			virtual bool PredictGrad(const Vector<double>&,
					     Matrix<double>&);
			virtual void SetRadius(const double);
			virtual void Save(std::ostream&) const;
			virtual bool Restore(std::istream&);
			virtual double TestError2();
		protected:
			//From GNGCluster
			virtual void Initialize();
			virtual double NeuronError( const int, const int);
			virtual void AdaptNormal (const int);
			virtual void AdaptExtended();
			virtual bool IsImproving();
			virtual void Adjust(const int, const int, const int);
			virtual GNGNeuron* GetNeuron(const int);
			virtual int* GetActiveRedirector();
			virtual int& GetActiveCounter();
		protected:
			virtual bool BuildInterpolate(TrainingPattern*,	
						const int, const int,
						const int);
			virtual bool BuildApproximate(TrainingPattern*,	
						const int, const int,
						const int);
			virtual bool TrainInterpolate();
			virtual bool TrainApproximate();
			//
			virtual void SetRadius();
			virtual bool ComputeWeights();
			virtual bool ComputeWeightsNOIF();
			virtual bool ComputeWeightsIF();
			virtual double TestError();
			virtual double norm2var(const Vector<double>&,
						const Vector<double>&) const;
			void SaveTP(std::ostream&) const;
			void SaveTestP(std::ostream&) const;
			void SaveIFs(std::ostream&) const;
			void SaveW(std::ostream&) const;
			void SaveInterp(std::ostream&) const;
			void SaveApprx(std::ostream&) const;
			bool RestoreTP(std::istream&, TrainingPattern*,
					const int, const int, const int);
			bool RestoreTestP(std::istream&);
			bool RestoreIFs(std::istream&);
			bool RestoreW(std::istream&);
			bool RestoreInterp(std::istream&);
			bool RestoreApprx(std::istream&);

	};
};


#endif

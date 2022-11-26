//////////////////////////////////////////////////////////////////////////////
////////////  Growing Neural Gas (GNG) 			//////////////////////
//
//
//
//	21 May 2005 Marios K. Karakasis  GNG Core
//	29 Sep 2006 Ioannis C. Kampolis  Compatible with EASY
//
//
//////////////////////////////////////////////////////////////////////////////


#ifndef GNG_H
#define GNG_H

#include "util.h"
#include "metamodel.h"
#include "gconfig.h"
#include <vector>
#include "rng.h"
#include "adv_ann.h"

namespace METAMODELS
{
	typedef NCONFIG::GConfig GConfig;
	//The exception
	class GNG_exception : public Default_exception
	{
		public:
			explicit GNG_exception(const estring& _msg)
				: Default_exception(""), msg(_msg) {}
			virtual ~GNG_exception() {}
			virtual const char* what() const 
				{ return join_str(label,msg); }
		private:
			static estring label;
			estring msg;
	};


	class GNGConfig : public GConfig
	{
		public:
			double mvc;	//Center movement factor
			double mvn;	//Neighbors movement factor
			int max_connect_age;	//Maximum age
			double addstep_dimmult;	//Factor for adding step
						//according to dimensionality
			double errdecr_interp;	//Error decrease factor
						//after interpolation
			double errdecr_gen;	//Gen. Error decrease factor
			double tp_close_frac;	//perc of tp close to a neuron
						//to be singled out
			double tp_close_prob;	//Probability to select a tp 
						//close to a neuron after 
						//adding it
			double tp_close_iter;	//Iterations to select tp close
						//to a neuron after adding it
						//with prob. tp_close_prob
		public:
			GNGConfig(std::ostream&);
			virtual ~GNGConfig();
		public:
			virtual bool ReadConfig(std::istream&);
			virtual bool CheckConfig();
	};

	class GNG;
	
	class GNGCluster 
	{
		friend class GNG;
		protected:
		virtual void Initialize() =0;
					// which neuron   which tp
		virtual double NeuronError( const int,     const int)=0;
					 // which pattern
		virtual void AdaptNormal( const int )=0;
		virtual void AdaptExtended()=0;
		virtual bool IsImproving()=0;
		virtual void Adjust(const int, const int, const int)=0;
		virtual GNGNeuron* GetNeuron(const int)=0;
		virtual int* GetActiveRedirector()=0;
		virtual int& GetActiveCounter()=0;
	};


	class GNG
	{
		public:
			//Constructor
		private:
			GNGCluster* cluster;	//The clustering info
			const GNGConfig* const config;	//The configuration

			struct tp_dist_t
			{
				double d;
				int i;
				bool operator< (const tp_dist_t& tp) const
				{
					return d<tp.d;
				}
			};

			RandUnif myrand;
		private:
			GNG(const GNG&);		//no copy constructor
			GNG& operator=(const GNG&);	//no assignment op.
		public:
			GNG(const GNGConfig* const _config,
			    GNGCluster* _cluster);	
		public:
			void Grow(const int Nneuron_max,
				  const int tp_Ntrain,
				  const int tp_Ntrain_max,
				  const int Nin,
				  const int Nneuron_start,
				  const int Nneuron_stop,
				  TrainingPattern* tp_train
				  );
		protected:
			const double norm2var(const int,
					      const Vector<double>&,
					      const Vector<double>&) const;
			
	};


				
};



#endif

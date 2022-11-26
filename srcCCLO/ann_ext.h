///////////////////////////////////////////////////////////////////////////////
//
//                   External Metamodel
//
//
//    Dec. 2006  I.C.Kampolis
//
//
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ANN_EXT_H
#define ANN_EXT_H

#include "adv_ann.h"
#include "filesystem.h"

namespace METAMODELS
{
	const char  dbfile[] = "meta.db";
	const char patfile[] = "meta.dat";
	const char resfile[] = "meta.res";
	const char  trfile[] = "meta_train.bat";
	const char usefile[] = "meta_use.bat";

	class ANNExternal : public AdvancedANN
	{
		public:
			virtual void Reset();
			virtual bool Build(TrainingPattern*, const int,
					   const int, const int);
			virtual bool CheckTrainingSet() const;
			virtual bool Predict(const Vector<double>&,
					     Vector<double>&);
			virtual bool PredictGrad(const Vector<double>&,
					     Matrix<double>&);
			virtual bool Train();
	};
};
#endif


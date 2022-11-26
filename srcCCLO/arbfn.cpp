#include "arbfn.h"
#include "rng.h"
#include "myblas.h"
#include "matrix.h"
#include "svd.h"
#include "filesystem.h"


using namespace METAMODELS;

ARBFN::ARBFN() 
{
	tp=0;
	npat=varN=objN=0;
	is_built=is_trained=0;
	//
	n_neuron=0;
	neurons=0;
	//
	b_interp=true;
	active=0;
	nactive=0;
	//
	w=0;
	//
	testpat=0;
	ntest=0;
	//
	test_to_tot_ratio=0.0;
	lms_corr_trace=0;
	ls_h=0;
	best_neuron=0;
	gngconfig=0;
	bIF=false;
	gngconfig=new GNGConfig(gngout);
}


ARBFN::~ARBFN()
{
	delete gngconfig; gngconfig=0;
	Reset();
}


void ARBFN::SetTrainingApproximation(const int _min_cen, const int _max_cen,
				     const double _ratio, 
				     const GNGConfig *  const _gng,
				     const double _rfac, 
				     const double _lrfac, const int idlegr)
{
	is_built=false;
	min_cen=_min_cen;
	max_cen=_max_cen;
	test_to_tot_ratio=_ratio;
	b_interp=false;
	//gngconfig=_gng;
	r_multfac=_rfac;
	lms_lrfac=_lrfac;
	idle_n_growth=idlegr;
	
	gngconfig->mvc=_gng->mvc;
	gngconfig->mvn=_gng->mvn;
	gngconfig->max_connect_age=_gng->max_connect_age;
	gngconfig->addstep_dimmult=_gng->addstep_dimmult;
	gngconfig->errdecr_interp=_gng->errdecr_interp;
	gngconfig->errdecr_gen=_gng->errdecr_gen;
	gngconfig->tp_close_frac=_gng->tp_close_frac;
	gngconfig->tp_close_prob=_gng->tp_close_prob;
	gngconfig->tp_close_iter=_gng->tp_close_iter;

}

void ARBFN::SetTrainingInterpolation()
{
	is_built=false;
	b_interp=true;
}

void ARBFN::SetIFUsage(const bool _bif)
{
	bIF=_bif;
}

bool ARBFN::SetIFs(const Matrix<double>& _IF)
{
	if (_IF.dim1()!=IF.dim1() || _IF.dim2()!=IF.dim2()) return false;
	IF=_IF;
	return true;
}
				     
void ARBFN::Reset()
{
	if (tp)
	{
		delete[] tp;
		tp=0;
	}
	npat=varN=objN=0;
	is_built=is_trained=false;
	//
	if (neurons)
	{
		for (int i=0;i<n_neuron;i++) delete neurons[i];
		delete[] neurons;
		neurons=0;
	}
	//
	if (active)
	{
		delete[] active;
		active=0;
	}
	nactive=0;
	//
	if (w)
	{
		for (int k=0;k<n_neuron;k++) delete[] w[k];
		delete[] w; w=0;
	}
	//
	if (testpat)
	{
		delete[] testpat;
		testpat=0;
	}
	ntest=0;
	if (lms_corr_trace)
	{
		delete[] lms_corr_trace; lms_corr_trace=0;
	}
	if (ls_h)
	{
		for (int i=0;i<n_neuron;i++) delete[] ls_h[i];
		delete[] ls_h; ls_h=0;
	}
	//
	n_neuron=0;
	//
	if (best_neuron)
	{
		delete[] best_neuron;
		best_neuron=0;
	}
	bIF=false;
	gngout.str("");
}
	
		
bool ARBFN::Build(TrainingPattern* _tp, const int _npat, const int _varN,
		  const int _objN)
{
	IF.Resize(_varN,_objN,1.0/double(_varN));
	if (b_interp) return BuildInterpolate(_tp,_npat,_varN,_objN);
	return BuildApproximate(_tp,_npat,_varN,_objN);
}


bool ARBFN::CheckTrainingSet() const
{
	//Check for sizes
	for (int i=0;i<npat;i++)
		if (tp[i].var.size()!=varN ||
		    tp[i].obj.size()!=objN   )
			return false;
	//Check duplicate entries
	for (int i=0;i<npat-1;i++)
		for (int j=i+1;j<npat;j++)
		{
			double sum=0;
			for (int k=0;k<varN;k++)
			{
				const double e=tp[i].var[k]-tp[j].var[k];
				sum+=e*e;
			}
			if (sum<deps) return false;
		}
	return true;
}

bool ARBFN::Train()
{
	if (!is_built) return false;
	if (b_interp) return TrainInterpolate();
	return TrainApproximate();
}


bool ARBFN::Predict(const Vector<double>& v, Vector<double>& o)
{
	if (o.size()!=objN) return false;
	o=0.0;
	for (int i=0;i<nactive;i++)
	{
		const int id=active[i];
		for (int k=0;k<objN;k++)
		{
			const double h = 
				neurons[ id ]->ActivateIF(v,IF,k);
			o[k]+=h*w[id][k];
		}
	}
	return true;
}

bool ARBFN::PredictGrad(const Vector<double>& v, Matrix<double>& g)
{
	if (!is_built) return false;
	g.Resize(varN,objN,0.0);
	Vector<double> der(varN);
	for (int i=0;i<nactive;i++)
	{
		const int id=active[i];
		der=0.0;
		for (int k=0;k<objN;k++)
		{
			neurons[ id ]->ActivateDerIF(v,IF,k,der);
			for (int j=0;j<varN;j++)
				g[j][k]+=der[j]*w[id][k];
		}
	}
	return true;
}



void ARBFN::SetRadius(const double _r)
{
	for (int i=0;i<n_neuron;i++) neurons[i]->SetR(_r);
}

void ARBFN::Save(std::ostream& out) const
{
	
	bin_write(out,is_built);
	if (!is_built) return;
	bin_write(out,b_interp);
	if (b_interp) SaveInterp(out);
	else SaveApprx(out);
}

bool ARBFN::Restore(std::istream& in)
{
	Reset();
	bin_read(in,is_built);
	if (!is_built) return istream_operational(in);
	//
	bin_read(in,b_interp);
	if (b_interp) RestoreInterp(in);
	else RestoreApprx(in);
	return istream_operational(in);
}

double ARBFN::TestError2()
{
	double sum=0;
	Vector<double> pred(objN); 
	int inum=0;
	for (int i=0;i<ntest;i++)
	{
			pred=0;
			Predict(testpat[i].var,pred);
			for (int k=0;k<objN;k++)
			{
				const double e=testpat[i].obj[k]-pred[k];
				if (dabs(testpat[i].obj[k])>deps)
				{
					sum+=dabs(e/testpat[i].obj[k]);
					inum++;
				}
			}
	}
	return sum/double(inum);
}


void ARBFN::Initialize()
{
	SetRadius();
	ComputeWeights();
	errtest_min=TestError();	
	nbest=nactive;
	for (int i=0;i<nactive;i++)
	{
		const int id= active[i];
		best_neuron[i].SetEqualTo( *neurons[id] );
	}

}

double ARBFN::NeuronError(const int ineuron, const int itp)
{
	Vector<double> oo; oo.resize(objN);
	const bool b1=Predict(tp[itp].var,oo);
	double sum=0;
	for (int i=0;i<objN;i++)
	{
		if (!b1) oo[i]=.0;
		const double e=oo[i]-tp[itp].obj[i];
		sum+=e*e;
	}
	return sum;
}
	
void ARBFN::AdaptNormal(const int itp)
{
	const Vector<double>& v=tp[itp].var;
	const Vector<double>& o=tp[itp].obj;
	double *h=ls_h[0];
	for (int k=0;k<objN;k++)
	{
		double& corr_trace	= lms_corr_trace[k];
		double corr_trace_tp 	= 0;
		double out		= 0;
		for (int i=0;i<nactive;i++)
		{
			const int id		= active[i];
			const RBFNeuron* neuron = neurons[ id ];
			double& h_i		= h[id];
			h_i			= neuron->ActivateIF(v,IF,k);
			corr_trace_tp	       += h_i*h_i;
			out		       += h[id]*w[id][k];
		}
		corr_trace  *= 	double(lms_Ntrain);
		corr_trace   =	(corr_trace+corr_trace_tp)/
					double(++lms_Ntrain);
		//Set the learning rate
		const double lr = (2.0/corr_trace)*lms_lrfac;
		const double out_diff = o[k]-out;

		for (int i=0;i<nactive;i++)
		{
			const int id = active[i];
			w[id][k]+=lr*h[id]*out_diff;
		}
	}
}


void ARBFN::AdaptExtended()
{
	SetRadius();
	ComputeWeights();
}

bool ARBFN::IsImproving()
{
	const double err_test=TestError();
	if (err_test> ERR_BLOWUP) return false;
	if (err_test <errtest_min)
	{
		nbest=nactive;
		for (int i=0;i<nactive;i++)
		{
			const int id= active[i];
			best_neuron[i].SetEqualTo( *neurons[id] );
		}
		errtest_min = err_test;
		return true;
	}
	if (nactive > nbest + idle_n_growth) return false;
	else return true;
}

			
	
void ARBFN::Adjust(const int inew, const int in1, const int in2)
{
	return;
}


GNGNeuron* ARBFN::GetNeuron(const int i)
{
	return neurons[i];
}

int* ARBFN::GetActiveRedirector()
{
	return active;
}

int& ARBFN::GetActiveCounter()
{
	return nactive;
}


//
bool ARBFN::BuildInterpolate(TrainingPattern* tpin, const int _npat,
			     const int _varN, const int _objN)
{
	Reset();
	npat=_npat;
	varN=_varN;
	objN=_objN;
	min_cen=max_cen=npat;
	n_neuron=nactive=npat;
	//Allocate
	tp=new TrainingPattern[npat];
	neurons=new ARBFN_Neuron*[npat];
	active=new int[npat];
	w=new double*[npat];
	ls_h=new double*[npat];
	//
	ntest=0;
	
	const double rdefault=1.0;
	for (int i=0;i<npat;i++)
	{
		tp[i].var.assign ( tpin[i].var);
		tp[i].obj.assign ( tpin[i].obj);
		//
		neurons[i]=new ARBFN_Neuron(tp[i].var,rdefault);
		active[i]=i;
		w[i]=new double[objN];
		ls_h[i]=new double[npat];
	}
	const bool b1= CheckTrainingSet();
	if (!b1)
	{
		Reset();
		return false;
	}
	is_built=true;
	return true;
}
	
bool ARBFN::BuildApproximate(TrainingPattern* tpin, const int _npat,
			     const int _varN, const int _objN)
{
	Reset();
	varN=_varN;
	objN=_objN;
	ntest=int(test_to_tot_ratio*double(_npat));
	npat=_npat-ntest;
	if (min_cen>=npat) return false;
	max_cen=min1(max_cen,npat-1);
	if (max_cen>=npat || max_cen <min_cen) return false;
	////////////////////////////////////////////
	tp=new TrainingPattern[npat];
	testpat=new TrainingPattern[ntest];
	//
	//Randomly populate training and testing sets
	bool *bchose=new bool[_npat];
	for (int i=0;i<_npat;i++) bchose[i]=false;
	int ich=0;
	RandUnif myrand;
	while (ich<npat)
	{
		const int id=min1( myrand(0,_npat), _npat-1);
		if (bchose[id]==false)
		{
			bchose[id]=true; ich++;
		}
	}
	int itrain,itest;
	itrain=itest=0;
	for (int i=0;i<_npat;i++)
	{
		if (bchose[i]==true)
		{
			tp[itrain].var.assign(tpin[i].var);
			tp[itrain].obj.assign(tpin[i].obj);
			itrain++;
		}
		else
		{
			testpat[itest].var.assign(tpin[i].var);
			testpat[itest].obj.assign(tpin[i].obj);
			itest++;
		}
	}
	delete[] bchose; bchose=0;
	///////////////////////////////////////////////////////
	const bool b1=CheckTrainingSet();
	if (!b1)
	{
		Reset();
		return false;
	}
	///////////////////////////////////////////////////////
	n_neuron=max_cen;
	neurons=new ARBFN_Neuron*[n_neuron];
	active =new int[n_neuron];
	w=new double*[n_neuron];
	lms_corr_trace=new double[objN];
	ls_h=new double*[n_neuron];
	best_neuron=new ARBFN_Neuron[n_neuron];
	///////////////////////////////////////////////////////
	for (int i=0;i<n_neuron;i++)
	{
		neurons[i]=new ARBFN_Neuron();
		neurons[i]->SetDimension(varN);
		active[i]=i;
		w[i]=new double[objN];
		ls_h[i]=new double[n_neuron];
	}
	nactive=min_cen;
	errtest_min=dmax;
	for (int i=0;i<objN;i++)
		lms_corr_trace[i]=0.0;
	nbest=min_cen;
	////////////////////////////////////////////////////////
	is_built=true;
	return true;
}


bool ARBFN::TrainInterpolate()
{
	is_trained=ComputeWeights();
	return is_trained;
}

bool ARBFN::TrainApproximate()
{
	//Grow the network using GNG
	GNG gng(gngconfig,this);
	gng.Grow(max_cen,npat,npat,varN,min_cen,max_cen,tp);
	/////////////////////////////////////////////////
	//
	//Now choose the best neurons and update the redirection matrices
	nactive=nbest;
	for (int i=0;i<nactive;i++)
	{
		active[i]=i;
		neurons[i]->SetEqualTo(best_neuron[i]);
	}
	is_trained=ComputeWeights();
	return is_trained;
}


void ARBFN::SetRadius()
{
	if (b_interp) {std::cout <<"ERROR 1\n"; exit(0);}
	const double rdefault=.1*sqrt(double(varN));
	for (int i=0;i<nactive;i++)
	{
		const int id            = active[i];
		std::set<int>& connect  = neurons[id]->GetConnectivity();
		std::set<int>::const_iterator it;
		double r=0;
		int n_nb=0;
		const Vector<double>& c = neurons[id]->GetCenterConst();
		for (it=connect.begin();it!=connect.end();it++)
		{
			//Not redirection (connectivity uses global numb.);
			const int idn    = *it; 
			ARBFN_Neuron* nb = neurons[idn];
			r+=norm2var(c,nb->GetCenterConst());
			n_nb++;
		}
		if (n_nb>0)
		{
			r/=double(n_nb)*r_multfac;
			neurons[id]->SetR(r);
		}
		else
			neurons[id]->SetR(rdefault);
	}
}
		
bool ARBFN::ComputeWeights()
{
	if (bIF) return ComputeWeightsIF();
	return ComputeWeightsNOIF();
}

bool ARBFN::ComputeWeightsNOIF()
{
	Matrix<double> h_ext(npat,nactive);
	Matrix<double> w_ext(nactive,objN);
	Matrix<double> b_ext(npat,objN);
	double corr_trace1 = 0;
	int    ncortrace   = npat;
	for (int i=0;i<npat;i++)
	{
		for (int j=0;j<nactive;j++)
		{
			const int id = active[j];
			h_ext[i][j]=neurons[id]->Activate(tp[i].var);
			corr_trace1 += h_ext[i][j]*h_ext[i][j];
		}
		for (int j=0;j<objN;j++)
			b_ext[i][j]=tp[i].obj[j];
	}
	///////////////////////////////////////////////////////////////
	//////////////////////
	bool b1=false;
	const double wtol=deps*1000.0;
	const int svditer=100;
	if (b_interp)
	{
		//Tote npat=active neurons
		//b1=Gauss_Elim(h_ext,w_ext,b_ext);
		b1=SVDSolve(h_ext,w_ext,b_ext,svditer,wtol);
		ofstream out("weights");
		for (int i=0;i<nactive;i++)
		{
			for (int j=0;j<objN;j++)
			{
				w[i][j]=w_ext[i][j];
				out<<w[i][j]<<endl;
			}
		}
		out.close();
	}
	else
	{
		for (int k=0;k<objN;k++)
			lms_corr_trace[k]=corr_trace1/double(npat);
		lms_Ntrain=npat;
		/////////////////////////////////////////////////////
		//Form the least squares matrix
		Matrix<double> Ht=h_ext.Transpose();
		Matrix<double> Blsq=Ht*b_ext;
		Matrix<double> Hlsq=Ht*h_ext;
		//b1=Gauss_Elim(Hlsq,w_ext,Blsq);
		b1=SVDSolve(h_ext,w_ext,b_ext,svditer,wtol);
		for (int i=0;i<nactive;i++)
		{
			const int id= active[i];
			for (int k=0;k<objN;k++)
				w[id][k]=w_ext[i][k];
		}
	}
	return b1;
}

bool ARBFN::ComputeWeightsIF()
{
	int    ncortrace   = npat;
	bool btot=true;
	for (int k=0;k<objN;k++)
	{
		Matrix<double> h_ext(npat,nactive);
		Matrix<double> w_ext(nactive,1);
		Matrix<double> b_ext(npat,1);
		double corr_trace1 = 0;
		for (int i=0;i<npat;i++)
		{
			for (int j=0;j<nactive;j++)
			{
				const int id = active[j];
				h_ext[i][j]  =
					neurons[id]->ActivateIF(tp[i].var,
							        IF,k);
				corr_trace1 += h_ext[i][j]*h_ext[i][j];
			}
			b_ext[i][0]=tp[i].obj[k];
		}
		///////////////////////////////////////////////////////
		//////////////////////
		bool b1=false;
		const double wtol=deps*1000.0;
		const int svditer=100;
		if (b_interp)
		{
			//Tote npat=active neurons
	//		b1=Gauss_Elim(h_ext,w_ext,b_ext);
			b1=SVDSolve(h_ext,w_ext,b_ext,svditer,wtol);
			for (int i=0;i<nactive;i++)
					w[i][k]=w_ext[i][0];
		}
		else
		{
			lms_corr_trace[k]=corr_trace1/double(npat);
			lms_Ntrain=npat;
			/////////////////////////////////////////////////////
			//Form the least squares matrix
			Matrix<double> Ht=h_ext.Transpose();
			Matrix<double> Blsq=Ht*b_ext;
			Matrix<double> Hlsq=Ht*h_ext;
			//b1=Gauss_Elim(Hlsq,w_ext,Blsq);
			b1=SVDSolve(h_ext,w_ext,b_ext,svditer,wtol);
			for (int i=0;i<nactive;i++)
			{
				const int id= active[i];
					w[id][k]=w_ext[i][0];
			}
		}
		btot= (btot && b1);
	}
	return btot;
}
		
double ARBFN::TestError()
{
	double sum=0;
	Vector<double> pred(objN); 
	for (int i=0;i<ntest;i++)
	{
			pred=0;
			Predict(testpat[i].var,pred);
			for (int k=0;k<objN;k++)
			{
				const double e=testpat[i].obj[k]-pred[k];
				sum+=e*e;
			}
	}
	return sum;
}

			
double ARBFN::norm2var(const Vector<double>& a, const Vector<double>& b) const
{
	double sum=0;
	for (int i=0;i<varN;i++)
	{
		const double e= a[i]-b[i];
		sum+=e*e;
	}
	return sqrt(sum);
}

void ARBFN::SaveTP(std::ostream& out) const
{
	for (int i=0;i<npat;i++)
	{
		for (int j=0;j<varN;j++) bin_write(out,tp[i].var[j]);
		for (int j=0;j<objN;j++) bin_write(out,tp[i].obj[j]);
	}
}

void ARBFN::SaveTestP(std::ostream& out) const
{
	for (int i=0;i<ntest;i++)
	{
		for (int j=0;j<varN;j++) bin_write(out,testpat[i].var[j]);
		for (int j=0;j<objN;j++) bin_write(out,testpat[i].obj[j]);
	}
}

void ARBFN::SaveIFs(std::ostream& out) const
{
	bin_write(out,bIF);
	const int d1=IF.dim1();
	const int d2=IF.dim2();
	bin_write(out,d1);
	bin_write(out,d2);
	for (int i=0;i<d1;i++)
		for (int j=0;j<d2;j++)
			bin_write(out,IF[i][j]);
	return;
}

void ARBFN::SaveW(std::ostream& out) const
{
	bin_write(out,is_trained);
	for (int i=0;i<n_neuron;i++)
		for (int j=0;j<objN;j++)
			bin_write(out,w[i][j]);
}


void ARBFN::SaveInterp(std::ostream& out) const
{
	bin_write(out,npat);
	bin_write(out,varN);
	bin_write(out,objN);
	SaveTP(out);
	SaveIFs(out);
	SaveW(out);
}
		
	
void ARBFN::SaveApprx(std::ostream& out) const
{
	bin_write(out,npat);
	bin_write(out,varN);
	bin_write(out,objN);
	//Config
	bin_write(out,ntest);
	bin_write(out,n_neuron);
	bin_write(out,test_to_tot_ratio);
	bin_write(out,min_cen);
	bin_write(out,max_cen);
	bin_write(out,r_multfac);
	bin_write(out,lms_lrfac);
	bin_write(out,idle_n_growth);
	//GNGConfig
	//
	bin_write(out,gngconfig->mvc);
	bin_write(out,gngconfig->mvn);
	bin_write(out,gngconfig->max_connect_age);
	bin_write(out,gngconfig->addstep_dimmult);
	bin_write(out,gngconfig->errdecr_interp);
	bin_write(out,gngconfig->errdecr_gen);
	bin_write(out,gngconfig->tp_close_frac);
	bin_write(out,gngconfig->tp_close_prob);
	bin_write(out,gngconfig->tp_close_iter);
	//Data
	SaveTP(out);
	SaveIFs(out);
	SaveW(out);
	SaveTestP(out);
	//
	for (int i=0;i<n_neuron;i++)
	{
		neurons[i]->Save(out);
		best_neuron[i].Save(out);
		bin_write(out,active[i]);
		for (int j=0;j<n_neuron;j++)
			bin_write(out,ls_h[i][j]);
	}
	for (int i=0;i<objN;i++) bin_write(out,lms_corr_trace[i]);
	bin_write(out,nactive);
	bin_write(out,errtest_min);
	bin_write(out,nbest);
	bin_write(out,lms_Ntrain);
}
	
	

bool ARBFN::RestoreTP(std::istream& in,TrainingPattern* _tp,
		      const int _npat,
	              const int _varN, const int _objN)
{
	for (int i=0;i<_npat;i++)
	{
		_tp[i].var.resize(_varN);
		_tp[i].obj.resize(_objN);
		for (int j=0;j<_varN;j++) bin_read(in,_tp[i].var[j]);
		for (int j=0;j<_objN;j++) bin_read(in,_tp[i].obj[j]);
	}
	return istream_operational(in);
}

bool ARBFN::RestoreTestP(std::istream& in) 
{
	for (int i=0;i<ntest;i++)
	{
		testpat[i].var.resize(varN);
		testpat[i].obj.resize(objN);
		for (int j=0;j<varN;j++) bin_read(in,testpat[i].var[j]);
		for (int j=0;j<objN;j++) bin_read(in,testpat[i].obj[j]);
	}
	return istream_operational(in);
}

bool ARBFN::RestoreIFs(std::istream& in) 
{
	bin_read(in,bIF);
	int d1,d2;
	d1=d2=0;
	bin_read(in,d1);
	bin_read(in,d2);
	if (d1!=varN || d2!=objN) return false;
	for (int i=0;i<d1;i++)
		for (int j=0;j<d2;j++)
			bin_read(in,IF[i][j]);
	return istream_operational(in);
}

bool ARBFN::RestoreW(std::istream& in) 
{
	bin_read(in,is_trained);
	for (int i=0;i<n_neuron;i++)
		for (int j=0;j<objN;j++)
			bin_read(in,w[i][j]);
	return istream_operational(in);
}

bool ARBFN::RestoreInterp(std::istream& in) 
{
	int _npat,_varN,_objN;
	bin_read(in,_npat);
	bin_read(in,_varN);
	bin_read(in,_objN);
	TrainingPattern* _tp=new TrainingPattern[_npat];
	RestoreTP(in,_tp,_npat,_varN,_objN);
	//
	const bool b1=BuildInterpolate(_tp,_npat,_varN,_objN);
	delete[] _tp; _tp=0;
	//
	if (!b1) return false;
	IF.Resize(_varN,_objN,1.0/double(_varN));
	if (!RestoreIFs(in)) return false;
	RestoreW(in);
	return istream_operational(in);
}
		

bool ARBFN::RestoreApprx(std::istream& in) 
{
	Reset();
	int _npat,_varN,_objN;
	bin_read(in,_npat);
	bin_read(in,_varN);
	bin_read(in,_objN);
	//Config
	bin_read(in,ntest);
	varN=_varN;
	objN=_objN;
	npat=_npat;
	//
	bin_read(in,n_neuron);
	bin_read(in,test_to_tot_ratio);
  	bin_read(in,min_cen);
  	bin_read(in,max_cen);
  	bin_read(in,r_multfac);
  	bin_read(in,lms_lrfac);
  	bin_read(in,idle_n_growth);
	//GNGConfig
	//
	bin_read(in,gngconfig->mvc);
	bin_read(in,gngconfig->mvn);
	bin_read(in,gngconfig->max_connect_age);
	bin_read(in,gngconfig->addstep_dimmult);
	bin_read(in,gngconfig->errdecr_interp);
	bin_read(in,gngconfig->errdecr_gen);
	bin_read(in,gngconfig->tp_close_frac);
	bin_read(in,gngconfig->tp_close_prob);
	bin_read(in,gngconfig->tp_close_iter);
	//Data
	tp=new TrainingPattern[npat];
	testpat=new TrainingPattern[ntest];
	//
	//
	max_cen=min1(max_cen,npat-1);
	n_neuron=max_cen;
	neurons=new ARBFN_Neuron*[n_neuron];
	active =new int[n_neuron];
	w=new double*[n_neuron];
	lms_corr_trace=new double[objN];
	ls_h=new double*[n_neuron];
	best_neuron=new ARBFN_Neuron[n_neuron];
	///////////////////////////////////////////////////////
	for (int i=0;i<n_neuron;i++)
	{
		neurons[i]=new ARBFN_Neuron();
		neurons[i]->SetDimension(varN);
		active[i]=i;
		w[i]=new double[objN];
		ls_h[i]=new double[n_neuron];
	}
	nactive=min_cen;
	errtest_min=dmax;
	for (int i=0;i<objN;i++)
		lms_corr_trace[i]=0.0;
	nbest=min_cen;
	////////////////////////////////////////////////////////
	is_built=true;
	//
	RestoreTP(in,tp,npat,varN,objN);
	IF.Resize(_varN,_objN,1.0/double(_varN));
	RestoreIFs(in);
	RestoreW(in);
	RestoreTestP(in);
	//
	for (int i=0;i<n_neuron;i++)
	{
		neurons[i]->Restore(in);
		best_neuron[i].Restore(in);
		bin_read(in,active[i]);
		for (int j=0;j<n_neuron;j++)
			bin_read(in,ls_h[i][j]);
	}
	for (int i=0;i<objN;i++) bin_read(in,lms_corr_trace[i]);
	bin_read(in,nactive);
	bin_read(in,errtest_min);
	bin_read(in,nbest);
	bin_read(in,lms_Ntrain);
	return istream_operational(in);
}



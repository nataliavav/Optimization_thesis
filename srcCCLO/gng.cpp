#include "gng.h"
#include "filesystem.h"
#include <set>
#include <list>
#include <cmath>

using namespace METAMODELS;

estring GNG_exception::label= "GNG exception:\n";


GNGConfig::GNGConfig(std::ostream& _out) : GConfig(_out)
{
	mvc=mvn=addstep_dimmult=errdecr_interp=errdecr_gen=0;
	tp_close_frac=tp_close_prob=tp_close_iter=0;
	max_connect_age=0;
}

GNGConfig::~GNGConfig()
{
}

bool GNGConfig::ReadConfig(std::istream& in)
{
	in >> mvc;		istream_nl(in);
	in >> mvn;		istream_nl(in);
	in >> max_connect_age;	istream_nl(in);
	in >> addstep_dimmult;	istream_nl(in);
	in >> errdecr_interp;	istream_nl(in);
	in >> errdecr_gen;	istream_nl(in);
	in >> tp_close_frac;	istream_nl(in);
	in >> tp_close_prob;	istream_nl(in);
	in >> tp_close_iter;	istream_nl(in);

	const bool b1=istream_operational(in);
	if (!b1) return false;
	const bool b2=CheckConfig();

	return b1&&b2;
}

bool GNGConfig::CheckConfig()
{
	if (mvc<=0) throw GNG_exception
		("Center movement factor should be > 0 \n");
	if (mvn<=0) throw GNG_exception
		("Center neighbours movement factor should be > 0 \n");
	if (max_connect_age<=0) throw GNG_exception
		("Maximum connection age should be >0\n");
	if (addstep_dimmult <=0 ) throw GNG_exception
		("Center insertion coef. for dimensionaliry should be >0\n");
	if (errdecr_interp <=0 || errdecr_interp>=1) throw GNG_exception
		("Error decrease factor should lie in (0,1)\n");
	if (errdecr_gen <=0 || errdecr_gen>=1) throw GNG_exception
		("Error decrease factor should lie in (0,1)\n");
	if (tp_close_frac <=0 || tp_close_frac>=1) throw GNG_exception
	    ("Patterns % to be regarded as neighbours should lie in (0,1)\n");
	if (tp_close_prob <0 || tp_close_prob>1) throw GNG_exception
		("Neighbor selection probability should lie in [0,1]\n");
	if (tp_close_iter <0 || tp_close_iter>1) throw GNG_exception
  	 ("Iterations % to select close to new node tp should lie in [0,1]\n");

	return true;
	
} 


GNG::GNG(const GNGConfig* const _config, GNGCluster* _cluster):
		config(_config), cluster(_cluster)
{
}

void GNG::Grow(const int Nneuron_max, const int tp_Ntrain, 
	       const int tp_Ntrain_max, const int Nin,
	       const int Nneuron_start, const int Nneuron_stop,
	       TrainingPattern* tp_train
	       )
{
	///////// Initialization
	const double mvc=config->mvc;
	const double mvn=config->mvn;
	const int max_connect_age=config->max_connect_age;
	const double addstep_dimmult=config->addstep_dimmult;
	const double errdecr_interp=config->errdecr_interp;
	const double errdecr_gen=config->errdecr_gen;
	const double tp_close_frac=config->tp_close_frac;
	const double tp_close_prob=config->tp_close_prob;
	const int    addstep= int(addstep_dimmult*double(Nin)); 
	const double tp_close_niter=max(0, int (config->tp_close_iter*
					        addstep_dimmult*double(Nin)));
	const int tp_Nclose = int(tp_close_frac*double(tp_Ntrain));
	//
	//Checks
	//
	if (Nneuron_max < 2)
		throw GNG_exception
		("Invalid maximum number of neurons\n");
	if (tp_Ntrain_max < 1)
		throw GNG_exception
		("Invalid maximum number of training patterns\n");
	/////
	//Memory allocation
	Vector<int> neurons_elim(Nneuron_max);
	int** connect_age=new int*[Nneuron_max];
	for (int i=0;i<Nneuron_max;i++) connect_age[i]=new int[Nneuron_max];
	double* err= new double[Nneuron_max];
	const int tp_Nclose_max=int(tp_Ntrain_max*double(tp_close_frac));
	int* tp_close=new int[tp_Nclose_max];
	tp_dist_t* tp_dist= new tp_dist_t[tp_Ntrain_max];
	double* tp_train_min=new double[Nin];
	double* tp_train_max=new double[Nin];
	std::set<int> neurons_active;
	std::list<int> neurons_avail;
	///
	//	Redirector and counter from cluster class
	int& Nneuron=cluster->GetActiveCounter();
	int *neurons_active_indx = cluster->GetActiveRedirector();
	//Main algorithm
	//
	//- Calculate min-max
	for (int i=0;i<Nin;i++)
		tp_train_min[i]=tp_train_max[i]=tp_train[0].var[i];
	for (int i=1;i<tp_Ntrain_max;i++)
		for (int k=0;k<Nin;k++)
		{
			tp_train_min[k]=min1(tp_train_min[k],
					     tp_train[i].var[k]);
			tp_train_max[k]=max1(tp_train_max[k],
					     tp_train[i].var[k]);
		}
	//- Get the neurons
	GNGNeuron** neurons=new GNGNeuron*[Nneuron_stop];
	for (int i=0;i<Nneuron_stop;i++) neurons[i]= cluster->GetNeuron(i);
	//- Set initial neuron data
	for (int i=0;i<Nneuron_start; i++)
	{
		GNGNeuron* neuron= neurons[i];
		Vector<double>& c = neuron->GetCenter();
		std::set<int>& connect=neuron->GetConnectivity();
		//- coords
		for (int k=0;k<Nin;k++)
			c[k]=myrand()*(tp_train_max[k]-tp_train_min[k])
				+tp_train_min[k];
		//-connectivity
		connect.clear();
		for (int j=0;j<Nneuron_start;j++)
		{
			if (i!=j) connect.insert(j);
			connect_age[i][j]=0;
		}

	}
	//Classify neurons to active & available, reset error to active,
	//and finally mark their index
	neurons_active.clear();
	neurons_avail.clear();
	Nneuron=0;
	for (int i=0;i<Nneuron_stop;i++) 
	{
		if (i<Nneuron_start)
		{
			neurons_active.insert(i);
			err[i]=0.0;
			neurons_active_indx[i]=i;
			Nneuron++;
		}
		else
		{
			neurons_avail.push_back(i);
			neurons[i]->GetConnectivity().clear();
		}
	}
	//Initialize the clustering class
	cluster->Initialize();
	//Set some counters
	int niter=0;
	int niter_add=0;
	RandUnif& rand=myrand;
	bool added=false;
	bool fracok=(tp_Nclose>0);
	while (neurons_active.size() < Nneuron_stop)
	{
		niter++;
		//Pick a training pattern
		int itp;
		if ((niter-niter_add <tp_close_niter&&niter) 
				&& added && fracok)
		{
			if (rand() < tp_close_prob)
				itp=tp_close[rand(0,tp_Nclose)];
			else
				itp=rand(0,tp_Ntrain);
		}
		else
			itp = rand(0, tp_Ntrain);
		const TrainingPattern& tp= tp_train[itp];
		//
		// Find the two closest neurons
		set<int>::const_iterator it=neurons_active.begin();
		int icl1,icl2;
		double dist_min1,dist_min2;
		icl1=*(it++);
		icl2=*(it++);
		dist_min1=norm2var(Nin,tp.var,neurons[icl1]->GetCenterConst());
		dist_min2=norm2var(Nin,tp.var,neurons[icl2]->GetCenterConst());
		if (dist_min2<dist_min1)
		{
			const double dt1=dist_min1; 
				dist_min1=dist_min2; dist_min2=dt1;
			const int it1=icl1;
				icl1=icl2;	icl2=it1;
		}
		for (;it!=neurons_active.end();it++)
		{
			const double dist=
				norm2var(Nin,tp.var,
					 neurons[*it]->GetCenterConst());
			if (dist<dist_min1)
			{
				dist_min2 = dist_min1;
				icl2=icl1;
				dist_min1=dist;
				icl1=*it;
			}
			else if (dist<dist_min2)
			{
				dist_min2 =  dist;
				icl2=*it;
			}
		}
		//Insert a connection between the closest neurons
		GNGNeuron* neuron_cl=neurons[icl1];
		set<int>& connect_cl=neuron_cl->GetConnectivity();
		connect_cl.insert(icl2);
		neurons[icl2]->GetConnectivity().insert(icl1);
		//Refresh the connection
		int* connect_age_cl 	= connect_age[icl1];
		connect_age_cl[icl2] 	= 0;
		connect_age[icl2][icl1] = 0;
		//
		Vector<double>& c=neuron_cl->GetCenter();
		//Move the best-matching neuron (icl1) towards the pattern
		for (int i=0;i<Nin;i++)
			c[i]+=mvc*(tp.var[i]-c[i]);
		// Adjust the neighborhood:
		int nelim=0;
		for (it= connect_cl.begin();it!=connect_cl.end();it++)
		{
			const int inb=*it;
			if (inb<0 || inb >=Nneuron_stop)
				throw GNG_exception("Invalid inb\n");
			GNGNeuron* neuron_nb = neurons[inb];
			Vector<double>& c =neuron_nb->GetCenter();
			for (int i=0;i<Nin;i++)
				c[i]+=mvn*(tp.var[i]-c[i]);
			//Increase the age of the connection
			const int age = ++connect_age_cl[inb];
			++connect_age[inb][icl1];
			if (age> max_connect_age)
				neurons_elim[nelim++] = inb;
		}
		for (int i=0;i<nelim;i++)
		{
			const int inb= neurons_elim[i];
			GNGNeuron* neuron_nb=neurons[inb];
			connect_cl.erase(inb);
			neuron_nb->GetConnectivity().erase(icl1);
			/////neurons[icl1] has connection. Just check it's
			//neighbour
			if (neuron_nb->GetConnectivity().empty())
			{
				neurons_active.erase(inb);
				neurons_avail.push_back(inb);
				if (neurons_active.size() <Nneuron_start)
					throw GNG_exception
					 ("Neurons_active < Neurons_start\n");
			}
		}
		// Set the indices of the active neurons
		Nneuron=0;
		for (it = neurons_active.begin();it!=neurons_active.end();
				it++,Nneuron++)
			neurons_active_indx[Nneuron] = *it;
		//Increase the error associated with the best-matching neuron
		err[icl1] += cluster->NeuronError( icl1,itp);
		cluster->AdaptNormal(itp);
		//Add a new neuron?
		if (niter%addstep == 0 )
		{
			niter_add = niter;
			if (!cluster->IsImproving())
				break;
			//Find the two neurons to interpolate a new one from:
			//1st Neuron: the one with the max associated error
			double err_max;
			it= neurons_active.begin();
			err_max= err[*it];
			int ineuron1 = *it; it++;
			for (; it!= neurons_active.end(); it ++)
				if (err[*it] >err_max)
				{
					err_max =err[*it];
					ineuron1= *it;
				}
			//2nd Neuron: Neighbour with largest associated error
			const set<int>& connect = neurons[ineuron1]->
							GetConnectivityConst();
			it =connect.begin();
			err_max= err[*it];
			int ineuron2 = *it; it++;
			for (; it!=connect.end();it++)
				if (err[*it]> err_max)
				{
					err_max =err[*it];
					ineuron2 = *it;
				}
			//Retrieve a neuron from the pool
			const int ineuron_new= neurons_avail.back();
			neurons_avail.pop_back();
			neurons_active.insert(ineuron_new);
			GNGNeuron* neuron_new= neurons[ineuron_new];
			//
			GNGNeuron* neuron1 = neurons[ineuron1];
			GNGNeuron* neuron2 = neurons[ineuron2];
			Vector<double>& c_new	 = neuron_new->GetCenter();
			const Vector<double>& c1 = neuron1->GetCenter();
			const Vector<double>& c2 = neuron2->GetCenter();
			for (int i=0;i<Nin;i++)
				c_new[i]=.5*(c1[i]+c2[i]);
			//
			cluster->Adjust(ineuron_new,ineuron1, ineuron2);
			//
			//Adjust connections
			neuron1->GetConnectivity().erase(ineuron2);
			neuron1->GetConnectivity().insert(ineuron_new);
			neuron2->GetConnectivity().erase(ineuron1);
			neuron2->GetConnectivity().insert(ineuron_new);
			//
			neuron_new->GetConnectivity().insert(ineuron1);
			neuron_new->GetConnectivity().insert(ineuron2);
			connect_age[ineuron_new][ineuron1]	= 0;
			connect_age[ineuron1][ineuron_new]	= 0;
			connect_age[ineuron_new][ineuron2]	= 0;
			connect_age[ineuron2][ineuron_new]	= 0;
			//Decrease the error of the initial neurons
			err[ineuron1] *= errdecr_interp;
			err[ineuron2] *= errdecr_interp;
			//Interpolate the error of the new neuron
			err[ineuron_new] = .5*(err[ineuron1]+err[ineuron2]);
			// Locate the training patterns closest to the
			// 	newly added neuron
			if (tp_Nclose)
			{
				for (int i=0;i<tp_Ntrain;i++)
				{
					tp_dist[i].d= norm2var(Nin,
							c_new,tp_train[i].var);
					tp_dist[i].i = i;
				}
				partial_sort(tp_dist,tp_dist+tp_Nclose,
						     tp_dist+tp_Ntrain);
				for (int i=0;i<tp_Nclose;i++)
					tp_close[i]=tp_dist[i].i;
			}
			added=true;
			neurons_active_indx[Nneuron++] = ineuron_new;
			cluster->AdaptExtended();
		}
		//
		//Decrease the error variables of all the neurons
		for (it =neurons_active.begin();it!=neurons_active.end();it++)
			err[*it]*=errdecr_gen;
	}
	//
	cluster->IsImproving();
	//////
	//Memory de-allocation
	for (int i=0;i<Nneuron_max;i++) delete[] connect_age[i];
	delete[] connect_age;
	delete[] err;
	delete[] tp_close;
	delete[] tp_dist;
	delete[] tp_train_min;
	delete[] tp_train_max;
	//////////////////////
		
}


const double GNG::norm2var(const int Nin, const Vector<double>& v1,
					  const Vector<double>& v2) const
{
	double sum=0;
	for (int i=0;i<Nin;i++)
	{
		const double e=v1[i]-v2[i];
		sum+=e*e;
	}
	return sqrt(sum);
}

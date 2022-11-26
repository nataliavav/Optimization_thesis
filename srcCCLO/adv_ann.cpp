
#include "adv_ann.h"
#include "metamodel.h"
#include "eigen.h"
#include "filesystem.h"

using namespace METAMODELS;

///////////////////////////////ACTIVATION

const double GaussianActivation::activ(const double u, const double r) const
{
	return exp(-u*u/r/r);
}

const double GaussianActivation::dactiv(const double u, const double r) const
{
	return -2.0*u*r/r*activ(u,r);
}

const double GaussianActivation::ddactiv(const double u,const double r) const
{
	return -2.0/r/r * activ(u,r) -2.0*u/r/r*dactiv(u,r);
}

///////////
const double MultiQuadActivation::activ(const double u, const double r) const
{
	return sqrt(u*u+r*r);
}

const double MultiQuadActivation::dactiv(const double u, const double r) const
{
	return u/activ(u,r);
}

const double MultiQuadActivation::ddactiv(const double u,const double r) const
{
	const double aa=activ(u,r);
	return 1.0/aa-u*dactiv(u,r)/aa/aa;
}

///////////
const double InvMultiQuadActivation::activ(const double u, const double r) const
{
	return 1.0/sqrt(u*u+r*r);
}

const double InvMultiQuadActivation::dactiv(const double u,const double r)const
{
	return -u/(u*u+r*r)*activ(u,r);
}

const double InvMultiQuadActivation::ddactiv(const double u,const double r)const
{
	const double aa=activ(u,r);
	const double da=dactiv(u,r);
	return -(u*da+aa*(r*r-u*u)/(u*u+r*r)) /(u*u+r*r);
}





///////////////////////////////Neuron
Neuron::Neuron()
{
	c.resize(0);
}

Neuron::Neuron(const Neuron& _c)
{
	SetEqualTo(_c);
}

Neuron::Neuron(const Vector<double>& _c)
{
	c.assign(_c);
}

Neuron::~Neuron()
{
}

void Neuron::SetEqualTo(const Neuron& _c)
{
	if (this!=&_c) c.assign(_c.c);
}

const int Neuron::GetDimension() const
{
	return c.size();
}

void Neuron::SetDimension(const int idim)
{
	c.resize(idim);
}

Vector<double>& Neuron::GetCenter()
{
	return c;
}

const Vector<double>& Neuron::GetCenterConst() const
{
	return c;
}

///////////////////////////////////  RBF Neuron

RBFNeuron::RBFNeuron() : Neuron()
{
//	a=new MultiQuadActivation();
	a=new GaussianActivation();
//	a=new InvMultiQuadActivation();	// vera 03/06/2014 for titan
	r=1.0;
}

RBFNeuron::RBFNeuron(const RBFNeuron& _c) : Neuron(_c)
{
//	a=new MultiQuadActivation();
	a=new GaussianActivation(); //if add other activs, just check
				    //with casting 
//	a=new InvMultiQuadActivation();	// vera 03/06/2014 for titan
	r=_c.r;
}

RBFNeuron::RBFNeuron(const Vector<double>& _c) : Neuron(_c)
{
//	a=new MultiQuadActivation();
	a=new GaussianActivation();
//	a=new InvMultiQuadActivation();	// vera 03/06/2014 for titan
	r=1.0;
}

RBFNeuron::RBFNeuron(const Vector<double>& _c, const double _r) :
			Neuron(_c)
{
//	a=new MultiQuadActivation();
	a=new GaussianActivation();
//	a=new InvMultiQuadActivation();	// vera 03/06/2014 for titan
	r=_r;
}

RBFNeuron::~RBFNeuron()
{
	delete a; a=0;
}

const double RBFNeuron::Activate(const Vector<double>& x) const
{
	if (x.size()!=c.size()) throw Mmodel_exception
		("RBFNeuron:Invalid pattern input size!\n");
		
	double sum=0;
	for (int i=0;i<c.size();i++)
		sum+=(x[i]-c[i])*(x[i]-c[i]);
	sum=sqrt(sum);

	return a->activ(sum,r);
}

const double RBFNeuron::ActivateDer(const Vector<double>& x,
				          Vector<double>& der) const
{
	if (x.size()!=c.size()) throw Mmodel_exception
		("RBFNeuron:Invalid pattern input size!\n");
	if (x.size()!=der.size()) throw Mmodel_exception
		("RBFNeuron:Invalid derivative vector input size!\n");
	double sum=0;
	for (int i=0;i<c.size();i++)
	{
		sum+=(x[i]-c[i])*(x[i]-c[i]);
		der[i]=2.0*(x[i]-c[i]);
	}
	sum=sqrt(sum);
	//
	double invsum=.5/sum;
	if (sum<deps) invsum=0;
	//
	const double dd=a->dactiv(sum,r);
	for (int i=0;i<c.size();i++)
		der[i]*=dd*invsum;
	return a->activ(sum,r);
}

const double RBFNeuron::ActivateIF(const Vector<double>& x,
				   const Matrix<double>& IF, 
				   const int io) const
{
	if (x.size()!=c.size()) throw Mmodel_exception
		("RBFNeuron:Invalid pattern input size!\n");
		
	double sum=0;
	for (int i=0;i<c.size();i++)
		sum+=IF[i][io]*(x[i]-c[i])*(x[i]-c[i]);
	sum*=double(c.size());
	sum=sqrt(sum);

	return a->activ(sum,r);
}

const double RBFNeuron::ActivateDerIF(const Vector<double>& x,
				      const Matrix<double>& IF,
				      const int io,
				      Vector<double>& der) const
{
	/////
	if (x.size()!=c.size()) throw Mmodel_exception
		("RBFNeuron:Invalid pattern input size!\n");
	if (x.size()!=der.size()) throw Mmodel_exception
		("RBFNeuron:Invalid derivative vector input size!\n");
	double sum=0;
	for (int i=0;i<c.size();i++)
	{
		sum+=IF[i][io]*(x[i]-c[i])*(x[i]-c[i]);
		der[i]=2.0*(x[i]-c[i])*IF[i][io]*double(c.size());
	}
	sum*=double(c.size());
	sum=sqrt(sum);
	//
	double invsum=.5/sum;
	if (sum<deps) invsum=0;
	//
	const double dd=a->dactiv(sum,r);
	for (int i=0;i<c.size();i++)
		der[i]*=dd*invsum;
	return a->activ(sum,r);
}



const double RBFNeuron::GetR() const
{
	return r;
}

void RBFNeuron::SetR(const double _r) 
{
	r=_r;
}


GNGNeuron::GNGNeuron() : Neuron()
{
	connect.clear();
}

GNGNeuron::GNGNeuron(const GNGNeuron& _c) : Neuron(_c)
{
	connect.clear();
	std::set<int>::const_iterator i;
	for (i=_c.connect.begin();i!=_c.connect.end();i++)
		connect.insert(*i);
}

GNGNeuron::GNGNeuron(const Vector<double>& _c) : Neuron(_c)
{
	connect.clear();
}

GNGNeuron::~GNGNeuron()
{
}

std::set<int>& GNGNeuron::GetConnectivity()
{
	return connect;
}
	

const std::set<int>& GNGNeuron::GetConnectivityConst() const
{
	return connect;
}
	

AdvancedANN::~AdvancedANN()
{
	//Do nothing
}

bool AdvancedANN::CheckTrainingSet() const
{
	for (int i=0;i<npat;i++)
		if (tp[i].var.size()!=varN ||
		    tp[i].obj.size()!=objN   )
			return false;
return true;
}



///////////////////////////////  MST

CMST::CMST(TrainingPattern* _tp, const int _npat, const int _varN, 
		const int _objN) : npat(_npat), varN(_varN), objN(_objN)
{
	//copy the training patterns to the local storage
	tp=new TrainingPattern[npat];
	dist=new double*[npat];
	for (int i=0;i<npat;i++)
	{
		dist[i]=new double[npat];
		tp[i].var.assign(  _tp[i].var  );
		tp[i].obj.assign(  _tp[i].obj  );
		tp[i].bfailed = tp[i].bfailed;
	}
	//
	CalcDistances();
	ConstructMST();
}

CMST::~CMST()
{
	delete[] tp; tp=0;
	for (int i=0;i<npat;i++)
		delete[] dist[i];
	delete[] dist; dist=0;
}

const double CMST::Norm2Var(const int i, const int j) const
{
	double sum=0;
	for (int k=0;k<varN;k++)
	{
		const double e= tp[i].var[k]-tp[j].var[k] ;
		sum+=e*e;
	}
	return sqrt(sum);
}

const double CMST::Norm2Obj(const int i, const int j) const
{
	double sum=0;
	for (int k=0;k<objN;k++)
	{
		const double e= tp[i].obj[k]-tp[j].obj[k] ;
		sum+=e*e;
	}
	return sqrt(sum);
}


void CMST::CalcDistances()
{
	for (int i=0;i<npat;i++)
	{
		if (tp[i].var.size()!=varN) throw Mmodel_exception
			("CMST:Invalid variables vector size!\n");
		if (tp[i].obj.size()!=objN) throw Mmodel_exception
			("CMST:Invalid objectives vector size!\n");
		dist[i][i]=-1;  //Do not allow self-linking in MST
		for (int j=i+1;j<npat;j++)
			dist[i][j]=dist[j][i]=Norm2Var(i,j);
	}
}
	
void CMST::ConstructMST()
{
	tree.resize(npat);
	const int n=npat;
	std::vector<bool> visited (npat,false);
	int src  = 0;	//current/active source node
	int dest = 0;	//current/active destination node
	std::priority_queue<edge_t> pq;
	pq.push(edge_t(src,dest,0.0));
	int Nunvisited=n;
	d_ave=0.0;
	double d2_sum=0.0;
	int nbranch=0;
	int iii=0;
	while (Nunvisited)
	{
		const edge_t& top_edge=pq.top();
		src  = top_edge.src;
		dest = top_edge.dest;
		pq.pop();
		if (!visited[dest])
		{
			if (dest!=src)
			{
				tree[src].push_back(dest);
				tree[dest].push_back(src);
				const double d = dist[src][dest];
				d_ave   += d;
				d2_sum  += d*d;
				++nbranch;
			}
			visited[dest] = true;
			--Nunvisited;
			for (int i=0;i<n;i++)
				if (dist[dest][i]> 0.0 && !visited[i] )
					pq.push(edge_t(dest,i,dist[dest][i]));
		}
	}
	const double inv_nbranch=1.0/double(nbranch);
	d_ave *= inv_nbranch;
	d_dev  = sqrt((d2_sum - double(nbranch)*d_ave*d_ave)*inv_nbranch);
}

void CMST::WritePath(std::ostream& out, const int id1, const int id2) const
{
	for (int i=0;i<objN;i++)
		out <<tp[id1].obj[i]<<"  ";
	out<<"\n";
	for (int i=0;i<objN;i++)
		out <<tp[id2].obj[i]<<"  ";
	out<<"\n";
	out<<"\n\n";
}

const double CMST::GetAveDist() const
{
	return d_ave;
}

const double CMST::GetStDev() const
{
	return d_dev;
}
	

void CMST::ConditionedWalk(const int start_node, const double minkeep_per,
			   const double TOL, int& Nkeep,
			   std::vector<int>& keep,
			   double& t_ave, double& t_dev) const
{
	typedef std::list<int> list_t;
	const int minNkeep=int(minkeep_per*double(npat));
	const list_t& t=tree[start_node];
	std::vector<bool> bvisited(npat,false);
	Nkeep=0;
	bvisited[start_node] = true ;
	//
	keep.resize(npat);
	for (list_t::const_iterator it= t.begin(); it!=t.end();it++)
		keep[Nkeep++] = *it;
	//
	bool bCont;
	int start=0;
	int end=Nkeep;
	int Nedges=0;
	t_ave=0.0;
	double d2_sum=0.0;
	do
	{
		bCont=false;
		bool bNoChk = Nkeep<minNkeep ? true :false;
		for (int i=start;i<end;i++)
		{
			const int curr_id = keep[i];
			bvisited [curr_id] = true;
			const list_t& neighb = tree[curr_id];
			const double* curr_dist = dist[curr_id];
			for (list_t::const_iterator it=neighb.begin();
					it!=neighb.end();it++)
			{
				const int id= *it;
				const double& d=curr_dist[id];
				if (!bvisited[id])
				{
					if (bNoChk || d<t_ave+TOL*t_dev)
					{
						keep[Nkeep++] = id;
						bCont=true;
						t_ave*=double(Nedges);
						t_ave+=d;
						const double invedges=1.0/
							double(++Nedges);
						t_ave*=invedges;
						d2_sum+=d*d;
						t_dev=sqrt((d2_sum-Nedges*
								d_ave*d_ave)
								*invedges);
					}
				}
			}
		}
		start=end;
		end=Nkeep;
	}while (bCont);
}


void CMST::WriteTree(const int start_node,std::ostream& out) const
{

	typedef std::list<int> list_t;
	const list_t& t=tree[start_node];
	std::vector<bool> bvisited(npat,false);
	int Nkeep=0;
	bvisited[start_node] = true ;
	//
	std::vector<int> keep;
	keep.resize(npat);
	for (list_t::const_iterator it= t.begin(); it!=t.end();it++)
	{
		WritePath(out,start_node,*it);
		keep[Nkeep++] = *it;
	}
	//
	bool bCont;
	int start=0;
	int end=Nkeep;
	int Nedges=0;
	double d2_sum=0.0;
	do
	{
		bCont=false;
		for (int i=start;i<end;i++)
		{
			const int curr_id = keep[i];
			bvisited [curr_id] = true;
			const list_t& neighb = tree[curr_id];
			const double* curr_dist = dist[curr_id];
			for (list_t::const_iterator it=neighb.begin();
					it!=neighb.end();it++)
			{
				const int id= *it;
				const double& d=curr_dist[id];
				if (!bvisited[id])
				{
					keep[Nkeep++] = id;
					bCont=true;
					WritePath(out,curr_id,id);
					
				}
			}
		}
		start=end;
		end=Nkeep;
	}while (bCont);
}
	
								
void CMST::WriteTree(const int start_node,std::string fname) const
{
	std::ofstream out(fname.c_str(),std::ios::out);
	if (!out) return;
	WriteTree(start_node,out);
	out.close();
}
	
						
ARBFN_Neuron::ARBFN_Neuron() : RBFNeuron() , GNGNeuron()
{
}

ARBFN_Neuron::ARBFN_Neuron(const ARBFN_Neuron& _c) : RBFNeuron(_c),
	GNGNeuron(_c), Neuron(_c)
{
}
ARBFN_Neuron::ARBFN_Neuron(const Vector<double>& _c) : RBFNeuron(_c),
	GNGNeuron(_c), Neuron(_c)
{
}

ARBFN_Neuron::ARBFN_Neuron(const Vector<double>& _c, const double _r):
		RBFNeuron(_c,_r), GNGNeuron(_c), Neuron(_c)
{
//	RBFNeuron::c.assign(_c);
}

ARBFN_Neuron::~ARBFN_Neuron()
{
}


void ARBFN_Neuron::SetEqualTo( const ARBFN_Neuron& _c)
{
	r=_c.r;
	RBFNeuron::c.assign( _c.GetCenterConst() );
	connect.clear();
	for (std::set<int>::const_iterator it=_c.connect.begin();
			it!=_c.connect.end();it++)
		connect.insert(*it);
}

	
Vector<double>& ARBFN_Neuron::GetCenter()
{
	return RBFNeuron::GetCenter();
}

const Vector<double>& ARBFN_Neuron::GetCenterConst() const
{
	return RBFNeuron::GetCenterConst();
}

void ARBFN_Neuron::SetDimension(const int idim)
{
	RBFNeuron::SetDimension(idim);
}

void ARBFN_Neuron::Save(std::ostream& out) const
{
	//Neuron part
	const int varN=c.size();
	bin_write(out,varN);
	///////////////////////////////////////////////
	for (int i=0;i<varN;i++) bin_write(out,c[i]);
	//GNGPart
	const int conN=connect.size();
	bin_write(out,conN);
	for (std::set<int>::const_iterator i=connect.begin();
		i!=connect.end();i++)
		bin_write(out,*i);
	//RBF Part
	bin_write(out,r);
}
		
bool ARBFN_Neuron::Restore(std::istream& in) 
{
	int varN,conN;
	varN=conN=0;
	//Neuron part
	bin_read(in,varN); 
	c.resize(varN);
	for (int i=0;i<varN;i++) bin_read(in,c[i]);
	//GNGPart
	bin_read(in,conN); connect.clear();
	for (int i=0;i<conN;i++)
	{
		int temp;
		bin_read(in,temp);
		connect.insert(temp);
	}
	//RBF Part
	bin_read(in,r);
	return istream_operational(in);
}

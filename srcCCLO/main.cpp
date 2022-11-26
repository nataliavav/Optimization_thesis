//CCLO
//Exw afhsei mesa elegxous ws sxolia

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <windows.h>

#include "adv_ann.h"
#include "arbfn.h"
#include "eigen.h"
#include "filesystem.h"
#include "gng.h"
#include "indiv.h"
#include "matrix.h"
#include "mymath.h"
#include "param.h"
#include "rng.h"
#include "svd.h"
#include "util.h"

using namespace std;

//Gia ta metamodels
METAMODELS::GNGConfig gngconfig1(std::cout);
int grbfn_min_cen;	//Minimum number of centers
int grbfn_max_cen;	//Maximum number of centers
double grbfn_rad_mult;	//Radius multiplier
double grbfn_tp_testr;	//Percentage of neigh
int grbfn_idle;		//Idle growth iterations
double grbfn_learn_r;	//Learning rate
METAMODELS::TrainingPattern *tp;
int ipca;
double rbfrad;
int intapprox;
Matrix<double> IFs;
METAMODELS::ARBFN rbfnet;
bool binterp;
int inondim;
Vector<double> v1,v2,o1,o2;

int objN;
int maxEval;
int HMS;
double HMCR;
double PAR;
double bw;
int gen_pop;
int elitemax;
Parameterization param;
int varN;
int conN;
vector<int>bits;
vector<double>bmin;
vector<double>bmax;
int run_times;
vector<double> gener_num;
double ge;
int patmin;
int pool_s;
int use_m;
int evals;
double pers;




double AdomB(const vector<double> A, const vector<double> B, const std::vector<double> &ext) {
    double sum=0.5-double(objN);
    for (int k=0;k<objN;k++)
        sum+=Hc((B[k]-A[k])/ext[k],0);
    return Hc(sum,0);
}

void SPEA2(const int dbN, vector<Individual>& db){
    /*for (int i=0;i<dbN;i++){                            //N Now
        //if(db[i].raw==0) {
            cout << "SPEA " << i + 1 << "  ";
            for (int j = 0; j < objN; j++) {
                cout << db[i].obj[j] << "  ";
            }
            cout << endl;
        //}
    }*/
    if (objN==1)
    {
        for (int idC=0; idC<dbN; idC++)
        {
            db[idC].phi=db[idC].pobj[0];
            db[idC].raw=db[idC].pobj[0];
            if (db[idC].flagEval == 2)
            {
                db[idC].phi=dmax;
                db[idC].raw=dmax;
            }
        }
        return;
    }
    //
    // Find the distance
    // =================
    vector<double> ext; ext.resize(objN);
    vector<double> mini; mini.resize(objN);
    for (int k=0;k<objN;k++)
    { ext[k]=-1.e10; mini[k]=1.e10 ; }

    vector<bool> isFailed; isFailed.resize(dbN); 	// vera Sept2016
    for (int id=0;id<dbN;id++) isFailed[id] = false;

    for (int id=0;id<dbN;id++)
    {
        if (db[id].flagEval == 2) continue;	// skip failed
        for (int k=0;k<objN;k++)
        {
            const double pobj = db[id].pobj[k];
            if (fabs(pobj) >1.e10)		// vera Sept2016
            {
                isFailed[id]=true;
                continue;
            }

            if (pobj<mini[k]) mini[k]=pobj;
            if (pobj> ext[k])  ext[k]=pobj;

        }
    }
    // vera Sept2016
    //int nFailed=0;
    //for (int id=0;id<dbN;id++)
    //	if (isFailed[id]) nFailed++ ;
    //cout<<nFailed<<" failed out of "<<dbN<<"\n";


    for (int k=0;k<objN;k++)
    {
        ext[k]-=mini[k];
        if (dabs(ext[k])<deps) ext[k]=1.0;
        //cout<<k<<" "<<ext[k]<<endl;
    }
    const int kth = static_cast<int>(sqrt(static_cast<float>(1+dbN/3)));
    //=============================
    for (int idC=0; idC<dbN; idC++)
        //=============================
    {

        if ( (db[idC].flagEval == 2 || isFailed[idC]) && db[idC].flagEval!=0 )  	// vera Sept2016
        {
            //cout<<" Skip failed (A) "<<idC<<"\n";
            db[idC].phi=dmax;
            db[idC].raw=dmax;
            continue;		// skip failed
        }

        vector<double> dist(dbN);
        indexx my_index;
        for (int id=0;id<dbN;id++)
        {
            if (db[id].flagEval == 2 || isFailed[id]) 	// vera Sept2016
            {
                dist[id]=dmax;
                //cout<<" Skip failed (B) "<<id<<endl;
                continue;	// skip failed
            }

            if (idC!=id)
            {
                double sum=0;
                for (int k=0;k<objN;k++)
                {
                    const double d=(db[idC].pobj[k] - db[id].pobj[k])/ ext[k];
                    sum+=d*d;
                }
                dist[id]=sum;
            }
            else
                dist[id]=dmax;
        }
        vector<int>    wksp(dbN);
        my_index(dbN,dist,wksp);
        const int kthid=wksp[kth];
        const double densDenom= 2.0 + sqrt(dist[kthid]+deps);
        const double density = 1.0/densDenom;

        double fi=0;
        double strength=0;
        for (int i=0;i<dbN;i++)
        {
            //if (db[i].flagEval == 2) continue;	// skip failed
            if (db[i].flagEval == 2 || isFailed[i]) continue; 	// vera Sept2016

            if (i==idC) continue;

            //Calculate the strength of i
            for (int j=0;j<dbN;j++)
            {
                if (db[j].flagEval == 2 || isFailed[j]) continue;	// skip failed
                if (i!=j) strength+=AdomB(db[i].pobj,db[j].pobj,ext);
            }
            strength/=double(dbN);
            const double bdom=AdomB(db[i].pobj,db[idC].pobj,ext);

            fi+=strength*bdom;
        }

        db[idC].raw = fi;
        db[idC].phi = fi+density;
    }
}


void tempDb(vector<Individual> &db_temp, int &db_temp_s, vector<Individual> &off, int off_s, vector<Individual> &elite,
            int eliteN, vector<Individual> &db) {
    if (eliteN>HMS){
        db_temp_s=off_s+eliteN;
    }
    else{
        db_temp_s=off_s+HMS;
    }
    db_temp.clear();
    db_temp.resize(db_temp_s);
    //cout<<db_temp_s<<"  "<<off_s<<"  "<<eliteN<<"  "<<HMS<<endl;

    if (db_temp_s==off_s+eliteN){
        for (int i=0;i<off_s;i++){
            db_temp[i].copy(off[i]);
        }
        for (int i=0;i<eliteN;i++){
            db_temp[off_s+i].copy(elite[i]);
        }
    }
    else{
        for (int i=0;i<off_s;i++){
            db_temp[i].copy(off[i]);
        }
        for (int i=0;i<HMS;i++){
            db_temp[off_s+i].copy(db[i]);
        }
    }
}

void ANNPredict(Vector<double>& p, Vector<double>& o)
{
    //NonDim vector
    for (int i=0;i<varN;i++)
        p[i]=(p[i]-v1[i])/(v2[i]-v1[i]);
    rbfnet.Predict(p,o);
    //Redim prediction
    for (int i=0;i<objN;i++)
        o[i]=o[i]*(o2[i]-o1[i])+o1[i];
}

bool CheckConfig()
{
    if (grbfn_min_cen>grbfn_max_cen) throw Default_exception
                ("GRBFN: Min.centers should be <= max.centers\n");
    if (grbfn_min_cen>=int(patmin*(1.0-grbfn_tp_testr)))
        throw Default_exception
                ("GRBFN: Training minus testing patterns should be more than minimum centers\n");
    if (grbfn_rad_mult<=0) throw Default_exception
                ("GRBFN: Radius multiplier should be>0\n");
    if (grbfn_tp_testr<=0 || grbfn_tp_testr>=1) throw Default_exception
                ("GRBFN: Testing-to-total patterns ratio should be between 0,1\n");
    if (grbfn_idle<=1) throw Default_exception
                ("GRBFN: Idle growth iterations should be >1\n");
    if (grbfn_learn_r<=0 ||grbfn_learn_r>1) throw Default_exception
                ("GRBFN: Learning rate fraction does not lie in (0,1]\n");
    const bool b1=gngconfig1.CheckConfig();
    if (!b1) return false;
    return true;
}

double RBF_RAD()
{
    std::cout<<"Computing Radius ..."; std::cout.flush();
    double dmax_N=-1.0; double sum=0.0;
    for (int i=0;i<patmin-1;i++)
        for (int j=i+1;j<patmin;j++)
        {

            sum=0.0;
            for (int k=0;k<varN;k++){
                sum+= (tp[i].var[k]-tp[j].var[k])* (tp[i].var[k]-tp[j].var[k]);
            }
            dmax_N=max1(dmax_N,(sum));
        }
    std::cout<<" OK "<<std::endl;
    return sqrt(dmax_N/(2.0*double(varN)));
}

void IFsPCA()
{
    Matrix<double> xx; xx.Resize(varN,varN,0);
    Vector<double> mm; mm.resize(varN); for(int i=0;i<varN;i++) mm[i]=0;
    //Compute mean
    for (int i=0;i<patmin;i++)
        for (int j=0;j<varN;j++)
            mm[j]+=tp[i].var[j];
    for (int j=0;j<varN;j++) mm[j]/=double(patmin);
    //Compute covariance matrix
    for (int i=0;i<patmin;i++)
    {
        for (int j=0;j<varN;j++)
            for(int k=0;k<varN;k++)
                xx[j][k]+=(tp[i].var[j]-mm[j])*
                          (tp[i].var[k]-mm[k]);
    }
    for (int i=0;i<varN;i++)
        for (int j=0;j<varN;j++)
            xx[i][j]/=double(patmin-1);
    Vector<double> ls;
    Matrix<double> eig;
    const bool b1=RealSymEigen(xx,ls,eig,100);
    eig.Dump(std::cout);
    if (!b1)
    {
        IFs.Resize(varN,1.0/double(varN));
        return;
    }
    //Find the largest eigenvalue
    double lmax=dabs(ls[0]);
    int ibest=0;
    for (int i=0;i<varN;i++) std::cout <<"eig("<<i<<")="<<ls[i]<<"\n";
    for (int i=1;i<varN;i++)
    {
        if (dabs(ls[i])>lmax)
        {
            lmax=dabs(ls[i]);
            ibest=i;
        }
    }
    lmax=0;
    //Compute IFs
    for (int i=0;i<varN;i++)
    {
        const double comp=dabs(eig[ibest][i]);
        for (int j=0;j<objN;j++)
            IFs[i][j]=comp;
        lmax+=comp;
    }
    if (lmax<deps)
    {
        IFs.Resize(varN,1.0/double(varN));
        return;
    }
    for (int i=0;i<varN;i++)
        for (int j=0;j<objN;j++)
            IFs[i][j]/=lmax;
}

void CalcBoundaries()
{
    //Step 1: Initialization
    v1.assign(tp[0].var);
    v2.assign(tp[0].var);
    o1.assign(tp[0].obj);
    o2.assign(tp[0].obj);
    //Step 2: Determine Limits
    for (int ipat=0;ipat<patmin;ipat++)
    {
        //a: Vars
        for (int ivar=0;ivar<varN;ivar++)
        {
            v1[ivar]=min1(v1[ivar],tp[ipat].var[ivar]);
            v2[ivar]=max1(v2[ivar],tp[ipat].var[ivar]);
        }
        //b: Objs
        for (int iobj=0;iobj<objN;iobj++)
        {
            o1[iobj]=min1(o1[iobj],tp[ipat].obj[iobj]);
            o2[iobj]=max1(o2[iobj],tp[ipat].obj[iobj]);
        }
    }

    //Step 3: Check limits
    //a: Vars
    for (int ivar=0;ivar<varN;ivar++)
        if (dabs( v1[ivar]-v2[ivar]) <deps)
        {
            v1[ivar]=0.0;
            v2[ivar]=1.0;
        }
    //b: Objs
    for (int iobj=0;iobj<objN;iobj++)
        if (dabs( o1[iobj]-o2[iobj]) <deps)
        {
            o1[iobj]=0.0;
            o2[iobj]=1.0;
        }
    for (int ivar=0;ivar<varN;ivar++)
        std::cout <<"id="<<ivar+1<<"   min="<<v1[ivar]
                  <<"   max="<<v2[ivar]<<"\n";
    for (int iobj=0;iobj<objN;iobj++)
        std::cout <<"id="<<iobj+1<<"   min="<<o1[iobj]
                  <<"   max="<<o2[iobj]<<"\n";

}

bool NonDim()
{
    //Step 1: Non-dimensionalize
    //--Patterns
    for (int ipat=0;ipat<patmin;ipat++)
    {
        //a: Vars
        for (int ivar=0;ivar<varN;ivar++)
        {
            const double e=v2[ivar]-v1[ivar];
            tp[ipat].var[ivar]=(tp[ipat].var[ivar]-v1[ivar])/e;
        }
        //b: Objs
        for (int iobj=0;iobj<objN;iobj++)
        {
            const double e=o2[iobj]-o1[iobj];
            tp[ipat].obj[iobj]=(tp[ipat].obj[iobj]-o1[iobj])/e;
        }
    }
    return true;
}

void getNonDominated(vector <Individual> & nDom,const int nonDom, const int nData, vector<Individual> indVec)
{
    nDom.clear();
    nDom.resize(nonDom);
    int ind=0;
    for (int id=0; id<nData; id++)
    {
        if (fabs(indVec[id].raw) < 1.e-13)
        {
            nDom[ind].copy(indVec[id]);
            ind++;
        }
    }
}

int countNonDominated(const int nData, vector<Individual>& indVec)
{
    int nonDom = 0;
    for (int id=0; id<nData; id++)
    {
        if (fabs(indVec[id].raw) < 1.e-13) nonDom++;
    }
    return nonDom;
}

void best_ind(vector<Individual>& db,vector<Individual>& pool,int pool_s1){
    int maxi=-1;
    double maxdi=-numeric_limits<double>::max();
    vector <int> mindbi;
    mindbi.clear();
    for (int i=0;i<HMS;i++){
        mindbi.push_back(i);
        if (pool[i].phi>maxdi){
            maxdi=pool[i].phi;
            maxi=i;
        }
    }
    for (int i=HMS;i<pool_s1;i++){
        if (pool[i].phi<maxdi){
            mindbi[maxi]=i;
            maxdi=-numeric_limits<double>::max();
            for (int l=0;l<HMS;l++){
                if (pool[mindbi[l]].phi>maxdi){
                    maxdi=pool[mindbi[l]].phi;
                    maxi=l;
                }
            }
        }
    }
    db.clear();
    db.resize(HMS);
    for (int i=0;i<HMS;i++){
        db[i].copy(pool[mindbi[i]]);
    }
}

vector<int> min_i(int patmin1, vector<double> distance,int pool_size ){

    int maxi=-1;
    double maxdi=-numeric_limits<double>::max();
    vector<int> mini;
    mini.clear();

    for (int i=0;i<patmin1;i++){
        mini.push_back(i);
        if (distance[i]>maxdi){
            maxdi=distance[i];
            maxi=i;
        }
    }

    for (int i=patmin1;i<pool_size;i++){
        if (distance[i]<maxdi){
            mini[maxi]=i;
            maxdi=distance[i];
            for (int l=0;l<patmin1;l++){
                if (distance[mini[l]]>maxdi){
                    maxdi=distance[mini[l]];
                    maxi=l;
                }
            }
        }
    }
    return mini;
}

vector<int> min_i_ind(int patmin1, vector<Individual>& off,int pool_size ){
    int maxi=-1;
    double maxdi=-numeric_limits<double>::max();
    vector<int> mini;
    mini.clear();

    for (int i=0;i<pool_size;i++){
        off[i].flagEval=2;
    }

    for (int i=0;i<patmin1;i++){
        mini.push_back(i);
        if (off[i].phi>maxdi){
            maxdi=off[i].phi;
            maxi=i;
        }
    }

    for (int i=patmin1;i<pool_size;i++){
        if (off[i].phi<maxdi){
            mini[maxi]=i;
            maxdi=-numeric_limits<double>::max();
            for (int l=0;l<patmin1;l++){
                if (off[mini[l]].phi>maxdi){
                    maxdi=off[mini[l]].phi;
                    maxi=l;
                }
            }
        }
    }

    for (int i=0;i<patmin1;i++){
        off[mini[i]].flagEval=0;
    }

    return mini;
}

void ExactEvaluations(vector<Individual>& db, int dbN ,int& ev1)
{
    // Perform exact evaluations
    ev1+=dbN;
    db.resize(dbN);
    for (int id=0; id<dbN; id++)
    {
        db[id].obj.resize(objN);
        db[id].con.resize(conN);
        //write task.dat
        std::ofstream out("task.dat",std::ios::out);
        out<<varN<<endl;
        for (int iv=0; iv<varN; iv++)
        {
            out<<db[id].var[iv]<<endl;
        }
        out.close();
        //
        //call task.bat
        //std::string scriptEval=".\objective.exe";
        system("objective.exe");   // bale ./task.bat prin to anebaseis
        //
        //read task.res
        //
        std::ifstream res("task.res",std::ios::in);
        db[id].flagEval = 0;
        db[id].flag = 0;

        if (!res)
        {
            res.close();
            db[id].flagEval = 1;
            for (int ib=0; ib<objN; ib++) db[id].obj.push_back(1.e50);
        }
        else
        {
            for (int ib=0; ib<objN; ib++)
            {
                double dummyR;
                res >> dummyR;
                db[id].obj[ib] = dummyR;
            }
            res.close();
        }
        //read task.cns
        std::ifstream cns("task.cns",std::ios::in);
        if(conN>0)
        {
            if (!cns)
            {
                cns.close();
                db[id].flagEval = 1;
                for (int ic=0; ic<conN; ic++) db[id].con.push_back(-1.e50);
            }
            else
            {
                for(int ic=0; ic<conN; ic++)
                {
                    double dummyC;
                    cns >>dummyC;
                    db[id].con[ic] = dummyC;
                }
                cns.close();
            }
        }
        //
        // apply penalty to compute penalized objective
        bool pen,death;

        param.Penalize (db[id].con,db[id].obj,db[id].pobj,pen,death);
    }



}

vector <int> trimElites(const int nonDom, vector<Individual>& nDm)
{
    //int* age = new int[nonDom];
    vector <int> age(nonDom);
    for (int id=0;id<nonDom;id++) age[id] = 0;

    // Trim extra members....
    int el2Trim = nonDom - elitemax;
    while (el2Trim > 0)
    {
        int im1=-1; int im2=-1;
        double distMin=1.e10;
        for (int ie1=0;     ie1<nonDom; ie1++)
        {
            if (age[ie1] < 0) continue;
            for (int ie2=ie1+1; ie2<nonDom; ie2++)
            {
                if (age[ie2] < 0) continue;
                double d=0;
                for (int k=0; k<objN; k++)
                {
                    const double ee = ( nDm[ie1].obj[k]-
                                        nDm[ie2].obj[k] );
                    d +=ee*ee;
                }
                if (d < distMin)
                {
                    distMin=d;  im1=ie1; im2=ie2;
                }
            }
        }

        double dist1=distMin;
        double dist2=distMin;
        for (int ie=0; ie<nonDom; ie++)
        {
            if (ie==im1 || ie==im2 || age[ie]<0) continue;
            double d1=0;
            double d2=0;
            for (int k=0; k<objN; k++)
            {
                const double ee1 = ( nDm[ie ].obj[k]-
                                     nDm[im1].obj[k] );
                d1 +=ee1*ee1;
                const double ee2 = ( nDm[ie ].obj[k]-
                                     nDm[im2].obj[k] );
                d2 +=ee2*ee2;
            }
            if (d1 < dist1) dist1=d1;
            if (d2 < dist2) dist2=d2;
        }
        int irem=-1;
        if (dist1 < dist2) irem = im1;
        else               irem = im2;
        age[irem] = -10;
        el2Trim--;
    }

    return age;
}

void Readcfg(const string& eaConfFile)
{
    cout<<"Reading configuration file...";

    ifstream conf(eaConfFile.c_str(),ios::in);

    skipLines(conf,4);

    conf>>objN;
    istream_nl(conf);

    conf>>maxEval;
    istream_nl(conf);

    skipLines(conf,3);

    conf>>HMS;
    istream_nl(conf);

    conf>>HMCR;
    istream_nl(conf);

    conf>>PAR;
    istream_nl(conf);

    conf>>bw;
    istream_nl(conf);

    conf>>gen_pop;
    istream_nl(conf);

    conf>>elitemax;
    istream_nl(conf);

    skipLines(conf,4);

    if (!param.ReadParameterization(conf)) cout<<"\n Error while reading Parameterization data\n"<<endl;
    varN=param.getVarN();
    istream_nl(conf);	// vera

    if (!param.ReadConstraints(conf)) cout<<"\n Error while reading Constraint data\n"<<endl;
    conN=param.getConN();

    for (int iv=0; iv<varN; iv++){
        int bit     = param.dvars[iv].bits;
        double vmin = param.dvars[iv].min;
        double vmax = param.dvars[iv].max;

        if (bit == 1) bit=0;
        if (bit == 0) vmax=vmin;

        bits.push_back(bit);
        bmin.push_back(vmin);
        bmax.push_back(vmax);
    }

    skipLines(conf,3);

    conf>>run_times;
    istream_nl(conf);

    for (int n=0;n<run_times;n++){
        conf>>ge;
        istream_nl(conf);
        gener_num.push_back(ge);
    }

    skipLines(conf,3);

    conf>>patmin;
    istream_nl(conf);

    conf>>pool_s;
    istream_nl(conf);

    conf>>pers;
    istream_nl(conf);

    conf>>use_m;
    istream_nl(conf);

    if (use_m==1) pool_s=HMS;

    conf.close();

    cout<<"OK"<<endl;

}


void Readmetamodels(){
    cout<<"Readnig metamodels configuration...";
    std::ifstream conf("ann.conf",std::ios::in);

    conf>>ipca;
    istream_nl(conf);

    conf>>rbfrad;
    istream_nl(conf);

    conf>>inondim;
    istream_nl(conf);

    conf>>intapprox;
    istream_nl(conf);

    conf>>grbfn_min_cen;
    istream_nl(conf);

    conf>>grbfn_max_cen;
    istream_nl(conf);

    conf>>grbfn_rad_mult;
    istream_nl(conf);

    conf>>grbfn_tp_testr;
    istream_nl(conf);

    conf>>grbfn_idle;
    istream_nl(conf);

    conf>>grbfn_learn_r;
    istream_nl(conf);

    gngconfig1.ReadConfig(conf);
    if (!CheckConfig() || !istream_operational(conf)){
        std::cout <<"Error in configuration!\n";
    }
    conf.close();

    cout<<"OK"<<endl;
}



int main(int argc, char* argv[])
{
    if (argc<2){
        cout<<"Expected: <program.exe> <input_file>"<<endl;
        exit(-1);
    }

    string eaConfFile = argv[1];
    Readcfg(eaConfFile);


    vector <vector <double> > objective2;
    objective2.clear();
    vector <double> temp;
    int pool_st=pool_s;
    int HMS_st=HMS;
    int gen_pop_st=gen_pop;

    for (int g=0;g<run_times;g++){
        pool_s=pool_st;
        HMS=HMS_st;
        gen_pop=gen_pop_st;
        //Prepare Random Number Generator
        int generator=gener_num[g];
        RandUnif myrand;
        for (int i=0;i<generator;i++) myrand(); //Warm it up

        //Dhmiourgia neas database kai arxikopoihsh metablhtwn
        cout<<"Creating random variables...";
        vector<Individual> pool;
        pool.clear();
        pool.resize(pool_s);

        for (int i=0;i<pool_s;i++){
            pool[i].var.resize(varN);
            for (int j=0;j<varN;j++){
                pool[i].var[j]=bmin[j]+myrand()*(bmax[j]-bmin[j]);
                //cout<<endl<< "stoixeio "<<i <<" metablhth "<<j<<" timh "<<pool[i].var[j]<<endl;
            }
        }
        cout<<"OK"<<endl;

        //Ypologismos antikeimenikhs (obj,pobj,phi)
        cout<<"Calculating objectve...";
        evals=0;
        ExactEvaluations(pool, pool_s, evals);
        SPEA2(pool_s,pool);
        /*cout<<endl<<pool_s<<endl;
        for (int i=0;i<pool_s;i++){
            for (int j=0;j<objN;j++){
                cout<<i+1<<"  "<< pool[i].obj[j] <<"   "<<pool[i].pobj[j] <<"   "<<pool[i].phi <<endl;
            }
        }*/
        cout<<"ok"<<endl;

        //Dhmioyrgia arxeiou, euresh min, 1 antikeimenikh
        stringstream result;
        result<<"RNG="<<generator<<".txt";
        ofstream myfile(result.str().c_str(),ios::out);

        double minN=numeric_limits<double>::max();
        int min_pos=-1;
        vector <double> min_var;
        double min_obj;
        if (objN==1){
            for (int id=0;id<pool_s;id++){
                if (pool[id].phi<minN){
                    minN = pool[id].phi;
                    min_pos=id;
                    min_var.clear();
                    for (int v=0;v<varN;v++){
                        min_var.push_back(pool[id].var[v]);
                    }
                }
                myfile<<id+1<<"   "<<pool[min_pos].phi<<"     ";
                for (int v=0;v<varN;v++){
                    myfile<<min_var[v]<<"  ";
                }
                myfile<<endl;
                /*"   ";
                for (int v=0;v<varN;v++){
                    myfile<<pool[min_pos].var[v]<<"   ";
                }
                myfile<<endl;*/
                temp.push_back(pool[min_pos].phi);
            }
            min_obj=pool[min_pos].obj[0];

        }

        //Find Elites and create db for HS
        vector<Individual> db(HMS);
        int eliteN;
        int nDmN;
        vector <Individual> nDm;
        vector <Individual> elite;
        vector <int> age;
        if (objN>1){
            cout<<"finding elites...";
            /*for (int i=0;i<pool_s;i++){
                cout<<"pool "<<i+1<<"  "<<pool[i].raw<<endl;
            }*/
            nDmN = countNonDominated(pool_s,pool);
            cout<<"nDmN= "<<nDmN<<endl;
            getNonDominated(nDm, nDmN,pool_s,pool);
            age = trimElites(nDmN,nDm);
            eliteN=0;
            for (int ie=0;ie<nDmN; ie++){
                if (age[ie] < 0) continue;
                eliteN++;
                elite.push_back(nDm[ie]);
            }
            /*for (int i=0;i<eliteN;i++){
                cout<<"elite "<<i+1<<"  "<<elite[i].raw<<endl;
            }*/
            cout<<"ok"<<endl;
        }
        if (use_m==0){
            //kalw synarthsh poy briskei ta hms kalytera
            cout<<"ypologismos arxikhs db...";
            best_ind(db, pool, pool_s);
            /*for (int i=0;i<HMS;i++){
                cout<<"db "<<i+1<<"  "<<db[i].raw<<endl;
            }*/
            cout<<"ok"<<endl;

            //Diabase arxeio dedomenwn gia metamodela
            Readmetamodels();
        }
        else{
            for (int i=0;i<HMS;i++){
                db[i].copy(pool[i]);
            }
            /*for (int i=0;i<HMS;i++){
                cout<<"db "<<i+1<<"  "<<db[i].raw<<endl;
            }*/
        }

        ////////////////////////////////////////////////////////
        //////////////////////EPANALHPSEIS//////////////////////
        ////////////////////////////////////////////////////////

        vector<Individual> off;
        vector<Individual> db_temp;
        vector<Individual> off2;
        vector<Individual> pool_temp;
        int db_temp_s;
        int a_spea;

        while (evals<maxEval){

            off.clear();
            off.resize(gen_pop);

            //Dhmioyrgia newn apogwnwn me thn methodo CCLO
            for (int gp=0;gp<gen_pop;gp++){
                off[gp].var.resize(varN);
                for (int j=0; j<varN; j++){
                    if (myrand()<=HMCR){
                        int i=myrand(0,HMS-1);
                        int i_1=myrand(0,HMS-1);
                        double i_2=myrand();
                        off[gp].var[j]=i_2*db[i].var[j]+(1.0-i_2)*db[i_1].var[j];
                        if((myrand()<=PAR)){
                            off[gp].var[j]=off[gp].var[j]=off[gp].var[j]+(-(bw/2.0)+myrand()*bw)*exp(-1.0 *evals/(1.0 *maxEval));;
                            while (off[gp].var[j]>bmax[j] || off[gp].var[j]<bmin[j]){
                                cout<<bmin[j]<<"  "<<off[gp].var[j]<<"  "<<bmax[j]<<endl;
                                if (off[gp].var[j]>bmax[j]) {
                                    off[gp].var[j]=off[gp].var[j]-(off[gp].var[j]-bmax[j])-(off[gp].var[j]-bmax[j])*myrand();
                                }
                                if(off[gp].var[j]<bmin[j]) {
                                    off[gp].var[j]=off[gp].var[j]+(bmin[j]-off[gp].var[j])+(bmin[j]-off[gp].var[j])*myrand();
                                }
                                cout<<off[gp].var[j]<<endl;
                            }
                        }
                    }
                    else{
                        off[gp].var[j]=bmin[j]+myrand()*(bmax[j]-bmin[j]);
                    }

                }

                //An thelw neo atomo me HS
                /*for (int gp=0;gp<gen_pop;gp++){
				off[gp].var.resize(varN);
				for (int j=0; j<varN; j++){
	 	   			if (myrand()<=HMCR){
	 	       			int i=myrand(0,HMS-1);
	            		off[gp].var[j]=db[i].var[j];
	           			if((myrand()<=PAR)){
	                		off[gp].var[j]=off[gp].var[j]-(bw/2.0)+myrand()*bw;
	                		if (off[gp].var[j]>bmax[j]) off[gp].var[j]=bmax[j];
	                		if(off[gp].var[j]<bmin[j]) off[gp].var[j]=bmin[j];
	            		}
	        		}
	        		else{
	        			off[gp].var[j]=bmin[j]+myrand()*(bmax[j]-bmin[j]);
	        		}

	   			}*/

                //Xrhsh twn metamodelwn
                if (use_m==0){

                    //Eyresh apostashs neoy atomoy apo kathe stoixeio toy pool
                    vector <double> distance;
                    distance.resize(pool_s);
                    for (int i=0;i<pool_s;i++){
                        distance[i]=0;
                        for (int j=0;j<varN;j++){
                            distance[i]+= (off[gp].var[j]-pool[i].var[j])*(off[gp].var[j]-pool[i].var[j]);
                        }
                    }
                    vector<int> mini=min_i(patmin,distance,pool_s);

                    //METAMODELS TRAIN
                    tp=new METAMODELS::TrainingPattern[patmin];
                    for (int i=0;i<patmin;i++){
                        tp[i].var.resize(varN);
                        tp[i].obj.resize(objN);
                        for (int k=0;k<varN;k++){
                            tp[i].var[k]=pool[mini[i]].var[k];
                        }
                        for (int k=0;k<objN;k++) tp[i].obj[k]=pool[mini[i]].pobj[k];
                        tp[i].bfailed=false;
                    }

                    IFs.Resize(varN,1.0/double(varN));
                    o1.resize(objN);  o1=0;
                    o2.resize(objN);  o2=1;
                    v1.resize(varN);  v1=0;
                    v2.resize(varN);  v2=1;

                    if (inondim==1){
                        CalcBoundaries();
                        NonDim();
                    }

                    if (rbfrad<0) rbfrad=RBF_RAD();

                    if (ipca==1) IFsPCA();
                    binterp=(intapprox==0);

                    if (binterp) rbfnet.SetTrainingInterpolation();

                    else rbfnet.SetTrainingApproximation (grbfn_min_cen, grbfn_max_cen, grbfn_tp_testr, &gngconfig1, grbfn_rad_mult, grbfn_learn_r, grbfn_idle);
                    if (binterp) rbfnet.SetRadius(rbfrad);
                    if (ipca==1){
                        rbfnet.SetIFs(IFs);
                        rbfnet.SetIFUsage(true);
                    }
                    rbfnet.Build(tp,patmin,varN,objN);
                    rbfnet.Train();
                    int vdat=varN;
                    Vector<double> p,o;
                    p.resize(varN); o.resize(objN);
                    for (int i=0;i<vdat;i++) p[i]=off[gp].var[i];
                    ANNPredict(p,o);

                    off[gp].obj.resize(objN);
                    off[gp].pobj.resize(objN);
                    off[gp].flagEval=0;
                    for (int i=0;i<objN;i++){
                        off[gp].pobj[i]=o[i];
                    }
                }
            }


            if (use_m==0){
                /*for (int i=0;i<gen_pop;i++){
                    for (int j=0;j<objN;j++){
                        cout<<off[i].pobj[j]<<"  ";
                    }
                    cout<<endl;
                }*/
                SPEA2(gen_pop,off);

                //Briskw ta kalytera atoma
                if ((evals+(int)round(pers/100.0*gen_pop))<=maxEval){
                    a_spea=(int)round(pers/100.0*gen_pop);
                }
                else {
                    a_spea=maxEval-evals;
                }

                /*for (int i=0;i<gen_pop;i++){
                    cout<<"off "<<i+1<<"  "<<off[i].raw<<endl;
                }*/
                vector<int> best_met=min_i_ind(a_spea,off,gen_pop);

                //Ypologizw pragmatika ta kalytera
                off2.resize(a_spea);

                for (int i=0;i<a_spea;i++){
                    off2[i].copy(off[best_met[i]]);
                }

                /*for (int i=0;i<a_spea;i++){
                    cout<<"off2 "<<i+1<<"  "<<off2[i].raw<<endl;
                }*/


                /*cout<<"before"<<endl;
                for (int i=0;i<a_spea;i++){
                    cout<<i+1<< "  ";
                    for (int j=0;j<objN;j++){
                        cout<<off2[i].pobj[j]<<"  ";
                    }
                    cout<<endl;
                }*/

                ExactEvaluations(off2,a_spea,evals) ;

                /*cout<<"after"<<endl;
                for (int i=0;i<a_spea;i++){
                    cout<<i+1<< "  ";
                    for (int j=0;j<objN;j++){
                        cout<<off2[i].pobj[j] <<"  ";
                    }
                    cout<<endl;
                }*/

                //ftiaxnw dbtemp me atoma poy stelnw gia aksiologhsh apo thn SPEA
                pool_temp.clear();
                pool_temp.resize(pool_s);
                for (int i=0;i<pool_s;i++){
                    pool_temp[i].copy(pool[i]);
                }
                pool.clear();
                pool.resize(pool_s+a_spea);
                for (int i=0;i<pool_s;i++){
                    pool[i].copy(pool_temp[i]);
                }
                for (int i=0;i<a_spea;i++){
                    pool[pool_s+i].copy(off2[i]);
                }
                pool_s+=a_spea;

                if (objN>1) {
                    for (int i = 0; i < a_spea; i++) {
                        //if(db_temp[i].raw==0){
                        /*cout << "off2 " << i + 1 << "  ";
                        for (int j = 0; j < objN; j++) {
                            cout << off2[i].obj[j] << "  ";
                        }
                        cout << endl;*/
                        //}
                    }
                    for (int i = 0; i < HMS; i++) {
                        //if(db_temp[i].raw==0){
                        /*cout << "db " << a_spea+ i + 1 << "  ";
                        for (int j = 0; j < objN; j++) {
                            cout << db[i].obj[j] << "  ";
                        }
                        cout << endl; */
                        //}
                    }

                    tempDb(db_temp, db_temp_s, off2, a_spea, elite, eliteN, db);
                }
                else{
                    db_temp_s=a_spea+HMS;
                    db_temp.clear();
                    db_temp.resize(db_temp_s);
                    for (int i=0;i<a_spea;i++){

                        db_temp[i].copy(off2[i]);
                    }
                    for (int i=0;i<HMS;i++){
                        db_temp[i+a_spea].copy(db[i]);

                    }

                }

            }
            else{

                if (evals+gen_pop>maxEval){
                    gen_pop=maxEval-evals;
                }
                ExactEvaluations(off,gen_pop,evals);
                /*for (int i=0;i<gen_pop;i++){
                    cout<<"off "<<i+1<<"  "<<off[i].raw<<endl;
                }*/

                //ftiaxnw dbtemp me atoma poy stelnw gia aksiologhsh apo thn SPEA
                if (objN>1) {
                    tempDb(db_temp, db_temp_s, off, gen_pop, elite, eliteN, db);
                    /*for (int i=0;i<db_temp_s;i++){
                        cout<<"db_temp "<<i+1<<"  "<<db_temp[i].raw<<endl;
                    }*/
                }
                else{
                    db_temp_s=HMS+gen_pop;
                    db_temp.clear();
                    db_temp.resize(db_temp_s);
                    for (int i=0;i<gen_pop;i++){
                        db_temp[i].copy(off[i]);

                    }
                    for (int i=0;i<HMS;i++) {
                        db_temp[i + gen_pop].copy(db[i]);

                    }
                }
            }

            //Aksiologhsh ths db_temp

            SPEA2(db_temp_s,db_temp);

            /*for (int i=0;i<db_temp_s;i++){
                if(db_temp[i].raw==0){
                    cout<<"db_temp "<<i+1<<"  ";
                    for (int j=0;j<objN;j++){
                        cout<<db_temp[i].obj[j]<<"  ";
                    }
                    cout<<endl;
                }
            }*/



            if (objN>1){
                //briskw nonDominated, kanw trim osa den thelw kai briskw elites
                elite.clear();
                nDm.clear();
                age.clear();
                eliteN=0;
                nDmN = countNonDominated(db_temp_s, db_temp);
                cout<<"nDmN= "<<nDmN<<endl;
                getNonDominated(nDm, nDmN, db_temp_s, db_temp);
                /*for (int i=0;i<db_temp_s;i++){
                    cout<<i<<"  "<<db_temp[i].raw<<endl;
                }*/
                age = trimElites(nDmN,nDm);
                for (int ie=0;ie<nDmN; ie++){
                    if (age[ie] < 0) continue;
                    eliteN++;
                    elite.push_back(nDm[ie]);
                }
                /*for (int i=0;i<eliteN;i++){
                    if(elite[i].raw==0){
                        cout<<"elite "<<i+1<<"  ";
                        for (int j=0;j<objN;j++){
                            cout<<elite[i].obj[j]<<"  ";
                        }
                        cout<<endl;
                    }
                }*/


                /*cout<< "PRIN"<<endl;
                for (int i=0;i<db_temp_s;i++){
                    if(db_temp[i].raw==0){ cout<<"  ";}
                        cout<<"db_temp "<<i+1<<"  "<< db_temp[i].phi<<"  ";
                        for (int j=0;j<objN;j++){
                            cout<<db_temp[i].obj[j]<<"  ";
                        }
                        cout<<endl;

                }*/
                best_ind(db,db_temp,db_temp_s);
                /*cout<<"TELIKH DB"<<endl;
                for (int i = 0; i < HMS; i++) {
                    //if(db_temp[i].raw==0){
                    cout << "db " << i + 1 << "  ";
                    for (int j = 0; j < objN; j++) {
                        cout << db[i].obj[j] << "  ";
                    }
                    cout << endl;
                    //}
                }*/

                db_temp.clear();

                //Grafw nDm se arxeio
                stringstream ss;
                if (use_m==0){
                    ss<<"m"<<evals<<"_RNG="<<generator<<".txt";
                }
                else{
                    ss<<"n"<<evals<<"_RNG="<<generator<<".txt";
                }

                ofstream myfile1(ss.str().c_str());
                for (int i=0;i<eliteN;i++){
                    for (int j=0;j<objN;j++){
                        myfile1 <<elite[i].obj[j]<<"  ";
                    }
                    myfile1 <<"  ";
                    for (int j=0;j<varN;j++){
                        myfile1 <<elite[i].var[j]<<"  ";
                    }
                    myfile1<<endl;
                }
                myfile1.close();
                //if (evals==174){ exit(0);}                      //Elegxos gia na stamataw ton kwdika ekei poy ton thelw
            }
            else{
                //Grafw areia eksodou
                if (use_m==0){
                    for (int id=0;id<round(pers/100.0*gen_pop);id++){
                        if (db_temp[id].phi<minN){
                            minN = db_temp[id].phi;
                            min_pos=id;
                            min_obj=db_temp[id].obj[0];
                            min_var.clear();
                            for (int v=0;v<varN;v++){
                                min_var.push_back(db_temp[id].var[v]);
                            }
                        }

                        myfile<<evals-round(pers/100.0*gen_pop)+id+1<<"   "<<min_obj<<"     ";
                        for (int v=0;v<varN;v++){
                            myfile<<min_var[v]<<"  ";
                        }
                        myfile<<endl;
                        temp.push_back(min_obj);
                    }
                }
                else{
                    for (int id=0;id<gen_pop;id++){
                        if (db_temp[id].phi<minN){
                            minN = db_temp[id].phi;
                            min_pos=id;
                            min_obj=db_temp[id].obj[0];
                            min_var.clear();
                            for (int v=0;v<varN;v++){
                                min_var.push_back(db_temp[id].var[v]);
                            }
                        }
                        myfile<<(evals-gen_pop+id+1)<<"   "<<min_obj<<"     ";
                        for (int v=0;v<varN;v++){
                            myfile<<min_var[v]<<"  ";
                        }
                        myfile<<endl;

                        /*for (int v=0;v<varN;v++){
                            myfile<<db_temp[min_pos].var[v]<<"   ";
                        }
                        myfile<<endl;*/
                        temp.push_back(min_obj);
                    }

                }
                cout<<"minimum "<<minN<<endl;
                best_ind(db,db_temp,db_temp_s);
                db_temp.clear();

            }

        }
        myfile.close();
        objective2.push_back(temp);
        temp.clear();
    }

    //Briskw mesh timh kai 3s gia monokrithriakh beltistopoihsh
    if (objN==1){
        vector <double> mean;
        mean.resize(maxEval);
        vector <double> sdev;
        sdev.resize(maxEval);

        //Ypologismos meshs timhs
        for (int plo=0;plo<maxEval;plo++){
            mean[plo]=0;
            for (int gene=0;gene<run_times;gene++){
                mean[plo]+=objective2[gene][plo];
            }
            mean[plo]/=run_times;
        }

        stringstream ss;
        if (use_m==0){
            ss<<"avgM"<<run_times<<".txt";
        }
        else{
            ss<<"avgN"<<run_times<<".txt";
        }
        ofstream finfile(ss.str().c_str(),ios::app);

        //Ypologismos typikhs apoklishs kai perasma sto arxeio a
        for (int plo=0;plo<maxEval;plo++){
            sdev[plo]=0;
            for (int gene=0;gene<run_times;gene++){
                sdev[plo]+=(mean[plo]-objective2[gene][plo])*(mean[plo]-objective2[gene][plo]);
            }
            sdev[plo]=3*sqrt(sdev[plo]/(run_times-1));
            finfile<< plo+1 <<"   "<<mean[plo]<<"   "<<sdev[plo]<<"   "<<endl;
        }

        finfile.close();
    }


    Beep(1000,1000);
    return 0;
}

#include<iostream>
#include<cmath>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<map>
#include<ctime>
#include<new>
#include "CNF.h"
#ifdef _MPI
#include <mpi.h>
#endif
using namespace std;
int num_procs;
int proc_id;
int num_tst = 309318;//332410;
int num_trn = 309318;//2410;//154096;
int num_all_trn,num_all_tst;
int dimX = 93+1;
int num_gates = 10;
int num_params = 0;
int train_num=20000;
int numModel;
int mModelTrain;
int mModelPred;
double* mWeight;
double* mGrad;
double* mGradSum;
//double grad[num_params];
//double weights[num_params];
//double grad_sum[num_params];
double *grad;
double *weights;
double *grad_sum;
double **trnX=NULL;
double ** tstX=NULL;
double *trnY=NULL;
double *tstY=NULL;
string datafile;
//double trnX[num_trn][dimX], tstX[num_tst][dimX];
//double trnY[num_trn], tstY[num_tst];
double reg = 0.00001;
//double W[num_gates][dimX];
//double gW[num_gates][dimX];
//double V[num_gates+1];
//double gV[num_gates+1];
double**  W;
double** gW;
double* V;
double* gV;
map<string,string> params;
void SetSeed()
{
	unsigned int randomSeed=0;
	ifstream in("/dev/urandom",ios::in);
	in.read((char*)&randomSeed, sizeof(unsigned)/sizeof(char));
	in.close();
	unsigned id=time(NULL);
	randomSeed=randomSeed*randomSeed+id*id;
	srand48(randomSeed);
	//we can set the random seed at only the main function
	srand(randomSeed);
}

//------ modified to Bug-Fixed version -------//
double Gates(double x)
{
	double sum=x;
	if(sum>=0)
	{
		double value=exp(-1.0*sum);
		return (double)( 1.0/(1.0+value) );
	}
	else
	{
		double value=exp(1.0*sum);
		return (double) (1.0*value/(1.0+value) );
	}

}
void MakeFeature()
{
	//normalize
	int normFlag=0;
	vector<double> mean(dimX,0);
	vector<double> std(dimX,0);
	if(params["-norm"]!="null") {
		normFlag=1;
		ifstream fpNorm(params["-norm"].c_str());
		if(!fpNorm.is_open())
		{
			cerr<<"Normalization file opened error\n";
			exit(-1);
		}
		while(fpNorm.good()) {
			for(int i=0; i<dimX; i++)
				fpNorm>>mean[i];
			for(int j=0; j<dimX; j++)
				fpNorm>>std[j];
		}
		fpNorm.close();
	}
	if(proc_id==0 && normFlag==1)cerr<<"norm loaded\n";
	//-> mask init
	int *mask=new int[dimX-1];
	for(int fi=0; fi<dimX-1; fi++)mask[fi]=1;
	//-> load training
	for(int i=0; i<num_trn; i++) 
	{
		for(int j=dimX-2; j>=0; j--) 
		{
			if(normFlag)trnX[i][j]=(trnX[i][j]-mean[j+1])/std[j+1];
			if(mask[j]==0)trnX[i][j]=0;			
		}
	}
	if(proc_id==0)cerr<<"train make feature done\n";
	//-> load testing
	for(int i=0; i<num_tst; i++) 
	{
		for(int j=dimX-2; j>=0; j--) 
		{
			if(normFlag)tstX[i][j]=(tstX[i][j]-mean[j+1])/std[j+1];
			if(mask[j]==0)tstX[i][j]=0;
		}
	}
	if(proc_id==0)cerr<<"test make feature done\n";
}
void LoadData()//load mpi ok
{
	ifstream fp(datafile.c_str());
	if(!fp.is_open())cerr<<"open file error"<<datafile<<endl;
	vector<vector< double > > vData;
	int nRead=0;
	while(fp.good()) {
		vector<double> oneData(dimX,0);
		int i=0;
		for(i=0; i<dimX; i++) {
			if(fp.eof())break;
			fp>>oneData[i];
		}
		if(i==dimX && (nRead % num_procs) ==proc_id)
			vData.push_back(oneData);
		nRead++;
	}
	num_trn=vData.size();
	cerr<<"data loaded "<<num_trn<<endl;
//	num_tst=num_trn/10;
	num_tst=num_trn/4;
	if(params["-act"]=="predict")
		num_tst=num_trn;
	num_trn=num_trn-num_tst;

	//--- load training data ----//
	trnY=new double[num_trn];
	trnX=new double*[num_trn];
	for(int i=0; i<num_trn; i++) 
	{
		trnY[i]=vData[i][0];
		trnX[i]=new double[dimX];
		for(int j=0; j<dimX-1; j++) 
		{
			trnX[i][j]=vData[i][j+1];
		}
		trnX[i][dimX-1]=1;
	}
	//--- load testing data ----//
	tstY=new double[num_tst];
	tstX=new double*[num_tst];
	for(int i=num_trn; i<num_trn+num_tst; i++) 
	{
		tstY[i-num_trn]=vData[i][0];
		tstX[i-num_trn]=new double[dimX];
		for(int j=0; j<dimX-1; j++) 
		{
			tstX[i-num_trn][j]=vData[i][j+1];
		}
		tstX[i-num_trn][dimX-1]=1;
	}
	//screen out
	if(proc_id==0)cerr << trnY[0] << " " << trnY[1] << " " << trnY[2] << endl;
	if(proc_id==0)cerr << "Read Data" << endl;
	MakeFeature();
}
void LoadDataListMpi(string listFile)
{
	vector<string> allDataFile;
	ifstream ifDataList(listFile.c_str());
	int nData=0;
	while(ifDataList.good()) {
		string s;
		ifDataList>>s;
		if( (nData % num_procs) == proc_id && s.length()>0)
			allDataFile.push_back(s);
		nData++;
	}
	if(proc_id==0) {
		cerr<<allDataFile.size();
	}
	ifDataList.close();
	vector<vector<double> > vData;
	for(int i=0; i<allDataFile.size(); i++) {
		string datafile=allDataFile[i];
		ifstream cdata(datafile.c_str());
		if(!cdata.is_open()) {
			cerr<<"Open datafile err: "<<datafile<<"\n";
		}
//		if(proc_id==0)cerr<<"Loading "<<datafile<<",";
		while(cdata.good()) {
			vector<double> data(dimX,0);
			for(int j=0; j<dimX; j++) {
				string x;
				cdata>>x;
				data[j]=atof(x.c_str());
			}
			if(params["-Bin"]!="") {
				if(GetModelType(data)==atoi(params["-Bin"].c_str()) ) {
					vData.push_back(data);
				}
			} else {
				vData.push_back(data);
			}
			if(vData.size()==2 && proc_id==0) {
				cerr<<"data: ";
				cerr<<vData[0].size()<<" "<<vData[0][0]<<" "<<vData[0][1]<<endl;
				cerr<<vData[1].size()<<" "<<vData[1][0]<<" "<<vData[1][1]<<endl;
			}
		}
		cdata.close();
		cdata.clear();
	}
	if(proc_id==0)cerr<<"\nFinish loading step 1\n";
	//copy vector data to pointer
	num_trn=vData.size();
	num_tst=num_trn/10;
	num_trn=num_trn-num_tst;
	if(proc_id==0)
		cerr<<"num_trn "<<num_trn<<" num_tst "<<num_tst<<endl;
	trnX=new double* [num_trn];
	for(int i=0; i<num_trn; i++)
		trnX[i]=new double[dimX];
	tstX=new double* [num_tst];
	for(int i=0; i<num_tst; i++)
		tstX[i]=new double[dimX];
	trnY=new double[num_trn];
	tstY=new double[num_tst];
	for(int h=0; h<num_trn; h++) {
		trnY[h]=vData[h][0];
		for(int j=0; j<dimX-1; j++) {
			trnX[h][j]=vData[h][j+1];
		}
		trnX[h][dimX-1]=1;
	}
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	if(proc_id==0)	cerr<<"\nFinish loading step 2\n";
#endif
	for(int i=num_trn; i<num_tst+num_trn; i++) {
		tstY[i-num_trn]=vData[i][0];
		tstX[i-num_trn][dimX-1]=1;
		for(int j=0; j<dimX-1; j++) {
			tstX[i-num_trn][j]=vData[i][j+1];
		}
	}
	if(proc_id<2)
		cerr << "Read Data List mpi proc id "<<proc_id<<" "<<trnY[0] << " " << trnY[1] << " " << trnY[2] << " num "<<num_trn<<" "<<num_tst<<endl;
	//cerr << "Read Data List mpi version" << endl;
	MakeFeature();
	if(proc_id<2)
		cerr << "After makefeature Data List mpi proc id "<<proc_id<<" "<<trnY[0] << " " << trnY[1] << " " << trnY[2] << endl;
}
void LoadDataMpi()  //old ,not use it again
{
	ifstream cdata(datafile.c_str());
	int h=0;
	for(int i=0; i<num_trn; i++) {
		double tmp;
		cdata >> trnY[h];
		trnY[h]=abs(trnY[h]);
		for(int j=0; j<dimX-1; j++) {
			cdata >> trnX[h][j];
		}
		trnX[h][dimX-1]=1;
		if(i % num_procs ==proc_id)
			h++;
	}
	num_all_trn=num_trn-num_tst;
	num_trn=h;
	num_all_tst=num_tst;
	num_tst=num_tst/num_procs;
	for(int i=0; i<num_tst; i++) {
		tstY[i]=trnY[num_trn-1-i];
		tstX[i][dimX-1]=1;
		for(int j=0; j<dimX-1; j++) {
			tstX[i][j]=trnX[num_trn-1-i][j];
		}
	}
	num_trn=num_trn-num_tst;
	cdata.close();
	if(proc_id<2)
		cerr << "Read Data mpi proc id "<<proc_id<<" "<<trnY[0] << " " << trnY[1] << " " << trnY[2] << endl;
	//cerr << "Read Data mpi version" << endl;
	MakeFeature();
}
class _LBFGS: public Minimizer
{
	void ComputeGradient(vector<double> &g, const vector<double> &x);
	double ComputeFunction(const vector<double> &x);
	void Report (const string &s);
public:
	_LBFGS();
	void Report (const vector<double> &theta, int iteration, double objective, double step_length);
};
_LBFGS::_LBFGS()
	: Minimizer(false)
{
}
void _LBFGS::Report (const string &s)
{
	if(proc_id==0)cerr << s << endl;
}
void assignW(const vector<double> x)
{
	int pivot=0;
	for(int i=0; i<num_gates; i++)
		for(int j=0; j<dimX; j++)
			W[i][j]=x[pivot],pivot++;
	for(int i=0; i<=num_gates; i++)
		V[i]=x[pivot],pivot++;
}
void assignW2(double *x)
{
	int pivot=0;
	for(int i=0; i<num_gates; i++)
		for(int j=0; j<dimX; j++)
			W[i][j]=x[pivot],pivot++;
	for(int i=0; i<=num_gates; i++)
		V[i]=x[pivot],pivot++;
}
void assignG(vector<double> &g)
{
	int pivot=0;
	for(int i=0; i<num_gates; i++)
		for(int j=0; j<dimX; j++)
			g[pivot]=gW[i][j],pivot++;
	for(int i=0; i<=num_gates; i++)
		g[pivot]=gV[i],pivot++;
}
void _LBFGS::Report (const vector<double> &theta, int iteration, double objective, double step_length)
{
	double err1_sum = 0, err2_sum = 0, obj_sum=0;
	double err1 = 0, err2 = 0;
	double maxErr=0;
	double maxAllErr=0;
	//cerr << num_procs << " " << proc_id << endl;
	//objective is 1/2(square of difference of all data)
	if((iteration%50) && iteration > 10) return;
	for(int i=0; i<num_params; i++)
		weights[i] = theta[i];
	if(proc_id==0)cerr<<"Report "<<iteration<<endl;
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(weights, num_params, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
	assignW2(weights);
	double gates[num_gates];
	if(proc_id==0)
		cerr<<"weights summary "<<weights[0]<<" "<<weights[1]<<" "<<weights[2]<<endl;
	vector<double> predY(num_tst,0);
	if(proc_id==0)cerr<<"predY size is "<<num_tst<<" num of model "<<numModel<<endl;
	if(params["-act"]=="predict") {
		
		printf("#-> %d \n",num_tst);
		for(int i=0; i<num_tst; i++) 
		{
			assignW2(weights);
			double pred = V[num_gates];
			for(int j=0; j<num_gates; j++) 
			{
				double y=0;
				for(int k=0; k<dimX; k++)
					y+=W[j][k]*tstX[i][k];
				gates[j]=Gates(y);
				pred+=V[j]*gates[j];
			}
			predY[i]=pred;
			printf("%lf\n",pred);
		}
		err1 = 0, err2 = 0;
		ofstream fout(params["-o"].c_str());
		int nError=0;
		for(int i=0; i<predY.size(); i++) 
		{
			fout<<predY[i]<<endl;
			{
				nError++;
				double ee=(fabs(predY[i])-fabs(tstY[i]))*(fabs(predY[i])-fabs(tstY[i]));
				err1 = err1+ fabs(predY[i]-tstY[i]);
				err2 = err2+ ee;
				cerr<<ee<<" ";
			}
		}
		fout.close();
		cerr<<endl;
		cerr<<"pred "<<predY[num_tst-1]<<" "<<tstY[num_tst-1]<<endl;
		cerr<<"pred "<<predY[num_tst-2]<<" "<<tstY[num_tst-2]<<endl;

		err1_sum=err1;
		err2_sum=err2;
		num_all_tst = num_tst;
		err1_sum /= num_all_tst;
		err2_sum = err2_sum / num_all_tst;
		if(proc_id==0) {
			cerr << "test ABS err = " << err1_sum << endl;
			cerr << "[result] test RMS-err " << err2_sum << " on "<< nError<<" "<<num_all_tst<< endl;
		}
		return;//End of prediction
	}
	//if(proc_id==0)cerr<<"error ";
	for(int i=0; i<num_tst; i++) {
		double pred = V[num_gates];
		for(int j=0; j<num_gates; j++) {
			double y=0;
			for(int k=0; k<dimX; k++)
				y+=W[j][k]*tstX[i][k];
			gates[j]=Gates(y);
			pred+=V[j]*gates[j];
		}
		predY[i]=pred;
		err1 += fabs(pred-tstY[i]);
		double ee= (pred-tstY[i])*(pred-tstY[i]);
		//if(proc_id==0)cerr<<ee<<" ";
		err2+=ee;
		if(maxErr<ee)maxErr=ee;
	}
	//cerr<<endl;
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&err1,&err1_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&err2,&err2_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&num_tst,&num_all_tst, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&num_trn,&num_all_trn, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&maxErr,&maxAllErr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	//	MPI_Reduce(&objective,&obj_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#else
	err1_sum=err1;
	err2_sum=err2;
	num_all_tst=num_tst;   //mod by WS//__111010__//
	num_all_trn=num_trn;   //mod by WS//__111010__//
	maxAllErr=maxErr;      //mod by WS//__111010__//
#endif
//save model
	if(proc_id==0) {
		cerr<<"err2 "<<err2<<endl;
		cerr << err1 << " " << err1_sum << endl;
		cerr << err2 << " " << err2_sum << endl;
		cerr << num_tst << " " << num_trn << endl;
		cerr << num_all_tst << " " << num_all_trn << endl;
		err1_sum /= num_all_tst;
		err2_sum = err2_sum / num_all_tst;
		cout << "Iter " <<  iteration << endl;
		cout << "Objective "<< objective<<endl;
		cout << "Average training error "<<objective/num_all_trn<<endl;
		cout << "num_all_tst " <<  num_all_tst << endl;
		cout << " ABS err = " << err1_sum << endl;
		cout << " RMS err = " << err2_sum << endl;
		cout << " Max err = " << maxAllErr << endl;
		string model_file = "./model/model_";
		model_file=model_file+params["-modeltag"]+"_"+params["-r"]+"_"+params["-nGates"]+"_"+params["-mask"]+"_";
		char buf[100];
		sprintf(buf,"%d",iteration);
		model_file+=buf;
		ofstream fout(model_file.c_str());
		// fout << "num_params: " << num_params << endl;
		// fout << "num_gates: " << num_gates << endl;
		// fout << "dimX: " << dimX << endl;
		for(map<string,string>::iterator it=params.begin(); it!=params.end(); ++it) {
			if(it->first=="-act") continue;
			if(it->second=="")it->second="null";
			fout<< it->first << " "<<it->second<<endl;
		}
		fout<<"-iterStart "<<iteration <<endl;
		fout << " weights: " << endl;
		for(int i=0; i<num_params; i++)
			fout << weights[i] << " ";
		fout << endl;
		fout.close();
	}
}
void _LBFGS::ComputeGradient(vector<double>&g, const vector<double> &x)
{
	for(int i=0; i<num_params; i++)
		weights[i] = x[i];
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(weights, num_params, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
	assignW2(weights);
	for(int i=0; i<num_gates; i++)
		for(int j=0; j<dimX; j++)
			gW[i][j]=0;
	for(int i=0; i<num_gates+1; i++)
		gV[i]=0;
	double gates[num_gates];
	for(int i=0; i<num_trn; i++) {
		double pred = V[num_gates];
		for(int j=0; j<num_gates; j++) {
			double y=0;
			for(int k=0; k<dimX; k++)
				y+=W[j][k]*trnX[i][k];
			gates[j]=Gates(y);
			pred+=V[j]*gates[j];
		}
		double err = pred-trnY[i];
		gV[num_gates] += err ;
		for(int j=0; j<num_gates; j++) {
			gV[j] += err*gates[j] ;
			for(int k=0; k<dimX; k++)
				gW[j][k] += err*V[j]*gates[j]*(1.-gates[j])*trnX[i][k] ;
		}
	}
	assignG(g);
	for(int i=0; i<num_params; i++)
		grad[i] = g[i],grad_sum[i]=0;
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(grad, grad_sum, num_params, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(grad_sum, num_params, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for(int i=0; i<num_params; i++) {
		g[i] = grad_sum[i]+2*reg*weights[i];
	}
#else
	for(int i=0; i<num_params; i++)
		g[i]=grad[i]+2*reg*weights[i]; //mod by WS//__111010__//
#endif
}
double _LBFGS::ComputeFunction(const vector<double> &x)
{
	for(int i=0; i<num_params; i++)
		weights[i] = x[i];
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(weights, num_params, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
	assignW2(weights);
	double obj=0,obj_sum=0;
	for(int i=0; i<num_trn; i++) {
		double pred = V[num_gates];
		for(int j=0; j<num_gates; j++) {
			double y=0;
			for(int k=0; k<dimX; k++)
				y+=W[j][k]*trnX[i][k];
			pred+=V[j]*Gates(y);
		}
		obj += (trnY[i]-pred)*(trnY[i]-pred);
	}
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&obj, &obj_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&obj_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
	obj_sum = obj;
#endif
//	obj_sum=obj_sum/2;
	//  for(int i=0;i<num_params;i++)
	//   obj_sum=obj_sum+reg*weights[i]*weights[i];
	return obj_sum;
}
int GetModelType(double* feature)
{
	float dist=feature[0];
	int mtype;
	if(dist<7) {
		mtype=0;
	} else {
		if(dist<12)
			mtype=1;
		else
			mtype=2;
	}
	return mtype;
	//model 0, templ geo dist in [-,7]
	//model ,templ geo dist in [7,12]
	//model , templ geo dist in [12,-]
	//
}
int GetModelType(vector<double> feature)  //For read data, data[0] is actually the label
{
	float dist=fabs(feature[1]-feature[2]);
	int mtype;
	if(dist<30) {
		mtype=0;
	} else {
		if(dist<60)
			mtype=1;
		else
			mtype=2;
	}
	return mtype;
	//model 0, templ geo dist in [-,7]
	//model ,templ geo dist in [7,12]
	//model , templ geo dist in [12,-]
	//
}
void Train()
{
	weights=NULL;
	dimX=atoi(params["-dim"].c_str());
	numModel=atoi(params["-mModel"].c_str());
	num_gates=atoi(params["-nGates"].c_str());
	datafile=params["-iFile"];
	reg=atof(params["-r"].c_str());
	if(params["-mModel"].length()>0)LoadModel();
	num_params = num_gates*dimX + num_gates+1;
	grad_sum=new double[num_params];
	W=new double*[num_gates];
	for(int i=0; i<num_gates; i++)
		W[i]=new double[dimX];
	gW=new double*[num_gates];
	for(int i=0; i<num_gates; i++)
		gW[i]=new double[dimX];
	V=new double[num_gates+1];
	gV=new double[num_gates+1];
	_LBFGS* lbfgs = new _LBFGS;
	vector<double> w_init(num_params,0);
	if((params["-mModel"].length()==0) || (weights==NULL) ) {
		weights=new double[num_params];//after load parameters and models
		for(int i=0; i<num_params; i++) {
			weights[i]=w_init[i] = 1.0*drand48()/10.0 - 0.05;
		}
	}
	else
	{
	  for(int i=0; i<num_params; i++) {
	    w_init[i]=weights[i];
	  }
	}
	//multi model
	if(proc_id==0) {
		cerr<<"params "<<num_params<<" "<<datafile<<" "<<endl;
		for(map<string,string>::iterator it=params.begin(); it!=params.end(); ++it) {
			cerr<< it->first << " "<<it->second<<endl;
		}
	}
	grad=new double[num_params];
	if(proc_id==0) {
		cerr<<"begin one model training\n";
		cerr<<"weights ... "<<w_init[0]<<" "<<w_init[1]<<endl;
	}
#ifdef _MPI
	if(params["-iFileList"]!="" && params["-iFileList"]!="null")
	LoadDataListMpi(params["-iFileList"]);	
	else if(params["-iFile"]!="")
		LoadData();	
		else
			cerr<<"No data loaded\n";
#else
	LoadData();
#endif
	lbfgs->LBFGS(w_init,train_num,atoi(params["-iterStart"].c_str()));
}
void Predict()
{
	datafile=params["-iFile"];
	LoadModel();
	cerr << "predict" << endl;
	num_tst=atoi(params["-nTst"].c_str());
	num_trn=atoi(params["-nTrn"].c_str());
	dimX=atoi(params["-dim"].c_str());
	num_gates=atoi(params["-nGates"].c_str());
	reg=atof(params["-r"].c_str());
	train_num=atoi(params["-nIter"].c_str());
	num_params = num_gates*dimX + num_gates+1;
	W=new double*[num_gates];
	for(int i=0; i<num_gates; i++)
		W[i]=new double[dimX];
	gW=new double*[num_gates];
	for(int i=0; i<num_gates; i++)
		gW[i]=new double[dimX];
	V=new double[num_gates+1];
	gV=new double[num_gates+1];
	//only predict
	if(proc_id==0) {
		cerr<<"params "<<num_params<<" "<<datafile<<" "<<endl;
		for(map<string,string>::iterator it=params.begin(); it!=params.end(); ++it) {
			cerr<< it->first << " "<<it->second<<endl;
		}
	}
	LoadData();
	vector<double> w_init(num_params*numModel,0);
	for(int i=0; i<num_params*numModel; i++)
		w_init[i]=weights[i];
	_LBFGS* lbfgs = new _LBFGS;
	lbfgs->Report(w_init,0,0,0);
}
void LoadModel()
{
	vector<string> mModelFile=StringSplit(params["-model"],":");
	if(proc_id==0)
		for(int i=0;i<mModelFile.size();i++){
			cerr<<"Model "<<i<<" : "<<mModelFile[i]<<endl;
		}
		
	numModel=mModelFile.size();//always > 0 ;check other places
	if(mModelFile.size()>=1) {
		weights=NULL;
		for(int mi=0; mi<mModelFile.size(); mi++) {
			string s;
			ifstream fin(mModelFile[mi].c_str(),ios::in);
			while(fin.good()) {
				//override params
				fin>>s;
				if(s=="weights:")break;
				string b;
				fin>>b;
				if(s=="-nTrn")continue;
				if(s=="-nTst")continue;
				if(s=="-iFile")continue;
				params[s]=b;
			}
			dimX=atoi(params["-dim"].c_str());//dimX and -dim should be true input dim +1
			num_gates=atoi(params["-nGates"].c_str());
			num_params = num_gates*dimX + num_gates+1;
			if(weights==NULL)weights=new double[num_params*mModelFile.size()];
			if(s=="weights:") {
				//continue to load all the weights.
				for(int i=0; i<num_params; i++) {
					fin>>weights[i+mi*num_params];
				}
			}
		}
	} else {
		cerr<<"no model is given\n";
	}
}
int main(int argc,char **argv)
{
	//initial,parse the argc and argv.
	num_procs=1;
	proc_id=0;
#ifdef _MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	if(proc_id==0)cerr<<"Starting mpi nnreg...\n";
#endif
	for(int i=1; i<argc; i=i+2)  {
		string a,b;
		a.assign(argv[i]);
		b.assign(argv[i+1]);
		params[a]=b;
	}
	if(params["-act"]!="predict") {
		Train();
	} else {
		Predict();
	}
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
	return 0;
}
/*
	string s;
	vector<string> mModelFile;
	check mModelPred
*/

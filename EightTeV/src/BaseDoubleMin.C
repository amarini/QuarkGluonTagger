#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "TROOT.h"
#include "TCanvas.h"
#include "THStack.h"


using namespace std;

double sqr(double a) {return a*a;}
#define EPSILON 0.001
bool CompareFloats (float A, float B) 
{
   float diff = A - B;
   return (diff < EPSILON) && (-diff < EPSILON);
}

class TRIO {
public:
	TRIO(){};
	~TRIO(){};
	TRIO(int id,float x,float w){pdgId=id;value=x;weight=w;};
	int pdgId;
	float value;
	float weight;
};
const inline bool operator<(TRIO&x,TRIO&y){return x.value<y.value;}

class BaseAnalyzer{
public:
	BaseAnalyzer(){
		PtMin=30;PtMax=80;RhoMin=0;RhoMax=15;EtaMin=0;EtaMax=2.0;
		t_mc=NULL;t_data=NULL;
		varName="QGLHisto";
		nBins=30;xMin=0;xMax=1.000001;
		nstep=10;stp0=.01;stp1=0.01;
		lmin=0;lmax=1.0;
		opt="CHI2 WW";
		aMin=0.5;aMax=1.3;bMin=0.5;bMax=1.5;
		iteration=0;
		readStartPoints=0;
		filenameStartPoints="";
		compress=0;
		}
	string varName;//QGL HISTO
	int iteration; //internal counter
	int nstep;
	float stp0;
	float stp1;
	int readStartPoints;
	string filenameStartPoints;
	int compress;

	float Chi2(TH1F*h1,TH1F*h2);	
	float LogL(float bw=1.0/30.);
	int ResetFast();	
	float function(float x0, float a ,float b,float min=0,float max=1);
	void ComputeMin();
	virtual void Loop(TChain *t,int type); //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc
	void FillHisto(TH1F *h,string var);
	pair<float, float> MinG(TGraph2D *g,double *min0=NULL,double*min1=NULL); //min0 minimum, min1=minlast point - 0-1
	void SpanMin();
	void CreateHisto(int type=3);
	void SetTrees(TChain *mc,TChain*data){t_data=data;t_mc=mc;}
	void Reset(TH1F *h);
	pair<float,float> SmearDoubleMin(float a0_q,float b0_q , float a0_g,float b0_g,int type); //type = 0 Q, 1 G
	void ComputeDoubleMin();
	void ComputeMinFast();

	pair<float,float> SmearDoubleMinFast(float a0_q,float b0_q , float a0_g,float b0_g,int type,int WriteOut=0); //type = 0 Q, 1 G
	void ComputeDoubleMinFast();
	void LoopFast();
	virtual void LoadBins();
	void ResetBins();
	TCanvas *DrawValidation();
	virtual bool ExtraCuts(){return true;}
//private:
	TChain *t_mc;
	TChain *t_data;
	float PtMin;
	float PtMax;
	float RhoMin;
	float RhoMax;
	float EtaMin;
	float EtaMax;
	
	float lmin;
	float lmax;
	float a_q,a_g,b_q,b_g;
	float alpha,beta;

	float aMin;
	float aMax;
	float bMin;
	float bMax;
	
	vector<float> alphaFast;
	vector<float> betaFast;
	vector<TH1F*> h_mcFast;

	Int_t nBins;
	string opt;
	Double_t xMin,xMax;
//SPAN
	vector<pair<float,float> > PtBins;
	vector<pair<float,float> > EtaBins;
	vector<pair<float,float> > RhoBins;
	
//---- internal not modify
	TH1F*h_mc,*h_data;
	TGraph2D *g2;
	
	map<string,float> treeVar;
	map<string,int> treeVarInt;
//----- vector with all the likelihood /MLP results - for a given selection /pdgId/value
	//vector<pair<int,float> > varAll;
	vector<TRIO> varAll;
	vector<float> varAllData;
	void CompressVector(int N=1000);
};


void BaseAnalyzer::CompressVector(int N){
vector<TRIO> *varNew=new vector<TRIO>;
float max=varAll[0].value;
float min=varAll[0].value;
//need to check pdgid too
for(int z=0;z<varAll.size();z++)
	{	
	if(varAll[z].value>max)max=varAll[z].value;
	if(varAll[z].value<min)min=varAll[z].value;
	}
//bin divide
max+=0.00001;
//(max-min)/N+min;
for(int i=0;i<N;i++)
	{
	float x=(max-min)/N*(i+0.5) + min ;
	varNew->push_back(TRIO( 1, x ,0) );
	varNew->push_back(TRIO( 21, x ,0) );
	}
float dx=(max-min)/N;
for(int z=0;z<varAll.size();z++)
	{
	float x=varAll[z].value;
	float w=varAll[z].weight;
	float i=varAll[z].pdgId;
	if( abs(i)<5 && abs(i)!=0 ) i=1;
	else i=21;//compress including everything else in gluon
	int n=floor( (x-min)/dx );
	if(i==1){(*varNew)[2*n].weight+=w;	
		if( fabs((*varNew)[2*n].value - varAll[z].value )> 0.01 ) printf("ERROR\n");
		}
	else { (*varNew)[2*n+1].weight+=w;	
		if( fabs((*varNew)[2*n+1].value - varAll[z].value )> 0.01 ) printf("ERROR\n");}
	}
varAll.clear();
varAll=*varNew;
}


float BaseAnalyzer::function(float x0, float a ,float b,float min,float max)
{
using namespace TMath;
float x=(x0-min)/(max-min); 
if(x<=0)x=0;
if(x>=1)x=1;
float x1= (TanH( a* ATanH(2*x-1)+b )/2+.5 ) ;
if(x<=0)x1=0; //prevent overflow and underflow bins
if(x>=1)x1=1;
return x1*(max-min)+min;
}

float BaseAnalyzer::Chi2(TH1F*h1,TH1F*h2)
{
	for(int i=1;i<=h1->GetNbinsX();i++)
		{
		if(h1->GetBinContent(i)==0)h1->SetBinError(i,1);
		if(h2->GetBinContent(i)==0)h2->SetBinError(i,1);
		}
	return h1->Chi2Test(h2,opt.c_str());
}
float BaseAnalyzer::LogL(float bw)
{
	double L=1;
	for(int i=0;i<int(varAllData.size());++i)
	{
	double densityEstimator=0;
	double Nj= varAll.size();
	for(int j=0;j<int(varAll.size());++j)
		{
		//weights, smearing
		float a=1,b=0;
		if(abs(varAll[j].pdgId)<5 && varAll[j].pdgId!=0) {a=a_q;b=b_q;}
		else // also the unmatched
			{a=a_g;b=b_g;}
		float value=function(varAll[j].value,a,b,lmin,lmax);
		densityEstimator+=varAll[j].weight *  TMath::Exp(- sqr(value-varAllData[i])/sqr(bw) );
		}
	L*=densityEstimator;
	}
	L=-TMath::Log(L); //- sign so -> minimize as chi2
	return float(L);
	
}

void BaseAnalyzer::Loop(TChain *t,int type){ //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc

		
		if(type&4) {lmin=1.0;lmax=0;} //reset lmin-lmax
		if(type&1) {delete h_data; CreateHisto(1);}
		if(type&10) {delete h_mc; CreateHisto(2);} //8+2
		if(type&32) {varAll.clear();} //reset varAll

}

void BaseAnalyzer::FillHisto(TH1F *h,string var){
	h->Fill(  function(treeVar[var],alpha,beta,lmin,lmax ), treeVar["eventWeight"]*treeVar["PUReWeight"] );
	}

void BaseAnalyzer::CreateHisto(int type)
	{	
	if(type&2){
	h_mc=new TH1F("hmc","hmc",nBins,xMin,xMax);
	h_mc->Sumw2();
	}
	if(type&1){
	h_data=new TH1F("hdata","hdata",nBins,xMin,xMax);
	h_data->Sumw2();
	}
	}
void BaseAnalyzer::Reset(TH1F *h)
	{
	for(int i=0;i<=h->GetNbinsX()+1;i++) {h->SetBinContent(i,0);h->SetBinError(i,0);}
	h->SetEntries(0);
	//h->Sumw2();
	}

pair<float, float> BaseAnalyzer::MinG(TGraph2D *g,double *min0,double*min1){
	pair<float,float> R(-99,-99);
	if(g==NULL) return R;
	if(g->GetN()==0)return R;
	double *x1,*y1,*z1;
	x1=g->GetX();
	y1=g->GetY();
	z1=g->GetZ();
	
	float a=z1[0];R=pair<float,float>(x1[0],y1[0]);
	for(int i=0;i<g->GetN();i++){if((z1[i]<a)||(a<0)){a=z1[i];    R=pair<float,float>(x1[i],y1[i]); }}
	
	if(min0!=NULL) *min0=a;
	if(min1!=NULL) *min1=z1[g2->GetN()-1]	;
	return R;
}

void BaseAnalyzer::ComputeMin(){
	g2=new TGraph2D(); //TODO
	
	alpha=1.0;beta=0;
	Loop(t_data,1);
	if(varName=="QGLMLP")
		Loop(t_mc,4);
	//scan
	alpha=1.0;beta=0;
	for(float ai=0.7; ai<=1.1; ai+=0.02)
		{
		Reset(h_mc);	
		alpha=ai;
		Loop(t_mc,2);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		g2->SetPoint(g2->GetN(),alpha,beta, Chi2(h_data,h_mc)  );	
		}
	alpha=1.0;beta=0;
	for(float bi=-0.5; bi<=0.5; bi+=0.01)
		{
		Reset(h_mc);	
		beta=bi;
		Loop(t_mc,2);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		g2->SetPoint(g2->GetN(),alpha,beta, Chi2(h_data,h_mc)  );	
		}
	
	//Find min0;min1
	float min0=1,min1=0;
	pair<float,float> R=MinG(g2);
	min0=R.first;min1=R.second;

	for(int i=-nstep;i<=nstep;i++)
	for(int j=-nstep;j<=nstep;j++)
        	{
        	alpha=min0+i*stp0;
        	beta=min1+j*stp1;
		Loop(t_mc,2);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		g2->SetPoint(g2->GetN(),alpha,beta, Chi2(h_data,h_mc)  );
		}
	//double min0,min1;	
	R=MinG(g2);
	printf("a=%.3f;b=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R.first,R.second,lmin,lmax);
	return;
	
}

void BaseAnalyzer::LoadBins(){
		PtBins.push_back(  pair<float,float>(30,50) );
		PtBins.push_back(  pair<float,float>(50,80) );
		PtBins.push_back(  pair<float,float>(80,120) );
		PtBins.push_back(  pair<float,float>(120,250) );
		
		RhoBins.push_back(  pair<float,float>(0,15) );
		RhoBins.push_back(  pair<float,float>(15,40) );
		
		EtaBins.push_back(  pair<float,float>(0,2) );
		EtaBins.push_back(  pair<float,float>(3,4.7) );
	return;
	}
void BaseAnalyzer::ResetBins(){
	PtBins.clear();
	EtaBins.clear();
	RhoBins.clear();
	return;
	}

void BaseAnalyzer::SpanMin(){

		ResetBins();
		LoadBins();

	for ( int e=0; e< int(EtaBins.size());e++)
	for ( int p=0; p< int(PtBins.size()) ;p++)
	for ( int r=0; r< int(RhoBins.size());r++)
		{
			//save temporary values
			float stp0_=stp0;
			float stp1_=stp1;
			int nstep_=nstep;
			if(e>0) {nstep*=2; fprintf(stderr,"Double*Double stps at high eta\n");};
		fprintf(stderr,"Bins: %d %d %d\n",e,p,r);
		PtMin=PtBins[p].first;PtMax=PtBins[p].second;
		EtaMin=EtaBins[e].first;EtaMax=EtaBins[e].second;
		RhoMin=RhoBins[r].first;RhoMax=RhoBins[r].second;
	
		//if(g2!=NULL)delete g2;	
		int t=0;
		if(varName=="QGLMLP")t=2;
		if(varName=="QGLHisto")t=3;
		int bin=(p+1)+(r+1)*10+(e+1)*100 + t*1000;
		//OLD printf("//%s: Pt=%.0f_%.0f Rho=%.0f_%.0f Eta=%.0f_%.0f\n",varName.c_str(),PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax);
		//OLD printf("case %d:",bin);
		printf("%s %.0f %.0f %.0f %.0f %.1f %.1f ",varName.c_str(),PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax);
		//ComputeMinFast(); //Be Fast!
		//ComputeDoubleMin(); //Be Slow!
		ComputeDoubleMinFast(); //Be Fast!
			//restore
			stp0=stp0_;
			stp1=stp1_;
			nstep=nstep_;
		}
	fprintf(stderr,"DONE\n");	
	}
void BaseAnalyzer::ComputeDoubleMin(){
	fprintf(stderr,"Compute DoubleMin\n");
	nstep=5; //otherwise too slow?
	pair<float,float> R_q,R_g;
	fprintf(stderr,"First Smear\n");
	R_q=SmearDoubleMin(1,0,1,0,0); //
	R_g=SmearDoubleMin(R_q.first,R_q.second,1,0,1); //
	R_q=SmearDoubleMin(R_q.first,R_q.second,R_g.first,R_g.second,0); //
	R_g=SmearDoubleMin(R_q.first,R_q.second,R_g.first,R_g.second,1); //

	printf("a_q=%.3f;b_q=%.3f;a_g=%.3f;b_g=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R_q.first,R_q.second,R_g.first,R_g.second,lmin,lmax);
	}
pair<float,float> BaseAnalyzer::SmearDoubleMin(float a0_q,float b0_q , float a0_g,float b0_g,int type){ //type = 0 Q, 1 G
	fprintf(stderr,"SmearDoubleMin\n");
	TGraph2D *g2_q=new TGraph2D(); 
	TGraph2D *g2_g=new TGraph2D(); 
	
	alpha=1.0;beta=0;

	//if(h_data!=NULL)delete h_data;
	//if(h_mc!=NULL)delete h_mc;
	fprintf(stderr,"Creating Histos\n");
	CreateHisto(3);
	
	fprintf(stderr,"Going to do Data Loop\n");
	Loop(t_data,1);

	if(varName=="QGLMLP")
		Loop(t_mc,4);
	//scan
	a_q=a0_q;b_q=b0_q;
	a_g=a0_g;b_g=b0_g;

	alpha=1.0;beta=0;
	fprintf(stderr,"Going to do span ai\n");
	for(float ai=0.7; ai<=1.1; ai+=0.02)
		{
		Reset(h_mc);	
		if(type==0)a_q=ai;
		if(type==1)a_g=ai;
		Loop(t_mc,8);
		for(int j=0;j<=h_mc->GetNbinsX()+1;j++)if(h_mc->GetBinError(j)==0)h_mc->SetBinError(j,1);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());

		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, Chi2(h_data,h_mc)  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, Chi2(h_data,h_mc)  );	
		}
	alpha=1.0;beta=0;
	a_q=a0_q;b_q=b0_q;
	a_g=a0_g;b_g=b0_g;

	fprintf(stderr,"Going to do span bi\n");
	for(float bi=-0.5; bi<=0.5; bi+=0.01)
		{
		Reset(h_mc);	
		if(type==0)b_q=bi;
		if(type==1)b_g=bi;
		Loop(t_mc,8);
		for(int j=0;j<=h_mc->GetNbinsX()+1;j++)if(h_mc->GetBinError(j)==0)h_mc->SetBinError(j,1);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, Chi2(h_data,h_mc)  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, Chi2(h_data,h_mc)  );	
		}
	
	//Find min0;min1
	float min0=1,min1=0;
	pair<float,float> R;
		if(type==0)R=MinG(g2_q);
		if(type==1)R=MinG(g2_g);
	min0=R.first;min1=R.second;

	for(int i=-nstep;i<=nstep;i++)
	for(int j=-nstep;j<=nstep;j++)
        	{
        	if(type==0)a_q=min0+i*stp0;
        	if(type==1)a_g=min0+i*stp0;
        	if(type==0)b_q=min1+j*stp1;
        	if(type==1)b_g=min1+j*stp1;
		Loop(t_mc,8);
		for(int k=0;k<=h_mc->GetNbinsX()+1;k++)if(h_mc->GetBinError(k)==0)h_mc->SetBinError(k,1);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, Chi2(h_data,h_mc)  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, Chi2(h_data,h_mc)  );	
		}
	
		if(type==0)R=MinG(g2_q);
		if(type==1)R=MinG(g2_g);
	//SAME ON G,& REDO
	//printf("a=%.3f;b=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R.first,R.second,lmin,lmax);
	return R;
}

int BaseAnalyzer::ResetFast()
{
	alphaFast.clear();
	betaFast.clear();
	for(int i=0;i<int(h_mcFast.size());i++){delete h_mcFast[i]; };
	h_mcFast.clear();
	return 0;

}
void BaseAnalyzer::ComputeMinFast(){
	g2=new TGraph2D(); 
	
	alpha=1.0;beta=0;
	Loop(t_data,1);
	for(int j=0;j<=h_data->GetNbinsX()+1;j++)if(h_data->GetBinError(j)==0)h_data->SetBinError(j,1);
	if(varName=="QGLMLP")
		Loop(t_mc,4);
	//scan
	//reset Fast
	ResetFast();
	
	alpha=1.0;beta=0;
	for(float ai=aMin; ai<=aMax; ai+=0.02)
		{
		Reset(h_mc);	
		alphaFast.push_back(ai);	
		betaFast.push_back(0);	
		h_mcFast.push_back(new TH1F( Form("hmc_%d",int(h_mcFast.size())),"hmc",nBins,xMin,xMax) );
		}
	alpha=1.0;beta=0;
	for(float bi=bMin; bi<=bMax; bi+=0.01)
		{
		Reset(h_mc);	
		alphaFast.push_back(1.0);	
		betaFast.push_back(bi);	
		h_mcFast.push_back(new TH1F( Form("hmc_%d",int(h_mcFast.size())),"hmc",nBins,xMin,xMax) );
		}

	for(int j=0;j< int(h_mcFast.size());j++) h_mcFast[j]->Sumw2();	
	Loop(t_mc,16);
		for(int i=0 ;i<int(alphaFast.size());i++)
			{
			for(int j=0;j<=h_mcFast[i]->GetNbinsX()+1;j++)if(h_mcFast[i]->GetBinError(j)==0)h_mcFast[i]->SetBinError(j,1);
			h_mcFast[i]->Scale(h_data->Integral()/h_mcFast[i]->Integral());
			g2->SetPoint(g2->GetN(),alphaFast[i],betaFast[i], Chi2(h_data,h_mcFast[i])  );

			}
	//Find min0;min1
	float min0=1,min1=0;
	pair<float,float> R=MinG(g2);
	min0=R.first;min1=R.second;
		
	ResetFast();

	delete g2;
	g2=new TGraph2D();

	for(int i=-nstep;i<=nstep;i++)
	for(int j=-nstep;j<=nstep;j++)
        	{
        	alphaFast.push_back(min0+i*stp0 );
        	betaFast.push_back(min1+j*stp1 );
		h_mcFast.push_back(new TH1F( Form("hmc_%d",int(h_mcFast.size())),"hmc",nBins,xMin,xMax) );
		//g2->SetPoint(g2->GetN(),alpha,beta, Chi2(h_data,h_mc)  );
		}
	alphaFast.push_back(min0);
	betaFast.push_back(0);
	h_mcFast.push_back(new TH1F( Form("hmc_%d",int(h_mcFast.size())),"hmc",nBins,xMin,xMax) );

	alphaFast.push_back(min1);
	betaFast.push_back(bMax);
	h_mcFast.push_back(new TH1F( Form("hmc_%d",int(h_mcFast.size())),"hmc",nBins,xMin,xMax) );
	
	alphaFast.push_back(1);
	betaFast.push_back(0);
	h_mcFast.push_back(new TH1F( Form("hmc_%d",int(h_mcFast.size())),"hmc",nBins,xMin,xMax) );
	
	for(int j=0;j<int(h_mcFast.size());j++) h_mcFast[j]->Sumw2();	
	Loop(t_mc,16);
		for(int i=0 ;i<int(alphaFast.size());i++)
			{
			for(int j=0;j<=h_mcFast[i]->GetNbinsX()+1;j++)if(h_mcFast[i]->GetBinError(j)==0)h_mcFast[i]->SetBinError(j,1);
			h_mcFast[i]->Scale(h_data->Integral()/h_mcFast[i]->Integral());
			g2->SetPoint(g2->GetN(),alphaFast[i],betaFast[i], Chi2(h_data,h_mcFast[i])  );
			}
	double m0,m1;
	R=MinG(g2,&m0,&m1);
	//OLD printf("a=%.3f;b=%.3f;lmin=%.3f;lmax=%.3f;break;//chi2=%.3lf; chi2_0=%.3lf\n",R.first,R.second,lmin,lmax,m0,m1);
	printf(" %.3f %.3f %.3f %.3f \n#chi2=%.3lf; chi2_0=%.3lf\n",R.first,R.second,lmin,lmax,m0,m1);
	{
	TFile *out=TFile::Open("output.root","UPDATE");out->cd();
	for(int i=0;i<int(h_mcFast.size());i++)
		{
		h_mcFast[i]->SetName(Form("%s_alpha%.2f_beta%.2f_lmin%.3f_lmax%.3f_pt%0f_%.0f_rho%.0f_%.0f_eta%.0f_%.0f",varName.c_str(),alphaFast[i],betaFast[i],lmin,lmax,PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax));
		h_mcFast[i]->Write();
		}
	h_data->SetName(Form("%s_data_pt%0f_%.0f_rho%.0f_%.0f_eta%.0f_%.0f",varName.c_str(),PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax));
	h_data->Write();
	}
	return;
}
//-----------------------------------DOUBLE MIN FAST ---------------------------------------------------------
void BaseAnalyzer::ComputeDoubleMinFast(){
	//nstep=5; //otherwise too slow?
	alpha=1.0;beta=0;
	iteration=0;
	Loop(t_data,1);
	if(varName=="QGLMLP")
		Loop(t_mc,4);
	Loop(t_mc,32);
	if(compress) CompressVector(1000);
	pair<float,float> R_q,R_g;
	if(!readStartPoints){
	fprintf(stderr,"Not Read Start points\n");
	R_g=SmearDoubleMinFast(1,0,1,0,1); //
	R_q=SmearDoubleMinFast(1,0,R_g.first,R_g.second,0); //
	}else{
	fprintf(stderr,"Read Start points\n");
	//read the points
	FILE *fr=fopen(filenameStartPoints.c_str(),"r");
	char _var_[1023];
	float _ptmin_,_ptmax_,_rhomin_,_rhomax_,_etamin_,_etamax_,_aq_,_bq_,_ag_,_bg_,_lmin_,_lmax_;
	float _aq_mean_=0,_bq_mean_=0,_ag_mean_=0,_bg_mean_=0;int _n_=0;
	while (fscanf(fr,"%s %f %f %f %f %f %f %f %f %f %f %f %f",_var_,&_ptmin_,&_ptmax_,&_rhomin_,&_rhomax_,&_etamin_,&_etamax_,&_aq_,&_bq_,&_ag_,&_bg_,&_lmin_,&_lmax_) != EOF) {
	if(varName==string(_var_) && (CompareFloats(_etamin_,EtaMin)) && ( CompareFloats(_etamax_,EtaMax)))
		{
		//cross comparing of the limit values. Adiacent bins
		if( CompareFloats(_ptmin_,PtMax) || CompareFloats(_ptmax_,PtMin) )	
			{
			_aq_mean_+=_aq_;
			_ag_mean_+=_ag_;
			_bq_mean_+=_bq_;
			_bg_mean_+=_bg_;
			_n_++; //can be 1 or 2
			}
		}
	}
	fclose(fr);
	R_q=pair<float,float>(_aq_mean_/_n_,_bq_mean_/_n_);
	R_g=pair<float,float>(_ag_mean_/_n_,_bg_mean_/_n_);
	R_g=SmearDoubleMinFast(R_q.first,R_q.second,R_g.first,R_g.second,1); //
	R_q=SmearDoubleMinFast(R_q.first,R_g.second,R_g.first,R_g.second,0); //
	}
	fprintf(stderr,"Common swap\n");
	R_g=SmearDoubleMinFast(R_q.first,R_q.second,R_g.first,R_g.second,1); //
	R_q=SmearDoubleMinFast(R_q.first,R_q.second,R_g.first,R_g.second,0,1); //

	printf("%.3f %.3f %.3f %.3f %.3f %.3f \n",R_q.first,R_q.second,R_g.first,R_g.second,lmin,lmax);
	fprintf(stderr,"PRINTED\n");
	}

void BaseAnalyzer::LoopFast()
{
//delete and recreate h_mc (empty)
delete  h_mc;
CreateHisto(2); 
fprintf(stderr,"LoopFast\n");
for(int z=0;z<int(varAll.size());++z){ 
	
	{alpha=1;beta=0;}	
	if(varAll[z].pdgId == 21){alpha=a_g; beta=b_g;} 
	if(fabs(varAll[z].pdgId) <5 ) {alpha=a_q;beta = b_q;}
//	if( fabs(treeVarInt["pdgIdPartJet0"])== 0) {alpha=1;beta=0;}
	if(fabs(varAll[z].pdgId) == 0) {alpha=a_g;beta=b_g;}
	if(fabs(varAll[z].pdgId) == 999) {alpha=a_g;beta=b_g;}

	treeVar[varName]=varAll[z].value;
	if(varAll[z].weight>=0) {treeVar["eventWeight"]=varAll[z].weight;treeVar["PUReweight"]=1;}
	FillHisto(h_mc,varName);
	}

}

TCanvas *BaseAnalyzer::DrawValidation()
{

TCanvas *c=new TCanvas(Form("validation"),"",800,1000);
TPad *up=new TPad("upPad","upPad",0,0.3,1,1);
TPad *dn=new TPad("dnPad","dnPad",0,0,1,0.3);

up->Draw();
dn->Draw();

up->cd();
	TH1F * h_q=new TH1F("h_q","h_q",nBins,xMin,xMax); h_q->Sumw2();
	TH1F * h_g=new TH1F("h_g","h_g",nBins,xMin,xMax); h_g->Sumw2();
	TH1F * h_b=new TH1F("h_b","h_b",nBins,xMin,xMax); h_b->Sumw2();
	TH1F * h_u=new TH1F("h_u","h_u",nBins,xMin,xMax); h_u->Sumw2();

alpha=1.0;
beta=0.0;
//Loop(t_mc,32);	
for(int z=0;z<int(varAll.size());++z){ 
	switch( abs(varAll[z].pdgId ) )
		{
		case 21:
			h_g->Fill(varAll[z].value,varAll[z].weight);break;
		case 1: case 2: case 3: case 4:
			h_q->Fill(varAll[z].value,varAll[z].weight);break;
		case 5:
			h_b->Fill(varAll[z].value,varAll[z].weight);break;
		case 0: default:
			h_u->Fill(varAll[z].value,varAll[z].weight);break;
		}
	}
//Loop(t_data,1);	
	//Normalization	
	{
	float norm=h_data->Integral();
	float q=h_q->Integral();
	float g=h_g->Integral();
	float b=h_b->Integral();
	float u=h_u->Integral();
	float mc=q+g+b+u;
	
	h_q->Scale( norm/mc );
	h_g->Scale( norm/mc );
	h_b->Scale( norm/mc );
	h_u->Scale( norm/mc );
	}
	//Style
	h_q->SetLineColor(kBlack); h_q->SetFillColor(kBlue-4);
	h_g->SetLineColor(kBlack); h_g->SetFillColor(kRed-4);
	h_b->SetLineColor(kBlack); h_b->SetFillColor(kGreen-4);
	h_u->SetLineColor(kBlack); h_u->SetFillColor(kGray+1);
	h_data->SetMarkerStyle(20);h_data->SetMarkerColor(kBlack);h_data->SetLineColor(kBlack);

	THStack *S=new THStack("stack","stack");
		S->Add(h_u);
		S->Add(h_b);
		S->Add(h_g);
		S->Add(h_q);
	S->Draw("HIST");
	S->Draw("AXIS X+ Y+ SAME");
	//h_data->Draw("AXIS");
	//h_data->Draw("AXIS X+ Y+ SAME");
	h_data->Draw("P SAME");
dn->cd();
	TH1F* r=(TH1F*)h_data->Clone("ratio");
	TH1F* mc=(TH1F*)h_q->Clone("mc");
	mc->Add(h_g);mc->Add(h_b);mc->Add(h_u);
	r->Divide(mc);
		delete mc;
	r->Draw("P");
return c;
}

pair<float,float> BaseAnalyzer::SmearDoubleMinFast(float a0_q,float b0_q , float a0_g,float b0_g,int type,int WriteOut){ //type = 0 Q, 1 G

	TGraph2D *g2_q=new TGraph2D(); g2_q->SetName("g2_q");
	TGraph2D *g2_g=new TGraph2D(); g2_g->SetName("g2_g");
	//scan
	a_q=a0_q;b_q=b0_q;
	a_g=a0_g;b_g=b0_g;

	alpha=1.0;beta=0;
	fprintf(stderr,"ai\n");
	for(float ai=0.7; ai<=1.1; ai+=0.02)
		{
		Reset(h_mc);	
		if(type==0)a_q=ai;
		if(type==1)a_g=ai;
		LoopFast();	
		h_mc->Scale(h_data->Integral()/h_mc->Integral());

		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, Chi2(h_data,h_mc)  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, Chi2(h_data,h_mc)  );	
		//if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, LogL(1./nBins)  );	
		//if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, LogL(1./nBins)  );	
		}
	alpha=1.0;beta=0;
	a_q=a0_q;b_q=b0_q;
	a_g=a0_g;b_g=b0_g;

	fprintf(stderr,"bi\n");
	for(float bi=-0.5; bi<=0.5; bi+=0.01)
		{
		Reset(h_mc);	
		if(type==0)b_q=bi;
		if(type==1)b_g=bi;
		LoopFast();
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, Chi2(h_data,h_mc)  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, Chi2(h_data,h_mc)  );	
		//if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, LogL(1./nBins)  );	
		//if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, LogL(1./nBins)  );	
		}
	
	//Find min0;min1
	float min0=1,min1=0;
	pair<float,float> R;
		if(type==0)R=MinG(g2_q);
		if(type==1)R=MinG(g2_g);
	min0=R.first;min1=R.second;

	fprintf(stderr,"double loop\n");
	for(int i=-nstep;i<=nstep;i++)
	for(int j=-nstep;j<=nstep;j++)
        	{
        	if(type==0)a_q=min0+i*stp0;
        	if(type==1)a_g=min0+i*stp0;
        	if(type==0)b_q=min1+j*stp1;
        	if(type==1)b_g=min1+j*stp1;
		LoopFast();
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		if(type==0){g2_q->SetPoint(g2_q->GetN(),a_q,b_q, Chi2(h_data,h_mc)  ); }
		if(type==1){g2_g->SetPoint(g2_g->GetN(),a_g,b_g, Chi2(h_data,h_mc)  );	}
		//if(type==0){g2_q->SetPoint(g2_q->GetN(),a_q,b_q, LogL(1./nBins)  ); }
		//if(type==1){g2_g->SetPoint(g2_g->GetN(),a_g,b_g, LogL(1./nBins)  );	}
		}
	
		if(type==0)R=MinG(g2_q);
		if(type==1)R=MinG(g2_g);
	
	
	if(WriteOut){	
	string name=Form("Results/output_%s_pt%.0f_%.0f_rho%.0f_%.0f_eta%.0f_%.0f",varName.c_str(),PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax);
	if(type==0)g2_q->SaveAs((name+"g2_q.root").c_str());
	if(type==1)g2_g->SaveAs((name+"g2_g.root").c_str());
	}
	
	iteration++;	
	return R;
}
//----------------------------------------------------------------------------------------------------




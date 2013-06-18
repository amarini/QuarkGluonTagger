

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
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
#include "TStyle.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TCanvas.h"

//---compute QGL & QGLMPL on the fly
#include "/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/src/QGMLPCalculator.cc"
#include "/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/src/parameters.cc"
#include "/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/src/QGLikelihoodCalculator.cc"
#include "/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/src/Bins.cc"

using namespace std;

double deltaPhi(double phi1,double phi2){
	double result= phi1-phi2;
	while (result> M_PI) result -=2*M_PI;
	while (result < -M_PI) result+=2*M_PI;
	return result;
	}
double deltaR(double eta1,double phi1, double eta2, double phi2){
	double deta=eta1-eta2;
	double dphi=deltaPhi(phi1,phi2);
	return sqrt(deta*deta+ dphi*dphi);
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

class Analyzer{
public:
	Analyzer(){
		PtMin=40;PtMax=60;RhoMin=0;RhoMax=15;EtaMin=0;EtaMax=2.0;
		t_mc=NULL;t_data=NULL;
		varName="QGLHisto";
		nBins=30;xMin=0;xMax=1.000001;
		nstep=10;stp0=.01;stp1=0.01;
		lmin=0;lmax=1.0;
		opt="CHI2 WW";
		aMin=0.5;aMax=1.3;bMin=0.5;bMax=1.5;
		TFile *Fpuw=TFile::Open("/afs/cern.ch/work/s/sunil/public/forTom/PU_rewt_flatP6.root");
		TFile *Fptetaw=TFile::Open("/afs/cern.ch/work/s/sunil/public/forTom/Jetpteta_rewt2D_flatP6.root ");
		puw=(TH1F*)Fpuw->Get("hist_WT")->Clone("hPU_wt");
		ptetaw=(TH2F*)Fptetaw->Get("hist_WT")->Clone("hPtEta_wt");
		qgl=new QGLikelihoodCalculator("/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/data/");//ReducedHisto_2012.root");
		qgmlp=new QGMLPCalculator("MLP","/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/data/",true); //prob
		}
	string varName;//QGL HISTO
	int nstep;
	float stp0;
	float stp1;
	
	int ResetFast();	
	float function(float x0, float a ,float b,float min=0,float max=1);
	void ComputeMin();
	void Loop(TChain *t,int type); //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc
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
//private:
	//histo reweight
		TH1F* puw;
		TH2F* ptetaw;
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
	
//---- internal not modify
	TH1F*h_mc,*h_data;
	TGraph2D *g2;
	
	map<string,float> treeVar;
	map<string,int> treeVarInt;
//----- vector with all the likelihood /MLP results - for a given selection /pdgId/value
	//vector<pair<int,float> > varAll;
	vector<TRIO> varAll;
//--- QGL QGMLP
	QGLikelihoodCalculator *qgl;
	QGMLPCalculator *qgmlp;
};


float Analyzer::function(float x0, float a ,float b,float min,float max)
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

void Analyzer::Loop(TChain *t,int type){ //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc

		
		treeVarInt["nvtx"] = -999; t->SetBranchAddress("nvtx",&treeVarInt["nvtx"]);
		treeVar["rho"] = -999; t->SetBranchAddress("rho",&treeVar["rho"]);
		Float_t jetPt[4]; t->SetBranchAddress("jetPt",&jetPt);
		Float_t jetEnergy[4]; t->SetBranchAddress("jetEnergy",&jetEnergy);
		Float_t jetBtag[4]; t->SetBranchAddress("jetBtag",&jetBtag);
		Float_t jetBeta[4]; t->SetBranchAddress("jetBeta",&jetBeta);
		Float_t jetEta[4]; t->SetBranchAddress("jetEta",&jetEta);
		Float_t jetPhi[4]; t->SetBranchAddress("jetPhi",&jetPhi);
		Float_t jetAxis_QC[2][4]; t->SetBranchAddress("jetAxis_QC",&jetAxis_QC);
		Float_t jetAxis[2][4]; t->SetBranchAddress("jetAxis",&jetAxis);
		Float_t jetPtD[4]; t->SetBranchAddress("jetPtD",&jetPtD);
		Float_t jetPtD_QC[4]; t->SetBranchAddress("jetPtD_QC",&jetPtD_QC);
		Int_t jetChgPart_ptcut[4]; t->SetBranchAddress("jetChgPart_ptcut",&jetChgPart_ptcut);
		Int_t jetChgPart_QC[4]; t->SetBranchAddress("jetChgPart_QC",&jetChgPart_QC);
		Int_t jetNeutralPart_ptcut[4]; t->SetBranchAddress("jetNeutralPart_ptcut",&jetNeutralPart_ptcut);
		vector<int> *partonId=0;if(type >1)t->SetBranchAddress("partonId",&partonId);
		vector<int> *partonSt=0;if(type >1)t->SetBranchAddress("partonSt",&partonSt);
		vector<float> *partonPt=0;if(type >1)t->SetBranchAddress("partonPt",&partonPt);
		vector<float> *partonEta=0;if(type >1)t->SetBranchAddress("partonEta",&partonEta);
		vector<float> *partonPhi=0;if(type >1)t->SetBranchAddress("partonPhi",&partonPhi);
		vector<float> *partonE=0;if(type >1)t->SetBranchAddress("partonE",&partonE);
		
		vector<bool> *triggerResult=0;if(type ==1)t->SetBranchAddress("triggerResult",&triggerResult);
		
		
		if(type&4) {lmin=1.0;lmax=0;} //reset lmin-lmax
		if(type&1) {delete h_data; CreateHisto(1);}
		if(type&10) {delete h_mc; CreateHisto(2);} //8+2
		if(type&32) {varAll.clear();} //reset varAll

		for(int i=0;i<t->GetEntries() ;i++)
			{
			t->GetEntry(i);
			treeVar["ptJet0"]=jetPt[0];
			treeVar["etaJet0"]=jetEta[0];
			treeVar["rhoPF"]=treeVar["rho"];
			
			//fprintf(stderr,"A: Pt: %f<%f<%f - Eta: %f<%f<%f: Rho: %f<%f<%f\n",PtMin,treeVar["ptJet0"],PtMax,EtaMin,treeVar["etaJet0"],EtaMax,RhoMin,treeVar["rhoPF"],RhoMax);
			if((treeVar["ptJet0"]<PtMin)||(treeVar["ptJet0"]>PtMax)||(fabs(treeVar["etaJet0"])<EtaMin)||(fabs(treeVar["etaJet0"])>EtaMax)|| (treeVar["rhoPF"]<RhoMin)||(treeVar["rhoPF"]>RhoMax))continue;
			//fprintf(stderr,"-B\n");
			//selection
			double muJet_dphi=deltaPhi(jetPhi[0],jetPhi[1]);
			if(fabs(muJet_dphi)<2.5) continue;
			//fprintf(stderr,"--C\n");
			if( ! (2.0 *jetPt[2]/ (jetPt[0]+jetPt[1])<.3) )continue; 
			//fprintf(stderr,"---D\n");
			if( jetBtag[0] >0.244)continue;
			//fprintf(stderr,"----E\n");
			//trigger --only on data
			if( type==1 && !( triggerResult != NULL && triggerResult->size()>1 && triggerResult->at(1) )) continue;
			
			//parton Matching
			double dR_min=999;
			int pos_min=999;
			//int part_min=5;
			
			//fprintf(stderr,"_______not NULL: %ld = %ld = %ld\n",partonPt,partonEta,partonPhi);
			if(type>1){ //only on MC
			for(int iPart=0;iPart<partonPt->size();iPart++)
				{
				double dR_ipart= deltaR(partonEta->at(iPart),partonPhi->at(iPart),jetEta[0],jetPhi[0]);
				if(dR_ipart< dR_min){dR_min=dR_ipart;pos_min=iPart;}
				}
			}
			if(dR_min<.3){
				//fprintf(stderr,"_______%f<%f\n",pos_min,partonId->size());
				treeVarInt["pdgIdPartJet0"]=partonId->at(pos_min);
				} else treeVarInt["pdgIdPartJet0"]=0;
			
			//fprintf(stderr,"_______E2:pos_min=%d dR=%f\n",pos_min,dR_min);
			map<TString,float> variables_MLP;	
			map<TString,float> variables_corr_MLP;	
			//map<TString,float> variables_QGL;	
			
//			variables_QGL["axis1"]= jetAxis_QC[0][0];
//			variables_QGL["axis2"]= jetAxis_QC[1][0];
//			variables_QGL["ptd"] = jetPtD_QC[0];
//			variables_QGL["mult"] = jetChgPart_QC[0]+jetNeutralPart_ptcut[0];
//			variables_QGL["pt"] = jetPt[0];
//			variables_QGL["eta"] = jetEta[0];
//			variables_QGL["rho"] = rho;
	
			//fprintf(stderr,"_______E3:nvtx: %d\n",treeVarInt["nvtx"]);
			if(fabs(jetEta[0])<2.5 && (jetBeta[0]<(1.0 -0.2*TMath::Log(treeVarInt["nvtx"]-0.67)))) continue;	
			//fprintf(stderr,"-----F\n");
			float sub_data=0.0;
			if(fabs(jetEta[0])>2.5 && type==1)sub_data=1.0;
			//Variables as general variables
			treeVar["mult"]=float(jetChgPart_QC[0]+jetNeutralPart_ptcut[0])-sub_data;
			treeVar["axis1"]=jetAxis_QC[0][0];
			treeVar["axis2"]=jetAxis_QC[1][0];
			treeVar["ptD"]=jetPtD_QC[0];
			//Discriminators - only if needed -  save time
			if(varName=="QGLHisto")treeVar["QGLHisto"] = qgl->computeQGLikelihood2012(jetPt[0],jetEta[0],treeVar["rho"],jetChgPart_QC[0]+jetNeutralPart_ptcut[0]-sub_data,jetPtD_QC[0],jetAxis_QC[1][0]);
			if(varName=="QGLMLP"){	
				variables_MLP["axis1"]=jetAxis_QC[0][0];
				variables_MLP["axis2"]=jetAxis_QC[1][0];
				variables_MLP["ptD"]=jetPtD_QC[0];
				variables_MLP["mult"]=jetChgPart_QC[0];
				
				variables_MLP["pt"]=jetPt[0];
				variables_MLP["eta"]=jetEta[0];
				variables_MLP["rho"]=treeVar["rho"];
				
				if(fabs(jetEta[0])>2.5){
					variables_MLP["axis1"]=jetAxis[0][0];
					variables_MLP["axis2"]=jetAxis[1][0];
					variables_MLP["ptD"]=jetPtD[0];
					variables_MLP["mult"]=jetChgPart_QC[0]+jetNeutralPart_ptcut[0]-sub_data;
					
					}
				
				variables_corr_MLP["axis1"] = variables_MLP["axis1"];
				variables_corr_MLP["axis2"] = variables_MLP["axis2"];
				variables_corr_MLP["ptD"] = variables_MLP["ptD"];
				variables_corr_MLP["mult"] = variables_MLP["mult"];
			
				//variables_corr_MLP=qgmlp->TEST(variables_MLP,variables_corr_MLP);
				
				treeVar["QGLMLP"]=qgmlp->QGvalue(variables_MLP);
			}
			//treeVar["pdgIdPartJet0"];
			//---------------------------
		
			//fprintf(stderr,"------G\n");
			if(type&1){
				//printf("passed selection - type 1 --> %.3f - %.3f\n",treeVar[varName],treeVar[varName+"Fwd"]);
				string var=varName;
			//	if(EtaMin>2.5)var+="Fwd"; //only in data fwd
				alpha=1; beta=0;
					//int bin=puw->FindBin(treeVar["rho"]);
					//int bin2=ptetaw->FindBin(jetPt[0],fabs(jetEta[0]) );
					//float weight=puw->GetBinContent(bin) *  ptetaw->GetBinContent(bin2);
					treeVar["eventWeight"]=1.; //data
					treeVar["PUReWeight"]=1;
				FillHisto(h_data,var);
				}	
			if(type&2){
				//mc
			//fprintf(stderr,"________notNull:%ld %ld %f\n",puw,ptetaw,treeVar["rho"]);
					int bin = puw->FindBin(treeVar["rho"]);
					int bin2 = ptetaw->FindBin(jetPt[0],fabs(jetEta[0]) );
					float weight=puw->GetBinContent(bin) *  ptetaw->GetBinContent(bin2);
			//fprintf(stderr,"________bin:%d,%d, w=%f\n",bin,bin2,weight);
					treeVar["eventWeight"]=weight;treeVar["PUReWeight"]=1;
				FillHisto(h_mc,varName);
				}
			if(type&4){
				string var=varName;
				if(lmin>treeVar[var]) lmin=treeVar[var];
				if(lmax<treeVar[var]) lmax=treeVar[var];
				}
			if(type&8){
				alpha=1;beta=0;
				if( treeVarInt["pdgIdPartJet0"] ==21) {alpha=a_g;beta=b_g;}
				if( fabs(treeVarInt["pdgIdPartJet0"]) < 5) {alpha=a_q;beta=b_q;}
				if( fabs(treeVarInt["pdgIdPartJet0"])== 0) {alpha=1;beta=0;}

				FillHisto(h_mc,varName);
				}
			if(type&16){
				for(int j=0;j<int(alphaFast.size());j++)
					{
					alpha=alphaFast[j];
					beta=betaFast[j];
					FillHisto(h_mcFast[j],varName);
					}
				}
			if(type&32){ 
					int bin=puw->FindBin(treeVar["rho"]);
					int bin2=ptetaw->FindBin(jetPt[0],fabs(jetEta[0]) );
					float weight=puw->GetBinContent(bin) *  ptetaw->GetBinContent(bin2);
					varAll.push_back(TRIO(treeVarInt["pdgIdPartJet0"],treeVar[varName],weight)); //w=-1 default
					//work-around to make it smear all in one shot
					//varAll.push_back(1,treeVar[varName],weight); //w=-1 default
				}
			}
}

void Analyzer::FillHisto(TH1F *h,string var){
	//printf("Filling Histos with %f\n",function(treeVar[var],alpha,beta,lmin,lmax ));
	h->Fill(  function(treeVar[var],alpha,beta,lmin,lmax ), treeVar["eventWeight"]*treeVar["PUReWeight"] );
	}

void Analyzer::CreateHisto(int type)
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
void Analyzer::Reset(TH1F *h)
	{
	for(int i=0;i<=h->GetNbinsX()+1;i++) {h->SetBinContent(i,0);h->SetBinError(i,0);}
	h->SetEntries(0);
	//h->Sumw2();
	}

pair<float, float> Analyzer::MinG(TGraph2D *g,double *min0,double*min1){
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

void Analyzer::ComputeMin(){
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
		g2->SetPoint(g2->GetN(),alpha,beta, h_data->Chi2Test(h_mc,opt.c_str())  );	
		}
	alpha=1.0;beta=0;
	for(float bi=-0.5; bi<=0.5; bi+=0.01)
		{
		Reset(h_mc);	
		beta=bi;
		Loop(t_mc,2);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		g2->SetPoint(g2->GetN(),alpha,beta, h_data->Chi2Test(h_mc,opt.c_str())  );	
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
		g2->SetPoint(g2->GetN(),alpha,beta, h_data->Chi2Test(h_mc,opt.c_str())  );
		}
	//double min0,min1;	
	R=MinG(g2);
	printf("a=%.3f;b=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R.first,R.second,lmin,lmax);
	return;
	
}

void Analyzer::SpanMin(){
	vector<pair<float,float> > PtBins;
	vector<pair<float,float> > EtaBins;
	vector<pair<float,float> > RhoBins;
	
		//PtBins.push_back(  pair<float,float>(40,60) );
		//PtBins.push_back(  pair<float,float>(60,80) );
		PtBins.push_back(  pair<float,float>(80,120) );
		PtBins.push_back(  pair<float,float>(120,250) );
		PtBins.push_back(  pair<float,float>(250,500) );
		
		RhoBins.push_back(  pair<float,float>(0,15) );
		RhoBins.push_back(  pair<float,float>(15,40) );
		
		EtaBins.push_back(  pair<float,float>(0,2) );
		EtaBins.push_back(  pair<float,float>(3,4.7) );
	
	for ( int e=0; e< int(EtaBins.size());e++)
	for ( int p=0; p< int(PtBins.size()) ;p++)
	for ( int r=0; r< int(RhoBins.size());r++)
		{
		fprintf(stderr,"Bins: %d %d %d\n",e,p,r);
		PtMin=PtBins[p].first;PtMax=PtBins[p].second;
		EtaMin=EtaBins[e].first;EtaMax=EtaBins[e].second;
		RhoMin=RhoBins[r].first;RhoMax=RhoBins[r].second;
	
		//if(g2!=NULL)delete g2;	
		int t=0;
		if(varName=="QGLMLP")t=2;
		if(varName=="QGLHisto")t=3;
		int bin=(p+1)+(r+1)*10+(e+1)*100 + t*1000;
		printf("//%s: Pt=%.0f_%.0f Rho=%.0f_%.0f Eta=%.0f_%.0f\n",varName.c_str(),PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax);
		printf("case %d:",bin);
		ComputeMinFast(); //Be Fast!
		//ComputeDoubleMin(); //Be Slow!
		//ComputeDoubleMinFast(); //Be Fast!
		}
	printf("DONE\n");	
	}
void Analyzer::ComputeDoubleMin(){
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
pair<float,float> Analyzer::SmearDoubleMin(float a0_q,float b0_q , float a0_g,float b0_g,int type){ //type = 0 Q, 1 G
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

		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, h_data->Chi2Test(h_mc,opt.c_str())  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, h_data->Chi2Test(h_mc,opt.c_str())  );	
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
		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, h_data->Chi2Test(h_mc,opt.c_str())  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, h_data->Chi2Test(h_mc,opt.c_str())  );	
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
		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, h_data->Chi2Test(h_mc,opt.c_str())  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, h_data->Chi2Test(h_mc,opt.c_str())  );	
		}
	
		if(type==0)R=MinG(g2_q);
		if(type==1)R=MinG(g2_g);
	//SAME ON G,& REDO
	//printf("a=%.3f;b=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R.first,R.second,lmin,lmax);
	return R;
}

int Analyzer::ResetFast()
{
	alphaFast.clear();
	betaFast.clear();
	for(int i=0;i<int(h_mcFast.size());i++){delete h_mcFast[i]; };
	h_mcFast.clear();
	return 0;

}
void Analyzer::ComputeMinFast(){
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
			g2->SetPoint(g2->GetN(),alphaFast[i],betaFast[i], h_data->Chi2Test(h_mcFast[i],opt.c_str())  );
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
		//g2->SetPoint(g2->GetN(),alpha,beta, h_data->Chi2Test(h_mc,opt.c_str())  );
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
			g2->SetPoint(g2->GetN(),alphaFast[i],betaFast[i], h_data->Chi2Test(h_mcFast[i],opt.c_str())  );
			}
	double m0,m1;
	R=MinG(g2,&m0,&m1);
	printf("a=%.3f;b=%.3f;lmin=%.3f;lmax=%.3f;break;//chi2=%.3lf; chi2_0=%.3lf\n",R.first,R.second,lmin,lmax,m0,m1);
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
void Analyzer::ComputeDoubleMinFast(){
	//nstep=5; //otherwise too slow?
	alpha=1.0;beta=0;
	Loop(t_data,1);
	if(varName=="QGLMLP")
		Loop(t_mc,4);
	Loop(t_mc,32);

	pair<float,float> R_q,R_g;
	R_g=SmearDoubleMinFast(1,0,1,0,1); //
	R_q=SmearDoubleMinFast(R_q.first,R_q.second,1,0,0); //
	R_g=SmearDoubleMinFast(R_q.first,R_q.second,R_g.first,R_g.second,1,1); //
	R_q=SmearDoubleMinFast(R_q.first,R_q.second,R_g.first,R_g.second,0,1); //

	printf("a_q=%.3f;b_q=%.3f;a_g=%.3f;b_g=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R_q.first,R_q.second,R_g.first,R_g.second,lmin,lmax);
	}

void Analyzer::LoopFast()
{
//delete and recreate h_mc (empty)
delete  h_mc;
CreateHisto(2); 
for(int z=0;z<int(varAll.size());++z){ 
	{alpha=1;beta=0;}
	if(varAll[z].pdgId == 21){alpha=a_g; beta=b_g; } 
	if(fabs(varAll[z].pdgId) <5 ) {alpha=a_q;beta = b_q;}if( fabs(treeVarInt["pdgIdPartJet0"])== 0) {alpha=1;beta=0;}
	if(fabs(varAll[z].pdgId) == 0) {alpha=1;beta=0;}

	treeVar[varName]=varAll[z].value;
	if(varAll[z].weight>=0) {treeVar["eventWeight"]=varAll[z].weight;treeVar["PUReweight"]=1;}
	FillHisto(h_mc,varName);
	}

}

pair<float,float> Analyzer::SmearDoubleMinFast(float a0_q,float b0_q , float a0_g,float b0_g,int type,int WriteOut){ //type = 0 Q, 1 G

	TGraph2D *g2_q=new TGraph2D(); g2_q->SetName("g2_q");
	TGraph2D *g2_g=new TGraph2D(); g2_g->SetName("g2_g");
	
	//scan
	a_q=a0_q;b_q=b0_q;
	a_g=a0_g;b_g=b0_g;

	alpha=1.0;beta=0;
	for(float ai=0.7; ai<=1.1; ai+=0.02)
		{
		Reset(h_mc);	
		if(type==0)a_q=ai;
		if(type==1)a_g=ai;
		LoopFast();	
		h_mc->Scale(h_data->Integral()/h_mc->Integral());

		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, h_data->Chi2Test(h_mc,opt.c_str())  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, h_data->Chi2Test(h_mc,opt.c_str())  );	
		}
	alpha=1.0;beta=0;
	a_q=a0_q;b_q=b0_q;
	a_g=a0_g;b_g=b0_g;

	for(float bi=-0.5; bi<=0.5; bi+=0.01)
		{
		Reset(h_mc);	
		if(type==0)b_q=bi;
		if(type==1)b_g=bi;
		LoopFast();
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, h_data->Chi2Test(h_mc,opt.c_str())  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, h_data->Chi2Test(h_mc,opt.c_str())  );	
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
		LoopFast();
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, h_data->Chi2Test(h_mc,opt.c_str())  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, h_data->Chi2Test(h_mc,opt.c_str())  );	
		}
	
		if(type==0)R=MinG(g2_q);
		if(type==1)R=MinG(g2_g);
	//SAME ON G,& REDO
	//printf("a=%.3f;b=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R.first,R.second,lmin,lmax);
	
	if(WriteOut){	
	string name=Form("Results/outputZJet2_%s_pt%.0f_%.0f_rho%.0f_%.0f_eta%.0f_%.0f",varName.c_str(),PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax);
	if(type==0)g2_q->SaveAs((name+"g2_q.root").c_str());
	if(type==1)g2_g->SaveAs((name+"g2_g.root").c_str());
	}
	
	return R;
}
//----------------------------------------------------------------------------------------------------



int ComputeDoubleMinDiJet(){
	system("[ -f output.root ] && rm output.root");
	TChain *mc=new TChain("Hbb/events");
	TChain *data=new TChain("Hbb/events");
		mc->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_flatQCD_P6_Dijets.root");

		data  ->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_data2012_JetMon_Dijets.root");
	Analyzer A;
	A.nstep=20;
	A.varName="QGLHisto";
//	A.varName="QGLMLP";
	A.CreateHisto();
	A.SetTrees(mc,data);
		freopen("/dev/null","w",stderr);
	//A.ComputeMinFast(); //A.ComputeDoubleMin;
	//A.ComputeDoubleMin();
		/*
		A.alpha=1;
		A.beta=0;
		A.Loop(mc,2)
		A.Loop(data,1)
		*/
	fprintf(stderr,"Going to do Span\n");
	A.SpanMin();
	A.varName="QGLMLP";
	A.SpanMin();
	return 0;
	}
/*
.L ComputeDoubleMinDiJet.C+
TChain *mc=new TChain("Hbb/events");
TChain *data=new TChain("Hbb/events");
mc->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_flatQCD_P6_Dijets.root");
data  ->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_data2012_JetMon_Dijets.root");

Analyzer A;
A.nstep=20;
A.varName="QGLHisto";
//A.varName="mult";A.nBins=50;A.xMin=0;A.xMax=50;A.lmin=0;A.lmax=100;
//A.varName="axis1";A.nBins=50;A.xMin=0;A.xMax=.5;
//A.varName="axis2";A.nBins=50;A.xMin=0;A.xMax=.5;
//A.varName="ptD";A.nBins=50;A.xMin=0;A.xMax=1.0001;
A.varName="rho";A.nBins=50;A.xMin=0;A.xMax=50;A.lmin=0;A.lmax=100;
A.RhoMin=0; A.RhoMax=100;A.PtMin=40;A.PtMax=60; A.EtaMin=0;A.EtaMax=2.0;
A.CreateHisto();
A.SetTrees(mc,data);
A.alpha=1;
A.beta=0;
A.Loop(mc,2);
A.Loop(data,1);
A.h_data->SetMarkerStyle(20);
A.h_data->DrawNormalized("P");
A.h_mc->DrawNormalized("HIST SAME");

*/

Analyzer *Check(float PtMin,float PtMax,float RhoMin,float RhoMax, float EtaMin,float EtaMax, float alpha, float beta , const char * varName="QGLHisto",float lmin=0,float lmax=1,TCanvas **pC=NULL){
Analyzer *A=new Analyzer();
TChain *mc=new TChain("Hbb/events");
TChain *data=new TChain("Hbb/events");
mc->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_flatQCD_P6_Dijets.root");
data  ->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_data2012_JetMon_Dijets.root");
A->nstep=20; A->varName=varName;
A->RhoMin=RhoMin; A->RhoMax=RhoMax;A->PtMin=PtMin;A->PtMax=PtMax; A->EtaMin=EtaMin;A->EtaMax=EtaMax;
A->CreateHisto();
A->SetTrees(mc,data);
A->alpha=1;
A->beta=0;
A->Loop(mc,2);
A->Loop(data,1);
TH1F* h_mc0=(TH1F*)A->h_mc->Clone("h_mc0");h_mc0->SetLineColor(kGreen+2);
A->lmin=lmin;
A->lmax=lmax;
A->alpha=alpha;
A->beta=beta;
A->Loop(mc,2);
//DoubleMin
//A->Loop(mc,32);
//A->a_q=0.91;A->b_q=0.15;A->a_g=1.02;A->b_g=0.13;
//A->LoopFast();

gStyle->SetOptStat(0);
TCanvas *c=new TCanvas("c","c",800,800);
A->h_mc->SetLineColor(kBlue);
	A->h_mc->GetXaxis()->SetTitle(varName);
	A->h_mc->SetMinimum(0);
TH1F * hmc1=(TH1F*)A->h_mc->DrawNormalized("HIST");
	hmc1->GetYaxis()->SetRangeUser(0,0.20);

h_mc0->SetLineWidth(2); h_mc0->SetLineStyle(2);
h_mc0->DrawNormalized("HIST SAME");
A->h_data->SetMarkerStyle(20);
A->h_data->DrawNormalized("P SAME");
	TLegend *L=new TLegend(0.4,.7,.6,.89);
	L->AddEntry(A->h_data,"data");
	L->AddEntry(h_mc0,"mc");
	L->AddEntry(A->h_mc,"mc+syst");
	L->Draw();
if(pC!=NULL) (*pC)=c;

return A;
}

